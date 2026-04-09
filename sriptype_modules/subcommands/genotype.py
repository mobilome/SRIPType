"""Genotype analysis subcommand.

RIP-seq genotyping pipeline (17 steps per sample):
  1. Initial BLAST alignment
  2. Extract unique matches
  3. Convert to BED format
  4. Extract 70bp flanking sequences
  5. Filter short sequences (>30bp)
  6. Get FASTA sequences
  7. Trim sequences (remove first 7bp)
  8. SINE head alignment
  9. Insertion site extraction
 10. Non-insertion site filtering
 11. Extract 100bp flanking sequences
 12. Filter short sequences (>50bp)
 13. Get FASTA sequences
 14. Trim sequences (remove first 20bp)
 15. RIP region alignment
 16. Deletion site extraction
 17. Genotype generation
"""

import os
import shutil
import subprocess
import time
import traceback
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

from sriptype_modules.utils import setup_logger

logger = setup_logger(__name__)

DESCRIPTION = "RIP-seq genotyping analysis pipeline"
ORDER = 2

# Bundled data directory (query fasta + reference databases)
_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

# External tools required by this subcommand
REQUIRED_TOOLS = ["blastn", "bedtools"]


def add_arguments(parser):
    """Add subcommand-specific arguments to the parser."""
    # Input
    parser.add_argument(
        "-i", "--mkdb-dir",
        required=True,
        help="Directory produced by 'sriptype mkdb' (contains FASTA, BLAST DB, and stats files)",
    )

    # Output
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Main output directory for results",
    )

    # Parallelism
    parser.add_argument(
        "-j", "--jobs",
        type=int,
        default=1,
        help="Number of samples to process in parallel (default: 1)",
    )





def _check_tools(tools):
    """Check that required external tools are available."""
    missing = []
    for tool in tools:
        if subprocess.run(
            ["which", tool],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).returncode != 0:
            missing.append(tool)
    if missing:
        raise RuntimeError(
            f"Required tools not found in PATH: {', '.join(missing)}. "
            "Please install them before running this subcommand."
        )


# ---------------------------------------------------------------------------
# Internal analysis functions (replace external scripts)
# ---------------------------------------------------------------------------

def _convert_blast_to_bed(input_path, output_path):
    """Convert BLAST outfmt-6 to BED format, keeping best hit per read.

    Replaces 5-30-converbed.awk.
    """
    max_score = {}
    best_len = {}
    bed_data = {}

    with open(input_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 12:
                continue

            read = parts[1]
            sstart = int(parts[8])
            send = int(parts[9])
            score = float(parts[11])
            aln_len = abs(send - sstart)

            if (read not in max_score
                    or score > max_score[read]
                    or (score == max_score[read] and aln_len > best_len[read])):
                max_score[read] = score
                best_len[read] = aln_len

                if sstart < send:
                    bed_start = sstart
                    bed_end = send + 1
                    strand = "+"
                else:
                    bed_start = send - 1
                    bed_end = sstart
                    strand = "-"

                bed_data[read] = f"{read}\t{bed_start}\t{bed_end}\t.\t.\t{strand}"

    with open(output_path, "w") as fh:
        for read in bed_data:
            fh.write(bed_data[read] + "\n")


def _extract_sites(file1_path, file2_path, output_path):
    """Filter file2 by read names found in file1 (colon-split col1).

    Replaces 5-30-extract_sites-fan.py.
    """
    valid_reads = set()
    with open(file1_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if not parts:
                continue
            read_name = parts[0].split(":", 1)[0]
            valid_reads.add(read_name)

    with open(file2_path) as f_in, open(output_path, "w") as f_out:
        for line in f_in:
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split("\t")
            if len(parts) < 2:
                continue
            if parts[1].strip() in valid_reads:
                f_out.write(stripped + "\n")


def _generate_deletions(file1_path, file2_path, output_path):
    """Cross-validate RIP/read relationships and generate deletion table.

    Replaces 5-30-deletion_generate.py.
    """
    read_rip = {}
    rip_counter = defaultdict(int)
    seen_reads = set()

    with open(file1_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            read_base = parts[0].split(":", 1)[0].strip()
            rip_name = parts[1].strip()
            if read_base not in seen_reads:
                seen_reads.add(read_base)
                read_rip[read_base] = rip_name
                rip_counter[rip_name] += 1

    with open(file2_path) as f_in, open(output_path, "w") as f_out:
        for line in f_in:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rip_file2 = parts[0].strip()
            read_file2 = parts[1].strip()
            if read_file2 in read_rip and rip_file2 == read_rip[read_file2]:
                count = rip_counter.get(rip_file2, 0)
                f_out.write(f"{rip_file2}\t{read_file2}\t{count}\t-\n")


def _generate_genetype(file1_path, file2_path, output_path):
    """Classify RIPs by insertion/deletion presence and generate genotype.

    Replaces 5-30-genetype_generate.py.
    Markers: "++" insertion only, "--" deletion only, "+-" both, "na" neither.
    """
    file1_data = defaultdict(int)
    with open(file1_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if parts:
                file1_data[parts[0]] += 1

    file2_data = {}
    with open(file2_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                file2_data[parts[0]] = parts[2]

    all_rips = set(file1_data.keys()) | set(file2_data.keys())
    with open(output_path, "w") as fh:
        for rip in all_rips:
            in_f1 = rip in file1_data
            in_f2 = rip in file2_data
            if in_f1 and in_f2:
                fh.write(f"{rip}\t{file1_data[rip]}\t{file2_data[rip]}\t+-\n")
            elif in_f1:
                fh.write(f"{rip}\t{file1_data[rip]}\t0\t++\n")
            elif in_f2:
                fh.write(f"{rip}\t0\t{file2_data[rip]}\t--\n")
            else:
                fh.write(f"{rip}\t0\t0\tna\n")


def _run_pipeline(sample, cfg):
    """Execute the full 17-step genotyping pipeline for one sample."""
    output_dir = cfg["output_dir"]
    threads = cfg["threads_per_sample"]

    # Create sample directory structure
    sample_dir = os.path.join(output_dir, sample)
    intermediate_dir = os.path.join(sample_dir, "intermediate")
    log_dir = os.path.join(sample_dir, "logs")
    for d in (intermediate_dir, log_dir):
        os.makedirs(d, exist_ok=True)

    # Set up a per-sample file logger
    import logging
    sample_logger = logging.getLogger(f"genotype.{sample}")
    sample_logger.setLevel(logging.INFO)
    fh = logging.FileHandler(os.path.join(log_dir, f"{sample}.log"))
    fh.setFormatter(logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    ))
    sample_logger.addHandler(fh)

    pfx = intermediate_dir  # shorthand for intermediate path prefix
    total_steps = 17
    result_file = f"{sample}_genetype_result.tsv"
    final_result_path = os.path.join(output_dir, result_file)

    try:
        start_time = time.time()
        sample_logger.info("Starting sample %s", sample)

        # Step 1: Initial BLAST alignment
        sample_logger.info("Step 1/%d: Initial BLAST alignment", total_steps)
        subprocess.run([
            "blastn",
            "-query", cfg["query_fasta"],
            "-db", os.path.join(cfg["mkdb_dir"], sample),
            "-dust", "no",
            "-max_target_seqs", "10000000",
            "-max_hsps", "1",
            "-word_size", "15",
            "-perc_identity", "90",
            "-qcov_hsp_perc", "90",
            "-num_threads", str(threads),
            "-outfmt", "6",
            "-out", os.path.join(pfx, f"1-ripL5l50-{sample}.tbl"),
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Step 2: Extract unique matches
        sample_logger.info("Step 2/%d: Extract unique matches", total_steps)
        awk_cmd = (
            f"awk '{{count[$2]++; lines[$2] = $0}} END "
            f"{{for (id in count) if (count[id] == 1) print lines[id]}}' "
            f"{pfx}/1-ripL5l50-{sample}.tbl "
            f"> {pfx}/2-ripL5l50-{sample}-u.tbl"
        )
        subprocess.run(awk_cmd, shell=True, check=True)

        # Step 3: Convert to BED format
        sample_logger.info("Step 3/%d: Convert to BED format", total_steps)
        _convert_blast_to_bed(
            os.path.join(pfx, f"2-ripL5l50-{sample}-u.tbl"),
            os.path.join(pfx, f"3-ripL5l50-{sample}-u.bed"),
        )

        # Step 4: Extract 70bp flanking sequences
        sample_logger.info("Step 4/%d: Extract 70bp flanking sequences", total_steps)
        subprocess.run([
            "bedtools", "flank",
            "-s", "-l", "0", "-r", "70",
            "-g", os.path.join(cfg["mkdb_dir"], f"{sample}.fasta.fai"),
            "-i", os.path.join(pfx, f"3-ripL5l50-{sample}-u.bed"),
        ], check=True, stdout=open(os.path.join(pfx, f"4-ripL5l50-{sample}-u-r70.bed"), "w"),
           stderr=subprocess.PIPE)

        # Step 5: Filter short sequences (>30bp)
        sample_logger.info("Step 5/%d: Filter short sequences (>30bp)", total_steps)
        subprocess.run(
            f"awk -F '\\t' -v OFS='\\t' '$3 - $2 >30 {{print}}' "
            f"{pfx}/4-ripL5l50-{sample}-u-r70.bed "
            f"> {pfx}/5-ripL5l50-{sample}-u-r70-above30.bed",
            shell=True, check=True,
        )

        # Step 6: Get FASTA sequences
        sample_logger.info("Step 6/%d: Get FASTA sequences", total_steps)
        subprocess.run([
            "bedtools", "getfasta",
            "-s",
            "-fi", os.path.join(cfg["mkdb_dir"], f"{sample}.fasta"),
            "-bed", os.path.join(pfx, f"5-ripL5l50-{sample}-u-r70-above30.bed"),
            "-fo", os.path.join(pfx, f"5-ripL5l50-{sample}-u-r70-above30.fasta"),
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Step 7: Trim sequences (remove first 7bp)
        sample_logger.info("Step 7/%d: Trim sequences (remove first 7bp)", total_steps)
        awk_trim7 = (
            f"awk '/^>/ {{if (seq) print substr(seq, 8); print; seq=\"\"; next}} "
            f"{{seq = seq $0}} END {{if (seq) print substr(seq, 8)}}' "
            f"{pfx}/5-ripL5l50-{sample}-u-r70-above30.fasta "
            f"> {pfx}/5-ripL5l50-{sample}-u-r70-above30-trm7.fasta"
        )
        subprocess.run(awk_trim7, shell=True, check=True)

        # Step 8: SINE head alignment
        sample_logger.info("Step 8/%d: SINE head alignment", total_steps)
        subprocess.run([
            "blastn",
            "-dust", "no",
            "-max_target_seqs", "1",
            "-max_hsps", "1",
            "-word_size", "15",
            "-perc_identity", "95",
            "-qcov_hsp_perc", "90",
            "-num_threads", str(threads),
            "-outfmt", "6",
            "-db", os.path.join(cfg["db_path"], "SINE-head100"),
            "-query", os.path.join(pfx, f"5-ripL5l50-{sample}-u-r70-above30-trm7.fasta"),
            "-out", os.path.join(pfx, f"6-ripL5l50-{sample}-u-r70-above30-trm7-SINEhead100.tbl"),
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Step 9: Insertion site extraction
        sample_logger.info("Step 9/%d: Insertion site extraction", total_steps)
        _extract_sites(
            os.path.join(pfx, f"6-ripL5l50-{sample}-u-r70-above30-trm7-SINEhead100.tbl"),
            os.path.join(pfx, f"2-ripL5l50-{sample}-u.tbl"),
            os.path.join(pfx, f"7-ripL5l50-{sample}-insertion.tbl"),
        )

        # Step 10: Non-insertion site filtering
        sample_logger.info("Step 10/%d: Non-insertion site filtering", total_steps)
        subprocess.run(
            f"awk 'NR==FNR {{seen[$2]++; next}} !($1 in seen)' "
            f"{pfx}/7-ripL5l50-{sample}-insertion.tbl "
            f"{pfx}/3-ripL5l50-{sample}-u.bed "
            f"> {pfx}/8-ripL5l50-{sample}-noin.bed",
            shell=True, check=True,
        )

        # Step 11: Extract 100bp flanking sequences
        sample_logger.info("Step 11/%d: Extract 100bp flanking sequences", total_steps)
        subprocess.run([
            "bedtools", "flank",
            "-s", "-l", "0", "-r", "100",
            "-g", os.path.join(cfg["mkdb_dir"], f"{sample}.fasta.fai"),
            "-i", os.path.join(pfx, f"8-ripL5l50-{sample}-noin.bed"),
        ], check=True, stdout=open(os.path.join(pfx, f"9-ripL5l50-{sample}-noin-r100.bed"), "w"),
           stderr=subprocess.PIPE)

        # Step 12: Filter short sequences (>50bp)
        sample_logger.info("Step 12/%d: Filter short sequences (>50bp)", total_steps)
        subprocess.run(
            f"awk -F '\\t' -v OFS='\\t' '$3 - $2 >50 {{print}}' "
            f"{pfx}/9-ripL5l50-{sample}-noin-r100.bed "
            f"> {pfx}/10-ripL5l50-{sample}-noin-r100-above50.bed",
            shell=True, check=True,
        )

        # Step 13: Get FASTA sequences
        sample_logger.info("Step 13/%d: Get FASTA sequences", total_steps)
        subprocess.run([
            "bedtools", "getfasta",
            "-s",
            "-fi", os.path.join(cfg["mkdb_dir"], f"{sample}.fasta"),
            "-bed", os.path.join(pfx, f"10-ripL5l50-{sample}-noin-r100-above50.bed"),
            "-fo", os.path.join(pfx, f"10-ripL5l50-{sample}-noin-r100-above50.fasta"),
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Step 14: Trim sequences (remove first 20bp)
        sample_logger.info("Step 14/%d: Trim sequences (remove first 20bp)", total_steps)
        awk_trim20 = (
            f"awk '/^>/ {{if (seq) print substr(seq, 21); print; seq=\"\"; next}} "
            f"{{seq = seq $0}} END {{if (seq) print substr(seq, 21)}}' "
            f"{pfx}/10-ripL5l50-{sample}-noin-r100-above50.fasta "
            f"> {pfx}/11-ripL5l50-{sample}-noin-r100-above50-trim20.fasta"
        )
        subprocess.run(awk_trim20, shell=True, check=True)

        # Step 15: RIP region alignment
        sample_logger.info("Step 15/%d: RIP region alignment", total_steps)
        subprocess.run([
            "blastn",
            "-query", os.path.join(pfx, f"11-ripL5l50-{sample}-noin-r100-above50-trim20.fasta"),
            "-db", os.path.join(cfg["db_path"], "rip-r110"),
            "-dust", "no",
            "-max_target_seqs", "10",
            "-max_hsps", "1",
            "-word_size", "15",
            "-perc_identity", "90",
            "-qcov_hsp_perc", "90",
            "-num_threads", str(threads),
            "-outfmt", "6",
            "-out", os.path.join(pfx, f"12-ripL5l50-{sample}-noin-r100-above50-trim20-ripr110.tbl"),
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Step 16: Deletion site extraction
        sample_logger.info("Step 16/%d: Deletion site extraction", total_steps)
        _generate_deletions(
            os.path.join(pfx, f"12-ripL5l50-{sample}-noin-r100-above50-trim20-ripr110.tbl"),
            os.path.join(pfx, f"2-ripL5l50-{sample}-u.tbl"),
            os.path.join(pfx, f"13-ripL5l50-{sample}-deletion.tbl"),
        )

        # Step 17: Genotype generation
        sample_logger.info("Step 17/%d: Genotype generation", total_steps)
        _generate_genetype(
            os.path.join(pfx, f"7-ripL5l50-{sample}-insertion.tbl"),
            os.path.join(pfx, f"13-ripL5l50-{sample}-deletion.tbl"),
            os.path.join(pfx, result_file),
        )

        # Move final result to main output directory
        shutil.move(os.path.join(pfx, result_file), final_result_path)

        elapsed = time.time() - start_time
        sample_logger.info(
            "Sample %s completed successfully in %.2f seconds -> %s",
            sample, elapsed, final_result_path,
        )
        return (sample, True, elapsed, None)

    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        stderr_msg = e.stderr.decode("utf-8", errors="replace") if e.stderr else ""
        err_msg = f"Command failed: {e.cmd}\n{stderr_msg}"
        sample_logger.error("Sample %s FAILED: %s", sample, err_msg)
        return (sample, False, elapsed, err_msg)
    except Exception as e:
        elapsed = time.time() - start_time
        err_msg = f"{e}\n{traceback.format_exc()}"
        sample_logger.error("Sample %s FAILED: %s", sample, err_msg)
        return (sample, False, elapsed, err_msg)
    finally:
        for handler in sample_logger.handlers[:]:
            handler.close()
            sample_logger.removeHandler(handler)


def run(args):
    """Execute the genotyping pipeline."""
    # Read sample list from mkdb output
    mkdb_dir = os.path.abspath(args.mkdb_dir)
    if not os.path.isdir(mkdb_dir):
        raise FileNotFoundError(f"mkdb output directory not found: {mkdb_dir}")

    sample_list_file = os.path.join(mkdb_dir, "sample_list.txt")
    if not os.path.isfile(sample_list_file):
        raise FileNotFoundError(
            f"sample_list.txt not found in {mkdb_dir}. "
            "Please run 'sriptype mkdb' first."
        )
    with open(sample_list_file) as fh:
        samples = [line.strip() for line in fh if line.strip()]
    if not samples:
        raise ValueError(f"No samples found in {sample_list_file}")

    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Build config dict (query_fasta and db_path use bundled data)
    cfg = {
        "mkdb_dir": mkdb_dir,
        "db_path": _DATA_DIR,
        "query_fasta": os.path.join(_DATA_DIR, "rip-shiftL5-l50.fasta"),
        "threads_per_sample": args.threads,
        "output_dir": output_dir,
    }

    jobs = args.jobs

    # Validate bundled data
    if not os.path.isfile(cfg["query_fasta"]):
        raise FileNotFoundError(f"Bundled query FASTA not found: {cfg['query_fasta']}")

    # Check external tools
    _check_tools(REQUIRED_TOOLS)

    logger.info("Starting genotype analysis for %d sample(s)", len(samples))
    logger.info("  Samples:           %s", ", ".join(samples))
    logger.info("  Parallel jobs:     %d", jobs)
    logger.info("  Threads/sample:    %d", cfg["threads_per_sample"])
    logger.info("  Output directory:  %s", output_dir)

    overall_start = time.time()
    success_count = 0
    failed_samples = []

    with ProcessPoolExecutor(max_workers=jobs) as executor:
        futures = {
            executor.submit(_run_pipeline, sample, cfg): sample
            for sample in samples
        }
        for future in as_completed(futures):
            sample_name = futures[future]
            try:
                sample, ok, elapsed, err = future.result()
                if ok:
                    success_count += 1
                    logger.info("[OK]   %s (%.1f s)", sample, elapsed)
                else:
                    failed_samples.append(sample)
                    logger.error("[FAIL] %s (%.1f s): %s", sample, elapsed, err)
            except Exception as e:
                failed_samples.append(sample_name)
                logger.error("[FAIL] %s: %s", sample_name, e)

    total_time = time.time() - overall_start

    # Summary

    logger.info("Genotype analysis completed")
    logger.info("  Total samples:  %d", len(samples))
    logger.info("  Succeeded:      %d", success_count)
    logger.info("  Failed:         %d", len(failed_samples))
    if failed_samples:
        logger.info("  Failed samples: %s", ", ".join(failed_samples))
    logger.info("  Total time:     %.2f s", total_time)
    logger.info("  Avg per sample: %.2f s", total_time / len(samples))


    if failed_samples:
        raise RuntimeError(
            f"{len(failed_samples)} sample(s) failed: {', '.join(failed_samples)}. "
            "Check per-sample logs in <output-dir>/<sample>/logs/ for details."
        )
