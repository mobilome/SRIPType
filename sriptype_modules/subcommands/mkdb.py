"""Data preparation subcommand.

Paired-end sequencing data processing and BLAST database construction pipeline:
  1. FLASH merge paired-end reads
  2. Merge FLASH output files (extendedFrags + notCombined)
  3. Convert merged FASTQ to FASTA with renamed headers
  4. Clean up intermediate files
  5. Extract sample IDs
  6. Build BLAST databases
  7. Compute sequence length statistics
"""

import glob
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

from sriptype_modules.core import validate_threads
from sriptype_modules.utils import setup_logger

logger = setup_logger(__name__)

DESCRIPTION = "Paired-end data processing and BLAST database construction"
ORDER = 1

# All recognized paired-end suffix patterns: (R1_suffix, R2_suffix)
PAIRED_SUFFIXES = [
    # Compressed pairs
    (".R1.fq.gz", ".R2.fq.gz"),
    ("_R1.fq.gz", "_R2.fq.gz"),
    (".R1.fastq.gz", ".R2.fastq.gz"),
    ("_R1.fastq.gz", "_R2.fastq.gz"),
    (".r1.fq.gz", ".r2.fq.gz"),
    ("_r1.fq.gz", "_r2.fq.gz"),
    (".r1.fastq.gz", ".r2.fastq.gz"),
    ("_r1.fastq.gz", "_r2.fastq.gz"),
    ("_1.fq.gz", "_2.fq.gz"),
    (".1.fq.gz", ".2.fq.gz"),
    ("_1.fastq.gz", "_2.fastq.gz"),
    (".1.fastq.gz", ".2.fastq.gz"),
    (".f1.fq.gz", ".r2.fq.gz"),
    ("_f1.fq.gz", "_r2.fq.gz"),
    (".f1.fastq.gz", ".r2.fastq.gz"),
    ("_f1.fastq.gz", "_r2.fastq.gz"),
    # Uncompressed pairs
    (".R1.fq", ".R2.fq"),
    ("_R1.fq", "_R2.fq"),
    (".R1.fastq", ".R2.fastq"),
    ("_R1.fastq", "_R2.fastq"),
    (".r1.fq", ".r2.fq"),
    ("_r1.fq", "_r2.fq"),
    (".r1.fastq", ".r2.fastq"),
    ("_r1.fastq", "_r2.fastq"),
    ("_1.fq", "_2.fq"),
    (".1.fq", ".2.fq"),
    ("_1.fastq", "_2.fastq"),
    (".1.fastq", ".2.fastq"),
    (".f1.fq", ".r2.fq"),
    ("_f1.fq", "_r2.fq"),
    (".f1.fastq", ".r2.fastq"),
    ("_f1.fastq", "_r2.fastq"),
]


def add_arguments(parser):
    """Add subcommand-specific arguments to the parser."""
    parser.add_argument(
        "-i", "--input-dir",
        required=True,
        help="Directory containing paired-end .R1.fq.gz / .R2.fq.gz files",
    )
    parser.add_argument(
        "-o", "--output-dir",
        default=None,
        help="Main output directory (default: <input-dir>/../blast_db)",
    )
    parser.add_argument(
        "-j", "--jobs",
        type=int,
        default=1,
        help="Number of parallel jobs (default: 1)",
    )


def _logged_run(cmd, shell=False, **kwargs):
    """Run subprocess with debug command logging."""
    if shell:
        logger.debug("CMD: %s", cmd)
    else:
        logger.debug("CMD: %s", " ".join(str(a) for a in cmd))
    return subprocess.run(cmd, shell=shell, **kwargs)


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


def _find_samples(input_dir):
    """Auto-detect paired-end samples with various naming conventions.

    Returns a list of (prefix, r1_suffix, r2_suffix) tuples.
    """
    # Collect all files under input_dir
    all_files = set()
    for root, _dirs, files in os.walk(input_dir):
        for f in files:
            all_files.add(os.path.join(root, f))

    # Sort suffix pairs by R1 suffix length (longest first) to avoid
    # shorter suffixes matching a substring of a longer one.
    sorted_suffixes = sorted(PAIRED_SUFFIXES, key=lambda x: len(x[0]), reverse=True)

    samples = []  # (prefix, r1_suffix, r2_suffix)
    matched_r1 = set()

    for fpath in sorted(all_files):
        if fpath in matched_r1:
            continue
        for r1_sfx, r2_sfx in sorted_suffixes:
            if fpath.endswith(r1_sfx):
                prefix = fpath[: -len(r1_sfx)]
                r2_path = prefix + r2_sfx
                if r2_path in all_files:
                    samples.append((prefix, r1_sfx, r2_sfx))
                    matched_r1.add(fpath)
                    break

    if not samples:
        raise FileNotFoundError(
            f"No paired-end sequencing files found in {input_dir}. "
            "Supported naming patterns include: .R1.fq.gz/.R2.fq.gz, "
            "_1.fastq/_2.fastq, etc."
        )

    # Log detected naming patterns
    patterns_found = sorted(set((r1, r2) for _, r1, r2 in samples))
    for r1, r2 in patterns_found:
        count = sum(1 for _, a, b in samples if a == r1 and b == r2)
        logger.debug("  Detected pattern %s / %s  (%d sample(s))", r1, r2, count)

    return samples


# --------------- Step-level worker functions (picklable for ProcessPoolExecutor) ---------------

def _worker_flash(args):
    """Worker: decompress (if needed) + FLASH merge for one sample."""
    prefix, r1_suffix, r2_suffix, flash_output_dir, threads = args
    sample = os.path.basename(prefix)
    r1_file = prefix + r1_suffix
    r2_file = prefix + r2_suffix

    r1_is_gz = r1_suffix.endswith(".gz")
    r2_is_gz = r2_suffix.endswith(".gz")
    temp_files = []

    try:
        # Prepare R1: decompress with pigz if gzipped
        if r1_is_gz:
            r1_fq = prefix + ".R1.tmp.fq"
            temp_files.append(r1_fq)
            with open(r1_fq, "wb") as f:
                _logged_run(
                    ["pigz", "-dc", "-p", str(threads), r1_file],
                    stdout=f, stderr=subprocess.PIPE, check=True,
                )
        else:
            r1_fq = r1_file

        # Prepare R2: decompress with pigz if gzipped
        if r2_is_gz:
            r2_fq = prefix + ".R2.tmp.fq"
            temp_files.append(r2_fq)
            with open(r2_fq, "wb") as f:
                _logged_run(
                    ["pigz", "-dc", "-p", str(threads), r2_file],
                    stdout=f, stderr=subprocess.PIPE, check=True,
                )
        else:
            r2_fq = r2_file

        # Run FLASH
        _logged_run(
            ["flash", r1_fq, r2_fq, "-t", str(threads),
             "-d", flash_output_dir, "-o", sample],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
        )
        return (sample, True, None)
    except Exception as e:
        return (sample, False, str(e))
    finally:
        # Clean temp decompressed files only
        for f in temp_files:
            if os.path.exists(f):
                os.remove(f)


def _worker_merge(args):
    """Worker: merge FLASH output files for one sample."""
    prefix = args
    sample = os.path.basename(prefix)
    frag_file = prefix + ".extendedFrags.fastq"
    not1_file = prefix + ".notCombined_1.fastq"
    not2_file = prefix + ".notCombined_2.fastq"
    output_fastq = prefix + ".merged.fastq"

    for f in (frag_file, not1_file, not2_file):
        if not os.path.isfile(f):
            return (sample, False, f"Missing file: {f}")

    try:
        with open(output_fastq, "wb") as out:
            for f in (frag_file, not1_file, not2_file):
                with open(f, "rb") as inp:
                    while True:
                        chunk = inp.read(1024 * 1024)
                        if not chunk:
                            break
                        out.write(chunk)
        return (sample, True, None)
    except Exception as e:
        return (sample, False, str(e))


def _worker_convert(args):
    """Worker: convert merged FASTQ to FASTA with renamed headers."""
    prefix = args
    sample = os.path.basename(prefix)
    input_fastq = prefix + ".merged.fastq"
    output_fasta = prefix + ".fasta"

    if not os.path.isfile(input_fastq):
        return (sample, False, f"Missing file: {input_fastq}")

    try:
        # Run seqkit fq2fa and pipe through awk for header renaming
        awk_script = 'BEGIN {i=0} /^>/ {print ">" "r" ++i; next} {print}'
        logger.debug("CMD: seqkit fq2fa %s | awk '%s'", input_fastq, awk_script)
        seqkit_proc = subprocess.Popen(
            ["seqkit", "fq2fa", input_fastq],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        awk_proc = subprocess.Popen(
            ["awk", awk_script],
            stdin=seqkit_proc.stdout,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        seqkit_proc.stdout.close()
        awk_out, awk_err = awk_proc.communicate()
        seqkit_proc.wait()

        if seqkit_proc.returncode != 0:
            stderr_msg = seqkit_proc.stderr.read().decode("utf-8", errors="replace")
            return (sample, False, f"seqkit failed: {stderr_msg}")
        if awk_proc.returncode != 0:
            return (sample, False, f"awk failed: {awk_err.decode('utf-8', errors='replace')}")

        with open(output_fasta, "wb") as f:
            f.write(awk_out)
        return (sample, True, None)
    except Exception as e:
        return (sample, False, str(e))


def _worker_makeblastdb(args):
    """Worker: build BLAST database for one sample."""
    sample_path, blast_db_dir = args
    sample = os.path.basename(sample_path).replace(".fasta", "")
    db_out = os.path.join(blast_db_dir, sample)

    try:
        _logged_run(
            ["makeblastdb", "-in", sample_path,
             "-input_type", "fasta", "-dbtype", "nucl",
             "-out", db_out],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
        )
        return (sample, True, None)
    except Exception as e:
        return (sample, False, str(e))


def _worker_seq_stats(args):
    """Worker: compute sequence length statistics for one sample."""
    sample_path, stats_dir = args
    sample = os.path.basename(sample_path).replace(".fasta", "")

    try:
        _logged_run(
            ["samtools", "faidx", sample_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
        )
        return (sample, True, None)
    except Exception as e:
        return (sample, False, str(e))


def _run_parallel(worker_func, tasks, max_workers, step_name):
    """Run a worker function in parallel over a list of tasks."""
    succeeded = 0
    failed_items = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(worker_func, t): t for t in tasks}
        for future in as_completed(futures):
            sample, ok, err = future.result()
            if ok:
                succeeded += 1
                logger.debug("  [OK] %s", sample)
            else:
                failed_items.append((sample, err))
    if failed_items:
        for sample, err in failed_items:
            logger.error("  [FAIL] %s: %s", sample, err)
        logger.debug("%s: %d succeeded, %d failed", step_name, succeeded, len(failed_items))
        raise RuntimeError(f"{step_name}: {len(failed_items)} sample(s) failed. Check logs above.")
    logger.debug("%s: all %d sample(s) succeeded", step_name, succeeded)


def run(args):
    """Execute the data-prepare pipeline."""
    input_dir = os.path.abspath(args.input_dir)
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    # Resolve output directories
    if args.output_dir:
        main_output = os.path.abspath(args.output_dir)
    else:
        main_output = os.path.join(os.path.dirname(input_dir), "blast_db")

    jobs = args.jobs
    threads = args.threads

    # Create output directory
    logger.info("Starting database creation...")
    os.makedirs(main_output, exist_ok=True)
    logger.debug("Output: %s", main_output)

    # =================== Part 1: FLASH merge & convert ===================

    # Find samples (auto-detect naming convention)
    samples = _find_samples(input_dir)
    logger.info("Found %d paired-end sample(s)", len(samples))

    # Check required tools (pigz only needed when gzipped files exist)
    has_compressed = any(
        r1.endswith(".gz") or r2.endswith(".gz")
        for _, r1, r2 in samples
    )
    required_tools = ["flash", "seqkit"]
    if has_compressed:
        required_tools.append("pigz")
    _check_tools(required_tools)

    # Step 1: FLASH merge
    logger.debug("Step 1/7: FLASH merging paired-end reads...")
    flash_tasks = [
        (prefix, r1_sfx, r2_sfx, main_output, threads)
        for prefix, r1_sfx, r2_sfx in samples
    ]
    _run_parallel(_worker_flash, flash_tasks, jobs, "FLASH merge")

    # Step 2: Merge FLASH output files
    logger.debug("Step 2/7: Merging FLASH output files...")
    frag_files = sorted(glob.glob(os.path.join(main_output, "*.extendedFrags.fastq")))
    merge_prefixes = [f[: -len(".extendedFrags.fastq")] for f in frag_files]
    _run_parallel(_worker_merge, merge_prefixes, jobs, "File merge")

    # Step 3: Convert to FASTA
    logger.debug("Step 3/7: Converting FASTQ to FASTA...")
    merged_files = sorted(glob.glob(os.path.join(main_output, "*.merged.fastq")))
    convert_prefixes = [f[: -len(".merged.fastq")] for f in merged_files]
    _run_parallel(_worker_convert, convert_prefixes, jobs, "FASTQ to FASTA")

    # Step 4: Clean intermediate files
    logger.debug("Step 4/7: Cleaning intermediate files...")
    patterns = [
        "*.extendedFrags.fastq",
        "*.notCombined_1.fastq",
        "*.notCombined_2.fastq",
        "*.merged.fastq",
        "*.hist",
        "*.histogram",
    ]
    removed = 0
    for pat in patterns:
        for f in glob.glob(os.path.join(main_output, pat)):
            os.remove(f)
            removed += 1
    logger.debug("  Removed %d intermediate file(s)", removed)

    # =================== Part 2: BLAST DB & seq stats ===================

    # Check required tools for Part 2
    _check_tools(["makeblastdb", "samtools"])

    # Find FASTA files
    fasta_files = sorted(glob.glob(os.path.join(main_output, "*.fasta")))
    if not fasta_files:
        raise FileNotFoundError(
            f"No *.fasta files found in {main_output}. "
            "Ensure FLASH processing was completed first."
        )

    # Step 5: Extract sample IDs
    logger.debug("Step 5/7: Extracting sample IDs...")
    sample_id_file = os.path.join(main_output, "sample_list.txt")
    sample_names = [
        os.path.basename(f).replace(".fasta", "") for f in fasta_files
    ]
    with open(sample_id_file, "w") as fh:
        fh.write("\n".join(sample_names) + "\n")
    logger.debug("  Wrote %d sample ID(s) to %s", len(sample_names), sample_id_file)

    # Step 6: Build BLAST databases
    logger.debug("Step 6/7: Building BLAST databases...")
    blast_tasks = [(f, main_output) for f in fasta_files]
    _run_parallel(_worker_makeblastdb, blast_tasks, jobs, "BLAST DB build")

    # Step 7: Sequence length statistics
    logger.debug("Step 7/7: Computing sequence length statistics...")
    stats_tasks = [(f, main_output) for f in fasta_files]
    _run_parallel(_worker_seq_stats, stats_tasks, jobs, "Sequence stats")

    # Final report
    n_fasta = len(fasta_files)
    n_blastdb = len(glob.glob(os.path.join(main_output, "*.n*")))
    n_stats = len(glob.glob(os.path.join(main_output, "*.fasta.fai")))

    logger.info("make database completed successfully!")
    logger.info("Output directory: %s", main_output)
    logger.info("Summary:")
    logger.info("  - Samples processed: %d", len(sample_names))
    logger.info("  - FASTA files:       %d", n_fasta)
    logger.info("  - BLAST DB files:    %d", n_blastdb)
    logger.info("  - Stats files:       %d", n_stats)