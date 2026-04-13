"""Merge and summarize genotype results.

Combines per-sample genotype results into a single matrix, filters sites
by minimum non-NA threshold (25% of samples), and generates summary
statistics including allele frequencies and genotype counts.

Outputs:
  - Merged genotype matrix (merged_genotype.tsv)
  - Per-site summary statistics (merged_genotype_sum.tsv)
  - Summary report (merged_genotype_report.txt)
"""

import glob
import math
import os
from collections import defaultdict

from sriptype_modules.utils import setup_logger

logger = setup_logger(__name__)

DESCRIPTION = "Merge per-sample genotype results and generate summary statistics"
ORDER = 3
USES_THREADS = False

# Bundled data directory
_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")


def add_arguments(parser):
    """Add subcommand-specific arguments to the parser."""
    parser.add_argument(
        "-i", "--input-dir",
        required=True,
        help="Directory containing *_genetype_result.tsv files (output of 'sriptype genotype')",
    )
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Output directory for merged results",
    )
    parser.add_argument(
        "--min-rate",
        type=float,
        default=0,
        help="Maximum NA rate per site; sites with NA/total >= this value are excluded from the merged matrix (default: 0, no filtering)",
    )


def _read_loci(locus_file):
    """Read reference locus BED file.

    Returns list of dicts with keys: chr, start, end, locus_id, strand.
    """
    loci = []
    with open(locus_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            loci.append({
                "chr": parts[0],
                "start": parts[1],
                "end": parts[2],
                "locus_id": parts[3],
                "strand": parts[4],
            })
    return loci


def _read_genotype_file(file_path):
    """Read a single genotype result file.

    Returns dict: locus_id -> genotype string.
    """
    genotype_dict = {}
    count = 0
    with open(file_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            genotype_dict[parts[0]] = parts[3]
            count += 1
    return genotype_dict, count


def _median(values):
    """Compute median of a list of numbers."""
    if not values:
        return 0.0
    s = sorted(values)
    n = len(s)
    if n % 2 == 1:
        return float(s[n // 2])
    return (s[n // 2 - 1] + s[n // 2]) / 2.0


def _mean(values):
    """Compute mean of a list of numbers."""
    if not values:
        return 0.0
    return sum(values) / len(values)


def run(args):
    """Execute the merge pipeline."""
    input_dir = os.path.abspath(args.input_dir)
    locus_file = os.path.join(_DATA_DIR, "51krip-array-ref.txt")
    output_dir = os.path.abspath(args.output_dir)
    min_rate = args.min_rate

    # Validate inputs
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if not os.path.isfile(locus_file):
        raise FileNotFoundError(f"Bundled locus file not found: {locus_file}")

    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)

    # Output file paths
    output_file = os.path.join(output_dir, "merged_genotype.tsv")
    summary_file = os.path.join(output_dir, "merged_genotype_sum.tsv")
    report_file = os.path.join(output_dir, "merged_genotype_report.txt")

    # Auto-discover genotype result files
    pattern = os.path.join(input_dir, "*_genetype_result.tsv")
    file_names = sorted(
        os.path.basename(f) for f in glob.glob(pattern)
    )

    if not file_names:
        raise FileNotFoundError(
            f"No genotype result files found in {input_dir}. "
            "Run 'sriptype genotype' first."
        )

    # Validate and read genotype files
    logger.info("Running merge pipeline...")
    valid_file_names = []
    genotype_data = {}

    for fname in file_names:
        fpath = os.path.join(input_dir, fname)
        if not os.path.isfile(fpath):
            logger.warning("File not found, skipping: %s", fpath)
            continue
        gdict, count = _read_genotype_file(fpath)
        genotype_data[fname] = gdict
        valid_file_names.append(fname)
        logger.debug("Read %s (%d loci)", fname, count)

    if not valid_file_names:
        raise FileNotFoundError("No valid genotype files to process")

    num_individuals = len(valid_file_names)
    logger.info("Read %d genotype file(s)", num_individuals)

    # Read locus file
    loci = _read_loci(locus_file)
    if not loci:
        raise ValueError(f"No loci found in {locus_file}")
    logger.info("Read %d loci from reference file", len(loci))

    # Build merged matrix: for each locus, look up genotype in each sample
    # Build data rows: base_info (5 cols) + genotypes (num_individuals cols)
    base_fields = ["chr", "start", "end", "locus_id", "strand"]
    data_rows = []
    for locus in loci:
        row = [locus[k] for k in base_fields]
        for fname in valid_file_names:
            gt = genotype_data[fname].get(locus["locus_id"], "NA")
            row.append(gt)
        data_rows.append(row)

    total_loci = len(data_rows)
    logger.info("Merged %d loci x %d samples", total_loci, num_individuals)

    # ---- Filter and compute stats ----
    sample_start = 5
    sample_end = sample_start + num_individuals

    all_site_stats = []
    all_sample_non_na_counts = [0] * num_individuals

    filtered_rows = []
    filtered_sample_non_na_counts = [0] * num_individuals

    for row in data_rows:
        genotypes = row[sample_start:sample_end]

        count_pp = genotypes.count("++")
        count_pm = genotypes.count("+-")
        count_mm = genotypes.count("--")
        count_na = genotypes.count("NA")

        total_ind = len(genotypes)
        perc_pp = round((count_pp / total_ind) * 100, 2) if total_ind else 0.0
        perc_pm = round((count_pm / total_ind) * 100, 2) if total_ind else 0.0
        perc_mm = round((count_mm / total_ind) * 100, 2) if total_ind else 0.0
        perc_na = round((count_na / total_ind) * 100, 2) if total_ind else 0.0

        plus_alleles = 2 * count_pp + count_pm
        minus_alleles = 2 * count_mm + count_pm
        total_alleles = plus_alleles + minus_alleles
        if total_alleles > 0:
            perc_plus = round((plus_alleles / total_alleles) * 100, 2)
            perc_minus = round((minus_alleles / total_alleles) * 100, 2)
        else:
            perc_plus = perc_minus = 0.0

        for idx, gt in enumerate(genotypes):
            if gt != "NA":
                all_sample_non_na_counts[idx] += 1

        stat = {
            "base_info": row[:5],
            "count_pp": count_pp,
            "count_pm": count_pm,
            "count_mm": count_mm,
            "count_na": count_na,
            "perc_pp": perc_pp,
            "perc_pm": perc_pm,
            "perc_mm": perc_mm,
            "perc_na": perc_na,
            "plus_alleles": plus_alleles,
            "minus_alleles": minus_alleles,
            "perc_plus": perc_plus,
            "perc_minus": perc_minus,
        }
        all_site_stats.append(stat)

        # Apply filtering for merged matrix
        if min_rate > 0 and count_na / num_individuals >= min_rate:
            continue
        filtered_rows.append(row)
        for idx, gt in enumerate(genotypes):
            if gt != "NA":
                filtered_sample_non_na_counts[idx] += 1

    total_retained = len(filtered_rows)
    if min_rate > 0:
        logger.info(
            "Filtered: %d / %d sites retained (NA rate < %.0f%%)",
            total_retained, total_loci, min_rate * 100,
        )
    else:
        logger.info("All %d sites included (no filtering)", total_loci)

    # ---- Write merged genotype matrix ----
    header1 = base_fields + valid_file_names
    header2 = ["" ] * 5 + [str(c) for c in filtered_sample_non_na_counts]

    with open(output_file, "w") as fout:
        fout.write("\t".join(header1) + "\n")
        fout.write("\t".join(header2) + "\n")
        for row in filtered_rows:
            fout.write("\t".join(str(x) for x in row[:sample_end]) + "\n")

    logger.info("Wrote merged matrix: %s", output_file)

    # ---- Write summary statistics (all sites, unfiltered) ----
    summary_header = base_fields + [
        "++", "+-", "--", "NA", "++%", "+-%", "--%", "NA%",
        "+", "-", "+%", "-%",
    ]

    with open(summary_file, "w") as fsum:
        fsum.write("\t".join(summary_header) + "\n")
        for stat in all_site_stats:
            parts = stat["base_info"] + [
                str(stat["count_pp"]), str(stat["count_pm"]),
                str(stat["count_mm"]), str(stat["count_na"]),
                str(stat["perc_pp"]), str(stat["perc_pm"]),
                str(stat["perc_mm"]), str(stat["perc_na"]),
                str(stat["plus_alleles"]), str(stat["minus_alleles"]),
                str(stat["perc_plus"]), str(stat["perc_minus"]),
            ]
            fsum.write("\t".join(parts) + "\n")

    logger.info("Wrote summary: %s", summary_file)

    # ---- Write report ----
    max_s = max(all_sample_non_na_counts) if all_sample_non_na_counts else 0
    min_s = min(all_sample_non_na_counts) if all_sample_non_na_counts else 0
    median_s = _median(all_sample_non_na_counts)
    mean_s = _mean(all_sample_non_na_counts)

    with open(report_file, "w") as frpt:
        frpt.write("# SRIPType Merge Summary Report\n\n")
        frpt.write(f"Total samples\t{num_individuals}\n")
        frpt.write(f"Initial loci\t{total_loci}\n")
        if min_rate > 0:
            frpt.write(f"Loci with NA rate < {min_rate:.0%} (retained)\t{total_retained}\n")
        frpt.write(f"Max successfully genotyped loci per sample\t{max_s}\n")
        frpt.write(f"Min successfully genotyped loci per sample\t{min_s}\n")
        frpt.write(f"Median successfully genotyped loci per sample\t{median_s:.2f}\n")
        frpt.write(f"Mean successfully genotyped loci per sample\t{mean_s:.2f}\n")

    logger.info("Wrote report: %s", report_file)

    # ---- Print report ----
    logger.info("Pipeline completed successfully!")
    logger.info("Summary:")
    logger.info("  Total samples:              %d", num_individuals)
    logger.info("  Initial loci:               %d", total_loci)
    if min_rate > 0:
        logger.info("  Retained loci (NA rate < %.0f%%): %d (%.2f%%)",
                    min_rate * 100, total_retained,
                    total_retained / total_loci * 100 if total_loci else 0)
    logger.info("  Max genotyped loci/sample:  %d", max_s)
    logger.info("  Min genotyped loci/sample:  %d", min_s)
    logger.info("  Median genotyped loci/sample: %.2f", median_s)
    logger.info("  Mean genotyped loci/sample: %.2f", mean_s)
