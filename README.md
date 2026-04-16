# SRIPType

A bioinformatics tool for 51k-SINE-RIP chip genotyping.

## Dependencies

SRIPType requires the following external tools to be installed and available in `PATH`:

- [FLASh](https://ccb.jhu.edu/software/FLASH/) — paired-end read merging
- [pigz](https://zlib.net/pigz/) — parallel gzip compression/decompression
- [SeqKit](https://bioinf.shenwei.me/seqkit/) — FASTA/FASTQ manipulation
- [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) — sequence alignment (`blastn`, `makeblastdb`)
- [BEDTools](https://bedtools.readthedocs.io/) — genomic interval operations
- [SAMtools](http://www.htslib.org/) — FASTA indexing (`samtools faidx`)

## Installation

### Conda (recommended)

```bash
conda install -c bioconda sriptype
```

### pip

```bash
pip install sriptype
```

> **Note:** pip only installs SRIPType itself; the external tools listed above must be installed separately. Installing via Conda handles all dependencies automatically.

### From source

```bash
git clone https://github.com/mobilome/SRIPType.git
cd SRIPType
pip install .
```

### Docker

```bash
docker build -t sriptype .
docker run --rm -v $(pwd):/data sriptype mkdb -i /data/input_dir -o /data/blast_db -j 4 -t 4
```

## Quick start

```bash
# Step 1: Build BLAST databases from paired-end data
sriptype mkdb -i raw_data/ -o blast_db/ -j 4 -t 4

# Step 2: Run genotyping analysis
sriptype genotype -i blast_db/ -o genotype_results/ -j 4 -t 4

# Step 3: Merge per-sample results into a summary matrix
sriptype merge -i genotype_results/ -o merged_output/
```

## Usage

```
sriptype [-v] <subcommand> [options]
```

### Subcommands

| Subcommand | Description |
|------------|-------------|
| `mkdb` | Paired-end data processing and BLAST database construction |
| `genotype` | RIP-seq genotyping analysis pipeline |
| `merge` | Merge per-sample genotype results and generate summary statistics |

### `mkdb` — Build BLAST databases

Process paired-end FASTQ files and create per-sample BLAST databases.

```bash
sriptype mkdb -i <input_dir> -o <output_dir> -j <jobs> -t <threads>
```

| Option | Description |
|--------|-------------|
| `-i, --input-dir` | Directory containing paired-end `.R1.fq.gz` / `.R2.fq.gz` files (**required**) |
| `-o, --output-dir` | Output directory (default: `<input-dir>/../blast_db`) |
| `-j, --jobs INT` | Number of parallel jobs (default: 1) |
| `-t, --threads INT` | Number of threads per job (default: 1) |
| `--verbose` | Enable verbose output |

### `genotype` — Genotyping analysis

Run the RIP-seq genotyping pipeline on each sample.

```bash
sriptype genotype -i <mkdb_dir> -o <output_dir> -j <jobs> -t <threads>
```

| Option | Description |
|--------|-------------|
| `-i, --mkdb-dir` | Directory produced by `sriptype mkdb` (**required**) |
| `-o, --output-dir` | Output directory for results (**required**) |
| `-j, --jobs INT` | Number of samples to process in parallel (default: 1) |
| `-t, --threads INT` | Number of threads per sample (default: 1) |
| `--verbose` | Enable verbose output |

### `merge` — Merge results

Merge per-sample genotype results into a unified matrix with summary statistics.

```bash
sriptype merge -i <input_dir> -o <output_dir> [--min-rate FLOAT]
```

| Option | Description |
|--------|-------------|
| `-i, --input-dir` | Directory containing `*_genetype_result.tsv` files (**required**) |
| `-o, --output-dir` | Output directory for merged results (**required**) |
| `--min-rate FLOAT` | Max NA rate per site; sites with NA/total ≥ this value are excluded (default: 0, no filtering) |
| `--verbose` | Enable verbose output |

Output files:

| File | Description |
|------|-------------|
| `merged_genotype.tsv` | Merged genotype matrix (filtered by `--min-rate`) |
| `merged_genotype_sum.tsv` | Per-site summary statistics (all sites, unfiltered) |
| `merged_genotype_report.txt` | Run summary report |

## License

MIT License
