# SRIPType

A bioinformatics toolkit for sequence typing and analysis.

## Installation

### Method 1: Direct execution (clone and run)

```bash
git clone https://github.com/your-org/SRIPType.git
cd SRIPType

# Run via absolute path
/path/to/SRIPType/sriptype --help

# Or add to PATH
export PATH=/path/to/SRIPType:$PATH
sriptype --help
```

### Method 2: pip install

```bash
git clone https://github.com/your-org/SRIPType.git
cd SRIPType
pip install .

sriptype --help
```

### Method 3: Docker

```bash
# Build
docker build -t sriptype .

# Run
docker run --rm -v $(pwd):/data sriptype example -i /data/input.fasta -o /data/output.txt -t 4
```

### Method 4: Conda (coming soon)

```bash
conda install -c your-channel sriptype
```

## Usage

```bash
# Show help
sriptype --help

# Show version
sriptype --version

# Run a subcommand
sriptype example -i input.fasta -o output.txt -t 4
```

### Global options

| Option | Description |
|--------|-------------|
| `-v, --version` | Show version and exit |
| `-t, --threads INT` | Number of threads (default: 1) |
| `--verbose` | Enable verbose output |

### Available subcommands

| Subcommand | Description |
|------------|-------------|
| `example` | Example subcommand (template) |

## Adding a new subcommand

1. Create a new Python file in `sriptype_modules/subcommands/`, e.g. `mycommand.py`
2. Define the following in the module:

```python
"""My new subcommand."""

from sriptype_modules.utils import check_file_exists, ensure_output_dir, setup_logger

logger = setup_logger(__name__)

DESCRIPTION = "Description of my subcommand"


def add_arguments(parser):
    """Add subcommand-specific arguments."""
    parser.add_argument("-i", "--input", required=True, help="Input file")
    parser.add_argument("-o", "--output", required=True, help="Output file")


def run(args):
    """Execute the subcommand."""
    input_file = check_file_exists(args.input)
    output_file = ensure_output_dir(args.output)
    threads = args.threads

    # Your logic here
    logger.info("Done!")
```

3. The subcommand will be automatically discovered and available as `sriptype mycommand`.

## Project Structure

```
SRIPType/
├── sriptype                      # Executable entry point
├── sriptype_modules/             # Core Python package
│   ├── __init__.py               # Version info
│   ├── cli.py                    # CLI argument parsing
│   ├── core.py                   # Core shared logic
│   ├── utils.py                  # Utility functions
│   └── subcommands/              # Subcommand modules
│       ├── __init__.py
│       └── example.py            # Example/template subcommand
├── tests/                        # Test directory
├── setup.py                      # Package setup
├── pyproject.toml                # Modern Python packaging config
├── requirements.txt              # Dependencies
├── Dockerfile                    # Docker build file
├── conda/meta.yaml               # Conda recipe
├── LICENSE
└── README.md
```

## License

MIT License
