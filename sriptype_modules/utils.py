"""Common utility functions for SRIPType."""

import logging
import os
import sys


def setup_logger(name="sriptype", level=logging.INFO):
    """Set up and return a logger instance."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        formatter = logging.Formatter(
            "[%(asctime)s] %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logger.setLevel(level)
    return logger


def check_file_exists(filepath, label="Input"):
    """Check if a file exists and is readable."""
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"{label} file not found: {filepath}")
    if not os.access(filepath, os.R_OK):
        raise PermissionError(f"{label} file is not readable: {filepath}")
    return os.path.abspath(filepath)


def ensure_output_dir(filepath):
    """Ensure the output directory exists, create if necessary."""
    outdir = os.path.dirname(os.path.abspath(filepath))
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    return os.path.abspath(filepath)
