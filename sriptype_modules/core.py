"""Core logic and shared processing for SRIPType."""

from sriptype_modules.utils import setup_logger

logger = setup_logger(__name__)


def validate_threads(threads):
    """Validate and return the number of threads to use."""
    import os
    max_threads = os.cpu_count() or 1
    if threads < 1:
        logger.warning("Thread count must be >= 1, using 1")
        return 1
    if threads > max_threads:
        logger.warning(
            "Requested %d threads but only %d CPUs available, using %d",
            threads, max_threads, max_threads,
        )
        return max_threads
    return threads
