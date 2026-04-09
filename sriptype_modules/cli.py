"""CLI argument parsing for SRIPType."""

import argparse
import importlib
import pkgutil
import sys

from sriptype_modules import __version__
import sriptype_modules.subcommands as subcommands_pkg


def discover_subcommands():
    """Dynamically discover all subcommand modules in the subcommands package."""
    subcommands = {}
    package_path = subcommands_pkg.__path__
    for importer, modname, ispkg in pkgutil.iter_modules(package_path):
        if modname.startswith("_"):
            continue
        module = importlib.import_module(f"sriptype_modules.subcommands.{modname}")
        if hasattr(module, "add_arguments") and hasattr(module, "run"):
            subcommands[modname] = module
    return subcommands


def build_parser():
    """Build the main argument parser with all subcommands."""
    parser = argparse.ArgumentParser(
        prog="sriptype",
        description="SRIPType - A bioinformatics tool for 51k-SINE-RIP chip genotyping",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  sriptype mkdb -i input_dir -o blast_db -j 4 -t 4\n"
            "  sriptype genotype -i blast_db -o genotype_output -j 4 -t 4\n"
            "\n"
            "Use 'sriptype <subcommand> -h' for subcommand-specific help.\n"
        ),
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    # Subcommands
    subparsers = parser.add_subparsers(
        title="Available subcommands",
        dest="subcommand",
        metavar="<subcommand>",
    )

    subcommands = discover_subcommands()
    for name, module in sorted(subcommands.items(),
                                key=lambda x: getattr(x[1], "ORDER", 99)):
        desc = getattr(module, "DESCRIPTION", f"{name} subcommand")
        sub_parser = subparsers.add_parser(
            name,
            help=desc,
            description=desc,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        module.add_arguments(sub_parser)
        # Add global arguments after subcommand-specific ones
        global_group = sub_parser.add_argument_group("global options")
        if getattr(module, "USES_THREADS", True):
            global_group.add_argument(
                "-t", "--threads",
                type=int,
                default=1,
                metavar="INT",
                help="Number of threads to use (default: 1)",
            )
        global_group.add_argument(
            "--verbose",
            action="store_true",
            default=False,
            help="Enable verbose/debug output",
        )
        sub_parser.set_defaults(func=module.run)

    return parser, subcommands


def main():
    """Main entry point for the CLI."""
    parser, subcommands = build_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.subcommand is None:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Set up logging level based on --verbose
    if args.verbose:
        import logging
        logging.getLogger().setLevel(logging.DEBUG)
        # Also update all sriptype_modules loggers so their own level
        # check doesn't filter out DEBUG messages
        for name in list(logging.Logger.manager.loggerDict):
            if name.startswith("sriptype_modules"):
                logging.getLogger(name).setLevel(logging.DEBUG)

    try:
        args.func(args)
    except KeyboardInterrupt:
        sys.stderr.write("\nInterrupted by user.\n")
        sys.exit(130)
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)
