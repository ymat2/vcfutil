import argparse
from importlib import metadata


__version__ = metadata.version("vcfutil")


def command_hello(name: str) -> None:
    print(f"Hello {name}!")


def command_version() -> None:
    print(f"vcfutil v{__version__}")


def main():
    parser = argparse.ArgumentParser(
        description = "A helper package for VCF manipulation written in python.",
        usage = "vcfutil <commands> [-h/--help] [Options]"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"vcfutil v{__version__}",
        help="Show the version and exit."
    )

    subparsers = parser.add_subparsers(
        dest="command",
        title="Commands",
        metavar=""
    )

    # hello
    help_txt = "Just say hello."
    help_txt += " See `vcfutil help -h`."
    parser_hello = subparsers.add_parser("hello", help=help_txt)
    parser_hello.add_argument("name")
    parser_hello.set_defaults(handler = command_hello)

    # version
    help_txt = "Show the version and exit."
    parser_version = subparsers.add_parser("version", help=help_txt)
    parser_version.set_defaults(handler = command_version)

    args = parser.parse_args()
    if hasattr(args, "handler"):
        args.handler(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
