import argparse
from importlib import metadata


__version__ = metadata.version("vcfutil")


def command_daf(args):
    if args.window_size and args.extract_daf:
        raise Exception("`--window_size` and `--extract_daf` cannot be used together.")
    elif args.window_size:
        from vcfutil.calc_wind_daf import main
        main(args)
    elif args.extract_daf:
        from vcfutil.extract_high_daf import main
        main(args)


def command_dataset(args):
    from vcfutil.dataset import main
    main(args)


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

    # daf
    help_txt = "Calculate delta allele frequency between populations."
    parser_daf = subparsers.add_parser("daf", help=help_txt)
    parser_daf.add_argument("--vcf", help="VCF file to calculate daf. Can be gziped.")
    #parser_daf.add_argument("--out", help="Prefix of output file.")
    parser_daf.add_argument("--chrom", help="Restrict analysis in a selected chromosome.")
    parser_daf.add_argument("--window_size", type=int, help="Number of SNPs within one window. Cannot be collocated with `--extract`.")
    parser_daf.add_argument("--extract_daf", type=float, help="Minimum delta allele frequenciy to be printed. Cannot be collocated with `--window_size`.")
    parser_daf.add_argument("--stats", choices=["mean", "median"], default="median", help="Statistics to calculate window-based frequencies.")
    parser_daf.add_argument("--pop1", help="Text file containing the target population. One sample for one line.")
    parser_daf.add_argument("--pop2", help="Text file containing the background population. One sample for one line.")
    parser_daf.set_defaults(handler = command_daf)

    # dataset
    help_txt = "Genetrate simple VCF for example."
    parser_data = subparsers.add_parser("dataset", help=help_txt)
    parser_data.add_argument("--out", help="Path for generated VCF.")
    parser_data.add_argument("--nchrom", type=int, help="Number of chromosomes. 2 by default.", default=2)
    parser_data.add_argument("--nsites_per_chrom", type=int, help="Number of sites per chromosome. 6 by default.", default=10)
    parser_data.add_argument("--nsamples_per_pop", type=int, help="Number of samples per population. 3 by default.", default=3)
    parser_data.add_argument("--npop", type=int, help="Number of populations. 2 by default.", default=2)
    parser_data.set_defaults(handler = command_dataset)

    # hello
    help_txt = "Just say hello."
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
