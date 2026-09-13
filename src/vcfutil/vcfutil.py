import argparse
from importlib import metadata


__version__ = metadata.version("vcfutil")


def command_daf(args):
    if args.window_size:
        if not args.window_step:
            args.window_step = args.window_size
        from vcfutil.calc_wind_daf import main
        main(args)
    if args.site is not None or args.site2 is not None:
        if args.site is not None and args.site2 is not None:
            raise Exception("Use one of `--site` or `--site2`.")
        else:
            from vcfutil.extract_high_daf import main
            main(args)


def command_dataset(args):
    from vcfutil.dataset import main
    main(args)


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
    help_txt = "Calculate delta allele frequency (DAF) between populations."
    parser_daf = subparsers.add_parser("daf", help=help_txt)
    parser_daf.add_argument("--vcf", help="VCF file to calculate daf. Can be gziped.")
    parser_daf.add_argument("--out", help="Prefix of output file.")
    parser_daf.add_argument("--window_size", type=int, help="Window size in bp.")
    parser_daf.add_argument("--window_step", type=int, help="Window step in bp.")
    parser_daf.add_argument("--site", type=float, default=None, help="Print site level DAF. Pass minimum DAF to be printed.")
    parser_daf.add_argument("--site2", type=float, default=None, help="Similar to `--site`. Supress output of allele counts.")
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
