import gzip
import sys
from . import vcf


def main(args):
    pop1 = vcf.get_samples_from_file(args.pop1)
    pop2 = vcf.get_samples_from_file(args.pop2)
    pop1_idx = vcf.extract_sample_index(args.vcf, pop1)
    pop2_idx = vcf.extract_sample_index(args.vcf, pop2)

    hdr = ["CHROM", "pos", "daf"]
    print("\t".join(hdr))

    handler = gzip.open if args.vcf.lower().endswith(".gz") else open
    with handler(args.vcf, "rt") as f:
        tmp_chrom = ""
        for line in f:
            if line.startswith("#"):
                continue
            else:
                fields = line.rstrip("\n").split("\t")
                chrom, pos = fields[0], fields[1]
                if chrom != tmp_chrom:
                    tmp_chrom = chrom
                    print(f"\tAnalyzing {tmp_chrom}", file=sys.stderr)
                else:
                    pop1_daf = vcf.calc_derived_allele_frequency(fields, pop1_idx)
                    pop2_daf = vcf.calc_derived_allele_frequency(fields, pop2_idx)
                    if (pop1_daf is not None) and (pop2_daf is not None):
                        delta_allele_freqency = pop1_daf - pop2_daf
                        if delta_allele_freqency >= args.extract_daf:
                            print("\t".join([chrom, pos, f"{delta_allele_freqency:.6f}"]))
                        else:
                            continue


if __name__ == "__main__":
    main()
