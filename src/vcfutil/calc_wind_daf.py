import gzip
import sys
from math import ceil
from . import vcf


def main(args):
    pop1 = vcf.get_samples_from_file(args.pop1)
    pop2 = vcf.get_samples_from_file(args.pop2)
    pop1_idx = vcf.extract_sample_index(args.vcf, pop1)
    pop2_idx = vcf.extract_sample_index(args.vcf, pop2)

    chrom2windows = {}  # chrom2windows[chrom] = {[0]: [number_of_SNPs, sum_of_daf], [1]: ...}
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
                    chrom2windows[tmp_chrom] = {}
                    first_idx = ceil((int(pos) - args.window_size) / args.window_step) if int(pos) > args.window_size else 0
                    last_idx = ceil(int(pos) / args.window_step)
                    pop1_daf = vcf.calc_derived_allele_frequency(fields, pop1_idx)
                    pop2_daf = vcf.calc_derived_allele_frequency(fields, pop2_idx)
                    if (pop1_daf is not None) and (pop2_daf is not None):
                        for idx in range(first_idx, last_idx):
                            if idx not in chrom2windows[tmp_chrom].keys():
                                chrom2windows[tmp_chrom][idx] = [0, 0]  # [number of SNPs, sum of daf]
                            chrom2windows[tmp_chrom][idx][0] += 1
                            chrom2windows[tmp_chrom][idx][1] += abs(pop1_daf - pop2_daf)
                else:
                    first_idx = ceil((int(pos) - args.window_size) / args.window_step) if int(pos) > args.window_size else 0
                    last_idx = ceil(int(pos) / args.window_step)
                    pop1_daf = vcf.calc_derived_allele_frequency(fields, pop1_idx)
                    pop2_daf = vcf.calc_derived_allele_frequency(fields, pop2_idx)
                    if (pop1_daf is not None) and (pop2_daf is not None):
                        for idx in range(first_idx, last_idx):
                            if idx not in chrom2windows[tmp_chrom].keys():
                                chrom2windows[tmp_chrom][idx] = [0, 0]  # [number of SNPs, sum of daf]
                            chrom2windows[tmp_chrom][idx][0] += 1
                            chrom2windows[tmp_chrom][idx][1] += abs(pop1_daf - pop2_daf)

    outfile = args.out + ".windowed.daf"
    with open(outfile, "w") as f:
        hdr = ["CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "MEAN_DAF"]
        f.write("\t".join(hdr)+"\n")
        for chr in chrom2windows:
            for idx,win in chrom2windows[chr].items():
                bin_start = args.window_step * int(idx) + 1
                bin_end = bin_start + args.window_size - 1
                n_variants = win[0]
                if args.window_mean:
                    mean_daf = f"{win[1]/args.window_size:.3f}"
                else:
                    mean_daf = f"{win[1]/win[0]:.3f}"
                f.write("\t".join([chr, str(bin_start), str(bin_end), str(n_variants), str(mean_daf)])+"\n")


if __name__ == "__main__":
    main()
