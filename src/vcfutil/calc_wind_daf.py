import gzip
import sys
from pathlib import Path
from statistics import mean, median


GT2AC = {
    "0/0":0,
    "0|0":0,
    "0/1":1,
    "0|1":1,
    "1/0":1,
    "1|0":1,
    "1/1":2,
    "1|1":2
}


def main(args):
    pop1 = get_samples_from_file(args.pop1)
    pop2 = get_samples_from_file(args.pop2)
    pop1_idx = extract_sample_index(args.vcf, pop1)
    pop2_idx = extract_sample_index(args.vcf, pop2)

    hdr = ["CHROM", "from", "to", f"{args.stats}_daf"]
    print("\t".join(hdr))

    handler = gzip.open if args.vcf.lower().endswith(".gz") else open
    with handler(args.vcf, "rt") as f:
        snp_num = 0
        tmp_chrom = ""
        pos_in_window = []
        daf_in_window = []
        for line in f:
            if line.startswith("#"):
                continue
            else:
                fields = line.rstrip("\n").split("\t")
                chrom, pos = fields[0], fields[1]
                if args.chrom:
                    if chrom == args.chrom:
                        pop1_daf = calc_derived_allele_frequency(fields, pop1_idx)
                        pop2_daf = calc_derived_allele_frequency(fields, pop2_idx)
                        if (pop1_daf is not None) and (pop2_daf is not None):
                            snp_num += 1
                            pos_in_window.append(pos)
                            daf_in_window.append(abs(pop1_daf - pop2_daf))
                            if snp_num >= args.window_size:
                                window_daf = median(daf_in_window) if args.stats == "median" else mean(daf_in_window)
                                print("\t".join([chrom, pos_in_window[0], pos_in_window[-1], f"{window_daf:.6f}"]))
                                snp_num = 0
                                pos_in_window.clear()
                                daf_in_window.clear()
                    else:
                        continue
                else:
                    if chrom != tmp_chrom:
                        tmp_chrom = chrom
                        print(f"\tAnalyzing {tmp_chrom}", file=sys.stderr)
                        snp_num = 0
                        pos_in_window.clear()
                        daf_in_window.clear()

                        pop1_daf = calc_derived_allele_frequency(fields, pop1_idx)
                        pop2_daf = calc_derived_allele_frequency(fields, pop2_idx)
                        if (pop1_daf is not None) and (pop2_daf is not None):
                            snp_num += 1
                            pos_in_window.append(pos)
                            daf_in_window.append(abs(pop1_daf - pop2_daf))
                            #daf_in_window.append(pop1_daf)
                    else:
                        pop1_daf = calc_derived_allele_frequency(fields, pop1_idx)
                        pop2_daf = calc_derived_allele_frequency(fields, pop2_idx)
                        if (pop1_daf is not None) and (pop2_daf is not None):
                            snp_num += 1
                            pos_in_window.append(pos)
                            daf_in_window.append(abs(pop1_daf - pop2_daf))
                            #daf_in_window.append(pop1_daf)
                            if snp_num >= args.window_size:
                                window_daf = median(daf_in_window) if args.stats == "median" else mean(daf_in_window)
                                print("\t".join([chrom, pos_in_window[0], pos_in_window[-1], f"{window_daf:.6f}"]))
                                snp_num = 0
                                pos_in_window.clear()
                                daf_in_window.clear()


def get_samples_from_file(file: Path) -> list:
    samples = []
    with open(file) as f:
        samples = [s.rstrip("\n") for s in f.readlines()]
    return samples


def extract_sample_index(vcf_file: Path, pop: list) -> list:
    sample_names = []
    sample_idx = []
    handler = gzip.open if vcf_file.lower().endswith(".gz") else open
    with handler(vcf_file, "rt") as vcf:
        for line in vcf:
            line = line.strip("\n")
            if line.startswith("#CHROM"):
                record = line.split("\t")
                sample_names = [record[i] for i in range(9, len(record))]
                break
    sample_idx = [i for i, sample in enumerate(sample_names) if sample in pop]
    return sample_idx


def get_genotypes(samples: list) -> list:
    """
    Takes each row of the VCF.
    Returns a list of GTs.
    """
    gts = []
    gts = [samples[i].split(":")[0] for i in range(9, len(samples))]
    return gts


def calc_derived_allele_frequency(sample_field: list, pop_idx: list) -> float | None:
    daf = 0
    gts = get_genotypes(sample_field)
    gts_of_pop = [gts[i] for i in pop_idx if gts[i] not in {"./.", ".|."}]
    if len(gts_of_pop) == 0:
        daf = None
    else:
        daf = sum([GT2AC[gt] for gt in gts_of_pop]) / (2*len(gts_of_pop))
    return daf


if __name__ == "__main__":
    main()
