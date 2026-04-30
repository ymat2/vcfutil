import gzip
import sys
from pathlib import Path


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
                    pop1_daf = calc_derived_allele_frequency(fields, pop1_idx)
                    pop2_daf = calc_derived_allele_frequency(fields, pop2_idx)
                    if (pop1_daf is not None) and (pop2_daf is not None):
                        delta_allele_freqency = pop1_daf - pop2_daf
                        if delta_allele_freqency >= args.extract_daf:
                            print("\t".join([chrom, pos, f"{delta_allele_freqency:.6f}"]))
                        else:
                            continue


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
