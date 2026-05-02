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
