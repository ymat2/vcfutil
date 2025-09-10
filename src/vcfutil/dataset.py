import random
from datetime import date

GTs = [
    "0/0",
    "0|0",
    "0/1",
    "0|1",
    "1/0",
    "1|0",
    "1/1",
    "1|1",
    "./.",
    ".|."
]

hdr = [
    "#CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT"
]


def main(args):
    lines = []
    today = date.today()
    lines.append('##fileformat=VCFv4.2')
    lines.append(f'##fileDate={today:%Y%m%d}')
    lines.append('##phasing=partial')
    lines.append('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">')
    lines.append('##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele Count">')
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    lines.append('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">')

    for pop in range(args.npop):
        for sample in range(args.nsamples_per_pop):
            hdr.append(f"pop{pop+1}_{sample+1}")
    lines.append("\t".join(hdr))

    for chrom in range(args.nchrom):
        sites = random.sample(range(100, 1000), k=args.nsites_per_chrom)
        sites = sorted(sites)
        for site in sites:
            lines.append(generate_record(chrom+1, site, args.nsamples_per_pop*args.npop))

    with open(args.out, "w") as f:
        for line in lines:
            f.write(line+"\n")


def generate_record(chrom: int, site: int, nsample: int) -> str:
    chrom = f"chr{chrom}"
    pos = str(site)
    id = "."
    ref, alt = random.sample(["A", "T", "G", "C"], k=2)
    qual = str(random.sample(range(20, 255), k=1)[0])
    filter = "PASS"
    format = "GT:DP"
    gts = random.choices(GTs, k=nsample)
    dps = random.choices(range(30, 100), k=nsample)
    info_dp = str(sum(dps))
    info_ac = str("".join(gts).count("1"))
    info = f"AC={info_ac};DP={info_dp}"
    recs = [gts[i]+":"+str(dps[i]) for i in range(nsample)]
    recs = [chrom, pos, id, ref, alt, qual, filter, info, format] + recs
    return "\t".join(recs)


if __name__ == "__main__":
    main()
