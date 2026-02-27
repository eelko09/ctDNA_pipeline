#!/usr/bin/env python3
import argparse
import gzip


def parse_info_field(info):
    mapping = {}
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            mapping[key] = value
        else:
            mapping[item] = True
    return mapping


def select_ann_for_alt(ann_entries, alt):
    if not ann_entries:
        return "", "", ""
    matching = []
    for ann in ann_entries:
        fields = ann.split("|")
        if not fields:
            continue
        allele = fields[0]
        if allele == alt:
            matching.append(fields)
    selected = matching[0] if matching else ann_entries[0].split("|")
    effect = selected[1] if len(selected) > 1 else ""
    impact = selected[2] if len(selected) > 2 else ""
    gene = selected[3] if len(selected) > 3 else ""
    return effect, impact, gene


def main():
    parser = argparse.ArgumentParser(description="Extract SnpEff ANN fields into TSV.")
    parser.add_argument("--vcf", required=True, help="Input VCF/VCF.GZ with ANN field")
    parser.add_argument("--out", required=True, help="Output TSV path")
    args = parser.parse_args()

    with gzip.open(args.vcf, "rt") if args.vcf.endswith(".gz") else open(args.vcf, "r") as src, open(args.out, "w") as out:
        out.write("CHROM\tPOS\tREF\tALT\tSNPEFF_EFFECT\tSNPEFF_IMPACT\tSNPEFF_GENE\n")
        for line in src:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _vid, ref, alts, _qual, _flt, info = parts[:8]
            info_map = parse_info_field(info)
            ann_raw = info_map.get("ANN", "")
            ann_entries = [x for x in ann_raw.split(",") if x] if ann_raw else []
            for alt in alts.split(","):
                effect, impact, gene = select_ann_for_alt(ann_entries, alt)
                out.write(
                    f"{chrom}\t{pos}\t{ref}\t{alt}\t{effect}\t{impact}\t{gene}\n"
                )


if __name__ == "__main__":
    main()
