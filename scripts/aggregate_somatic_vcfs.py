"""Aggregate per-sample Mutect2 VCFs into a long + matrix TSV.

This is intentionally simple and robust.
For real-world use, you may want to parse INFO/FORMAT more deeply.
"""

import argparse
import os
import pandas as pd


def parse_vcf(vcf_path, sample):
    rows = []
    opener = open
    if vcf_path.endswith(".gz"):
        import gzip
        opener = gzip.open

    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("
").split("	")
            chrom, pos, vid, ref, alt, qual, flt, info = parts[:8]
            rows.append(
                {
                    "sample": sample,
                    "chrom": chrom,
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "filter": flt,
                    "qual": qual,
                    "info": info,
                    "variant_id": f"{chrom}:{pos}:{ref}>{alt}",
                }
            )

    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcfs", nargs="+", required=True)
    ap.add_argument("--out-long", required=True)
    ap.add_argument("--out-matrix", required=True)
    args = ap.parse_args()

    dfs = []
    for vcf in args.vcfs:
        sample = os.path.basename(vcf).replace(".filtered.vcf.gz", "")
        df = parse_vcf(vcf, sample)
        dfs.append(df)

    long_df = pd.concat(dfs, axis=0, ignore_index=True) if dfs else pd.DataFrame()

    os.makedirs(os.path.dirname(args.out_long), exist_ok=True)
    long_df.to_csv(args.out_long, sep="	", index=False)

    if long_df.shape[0] == 0:
        pd.DataFrame().to_csv(args.out_matrix, sep="	", index=False)
        return

    # Matrix: 1 if variant present and PASS, else 0
    long_df["is_pass"] = (long_df["filter"] == "PASS").astype(int)

    matrix = (
        long_df.pivot_table(
            index="variant_id",
            columns="sample",
            values="is_pass",
            aggfunc="max",
            fill_value=0,
        )
        .reset_index()
    )

    matrix.to_csv(args.out_matrix, sep="	", index=False)


if __name__ == "__main__":
    main()