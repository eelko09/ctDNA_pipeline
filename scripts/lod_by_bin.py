#!/usr/bin/env python3
import argparse
import json
import math
import os
import pandas as pd


def parse_mean_coverage(path):
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t")
    if df.empty or "mean" not in df.columns:
        return None
    return float(df.loc[0, "mean"])


def parse_contamination(path):
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t")
    if df.empty or df.shape[1] < 2:
        return None
    try:
        return float(df.iloc[0, 1])
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser(description="Generate sample-level LOD/callable bin summary.")
    parser.add_argument("--samples", required=True, help="Comma-separated sample list")
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--bins-json", required=True, help="JSON list of bins with name/min_af/max_af")
    parser.add_argument("--min-alt-reads", required=True, type=float)
    parser.add_argument("--max-contamination", required=True, type=float)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    bins = json.loads(args.bins_json)
    samples = [s for s in args.samples.split(",") if s]
    records = []

    for sample in samples:
        cov_path = os.path.join(
            args.results_dir, "coverage", sample, f"{sample}.mosdepth.summary.txt"
        )
        contam_path = os.path.join(
            args.results_dir, "mutect2", f"{sample}.contamination.table"
        )
        mean_cov = parse_mean_coverage(cov_path)
        contamination = parse_contamination(contam_path)
        contam_ok = (
            contamination is not None and contamination <= args.max_contamination
        )

        for bin_cfg in bins:
            bin_name = str(bin_cfg["name"])
            min_af = float(bin_cfg["min_af"])
            max_af = float(bin_cfg["max_af"])
            required_depth = math.ceil(args.min_alt_reads / max(min_af, 1e-9))
            callable_flag = (
                mean_cov is not None and mean_cov >= required_depth and contam_ok
            )
            records.append(
                {
                    "sample": sample,
                    "bin_name": bin_name,
                    "min_af": min_af,
                    "max_af": max_af,
                    "mean_coverage": mean_cov,
                    "required_depth": required_depth,
                    "contamination": contamination,
                    "contamination_gate": contam_ok,
                    "callable": callable_flag,
                }
            )

    pd.DataFrame(records).to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
