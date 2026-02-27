#!/usr/bin/env python
import argparse
import os
import re
import pandas as pd


def parse_flagstat_mapped_pct(path):
    if not os.path.exists(path):
        return None
    pattern = re.compile(r"\(([\d.]+)%\s*:\s*N/A\)")
    with open(path) as handle:
        for line in handle:
            if " mapped (" in line:
                match = pattern.search(line)
                if match:
                    return float(match.group(1))
    return None


def parse_dup_fraction(path):
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[0] == "LIBRARY":
                continue
            try:
                return float(parts[8])
            except ValueError:
                return None
    return None


def parse_mean_coverage(path):
    if not os.path.exists(path):
        return None
    cov_df = pd.read_csv(path, sep="\t")
    if "mean" in cov_df.columns and not cov_df.empty:
        return float(cov_df.loc[0, "mean"])
    return None


def parse_contamination(path):
    if not os.path.exists(path):
        return None
    cont_df = pd.read_csv(path, sep="\t")
    if cont_df.empty or cont_df.shape[1] < 2:
        return None
    try:
        return float(cont_df.iloc[0, 1])
    except ValueError:
        return None


def pass_if_present(value, comparator):
    if value is None:
        return False
    return comparator(value)


def main():
    parser = argparse.ArgumentParser(description="Evaluate ctDNA sample QC gates.")
    parser.add_argument("--samples", required=True, help="Comma-separated sample list.")
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--min-mapped-pct", required=True, type=float)
    parser.add_argument("--min-mean-coverage", required=True, type=float)
    parser.add_argument("--max-dup-fraction", required=True, type=float)
    parser.add_argument("--max-contamination", required=True, type=float)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    samples = [s for s in args.samples.split(",") if s]
    records = []

    for sample in samples:
        flagstat = os.path.join(args.results_dir, "qc", sample, f"{sample}.flagstat.txt")
        dup_metrics = os.path.join(args.results_dir, "qc", sample, f"{sample}.dup_metrics.txt")
        mosdepth = os.path.join(args.results_dir, "coverage", sample, f"{sample}.mosdepth.summary.txt")
        contamination = os.path.join(args.results_dir, "mutect2", f"{sample}.contamination.table")

        mapped_pct = parse_flagstat_mapped_pct(flagstat)
        dup_fraction = parse_dup_fraction(dup_metrics)
        mean_coverage = parse_mean_coverage(mosdepth)
        contam = parse_contamination(contamination)

        mapped_pass = pass_if_present(mapped_pct, lambda x: x >= args.min_mapped_pct)
        coverage_pass = pass_if_present(mean_coverage, lambda x: x >= args.min_mean_coverage)
        dup_pass = pass_if_present(dup_fraction, lambda x: x <= args.max_dup_fraction)
        contam_pass = pass_if_present(contam, lambda x: x <= args.max_contamination)
        overall_pass = mapped_pass and coverage_pass and dup_pass and contam_pass

        records.append(
            {
                "sample": sample,
                "mapped_pct": mapped_pct,
                "mapped_gate": mapped_pass,
                "mean_coverage": mean_coverage,
                "coverage_gate": coverage_pass,
                "dup_fraction": dup_fraction,
                "dup_gate": dup_pass,
                "contamination": contam,
                "contamination_gate": contam_pass,
                "qc_pass": overall_pass,
            }
        )

    out_df = pd.DataFrame(records)
    out_df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
