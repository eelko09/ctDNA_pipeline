#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

import pandas as pd


def as_bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def load_calls(path, max_vaf):
    df = pd.read_csv(path, sep="\t")
    required = {"CHROM", "POS", "REF", "ALT"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")

    vaf_col = None
    for candidate in ("AF", "VAF", "af", "vaf"):
        if candidate in df.columns:
            vaf_col = candidate
            break
    if vaf_col is not None:
        df[vaf_col] = pd.to_numeric(df[vaf_col], errors="coerce")
        df = df[df[vaf_col].notna() & (df[vaf_col] <= max_vaf)]

    if df.empty:
        return set()

    return set(
        df["CHROM"].astype(str)
        + ":"
        + df["POS"].astype(str)
        + ":"
        + df["REF"].astype(str)
        + ":"
        + df["ALT"].astype(str)
    )


def split_variant_id(variant_id):
    chrom, pos, ref, alt = variant_id.split(":", 3)
    return chrom, pos, ref, alt


def main():
    parser = argparse.ArgumentParser(description="Build pooled PBMC blacklist from per-sample TSV calls.")
    parser.add_argument("--enabled", required=True)
    parser.add_argument("--calls-dir", default="")
    parser.add_argument("--max-vaf", required=True, type=float)
    parser.add_argument("--min-recurrence", required=True, type=int)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    enabled = as_bool(args.enabled)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not enabled:
        pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT", "RECURRENCE", "SAMPLES"]).to_csv(
            out_path, sep="\t", index=False
        )
        return

    calls_dir = Path(args.calls_dir)
    if not calls_dir.exists() or not calls_dir.is_dir():
        raise FileNotFoundError(f"PBMC calls directory not found: {calls_dir}")

    files = sorted(calls_dir.glob("*.tsv"))
    counts = {}
    samples = {}
    for file_path in files:
        sample_id = file_path.stem
        variant_ids = load_calls(file_path, args.max_vaf)
        for variant_id in variant_ids:
            counts[variant_id] = counts.get(variant_id, 0) + 1
            samples.setdefault(variant_id, []).append(sample_id)

    rows = []
    for variant_id, recurrence in counts.items():
        if recurrence < args.min_recurrence:
            continue
        chrom, pos, ref, alt = split_variant_id(variant_id)
        rows.append(
            {
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "RECURRENCE": recurrence,
                "SAMPLES": ",".join(sorted(samples.get(variant_id, []))),
            }
        )

    output_df = pd.DataFrame(rows, columns=["CHROM", "POS", "REF", "ALT", "RECURRENCE", "SAMPLES"])
    output_df.to_csv(out_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
