#!/usr/bin/env python3
import argparse
import os

import pandas as pd
from pandas.errors import EmptyDataError


def as_bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def load_tsv(path):
    try:
        return pd.read_csv(path, sep="\t")
    except EmptyDataError:
        return pd.DataFrame()


def variant_ids(df):
    required = {"CHROM", "POS", "REF", "ALT"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns {sorted(missing)}")
    return set(
        df["CHROM"].astype(str)
        + ":"
        + df["POS"].astype(str)
        + ":"
        + df["REF"].astype(str)
        + ":"
        + df["ALT"].astype(str)
    )


def main():
    parser = argparse.ArgumentParser(
        description="Optional tumor-informed filter against known tissue variants."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--enabled", required=True)
    parser.add_argument("--known-dir", default="")
    parser.add_argument("--require-known", required=True)
    parser.add_argument("--fail-on-missing-known", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    enabled = as_bool(args.enabled)
    require_known = as_bool(args.require_known)
    fail_on_missing_known = as_bool(args.fail_on_missing_known)

    df = load_tsv(args.input)
    if df.empty:
        df.to_csv(args.out, sep="\t", index=False)
        return

    if not enabled:
        df["tumor_informed_match"] = False
        df.to_csv(args.out, sep="\t", index=False)
        return

    known_path = os.path.join(args.known_dir, f"{args.sample}.tsv")
    if not os.path.exists(known_path):
        if fail_on_missing_known:
            raise FileNotFoundError(
                f"Tumor-informed mode enabled but known-variant file missing: {known_path}"
            )
        df["tumor_informed_match"] = False
        df.to_csv(args.out, sep="\t", index=False)
        return

    known_df = load_tsv(known_path)
    if known_df.empty:
        known_set = set()
    else:
        known_set = variant_ids(known_df)

    # Recompute row-wise IDs for filtering.
    row_ids = (
        df["CHROM"].astype(str)
        + ":"
        + df["POS"].astype(str)
        + ":"
        + df["REF"].astype(str)
        + ":"
        + df["ALT"].astype(str)
    )
    df["tumor_informed_match"] = row_ids.isin(known_set)

    if require_known:
        df = df[df["tumor_informed_match"]].copy()

    df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
