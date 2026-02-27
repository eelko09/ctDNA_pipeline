#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd
from pandas.errors import EmptyDataError


def as_bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def load_table(path):
    try:
        return pd.read_csv(path, sep="\t")
    except EmptyDataError:
        return pd.DataFrame()


def add_variant_id(df):
    required = {"CHROM", "POS", "REF", "ALT"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns {sorted(missing)}")
    df["variant_id"] = (
        df["CHROM"].astype(str)
        + ":"
        + df["POS"].astype(str)
        + ":"
        + df["REF"].astype(str)
        + ":"
        + df["ALT"].astype(str)
    )
    return df


def main():
    parser = argparse.ArgumentParser(description="Apply PBMC blacklist to clinical variant table.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--blacklist", required=True)
    parser.add_argument("--enabled", required=True)
    parser.add_argument("--fail-on-match", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    enabled = as_bool(args.enabled)
    fail_on_match = as_bool(args.fail_on_match)

    input_df = load_table(args.input)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if input_df.empty:
        input_df.to_csv(out_path, sep="\t", index=False)
        return

    input_df = add_variant_id(input_df)

    if not enabled:
        input_df.drop(columns=["variant_id"]).to_csv(out_path, sep="\t", index=False)
        return

    blacklist_df = load_table(args.blacklist)
    if blacklist_df.empty:
        input_df["pbmc_blacklist_match"] = False
        input_df.drop(columns=["variant_id"]).to_csv(out_path, sep="\t", index=False)
        return

    blacklist_df = add_variant_id(blacklist_df)
    blocked = set(blacklist_df["variant_id"].tolist())
    input_df["pbmc_blacklist_match"] = input_df["variant_id"].isin(blocked)

    if fail_on_match:
        input_df = input_df[~input_df["pbmc_blacklist_match"]].copy()

    input_df.drop(columns=["variant_id"]).to_csv(out_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
