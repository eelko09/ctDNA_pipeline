#!/usr/bin/env python3
import argparse
import pandas as pd
from pandas.errors import EmptyDataError


def as_bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def main():
    parser = argparse.ArgumentParser(description="Gate clinical output table.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--enabled", required=True)
    parser.add_argument("--accepted-gates", required=True, help="Comma-separated gates")
    parser.add_argument("--include-only-annotated", required=True)
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input, sep="\t")
    except EmptyDataError:
        pd.DataFrame().to_csv(args.output, sep="\t", index=False)
        return
    enabled = as_bool(args.enabled)
    include_only_annotated = as_bool(args.include_only_annotated)
    accepted = {g.strip() for g in args.accepted_gates.split(",") if g.strip()}

    if enabled:
        if "support_gate" in df.columns:
            df = df[df["support_gate"].isin(accepted)].copy()
        if include_only_annotated and "clinical_tier" in df.columns:
            df = df[df["clinical_tier"].astype(str).str.len() > 0].copy()

    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
