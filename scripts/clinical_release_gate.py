#!/usr/bin/env python3
import argparse
import json
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


def main():
    parser = argparse.ArgumentParser(description="Generate per-sample clinical release gating.")
    parser.add_argument("--samples", required=True)
    parser.add_argument("--qc-gates", required=True)
    parser.add_argument("--lod", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--variant-tables", required=True)
    parser.add_argument("--require-qc-pass", required=True)
    parser.add_argument("--require-manifest-git-sha", required=True)
    parser.add_argument("--require-variants", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    samples = [s for s in args.samples.split(",") if s]
    require_qc_pass = as_bool(args.require_qc_pass)
    require_manifest_git_sha = as_bool(args.require_manifest_git_sha)
    require_variants = as_bool(args.require_variants)
    variant_table_paths = [p for p in args.variant_tables.split(",") if p]

    qc_df = load_tsv(args.qc_gates)
    qc_pass_map = {}
    if not qc_df.empty and {"sample", "qc_pass"}.issubset(qc_df.columns):
        qc_pass_map = dict(zip(qc_df["sample"].astype(str), qc_df["qc_pass"].astype(bool)))

    lod_ok = os.path.exists(args.lod) and os.path.getsize(args.lod) > 0

    manifest_ok = False
    if os.path.exists(args.manifest):
        with open(args.manifest) as handle:
            manifest = json.load(handle)
        if require_manifest_git_sha:
            git_sha = str(manifest.get("git_sha", "")).strip().lower()
            manifest_ok = git_sha not in {"", "unknown", "none", "null"}
        else:
            manifest_ok = True

    table_by_sample = {}
    for path in variant_table_paths:
        sample = os.path.basename(path).split(".")[0]
        table_by_sample[sample] = path

    rows = []
    for sample in samples:
        reasons = []

        qc_pass = bool(qc_pass_map.get(sample, False))
        if require_qc_pass and not qc_pass:
            reasons.append("qc_fail")

        if not lod_ok:
            reasons.append("lod_missing")

        if not manifest_ok:
            reasons.append("manifest_missing_or_unversioned")

        variant_count = 0
        vt_path = table_by_sample.get(sample, "")
        if vt_path and os.path.exists(vt_path):
            vt_df = load_tsv(vt_path)
            variant_count = 0 if vt_df.empty else int(vt_df.shape[0])
        else:
            reasons.append("variant_table_missing")

        if require_variants and variant_count == 0:
            reasons.append("no_variants")

        release_status = "PASS" if not reasons else "FAIL"
        rows.append(
            {
                "sample": sample,
                "qc_pass": qc_pass,
                "lod_report_present": lod_ok,
                "manifest_ok": manifest_ok,
                "variant_count": variant_count,
                "release_status": release_status,
                "release_reasons": ",".join(reasons),
            }
        )

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
