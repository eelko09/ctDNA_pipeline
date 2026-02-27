#!/usr/bin/env python3
import argparse
import os
import pandas as pd
from pandas.errors import EmptyDataError


def as_bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def load_variant_keys(path):
    try:
        df = pd.read_csv(path, sep="\t")
    except EmptyDataError:
        cols = [
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "DP",
            "AF",
        ]
        return pd.DataFrame(columns=cols + ["variant_id"])

    # Backward compatibility: older outputs without headers.
    if not {"CHROM", "POS", "REF", "ALT"}.issubset(df.columns) and df.shape[1] >= 8:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "DP", "AF"],
        )

    required = {"CHROM", "POS", "REF", "ALT"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Variant table missing required columns {sorted(missing)}: {path}")
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


def load_orthogonal_calls(path):
    calls = pd.read_csv(path, sep="\t")
    required = {"CHROM", "POS", "REF", "ALT"}
    missing = required - set(calls.columns)
    if missing:
        raise ValueError(
            f"Orthogonal calls missing required columns {sorted(missing)}: {path}"
        )
    return set(
        calls["CHROM"].astype(str)
        + ":"
        + calls["POS"].astype(str)
        + ":"
        + calls["REF"].astype(str)
        + ":"
        + calls["ALT"].astype(str)
    )


def load_chip_panel(path):
    panel = pd.read_csv(path, sep="\t")
    required = {"CHROM", "POS", "REF", "ALT", "GENE"}
    missing = required - set(panel.columns)
    if missing:
        raise ValueError(f"CHIP panel missing required columns {sorted(missing)}: {path}")
    panel["variant_id"] = (
        panel["CHROM"].astype(str)
        + ":"
        + panel["POS"].astype(str)
        + ":"
        + panel["REF"].astype(str)
        + ":"
        + panel["ALT"].astype(str)
    )
    gene_map = {}
    for _, row in panel.iterrows():
        gene_map[row["variant_id"]] = str(row["GENE"])
    return gene_map


def load_annotation_panel(path):
    panel = pd.read_csv(path, sep="\t")
    required = {"CHROM", "POS", "REF", "ALT", "GENE", "CLINICAL_TIER", "ACTIONABILITY"}
    missing = required - set(panel.columns)
    if missing:
        raise ValueError(
            f"Clinical annotation panel missing required columns {sorted(missing)}: {path}"
        )
    panel["variant_id"] = (
        panel["CHROM"].astype(str)
        + ":"
        + panel["POS"].astype(str)
        + ":"
        + panel["REF"].astype(str)
        + ":"
        + panel["ALT"].astype(str)
    )
    ann = {}
    for _, row in panel.iterrows():
        ann[row["variant_id"]] = {
            "clinical_gene": str(row["GENE"]),
            "clinical_tier": str(row["CLINICAL_TIER"]),
            "actionability": str(row["ACTIONABILITY"]),
        }
    return ann


def load_clinvar_cosmic(path):
    panel = pd.read_csv(path, sep="\t")
    required = {"CHROM", "POS", "REF", "ALT", "CLINVAR", "COSMIC"}
    missing = required - set(panel.columns)
    if missing:
        raise ValueError(
            f"ClinVar/COSMIC table missing required columns {sorted(missing)}: {path}"
        )
    panel["variant_id"] = (
        panel["CHROM"].astype(str)
        + ":"
        + panel["POS"].astype(str)
        + ":"
        + panel["REF"].astype(str)
        + ":"
        + panel["ALT"].astype(str)
    )
    ann = {}
    for _, row in panel.iterrows():
        ann[row["variant_id"]] = {
            "clinvar": str(row["CLINVAR"]),
            "cosmic": str(row["COSMIC"]),
        }
    return ann


def load_simple_varset(path):
    if not path or not os.path.exists(path):
        return set()
    calls = pd.read_csv(path, sep="\t")
    required = {"CHROM", "POS", "REF", "ALT"}
    missing = required - set(calls.columns)
    if missing:
        raise ValueError(f"Missing {sorted(missing)} in {path}")
    return set(
        calls["CHROM"].astype(str)
        + ":"
        + calls["POS"].astype(str)
        + ":"
        + calls["REF"].astype(str)
        + ":"
        + calls["ALT"].astype(str)
    )


def load_snpeff_map(path):
    if not path or not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    required = {"CHROM", "POS", "REF", "ALT", "SNPEFF_EFFECT", "SNPEFF_IMPACT", "SNPEFF_GENE"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"SnpEff table missing required columns {sorted(missing)}: {path}")
    df["variant_id"] = (
        df["CHROM"].astype(str)
        + ":"
        + df["POS"].astype(str)
        + ":"
        + df["REF"].astype(str)
        + ":"
        + df["ALT"].astype(str)
    )
    mapping = {}
    for _, row in df.iterrows():
        mapping[row["variant_id"]] = {
            "snpeff_effect": str(row["SNPEFF_EFFECT"]),
            "snpeff_impact": str(row["SNPEFF_IMPACT"]),
            "snpeff_gene": str(row["SNPEFF_GENE"]),
        }
    return mapping


def main():
    parser = argparse.ArgumentParser(description="Annotate consensus and CHIP flags.")
    parser.add_argument("--input", required=True, help="Mutect-derived variant TSV")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--orth-enabled", required=True)
    parser.add_argument("--orth-calls-dir", default="")
    parser.add_argument("--chip-enabled", required=True)
    parser.add_argument("--chip-panel", default="")
    parser.add_argument("--wbc-enabled", required=True)
    parser.add_argument("--wbc-calls-dir", default="")
    parser.add_argument("--wbc-fail-on-support", required=True)
    parser.add_argument("--normal-sample", default="")
    parser.add_argument("--clinical-annotations-enabled", required=True)
    parser.add_argument("--clinical-annotations-panel", default="")
    parser.add_argument("--clinvar-cosmic-tsv", default="")
    parser.add_argument("--varscan-enabled", required=True)
    parser.add_argument("--varscan-tsv", default="")
    parser.add_argument("--snpeff-enabled", required=True)
    parser.add_argument("--snpeff-tsv", default="")
    parser.add_argument("--min-dp", required=True, type=float)
    parser.add_argument("--min-alt-reads", required=True, type=float)
    parser.add_argument("--min-af", required=True, type=float)
    parser.add_argument("--low-vaf-threshold", required=True, type=float)
    parser.add_argument("--require-orthogonal-low-vaf", required=True)
    parser.add_argument("--chip-flag-action", default="review")
    args = parser.parse_args()

    df = load_variant_keys(args.input)

    orth_enabled = as_bool(args.orth_enabled)
    chip_enabled = as_bool(args.chip_enabled)
    wbc_enabled = as_bool(args.wbc_enabled)
    wbc_fail_on_support = as_bool(args.wbc_fail_on_support)
    clinical_ann_enabled = as_bool(args.clinical_annotations_enabled)
    varscan_enabled = as_bool(args.varscan_enabled)
    snpeff_enabled = as_bool(args.snpeff_enabled)
    require_orthogonal_low_vaf = as_bool(args.require_orthogonal_low_vaf)

    orth_supported = set()
    if orth_enabled:
        orth_path = os.path.join(args.orth_calls_dir, f"{args.sample}.tsv")
        if not os.path.exists(orth_path):
            raise FileNotFoundError(
                f"Orthogonal cross-check enabled but file missing: {orth_path}"
            )
        orth_supported = load_orthogonal_calls(orth_path)

    chip_genes = {}
    if chip_enabled:
        if not os.path.exists(args.chip_panel):
            raise FileNotFoundError(
                f"CHIP panel enabled but file missing: {args.chip_panel}"
            )
        chip_genes = load_chip_panel(args.chip_panel)

    wbc_supported = set()
    if wbc_enabled and args.normal_sample:
        wbc_path = os.path.join(args.wbc_calls_dir, f"{args.normal_sample}.tsv")
        if os.path.exists(wbc_path):
            wbc_supported = load_orthogonal_calls(wbc_path)

    ann_map = {}
    if clinical_ann_enabled:
        if not os.path.exists(args.clinical_annotations_panel):
            raise FileNotFoundError(
                "Clinical annotations enabled but panel missing: "
                f"{args.clinical_annotations_panel}"
            )
        ann_map = load_annotation_panel(args.clinical_annotations_panel)

    clinvar_cosmic = {}
    if args.clinvar_cosmic_tsv and os.path.exists(args.clinvar_cosmic_tsv):
        clinvar_cosmic = load_clinvar_cosmic(args.clinvar_cosmic_tsv)

    varscan_supported = set()
    if varscan_enabled:
        if not os.path.exists(args.varscan_tsv):
            raise FileNotFoundError(f"VarScan enabled but TSV missing: {args.varscan_tsv}")
        varscan_supported = load_simple_varset(args.varscan_tsv)

    snpeff_map = {}
    if snpeff_enabled:
        if not os.path.exists(args.snpeff_tsv):
            raise FileNotFoundError(f"SnpEff enabled but TSV missing: {args.snpeff_tsv}")
        snpeff_map = load_snpeff_map(args.snpeff_tsv)

    if orth_enabled:
        df["orthogonal_support"] = df["variant_id"].isin(orth_supported)
        df["consensus_flag"] = df["orthogonal_support"].map(
            lambda x: "consensus" if x else "mutect_only"
        )
    else:
        df["orthogonal_support"] = False
        df["consensus_flag"] = "not_evaluated"

    if varscan_enabled:
        df["varscan_support"] = df["variant_id"].isin(varscan_supported)
        if orth_enabled:
            df["consensus_flag"] = (
                df["orthogonal_support"] | df["varscan_support"]
            ).map(lambda x: "consensus" if x else "mutect_only")
        else:
            df["consensus_flag"] = df["varscan_support"].map(
                lambda x: "consensus" if x else "mutect_only"
            )
    else:
        df["varscan_support"] = False

    if chip_enabled:
        df["chip_gene"] = df["variant_id"].map(chip_genes).fillna("")
        df["chip_flag"] = df["chip_gene"] != ""
    else:
        df["chip_gene"] = ""
        df["chip_flag"] = False

    if wbc_enabled and args.normal_sample:
        df["matched_wbc_support"] = df["variant_id"].isin(wbc_supported)
    else:
        df["matched_wbc_support"] = False

    if clinical_ann_enabled:
        df["clinical_gene"] = df["variant_id"].map(
            lambda vid: ann_map.get(vid, {}).get("clinical_gene", "")
        )
        df["clinical_tier"] = df["variant_id"].map(
            lambda vid: ann_map.get(vid, {}).get("clinical_tier", "")
        )
        df["actionability"] = df["variant_id"].map(
            lambda vid: ann_map.get(vid, {}).get("actionability", "")
        )
    else:
        df["clinical_gene"] = ""
        df["clinical_tier"] = ""
        df["actionability"] = ""

    if clinvar_cosmic:
        df["clinvar"] = df["variant_id"].map(
            lambda vid: clinvar_cosmic.get(vid, {}).get("clinvar", "")
        )
        df["cosmic"] = df["variant_id"].map(
            lambda vid: clinvar_cosmic.get(vid, {}).get("cosmic", "")
        )
    else:
        df["clinvar"] = ""
        df["cosmic"] = ""

    if snpeff_enabled:
        df["snpeff_effect"] = df["variant_id"].map(
            lambda vid: snpeff_map.get(vid, {}).get("snpeff_effect", "")
        )
        df["snpeff_impact"] = df["variant_id"].map(
            lambda vid: snpeff_map.get(vid, {}).get("snpeff_impact", "")
        )
        df["snpeff_gene"] = df["variant_id"].map(
            lambda vid: snpeff_map.get(vid, {}).get("snpeff_gene", "")
        )
    else:
        df["snpeff_effect"] = ""
        df["snpeff_impact"] = ""
        df["snpeff_gene"] = ""

    # Ensure numeric support fields are available for gating.
    df["DP"] = pd.to_numeric(df.get("DP", 0), errors="coerce").fillna(0.0)
    df["AF"] = pd.to_numeric(df.get("AF", 0), errors="coerce").fillna(0.0)
    df["ALT_COUNT"] = (df["DP"] * df["AF"]).round(0)

    def support_gate(row):
        reasons = []
        if row["DP"] < args.min_dp:
            reasons.append("low_dp")
        if row["AF"] < args.min_af:
            reasons.append("low_af")
        if row["ALT_COUNT"] < args.min_alt_reads:
            reasons.append("low_alt_reads")

        if (
            require_orthogonal_low_vaf
            and orth_enabled
            and row["AF"] < args.low_vaf_threshold
            and not row["orthogonal_support"]
        ):
            reasons.append("needs_orthogonal_support")

        if wbc_enabled and wbc_fail_on_support and row["matched_wbc_support"]:
            reasons.append("matched_wbc")

        chip_action = args.chip_flag_action.lower()
        if row["chip_flag"] and chip_action == "review":
            reasons.append("chip_review")
        if row["chip_flag"] and chip_action == "fail":
            reasons.append("chip_fail")

        if not reasons:
            return "PASS", ""
        if "chip_review" in reasons and len(reasons) == 1:
            return "REVIEW", "chip_review"
        return "FAIL", ",".join(reasons)

    gate = df.apply(support_gate, axis=1, result_type="expand")
    df["support_gate"] = gate[0]
    df["support_reasons"] = gate[1]

    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
