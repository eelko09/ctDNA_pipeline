#!/usr/bin/env python3
import sys
from pathlib import Path

import yaml


def fail(message):
    raise ValueError(message)


def require_key(dct, key, where):
    if key not in dct:
        fail(f"Missing key '{key}' in {where}")


def require_positive_number(value, name, allow_zero=False):
    if not isinstance(value, (int, float)):
        fail(f"{name} must be numeric, got {type(value).__name__}")
    if allow_zero:
        if value < 0:
            fail(f"{name} must be >= 0, got {value}")
    else:
        if value <= 0:
            fail(f"{name} must be > 0, got {value}")


def validate_resources(cfg):
    require_key(cfg, "resources", "config")
    resources = cfg["resources"]
    expected_rules = ["fastqc", "fastp", "bwa_mem", "samtools_sort", "gatk", "mutect2"]
    for rule in expected_rules:
        require_key(resources, rule, "resources")
        require_key(resources[rule], "threads", f"resources.{rule}")
        require_key(resources[rule], "mem_mb", f"resources.{rule}")
        require_positive_number(resources[rule]["threads"], f"resources.{rule}.threads")
        require_positive_number(resources[rule]["mem_mb"], f"resources.{rule}.mem_mb")


def validate_references(cfg):
    require_key(cfg, "references", "config")
    refs = cfg["references"]
    for key in ["reference_fasta", "panel_bed", "germline_resource", "common_variants", "pon_vcf"]:
        require_key(refs, key, "references")
    require_key(refs, "known_sites", "references")
    require_key(refs["known_sites"], "dbsnp", "references.known_sites")
    require_key(refs["known_sites"], "mills_indels", "references.known_sites")


def validate_variant_calling(cfg):
    vc = cfg.get("variant_calling", {})
    mode = vc.get("mode", "tumor_only")
    if mode not in {"tumor_only", "tumor_normal", "auto"}:
        fail(f"variant_calling.mode must be one of tumor_only/tumor_normal/auto, got {mode}")

    mutect2 = vc.get("mutect2", {})
    maf = mutect2.get("min_allele_fraction", 0.001)
    require_positive_number(maf, "variant_calling.mutect2.min_allele_fraction")
    if maf >= 1:
        fail("variant_calling.mutect2.min_allele_fraction must be < 1")

    post = vc.get("postfilter", {})
    pass_only = post.get("pass_only", True)
    if not isinstance(pass_only, bool):
        fail("variant_calling.postfilter.pass_only must be boolean")
    require_positive_number(post.get("min_dp", 1), "variant_calling.postfilter.min_dp")
    require_positive_number(post.get("min_alt_reads", 1), "variant_calling.postfilter.min_alt_reads")
    min_af = post.get("min_af", 0.001)
    require_positive_number(min_af, "variant_calling.postfilter.min_af")
    if min_af >= 1:
        fail("variant_calling.postfilter.min_af must be < 1")


def validate_qc_gates(cfg):
    gates = cfg.get("qc_gates", {})
    require_positive_number(gates.get("min_mapped_pct", 1), "qc_gates.min_mapped_pct")
    require_positive_number(gates.get("min_mean_coverage", 1), "qc_gates.min_mean_coverage")
    require_positive_number(gates.get("max_dup_fraction", 1), "qc_gates.max_dup_fraction", allow_zero=True)
    require_positive_number(gates.get("max_contamination", 1), "qc_gates.max_contamination", allow_zero=True)
    if gates.get("min_mapped_pct", 100) > 100:
        fail("qc_gates.min_mapped_pct must be <= 100")
    if gates.get("max_dup_fraction", 1) > 1:
        fail("qc_gates.max_dup_fraction must be <= 1")
    if gates.get("max_contamination", 1) > 1:
        fail("qc_gates.max_contamination must be <= 1")


def validate_assay(cfg):
    assay = cfg.get("assay", {})
    if not assay:
        return

    umi = assay.get("umi", {})
    if umi:
        enabled = umi.get("enabled", False)
        if not isinstance(enabled, bool):
            fail("assay.umi.enabled must be boolean")
        method = umi.get("method", "umi_tools")
        if method not in {"umi_tools"}:
            fail("assay.umi.method currently supports only: umi_tools")
        if "bc_pattern" in umi and not isinstance(umi["bc_pattern"], str):
            fail("assay.umi.bc_pattern must be a string")

    orth = assay.get("orthogonal", {})
    if orth:
        enabled = orth.get("enabled", False)
        if not isinstance(enabled, bool):
            fail("assay.orthogonal.enabled must be boolean")
        if "calls_dir" in orth and not isinstance(orth["calls_dir"], str):
            fail("assay.orthogonal.calls_dir must be a string path")

    chip = assay.get("chip", {})
    if chip:
        enabled = chip.get("enabled", False)
        if not isinstance(enabled, bool):
            fail("assay.chip.enabled must be boolean")
        if "panel_tsv" in chip and not isinstance(chip["panel_tsv"], str):
            fail("assay.chip.panel_tsv must be a string path")

    wbc = assay.get("wbc_filter", {})
    if wbc:
        enabled = wbc.get("enabled", False)
        if not isinstance(enabled, bool):
            fail("assay.wbc_filter.enabled must be boolean")
        if "calls_dir" in wbc and not isinstance(wbc["calls_dir"], str):
            fail("assay.wbc_filter.calls_dir must be a string path")
        if "fail_on_support" in wbc and not isinstance(wbc["fail_on_support"], bool):
            fail("assay.wbc_filter.fail_on_support must be boolean")


def validate_clinical_annotations(cfg):
    ann = cfg.get("clinical_annotations", {})
    if not ann:
        return
    if "enabled" in ann and not isinstance(ann["enabled"], bool):
        fail("clinical_annotations.enabled must be boolean")
    if "panel_tsv" in ann and not isinstance(ann["panel_tsv"], str):
        fail("clinical_annotations.panel_tsv must be a string path")


def validate_clinical_output(cfg):
    out = cfg.get("clinical_output", {})
    if not out:
        return
    if "enabled" in out and not isinstance(out["enabled"], bool):
        fail("clinical_output.enabled must be boolean")
    accepted = out.get("accepted_support_gates", [])
    if not isinstance(accepted, list):
        fail("clinical_output.accepted_support_gates must be a list")
    valid = {"PASS", "REVIEW", "FAIL"}
    for gate in accepted:
        if gate not in valid:
            fail(f"clinical_output.accepted_support_gates contains invalid value: {gate}")
    if "include_only_annotated" in out and not isinstance(out["include_only_annotated"], bool):
        fail("clinical_output.include_only_annotated must be boolean")


def validate_pair_repair(cfg):
    rep = cfg.get("pair_repair", {})
    if not rep:
        return
    if "enabled" in rep and not isinstance(rep["enabled"], bool):
        fail("pair_repair.enabled must be boolean")
    if "max_singleton_fraction" in rep:
        frac = rep["max_singleton_fraction"]
        require_positive_number(frac, "pair_repair.max_singleton_fraction", allow_zero=True)
        if frac > 1:
            fail("pair_repair.max_singleton_fraction must be <= 1")


def validate_clinical_gates(cfg):
    gates = cfg.get("clinical_support_gates", {})
    if not gates:
        return
    require_positive_number(gates.get("min_dp", 1), "clinical_support_gates.min_dp")
    require_positive_number(
        gates.get("min_alt_reads", 1), "clinical_support_gates.min_alt_reads"
    )
    min_af = gates.get("min_af", 0.001)
    require_positive_number(min_af, "clinical_support_gates.min_af")
    if min_af >= 1:
        fail("clinical_support_gates.min_af must be < 1")
    low_vaf = gates.get("low_vaf_threshold", 0.01)
    require_positive_number(low_vaf, "clinical_support_gates.low_vaf_threshold")
    if low_vaf >= 1:
        fail("clinical_support_gates.low_vaf_threshold must be < 1")
    if not isinstance(gates.get("require_orthogonal_low_vaf", True), bool):
        fail("clinical_support_gates.require_orthogonal_low_vaf must be boolean")
    if gates.get("chip_flag_action", "review") not in {"review", "fail", "ignore"}:
        fail("clinical_support_gates.chip_flag_action must be review/fail/ignore")


def validate_lod(cfg):
    lod = cfg.get("lod", {})
    if not lod:
        return
    bins = lod.get("bins", [])
    if not isinstance(bins, list) or not bins:
        fail("lod.bins must be a non-empty list")
    for idx, item in enumerate(bins):
        if not isinstance(item, dict):
            fail(f"lod.bins[{idx}] must be a mapping")
        for key in ["name", "min_af", "max_af"]:
            require_key(item, key, f"lod.bins[{idx}]")
        require_positive_number(item["min_af"], f"lod.bins[{idx}].min_af")
        require_positive_number(item["max_af"], f"lod.bins[{idx}].max_af")
        if item["min_af"] >= item["max_af"]:
            fail(f"lod.bins[{idx}] min_af must be < max_af")


def validate_config(path):
    cfg_path = Path(path)
    if not cfg_path.exists():
        fail(f"Config does not exist: {path}")
    with cfg_path.open() as handle:
        cfg = yaml.safe_load(handle)
    if not isinstance(cfg, dict):
        fail(f"Config must parse to a mapping: {path}")

    require_key(cfg, "samples_tsv", str(cfg_path))
    require_key(cfg, "testing", str(cfg_path))
    validate_references(cfg)
    validate_resources(cfg)
    validate_variant_calling(cfg)
    validate_qc_gates(cfg)
    validate_assay(cfg)
    validate_clinical_gates(cfg)
    validate_lod(cfg)
    validate_clinical_annotations(cfg)
    validate_clinical_output(cfg)
    validate_pair_repair(cfg)
    print(f"OK: {path}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python tests/validate_config.py <config.yaml> [<config2.yaml> ...]")
        sys.exit(2)
    for config_path in sys.argv[1:]:
        validate_config(config_path)


if __name__ == "__main__":
    main()
