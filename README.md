# ctDNA / cfDNA Snakemake Pipeline

Production-oriented Snakemake workflow for ctDNA/cfDNA short-read analysis:
FastQC -> fastp -> BWA -> duplicate marking -> BQSR -> Mutect2 -> filtering -> report generation.

## Repository layout

```
.
├── Snakefile                  # Entrypoint (includes workflow/Snakefile)
├── workflow/
│   ├── Snakefile              # Main workflow
│   ├── config.yaml            # Pipeline configuration
│   ├── samples.tsv            # Sample sheet
│   └── rules/
├── envs/                      # Conda env specs per tool
├── scripts/                   # Reporting and helper scripts
├── profiles/                  # laptop/slurm Snakemake profiles
├── tests/                     # small test fixtures
└── ROADMAP_V2.md
```

## Sample sheet format

Required columns:
- `sample`
- `R1_fastq`
- `R2_fastq`

Optional columns:
- `type` (`tumor` or `normal`)
- `normal_sample` (matched normal sample ID for tumor rows)

Example:

```tsv
sample	type	R1_fastq	R2_fastq	normal_sample
TUMOR001	tumor	TUMOR001_R1.fastq.gz	TUMOR001_R2.fastq.gz	NORMAL001
NORMAL001	normal	NORMAL001_R1.fastq.gz	NORMAL001_R2.fastq.gz
```

## Quick start

1) Create a Snakemake environment:

```bash
mamba create -n snakemake -c conda-forge -c bioconda snakemake
mamba activate snakemake
```

2) Edit:
- `workflow/config.yaml`
- `workflow/samples.tsv`

3) Run lint:

```bash
snakemake --lint -s workflow/Snakefile --configfile workflow/config.yaml
```

4) Run pipeline:

```bash
snakemake -s workflow/Snakefile --configfile workflow/config.yaml --profile profiles/laptop --use-conda --rerun-incomplete
```

## Calling mode (v2 scaffold)

Configured in `workflow/config.yaml`:
- `variant_calling.mode`: `tumor_only`, `tumor_normal`, or `auto`
- `variant_calling.mutect2.min_allele_fraction`
- `variant_calling.postfilter.*` thresholds
- `qc_gates.*` thresholds for `results/reports/qc_gates.tsv`

## Assay-specific enhancements (Phase 5 scaffold)

- UMI branch:
  - enable `assay.umi.enabled: true`
  - configure `assay.umi.bc_pattern` for your assay
  - pipeline will run `umi_tools extract` before `fastp`
- Orthogonal cross-check:
  - enable `assay.orthogonal.enabled: true`
  - place per-sample TSV files in `assay.orthogonal.calls_dir` named `<sample>.tsv`
  - required columns: `CHROM`, `POS`, `REF`, `ALT`
- CHIP panel flagging:
  - enable `assay.chip.enabled: true`
  - provide `assay.chip.panel_tsv` with columns: `CHROM`, `POS`, `REF`, `ALT`, `GENE`
- Output:
  - `results/variants/{sample}.variants.flagged.tsv`
  - includes `consensus_flag`, `orthogonal_support`, `chip_flag`, `chip_gene`

## Clinical support and audit outputs

- Variant-level support gate fields are added in flagged TSVs:
  - `support_gate` (`PASS` / `REVIEW` / `FAIL`)
  - `support_reasons`
- LOD/callable summary:
  - `results/reports/lod_by_bin.tsv`
- Run audit manifest:
  - `results/reports/run_manifest.json` (includes config hash, sample lists, git SHA if available)

## Optional matched-WBC and clinical outputs

- Matched-WBC filtering:
  - set `assay.wbc_filter.enabled: true`
  - provide `assay.wbc_filter.calls_dir/<normal_sample>.tsv` with `CHROM,POS,REF,ALT`
  - variants seen in matched WBC can be hard-failed (`fail_on_support: true`)
- Clinical annotations:
  - set `clinical_annotations.enabled: true`
  - provide `clinical_annotations.panel_tsv` with:
    - `CHROM, POS, REF, ALT, GENE, CLINICAL_TIER, ACTIONABILITY`
- Clinical output gating:
  - configure `clinical_output.accepted_support_gates` (e.g. `["PASS","REVIEW"]`)
  - optional `clinical_output.include_only_annotated: true`
  - output table: `results/variants/{sample}.clinical.tsv`

## Optional tumor-informed mode (P0 scaffold)

- Enable `tumor_informed.enabled: true`
- Provide `tumor_informed.known_variants_dir/<sample>.tsv` with:
  - `CHROM`, `POS`, `REF`, `ALT`
- Controls:
  - `tumor_informed.require_known` (keep only known tissue variants)
  - `tumor_informed.fail_on_missing_known`
- Output:
  - `results/variants/{sample}.clinical.tumor_informed.tsv`

## Clinical release gating (P0 scaffold)

- Configure `clinical_release.*`:
  - `require_qc_pass`
  - `require_manifest_git_sha`
  - `require_variants`
- Output:
  - `results/reports/clinical_release_gate.tsv`

## Validation and regression checks

- Config/schema validation:
  - `python tests/validate_config.py workflow/config.yaml tests/config.test.yaml tests/config.ci.yaml`
- Synthetic low-VAF sensitivity bins:
  - `python tests/assert_sensitivity_bins.py --truth tests/truthset/synthetic_truth.tsv --calls tests/truthset/synthetic_calls.tsv --bin 0.001:0.005=0.60 --bin 0.005:0.010=0.80 --bin 0.010:0.050=0.95`
- CI dry-run (mock references):
  - `snakemake -n -s workflow/Snakefile --configfile tests/config.ci.yaml --cores 1`

## Notes for publication

- Keep large data and runtime outputs out of git (see `.gitignore`).
- Do not publish patient/sample-derived FASTQ/BAM/VCF files.
- Ensure reference files are distributed via documented download/build steps, not committed binaries.

## License

MIT. See `LICENSE`.
