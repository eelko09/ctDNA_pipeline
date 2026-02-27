# ctDNA Pipeline v2 Roadmap

## Phase 1: Robust Calling Modes (Now Scaffolded)
- Add configurable `variant_calling.mode`:
  - `tumor_only`
  - `tumor_normal`
  - `auto`
- Use `samples.tsv` metadata (`type`, `normal_sample`) to drive called sample selection.
- Enforce validation for `tumor_normal` mode (no missing/unknown matched normals).

## Phase 2: ctDNA-Oriented Postfilters (Now Scaffolded)
- Keep `FilterMutectCalls` as probabilistic filter.
- Add configurable hard postfilter (`min_dp`, `min_alt_reads`, `min_af`, `pass_only`).
- Emit a final VCF target (`*.filtered.final.vcf.gz`) for downstream reporting.

## Phase 3: Objective QC Gates (Now Scaffolded)
- Add configurable QC gates:
  - minimum mapped percentage
  - minimum mean panel coverage
  - maximum duplication fraction
  - maximum contamination
- Emit `results/reports/qc_gates.tsv` with per-sample gate outcomes and overall pass/fail.

## Phase 4: Validation and Regression Testing (Next)
- Add CI dry run:
  - `snakemake -n -s workflow/Snakefile --configfile workflow/config.yaml`
- Add config/schema tests for expected keys and value ranges.
- Add one synthetic low-VAF truth set and assert expected sensitivity bins.

## Phase 5: Assay-Specific Enhancements (Next)
- Add UMI-aware branch (`fgbio`/`umi-tools`) for UMI assays.
- Add optional orthogonal caller cross-check and consensus flagging.
- Add CHIP-focused postfilter/annotation panel for common hematopoietic genes.
