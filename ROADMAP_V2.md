# ctDNA Pipeline v2 Roadmap

## P0: Clinical-readiness baseline (In Progress)
- Implemented:
  - Calling modes (`tumor_only`, `tumor_normal`, `auto`)
  - ctDNA hard postfilter over `FilterMutectCalls`
  - QC gates and LOD-by-bin report scaffolding
  - Optional WBC/PBMC filtering
  - Optional VarScan cross-check and SnpEff annotation
  - Optional tumor-informed filter (`tumor_informed.*`)
  - Clinical release gate report (`results/reports/clinical_release_gate.tsv`)
- Remaining P0:
  - Add explicit sensitivity/precision regression assertions from workflow outputs
  - Add strict CI smoke run for non-dry-run test fixture

## P1: Molecule-aware error suppression
- Add SSCS consensus path for UMI assays (fgbio/UMI-tools consensus stage).
- Add duplex consensus (DCS) path where assay supports it.
- Gate on consensus molecule metrics (not just read-level DP/AF).
- Add strand/read-position/bias artifacts into release gate.

## P2: Clinical interpretation depth
- Add optional OncoKB/CIViC/AMP tier mapping inputs.
- Add transcript/HGVS normalization and canonical transcript policy.
- Add tumor-informed MRD mode with patient-specific panel generation support.
- Add longitudinal ctDNA trend reporting and reproducibility dashboards.
