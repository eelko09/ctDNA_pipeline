Place per-sample PBMC call TSV files in this directory as `<sample>.tsv`.

Required columns:
- `CHROM`
- `POS`
- `REF`
- `ALT`

Optional frequency columns used for filtering (`pbmc_blacklist.max_vaf`):
- `AF` or `VAF`
