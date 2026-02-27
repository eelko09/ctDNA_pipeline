#!/usr/bin/env python
import argparse
import pandas as pd
import os
from jinja2 import Environment, FileSystemLoader

def read_variant_tables(files):
    dfs = []
    for f in files:
        if os.path.exists(f) and os.path.getsize(f) > 0:
            df = pd.read_csv(f, sep="\t")
            df['sample'] = os.path.basename(f).split(".")[0]
            dfs.append(df)
    if dfs:
        return pd.concat(dfs, ignore_index=True)
    else:
        return pd.DataFrame()

def read_qc_tables(flagstats, samtools_stats, dup_metrics, coverage, contamination):
    records = []
    for f_flag, f_stats, f_dup, f_cov, f_contam in zip(flagstats, samtools_stats, dup_metrics, coverage, contamination):
        sample = os.path.basename(f_flag).split(".")[0]
        record = {'sample': sample}

        # Flagstat
        #if os.path.exists(f_flag) and os.path.getsize(f_flag) > 0:
        #    with open(f_flag) as fh:
        #        record['flagstat'] = fh.read().strip()

        # Samtools stats
        #if os.path.exists(f_stats) and os.path.getsize(f_stats) > 0:
        #    with open(f_stats) as fh:
        #        record['samtools_stats'] = fh.read().strip()

        # Duplication metrics
        #if os.path.exists(f_dup) and os.path.getsize(f_dup) > 0:
        #    with open(f_dup) as fh:
        #        record['dup_metrics'] = fh.read().strip()

        # Coverage
        if os.path.exists(f_cov) and os.path.getsize(f_cov) > 0:
            cov_df = pd.read_csv(f_cov, sep="\t")
            record['mean_coverage'] = cov_df['mean'][0] if 'mean' in cov_df.columns else None

        # Contamination
        if os.path.exists(f_contam) and os.path.getsize(f_contam) > 0:
            cont_df = pd.read_csv(f_contam, sep="\t")
            record['contamination'] = cont_df.iloc[0,1] if not cont_df.empty else None

        records.append(record)
    return pd.DataFrame(records)

def generate_html_report(qc_df, variants_df, out_html):
    script_dir = os.path.dirname(__file__)
    template_name = "ctdna_report.html.j2"
    template_dirs = [
        os.path.join(script_dir, "templates"),
        script_dir,
    ]
    env = Environment(loader=FileSystemLoader(searchpath=template_dirs))
    template = env.get_template(template_name)

    html_content = template.render(qc=qc_df.to_dict(orient="records"),
                                   variants=variants_df.to_dict(orient="records"))
    with open(out_html, "w") as fh:
        fh.write(html_content)

def main():
    parser = argparse.ArgumentParser(description="Generate ctDNA summary report")
    parser.add_argument("--variants-matrix", required=True, nargs="+", help="Variant tables per sample")
    parser.add_argument("--qc", required=True, nargs="+", help="Flagstat files per sample")
    parser.add_argument("--samtools-stats", required=True, nargs="+", help="Samtools stats per sample")
    parser.add_argument("--dup-metrics", required=True, nargs="+", help="Duplication metrics per sample")
    parser.add_argument("--coverage", required=True, nargs="+", help="Mosdepth coverage summaries per sample")
    parser.add_argument("--contamination", required=True, nargs="+", help="Mutect2 contamination tables per sample")
    parser.add_argument("--out", required=True, help="Output HTML report path")
    parser.add_argument("--qc-out", required=False, help="Optional QC summary TSV output")
    parser.add_argument("--variants-out", required=False, help="Optional variant summary TSV output")
    parser.add_argument("--samples-tsv", required=False, help="Optional samples TSV")
    parser.add_argument("--results-dir", required=False, help="Optional results directory")
    args = parser.parse_args()

    # Read variants
    variants_df = read_variant_tables(args.variants_matrix)

    # Read QC
    qc_df = read_qc_tables(args.qc, args.samtools_stats, args.dup_metrics, args.coverage, args.contamination)

    # Optional TSV outputs
    if args.qc_out:
        qc_df.to_csv(args.qc_out, sep="\t", index=False)
    if args.variants_out:
        variants_df.to_csv(args.variants_out, sep="\t", index=False)

    # Generate HTML
    generate_html_report(qc_df, variants_df, args.out)

if __name__ == "__main__":
    main()
