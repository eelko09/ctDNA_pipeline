import argparse
import os
import sys
import pandas as pd


def main():
ap = argparse.ArgumentParser()
ap.add_argument("--samples", required=True)
ap.add_argument("--data-dir", required=True)
ap.add_argument("--require-fastq", action="store_true")
args = ap.parse_args()

if not os.path.exists(args.samples):
raise FileNotFoundError(args.samples)

df = pd.read_csv(args.samples, sep=" ").fillna("")

for col in ["sample", "r1", "r2"]:
if col not in df.columns:
raise ValueError(f"samples.tsv missing column: {col}")

if df["sample"].duplicated().any():
dups = df.loc[df["sample"].duplicated(), "sample"].tolist()
raise ValueError(f"Duplicate sample IDs: {dups}")

missing = []
for _, row in df.iterrows():
s = str(row["sample"])
r1 = os.path.join(args.data_dir, str(row["r1"])) if row["r1"] else os.path.join(args.data_dir, f"{s}_R1.fastq.gz")
r2 = os.path.join(args.data_dir, str(row["r2"])) if row["r2"] else os.path.join(args.data_dir, f"{s}_R2.fastq.gz")
if args.require_fastq:
if not os.path.exists(r1):
missing.append(r1)
if not os.path.exists(r2):
missing.append(r2)

if missing:
print("[ERROR] Missing FASTQs:")
for f in missing:
print(" -", f)
sys.exit(2)

print("[OK] Inputs validated")


if __name__ == "__main__":
main()