#!/usr/bin/env python3
import argparse
import csv
import sys


def parse_bin_spec(spec):
    # Example: 0.001:0.005=0.60
    try:
        span, threshold = spec.split("=")
        lower, upper = span.split(":")
        return (float(lower), float(upper), float(threshold))
    except Exception as exc:
        raise ValueError(f"Invalid bin spec '{spec}': {exc}") from exc


def load_truth(path):
    truth = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            truth[row["variant_id"]] = float(row["truth_af"])
    return truth


def load_calls(path):
    calls = set()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            calls.add(row["variant_id"])
    return calls


def main():
    parser = argparse.ArgumentParser(description="Assert sensitivity across low-VAF bins.")
    parser.add_argument("--truth", required=True, help="TSV with columns: variant_id, truth_af")
    parser.add_argument("--calls", required=True, help="TSV with column: variant_id")
    parser.add_argument(
        "--bin",
        action="append",
        required=True,
        help="Bin spec: lower:upper=min_recall, e.g. 0.001:0.005=0.60",
    )
    args = parser.parse_args()

    bins = [parse_bin_spec(item) for item in args.bin]
    truth = load_truth(args.truth)
    calls = load_calls(args.calls)

    failed = False
    for lower, upper, min_recall in bins:
        in_bin = [vid for vid, af in truth.items() if lower <= af < upper]
        if not in_bin:
            print(f"WARN: no truth variants in bin [{lower}, {upper})")
            continue
        detected = sum(1 for vid in in_bin if vid in calls)
        recall = detected / len(in_bin)
        print(
            f"bin [{lower:.4f}, {upper:.4f}) "
            f"detected={detected}/{len(in_bin)} recall={recall:.3f} threshold={min_recall:.3f}"
        )
        if recall < min_recall:
            failed = True
            print(
                f"FAIL: sensitivity below threshold in bin [{lower}, {upper}): "
                f"{recall:.3f} < {min_recall:.3f}"
            )

    if failed:
        sys.exit(1)
    print("PASS: all sensitivity bins meet thresholds")


if __name__ == "__main__":
    main()
