#!/usr/bin/env python3
import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone


def git_sha():
    try:
        return (
            subprocess.check_output(
                ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL, text=True
            )
            .strip()
        )
    except Exception:
        return "unknown"


def main():
    parser = argparse.ArgumentParser(description="Emit run manifest JSON.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--called-samples", required=True)
    parser.add_argument("--all-samples", required=True)
    parser.add_argument("--config-json", required=True)
    args = parser.parse_args()

    config_obj = json.loads(args.config_json)
    config_canonical = json.dumps(config_obj, sort_keys=True, separators=(",", ":"))
    config_sha256 = hashlib.sha256(config_canonical.encode("utf-8")).hexdigest()

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_sha": git_sha(),
        "called_samples": [s for s in args.called_samples.split(",") if s],
        "all_samples": [s for s in args.all_samples.split(",") if s],
        "config_sha256": config_sha256,
        "pipeline": {
            "name": "ctDNA_pipeline",
            "supports": {
                "calling_mode": True,
                "hard_postfilter": True,
                "qc_gates": True,
                "lod_by_bin": True,
                "variant_support_gates": True,
            },
        },
    }

    with open(args.output, "w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
