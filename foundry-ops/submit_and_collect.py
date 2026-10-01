#!/usr/bin/env python3
"""Collect finished Nextflow predict-workflow outputs into preds/<arm>/seed_<seed>/
and append rows to a job tracker CSV. Idempotent - safe to re-run as more
group output dirs finish throughout the day.

Usage:
  python3 submit_and_collect.py --arm-manifest arm_manifest.csv \
      --job-tracker job_tracker.csv --preds-dir preds GROUP_OUTPUT_DIR [GROUP_OUTPUT_DIR ...]
"""
import argparse
import csv
import datetime
import json
import os
import shutil
import sys

TRACKER_FIELDS = [
    "arm_id", "query_id", "seed", "sample", "status",
    "cif_path", "scores_path", "timestamp", "group_dir",
]


def load_arm_manifest(path):
    mapping = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            mapping[row["query_id"]] = row["arm_id"]
    return mapping


def load_existing_keys(tracker_path):
    keys = set()
    if os.path.exists(tracker_path):
        with open(tracker_path, newline="") as f:
            for row in csv.DictReader(f):
                keys.add((row["arm_id"], row["query_id"], row["seed"], row["sample"]))
    return keys


def collect_group(group_dir, arm_manifest, preds_dir, existing_keys, tracker_rows):
    result_path = os.path.join(group_dir, "result.json")
    if not os.path.exists(result_path):
        print(f"skip {group_dir}: no result.json (not finished yet)", file=sys.stderr)
        return
    with open(result_path) as f:
        result = json.load(f)

    for pred in result["result"]["predictions"]:
        query_id = pred["query_id"]
        seed = str(pred["seed"])
        sample = str(pred["sample"])
        arm_id = arm_manifest.get(query_id)
        if arm_id is None:
            print(f"warn: query_id {query_id!r} not in arm manifest, skipping", file=sys.stderr)
            continue

        key = (arm_id, query_id, seed, sample)

        dest_dir = os.path.join(preds_dir, arm_id, f"seed_{seed}")
        os.makedirs(dest_dir, exist_ok=True)
        src_cif = os.path.join(group_dir, pred["path"])
        src_scores = os.path.join(group_dir, pred["scores_path"])
        dest_cif = os.path.join(dest_dir, os.path.basename(pred["path"]))
        dest_scores = os.path.join(dest_dir, os.path.basename(pred["scores_path"]))

        if not os.path.exists(dest_cif):
            shutil.copy2(src_cif, dest_cif)
        if not os.path.exists(dest_scores):
            shutil.copy2(src_scores, dest_scores)

        if key in existing_keys:
            continue

        tracker_rows.append({
            "arm_id": arm_id,
            "query_id": query_id,
            "seed": seed,
            "sample": sample,
            "status": "collected",
            "cif_path": dest_cif,
            "scores_path": dest_scores,
            "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "group_dir": group_dir,
        })
        existing_keys.add(key)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm-manifest", required=True)
    ap.add_argument("--job-tracker", required=True)
    ap.add_argument("--preds-dir", required=True)
    ap.add_argument("group_dirs", nargs="+")
    args = ap.parse_args()

    arm_manifest = load_arm_manifest(args.arm_manifest)
    existing_keys = load_existing_keys(args.job_tracker)
    tracker_rows = []

    for group_dir in args.group_dirs:
        collect_group(group_dir, arm_manifest, args.preds_dir, existing_keys, tracker_rows)

    if tracker_rows:
        write_header = not os.path.exists(args.job_tracker)
        with open(args.job_tracker, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=TRACKER_FIELDS)
            if write_header:
                writer.writeheader()
            writer.writerows(tracker_rows)

    print(f"collected {len(tracker_rows)} new prediction(s) into {args.preds_dir}")


if __name__ == "__main__":
    main()
