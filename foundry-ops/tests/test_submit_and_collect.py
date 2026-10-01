import csv
import json
import os

import pytest

import submit_and_collect as sac


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def make_group(tmp_path, name, predictions, extra_files=True):
    """Create a fake Nextflow --output dir with a result.json and the
    .cif/.scores.json files it references."""
    group_dir = tmp_path / name
    group_dir.mkdir()
    result = {"result": {"predictions": predictions}}
    (group_dir / "result.json").write_text(json.dumps(result))

    if extra_files:
        for pred in predictions:
            cif_path = group_dir / pred["path"]
            cif_path.parent.mkdir(parents=True, exist_ok=True)
            cif_path.write_text(f"CIF DATA for {pred['query_id']}")
            scores_path = group_dir / pred["scores_path"]
            scores_path.parent.mkdir(parents=True, exist_ok=True)
            scores_path.write_text(json.dumps({"plddt": 90.0}))

    return group_dir


PRED_A0 = {
    "path": "a0_baseline_agonist/seed_42/a0_baseline_agonist_seed_42_sample_1_model.cif",
    "query_id": "a0_baseline_agonist",
    "query_index": 0,
    "sample": 1,
    "scores_path": "a0_baseline_agonist/seed_42/a0_baseline_agonist_seed_42_sample_1_scores.json",
    "seed": 42,
}


# ---------------------------------------------------------------------------
# load_arm_manifest / load_existing_keys
# ---------------------------------------------------------------------------

def test_load_arm_manifest(tmp_path):
    manifest_path = tmp_path / "arm_manifest.csv"
    write_csv(
        manifest_path,
        [{"query_id": "a0_baseline_agonist", "arm_id": "A0", "seed_count": "2",
          "conditioning": "none", "ligand": "25-HC", "group_dir": "group1"}],
        ["query_id", "arm_id", "seed_count", "conditioning", "ligand", "group_dir"],
    )
    mapping = sac.load_arm_manifest(manifest_path)
    assert mapping == {"a0_baseline_agonist": "A0"}


def test_load_existing_keys_missing_file_returns_empty_set(tmp_path):
    assert sac.load_existing_keys(tmp_path / "does_not_exist.csv") == set()


def test_load_existing_keys_parses_rows(tmp_path):
    tracker_path = tmp_path / "job_tracker.csv"
    write_csv(
        tracker_path,
        [{"arm_id": "A0", "query_id": "a0_baseline_agonist", "seed": "42", "sample": "1",
          "status": "collected", "cif_path": "x.cif", "scores_path": "x.json",
          "timestamp": "t", "group_dir": "g"}],
        sac.TRACKER_FIELDS,
    )
    keys = sac.load_existing_keys(tracker_path)
    assert keys == {("A0", "a0_baseline_agonist", "42", "1")}


# ---------------------------------------------------------------------------
# collect_group
# ---------------------------------------------------------------------------

def test_collect_group_skips_unfinished_group(tmp_path, capsys):
    group_dir = tmp_path / "not_done_yet"
    group_dir.mkdir()

    tracker_rows = []
    sac.collect_group(str(group_dir), {}, str(tmp_path / "preds"), set(), tracker_rows)

    assert tracker_rows == []
    assert "no result.json" in capsys.readouterr().err


def test_collect_group_copies_files_and_appends_tracker_row(tmp_path):
    group_dir = make_group(tmp_path, "group1", [PRED_A0])
    preds_dir = tmp_path / "preds"
    arm_manifest = {"a0_baseline_agonist": "A0"}
    tracker_rows = []
    existing_keys = set()

    sac.collect_group(str(group_dir), arm_manifest, str(preds_dir), existing_keys, tracker_rows)

    dest_cif = preds_dir / "A0" / "seed_42" / "a0_baseline_agonist_seed_42_sample_1_model.cif"
    dest_scores = preds_dir / "A0" / "seed_42" / "a0_baseline_agonist_seed_42_sample_1_scores.json"
    assert dest_cif.exists()
    assert dest_scores.exists()

    assert len(tracker_rows) == 1
    row = tracker_rows[0]
    assert row["arm_id"] == "A0"
    assert row["query_id"] == "a0_baseline_agonist"
    assert row["seed"] == "42"
    assert row["sample"] == "1"
    assert row["status"] == "collected"
    assert ("A0", "a0_baseline_agonist", "42", "1") in existing_keys


def test_collect_group_warns_and_skips_unknown_query_id(tmp_path, capsys):
    group_dir = make_group(tmp_path, "group1", [PRED_A0])
    tracker_rows = []

    sac.collect_group(str(group_dir), {}, str(tmp_path / "preds"), set(), tracker_rows)

    assert tracker_rows == []
    assert "not in arm manifest" in capsys.readouterr().err


def test_collect_group_is_idempotent_on_rerun(tmp_path):
    group_dir = make_group(tmp_path, "group1", [PRED_A0])
    preds_dir = tmp_path / "preds"
    arm_manifest = {"a0_baseline_agonist": "A0"}

    tracker_rows = []
    existing_keys = set()
    sac.collect_group(str(group_dir), arm_manifest, str(preds_dir), existing_keys, tracker_rows)
    assert len(tracker_rows) == 1

    # Simulate a second invocation with a fresh tracker_rows list but the
    # same existing_keys (as main() would reconstruct from job_tracker.csv).
    tracker_rows_second = []
    sac.collect_group(str(group_dir), arm_manifest, str(preds_dir), existing_keys, tracker_rows_second)
    assert tracker_rows_second == []


def test_collect_group_does_not_overwrite_existing_dest_files(tmp_path):
    group_dir = make_group(tmp_path, "group1", [PRED_A0])
    preds_dir = tmp_path / "preds"
    dest_cif = preds_dir / "A0" / "seed_42" / "a0_baseline_agonist_seed_42_sample_1_model.cif"
    dest_cif.parent.mkdir(parents=True)
    dest_cif.write_text("already here, do not clobber")

    sac.collect_group(str(group_dir), {"a0_baseline_agonist": "A0"}, str(preds_dir), set(), [])

    assert dest_cif.read_text() == "already here, do not clobber"


# ---------------------------------------------------------------------------
# main() - integration tests
# ---------------------------------------------------------------------------

def run_main(monkeypatch, args):
    monkeypatch.setattr("sys.argv", ["submit_and_collect.py"] + args)
    sac.main()


def test_main_end_to_end_writes_tracker_with_header(tmp_path, monkeypatch, capsys):
    group_dir = make_group(tmp_path, "group1", [PRED_A0])
    manifest_path = tmp_path / "arm_manifest.csv"
    write_csv(
        manifest_path,
        [{"query_id": "a0_baseline_agonist", "arm_id": "A0", "seed_count": "2",
          "conditioning": "none", "ligand": "25-HC", "group_dir": "group1"}],
        ["query_id", "arm_id", "seed_count", "conditioning", "ligand", "group_dir"],
    )
    tracker_path = tmp_path / "job_tracker.csv"
    preds_dir = tmp_path / "preds"

    run_main(monkeypatch, [
        "--arm-manifest", str(manifest_path),
        "--job-tracker", str(tracker_path),
        "--preds-dir", str(preds_dir),
        str(group_dir),
    ])

    assert "collected 1 new prediction(s)" in capsys.readouterr().out
    rows = list(csv.DictReader(open(tracker_path)))
    assert len(rows) == 1
    assert rows[0]["arm_id"] == "A0"


def test_main_rerun_collects_zero_new_predictions(tmp_path, monkeypatch, capsys):
    group_dir = make_group(tmp_path, "group1", [PRED_A0])
    manifest_path = tmp_path / "arm_manifest.csv"
    write_csv(
        manifest_path,
        [{"query_id": "a0_baseline_agonist", "arm_id": "A0", "seed_count": "2",
          "conditioning": "none", "ligand": "25-HC", "group_dir": "group1"}],
        ["query_id", "arm_id", "seed_count", "conditioning", "ligand", "group_dir"],
    )
    tracker_path = tmp_path / "job_tracker.csv"
    preds_dir = tmp_path / "preds"
    args = [
        "--arm-manifest", str(manifest_path),
        "--job-tracker", str(tracker_path),
        "--preds-dir", str(preds_dir),
        str(group_dir),
    ]

    run_main(monkeypatch, args)
    capsys.readouterr()  # discard first run's output
    run_main(monkeypatch, args)

    assert "collected 0 new prediction(s)" in capsys.readouterr().out
    rows = list(csv.DictReader(open(tracker_path)))
    assert len(rows) == 1  # still just the one row, not duplicated


def test_main_multiple_group_dirs_in_one_invocation(tmp_path, monkeypatch):
    pred_b = dict(PRED_A0)
    pred_b["path"] = "a1_active_steer/seed_42/a1_active_steer_seed_42_sample_1_model.cif"
    pred_b["scores_path"] = "a1_active_steer/seed_42/a1_active_steer_seed_42_sample_1_scores.json"
    pred_b["query_id"] = "a1_active_steer"

    group1 = make_group(tmp_path, "group1", [PRED_A0])
    group2 = make_group(tmp_path, "group2", [pred_b])

    manifest_path = tmp_path / "arm_manifest.csv"
    write_csv(
        manifest_path,
        [
            {"query_id": "a0_baseline_agonist", "arm_id": "A0", "seed_count": "2",
             "conditioning": "none", "ligand": "25-HC", "group_dir": "group1"},
            {"query_id": "a1_active_steer", "arm_id": "A1", "seed_count": "2",
             "conditioning": "active_pocket", "ligand": "25-HC", "group_dir": "group2"},
        ],
        ["query_id", "arm_id", "seed_count", "conditioning", "ligand", "group_dir"],
    )
    tracker_path = tmp_path / "job_tracker.csv"
    preds_dir = tmp_path / "preds"

    run_main(monkeypatch, [
        "--arm-manifest", str(manifest_path),
        "--job-tracker", str(tracker_path),
        "--preds-dir", str(preds_dir),
        str(group1), str(group2),
    ])

    rows = list(csv.DictReader(open(tracker_path)))
    assert {r["arm_id"] for r in rows} == {"A0", "A1"}
