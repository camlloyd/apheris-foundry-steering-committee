"""
The Day-1 figure, auto-generated from the scorer's CSV. Run it as results
land; the figure must exist before you leave on Day 1 (demos are 09:30).

One figure, two panels, one message:

    left : state recovery rate by arm -- the headline. Negative controls
           (shared pocket, wrong template) are visually flagged so the
           audience sees the steering-vs-perturbation argument at a glance.
           Diamond markers show how often the independent DFG order
           parameter agrees with the RMSD verdict (low agreement = shaky
           state call = say so out loud).
    right: mean ligand centroid distance by arm -- the pose metric, shown
           NEXT TO state recovery, never instead of it. An arm that improves
           pose while leaving the state wrong is a result, not a success
           (KinConfBench, npj Drug Discovery 2026).

Usage:
    python figure.py results.csv --out-prefix figure_state_recovery
    python figure.py results.csv --target-state DFG-out --md summary_table.md
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt

# colorblind-safe: blue levers, orange controls, grey baseline, green combos
COLORS = {"baseline": "#999999", "control": "#FF9400",
          "lever": "#0279EE", "combo": "#75A025", "other_model": "#7E57C2"}


def load_rows(path: str) -> list:
    with open(path) as fh:
        return list(csv.DictReader(fh))


def _arm_kind(arm_id: str) -> str:
    a = arm_id.lower()
    if "baseline" in a or a.startswith("b0"):
        return "baseline"
    if "shared" in a or "wrong" in a:
        return "control"
    if "conflict" in a or "plus" in a or "msa_free" in a and "pocket" in a:
        return "combo"
    if a.startswith("b"):
        return "other_model"
    return "lever"


def _arm_sort_key(arm_id: str):
    m = re.match(r"([A-Z])(\d+)", arm_id)
    return (m.group(1), int(m.group(2))) if m else ("Z", 99)


def aggregate(rows: list) -> dict:
    by_arm = defaultdict(list)
    for r in rows:
        by_arm[r["arm"]].append(r)

    out = {}
    for arm, rs in by_arm.items():
        n = len(rs)
        def _rate(key):
            vals = [r[key] for r in rs if r.get(key) not in ("", None)]
            return sum(float(v) for v in vals) / len(vals) if vals else None
        dists = [float(r["ligand_centroid_dist_A"]) for r in rs
                 if r.get("ligand_centroid_dist_A") not in ("", None)]
        out[arm] = {
            "n": n,
            "state_rate": _rate("target_state_recovered"),
            "dfg_agree": _rate("dfg_agrees_with_target"),
            "decisive": _rate("decisive"),
            "lig_dist": sum(dists) / len(dists) if dists else None,
        }
    return out


def make_figure(results_csv: str, out_prefix: str = "figure_state_recovery",
                target_state: str = "target state") -> str:
    rows = load_rows(results_csv)
    if not rows:
        raise SystemExit(f"no rows in {results_csv}")
    agg = aggregate(rows)
    arms = sorted(agg, key=_arm_sort_key)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2),
                                   gridspec_kw={"width_ratios": [3, 2]})

    x = range(len(arms))
    kinds = [_arm_kind(a) for a in arms]
    colors = [COLORS[k] for k in kinds]
    rates = [agg[a]["state_rate"] or 0 for a in arms]

    bars = ax1.bar(x, rates, color=colors, edgecolor="black", linewidth=0.6)
    for b, k in zip(bars, kinds):
        if k == "control":
            b.set_hatch("//")
    for xi, a in zip(x, arms):
        n_hit = int(round((agg[a]["state_rate"] or 0) * agg[a]["n"]))
        ax1.text(xi, (agg[a]["state_rate"] or 0) + 0.03,
                 f"{n_hit}/{agg[a]['n']}", ha="center", fontsize=8)
        dagree = agg[a]["dfg_agree"]
        if dagree is not None:
            ax1.plot(xi, dagree, marker="D", color="black", markersize=5,
                     linestyle="none")

    ax1.set_xticks(list(x))
    ax1.set_xticklabels([a.replace("_", "\n") for a in arms], fontsize=7)
    ax1.set_ylim(0, 1.12)
    ax1.set_ylabel(f"Fraction landing in {target_state}")
    ax1.set_title("State recovery by arm")

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    handles = [
        Patch(facecolor=COLORS["baseline"], edgecolor="black", label="baseline"),
        Patch(facecolor=COLORS["lever"], edgecolor="black", label="single lever"),
        Patch(facecolor=COLORS["combo"], edgecolor="black", label="combination"),
        Patch(facecolor=COLORS["control"], edgecolor="black", hatch="//",
              label="negative control"),
        Patch(facecolor=COLORS["other_model"], edgecolor="black",
              label="OpenFold3"),
        Line2D([0], [0], marker="D", color="black", linestyle="none",
               markersize=5, label="DFG metric agrees"),
    ]
    ax1.legend(handles=handles, fontsize=7, frameon=False, loc="upper left")

    dists = [agg[a]["lig_dist"] for a in arms]
    ax2.bar(x, [d or 0 for d in dists], color=colors, edgecolor="black",
            linewidth=0.6)
    ax2.set_xticks(list(x))
    ax2.set_xticklabels([a.split("_")[0] for a in arms], fontsize=7)
    ax2.set_ylabel("Mean ligand centroid distance (A)")
    ax2.set_title("Ligand pose accuracy")

    fig.tight_layout()
    png = f"{out_prefix}.png"
    fig.savefig(png, dpi=300)
    fig.savefig(f"{out_prefix}.svg")
    print(f"wrote {png} and {out_prefix}.svg")
    return png


def summary_markdown(results_csv: str, target_state: str = "target state") -> str:
    """The results table as markdown, ready to paste into the slide notes."""
    agg = aggregate(load_rows(results_csv))
    lines = [
        f"| arm | n | state recovery ({target_state}) | DFG agrees | decisive | ligand dist (A) |",
        "|---|---|---|---|---|---|",
    ]
    for arm in sorted(agg, key=_arm_sort_key):
        a = agg[arm]
        def f(v, pct=True):
            if v is None:
                return "-"
            return f"{v:.0%}" if pct else f"{v:.2f}"
        lines.append(
            f"| {arm} | {a['n']} | {f(a['state_rate'])} | "
            f"{f(a['dfg_agree'])} | {f(a['decisive'])} | {f(a['lig_dist'], False)} |"
        )
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results_csv")
    ap.add_argument("--out-prefix", default="figure_state_recovery")
    ap.add_argument("--target-state", default="target state")
    ap.add_argument("--md", default=None, help="also write a markdown table here")
    args = ap.parse_args()

    make_figure(args.results_csv, args.out_prefix, args.target_state)
    if args.md:
        with open(args.md, "w") as fh:
            fh.write(summary_markdown(args.results_csv, args.target_state) + "\n")
        print(f"wrote {args.md}")


if __name__ == "__main__":
    main()
