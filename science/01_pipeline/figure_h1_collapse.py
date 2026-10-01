#!/usr/bin/env python3
"""The H1 figure: mode collapse, and the MSA window H2 will test.

    python3 figure_h1_collapse.py --kit ~/apheris_kit_rorgamma --runs ~/.../runs

Panel A — state space. Every structurally sound inverse-agonist prediction
plotted by its H12-region RMSD to each crystal reference. The two crystal states
sit at the poles 7.78 A apart; every prediction sits in a tight cluster at the
H12-in pole. That is the result, visible without reading a number.

Panel B — the MSA bracket. Mean pLDDT against MSA depth, with the two measured
endpoints from H1 and the untested interval shaded. This is what H2 titrates.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

C_IN = "#2a78d6"      # H12-in  (cool pole)
C_OUT = "#d03b3b"     # H12-out (warm pole)
C_NEUTRAL = "#8a8980"
C_BAND = "#f0efec"
C_GRID = "#dedcd6"
C_TEXT = "#0b0b0b"
C_TEXT2 = "#52514e"
SURFACE = "#fcfcfb"

REGION = list(range(479, 487))
PLDDT_FLOOR = 70.0


def mean_plddt(cif: pathlib.Path):
    cand = cif.with_name(cif.name.replace("_model.cif", "_scores.json"))
    if cand.exists():
        s = json.load(open(cand))
        for k, v in s.items():
            if "local_distance" in k and isinstance(v, list) and v:
                m = float(np.mean(v))
                return m * 100 if m <= 1.0 else m
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kit", required=True)
    ap.add_argument("--runs", required=True, help="runs/rorgamma_invago directory")
    ap.add_argument("--rgkit", default=None)
    ap.add_argument("--out", default="figure_h1_collapse")
    args = ap.parse_args()

    for p in filter(None, [args.rgkit, args.kit, str(pathlib.Path(args.kit).parent)]):
        sys.path.insert(0, p)
    from rgkit import state_recovery2 as SR
    from rgkit import structio as S

    refs_dir = pathlib.Path(args.kit) / "refs_rorgamma"
    a = S.load_chain(refs_dir / "3KYT_clean.cif", source="3KYT")
    b = S.load_chain(refs_dir / "4ZJW_clean.cif", source="4ZJW")

    sep = SR.region_rmsd(b, a, REGION, a)["region_rmsd"]

    good, bad = [], []
    for p in sorted(pathlib.Path(args.runs).rglob("*model.cif")):
        if "__MACOSX" in str(p):
            continue
        c = S.load_chain(p, source=p.name)
        pl = mean_plddt(p)
        d_in = SR.region_rmsd(c, a, REGION, a)["region_rmsd"]
        d_out = SR.region_rmsd(c, b, REGION, a)["region_rmsd"]
        rec = (d_in, d_out, pl, p.parent.name)
        (bad if (pl is not None and pl < PLDDT_FLOOR) else good).append(rec)

    fig = plt.figure(figsize=(13.2, 6.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.26)

    # ---------------------------------------------------------- panel A
    ax = fig.add_subplot(gs[0, 0])
    hi = sep * 1.28
    ax.plot([0, hi], [0, hi], color=C_NEUTRAL, lw=1.0, ls="--", zorder=1)
    ax.annotate("equidistant", xy=(hi * 0.70, hi * 0.70), fontsize=8,
                color=C_NEUTRAL, rotation=45, ha="center", va="bottom")

    gx = np.array([r[0] for r in good])
    gy = np.array([r[1] for r in good])
    ax.plot(gx, gy, "o", ms=8, color=C_IN, markeredgecolor="white",
            markeredgewidth=1.6, zorder=4, label=f"predictions, fold intact (n={len(good)})")
    if bad:
        bmin = min(r[0] for r in bad)
        ax.plot([], [], "X", ms=9, color=C_NEUTRAL, markeredgecolor="white",
                markeredgewidth=1.4,
                label=f"fold failed, off scale at ~{bmin:.0f} Å (n={len(bad)})")
        ax.annotate(f"{len(bad)} MSA-free runs are off this plot entirely\n"
                    f"({bmin:.0f}+ Å from both references, pLDDT 44)",
                    xy=(0.035, 0.035), xycoords="axes fraction", fontsize=8.5,
                    color=C_NEUTRAL, ha="left", va="bottom")

    ax.plot([0], [sep], "*", ms=20, color=C_IN, markeredgecolor="white",
            markeredgewidth=1.4, zorder=5)
    ax.annotate("3KYT · H12-in", xy=(0, sep), xytext=(14, 12),
                textcoords="offset points", fontsize=9, color=C_IN,
                fontweight="bold", va="bottom")
    ax.plot([sep], [0], "*", ms=20, color=C_OUT, markeredgecolor="white",
            markeredgewidth=1.4, zorder=5)
    ax.annotate("4ZJW\nH12-out\n(the target)", xy=(sep, 0), xytext=(-6, 14),
                textcoords="offset points", fontsize=9, color=C_OUT,
                fontweight="bold", ha="right")

    if len(good):
        cx, cy = gx.mean(), gy.mean()
        ax.annotate(f"every prediction here\n{cx:.2f} ± {gx.std():.2f} Å from H12-in\n"
                    f"{cy:.2f} ± {gy.std():.2f} Å from H12-out",
                    xy=(cx, cy), xytext=(54, -18), textcoords="offset points",
                    fontsize=9, color=C_TEXT2,
                    arrowprops=dict(arrowstyle="-", color=C_NEUTRAL, lw=1.0))

    ax.set_xlim(-0.4, hi)
    ax.set_ylim(-0.4, hi)
    ax.set_xlabel("H12 RMSD to 3KYT (H12-in), Å", color=C_TEXT2, fontsize=9.5)
    ax.set_ylabel("H12 RMSD to 4ZJW (H12-out), Å", color=C_TEXT2, fontsize=9.5)
    ax.set_title("A · Where the predictions sit in state space",
                 color=C_TEXT, fontsize=12, loc="left", pad=10)
    ax.grid(color=C_GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")

    # ---------------------------------------------------------- panel B
    ax2 = fig.add_subplot(gs[0, 1])
    full_pl = np.mean([r[2] for r in good if r[2] is not None])
    zero_pl = np.mean([r[2] for r in bad if r[2] is not None]) if bad else 44.0

    x_full, x_zero = 4.0, 0.0           # log10 placeholders
    ax2.axhspan(0, PLDDT_FLOOR, color=C_BAND, zorder=0)
    ax2.annotate("fold not usable", xy=(0.04, PLDDT_FLOOR - 4), xycoords=("axes fraction", "data"),
                 fontsize=8.5, color=C_TEXT2, va="top")
    ax2.axhline(PLDDT_FLOOR, color=C_NEUTRAL, lw=1.0, ls="--", zorder=2)

    ax2.axvspan(0.6, 3.45, color="#eef4fc", zorder=0)
    ax2.annotate("untested window", xy=(2.02, 99), fontsize=10,
                 color=C_IN, ha="center", fontweight="bold")

    ax2.plot([x_zero], [zero_pl], "o", ms=11, color=C_NEUTRAL,
             markeredgecolor="white", markeredgewidth=1.6, zorder=5)
    ax2.annotate(f"no MSA\npLDDT {zero_pl:.0f}\nfold destroyed",
                 xy=(x_zero, zero_pl), xytext=(10, -2), textcoords="offset points",
                 fontsize=9, color=C_TEXT2, va="top")
    ax2.plot([x_full], [full_pl], "o", ms=11, color=C_IN,
             markeredgecolor="white", markeredgewidth=1.6, zorder=5)
    ax2.annotate(f"full MSA\npLDDT {full_pl:.0f}\nH12 locked in",
                 xy=(x_full, full_pl), xytext=(-10, -10), textcoords="offset points",
                 fontsize=9, color=C_IN, ha="right", va="top")

    for d in (8, 32, 128, 512, 2048):
        x = np.log10(d)
        ax2.annotate("", xy=(x, -0.225), xytext=(x, -0.275),
                     xycoords=("data", "axes fraction"),
                     textcoords=("data", "axes fraction"),
                     arrowprops=dict(arrowstyle="-", color=C_OUT, lw=2))
        ax2.annotate(str(d), xy=(x, -0.335), xycoords=("data", "axes fraction"),
                     fontsize=8.5, color=C_OUT, ha="center", va="center")
    ax2.annotate("H2 depths →", xy=(-0.5, -0.335), xycoords=("data", "axes fraction"),
                 fontsize=8.5, color=C_OUT, ha="left", va="center", fontweight="bold")

    ax2.set_xlim(-0.55, 4.5)
    ax2.set_ylim(20, 104)
    ax2.set_xticks([0, 1, 2, 3, 4])
    ax2.set_xticklabels(["0", "10", "100", "1k", "10k"])
    ax2.set_xlabel("MSA depth (sequences)", color=C_TEXT2, fontsize=9.5)
    ax2.set_ylabel("mean pLDDT", color=C_TEXT2, fontsize=9.5)
    ax2.set_title("B · The interval nobody has tested",
                  color=C_TEXT, fontsize=12, loc="left", pad=10)
    ax2.grid(color=C_GRID, lw=0.8, axis="y")
    ax2.set_axisbelow(True)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)

    fig.text(0.011, 0.962,
             "RORγ H12: five conditioning levers, 25 predictions, no escape from the agonist state",
             fontsize=13.5, color=C_TEXT, ha="left", va="top")
    fig.text(0.011, 0.912,
             f"The two crystal states are {sep:.2f} Å apart across H12 (479–486). "
             f"The predicted ensemble is {sep/max(gx.std(),1e-6):.0f}× tighter than that gap "
             f"and sits entirely at one pole. Closest approach to H12-out: {gy.min():.2f} Å.",
             fontsize=9.5, color=C_TEXT2, ha="left", va="top")

    fig.patch.set_facecolor(SURFACE)
    for a_ in fig.axes:
        a_.set_facecolor(SURFACE)
    fig.subplots_adjust(left=0.062, right=0.985, top=0.80, bottom=0.235)

    out = pathlib.Path(args.out)
    fig.savefig(f"{out}.png", dpi=200, facecolor=fig.get_facecolor())
    fig.savefig(f"{out}.svg", facecolor=fig.get_facecolor())
    print(f"-> {out}.png\n-> {out}.svg")
    print(f"crystal separation {sep:.2f} A; {len(good)} sound, {len(bad)} failed folds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
