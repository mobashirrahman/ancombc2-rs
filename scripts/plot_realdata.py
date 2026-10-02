#!/usr/bin/env python3
"""Layer 4 -- the ``y = x`` scatters the definition of done asks for.

PLAN.md §17 requires that real-dataset agreement be inspected "visually on
``y = x``". ``scripts/compare_realdata.py`` computes the numbers; this draws them.

    python3 scripts/plot_realdata.py --data validation/realdata/results \\
        --out validation/realdata/results/scatter

One SVG per dataset, four panels -- ``beta``, ``se``, ``p``, ``q`` -- each with the
identity line drawn through it, because a correlation of 1.0 and a scatter that is
obviously not the identity are different claims and only the second one shows which
it is. A fifth panel plots the significant-call indicator.

SVG rather than PNG, and matplotlib rather than a hand-rolled renderer: SVG is text,
so the figures are diffable in git and a reader can see what was committed without a
renderer, and matplotlib handles the axis scaling that would otherwise be the fiddly
part.

Both arms are joined on the taxon name, matching ``compare_realdata.py``: the two
arms emit their own taxon order. Points are drawn semi-transparent so the density is
visible, and points where the two arms disagree on significance are drawn last and
in a different colour, so a scatter that looks perfect but disagrees on 3 calls is
not silently reassuring.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
# Keep text as text in the SVG. The default (`svg.fonttype: path`) converts every
# label to a glyph outline, which is why the first version of these figures was
# 127 kB of unsearchable path data and could not be reviewed in a diff. With
# `none`, the axis labels and the `max|d|` annotations are readable in the file
# itself, which is the entire reason for choosing SVG over PNG.
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt  # noqa: E402

PANELS = ("beta", "se", "p", "q")
AXIS_LABEL = {"beta": "beta", "se": "SE", "p": "p-value", "q": "q-value"}


def load_arm(path: Path) -> dict[str, dict]:
    """One record per dataset from an arm's jsonl.

    Each quantity is a *list parallel to* ``taxa``, already reduced to the one
    coefficient this comparison uses -- ``compare_realdata.py`` writes it that way,
    so the taxon name is what pairs the two arms and the list index alone is not
    safe to zip: the arms emit their own order and it is not the same one.
    """
    out: dict[str, dict] = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            taxa = rec["taxa"]
            entry = {"coefficient": rec["coefficient"], "vals": {}}
            for field in PANELS:
                vals = rec.get(field) or []
                entry["vals"][field] = dict(zip(taxa, vals))
            entry["sig"] = dict(zip(taxa, rec.get("diff_abn") or []))
            out[rec["dataset"]] = entry
    return out


def pairs(rust: dict, r: dict, field: str) -> tuple[list, list, list]:
    xs, ys, disagree = [], [], []
    ru = rust["vals"].get(field, {})
    rv = r["vals"].get(field, {})
    for taxon in sorted(set(ru) & set(rv)):
        a, b = ru[taxon], rv[taxon]
        if a is None or b is None or not _finite(a) or not _finite(b):
            continue
        sa = bool(rust["sig"].get(taxon))
        sb = bool(r["sig"].get(taxon))
        xs.append(a)
        ys.append(b)
        disagree.append(sa != sb)
    return xs, ys, disagree


def _finite(v) -> bool:
    try:
        return not math.isnan(float(v)) and not math.isinf(float(v))
    except (TypeError, ValueError):
        return False


def draw_axis(ax, xs, ys, disagree, label: str, title: str) -> None:
    finite = [v for v in xs + ys if _finite(v)]
    if not finite:
        ax.set_title(f"{title}: no finite pairs")
        ax.set_xlabel("R")
        ax.set_ylabel("Rust")
        return
    lo, hi = min(finite), max(finite)
    # `q` and `p` live on [0, 1]; a linear identity line is what the plan asks for
    # and is kept rather than switching axes, because switching axes is exactly the
    # thing that can make a disagreement look like agreement.
    ax.plot([lo, hi], [lo, hi], color="0.55", lw=1.0, ls="--", zorder=1)
    agree_x = [x for x, d in zip(xs, disagree) if not d]
    agree_y = [y for y, d in zip(ys, disagree) if not d]
    dis_x = [x for x, d in zip(xs, disagree) if d]
    dis_y = [y for y, d in zip(ys, disagree) if d]
    if agree_x:
        ax.scatter(agree_x, agree_y, s=8, alpha=0.35, color="#2166ac", zorder=2,
                   linewidths=0)
    if dis_x:
        ax.scatter(dis_x, dis_y, s=22, alpha=0.95, color="#b2182b", zorder=3,
                   marker="D", linewidths=0,
                   label=f"differing call ({len(dis_x)})")
        ax.legend(loc="upper left", fontsize=7, frameon=False)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(f"R {label}", fontsize=8)
    ax.set_ylabel(f"Rust {label}", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    # The number that actually matters, on the figure.
    worst = max((abs(a - b) for a, b in zip(xs, ys)), default=0.0)
    ax.text(0.97, 0.04, f"max|d| = {worst:.3g}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7, color="0.3")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", default="validation/realdata/results")
    ap.add_argument("--out", default="validation/realdata/results/scatter")
    args = ap.parse_args()

    data = Path(args.data)
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    rust_all = load_arm(data / "rust.jsonl")
    r_all = load_arm(data / "r.jsonl")
    summary = json.loads((data / "summary.json").read_text())
    agree_by_ds = {
        d["dataset"]: d.get("diff_abn", {}).get("agreement")
        for d in summary["datasets"]
    }

    written = []
    for name in sorted(set(rust_all) & set(r_all)):
        rust, r = rust_all[name], r_all[name]
        fig, axes = plt.subplots(1, len(PANELS) + 1, figsize=(4.0 * (len(PANELS) + 1), 3.8))
        for ax, field in zip(axes, PANELS):
            xs, ys, dis = pairs(rust, r, field)
            draw_axis(ax, xs, ys, dis, AXIS_LABEL[field], field)
        # Fifth panel: the significance indicator itself, not a continuous proxy.
        ax = axes[-1]
        ru = rust["sig"]
        rv = r["sig"]
        shared = sorted(set(ru) & set(rv))
        ax.scatter([1.0 if rv[t] else 0.0 for t in shared],
                   [1.0 if ru[t] else 0.0 for t in shared],
                   s=10, alpha=0.4, color="#2166ac", linewidths=0)
        ax.set_xlim(-0.4, 1.4)
        ax.set_ylim(-0.4, 1.4)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["not sig", "sig"], fontsize=7)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["not sig", "sig"], fontsize=7)
        ax.set_xlabel("R call", fontsize=8)
        ax.set_ylabel("Rust call", fontsize=8)
        # The agreement in the title is computed from the plotted points, not read
        # from `summary.json`. Reading it from the summary meant a figure could
        # caption itself "agree 1.0" while showing three mismatched calls, because
        # the summary and the arms behind it were free to disagree -- which is
        # exactly what a mutation test does.
        n_agree = sum(1 for t in shared if ru[t] == rv[t])
        frac = n_agree / len(shared) if shared else float("nan")
        title = f"significance  ({n_agree}/{len(shared)} agree = {frac:.4f})"
        recorded = agree_by_ds.get(name)
        if recorded is not None and abs(recorded - frac) > 1e-12:
            # Not a warning: the figure and the recorded summary disagree, and
            # saying so on the figure is more use than picking one silently.
            title += f"\nsummary.json says {recorded}"
        ax.set_title(title, fontsize=9)
        ax.tick_params(labelsize=7)

        fig.suptitle(
            f"{name}: Rust vs R on {rust['coefficient']}  "
            f"(n = {len(shared)} taxa; every point on the dashed line is an exact match)",
            fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        path = outdir / f"{name}.svg"
        fig.savefig(path, format="svg")
        plt.close(fig)
        written.append(path)
        print(f"wrote {path}")

    if not written:
        print("no datasets in both arms; nothing drawn", file=sys.stderr)
        return 1
    return 0


import sys  # noqa: E402  (used only in main's error path)

if __name__ == "__main__":
    raise SystemExit(main())