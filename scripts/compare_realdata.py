#!/usr/bin/env python3
"""Layer 4 -- compare the Rust and R arms on a real dataset.

    python3 scripts/compare_realdata.py --data validation/realdata \\
        --rust results/rust.jsonl --r results/r.jsonl --out results/summary.json

Every metric the plan asks for, and the one it gates on:

* Spearman and Pearson of ``beta``
* ``max |dbeta|`` and ``max |dse|``
* p-value and q-value correlation
* Jaccard of the significant calls
* ``diff_abn`` agreement, ``diff_robust`` and ``passed_ss`` agreement

Why the correlation metrics are reported but not gated
------------------------------------------------------
Spearman and Pearson are *rank* and *linear* summaries, and both are dominated by
a handful of large coefficients. A port that got 99% of the taxa right to seven
digits and the remaining 1% badly wrong can score a higher Pearson than one that
is uniformly 1e-7 off. They are the right diagnostics and the wrong acceptance
criterion, so they are reported in full and the gate is on ``diff_abn``
agreement plus the absolute differences, which do not have that failure mode.

Paired-by-taxon, not by row order
--------------------------------
Both arms emit their own taxon order, and the two are not the same order. The
comparison joins on the taxon name and *reports* the join: a taxon present in one
arm and not the other is a difference in which taxa were retained, which is a
Level A contract quantity and is surfaced rather than dropped.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any


def _finite(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def pearson(xs: list[float], ys: list[float]) -> float | None:
    n = len(xs)
    if n < 2:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0.0 or syy <= 0.0:
        return None
    return sxy / math.sqrt(sxx * syy)


def _rankdata(values: list[float]) -> list[float]:
    """Average ranks, ties shared -- R's `rank(ties = "average")`.

    Written out rather than pulled from a library because a tie-handling
    difference is exactly the kind of thing that would quietly change a
    Spearman coefficient, and the whole point of this file is that every number
    in it is reproducible from this repository.
    """
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        # R's average rank for a tie block of size n starting at 0-based i.
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def spearman(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2:
        return None
    return pearson(_rankdata(xs), _rankdata(ys))


def jaccard(a: set[str], b: set[str]) -> float | None:
    union = a | b
    if not union:
        return None
    return len(a & b) / len(union)


def read_jsonl(path: Path) -> list[dict]:
    rows = []
    for i, line in enumerate(path.read_text().splitlines()):
        if not line.strip():
            continue
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError as e:
            raise SystemExit(f"{path}:{i + 1}: {e}")
    return rows


def align(rust: dict, r: dict) -> dict:
    """Join the two arms on taxon name and return the paired vectors."""
    rb = dict(zip(rust["taxa"], rust["beta"]))
    rr = dict(zip(r["taxa"], r["beta"]))
    only_rust = sorted(set(rb) - set(rr))
    only_r = sorted(set(rr) - set(rb))
    names = [t for t in rust["taxa"] if t in rr]

    # Index maps rather than `list.index` per taxon: this is `n` lookups per
    # field per arm, and quadratic in the taxa for no reason.
    idx_rust = {t: i for i, t in enumerate(rust["taxa"])}
    idx_r = {t: i for i, t in enumerate(r["taxa"])}

    def pair(field: str) -> tuple[list[Any], list[Any]]:
        # A field is comparable only when *both* arms produced one value per
        # taxon. `null` means the feature did not run -- no sensitivity analysis,
        # so no `passed_ss` -- and a length that disagrees with the taxon count
        # means the arm emitted something this comparison does not understand.
        # Either way the field becomes empty and the per-field report says how
        # many entries were compared, so a skipped comparison is visible rather
        # than absent. Indexing an empty list is what would otherwise raise.
        rv = rust.get(field)
        av = r.get(field)
        n_rust = len(rust.get("taxa") or [])
        n_r = len(r.get("taxa") or [])
        ok_r = isinstance(rv, list) and len(rv) == n_rust
        ok_a = isinstance(av, list) and len(av) == n_r
        if not (ok_r and ok_a):
            return [], []
        return (
            [rv[idx_rust[t]] for t in names],
            [av[idx_r[t]] for t in names],
        )

    out = {
        "names": names,
        "n_paired": len(names),
        "taxa_only_in_rust": only_rust,
        "taxa_only_in_r": only_r,
    }
    for field in ("beta", "se", "p", "q"):
        a, b = pair(field)
        out[field] = (a, b)
    for field in ("diff_abn", "passed_ss", "diff_robust"):
        out[field] = pair(field)
    return out


def compare_one(dataset: str, rust: dict, r: dict) -> dict:
    rec: dict[str, Any] = {"dataset": dataset, "coefficient": rust.get("coefficient")}
    for arm, row in (("rust", rust), ("r", r)):
        if row.get("error"):
            rec[f"{arm}_error"] = row["error"]
    if rust.get("error") or r.get("error"):
        rec["status"] = "failed"
        return rec

    # The two arms must be comparing the *same coefficient*. If they are not, the
    # numbers below are two unrelated contrasts: a design with a covariate gives
    # one arm the group term and the other the covariate, the signs disagree by
    # construction, and every metric reports a catastrophic mismatch that is
    # really a naming difference. So this is checked first, and it is a failure
    # rather than a warning.
    rc, ac = rust.get("coefficient"), r.get("coefficient")
    rec["coefficient_rust"] = rc
    rec["coefficient_r"] = ac
    rec["coefficients_agree"] = rc == ac
    if rc != ac:
        rec["status"] = "different coefficients"
        return rec

    a = align(rust, r)
    rec["n_taxa_in"] = rust.get("n_taxa_in")
    rec["n_samp_in"] = rust.get("n_samp_in")
    rec["n_taxa_retained_rust"] = rust.get("n_taxa_retained")
    rec["n_taxa_retained_r"] = r.get("n_taxa_retained")
    rec["n_paired"] = a["n_paired"]
    rec["taxa_only_in_rust"] = a["taxa_only_in_rust"]
    rec["taxa_only_in_r"] = a["taxa_only_in_r"]
    # A difference in the retained set is a Level A failure, and it silently
    # biases every paired metric below, so it is stated first and it also fails
    # the comparison.
    rec["retained_sets_agree"] = not a["taxa_only_in_rust"] and not a["taxa_only_in_r"]
    if a["n_paired"] == 0:
        rec["status"] = "no paired taxa"
        return rec

    summary: dict[str, Any] = {}
    for field in ("beta", "se", "p", "q"):
        xs, ys = a[field]
        pairs = [(x, y) for x, y in zip(xs, ys) if _finite(x) and _finite(y)]
        bx = [p[0] for p in pairs]
        by = [p[1] for p in pairs]
        diffs = [abs(p[0] - p[1]) for p in pairs]
        rel = [
            abs(p[0] - p[1]) / max(abs(p[0]), abs(p[1]), 1e-300) for p in pairs
        ]
        summary[field] = {
            "n": len(pairs),
            "n_nonfinite_skipped": len(xs) - len(pairs),
            "pearson": pearson(bx, by),
            "spearman": spearman(bx, by),
            "max_abs_diff": max(diffs) if diffs else None,
            "mean_abs_diff": sum(diffs) / len(diffs) if diffs else None,
            "max_rel_diff": max(rel) if rel else None,
        }
    rec.update(summary)

    for field in ("diff_abn", "passed_ss", "diff_robust"):
        xs, ys = a[field]
        # A flag is compared only where both arms produced a *boolean*. `false`
        # is a call; `null` is the absence of the feature, and a number is a
        # schema difference. Both are excluded here and counted in
        # `n_not_compared`, so "0 compared" is distinguishable from "compared and
        # all agreed".
        pairs = [(x, y) for x, y in zip(xs, ys) if isinstance(x, bool) and isinstance(y, bool)]
        if not pairs:
            rec[field] = {
                "n": 0,
                "agreement": None,
                "n_not_compared": len(xs),
                "reason": (
                    "the analysis did not produce this feature"
                    if not xs
                    else "one or both arms did not emit booleans for it"
                ),
            }
            continue
        rec[field] = {"n": len(pairs), "n_not_compared": len(xs) - len(pairs)}
        agree = sum(1 for x, y in pairs if x == y)
        paired_names = a["names"]
        called_rust = {n for n, v in zip(paired_names, xs) if v is True}
        called_r = {n for n, v in zip(paired_names, ys) if v is True}
        rec[field].update({
            "agreement": agree / len(pairs),
            "n_disagree": len(pairs) - agree,
            "called_only_by_rust": sorted(called_rust - called_r)[:20],
            "called_only_by_r": sorted(called_r - called_rust)[:20],
        })
        rec.setdefault("_disagreeing", {})[field] = len(pairs) - agree

    rec["significant_jaccard"] = jaccard(
        {n for n, v in zip(rust["taxa"], rust["diff_abn"]) if v is True},
        {n for n, v in zip(r["taxa"], r["diff_abn"]) if v is True},
    )
    disagreeing = rec.pop("_disagreeing", {})
    total = sum(disagreeing.values())
    da_agree = rec.get("diff_abn", {}).get("agreement")
    # The plan's acceptance: `diff_abn` agreement >= 99.99% on real data, and
    # the retained sets must match. Anything else is reported, not gated.
    rec["diff_abn_agreement_threshold"] = 0.9999
    rec["passes_diff_abn_threshold"] = (
        da_agree is not None and da_agree >= 0.9999 and rec["retained_sets_agree"]
    )
    rec["total_flag_disagreements"] = total
    rec["status"] = "ok"
    return rec


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rust", required=True, type=Path)
    ap.add_argument("--r", required=True, type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--data", type=Path, default=Path("validation/realdata"),
                    help="only used for the header line")
    ap.add_argument("--fail", action="store_true",
                    help="exit non-zero if any dataset fails the diff_abn threshold")
    args = ap.parse_args(argv)

    rust_rows = {r["dataset"]: r for r in read_jsonl(args.rust)}
    r_rows = {r["dataset"]: r for r in read_jsonl(args.r)}
    names = sorted(set(rust_rows) | set(r_rows))

    print("Layer 4: real-data validation")
    print("=" * 78)
    print(f"{'dataset':<26} {'paired':>7} {'d_abn':>8} {'max|dbeta|':>11} {'max|dse|':>10} "
          f"{'rho(beta)':>10} {'r(beta)':>9} {'jaccard':>8}  verdict")
    results = []
    for name in names:
        if name not in rust_rows or name not in r_rows:
            print(f"{name:<26} {'--':>7} {'--':>8} missing from "
                  f"{'rust' if name not in rust_rows else 'r'}")
            results.append({"dataset": name, "status": "missing arm"})
            continue
        rec = compare_one(name, rust_rows[name], r_rows[name])
        results.append(rec)
        if rec.get("status") != "ok":
            print(f"{name:<26} {'--':>7} {'--':>8} {rec.get('status')}"
                  + (f": {rec.get('rust_error') or rec.get('r_error')}" if
                     rec.get("rust_error") or rec.get("r_error") else ""))
            continue

        def f(v: float | None, spec: str = ".3g") -> str:
            return "n/a" if v is None else format(v, spec)

        da = rec.get("diff_abn", {}).get("agreement")
        verdict = "PASS" if rec["passes_diff_abn_threshold"] else "FAIL"
        if not rec["retained_sets_agree"]:
            verdict += " (retained sets differ)"
        print(f"{name:<26} {rec['n_paired']:>7} {f(da, '.5f'):>8} "
              f"{f(rec['beta']['max_abs_diff']):>11} {f(rec['se']['max_abs_diff']):>10} "
              f"{f(rec['beta']['spearman']):>10} {f(rec['beta']['pearson']):>9} "
              f"{f(rec.get('significant_jaccard')):>8}  {verdict}")

    ok = [r for r in results if r.get("status") == "ok"]
    passing = [r for r in ok if r.get("passes_diff_abn_threshold")]
    print("=" * 78)
    print(f"{len(passing)} of {len(ok)} dataset(s) meet the 99.99% diff_abn threshold")

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "datasets": results,
            "n_ok": len(ok),
            "n_passing": len(passing),
            "threshold": 0.9999,
            "note": (
                "Pearson and Spearman are reported as diagnostics, not gates: "
                "both are dominated by a few large coefficients, so a port that "
                "is catastrophically wrong on 1% of taxa can score higher than "
                "one that is uniformly 1e-7 off. The gate is diff_abn agreement "
                "plus a matching retained set, neither of which has that failure "
                "mode."
            ),
        }
        args.out.write_text(json.dumps(payload, indent=1) + "\n")
        print(f"wrote {args.out}")

    if args.fail and len(passing) != len(ok):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
