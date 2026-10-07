#!/usr/bin/env python3
"""Derive validation/exact/cases.json from validation/exact/fixtures.json.

One source of truth. `fixtures.json` says how to build each input and what the
installed original actually did; `cases.json` says which of those the release
gates require. Writing them separately is how a case ends up in the manifest but
not in the suite, or in both with different outcomes.

The `small` case set is chosen by tag, not by a hand-written list, so adding a
case to `fixtures.json` with a `small` tag puts it in CI without anyone editing a
second file. The `required` set is every case.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SPEC = REPO / "validation" / "exact" / "fixtures.json"
OUT = REPO / "validation" / "exact" / "cases.json"

# Tags that make a case cheap enough for a push. Chosen for cost, not for
# convenience: the lmer paths and the fx03-scale cases are deliberately absent.
SMALL_TAGS = {"small"}


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spec", default=str(SPEC))
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--check", action="store_true",
                    help="fail if the existing cases.json is stale instead of rewriting it")
    args = ap.parse_args(argv)

    spec = json.loads(Path(args.spec).read_text())
    cases = spec["cases"]

    ids = [c["id"] for c in cases]
    if len(set(ids)) != len(ids):
        dupes = sorted({i for i in ids if ids.count(i) > 1})
        print(f"duplicate case ids in the fixture manifest: {', '.join(dupes)}",
              file=sys.stderr)
        return 1
    for c in cases:
        if "expect" not in c or "outcome" not in c["expect"]:
            print(f"case {c['id']} has no expect.outcome; every case must declare one",
                  file=sys.stderr)
            return 1

    def entry(c: dict, inputs_rel: str) -> dict:
        return {
            "id": c["id"],
            "input": f"inputs/{c['id']}.rds",
            "fn": c.get("fn", "ancombc2"),
            "call_style": c.get("call_style", "named"),
            "expect": {"outcome": c["expect"]["outcome"]},
            "tags": sorted(c.get("tags") or []),
        }

    small = [entry(c, "inputs") for c in cases
             if SMALL_TAGS & set(c.get("tags") or [])]
    required = [entry(c, "inputs") for c in cases]

    doc = {
        "schema": "ancombc2-exact-cases/1",
        "generated_by": "scripts/make_case_manifest.py",
        "generated_from": str(Path(args.spec).relative_to(REPO)),
        "note": (
            "Every case both arms must execute, and the outcome each must produce. "
            "A case absent from an observed run is a failed gate, not a skip: "
            "reference errors are comparable outcomes and are declared with "
            "expect.outcome = error. This file is generated from fixtures.json so the "
            "manifest and the suite cannot describe different case sets; edit the "
            "fixture manifest, not this file. Removing a case requires a recorded "
            "rationale in IMPLEMENTATION_STATUS.md and a new campaign."),
        "case_sets": {
            "small": {
                "description": ("Strict cases for normal CI. Selected by the `small` "
                                "tag in fixtures.json: no lmer path, no fx03-scale "
                                "case. Broad enough that a schema or default change "
                                "breaks them."),
                "allow_extra": False,
                "cases": small,
            },
            "required": {
                "description": "Every case in the fixture manifest.",
                "allow_extra": False,
                "cases": required,
            },
        },
    }

    text = json.dumps(doc, indent=2) + "\n"
    out = Path(args.out)
    if args.check:
        if not out.is_file():
            print(f"{out} does not exist; run without --check", file=sys.stderr)
            return 1
        if out.read_text() != text:
            print(f"{out} is stale; run scripts/make_case_manifest.py without --check",
                  file=sys.stderr)
            return 1
        print(f"{out.relative_to(REPO)} is up to date ({len(small)} small, "
              f"{len(required)} required)")
        return 0

    out.write_text(text)
    print(f"wrote {out.relative_to(REPO)} ({len(small)} small, {len(required)} required)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
