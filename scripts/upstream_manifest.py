#!/usr/bin/env python3
"""Write and verify r/ANCOMBC/inst/UPSTREAM_MANIFEST.sha256 (S05).

Three facts per carried file, and the point of the file is that they can be
re-derived rather than trusted:

  * what upstream pinned (from reference/env/oracle.sha256, where it is listed),
  * what the pinned checkout at reference/ANCOMBC holds,
  * what this build holds.

Files that upstream's manifest does not cover -- man/, tests/, vignettes/,
data/, inst/CITATION -- are marked as such. Calling them "added by this build"
would be false and would quietly weaken the provenance record: they came from
upstream, the existing oracle digest list simply does not mention them.

A local file that differs from the pinned checkout is reported as a deliberate
modification and must also be named in inst/REPLACEMENT_PROVENANCE.md. This
script does not decide whether a difference is acceptable; it makes the
difference impossible to miss.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PKG = REPO / "r" / "ANCOMBC"
PINNED = REPO / "reference" / "ANCOMBC"
ORACLE_MANIFEST = REPO / "reference" / "env" / "oracle.sha256"
OUT = PKG / "inst" / "UPSTREAM_MANIFEST.sha256"
PROVENANCE = PKG / "inst" / "REPLACEMENT_PROVENANCE.md"

# Written by this build, not carried from upstream.
BUILT_BY_US = {"DESCRIPTION", "inst/REPLACEMENT_PROVENANCE.md",
               "inst/UPSTREAM_MANIFEST.sha256"}

# Compiler output. Not carried from upstream and not written by hand either: the
# package compiles C and links a Rust library *inside the source tree*, so
# `r/ANCOMBC/src/` fills with `.o` and `.so` files on every `R CMD INSTALL`.
#
# A manifest that pins them would (a) fail `--verify` after every rebuild, for
# reasons that have nothing to do with provenance, and (b) record a build
# artifact's digest as if it meant something. Both are worse than not listing
# them: the record is about which *sources* this package is made of.
#
# They are also not gitignored, deliberately. A `.gitignore` rule that hides
# build output from `git status` is convenient; a provenance script that can see
# everything is more useful, and the exclusion here is explicit and auditable.
COMPILED_SUFFIXES = (".o", ".so", ".dll", ".dylib", ".a", ".mod")

HEADER = """\
# Upstream file manifest for this replacement package.
#
# ANCOMBC 2.15.2 @ dc4febdf59badb3a8dfe0c767ef2186323c2199a
# https://github.com/FrederickHuangLin/ANCOMBC -- Artistic-2.0.
#
# Fields, tab separated:
#   path  sha256_of_this_build  sha256_of_pinned_checkout  upstream_pinned_sha256
#
# `upstream_pinned_sha256` is the value from reference/env/oracle.sha256 for the
# files that oracle manifest covers -- DESCRIPTION, NAMESPACE, NEWS, R/* and
# README.md. It reads `not-covered` for man/, tests/, vignettes/, data/ and
# inst/CITATION: upstream ships those, the existing digest list simply does not
# mention them, and calling them new would be false. For those files the pinned
# checkout is the reference and the two digests are what must agree.
#
# Compiler output (.o, .so, ...) is excluded: this package builds in the source
# tree, so those files appear and change on every install and say nothing about
# which sources the package is made of.
#
# A row where column 2 and column 3 differ is a deliberate local modification.
# Every such row must also be named in inst/REPLACEMENT_PROVENANCE.md with a
# reason. Regenerate and re-check with:
#   python3 scripts/upstream_manifest.py --write
#   python3 scripts/upstream_manifest.py --verify
#
"""


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def oracle_pins() -> dict[str, str]:
    out: dict[str, str] = {}
    for line in ORACLE_MANIFEST.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        digest, name = line.split("  ", 1)
        out[name.strip()] = digest.strip()
    return out


def rows() -> list[tuple[str, str, str, str]]:
    pins = oracle_pins()
    out = []
    for p in sorted(PKG.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(PKG).as_posix()
        if rel == OUT.relative_to(PKG).as_posix():
            continue
        if p.name.endswith(COMPILED_SUFFIXES):
            continue
        local = sha256_file(p)
        pinned_path = PINNED / rel
        pinned = sha256_file(pinned_path) if pinned_path.is_file() else "<absent>"
        upstream = pins.get(rel, "not-covered")
        if rel in BUILT_BY_US:
            upstream = "written-by-this-build"
            pinned = "n/a"
        out.append((rel, local, pinned, upstream))
    return out


def render(rs: list[tuple[str, str, str, str]]) -> str:
    lines = [HEADER.rstrip("\n"), ""]
    for rel, local, pinned, upstream in rs:
        lines.append(f"{rel}\t{local}\t{pinned}\t{upstream}")
    return "\n".join(lines) + "\n"


def verify() -> int:
    if not OUT.is_file():
        print(f"no manifest: {OUT}", file=sys.stderr)
        return 2
    recorded: dict[str, str] = {}
    for line in OUT.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) == 4:
            recorded[parts[0]] = parts[1]

    problems: list[str] = []
    current = {rel: local for rel, local, _, _ in rows()}

    for rel, local in sorted(current.items()):
        want = recorded.get(rel)
        if want is None:
            problems.append(f"{rel}: present in the package but absent from the manifest")
        elif want != local:
            problems.append(f"{rel}: manifest says {want[:12]}..., file is {local[:12]}...")

    # The pinned checkout must still be what the digests were taken from. This is
    # S01's oracle pin restated where the package is built, and it is what makes
    # "byte-identical copy" a checkable claim rather than an intention.
    pins = oracle_pins()
    for rel, digest in sorted(pins.items()):
        p = PINNED / rel
        if not p.is_file():
            problems.append(f"pinned source {rel} is missing")
        elif sha256_file(p) != digest:
            problems.append(f"pinned source {rel} no longer matches reference/env/oracle.sha256")

    modified = [rel for rel, local, pinned, _ in rows()
                if pinned not in ("<absent>", "n/a") and local != pinned]
    undocumented = [rel for rel in modified
                    if PROVENANCE.is_file()
                    and rel not in PROVENANCE.read_text()]
    if undocumented:
        problems.append(
            f"{len(undocumented)} locally modified file(s) not named in "
            f"inst/REPLACEMENT_PROVENANCE.md: {', '.join(undocumented)}")

    print(f"files in package: {len(current)}   manifest rows: {len(recorded)}   "
          f"locally modified: {len(modified)}")
    for rel in modified:
        print(f"  modified: {rel}")
    if problems:
        print("UPSTREAM MANIFEST: INVALID")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("UPSTREAM MANIFEST: consistent")
    return 0


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--verify", action="store_true")
    args = ap.parse_args(argv)
    if args.verify or not args.write:
        return verify()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(render(rows()))
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
