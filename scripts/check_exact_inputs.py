#!/usr/bin/env python3
"""Hash and verify the exact inputs (IMPROVED_PLAN.md S04).

Two jobs, one place to do the hashing:

  --write     rebuild nothing; digest every input RDS and every committed source
              file it was read from, and write validation/exact/inputs/manifest.json
  --verify    (default) rebuild every input in a temporary directory from the
              same spec and require every digest to match the manifest

Why R does not do this: R has no base sha256, so a digest computed in R and one
computed in Python are two implementations that can disagree. This script owns
hashing; scripts/make_exact_inputs.R owns construction.

Why it matters at all: the legacy goldens under validation/golden were computed
before the fixture files round-tripped through lossy text. An exact contract
cannot compare a Rust run over rounded TSV with an original run over earlier
unrounded values, so the inputs are digested and the digest is part of the case's
identity. A changed input is a different case, not a new result for the old one.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BUILDER = REPO / "scripts" / "make_exact_inputs.R"
PROFILES = REPO / "validation" / "exact" / "profiles"

DEFAULT_SPEC = REPO / "validation" / "exact" / "fixtures.json"
DEFAULT_INPUTS = REPO / "validation" / "exact" / "inputs"
DEFAULT_MANIFEST = DEFAULT_INPUTS / "manifest.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def profile_interpreter(profile_id: str | None) -> list[str]:
    """The interpreter that builds the inputs, taken from the profile.

    The input RDS is a serialized R object. A file written by one R and read by
    another is usually fine but is not guaranteed to be, so the builder's
    interpreter is recorded and the manifest records which one it was.
    """
    if not profile_id:
        return []
    p = PROFILES / f"{profile_id}.json"
    if not p.is_file():
        return []
    return list(json.loads(p.read_text())["runtime"]["interpreter_argv"])


def build(spec: Path, out_dir: Path, index: Path, profile: str | None,
          only: list[str] | None = None) -> None:
    argv = profile_interpreter(profile) + [
        str(BUILDER),
        "--spec", str(spec),
        "--out-dir", str(out_dir),
        "--index", str(index),
        "--profile", profile or "unspecified",
    ]
    if only:
        argv += ["--only", ",".join(only)]
    env = dict(os.environ)
    env.pop("R_PROFILE_USER", None)
    if profile:
        p = PROFILES / f"{profile}.json"
        if p.is_file():
            env.update({str(k): str(v)
                        for k, v in json.loads(p.read_text())["runtime"].get("env", {}).items()})
    proc = subprocess.run(argv, cwd=str(REPO), env=env, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"make_exact_inputs.R failed:\n{proc.stdout}\n{proc.stderr}")


def load_index(path: Path) -> list[dict]:
    doc = json.loads(path.read_text())
    if doc.get("schema") != "ancombc2-exact-input-index/1":
        raise RuntimeError(f"{path}: unexpected schema {doc.get('schema')!r}")
    return doc["cases"]


def digest_inputs(inputs_dir: Path, cases: list[dict]) -> dict:
    """Digest every input and every committed file it was read from."""
    out: dict[str, dict] = {}
    for cs in cases:
        p = Path(cs["input"])
        if not p.is_absolute():
            p = inputs_dir / p.name
        if not p.is_file():
            raise RuntimeError(f"case {cs['id']}: input not found at {p}")
        src = cs["source"]
        digests = {}
        for label in ("counts", "meta", "config"):
            rel = src.get(label)
            if not rel or (isinstance(rel, str) and rel == "null"):
                digests[label] = None
                continue
            f = Path(rel)
            if not f.is_absolute():
                f = REPO / f
            digests[label] = sha256_file(f) if f.is_file() else None
        out[cs["id"]] = {
            "input": str(p.relative_to(REPO)) if p.is_relative_to(REPO) else str(p),
            "input_sha256": sha256_file(p),
            "input_bytes": p.stat().st_size,
            "source": src,
            "source_sha256": digests,
            "expect": cs["expect"],
            "fn": cs["fn"],
            "call_style": cs["call_style"],
            "tags": cs["tags"],
            "n_tax": cs["n_tax"],
            "n_samp": cs["n_samp"],
            "threads": cs["threads"],
            "seed": cs["seed"],
        }
    return out


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spec", default=str(DEFAULT_SPEC))
    ap.add_argument("--inputs", default=str(DEFAULT_INPUTS))
    ap.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    ap.add_argument("--profile", default="linux-r453-openblas")
    ap.add_argument("--write", action="store_true", help="write the manifest and stop")
    ap.add_argument("--verify", action="store_true",
                    help="rebuild in a temporary directory and compare (default)")
    ap.add_argument("--only", action="append", default=[])
    ap.add_argument("--json-out")
    args = ap.parse_args(argv)

    spec = Path(args.spec)
    inputs = Path(args.inputs)
    manifest = Path(args.manifest)
    if not spec.is_file():
        print(f"no such spec: {spec}", file=sys.stderr)
        return 2

    index_path = inputs / "built.json"
    if args.write or not args.verify:
        if not index_path.is_file():
            build(spec, inputs, index_path, args.profile)
        cases = load_index(index_path)
        digests = digest_inputs(inputs, cases)
        doc = {
            "schema": "ancombc2-exact-input-manifest/1",
            "spec": str(spec.relative_to(REPO)),
            "spec_sha256": sha256_file(spec),
            "profile": args.profile,
            "builder": "scripts/make_exact_inputs.R",
            "hashed_by": "scripts/check_exact_inputs.py",
            "digest_algorithm": "sha256",
            "n_cases": len(digests),
            "cases": digests,
        }
        manifest.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n")
        print(f"wrote {manifest} ({len(digests)} cases)")
        return 0

    # ---- verify: rebuild in a temporary directory and compare -----------
    if not manifest.is_file():
        print(f"no manifest: {manifest}; run with --write first", file=sys.stderr)
        return 2
    man = json.loads(manifest.read_text())
    if man.get("schema") != "ancombc2-exact-input-manifest/1":
        print(f"{manifest}: unexpected schema {man.get('schema')!r}", file=sys.stderr)
        return 2

    problems: list[str] = []

    if sha256_file(spec) != man["spec_sha256"]:
        problems.append(
            f"{spec.name} changed since the manifest was written "
            f"({man['spec_sha256'][:12]} -> {sha256_file(spec)[:12]}); "
            f"the case set changed, so this is a different suite, not a re-run")

    # The installed original's signature is asserted here rather than trusted:
    # a positional case is pinned to the documented formal order, and if the
    # installed package's order ever differs the positional cases mean something
    # else.
    declared = json.loads(spec.read_text()).get("formal_order")
    if declared:
        argv = profile_interpreter(args.profile) + [
            "-e", "cat(paste(names(formals(ANCOMBC::ancombc2)), collapse = '\\n'), '\\n')"]
        env = dict(os.environ)
        pj = PROFILES / f"{args.profile}.json"
        if pj.is_file():
            env.update({str(k): str(v)
                        for k, v in json.loads(pj.read_text())["runtime"].get("env", {}).items()})
        got = subprocess.run(argv, capture_output=True, text=True, env=env, cwd=str(REPO))
        if got.returncode != 0:
            problems.append(f"could not read the installed ancombc2 formals: {got.stderr.strip()}")
        else:
            actual = got.stdout.split()
            if actual != list(declared):
                first = next((i for i, (a, b) in enumerate(zip(declared, actual)) if a != b), None)
                problems.append(
                    "the manifest's formal_order does not match the installed "
                    f"ancombc2: declared {len(declared)} formals, installed {len(actual)}"
                    + (f"; first difference at index {first}" if first is not None else ""))

    with tempfile.TemporaryDirectory(prefix="ancombc2-exact-regen-") as td:
        tdir = Path(td) / "inputs"
        tindex = Path(td) / "built.json"
        try:
            build(spec, tdir, tindex, args.profile)
        except RuntimeError as exc:
            print(f"regeneration failed:\n{exc}", file=sys.stderr)
            return 1
        tcases = load_index(tindex)
        tdigests = digest_inputs(tdir, tcases)

        want = man["cases"]
        for cid in sorted(set(want) | set(tdigests)):
            if cid not in want:
                problems.append(f"case {cid} exists but is not in the manifest")
                continue
            if cid not in tdigests:
                problems.append(f"case {cid} is in the manifest but did not regenerate")
                continue
            for field in ("input_sha256", "n_tax", "n_samp", "call_style", "fn"):
                if want[cid][field] != tdigests[cid][field]:
                    problems.append(
                        f"case {cid}: {field} changed on regeneration "
                        f"({want[cid][field]!r} -> {tdigests[cid][field]!r})")
            for label in ("counts", "meta", "config"):
                a = (want[cid]["source_sha256"] or {}).get(label)
                b = (tdigests[cid]["source_sha256"] or {}).get(label)
                if a != b:
                    problems.append(
                        f"case {cid}: source {label} changed since the manifest "
                        f"({str(a)[:12]} -> {str(b)[:12]})")

    print(f"profile: {args.profile}   cases: {len(man['cases'])}")
    if problems:
        print("REGENERATION: NOT REPRODUCIBLE")
        for p in problems:
            print(f"  - {p}")
        rc = 1
    else:
        print("REGENERATION: byte-identical for every input and every source file")
        rc = 0

    if args.json_out:
        Path(args.json_out).write_text(json.dumps(
            {"profile": args.profile, "n_cases": len(man["cases"]),
             "problems": problems, "exit_code": rc}, indent=2) + "\n")
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
