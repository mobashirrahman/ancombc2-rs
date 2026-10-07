#!/usr/bin/env python3
"""Launch exact-comparison arms and compare their captured bytes (S02/S03).

Each arm is a *fresh R process* with an isolated library containing exactly one
ANCOMBC. Both arms read the same input RDS and call the same exported function.
This driver never loads either package itself, so it cannot become an accidental
third arm.

    python3 scripts/run_exact.py --profile <profile.json> --case <input.rds> \
        --repeats 2 [--out DIR]

Exit codes
    0  every requested arm ran and every comparison agreed
    1  a comparison disagreed, or an arm failed unexpectedly
    2  the invocation or the environment is wrong (no profile, no input)
    3  the analysis itself raised an error (a comparable outcome, not a failure
       of the runner) -- only when --allow-error is given
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_profile import REPO, package_identity_digest, run_probe  # noqa: E402

RUNNER = REPO / "scripts" / "exact_runner.R"


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


@dataclass
class ArmRun:
    arm: str
    repeat: int
    argv: list[str]
    returncode: int
    out_dir: Path
    stdout: str
    stderr: str
    outcome: str = "<none>"
    pkg_path: str = "<none>"
    pkg_version: str = "<none>"
    pkg_identity: str = "<none>"
    pkg_is_replacement: str = "<none>"
    result_sha: str | None = None
    result_len: int = 0
    error_sha: str | None = None
    error_len: int = 0
    n_conditions: int = 0

    def to_json(self) -> dict:
        return {
            "arm": self.arm,
            "repeat": self.repeat,
            "argv": self.argv,
            "returncode": self.returncode,
            "out_dir": str(self.out_dir),
            "outcome": self.outcome,
            "pkg_path": self.pkg_path,
            "pkg_version": self.pkg_version,
            "pkg_identity": self.pkg_identity,
            "pkg_is_replacement": self.pkg_is_replacement,
            "result_sha256": self.result_sha,
            "result_bytes": self.result_len,
            "error_sha256": self.error_sha,
            "error_bytes": self.error_len,
            "n_conditions": self.n_conditions,
            "stdout_tail": self.stdout[-4000:],
            "stderr_tail": self.stderr[-4000:],
        }


def read_kv(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path.is_file():
        return out
    for line in path.read_text(errors="replace").splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            out[k.strip()] = v.strip()
    return out


def identity_digest_of(pkg_path: str) -> str:
    """Digest the directory the R process actually resolved ANCOMBC from.

    One definition, in verify_profile.py. The R side only reports the path, so
    there is no second notion of "the same package" to drift.
    """
    if not pkg_path or pkg_path == "<none>":
        return "<absent>"
    p = Path(pkg_path)
    if not p.is_dir():
        return "<absent>"
    return package_identity_digest(p)


def run_arm(profile: dict, arm: str, case: Path, out_root: Path, repeat: int,
            extra_env: dict[str, str] | None = None,
            timeout: int = 7200) -> ArmRun:
    lib = Path(profile["oracle"]["installed_library"] if arm == "original"
               else profile["replacement"]["installed_library"]).resolve()
    out_dir = out_root / f"{arm}-r{repeat}"
    out_dir.mkdir(parents=True, exist_ok=True)

    rscript = profile["runtime"]["interpreter_argv"]
    argv = list(rscript) + [
        str(RUNNER),
        "--library", str(lib),
        "--arm", arm,
        "--input", str(case),
        "--out", str(out_dir),
    ]
    if arm == "candidate" and profile.get("forbid_path"):
        argv += ["--forbid-path", str(Path(profile["forbid_path"]).resolve())]

    env = dict(os.environ)
    env.update({str(k): str(v) for k, v in (profile["runtime"].get("env") or {}).items()})
    # Each arm gets its OWN library path: its own package, then the profile's
    # site library for the dependencies, and nothing else.
    #
    # This is not tidiness. The two arms are the same package name, so a shared
    # R_LIBS makes `find.package("ANCOMBC")` return whichever library comes
    # first -- and `scripts/exact_runner.R` refuses to run in that case. That
    # refusal was observed here for real when the profile's R_LIBS briefly held
    # both libraries: every candidate arm reported "the candidate resolved a
    # package that does not declare itself a replacement". The guard works; the
    # shared path was the bug.
    site = [p for p in (profile["runtime"].get("library_path") or [])
            if p not in (profile["oracle"]["installed_library"],
                         (profile.get("replacement") or {}).get("installed_library"))]
    own = str(lib)
    env["R_LIBS"] = ":".join([own, *site])
    env["R_LIBS_USER"] = own
    env["R_LIBS_SITE"] = ":".join(site)
    if extra_env:
        env.update({str(k): str(v) for k, v in extra_env.items()})

    proc = subprocess.run(argv, capture_output=True, text=True, env=env,
                          cwd=str(REPO), timeout=timeout)
    run = ArmRun(arm=arm, repeat=repeat, argv=argv, returncode=proc.returncode,
                 out_dir=out_dir, stdout=proc.stdout, stderr=proc.stderr)
    meta = read_kv(out_dir / "meta.tsv")
    run.outcome = meta.get("outcome", "<no-meta>")
    run.pkg_path = meta.get("pkg_path", "<none>")
    run.pkg_version = meta.get("pkg_version", "<none>")
    run.pkg_is_replacement = meta.get("pkg_is_replacement", "<none>")
    run.n_conditions = int(meta.get("n_conditions", "0") or 0)
    run.pkg_identity = identity_digest_of(run.pkg_path)
    rb = out_dir / "result.bytes"
    if rb.is_file() and rb.stat().st_size > 0:
        run.result_sha = sha256_file(rb)
        run.result_len = rb.stat().st_size
    eb = out_dir / "error.bytes"
    if eb.is_file() and eb.stat().st_size > 0:
        run.error_sha = sha256_file(eb)
        run.error_len = eb.stat().st_size
    return run


# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# negative checks for the runner's guards (S02 acceptance)
# ---------------------------------------------------------------------------


def _invocation(profile: dict, lib: str, arm: str, case: Path, out: Path,
                extra: list[str] | None = None) -> tuple[int, str]:
    argv = list(profile["runtime"]["interpreter_argv"]) + [
        str(RUNNER), "--library", lib, "--arm", arm,
        "--input", str(case), "--out", str(out),
    ] + (extra or [])
    env = dict(os.environ)
    env.update({str(k): str(v) for k, v in (profile["runtime"].get("env") or {}).items()})
    p = subprocess.run(argv, capture_output=True, text=True, env=env,
                       cwd=str(REPO), timeout=600)
    return p.returncode, (p.stderr or p.stdout)


def _r_eval(profile: dict, code: str, *args: str) -> None:
    """Run one R expression in the profile's interpreter.

    The expression is a single argv element: splitting it would hand R several
    separate statements and it would only run the first.
    """
    argv = list(profile["runtime"]["interpreter_argv"]) + ["-e", code, *args]
    env = dict(os.environ)
    env.update({str(k): str(v) for k, v in (profile["runtime"].get("env") or {}).items()})
    p = subprocess.run(argv, capture_output=True, text=True, cwd=str(REPO), env=env)
    if p.returncode != 0:
        raise RuntimeError(f"R failed: {p.stderr or p.stdout}")


def selftest(profile: dict, case: Path, work: Path) -> int:
    import tempfile

    lib = profile["oracle"]["installed_library"]
    bad = 0

    def check(name: str, ok: bool, detail: str = "") -> None:
        nonlocal bad
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"\n         {detail}" if detail and not ok else ""))
        if not ok:
            bad += 1

    # A library that does not exist at all.
    rc, msg = _invocation(profile, str(work / "no-such-lib"), "original", case, work / "o1")
    check("absent_library_rejected", rc != 0, f"rc={rc}: {msg[:400]}")

    # A library that exists but holds no ANCOMBC. find.package() would still
    # search the rest of .libPaths(), so this is the case where an arm can
    # silently measure a package it was not given.
    empty = work / "empty-lib"
    empty.mkdir(parents=True, exist_ok=True)
    saved = dict(profile["runtime"]["env"])
    profile["runtime"]["env"] = {**saved, "R_LIBS": str(empty), "R_LIBS_USER": str(empty),
                                 "R_LIBS_SITE": str(empty)}
    rc, msg = _invocation(profile, str(empty), "original", case, work / "o2")
    profile["runtime"]["env"] = saved
    check("library_without_ancombc_rejected", rc != 0, f"rc={rc}: {msg[:400]}")

    # The candidate arm pointed at the original: the marker guard must refuse.
    rc, msg = _invocation(profile, lib, "candidate", case, work / "o3")
    check("original_refused_as_candidate", rc != 0 and "REPLACEMENT_PROVENANCE" in msg,
          f"rc={rc}: {msg[:400]}")

    # The original arm pointed at a directory that claims to be a replacement.
    fake = work / "fake-replacement"
    (fake / "inst").mkdir(parents=True, exist_ok=True)
    (fake / "inst" / "REPLACEMENT_PROVENANCE.md").write_text("not really a replacement\n")
    rc, msg = _invocation(profile, str(fake), "original", case, work / "o4")
    check("replacement_refused_as_original", rc != 0, f"rc={rc}: {msg[:400]}")

    # A missing input file, and an input that is not the right schema.
    rc, msg = _invocation(profile, lib, "original", work / "nope.rds", work / "o5")
    check("absent_input_rejected", rc != 0, f"rc={rc}: {msg[:400]}")

    bogus = work / "bogus.rds"
    _r_eval(profile,
            "saveRDS(list(schema = 'wrong/9', case_id = 'x', fn = 'ancombc2', "
            "args = list()), commandArgs(TRUE)[1], version = 3)",
            str(bogus))
    rc, msg = _invocation(profile, lib, "original", bogus, work / "o6")
    check("wrong_input_schema_rejected", rc != 0, f"rc={rc}: {msg[:400]}")

    # An unknown arm name is not silently treated as one of the two.
    rc, msg = _invocation(profile, lib, "oracle", case, work / "o7")
    check("unknown_arm_rejected", rc != 0, f"rc={rc}: {msg[:400]}")

    # An input naming a function the package does not have.
    nofn = work / "nofn.rds"
    _r_eval(profile,
            "a <- commandArgs(TRUE); i <- readRDS(a[1]); i$fn <- 'no_such_function'; "
            "saveRDS(i, a[2], version = 3)",
            str(case), str(nofn))
    rc, msg = _invocation(profile, lib, "original", nofn, work / "o8")
    check("unknown_function_rejected", rc != 0, f"rc={rc}: {msg[:400]}")

    return 1 if bad else 0


# ---------------------------------------------------------------------------


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profile", required=True)
    ap.add_argument("--case", action="append", default=[],
                    help="input RDS; repeat for several cases")
    ap.add_argument("--input-index", default=None,
                    help="built.json from scripts/make_exact_inputs.R; every case "
                         "in it is run. Pairs with --case-set so the manifest and "
                         "the run cannot describe different suites.")
    ap.add_argument("--case-set", default=None,
                    help="only cases carrying this tag are run")
    ap.add_argument("--observed-out", default=None,
                    help="write the observed case list in the shape "
                         "scripts/check_exact.R --cases reads")
    ap.add_argument("--goldens", default=None,
                    help="directory to store the original's captures in. The "
                         "candidate arm is never given this path; see S04.")
    ap.add_argument("--original-only", action="store_true",
                    help="run only the original (goldens, baseline captures)")
    ap.add_argument("--repeats", type=int, default=2,
                    help="runs per arm per case; >=2 exposes reference "
                         "nondeterminism rather than assuming it away")
    ap.add_argument("--out", default=None, help="output root (default: a temp dir)")
    ap.add_argument("--allow-error", action="store_true",
                    help="an analysis error is a comparable outcome, not a "
                         "runner failure; compare the serialized condition")
    ap.add_argument("--selftest", action="store_true",
                    help="only run the runner's guard checks, no analysis")
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args(argv)

    profile_path = Path(args.profile)
    if not profile_path.is_file():
        print(f"no such profile: {profile_path}", file=sys.stderr)
        return 2
    profile = json.loads(profile_path.read_text())

    cases = list(args.case)
    if args.input_index:
        idx_path = Path(args.input_index)
        if not idx_path.is_file():
            print(f"no such input index: {idx_path}", file=sys.stderr)
            return 2
        doc = json.loads(idx_path.read_text())
        if doc.get("schema") != "ancombc2-exact-input-index/1":
            print(f"{idx_path}: unexpected schema {doc.get('schema')!r}", file=sys.stderr)
            return 2
        for cs in doc["cases"]:
            if args.case_set and args.case_set not in (cs.get("tags") or []):
                continue
            cases.append(cs["input"])
    if not cases:
        print("no cases selected: pass --case or --input-index", file=sys.stderr)
        return 2

    for case in cases:
        if not Path(case).is_file():
            print(f"no such case input: {case}", file=sys.stderr)
            return 2

    # Refuse to score anything in a runtime that is not the pinned one. This is
    # the same preflight S01 runs, repeated here because the arms are launched
    # separately from it and a later edit could drift.
    try:
        probe, _ = run_probe(profile)
    except RuntimeError as exc:
        print(f"preflight failed: {exc}", file=sys.stderr)
        return 2
    want = profile["runtime"]["r_version"]
    if probe.get("r_version") != want:
        print(f"preflight failed: R is {probe.get('r_version')}, profile pins {want}",
              file=sys.stderr)
        return 2

    import tempfile

    if args.selftest:
        print(f"runner guard checks for profile {profile['id']}")
        with tempfile.TemporaryDirectory(prefix="ancombc2-exact-selftest-") as td:
            return selftest(profile, Path(args.case[0]), Path(td))

    tmp = None
    if args.out:
        out_root = Path(args.out)
    else:
        tmp = tempfile.TemporaryDirectory(prefix="ancombc2-exact-")
        out_root = Path(tmp.name)
    out_root.mkdir(parents=True, exist_ok=True)

    report: dict = {"profile": profile["id"], "cases": []}
    problems: list[str] = []
    nondeterminism: list[str] = []

    observed_ids: list[str] = []
    goldens_root = Path(args.goldens) if args.goldens else None
    if goldens_root:
        goldens_root.mkdir(parents=True, exist_ok=True)

    for case in cases:
        case_p = Path(case)
        entry: dict = {"case": str(case_p),
                       "input_sha256": sha256_file(case_p),
                       "runs": [], "repeats": {}}
        per_arm: dict[str, list[ArmRun]] = {}
        arms = ["original"] if args.original_only else ["original", "candidate"]
        for arm in arms:
            if arm == "candidate" and "replacement" not in profile:
                entry["runs"].append({"arm": "candidate",
                                      "status": "not-installed-yet",
                                      "note": "S05 builds r/ANCOMBC; until then "
                                              "there is no candidate to compare"})
                continue
            runs = []
            for rep in range(1, args.repeats + 1):
                if goldens_root and arm == "original":
                    r = run_arm(profile, arm, case_p, goldens_root / case_p.stem, rep)
                else:
                    r = run_arm(profile, arm, case_p, out_root / case_p.stem, rep)
                runs.append(r)
                entry["runs"].append(r.to_json())
                if arm == "original":
                    expect = profile["oracle"]["installed_identity_sha256"]
                    if r.pkg_identity != expect:
                        problems.append(
                            f"{case_p.stem}/{arm} r{rep}: package identity "
                            f"{r.pkg_identity} != pinned {expect}")
                else:
                    if r.pkg_is_replacement != "TRUE":
                        problems.append(
                            f"{case_p.stem}/{arm} r{rep}: the candidate resolved a "
                            f"package that does not declare itself a replacement")
            per_arm[arm] = runs
            if runs:
                observed_ids.append(case_p.stem)

        # A build that changes in the middle of a campaign invalidates it. This
        # is not hypothetical: reinstalling the candidate while a 95-case
        # campaign was running removed .rlib/replacement/ANCOMBC for a moment,
        # and one candidate arm correctly refused to run with "ANCOMBC is not
        # installed". That case's row is not a result and cannot be scored.
        for arm, runs in per_arm.items():
            digests = {r.pkg_identity for r in runs if r.pkg_identity not in ("<absent>", "")}
            if len(digests) > 1:
                problems.append(
                    f"{case_p.stem}/{arm}: the installed package changed during the "
                    f"campaign ({sorted(d[:12] for d in digests)}); this row is not a "
                    f"result. Reinstalls and campaign runs must not overlap.")

        # -- reference repeatability, per arm ----------------------------
        for arm, runs in per_arm.items():
            shas = {r.result_sha or r.error_sha for r in runs}
            outcomes = {r.outcome for r in runs}
            entry["repeats"][arm] = {
                "n": len(runs),
                "outcomes": sorted(outcomes),
                "distinct_payload_sha256": sorted(x for x in shas if x),
                "self_deterministic": len(shas) <= 1 and len(outcomes) <= 1,
            }
            if len(shas) > 1 or len(outcomes) > 1:
                # The original is the reference. If it is not reproducible, every
                # later comparison is uninterpretable, so this is a diagnosed
                # failure and never a skipped case.
                msg = (f"{case_p.stem}/{arm}: {len(shas)} distinct payloads across "
                       f"{len(runs)} runs -- the ARM IS NOT REPEATABLE, so this case "
                       f"cannot certify equality")
                (nondeterminism if arm == "original" else problems).append(msg)

        # A run whose package could not be identified, or whose outcome could not
        # be read, is infrastructure failure, not a mismatch. It is reported
        # separately so it can never be mistaken for a numerical difference.
        for arm, runs in per_arm.items():
            for r in runs:
                if r.outcome == "<no-meta>":
                    entry.setdefault("infrastructure_failures", []).append({
                        "arm": arm, "repeat": r.repeat, "returncode": r.returncode,
                        "stderr_tail": r.stderr[-1000:],
                    })
                    problems.append(
                        f"{case_p.stem}/{arm} r{r.repeat}: no result was produced "
                        f"(rc={r.returncode}); see infrastructure_failures")

        # -- original vs candidate ---------------------------------------
        o, c = per_arm.get("original"), per_arm.get("candidate")
        usable = lambda runs: bool(runs) and all(r.outcome != "<no-meta>" for r in runs)
        if o and c and usable(o) and usable(c):
            op, cp = o[0].outcome, c[0].outcome
            if op != cp:
                problems.append(f"{case_p.stem}: outcome {o[0].arm}={op} but "
                                f"{c[0].arm}={cp}")
            elif op == "error":
                entry["error_bytes_match"] = (o[0].error_sha == c[0].error_sha)
                if not entry["error_bytes_match"]:
                    problems.append(f"{case_p.stem}: error conditions differ")
                else:
                    entry["outcome"] = "error-match"
            else:
                entry["result_bytes_match"] = (o[0].result_sha == c[0].result_sha)
                entry["result_bytes"] = o[0].result_len
                if not entry["result_bytes_match"]:
                    problems.append(f"{case_p.stem}: result bytes differ "
                                    f"({o[0].result_sha} vs {c[0].result_sha})")
                else:
                    entry["outcome"] = "match"

        report["cases"].append(entry)

    report["problems"] = problems
    report["nondeterminism"] = nondeterminism
    report["out_root"] = str(out_root)
    if args.goldens:
        report["golden_root"] = str(goldens_root)
    if args.json_out:
        Path(args.json_out).write_text(json.dumps(report, indent=2) + "\n")
    if args.observed_out:
        # The shape scripts/check_exact.R --cases reads. Written only from cases
        # that actually ran, so the manifest check cannot be satisfied by a list
        # somebody typed.
        Path(args.observed_out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.observed_out).write_text(json.dumps(
            {"schema": "ancombc2-exact-observed/1",
             "profile": profile["id"],
             "observed_by": "scripts/run_exact.py",
             "cases": [{"id": i} for i in dict.fromkeys(observed_ids)]},
            indent=2) + "\n")

    print(f"profile: {profile['id']}   cases: {len(cases)}   "
          f"repeats: {args.repeats}")
    for entry in report["cases"]:
        rp = entry.get("repeats", {})
        rep_txt = ", ".join(f"{a}:{'ok' if v['self_deterministic'] else 'NOT-REPEATABLE'}"
                            for a, v in rp.items())
        print(f"  {Path(entry['case']).stem}: {entry.get('outcome', 'no-comparison')} "
              f"[{rep_txt}] input={entry['input_sha256'][:12]}")
    if nondeterminism:
        print("\nNONDETERMINISM (diagnosed, not ignored):")
        for m in nondeterminism:
            print(f"  {m}")
    if problems:
        print("\nPROBLEMS:")
        for m in problems:
            print(f"  {m}")

    if tmp is not None and not problems and not nondeterminism:
        tmp.cleanup()
    elif tmp is not None:
        print(f"\ncaptures kept at {out_root}")

    return 1 if (problems or nondeterminism) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
