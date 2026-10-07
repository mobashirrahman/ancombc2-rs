#!/usr/bin/env python3
"""Run an R script in a certified profile's interpreter, with its environment.

    python3 scripts/with_profile_r.py --profile <profile.json> scripts/foo.R --bar

The point is that no Makefile target, log line or evidence file has to name the
interpreter, the library path or the thread environment by hand. If the profile
changes, every R invocation follows it, and a target cannot silently run in a
different runtime than the one that was certified.

A profile whose declared R version is not what actually answers is a failure
here, not a warning: an evidence file produced in the wrong runtime is worse
than no evidence file.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_profile import REPO, run_probe  # noqa: E402


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profile", required=True)
    ap.add_argument("--r-cmd", action="store_true",
                    help="the command is an `R CMD ...` invocation. `R CMD INSTALL` "
                         "must be run by the profile's own R: the R on PATH is a "
                         "different interpreter, and installing this package with "
                         "it fails with 'this R is version 4.3.3, package ANCOMBC "
                         "requires R >= 4.5.0'.")
    ap.add_argument("--set-env", action="append", default=[],
                    metavar="K=V",
                    help="override one environment variable for this run. Needed "
                         "when the command must resolve ANCOMBC from a different "
                         "library than the profile's preflight uses: the profile's "
                         "R_LIBS is the *original* arm's, so a test that exercises "
                         "the replacement has to say so explicitly rather than "
                         "silently loading the original.")
    ap.add_argument("--exec", action="store_true",
                    help="run the command directly in the profile's environment "
                         "instead of under Rscript. `R CMD INSTALL` is not an R "
                         "expression: passing it to Rscript produces "
                         "`Fatal error: cannot open file 'R'`, which looks like a "
                         "missing interpreter and is not one.")
    ap.add_argument("command", nargs=argparse.REMAINDER,
                    help="R script and its arguments")
    args = ap.parse_args(argv)

    cmd = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not cmd:
        ap.error("no R command given")
    # `R CMD INSTALL` is not an R expression and must not be prefixed with
    # Rscript. The caller says so with a leading `--`, which is consumed here.
    if args.command[:1] == ["--"]:
        cmd = args.command[1:]

    profile_path = Path(args.profile)
    if not profile_path.is_file():
        print(f"no such profile: {profile_path}", file=sys.stderr)
        return 2
    profile = json.loads(profile_path.read_text())

    try:
        probe, _ = run_probe(profile)
    except RuntimeError as exc:
        print(f"preflight failed: {exc}", file=sys.stderr)
        return 2
    want = profile["runtime"]["r_version"]
    if probe.get("r_version") != want:
        print(f"preflight failed: profile {profile['id']} pins R {want}, "
              f"the interpreter answered {probe.get('r_version')}", file=sys.stderr)
        return 2

    env = dict(os.environ)
    env.update({str(k): str(v) for k, v in (profile["runtime"].get("env") or {}).items()})
    env.pop("R_PROFILE_USER", None)
    for kv in args.set_env:
        k, _, v = kv.partition("=")
        if not k:
            ap.error(f"--set-env needs K=V, got {kv!r}")
        env[k] = v

    if args.r_cmd:
        if cmd[:1] != ["R"]:
            ap.error("--r-cmd expects the command to start with R")
        argv = list(profile["runtime"]["r_cmd_argv"]) + cmd[1:]
    elif args.exec:
        argv = cmd
    else:
        argv = list(profile["runtime"]["interpreter_argv"]) + cmd
    proc = subprocess.run(argv, cwd=str(REPO), env=env)
    return proc.returncode


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
