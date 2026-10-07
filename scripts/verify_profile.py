#!/usr/bin/env python3
"""Preflight verifier for a certified execution profile (IMPROVED_PLAN.md S01).

A profile pins one runtime in which the original and the replacement are
compared. This script proves that the runtime in front of it *is* that runtime,
before any analysis runs. It exits non-zero if any pinned fact disagrees, if the
oracle digests do not match, or if a declared dependency is absent -- so a
mismatch is a failed preflight rather than a silently different measurement.

Nothing here reads a golden, an oracle result, or the candidate's output. It
reads the profile, the pinned source, the installed package on disk, and the
running interpreter.

    python3 scripts/verify_profile.py --profile validation/exact/profiles/linux-r453-openblas.json
    python3 scripts/verify_profile.py --profile <p> --selftest
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PROBE = REPO / "scripts" / "probe_runtime.R"


# --------------------------------------------------------------------------
# result bookkeeping
# --------------------------------------------------------------------------


@dataclass
class Check:
    name: str
    ok: bool
    expected: str
    actual: str
    detail: str = ""

    def line(self) -> str:
        mark = "PASS" if self.ok else "FAIL"
        out = f"  [{mark}] {self.name}"
        if not self.ok:
            out += f"\n         expected: {self.expected}\n         actual:   {self.actual}"
        if self.detail:
            out += f"\n         {self.detail}"
        return out


@dataclass
class Report:
    profile_id: str
    checks: list[Check] = field(default_factory=list)

    def add(self, name: str, ok: bool, expected, actual, detail: str = "") -> None:
        self.checks.append(Check(name, bool(ok), _s(expected), _s(actual), detail))

    @property
    def failed(self) -> list[Check]:
        return [c for c in self.checks if not c.ok]

    def exit_code(self) -> int:
        return 1 if self.failed else 0

    def render(self) -> str:
        head = f"profile: {self.profile_id}"
        tail = f"  {len(self.checks) - len(self.failed)}/{len(self.checks)} checks passed"
        return "\n".join([head, ""] + [c.line() for c in self.checks] + ["", tail])


def _s(v) -> str:
    if isinstance(v, (list, tuple)):
        return ",".join(str(x) for x in v)
    return str(v)


# --------------------------------------------------------------------------
# digests
# --------------------------------------------------------------------------


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# A DESCRIPTION carries build-time fields that change on every reinstall. They
# say nothing about behaviour, so they are excluded from the identity digest;
# every field that can change a result is kept.
_VOLATILE_DCF = re.compile(
    r"^(Built|Packaged|NeedsCompilation|Repository|Repository|Source|"
    r"Installation|Multiarch|LinkedTo)\s*:", re.M
)


def canonical_dcf(text: str) -> str:
    return _VOLATILE_DCF.sub("", text)


def package_identity_digest(pkg_dir: Path) -> str:
    """Behavioural identity of an installed pure-R package.

    DESCRIPTION (minus build metadata), NAMESPACE and every R/*.R, hashed with
    the path recorded alongside each digest so a rename cannot be absorbed.
    """
    parts: list[str] = []
    files: list[Path] = []
    for name in ("DESCRIPTION", "NAMESPACE"):
        p = pkg_dir / name
        if p.exists():
            files.append(p)
    rd = pkg_dir / "R"
    if rd.is_dir():
        files.extend(sorted(rd.rglob("*.R")))
    for p in sorted(files):
        raw = p.read_text(encoding="utf-8", errors="replace")
        if p.name == "DESCRIPTION":
            raw = canonical_dcf(raw)
        parts.append(f"{p.relative_to(pkg_dir).as_posix()} {hashlib.sha256(raw.encode()).hexdigest()}")
    if not parts:
        return "<no files>"
    return hashlib.sha256("\n".join(parts).encode()).hexdigest()


def parse_oracle_manifest(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        digest, _, name = line.partition("  ")
        out[name.strip()] = digest.strip()
    return out


# --------------------------------------------------------------------------
# runtime probe
# --------------------------------------------------------------------------


def parse_probe(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in text.splitlines():
        if "=" not in line:
            continue
        key, _, raw = line.partition("=")
        val = raw.strip()
        if val.startswith("'") and val.endswith("'") and len(val) >= 2:
            val = val[1:-1].replace("'\\''", "'")
        out[key.strip()] = val
    return out


def run_probe(profile: dict) -> tuple[dict[str, str], str]:
    argv = list(profile["runtime"]["probe_argv"])
    env = dict(os.environ)
    env.update({str(k): str(v) for k, v in (profile["runtime"].get("env") or {}).items()})
    env.pop("R_PROFILE_USER", None)
    proc = subprocess.run(argv, capture_output=True, text=True, env=env, cwd=str(REPO))
    if proc.returncode != 0:
        raise RuntimeError(
            f"runtime probe failed ({proc.returncode}): {argv}\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    return parse_probe(proc.stdout), proc.stderr


# --------------------------------------------------------------------------
# the checks
# --------------------------------------------------------------------------


def verify(profile: dict) -> Report:
    rep = Report(profile_id=profile["id"])

    # ---- 1. pinned oracle source digests -------------------------------
    oracle = profile["oracle"]
    src = (REPO / oracle["source_dir"]).resolve()
    manifest_path = (REPO / oracle["sha256_manifest"]).resolve()
    if not src.is_dir():
        rep.add("oracle.source_dir_exists", False, oracle["source_dir"], str(src), "pinned source directory is absent")
        return rep
    if not manifest_path.is_file():
        rep.add("oracle.manifest_exists", False, oracle["sha256_manifest"], str(manifest_path))
        return rep

    want = parse_oracle_manifest(manifest_path)
    bad: list[str] = []
    missing: list[str] = []
    for name, digest in sorted(want.items()):
        p = src / name
        if not p.is_file():
            missing.append(name)
        elif sha256_file(p) != digest:
            bad.append(name)
    rep.add(
        "oracle.file_digests",
        not bad and not missing,
        f"{len(want)} files matching {oracle['sha256_manifest']}",
        f"{len(want) - len(bad) - len(missing)} match"
        + (f"; mismatched={','.join(bad)}" if bad else "")
        + (f"; absent={','.join(missing)}" if missing else ""),
        "reference/ANCOMBC must stay byte-identical to the pinned commit",
    )

    desc = (src / "DESCRIPTION").read_text()
    m = re.search(r"^Version:\s*(\S+)\s*$", desc, re.M)
    rep.add("oracle.version", bool(m) and m.group(1) == oracle["version"],
             oracle["version"], m.group(1) if m else "<absent>")

    m = re.search(r"^License:\s*(.+?)\s*$", desc, re.M)
    rep.add("oracle.license", bool(m) and m.group(1) == oracle["license"],
             oracle["license"], m.group(1) if m else "<absent>",
             "bundled upstream code must keep its licence")

    m = re.search(r"^Depends:\s*(.+?)\s*$", desc, re.M)
    rep.add("oracle.depends", bool(m) and m.group(1) == oracle["depends"],
             oracle["depends"], m.group(1) if m else "<absent>")

    # ---- 2. the running interpreter ------------------------------------
    try:
        probe, stderr = run_probe(profile)
    except RuntimeError as exc:
        rep.add("runtime.probe_runs", False, profile["runtime"]["probe_argv"][0], str(exc))
        return rep
    rep.add("runtime.probe_runs", True, "exit 0", "exit 0")

    want_r = profile["runtime"]
    rep.add("runtime.r_version", probe.get("r_version") == want_r["r_version"],
             want_r["r_version"], probe.get("r_version", "<absent>"),
             "exactness is per runtime; a different R is a different profile")
    rep.add("runtime.r_platform", probe.get("r_platform") == want_r["r_platform"],
             want_r["r_platform"], probe.get("r_platform", "<absent>"))
    rep.add("runtime.r_arch", probe.get("r_arch") == want_r["r_arch"],
             want_r["r_arch"], probe.get("r_arch", "<absent>"))

    # The pinned Depends is a floor on the interpreter, not a preference. A
    # runtime below it cannot host the installed original at all, so it must be
    # identified as unsupported here rather than at the first analysis.
    floor = _r_version_floor(oracle["depends"])
    if floor is not None:
        have = tuple(int(x) for x in probe.get("r_version", "0").split(".")[:2])
        rep.add("runtime.satisfies_pinned_depends", have >= floor[:2],
                f"R >= {floor[0]}.{floor[1]} (pinned Depends)", probe.get("r_version", "<absent>"),
                "below the floor the original cannot be installed, so this profile cannot be exact-certified")

    # ---- 3. BLAS / LAPACK identity -------------------------------------
    blas = want_r["blas"]
    actual_blas = probe.get("blas_vendor", "<absent>")
    ok_blas = blas["kind"] in Path(actual_blas).name
    if ok_blas:
        ok_blas = blas["version"] in actual_blas
    rep.add("runtime.blas", ok_blas,
             f"{blas['kind']} {blas['version']} (substring of the linked .so)",
             actual_blas,
             "P4 is a ratio between two programs; swapping BLAS would be choosing the ruler")
    lap = want_r.get("laplack")
    if lap:
        rep.add("runtime.lapack", lap["resolves_to"] in actual_blas,
                lap["resolves_to"], actual_blas,
                "LAPACK is the same OpenBLAS shared object")

    # ---- 4. controlled library path -------------------------------------
    # Entries may be repo-relative so a committed profile carries no host path.
    want_libs = [
        str(REPO / lib) if not Path(lib).is_absolute() else lib
        for lib in want_r["library_path"]
    ]
    rep.add("runtime.library_path", probe.get("libpaths", "").split("|") == want_libs,
            want_libs, probe.get("libpaths", "<absent>").split("|"),
            "an unlisted library can shadow a dependency and silently change a number")
    rep.add("runtime.rng_kind", probe.get("rng_kind") == want_r["rng_kind"],
            want_r["rng_kind"], probe.get("rng_kind", "<absent>"))
    rep.add("runtime.locale", probe.get("locale") == want_r["locale"],
            want_r["locale"], probe.get("locale", "<absent>"),
            "locale changes collation, and collation changes p.adjust ordering")

    # ---- 5. installed original ------------------------------------------
    found = probe.get("ancombc_found", "<not-found>")
    lib_dir = Path(found) if found != "<not-found>" else None
    if lib_dir is None or not lib_dir.is_dir():
        rep.add("original.installed", False, profile["oracle"]["installed_library"],
                found, "the pinned original is not installed in this profile's library")
        return rep

    rep.add("original.resolved_from_pinned_library",
            str(Path(found).resolve()).startswith(str(Path(profile["oracle"]["installed_library"]).resolve())),
            profile["oracle"]["installed_library"], found,
            "the runner must not be able to resolve the original as the candidate")
    rep.add("original.version", probe.get("ancombc_version") == oracle["version"],
            oracle["version"], probe.get("ancombc_version", "<absent>"))

    digest = package_identity_digest(lib_dir)
    rep.add("original.identity_digest", digest == oracle["installed_identity_sha256"],
            oracle["installed_identity_sha256"], digest,
            "DESCRIPTION minus build metadata, NAMESPACE and R/*.R")

    rep.add("original.depends", probe.get("ancombc_depends") == oracle["depends"],
            oracle["depends"], probe.get("ancombc_depends", "<absent>"))

    # ---- 6. declared dependencies resolve, at the recorded versions -----
    imports = [re.sub(r"\s*\(.*", "", p).strip()
               for p in re.split(r"[,\n]", oracle["imports"] or probe.get("ancombc_imports", ""))]
    imports = [p for p in imports if p]
    installed: dict[str, list[str]] = {}
    for entry in probe.get("installed_packages", "").split("|"):
        if "@" in entry:
            name, _, ver = entry.partition("@")
            installed.setdefault(name, set()).add(ver)  # type: ignore[arg-type]
    absent = sorted(p for p in imports if p not in installed)
    rep.add("original.all_imports_resolve", not absent,
            f"all {len(imports)} Imports present", f"absent={','.join(absent) or '<none>'}",
            "a missing dependency must fail before analysis, not during it")

    want_deps = oracle.get("required_dependency_versions", {})
    wrong = []
    for name, versions in sorted(want_deps.items()):
        got = installed.get(name)
        if not got:
            wrong.append(f"{name}:<absent>")
        elif not (set(versions) & got):
            wrong.append(f"{name}:{'/'.join(sorted(got))} not in {versions}")
    rep.add("original.dependency_versions", not wrong,
            f"{len(want_deps)} declared versions", "; ".join(wrong) or "all match")

    # ---- 7. host -------------------------------------------------------
    host = profile["host"]
    try:
        physical, logical = _host_cores()
    except Exception as exc:  # pragma: no cover - diagnostics only
        physical = logical = -1
        rep.add("host.cores_readable", False, "lscpu", str(exc))
    else:
        rep.add("host.cores_readable", True, "lscpu", "ok")
    if physical >= 0:
        rep.add("host.physical_cores_at_least", physical >= host["min_physical_cores"],
                f">= {host['min_physical_cores']} physical cores", f"{physical} physical / {logical} logical",
                host["core_requirement_note"])

    return rep


def _r_version_floor(depends: str) -> tuple[int, ...] | None:
    """Extract the `R (>= a.b.c)` floor from a Depends field."""
    m = re.search(r"R\s*\(\s*>=\s*([0-9]+(?:\.[0-9]+)*)\s*\)", depends or "")
    if not m:
        return None
    return tuple(int(x) for x in m.group(1).split("."))


def _host_cores() -> tuple[int, int]:
    out = subprocess.run(["lscpu"], capture_output=True, text=True, check=True).stdout
    cores = sockets = 0
    for line in out.splitlines():
        m = re.match(r"^Core\(s\) per socket:\s*(\d+)", line)
        if m:
            cores = int(m.group(1))
        m = re.match(r"^Socket\(s\):\s*(\d+)", line)
        if m:
            sockets = int(m.group(1))
    m = re.search(r"^CPU\(s\):\s*(\d+)", out, re.M)
    logical = int(m.group(1)) if m else os.cpu_count() or 0
    physical = cores * sockets if cores and sockets else logical
    return physical, logical


# --------------------------------------------------------------------------
# negative checks
# --------------------------------------------------------------------------


MUTATIONS: list[str] = [
    "wrong_oracle_digest",
    "wrong_oracle_version",
    "wrong_r_version",
    "wrong_blas",
    "wrong_package_hash",
    "wrong_dependency_version",
    "absent_dependency",
    "wrong_library_path",
    "wrong_locale",
    "absent_original",
    "insufficient_physical_cores",
]

# Which check each mutation must trip. A mutation that fails for some *other*
# reason is not evidence: a non-zero exit caused by, say, a missing manifest
# file says nothing about digest verification.
EXPECTED_FAILING: dict[str, str] = {
    "wrong_oracle_digest": "oracle.file_digests",
    "wrong_oracle_version": "oracle.version",
    "wrong_r_version": "runtime.r_version",
    "wrong_blas": "runtime.blas",
    "wrong_package_hash": "original.identity_digest",
    "wrong_dependency_version": "original.dependency_versions",
    "absent_dependency": "original.all_imports_resolve",
    "wrong_library_path": "runtime.library_path",
    "wrong_locale": "runtime.locale",
    "absent_original": "original.installed",
    "insufficient_physical_cores": "host.physical_cores_at_least",
}


def _mutate(profile: dict, kind: str, scratch: Path) -> dict:
    p = json.loads(json.dumps(profile))
    if kind == "wrong_oracle_digest":
        # A real manifest with one byte changed, not a missing file: this tests
        # digest verification, not the existence check.
        src = Path(REPO / profile["oracle"]["sha256_manifest"])
        lines = src.read_text().splitlines()
        for i, line in enumerate(lines):
            if line.strip() and not line.startswith("#"):
                digest, name = line.split("  ", 1)
                flipped = ("0" if digest[0] != "0" else "1") + digest[1:]
                lines[i] = f"{flipped}  {name}"
                break
        bad = scratch / "oracle-tampered.sha256"
        bad.write_text("\n".join(lines) + "\n")
        p["oracle"]["sha256_manifest"] = str(bad)
    elif kind == "wrong_oracle_version":
        p["oracle"]["version"] = "2.15.1"
    elif kind == "wrong_r_version":
        p["runtime"]["r_version"] = "9.9.9"
    elif kind == "wrong_blas":
        p["runtime"]["blas"] = {"kind": "mkl", "version": "2.0.0"}
    elif kind == "wrong_package_hash":
        p["oracle"]["installed_identity_sha256"] = "0" * 64
    elif kind == "wrong_dependency_version":
        p["oracle"]["required_dependency_versions"]["MASS"] = ["0.0-not-a-version"]
    elif kind == "absent_dependency":
        p["oracle"]["imports"] = p["oracle"]["imports"] + ", notAPackageXYZ"
    elif kind == "wrong_library_path":
        p["runtime"]["library_path"] = ["/nonexistent/library/a", "/nonexistent/library/b"]
    elif kind == "wrong_locale":
        p["runtime"]["locale"] = "LC_COLLATE=C;LC_NUMERIC=de_DE.UTF-8"
    elif kind == "absent_original":
        p["oracle"]["installed_library"] = "/nonexistent/library/ANCOMBC"
        for var in ("R_LIBS", "R_LIBS_USER", "R_LIBS_SITE"):
            p["runtime"]["env"][var] = "/nonexistent/library"
    elif kind == "insufficient_physical_cores":
        p["host"]["min_physical_cores"] = 1024
    else:  # pragma: no cover
        raise ValueError(kind)
    return p


def selftest(profile: dict) -> int:
    """Each mutation of the profile must produce a non-zero preflight, and must
    trip the check it was designed to exercise."""
    bad = 0
    print(f"negative checks for profile {profile['id']}")
    with tempfile.TemporaryDirectory() as td:
        scratch = Path(td)
        for kind in MUTATIONS:
            bad_p = _mutate(profile, kind, scratch)
            bad_p.pop("_drop_manifest", None)
            want = EXPECTED_FAILING[kind]
            rep = verify(bad_p)
            code = rep.exit_code()
            names = {c.name for c in rep.failed}
            ok = code != 0 and want in names
            extra = sorted(names - {want})
            print(f"  [{'PASS' if ok else 'FAIL'}] mutation={kind}: exit={code}, "
                  f"expected={want}, failing={sorted(names)}"
                  + (f"  (also {extra})" if extra else ""))
            if not ok:
                bad += 1

    rep = verify(profile)
    code = rep.exit_code()
    ok = code == 0
    print(f"  [{'PASS' if ok else 'FAIL'}] unmutated profile: exit={code}, "
          f"{len(rep.checks) - len(rep.failed)}/{len(rep.checks)} checks passed")
    if not ok:
        print(rep.render())
        bad += 1
    return 1 if bad else 0


# --------------------------------------------------------------------------


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profile", required=True, help="path to a profile JSON")
    ap.add_argument("--json-out", help="write the machine-readable report here")
    ap.add_argument("--selftest", action="store_true", help="also run the negative checks")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    profile = json.loads(Path(args.profile).read_text())
    profile.pop("_drop_manifest", None)

    if args.selftest:
        rc = selftest(profile)
    else:
        rep = verify(profile)
        if not args.quiet:
            print(rep.render())
        rc = rep.exit_code()
        if args.json_out:
            Path(args.json_out).write_text(
                json.dumps(
                    {
                        "profile": profile["id"],
                        "checks": [
                            {"name": c.name, "ok": c.ok, "expected": c.expected,
                             "actual": c.actual, "detail": c.detail}
                            for c in rep.checks
                        ],
                        "exit_code": rc,
                    },
                    indent=2,
                )
                + "\n"
            )
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
