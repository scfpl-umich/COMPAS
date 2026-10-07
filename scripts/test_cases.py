#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""
COMPAS Test Builder (Serial / MPI Version - HPC Portable)

Features:
  - Verbose by default: print the commands, log files and timings of each step
    (--quiet prints only the progress and the summary)
  - Auto log tail on failure
  - Filter test cases using --only
  - Run-time variants of each case from <case>/variants.txt (skip with --no-variants)
  - CSV summary output
  - Pass-through args to make using `--` (after -j, so that -jN limits the build jobs)
  - MPI ranks via --ranks (default 1, which runs without mpirun)
  - Timeout (portable, no external `timeout`)
  - GitHub Actions annotations
  - Keep the output of the last run with --dont-clean

scripts/test_cases.sh takes the same options and has the same defaults.

Every directory in exec/_Tests that contains a GNUmakefile is tested. The cases listed in
exec/_Tests/test_list.txt run first, in that order, and any other case directory runs after
them with a warning. A run passes if it finishes, or if it is still running without error
when the timeout is reached. A run that prints NaN fails in either case.

A case may contain variants.txt. Each non-comment line is a variant name followed by
ParmParse overrides appended to the command line of the base run. Output is cleaned before
every run, except that a variant whose name starts with "+" continues from the output of
the previous run, e.g. to restart from its checkpoint. See scripts/test_cases.sh.

A case directory may hold several builds and several base runs. In variants.txt a line
"[LABEL] ARGS" starts a section: the base run LABEL, built with "make ARGS" (e.g. BUILD=6Eq)
and run with the inputs file inputs=FILE among ARGS (default prob/inputs); the variant lines
below it use the same executable and inputs. Sections with the same make arguments share one
build, and each build must have its own executable (USERSuffix in the GNUmakefile). Without
section lines a case is one build, run with prob/inputs, as before.
"""

import os
import re
import sys
import csv
import time
import shutil
import argparse
import subprocess
from pathlib import Path
from contextlib import contextmanager

# ============================================================
# Colors (auto-disabled if not a terminal)
# ============================================================

IS_TERMINAL = sys.stdout.isatty()

def color(code, text):
    if IS_TERMINAL:
        return f"\033[{code}m{text}\033[0m"
    return text

def header(msg):  print(color("1;35", f"=== {msg} ==="))
def info(msg):    print(color("1;36", msg))
def success(msg): print(color("1;32", msg))
def warning(msg): print(color("1;33", msg))
def error(msg):   print(color("1;31", msg))

def vprint(msg, verbose):
    if verbose:
        print(color("0;36", f"[] {msg}"))

# ============================================================
# GitHub Actions helpers
# ============================================================

def gh_active():
    return os.environ.get("GITHUB_ACTIONS") == "true"

def gh_error(msg, file=None):
    if gh_active():
        if file:
            print(f"::error file={file}::{msg}")
        else:
            print(f"::error::{msg}")

def gh_warning(msg):
    if gh_active():
        print(f"::warning::{msg}")

def gh_group(title):
    if gh_active():
        print(f"::group::{title}")
    else:
        print(color("1;35", f"--- {title} ---"))

def gh_endgroup():
    if gh_active():
        print("::endgroup::")

# ============================================================
# Portable timeout wrapper
# ============================================================

def run_with_timeout(cmd, timeout, stdout=None, stderr=None):
    try:
        proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=stdout, stderr=stderr)
        try:
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.terminate()
            try:
                proc.wait(5)
            except subprocess.TimeoutExpired:
                proc.kill()
            return 124
        return proc.returncode
    except FileNotFoundError:
        return 127

# ============================================================
# Directory isolation (like bash subshell cd)
# ============================================================

@contextmanager
def pushd(path: Path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)

# ============================================================
# Find executable (*.ex) in the case directory
# ============================================================

def find_executable(path: Path):
    for p in sorted(path.glob("*.ex")):
        if p.is_file() and os.access(p, os.X_OK):
            return p
    return None

# ============================================================
# Test list and variants
# ============================================================

def load_cases(test_list: Path, test_base: Path):
    cases = []
    with open(test_list) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cases.append(line)
    for d in sorted(test_base.iterdir()):
        if d.is_dir() and (d / "GNUmakefile").exists() and d.name not in cases:
            warning(f"Case {d.name} is not in test_list.txt; running it after the listed cases")
            cases.append(d.name)
    return cases

SECTION_RE = re.compile(r"^\[([^\]\s]+)\]\s*(.*)$")

def load_sections(case_dir: Path, case: str):
    """The base runs of a case directory, in file order: a list of dicts with the label, the
    make arguments, the inputs file and the variants [(name, overrides, keep)]. Without a
    section line in variants.txt: one base run named after the case, plain make, prob/inputs."""
    path = case_dir / "variants.txt"
    sections, current = [], None
    lines = path.read_text().splitlines() if path.exists() else []
    has_sections = any(SECTION_RE.match(l.split("#", 1)[0].strip()) for l in lines)
    if not has_sections:
        current = {"label": case, "make_args": [], "inputs": "./prob/inputs", "variants": []}
        sections.append(current)
    for raw in lines:
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        m = SECTION_RE.match(line)
        if m:
            make_args, inputs = [], "./prob/inputs"
            for word in m.group(2).split():
                if word.startswith("inputs="):
                    inputs = word.split("=", 1)[1]
                else:
                    make_args.append(word)
            current = {"label": m.group(1), "make_args": make_args, "inputs": inputs, "variants": []}
            sections.append(current)
            continue
        if current is None:
            raise ValueError(f"{path}: variant line before the first [LABEL] line: {raw.strip()}")
        words = line.split()
        name, overrides = words[0], words[1:]
        keep = name.startswith("+")
        current["variants"].append((name.lstrip("+"), overrides, keep))
    return sections

EXE_RE = re.compile(r"^executable is (\S+)\s*$", re.MULTILINE)

def query_executable(make_args):
    """The executable that "make ARGS" builds, from AMReX's print-executable rule."""
    out = subprocess.run(["make", "--no-print-directory", *make_args, "print-executable"],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout
    m = EXE_RE.search(out)
    return Path(m.group(1)) if m else None

# ============================================================
# Run one executable and classify the result
# ============================================================

NAN_RE = re.compile(r"(^|[^A-Za-z])nan([^A-Za-z]|$)", re.IGNORECASE | re.MULTILINE)

def unused_overrides(contents, overrides):
    idx = contents.find("Unused ParmParse Variables")
    if idx < 0:
        return []
    unused = contents[idx:]
    keys = [o.split("=", 1)[0] for o in overrides if "=" in o]
    return [k for k in keys if f"::{k}(" in unused]

def run_one(label, executable, inputs, overrides, run_log, mpi_ranks, timeout_sec, tail_lines, verbose):
    if mpi_ranks == 1:
        cmd = [str(executable), inputs, *overrides]
    else:
        cmd = ["mpirun", "-n", str(mpi_ranks), str(executable), inputs, *overrides]

    vprint(f"Run command: {' '.join(cmd)}", verbose)
    vprint(f"Run log: {run_log}", verbose)

    start_run = time.time()
    with open(run_log, "w") as rf:
        exit_code = run_with_timeout(cmd, timeout_sec, stdout=rf, stderr=rf)
    run_time = int(time.time() - start_run)

    with open(run_log, errors="replace") as rf:
        contents = rf.read()
    unused = unused_overrides(contents, overrides)

    if NAN_RE.search(contents):
        status = "NAN"
        error(f"✖ {label}: NaN detected in output (exit={exit_code})")
    elif exit_code == 124:
        status = "TIMEOUT"
        success(f"✔ {label} reached timeout")
    elif exit_code != 0:
        status = "RUN_FAIL"
        error(f"✖ {label}: runtime failure (exit={exit_code})")
    elif unused:
        status = "UNUSED_PARAM"
        error(f"✖ {label}: override(s) never read by the code: {' '.join(unused)}")
    else:
        status = "OK"
        success(f"✔ {label} ran successfully in {run_time}s")

    if status in ("NAN", "RUN_FAIL"):
        gh_error("Runtime failure", str(run_log))
        print("--- Run log (tail) ---")
        for line in contents.splitlines()[-tail_lines:]:
            print(line)

    return status, run_time

# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    # Verbose by default
    parser.add_argument("--verbose", "-v", action="store_true", default=True,
                        help="Print the commands, log files and timings of each step (default)")
    parser.add_argument("--quiet", "-q", action="store_false", dest="verbose",
                        help="Print only the progress and the summary")
    parser.add_argument("--tail", type=int, default=40,
                        help="Lines of the log printed when a build or run fails")
    parser.add_argument("--only", type=str, default="",
                        help="Run only the cases whose directory name or run label contains this pattern")
    parser.add_argument("--ranks", type=int, default=1,
                        help="MPI ranks per run; with 1 the executable runs without mpirun")
    parser.add_argument("--timeout", type=int, default=10,
                        help="Time limit of each run, in seconds")
    parser.add_argument("--no-variants", action="store_true",
                        help="Run only the base inputs of each case", default=False)
    parser.add_argument("--dont-clean", action="store_true",
                        help="Keep the output of the last run of each case",
                        default=False)

    parser.add_argument("make_opts", nargs=argparse.REMAINDER)

    args = parser.parse_args()

    verbose = args.verbose
    tail_lines = args.tail
    only_pattern = args.only
    mpi_ranks = args.ranks
    timeout_sec = args.timeout
    dont_clean = args.dont_clean
    run_variants = not args.no_variants

    make_opts = args.make_opts
    if make_opts and make_opts[0] == "--":
        make_opts = make_opts[1:]

    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent
    test_base = repo_root / "exec" / "_Tests"
    test_list = test_base / "test_list.txt"
    clean_script = repo_root / "exec" / "clean.sh"
    log_dir = script_dir / "logs"

    log_dir.mkdir(exist_ok=True)

    # Cleanup old logs
    for f in log_dir.glob("*.log"):
        f.unlink()
    for f in script_dir.glob("*.tmp"):
        f.unlink()

    csv_file = script_dir / "test_summary.csv"
    if csv_file.exists():
        csv_file.unlink()

    if not test_list.exists():
        error(f"Test list file not found: {test_list}")
        sys.exit(1)

    def clean_output():
        subprocess.run(["bash", str(clean_script)],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # ========================================================
    # Load test list
    # ========================================================

    cases = load_cases(test_list, test_base)

    # Base runs of each case directory; --only keeps a whole directory whose name matches,
    # or the base runs whose label matches
    case_sections, owner = {}, {}
    for c in cases:
        try:
            secs = load_sections(test_base / c, c) if (test_base / c).is_dir() else \
                   [{"label": c, "make_args": [], "inputs": "./prob/inputs", "variants": []}]
        except ValueError as e:
            error(str(e))
            sys.exit(2)
        # Logs and summary rows are named by label, so a label may appear only once
        for sec in secs:
            if sec["label"] in owner:
                error(f"Label {sec['label']} is in both {owner[sec['label']]} and {c}")
                sys.exit(2)
            owner[sec["label"]] = c
        if only_pattern and only_pattern not in c:
            secs = [s for s in secs if only_pattern in s["label"]]
        if secs:
            case_sections[c] = secs
    cases = [c for c in cases if c in case_sections]

    if not cases:
        error(f"No test cases match: {only_pattern or '<none>'}")
        sys.exit(2)

    header("Test Cases to Build + Run")
    for c in cases:
        labels = [s["label"] for s in case_sections[c]]
        print(f"  - {c}" + ("" if labels == [c] else f": {' '.join(labels)}"))
    print()

    header("Configuration")
    info(
        f"Verbose={verbose} | Ranks={mpi_ranks} | "
        f"Timeout={timeout_sec}s | Only={only_pattern or '<all>'} | "
        f"Variants={run_variants} | CleanScript={not dont_clean}"
    )
    print()

    with open(csv_file, "w", newline="") as csvf:
        writer = csv.writer(csvf)
        writer.writerow([
            "case","build_status","run_status",
            "build_seconds","run_seconds","total_seconds"
        ])

    def record(row):
        with open(csv_file, "a", newline="") as csvf:
            csv.writer(csvf).writerow(row)

    # ========================================================
    # Build + Run loop
    # ========================================================

    for idx, case in enumerate(cases, start=1):

        case_dir = test_base / case
        sections = case_sections[case]

        print()
        print(color("1;34", "=" * 70))
        print(color("1;34", f"[{idx}/{len(cases)}] Test Case: {case}"))
        print(color("1;34", "=" * 70))

        if not case_dir.exists():
            error(f"Missing directory: {case_dir}")
            gh_error("Missing test directory", str(log_dir / f"{case}.log"))
            record([case,"FAILED","SKIP",0,0,0])
            continue

        gh_group(f"Building {case}")

        with pushd(case_dir):

            vprint(f"Case directory: {case_dir}", verbose)
            vprint("Running: make realclean", verbose)

            # realclean removes every build of the directory
            subprocess.run(
                ["make", *sections[0]["make_args"], "realclean"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )

            # One build job per CPU, unless the make arguments after -- set their own -j
            ncpu = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 4)
            several_builds = len({tuple(s["make_args"]) for s in sections}) > 1

            builds = {}      # make arguments -> (build status, executable, build seconds)
            exe_owner = {}   # executable -> make arguments that built it

            # Each build must write its own executable. Ask make before building anything:
            # checking after a build is too late, the build has already overwritten the
            # executable of an earlier build that later sections still run.
            planned = {}
            for k in dict.fromkeys(tuple(s["make_args"]) for s in sections):
                e = query_executable([*k, *make_opts])
                if e is not None:
                    planned.setdefault(e.name, []).append(k)
            for name, keys in planned.items():
                if len(keys) > 1:
                    error(f"✖ {' and '.join('make ' + ' '.join(k) for k in keys)} write the same "
                          f"executable {name}; give each build its own USERSuffix in the GNUmakefile")
                    for k in keys:
                        builds[k] = ("NOEXE", None, 0)

            for sec in sections:

                label = sec["label"]
                key = tuple(sec["make_args"])
                build_time = 0

                # ================= BUILD (once per make arguments) =================

                if key not in builds:
                    log_file = log_dir / f"{label}.log"
                    start_build = time.time()
                    build_cmd = ["make", f"-j{ncpu}", *sec["make_args"], *make_opts]

                    vprint(f"Build command: {' '.join(build_cmd)}", verbose)
                    vprint(f"Build log: {log_file}", verbose)

                    with open(log_file, "w") as lf:
                        ret = subprocess.run(
                            build_cmd,
                            stdout=lf,
                            stderr=lf
                        )

                    build_time = int(time.time() - start_build)

                    vprint(f"Build return code: {ret.returncode}", verbose)
                    vprint(f"Build duration: {build_time}s", verbose)

                    executable = None
                    if ret.returncode != 0:
                        error(f"✖ Build failed ({build_time}s)")
                        gh_error("Build failed", str(log_file))

                        print("--- Build log (tail) ---")
                        with open(log_file) as lf:
                            for line in lf.readlines()[-tail_lines:]:
                                print(line.rstrip())
                        builds[key] = ("FAILED", None, build_time)
                    else:
                        success(f"✔ Built in {build_time}s")

                        # ================= LOCATE EXECUTABLE =================

                        vprint("Asking make for the executable (print-executable)...", verbose)
                        executable = query_executable([*sec["make_args"], *make_opts])
                        if executable is None or not executable.is_file():
                            vprint("Searching for executable (*.ex)...", verbose)
                            executable = None if several_builds else find_executable(Path("."))

                        if not executable:
                            error("No executable (*.ex) found")
                            builds[key] = ("NOEXE", None, build_time)
                        else:
                            executable = executable.resolve()
                            vprint(f"Executable found: {executable}", verbose)
                            if executable in exe_owner and exe_owner[executable] != key:
                                error(f"✖ make {' '.join(key)} and make {' '.join(exe_owner[executable])} "
                                      f"write the same executable {executable.name}; give each build "
                                      "its own USERSuffix in the GNUmakefile")
                                builds[key] = ("NOEXE", None, build_time)
                            else:
                                exe_owner[executable] = key
                                builds[key] = ("OK", executable, build_time)

                    # Keep the disk use to one build's object files
                    if several_builds:
                        shutil.rmtree("tmp_build_dir", ignore_errors=True)

                bstatus, executable, _ = builds[key]
                if bstatus == "FAILED":
                    record([label,"FAILED","SKIP",build_time,0,build_time])
                    continue
                if bstatus == "NOEXE":
                    record([label,"OK","NOEXE",build_time,0,build_time])
                    continue

                # ================= RUN (base inputs) =================

                clean_output()
                status, run_time = run_one(label, executable, sec["inputs"], [],
                                           log_dir / f"{label}_run.log",
                                           mpi_ranks, timeout_sec, tail_lines, verbose)
                record([label,"OK",status,build_time,run_time,build_time+run_time])

                # ================= RUN (variants) =================

                if run_variants:
                    for name, overrides, keep in sec["variants"]:
                        if not keep:
                            clean_output()
                        vlabel = f"{label}:{name}"
                        vlog = log_dir / f"{label}__{name}_run.log"
                        status, run_time = run_one(vlabel, executable, sec["inputs"], overrides, vlog,
                                                   mpi_ranks, timeout_sec, tail_lines, verbose)
                        record([vlabel,"OK",status,0,run_time,run_time])

            # ================= CLEANUP =================

            vprint("Starting cleanup phase...", verbose)

            if not dont_clean:
                vprint(f"Executing {clean_script}", verbose)
                clean_output()

            vprint("Running final make realclean", verbose)
            subprocess.run(
                ["make", *sections[0]["make_args"], "realclean"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )

            vprint("Cleanup complete", verbose)

        gh_endgroup()

    # ========================================================
    # Summary
    # ========================================================

    header("Summary")

    passed = []
    failed = []
    total = 0

    with open(csv_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            total += int(row["total_seconds"])
            if (row["build_status"] == "OK" and
                row["run_status"] in ("OK","TIMEOUT")):
                passed.append(
                    f"{row['case']} "
                    f"({row['run_status']}, {row['total_seconds']}s)"
                )
            else:
                failed.append(
                    f"{row['case']} "
                    f"({row['run_status']}, {row['total_seconds']}s)"
                )

    if verbose:
        info("Detailed timing breakdown:")
        with open(csv_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                print(color("0;38",
                    f"{row['case']:40} "
                    f"Build={row['build_seconds']:>4}s  "
                    f"Run={row['run_seconds']:>4}s  "
                    f"Total={row['total_seconds']:>4}s  "
                    f"Status={row['run_status']}"
                ))
        print()

    success(f"Passed: {len(passed)}")
    for p in passed:
        print(color("1;32", f"  {p}"))
    print()

    if failed:
        error(f"Failed: {len(failed)}")
        for f in failed:
            print(color("1;31", f"  {f}"))
        gh_warning("Some test cases failed.")
        print()
        warning(f"Total time: {total}s")
        warning(f"Summary CSV: {csv_file}")
        sys.exit(1)

    success(f"Total time: {total}s")
    success(f"Summary CSV: {csv_file}")

if __name__ == "__main__":
    main()
