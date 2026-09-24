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
"""

import os
import re
import sys
import csv
import time
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

def load_variants(path: Path):
    variants = []
    if not path.exists():
        return variants
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            words = line.split()
            name, overrides = words[0], words[1:]
            keep = name.startswith("+")
            variants.append((name.lstrip("+"), overrides, keep))
    return variants

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

def run_one(label, executable, overrides, run_log, mpi_ranks, timeout_sec, tail_lines, verbose):
    if mpi_ranks == 1:
        cmd = [str(executable), "./prob/inputs", *overrides]
    else:
        cmd = ["mpirun", "-n", str(mpi_ranks), str(executable), "./prob/inputs", *overrides]

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
                        help="Run only the cases whose name contains this pattern")
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

    if only_pattern:
        cases = [c for c in cases if only_pattern in c]

    if not cases:
        error(f"No test cases match: {only_pattern or '<none>'}")
        sys.exit(2)

    header("Test Cases to Build + Run")
    for c in cases:
        print(f"  - {c}")
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
        log_file = log_dir / f"{case}.log"
        run_log = log_dir / f"{case}_run.log"

        print()
        print(color("1;34", "=" * 70))
        print(color("1;34", f"[{idx}/{len(cases)}] Test Case: {case}"))
        print(color("1;34", "=" * 70))

        if not case_dir.exists():
            error(f"Missing directory: {case_dir}")
            gh_error("Missing test directory", str(log_file))
            record([case,"FAILED","SKIP",0,0,0])
            continue

        gh_group(f"Building {case}")

        with pushd(case_dir):

            # ================= BUILD =================

            start_build = time.time()

            vprint(f"Case directory: {case_dir}", verbose)
            vprint("Running: make realclean", verbose)

            subprocess.run(
                ["make", "realclean"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )

            # One build job per CPU, unless the make arguments after -- set their own -j
            ncpu = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 4)
            build_cmd = ["make", f"-j{ncpu}", *make_opts]

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

            if ret.returncode != 0:
                error(f"✖ Build failed ({build_time}s)")
                gh_error("Build failed", str(log_file))

                print("--- Build log (tail) ---")
                with open(log_file) as lf:
                    for line in lf.readlines()[-tail_lines:]:
                        print(line.rstrip())

                record([case,"FAILED","SKIP",build_time,0,build_time])
                gh_endgroup()
                continue

            success(f"✔ Built in {build_time}s")

            # ================= LOCATE EXECUTABLE =================

            vprint("Searching for executable (*.ex)...", verbose)
            executable = find_executable(Path("."))

            if not executable:
                error("No executable (*.ex) found")
                record([case,"OK","NOEXE",build_time,0,build_time])
                gh_endgroup()
                continue

            executable = executable.resolve()
            vprint(f"Executable found: {executable}", verbose)

            # ================= RUN (base inputs) =================

            clean_output()
            status, run_time = run_one(case, executable, [], run_log,
                                       mpi_ranks, timeout_sec, tail_lines, verbose)
            record([case,"OK",status,build_time,run_time,build_time+run_time])

            # ================= RUN (variants) =================

            if run_variants:
                for name, overrides, keep in load_variants(Path("variants.txt")):
                    if not keep:
                        clean_output()
                    label = f"{case}:{name}"
                    vlog = log_dir / f"{case}__{name}_run.log"
                    status, run_time = run_one(label, executable, overrides, vlog,
                                               mpi_ranks, timeout_sec, tail_lines, verbose)
                    record([label,"OK",status,0,run_time,run_time])

            # ================= CLEANUP =================

            vprint("Starting cleanup phase...", verbose)

            if not dont_clean:
                vprint(f"Executing {clean_script}", verbose)
                clean_output()

            vprint("Running final make realclean", verbose)
            subprocess.run(
                ["make", "realclean"],
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
