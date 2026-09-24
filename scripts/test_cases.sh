#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#===============================================================================
# COMPAS Test Builder (Serial / MPI Version)
#
# Features:
#   - Verbose by default: print the commands, log files and timings of each step
#     (--quiet prints only the progress and the summary)
#   - Auto log tail on failure
#   - Filter test cases using --only <pattern>
#   - Run-time variants of each case from <case>/variants.txt (skip with --no-variants)
#   - CSV summary output
#   - Pass-through args to GNUmakefile using `--`
#   - MPI ranks via --ranks (default 1, which runs without mpirun)
#   - Keep the output of the last run with --dont-clean
#
# scripts/test_cases.py takes the same options and has the same defaults.
#
# Every directory in exec/_Tests that contains a GNUmakefile is tested. The cases
# listed in exec/_Tests/test_list.txt run first, in that order, and any other case
# directory runs after them with a warning so that it can be added to the list.
#
# A run passes if it finishes, or if it is still running without error when the
# timeout is reached. A run that prints NaN fails in either case.
#
#===============================================================================
# USAGE:
#
# Build all tests:
#     ./test_cases.sh
#
# Quiet, without the commands, log files and timings:
#     ./test_cases.sh --quiet
#
# Tail more log lines on failure:
#     ./test_cases.sh --tail=100
#
# Options that take a value accept both --opt=VALUE and --opt VALUE.
#
# Filter test cases with patterns:
#     ./test_cases.sh --only Advection
#
# Parallel runs and a longer timeout:
#     ./test_cases.sh --ranks=4 --timeout=30
#
# Base runs only, without the variants:
#     ./test_cases.sh --no-variants
#
# Keep the output of the last run of each case:
#     ./test_cases.sh --only Advection-5Eq --dont-clean
#
# Pass-through arguments to GNUmakefile (everything after -- goes to make, after -j,
# so that -jN limits the build jobs):
#     ./test_cases.sh -- CXX=clang++ OPT=-O0
#     ./test_cases.sh --quiet -- DEBUG=TRUE -j4
#
#===============================================================================
# VARIANTS
#
# A case may contain a file variants.txt. Each non-comment line is a variant name
# followed by ParmParse overrides that are appended to the command line of the base
# run, for example
#
#     weno5-hll   FiniteVolume.Scheme=WENO5 Physics.RiemannSolver=HLL
#
# The executable is built once and each variant is run with it. Output is cleaned
# before every run, except that a variant whose name starts with "+" continues from
# the output of the previous run, e.g. to restart from its checkpoint.
#
#===============================================================================

set -o pipefail

#############################################
# Colors
#############################################
RED="\033[1;31m"; GREEN="\033[1;32m"; YELLOW="\033[1;33m"
BLUE="\033[1;34m"; CYAN="\033[1;36m"; MAGENTA="\033[1;35m"
RESET="\033[0m"

#############################################
# Paths
#############################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
TEST_BASE="${REPO_ROOT}/exec/_Tests"
TEST_LIST="${TEST_BASE}/test_list.txt"
CLEAN_SCRIPT="${REPO_ROOT}/exec/clean.sh"
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"

#############################################
# CLI args
#############################################
VERBOSE=1     # verbose by default
TAIL_LINES=40
ONLY_PATTERN=""
MAKE_OPTS=()
MPI_RANKS=1
TIMEOUT_SEC=10
RUN_VARIANTS=1
DONT_CLEAN=0

while [ $# -gt 0 ]; do
    case "$1" in
        # The help is the comment block at the top of this file, without the license header.
        -h|--help)     sed -n '2,/^[^#]/{2,/ABSOLUTELY NO WARRANTY/d;/^#/!d;/^#====/d;s/^# \{0,1\}//;p;}' "${BASH_SOURCE[0]}"; exit 0 ;;
        -v|--verbose)  VERBOSE=1 ;;
        -q|--quiet)    VERBOSE=0 ;;
        --tail=*)      TAIL_LINES="${1#*=}" ;;
        --tail)        TAIL_LINES="$2"; shift ;;
        --only=*)      ONLY_PATTERN="${1#*=}" ;;
        --only)        ONLY_PATTERN="$2"; shift ;;
        --ranks=*)     MPI_RANKS="${1#*=}" ;;
        --ranks)       MPI_RANKS="$2"; shift ;;
        --timeout=*)   TIMEOUT_SEC="${1#*=}" ;;
        --timeout)     TIMEOUT_SEC="$2"; shift ;;
        --no-variants) RUN_VARIANTS=0 ;;
        --dont-clean)  DONT_CLEAN=1 ;;
        --) shift; MAKE_OPTS=("$@"); break ;;
        *) echo "Unknown option: $1 (see --help)"; exit 2 ;;
    esac
    shift
done

# Build jobs: one per CPU, unless the make arguments after -- set their own -j.
DEFAULT_JOBS=$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)

#############################################
# Cleanup leftover logs
#############################################
rm -f "${LOG_DIR}"/*.log 2>/dev/null
rm -f "${SCRIPT_DIR}"/*.tmp 2>/dev/null
rm -f "${SCRIPT_DIR}/test_summary.csv" 2>/dev/null

#############################################
# Output files
#############################################
CSV_FILE="${SCRIPT_DIR}/test_summary.csv"  # Merged summary file
: > "${CSV_FILE}"
echo "case,build_status,run_status,build_seconds,run_seconds,total_seconds" > "${CSV_FILE}"

#############################################
# Check if output is a terminal (interactive) or redirected to a file
#############################################
is_terminal() {
    [ -t 1 ] && return 0 || return 1
}

#############################################
# Pretty printing
#############################################
header()  { echo -e "${MAGENTA}=== $1 ===${RESET}"; }
info() {
    if is_terminal; then
        echo -e "${CYAN}$1${RESET}"
    else
        echo "$1"  # No colors if not terminal
    fi
}
success() {
    if is_terminal; then
        echo -e "${GREEN}$1${RESET}"
    else
        echo "$1"  # No colors if not terminal
    fi
}
warning() {
    if is_terminal; then
        echo -e "${YELLOW}$1${RESET}"
    else
        echo "$1"  # No colors if not terminal
    fi
}
error() {
    if is_terminal; then
        echo -e "${RED}$1${RESET}"
    else
        echo "$1"  # No colors if not terminal
    fi
}
vprint() { [ "$VERBOSE" -eq 1 ] && info "[] $1"; return 0; }

#############################################
# GitHub Actions helpers
#############################################
gh_active() { [[ "${GITHUB_ACTIONS}" == "true" ]]; }

gh_error() { gh_active && echo "::error file=$2::$1"; }
gh_warning() { gh_active && echo "::warning::$1"; }

gh_group() {
    if gh_active; then
        echo "::group::$1"
    else
        echo -e "${MAGENTA}--- $1 ---${RESET}"
    fi
}
gh_endgroup() { gh_active && echo "::endgroup::"; }

#############################################
# Validate test list
#############################################
if [ ! -f "${TEST_LIST}" ]; then
    error "Test list file not found: ${TEST_LIST}"
    exit 1
fi

#############################################
# Load test list, then add any case directory it does not name
#############################################
CASES=()
while IFS= read -r raw || [ -n "$raw" ]; do
    line="$(printf '%s' "$raw" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
    case "$line" in ""|\#*) continue ;; esac
    CASES+=("$line")
done < "${TEST_LIST}"

for dir in "${TEST_BASE}"/*/; do
    name="$(basename "${dir}")"
    [ -f "${dir}/GNUmakefile" ] || continue
    listed=0
    for c in "${CASES[@]}"; do [ "$c" = "$name" ] && listed=1 && break; done
    if [ $listed -eq 0 ]; then
        warning "Case ${name} is not in test_list.txt; running it after the listed cases"
        CASES+=("$name")
    fi
done

#############################################
# Apply --only filter (safe substring match)
#############################################
if [ -n "${ONLY_PATTERN}" ]; then
    FILTERED=()
    for c in "${CASES[@]}"; do
        if [[ "$c" == *"$ONLY_PATTERN"* ]]; then
            FILTERED+=("$c")
        fi
    done
    CASES=("${FILTERED[@]}")
fi

N=${#CASES[@]}
if [ $N -eq 0 ]; then
    error "No test cases match: ${ONLY_PATTERN:-<none>}"
    exit 2
fi

header "Test Cases to Build + Run"
for c in "${CASES[@]}"; do echo "  - $c"; done
echo ""

#############################################
# Timeout wrapper (silent unless verbose)
#############################################
timeout_cmd() {
    if command -v timeout >/dev/null 2>&1; then
        timeout "$@"
    elif command -v gtimeout >/dev/null 2>&1; then
        gtimeout "$@"
    else
        warning "No timeout command found; running without timeout"
        shift
        "$@"   # direct execution
    fi
}

# Warn once on the terminal, since run_one sends timeout_cmd's output to the run log.
if ! command -v timeout >/dev/null 2>&1 && ! command -v gtimeout >/dev/null 2>&1; then
    warning "No timeout or gtimeout found (on macOS: brew install coreutils); runs will not be time-limited"
fi

#############################################
# Portable executable finder (*.ex)
#############################################
find_executable() {
    find . -maxdepth 1 -type f -name "*.ex" -perm -u+x 2>/dev/null | head -n 1
}

#############################################
# Run the executable once and record the result
#   run_one <label> <run_log> <build_seconds> [overrides...]
# Appends one row to the CSV file.
#############################################
run_one() {
    local label="$1" run_log="$2" btime="$3"
    shift 3

    local rstart rend rtime exit_code status cmd
    if [ "$MPI_RANKS" -eq 1 ]; then
        cmd=("${executable}" ./prob/inputs "$@")
    else
        cmd=(mpirun -n "${MPI_RANKS}" "${executable}" ./prob/inputs "$@")
    fi
    vprint "Run command: ${cmd[*]}"
    vprint "Run log: ${run_log}"

    rstart=$(date +%s)
    timeout_cmd "${TIMEOUT_SEC}" "${cmd[@]}" < /dev/null > "${run_log}" 2>&1
    exit_code=$?
    rend=$(date +%s)
    rtime=$((rend - rstart))

    # Overrides that AMReX reports as unused were never read, so the variant tested nothing
    local unused="" key
    for key in "$@"; do
        key="${key%%=*}"
        [ -z "$key" ] && continue
        if sed -n '/Unused ParmParse Variables/,$p' "${run_log}" | grep -q "::${key}("; then
            unused="${unused} ${key}"
        fi
    done

    if grep -Eqi '(^|[^A-Za-z])nan([^A-Za-z]|$)' "${run_log}"; then
        status="NAN"
    elif [[ "$exit_code" -eq 124 || "$exit_code" -eq 15 || "$exit_code" -eq 14 ]]; then
        status="TIMEOUT"
    elif [ "$exit_code" -ne 0 ]; then
        status="RUN_FAIL"
    elif [ -n "$unused" ]; then
        status="UNUSED_PARAM"
    else
        status="OK"
    fi

    case "$status" in
        OK)           success "✔ ${label} ran successfully in ${rtime}s" ;;
        TIMEOUT)      success "✔ ${label} reached timeout" ;;
        NAN)          error "✖ ${label}: NaN detected in output (exit=${exit_code})" ;;
        RUN_FAIL)     error "✖ ${label}: runtime failure (exit=${exit_code})" ;;
        UNUSED_PARAM) error "✖ ${label}: override(s) never read by the code:${unused}" ;;
    esac
    if [ "$status" = "NAN" ] || [ "$status" = "RUN_FAIL" ]; then
        gh_error "Runtime failure" "${run_log}"
        gh_group "Run log (tail)"
        tail -n "${TAIL_LINES}" "${run_log}"
        gh_endgroup
    fi

    echo "${label},OK,${status},${btime},${rtime},$((btime+rtime))" >> "${CSV_FILE}"
}

clean_output() {
    bash "${CLEAN_SCRIPT}" &>/dev/null
}

#############################################
# Build + Run a single case and its variants
#############################################
build_case() {
    local case="$1" idx="$2" total="$3"
    local case_dir="${TEST_BASE}/${case}"
    local log="${LOG_DIR}/${case}.log"
    local run_log="${LOG_DIR}/${case}_run.log"

    echo -e "${BLUE}[${idx}/${total}]${RESET} ${CYAN}${case}${RESET}"

    if [ ! -d "${case_dir}" ]; then
        error "Missing directory: ${case_dir}"
        gh_error "Missing test directory: ${case_dir}" "${log}"
        echo "${case},FAILED,SKIP,0,0,0" >> "${CSV_FILE}"
        return
    fi

    gh_group "Building ${case}"

    (
        cd "${case_dir}" || { error "cd failed"; exit 1; }

        ########################################
        # BUILD
        ########################################
        bstart=$(date +%s)

        vprint "Case directory: ${case_dir}"
        vprint "Running: make realclean"
        make realclean > /dev/null 2>&1 || true

        vprint "Build command: make -j ${DEFAULT_JOBS} ${MAKE_OPTS[*]}"
        vprint "Build log: ${log}"
        if make -j "${DEFAULT_JOBS}" "${MAKE_OPTS[@]}" > "${log}" 2>&1; then
            bend=$(date +%s)
            btime=$((bend - bstart))
            success "✔ ${case} built in ${btime}s"
        else
            bend=$(date +%s)
            btime=$((bend - bstart))
            error "✖ Build failed (${btime}s)"
            gh_error "Build failed" "${log}"

            gh_group "Build log (tail)"
            tail -n "${TAIL_LINES}" "${log}"
            gh_endgroup

            echo "${case},FAILED,SKIP,${btime},0,${btime}" >> "${CSV_FILE}"
            exit 0
        fi

        ########################################
        # LOCATE EXECUTABLE
        ########################################
        executable="$(find_executable)"
        if [ -z "$executable" ]; then
            error "No executable (*.ex) found"
            echo "${case},OK,NOEXE,${btime},0,${btime}" >> "${CSV_FILE}"
            exit 0
        fi
        vprint "Executable found: ${executable}"

        ########################################
        # RUN (base inputs)
        ########################################
        clean_output
        run_one "${case}" "${run_log}" "${btime}"

        ########################################
        # RUN (variants)
        ########################################
        if [ "$RUN_VARIANTS" -eq 1 ] && [ -f variants.txt ]; then
            while IFS= read -r raw || [ -n "$raw" ]; do
                line="$(printf '%s' "$raw" | tr -d '\r' | sed 's/#.*//;s/^[[:space:]]*//;s/[[:space:]]*$//')"
                [ -z "$line" ] && continue
                read -r -a words <<< "$line"
                vname="${words[0]}"
                overrides=("${words[@]:1}")
                if [[ "$vname" == +* ]]; then
                    vname="${vname#+}"
                else
                    clean_output
                fi
                run_one "${case}:${vname}" "${LOG_DIR}/${case}__${vname}_run.log" 0 "${overrides[@]}"
            done < variants.txt
        fi

        if [ "$DONT_CLEAN" -eq 0 ]; then
            vprint "Executing ${CLEAN_SCRIPT}"
            clean_output
        fi
        vprint "Running final make realclean"
        make realclean > /dev/null 2>&1

        exit 0
    )

    gh_endgroup
}

#############################################
# Build + Run all cases
#############################################
header "Building + Running Test Cases"
info "Verbose=${VERBOSE} | Ranks=${MPI_RANKS} | Timeout=${TIMEOUT_SEC}s | Only=${ONLY_PATTERN:-<all>} | Variants=${RUN_VARIANTS} | CleanScript=$((1 - DONT_CLEAN))"
echo ""

i=1
for case in "${CASES[@]}"; do
    build_case "${case}" "${i}" "${N}"
    i=$((i+1))
done

#############################################
# Summary
#############################################
header "Summary"

PASSED=()
FAILED=()
TOTAL=0

while IFS=',' read -r c bstat rstat btime rtime tot; do
    [[ "$c" == "case" ]] && continue
    TOTAL=$((TOTAL + tot))

    if [ "$bstat" = "OK" ] && [[ "$rstat" == "OK" || "$rstat" == "TIMEOUT" ]]; then
        PASSED+=("${c} (${rstat}, ${tot}s)")
    else
        FAILED+=("${c} (${rstat}, ${tot}s)")
    fi
done < "${CSV_FILE}"

if [ "$VERBOSE" -eq 1 ]; then
    info "Detailed timing breakdown:"
    while IFS=',' read -r c bstat rstat btime rtime tot; do
        [[ "$c" == "case" ]] && continue
        printf '%-40s Build=%4ss  Run=%4ss  Total=%4ss  Status=%s\n' "$c" "$btime" "$rtime" "$tot" "$rstat"
    done < "${CSV_FILE}"
    echo ""
fi

success "Passed: ${#PASSED[@]}"
for x in "${PASSED[@]}"; do echo -e "  ${GREEN}${x}${RESET}"; done
echo ""

if [ ${#FAILED[@]} -gt 0 ]; then
    error "Failed: ${#FAILED[@]}"
    for x in "${FAILED[@]}"; do echo -e "  ${RED}${x}${RESET}"; done
    gh_warning "Some test cases failed."

    echo ""
    warning "Total time: ${TOTAL}s"
    warning "Summary CSV: ${CSV_FILE}"
    exit 1
fi

success "Total time: ${TOTAL}s"
success "Summary CSV: ${CSV_FILE}"
