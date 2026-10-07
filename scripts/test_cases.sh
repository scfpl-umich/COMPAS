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
# SEVERAL BUILDS OR BASE RUNS IN ONE CASE
#
# A line "[LABEL] ARGS" in variants.txt starts a section: the base run LABEL, built
# with "make ARGS" and run with the inputs file given by inputs=FILE among ARGS
# (default prob/inputs). The variant lines below it belong to it, for example
#
#     [Couette2Layer-6Eq]   BUILD=6Eq   inputs=prob/inputs.Couette2Layer-6Eq
#     no-relaxation         max_step=20 Physics.pressure_relaxation=0
#
# Sections with the same make arguments share one build. Each build must write its
# own executable (USERSuffix in the GNUmakefile); the runner asks make for it
# ("make ARGS print-executable"). Without section lines a case is one build, run
# with prob/inputs and named after its directory. --only matches the directory
# name or the labels.
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
# Sections of a case: "[LABEL] ARGS" lines of variants.txt
#   parse_sections VARIANTS_FILE DEFAULT_LABEL fills
#   S_LABEL[s], S_ARGS[s] (make arguments), S_INPUTS[s], and per variant line
#   V_SEC[j] (its section), V_NAME[j], V_PLUS[j], V_OVR[j]
#############################################
parse_sections() {
    S_LABEL=(); S_ARGS=(); S_INPUTS=(); V_SEC=(); V_NAME=(); V_PLUS=(); V_OVR=()
    local f="$1" raw line w lbl args inp n plus words=()
    if ! { [ -f "$f" ] && grep -Eq '^[[:space:]]*\[[^]]+\]' "$f"; }; then
        S_LABEL=("$2"); S_ARGS=(""); S_INPUTS=("./prob/inputs")
    fi
    [ -f "$f" ] || return 0
    while IFS= read -r raw || [ -n "$raw" ]; do
        line="$(printf '%s' "$raw" | tr -d '\r' | sed 's/#.*//;s/^[[:space:]]*//;s/[[:space:]]*$//')"
        [ -z "$line" ] && continue
        case "$line" in
            \[*\]*)
                lbl="${line#\[}"; lbl="${lbl%%\]*}"
                read -r -a words <<< "${line#*\]}"
                args=""; inp="./prob/inputs"
                for w in "${words[@]}"; do
                    case "$w" in inputs=*) inp="${w#inputs=}" ;; *) args="${args:+${args} }$w" ;; esac
                done
                S_LABEL+=("${lbl}"); S_ARGS+=("${args}"); S_INPUTS+=("${inp}")
                continue ;;
        esac
        if [ ${#S_LABEL[@]} -eq 0 ]; then
            error "${f}: variant line before the first [LABEL] line: ${line}"
            return 1
        fi
        read -r -a words <<< "$line"
        n="${words[0]}"; plus=0
        case "$n" in +*) plus=1; n="${n#+}" ;; esac
        V_SEC+=($(( ${#S_LABEL[@]} - 1 ))); V_NAME+=("$n"); V_PLUS+=("$plus"); V_OVR+=("${words[*]:1}")
    done < "$f"
}

#############################################
# Apply --only filter (safe substring match on the directory name or on the labels);
# CASE_SEL[i] is "*" (every base run of case i) or the labels to run
#############################################
FILTERED=(); CASE_SEL=(); ALL_LABELS=()
for c in "${CASES[@]}"; do
    parse_sections "${TEST_BASE}/${c}/variants.txt" "${c}" || exit 2
    # Logs and summary rows are named by label, so a label may appear only once
    for l in "${S_LABEL[@]}"; do
        for x in "${ALL_LABELS[@]}"; do
            if [ "${x%%	*}" = "$l" ]; then
                error "Label ${l} is in both ${x#*	} and ${c}"
                exit 2
            fi
        done
        ALL_LABELS+=("${l}	${c}")
    done
    if [ -z "${ONLY_PATTERN}" ] || [[ "$c" == *"$ONLY_PATTERN"* ]]; then
        FILTERED+=("$c"); CASE_SEL+=("*")
        continue
    fi
    sel=""
    for l in "${S_LABEL[@]}"; do
        [[ "$l" == *"$ONLY_PATTERN"* ]] && sel="${sel:+${sel} }${l}"
    done
    if [ -n "${sel}" ]; then FILTERED+=("$c"); CASE_SEL+=("${sel}"); fi
done
CASES=("${FILTERED[@]}")

N=${#CASES[@]}
if [ $N -eq 0 ]; then
    error "No test cases match: ${ONLY_PATTERN:-<none>}"
    exit 2
fi

header "Test Cases to Build + Run"
i=0
for c in "${CASES[@]}"; do
    parse_sections "${TEST_BASE}/${c}/variants.txt" "${c}" 2>/dev/null
    if [ "${CASE_SEL[$i]}" != "*" ]; then echo "  - $c: ${CASE_SEL[$i]}"
    elif [ ${#S_LABEL[@]} -eq 1 ] && [ "${S_LABEL[0]}" = "$c" ]; then echo "  - $c"
    else echo "  - $c: ${S_LABEL[*]}"; fi
    i=$((i+1))
done
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
#   run_one <label> <run_log> <build_seconds> <inputs> [overrides...]
# Appends one row to the CSV file.
#############################################
run_one() {
    local label="$1" run_log="$2" btime="$3" inputs="$4"
    shift 4

    local rstart rend rtime exit_code status cmd
    if [ "$MPI_RANKS" -eq 1 ]; then
        cmd=("${executable}" "${inputs}" "$@")
    else
        cmd=(mpirun -n "${MPI_RANKS}" "${executable}" "${inputs}" "$@")
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
# Build + Run a single case: each of its builds once, then its base runs and variants
#   build_case <case> <index> <total> <selected labels or "*">
#############################################
build_case() {
    local case="$1" idx="$2" total="$3" sel="$4"
    local case_dir="${TEST_BASE}/${case}"

    echo -e "${BLUE}[${idx}/${total}]${RESET} ${CYAN}${case}${RESET}"

    if [ ! -d "${case_dir}" ]; then
        error "Missing directory: ${case_dir}"
        gh_error "Missing test directory: ${case_dir}" "${LOG_DIR}/${case}.log"
        echo "${case},FAILED,SKIP,0,0,0" >> "${CSV_FILE}"
        return
    fi

    gh_group "Building ${case}"

    (
        cd "${case_dir}" || { error "cd failed"; exit 1; }

        if ! parse_sections variants.txt "${case}"; then
            echo "${case},FAILED,SKIP,0,0,0" >> "${CSV_FILE}"
            exit 0
        fi
        local ns=${#S_LABEL[@]} nv=${#V_NAME[@]} s j b
        local nbuilds
        nbuilds=$(printf '%s\n' "${S_ARGS[@]}" | sort -u | wc -l | tr -d ' ')

        vprint "Case directory: ${case_dir}"
        vprint "Running: make realclean"
        read -r -a margv <<< "${S_ARGS[0]}"
        make "${margv[@]}" realclean > /dev/null 2>&1 || true

        BUILT_ARGS=(); BUILT_EXE=(); BUILT_STATUS=()
        # Each build must write its own executable: ask make before building anything (after
        # a build is too late, it has overwritten the executable of an earlier build)
        local pa=() pe=() e
        s=0
        while [ $s -lt $ns ]; do
            b=0; e=""
            while [ $b -lt ${#pa[@]} ]; do [ "${pa[$b]}" = "${S_ARGS[$s]}" ] && e=seen; b=$((b+1)); done
            if [ -z "$e" ]; then
                read -r -a margv <<< "${S_ARGS[$s]}"
                e="$(make --no-print-directory "${margv[@]}" "${MAKE_OPTS[@]}" print-executable 2>/dev/null \
                     | sed -n 's/^executable is \([^ ]*\) *$/\1/p' | head -1)"
                b=0
                while [ $b -lt ${#pe[@]} ]; do
                    if [ -n "$e" ] && [ "${pe[$b]}" = "$e" ]; then
                        error "✖ make ${pa[$b]} and make ${S_ARGS[$s]} write the same executable ${e}; give each build its own USERSuffix"
                        BUILT_ARGS+=("${pa[$b]}" "${S_ARGS[$s]}"); BUILT_EXE+=("" ""); BUILT_STATUS+=(NOEXE NOEXE)
                    fi
                    b=$((b+1))
                done
                pa+=("${S_ARGS[$s]}"); pe+=("$e")
            fi
            s=$((s+1))
        done
        s=0
        while [ $s -lt $ns ]; do
            local label="${S_LABEL[$s]}" margs="${S_ARGS[$s]}" inputs="${S_INPUTS[$s]}"
            if [ "${sel}" != "*" ] && [[ " ${sel} " != *" ${label} "* ]]; then s=$((s+1)); continue; fi

            ########################################
            # BUILD (once per set of make arguments)
            ########################################
            local btime=0 bstatus="" found=-1
            executable=""
            b=0
            while [ $b -lt ${#BUILT_ARGS[@]} ]; do
                [ "${BUILT_ARGS[$b]}" = "${margs}" ] && found=$b
                b=$((b+1))
            done
            if [ $found -ge 0 ]; then
                executable="${BUILT_EXE[$found]}"; bstatus="${BUILT_STATUS[$found]}"
            else
                local log="${LOG_DIR}/${label}.log"
                read -r -a margv <<< "${margs}"
                bstart=$(date +%s)
                vprint "Build command: make -j ${DEFAULT_JOBS} ${margs:+${margs} }${MAKE_OPTS[*]}"
                vprint "Build log: ${log}"
                if make -j "${DEFAULT_JOBS}" "${margv[@]}" "${MAKE_OPTS[@]}" > "${log}" 2>&1; then
                    bend=$(date +%s); btime=$((bend - bstart))
                    success "✔ ${label} built in ${btime}s"
                    bstatus=OK
                    ########################################
                    # LOCATE EXECUTABLE
                    ########################################
                    executable="$(make --no-print-directory "${margv[@]}" "${MAKE_OPTS[@]}" print-executable 2>/dev/null \
                                  | sed -n 's/^executable is \([^ ]*\) *$/\1/p' | head -1)"
                    [ -n "${executable}" ] && [ -f "${executable}" ] && executable="./${executable}" || executable=""
                    if [ -z "${executable}" ] && [ "${nbuilds}" -eq 1 ]; then executable="$(find_executable)"; fi
                    if [ -z "$executable" ]; then
                        error "No executable (*.ex) found"
                        bstatus=NOEXE
                    else
                        vprint "Executable found: ${executable}"
                        b=0
                        while [ $b -lt ${#BUILT_EXE[@]} ]; do
                            if [ "${BUILT_EXE[$b]}" = "${executable}" ]; then
                                error "✖ make ${BUILT_ARGS[$b]} and make ${margs} write the same executable ${executable}; give each build its own USERSuffix"
                                bstatus=NOEXE
                            fi
                            b=$((b+1))
                        done
                    fi
                else
                    bend=$(date +%s); btime=$((bend - bstart))
                    error "✖ Build failed (${btime}s)"
                    gh_error "Build failed" "${log}"

                    gh_group "Build log (tail)"
                    tail -n "${TAIL_LINES}" "${log}"
                    gh_endgroup
                    bstatus=FAILED
                fi
                # Keep the disk use to one build's object files
                [ "${nbuilds}" -gt 1 ] && rm -rf tmp_build_dir
                BUILT_ARGS+=("${margs}"); BUILT_EXE+=("${executable}"); BUILT_STATUS+=("${bstatus}")
            fi
            if [ "${bstatus}" = "FAILED" ]; then
                echo "${label},FAILED,SKIP,${btime},0,${btime}" >> "${CSV_FILE}"
                s=$((s+1)); continue
            fi
            if [ "${bstatus}" = "NOEXE" ]; then
                echo "${label},OK,NOEXE,${btime},0,${btime}" >> "${CSV_FILE}"
                s=$((s+1)); continue
            fi

            ########################################
            # RUN (base inputs)
            ########################################
            clean_output
            run_one "${label}" "${LOG_DIR}/${label}_run.log" "${btime}" "${inputs}"

            ########################################
            # RUN (variants)
            ########################################
            if [ "$RUN_VARIANTS" -eq 1 ]; then
                j=0
                while [ $j -lt $nv ]; do
                    if [ "${V_SEC[$j]}" = "$s" ]; then
                        [ "${V_PLUS[$j]}" = "1" ] || clean_output
                        read -r -a overrides <<< "${V_OVR[$j]}"
                        run_one "${label}:${V_NAME[$j]}" "${LOG_DIR}/${label}__${V_NAME[$j]}_run.log" 0 "${inputs}" "${overrides[@]}"
                    fi
                    j=$((j+1))
                done
            fi
            s=$((s+1))
        done

        if [ "$DONT_CLEAN" -eq 0 ]; then
            vprint "Executing ${CLEAN_SCRIPT}"
            clean_output
        fi
        vprint "Running final make realclean"
        read -r -a margv <<< "${S_ARGS[0]}"
        make "${margv[@]}" realclean > /dev/null 2>&1

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
    build_case "${case}" "${i}" "${N}" "${CASE_SEL[$((i-1))]}"
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
