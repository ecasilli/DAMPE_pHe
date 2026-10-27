#!/usr/bin/env bash
# Batch driver for unfolding_smooth_pHe_period.C
# Place in ANALYSIS_pHe2026_MONTHLY, together with the ROOT macro.
set -Eeuo pipefail
shopt -s nullglob

cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

ROOT_BIN="${ROOT_BIN:-root}"
DATA_DIR="${DATA_DIR:-ROOT_FILES}"
LIVETIME_FILE="${LIVETIME_FILE:-TXT_FILES/effective_livetime.txt}"
MACRO="${MACRO:-unfolding_smooth_pHe_monthly.cpp}"
TOYS="${TOYS:-10000}"
FORCE="${FORCE:-0}"
LOG_DIR="${LOG_DIR:-LOG_FILES/unfolding}"
SUFFIX="_3sigmaLow_6sigmaUp_PSDprogr_STKcharge450_comb_vert0e7_10TeV_EPOSLHC"

usage() {
cat <<'EOF'
Usage:
  ./run_unfolding.sh all                 # All annual and monthly Load outputs
  ./run_unfolding.sh monthly             # Available monthly Load outputs
  ./run_unfolding.sh annual              # Available annual Load outputs
  ./run_unfolding.sh 2016               # All available periods in 2016
  ./run_unfolding.sh 2016 01            # One month
  ./run_unfolding.sh 2016 0             # Full year
  ./run_unfolding.sh file 2016 01 PATH  # Explicit ROOT Load output path

Environment:
  LIVETIME_FILE=TXT_FILES/effective_livetime.txt  (6 columns, MM/YYYY)
  TOYS=10000        Number of toy-MC samples (0 for fast tests)
  FORCE=1           Overwrite existing complete outputs (default: skip)
  DATA_DIR=ROOT_FILES     Folder scanned for PHe_skim_*.root
  ROOT_BIN=root      ROOT executable

Automatically recognized data filenames:
  ROOT_FILES/PHe_skim_2016_...root
  ROOT_FILES/PHe_skim_2016_01_...root
The macro always takes the EXACT input path; no changes are made to Load files.
EOF
}

[[ "$TOYS" =~ ^[0-9]+$ ]] || { echo "Invalid TOYS: $TOYS" >&2; exit 2; }
[[ "$FORCE" == 0 || "$FORCE" == 1 ]] || { echo "FORCE must be 0 or 1" >&2; exit 2; }
if [[ $# -eq 0 || ${1:-} == -h || ${1:-} == --help ]]; then usage; exit 0; fi
[[ -f "$MACRO" ]] || { echo "ERROR: missing macro: $MACRO" >&2; exit 2; }
[[ -f "$LIVETIME_FILE" ]] || { echo "ERROR: missing livetime table: $LIVETIME_FILE" >&2; exit 2; }
command -v "$ROOT_BIN" >/dev/null 2>&1 || { echo "ERROR: cannot find ROOT: $ROOT_BIN" >&2; exit 2; }
mkdir -p "$LOG_DIR" ROOT_FILES TXT_FILES

make_month() {
  local text=$1
  [[ "$text" =~ ^[0-9]{1,2}$ ]] || return 1
  local num=$((10#$text))
  (( num >= 0 && num <= 12 )) || return 1
  printf '%d' "$num"
}

run_one() {
  local year=$1 month=$2 infile=$3 tag root_out dat_out log macro_call
  [[ "$year" =~ ^[0-9]{4}$ ]] || { echo "Invalid year: $year" >&2; return 2; }
  [[ -s "$infile" ]] || { echo "ERROR: missing/empty input: $infile" >&2; return 1; }
  if (( month == 0 )); then
    tag="$year"
  else
    printf -v tag '%s_%02d' "$year" "$month"
  fi
  root_out="ROOT_FILES/unfold_results_pHe_${tag}${SUFFIX}.root"
  dat_out="TXT_FILES/flux_spectrum_pHe_${tag}${SUFFIX}.dat"
  log="${LOG_DIR}/unfold_${tag}.log"

  if [[ "$FORCE" == 0 && -s "$root_out" && -s "$dat_out" ]]; then
    echo "SKIP ${tag} (outputs already exist; FORCE=1 to rerun)"
    return 0
  fi
  # ROOT's .C(args) syntax uses C++ string literals; reject problematic quotes.
  if [[ "$infile" == *\"* || "$infile" == *\\* || "$LIVETIME_FILE" == *\"* || "$LIVETIME_FILE" == *\\* ]]; then
    echo "ERROR: paths containing double quotes/backslashes are unsupported" >&2
    return 2
  fi
  macro_call="${MACRO}(${year},${month},\"${infile}\",\"${LIVETIME_FILE}\",${TOYS})"
  echo "RUN  ${tag}: ${infile} (toys=${TOYS})"
  if "$ROOT_BIN" -l -b -q "$macro_call" >"$log" 2>&1 && \
     grep -Fq "Results saved to ${root_out}" "$log" && \
     [[ -s "$root_out" && -s "$dat_out" ]]; then
    echo "OK   ${tag} -> ${dat_out}"
  else
    echo "FAIL ${tag} (see ${log})" >&2
    tail -n 12 "$log" >&2 || true
    return 1
  fi
}

# Explicit mode supports filenames other than PHe_skim_YYYY_MM_....root.
if [[ "$1" == file ]]; then
  [[ $# -eq 4 ]] || { usage >&2; exit 2; }
  [[ "$2" =~ ^[0-9]{4}$ ]] || { usage >&2; exit 2; }
  month=$(make_month "$3") || { usage >&2; exit 2; }
  run_one "$2" "$month" "$4"
  exit $?
fi

if [[ "$1" == all || "$1" == monthly || "$1" == annual ]]; then
  [[ $# -eq 1 ]] || { usage >&2; exit 2; }
  mode=$1
  year_filter=""
  month_filter=""
elif [[ "$1" =~ ^[0-9]{4}$ ]]; then
  (( $# <= 2 )) || { usage >&2; exit 2; }
  mode=all
  year_filter=$1
  month_filter=""
  if (( $# == 2 )); then
    month_filter=$(make_month "$2") || { usage >&2; exit 2; }
  fi
else
  usage >&2
  exit 2
fi

# Do not use a hard-coded list of months: discover whatever Load has produced.
declare -A per_period_input=()
while IFS= read -r -d '' path; do
  name=${path##*/}
  year=''; month=''
  if [[ "$name" =~ ^PHe_skim_([0-9]{4})_(0?[1-9]|1[0-2])_ ]]; then
    year=${BASH_REMATCH[1]}
    month=$((10#${BASH_REMATCH[2]}))
  elif [[ "$name" =~ ^PHe_skim_([0-9]{4})_ ]]; then
    year=${BASH_REMATCH[1]}
    month=0
  else
    continue
  fi
  [[ -z "$year_filter" || "$year" == "$year_filter" ]] || continue
  [[ -z "$month_filter" || "$month" == "$month_filter" ]] || continue
  [[ "$mode" != monthly || "$month" != 0 ]] || continue
  [[ "$mode" != annual || "$month" == 0 ]] || continue
  printf -v key '%s_%02d' "$year" "$month"
  if [[ -n ${per_period_input[$key]+x} ]]; then
    echo "ERROR: two Load files for period $key (would overwrite results):" >&2
    echo "  ${per_period_input[$key]}" >&2
    echo "  $path" >&2
    echo "Choose one explicitly via: ./run_unfolding.sh file YEAR MONTH PATH" >&2
    exit 2
  fi
  per_period_input[$key]=$path
done < <(find "$DATA_DIR" -maxdepth 1 -type f -name 'PHe_skim_*.root' -print0 | sort -z)

if (( ${#per_period_input[@]} == 0 )); then
  echo "No matching Load outputs in $DATA_DIR (see naming format in --help)." >&2
  exit 2
fi

ok=0
failed=0
while IFS= read -r key; do
  year=${key%_*}
  month=$((10#${key#*_}))
  if run_one "$year" "$month" "${per_period_input[$key]}"; then
    ((ok+=1))
  else
    ((failed+=1))
  fi
done < <(printf '%s\n' "${!per_period_input[@]}" | sort)
echo "Finished: ${ok} successful/skipped; ${failed} failed. Logs: $LOG_DIR"
(( failed == 0 ))

