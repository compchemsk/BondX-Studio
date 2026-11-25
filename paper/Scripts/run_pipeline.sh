#!/usr/bin/env bash
# run_pipeline.sh -- full re-run of the paper analysis in a fresh folder.
# Usage: bash run_pipeline.sh <new_run_folder> <master_paper_folder> <GMTKN55 HAL59 folder> <hIYD_for_BondX folder>
# The QM outputs (*.log, *.out) are taken from the master qm_package. Every input file is written again.
set -euo pipefail
RUN=$1; MASTER=$2; GMTKN=$3; HIYD=$4
if [ -e "$RUN" ]; then echo "folder $RUN exists; choose a new folder"; exit 1; fi
mkdir -p "$RUN/data" "$RUN/results" "$RUN/tables" "$RUN/Figures_V2_MS" "$RUN/Figures_V2_ESI"
cp -r "$MASTER/Scripts" "$RUN/Scripts"
cp -r "$MASTER/data/protein" "$RUN/data/protein"
cd "$RUN/Scripts"
python3 01_prepare_hal59.py "$GMTKN" ../data/hal59
python3 02_make_scans.py ../data/hal59 ../data/scans
python3 03_t1_benchmark.py ..
python3 06_t1_sensitivity.py ..
python3 14_ligand_check.py ..
python3 04_make_qm_inputs_small.py ../data ../qm_package
python3 05_protein_cases.py ..
# link the QM outputs of the master package into the new package
( cd "$MASTER/qm_package" && find . -type f \( -name '*.log' -o -name '*.out' \) ) | while read -r f; do
  mkdir -p "../qm_package/$(dirname "$f")"; ln "$MASTER/qm_package/$f" "../qm_package/$f"; done
python3 07_collect_results.py ../qm_package ../results
python3 13_orca_cp_analysis.py ..
python3 15_hiyd_analysis.py .. "$HIYD"
python3 17_qm_nbo_eda_analysis.py ..
python3 16_hiyd_figures_tables.py ..
cp "$MASTER/results/performance.csv" ../results/performance.csv   # timings come from script 09 (machine dependent)
python3 11_make_tables.py ..
echo "run finished: $RUN"
