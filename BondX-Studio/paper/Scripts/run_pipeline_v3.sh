#!/usr/bin/env bash
# run_pipeline_v3.sh -- re-runs the analysis of BondX Studio 3.0 (scripts 18-21, 27, 28, 11, and 26) in a fresh folder.
# Usage: bash run_pipeline_v3.sh <new_run_folder> <master_paper_folder> <hIYD_for_BondX folder> <app folder>
# Inputs: the QM outputs of the master qm_package, the protein files, and the V2 result files of the master
# (written by run_pipeline.sh). No QM job is run.
set -euo pipefail
RUN=$1; MASTER=$2; HIYD=$3; APP=$4
if [ -e "$RUN" ]; then echo "folder $RUN exists; choose a new folder"; exit 1; fi
mkdir -p "$RUN/results/v3" "$RUN/tables"
cp -r "$MASTER/Scripts" "$RUN/Scripts"; rm -rf "$RUN/Scripts/__pycache__"
ln -s "$(readlink -f "$MASTER/qm_package")" "$RUN/qm_package"
ln -s "$(readlink -f "$MASTER/data")" "$RUN/data"
ln -s "$(readlink -f "$MASTER/qm_extra")" "$RUN/qm_extra"
mkdir -p "$RUN/results/extra"; cp "$MASTER/results/extra/jobs_extra.csv" "$MASTER/results/extra/results_extra.csv" "$RUN/results/extra/"
find "$MASTER/results" -maxdepth 1 -type f -exec cp {} "$RUN/results/" \;
cp -r "$MASTER/results/hiyd" "$RUN/results/hiyd"
cd "$RUN/Scripts"
node 18_v3_nbo_reparse.js ../qm_package ../results/v3
node 19_v3_field.js ../data/protein ../results/v3
node 20_v3_hiyd_t2.js "$HIYD" ../results/v3
python3 21_v3_analysis.py ..
node 27_v3_extra_nbo.js ../qm_extra ../results/v3
python3 27_extra_parse.py ../qm_extra ../results/extra
python3 28_v3_round2_analysis.py ..
python3 11_make_tables.py .. "$APP"
python3 26_v3_tables.py .. "$APP"
echo "run finished: $RUN"
