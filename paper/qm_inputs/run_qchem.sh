#!/bin/bash
#SBATCH -J bx_qchem                 # job name
#SBATCH -o %j-out-qchem             # stdout
#SBATCH -e %j-err-qchem             # stderr
#SBATCH -p medium                   # partition
#SBATCH -A sabyashachicy_spl        # account
#SBATCH --nodes=1                   # Q-Chem OpenMP runs on one node
#SBATCH --ntasks=1                  # ONE process ...
#SBATCH --cpus-per-task=40          # ... with 40 threads
#SBATCH --mem=60000                 # MB per node
#SBATCH -t 30-00:00:00
# ------------------------------------------------------------------
# BondX - all Q-Chem jobs (NBO + ALMO-EDA): 624 inputs
#   cd qm_package
#   sbatch run_qchem.sh                  -> runs everything
#   sbatch run_qchem.sh F_protein        -> only inputs whose path contains F_protein
#   sbatch run_qchem.sh qchem_eda        -> only the EDA jobs
# Finished jobs are skipped, so you can resubmit after a time-out.
# You can submit the same script twice; the two jobs share the list.
# ------------------------------------------------------------------
PROG=qchem
SUBDIRS="qchem qchem_eda"
EXT=in
OUTEXT=out

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-40}
export MKL_NUM_THREADS=$OMP_NUM_THREADS
module purge
module load apps/qchem/6.0
export QCSCRATCH=/scratch/sabyashachicy/qchem/${SLURM_JOB_ID:-$$}   # same scratch as your job.sh
export QCLOCALSCR=$QCSCRATCH
mkdir -p "$QCSCRATCH"

finished ()    { grep -q "Thank you very much for using Q-Chem" "$1" 2>/dev/null; }
skip_reason () { return 1; }
run_one ()     { qchem -nt "$NT" "$1" "$2.out" > /dev/null 2>&1; }
cleanup ()     { rm -rf "$QCSCRATCH"/* 2>/dev/null; }
final_cleanup (){ rm -rf "$QCSCRATCH"; }

# ------------------------------------------------------------------
# Common part (same in all three scripts)
# ------------------------------------------------------------------
PKG=${SLURM_SUBMIT_DIR:-$PWD}
cd "$PKG" || exit 1
found=0; for d in A_HAL59 B_scan_distance C_scan_angle D_method_basis E_eda_settings F_protein; do [ -d "$d" ] && found=1; done
if [ $found = 0 ]; then
    echo "No input folders here. Submit this script from inside the qm_package folder."; exit 1
fi
NT=${SLURM_CPUS_PER_TASK:-${SLURM_NTASKS:-40}}
JOBID=${SLURM_JOB_ID:-$$}
LOGDIR=$PKG/logs; mkdir -p "$LOGDIR"
SUMMARY=$LOGDIR/${PROG}_$JOBID.txt

# Build the input list: small sets first, protein models last.
make_list () {
    local set sub f
    for set in A_HAL59 E_eda_settings B_scan_distance C_scan_angle D_method_basis F_protein; do
        for sub in $SUBDIRS; do
            [ -d "$set/$sub" ] || continue
            for f in "$set/$sub"/*.$EXT; do [ -e "$f" ] && echo "$f"; done
        done
    done
}

# A lock folder (<input>.lock) lets two copies of this script share the work.
# A lock left behind by a job that is no longer in the queue is taken over.
owner_gone () {      # true only when we are sure the owning job has ended
    local o=$1 msg rc
    if [ -n "$SLURM_JOB_ID" ] && command -v squeue >/dev/null; then
        msg=$(squeue -h -j "$o" -o %T 2>&1); rc=$?
        [ $rc = 0 ] && [ -z "$msg" ] && return 0
        echo "$msg" | grep -qi "invalid job id" && return 0
        return 1
    fi
    ! kill -0 "$o" 2>/dev/null   # outside Slurm: owner is a process id
}
claim () {
    local lk=$1.lock owner
    if mkdir "$lk" 2>/dev/null; then echo "$JOBID" > "$lk/owner"; return 0; fi
    owner=$(cat "$lk/owner" 2>/dev/null)
    if [ -z "$owner" ]; then
        [ -n "$(find "$lk" -maxdepth 0 -mmin +10 2>/dev/null)" ] || return 1
    elif ! owner_gone "$owner"; then
        return 1
    fi
    rm -rf "$lk"; mkdir "$lk" 2>/dev/null || return 1
    echo "$JOBID" > "$lk/owner"; return 0
}
release () { rm -rf "$1.lock"; }

# Optional filters: sbatch run_xxx.sh F_protein   or   sbatch run_xxx.sh _cp
FILTER=$(IFS='|'; echo "$*")
LIST=$LOGDIR/${PROG}_$JOBID.list
if [ -n "$FILTER" ]; then make_list | grep -E "$FILTER" > "$LIST"; else make_list > "$LIST"; fi

N=$(wc -l < "$LIST"); nok=0; nfail=0; nskip=0; i=0
echo "Job $JOBID on $(hostname) | $PROG | $NT cores | $N inputs | start $(date)" | tee "$SUMMARY"

while read -r inp <&3; do   # list on channel 3: programs cannot eat it
    i=$((i+1))
    dir=$(dirname "$inp"); base=$(basename "$inp"); stem=${base%.*}
    out=$dir/$stem.$OUTEXT
    if finished "$out"; then nskip=$((nskip+1)); continue; fi
    if skip_reason "$inp"; then nskip=$((nskip+1)); continue; fi
    claim "$PKG/$dir/$stem" || { echo "[$i/$N] busy  $inp (another job has it)"; continue; }
    if finished "$out"; then release "$PKG/$dir/$stem"; nskip=$((nskip+1)); continue; fi
    t0=$(date +%s)
    ( cd "$PKG/$dir" && run_one "$base" "$stem" ) < /dev/null
    dt=$(( $(date +%s) - t0 ))
    if finished "$out"; then
        nok=$((nok+1)); echo "[$i/$N] OK    $inp  (${dt} s)" | tee -a "$SUMMARY"
        cleanup "$PKG/$dir" "$stem"
    else
        nfail=$((nfail+1)); echo "[$i/$N] FAIL  $inp  (${dt} s) -> see $out" | tee -a "$SUMMARY"
    fi
    release "$PKG/$dir/$stem"
done 3< "$LIST"

echo "Done: $nok finished now, $nskip skipped, $nfail failed | end $(date)" | tee -a "$SUMMARY"
final_cleanup
