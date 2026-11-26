#!/bin/bash
#SBATCH -J bx_orca                  # job name
#SBATCH -o %j-out-orca              # stdout
#SBATCH -e %j-err-orca              # stderr
#SBATCH -p medium                   # partition
#SBATCH -A sabyashachicy_spl        # account
#SBATCH --nodes=1
#SBATCH -n 40                       # 40 MPI ranks (inputs use %pal nprocs 40)
#SBATCH --mem=60000                 # MB per node (inputs use %maxcore 937)
#SBATCH -t 30-00:00:00
# ------------------------------------------------------------------
# BondX - all ORCA jobs (NBO + counterpoise): 1000 inputs
#   cd qm_package
#   sbatch run_orca.sh                   -> runs everything
#   sbatch run_orca.sh orca_cp           -> only the counterpoise jobs
#   sbatch run_orca.sh F_protein         -> only inputs whose path contains F_protein
# Finished jobs are skipped, so you can resubmit after a time-out.
# You can submit the same script twice; the two jobs share the list.
#
# NBO: ORCA calls an external NBO 7 program. Put its path in NBOEXE below.
# Without NBO 7, the 382 NBO inputs are skipped and the 618 counterpoise
# inputs still run.
# ------------------------------------------------------------------
PROG=orca
SUBDIRS="orca orca_cp"
EXT=inp
OUTEXT=out

NBOEXE=""          # e.g. /home/apps/nbo7/bin/nbo7.i8.exe  (NBO 6 also works)

module load apps/orca/5.0.1
export ORCAPATH=/home/apps/orca/orca_5_0_1_linux_x86-64_shared_openmpi411/
export PATH=$ORCAPATH:$PATH
export LD_LIBRARY_PATH=$ORCAPATH:$LD_LIBRARY_PATH
exe=$ORCAPATH/orca                  # ORCA needs its full path to run in parallel
export OMP_NUM_THREADS=1

if [ -z "$NBOEXE" ]; then             # try to find NBO 7 (or NBO 6) by itself
    module load nbo 2>/dev/null || module load apps/nbo 2>/dev/null || module load apps/nbo7 2>/dev/null
    for n in nbo7.i8.exe nbo7.i4.exe nbo6.i8.exe nbo6.i4.exe; do
        NBOEXE=$(command -v $n 2>/dev/null) && break
    done
fi
if [ -z "$NBOEXE" ]; then
    NBOEXE=$(find /home/apps /home/sabyashachicy /opt /usr/local -maxdepth 6 \( -name 'nbo7.i[48].exe' -o -name 'nbo6.i[48].exe' \) 2>/dev/null | head -1)
fi
if [ -n "$NBOEXE" ] && [ -x "$NBOEXE" ]; then
    export NBOEXE; HAVE_NBO=1; echo "NBO program: $NBOEXE"
else
    HAVE_NBO=0; echo "NBO 7 not found: ORCA NBO inputs are skipped, counterpoise inputs run."
fi

finished ()    { grep -q "ORCA TERMINATED NORMALLY" "$1" 2>/dev/null; }
skip_reason () { [ $HAVE_NBO = 0 ] && grep -q "^!.*NBO" "$1"; }
run_one ()     { "$exe" "$1" > "$2.out" 2>&1; }
cleanup ()     { ( cd "$1" && rm -f "$2".gbw "$2".tmp "$2".*.tmp "$2".densities "$2".densitiesinfo "$2".bibtex "$2"_property.txt "$2".ges "$2".cpcm "$2".cpcm_corr ); }
final_cleanup (){ :; }

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
