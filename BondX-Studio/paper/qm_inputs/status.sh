#!/bin/bash
# Shows how many BondX jobs have finished.   Usage:  cd qm_package ; bash status.sh
cd "$(dirname "$0")"
count () {   # $1 label, $2 subfolders, $3 input ext, $4 output ext, $5 success text
    local tot=0 ok=0 bad=0 sub f o
    for sub in $2; do
        for f in */$sub/*.$3; do
            [ -e "$f" ] || continue; tot=$((tot+1)); o=${f%.*}.$4
            if grep -q "$5" "$o" 2>/dev/null; then ok=$((ok+1)); elif [ -s "$o" ] && [ ! -d "${f%.*}.lock" ]; then bad=$((bad+1)); fi
        done
    done
    printf "%-9s %5d inputs  %5d finished  %5d failed  %5d to do\n" "$1" $tot $ok $bad $((tot-ok-bad))
}
count Q-Chem   "qchem qchem_eda"      in  out "Thank you very much for using Q-Chem"
count Gaussian "gaussian gaussian_cp" gjf log "Normal termination of Gaussian"
count ORCA     "orca orca_cp"         inp out "ORCA TERMINATED NORMALLY"
echo "Failed jobs are retried when you submit the script again."
