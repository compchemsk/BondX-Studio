#!/bin/bash
# BondX protein-model update (23 Nov 2026).
# Ligand protonation was corrected for 1P5E (TBS), 2W0S (BVP) and 4GQL (R47); all ORCA counterpoise
# jobs now remove near-linearly-dependent functions (Sthresh 1e-6). This script finds outputs that were
# made from an older input: the output does not contain the model ID of the current input (or, for ORCA
# counterpoise jobs, the Sthresh line). Such outputs are moved to F_protein/_stale/ (nothing is deleted).
# The run scripts then compute these jobs again with the new inputs.
# Usage:  cd qm_package ; bash update_protein.sh
cd "$(dirname "$0")" || exit 1
n=0; busy=0
for inp in F_protein/gaussian/*.gjf F_protein/gaussian_cp/*.gjf F_protein/qchem/*.in F_protein/qchem_eda/*.in F_protein/orca/*.inp F_protein/orca_cp/*.inp; do
    [ -e "$inp" ] || continue
    case "$inp" in *.gjf) out=${inp%.gjf}.log ;; *) out=${inp%.*}.out ;; esac
    [ -e "$out" ] || continue
    if [ -d "${inp%.*}.lock" ]; then echo "running now, check again later: $inp"; busy=$((busy+1)); continue; fi
    id=$(grep -o "BX-[0-9a-f]\{10\}" "$inp" | head -1)
    stale=0
    grep -q "$id" "$out" || stale=1
    case "$inp" in */orca_cp/*) grep -q "Sthresh" "$out" || stale=1 ;; esac
    if [ $stale = 1 ]; then
        d=F_protein/_stale/$(basename "$(dirname "$inp")"); mkdir -p "$d"; mv "$out" "$d/"
        n=$((n+1)); echo "stale: $out"
    fi
done
echo "$n stale output(s) moved to F_protein/_stale/; $busy job(s) still running."
[ $busy -gt 0 ] && echo "Run this script again when those jobs have finished."
exit 0
