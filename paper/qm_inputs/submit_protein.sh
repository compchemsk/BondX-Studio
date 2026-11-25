#!/bin/bash
# Submits the BondX protein jobs over several Slurm accounts (one node each).
#   ORCA      -> sabyashachicy_spl  (your third spl node)
#   Gaussian  -> sabyashachicy_imf
#   Q-Chem    -> sabyashachicy_doe
#   Q-Chem    -> sabyashachicy_fdf  (second copy; the two copies share the inputs)
#   Q-Chem    -> sabyashachicy_spl  (third copy; starts on the spl node after ORCA ends)
# If an account is refused, the next one in its list is tried.
# Usage:  cd qm_package ; bash submit_protein.sh
# Change the lists below if you want another arrangement.
cd "$(dirname "$0")" || exit 1

ORCA_ACC="sabyashachicy_spl sabyashachicy_imf"
G16_ACC="sabyashachicy_imf sabyashachicy_fdf sabyashachicy_doe"
QC1_ACC="sabyashachicy_doe sabyashachicy_fdf sabyashachicy_imf"
QC2_ACC="sabyashachicy_fdf sabyashachicy_doe"

submit () {   # $1 = account list, $2 = script, $3.. = extra sbatch options
    local accs=$1 script=$2; shift 2
    for a in $accs; do
        out=$(sbatch --parsable -A "$a" "$@" "$script" 2>&1)
        if [[ $out =~ ^[0-9]+ ]]; then echo "${out%%;*}"; echo "  $script -> account $a, job ${out%%;*}" >&2; return 0; fi
        echo "  $script: account $a refused ($out)" >&2
    done
    echo "  $script: not submitted" >&2; return 1
}

echo "Submitting BondX protein jobs from $PWD"
J_ORCA=$(submit "$ORCA_ACC" run_orca.sh)
submit "$G16_ACC" run_gaussian.sh > /dev/null
submit "$QC1_ACC" run_qchem.sh > /dev/null
submit "$QC2_ACC" run_qchem.sh > /dev/null
if [ -n "$J_ORCA" ]; then submit "sabyashachicy_spl" run_qchem.sh --dependency=afterany:$J_ORCA > /dev/null; fi
echo "Check with:  squeue -u $USER   and   bash status.sh"
