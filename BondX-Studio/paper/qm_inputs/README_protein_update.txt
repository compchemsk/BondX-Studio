BondX protein update (23 November 2026)
=======================================

What changed
------------
BondX now builds ligand hydrogens and the ligand charge from one structure.
Three ligands were wrong before. They are correct now.
  1P5E  TBS  had 3 N-H hydrogens. It now has 1 N-H (C6HBr4N3, neutral).
  2W0S  BVP  had charge -2 with 2 phosphate O-H. Its phosphate now has no H (charge -2).
  4GQL  R47  kept 30 file hydrogens. Two polar H were missing. It now has 32 H (charge -3).
This changes 44 of the 89 protein models (all variants of 1P5E, 2W0S and 4GQL).
The other 45 models (1US0, 4AGL, 4DST) are unchanged.
All 27 ORCA counterpoise inputs of the protein clusters also have new SCF settings
(Sthresh 1e-6, DefGrid3, SlowConv). Two old protein runs ended in wrong SCF states.
The HAL59 and scan jobs (sets A-E) are not affected.

Jobs to run again
-----------------
  Gaussian  44 NBO + 5 counterpoise
  Q-Chem    44 NBO + 5 ALMO-EDA
  ORCA      27 counterpoise (all protein clusters)

Steps (in the folder where your Gaussian and Q-Chem jobs run)
--------------------------------------------------------------
  cd /scratch/19cy40038/IYD/BondX
  unzip -o BondX_protein_update.zip      # replaces the protein inputs and the run scripts
  cd qm_package
  bash update_protein.sh                 # moves outputs of old inputs to F_protein/_stale/
  sbatch run_orca.sh F_protein           # the 27 ORCA protein jobs (about 10 h)

Gaussian and Q-Chem
-------------------
If your Gaussian and Q-Chem jobs are still running, do nothing more.
They run the protein models last. They read each input when they reach it.
They will use the new inputs.
If they have already finished, submit the protein part again:
  sbatch run_gaussian.sh F_protein
  sbatch run_qchem.sh F_protein
When all jobs have finished, run update_protein.sh once more.
It catches a job that was running an old input during the update.
Then submit the protein part again if it moved any output.

Check
-----
  bash status.sh

New folder with several accounts (recommended)
----------------------------------------------
  cd /scratch/19cy40038/IYD/BondX/new-bugfix
  unzip -o BondX_protein_update.zip
  cd qm_package
  bash submit_protein.sh
This puts ORCA on sabyashachicy_spl, Gaussian on sabyashachicy_imf and two Q-Chem copies
on sabyashachicy_doe and sabyashachicy_fdf. A third Q-Chem copy starts on the spl node
when ORCA has finished. A refused account is replaced by the next one in the list.
