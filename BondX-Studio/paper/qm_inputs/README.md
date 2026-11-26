# BondX paper V2 — QM jobs to run on the cluster

> **Update 23 November 2026:** the protein inputs (set F) were corrected. Read `README_protein_update.txt` first.


BondX Studio made every input in this folder. The `Scripts/` folder of the paper package holds the generating scripts.
Each model has three NBO inputs: Gaussian 16, ORCA and Q-Chem. Each model also has a manifest. The manifest lists every QM atom in input order.
ALMO-EDA exists only in Q-Chem. Gaussian and ORCA have no ALMO-EDA. Each Q-Chem EDA job of sets A, B, C and F therefore has two partner jobs:
- a Gaussian counterpoise job (`counterpoise=2`, one input);
- an ORCA counterpoise set (three single points: dimer AB, fragment A with B as ghost atoms, fragment B with A as ghost atoms).

These jobs give the counterpoise-corrected interaction energy in Gaussian and ORCA. The paper compares it with the ALMO-EDA total.

## What is inside

| Folder | Content | G16 NBO | ORCA NBO | Q-Chem NBO | Q-Chem EDA | G16 CP | ORCA CP |
|---|---|---:|---:|---:|---:|---:|---:|
| `A_HAL59` | 59 benchmark dimers (XB51 + X40, CCSD(T)/CBS references) | 59 | 59 | 59 | 59 | 59 | 177 |
| `F_protein` | 6 protein–ligand systems × 14 model variants + 3 crystal copies + 2 His110 states | 89 | 89 | 89 | 9 | 9 | 27 |
| `B_scan_distance` | 12 dimers × 8 X···Y distances (0.90–1.50 × d_eq) | 96 | 96 | 96 | 96 | 96 | 288 |
| `C_scan_angle` | 6 dimers × 7 C–X···Y angles (180°–120°) | 42 | 42 | 42 | 42 | 42 | 126 |
| `D_method_basis` | 12 dimers × 5 functionals (def2-TZVP) + 3 extra basis sets (ωB97X-D) | 84 | 96 | 96 | – | – | – |
| `E_eda_settings` | 12 dimers: EDA with def2-SVPD, def2-QZVP, and def2-TZVPD without BSSE | – | – | – | 36 | – | – |
| **Total** | **2200 input files** | 370 | 382 | 382 | 242 | 206 | 618 |

Sub-folders of each set:
- `gaussian/` (`.gjf`) and `gaussian_cp/` (`_cp.gjf`);
- `orca/` (`.inp`, with `.pc` point-charge files for QM/MM) and `orca_cp/` (`_cp_AB.inp`, `_cp_A.inp`, `_cp_B.inp`);
- `qchem/` (`.in`) and `qchem_eda/` (`_eda.in`);
- `manifests/` (`.json`), `xyz/`, and `mm_charges/` (QM/MM point charges as CSV).

## Levels of theory

| Job | Gaussian 16 | ORCA | Q-Chem |
|---|---|---|---|
| NBO (main) | `wB97XD/def2TZVP pop=(nbo)` | `wB97X-D3 def2-TZVP RIJCOSX def2/J NBO` | `METHOD wB97X-D, BASIS def2-TZVP, NBO 1` |
| EDA | – | – | `wB97X-V/def2-TZVPD, EDA2 1, EDA_BSSE TRUE` |
| Counterpoise interaction energy | `wB97XD/def2TZVP counterpoise=2` | `wB97X-V def2-TZVPD RIJCOSX` (AB, A, B) | (EDA total) |

- I uses the def2 ECP in all three programs. Gaussian and ORCA add it by themselves. The Q-Chem inputs contain `ECP def2-ECP`.
- ORCA has no ωB97X-D. Its inputs use ωB97X-D3. The paper reports the ORCA numbers as a separate column.
- Protein cluster models use a dielectric of ε = 4: CPCM in Gaussian and ORCA, PCM in Q-Chem. The protein EDA inputs are gas phase.
- QM/MM models carry point charges: `charge` in Gaussian, `%pointcharges` in ORCA, `$external_charges` in Q-Chem.
- The method set D uses D3(BJ) for B3LYP and PBE0. ωB97X-D, M06-2X and HF carry no extra dispersion.
- Gaussian 16 has no built-in def2-TZVPD. The def2-TZVPD jobs in set D therefore exist for ORCA and Q-Chem only.

## How to run

There are three job scripts. Each script runs every input of one program, one after another, on one 40-core node.

| Script | Program | Inputs | Folders |
|---|---|---|---|
| `run_qchem.sh` | Q-Chem 6.0 | 624 | `qchem`, `qchem_eda` |
| `run_gaussian.sh` | Gaussian 16 | 576 | `gaussian`, `gaussian_cp` |
| `run_orca.sh` | ORCA 5.0.1 | 1000 | `orca`, `orca_cp` |

1. Copy the whole `qm_package` folder to your scratch space, for example `/scratch/19cy40038/`.
2. Open `run_orca.sh`. Put the path of the NBO 7 program in `NBOEXE`. ORCA needs it for the NBO step. The script also searches a few usual places. Without NBO 7, the 382 ORCA NBO inputs are skipped. The 618 ORCA counterpoise inputs still run.
3. Submit the three jobs from inside the folder:

   ```bash
   cd /scratch/19cy40038/qm_package
   sbatch run_qchem.sh
   sbatch run_gaussian.sh
   sbatch run_orca.sh
   ```

4. Check the progress at any time:

   ```bash
   bash status.sh
   ```

Notes:

- The scripts use your settings: partition `medium`, account `sabyashachicy_spl`, 40 cores, 60000 MB and 30 days.
- Q-Chem uses `module load apps/qchem/6.0`. Gaussian uses `g16root=/home/sabyashachicy`. ORCA uses `module load apps/orca/5.0.1`.
- Small dimers run first. The protein models run last.
- A finished input is skipped. A job that ran out of time can be submitted again. Failed inputs are tried again in the new job.
- You can submit the same script twice. The two jobs then share the inputs and finish sooner. A lock folder (`*.lock`) marks the input in use.
- You can run a part only. Give a word from the path: `sbatch run_qchem.sh F_protein` or `sbatch run_orca.sh orca_cp`.
- Each job writes a short report in `logs/`. Outputs appear next to the inputs. Q-Chem and ORCA write `*.out`. Gaussian writes `*.log`.
- The scripts delete Gaussian `.chk` files and ORCA `.gbw` and `.tmp` files after a finished job. This saves disk space.

## What to send back

Zip the whole `qm_package` folder with the outputs inside it. Keep the folder structure:

```bash
zip -r qm_package_outputs.zip qm_package -x "*.chk" "*.gbw" "*.tmp*"
```

Check for errors before you send it:

```bash
grep -L "Thank you very much for using Q-Chem" */qchem*/*.out
grep -L "Normal termination" */gaussian*/*.log
grep -L "ORCA TERMINATED NORMALLY" */orca*/*.out
```

These commands list the unfinished files.

The results are collected with:

```bash
python3 Scripts/07_collect_results.py qm_package results
python3 Scripts/08_make_figures.py .
```

## Notes

- Do not change atom order. Do not let the programs reorient the molecule: `nosymm` and `SYM_IGNORE` are set. The NBO gate checks the output geometry against the manifest.
- The NBO print threshold is the program default. Each output prints it. BondX records it. With a threshold of 0.5 kcal/mol (the gate value) the gate decision is not affected. On the long-distance scan points, E(2) values below the print threshold are not printed. The tables report them as below the threshold.
- Memory: 50 GB and 40 threads per job (Gaussian `%mem=50GB`, ORCA `%maxcore 937` × 40, Q-Chem `MEM_TOTAL 50000`). This matches a 60 GB node.
