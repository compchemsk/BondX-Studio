# Paper scripts and data

This folder reproduces every number, figure, and table of the paper
*BondX Studio: Screening Geometric Halogen Contacts with Force-Field-Consistent QM Models, NBO Analysis, and ALMO-EDA*.

## Requirements

- Python 3.10 or newer with numpy, pandas, scipy, matplotlib, and Pillow
- Node.js 18 or newer (the scripts call the BondX engine through `Scripts/bondx_cli.js` or `require('./engine.js')`)
- Playwright with Chromium (only for the screenshots of script 23)
- XeLaTeX with TikZ (only for Figure 1, `Fig1_workflow.tex`)
- The GMTKN55 HAL59 folder (only for script 01): <https://github.com/grimme-lab/GMTKN55>
- Figures use Arial; where Arial is not installed, the metric-identical Liberation Sans is used.

## Folders

| Folder | Content |
|---|---|
| `Scripts/` | All scripts. `engine.js` and `ff_data.js` are identical to the web application (version 3.0.0). |
| `data/hal59/` | The 59 HAL59 dimers (XYZ) and `hal59.csv` with reference energies and X, A, and Y atoms |
| `data/scans/` | Rigid distance and angle scans of 12 dimers |
| `data/protein/` | The 56 protein files of the test set |
| `qm_inputs/` | All 2200 QM input files with their manifests and job lists |
| `results/` | Parsed QM results and analysis tables (CSV and JSON) |
| `results/v3/` | Version-3.0 analysis: valence NBO screen, grey zone, BXQI, field tests, hIYD T2 and T5 |
| `results/extra/` | Job sets G (poor Lewis bases and controls), J (4DST with the intact amide), and K (capped published hIYD poses) |
| `BondX_paper_scripts.ipynb` | All Python scripts in one notebook, with the order in which to run them |
| `tables/` | Generated LaTeX tables of the ESI |
| `checks/` | Independent re-implementations of the QM parsing and of the hIYD frame analysis |

The QM output files (about 400 MB) are not in this repository.
The parsed results in `results/` are enough for scripts 11, 12, 16, 21, 22, 25, and 26.

## Order of the scripts

| Step | Script | Task |
|---|---|---|
| 0 | `00_pubmed_xb_counts.py` | PubMed counts (ESI Table S1); needs internet |
| 1 | `01_prepare_hal59.py` | HAL59 geometries and reference energies |
| 2 | `02_make_scans.py` | Rigid distance and angle scans |
| 3 | `03_t1_benchmark.py` | Tier T1 on HAL59 and the scans |
| 4 | `06_t1_sensitivity.py` | T1 parameter grid |
| 5 | `14_ligand_check.py` | Ligand hydrogens against the PDB FORMUL records |
| 6 | `04_make_qm_inputs_small.py` | QM inputs of sets A–E |
| 7 | `05_protein_cases.py` | Protein T1, 89 QM models, and set F inputs |
| 8 | (cluster) | Run the QM jobs (`qm_inputs/run_*.sh`) |
| 9 | `07_collect_results.py` | Read every QM output with the BondX parsers |
| 10 | `13_orca_cp_analysis.py` | ORCA counterpoise statistics |
| 11 | `15_hiyd_analysis.py` | hIYD case study (needs the hIYD trajectories) |
| 12 | `17_qm_nbo_eda_analysis.py`, `16_hiyd_figures_tables.py` | Rydberg counts, EDA and counterpoise columns, and hIYD tables |
| 13 | `18_v3_nbo_reparse.js` | Re-read all NBO outputs with the version-3.0 screen |
| 14 | `19_v3_field.js` | MM field tests (explicit H, RCD, radius, and force field) |
| 15 | `20_v3_hiyd_t2.js` | Tier T2 on the 108 hIYD windows |
| 16 | `21_v3_analysis.py` | Statistics, grey zone, BXQI, field versus ΣE(2), and hIYD T5 |
| 16a | `27_v3_extra_nbo.js`, `27_extra_parse.py` | Sets G, J, and K with the engine and with an independent Python parser (needs the outputs of these sets) |
| 16b | `28_v3_round2_analysis.py` | Acceptor classes, grey-zone extremes, competing hydrogen bonds, and hIYD contact persistence |
| 17 | `09_performance.py` | Timings |
| 18 | `11_make_tables.py`, `26_v3_tables.py` | ESI tables |
| 19 | `12_parameter_figure.py`, `22_v3_figures.py`, `25_v3_toc.py` | Figures 3–9, S1–S10, S12, and the TOC graphic |
| 20 | `23_v3_screenshots.py`, `24_v3_screen_figures.py` | Interface screenshots (Figure 2 and Figure S11) |
| 21 | `inline_tables.py` | Writes the tables into `ESI_V3.tex` |

`run_pipeline.sh` runs steps 1–12 in a fresh folder, and `run_pipeline_v3.sh` runs steps 13–16b and 18.
Each was run three times in fresh folders; all result files agree, and the timing columns are the only exception.

## Notes

- The hIYD trajectories come from our earlier study (Karmakar and Mishra, 2026) and are not part of this repository. `results/hiyd/` holds the parsed hIYD results.
- `bxstyle3.py` sets the common figure style (Arial, large labels, and legends without frames).
- The 89 protein models were built with BondX Studio 2.0, which differs from the defaults of version 3.0 in four settings (backbone amides, hydroxyl rotamers, explicit MM hydrogens, and the embedding radius); version 3.0 also adds hydrogen-bond donors within 3.6 Å of the halogen to fragment 2. Section S6 of the Supporting Information lists the effect on each system.
