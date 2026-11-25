# Paper scripts and data

This folder reproduces every number, figure and table of the paper
*BondX Studio: From Geometric Contacts to Verified Halogen Bonds with Force-Field-Consistent QM Models, NBO and ALMO-EDA*.

## Requirements

- Python 3.10 or newer with numpy, pandas, scipy, matplotlib and Pillow
- Node.js 18 or newer (the scripts call the BondX engine through `Scripts/bondx_cli.js`)
- Playwright with Chromium (only for the screenshots of script 10)
- The GMTKN55 HAL59 folder (only for script 01): <https://github.com/grimme-lab/GMTKN55>

## Folders

| Folder | Content |
|---|---|
| `Scripts/` | All scripts. `engine.js` and `ff_data.js` are identical to the web application. |
| `data/hal59/` | The 59 HAL59 dimers (XYZ) and `hal59.csv` with reference energies and X/A/Y atoms |
| `data/scans/` | Rigid distance and angle scans of 12 dimers |
| `data/protein/` | The 56 protein files of the test set |
| `qm_inputs/` | All 2200 QM input files with their manifests and job lists |
| `results/` | Parsed QM results and analysis tables (CSV and JSON) |
| `tables/` | Generated LaTeX tables of the ESI |
| `checks/` | Independent re-implementations of the QM parsing and of the hIYD frame analysis |

The QM output files (about 400 MB) are not in this repository.
The parsed results in `results/` are enough for scripts 08, 11, 12, 16 and 17.

## Order of the scripts

| Step | Script | Task |
|---|---|---|
| 0 | `00_pubmed_xb_counts.py` | PubMed counts (Figure 1); needs internet |
| 1 | `01_prepare_hal59.py` | HAL59 geometries and reference energies |
| 2 | `02_make_scans.py` | Rigid distance and angle scans |
| 3 | `03_t1_benchmark.py` | Tier T1 on HAL59 and the scans |
| 4 | `06_t1_sensitivity.py` | T1 parameter grid |
| 5 | `14_ligand_check.py` | Ligand hydrogens against the PDB FORMUL records |
| 6 | `04_make_qm_inputs_small.py` | QM inputs of sets A–E |
| 7 | `05_protein_cases.py` | Protein T1, 89 QM models, set F inputs, MM field |
| 8 | (cluster) | Run the QM jobs (`qm_inputs/run_*.sh`) |
| 9 | `07_collect_results.py` | Read every QM output with the BondX parsers |
| 10 | `13_orca_cp_analysis.py` | ORCA counterpoise statistics |
| 11 | `15_hiyd_analysis.py` | hIYD case study (needs the hIYD trajectories) |
| 12 | `17_qm_nbo_eda_analysis.py` | NBO gate, ALMO-EDA and protein analysis |
| 13 | `16_hiyd_figures_tables.py` | hIYD figures and tables |
| 14 | `09_performance.py` | Timings |
| 15 | `11_make_tables.py` | ESI tables |
| 16 | `08_make_figures.py`, `12_parameter_figure.py` | Figures |
| 17 | `10_screenshots.py` | Interface screenshots |
| 18 | `inline_tables.py` | Writes the tables into `ESI_V2.tex` |

`run_pipeline.sh` runs steps 1–15 in a fresh folder. We ran it three times. All result files agree. The timing columns are the only exception.

## Notes

- The hIYD trajectories come from our earlier study (Karmakar and Mishra, 2026). They are not part of this repository. `results/hiyd/` holds the parsed hIYD results.
- `bxstyle.py` sets the common figure style (large labels, legends without frames).
