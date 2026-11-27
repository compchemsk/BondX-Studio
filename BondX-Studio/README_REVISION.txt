BondX Studio — Revised pre-submission manuscript package

MAIN FILES
- Manuscript_V3.tex
- ESI_V3.tex
- MS_Figures_V3/ : available manuscript figures, including the three separately cropped Figure 2 panels
- ESI_Figures_V3/ : available ESI figures
- app/ : updated application source, regression tests, and generated app/index.html
- index.html : updated repository-root deployment file
- paper/ : scripts, input datasets, regenerated results, and tables
- Checks/validation_summary.txt : final automated-test and syntax-check summary

T1 IMPLEMENTATION
- The default donor set is C-bound Cl, Br, and I; non-C-bound HAL59 dimers remain exploratory QM cases and are excluded from T1 performance statistics.
- Atom-acceptor theta2 scoring now uses hybridization-dependent optima: sp 180 degrees, sp2 120 degrees, and sp3 109.5 degrees, with a heuristic hybridization assignment and an sp2 fallback when unresolved.
- New regression tests cover the C-bound donor scope and linear nitrile scoring.
- The existing distance cut-offs were retained and are now described transparently as empirical operational limits. A radius-only recalibration has NOT been performed.

VALIDATION
- 123 automated engine tests pass.
- Independent force-field reference comparisons pass.
- Updated HAL59 T1 benchmark, figure-generation and table-generation scripts were run with the revised engine.
- The complete end-to-end paper pipeline was NOT rerun after the T1 implementation changes.

MISSING INPUT ASSETS — NOT FABRICATED
The uploaded materials did not include these files, which are still referenced by the LaTeX sources:
- MS_Figures_V3/Fig1_workflow.pdf
- MS_Figures_V3/Fig3_T1_model.pdf
- MS_Figures_V3/TOC_graphic.pdf
- ESI_Figures_V3/FigS11_screens.pdf
- References_V3.bib

Consequently, the LaTeX sources pass syntax compilation using temporary, local placeholders for missing figures, but a final publication PDF with correct figures and references cannot be produced from the provided assets alone. The temporary placeholders are NOT included in this package. Add the original figure files and References_V3.bib before producing the final PDFs.

The revised working tree is based on public repository commit:
971f1fa3bd1adcf74b2b6b078bd7ec408d4f1bc1
The revised files in this package have not yet been committed or pushed to the public repository. No Zenodo DOI has been registered.
