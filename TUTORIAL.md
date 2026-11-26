# BondX Studio tutorial

This is the tutorial of BondX Studio 3.0. The same text opens with the **Tutorial** button inside the program.

---

## Page 1 · What BondX Studio does

BondX Studio screens halogen bonds (XBs) in proteins, ligands, and crystals. It runs in your browser, and your structures and outputs stay on your computer. Network calls are made only to download a PDB entry from RCSB and to load the 3Dmol.js viewer; the hosted website also counts anonymous page views.

The analysis is split into tiers, and each tier adds one kind of evidence. The sidebar uses the same tier labels.

| Tier | Page | Evidence |
|---|---|---|
| T0 | Load structure | Atoms, chains, hydrogens, charges, and virtual sites |
| T1 | Find halogen contacts | Geometry: distance, two angles, and the geometric score *G* |
| T2 | Dynamics | Occupancy and residence time of each contact over a trajectory |
| T3a | QM model & charges | Cluster or QM/MM input files with force-field protonation and charges |
| T3b | NBO screen | Orbital support from the NBO donor–acceptor energy ΣE(2) |
| T4 | ALMO-EDA | Frozen, polarization, and charge-transfer energies |
| T5 | Rank & report | Ranking index BXQI and the full run record |

**Quick start.** Open **T0 Load structure**, press **Load demo**, and then press the button at the bottom right of each page. The demo takes less than one minute.

You can stop after any tier. A result from T1 alone is valid, but it is marked as geometry only.

---

## Page 2 · Load a structure and find contacts

### T0 · Load structure

- Type a four-letter PDB ID and press **Fetch**, or drop a local `.pdb`, `.ent`, or `.pqr` file.
- A PQR file brings its own partial charges (for example, from PDB2PQR with PROPKA). A PDB file uses the force-field charges of the program.
- Pick the chains to keep and check the composition table.
- Extra-point atoms (EP, LP, and XX) of σ-hole force fields are removed on load, and the message says how many were removed.

### T1 · Find halogen contacts

- The program finds every C–X···Y contact with X = Cl, Br, or I.
- It measures the X···Y distance *d*, the C–X···Y angle θ1, and the X···Y–R angle θ2, and it reports the score *G* between 0 and 1.
- Poor Lewis bases are flagged: N atoms with a delocalized lone pair (backbone amide N, Asn ND2, Gln NE2, and Trp NE1) and protonated N atoms (Lys NZ and Arg). Such a contact can also be an N–H···X hydrogen bond.
- Contacts with *d* below 0.80 of the van der Waals sum are flagged as short.
- Sort or filter the table, and click a row to see the contact in the 3D viewer.

All cut-offs are on the page **All parameters**. A changed value is marked in orange and is written to the report.

---

## Page 3 · Dynamics and the QM model

### T2 · Dynamics (optional)

- Drop a topology `.pdb` whose atom order matches the trajectory, and then a `.dcd`, multi-model `.pdb`, or multi-frame `.xyz` file.
- Convert GROMACS or AMBER files once with `gmx trjconv` or `cpptraj`.
- Set the time between saved frames (ps). The residence time τ and the scale τ0 are in ps, so the score *D* does not depend on the saving stride.
- Frames of biased runs (umbrella windows) are counted without reweighting.

### T3a · QM model and charges

- Choose a **cluster** model or a **QM/MM** model with point charges.
- Cut bonds are capped with hydrogen atoms. A backbone acceptor keeps its amide: the neighbouring residue is added.
- Check the ligand charge. Aliphatic amines, tetrazoles, and acyl sulfonamides are flagged, because their pH-7 state is usually charged. A `.mol2` or `.pqr` file with a pKa-based state is the recommended source.
- Gaussian 16 with its built-in NBO is the tested route. ORCA and Q-Chem NBO inputs are untested and are hidden unless you allow untested routes.
- Download the input files. Each file carries a manifest ID that links the output back to the contact. BondX Studio does not run QM jobs.

---

## Page 4 · NBO screen, ALMO-EDA, and the ranking

### T3b · NBO screen

- Drop the Gaussian output. The program reads the second-order perturbation table.
- It adds the valence terms LP(Y) → σ*(C–X), or the ring bonds for a π acceptor. The Rydberg terms LP(Y) → RY*(X) are shown beside the sum but are not part of the decision.
- Intramolecular terms are excluded; only the T1 acceptor counts.
- ΣE(2) ≥ 0.8 kcal/mol: orbital-supported. 0.3–0.8: inconclusive. Below 0.3: not supported.
- Donation from the halogen into N–H or O–H antibonds of the acceptor is shown as a sign of a competing hydrogen bond.
- ΣE(2) is not a binding energy. It stays high for bent and compressed contacts, so the T1 limits remain part of the decision.

### T4 · ALMO-EDA

- Drop the Q-Chem EDA output. The program reads FRZ, POL, CT, and the total energy in kcal/mol, and it reports χCT = CT / (POL + CT).
- In a protein cluster the EDA covers every contact between the two fragments, not the halogen bond alone.

### T5 · Rank and report

- BXQI is a weighted geometric mean of the tiers that could be evaluated, from 0 to 100. It ranks candidates; it has no grades.
- The ranking first follows the NBO decision and then the index.
- Press **Report** for a single HTML file with every input, parameter, and result.

---

## Page 5 · What it can and cannot do

**It can**

- Read PDB, PQR, CIF, XYZ, DCD, and multi-model PDB files.
- Find Cl, Br, and I contacts to O, N, S, protein rings, and nucleobases (ligand and cofactor rings on request).
- Screen many PDB entries in one batch.
- Build cluster and QM/MM models with ff14SB, ff19SB, or CHARMM36 charges and the RC, RCD, or Z1–Z3 boundary schemes.
- Write Gaussian inputs (tested) and Q-Chem EDA inputs (tested).
- Read Gaussian NBO output and Q-Chem ALMO-EDA output.
- Save and reload a full session as JSON, and write a complete HTML report and CSV tables.

**It cannot**

- Run quantum-chemistry jobs. You run them yourself.
- Read GROMACS `.xtc` or AMBER `.nc` files directly.
- Give a tested result for fluorine. Fluorine is off by default.
- Predict pKa values. Use PROPKA, Epik, or a similar tool and load the result.
- Prove a halogen bond. The NBO screen shows orbital support; it does not test the energy.
- Replace a free-energy study of binding.

**Tested routes**

| Task | Program | Status |
|---|---|---|
| NBO screen | Gaussian 16 with its built-in NBO 3.1 | tested on 59 dimers, 138 scan points, and 89 protein models |
| ALMO-EDA | Q-Chem 6.2.2 | tested |
| Batch screening | local PDB files or PDB IDs | tested on three test-set proteins |
| Crystal contacts and QM model | CIF with symmetry operators | tested on a small test CIF |
| Counterpoise binding | ORCA 5 or Gaussian 16 | used in the paper; not read by the app |
| NBO screen | ORCA with NBO 7 | untested; hidden by default |
| NBO screen | Q-Chem with NBO 5 | untested; hidden by default |

Source code, tests, and paper scripts: <https://github.com/compchemsk/BondX-Studio>.
