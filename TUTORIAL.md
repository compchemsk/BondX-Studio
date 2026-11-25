# BondX Studio tutorial

This is the five-page tutorial of BondX Studio. The same text opens with the **Tutorial** button inside the program.

---

## Page 1 · What BondX Studio does

BondX Studio checks halogen bonds (XBs) in proteins, ligands and crystals. It runs in your browser. No file leaves your computer. The only network call is the download of a PDB entry from RCSB.

The analysis is split into tiers. Each tier adds one kind of evidence.

| Tier | Step | Evidence |
|---|---|---|
| T0 | 1 Load structure | Atoms, chains, hydrogens, charges, virtual sites |
| T1 | 2 Find halogen bonds | Geometry: distance, two angles, geometric score *G* |
| T2 | 3 Dynamics | Occupancy of each contact over a trajectory |
| T3 | 4–5 QM model, NBO gate | Cluster or QM/MM input files; NBO donor–acceptor energy ΣE(2) |
| T4 | 6 ALMO-EDA | Frozen, polarisation and charge-transfer energies |
| T5 | 7 Score & report | Quality index BXQI and grade A–D |

**Quick start.** Open **1 Load structure**. Press **Load demo**. Then press the button at the bottom right of each page. The demo takes less than one minute.

You can stop after any tier. A result from T1 alone is valid. It is then marked as geometry only.

---

## Page 2 · Load a structure and find contacts

### Step 1 · Load structure (T0)

- Type a four-letter PDB ID and press **Fetch**. Or drop a local `.pdb`, `.ent` or `.pqr` file.
- A PQR file brings its own partial charges. A PDB file uses the force-field charges of the program.
- Pick the chains to keep. Check the composition table.
- Extra-point atoms (EP, LP) of σ-hole force fields are removed on load. The message says how many were removed.

### Step 2 · Find halogen bonds (T1)

- The program finds every C–X···Y contact with X = Cl, Br or I.
- It measures the X···Y distance *d*, the C–X···Y angle θ1 and the X···Y–R angle θ2.
- It reports the score *G* between 0 and 1. A larger *G* means a better geometry.
- Each acceptor is classed as a Lewis base or not. A backbone N–H has no lone pair. Such a contact is flagged as weak.
- Sort or filter the table. Click a row to see the contact in the 3D viewer.

All cut-offs are on the page **All parameters**. A changed value is marked in orange. It is also written to the report.

---

## Page 3 · Dynamics and the QM model

### Step 3 · Dynamics (T2, optional)

- Drop a topology `.pdb`. Its atom order must match the trajectory.
- Drop a `.dcd`, a multi-model `.pdb` or a multi-frame `.xyz` file.
- Convert GROMACS or AMBER files once with `gmx trjconv` or `cpptraj`.
- The program applies the T1 test to every frame. It reports the occupancy of each contact.

### Step 4 · QM model and charges (T3)

- Choose a **cluster** model or a **QM/MM** model with point charges.
- Cut bonds are capped with hydrogen atoms. Ligand charges can come from a `.mol2` file.
- Pick the program: Gaussian, ORCA or Q-Chem. Pick the method, the basis set and the solvent model.
- Download the input files. Each file carries a manifest ID. The ID links the output back to the contact.
- Run the jobs on your own computer or cluster. BondX Studio does not run QM jobs.

---

## Page 4 · NBO gate, ALMO-EDA and the score

### Step 5 · NBO gate (T3)

- Drop the QM output (`.log` or `.out`). The program reads the second-order perturbation table.
- It adds the terms LP(Y) → σ*(C–X) and LP(Y) → RY*(X). For a π acceptor it adds the ring bonds.
- Intramolecular terms are excluded. Only the T1 acceptor counts.
- A contact passes the gate if ΣE(2) ≥ 0.5 kcal/mol.
- For an ionic NBO structure the halogen has a vacant LP*. Donation into LP*(X) is then counted.

### Step 6 · ALMO-EDA (T4)

- Drop the Q-Chem EDA output. The program reads FRZ, POL, CT and the total energy.
- It reports χCT = CT / (POL + CT).

### Step 7 · Score and report (T5)

- BXQI is a weighted geometric mean of the tier scores. It runs from 0 to 100.
- Grades: A ≥ 75, B ≥ 50, C ≥ 25, else D. Grade A needs a passed NBO gate.
- Press **Report** for a single HTML file. It lists every input, parameter and result.

---

## Page 5 · What it can and cannot do

**It can**

- Read PDB, PQR, CIF, XYZ, DCD and multi-model PDB files.
- Find Cl, Br and I halogen bonds to O, N, S and aromatic rings.
- Screen many PDB entries in one batch.
- Build cluster and QM/MM models with ff14SB, ff19SB or CHARMM36 charges.
- Write Gaussian, ORCA and Q-Chem inputs.
- Read Gaussian NBO output and Q-Chem ALMO-EDA output.
- Save and reload a full session as JSON.
- Write a complete HTML report and CSV tables.

**It cannot**

- Run quantum-chemistry jobs. You run them yourself.
- Read GROMACS `.xtc` or AMBER `.nc` files directly.
- Give a tested result for fluorine. Fluorine is off by default.
- Replace a full free-energy study of binding.
- Decide on a weak contact without the NBO tier.

**Tested routes**

| Task | Program | Status |
|---|---|---|
| NBO gate | Gaussian 16 with NBO 3.1 | tested on 59 dimers and 9 protein contacts |
| ALMO-EDA | Q-Chem 6.2 | tested |
| Counterpoise binding | ORCA 5 or Gaussian 16 | used in the paper; not read by the app |
| NBO gate | ORCA with NBO | input written; not tested |
| NBO gate | Q-Chem with NBO 5 | input written; not tested |

Source code, tests and paper scripts: <https://github.com/compchemsk/BondX-Studio>.
