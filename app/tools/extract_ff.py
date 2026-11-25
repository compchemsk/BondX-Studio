#!/usr/bin/env python3
"""
extract_ff.py -- builds src/ff_data.js for BondX Studio from the OpenMM force-field XML files
(openmm/wrappers/python/openmm/app/data, git checkout) and validates every residue.

Force fields written:
  ff14SB   : AMBER ff14SB protein + OL15 DNA + OL3 RNA + TIP3P/SPC-E/OPC3 water + Joung-Cheatham/Li-Merz ions
  ff19SB   : AMBER ff19SB protein + OL15/OL3 + same water/ions (charges compared with ff14SB)
  charmm36 : CHARMM36 protein (patches NTER/GLYP/PROP/CTER/ASPP/GLUP/LSN/DISU applied here) + CHARMM TIP3P + ions

For every residue variant we store, keyed by PDB (v3) heavy-atom names:
   atoms[name] = [q_heavy, type, [[H_name, q_H], ...]]
listing the force-field hydrogens bonded to that heavy atom (names as in the force field).
Net charge of every variant must be an integer within 1e-4, otherwise it aborts.
"""
import json, sys, copy, xml.etree.ElementTree as ET
from pathlib import Path

D = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/omm/wrappers/python/openmm/app/data")
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "src/ff_data.js")
log = []

def elem_of(type_name, types):
    return types.get(type_name, "")

def load_types(root):
    t = {}
    for ty in root.iter("Type"):
        t[ty.get("name")] = ty.get("element") or ""
    return t

def load_lj(root, name_attr="type"):
    lj = {}
    for nb in root.iter("NonbondedForce"):
        for a in nb.findall("Atom"):
            k = a.get("type") or a.get("class")
            if a.get("sigma") is not None:
                lj[k] = [round(float(a.get("sigma")) * 10, 5), round(float(a.get("epsilon")), 6)]  # A, kJ/mol
    return lj

def residues(root):
    out = {}
    for r in root.iter("Residue"):
        atoms = [(a.get("name"), float(a.get("charge")), a.get("type")) for a in r.findall("Atom")]
        bonds = [(b.get("atomName1"), b.get("atomName2")) for b in r.findall("Bond")]
        out[r.get("name")] = {"atoms": atoms, "bonds": bonds}
    return out

def fold(resname, atoms, bonds, types, rename=None, allow_frac=False):
    """heavy-atom table with hydrogen charges folded by parent."""
    rename = rename or {}
    el = {n: elem_of(t, types) for n, q, t in atoms}
    q = {n: c for n, c, t in atoms}
    ty = {n: t for n, c, t in atoms}
    parent = {}
    for a, b in bonds:
        if el.get(a) == "H" and el.get(b) not in ("H", None): parent[a] = b
        if el.get(b) == "H" and el.get(a) not in ("H", None): parent[b] = a
    table = {}
    for n, c, t in atoms:
        if el[n] == "H": continue
        table[n] = [c, t, []]
    for n, c, t in atoms:
        if el[n] != "H": continue
        if n not in parent: raise SystemExit(f"{resname}: hydrogen {n} without parent")
        table[parent[n]][2].append([n, round(q[n], 6)])
    out = {}
    for n, rec in table.items():
        out[rename.get(n, n)] = [round(rec[0], 6), rec[1].split("-", 1)[-1] if rec[1] else "", rec[2]]
    net = sum(c for n, c, t in atoms)
    if abs(net - round(net)) > 1e-4 and not allow_frac: raise SystemExit(f"{resname}: non-integer net charge {net}")
    return {"q": round(net, 4) if allow_frac else int(round(net)), "atoms": out, "sum_exact": round(net, 6)}

# ------------------------------------------------------------------ AMBER
def amber(protein_xml, label):
    root = ET.parse(D / protein_xml).getroot(); types = load_types(root); res = residues(root)
    out = {}
    for name, r in res.items():
        out[name] = fold(name, r["atoms"], r["bonds"], types)
    return out, load_lj(root), root.find("Info/Reference").text if root.find("Info/Reference") is not None else ""

def nucleic():
    out = {}; lj = {}
    for f in ["amber14/DNA.OL15.xml", "amber14/RNA.OL3.xml"]:
        root = ET.parse(D / f).getroot(); types = load_types(root)
        for name, r in residues(root).items():
            if name in ("OHE",): continue
            out[name] = fold(name, r["atoms"], r["bonds"], types, allow_frac=name[-1] in "35")
        lj.update(load_lj(root))
    # 5'- and 3'-terminal nucleotides carry fractional charges; any 5' + 3' pair must sum to -1
    for a in [k for k in out if k.endswith("5")]:
        for b in [k for k in out if k.endswith("3") and k.startswith("D") == a.startswith("D")]:
            s_ = out[a]["sum_exact"] + out[b]["sum_exact"]
            if abs(s_ - round(s_)) > 1e-4 or round(s_) != -1: raise SystemExit(f"nucleic termini {a}+{b} = {s_}")
    log.append("nucleic: every 5'-terminal + 3'-terminal pair sums to -1 (checked)")
    return out, lj

def water_ions(xml_file, water_key):
    root = ET.parse(D / xml_file).getroot(); types = load_types(root); out = {}; ions = {}
    for name, r in residues(root).items():
        if len(r["atoms"]) == 1:
            ions[name.upper()] = round(r["atoms"][0][1], 6)
        elif name == "HOH":
            out[water_key] = fold(water_key, r["atoms"], r["bonds"], types)
    return out, ions, load_lj(root)

# ------------------------------------------------------------------ CHARMM36 (patches applied here)
def charmm():
    root = ET.parse(D / "charmm36.xml").getroot(); types = load_types(root)
    res = residues(root)
    patches = {p.get("name"): p for p in root.iter("Patch")}
    def apply(base, *pnames, prefix=None):
        atoms = {n: [q, t] for n, q, t in base["atoms"]}; order = [n for n, q, t in base["atoms"]]
        bonds = set(tuple(b) for b in base["bonds"])
        for pn in pnames:
            p = patches[pn]
            for ch in p:
                nm = ch.get("name") or ""
                if prefix:
                    if not nm.startswith(prefix + ":") and ch.tag in ("ChangeAtom", "AddAtom", "RemoveAtom"): continue
                    nm = nm.split(":", 1)[-1]
                if ch.tag == "ChangeAtom": atoms[nm] = [float(ch.get("charge")), ch.get("type")]
                elif ch.tag == "AddAtom": atoms[nm] = [float(ch.get("charge")), ch.get("type")]; order.append(nm)
                elif ch.tag == "RemoveAtom": atoms.pop(nm, None); order = [o for o in order if o != nm]
                elif ch.tag == "AddBond":
                    a, b = ch.get("atomName1"), ch.get("atomName2")
                    if prefix:
                        if not (a.startswith(prefix + ":") and b.startswith(prefix + ":")): continue
                        a, b = a.split(":")[1], b.split(":")[1]
                    bonds.add((a, b))
                elif ch.tag == "RemoveBond":
                    a, b = ch.get("atomName1"), ch.get("atomName2"); bonds.discard((a, b)); bonds.discard((b, a))
        bonds = [b for b in bonds if b[0] in atoms and b[1] in atoms]
        return {"atoms": [(n, atoms[n][0], atoms[n][1]) for n in order if n in atoms], "bonds": bonds}
    ren_common = {"OT1": "O", "OT2": "OXT"}
    out = {}
    base_map = {"HSD": "HID", "HSE": "HIE", "HSP": "HIP"}
    std = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    for n in std:
        key = base_map.get(n, n); ren = dict(ren_common)
        if n == "ILE": ren["CD"] = "CD1"
        out[key] = fold(key, res[n]["atoms"], res[n]["bonds"], types, ren)
        nt = "GLYP" if n == "GLY" else "PROP" if n == "PRO" else "NTER"
        r = apply(res[n], nt); out["N" + key] = fold("N" + key, r["atoms"], r["bonds"], types, ren)
        r = apply(res[n], "CTER"); out["C" + key] = fold("C" + key, r["atoms"], r["bonds"], types, ren)
    for pn, base, key in [("ASPP", "ASP", "ASH"), ("GLUP", "GLU", "GLH"), ("LSN", "LYS", "LYN")]:
        r = apply(res[base], pn); out[key] = fold(key, r["atoms"], r["bonds"], types)
    out["CYM"] = fold("CYM", res["CYM"]["atoms"], res["CYM"]["bonds"], types)
    r = apply(res["CYS"], "DISU", prefix="1"); out["CYX"] = fold("CYX", r["atoms"], r["bonds"], types)
    for key in ["CYX"]:
        for t, p in [("N", "NTER"), ("C", "CTER")]:
            r = apply(apply(res["CYS"], "DISU", prefix="1"), p); out[t + key] = fold(t + key, r["atoms"], r["bonds"], types, ren_common)
    # water + ions
    wroot = ET.parse(D / "charmm36/water.xml").getroot(); wtypes = load_types(wroot); ions = {}
    for name, r in residues(wroot).items():
        if name == "TIP3":
            out_w = fold("TIP3", [( {"OH2": "O"}.get(a, a), q, t) for a, q, t in r["atoms"]], [({"OH2": "O"}.get(a, a), {"OH2": "O"}.get(b, b)) for a, b in r["bonds"]], wtypes)
        elif len(r["atoms"]) == 1: ions[name] = r["atoms"][0][1]
    lj = load_lj(root); lj.update(load_lj(wroot))
    return out, {"TIP3P (CHARMM)": out_w}, ions, lj

def main():
    FF = {}
    a14, lj14, ref14 = amber("amber14/protein.ff14SB.xml", "ff14SB")
    a19, lj19, ref19 = amber("amber19/protein.ff19SB.xml", "ff19SB")
    nuc, ljn = nucleic()
    w1, ions_a, ljw = water_ions("amber14/tip3p.xml", "TIP3P")
    w2, _, ljw2 = water_ions("amber14/spce.xml", "SPC/E")
    w3, _, ljw3 = water_ions("amber14/opc3.xml", "OPC3")
    waters_a = {**w1, **w2, **w3}
    for k, v in nuc.items(): a14.setdefault(k, copy.deepcopy(v)); a19.setdefault(k, copy.deepcopy(v))
    # compare ff14SB / ff19SB charges
    diffs = []
    for k in a14:
        if k in a19:
            for an, rec in a14[k]["atoms"].items():
                r2 = a19[k]["atoms"].get(an)
                if r2 is None or abs(r2[0] - rec[0]) > 1e-6 or [round(h[1],6) for h in r2[2]] != [round(h[1],6) for h in rec[2]]: diffs.append(f"{k}:{an}")
    log.append(f"ff14SB vs ff19SB: {len(diffs)} atom charge differences ({', '.join(diffs[:12])}{' ...' if len(diffs) > 12 else ''})")
    c36, wc, ions_c, ljc = charmm()
    ions_a = {k: v for k, v in ions_a.items()}
    FF["ff14SB"] = {"label": "AMBER ff14SB (+ OL15 DNA, OL3 RNA)", "ref": ref14, "residues": a14, "waters": waters_a, "ions": ions_a, "lj": {**lj14, **ljn, **ljw, **ljw2, **ljw3}}
    FF["ff19SB"] = {"label": "AMBER ff19SB (+ OL15 DNA, OL3 RNA)", "ref": ref19, "residues": a19, "waters": copy.deepcopy(waters_a), "ions": ions_a, "lj": {**lj19, **ljn, **ljw, **ljw2, **ljw3}}
    FF["charmm36"] = {"label": "CHARMM36 (protein; nucleic acids fall back to AMBER OL15/OL3)", "ref": "Huang, J. et al. CHARMM36m. Nat. Methods 2017, 14, 71-73; Best, R. B. et al. J. Chem. Theory Comput. 2012, 8, 3257-3273 (CHARMM36 protein).",
                      "residues": {**{k: copy.deepcopy(v) for k, v in nuc.items()}, **c36}, "waters": wc, "ions": ions_c, "lj": {**ljc, **ljn}}
    # validation summary
    for name, ff in FF.items():
        n = len(ff["residues"]); bad = [k for k, v in ff["residues"].items() if abs(v["sum_exact"] - v["q"]) > 1e-4 and not (k[-1] in "35" and k[:-1] in ("DA","DC","DG","DT","A","C","G","U"))]
        log.append(f"{name}: {n} residue variants, waters {list(ff['waters'])}, {len(ff['ions'])} ions, {len(ff['lj'])} LJ types, non-integer: {bad}")
        for v in ff["residues"].values(): v.pop("sum_exact")
        for v in ff["waters"].values(): v.pop("sum_exact", None)
    meta = {"generated_from": "OpenMM app/data force-field XML (github.com/openmm/openmm, master)", "units": {"charge": "e", "sigma": "Angstrom", "epsilon": "kJ/mol"}}
    js = "/* Generated by tools/extract_ff.py -- do not edit by hand.\n * " + "\n * ".join(log) + "\n */\n(function (root) {\n  const FF = " + json.dumps(FF, separators=(",", ":")) + ";\n  FF._meta = " + json.dumps(meta) + ";\n  if (typeof module === 'object' && module.exports) module.exports = FF; else root.BondXFF = FF;\n})(typeof self !== 'undefined' ? self : this);\n"
    OUT.parent.mkdir(parents=True, exist_ok=True); OUT.write_text(js)
    print("\n".join(log)); print(f"wrote {OUT} ({len(js)/1024:.0f} kB)")

if __name__ == "__main__":
    main()
