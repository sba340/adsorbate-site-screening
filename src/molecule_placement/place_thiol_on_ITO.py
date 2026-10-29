import itertools, os
import numpy as np
from ase.io import read, write
from ase.geometry import find_mic, get_distances
from ase.constraints import FixAtoms

# =====================================================================
# Thiols (PFBT, FBzS) on hydroxylated ITO(222). Same slab model and anchor
# search as place_PA_on_ITO.py, with the orientation scan of the Au script.
# Anchor types:
#   "M" : exposed, under-coordinated In/Sn. S placed D_M_S from the cation along
#         its open direction. H policies as in the Au script:
#           "thiolate" = S-H hydrogen removed, "thiol" = intact S-H.
#   "OH": surface hydroxyl, thiol condensation  R-SH + HO-M -> R-S-M + H2O.
#         S replaces the hydroxyl O (moved outward until the mean M-S distance is
#         D_M_S); the hydroxyl (O+H) and the S-H hydrogen are removed as water.
#         Switch off with USE_OH_CONDENSATION = False.
# =====================================================================

# ---------------- user settings ----------------
SLAB_FILE = "/content/CONTCAR"                  # ITO slab
MOLECULES = {
    "PFBT": "/content/CONTCAR_PFBT",   # adjust to your actual file names
    "FBzS": "/content/CONTCAR_FBzS",
}
OUT_ROOT   = "/content/POSCAR_thiol_on_ITO"
METALS     = ("In", "Sn")
D_M_S      = {"In": 2.55, "Sn": 2.45}        # target M-S distance (Angstrom)
USE_OH_CONDENSATION = True
TILTS      = (0.0, 30.0)                     # tilt of the S->C bond from the surface normal (degrees)
PHIS       = tuple(np.arange(0, 360, 60))    # azimuth of the tilt (degrees), used when tilt > 0
PSIS       = tuple(np.arange(0, 360, 30))    # spin about the S->C bond (degrees)
OH_CUTOFF   = 1.20    # O-H distance that defines a surface hydroxyl (A)
MO_CUTOFF   = 2.70    # M-O cutoff for coordination / open direction (A)
SURF_DEPTH  = 2.5     # metals within this depth of the highest metal are surface candidates (A)
MIN_OPEN    = 0.35
MIN_UP      = 0.30
CLEAR_HEAVY = 2.50    # open-site S must be this far from other slab heavy atoms (A)
CLEAR_H     = 2.00    # ... and this far from slab H atoms (A)
N_OH_ANCHORS = 8      # surface hydroxyls used as anchors (highest first)
N_M_ANCHORS  = 8      # open metal sites used as anchors (most exposed first)
VACUUM_TOP = 12.0
FREEZE_BELOW = 3.0    # freeze slab atoms with z below this value (A); 0 = none
# -----------------------------------------------

# ================= slab analysis (done once) =================
slab = read(SLAB_FILE)
slab.set_pbc(True)
pos = slab.positions
sym_s = np.array(slab.get_chemical_symbols())
ns = len(slab)
CELL, PBC = slab.cell, slab.pbc
metal_idx = np.where(np.isin(sym_s, METALS))[0]
O_idx_s = np.where(sym_s == "O")[0]
H_idx_s = np.where(sym_s == "H")[0]
zmetal = pos[metal_idx, 2].max()
print(f"Slab: {slab.get_chemical_formula()}  highest metal z = {zmetal:.2f}")
Dss = slab.get_all_distances(mic=True)

def mic_vec(a, b):
    return find_mic((b - a)[None, :], CELL, PBC)[0][0]

def unwrap_to(ref, p):
    return ref + mic_vec(ref, p)

# ---- anchor type "OH": surface hydroxyls (thiol condensation) ----
anchors = []
if USE_OH_CONDENSATION:
    hyd_O = {}
    for h in H_idx_s:
        o = int(O_idx_s[np.argmin(Dss[h, O_idx_s])])
        if Dss[h, o] < OH_CUTOFF:
            hyd_O[o] = int(h)
    cand = []
    for o, h in hyd_O.items():
        ms = [int(m) for m in metal_idx if Dss[o, m] < MO_CUTOFF]
        if not ms:
            continue
        mp = [unwrap_to(pos[o], pos[m]) for m in ms]
        cen = np.mean(mp, axis=0)
        d = pos[o] - cen
        d = d / np.linalg.norm(d) if np.linalg.norm(d) > 1e-6 else np.array([0, 0, 1.0])
        target = np.mean([D_M_S[sym_s[m]] for m in ms])
        best_t, best_e = 0.0, 1e9
        for t in np.arange(0.0, 2.01, 0.05):          # move S outward until mean M-S = target
            e = abs(np.mean([np.linalg.norm(pos[o] + t * d - q) for q in mp]) - target)
            if e < best_e:
                best_t, best_e = t, e
        p = pos[o] + best_t * d
        dq = np.linalg.norm(find_mic(pos - p, CELL, PBC)[0], axis=1)
        if all(dq[j] >= (1.8 if sym_s[j] == "H" else 2.0) for j in range(ns) if j not in (o, h)):
            cand.append(dict(kind="OH", tag=f"OH{o}", pos=p, dir=d, leave=[o, h], score=pos[o, 2]))
    cand.sort(key=lambda a: -a["score"])
    anchors += cand[:N_OH_ANCHORS]
    print(f"Surface hydroxyls found: {len(hyd_O)} (anchors kept: {min(len(cand), N_OH_ANCHORS)})")

# ---- anchor type "M": exposed, unobstructed under-coordinated cations ----
candM = []
for m in metal_idx:
    if pos[m, 2] < zmetal - SURF_DEPTH:
        continue
    vecs = slab.get_distances(int(m), range(ns), mic=True, vector=True)
    dist = np.linalg.norm(vecs, axis=1)
    nb = [j for j in range(ns) if j != m and sym_s[j] == "O" and dist[j] < MO_CUTOFF]
    if not nb:
        continue
    u = vecs[nb] / dist[nb][:, None]
    ov = -u.sum(axis=0)
    mag = np.linalg.norm(ov)
    if mag < MIN_OPEN:
        continue
    d = ov / mag
    if d[2] < MIN_UP:
        continue
    el = sym_s[m]
    p = pos[m] + D_M_S[el] * d
    dq = np.linalg.norm(find_mic(pos - p, CELL, PBC)[0], axis=1)
    if all(dq[j] >= (CLEAR_H if sym_s[j] == "H" else CLEAR_HEAVY) for j in range(ns) if j != m):
        candM.append(dict(kind="M", tag=f"{el}{int(m)}", pos=p, dir=d, leave=[], score=mag + 0.05 * pos[m, 2]))
candM.sort(key=lambda a: -a["score"])
anchors += candM[:N_M_ANCHORS]
print(f"Open metal sites kept: {min(len(candM), N_M_ANCHORS)}")
assert anchors, "No anchors found; relax the site-finding settings."
print("Anchors:", [a["tag"] for a in anchors])

fix_all = lambda arr: [i for i in range(len(arr)) if FREEZE_BELOW > 0 and arr[i, 2] < FREEZE_BELOW]

# ================= rotation helpers =================
def skew(v):
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

def rot_axis(axis, ang):
    axis = np.asarray(axis, float); axis = axis / np.linalg.norm(axis)
    K = skew(axis)
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * (K @ K)

def rot_align(a, b):
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); c = float(np.dot(a, b)); s = np.linalg.norm(v)
    if s < 1e-8:
        if c > 0:
            return np.eye(3)
        ax = np.cross(a, [1, 0, 0])
        if np.linalg.norm(ax) < 1e-6:
            ax = np.cross(a, [0, 1, 0])
        return rot_axis(ax, np.pi)
    K = skew(v)
    return np.eye(3) + K + (K @ K) * ((1 - c) / s**2)

# ================= per-molecule builder =================
def process(label, mol_file):
    out_dir = os.path.join(OUT_ROOT, label)
    os.makedirs(out_dir, exist_ok=True)
    mol = read(mol_file)
    mol.set_pbc(False)

    sym = np.array(mol.get_chemical_symbols())
    assert (sym == "S").sum() == 1, f"{label}: expected exactly one S"
    iS = int(np.where(sym == "S")[0][0])
    dS = np.linalg.norm(mol.positions - mol.positions[iS], axis=1)
    Cs = [int(i) for i in np.where(sym == "C")[0] if dS[i] < 2.0]
    assert len(Cs) == 1, f"{label}: expected one C bonded to S, found {len(Cs)}"
    iC = Cs[0]
    H_S = [int(h) for h in np.where(sym == "H")[0] if dS[h] < 1.6]
    print(f"\n[{label}] S: {iS}  C: {iC}  H on S: {H_S}")

    u0 = mol.positions[iC] - mol.positions[iS]
    summary = []

    for anc in anchors:
        spos = anc["pos"]
        leave = set(anc["leave"])
        keep_idx = [i for i in range(ns) if i not in leave]
        s_r = slab[keep_idx]
        sxyz = s_r.positions
        if anc["kind"] == "OH":
            policies = ("H2O",) if H_S else ()      # condensation needs the S-H hydrogen
        else:
            policies = ("thiolate", "thiol") if H_S else ("thiolate",)
        for tilt in TILTS:
            for policy in policies:
                drop = set(H_S) if policy in ("H2O", "thiolate") else set()
                excl = {iS} | drop
                rows = [i for i in range(len(mol)) if i not in excl]
                best = None
                for phi in (PHIS if tilt > 0 else (0.0,)):
                    t = np.array([np.sin(np.deg2rad(tilt)) * np.cos(np.deg2rad(phi)),
                                  np.sin(np.deg2rad(tilt)) * np.sin(np.deg2rad(phi)),
                                  np.cos(np.deg2rad(tilt))])
                    R1 = rot_align(u0, t)
                    for psi in PSIS:
                        R = rot_axis(t, np.deg2rad(psi)) @ R1
                        newpos = (mol.positions - mol.positions[iS]) @ R.T + spos
                        c = get_distances(newpos[rows], sxyz, cell=CELL, pbc=PBC)[1].min()
                        if best is None or c > best[0]:
                            best = (c, newpos)
                m = mol.copy()
                m.positions = best[1]
                if drop:
                    m = m[[i for i in range(len(m)) if i not in drop]]
                s2 = s_r.copy()
                s2 += m
                nsr = len(s_r)
                zmax = s2.positions[:, 2].max()
                c_new = max(slab.cell[2, 2], zmax + VACUUM_TOP)
                cell = s2.cell.array.copy(); cell[2] = [0, 0, c_new]
                s2.set_cell(cell, scale_atoms=False)
                s2.set_constraint(FixAtoms(indices=fix_all(s2.positions[:nsr])))
                Dm = s2.get_all_distances(mic=True)
                Dno = s2.get_all_distances(mic=False)
                sub, subn = Dm[nsr:, nsr:], Dno[nsr:, nsr:]
                img = sub[sub < subn - 0.01]
                dimg = img.min() if img.size else np.inf
                dMS = get_distances([spos], sxyz, cell=CELL, pbc=PBC)[1].min()
                nw = 1 if policy == "H2O" else 0
                name = f"POSCAR_{label}_{anc['tag']}_tilt{int(tilt)}_{policy}"
                write(os.path.join(out_dir, name), s2, format="vasp", direct=True, vasp5=True, sort=True)
                summary.append((name, nw, dMS, best[0], dimg))
    return summary

for label, path in MOLECULES.items():
    summ = process(label, path)
    print(f"{'file':50s} {'H2O':>3s} {'d(M-S)':>7s} {'min d(other-slab)':>18s} {'min d(image)':>13s}")
    for n, nw, d0, d1, d2 in summ:
        flag = "  <-- check" if (d1 < 2.0 or d2 < 2.0 or d0 < 2.0) else ""
        print(f"{n:50s} {nw:3d} {d0:7.2f} {d1:18.2f} {d2:13.2f}{flag}")
