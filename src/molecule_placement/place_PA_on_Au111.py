import itertools, os
import numpy as np
from ase.io import read, write
from ase.geometry import find_mic
from ase.constraints import FixAtoms

# ---------------- user settings ----------------
SLAB_FILE = "/content/CONTCAR_Au"
MOLECULES = {                                   # label: path to optimized molecule
    "F5BzPA": "/content/CONTCAR_F5BzPA",
    "FBzPA":  "/content/CONTCAR_FBzPA",
}
OUT_ROOT    = "/content/POSCAR_PA_on_Au111"
D_AU_O      = 2.20    # target Au-O distance for bound O (Angstrom)
D_PO        = 1.52    # P-O bond length used to build target for P
D_PC        = 1.82    # P-C bond length used to build target for C
H_INTERACT  = 2.90    # height above top Au plane for the non-bonded O in mode (b)
VACUUM_TOP  = 12.0    # vacuum above the highest atom (Angstrom)
N_FREEZE    = 2       # number of bottom Au layers to freeze (0 = none)
H_POLICIES  = ("remove", "keep")  # remove H from bound O (P-O-Au) or keep it
MIN_CONTACT = 1.9     # minimum distance of non-bound atoms to Au after lifting (Angstrom)
MAX_LIFT    = 1.0     # largest upward shift allowed (Angstrom)
LIFT_STEP   = 0.05
# -----------------------------------------------

# ================= slab analysis (done once) =================
slab = read(SLAB_FILE)
pos = slab.positions
ztop = pos[:, 2].max()
top = np.where(pos[:, 2] > ztop - 0.5)[0]
second = np.where((pos[:, 2] < ztop - 0.5) & (pos[:, 2] > ztop - 3.0))[0]
dall = slab.get_all_distances(mic=True)
dt = dall[np.ix_(top, top)]
nn = dt[dt > 0.1].min()
print(f"Top layer atoms: {len(top)}, nearest-neighbor distance: {nn:.3f} A")

xy_center = slab.cell[:2, :2].sum(axis=0) / 2
anchor = top[np.argmin(np.linalg.norm(pos[top, :2] - xy_center, axis=1))]
vec = {j: slab.get_distance(anchor, j, mic=True, vector=True)
       for j in top if j != anchor and slab.get_distance(anchor, j, mic=True) < 1.1 * nn}
nbrs = list(vec.keys())

def site_z(r_inplane):
    return ztop + np.sqrt(max(D_AU_O**2 - r_inplane**2, 0.5**2))

Z_ATOP, Z_BRIDGE, Z_HOLLOW = site_z(0.0), site_z(nn / 2), site_z(nn / np.sqrt(3))

tri = {}   # fcc / hcp triangles containing the anchor atom
for j, k in itertools.combinations(nbrs, 2):
    if slab.get_distance(j, k, mic=True) < 1.1 * nn:
        cen = pos[anchor] + (vec[j] + vec[k]) / 3
        d = pos[second] - cen
        d[:, 2] = 0.0
        dmin = find_mic(d, slab.cell, slab.pbc)[1].min()
        kind = "hcp" if dmin < 0.3 else "fcc"
        tri.setdefault(kind, (j, k, cen))

def atop(i):
    p = pos[i].copy(); p[2] = Z_ATOP; return p

A = atop(anchor)
B = pos[anchor] + vec[nbrs[0]]; B[2] = Z_ATOP
bridge = pos[anchor] + vec[nbrs[0]] / 2; bridge[2] = Z_BRIDGE

site_lookup = {"atop": A, "bridge": bridge}
for kind in tri:
    c = tri[kind][2].copy(); c[2] = Z_HOLLOW
    site_lookup[kind] = c

# layers to freeze
zs = np.sort(np.unique(np.round(pos[:, 2], 1)))
layer_z = []
for z in zs:
    if not layer_z or z - layer_z[-1] > 0.5:
        layer_z.append(z)
fix_idx = [i for i in range(len(slab)) if N_FREEZE > 0 and pos[i, 2] < layer_z[N_FREEZE - 1] + 0.5]

# ================= generic helpers =================
def kabsch(P, Q, w):
    w = np.asarray(w, float)
    cp = (P * w[:, None]).sum(0) / w.sum()
    cq = (Q * w[:, None]).sum(0) / w.sum()
    H = (P - cp).T @ ((Q - cq) * w[:, None])
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cp, cq

def P_C_targets(O_targets):
    """P above the centroid of the O targets, C above P along +z."""
    T = np.array(O_targets)
    cxy = T[:, :2].mean(0)
    r = np.linalg.norm(T[:, :2] - cxy, axis=1).mean()
    dz = np.sqrt(max(D_PO**2 - r**2, 0.5**2))
    Pt = np.array([cxy[0], cxy[1], T[:, 2].mean() + dz])
    Ct = Pt + np.array([0, 0, D_PC])
    return Pt, Ct

# ================= per-molecule builder =================
def process(label, mol_file):
    out_dir = os.path.join(OUT_ROOT, label)
    os.makedirs(out_dir, exist_ok=True)
    mol = read(mol_file)
    mol.set_pbc(False)

    # ---- molecule bookkeeping ----
    sym = np.array(mol.get_chemical_symbols())
    iP = int(np.where(sym == "P")[0][0])
    dP = np.linalg.norm(mol.positions - mol.positions[iP], axis=1)
    O_idx = [int(i) for i in np.where(sym == "O")[0] if dP[i] < 2.0]
    iC = int([i for i in np.where(sym == "C")[0] if dP[i] < 2.2][0])
    assert len(O_idx) == 3, f"{label}: expected 3 O on P, found {len(O_idx)} (is the molecule unwrapped?)"

    H_of = {}
    for o in O_idx:
        do = np.linalg.norm(mol.positions - mol.positions[o], axis=1)
        H_of[o] = [int(h) for h in np.where(sym == "H")[0] if do[h] < 1.2]
    oh_O = [o for o in O_idx if H_of[o]]
    po_O = [o for o in O_idx if not H_of[o]]
    assert len(oh_O) == 2 and len(po_O) == 1, f"{label}: expected 2 P-OH and 1 P=O"
    print(f"\n[{label}] P: {iP}  C: {iC}  P-OH O: {oh_O}  P=O O: {po_O}")

    # ---- fitting helpers that depend on this molecule ----
    def fit(idx, targets, weights):
        P = mol.positions[idx]
        Q = np.array(targets)
        R, cp, cq = kabsch(P, Q, weights)
        new = (mol.positions - cp) @ R.T + cq
        rmsd = np.sqrt(np.average(((new[idx] - Q) ** 2).sum(1), weights=weights))
        return new, rmsd

    def clash_metric(new_pos, exclude):
        """Smallest distance between molecule atoms (except those in `exclude`) and slab atoms."""
        tmp = slab.copy()
        m = mol.copy(); m.positions = new_pos
        tmp += m
        D = tmp.get_all_distances(mic=True)
        ns = len(slab)
        rows = [ns + i for i in range(len(mol)) if i not in exclude]
        return D[np.ix_(rows, range(ns))].min()

    def best_assignment(O_list, O_targets, weights_O):
        best = None
        Pt, Ct = P_C_targets(O_targets)
        for perm in itertools.permutations(range(len(O_targets))):
            tg = [O_targets[perm[p]] for p in range(len(O_list))]
            new, rmsd = fit(list(O_list) + [iP, iC], tg + [Pt, Ct], list(weights_O) + [1.0, 1.0])
            if best is None or rmsd < best[1]:
                best = (new, rmsd)
        return best

    # ---- mode builders: return (new_positions, bound_O_list, rmsd) ----
    def mono(site_pos):
        O_b = oh_O[0]
        best = None
        for phi in np.deg2rad(np.arange(0, 360, 30)):
            d_PO = np.array([np.sin(np.deg2rad(109.5)) * np.cos(phi),
                             np.sin(np.deg2rad(109.5)) * np.sin(phi),
                             np.cos(np.deg2rad(109.5))])       # P -> O direction
            Pt = site_pos - D_PO * d_PO
            Ct = Pt + np.array([0, 0, D_PC])
            new, rmsd = fit([O_b, iP, iC], [site_pos, Pt, Ct], [1, 1, 1])
            c = clash_metric(new, {O_b})
            if best is None or c > best[0]:
                best = (c, new, rmsd)
        return best[1], [O_b], best[2]

    def bidentate():
        new, rmsd = best_assignment(oh_O, [A, B], [1, 1])
        return new, list(oh_O), rmsd

    def bi_electro():
        O_b, O_i = oh_O[0], po_O[0]
        tgt_b = A
        tgt_i = B.copy(); tgt_i[2] = ztop + H_INTERACT
        Pt, Ct = P_C_targets([tgt_b, tgt_i])
        new, rmsd = fit([O_b, O_i, iP, iC], [tgt_b, tgt_i, Pt, Ct], [2, 1, 1, 1])
        return new, [O_b], rmsd

    def tridentate(kind):
        j, k, _ = tri[kind]
        T = [atop(anchor), atop(j), atop(k)]
        new, rmsd = best_assignment(O_idx, T, [1, 1, 1])
        return new, list(O_idx), rmsd

    jobs = [(f"mono_{s}", (lambda s=s: mono(site_lookup[s]))) for s in site_lookup]
    jobs.append(("bi_atop-atop", bidentate))
    jobs.append(("bi-electro_atop", bi_electro))
    jobs += [(f"tri_{kd}", (lambda kd=kd: tridentate(kd))) for kd in tri]

    # ---- build and write ----
    summary = []
    for mode, builder in jobs:
        new_pos, bound, rmsd = builder()
        for policy in H_POLICIES:
            drop = {h for o in bound for h in H_of[o]} if policy == "remove" else set()
            exclude = set(bound) | drop

            # lift the molecule until non-bound atoms are at least MIN_CONTACT from Au
            lift = 0.0
            while (clash_metric(new_pos + np.array([0, 0, lift]), exclude) < MIN_CONTACT
                   and lift < MAX_LIFT):
                lift += LIFT_STEP
            placed = new_pos + np.array([0, 0, lift])
            contact = clash_metric(placed, exclude)

            m = mol.copy()
            m.positions = placed
            if drop:
                m = m[[i for i in range(len(m)) if i not in drop]]
            s2 = slab.copy()
            s2 += m
            zmax = s2.positions[:, 2].max()
            c_new = max(slab.cell[2, 2], zmax + VACUUM_TOP)
            cell = s2.cell.array.copy(); cell[2] = [0, 0, c_new]
            s2.set_cell(cell, scale_atoms=False)
            s2.set_constraint(FixAtoms(indices=fix_idx))

            Dm = s2.get_all_distances(mic=True)
            Dno = s2.get_all_distances(mic=False)
            ns = len(slab)
            dmin_slab = Dm[ns:, :ns].min()
            sub, subn = Dm[ns:, ns:], Dno[ns:, ns:]
            img = sub[sub < subn - 0.01]
            dimg = img.min() if img.size else np.inf

            name = f"POSCAR_{label}_{mode}_H{policy}"
            write(os.path.join(out_dir, name), s2, format="vasp", direct=True, vasp5=True, sort=True)
            summary.append((name, rmsd, lift, contact, dmin_slab, dimg))
    return summary

# ================= run all molecules =================
for label, path in MOLECULES.items():
    summ = process(label, path)
    print(f"{'file':42s} {'fit rmsd':>9s} {'lift':>6s} {'min d(non-bound)':>17s} {'min d(all)':>11s} {'min d(image)':>13s}")
    for n, r, lf, c, d1, d2 in summ:
        flag = "  <-- check" if (c < MIN_CONTACT - 0.05 or d2 < 2.0) else ""
        print(f"{n:42s} {r:9.2f} {lf:6.2f} {c:17.2f} {d1:11.2f} {d2:13.2f}{flag}")
