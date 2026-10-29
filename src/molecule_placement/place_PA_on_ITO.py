import itertools, os
import numpy as np
from ase.io import read, write
from ase.geometry import find_mic, get_distances
from ase.constraints import FixAtoms

# =====================================================================
# Phosphonic acid (PA) on hydroxylated ITO(222), following the model of
# Paramonov et al., Chem. Mater. 2008, 20, 5131 (Scheme 1):
#   (a) monodentate         : 1 P-O-M bond
#   (b) bidentate + electro.: 1 P-O-M bond + P-OH...O(ITO) hydrogen bond
#   (c) bidentate           : 2 P-O-M bonds
#   (d) tridentate          : 3 P-O-M bonds
# Chemisorption proceeds by condensation with surface hydroxyls,
#   P-OH + HO-M  ->  P-O-M + H2O,
# eliminating up to two water molecules per PA. The two sampled anchor types:
#   "OH": a surface hydroxyl. The P-O takes the position of the hydroxyl O;
#         the hydroxyl (O+H) and the P-OH hydrogen are removed as water.
#   "M" : an under-coordinated metal ion. The P-O binds the open site, no
#         water is eliminated (the acid H stays on that O).
# =====================================================================

# ---------------- user settings ----------------
SLAB_FILE = "/content/CONTCAR"                  # ITO slab
MOLECULES = {                                   # label: path to optimized molecule
    "F5BzPA": "/content/CONTCAR_F5BzPA",
    "FBzPA":  "/content/CONTCAR_FBzPA",
}
OUT_ROOT    = "/content/POSCAR_PA_on_ITO"
METALS      = ("In", "Sn")
D_M_O       = {"In": 2.15, "Sn": 2.05}   # M-O distance for O placed on an open metal site (A)
D_PO        = 1.52    # P-O bond length
D_PC        = 1.82    # P-C bond length
D_HBOND     = 1.70    # target H...O distance for the P-OH...O(ITO) hydrogen bond in mode (b)
OH_CUTOFF   = 1.20    # O-H distance that defines a surface hydroxyl (A)
MO_CUTOFF   = 2.70    # M-O cutoff for coordination / open direction (A)
SURF_DEPTH  = 2.5     # metals within this depth of the highest metal are surface candidates (A)
MIN_OPEN    = 0.35    # minimum |sum of unit vectors to O neighbours| for an exposed cation
MIN_UP      = 0.30    # minimum z component of the open direction
CLEAR_HEAVY = 2.00    # open-site O must be this far from other slab heavy atoms (A)
CLEAR_H     = 1.50    # ... and this far from slab H atoms (A)
N_OH_ANCHORS = 12     # surface hydroxyls kept as anchors (highest first)
N_M_ANCHORS  = 8      # open metal sites kept as anchors (most exposed first)
N_MONO       = 4      # monodentate anchors per anchor type
N_HBOND      = 4      # mode (b) anchors per anchor type
N_PER_CLASS_BI  = 3   # bidentate pairs per water-count class (0, 1, 2 H2O)
N_PER_CLASS_TRI = 2   # tridentate triples per water-count class (0, 1, 2 H2O)
OO_MIN, OO_MAX = 2.3, 4.4   # allowed distance between anchor positions for multidentate modes (A)
VACUUM_TOP  = 12.0
FREEZE_BELOW = 3.0    # freeze slab atoms with z below this value (A); 0 = none
MIN_CONTACT = 1.8     # minimum distance of non-bound molecule atoms to slab after lifting (A)
MAX_LIFT    = 1.0
LIFT_STEP   = 0.05
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
heavy_idx = np.where(sym_s != "H")[0]
zmetal = pos[metal_idx, 2].max()
ztop = pos[heavy_idx, 2].max()
print(f"Slab: {slab.get_chemical_formula()}  highest metal z = {zmetal:.2f}  highest heavy atom z = {ztop:.2f}")

Dss = slab.get_all_distances(mic=True)

def mic_vec(a, b):
    return find_mic((b - a)[None, :], CELL, PBC)[0][0]

def unwrap_to(ref, p):
    return ref + mic_vec(ref, p)

# ---- anchor type "OH": surface hydroxyls ----
anchors_OH = []
hyd_O = {}
for h in H_idx_s:
    o = int(O_idx_s[np.argmin(Dss[h, O_idx_s])])
    if Dss[h, o] < OH_CUTOFF:
        hyd_O[o] = int(h)
for o, h in hyd_O.items():
    ms = [int(m) for m in metal_idx if Dss[o, m] < MO_CUTOFF]
    if ms:
        cen = np.mean([unwrap_to(pos[o], pos[m]) for m in ms], axis=0)
        d = pos[o] - cen
        d = d / np.linalg.norm(d) if np.linalg.norm(d) > 1e-6 else np.array([0, 0, 1.0])
    else:
        d = np.array([0, 0, 1.0])
    anchors_OH.append(dict(kind="OH", tag=f"OH{o}", pos=pos[o].copy(), dir=d, leave=[o, h], score=pos[o, 2], nmetal=len(ms)))
anchors_OH.sort(key=lambda a: -a["score"])
anchors_OH = anchors_OH[:N_OH_ANCHORS]

# ---- anchor type "M": exposed, unobstructed under-coordinated cations ----
anchors_M = []
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
    p = pos[m] + D_M_O[el] * d
    dq = np.linalg.norm(find_mic(pos - p, CELL, PBC)[0], axis=1)
    ok = all(dq[j] >= (CLEAR_H if sym_s[j] == "H" else CLEAR_HEAVY) for j in range(ns) if j != m)
    if ok:
        anchors_M.append(dict(kind="M", tag=f"{el}{int(m)}", pos=p, dir=d, leave=[], score=mag + 0.05 * pos[m, 2], el=el))
anchors_M.sort(key=lambda a: -a["score"])
anchors_M = anchors_M[:N_M_ANCHORS]

print(f"Surface hydroxyls found: {len(hyd_O)} (anchors kept: {len(anchors_OH)})")
print(f"Open metal sites kept: {len(anchors_M)}")
assert anchors_OH or anchors_M, "No anchors found; relax the site-finding settings."
pool = anchors_OH + anchors_M

# acceptor oxygens for the P-OH...O(ITO) hydrogen bond (surface oxygens near the top)
acceptors = [int(o) for o in O_idx_s if pos[o, 2] > ztop - 2.0]

def apos(a, ref):
    return unwrap_to(ref, a["pos"])

def adist(a, b):
    return np.linalg.norm(mic_vec(a["pos"], b["pos"]))

def leaves_overlap(*anc):
    s = []
    for a in anc:
        s += a["leave"]
    return len(s) != len(set(s))

def nwater(*anc):
    return sum(a["kind"] == "OH" for a in anc)

# ---- choose multidentate anchor sets, balanced over water-count classes ----
def best_per_class(combos, n_per, nanch):
    out = []
    for nw in range(0, 3):
        sel = sorted([c for c in combos if nwater(*c[1:]) == nw], key=lambda c: c[0])[:n_per]
        out += sel
    return out

pairs = []
for a, b in itertools.combinations(pool, 2):
    d = adist(a, b)
    if OO_MIN <= d <= OO_MAX and not leaves_overlap(a, b):
        pairs.append((abs(d - 2.9), a, b))
pairs = best_per_class(pairs, N_PER_CLASS_BI, 2)

triples = []
for a, b, c in itertools.combinations(pool, 3):
    if nwater(a, b, c) > 2 or all(x["kind"] == "OH" for x in (a, b, c)):
        continue            # at most two P-OH can condense; the P=O oxygen needs an open M site
    ds = [adist(a, b), adist(a, c), adist(b, c)]
    if all(OO_MIN <= d <= OO_MAX for d in ds) and not leaves_overlap(a, b, c):
        triples.append((sum(abs(d - 3.0) for d in ds), a, b, c))
triples = best_per_class(triples, N_PER_CLASS_TRI, 3)

mono_anch = anchors_OH[:N_MONO] + anchors_M[:N_MONO]
hb_anch = anchors_OH[:N_HBOND] + anchors_M[:N_HBOND]
print("Monodentate anchors :", [a["tag"] for a in mono_anch])
print("Bidentate pairs     :", [(a["tag"], b["tag"], nwater(a, b)) for _, a, b in pairs])
print("Tridentate triples  :", [(a["tag"], b["tag"], c["tag"], nwater(a, b, c)) for _, a, b, c in triples])

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
    T = np.array(O_targets)
    cxy = T[:, :2].mean(0)
    r = np.linalg.norm(T[:, :2] - cxy, axis=1).mean()
    dz = np.sqrt(max(D_PO**2 - r**2, 0.5**2))
    Pt = np.array([cxy[0], cxy[1], T[:, 2].mean() + dz])
    Ct = Pt + np.array([0, 0, D_PC])
    return Pt, Ct

def slab_xyz(leave):
    mask = np.ones(ns, bool)
    mask[list(leave)] = False
    return pos[mask]

def min_d(A, B):
    return get_distances(A, B, cell=CELL, pbc=PBC)[1].min()

# ================= per-molecule builder =================
def process(label, mol_file):
    out_dir = os.path.join(OUT_ROOT, label)
    os.makedirs(out_dir, exist_ok=True)
    mol = read(mol_file)
    mol.set_pbc(False)
    nm = len(mol)

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

    def fit(idx, targets, weights):
        P = mol.positions[idx]
        Q = np.array(targets)
        R, cp, cq = kabsch(P, Q, weights)
        new = (mol.positions - cp) @ R.T + cq
        rmsd = np.sqrt(np.average(((new[idx] - Q) ** 2).sum(1), weights=weights))
        return new, rmsd

    def clash(new_pos, exclude, sxyz):
        rows = [i for i in range(nm) if i not in exclude]
        return min_d(new_pos[rows], sxyz)

    def acid_H_removed(bound, anc_of):
        """P-OH hydrogens lost as water: only for O condensed onto a surface hydroxyl."""
        return {h for o in bound if anc_of[o]["kind"] == "OH" for h in H_of[o]}

    def leaving(anc_of):
        s = set()
        for a in anc_of.values():
            s |= set(a["leave"])
        return s

    def assign(O_list, anchors, weights, constrain_PO=False):
        """Best O -> anchor assignment by Kabsch rmsd. P=O oxygen must go to an open M site."""
        best = None
        ref = anchors[0]["pos"]
        T = [apos(a, ref) for a in anchors]
        Pt, Ct = P_C_targets(T)
        for perm in itertools.permutations(range(len(anchors))):
            if constrain_PO and any(O_list[p] in po_O and anchors[perm[p]]["kind"] != "M" for p in range(len(O_list))):
                continue
            tg = [T[perm[p]] for p in range(len(O_list))]
            new, rmsd = fit(list(O_list) + [iP, iC], tg + [Pt, Ct], list(weights) + [1.0, 1.0])
            if best is None or rmsd < best[1]:
                best = (new, rmsd, {O_list[p]: anchors[perm[p]] for p in range(len(O_list))})
        return best

    def mono_scan(anchor, tilts=(0.0,), psis=(0.0,), extra_score=None):
        """Place O_b on the anchor; scan P-O azimuth (and optional tilt) for the best geometry."""
        O_b = oh_O[0]
        sp = anchor["pos"]
        leave = set(anchor["leave"])
        sxyz = slab_xyz(leave)
        anc_of = {O_b: anchor}
        drop = acid_H_removed([O_b], anc_of)
        excl = {O_b} | drop
        best = None
        for phi in np.deg2rad(np.arange(0, 360, 30)):
            d_PO = np.array([np.sin(np.deg2rad(109.5)) * np.cos(phi),
                             np.sin(np.deg2rad(109.5)) * np.sin(phi),
                             np.cos(np.deg2rad(109.5))])
            Pt = sp - D_PO * d_PO
            for tilt in tilts:
                for psi in (psis if tilt > 0 else (0.0,)):
                    t, ps = np.deg2rad(tilt), np.deg2rad(psi)
                    cdir = np.array([np.sin(t) * np.cos(ps), np.sin(t) * np.sin(ps), np.cos(t)])
                    Ct = Pt + D_PC * cdir
                    new, rmsd = fit([O_b, iP, iC], [sp, Pt, Ct], [1, 1, 1])
                    c = clash(new, excl, sxyz)
                    sc = c if extra_score is None else extra_score(new, c, anc_of, leave)
                    if best is None or sc > best[0]:
                        best = (sc, new, rmsd)
        return best[1], [O_b], anc_of, best[2]

    # ---- mode (a): monodentate ----
    def mono(anchor):
        return mono_scan(anchor)

    # ---- mode (b): one P-O-M bond + P-OH...O(ITO) hydrogen bond ----
    def hbond(anchor):
        O_h = oh_O[1]
        H_i = H_of[O_h][0]
        def score(new, c, anc_of, leave):
            acc = [j for j in acceptors if j not in leave]
            dHA = min_d(new[[H_i]], pos[acc])
            pen = 0.0 if (c >= 1.5 and dHA >= 1.4) else 10.0
            return -abs(dHA - D_HBOND) - pen
        new, bound, anc_of, rmsd = mono_scan(anchor, tilts=(0.0, 25.0, 50.0), psis=np.arange(0, 360, 45), extra_score=score)
        return new, bound, anc_of, rmsd

    # ---- mode (c): bidentate ----
    def bidentate(a, b):
        new, rmsd, anc_of = assign(oh_O, [a, b], [1, 1])
        return new, list(oh_O), anc_of, rmsd

    # ---- mode (d): tridentate ----
    def tridentate(a, b, c):
        Ol = [oh_O[0], oh_O[1], po_O[0]]
        new, rmsd, anc_of = assign(Ol, [a, b, c], [1, 1, 1], constrain_PO=True)
        return new, Ol, anc_of, rmsd

    jobs = [(f"a-mono_{a['tag']}", (lambda a=a: mono(a))) for a in mono_anch]
    jobs += [(f"b-bi-Hbond_{a['tag']}", (lambda a=a: hbond(a))) for a in hb_anch]
    jobs += [(f"c-bi_{a['tag']}-{b['tag']}", (lambda a=a, b=b: bidentate(a, b))) for _, a, b in pairs]
    jobs += [(f"d-tri_{a['tag']}-{b['tag']}-{c['tag']}", (lambda a=a, b=b, c=c: tridentate(a, b, c))) for _, a, b, c in triples]

    summary = []
    for mode, builder in jobs:
        new_pos, bound, anc_of, rmsd = builder()
        leave = leaving(anc_of)
        drop = acid_H_removed(bound, anc_of)
        nw = int(sum(anc_of[o]["kind"] == "OH" for o in bound))
        exclude = set(bound) | drop
        sxyz = slab_xyz(leave)

        lift = 0.0
        while clash(new_pos + np.array([0, 0, lift]), exclude, sxyz) < MIN_CONTACT and lift < MAX_LIFT:
            lift += LIFT_STEP
        placed = new_pos + np.array([0, 0, lift])
        contact = clash(placed, exclude, sxyz)

        # H...O distance for mode (b)
        dHA = np.nan
        if mode.startswith("b-"):
            H_i = H_of[oh_O[1]][0]
            acc = [j for j in acceptors if j not in leave]
            dHA = min_d(placed[[H_i]], pos[acc])

        m = mol.copy()
        m.positions = placed
        if drop:
            m = m[[i for i in range(nm) if i not in drop]]
        s_r = slab[[i for i in range(ns) if i not in leave]]
        s2 = s_r.copy()
        s2 += m
        nsr = len(s_r)
        zmax = s2.positions[:, 2].max()
        c_new = max(slab.cell[2, 2], zmax + VACUUM_TOP)
        cell = s2.cell.array.copy(); cell[2] = [0, 0, c_new]
        s2.set_cell(cell, scale_atoms=False)
        fix_idx = [i for i in range(nsr) if FREEZE_BELOW > 0 and s2.positions[i, 2] < FREEZE_BELOW]
        s2.set_constraint(FixAtoms(indices=fix_idx))

        Dm = s2.get_all_distances(mic=True)
        Dno = s2.get_all_distances(mic=False)
        dmin_slab = Dm[nsr:, :nsr].min()
        sub, subn = Dm[nsr:, nsr:], Dno[nsr:, nsr:]
        img = sub[sub < subn - 0.01]
        dimg = img.min() if img.size else np.inf

        name = f"POSCAR_{label}_{mode}_{nw}H2O"
        write(os.path.join(out_dir, name), s2, format="vasp", direct=True, vasp5=True, sort=True)
        summary.append((name, nw, rmsd, lift, contact, dmin_slab, dimg, dHA))
    return summary

# ================= run all molecules =================
for label, path in MOLECULES.items():
    summ = process(label, path)
    print(f"{'file':62s} {'H2O':>3s} {'rmsd':>5s} {'lift':>5s} {'d(nonbd)':>8s} {'d(all)':>7s} {'d(img)':>7s} {'d(H..O)':>8s}")
    for n, nw, r, lf, c, d1, d2, dh in summ:
        bad = (c < MIN_CONTACT - 0.05 or d1 < 1.6 or d2 < 2.0 or (not np.isnan(dh) and dh > 2.3))
        dhs = "" if np.isnan(dh) else f"{dh:8.2f}"
        print(f"{n:62s} {nw:3d} {r:5.2f} {lf:5.2f} {c:8.2f} {d1:7.2f} {d2:7.2f} {dhs:>8s}{'  <-- check' if bad else ''}")
