import itertools, os
import numpy as np
from ase.io import read, write
from ase.geometry import find_mic, get_distances
from ase.constraints import FixAtoms

# ---------------- user settings ----------------
SLAB_FILE = "/content/CONTCAR_Au"
MOLECULES = {
    "PFBT": "/content/CONTCAR_PFBT",   # adjust to your actual file names
    "FBzS": "/content/CONTCAR_FBzS",
}
OUT_ROOT   = "/content/POSCAR_thiol_on_Au111"
D_AU_S     = 2.40    # target Au-S distance (Angstrom)
TILTS      = (0.0, 30.0)                     # tilt of the S->C bond from the surface normal (degrees)
PHIS       = tuple(np.arange(0, 360, 60))    # azimuth of the tilt (degrees), used when tilt > 0
PSIS       = tuple(np.arange(0, 360, 30))    # spin about the S->C bond (degrees)
VACUUM_TOP = 12.0
N_FREEZE   = 2
H_POLICIES = ("remove", "keep")  # remove = thiolate (S-Au), keep = intact S-H thiol
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

def z_for_xy(p):
    """Height of S above the top plane so that the nearest top-layer Au is D_AU_S away."""
    d = pos[top] - p
    d[:, 2] = 0.0
    r = find_mic(d, slab.cell, slab.pbc)[1].min()
    return ztop + np.sqrt(max(D_AU_S**2 - r**2, 0.5**2))

def with_z(p):
    q = np.array(p, float); q[2] = z_for_xy(q); return q

tri = {}
for j, k in itertools.combinations(nbrs, 2):
    if slab.get_distance(j, k, mic=True) < 1.1 * nn:
        cen = pos[anchor] + (vec[j] + vec[k]) / 3
        d = pos[second] - cen
        d[:, 2] = 0.0
        dmin = find_mic(d, slab.cell, slab.pbc)[1].min()
        kind = "hcp" if dmin < 0.3 else "fcc"
        tri.setdefault(kind, (j, k, cen))

sites = {"atop": with_z(pos[anchor]),
         "bridge": with_z(pos[anchor] + vec[nbrs[0]] / 2)}
for kind, (j, k, cen) in tri.items():
    sites[kind] = with_z(cen)
    mid = pos[anchor] + vec[j] / 2
    sites[f"bridge-{kind}"] = with_z(mid + 0.5 * (cen - mid))   # bridge shifted toward that hollow

# layers to freeze
zs = np.sort(np.unique(np.round(pos[:, 2], 1)))
layer_z = []
for z in zs:
    if not layer_z or z - layer_z[-1] > 0.5:
        layer_z.append(z)
fix_idx = [i for i in range(len(slab)) if N_FREEZE > 0 and pos[i, 2] < layer_z[N_FREEZE - 1] + 0.5]

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

    policies = H_POLICIES if H_S else ("keep",)
    u0 = mol.positions[iC] - mol.positions[iS]
    summary = []

    for site, spos in sites.items():
        for tilt in TILTS:
            for policy in policies:
                excl = {iS} | (set(H_S) if policy == "remove" else set())
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
                        D = get_distances(newpos[rows], pos, cell=slab.cell, pbc=slab.pbc)[1]
                        c = D.min()
                        if best is None or c > best[0]:
                            best = (c, newpos)
                m = mol.copy()
                m.positions = best[1]
                if policy == "remove":
                    m = m[[i for i in range(len(m)) if i not in set(H_S)]]
                s2 = slab.copy()
                s2 += m
                zmax = s2.positions[:, 2].max()
                c_new = max(slab.cell[2, 2], zmax + VACUUM_TOP)
                cell = s2.cell.array.copy(); cell[2] = [0, 0, c_new]
                s2.set_cell(cell, scale_atoms=False)
                s2.set_constraint(FixAtoms(indices=fix_idx))
                ns = len(slab)
                Dm = s2.get_all_distances(mic=True)
                Dno = s2.get_all_distances(mic=False)
                sub, subn = Dm[ns:, ns:], Dno[ns:, ns:]
                img = sub[sub < subn - 0.01]
                dimg = img.min() if img.size else np.inf
                dAuS = get_distances([spos], pos, cell=slab.cell, pbc=slab.pbc)[1].min()
                name = f"POSCAR_{label}_{site}_tilt{int(tilt)}_H{policy}"
                write(os.path.join(out_dir, name), s2, format="vasp", direct=True, vasp5=True, sort=True)
                summary.append((name, dAuS, best[0], dimg))
    return summary

for label, path in MOLECULES.items():
    summ = process(label, path)
    print(f"{'file':44s} {'d(Au-S)':>8s} {'min d(other-Au)':>16s} {'min d(image)':>13s}")
    for n, d0, d1, d2 in summ:
        flag = "  <-- check" if (d1 < 2.0 or d2 < 2.0) else ""
        print(f"{n:44s} {d0:8.2f} {d1:16.2f} {d2:13.2f}{flag}")
