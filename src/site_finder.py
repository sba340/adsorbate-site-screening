"""
Symmetry check and adsorption-site enumeration on a relaxed metal slab.

Site types identified: ontop, bridge, hollow_fcc, hollow_hcp
Method: KMeans layer separation + Delaunay triangulation of the top layer,
with 3x3 periodic tiling so edge sites near cell boundaries aren't missed.
"""
import numpy as np
from sklearn.cluster import KMeans
from scipy.spatial import Delaunay


def check_slab_symmetry(atoms, symprec=0.1):
    """Return (space_group_symbol, n_symmetry_ops) for a relaxed slab.

    A high-symmetry result (e.g. many ops, simple space group) means
    candidate sites of the same type are physically equivalent, so DFT
    only needs to be run on one representative site per type.
    """
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    pmg_slab = AseAtomsAdaptor.get_structure(atoms)
    sga = SpacegroupAnalyzer(pmg_slab, symprec=symprec)
    return sga.get_space_group_symbol(), len(sga.get_symmetry_operations())


def get_layer_indices(atoms, n_layers, seed=42):
    """Cluster atoms into n_layers by z-coordinate, ordered top -> bottom."""
    z = atoms.positions[:, 2].reshape(-1, 1)
    km = KMeans(n_clusters=n_layers, random_state=seed, n_init=10).fit(z)
    means = [atoms.positions[km.labels_ == i, 2].mean() for i in range(n_layers)]
    order = np.argsort(means)[::-1]
    return [np.where(km.labels_ == order[k])[0] for k in range(n_layers)]


def _periodic_xy(xy, cell, reps=(-1, 0, 1)):
    a, b = cell[0][:2], cell[1][:2]
    return np.vstack([xy + i * np.array(a) + j * np.array(b) for i in reps for j in reps])


def _in_central_cell(xy, cell, tol=1e-3):
    A = np.array([cell[0][:2], cell[1][:2]]).T
    frac = np.linalg.solve(A, xy)
    return np.all(frac >= -tol) and np.all(frac < 1 - tol)


def find_sites(atoms, layers, nn_cutoff=3.0, hcp_thresh=0.75):
    """Enumerate ontop/bridge/hollow_fcc/hollow_hcp candidate sites."""
    top_idx, second_idx = layers[0], layers[1]
    z_surf = atoms.positions[top_idx, 2].mean()
    top_xy = atoms.positions[top_idx, :2]
    second_xy = _periodic_xy(atoms.positions[second_idx, :2], atoms.get_cell())

    sites = [dict(type="ontop", xy=xy) for xy in top_xy]

    tiled = _periodic_xy(top_xy, atoms.get_cell())
    tri = Delaunay(tiled)

    seen = set()
    for simplex in tri.simplices:
        for i in range(3):
            a, b = simplex[i], simplex[(i + 1) % 3]
            key = tuple(sorted((a, b)))
            if key in seen:
                continue
            seen.add(key)
            p1, p2 = tiled[a], tiled[b]
            if np.linalg.norm(p1 - p2) > nn_cutoff:
                continue
            mid = (p1 + p2) / 2
            if _in_central_cell(mid, atoms.get_cell()):
                sites.append(dict(type="bridge", xy=mid))

    for simplex in tri.simplices:
        pts = tiled[simplex]
        edges = [np.linalg.norm(pts[i] - pts[(i + 1) % 3]) for i in range(3)]
        if max(edges) > nn_cutoff:
            continue
        c = pts.mean(axis=0)
        if not _in_central_cell(c, atoms.get_cell()):
            continue
        d = np.linalg.norm(second_xy - c, axis=1)
        t = "hollow_hcp" if d.min() < hcp_thresh else "hollow_fcc"
        sites.append(dict(type=t, xy=c))

    for s in sites:
        s["z"] = z_surf
    return sites


def dedupe(sites, tol=0.3):
    """Drop symmetry-equivalent duplicates of the same site type."""
    kept = []
    for s in sites:
        if not any(
            k["type"] == s["type"] and np.linalg.norm(np.array(k["xy"]) - np.array(s["xy"])) < tol
            for k in kept
        ):
            kept.append(s)
    return kept
