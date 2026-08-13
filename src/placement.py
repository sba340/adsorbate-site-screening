"""
Orient an adsorbate molecule by its anchor atom and place it on a slab site,
nudging upward until no unphysical clash with the slab remains.
"""
import numpy as np


def find_anchor_index(atoms, symbol):
    idx = [i for i, a in enumerate(atoms) if a.symbol == symbol]
    if len(idx) != 1:
        raise ValueError(f"Expected exactly 1 {symbol} atom, found {len(idx)}")
    return idx[0]


def place_adsorbate(slab, ads, anchor_idx, xy, z, height):
    """Return slab+adsorbate combined Atoms with the adsorbate's anchor atom
    positioned at (xy, z + height), tail atoms pointed away from the surface,
    and a small vertical clash-avoidance nudge applied if needed.
    """
    mol = ads.copy()
    n = len(mol)
    if n > 1:
        anchor_pos = mol.positions[anchor_idx].copy()
        others = [i for i in range(n) if i != anchor_idx]
        tail = mol.positions[others].mean(axis=0) - anchor_pos
        tn = np.linalg.norm(tail)
        if tn > 1e-6:
            tail /= tn
            target = np.array([0.0, 0.0, 1.0])
            axis = np.cross(tail, target)
            an = np.linalg.norm(axis)
            if an > 1e-6:
                axis /= an
                angle = np.degrees(np.arccos(np.clip(np.dot(tail, target), -1, 1)))
                mol.translate(-anchor_pos)
                mol.rotate(angle, axis)
                mol.translate(anchor_pos)

    mol.translate(np.array([xy[0], xy[1], z + height]) - mol.positions[anchor_idx])
    combined = slab.copy() + mol
    n_slab = len(slab)

    guard = 0
    while guard < 50:
        d = combined.get_all_distances()
        if np.min(d[n_slab:, :n_slab]) >= 1.0:
            break
        combined.positions[n_slab:, 2] += 0.1
        guard += 1
    return combined
