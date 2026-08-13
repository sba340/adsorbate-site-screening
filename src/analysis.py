"""
Post-DFT checks on relaxed adsorbate/slab structures: dissociation,
desorption, and final binding site classification by coordination number.
"""
import re
from ase.neighborlist import NeighborList, natural_cutoffs


def check_dissociation(initial_ads, final_ads):
    try:
        inl = NeighborList(natural_cutoffs(initial_ads), self_interaction=False, bothways=True)
        inl.update(initial_ads)
        fnl = NeighborList(natural_cutoffs(final_ads), self_interaction=False, bothways=True)
        fnl.update(final_ads)
        for i in range(len(initial_ads)):
            if set(inl.get_neighbors(i)[0]) != set(fnl.get_neighbors(i)[0]):
                return True
        return False
    except Exception:
        return True


def check_desorption(system, n_slab_atoms, cushion=1.5):
    try:
        cutoffs = [c * cushion for c in natural_cutoffs(system)]
        nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
        nl.update(system)
        for ads_idx in range(n_slab_atoms, len(system)):
            if any(n < n_slab_atoms for n in nl.get_neighbors(ads_idx)[0]):
                return False
        return True
    except Exception:
        return True


def parse_label_site(path, metal_symbol):
    fname = path.split("/")[-1]
    pattern = rf"CONTCAR_{metal_symbol}_(.+)_(ontop|bridge|hollow_fcc|hollow_hcp)$"
    m = re.match(pattern, fname)
    if not m:
        raise ValueError(f"Filename doesn't match expected pattern: {fname}")
    return m.group(1), m.group(2)


def analyze_final_structure(contcar_path, poscar_root, metal_symbol, anchor_symbol, coord_cutoff=3.2):
    from ase.io import read

    ads_label, site_label = parse_label_site(contcar_path, metal_symbol)
    final_system = read(contcar_path)

    poscar_path = f"{poscar_root}/{ads_label}_{site_label}/POSCAR"
    initial_system = read(poscar_path)
    n_slab = sum(1 for a in initial_system if a.symbol == metal_symbol)

    initial_ads = initial_system[n_slab:].copy()
    final_ads = final_system[n_slab:].copy()

    dissociated = check_dissociation(initial_ads, final_ads) if len(initial_ads) > 1 else False
    desorbed = check_desorption(final_system, n_slab)

    anchor_idx = [i for i, a in enumerate(final_system) if a.symbol == anchor_symbol][0]
    metal_idx = [i for i, a in enumerate(final_system) if a.symbol == metal_symbol]

    top_metal_z = final_system.positions[metal_idx, 2].max()
    height = final_system.positions[anchor_idx, 2] - top_metal_z
    dists = sorted(final_system.get_distances(anchor_idx, metal_idx, mic=True))
    n_coord = sum(d < coord_cutoff for d in dists)

    if n_coord == 1:
        final_site_type = "ontop"
    elif n_coord == 2:
        final_site_type = "bridge"
    elif n_coord >= 3:
        final_site_type = "hollow (fcc/hcp not distinguished here)"
    else:
        final_site_type = "desorbed / undercoordinated"

    return {
        "file": contcar_path.split("/")[-1],
        "intended_site": site_label,
        "final_site_by_coordination": final_site_type,
        "height_above_surface": round(height, 3),
        "anchor_metal_dist": round(dists[0], 3),
        "coordination_number": n_coord,
        "dissociated": dissociated,
        "desorbed": desorbed,
        "site_moved": site_label not in final_site_type,
    }
