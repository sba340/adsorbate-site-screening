"""
Generate one POSCAR per symmetry-inequivalent adsorption site type
(ontop, bridge, hollow_fcc, hollow_hcp) for a given slab + adsorbate pair.

Usage:
    python src/generate_sites.py \
        --slab CONTCAR_Co --ads CONTCAR_SH2 --anchor S \
        --out-root adsorption_sites

Replaces the Colab-specific cells in notebooks/Adsorbate_placement_on_slab.ipynb
(no hardcoded /content/ paths, no google.colab.files.download call).
"""
import argparse
import os

import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms

from site_finder import check_slab_symmetry, get_layer_indices, find_sites, dedupe
from placement import find_anchor_index, place_adsorbate


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--slab", required=True, help="Path to relaxed, bare slab CONTCAR")
    p.add_argument("--ads", required=True, help="Path to adsorbate CONTCAR/geometry file")
    p.add_argument("--anchor", required=True, help="Element symbol of the adsorbate's binding atom, e.g. S")
    p.add_argument("--out-root", default="adsorption_sites", help="Output directory for generated POSCARs")
    p.add_argument("--n-layers", type=int, default=5)
    p.add_argument("--n-fixed-layers", type=int, default=2)
    p.add_argument("--height", type=float, default=1.8, help="Anchor height above surface, Angstrom")
    p.add_argument("--nn-cutoff", type=float, default=3.0)
    p.add_argument("--hcp-thresh", type=float, default=0.75)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--skip-symmetry-check", action="store_true",
                    help="Skip the SpacegroupAnalyzer sanity check on the slab")
    return p.parse_args()


def main():
    args = parse_args()
    ads_label = os.path.basename(args.ads).split("_", 1)[-1] if "_" in os.path.basename(args.ads) else os.path.basename(args.ads)
    os.makedirs(args.out_root, exist_ok=True)

    slab = read(args.slab)
    adsorbate = read(args.ads)

    if not args.skip_symmetry_check:
        sg, n_ops = check_slab_symmetry(slab)
        print(f"Slab space group: {sg} ({n_ops} symmetry ops)")
        print("Review this before trusting one representative site per type.\n")

    anchor_idx = find_anchor_index(adsorbate, args.anchor)
    layers = get_layer_indices(slab, args.n_layers, args.seed)
    fixed = np.concatenate(layers[-args.n_fixed_layers:]).tolist()

    sites = dedupe(find_sites(slab, layers, args.nn_cutoff, args.hcp_thresh))
    by_type = {}
    for s in sites:
        by_type.setdefault(s["type"], []).append(s)

    for site_type, site_list in by_type.items():
        print(f"  {site_type}: {len(site_list)} symmetry-equivalent candidate(s) found")

    chosen_sites = {t: s[0] for t, s in by_type.items()}

    for site_type, site in chosen_sites.items():
        system = place_adsorbate(slab, adsorbate, anchor_idx, site["xy"], site["z"], args.height)
        system.set_constraint(FixAtoms(indices=fixed))
        folder = os.path.join(args.out_root, f"{ads_label}_{site_type}")
        os.makedirs(folder, exist_ok=True)
        write(os.path.join(folder, "POSCAR"), system, format="vasp", vasp5=True, direct=True)
        print(f"wrote {folder}/POSCAR")


if __name__ == "__main__":
    main()
