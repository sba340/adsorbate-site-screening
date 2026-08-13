"""
Classify each relaxed (post-DFT) CONTCAR by final adsorption site,
flagging dissociation or desorption vs. the intended site.

Usage:
    python src/analyze_results.py \
        --pattern "CONTCAR_Co_*" --poscar-root adsorption_sites \
        --metal Co --anchor S
"""
import argparse
import glob

from analysis import analyze_final_structure


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pattern", required=True, help="Glob pattern for relaxed CONTCARs, e.g. 'CONTCAR_Co_*'")
    p.add_argument("--poscar-root", required=True, help="Directory of pre-DFT POSCARs from generate_sites.py")
    p.add_argument("--metal", required=True, help="Slab element symbol, e.g. Co")
    p.add_argument("--anchor", required=True, help="Adsorbate anchor element symbol, e.g. S")
    p.add_argument("--coord-cutoff", type=float, default=3.2)
    return p.parse_args()


def main():
    args = parse_args()
    results = []
    for path in sorted(glob.glob(args.pattern)):
        try:
            results.append(
                analyze_final_structure(path, args.poscar_root, args.metal, args.anchor, args.coord_cutoff)
            )
        except Exception as e:
            print(f"Skipping {path}: {e}")

    for r in results:
        print(r)


if __name__ == "__main__":
    main()
