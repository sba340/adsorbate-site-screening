"""
Build a synthetic Co(111) slab + a simple H2S molecule and write them as
CONTCAR-format files, purely so the repo has runnable example inputs that
are NOT real, unpublished research structures.

Usage:
    python examples/make_example_structures.py
"""
from ase import Atoms
from ase.build import fcc111
from ase.io import write

if __name__ == "__main__":
    slab = fcc111("Co", size=(3, 3, 5), vacuum=12.0, a=2.51)
    write("examples/CONTCAR_Co_example", slab, format="vasp", vasp5=True, direct=True)

    # Simple bent H2S geometry (not built from the ASE g2 database, which
    # doesn't include H2S) with a 15 A vacuum cell so `read`/`write` as a
    # standalone VASP-format file behaves the same as the CONTCAR workflow.
    ads = Atoms(
        "SH2",
        positions=[
            (0.0, 0.0, 0.0),
            (0.962, 0.0, 0.932),
            (-0.962, 0.0, 0.932),
        ],
    )
    ads.center(vacuum=10.0)
    write("examples/CONTCAR_SH2_example", ads, format="vasp", vasp5=True, direct=True)

    print("Wrote examples/CONTCAR_Co_example and examples/CONTCAR_SH2_example")
