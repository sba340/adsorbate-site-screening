# Notebooks

| Notebook | Purpose |
|---|---|
| `Adsorbate_placement_on_slab.ipynb` | Original Colab version of the site-enumeration and placement workflow (S, SH, SH2, SCH3, SHCH3 on Co(111)). The portable version lives in `src/`. |
| `SAM_placement_on_Au111.ipynb` | Anchoring-group placement for phosphonic acid and thiol molecules on Au(111). Writes one POSCAR per binding mode, site and H policy. |

## SAM_placement_on_Au111.ipynb

Builds starting geometries for VASP relaxations of F5BzPA, FBzPA (phosphonic acids) and PFBT, FBzS (thiols) on a relaxed Au(111) slab.

### Inputs

| File | Description |
|---|---|
| `CONTCAR_Au` | Relaxed Au(111) slab |
| `CONTCAR_F5BzPA`, `CONTCAR_FBzPA` | Relaxed gas-phase phosphonic acids (one P, three O: two P-OH, one P=O) |
| `CONTCAR_PFBT`, `CONTCAR_FBzS` | Relaxed gas-phase thiols (exactly one S bonded to one C) |

Molecules must be unwrapped (not split across periodic boundaries). The notebook asserts the expected P/O/H and S/C connectivity and stops with a message if it does not match.

### Part 1: phosphonic acids

Binding modes: `mono` (atop, bridge, fcc, hcp), `bi_atop-atop`, `bi-electro_atop`, `tri` (fcc, hcp). Each mode is written with the bound O-H hydrogen removed (`Hremove`) or kept (`Hkeep`).

Placement fits the anchoring O, P and C atoms to target positions with a weighted Kabsch fit. The molecule is then lifted in `LIFT_STEP` increments (up to `MAX_LIFT`) until non-bound atoms are at least `MIN_CONTACT` from Au.

| Setting | Default | Meaning |
|---|---|---|
| `D_AU_O` | 2.20 A | Target Au-O distance for bound O |
| `D_PO`, `D_PC` | 1.52, 1.82 A | P-O and P-C bond lengths used to build P and C targets |
| `H_INTERACT` | 2.90 A | Height of the non-bonded O above the top Au plane in bidentate-electrostatic mode |
| `MIN_CONTACT` | 1.9 A | Minimum non-bound atom to Au distance after lifting |
| `MAX_LIFT`, `LIFT_STEP` | 1.0, 0.05 A | Upward shift limit and increment |

### Part 2: thiols

Sites: atop, bridge, fcc, hcp, and bridge shifted toward each hollow. For every site, the S-C axis is tilted (`TILTS`), rotated in azimuth (`PHIS`), and spun about the S-C bond (`PSIS`). The orientation with the largest minimum distance to the slab is kept. `Hremove` gives the thiolate (S-Au), `Hkeep` the intact S-H thiol.

| Setting | Default | Meaning |
|---|---|---|
| `D_AU_S` | 2.40 A | Target Au-S distance |
| `TILTS` | 0, 30 deg | Tilt of S-C from the surface normal |
| `PHIS` | 0 to 300, step 60 deg | Tilt azimuth (used when tilt > 0) |
| `PSIS` | 0 to 330, step 30 deg | Spin about S-C |

### Shared settings

`VACUUM_TOP` (12 A) sets vacuum above the highest atom. `N_FREEZE` (2) fixes the bottom Au layers with `FixAtoms`; set to 0 to relax everything. `H_POLICIES` controls which H treatments are generated.

### Output

```
POSCAR_PA_on_Au111/<molecule>/POSCAR_<molecule>_<mode>_H<policy>
POSCAR_thiol_on_Au111/<molecule>/POSCAR_<molecule>_<site>_tilt<deg>_H<policy>
```

Each run prints a summary table. Rows marked `<-- check` have a contact below `MIN_CONTACT` or an image distance below 2.0 A, and should be inspected before submitting.

### Running

Colab: upload the CONTCARs to `/content/` and run all cells.

Local or cluster: `pip install ase numpy`, edit `SLAB_FILE`, `MOLECULES` and `OUT_ROOT` in each settings block, and skip the `google.colab` zip/download cells.

### Notes

- Developed with an Au(111) slab with 16 atoms per layer and a nearest-neighbor distance of 2.88 A. Check the layer detection and site geometry for other slabs.
- Generated structures are starting guesses, not relaxed geometries. Run VASP relaxations before drawing conclusions about preferred modes.
- Cleanup cells at the end of each part delete the output folder and zip. Run them only after downloading.
