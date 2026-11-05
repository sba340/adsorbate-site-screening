# Notebooks

| Notebook | Purpose |
|---|---|
| `Adsorbate_placement_on_slab.ipynb` | Original Colab version of the site-enumeration and placement workflow (S, SH, SH2, SCH3, SHCH3 on Co(111)). The portable version lives in `src/`. |
| `ITO_Au.ipynb` | Placement of phosphonic acid and thiol molecules on Au(111) and on hydroxylated ITO(222). Four workflows, each writing one POSCAR per binding mode, site and H treatment. |

## ITO_Au.ipynb

Builds starting geometries for VASP relaxations of F5BzPA and FBzPA (phosphonic acids) and PFBT and FBzS (thiols). The notebook contains four workflows:

| Part | Molecules | Surface |
|---|---|---|
| 1 | F5BzPA, FBzPA | Au(111) |
| 2 | PFBT, FBzS | Au(111) |
| 3 | F5BzPA, FBzPA | Hydroxylated ITO(222) |
| 4 | PFBT, FBzS | Hydroxylated ITO(222) |

The same code is available as standalone scripts in `src/molecule_placement/`. The ITO binding modes follow Paramonov et al., Chem. Mater. 2008, 20, 5131 (Scheme 1).

### Inputs

| File | Used in | Description |
|---|---|---|
| `CONTCAR_Au` | Parts 1, 2 | Relaxed Au(111) slab |
| `CONTCAR` | Parts 3, 4 | Hydroxylated ITO(222) slab |
| `CONTCAR_F5BzPA`, `CONTCAR_FBzPA` | Parts 1, 3 | Relaxed gas-phase phosphonic acids (one P, three O: two P-OH, one P=O) |
| `CONTCAR_PFBT`, `CONTCAR_FBzS` | Parts 2, 4 | Relaxed gas-phase thiols (exactly one S bonded to one C) |

Molecules must be unwrapped (not split across periodic boundaries). Each part asserts the expected P/O/H or S/C connectivity and stops with a message if it does not match.

### Part 1: phosphonic acids on Au(111)

Binding modes: `mono` (atop, bridge, fcc, hcp), `bi_atop-atop`, `bi-electro_atop`, `tri` (fcc, hcp). Each mode is written with the bound O-H hydrogen removed (`Hremove`) or kept (`Hkeep`).

Placement fits the anchoring O, P and C atoms to target positions with a weighted Kabsch fit. The molecule is then lifted in `LIFT_STEP` increments (up to `MAX_LIFT`) until non-bound atoms are at least `MIN_CONTACT` from Au.

| Setting | Default | Meaning |
|---|---|---|
| `D_AU_O` | 2.20 A | Target Au-O distance for bound O |
| `D_PO`, `D_PC` | 1.52, 1.82 A | P-O and P-C bond lengths used to build P and C targets |
| `H_INTERACT` | 2.90 A | Height of the non-bonded O above the top Au plane in bidentate-electrostatic mode |
| `MIN_CONTACT` | 1.9 A | Minimum non-bound atom to Au distance after lifting |
| `MAX_LIFT`, `LIFT_STEP` | 1.0, 0.05 A | Upward shift limit and increment |

### Part 2: thiols on Au(111)

Sites: atop, bridge, fcc, hcp, and bridge shifted toward each hollow. For every site, the S-C axis is tilted (`TILTS`), rotated in azimuth (`PHIS`), and spun about the S-C bond (`PSIS`). The orientation with the largest minimum distance to the slab is kept. `Hremove` gives the thiolate (S-Au), `Hkeep` the intact S-H thiol.

| Setting | Default | Meaning |
|---|---|---|
| `D_AU_S` | 2.40 A | Target Au-S distance |
| `TILTS` | 0, 30 deg | Tilt of S-C from the surface normal |
| `PHIS` | 0 to 300, step 60 deg | Tilt azimuth (used when tilt > 0) |
| `PSIS` | 0 to 330, step 30 deg | Spin about S-C |

### Part 3: phosphonic acids on ITO(222)

Binding modes follow Paramonov et al.:

| Mode | Label | Description |
|---|---|---|
| a | `a-mono` | Monodentate, one P-O-M bond |
| b | `b-bi-Hbond` | One P-O-M bond plus a P-OH...O(ITO) hydrogen bond |
| c | `c-bi` | Bidentate, two P-O-M bonds |
| d | `d-tri` | Tridentate, three P-O-M bonds |

Two anchor types are sampled. For `OH` anchors, the P-O takes the position of a surface hydroxyl, and the hydroxyl and the P-OH hydrogen leave as water (P-OH + HO-M gives P-O-M + H2O). For `M` anchors, the P-O binds an under-coordinated In or Sn site and no water is eliminated. The number of water molecules removed is written into each file name (`0H2O`, `1H2O`, `2H2O`).

| Setting | Default | Meaning |
|---|---|---|
| `METALS` | In, Sn | Metals considered as surface cation sites |
| `D_M_O` | In 2.15, Sn 2.05 A | M-O distance for O placed on an open metal site |
| `D_PO`, `D_PC` | 1.52, 1.82 A | P-O and P-C bond lengths |
| `D_HBOND` | 1.70 A | Target H...O distance for the hydrogen bond in mode b |
| `OH_CUTOFF` | 1.20 A | O-H distance that defines a surface hydroxyl |
| `MO_CUTOFF` | 2.70 A | M-O cutoff for coordination and open direction |
| `SURF_DEPTH` | 2.5 A | Metals within this depth of the highest metal are surface candidates |
| `N_OH_ANCHORS`, `N_M_ANCHORS` | 12, 8 | Hydroxyl and open-metal anchors kept |
| `N_MONO`, `N_HBOND` | 4, 4 | Anchors per type for modes a and b |
| `OO_MIN`, `OO_MAX` | 2.3, 4.4 A | Allowed spacing between anchors in multidentate modes |
| `MIN_CONTACT` | 1.8 A | Minimum non-bound atom to slab distance after lifting |

### Part 4: thiols on ITO(222)

Uses the same slab model and anchor search as Part 3, with the orientation scan from Part 2 (`TILTS`, `PHIS`, `PSIS`). Anchor types:

| Anchor | Description | H treatment |
|---|---|---|
| `M` | S placed `D_M_S` from an exposed In or Sn cation along its open direction | `thiolate` (S-H hydrogen removed) or `thiol` (intact S-H) |
| `OH` | Condensation R-SH + HO-M gives R-S-M + H2O. S replaces the hydroxyl O, and the hydroxyl and S-H hydrogen leave as water | `H2O` |

Set `USE_OH_CONDENSATION = False` to skip the hydroxyl anchors.

| Setting | Default | Meaning |
|---|---|---|
| `D_M_S` | In 2.55, Sn 2.45 A | Target M-S distance |
| `CLEAR_HEAVY`, `CLEAR_H` | 2.50, 2.00 A | Minimum distance of an open-site S to slab heavy atoms and H atoms |
| `N_OH_ANCHORS`, `N_M_ANCHORS` | 8, 8 | Hydroxyl and open-metal anchors used |

### Shared settings

`VACUUM_TOP` (12 A) sets vacuum above the highest atom. On Au(111), `N_FREEZE` (2) fixes the bottom Au layers with `FixAtoms`; set to 0 to relax everything. On ITO, `FREEZE_BELOW` (3.0 A) fixes slab atoms with z below that value; set to 0 to relax everything. On Au(111), `H_POLICIES` controls which H treatments are generated.

### Output

```
POSCAR_PA_on_Au111/<molecule>/POSCAR_<molecule>_<mode>_H<policy>
POSCAR_thiol_on_Au111/<molecule>/POSCAR_<molecule>_<site>_tilt<deg>_H<policy>
POSCAR_PA_on_ITO/<molecule>/POSCAR_<molecule>_<mode>_<n>H2O
POSCAR_thiol_on_ITO/<molecule>/POSCAR_<molecule>_<anchor>_tilt<deg>_<policy>
```

Each run prints a summary table. Rows marked `<-- check` have a contact or image distance below the acceptable limit (and, for mode b on ITO, an H...O distance above 2.3 A) and should be inspected before submitting.

### Running

Colab: the first cell installs ASE. Upload the CONTCARs to `/content/` and run the cells of the part you need.

Local or cluster: `pip install ase numpy`, edit `SLAB_FILE`, `MOLECULES` and `OUT_ROOT` in each settings block, and skip the `google.colab` zip/download cells.

### Notes

- The Au(111) parts were developed with a slab of 16 atoms per layer and a nearest-neighbor distance of 2.88 A. Check the layer detection and site geometry for other slabs.
- The ITO parts search the slab for surface hydroxyls and under-coordinated In or Sn sites. If none are found, the notebook stops with a message to relax the site-finding settings.
- Generated structures are starting guesses, not relaxed geometries. Run VASP relaxations before drawing conclusions about preferred modes.
- Cleanup cells at the end of each part delete the output folder and zip. Run them only after downloading.
