# Notebooks

| Notebook | Purpose |
|---|---|
| `Adsorbate_placement_on_slab.ipynb` | Placement of small adsorbates (S, SH, SH2, SCH3, SHCH3 on Co(111)) at on-top, bridge, hcp and fcc sites, with a slab symmetry check before placement and a site analysis after DFT. The portable version lives in `src/`. |
| `ITO_Au.ipynb` | Placement of phosphonic acid and thiol molecules on Au(111) and on hydroxylated ITO(222). Four workflows, each writing one POSCAR per binding mode, site and H treatment. |
## Adsorbate_placement_on_slab.ipynb

Builds starting geometries for VASP relaxations of one or more adsorbates on a clean metal slab, and analyzes the relaxed structures afterward. The defaults are set up for a five-layer Co(111) slab, but the site finder works for other HCP or FCC close-packed surfaces.

The notebook has three stages:

1. **Slab symmetry check.** Reads the relaxed bare slab and reports its space group and the number of symmetry operations (pymatgen, `symprec` 0.1 A). A high-symmetry result supports running VASP on one representative site per type. A low-symmetry result (for example P1) means the sites are not equivalent and more should be kept.
2. **Adsorbate placement.** Writes one POSCAR per adsorbate and site type.
3. **Post-DFT site analysis.** Compares a relaxed CONTCAR with its starting POSCAR and reports where the adsorbate ended up.

### Inputs

| File | Description |
|---|---|
| `CONTCAR_Co` | Relaxed bare slab |
| `CONTCAR_<adsorbate>` | Relaxed gas-phase adsorbate, one file per species (for example `CONTCAR_SH2`). The label is taken from the file name. |
| `POSCAR`, `CONTCAR` | Starting and relaxed structures of one adsorbate on the slab, used by the post-DFT analysis |

### Placement method

Sites are found from a Delaunay triangulation of the top-layer atoms, including periodic images. Layers are identified by clustering atoms on their z coordinate (KMeans). The sites are on-top, bridge, hollow_hcp and hollow_fcc. A hollow is labeled hcp if a second-layer atom lies below its centroid and fcc otherwise. One representative site per type is written.

The binding atom is chosen automatically from `ANCHOR_PRIORITY` (S, O, N, C, H, first match wins), unless `ANCHOR_OVERRIDE` names a different element for that adsorbate. If several atoms of that element exist, the one with the lowest z in the input file is used. The molecule is rotated so that the anchor faces the surface and the rest of the molecule points away from it. This includes the case where the CONTCAR stores the tail pointing straight down (SH2), which is handled with an explicit 180 degree flip. The anchor is then placed at `HEIGHT` above the top layer, and the adsorbate is moved up in 0.1 A steps until no adsorbate to slab distance is below 1.0 A.

| Setting | Default | Meaning |
|---|---|---|
| `SLAB_FILE` | `/content/CONTCAR_Co` | Clean slab |
| `ADS_FILES` | S, SH, SH2, SHCH3, SCH3 | List of adsorbate CONTCARs |
| `ANCHOR_PRIORITY` | S, O, N, C, H | Order in which the binding atom is chosen |
| `ANCHOR_OVERRIDE` | empty | Per-adsorbate binding element, for example `"CO": "C"` |
| `N_LAYERS` | 5 | Total layers in the slab |
| `N_FIXED_LAYERS` | 2 | Bottom layers fixed with `FixAtoms` |
| `HEIGHT` | 1.8 A | Height of the anchor above the mean z of the top layer |
| `NN_CUTOFF` | 3.0 A | Maximum bond length for bridge and hollow detection |
| `HCP_THRESH` | 0.75 A | Distance threshold for hcp versus fcc hollows |
| `SEED` | 42 | KMeans random seed |
| `OUT_ROOT` | `/content/adsorption_sites` | Output folder |

### Output

```
adsorption_sites/<adsorbate>/<site>/POSCAR
```

where `<site>` is `ontop`, `bridge`, `hollow_hcp` or `hollow_fcc`. The output folder is zipped and downloaded at the end of the run when running in Colab.

### Post-DFT site analysis

Set `ads_label` and `site_label` to the adsorbate and the intended site, then run the last cell on a `POSCAR` and `CONTCAR` pair. It prints one dictionary with:

| Key | Meaning |
|---|---|
| `final_site_by_coordination` | Site type from the number of metal atoms within `COORD_CUTOFF` of the anchor (1 ontop, 2 bridge, 3 or more hollow, classified as hcp or fcc) |
| `site_moved` | True if the final site differs from the intended one (hcp and fcc both count as hollow) |
| `height_above_surface` | Anchor height above the top metal layer |
| `anchor_metal_dist` | Shortest anchor to metal distance |
| `coordination_number` | Number of metal neighbors of the anchor |
| `dissociated` | True if the bonding within the adsorbate changed between POSCAR and CONTCAR |
| `desorbed` | True if no adsorbate atom has a slab neighbor within 1.5 times the natural bond cutoffs |

Settings: `METAL_SYMBOL` (Co), `ANCHOR_SYMBOL` (S), `COORD_CUTOFF` (3.2 A) and `LAYER_TOL` (0.5 A). The metal atoms are assumed to come first in the file, followed by the adsorbate atoms.

### Running

Colab: the first cells install `ase`, `scikit-learn`, `scipy`, `numpy` and `pymatgen`. Upload the CONTCARs to `/content/` and run the cells in order. Run the cleaning cell before each new placement run. It deletes `/content/adsorption_sites` and everything in it.

Local or cluster: `pip install ase scikit-learn scipy numpy pymatgen`, edit the paths in the settings block, and skip the Colab download step (the notebook skips it automatically outside Colab).

### Notes

- Only one site per type is written. This is appropriate when the symmetry check shows that all sites of a type are equivalent.
- Generated structures are starting guesses, not relaxed geometries.
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
