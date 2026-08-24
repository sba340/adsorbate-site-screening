# adsorbate-site-screening
Workflow for enumerating symmetry-inequivalent adsorption sites (ontop,
bridge, hollow_fcc, hollow_hcp) on a relaxed metal slab, placing an adsorbate
on each, and screening the resulting VASP CONTCARs for dissociation,
desorption, and site migration.

Built for bimetallic/transition-metal catalyst adsorbate screening
(e.g. S, SH, SH2, SCH3, SHCH3 on Co(111)), but the site-finding and
placement logic is general to any close-packed metal slab + small-molecule
adsorbate.

## How it works
1. **Symmetry check** (`src/site_finder.py::check_slab_symmetry`) — confirms
   the clean, relaxed slab is symmetric enough that candidate sites of a
   given type are physically equivalent, so DFT only needs to run on one
   representative site per type rather than every candidate found.
2. **Site enumeration** (`src/site_finder.py`) — KMeans-separates atomic
   layers by z, then Delaunay-triangulates the top layer (with 3x3 periodic
   tiling) to find ontop, bridge, and hollow sites, classifying hollows as
   fcc or hcp by proximity to the second layer.
3. **Placement** (`src/placement.py`) — orients the adsorbate by its anchor
   atom, places it above each chosen site, and nudges it upward if the
   initial placement clashes with the slab.
4. **Post-DFT analysis** (`src/analysis.py`) — after VASP relaxation, checks
   each CONTCAR for bond dissociation, desorption, and reclassifies the
   final site by anchor-atom coordination number.

### Final site classification
After DFT, the adsorbate's final site is reclassified independently of
its intended starting site, using the anchor atom's coordination number
to the metal slab (distance < `COORD_CUTOFF`):

- 1 metal neighbor → ontop
- 2 → bridge
- ≥3 → hollow, then split into fcc/hcp by the anchor's lateral proximity
  (minimum-image, periodic-boundary-aware) to the second metal layer

`site_moved` flags when the final classification differs from the
intended one, which is common when a shallow starting site (e.g. bridge)
relaxes into a lower-energy hollow.

`dissociated` is only meaningful for multi-atom adsorbates (SH, SH2,
SCH3, SHCH3); for monatomic adsorbates (S) it always returns `False`
by construction, not because dissociation was checked and ruled out.

## Usage
```bash
pip install -r requirements.txt
# 1. Generate one POSCAR per site type
python src/generate_sites.py \
    --slab CONTCAR_Co --ads CONTCAR_SH2 --anchor S \
    --out-root adsorption_sites
# 2. Run VASP relaxations on each adsorption_sites/*/POSCAR yourself (SLURM, etc.)
# 3. Classify the relaxed CONTCARs
python src/analyze_results.py \
    --pattern "CONTCAR_Co_*" --poscar-root adsorption_sites \
    --metal Co --anchor S
```

Try it on synthetic example structures first (no real research data
required):
```bash
python examples/make_example_structures.py
python src/generate_sites.py \
    --slab examples/CONTCAR_Co_example --ads examples/CONTCAR_SH2_example \
    --anchor S --out-root examples/adsorption_sites
```

## Notebook
`notebooks/Adsorbate_placement_on_slab.ipynb` is the original Colab
notebook this workflow was developed in (uses `/content/` paths and
`google.colab.files.download`). The `src/` scripts are the portable,
cluster/local-friendly version of the same logic and are what this repo
is meant to be used from going forward.

## Notes
- Site-finding parameters (`--nn-cutoff`, `--hcp-thresh`, `--n-layers`,
  `--n-fixed-layers`) were tuned for a Co(111) 5-layer slab; check they
  make sense for other metals/facets/slab thicknesses.
- `--anchor` must match exactly one atom of that element in the adsorbate
  file, since orientation and site placement are both anchored on it.
- `analyze_results.py` requires exactly one anchor atom in the final
  CONTCAR and raises an error otherwise. This is intentional: more than
  one matching atom usually means a dissociation fragment or a second
  adsorbate landed in the cell, and silently picking one would give a
  misleading result.
