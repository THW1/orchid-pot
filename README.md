# Voronoi orchid basket generator

`orchid_pot.py` generates a printable STL of an orchid basket or pot (Vanda style), with
planar Voronoi walls and floor, a solid bottom band, a top rim and optional hanging
eyelets. The shape is chosen as a subcommand:

| Shape | Description | Eyelets (default) |
|---|---|---|
| `box` | rounded-square box | 4, at the corners |
| `cylinder` | cylinder, or with `--taper` a cone frustum (flowerpot) | 3 |
| `orb` | **experimental**: "onion" (stretched ellipsoid) with flat floor and top opening | 3 |

It is tuned for PETG with a 0.4 mm nozzle, printed upright without supports.

Design rationale and rejected alternatives are in the context header of `orchid_pot.py`.

## Requirements

| Component | Needed for | Tested with |
|---|---|---|
| Python 3 | everything | 3.14.8 |
| Python packages (`requirements.txt`) | generating STLs | see below |
| Blender | preview renders only (`render.py`), optional | 5.2.2 LTS |
| pytest, pyflakes (`requirements-dev.txt`) | tests and lint, optional | 8.4.2, 3.4.0 |

Install the Python packages:

```bash
pip3 install -r requirements.txt          # numpy scipy scikit-image trimesh fast-simplification
                                          # shapely rtree networkx manifold3d
pip3 install -r requirements-dev.txt      # adds pytest and pyflakes
```

`rtree` and `networkx` are not imported directly. trimesh needs them for the slice
polygons that the island and cantilever checks use. `manifold3d` is trimesh's boolean
engine for the tray fit check. The lower bounds in
`requirements.txt` are the tested versions; older versions are untested.

Blender comes from the distribution (for example `dnf install blender`) or from
blender.org. `render.py` uses Blender's own Python and the built-in STL importer, so no
pip packages are needed for it.

A full run takes about 30 s and peaks at about 1.3 GB RAM (default settings, 6 cores).

## Usage

```bash
python3 orchid_pot.py box                              # → orchid_box.stl + orchid_box_report.json
python3 orchid_pot.py cylinder --taper 8               # flowerpot, 150 mm at the rim
python3 orchid_pot.py cylinder --taper 12 --eyelets 0  # standing pot, no eyelets
python3 orchid_pot.py orb                              # onion, 150 mm belly (~95 s, ~2.4 GB)
python3 orchid_pot.py box --seed 3 --out pot_s3        # another random pattern, other file names
python3 orchid_pot.py cylinder --graph-only --seed 3   # pattern metrics only (~2 s), no mesh
python3 orchid_pot.py cylinder --parts pot,tray        # pot + drip tray (saucer)
python3 orchid_pot.py box --parts tray                 # only the tray (~5 s), pot unchanged
python3 orchid_pot.py orb --parts tray --tray clip     # clip-on tray for hanging pots
python3 orchid_pot.py box --from-report pot_s3_report.json --parts tray --tray clip  # tray for that pot
python3 orchid_pot.py cylinder --taper 8 --parts pot,planter --eyelets 0 --planter-eyelets 3  # hanging planter
python3 orchid_pot.py cylinder --help                  # options for that shape
blender -b --factory-startup --python render.py -- $PWD/orchid_box.stl $PWD/preview_box
```

`render.py` needs **absolute** paths. It writes `<prefix>_{persp,side,top,detail}.png`.

**Workflow:** compare seeds and parameters with `--graph-only` first, then build the mesh,
then check the report and previews.

## Options

Lengths are in mm. Recommended ranges are based on the runs in this project. Values
outside them work, but check the report.

### Shape options

| Shape | Option | Default | Recommended | Meaning |
|---|---|---|---|---|
| `box` | `--size` | 150 | – | Outer width and depth |
| `box` | `--corner-r` | 20 | 15–25 | Outer radius of the vertical edges; must exceed `--wall`/2 |
| `cylinder` | `--diameter` | 150 | – | Outer diameter at the rim |
| `cylinder` | `--taper` | 0 | −15…+15 | Wall lean from vertical in degrees; > 0 widens upward. Limit ±20 |
| `orb` | `--diameter` | 150 | – | Outer diameter at the belly |
| `orb` | `--stretch` | 1.2 | 1.0–1.4 | Height/width ratio of the ellipsoid; 1 = sphere. Allowed 0.6–2.0 |
| `orb` | `--min-wall-angle` | 55 | 50–60 | Wall angle from horizontal where floor and opening are cut. Allowed 45–80 |

At `--taper 8` the floor diameter is 121.9 mm; at `--taper 15` it is 96.4 mm (for a
150 mm rim, 100 mm height).

The orb has no `--height` or `--zstretch`: floor, opening and height follow from
`--diameter`, `--stretch` and `--min-wall-angle`, and the stretch elongates the cells.
Values at 150 mm diameter (mid-wall floor/opening Ø, rim height):

| `--stretch` | `--min-wall-angle 50` | `55` (default) | `60` |
|---|---|---|---|
| 1.2 (default) | 101.9 / 124.5 | 110.7 / 113.0 | 119.0 / 100.0 |

`--stretch 1.0` gives 118.8 / 83.2 and `--stretch 1.4` gives 103.0 / 144.2 (both at 55°).
Below about 50°, few strut directions can still be printed near the floor and opening.

### Common options (all shapes)

| Option | Default | Recommended | Meaning |
|---|---|---|---|
| `--height` | 100 | – | Body height; the eyelets add about 11 mm (box, cylinder) |
| `--eyelets` | box 4, cylinder/orb 3 | 0, 3, 4 | Hanging eyelets; 0 = standing pot. Box: 0, 2 or 4 (corners) |
| `--wall` | 5 | 4–6 | Wall thickness |
| `--strut` | 4.5 | 4–5 | Wall strut width (about 10 lines at 0.42 mm) |
| `--edge-round` | 1.0 | 0.5–1.5 | Rounding of the strut face edges |
| `--cell` | 17 | 14–20 | Wall cell size (seed spacing). Smaller means smaller openings, more weight |
| `--zstretch` | 1.2 | 1.2–1.4 | Vertical cell elongation (box, cylinder). Lower is more isotropic but needs more fixes (see below) |
| `--min-angle` | 35 | 35–40 | Minimum strut angle from horizontal that counts as printable |
| `--max-bridge` | 10 | 6–10 | Flat struts up to this horizontal span stay flat; longer ones become pointed arches |
| `--max-arch` | 25 | 20–25 | Wider shallow struts are not arched (kept as flat bridges) |
| `--grade-top` | 0.25 | 0.2–0.3 | Top fraction of the height where cells shrink toward the rim; 0 = off (fixed in v0.6.1) |
| `--grade-scale` | 0.6 | 0.5–0.7 | Cell size factor at the rim; 1 = off |
| `--merge` | 3 | 2–4 | Voronoi nodes closer than this are merged (avoids blobs) |
| `--band` | 5 | 4–6 | Bottom band height |
| `--rim` | 6 | 5–8 | Top rim height |
| `--floor` | 5 | 4–6 | Floor thickness |
| `--floor-cell` | 16 | 12–18 | Floor cell size. Use the smaller end for a bark medium |
| `--floor-strut` | 4.5 | 4–5 | Floor strut width |
| `--eyelet-hole` | 6 | 5–8 | Eyelet hole diameter (S-hook or wire) |
| `--eyelet-wire` | 4.5 | 4–5 | Eyelet ring width |
| `--fillet` | 1.5 | 1–2 | Fillet radius at strut junctions |
| `--voxel` | 0.4 | 0.3–0.5 | Mesh resolution. 0.3 is finer but about 2.4× the grid memory and slower |
| `--decimate` | 0.75 | 0.6–0.8 | Fraction of triangles removed; 0 = off. Retries once milder if the result is not watertight |
| `--layer` | 0.2 | your layer height | Layer height used by the island and cantilever checks |
| `--seed` | 7 | any | Random pattern |
| `--out` | `orchid_<shape>` | – | Output file prefix |
| `--graph-only` | off | – | Stop after the pattern and its metrics |

### Accessories (separate prints, all shapes)

| Option | Default | Recommended | Meaning |
|---|---|---|---|
| `--from-report` | – | a pot's `_report.json` | Take all pot parameters (and `--out`) from that pot, so the accessory matches it. Only accessory options come from the command line |
| `--parts` | `pot` | `pot`, `tray`, `planter`, combinations | Parts to build. Each writes its own STL + report: `orchid_<shape>_tray.stl` (saucer), `orchid_<shape>_tray_clip.stl` (clip), `orchid_<shape>_planter.stl` |
| `--tray` | `saucer` | `saucer`, `clip` | `saucer` for standing pots; `clip` snaps onto the pot's bottom band, for hanging pots |
| `--tray-gap` | 3 | 2–5 | Gap between the pot's outline (widest at any height) and the tray wall |
| `--tray-wall` | 2.4 | 2.0–3.0 | Tray wall thickness (about 6 lines for watertight walls) |
| `--tray-floor` | 2.0 | 1.6–2.4 | Tray floor thickness |
| `--standoff` | saucer 6, clip auto (~10–12) | 4–8 (saucer) | Height of the pot bottom above the tray floor. Water depth = standoff − 1 mm. Clip: chosen so the arm stays ≤ 2 % strain |
| `--tray-rise` | 4 | 2–8 | Tray wall height above the pot bottom |
| `--planter-gap` | 8 | 6–12 | Gap between the pot and the planter's inner wall |
| `--planter-wall` | 2.0 | 1.7–2.4 | Planter wall thickness, measured normal to the wall. 2.0 mm is about 5 lines at 0.42 mm, the upper end of the 4–5 perimeters for watertight walls |
| `--planter-floor` | 2.0 | 1.6–2.4 | Planter floor thickness |
| `--planter-ribs` | 4 | 3–4 | Vertical centring ribs inside (45° lead-in at the top) |
| `--planter-rib-clearance` | 1.0 | 0.5–1.5 | Clearance between ribs and pot |
| `--planter-eyelets` | 0 | 0, 3, 4 | Eyelets on the planter rim for a hanging planter; then use `--eyelets 0` for the pot. Box: 0, 2 or 4 |
| `--clips` | box 4, round 3 | 3–4 | Number of clips (box: side centres; round: 120° apart, each moved into a free lattice gap) |
| `--clip-width` | 8 | 6–10 | Clip width along the band |
| `--clip-thickness` | 1.6 | 1.2–2.0 | Clip arm thickness (4 lines). Thicker = stiffer, needs a longer arm |
| `--clip-hook` | 0.6 | 0.4–1.0 | Hook reach onto the flat band top, past the band's 1.5 mm edge rounding |
| `--clip-clearance` | 0.3 | 0.2–0.5 | Clearance between clip and pot; adjust after a test print |

**Clip tray:** flexible arms rise beside the bottom band. Their hooks snap over the band top
into lattice openings and sit on the flat band top.
- **Placement:** clip positions are checked against the strut paths themselves, over the
  whole hook height: at least 1.5 mm to every strut surface (`clip_strut_clearance_mm`).
- **Negative taper:** the band is widest at the bottom, so the arm stands further out and
  the snap deflection grows (`clip_snap_deflection_mm`); the standoff follows
  automatically. When seated, the pot only touches the
feet; when lifted 0.6–1.2 mm, all hooks engage. Arm strain is about 2 % (PETG snap-fit
guideline 1.5–2 %). To remove the tray, bend the arms outward one by one. The clearance
is an estimate: print a test tray first [UNVERIFIED] and adjust `--clip-clearance`.
Defaults: box about 218 ml / 98 g, cylinder about 171 ml / 84 g, orb about 187 ml / 84 g.

**Planter (`--parts planter`, decorative):** a closed outer pot. Box and cylinder planters
are flush with the pot rim; the orb gets a bowl up to its belly, because a closed onion
could not be assembled.
- **Insertion rule:** at every height, the inner wall follows the widest part of the pot
  *below* it, plus the gap, because the pot slides in from above. So the box planter is
  straight, a positive taper is followed, and a negative taper (pot widest at the floor)
  gives a straight cylinder.
- **Inside:** the pot stands on the same feet as in the tray (6 mm standoff, 5 mm water
  reservoir) and is centred by ribs with 1 mm clearance.
- **Defaults:** box about 256 g / 132 ml reservoir, cylinder about 215 g / 105 ml, orb bowl
  about 130 g / 66 ml.
- **Fit check:** like the tray, plus a sliding-in test (no overlap 10, 30 and 60 mm above
  the end position). A closed planter keeps air from the roots; use it only for less
  demanding plants.

**Matching tray and pot:**
- **Fingerprints:** pot and tray reports carry a `fingerprint`. `shape` hashes outline and
  band parameters; `pattern` hashes the strut paths themselves, which is what the clips
  depend on.
- **Fit check:** if the pot's STL and report sit next to the tray (same `--out` prefix,
  e.g. after `--parts pot,tray` or with `--from-report`), the tray run checks it
  automatically by mesh intersection. Seated, the pot must touch only the feet; for the
  clip tray, lifted 1.5 mm, every hook must catch the band.
- **On mismatch:** the run stops with an error, e.g. a tray built with another `--seed`.
  The result is in the report as `fit_check`. Without pot files next to the tray, the
  check is skipped and noted.

The pot stands on 8 feet under its solid bottom band, not under the floor lattice, so
drainage stays open. Each foot has an outer stop, which leaves about 1 mm side play. The
water line stays 1 mm below the pot bottom. Capacities with defaults: box about 117 ml,
cylinder/orb about 92–93 ml. For the orb the tray is belly-wide, because drips also fall
from the belly. The report gives `tray.pot_bottom_z_mm` and `tray.water_capacity_ml`.

### Effect of `--zstretch` (box, seed 7, other options default)

| zstretch | Nodes lifted | Struts kept flat | Opening height p90 / max | Largest opening (inscribed) |
|---|---|---|---|---|
| 1.0 | 19 | 25 | 24.7 / 32.5 | 23.1 |
| 1.2 | 9 | 14 | 30.5 / 35.5 | 22.6 |
| 1.5 | 2 | 4 | 36.5 / 42.0 | 24.1 |

Higher values give a taller, more "upward" pattern that needs fewer printability fixes.
Lower values give rounder cells but more flat bridges.

## How printability is ensured

- **Nodes:** every strut junction needs a strut below it steeper than `--min-angle`.
  Junctions that miss this are lifted slightly; no strut is added or removed.
- **Shallow struts:** short ones stay flat bridges. Longer ones become pointed arches,
  unless an arch would hit a neighbouring strut.
- **Rim:** the rim gets gussets in each gap, grown layer by layer from the struts. The
  rim then never bridges more than a few mm.

## Tests

```bash
python3 -m pytest -q tests/     # ~100 s
```

The tests check that box, cylinder and cylinder `--taper 8` (seed 7) give byte-identical
STLs and reports compared with `tests/golden/`. The box hash is the v0.4.0 release.
Regenerate the golden files only for an intended output change. The experimental orb has
only a smoke test of its strut graph; the trays (box, orb) are checked for one watertight
body, no islands and at least 80 ml capacity, and the clip trays for all clips placed and
≤ 2 % arm strain, and the box planter must fit (seated and sliding in) around the repo's
box pot. A clip tray for a pot with another pattern must be refused, and the
same tray via `--from-report` must pass the fit check. Every STL must stay watertight
after export and reload (`stl_roundtrip_watertight`). Where float32 export would split
off zero-volume slivers, they are removed (`stl_slivers_removed`).

## Report (`<out>_report.json`)

Key fields to check after each run:

| Field | Target |
|---|---|
| `watertight`, `bodies` | `true`, `1` |
| `unsupported_islands` | 0 |
| `cantilevers` | empty, or only short flat wall struts (reach ≤ about 6 mm) |
| `rim_unfilled_max_mm` | ≤ 6 (reach ≤ 3 mm) |
| `stability.tip_over_angle_deg` | orb only: tilt at which the empty pot tips over |
| `strut_crossings` | 0 |
| `wall_metrics.inscribed_dia_mm_p50/p90/max` | Opening size; fits roots or medium? |
| `mass_g_petg_1.27` | Weight estimate (PETG 1.23–1.29 g/cm³) |

`cantilevers` approximates OrcaSlicer's rule (TreeSupport: an overhang point more than
3 mm from its support), assuming a 30° threshold and 0.42 mm lines. ElegooSlicer is an
OrcaSlicer fork. With the actual Elegoo profile the slicer may be more lenient: seed 7
showed no warning there.

## Printing notes (PETG, 0.4 mm nozzle)

- Print upright on the floor, without supports. Optionally paint support on any strut
  listed under `cantilevers`.
- Lattices string easily with PETG, so tune retraction and travel.
- Trays need watertight settings, separate from the lattice: 4–5 perimeters, about
  102–105 % flow, about 240–245 °C, dried PETG
  ([Prusa](https://blog.prusa3d.com/watertight-3d-printing-part-2_53638/),
  [All3DP](https://all3dp.com/2/watertight-3d-print-tutorial/)).
- The weight with defaults is about 210–220 g (box) and about 170 g (cylinder, orb);
  the slicer gives the exact value.

## Limitations

- Box: on the curved corners, struts are about 14% wider on the outer face than on the
  inner face.
- Tapered cylinder: the pattern is computed in an approximate unrolled layout, so cells
  get slightly larger toward the wider end. Opening sizes in the report are in those
  layout units.
- Struts where an arch would collide stay as flat bridges, with no gusset underneath.
  Long ones show up in `cantilevers`.
- Since v0.6.1 the gussets no longer cause a staircase texture on the struts under the rim;
  only the gusset faces themselves come from a voxel field, smoothed by a 1-voxel blur.
- Orb (experimental): the opening metrics use a layout at the belly radius, so sizes near
  the floor and rim are overstated by up to about 25%. Mirror-snapping moves cells at
  the band and rim lines slightly. Not yet test-printed.
