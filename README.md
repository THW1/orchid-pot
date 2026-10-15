# Voronoi orchid basket generator

`orchid_pot.py` generates a printable STL of a hanging orchid basket (Vanda style). It is
a rounded-square box with planar Voronoi walls and floor, a solid bottom band, a top rim
and four hanging eyelets. It is tuned for PETG with a 0.4 mm nozzle, printed upright
without supports.

Design rationale and rejected alternatives are in the context header of `orchid_pot.py`.

## Requirements

Python 3 with `numpy`, `scipy`, `scikit-image`, `trimesh`, `fast_simplification` and
`shapely`. The previews also need Blender (tested with 5.2.2).

A full run takes about 30 s and peaks at about 1.3 GB RAM (default settings, 6 cores).

## Usage

```bash
python3 orchid_pot.py                          # defaults → orchid_pot.stl + orchid_pot_report.json
python3 orchid_pot.py --seed 3 --out pot_s3    # another random pattern, other file names
python3 orchid_pot.py --graph-only --seed 3    # pattern metrics only (~2 s), no mesh
blender -b --factory-startup --python render.py -- $PWD/orchid_pot.stl $PWD/preview
```

`render.py` needs **absolute** paths. It writes `preview_{persp,side,top,detail}.png`.

**Workflow:** compare seeds and parameters with `--graph-only` first, then build the mesh,
then check the report and previews.

## Options

Lengths are in mm. Recommended ranges are based on the runs in this project. Values
outside them work, but check the report.

| Option | Default | Recommended | Meaning |
|---|---|---|---|
| `--size` | 150 | – | Outer width and depth |
| `--height` | 100 | – | Body height; the eyelets add about 11 mm |
| `--corner-r` | 20 | 15–25 | Outer radius of the vertical edges; must exceed `--wall`/2 |
| `--wall` | 5 | 4–6 | Wall thickness |
| `--strut` | 4.5 | 4–5 | Wall strut width (about 10 lines at 0.42 mm) |
| `--edge-round` | 1.0 | 0.5–1.5 | Rounding of the strut face edges |
| `--cell` | 17 | 14–20 | Wall cell size (seed spacing). Smaller means smaller openings, more weight |
| `--zstretch` | 1.2 | 1.2–1.4 | Vertical cell elongation. Lower is more isotropic but needs more fixes (see below) |
| `--min-angle` | 35 | 35–40 | Minimum strut angle from horizontal that counts as printable |
| `--max-bridge` | 10 | 6–10 | Flat struts up to this horizontal span stay flat; longer ones become pointed arches |
| `--max-arch` | 25 | 20–25 | Wider shallow struts are not arched (kept as flat bridges) |
| `--grade-top` | 0.25 | 0.2–0.3 | Top fraction of the height where cells shrink toward the rim; 0 = off |
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
| `--out` | `orchid_pot` | – | Output file prefix |
| `--graph-only` | off | – | Stop after the pattern and its metrics |

### Effect of `--zstretch` (seed 7, other options default)

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

## Report (`<out>_report.json`)

Key fields to check after each run:

| Field | Target |
|---|---|
| `watertight`, `bodies` | `true`, `1` |
| `unsupported_islands` | 0 |
| `cantilevers` | empty, or only short flat wall struts (reach ≤ about 6 mm) |
| `rim_unfilled_max_mm` | ≤ 6 (reach ≤ 3 mm) |
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
- The weight is about 210–220 g with defaults; the slicer gives the exact value.

## Limitations

- On the curved corners, struts are about 14% wider on the outer face than on the inner
  face.
- Struts where an arch would collide stay as flat bridges, with no gusset underneath.
  Long ones show up in `cantilevers`.
- The gussets show a slight 0.4 mm staircase texture (the voxel pitch).
