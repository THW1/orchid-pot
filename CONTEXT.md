# Project context: orchid-pot
Verified-at: uncommitted, 2026-10-15

## Purpose
This project generates a printable STL of a Voronoi hanging orchid basket for PETG with a
0.4 mm nozzle, printed upright without supports. The why of each design choice lives in
the context header of `orchid_pot.py`. User documentation is in `README.md`.

## Layout (load-bearing files only)
```
orchid_pot.py   generator CLI: strut graph → SDF → mesh → checks → STL + JSON report
render.py       Blender (headless, Cycles) preview renders of an STL
README.md       usage, options table, recommended ranges, report targets
v2-3d/          archived v2 (3D-lattice wall): code, STL, report, renders; not maintained
```
Generated (not source): `orchid_pot.stl`, `orchid_pot_report.json`, `preview_*.png`, `*.log`.

## Flow
1. `orchid_pot.py::main` → `parse_args`, then `WallFrame` (geometry of the unrolled wall).
2. `wall_graph`: 2D Voronoi on the unrolled wall, lift, classify, arches. Its output feeds
   `wall_metrics` and `rim_gussets` (the gap list).
3. `floor_segments`: 2D Voronoi of the floor.
4. `build_sdf` calls `wall_pattern` (2D field with struts, grown gussets and eyelets),
   then composes each z slice with band, rim and floor.
5. `to_mesh` (marching cubes plus `clean`), decimation (one milder retry), `check_mesh`
   (with `cantilevers`), then the STL and report are written.
6. `--graph-only` stops after steps 2–3 and prints the metrics.

## Relationships
- `orchid_pot.py::parse_args` → `README.md:Options`: the table mirrors flags, defaults
  and ranges.
- `orchid_pot.py::wall_graph` → `orchid_pot.py::wall_pattern`: polylines in unwrapped u;
  the consumer shifts by ±P.
- `orchid_pot.py::rim_gussets` → `orchid_pot.py::wall_pattern`: list of
  (gap midpoint u, gap width), only gaps wider than one strut.
- `orchid_pot.py::wall_pattern` → `orchid_pot.py::build_sdf`: F rows equal the voxel
  z axis.
- `orchid_pot.py::WallFrame` → `orchid_pot.py::build_sdf` via `unmap_u` (the x,y → u
  lookup).
- `orchid_pot.py::cantilevers` → OrcaSlicer TreeSupport rule (external; source URL in
  the README).
- `orchid_pot.py::main` → `render.py`: the STL path, passed as an absolute argument.

## Run and test
```bash
python3 orchid_pot.py --graph-only            # ~2 s: pattern metrics
python3 orchid_pot.py                         # ~30 s, ~1.3 GB: STL + report
blender -b --factory-startup --python render.py -- $PWD/orchid_pot.stl $PWD/preview
python3 -m pyflakes orchid_pot.py render.py
```
There is no test suite. Verification is the report:
- `watertight` true, `bodies` 1, `unsupported_islands` 0, `strut_crossings` 0
- `cantilevers` without a `rim_underside` entry
- a look at the previews

## Conventions
- All lengths are in mm. Wall coordinates are (u, z); u is arc length along the wall
  mid-curve and periodic with P.
- Every check result goes into the JSON report (`stats` dict), never only to stdout.
- The zsh shell does not split `$var`: use `${=cfg}` when looping over option strings.

## Known limitations
See `Limitations:` in the `orchid_pot.py` header, and the README "Limitations" section.
