# Project context: orchid-pot
Verified-at: v0.6.0 (git tag)

## Purpose
This project generates printable STLs of Voronoi orchid baskets and pots for PETG with a
0.4 mm nozzle, printed upright without supports. Shapes are subcommands: `box`,
`cylinder` (optionally tapered) and the experimental `orb` (onion). The why of each
design choice lives in the context headers of `orchid_pot.py`, `potlib.py` and
`orb.py`. User documentation is in `README.md`.

## Layout (load-bearing files only)
```
orchid_pot.py           CLI + box/cylinder pipeline: frame → strut graph → SDF → potlib.finish
potlib.py               shape-independent helpers: SDF primitives, Poisson, Grid, floor, metrics,
                        meshing, checks, finish() (mesh → decimate → check → write)
orb.py                  EXPERIMENTAL orb pipeline (sphere Voronoi, 3D printability, shell SDF)
render.py               Blender (headless, Cycles) preview renders of an STL
tests/test_regression.py  byte-identical STL + report for box, cylinder, cylinder --taper 8;
                        orb strut-graph smoke test
tests/golden/           reference reports and STL hashes (box hash = v0.4.0 release)
README.md               usage, options tables, recommended ranges, report targets
v2-3d/                  archived v2 (3D-lattice wall): code, STL, report, renders; not maintained
```
Generated (not source): `orchid_<shape>.stl`, `orchid_<shape>_report.json`,
`preview_<shape>_*.png`.

## Flow
1. `orchid_pot.py::main` → `parse_args` (subcommand = shape). For `orb` it imports
   `orb` and hands over to `orb.py::run`. Otherwise `make_frame` returns a `BoxFrame`
   or `ConeFrame`.
2. `wall_graph`: 2D Voronoi on the unrolled wall, lift, classify, arches. Its output feeds
   `wall_metrics` and `rim_gussets` (the gap list).
3. `floor_segments`: 2D Voronoi of the floor, clipped by `footprint_sdf`.
4. `build_sdf` calls `wall_pattern` (2D field with struts, grown gussets and eyelets),
   then composes each z slice using `wall_fields`, band, rim and floor.
5. `potlib.py::finish`: `to_mesh` (marching cubes, `clean`, micro-void fill),
   decimation (one milder retry), `check_mesh` (with `cantilevers`), optional
   `extra_checks` (orb: stability), then the STL and report are written.
6. `--graph-only` stops after steps 2–3 and prints the metrics.
7. Orb: `orb.py::orb_graph` → `orb_metrics`, `rim_gaps` → `build_orb_sdf` (shell
   projection, struts, band/rim/floor, eyelets, `grow_rim_gussets`) →
   `potlib.py::finish`.

## Relationships
- `orchid_pot.py::parse_args` → `README.md:Options`: the tables mirror flags, defaults
  and ranges.
- Frame interface used by the pipeline: `rise`, `hor`, `scale`, `unmap_u`,
  `wall_fields`, `footprint_sdf`, `eyelet_u`, `P`, `ku_min`, `extent`. It is implemented
  by `orchid_pot.py::BoxFrame` (via `WallFrame`) and `orchid_pot.py::ConeFrame`.
- `orchid_pot.py::wall_graph` → `orchid_pot.py::wall_pattern`: polylines in unwrapped u;
  the consumer shifts by ±P.
- `orchid_pot.py::rim_gussets` → `orchid_pot.py::wall_pattern`: list of
  (gap midpoint u, gap width), only gaps wider than one strut.
- `orchid_pot.py::wall_pattern` → `orchid_pot.py::build_sdf`: F rows equal the voxel
  z axis.
- `orb.py` → `potlib.py` only. The orb never imports `orchid_pot.py`, and
  `orchid_pot.py` imports `orb` only for the `orb` subcommand.
- `orb.py::run` sets `args.height` to the derived rim height before the potlib checks,
  which use it for the rim zone.
- `tests/test_regression.py` → `tests/golden/`: golden files change only together with
  an intended output change.
- `potlib.py::cantilevers` → OrcaSlicer TreeSupport rule (external; source URL in
  the README).
- `orchid_pot.py::main` → `render.py`: the STL path, passed as an absolute argument.

## Run and test
```bash
python3 orchid_pot.py cylinder --graph-only   # ~2 s: pattern metrics
python3 orchid_pot.py box                     # ~30 s, ~1.3 GB: STL + report
python3 -m pytest -q tests/                   # ~100 s regression
blender -b --factory-startup --python render.py -- $PWD/orchid_box.stl $PWD/preview_box
python3 orchid_pot.py orb                     # ~95 s, ~2.4 GB (experimental)
python3 -m pyflakes *.py tests/*.py
```

## Conventions
- All lengths are in mm. Wall coordinates are (u, z); u is chart arc length along the
  wall mid-curve and periodic with P.
- New shapes must not change box or cylinder output: run the regression tests. The
  experimental orb lives in `orb.py`; changes to `potlib.py` affect all shapes.
- Every check result goes into the JSON report (`stats` dict), never only to stdout.
- Releases are annotated git tags `vX.Y.Z`. The zsh shell does not split `$var`: use
  `${=cfg}`.

## Known limitations
See `Limitations:` in the `orchid_pot.py` header, and the README "Limitations" section.
