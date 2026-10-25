# Example designs

These files are licensed under the Creative Commons Attribution 4.0 International License
(CC BY 4.0, full text in `LICENSE` here, <https://creativecommons.org/licenses/by/4.0/>).
Attribution: "Orchid pot design by Thomas Hartmann (orchid-pot)". The generator code in
the repository root is MIT-licensed (see `../LICENSE`); designs you generate with it are
yours.

| File | Design | Command (repository root) |
|---|---|---|
| `orchid_box.stl` (+ report) | 150 × 150 × 100 mm Voronoi box for Vanda orchids, 4 eyelets | `python3 orchid_pot.py box --out examples/orchid_box` |
| `orchid_cylinder_membrane.stl` (+ report) | Tapered cylinder, Voronoi relief with a 1.2 mm membrane, 25 % of the cells open, snap groove, 3 eyelets | `python3 orchid_pot.py cylinder --taper 8 --wall 3 --strut 3 --cell 10 --membrane --membrane-open 0.25 --clip-groove --groove-depth 2.0 --parts pot,tray --tray clip --out examples/orchid_cylinder_membrane` |
| `orchid_cylinder_membrane_tray_clip.stl` (+ report) | Clip-on drip tray for the cylinder (3 groove clips) | same command |
| `preview_*.png` | Renders of the box and the cylinder | `render.py` (Blender) |

Print: PETG, 0.4 mm nozzle, upright, no supports (see the main README).
