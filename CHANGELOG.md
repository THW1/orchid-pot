# Changelog

All releases are annotated git tags. Details and verification results are in the commit
messages.

## v0.13.0 (2026-10-24)
- **Groove band grows inward** (`--clip-groove`).
  - The bottom band thickens inward to `--groove-depth` + 1.5 when the wall is thinner;
    `--band-wall` sets it explicitly, within [`--wall`, `--wall` + 2].
  - Thin walls keep a full-depth groove: `--wall 2.5 --groove-depth 2.0` gives a 3.5 mm
    band and a clip reach of 1.8 mm (was capped at 0.8 mm).
  - The outer face stays flush, so trays and planters fit unchanged.
  - Box, cylinder and orb.
- **`--membrane-open F`** (box, cylinder): a share F of the wall cells stays open, spread
  over the wall and picked by `--seed`; the rest is closed by the membrane. Combines with
  `--membrane-top`. The report lists `open_cells`.
- **Decimation on dense lattices.**
  - Fixed: dense lattices (`--cell 3`) were written undecimated (2.6M triangles).
  - Duplicate triangles from the decimator are dropped.
  - Retry ladder 0.75 / 0.6 / 0.5 / 0.3.
  - Example: a 90 mm `--cell 3` box now has 616k triangles instead of 1.23M.
- Existing outputs unchanged (goldens byte-identical); 3 new tests.

## v0.12.1 (2026-10-23)
- Fix: the groove clip tray did not fit a `--membrane` pot (fit check: seated overlap
  1.2 mm³). The floor plate shares its edge with the lattice floor; their smooth union
  pushed the band's outer face out by more than the 0.2 mm clip clearance. Now a hard
  union (box, cylinder, orb); band face as without membrane (±0.01 mm). New test.
- Mesh tail: harmless trimesh divide-by-zero warnings on zero-volume slivers are silenced.
  `voids_filled` no longer counts zero-volume slivers. Stats of rejected decimation
  candidates no longer leak into the report.

## v0.12.0 (2026-10-23)
- `--membrane [T]` for Voronoi pots (box, cylinder, orb): closes the cells with a thin wall
  (default 1.2 mm = 3 lines) on the inner face, plus a plate on the floor underside with
  drain holes (`--drain-*`, default 5).
  - The struts stay load-bearing relief outside; the inside is smooth.
  - Measured membrane 1.20–1.22 mm.
  - Box 283 g vs 214 g Voronoi.
- `--membrane-top Z`: the membrane ends at Z, so the cells above stay open for ventilation.
- Membrane is cut flat at the inner face after the smooth join, so the struts leave no
  inward bulge.
- Errors for misuse:
  - `--membrane*` with solid pots;
  - `--membrane-top` without `--membrane`;
  - `--drain-*` on open lattice pots;
  - a membrane thicker than wall − 1 mm;
  - lattice clip trays with `--membrane` (needs `--clip-groove`).
- Mesh tail: enclosed voids that decimation closes off are filled (`voids_filled`;
  orb `--membrane` had one of 1.3 mm³). Single-shell meshes untouched.
- Pots without membrane unchanged (goldens byte-identical); new test.

## v0.11.0 (2026-10-19)
- `--vent-slots N` for solid pots (box, cylinder, orb): vertical ventilation slots under the
  rim, for orchids that need more air than house plants.
  - `--vent-width` 4 mm and `--vent-height` 25 mm by default; 45° pointed top, so no
    support is needed. The apex sits at the rim bottom.
  - The slots sit between the eyelets, with at least 4 mm of wall between two slots
    (otherwise an error).
  - The report gains a `vent_slots` entry.
- `--vent-*` without `--wall-style solid`, and width/height without `--vent-slots`, are errors.
- Pots without slots unchanged; new test.

## v0.10.1 (2026-10-19)
- `--clip-hook` together with `--clip-groove` (also from `--from-report`) is now an error.
  It was silently ignored: in groove mode the hook reach is `--groove-depth` −
  `--clip-clearance`. New test; pots and trays unchanged.

## v0.10.0 (2026-10-18)
- `--wall-style solid`: plain closed pots for ordinary plants (box, cylinder, orb).
  - Solid wall and floor, default 2.4 mm (3.2 mm with `--clip-groove`).
  - Drainage holes `--drain-holes` (default 5) / `--drain-dia` (8 mm).
  - No strut graph.
  - Wall and band/rim are joined with a hard union: a smooth union of the coinciding
    surfaces had produced a 0.4 mm ridge at band and rim, which a groove clip touched.
  - Rim, eyelets, snap groove, saucer, groove clip tray and planter work unchanged.
- Voronoi pots unchanged (goldens byte-identical); new test for a solid box + saucer.

## v0.9.0 (2026-10-18)
- **Snap groove** (`--clip-groove`, pot option): a continuous groove in the bottom band
  (default band 8 mm, groove 1.5 mm deep, 3 mm opening). Flat floor on a layer boundary
  as the catch face; 45° ceiling as the lead-in. Prints without support on all shapes.
- **Clip tray in groove mode:**
  - clips exactly evenly spaced, independent of the lattice, default width 10 mm;
  - the arm follows the band (also on the orb's lean);
  - vertical play = clearance only (0.3 mm, was about 1.05 mm);
  - fingerprint check on shape incl. groove only.
- **Clip robustness (both modes):** 3 mm root fillet (a printed clip broke at its 0.5 mm
  root) and `--clip-strain` (default 1.5 %, was a fixed 2 %). Lattice clip arms therefore
  get longer, so the standoff rises.
- Verified with fit checks: box, cylinder, taper +8/−8, orb and a 125 mm orb (wall 4,
  4 clips); pots without support islands. New test for groove pot + tray.

## v0.8.3 (2026-10-18)
- Clips are placed as an evenly spaced set, rotated as a whole to fit the lattice, with a
  guaranteed minimum spacing of 75 % of the even spacing. Before, each clip took its
  nearest free spot on its own, and two clips could end up about 10° apart.
- Measured over 7 seeds: round shapes 115–120° (was 92–119°), box 86–89° (was 81–86°).
- New report fields `clip_spacing_min_pct_of_even` and `clip_shift_from_targets_mm`.

## v0.8.2 (2026-10-18)
- Accessory options that the requested parts would ignore are now an error: `--clip-*`
  without `--tray clip`, `--tray-*` without `tray`, `--planter-*` without `planter` in
  `--parts`. Before, a run with `--clip-*` but the default saucer silently built a tray
  without clips.
- `--clip-hook` longer than the flat band top (band width minus the edge rounding) is
  refused. The limit names the maximum, e.g. 1.0 mm for `--wall 4`.
- Fit check: the clip lift test follows the hook's vertical play (band-top allowance +
  clearance + 0.75 mm) instead of a fixed 1.5 mm, which false-failed with
  `--clip-clearance 1`. New report field `clip_vertical_play_mm`.
- A part that fails the fit check no longer stops the other parts; all parts are built,
  and the run fails at the end.

## v0.8.1 (2026-10-17)
- Review fixes for accessories:
  - exact pot height in the fit check;
  - signed band lean and separate arm/tip radii, so negative tapers work;
  - clip placement checked against the strut paths over the whole hook height;
  - per-clip engagement count;
  - `--fillet`-derived band-top allowance;
  - `--from-report` hardening;
  - fewer STL round trips;
  - strut paths reused.
- Planter wall default 2.4 → 2.0 mm.

## v0.8.0 (2026-10-17)
- Planter (`--parts planter`): closed outer pot; the orb gets a bowl. The inner wall
  follows the insertion envelope; feet, centring ribs and optional rim eyelets.

## v0.7.2 (2026-10-17)
- Pot/accessory fingerprints, automatic mesh fit check, `--from-report`, STL sliver
  repair (`export_safe`).

## v0.7.1 (2026-10-17)
- Clip-on tray for hanging pots; STL round-trip check of every decimated mesh.

## v0.7.0 (2026-10-17)
- `accessories.py` with the drip tray (saucer) for all shapes.

## v0.6.1 (2026-10-17)
- Review fixes: orb mirror wrap, `--grade-top 0`, orb floor joint, arch clash check, input
  validation. Gussets merged without shifting existing surfaces.

## v0.6.0 (2026-10-17)
- Experimental `orb` (onion) shape in `orb.py`; shared helpers in `potlib.py`;
  `requirements.txt`.

## v0.5.0 (2026-10-17)
- Shape subcommands; cylinder with optional `--taper`; optional eyelets; regression tests.

## v0.4.0 (2026-10-15)
- Box with planar Voronoi wall, graded cells, rim gussets.
