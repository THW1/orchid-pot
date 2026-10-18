# Changelog

All releases are annotated git tags. Details and verification results are in the commit
messages.

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
