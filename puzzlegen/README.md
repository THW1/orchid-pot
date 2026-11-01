# puzzlegen

Generator for irregular interlocking puzzles (Voronoi, Penrose, hat/Spectre tilings)
for 3D printing (STL/3MF) and laser or plotter cutting (SVG/DXF), with an optional
matching frame and a verified unique solution.

**Status:** design phase. See [`DESIGN.md`](DESIGN.md). No code yet.

This folder is self-contained and independent of the orchid-pot code. It lives here
temporarily. To move it into its own repository with history:

```bash
git subtree split --prefix=puzzlegen -b puzzlegen-only
# then push branch puzzlegen-only to the new repository's main
```
