#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# === CONTEXT ===
# Purpose:     Generate a printable STL of an orchid basket/pot (Vanda style) with planar
#              Voronoi walls and floor, solid bottom band, top rim and optional flat hanging
#              eyelets. Shapes are subcommands: box (rounded square) and cylinder (optionally
#              tapered). Target: PETG, 0.4 mm nozzle, printed upright without supports.
# Decisions:   Walls are a planar (2D) Voronoi on the unrolled wall (u = arc length
#              along the wall mid-curve, z), extruded through the wall thickness like
#              the floor. Rejected: a 3D Voronoi lattice in a 12 mm shell (v2, archived
#              in v2-3d/). With ~17 mm cells it was only one cell deep, so it read as the
#              same 2D pattern on both faces plus through-struts.
#              Seeds are mirrored at the band/rim lines, so the pattern closes exactly
#              into them. They are copied by +-P for periodicity; only edges whose
#              midpoint lies in the central copy are kept, because the outer copies have
#              wrong cells at their far ends (this leaked 33 spurious faces into v2).
#              z is compressed by --zstretch before the Voronoi: taller cells give more
#              steep struts.
#              Printability is enforced on the strut graph. Every node must reach a lower
#              node via a steep strut (no node starts in mid-air). Nodes that miss this
#              are lifted (bottom-up) until their best lower strut is steep enough, so
#              the Voronoi topology is untouched (details at the lift step).
#              Long shallow struts between supported nodes become pointed arches (legs at
#              >= min-angle+5 deg). An arch must keep one strut width of clearance from
#              other struts (the plane has no depth to dodge into); otherwise it is
#              retried flatter, then kept as a flat bridge. It is never dropped, since a
#              dropped strut merges two cells. Under the rim, pointed-arch gussets
#              replace bridges (rim_gussets, wall_pattern). Cells shrink over the top 25%,
#              so more struts reach the rim and the gussets stay small.
#              Rejected: random vs jittered-lattice seeding made little difference once
#              lifting was in place, so the plain random (Poisson) Voronoi stays.
#              Shapes are frames (BoxFrame, ConeFrame) behind one interface: rise, hor,
#              scale, unmap_u, wall_fields, footprint_sdf, eyelet_u. The pipeline never
#              branches on the shape. BoxFrame returns the v0.4.0 expressions exactly, so
#              the box STL is byte-identical to v0.4.0 (tests/test_regression.py). The
#              experimental orb is meant to live in its own module, not as a third frame
#              here (its walls lean and are not developable).
#              A taper's metric factor is always taken at a strut's lower end. Lift and
#              classify must agree on it; with mismatched heights, the -8 deg taper lost
#              34 nodes to the safety net.
#              Mesh = SDF + marching cubes. The wall pattern is a 2D smooth-min distance
#              field in (u, z), looked up per voxel via the inverse wall mapping. This
#              gives node fillets and rounded face edges, with no booleans.
# Constraints: numpy, scipy, scikit-image, trimesh, fast_simplification, shapely.
#              The voxel grid at 0.4 mm is ~45M float32 (~180 MB); peak RSS ~1.3 GB,
#              ~30 s on 6 cores. The grid grows with 1/voxel^3.
# Related:     render.py renders the STL. README.md documents the CLI; its options
#              table mirrors parse_args, so the two change together.
# Invariants:  The band/rim mirror lines sit inside the band/rim solids; edges on them
#              are anchors only and are not drawn. Grid origin has an odd offset, so no
#              sample lands on a flat face (zero-area marching-cubes slivers otherwise).
# Tests:       tests/test_regression.py (byte-identical STL + report for box, cylinder,
#              cylinder --taper 8). Per run, the report checks: watertight, one body,
#              islands per 0.2 mm layer, overhang area, opening sizes, approximate
#              OrcaSlicer cantilever check. Preview with render.py.
# Limitations: In-plane distances are measured on the mid-curve, so at the corners
#              struts are ~14% wider on the outer face than on the inner face (box only).
#              Tapered cylinders use an approximate chart (see ConeFrame): cells and struts
#              follow the local radius only to first order, and the opening metrics are in
#              chart units. Wall arches that clash stay flat bridges; the cantilever check
#              reports them.
# === END CONTEXT ===
import argparse
import json
import math
import time

import numpy as np
import shapely
import trimesh
import fast_simplification
from scipy import ndimage
from scipy.spatial import Voronoi, cKDTree
from shapely.ops import unary_union
from skimage.measure import marching_cubes


# CTX: CLI contract: one subcommand ("noun") per shape; options shared by all shapes live in
#      a parent parser, so `orchid_pot.py <shape> --help` lists only what applies.
#      README.md's options tables mirror every flag, default and range.
# Used-by: orchid_pot.py::main, README.md:Options
def parse_args(argv=None):
    common = argparse.ArgumentParser(add_help=False)
    g = common.add_argument_group("common")
    g.add_argument("--height", type=float, default=100.0, help="body height without eyelets [mm]")
    g.add_argument("--wall", type=float, default=5.0, help="wall thickness [mm]")
    g.add_argument("--strut", type=float, default=4.5, help="wall strut width in the wall plane [mm]")
    g.add_argument("--cell", type=float, default=17.0, help="seed spacing of the wall Voronoi [mm]")
    g.add_argument("--zstretch", type=float, default=1.2, help="vertical cell elongation factor")
    g.add_argument("--eyelets", type=int, default=None,
                   help="hanging eyelets: 0 = none (standing pot), else >= 2 (default: shape)")
    g.add_argument("--seed", type=int, default=7, help="random pattern")
    g.add_argument("--out", default=None, help="output prefix (default: orchid_<shape>)")
    g.add_argument("--graph-only", action="store_true", help="stop after the strut graph and its metrics")
    g = common.add_argument_group("details (rarely changed)")
    g.add_argument("--edge-round", type=float, default=1.0, help="rounding of strut face edges [mm]")
    g.add_argument("--band", type=float, default=5.0, help="bottom band height [mm]")
    g.add_argument("--rim", type=float, default=6.0, help="top rim height [mm]")
    g.add_argument("--floor", type=float, default=5.0, help="floor lattice thickness [mm]")
    g.add_argument("--floor-cell", type=float, default=16.0, help="floor Voronoi seed spacing [mm]")
    g.add_argument("--floor-strut", type=float, default=4.5, help="floor strut width [mm]")
    g.add_argument("--eyelet-hole", type=float, default=6.0, help="eyelet hole diameter [mm]")
    g.add_argument("--eyelet-wire", type=float, default=4.5, help="eyelet ring width [mm]")
    g.add_argument("--fillet", type=float, default=1.5, help="smooth-union radius at nodes [mm]")
    g = common.add_argument_group("printability (advanced)")
    g.add_argument("--min-angle", type=float, default=35.0, help="min strut elevation [deg]")
    g.add_argument("--max-bridge", type=float, default=10.0, help="max horizontal span of flat bridges [mm]")
    g.add_argument("--max-arch", type=float, default=25.0, help="max horizontal span of an arch [mm]")
    g.add_argument("--grade-top", type=float, default=0.25, help="top fraction of height with graded cells")
    g.add_argument("--grade-scale", type=float, default=0.6, help="cell size factor at the rim")
    g.add_argument("--merge", type=float, default=3.0, help="merge Voronoi nodes closer than this [mm]")
    g.add_argument("--layer", type=float, default=0.2, help="layer height for island check [mm]")
    g = common.add_argument_group("mesh")
    g.add_argument("--voxel", type=float, default=0.4, help="SDF grid pitch [mm]")
    g.add_argument("--decimate", type=float, default=0.75, help="fraction of triangles to remove")

    p = argparse.ArgumentParser(description="Voronoi orchid basket STL generator")
    sub = p.add_subparsers(dest="shape", required=True, metavar="SHAPE")
    b = sub.add_parser("box", parents=[common], help="rounded-square box (default 4 eyelets)",
                       description="Rounded-square box with planar Voronoi walls.")
    g = b.add_argument_group("box")
    g.add_argument("--size", type=float, default=150.0, help="outer width/depth [mm]")
    g.add_argument("--corner-r", type=float, default=20.0, help="outer vertical edge radius [mm]")
    b.set_defaults(eyelets_default=4)
    c = sub.add_parser("cylinder", parents=[common], help="round pot, optionally tapered (default 3 eyelets)",
                       description="Cylinder or cone frustum with planar Voronoi walls.")
    g = c.add_argument_group("cylinder")
    g.add_argument("--diameter", type=float, default=150.0, help="outer diameter at the rim [mm]")
    g.add_argument("--taper", type=float, default=0.0,
                   help="wall lean from vertical [deg]; > 0 widens upward (flowerpot), max +-20")
    c.set_defaults(eyelets_default=3)
    args = p.parse_args(argv)
    if args.eyelets is None:
        args.eyelets = args.eyelets_default
    del args.eyelets_default
    if args.eyelets == 1 or args.eyelets < 0:
        p.error("--eyelets must be 0 or >= 2")
    if args.shape == "cylinder" and abs(args.taper) > 20:
        p.error("--taper is limited to +-20 deg (wall lean and the chart approximation)")
    if args.out is None:
        args.out = f"orchid_{args.shape}"
    return args


# ---------- geometry helpers ----------

# CTX: Single source of the unrolled wall coordinate (u = arc length on the mid-curve).
# Used-by: orchid_pot.py::wall_graph, orchid_pot.py::wall_pattern, orchid_pot.py::build_sdf
# Invariant: unmap_u inverts map outside the core square (round trip ~1e-13 mm). The build
#            no longer calls map; it is kept as the reference that unmap_u is checked against.
class WallFrame:
    """Rounded-square mid-curve of the wall; maps (u, z, w) <-> (x, y, z)."""

    def __init__(self, half, corner_r, wall):
        self.t = wall
        self.a = half - wall / 2                 # mid-curve half size
        self.r = corner_r - wall / 2             # mid-curve corner radius
        assert self.r > 0, "corner radius must exceed wall/2"
        self.c = self.a - self.r
        self.side = 2 * self.c
        self.quad = self.side + math.pi * self.r / 2
        self.P = 4 * self.quad

    def map(self, u, z, w):
        u = np.mod(u, self.P)
        k = np.floor(u / self.quad).astype(int)
        ul = u - k * self.quad
        on_side = ul < self.side
        th = np.where(on_side, 0.0, (ul - self.side) / self.r)
        px = np.where(on_side, self.a, self.c + self.r * np.cos(th))
        py = np.where(on_side, -self.c + ul, self.c + self.r * np.sin(th))
        nx, ny = np.cos(th), np.sin(th)
        ang = k * (math.pi / 2)
        ca, sa = np.cos(ang), np.sin(ang)
        off = w - self.t / 2
        x = ca * px - sa * py + off * (ca * nx - sa * ny)
        y = sa * px + ca * py + off * (sa * nx + ca * ny)
        return np.stack([x, y, z], axis=-1)

    def unmap_u(self, x, y):
        """Arc length u of the nearest mid-curve point (valid outside the core square)."""
        qx, qy = np.clip(x, -self.c, self.c), np.clip(y, -self.c, self.c)
        vx, vy = x - qx, y - qy
        psi = np.arctan2(vy, vx) % (2 * math.pi)
        side = (vx == 0) | (vy == 0)             # exact: clip returns x itself on a side
        k = np.where(side, np.round(psi / (math.pi / 2)),
                     np.floor(psi / (math.pi / 2))).astype(int) % 4
        ang = -k * (math.pi / 2)
        ry = np.sin(ang) * qx + np.cos(ang) * qy  # core point rotated back onto side 0
        ul = np.where(side, ry + self.c, self.side + self.r * (psi - k * (math.pi / 2)))
        return (k * self.quad + ul) % self.P

    # --- frame interface shared with ConeFrame: vertical walls, chart == real metric ---
    ku_min = 1.0

    def rise(self, h, tan, z=None):
        """Vertical rise a strut of chart-horizontal length h needs to reach atan(tan)."""
        return tan * h

    def hor(self, h, z=None):
        """Real horizontal length of a chart-horizontal length h at height z."""
        return h

    def scale(self, z):
        """(u, z) chart -> local real metric factors at height(s) z."""
        return np.ones_like(np.asarray(z, float)), 1.0

    def eyelet_u(self, n):
        # middle of the corner arcs for n = 4; same float expression as before (v0.4.0)
        return [i * (self.P / n) + self.side + self.r * math.pi / 4 for i in range(n)]


class BoxFrame(WallFrame):
    def __init__(self, args):
        super().__init__(args.size / 2, args.corner_r, args.wall)
        self.half_out, self.corner_out = args.size / 2, args.corner_r
        self.extent = args.size / 2
        self._dout = None

    def footprint_sdf(self, x, y):
        return rounded_square_sdf(x, y, self.half_out, self.corner_out)

    def wall_fields(self, X2, Y2, z):
        """(outer-surface SDF, |distance from the mid-surface|) of the slice at height z."""
        if self._dout is None:                  # vertical walls: the same for every slice
            self._dout = rounded_square_sdf(X2, Y2, self.half_out, self.corner_out)
            self._pw = np.abs(self._dout + self.t / 2)
        return self._dout, self._pw

    def describe(self):
        return {"outer_mm": [2 * self.half_out, 2 * self.half_out], "corner_r_mm": self.corner_out}


# CTX: Cylinder / cone frustum. Chart u = theta * r_ref (r_ref = mid-surface radius at mid
#      lattice height); with a taper the chart is not isometric, so lengths use the local
#      factor ku(z) = r(z) / r_ref and the slant factor 1/cos(taper).
# Used-by: orchid_pot.py::make_frame
# Rejected: exact cone development (annular sector): periodicity becomes a rotation and the
#           band/rim mirror lines become arcs; for |taper| <= 20 deg the chart error is a few %.
class ConeFrame:
    def __init__(self, args, z_top):
        self.t = args.wall
        self.H = args.height
        tt = math.radians(args.taper)
        self.tan_t, self.cos_t = math.tan(tt), math.cos(tt)
        self.r_top_out = args.diameter / 2
        self.r_ref = self.r_mid((args.band / 2 + args.height - args.rim / 2) / 2)
        assert self.r_mid(0) > self.t and self.r_mid(z_top) > self.t, "taper too strong"
        self.P = 2 * math.pi * self.r_ref
        self.ku_min = min(self.ku(0.0), self.ku(z_top))
        self.extent = max(self.r_out(0.0), self.r_out(z_top)) + self.t

    def r_out(self, z):
        return self.r_top_out - (self.H - np.asarray(z, float)) * self.tan_t

    def r_mid(self, z):
        return self.r_out(z) - self.t / 2 / self.cos_t

    def ku(self, z):
        return self.r_mid(z) / self.r_ref

    def rise(self, h, tan, z=0.0):
        # rise r over tangential run h*ku on a wall leaning by the taper: the strut also
        # moves radially by r*tan_t, so tan(e) = r / hypot(h*ku, r*tan_t) >= tan
        return tan * h * float(self.ku(z)) / math.sqrt(1 - (tan * self.tan_t) ** 2)

    def hor(self, h, z=0.0):
        return h * float(self.ku(z))

    def scale(self, z):
        return self.ku(z), 1 / self.cos_t

    def eyelet_u(self, n):
        return [i * (self.P / n) for i in range(n)]

    def unmap_u(self, x, y):
        return (np.arctan2(y, x) % (2 * math.pi)) * self.r_ref

    def footprint_sdf(self, x, y):
        return np.hypot(x, y) - float(self.r_out(0.0))

    def wall_fields(self, X2, Y2, z):
        wn = (np.hypot(X2, Y2) - float(self.r_mid(z))) * self.cos_t   # normal distance to mid-surface
        return wn - self.t / 2, np.abs(wn)

    def describe(self):
        return {"outer_dia_bottom/top_mm": [round(2 * float(self.r_out(0.0)), 1), 2 * self.r_top_out],
                "taper_deg": round(math.degrees(math.atan(self.tan_t)), 2)}


def make_frame(args):
    if args.shape == "box":
        return BoxFrame(args)
    z_top = args.height + (args.eyelet_hole + 2 * args.eyelet_wire + 1 if args.eyelets else 0)
    return ConeFrame(args, z_top)


def rounded_square_sdf(x, y, half, r):
    qx = np.abs(x) - (half - r)
    qy = np.abs(y) - (half - r)
    out = np.hypot(np.maximum(qx, 0), np.maximum(qy, 0))
    return out + np.minimum(np.maximum(qx, qy), 0) - r


def rounded_box2(px, pz, hx, hz, rr):
    qx = px - (hx - rr)
    qz = pz - (hz - rr)
    out = np.hypot(np.maximum(qx, 0), np.maximum(qz, 0))
    return out + np.minimum(np.maximum(qx, qz), 0) - rr


def round_isect(a, b, rr):
    """Intersection of two SDFs with the shared edge rounded by rr."""
    qa, qb = a + rr, b + rr
    return np.hypot(np.maximum(qa, 0), np.maximum(qb, 0)) + np.minimum(np.maximum(qa, qb), 0) - rr


def smin(a, b, k):
    h = np.clip(0.5 + 0.5 * (b - a) / k, 0.0, 1.0)
    return b * (1 - h) + a * h - k * h * (1 - h)


def poisson(rng, lo, hi, spacing, period=None, tries=60000):
    """Dart throwing. spacing: float, or callable(point) -> local spacing (pairs use the
    mean of both). period: index of a periodic axis -> length."""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    sp = spacing if callable(spacing) else (lambda p: spacing)
    pts, rad = [], []
    for _ in range(tries):
        p = lo + rng.random(len(lo)) * (hi - lo)
        r = sp(p)
        if pts:
            d = np.abs(np.asarray(pts) - p)
            if period:
                for ax, L in period.items():
                    d[:, ax] = np.minimum(d[:, ax], L - d[:, ax])
            if np.any(np.sqrt((d ** 2).sum(1)) < (np.asarray(rad) + r) / 2):
                continue
        pts.append(p); rad.append(r)
    return np.asarray(pts)


def pt_seg_dist(P, A, B):
    """Distance of points P to segments A-B (broadcast over rows)."""
    ab = B - A
    t = np.clip(((P - A) * ab).sum(-1) / np.maximum((ab * ab).sum(-1), 1e-12), 0, 1)
    return np.linalg.norm(A + t[..., None] * ab - P, axis=-1)


def seg_seg_dist(a, b, A, B):
    """Min distance between segment a-b and each segment A[i]-B[i] in 2D."""
    d = np.minimum.reduce([pt_seg_dist(a[None], A, B), pt_seg_dist(b[None], A, B),
                           pt_seg_dist(A, a[None], b[None]), pt_seg_dist(B, a[None], b[None])])
    cr = lambda o, p, q: (p[..., 0] - o[..., 0]) * (q[..., 1] - o[..., 1]) - \
        (p[..., 1] - o[..., 1]) * (q[..., 0] - o[..., 0])
    hit = (cr(a[None], b[None], A) * cr(a[None], b[None], B) < 0) & \
          (cr(A, B, a[None]) * cr(A, B, b[None]) < 0)
    d[hit] = 0.0
    return d


# ---------- wall graph (planar Voronoi on the unrolled wall) ----------

# CTX: Returns wall struts as polylines in UNWRAPPED (u, z); u may leave [0, P), so every
#      consumer wraps or shifts by +-P. Band/rim anchor edges are not returned.
# Used-by: orchid_pot.py::main, orchid_pot.py::wall_metrics, orchid_pot.py::rim_gussets, orchid_pot.py::wall_pattern
def wall_graph(args, fr, rng, stats):
    s = args.zstretch
    z0, z1 = args.band / 2 + 0.5, args.height - args.rim / 2
    eps = 1e-4
    # graded cells: spacing shrinks linearly over the top --grade-top fraction of the
    # height to --grade-scale x cell at the rim -> more struts reach the rim
    zg = args.height * (1 - args.grade_top)
    cell_at = lambda p: args.cell * (1 - (1 - args.grade_scale) *
                                     min(max((p[1] * s - zg) / (z1 - zg), 0.0), 1.0))
    seeds = poisson(rng, [0, z0 / s], [fr.P, z1 / s], cell_at, period={0: fr.P})
    B = np.vstack([seeds, seeds * [1, -1] + [0, 2 * z0 / s], seeds * [1, -1] + [0, 2 * z1 / s]])
    C = np.vstack([B, B + [fr.P, 0], B - [fr.P, 0]])
    vor = Voronoi(C)
    V = vor.vertices * [1, s]
    inside = (V[:, 1] > z0 - eps) & (V[:, 1] < z1 + eps)
    raw = [(a, b) for a, b in vor.ridge_vertices
           if a >= 0 and b >= 0 and inside[a] and inside[b]
           and 0 <= (V[a, 0] + V[b, 0]) / 2 < fr.P]

    used = sorted({v for e in raw for v in e})
    idx = {v: n for n, v in enumerate(used)}
    nodes = np.column_stack([V[used, 0] % fr.P, V[used, 1]])

    def dvec(i, j):
        d = nodes[j] - nodes[i]
        d[0] = (d[0] + fr.P / 2) % fr.P - fr.P / 2
        return d

    # union-find: seam duplicates (same point via two copies) and short edges
    parent = list(range(len(nodes)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i

    tree = cKDTree(np.column_stack([nodes[:, 0] % fr.P, nodes[:, 1]]), boxsize=[fr.P, 1e6])
    for i, j in tree.query_pairs(1e-3):
        parent[find(i)] = find(j)
    for a, b in raw:
        if np.linalg.norm(dvec(idx[a], idx[b])) < args.merge:
            parent[find(idx[a])] = find(idx[b])
    groups = {}
    for i in range(len(nodes)):
        groups.setdefault(find(i), []).append(i)
    new_id, merged = {}, []
    for g in groups.values():
        p = nodes[g[0]] + np.mean([dvec(g[0], i) for i in g], axis=0)
        p[0] %= fr.P
        zs = nodes[g, 1]
        if np.any(np.abs(zs - z0) < eps): p[1] = z0      # keep anchors on their lines
        elif np.any(np.abs(zs - z1) < eps): p[1] = z1
        for i in g:
            new_id[i] = len(merged)
        merged.append(p)
    nodes = np.asarray(merged)
    E = set()
    for a, b in raw:
        i, j = new_id[idx[a]], new_id[idx[b]]
        if i != j:
            E.add((min(i, j), max(i, j)))
    stats["voronoi_nodes"], stats["voronoi_edges"] = len(nodes), len(E)

    # ---- printability ----
    bottom = np.abs(nodes[:, 1] - z0) < eps
    top = np.abs(nodes[:, 1] - z1) < eps
    tan_min = math.tan(math.radians(args.min_angle))
    tan_arch = math.tan(math.radians(args.min_angle + 5))
    anchor = lambda i, j: (bottom[i] and bottom[j]) or (top[i] and top[j])

    # Lift: bottom-up, a node without a steep strut to a lower node is raised until its
    # cheapest lower strut reaches min-angle (capped at the rim line, where it anchors).
    # The topology stays a pure Voronoi: no strut is added or removed, so cells cannot
    # merge into tall holes. Rejected: deleting such nodes (v1) and adding rescue struts
    # (v2/v3 draft). In one plane the extra struts collide with neighbours, and the
    # deletions left 70-90 mm holes.
    adj = {i: [j for e in E if i in e and not anchor(*e) for j in e if j != i]
           for i in range(len(nodes))}
    lifts = []
    for n in sorted(range(len(nodes)), key=lambda n: nodes[n, 1]):
        if bottom[n] or top[n] or not adj[n]:
            continue
        hl = {m: abs(dvec(n, m)[0]) for m in adj[n]}
        lower = [m for m in adj[n] if nodes[m, 1] <= nodes[n, 1]] or adj[n]
        # taper metric is always taken at a strut's lower end, here and in classify
        zn = nodes[n, 1]
        if any(zn - nodes[m, 1] >= fr.rise(hl[m], tan_min, nodes[m, 1]) for m in lower):
            continue
        need = min(nodes[m, 1] + fr.rise(hl[m], tan_min, nodes[m, 1]) * (1 + 1e-6) for m in lower)
        z_new = min(need, z1)
        lifts.append(z_new - nodes[n, 1])
        nodes[n, 1] = z_new
        top[n] = z_new >= z1
    stats["nodes_lifted"] = len(lifts)
    stats["lift_mm_median/max"] = [round(float(np.median(lifts)), 1) if lifts else 0,
                                   round(float(max(lifts)), 1) if lifts else 0]

    def classify(i, j):
        d = dvec(i, j)
        hl = abs(d[0])
        zm = min(nodes[i, 1], nodes[j, 1])
        if anchor(i, j):
            return "plane"
        if abs(d[1]) >= fr.rise(hl, tan_min, zm):
            return "steep"
        return "short" if fr.hor(hl, zm) <= args.max_bridge else "long"

    arch_tan = {}                                 # per-edge override after a clash

    def apex(i, j):
        d = dvec(i, j)
        p = nodes[i] + d / 2
        zt = max(nodes[i, 1], nodes[j, 1])
        p[1] = zt + fr.rise(abs(d[0]) / 2, arch_tan.get((i, j), tan_arch), zt)
        return p

    def shape(e):
        """Strut centre line in unwrapped (u, z), starting at node e[0]."""
        a = nodes[e[0]]
        b = a + dvec(*e)
        return [a, apex(*e), b] if cls[e] == "long" else [a, b]

    def segments(E):
        """All drawn strut segments with their endpoint ids (apex ids are unique < 0)."""
        A, B_, ids, keys = [], [], [], []
        for n, e in enumerate(E):
            if cls[e] == "plane":
                continue
            pts = shape(e)
            pid = [e[0], -1 - n, e[1]] if len(pts) == 3 else [e[0], e[1]]
            for k_ in range(len(pts) - 1):
                A.append(pts[k_]); B_.append(pts[k_ + 1])
                ids.append((pid[k_], pid[k_ + 1])); keys.append(e)
        return np.asarray(A), np.asarray(B_), ids, keys

    def clashes(a, b, own_ids, own_edge, S, clear):
        """True if segment a-b comes closer than `clear` to a non-adjacent strut."""
        A, B_, ids, keys = S
        keep = np.array([k_ != own_edge and not (set(p) & own_ids) for p, k_ in zip(ids, keys)])
        if not keep.any():
            return False
        for sh in (-fr.P, 0.0, fr.P):
            if np.any(seg_seg_dist(a, b, A[keep] + [sh, 0], B_[keep] + [sh, 0]) < clear):
                return True
        return False

    # Long shallow struts become arches unless the apex would leave the lattice zone or
    # the span is too wide; those stay flat bridges (both ends are supported).
    cls = {e: classify(*e) for e in E}
    n_flat = 0
    for e in E:
        if cls[e] == "long" and (apex(*e)[1] > z1 or
                                 fr.hor(abs(dvec(*e)[0]), min(nodes[e[0], 1], nodes[e[1], 1])) > args.max_arch):
            cls[e] = "short"; n_flat += 1

    removed_nodes = 0
    while True:
        adj = {i: [] for i in range(len(nodes))}
        for i, j in E:
            adj[i].append(j); adj[j].append(i)
        alive = {i for i in adj if adj[i]}
        sup = {}
        for i in sorted(alive, key=lambda n: nodes[n, 1]):
            sup[i] = bottom[i] or any(
                sup.get(j, False) and nodes[j, 1] < nodes[i, 1] - eps
                and cls[(min(i, j), max(i, j))] == "steep" for j in adj[i])
        # safety net only: after lifting, every node should be supported
        bad = {i for i in alive if not sup[i] and not top[i]}
        bad |= {i for i in alive if len(adj[i]) < 2 and not (bottom[i] or top[i])}
        bad_e = {e for e in E if cls[e] in ("short", "long") and not (sup[e[0]] and sup[e[1]])}
        if bad or bad_e:
            removed_nodes += len(bad)
            E = {e for e in E if e not in bad_e and e[0] not in bad and e[1] not in bad}
            continue
        # arches: the widest clashing one first gets a flatter apex (min-angle, no
        # margin); if that clashes too it stays a flat bridge. Then re-check.
        S = segments(E)
        clash = None
        for e in sorted((e for e in E if cls[e] == "long"), key=lambda e: -abs(dvec(*e)[0])):
            a, m, b = shape(e)
            clear = args.strut / float(fr.scale(min(nodes[e[0], 1], nodes[e[1], 1]))[0])  # chart units
            if clashes(a, m, {e[0]}, e, S, clear) or clashes(m, b, {e[1]}, e, S, clear):
                clash = e
                break
        if clash is None:
            break
        if clash not in arch_tan:
            arch_tan[clash] = tan_min
        else:
            cls[clash] = "short"; n_flat += 1

    stats["removed_unsupported_or_leaf_nodes"] = removed_nodes
    stats["shallow_kept_flat"] = n_flat
    # planarity: lifting must not make struts cross (0 expected)
    A, B_, ids, keys = segments(E)
    stats["strut_crossings"] = int(sum(
        np.sum(seg_seg_dist(A[k_], B_[k_], A[k_ + 1:] + [sh, 0], B_[k_ + 1:] + [sh, 0])[
            [not (set(ids[k_]) & set(p)) for p in ids[k_ + 1:]]] == 0)
        for k_ in range(len(A) - 1) for sh in (-fr.P, 0.0, fr.P)))
    stats["edges_final"] = len(E)
    stats["edge_classes"] = {c: sum(cls[e] == c for e in E) for c in ("plane", "steep", "short", "long")}
    lens = [np.linalg.norm(dvec(i, j)) for i, j in E if cls[(i, j)] != "plane"]
    stats["edge_len_mm"] = {"min": round(min(lens), 1), "median": round(float(np.median(lens)), 1),
                            "max": round(max(lens), 1)}

    # anchor edges on the band/rim lines lie inside those solids and are not drawn
    polys = [np.asarray(shape(e)) for e in E if cls[e] != "plane"]
    wall_metrics(args, fr, polys, stats)
    return polys, rim_gussets(args, fr, polys, stats)


def rim_gussets(args, fr, polys, stats):
    """Gaps between struts arriving at the rim bottom, as (gap midpoint u, gap width).
    Each gap gets a pointed-arch gusset in wall_pattern, so the rim grows out of the
    struts layer by layer instead of bridging. A slicer flags a bridge as a 'floating
    cantilever' once any point is > 3 mm from its support (OrcaSlicer/ElegooSlicer
    TreeSupport), i.e. any rim gap > ~6 mm. Rejected: rim posts, which shorten gaps
    but cannot get them below 6 mm without clashing with the lattice."""
    zr = args.height - args.rim
    us = sorted((a[0] + (zr - a[1]) / (b[1] - a[1]) * (b[0] - a[0])) % fr.P
                for p in polys for a, b in zip(p[:-1], p[1:]) if (a[1] - zr) * (b[1] - zr) < 0)
    gaps = [((us[(k + 1) % len(us)] - u) % fr.P, u) for k, u in enumerate(us)]
    stats["rim_gaps_mm_p50/max"] = [round(float(np.median([g for g, _ in gaps])), 1),
                                    round(float(max(g for g, _ in gaps)), 1)]
    return [(u + g / 2, g) for g, u in gaps if fr.hor(g, zr) > args.strut]


def wall_metrics(args, fr, polys, stats, res=0.5):
    """Orientation mix and opening sizes of the unrolled wall."""
    r = args.strut / 2
    zlo, zhi = args.band, args.height - args.rim
    nu, nz = int(fr.P / res), int((zhi - zlo) / res)
    U = (np.arange(nu) + 0.5) * res
    Z = zlo + (np.arange(nz) + 0.5) * res
    # strut length per angle from vertical
    hist = np.zeros(4)
    for p in polys:
        for d in np.diff(p, axis=0):
            ang = math.degrees(math.atan2(abs(d[0]), abs(d[1])))
            hist[min(int(ang // 22.5), 3)] += math.hypot(d[0], d[1])
    solid = np.zeros((nz, nu), bool)
    for p in polys:
        for a, b in zip(p[:-1], p[1:]):
            for shift in (-fr.P, 0.0, fr.P):
                a2, b2 = a + [shift, 0], b + [shift, 0]
                lo_ = np.minimum(a2, b2) - r; hi_ = np.maximum(a2, b2) + r
                if hi_[0] < 0 or lo_[0] > fr.P:
                    continue
                iu = slice(max(int(lo_[0] / res), 0), min(int(hi_[0] / res) + 2, nu))
                iz = slice(max(int((lo_[1] - zlo) / res), 0), min(int((hi_[1] - zlo) / res) + 2, nz))
                uu, zz = U[iu][None, :], Z[iz][:, None]
                ab = b2 - a2
                t = np.clip(((uu - a2[0]) * ab[0] + (zz - a2[1]) * ab[1]) / max(ab @ ab, 1e-9), 0, 1)
                solid[iz, iu] |= np.hypot(uu - a2[0] - t * ab[0], zz - a2[1] - t * ab[1]) <= r
    empty = ~solid
    lab, _ = ndimage.label(empty)
    for row in range(nz):                         # wall is closed around: join the seam
        a_, b_ = lab[row, 0], lab[row, -1]
        if a_ and b_ and a_ != b_:
            lab[lab == b_] = a_
    edt = ndimage.distance_transform_edt(empty) * res
    ids = [i for i in np.unique(lab) if i]
    incl = np.asarray(ndimage.maximum(edt, lab, ids)) * 2
    objs = ndimage.find_objects(lab)
    heights = np.array([(objs[i - 1][0].stop - objs[i - 1][0].start) * res for i in ids])
    widths = np.array([(objs[i - 1][1].stop - objs[i - 1][1].start) * res for i in ids])
    big = incl >= 2.0                             # ignore slivers between touching struts
    stats["wall_metrics"] = {
        "struts": len(polys),
        "len_share_deg_0-22/22-45/45-67/67-90": [round(x, 2) for x in hist / max(hist.sum(), 1e-9)],
        "openings": int(big.sum()),
        "open_area_share": round(float(empty.mean()), 2),
        "inscribed_dia_mm_p50/p90/max": [round(float(np.percentile(incl[big], q)), 1) for q in (50, 90, 100)],
        "height_mm_p50/p90/max": [round(float(np.percentile(heights[big], q)), 1) for q in (50, 90, 100)],
        "height_to_width_p90": round(float(np.percentile(heights[big] / np.maximum(widths[big], res), 90)), 1)}


def floor_segments(args, fr, rng):
    h = fr.extent
    pts = poisson(rng, [-h - 15, -h - 15], [h + 15, h + 15], args.floor_cell, tries=40000)
    vor = Voronoi(pts)
    segs = []
    for rv in vor.ridge_vertices:
        if -1 in rv:
            continue
        a, b = vor.vertices[rv[0]], vor.vertices[rv[1]]
        if fr.footprint_sdf(*((a + b) / 2)) < 0:
            segs.append((a, b))
    return segs


# ---------- SDF ----------

class Grid:
    def __init__(self, lo, hi, h):
        self.h = h
        self.lo = np.asarray(lo, float)
        self.n = np.ceil((np.asarray(hi) - self.lo) / h).astype(int) + 1
        self.ax = [self.lo[i] + h * np.arange(self.n[i]) for i in range(3)]
        self.D = np.empty(self.n, np.float32)

    def box(self, pmin, pmax):
        i0 = np.clip(np.floor((np.asarray(pmin) - self.lo) / self.h).astype(int), 0, self.n - 1)
        i1 = np.clip(np.ceil((np.asarray(pmax) - self.lo) / self.h).astype(int) + 1, 1, self.n)
        sl = tuple(slice(a, b) for a, b in zip(i0, i1))
        X = self.ax[0][sl[0]][:, None, None]
        Y = self.ax[1][sl[1]][None, :, None]
        Z = self.ax[2][sl[2]][None, None, :]
        return sl, X, Y, Z


# CTX: 2D field F[z_row, u_col] of the whole wall; rows are exactly the voxel z axis, so
#      build_sdf reads F[iz] per slice and interpolates over u with period P.
# Used-by: orchid_pot.py::build_sdf | Depends: orchid_pot.py::rim_gussets (gap list contract)
# Invariant: gussets grow from S = F <= 0, so all struts must be in F before the growth step.
def wall_pattern(args, fr, polys, gussets, Z, ez):
    """2D smooth-min SDF of wall struts, rim gussets and eyelet rings on the unrolled
    wall, rows = Z."""
    nu = int(fr.P / 0.2)
    U = (np.arange(nu) + 0.5) * fr.P / nu
    du = fr.P / nu
    dz = Z[1] - Z[0]
    F = np.full((len(Z), nu), 10.0, np.float32)

    def add(lo, hi, dist):
        # a grid column U stands for the unwrapped coordinate U + sh
        for sh in (-fr.P, 0.0, fr.P):
            iu0 = max(int(np.floor((lo[0] - sh - U[0]) / du)), 0)
            iu1 = min(int(np.ceil((hi[0] - sh - U[0]) / du)) + 1, nu)
            iz0 = max(int(np.floor((lo[1] - Z[0]) / dz)), 0)
            iz1 = min(int(np.ceil((hi[1] - Z[0]) / dz)) + 1, len(Z))
            if iu1 <= iu0 or iz1 <= iz0:
                continue
            Q = np.stack(np.broadcast_arrays((U[iu0:iu1] + sh)[None, :], Z[iz0:iz1][:, None]), -1)
            # chart -> local real metric of each row (all ones for vertical walls); the
            # dist callbacks scale their own geometry with the same factors M
            ku, kz = fr.scale(Z[iz0:iz1])
            M = np.stack(np.broadcast_arrays(ku[:, None], np.full((iz1 - iz0, 1), kz)), -1)
            sl = (slice(iz0, iz1), slice(iu0, iu1))
            F[sl] = smin(F[sl], dist(Q * M, M), args.fillet)

    pad = (args.strut / 2 + args.fillet + 1) / fr.ku_min
    hw = args.strut / 2
    for p in polys:
        add(p.min(0) - pad, p.max(0) + pad,
            lambda Q, M, p=p: np.min([pt_seg_dist(Q, a * M, b * M) for a, b in zip(p[:-1], p[1:])], 0) - hw)
    # Rim gussets, grown row by row from the actual strut material. Each row may extend
    # the row below by at most dz / tan(min-angle) sideways (a printable overhang), and
    # only inside a pointed-arch envelope per gap (flanks at min-angle + 5 deg, apex 1 mm
    # inside the rim). Growing from real material means a gusset can never start in
    # mid-air. Rejected: analytic wedges bounded by the arriving strut's line; where the
    # strut bends, the extended line left a floating tip. Also rejected: the same gussets
    # under flat wall struts. The support strut below their end node often leans away,
    # so there is nothing to grow from; the partial gussets only moved the warnings.
    zr = args.height - args.rim
    te = math.tan(math.radians(args.min_angle + 5))
    ts = math.tan(math.radians(args.min_angle))
    ku_r, kz = fr.scale(zr)                        # gussets stay within ~10 mm of the rim
    ku_r = float(ku_r)
    env = np.zeros(F.shape, bool)

    def envelope(um, half, z_line):
        """Allow rows z with z_line(u) - te*|u - um| <= z <= z_line(u) over the span."""
        d_ = (U - um + fr.P / 2) % fr.P - fr.P / 2
        cols = np.abs(d_) < half
        zl = z_line(um + d_[cols])
        env[:, cols] |= (Z[:, None] >= zl - te * (np.abs(d_[cols]) * ku_r)) & (Z[:, None] <= zl)

    for um, g in gussets:
        envelope(um, g / 2, lambda u: np.full_like(u, zr + 1.0))
    # exact growth distance per row; a whole-pixel dilation (0.4 instead of 0.57 mm)
    # grows slower than the envelope flanks widen, and the gussets never close
    grow = dz / ts / ku_r
    S = F <= 0
    G = S.copy()
    for k_ in range(1, len(Z)):
        if env[k_].any() and G[k_ - 1].any():
            d1 = ndimage.distance_transform_edt(~np.tile(G[k_ - 1], 3))[nu:2 * nu] * du
            G[k_] |= (d1 <= grow) & env[k_]
    added = G & ~S
    rows = np.where(added.any(1))[0]
    # unfilled run directly under the rim: remaining bridge length (0 = fully gusseted)
    top_row = G[np.searchsorted(Z, zr) - 1]
    runs, cur = [], 0
    for v in np.concatenate([top_row, top_row]):
        cur = 0 if v else cur + 1; runs.append(cur)
    stats_g = {"rim_gussets": len(gussets),
               "gusset_area_cm2": round(float(added.sum()) * du * dz / 100, 1),
               "rim_unfilled_max_mm": round(min(max(runs), nu) * du, 1)}
    if len(rows):
        # signed distance of the grown mask (periodic in u), merged into F
        r0, r1 = max(rows[0] - 8, 0), min(rows[-1] + 9, len(Z))
        T = np.concatenate([G[r0:r1]] * 3, axis=1)
        smp = (dz * kz, du * ku_r)
        sd = (ndimage.distance_transform_edt(~T, sampling=smp)
              - ndimage.distance_transform_edt(T, sampling=smp))[:, nu:2 * nu]
        F[r0:r1] = np.minimum(F[r0:r1], sd)

    rt = args.eyelet_wire / 2
    Rt = args.eyelet_hole / 2 + rt
    for uc in fr.eyelet_u(args.eyelets):
        c = np.array([uc, ez])
        add(c - (Rt + rt + 1) / fr.ku_min, c + (Rt + rt + 1) / fr.ku_min,
            lambda Q, M, c=c: np.abs(np.linalg.norm(Q - c * M, axis=-1) - Rt) - rt)
    return U, F, stats_g


def build_sdf(args, fr, polys, gussets, fsegs, stats):
    half = fr.extent
    H = args.height
    k = args.fillet
    rt = args.eyelet_wire / 2
    Rt = args.eyelet_hole / 2 + rt
    ez = H + 0.5 + args.eyelet_hole / 2          # hole bottom 0.5 mm above the rim top
    m = 2.0
    # odd origin offset: no grid sample lands exactly on a flat face (z=0, band top, ...),
    # which would make marching cubes emit zero-area slivers
    o = 0.137 * args.voxel
    ztop = ez + Rt + rt + 1.0 if args.eyelets else H + 1.0
    g = Grid([-half - m + o, -half - m + o, -1.0 + o], [half + m, half + m, ztop], args.voxel)
    U, F, stats_g = wall_pattern(args, fr, polys, gussets, g.ax[2], ez)
    stats.update(stats_g)

    X2 = g.ax[0][:, None]; Y2 = g.ax[1][None, :]
    dfoot = fr.footprint_sdf(X2, Y2)
    uxy = fr.unmap_u(*np.broadcast_arrays(X2, Y2))
    fl = np.full(dfoot.shape, 10.0)
    fw = args.floor_strut / 2
    for a, b in fsegs:
        lo_ = np.minimum(a, b) - fw - 1; hi_ = np.maximum(a, b) + fw + 1
        sl, X, Y, _ = g.box([lo_[0], lo_[1], 0], [hi_[0], hi_[1], 0])
        ab = b - a
        px, py = X[:, :, 0] - a[0], Y[:, :, 0] - a[1]
        t = np.clip((px * ab[0] + py * ab[1]) / ab.dot(ab), 0, 1)
        fl[sl[:2]] = np.minimum(fl[sl[:2]], np.hypot(px - t * ab[0], py - t * ab[1]) - fw)
    fl = np.maximum(fl, dfoot + 1.0)

    zc_b = args.band / 2 - 1
    for iz, z in enumerate(g.ax[2]):
        dpat = np.interp(uxy, U, F[iz], period=fr.P)
        _, pw = fr.wall_fields(X2, Y2, z)
        wall = round_isect(dpat, pw - fr.t / 2, args.edge_round)
        band = rounded_box2(pw, abs(z - zc_b), fr.t / 2, args.band / 2 + 1, 1.5)
        rim = rounded_box2(pw, abs(z - (H - args.rim / 2)), fr.t / 2, args.rim / 2, 2.0)
        flo = np.maximum(fl, abs(z - (args.floor / 2 - 0.5)) - (args.floor / 2 + 0.5))
        solid = smin(smin(band, rim, 2.0), flo, 2.0)
        g.D[:, :, iz] = np.maximum(smin(wall, solid, k), -z)
    return g


# ---------- mesh + checks ----------

def clean(mesh):
    # marching cubes and quadric decimation both leave a few zero-area faces, which
    # show up as non-manifold edges / separate bodies
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.remove_unreferenced_vertices()
    mesh.merge_vertices()
    return mesh


def to_mesh(g, stats):
    v, f, _, _ = marching_cubes(g.D, level=0.0, spacing=(g.h,) * 3)
    mesh = clean(trimesh.Trimesh(v + g.lo, f, process=True))
    # Sub-mm^3 shells are voids enclosed in the material (seen with tapered walls: smooth-min
    # pockets between near-touching fields); dropping their shell fills them. They also
    # made the decimated mesh non-watertight.
    parts = mesh.split(only_watertight=False)
    if len(parts) > 1:
        keep = [q for q in parts if abs(q.volume) >= 1.0]
        stats["micro_voids_filled"] = len(parts) - len(keep)
        mesh = trimesh.util.concatenate(keep)
    if mesh.volume < 0:
        mesh.invert()
    return mesh


def check_mesh(mesh, args, stats):
    stats["triangles"] = int(len(mesh.faces))
    stats["watertight"] = bool(mesh.is_watertight)
    stats["winding_consistent"] = bool(mesh.is_winding_consistent)
    stats["volume_cm3"] = round(mesh.volume / 1000, 1)
    stats["mass_g_petg_1.27"] = round(mesh.volume / 1000 * 1.27, 0)
    stats["bbox_mm"] = [round(x, 1) for x in mesh.extents]
    stats["bodies"] = int(len(mesh.split(only_watertight=False)))
    # downward faces steeper than the allowed overhang, excluding the bed face
    cz = mesh.triangles_center[:, 2]
    nz = mesh.face_normals[:, 2]
    lim = -math.cos(math.radians(args.min_angle))
    ov = (nz < lim) & (cz > 0.3)
    rimz = args.height - args.rim
    stats["overhang_area_cm2"] = {
        "rim_underside": round(mesh.area_faces[ov & (np.abs(cz - rimz) < 1.5)].sum() / 100, 1),
        "elsewhere": round(mesh.area_faces[ov & (np.abs(cz - rimz) >= 1.5)].sum() / 100, 1),
        "total_surface": round(mesh.area / 100, 1)}
    # layer islands: polygons with no overlap to the layer below need support.
    # Specks smaller than one extrusion footprint (0.4 x 0.4 mm) are not printed by a
    # slicer; they are dropped from the support union as well, so an island that starts
    # as a speck is still caught on the next layer.
    zs = np.arange(args.layer / 2, mesh.bounds[1, 2], args.layer)
    secs = mesh.section_multiplane([0, 0, 0], [0, 0, 1], zs)
    prev, islands = None, []
    for z, s in zip(zs, secs):
        polys = [p for p in s.polygons_full if p.area >= 0.16] if s is not None else []
        if prev is not None:
            for p in polys:
                if not p.intersects(prev):
                    c = p.centroid
                    islands.append([round(c.x, 1), round(c.y, 1), round(float(z), 1)])
        prev = shapely.make_valid(unary_union(polys)) if polys else None
    stats["unsupported_islands"] = len(islands)
    stats["island_positions"] = islands[:20]
    stats["cantilevers"] = cantilevers(mesh, args, zs, secs)


def cantilevers(mesh, args, zs, secs, ew=0.42, thresh=30.0, limit=3.0):
    """Approximate replica of OrcaSlicer's cantilever test (TreeSupport detect_overhangs).
    Overhang = layer minus the layer below grown by layer/tan(threshold). It is a
    cantilever if a contour point lies > 3 mm from where the overhang touches the layer
    below grown by extrusion width + 0.1 mm. Defaults assumed: 30 deg threshold, 0.42 mm
    line width."""
    grow = args.layer / math.tan(math.radians(thresh))
    hits, prev = [], None
    for z, s in zip(zs, secs):
        cur = unary_union([p for p in s.polygons_full if p.area >= 0.16]) if s is not None else None
        if prev is not None and cur is not None:
            ovh = cur.difference(prev.buffer(grow))
            base = prev.buffer(ew + 0.1)
            for p in getattr(ovh, "geoms", [ovh]):
                if p.is_empty or p.area < 0.16:
                    continue
                contact = p.intersection(base)
                if contact.is_empty:
                    continue
                d = max(contact.distance(shapely.Point(c)) for c in p.exterior.coords)
                if d > limit:
                    c = p.centroid
                    hits.append((round(float(z), 1), round(d, 1), round(c.x, 1), round(c.y, 1)))
        prev = shapely.make_valid(cur) if cur is not None else None
    rim0 = args.height - args.rim
    zone = lambda z: "rim_underside" if abs(z - rim0) < 1.0 else ("eyelets" if z > args.height else "wall")
    out = {}
    for z, d, x, y in hits:
        o = out.setdefault(zone(z), {"regions": 0, "max_reach_mm": 0.0, "layers": set()})
        o["regions"] += 1; o["max_reach_mm"] = max(o["max_reach_mm"], d); o["layers"].add(z)
    for o in out.values():
        o["layers"] = sorted(o["layers"])[:12]
    return out


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    fr = make_frame(args)
    stats = {"params": vars(args), "shape": {"type": args.shape, **fr.describe()}}
    t0 = time.time()
    polys, gussets = wall_graph(args, fr, rng, stats)
    fsegs = floor_segments(args, fr, rng)
    print(f"graph {time.time() - t0:.1f}s: {len(polys)} wall struts, {len(fsegs)} floor struts")
    if args.graph_only:
        print(json.dumps({k: v for k, v in stats.items() if k != "params"}, indent=1))
        return
    g = build_sdf(args, fr, polys, gussets, fsegs, stats)
    print(f"sdf {time.time() - t0:.1f}s grid {g.n.tolist()}")
    mesh = to_mesh(g, stats)
    del g
    print(f"mesh {time.time() - t0:.1f}s {len(mesh.faces)} tris")
    # quadric decimation occasionally leaves a non-manifold spot; retry milder once
    for ratio in (args.decimate, args.decimate * 0.8) if args.decimate > 0 else ():
        v, f = fast_simplification.simplify(mesh.vertices, mesh.faces, ratio)
        dec = clean(trimesh.Trimesh(v, f, process=True))
        if dec.is_watertight:
            mesh = dec
            stats["decimation_ratio"] = ratio
            break
    else:
        stats["decimation_skipped"] = "result not watertight or disabled"
    check_mesh(mesh, args, stats)
    mesh.export(args.out + ".stl")
    with open(args.out + "_report.json", "w") as fh:
        json.dump(stats, fh, indent=1)
    print(f"done {time.time() - t0:.1f}s")
    print(json.dumps({k: v for k, v in stats.items() if k != "params"}, indent=1))


if __name__ == "__main__":
    main()
