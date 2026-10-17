# SPDX-License-Identifier: MIT
# === CONTEXT ===
# Purpose:     Shape-independent helpers shared by all shapes: SDF primitives, Poisson
#              seeding, segment distances, the voxel Grid, floor pattern, wall metrics,
#              meshing (marching cubes, cleanup, decimation) and all report checks.
# Decisions:   Split out of orchid_pot.py (v0.5.0, verbatim move) so orb.py can reuse the
#              helpers without importing or touching the box/cylinder pipeline. Later
#              additions shared by both pipelines: grade_fraction, floor_field,
#              require_rim_crossings (v0.6.1), finish(extra_checks=...) and the STL round
#              trip check of every decimation candidate (v0.7.1, stl_roundtrip_ok).
#              Behaviour changes here affect every shape: run tests/ (byte-identical
#              goldens for box and cylinder).
# Related:     Used by orchid_pot.py (box, cylinder) and orb.py.
# === END CONTEXT ===
import io
import json
import math
import time

import numpy as np
import shapely
import trimesh
import fast_simplification
from scipy import ndimage
from scipy.spatial import Voronoi
from shapely.ops import unary_union
from skimage.measure import marching_cubes


# ---------- SDF primitives and sampling ----------

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


# ---------- metrics and floor ----------

def grade_fraction(z, zg, z1):
    """Cell grading weight: 0 below zg, rising linearly to 1 at z1. If the range is empty
    (e.g. --grade-top 0) grading is off; the bare formula would invert and shrink every cell."""
    return min(max((z - zg) / (z1 - zg), 0.0), 1.0) if z1 > zg else 0.0


def require_rim_crossings(n):
    """Rim gaps need >= 2 struts crossing the rim line (0: max() of nothing; 1: a full
    turn would read as a 0 mm gap)."""
    if n < 2:
        raise SystemExit("error: fewer than 2 struts reach the rim line; "
                         "try a smaller --cell or another --seed")


def floor_field(g, fsegs, dfoot, args):
    """2D SDF of the floor lattice, clipped 1 mm inside the footprint (dfoot = SDF of the
    footprint on the grid's xy plane)."""
    fl = np.full(dfoot.shape, 10.0)
    fw = args.floor_strut / 2
    for a, b in fsegs:
        lo_ = np.minimum(a, b) - fw - 1; hi_ = np.maximum(a, b) + fw + 1
        sl, X, Y, _ = g.box([lo_[0], lo_[1], 0], [hi_[0], hi_[1], 0])
        ab = b - a
        px, py = X[:, :, 0] - a[0], Y[:, :, 0] - a[1]
        t = np.clip((px * ab[0] + py * ab[1]) / ab.dot(ab), 0, 1)
        fl[sl[:2]] = np.minimum(fl[sl[:2]], np.hypot(px - t * ab[0], py - t * ab[1]) - fw)
    return np.maximum(fl, dfoot + 1.0)

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


# ---------- voxel grid ----------

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


# CTX: Common tail of every pipeline: mesh -> decimate (one milder retry) -> checks ->
#      write. extra_checks runs after check_mesh and before writing, so its stats land in
#      the report.
# Used-by: orchid_pot.py::main, orb.py::run
# Invariant: g.D is released before decimation; callers must not use the grid afterwards.
def stl_roundtrip_ok(mesh):
    """True if the mesh is still watertight after binary STL export + reload."""
    back = trimesh.load(io.BytesIO(mesh.export(file_type="stl")), file_type="stl")
    return bool(back.is_watertight)


def finish(g, args, stats, t0, extra_checks=None):
    """Mesh the SDF grid, decimate, check, write <out>.stl and <out>_report.json.
    extra_checks(mesh, stats): optional shape-specific checks (orb: stability)."""
    mesh = to_mesh(g, stats)
    g.D = None                                 # free the voxel grid (~180 MB) before decimating
    print(f"mesh {time.time() - t0:.1f}s {len(mesh.faces)} tris")
    # quadric decimation occasionally leaves a non-manifold spot; retry milder once. The
    # candidate must also survive the STL round trip: float32 export can collapse a sliver
    # that is fine in memory (seen on a clip tray: watertight in memory, not after reload).
    for ratio in (args.decimate, args.decimate * 0.8) if args.decimate > 0 else ():
        v, f = fast_simplification.simplify(mesh.vertices, mesh.faces, ratio)
        dec = clean(trimesh.Trimesh(v, f, process=True))
        if dec.is_watertight and stl_roundtrip_ok(dec):
            mesh = dec
            stats["decimation_ratio"] = ratio
            break
    else:
        stats["decimation_skipped"] = "result not watertight (after STL round trip) or disabled"
    stats["stl_roundtrip_watertight"] = stl_roundtrip_ok(mesh)
    check_mesh(mesh, args, stats)
    if extra_checks:
        extra_checks(mesh, stats)
    mesh.export(args.out + ".stl")
    with open(args.out + "_report.json", "w") as fh:
        json.dump(stats, fh, indent=1)
    print(f"done {time.time() - t0:.1f}s")
    print(json.dumps({k: v for k, v in stats.items() if k != "params"}, indent=1))
