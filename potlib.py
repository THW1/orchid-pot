# SPDX-License-Identifier: MIT
# === CONTEXT ===
# Purpose:     Shape-independent helpers shared by all shapes: SDF primitives, Poisson
#              seeding, segment distances, the voxel Grid, floor pattern, wall metrics,
#              meshing (marching cubes, cleanup, decimation) and all report checks.
# Decisions:   Split out of orchid_pot.py (v0.5.0, verbatim move) so orb.py can reuse the
#              helpers without importing or touching the box/cylinder pipeline. Later
#              additions shared by both pipelines: grade_fraction, floor_field,
#              groove_geometry/groove_sdf (v0.9.0, the snap groove shared by pots and clip
#              tray), require_rim_crossings (v0.6.1), finish(extra_checks=...) and the STL round
#              trip check of every decimation candidate (v0.7.1, stl_roundtrip /
#              export_safe).
#              Behaviour changes here affect every shape: run tests/ (byte-identical
#              goldens for box and cylinder).
# Related:     Used by orchid_pot.py (box, cylinder) and orb.py.
# === END CONTEXT ===
import hashlib
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
        with np.errstate(divide="ignore", invalid="ignore"):   # zero-volume slivers
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


def drain_holes_sdf(X, Y, args, r_inner):
    """2D SDF of the drainage holes of a solid floor (negative inside a hole): one centre
    hole plus the rest evenly on a ring at 45 % of the inner floor radius r_inner."""
    n, r = args.drain_holes, args.drain_dia / 2
    pts = [(0.0, 0.0)] if n >= 1 else []
    pts += [(0.45 * r_inner * math.cos(2 * math.pi * k / (n - 1)),
             0.45 * r_inner * math.sin(2 * math.pi * k / (n - 1))) for k in range(n - 1)]
    d = np.full(np.broadcast(X, Y).shape, 1e3)
    for px, py in pts:
        d = np.minimum(d, np.hypot(X - px, Y - py) - r)
    return d


def membrane_top(args, H, stats):
    """Top height of the cell membrane (--membrane), or None when off. Default: into the
    rim (closed up to the top); --membrane-top lower leaves the cells above it open."""
    if not getattr(args, "membrane", 0):
        return None
    if args.membrane_top is None:
        zt = H - args.rim / 2
    elif args.band < args.membrane_top <= H - args.rim:
        zt = args.membrane_top
    else:
        raise SystemExit(f"--membrane-top must be in (--band, rim bottom] = ({args.band}, {H - args.rim:.1f}] mm")
    stats["membrane"] = {"thickness_mm": args.membrane, "top_mm": round(zt, 1),
                         "open_below_rim_mm": round(max(H - args.rim - zt, 0.0), 1)}
    return zt


def membrane_wall(wall, d_in, z, args, zt, k):
    """Lattice wall SDF plus the cell membrane: a slab of --membrane thickness on the inner
    face (d_in = distance from the inner face towards the outside), up to height zt.
    The membrane is a half-space smooth-joined to the struts (fillet k at the strut sides)
    and only then cut flat at the inner face: a smooth union of the coinciding inner faces
    would bulge them inwards by k/4 (cf. the solid-wall band ridge)."""
    return np.maximum(smin(wall, np.maximum(d_in - args.membrane, z - zt), k), -d_in)


VENT_MIN_PILLAR = 4.0                      # wall left between two vent slots [mm]


def vent_geometry(args, H, L_min, stats):
    """Vertical ventilation slots under the rim of a solid pot (--vent-slots), or None.
    Returns (z_bottom, z_apex). The apex sits at the rim bottom, so the full rim stays a
    closed ring. L_min: smallest mid-surface perimeter over the slot height (pillar check)."""
    if not args.vent_slots or args.wall_style != "solid":
        return None
    n, w, h = args.vent_slots, args.vent_width, args.vent_height
    zt = H - args.rim
    zb = zt - h
    if zb < args.band + 2.0:
        raise SystemExit(f"--vent-height {h}: slots would reach into the bottom band "
                         f"(max {zt - args.band - 2.0:.1f} mm)")
    pillar = L_min / n - w
    if pillar < VENT_MIN_PILLAR:
        raise SystemExit(f"--vent-slots {n}: only {pillar:.1f} mm wall between slots "
                         f"(need >= {VENT_MIN_PILLAR}); max {int(L_min // (w + VENT_MIN_PILLAR))} slots")
    stats["vent_slots"] = {"n": n, "width/height_mm": [w, h], "z_bottom/apex_mm": [round(zb, 1), round(zt, 1)],
                           "pillar_min_mm": round(pillar, 1),
                           "open_area_cm2": round(n * (w * (h - w / 2) + w * w / 4) / 100, 1)}
    return zb, zt


def vent_slots_sdf(frac, L, z, args, zt):
    """2D SDF (negative inside) of the vent slots at perimeter fraction frac in [0, 1)
    (0 = first eyelet position; slot k is centred at (k + 0.5)/n, i.e. between eyelets),
    perimeter L and height z. Flat bottom, 45 deg pointed roof with its apex at zt: prints
    without support. Applied as a horizontal cut through the wall."""
    n, w = args.vent_slots, args.vent_width / 2
    q = frac * n
    du = np.abs(q - np.floor(q) - 0.5) * L / n
    zb = zt - args.vent_height
    return np.maximum(np.maximum(du - w, zb - z), (z + du - zt) / math.sqrt(2))


def groove_geometry(args):
    """Continuous snap groove in the bottom band's outer face (--clip-groove), or None.
    Returns (z0, depth, opening) in pot coordinates: flat groove floor at z0 (the clip's
    catch face; it faces up, so it prints as a floor), opening height at the outer face,
    ceiling at 45 deg (self-supporting, and the clip's lead-in). >= 2.5 mm of solid band
    stay above the groove. The floor is snapped down to a layer boundary (multiple of
    --layer): the catch face is then a crisp layer surface, and it never lies on a slice
    plane of the island check (z0 = 2.5 did, and the slice fell apart into fake islands).
    Shared by the pots (cut) and the clip tray (hooks), so both agree."""
    if not getattr(args, "clip_groove", False):
        return None
    z0 = math.floor((args.band - args.groove_height - 2.5) / args.layer + 1e-9) * args.layer
    return round(z0, 6), args.groove_depth, args.groove_height


def groove_sdf(n, z, gg):
    """SDF of the groove region; n = horizontal distance inward from the band's outer face
    (negative outside), z = pot height. Negative inside the groove."""
    z0, d, h = gg
    return np.maximum(np.maximum(n - d, z0 - z), ((z - z0) + np.maximum(n, 0) - h) / math.sqrt(2))


# Parameters that set the pot's outline, band and band-top height (what a tray sits on).
SHAPE_KEYS = ("shape", "size", "corner_r", "diameter", "taper", "height", "stretch",
              "min_wall_angle", "wall", "band", "rim", "floor", "fillet")


# CTX: Pot identity for accessories. 'pattern' hashes the strut paths, not parameters, so it
#      also catches library-version changes. SHAPE_KEYS must list every parameter that moves
#      the outline, band or band top, or a mismatched tray would pass.
# Used-by: orchid_pot.py::build_pot, orb.py::run, accessories.py::run
def fingerprint(args, struts=None):
    """Identity of a pot as seen by its accessories: 'shape' hashes SHAPE_KEYS (orb: without
    the derived height); 'pattern' hashes the strut paths themselves (rounded to 1 um), so
    any parameter or library change that moves a strut changes it."""
    d = vars(args)
    keys = [k for k in SHAPE_KEYS if k in d and not (k == "height" and d.get("shape") == "orb")]
    if d.get("wall_style") == "solid":        # only then: Voronoi pots keep their hash
        keys += ["wall_style"]
    if d.get("clip_groove"):                  # only then: older/groove-less pots keep their hash
        keys += ["clip_groove", "groove_depth", "groove_height", "layer"]   # layer: groove floor
        if (d.get("band_wall") or 0) > d["wall"]:   # only then: thicker-band pots differ
            keys += ["band_wall"]
    out = {"shape": hashlib.sha256(json.dumps({k: d[k] for k in keys}, sort_keys=True).encode()).hexdigest()[:12]}
    if struts is not None:
        q = np.round(np.vstack([np.asarray(p, float).reshape(len(p), -1)[:, :3] for p in struts]) * 1000).astype(np.int64)
        out["pattern"] = hashlib.sha256(q.tobytes()).hexdigest()[:12]
    return out


def stl_roundtrip(mesh):
    """(ok, reloaded): the mesh after binary STL export + reload, and whether it is still
    watertight."""
    back = trimesh.load(io.BytesIO(mesh.export(file_type="stl")), file_type="stl")
    return bool(back.is_watertight), back


def fill_voids(mesh, stats):
    """Drop shells oriented against the main body: enclosed voids. Decimation can close a
    sub-voxel channel into a void (orb --membrane: 1.3 mm^3 where rim gussets meet the
    membrane). Zero-volume slivers are not voids; export_safe handles them. A single-shell
    mesh is returned untouched (outputs stay byte-identical)."""
    parts = mesh.split(only_watertight=False)
    if len(parts) < 2:
        return mesh
    with np.errstate(divide="ignore", invalid="ignore"):       # zero-volume slivers
        vols = [float(q.volume) for q in parts]
    s = math.copysign(1.0, max(vols, key=abs))
    keep = [q for q, v in zip(parts, vols) if v * s > -1e-6]
    if len(keep) == len(parts):
        return mesh
    stats["voids_filled"] = len(parts) - len(keep)
    return clean(trimesh.util.concatenate(keep))


def export_safe(mesh, stats):
    """The mesh as it will be in the STL file. If the float32 round trip breaks it (vertices
    that merge in float32 split off zero-volume two-triangle slivers; seen on box --seed 3),
    use the round-tripped mesh without those slivers. Meshes that pass stay untouched, so
    outputs that were fine do not change. None if no watertight STL can be had."""
    ok, back = stl_roundtrip(mesh)
    if ok:
        return mesh
    parts = back.split(only_watertight=False)
    with np.errstate(divide="ignore", invalid="ignore"):       # the slivers have zero volume
        keep = [q for q in parts if abs(q.volume) >= 1.0]
    fixed = clean(trimesh.util.concatenate(keep))
    if fixed.is_watertight and stl_roundtrip(fixed)[0]:
        stats["stl_slivers_removed"] = stats.get("stl_slivers_removed", 0) + len(parts) - len(keep)
        return fixed
    return None


# CTX: Common tail of every pipeline: mesh -> decimate (one milder retry) -> checks ->
#      write. extra_checks runs after check_mesh and before writing, so its stats land in
#      the report.
# Used-by: orchid_pot.py::main, orb.py::run
# Invariant: g.D is released before decimation; callers must not use the grid afterwards.
def finish(g, args, stats, t0, extra_checks=None):
    """Mesh the SDF grid, decimate, check, write <out>.stl and <out>_report.json.
    extra_checks(mesh, stats): optional shape-specific checks (orb: stability)."""
    mesh = to_mesh(g, stats)
    g.D = None                                 # free the voxel grid (~180 MB) before decimating
    print(f"mesh {time.time() - t0:.1f}s {len(mesh.faces)} tris")
    # quadric decimation occasionally leaves a non-manifold spot; retry milder once. The
    # candidate must also survive the STL round trip: float32 export can collapse a sliver
    # that is fine in memory (seen on a clip tray: watertight in memory, not after reload).
    cand = None
    for ratio in (args.decimate, args.decimate * 0.8) if args.decimate > 0 else ():
        v, f = fast_simplification.simplify(mesh.vertices, mesh.faces, ratio)
        st = {}                                # stats of a rejected candidate are dropped
        dec = fill_voids(clean(trimesh.Trimesh(v, f, process=True)), st)
        cand = export_safe(dec, st) if dec.is_watertight else None
        if cand is not None:
            stats.update(st)
            mesh = cand
            stats["decimation_ratio"] = ratio
            break
    else:
        stats["decimation_skipped"] = "result not watertight (after STL round trip) or disabled"
        cand = export_safe(mesh, stats)
        mesh = cand if cand is not None else mesh
    # export_safe verified every accepted mesh; only a total failure leaves False here
    stats["stl_roundtrip_watertight"] = cand is not None
    check_mesh(mesh, args, stats)
    if extra_checks:
        extra_checks(mesh, stats)
    mesh.export(args.out + ".stl")
    with open(args.out + "_report.json", "w") as fh:
        json.dump(stats, fh, indent=1)
    print(f"done {time.time() - t0:.1f}s")
    print(json.dumps({k: v for k, v in stats.items() if k != "params"}, indent=1))
