# SPDX-License-Identifier: MIT
# === CONTEXT ===
# Purpose:     EXPERIMENTAL `orb` shape: an "onion", i.e. a vertically stretched ellipsoid
#              with a flat floor and a top opening, Voronoi wall, optional eyelets. Called
#              lazily from orchid_pot.py::main.
# Decisions:   Own module, so the box/cylinder pipeline stays untouched (its frames assume
#              an unrollable wall; this wall is curved in two directions). Shared helpers
#              come from potlib.py only.
#              Floor and opening are cut where the wall angle reaches --min-wall-angle.
#              Below that, almost no strut direction can reach the 35 deg printing limit
#              (allowed spread at a 55 deg wall: +-45.6 deg around up-the-wall).
#              Pattern: SphericalVoronoi of Poisson seeds on a sphere, mapped to the
#              ellipsoid by stretching z. The stretch elongates cells vertically, like
#              --zstretch for box/cylinder. Seeds are mirrored in polar angle at the
#              band/rim anchor latitudes, so cells close into them. Only seeds whose image
#              stays beyond the latitude are mirrored (otherwise it wraps over the pole into
#              the band). That mirror is not an isometry, so boundary vertices are snapped
#              onto the latitude.
#              Printability as for box/cylinder but on true 3D angles. Nodes are lifted
#              along the meridian (up the curved wall), long shallow struts become arches
#              with the apex lifted along the meridian, and clashing arches stay flat.
#              SDF: each voxel near the shell is projected onto the mid-surface (Newton on
#              the meridian ellipse); struts are distances from that projection. Rim
#              gussets grow per voxel slice (z = print layers): each layer may extend the
#              one below by dz / tan(min-angle), inside a pointed-arch envelope per gap.
#              Rejected: a 2D chart field as for box/cylinder (non-developable wall).
# Constraints: Same packages as potlib.py; peak memory is above the box (shell arrays).
# Invariants:  args.height is set to the derived rim height before potlib checks run; they
#              use it for the rim zone.
# Tests:       tests/test_regression.py::test_orb_graph_smoke (no golden yet, experimental).
#              Per run: report checks plus stability and rim_unfilled_max_mm.
# Limitations: The anchor snapping moves boundary vertices slightly. Opening metrics use a
#              chart at the belly radius (sizes near floor/rim overstated by up to 1/0.8).
#              Strut width is measured along chords, not geodesics.
# === END CONTEXT ===
import json
import math
import time
from types import SimpleNamespace

import numpy as np
from scipy import ndimage
from scipy.spatial import SphericalVoronoi, cKDTree

from potlib import (Grid, fingerprint, finish, floor_field, floor_segments, grade_fraction, pt_seg_dist,
                    require_rim_crossings, round_isect, rounded_box2, smin, wall_metrics)


class OrbShape:
    """Mid-surface ellipsoid (a horizontal, c vertical semi-axis), floor at z = 0.
    Direction d on the unit sphere <-> point (a dx, a dy, zc + c dz); polar angle phi is
    measured from the bottom pole, so dz = -cos(phi)."""

    def __init__(self, args):
        self.t = args.wall
        self.a_out = args.diameter / 2
        self.a = self.a_out - self.t / 2
        self.c = args.stretch * self.a_out - self.t / 2
        al = math.radians(args.min_wall_angle)
        self.phi_f = math.atan(math.tan(al) * self.a / self.c)    # slope = atan(c/a tan phi)
        self.zc = self.c * math.cos(self.phi_f)
        self.H = 2 * self.zc
        self.r_floor = self.a * math.sin(self.phi_f)
        # outer face at z = 0: the wall (slope al) is t / sin(al) wide along a horizontal
        self.r_floor_out = self.r_floor + self.t / 2 / math.sin(al)
        self.extent = self.a_out + self.t

    def z_of_phi(self, phi):
        return self.zc - self.c * np.cos(phi)

    def phi_of_z(self, z):
        return np.arccos(np.clip((self.zc - np.asarray(z, float)) / self.c, -1, 1))

    def to3d(self, d):
        d = np.asarray(d, float)
        return np.stack([self.a * d[..., 0], self.a * d[..., 1], self.zc + self.c * d[..., 2]], -1)

    def footprint_sdf(self, x, y):
        # outer footprint, as for the box: the floor (clipped 1 mm inside it) then overlaps
        # the outward-leaning wall over its full thickness; the mid radius left a groove
        return np.hypot(x, y) - self.r_floor_out

    def describe(self):
        return {"belly_dia_mm": 2 * self.a_out, "floor_opening_dia_mid_mm": round(2 * self.r_floor, 1),
                "rim_height_mm": round(self.H, 1)}


def dir_of(phi, th):
    return np.stack([np.sin(phi) * np.cos(th), np.sin(phi) * np.sin(th), -np.cos(phi)], -1)


def phi_th(d):
    return np.arccos(np.clip(-d[..., 2], -1, 1)), np.arctan2(d[..., 1], d[..., 0])


def slerp_path(d0, d1, sh, step=3.0):
    """Points along the great circle d0 -> d1, mapped to the ellipsoid (3D polyline)."""
    ang = math.acos(min(1.0, max(-1.0, float(d0 @ d1))))
    n = max(1, math.ceil(ang * sh.a / step))
    if ang < 1e-9:
        return sh.to3d(np.array([d0, d1]))
    ts = np.linspace(0, 1, n + 1)[:, None]
    pts = (np.sin((1 - ts) * ang) * d0 + np.sin(ts * ang) * d1) / math.sin(ang)
    return sh.to3d(pts / np.linalg.norm(pts, axis=1, keepdims=True))


def densify(p, step=1.0):
    seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
    out = [a + (b - a) * np.linspace(0, 1, max(1, math.ceil(L / step)), endpoint=False)[:, None]
           for a, b, L in zip(p[:-1], p[1:], seg)]
    return np.vstack(out + [p[-1:]])


def steep_enough(p, q, tan):
    """Chord from the lower to the higher point rises at least atan(tan)."""
    lo, hi = (p, q) if p[2] <= q[2] else (q, p)
    return hi[2] - lo[2] >= tan * math.hypot(hi[0] - lo[0], hi[1] - lo[1])


# ---------- pattern ----------

# CTX: Returns struts as 3D polylines ON the mid-surface (anchor edges excluded); every
#      consumer (orb_metrics, rim_gaps, build_orb_sdf) relies on points lying on it.
# Used-by: orb.py::run
# Invariant: point index ranges [0, n) seeds, [n, n + nb) band mirrors, [n + nb, ...) rim
#            mirrors decide which vertices are anchors; keep them in sync with the masks.
def orb_graph(args, sh, rng, stats):
    z0, z1 = args.band / 2 + 0.5, sh.H - args.rim / 2           # anchor latitudes (as box)
    p0, p1 = float(sh.phi_of_z(z0)), float(sh.phi_of_z(z1))
    zg = sh.H * (1 - args.grade_top)
    cell_at = lambda z: args.cell * (1 - (1 - args.grade_scale) * grade_fraction(z, zg, z1))

    # Poisson seeds on the sphere of radius a, area-uniform in the band p0..p1
    seeds, rad = [], []
    cz0, cz1 = -math.cos(p0), -math.cos(p1)
    for _ in range(60000):
        dz = cz0 + rng.random() * (cz1 - cz0)
        th = rng.random() * 2 * math.pi
        s_ = math.sqrt(1 - dz * dz)
        d = np.array([s_ * math.cos(th), s_ * math.sin(th), dz])
        r = cell_at(sh.zc + sh.c * dz)
        if seeds and np.any(sh.a * np.linalg.norm(np.asarray(seeds) - d, axis=1) < (np.asarray(rad) + r) / 2):
            continue
        seeds.append(d); rad.append(r)
    seeds = np.asarray(seeds)
    n = len(seeds)
    ph, th = phi_th(seeds)
    # Mirror in polar angle only seeds whose image stays beyond its anchor latitude: for
    # ph > 2*p0 (or < 2*p1 - pi) the image wraps over the pole back into the band and
    # pins interior vertices to the anchor line (seen at --stretch 1.4 --min-wall-angle 50)
    mb, mr = ph < 2 * p0, ph > 2 * p1 - math.pi
    nb = int(mb.sum())
    pts = np.vstack([seeds, dir_of(2 * p0 - ph[mb], th[mb]), dir_of(2 * p1 - ph[mr], th[mr])])
    sv = SphericalVoronoi(pts, radius=1.0)
    sv.sort_vertices_of_regions()
    owner = {}                                     # vertex -> set of point indices
    for ip, reg in enumerate(sv.regions):
        for v in reg:
            owner.setdefault(v, set()).add(ip)
    raw = set()
    for ip in range(n):
        reg = sv.regions[ip]
        for k in range(len(reg)):
            a_, b_ = reg[k], reg[(k + 1) % len(reg)]
            raw.add((min(a_, b_), max(a_, b_)))
    used = sorted({v for e in raw for v in e})
    side = {}                                      # boundary vertices: 0 = band, 1 = rim
    D = {}
    for v in used:
        o = owner[v]
        d = sv.vertices[v].copy()
        if any(n <= q < n + nb for q in o):
            side[v] = 0
        elif any(q >= n + nb for q in o):
            side[v] = 1
        if v in side:                              # snap onto the anchor latitude
            _, t_ = phi_th(d)
            d = dir_of(p0 if side[v] == 0 else p1, t_)
        D[v] = d

    # union-find merge of short edges (3D distance on the ellipsoid)
    parent = {v: v for v in used}

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]; v = parent[v]
        return v

    for a_, b_ in raw:
        if np.linalg.norm(sh.to3d(D[a_]) - sh.to3d(D[b_])) < args.merge:
            parent[find(a_)] = find(b_)
    groups = {}
    for v in used:
        groups.setdefault(find(v), []).append(v)
    nid, dirs, bside = {}, [], []
    for g in groups.values():
        d = np.mean([D[v] for v in g], axis=0); d /= np.linalg.norm(d)
        sides = {side[v] for v in g if v in side}
        sd = min(sides) if sides else -1
        if sd >= 0:
            d = dir_of(p0 if sd == 0 else p1, phi_th(d)[1])
        for v in g:
            nid[v] = len(dirs)
        dirs.append(d); bside.append(sd)
    dirs, bside = np.asarray(dirs), np.asarray(bside)
    E = {(min(nid[a_], nid[b_]), max(nid[a_], nid[b_])) for a_, b_ in raw if nid[a_] != nid[b_]}
    stats["voronoi_nodes"], stats["voronoi_edges"] = len(dirs), len(E)

    bottom, top = bside == 0, bside == 1
    tan_min = math.tan(math.radians(args.min_angle))
    tan_arch = math.tan(math.radians(args.min_angle + 5))
    P3 = lambda i: sh.to3d(dirs[i])
    anchor = lambda i, j: (bottom[i] and bottom[j]) or (top[i] and top[j])

    # Lift (bottom-up): move a node up its meridian until one lower strut is steep enough
    adj = {i: [j for e in E if i in e and not anchor(*e) for j in e if j != i] for i in range(len(dirs))}
    lifts = []
    for i in sorted(range(len(dirs)), key=lambda i: P3(i)[2]):
        if bottom[i] or top[i] or not adj[i]:
            continue
        zi = P3(i)[2]
        lower = [m for m in adj[i] if P3(m)[2] <= zi] or adj[i]
        if any(steep_enough(P3(m), P3(i), tan_min) for m in lower):
            continue
        phi_i, th_i = phi_th(dirs[i])
        ok = lambda ph: any(steep_enough(P3(m), sh.to3d(dir_of(ph, th_i)), tan_min) for m in lower)
        lo_, hi_ = float(phi_i), p1
        if not ok(hi_):
            dirs[i] = dir_of(p1, th_i); top[i] = True
        else:
            for _ in range(40):
                mid = (lo_ + hi_) / 2
                lo_, hi_ = (lo_, mid) if ok(mid) else (mid, hi_)
            dirs[i] = dir_of(hi_, th_i)
        lifts.append(float(P3(i)[2] - zi))
    stats["nodes_lifted"] = len(lifts)
    stats["lift_mm_median/max"] = [round(float(np.median(lifts)), 1) if lifts else 0,
                                   round(float(max(lifts)), 1) if lifts else 0]

    def horiz(i, j):
        p, q = P3(i), P3(j)
        return math.hypot(p[0] - q[0], p[1] - q[1])

    def classify(i, j):
        if anchor(i, j):
            return "plane"
        if steep_enough(P3(i), P3(j), tan_min):
            return "steep"
        return "short" if horiz(i, j) <= args.max_bridge else "long"

    arch_tan, apex_cache = {}, {}

    def apex(e):
        """Apex direction: edge midpoint lifted along its meridian until both legs reach
        the arch angle; None if that needs to go above the rim anchor latitude."""
        key = (e, arch_tan.get(e, tan_arch))
        if key in apex_cache:
            return apex_cache[key]
        tn = key[1]
        m = dirs[e[0]] + dirs[e[1]]; m /= np.linalg.norm(m)
        ph_m, th_m = phi_th(m)
        ok = lambda ph: all(steep_enough(P3(v), sh.to3d(dir_of(ph, th_m)), tn) for v in e)
        res = None
        if ok(p1):
            lo_, hi_ = float(ph_m), p1
            for _ in range(40):
                mid = (lo_ + hi_) / 2
                lo_, hi_ = (lo_, mid) if ok(mid) else (mid, hi_)
            res = dir_of(hi_, th_m)
        apex_cache[key] = res
        return res

    def legs_of(e):
        """Arch legs (from e[0], from e[1]) resampled at <= 1 mm for the clash test."""
        ap = apex(e)
        return [densify(slerp_path(dirs[v], ap, sh)) for v in e]

    def shape_of(e):
        if cls[e] == "long":
            ap = apex(e)
            return np.vstack([slerp_path(dirs[e[0]], ap, sh), slerp_path(ap, dirs[e[1]], sh)[1:]])
        return slerp_path(dirs[e[0]], dirs[e[1]], sh)

    cls = {e: classify(*e) for e in E}
    n_flat = 0
    for e in E:
        if cls[e] == "long" and (apex(e) is None or horiz(*e) > args.max_arch):
            cls[e] = "short"; n_flat += 1

    def samples(E):
        pts, own = [], []
        for e in E:
            if cls[e] == "plane":
                continue
            p = shape_of(e)
            seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
            for a_, b_, L in zip(p[:-1], p[1:], seg):
                k_ = max(1, math.ceil(L))
                pts.append(a_ + (b_ - a_) * np.linspace(0, 1, k_, endpoint=False)[:, None]); own += [e] * k_
        return np.vstack(pts), own

    removed = 0
    while True:
        adj = {i: [] for i in range(len(dirs))}
        for i, j in E:
            adj[i].append(j); adj[j].append(i)
        alive = {i for i in adj if adj[i]}
        sup = {}
        for i in sorted(alive, key=lambda i: P3(i)[2]):
            sup[i] = bottom[i] or any(sup.get(j, False) and P3(j)[2] < P3(i)[2] - 1e-4
                                      and cls[(min(i, j), max(i, j))] == "steep" for j in adj[i])
        bad = {i for i in alive if not sup[i] and not top[i]}
        bad |= {i for i in alive if len(adj[i]) < 2 and not (bottom[i] or top[i])}
        bad_e = {e for e in E if cls[e] in ("short", "long") and not (sup[e[0]] and sup[e[1]])}
        if bad or bad_e:
            removed += len(bad)
            E = {e for e in E if e not in bad_e and e[0] not in bad and e[1] not in bad}
            continue
        # arches must keep one strut width from non-adjacent struts; else flatter, then flat
        pts, own = samples(E)
        tree = cKDTree(pts)
        clash = None
        for e in sorted((e for e in E if cls[e] == "long"), key=lambda e: -horiz(*e)):
            # per leg: only struts at that leg's own end node are exempt (as for the box)
            for v, leg in zip(e, legs_of(e)):
                near = {own[k] for ks in tree.query_ball_point(leg, args.strut) for k in ks}
                if any(o != e and v not in o for o in near):
                    clash = e
            if clash:
                break
        if clash is None:
            break
        if clash not in arch_tan:
            arch_tan[clash] = tan_min
            if apex(clash) is None:
                cls[clash] = "short"; n_flat += 1
        else:
            cls[clash] = "short"; n_flat += 1
    stats["removed_unsupported_or_leaf_nodes"] = removed
    stats["shallow_kept_flat"] = n_flat
    stats["edges_final"] = len(E)
    stats["edge_classes"] = {c: sum(cls[e] == c for e in E) for c in ("plane", "steep", "short", "long")}
    polys = [shape_of(e) for e in E if cls[e] != "plane"]
    # struts closer than one width without sharing a node (0 expected)
    pts, own = samples(E)
    tree = cKDTree(pts)
    stats["strut_clearance_violations"] = len({tuple(sorted((own[i], own[j]))) for i, j in
                                               tree.query_pairs(args.strut * 0.5) if not (set(own[i]) & set(own[j]))})
    return polys


def orb_metrics(args, sh, polys, stats):
    """Opening sizes etc. in a chart u = theta * belly radius (approximate, see header)."""
    chart = []
    for p in polys:
        th = np.unwrap(np.arctan2(p[:, 1], p[:, 0]))
        chart.append(np.column_stack([th * sh.a, p[:, 2]]))
    wall_metrics(args, SimpleNamespace(P=2 * math.pi * sh.a), chart, stats)


def rim_gaps(args, sh, polys, stats):
    zr = sh.H - args.rim
    ths = []
    for p in polys:
        for a, b in zip(p[:-1], p[1:]):
            if (a[2] - zr) * (b[2] - zr) < 0:
                q = a + (zr - a[2]) / (b[2] - a[2]) * (b - a)
                ths.append(math.atan2(q[1], q[0]) % (2 * math.pi))
    ths.sort()
    require_rim_crossings(len(ths))
    rr = sh.a * math.sin(float(sh.phi_of_z(zr)))
    gaps = [(ths[(k + 1) % len(ths)] - t) % (2 * math.pi) for k, t in enumerate(ths)]
    stats["rim_gaps_mm_p50/max"] = [round(float(np.median(gaps)) * rr, 1), round(float(max(gaps)) * rr, 1)]
    return [((t + g / 2) % (2 * math.pi), g / 2) for t, g in zip(ths, gaps) if g * rr > args.strut]


# ---------- SDF ----------

def build_orb_sdf(args, sh, polys, gaps, fsegs, stats):
    k, t = args.fillet, sh.t
    rt = args.eyelet_wire / 2
    Rt = args.eyelet_hole / 2 + rt
    H = sh.H
    ez = H + 0.5 + args.eyelet_hole / 2
    o = 0.137 * args.voxel
    E_ = sh.extent + 2.0
    ztop = ez + Rt + rt + 1.0 if args.eyelets else H + 1.0
    g = Grid([-E_ + o, -E_ + o, -1.0 + o], [E_, E_, ztop], args.voxel)
    X2, Y2 = np.meshgrid(g.ax[0], g.ax[1], indexing="ij")
    rho2, th2 = np.hypot(X2, Y2), np.arctan2(Y2, X2)

    fl = floor_field(g, fsegs, sh.footprint_sdf(X2, Y2), args)

    # shell voxels: project onto the mid-surface meridian ellipse (Newton on phi)
    a, c, zc = sh.a, sh.c, sh.zc
    sel_idx, Ep, W, Zs = [], [], [], []
    for iz, z in enumerate(g.ax[2]):
        q = np.hypot(rho2 / a, (z - zc) / c)
        m = np.abs(q - 1) * min(a, c) < t / 2 + 3.0
        if not m.any():
            continue
        rho, th = rho2[m], th2[m]
        phi = np.arctan2(rho / a, (zc - z) / c)
        for _ in range(5):
            ex, ezz = a * np.sin(phi), zc - c * np.cos(phi)
            dx, dzz = a * np.cos(phi), c * np.sin(phi)
            f = (ex - rho) * dx + (ezz - z) * dzz
            fp = dx * dx + dzz * dzz + (ex - rho) * (-a * np.sin(phi)) + (ezz - z) * (c * np.cos(phi))
            phi = phi - f / fp
        ex, ezz = a * np.sin(phi), zc - c * np.cos(phi)
        nr, nz = c * np.sin(phi), -a * np.cos(phi)
        w = np.hypot(rho - ex, z - ezz) * np.sign((rho - ex) * nr + (z - ezz) * nz)
        ii, jj = np.nonzero(m)
        sel_idx.append(np.column_stack([ii, jj, np.full(len(ii), iz)]))
        Ep.append(np.column_stack([ex * np.cos(th), ex * np.sin(th), ezz]))
        W.append(w); Zs.append(np.full(len(ii), z))
    sel_idx, Ep, W, Zs = np.vstack(sel_idx), np.vstack(Ep), np.concatenate(W), np.concatenate(Zs)
    stats["shell_voxels"] = int(len(W))

    # strut field on the projected points (hard-min within a strut, smooth-min between)
    order = np.argsort(Ep[:, 2]); Ez_sorted = Ep[order, 2]
    Dst = np.full(len(W), 10.0)
    pad = args.strut / 2 + k + 1
    for p in polys:
        lo_, hi_ = p.min(0) - pad, p.max(0) + pad
        i0, i1 = np.searchsorted(Ez_sorted, [lo_[2], hi_[2]])
        cand = order[i0:i1]
        cand = cand[np.all((Ep[cand, :2] >= lo_[:2]) & (Ep[cand, :2] <= hi_[:2]), axis=1)]
        if not len(cand):
            continue
        d = np.min([pt_seg_dist(Ep[cand], s0, s1) for s0, s1 in zip(p[:-1], p[1:])], 0)
        Dst[cand] = smin(Dst[cand], d - args.strut / 2, k)
    aw = np.abs(W)
    wall = round_isect(Dst, aw - t / 2, args.edge_round)
    band = rounded_box2(aw, np.abs(Zs - (args.band / 2 - 1)), t / 2, args.band / 2 + 1, 1.5)
    rim = rounded_box2(aw, np.abs(Zs - (H - args.rim / 2)), t / 2, args.rim / 2, 2.0)

    for iz, z in enumerate(g.ax[2]):
        g.D[:, :, iz] = np.maximum(fl, abs(z - (args.floor / 2 - 0.5)) - (args.floor / 2 + 0.5))
    ix, iy, izz = sel_idx.T
    flo = g.D[ix, iy, izz]
    g.D[ix, iy, izz] = smin(wall, smin(smin(band, rim, 2.0), flo, 2.0), k)

    # eyelets: flat vertical rings on the rim, in the rim's tangent plane
    rtop = a * math.sin(float(sh.phi_of_z(H)))
    for i in range(args.eyelets):
        thc = 2 * math.pi * i / args.eyelets
        cen = np.array([rtop * math.cos(thc), rtop * math.sin(thc), ez])
        sl, X, Y, Z = g.box(cen - Rt - rt - t, cen + Rt + rt + t)
        tx, ty = -math.sin(thc), math.cos(thc)
        a1 = (X - cen[0]) * tx + (Y - cen[1]) * ty
        ar = (X - cen[0]) * ty - (Y - cen[1]) * tx         # radial offset from the ring plane
        ring = round_isect(np.abs(np.hypot(a1, Z - cen[2]) - Rt) - rt, np.abs(ar) - t / 2, args.edge_round)
        g.D[sl] = smin(g.D[sl], ring, k)

    grow_rim_gussets(args, sh, g, gaps, sel_idx, aw, stats)
    for iz, z in enumerate(g.ax[2]):
        g.D[:, :, iz] = np.maximum(g.D[:, :, iz], -z)
    return g


# CTX: Layer-wise gusset growth in the voxel grid (z slices = print layers).
# Used-by: orb.py::build_orb_sdf | Depends: orb.py::rim_gaps (gap list contract)
# Invariant: runs after all struts, rim and eyelets are in g.D; merges the ADDED voxels
#            only (blur + smooth union), same rule as orchid_pot.py::wall_pattern.
# Rejected: SDF of all grown material (shifted every surface ~0.2 mm, v0.6.0).
def grow_rim_gussets(args, sh, g, gaps, sel_idx, aw, stats):
    """Per voxel slice (= print layers): material may extend the layer below by
    dz / tan(min-angle) sideways, inside the shell and inside a pointed-arch envelope per
    rim gap (flanks at min-angle + 5 deg, apex 1 mm inside the rim)."""
    zr = sh.H - args.rim
    te = math.tan(math.radians(args.min_angle + 5))
    ts = math.tan(math.radians(args.min_angle))
    dz = g.h
    rr = sh.a * math.sin(float(sh.phi_of_z(zr)))
    hmax = te * max([hw for _, hw in gaps], default=0) * rr
    rows = np.where((g.ax[2] > zr - hmax - 2 * dz) & (g.ax[2] < zr + 1.0))[0]
    if not gaps or not len(rows):
        return
    X2, Y2 = np.meshgrid(g.ax[0], g.ax[1], indexing="ij")
    rho2, th2 = np.hypot(X2, Y2), np.arctan2(Y2, X2) % (2 * math.pi)
    arc = np.full(X2.shape, np.inf)                # |arc distance to the gap midpoint|
    for thm, hw in gaps:
        dth = np.abs((th2 - thm + math.pi) % (2 * math.pi) - math.pi)
        arc = np.where(dth < hw, np.minimum(arc, dth * rho2), arc)
    shell = np.zeros(g.n, bool)
    inside = aw <= sh.t / 2
    shell[tuple(sel_idx[inside].T)] = True
    S = g.D[:, :, rows] <= 0
    G = S.copy()
    for k_ in range(1, len(rows)):
        z = g.ax[2][rows[k_]]
        env = (z >= zr + 1.0 - te * arc) & shell[:, :, rows[k_]]
        if env.any() and G[:, :, k_ - 1].any():
            d1 = ndimage.distance_transform_edt(~G[:, :, k_ - 1]) * g.h
            G[:, :, k_] |= (d1 <= dz / ts) & env
    added = G & ~S
    # remaining bridge under the rim: longest unfilled arc of the mid-surface circle in the
    # last layer below the rim (box: rim_unfilled_max_mm; <= 6 mm keeps the reach <= 3 mm)
    kb = int(np.searchsorted(g.ax[2][rows], zr)) - 1
    na = int(2 * math.pi * rr / 0.2)
    tha = np.arange(na) * 2 * math.pi / na
    ix = np.clip(np.round((rr * np.cos(tha) - g.lo[0]) / g.h).astype(int), 0, g.n[0] - 1)
    iy = np.clip(np.round((rr * np.sin(tha) - g.lo[1]) / g.h).astype(int), 0, g.n[1] - 1)
    filled = G[ix, iy, kb]
    runs, cur = [0], 0
    for v in np.concatenate([filled, filled]):
        cur = 0 if v else cur + 1; runs.append(cur)
    stats["rim_unfilled_max_mm"] = round(min(max(runs), na) * 2 * math.pi * rr / na, 1)
    stats["rim_gussets"] = len(gaps)
    stats["gusset_volume_cm3"] = round(float(added.sum()) * g.h ** 3 / 1000, 1)
    # SDF of the ADDED material only: inside the original solid it is >= +h, so min() keeps
    # the smooth field there. An SDF of all of G (<= -h inside) replaced it and pushed every
    # surface in these rows ~0.2 mm outward. Blended with the strut fillet like every other
    # union: a hard min (or np.where switch) leaves creases that break the decimation.
    sd = ndimage.distance_transform_edt(~added) * g.h - ndimage.distance_transform_edt(added) * g.h
    # 1-voxel blur: removes 1-voxel slivers/pinches where two gusset fronts meet near the
    # apex (they made a second body / non-manifold spots after export and decimation)
    sd = ndimage.gaussian_filter(sd, 1.0, mode="nearest")
    g.D[:, :, rows] = smin(g.D[:, :, rows], sd, args.fillet)


def stability(sh):
    def check(mesh, stats):
        com = mesh.center_mass
        sec = mesh.section(plane_origin=[0, 0, 0.3], plane_normal=[0, 0, 1])
        r_contact = float(np.max(np.hypot(sec.vertices[:, 0], sec.vertices[:, 1]))) if sec is not None else 0.0
        margin = r_contact - math.hypot(com[0], com[1])
        stats["stability"] = {"floor_contact_radius_mm": round(r_contact, 1),
                              "center_of_mass_z_mm": round(float(com[2]), 1),
                              "tip_over_angle_deg": round(math.degrees(math.atan2(margin, com[2])), 1)}
    return check


def run(args):
    sh = OrbShape(args)
    args.height = sh.H                             # potlib checks use it for the rim zone
    rng = np.random.default_rng(args.seed)
    stats = {"params": vars(args), "shape": {"type": "orb (experimental)", **sh.describe()}}
    t0 = time.time()
    polys = orb_graph(args, sh, rng, stats)
    stats["fingerprint"] = fingerprint(args, polys)      # accessories match against this
    orb_metrics(args, sh, polys, stats)
    gaps = rim_gaps(args, sh, polys, stats)
    fsegs = floor_segments(args, sh, rng)
    print(f"graph {time.time() - t0:.1f}s: {len(polys)} wall struts, {len(fsegs)} floor struts")
    if args.graph_only:
        print(json.dumps({k: v for k, v in stats.items() if k != "params"}, indent=1))
        return
    g = build_orb_sdf(args, sh, polys, gaps, fsegs, stats)
    print(f"sdf {time.time() - t0:.1f}s grid {g.n.tolist()}")
    finish(g, args, stats, t0, extra_checks=stability(sh))
    return polys                               # the clip tray reuses the strut paths
