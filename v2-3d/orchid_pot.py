#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# === CONTEXT ===
# Purpose:     Generate a printable STL of a hanging orchid basket (Vanda style):
#              rounded-square box whose walls are a 3D Voronoi strut lattice inside a
#              thick shell, with a solid bottom band, top rim, 2D Voronoi floor and
#              four hanging eyelets. Target: PETG, 0.4 mm nozzle, printed upright.
# Decisions:   3D Voronoi is computed in a flattened wall frame (u = arc length along
#              the wall mid-curve, z, w = depth through the wall) and mapped back onto
#              the rounded square. Seeds are mirrored across w and z planes, so cells
#              are clipped exactly to the shell: the inner and outer faces get closed
#              strut networks, with no dangling stubs. Rejected: clipping a free 3D
#              Voronoi to the shell, which leaves spikes on both surfaces.
#              z is compressed by --zstretch before the Voronoi and stretched back
#              afterwards. This elongates cells vertically, so more struts are steep
#              enough to print without supports.
#              Printability is enforced on the strut graph, not the mesh. Every node
#              must reach a lower node through a steep strut (no node starts in mid-air).
#              Shallow struts: short ones stay as flat bridges; long ones become pointed
#              arches (two legs >= min-angle+5 deg). An unsupported node gets one rescue
#              strut down to a supported node before it is removed. Rejected: dropping
#              shallow struts / unsupported nodes. At zstretch 1.5 that left the faces
#              with ~13% non-vertical strut length and openings up to the full wall
#              height. zstretch 1.2 is the sweet spot: at 1.0, ~125 nodes need rescue struts.
#              Mesh = smooth-min SDF of capsules + marching cubes. The smooth union
#              puts fillets at the nodes (stronger, no notch stress). Rejected: Blender
#              boolean unions (slow and fragile with ~1000 cylinders).
# Constraints: numpy, scipy, scikit-image, trimesh, fast_simplification, shapely.
#              The voxel grid at 0.4 mm is ~45M float32 (~180 MB) and peaks at a few GB
#              during marching cubes.
# Invariants:  Mirror planes sit at strut-radius distance inside the shell/band/rim, so
#              the in-plane networks stay inside the nominal envelope. Strut pieces of
#              one edge are hard-min'ed before the smooth-min; otherwise the joints bulge.
# Tests:       The report's mesh checks (watertight, islands per layer, overhang area);
#              preview with render.py.
# Limitations: The overhang angle uses (u, z, w) as a local isometry of the wall. That
#              is exact on the flat sides and approximate at the corners. Bridges under
#              the top rim are not removed; their max span is reported instead.
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
from scipy.spatial import Voronoi
from shapely.ops import unary_union
from skimage.measure import marching_cubes


def parse_args():
    p = argparse.ArgumentParser(description="Voronoi orchid basket STL generator")
    a = p.add_argument
    a("--size", type=float, default=150.0, help="outer width/depth [mm]")
    a("--height", type=float, default=100.0, help="body height without eyelets [mm]")
    a("--corner-r", type=float, default=20.0, help="outer vertical edge radius [mm]")
    a("--wall", type=float, default=12.0, help="nominal wall thickness [mm]")
    a("--strut", type=float, default=4.0, help="strut diameter [mm]")
    a("--cell", type=float, default=17.0, help="seed spacing of the wall Voronoi [mm]")
    a("--zstretch", type=float, default=1.2, help="vertical cell elongation factor")
    a("--min-angle", type=float, default=35.0, help="min strut elevation [deg]")
    a("--max-bridge", type=float, default=10.0, help="max horizontal span of shallow struts [mm]")
    a("--max-arch", type=float, default=25.0, help="max horizontal span of an arch [mm]")
    a("--rescue-len", type=float, default=25.0, help="max length of a rescue support strut [mm]")
    a("--merge", type=float, default=3.0, help="merge Voronoi nodes closer than this [mm]")
    a("--band", type=float, default=5.0, help="bottom band height [mm]")
    a("--rim", type=float, default=6.0, help="top rim height [mm]")
    a("--floor", type=float, default=5.0, help="floor lattice thickness [mm]")
    a("--floor-cell", type=float, default=16.0, help="floor Voronoi seed spacing [mm]")
    a("--floor-strut", type=float, default=4.5, help="floor strut width [mm]")
    a("--eyelet-hole", type=float, default=6.0, help="eyelet hole diameter [mm]")
    a("--eyelet-wire", type=float, default=4.0, help="eyelet ring thickness [mm]")
    a("--fillet", type=float, default=1.5, help="smooth-union radius at nodes [mm]")
    a("--voxel", type=float, default=0.4, help="SDF grid pitch [mm]")
    a("--decimate", type=float, default=0.75, help="fraction of triangles to remove")
    a("--layer", type=float, default=0.2, help="layer height for island check [mm]")
    a("--seed", type=int, default=7)
    a("--out", default="orchid_pot")
    a("--graph-only", action="store_true", help="stop after the strut graph and its metrics")
    return p.parse_args()


# ---------- geometry helpers ----------

class WallFrame:
    """Rounded-square mid-curve of the wall; maps (u, z, w) -> (x, y, z)."""

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


def smin(a, b, k):
    h = np.clip(0.5 + 0.5 * (b - a) / k, 0.0, 1.0)
    return b * (1 - h) + a * h - k * h * (1 - h)


def poisson(rng, lo, hi, spacing, period=None, tries=60000):
    """Dart throwing. period: index of a periodic axis -> length."""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    pts = []
    for _ in range(tries):
        p = lo + rng.random(len(lo)) * (hi - lo)
        if pts:
            d = np.abs(np.asarray(pts) - p)
            if period:
                for ax, L in period.items():
                    d[:, ax] = np.minimum(d[:, ax], L - d[:, ax])
            if np.min(np.sqrt((d ** 2).sum(1))) < spacing:
                continue
        pts.append(p)
    return np.asarray(pts)


# ---------- wall lattice graph ----------

def wall_graph(args, fr, rng, stats):
    r = args.strut / 2
    s = args.zstretch
    z0, z1 = args.band / 2 + 0.5, args.height - args.rim / 2
    w0, w1 = r, fr.t - r
    # Voronoi in compressed space: z' = z / s
    seeds = poisson(rng, [0, z0 / s, w0], [fr.P, z1 / s, w1], args.cell, period={0: fr.P})
    A = [seeds]
    for plane, ax in ((w0, 2), (w1, 2)):
        m = seeds.copy(); m[:, ax] = 2 * plane - m[:, ax]; A.append(m)
    A = np.vstack(A)
    B = [A]
    for plane in (z0 / s, z1 / s):
        m = A.copy(); m[:, 1] = 2 * plane - m[:, 1]; B.append(m)
    B = np.vstack(B)
    C = np.vstack([B, B + [fr.P, 0, 0], B - [fr.P, 0, 0]])
    vor = Voronoi(C)
    V = vor.vertices.copy()
    V[:, 1] *= s
    eps = 1e-4
    inside = ((V[:, 1] > z0 - eps) & (V[:, 1] < z1 + eps)
              & (V[:, 2] > w0 - eps) & (V[:, 2] < w1 + eps))

    # Edges = consecutive vertices of each ridge polygon, ordered by angle in its plane
    edges = set()
    for (pa, pb), rv in zip(vor.ridge_points, vor.ridge_vertices):
        if -1 in rv or len(rv) < 3:
            continue
        rv = np.asarray(rv)
        if not inside[rv].all():
            continue
        nrm = C[pb] - C[pa]; nrm /= np.linalg.norm(nrm)
        e1 = np.cross(nrm, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1)
        e2 = np.cross(nrm, e1)
        q = vor.vertices[rv] - vor.vertices[rv].mean(0)
        order = rv[np.argsort(np.arctan2(q @ e2, q @ e1))]
        for i in range(len(order)):
            a_, b_ = order[i], order[(i + 1) % len(order)]
            edges.add((min(a_, b_), max(a_, b_)))

    # Canonical nodes: wrap u, then union-find merge of near nodes
    used = sorted({v for e in edges for v in e})
    by_pos, key = {}, {}
    nodes = []
    for v in used:
        u = V[v, 0] % fr.P
        k = (round(u, 3) % round(fr.P, 3), round(V[v, 1], 3), round(V[v, 2], 3))
        if k not in by_pos:
            by_pos[k] = len(nodes); nodes.append([u, V[v, 1], V[v, 2]])
        key[v] = by_pos[k]
    nodes = np.asarray(nodes)
    E = {(min(key[a], key[b]), max(key[a], key[b])) for a, b in edges if key[a] != key[b]}

    def dvec(i, j):
        d = nodes[j] - nodes[i]
        d[0] = (d[0] + fr.P / 2) % fr.P - fr.P / 2
        return d

    parent = list(range(len(nodes)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i

    for i, j in E:
        if np.linalg.norm(dvec(i, j)) < args.merge:
            parent[find(i)] = find(j)
    groups = {}
    for i in range(len(nodes)):
        groups.setdefault(find(i), []).append(i)
    new_id, merged = {}, []
    for g in groups.values():
        ref = nodes[g[0]]
        ds = np.array([dvec(g[0], i) for i in g])
        p = ref + ds.mean(0)
        p[0] %= fr.P
        for ax, lo_, hi_ in ((1, z0, z1), (2, w0, w1)):   # keep plane nodes on their planes
            pl = nodes[g, ax]
            if np.any(np.abs(pl - lo_) < eps): p[ax] = lo_
            elif np.any(np.abs(pl - hi_) < eps): p[ax] = hi_
        for i in g:
            new_id[i] = len(merged)
        merged.append(p)
    nodes = np.asarray(merged)
    E = {(min(new_id[a], new_id[b]), max(new_id[a], new_id[b]))
         for a, b in E if new_id[a] != new_id[b]}
    stats["voronoi_nodes"], stats["voronoi_edges"] = len(nodes), len(E)

    # ---- printability filter ----
    bottom = np.abs(nodes[:, 1] - z0) < eps
    top = np.abs(nodes[:, 1] - z1) < eps
    tan_min = math.tan(math.radians(args.min_angle))

    def classify(i, j):
        d = dvec(i, j)
        hl = math.hypot(d[0], d[2])
        if (bottom[i] and bottom[j]) or (top[i] and top[j]):
            return "plane"
        if abs(d[1]) >= tan_min * hl:
            return "steep"
        return "short" if hl <= args.max_bridge else "long"

    # Long shallow struts become pointed arches: two legs rising at arch_angle to an apex
    # above the midpoint, like a gothic arch. This keeps the lateral link, which plain
    # removal would lose (removal turned the faces into tall vertical slits). Arches whose
    # apex would leave the lattice zone are dropped; the rim braces that height anyway.
    tan_arch = math.tan(math.radians(args.min_angle + 5))

    def apex(i, j):
        d = dvec(i, j)
        p = nodes[i] + d / 2
        p[1] = max(nodes[i, 1], nodes[j, 1]) + tan_arch * math.hypot(d[0], d[2]) / 2
        return p

    cls = {e: classify(*e) for e in E}
    for e in E:
        d = dvec(*e)
        if cls[e] == "long" and (apex(*e)[1] > z1 or math.hypot(d[0], d[2]) > args.max_arch):
            cls[e] = "drop"
    stats["dropped_arch_too_high_or_wide"] = sum(c == "drop" for c in cls.values())
    E = {e for e in E if cls[e] != "drop"}
    removed_nodes = 0
    rescued, n_rescue = set(), 0
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
        # Rescue an unsupported node with one extra steep strut down to the nearest
        # supported node (same face preferred), instead of deleting it: deletion merges
        # cells into tall holes. Rescued edges are tracked to cap one rescue per node.
        added = False
        for i in sorted(alive, key=lambda n: nodes[n, 1]):
            if sup[i] or top[i] or i in rescued:
                continue
            best = None
            for j in alive:
                if not sup[j] or nodes[j, 1] >= nodes[i, 1] - eps:
                    continue
                d = dvec(i, j)
                hl = math.hypot(d[0], d[2])
                L = np.linalg.norm(d)
                if -d[1] < tan_min * hl or L > args.rescue_len:
                    continue
                score = L + (0 if abs(d[2]) < eps else 5.0)
                if best is None or score < best[0]:
                    best = (score, j)
            rescued.add(i)
            if best:
                e = (min(i, best[1]), max(i, best[1]))
                E.add(e); cls[e] = "steep"; added = True
                n_rescue += 1
        if added:
            continue
        # top-plane nodes are carried by the rim bridge, not by struts
        bad = {i for i in alive if not sup[i] and not top[i]}
        leaf = {i for i in alive if len(adj[i]) < 2 and not (bottom[i] or top[i])}
        bad |= leaf
        # shallow struts (flat bridges and arches) need both ends supported
        bad_e = {e for e in E if cls[e] in ("short", "long") and not (sup[e[0]] and sup[e[1]])}
        if not bad and not bad_e:
            break
        removed_nodes += len(bad)
        E = {e for e in E if e not in bad_e and e[0] not in bad and e[1] not in bad}
    stats["removed_unsupported_or_leaf_nodes"] = removed_nodes
    stats["rescue_struts"] = n_rescue
    stats["edges_final"] = len(E)
    stats["edge_classes"] = {c: sum(cls[e] == c for e in E) for c in ("plane", "steep", "short", "long")}
    lens = [np.linalg.norm(dvec(i, j)) for i, j in E if cls[(i, j)] != "plane"]
    stats["edge_len_mm"] = {"min": round(min(lens), 1), "median": round(float(np.median(lens)), 1),
                            "max": round(max(lens), 1)}
    # top-rim bridging: horizontal gaps between struts arriving from below
    # (top nodes linked only by in-plane rim edges do not hold the rim up)
    tu = np.sort([nodes[i, 0] for e in E if cls[e] != "plane" for i in e if top[i]])
    gaps = np.diff(np.concatenate([tu, [tu[0] + fr.P]])) if len(tu) else [fr.P]
    stats["rim_support_points"] = int(len(tu))
    stats["rim_max_gap_mm"] = round(float(np.max(gaps)), 1)

    # polylines in (u, z, w), subdivided so struts follow the corner curvature.
    # Plane edges lie inside the band/rim solids; their nodes stay as anchors only.
    uzw_polys = []
    for i, j in E:
        if cls[(i, j)] == "plane":
            continue
        pts = [nodes[i], apex(i, j), nodes[i] + dvec(i, j)] if cls[(i, j)] == "long" \
            else [nodes[i], nodes[i] + dvec(i, j)]
        poly = [pts[0]]
        for a, b in zip(pts[:-1], pts[1:]):
            n = max(1, math.ceil(np.linalg.norm(b - a) / 4.0))
            poly += [a + (b - a) * t for t in np.linspace(0, 1, n + 1)[1:]]
        uzw_polys.append(np.asarray(poly))
    wall_metrics(args, fr, uzw_polys, w0, w1, stats)
    return [fr.map(p[:, 0], p[:, 1], p[:, 2]) for p in uzw_polys]


def wall_metrics(args, fr, uzw_polys, w0, w1, stats, res=0.5):
    """Orientation mix and opening sizes of the unrolled wall, per face and see-through."""
    r = args.strut / 2
    zlo, zhi = args.band, args.height - args.rim
    nu, nz = int(fr.P / res), int((zhi - zlo) / res)
    U = (np.arange(nu) + 0.5) * res
    Z = zlo + (np.arange(nz) + 0.5) * res
    sets = {"outer": [], "inner": [], "see_through": uzw_polys}
    for p in uzw_polys:
        if np.all(np.abs(p[:, 2] - w1) < 1e-3): sets["outer"].append(p)
        elif np.all(np.abs(p[:, 2] - w0) < 1e-3): sets["inner"].append(p)
    sets["through"] = [p for p in uzw_polys
                       if not any(p is q for q in sets["outer"] + sets["inner"])]
    res_ = {}
    for name, polys in sets.items():
        # strut length per in-plane angle from vertical (u-z plane, w ignored)
        hist = np.zeros(4)
        for p in polys:
            for d in np.diff(p, axis=0):
                ang = math.degrees(math.atan2(abs(d[0]), abs(d[1])))
                hist[min(int(ang // 22.5), 3)] += math.hypot(d[0], d[1])
        entry = {"struts": len(polys),
                 "len_share_deg_0-22/22-45/45-67/67-90":
                     [round(x, 2) for x in hist / max(hist.sum(), 1e-9)]}
        if name == "through":
            res_[name] = entry
            continue
        solid = np.zeros((nz, nu), bool)
        for p in polys:
            for a, b in zip(p[:-1], p[1:]):
                for shift in (-fr.P, 0.0, fr.P):
                    a2, b2 = a[:2] + [shift, 0], b[:2] + [shift, 0]
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
        lab, n = ndimage.label(empty)
        for row in range(nz):                     # wall is closed around: join the seam
            a_, b_ = lab[row, 0], lab[row, -1]
            if a_ and b_ and a_ != b_:
                lab[lab == b_] = a_
        edt = ndimage.distance_transform_edt(empty) * res
        ids = [i for i in np.unique(lab) if i]
        incl = np.asarray(ndimage.maximum(edt, lab, ids)) * 2
        objs = ndimage.find_objects(lab)
        heights = np.array([(objs[i - 1][0].stop - objs[i - 1][0].start) * res for i in ids])
        widths = np.array([(objs[i - 1][1].stop - objs[i - 1][1].start) * res for i in ids])
        big = incl >= 2.0                         # ignore slivers between touching struts
        entry.update({"openings": int(big.sum()),
                      "inscribed_dia_mm_p50/p90/max": [round(float(np.percentile(incl[big], q)), 1)
                                                       for q in (50, 90, 100)],
                      "height_mm_p50/p90/max": [round(float(np.percentile(heights[big], q)), 1)
                                                for q in (50, 90, 100)],
                      "height_to_width_p90": round(float(np.percentile(
                          heights[big] / np.maximum(widths[big], res), 90)), 1)})
        res_[name] = entry
    stats["wall_metrics"] = res_


def floor_segments(args, rng):
    h = args.size / 2
    pts = poisson(rng, [-h - 15, -h - 15], [h + 15, h + 15], args.floor_cell, tries=40000)
    vor = Voronoi(pts)
    segs = []
    for rv in vor.ridge_vertices:
        if -1 in rv:
            continue
        a, b = vor.vertices[rv[0]], vor.vertices[rv[1]]
        if rounded_square_sdf(*((a + b) / 2), h, args.corner_r) < 0:
            segs.append((a, b))
    return segs


# ---------- SDF ----------

class Grid:
    def __init__(self, lo, hi, h):
        self.h = h
        self.lo = np.asarray(lo, float)
        self.n = np.ceil((np.asarray(hi) - self.lo) / h).astype(int) + 1
        self.ax = [self.lo[i] + h * np.arange(self.n[i]) for i in range(3)]
        self.D = np.full(self.n, 10.0, np.float32)

    def box(self, pmin, pmax):
        i0 = np.clip(np.floor((np.asarray(pmin) - self.lo) / self.h).astype(int), 0, self.n - 1)
        i1 = np.clip(np.ceil((np.asarray(pmax) - self.lo) / self.h).astype(int) + 1, 1, self.n)
        sl = tuple(slice(a, b) for a, b in zip(i0, i1))
        X = self.ax[0][sl[0]][:, None, None]
        Y = self.ax[1][sl[1]][None, :, None]
        Z = self.ax[2][sl[2]][None, None, :]
        return sl, X, Y, Z


def seg_dist(X, Y, Z, a, b):
    ab = b - a
    px, py, pz = X - a[0], Y - a[1], Z - a[2]
    t = np.clip((px * ab[0] + py * ab[1] + pz * ab[2]) / ab.dot(ab), 0, 1)
    return np.sqrt((px - t * ab[0]) ** 2 + (py - t * ab[1]) ** 2 + (pz - t * ab[2]) ** 2)


def build_sdf(args, fr, polys, fsegs):
    half = args.size / 2
    H = args.height
    r = args.strut / 2
    k = args.fillet
    rt = args.eyelet_wire / 2
    Rt = args.eyelet_hole / 2 + rt
    ez = H + 0.5 + args.eyelet_hole / 2          # hole bottom 0.5 mm above the rim top
    m = 2.0
    # odd origin offset: no grid sample lands exactly on a flat face (z=0, band top, ...),
    # which would make marching cubes emit zero-area slivers
    o = 0.137 * args.voxel
    g = Grid([-half - m + o, -half - m + o, -1.0 + o], [half + m, half + m, ez + Rt + rt + 1.0],
             args.voxel)
    pad = r + k + 2 * args.voxel

    for pl in polys:
        sl, X, Y, Z = g.box(pl.min(0) - pad, pl.max(0) + pad)
        d = None
        for a, b in zip(pl[:-1], pl[1:]):
            di = seg_dist(X, Y, Z, a, b)
            d = di if d is None else np.minimum(d, di)
        g.D[sl] = smin(g.D[sl], d - r, k)

    # 2D fields for floor / band / rim
    X2 = g.ax[0][:, None]; Y2 = g.ax[1][None, :]
    dout = rounded_square_sdf(X2, Y2, half, args.corner_r)
    fl = np.full(dout.shape, 10.0)
    fw = args.floor_strut / 2
    for a, b in fsegs:
        lo_ = np.minimum(a, b) - fw - 1; hi_ = np.maximum(a, b) + fw + 1
        sl, X, Y, _ = g.box([lo_[0], lo_[1], 0], [hi_[0], hi_[1], 0])
        ab = b - a
        px, py = X[:, :, 0] - a[0], Y[:, :, 0] - a[1]
        t = np.clip((px * ab[0] + py * ab[1]) / ab.dot(ab), 0, 1)
        fl[sl[:2]] = np.minimum(fl[sl[:2]], np.hypot(px - t * ab[0], py - t * ab[1]) - fw)
    fl = np.maximum(fl, dout + 1.0)
    pw = np.abs(dout + fr.t / 2)

    zc_b = args.band / 2 - 1
    for iz, z in enumerate(g.ax[2]):
        band = rounded_box2(pw, abs(z - zc_b), fr.t / 2, args.band / 2 + 1, 1.5)
        rim = rounded_box2(pw, abs(z - (H - args.rim / 2)), fr.t / 2, args.rim / 2, 2.0)
        flo = np.maximum(fl, abs(z - (args.floor / 2 - 0.5)) - (args.floor / 2 + 0.5))
        solid = smin(smin(band, rim, 2.0), flo, 2.0)
        g.D[:, :, iz] = np.maximum(smin(g.D[:, :, iz], solid, k), -z)

    # eyelets: torus at each corner, ring plane = vertical + diagonal (in-plane cord pull)
    cc = fr.c + fr.r / math.sqrt(2)
    for sx in (-1, 1):
        for sy in (-1, 1):
            c = np.array([sx * cc, sy * cc, ez])
            e1 = np.array([sx, sy, 0]) / math.sqrt(2)
            e3 = np.array([-sy, sx, 0]) / math.sqrt(2)
            sl, X, Y, Z = g.box(c - Rt - rt - k, c + Rt + rt + k)
            px, py, pz = X - c[0], Y - c[1], Z - c[2]
            a1 = px * e1[0] + py * e1[1]
            a3 = px * e3[0] + py * e3[1]
            d = np.sqrt((np.sqrt(a1 ** 2 + pz ** 2) - Rt) ** 2 + a3 ** 2) - rt
            g.D[sl] = smin(g.D[sl], d, k)
    return g


# ---------- mesh + checks ----------

def clean(mesh):
    # marching cubes and quadric decimation both leave a few zero-area faces, which
    # show up as non-manifold edges / separate bodies
    mesh.update_faces(mesh.nondegenerate_faces())
    mesh.remove_unreferenced_vertices()
    mesh.merge_vertices()
    return mesh


def to_mesh(g):
    v, f, _, _ = marching_cubes(g.D, level=0.0, spacing=(g.h,) * 3)
    mesh = clean(trimesh.Trimesh(v + g.lo, f, process=True))
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


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    fr = WallFrame(args.size / 2, args.corner_r, args.wall)
    stats = {"params": vars(args)}
    t0 = time.time()
    polys = wall_graph(args, fr, rng, stats)
    fsegs = floor_segments(args, rng)
    print(f"graph {time.time() - t0:.1f}s: {len(polys)} wall struts, {len(fsegs)} floor struts")
    if args.graph_only:
        print(json.dumps({k: v for k, v in stats.items() if k != "params"}, indent=1))
        return
    g = build_sdf(args, fr, polys, fsegs)
    print(f"sdf {time.time() - t0:.1f}s grid {g.n.tolist()}")
    mesh = to_mesh(g)
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
