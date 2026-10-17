# SPDX-License-Identifier: MIT
# === CONTEXT ===
# Purpose:     Optional accessory parts printed separately from the pot: drip tray ("saucer"
#              now; "clip" and the planter follow). Decorative, for less demanding plants.
#              Called from orchid_pot.py::main for --parts tray.
# Decisions:   Own module, so the pot pipelines stay untouched. It reads the pot geometry
#              read-only through small adapters (outline, band, wall lean) instead of adding
#              methods to the pot's frames.
#              The tray outline is the pot's vertical projection, i.e. the widest outline at
#              any height, plus a gap. Drips can leave an outward-leaning wall anywhere
#              (taper > 0, orb belly), so the floor footprint alone is too small.
#              The pot stands on feet under its solid bottom band, not under the floor
#              lattice, so drainage stays open. The water line stays 1 mm below the pot
#              bottom: the floor is never in standing water and air reaches it from below.
#              Each foot has an outer stop (cleared by the band's lean), so the pot cannot
#              slide in a tray that is wider than its band.
# Constraints: Same packages as potlib.py. Watertightness is a print setting (4-5 walls,
#              ~102-105 % flow, dry PETG), not geometry; see README.
# Invariants:  Tray z = 0 is the bed; the pot bottom sits at floor + standoff (reported as
#              pot_bottom_z_mm).
# Tests:       tests/test_regression.py::test_tray_saucer_checks.
# === END CONTEXT ===
import copy
import math
import time

import numpy as np

from potlib import Grid, finish, round_isect, smin


class TrayGeometry:
    """Read-only view of a pot shape for accessories: outline SDF of the vertical
    projection, band centre line at z = 0 (points + outward normals), band width at z = 0
    and the outward lean of the band per mm of height."""

    def __init__(self, args, geom):
        self.shape = args.shape
        t = args.wall
        if args.shape == "box":
            self.outline = geom.footprint_sdf                    # vertical walls
            self.extent = geom.half_out
            us = [k * geom.quad + d for k in range(4) for d in (geom.c, geom.side + geom.r * math.pi / 4)]
            mid = geom.map(np.array(us), np.zeros(len(us)), np.full(len(us), t / 2))
            out = geom.map(np.array(us), np.zeros(len(us)), np.full(len(us), t))
            self.band_pts = mid[:, :2]
            n = out[:, :2] - mid[:, :2]
            self.band_nrm = n / np.linalg.norm(n, axis=1, keepdims=True)
            self.band_w, self.lean = t, 0.0
        else:
            if args.shape == "cylinder":
                r_proj = max(float(geom.r_out(0.0)), float(geom.r_out(geom.H)))
                r_mid0, self.band_w = float(geom.r_mid(0.0)), t / geom.cos_t
                self.lean = max(geom.tan_t, 0.0)
            else:                                                # orb
                r_proj = geom.a_out
                r_mid0 = geom.r_floor
                al = math.radians(args.min_wall_angle)
                self.band_w, self.lean = t / math.sin(al), 1 / math.tan(al)
            self.outline = lambda x, y, R=r_proj: np.hypot(x, y) - R
            self.extent = r_proj
            th = np.arange(8) * 2 * math.pi / 8
            self.band_nrm = np.column_stack([np.cos(th), np.sin(th)])
            self.band_pts = r_mid0 * self.band_nrm


def tray_sdf(args, tg, stats):
    gap, wall, ft, so = args.tray_gap, args.tray_wall, args.tray_floor, args.standoff
    h_tray = ft + so + args.tray_rise
    o = 0.137 * args.voxel
    E = tg.extent + gap + wall + 2.0
    g = Grid([-E + o, -E + o, -1.0 + o], [E, E, h_tray + 1.0], args.voxel)
    X2, Y2 = np.meshgrid(g.ax[0], g.ax[1], indexing="ij")
    d = tg.outline(X2, Y2)                                   # < 0 inside the pot projection
    wall2 = np.abs(d - gap - wall / 2) - wall / 2             # tray wall ring
    base2 = d - gap - wall                                    # tray footprint

    # feet under the band centre line (+ outer stop that clears the band's lean)
    hw_t = 5.0                                                # half width along the band
    # 1 mm side clearance: 0.5 mm measured only ~0.4 mm real play (fillets on the band), too
    # tight with the over-extrusion of watertight print settings
    stop_h, stop_t, clear = 3.0, 2.0, 1.0
    feet = []
    for p, n in zip(tg.band_pts, tg.band_nrm):
        tvec = np.array([-n[1], n[0]])
        r_ = (X2 - p[0]) * n[0] + (Y2 - p[1]) * n[1]
        s_ = (X2 - p[0]) * tvec[0] + (Y2 - p[1]) * tvec[1]
        foot2 = np.maximum(np.abs(r_) - tg.band_w / 2, np.abs(s_) - hw_t)
        r_stop = tg.band_w / 2 + clear + tg.lean * stop_h + stop_t / 2
        stop2 = np.maximum(np.abs(r_ - r_stop) - stop_t / 2, np.abs(s_) - hw_t)
        feet.append((foot2, stop2))

    for iz, z in enumerate(g.ax[2]):
        floor = np.maximum(base2, abs(z - ft / 2) - ft / 2)
        tw = round_isect(wall2, np.maximum(z - h_tray, -z), 0.8)
        D = smin(floor, tw, 1.0)
        for foot2, stop2 in feet:
            f = round_isect(foot2, np.maximum(z - (ft + so), -z), 0.5)
            s = round_isect(stop2, np.maximum(z - (ft + so + stop_h), -z), 0.5)
            D = smin(D, np.minimum(f, s), 0.8)
        g.D[:, :, iz] = np.maximum(D, -z)

    # water capacity: cavity above the floor up to 1 mm below the pot bottom, minus feet
    cav = (d < gap) & np.all([(f > 0) & (s > 0) for f, s in feet], axis=0)
    depth = so - 1.0
    stats["tray"] = {"type": args.tray, "outer_extent_mm": round(2 * (tg.extent + gap + wall), 1),
                     "height_mm": round(h_tray, 1), "pot_bottom_z_mm": round(ft + so, 1),
                     "water_depth_mm": depth,
                     "water_capacity_ml": round(float(cav.sum()) * g.h ** 2 * depth / 1000, 0),
                     "feet": len(feet)}
    return g


def run(args, parts, geom):
    t0 = time.time()
    tg = TrayGeometry(args, geom)
    if "tray" in parts:
        a = copy.copy(args)
        # height far above the part: potlib's rim/eyelet zones then never apply to a tray
        a.out, a.height, a.eyelets = args.out + "_tray", 1e6, 0
        stats = {"params": vars(a), "part": "tray"}
        g = tray_sdf(a, tg, stats)
        print(f"tray sdf {time.time() - t0:.1f}s grid {g.n.tolist()}")
        finish(g, a, stats, t0)
