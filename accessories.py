# SPDX-License-Identifier: MIT
# === CONTEXT ===
# Purpose:     Optional accessory parts printed separately from the pot: drip tray as
#              "saucer" (standing) or "clip" (snaps onto the pot, for hanging) and a closed
#              planter (orb: bowl). Decorative, for less demanding plants.
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
#              Clip tray (hanging): flexible arms rise beside the band; their hooks reach past
#              the band's edge rounding onto its flat top, inside lattice gaps found from the
#              pot's strut paths (passed in by orchid_pot.main). The catch face sits above the
#              band's real top (smooth-union bulge + clearance), so a seated pot only touches
#              the feet. The standoff is chosen so the arm stays at <= 2 % strain.
#              Rejected: hooks on the rounded edge only (a detent the load cams open; the orb's
#              leaning band was not gripped at all) and a 1 mm hook with a 6 mm standoff.
#              All verified by mesh intersection with the real pots (seated, lifted, pressed).
#              The same check runs on every tray build when the pot's STL/report sit next to
#              it (fit_checker), after a fingerprint match. Rejected: trusting the same
#              parameters were given (a tray for another seed silently collided).
#              Planter: inner wall = insertion envelope (widest pot outline below each
#              height) + gap, see planter_sdf. Rejected: offsetting the pot surface (a
#              negative taper or the orb could then not be inserted). Orb: bowl to the belly,
#              since a closed onion cannot be assembled; a two-part screw shell was deferred
#              (large threads need clearance trials and hide the pattern completely).
# Constraints: Same packages as potlib.py. Watertightness is a print setting (4-5 walls,
#              ~102-105 % flow, dry PETG), not geometry; see README.
# Invariants:  Tray z = 0 is the bed; the pot bottom sits at floor + standoff (reported as
#              pot_bottom_z_mm).
# Tests:       tests/test_regression.py::test_tray_saucer_checks, ::test_tray_clip_checks,
#              ::test_tray_refuses_foreign_pot_and_from_report, ::test_planter_box_fits_repo_pot.
# === END CONTEXT ===
import copy
import json
import math
import os
import time

import numpy as np
import trimesh

from potlib import Grid, fingerprint, finish, round_isect, smin


def band_top_allowance(args):
    """How far the band's real top rises above --band (the hook's catch face must clear it):
    the band/floor smooth union (k = 2.0 in orchid_pot.build_sdf / orb.build_orb_sdf) adds up
    to k/4 = 0.5 mm, strut fillets about --fillet / 6 (calibrated at the default 1.5: band top
    measured 5.46-5.72 mm for a 5 mm band). Larger fillets raise it; the fit check verifies."""
    return 0.5 + args.fillet / 6


CLIP_STRUT_MARGIN = 1.5        # mm between a clip and a strut surface (covers the strut fillet)
CLIP_LIP_H = 2.0               # hook (lip) height above its catch face: geometry and clearance check
# The band's outer top edge is rounded with r = 1.5 mm in the wall's normal direction
# (rounded_box2 in orchid_pot.build_sdf and orb.build_orb_sdf); horizontally that is
# 1.5 * band_w / wall (orb at 55 deg: 1.83 mm). A hook that only reaches the rounding acts as
# a detent (the load cams it outward), so --clip-hook is measured onto the flat band top.
BAND_EDGE_ROUND = 1.5
CLIP_STRAIN_MAX = 0.02         # PETG cantilever snap-fit guideline 1.5-2 %


# CTX: Read-only adapter from the three pot shapes to what accessories need; the only place
#      that knows the pots' internals. Band parameter: box = arc length u, round = angle.
# Used-by: accessories.py::place_clips, accessories.py::feet_fields, accessories.py::tray_sdf, accessories.py::planter_sdf
# Depends: orchid_pot.py::BoxFrame, orchid_pot.py::ConeFrame, orb.py::OrbShape (attributes)
class TrayGeometry:
    """Read-only view of a pot shape for accessories: outline SDF of the vertical
    projection, the band centre line at z = 0 as a periodic parameter (box: arc length u
    on the mid-curve, round shapes: angle theta), band width at z = 0 and the outward lean
    of the band per mm of height (lean: clamped >= 0, for the foot stops, which only need to
    clear an outward lean; lean_signed: real, for the hook's reach onto the band top)."""

    def __init__(self, args, geom):
        self.shape, self.geom, t = args.shape, geom, args.wall
        if args.shape == "box":
            self.outline = geom.footprint_sdf                    # vertical walls
            self.extent = geom.half_out
            self.period, self.arc = geom.P, 1.0
            self.foot_params = [k * geom.quad + d for k in range(4)
                                for d in (geom.c, geom.side + geom.r * math.pi / 4)]
            self.band_w, self.lean, self.lean_signed = t, 0.0, 0.0
        else:
            if args.shape == "cylinder":
                r_proj = max(float(geom.r_out(0.0)), float(geom.r_out(geom.H)))
                self.r_mid0, self.band_w = float(geom.r_mid(0.0)), t / geom.cos_t
                self.lean_signed = geom.tan_t
                self.lean = max(geom.tan_t, 0.0)
            else:                                                # orb
                r_proj = geom.a_out
                self.r_mid0 = geom.r_floor
                al = math.radians(args.min_wall_angle)
                self.band_w, self.lean = t / math.sin(al), 1 / math.tan(al)
                self.lean_signed = self.lean
            self.outline = lambda x, y, R=r_proj: np.hypot(x, y) - R
            self.extent = r_proj
            self.period, self.arc = 2 * math.pi, self.r_mid0
            self.foot_params = [k * math.pi / 4 for k in range(8)]

    def band_at(self, params):
        """Band centre-line points (x, y) at z = 0 and outward unit normals."""
        q = np.asarray(params, float)
        if self.shape == "box":
            t = self.geom.t
            mid = self.geom.map(q, np.zeros(len(q)), np.full(len(q), t / 2))
            out = self.geom.map(q, np.zeros(len(q)), np.full(len(q), t))
            n = out[:, :2] - mid[:, :2]
            return mid[:, :2], n / np.linalg.norm(n, axis=1, keepdims=True)
        n = np.column_stack([np.cos(q), np.sin(q)])
        return self.r_mid0 * n, n

    def clip_targets(self, n):
        if self.shape == "box":                                  # side centres
            return [k * self.period / n + self.geom.c for k in range(n)]
        return [k * self.period / n for k in range(n)]


def strut_chart(tg, struts, step=0.25):
    """Strut centre lines as points (arc length along the band line, pot height z),
    resampled every `step` mm; box/cylinder come in their (u, z) chart, the orb in 3D."""
    out = []
    for p in struts:
        p = np.asarray(p, float)
        seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
        pts = np.vstack([a_ + (b_ - a_) * np.linspace(0, 1, max(1, math.ceil(L / step)), endpoint=False)[:, None]
                         for a_, b_, L in zip(p[:-1], p[1:], seg)] + [p[-1:]])
        if p.shape[1] == 2:
            q = pts[:, 0] % tg.period if tg.shape == "box" else (pts[:, 0] / tg.geom.r_ref) % tg.period
            z = pts[:, 1]
        else:
            q, z = np.arctan2(pts[:, 1], pts[:, 0]) % tg.period, pts[:, 2]
        out.append(np.column_stack([q * tg.arc, z]))
    return np.vstack(out)


# CTX: Clip positions from the pot's own strut paths (same seed); must agree with the hook
#      region tray_sdf builds (catch height band_top_allowance + clearance, CLIP_LIP_H).
# Used-by: accessories.py::run | Depends: orchid_pot.py::wall_graph, orb.py::orb_graph (paths)
def place_clips(tg, args, struts, stats):
    """Clip parameters. A position is free if the clip's hook region (its width along the
    band x the hook's height range) keeps CLIP_STRUT_MARGIN from every strut surface, i.e.
    strut/2 + margin from the strut centre lines, checked over the whole hook height against
    the strut paths themselves (works without the pot STL). Each target (box: side centres,
    round: even angles) takes the nearest free position; clips stay 2 clip widths apart.
    Rejected: a gap width measured at one height (struts at 35 deg move 1.4 mm sideways per
    mm of height, so they could cut into the hook)."""
    n = args.clips if args.clips is not None else (4 if tg.shape == "box" else 3)
    z_lo = args.band + band_top_allowance(args) + args.clip_clearance
    z_hi = z_lo + CLIP_LIP_H
    need = args.strut / 2 + CLIP_STRUT_MARGIN
    pts = strut_chart(tg, struts)
    pts = pts[(pts[:, 1] > z_lo - need - 1) & (pts[:, 1] < z_hi + need + 1)]
    L = tg.period * tg.arc                                    # band line length
    cand = np.arange(0, L, 0.5)
    ds = np.abs((pts[None, :, 0] - cand[:, None] + L / 2) % L - L / 2)
    dx = np.maximum(ds - args.clip_width / 2, 0)
    dz = np.maximum(np.maximum(pts[None, :, 1] - z_hi, z_lo - pts[None, :, 1]), 0)
    clear = np.hypot(dx, dz).min(axis=1) - args.strut / 2     # to the strut surface
    free = clear >= CLIP_STRUT_MARGIN
    params, clears = [], []
    for t in tg.clip_targets(n):
        ok = free.copy()
        for c in params:
            ok &= np.abs((cand - c + L / 2) % L - L / 2) >= 2 * args.clip_width
        if not ok.any():
            raise SystemExit("error: not enough free lattice room above the band for the clips")
        d = np.abs((cand - t * tg.arc + L / 2) % L - L / 2)
        k = int(np.argmin(np.where(ok, d, np.inf)))
        params.append(cand[k]); clears.append(round(float(clear[k]), 2))
    stats["clip_strut_clearance_mm"] = clears
    return [c / tg.arc for c in params]


FOOT_HW, FOOT_STOP_H = 5.0, 3.0   # foot half width along the band; stop height above the foot


def feet_fields(tg, X2, Y2, params):
    """2D fields (foot, stop) of the pot's feet under its band centre line at the given
    band parameters. The stop sits outside the band, cleared by the band's lean, with 1 mm
    side clearance (0.5 mm measured only ~0.4 mm real play because of the band's fillets,
    too tight with the over-extrusion of watertight print settings)."""
    stop_t, clear = 2.0, 1.0
    out = []
    pts, nrm = tg.band_at(params)
    for p, n in zip(pts, nrm):
        r_ = (X2 - p[0]) * n[0] + (Y2 - p[1]) * n[1]
        s_ = (X2 - p[0]) * -n[1] + (Y2 - p[1]) * n[0]
        foot2 = np.maximum(np.abs(r_) - tg.band_w / 2, np.abs(s_) - FOOT_HW)
        r_stop = tg.band_w / 2 + clear + tg.lean * FOOT_STOP_H + stop_t / 2
        stop2 = np.maximum(np.abs(r_ - r_stop) - stop_t / 2, np.abs(s_) - FOOT_HW)
        out.append((foot2, stop2))
    return out


def clip_geometry(args, tg):
    """Radial clip geometry relative to the band centre line at z = 0:
    r_top  band outer face at the band top (signed lean: narrower for taper < 0),
    r_wide widest band face below the hook (top for an outward lean, bottom for taper < 0),
    reach  hook reach past r_top: edge rounding (horizontal) + --clip-hook onto the flat top,
    r_tip  hook tip, r_arm arm inner face (clears r_wide), defl snap deflection r_wide - r_tip.
    Rejected: one radius for both (taper -8: the arm then grazed the band's lower edge)."""
    r_top = tg.band_w / 2 + tg.lean_signed * args.band
    r_wide = tg.band_w / 2 + max(tg.lean_signed * args.band, 0.0)
    reach = BAND_EDGE_ROUND * tg.band_w / args.wall + args.clip_hook
    r_tip = r_top - reach
    return {"r_top": r_top, "r_wide": r_wide, "reach": reach, "r_tip": r_tip,
            "r_arm": r_wide + args.clip_clearance, "defl": r_wide - r_tip}


def clip_standoff(args, tg):
    """Standoff giving a clip arm long enough for CLIP_STRAIN_MAX at the snap deflection:
    strain = 1.5 t y / L^2 and L = standoff + band + band_top_allowance + clearance."""
    y = clip_geometry(args, tg)["defl"]
    L = math.sqrt(1.5 * args.clip_thickness * y / CLIP_STRAIN_MAX)
    return max(6.0, round(L - args.band - band_top_allowance(args) - args.clip_clearance + 0.05, 1))


def tray_sdf(args, tg, stats, clip_params=()):
    gap, wall, ft, so = args.tray_gap, args.tray_wall, args.tray_floor, args.standoff
    h_tray = ft + so + args.tray_rise
    ct, cc, cw = args.clip_thickness, args.clip_clearance, args.clip_width
    cg = clip_geometry(args, tg)
    cpts, cnrm = tg.band_at(clip_params) if len(clip_params) else (np.zeros((0, 2)), np.zeros((0, 2)))
    # the arm must be able to flex outward by the snap deflection without touching the tray wall
    need = cg["r_arm"] + ct + cg["defl"] + 0.5
    rs = np.arange(0, 60, 0.05)
    for p, n in zip(cpts, cnrm):
        o_ = tg.outline(p[0] + rs * n[0], p[1] + rs * n[1])
        r_wall = rs[np.argmax(o_ >= gap)]
        gap = max(gap, gap + need - r_wall)
    o = 0.137 * args.voxel
    E = tg.extent + gap + wall + 2.0
    # The band top is not flat at nominal --band: the pot's band/floor smooth union (k = 2.0
    # in orchid_pot.build_sdf / orb.build_orb_sdf) lifts it by up to k/4 = 0.5 mm, and strut
    # fillets add more. Measured under the hooks at defaults: 5.46-5.72 mm for a 5 mm band.
    z_catch = ft + so + args.band + band_top_allowance(args) + cc
    lip_h = CLIP_LIP_H
    g = Grid([-E + o, -E + o, -1.0 + o], [E, E, max(h_tray, z_catch + lip_h) + 1.0], args.voxel)
    X2, Y2 = np.meshgrid(g.ax[0], g.ax[1], indexing="ij")
    d = tg.outline(X2, Y2)                                   # < 0 inside the pot projection
    wall2 = np.abs(d - gap - wall / 2) - wall / 2             # tray wall ring
    base2 = d - gap - wall                                    # tray footprint

    # feet under the band centre line (+ outer stop that clears the band's lean); feet next
    # to a clip are dropped (the clip arm locates the pot there)
    keep = [q for q in tg.foot_params
            if all(abs(((q - c) + tg.period / 2) % tg.period - tg.period / 2) * tg.arc > FOOT_HW + cw / 2 + 2
                   for c in clip_params)]
    feet = feet_fields(tg, X2, Y2, keep)
    stop_h = FOOT_STOP_H
    local = lambda p, n: ((X2 - p[0]) * n[0] + (Y2 - p[1]) * n[1], (X2 - p[0]) * -n[1] + (Y2 - p[1]) * n[0])
    # clips: vertical arm beside the band, hook over the band top into a lattice gap. The
    # hook's top is chamfered 45 deg (lead-in when the pot is pushed down); its underside is
    # the catch face, a 1.3 mm overhang that prints without support.
    clips = []
    r_ai, r_tip = cg["r_arm"], cg["r_tip"]
    for p, n in zip(cpts, cnrm):
        r_, s_ = local(p, n)
        clips.append((np.maximum(np.abs(r_ - (r_ai + ct / 2)) - ct / 2, np.abs(s_) - cw / 2), r_, np.abs(s_) - cw / 2))

    for iz, z in enumerate(g.ax[2]):
        floor = np.maximum(base2, abs(z - ft / 2) - ft / 2)
        tw = round_isect(wall2, np.maximum(z - h_tray, -z), 0.8)
        D = smin(floor, tw, 1.0)
        for foot2, stop2 in feet:
            f = round_isect(foot2, np.maximum(z - (ft + so), -z), 0.5)
            s = round_isect(stop2, np.maximum(z - (ft + so + stop_h), -z), 0.5)
            D = smin(D, np.minimum(f, s), 0.8)
        for arm2, r_, sw in clips:
            arm = round_isect(arm2, np.maximum(z - (z_catch + lip_h), (ft - 0.5) - z), 0.4)
            lip = np.maximum.reduce([r_tip - r_, r_ - (r_ai + ct), sw, np.full_like(r_, z_catch - z),
                                     np.full_like(r_, z - (z_catch + lip_h)),
                                     ((z - z_catch) - 0.3 - (r_ - r_tip)) / math.sqrt(2)])
            D = smin(D, np.minimum(arm, lip), 0.5)
        g.D[:, :, iz] = np.maximum(D, -z)

    # water capacity: cavity above the floor up to 1 mm below the pot bottom, minus feet/clips
    cav = (d < gap) & np.all([(f > 0) & (s > 0) for f, s in feet] + [a > 0 for a, _, _ in clips], axis=0)
    depth = so - 1.0
    stats["tray"] = {"type": args.tray, "outer_extent_mm": round(2 * (tg.extent + gap + wall), 1),
                     "height_mm": round(h_tray, 1), "pot_bottom_z_mm": round(ft + so, 1),
                     "water_depth_mm": depth,
                     "water_capacity_ml": round(float(cav.sum()) * g.h ** 2 * depth / 1000, 0),
                     "feet": len(feet), "tray_gap_used_mm": round(gap, 2)}
    if len(clip_params):
        L = z_catch - ft                                      # free arm length (floor top -> catch)
        stats["tray"].update({"clips": len(clip_params), "clip_catch_z_mm": round(z_catch, 2),
                              "clip_arm_length_mm": round(L, 1),
                              "clip_hook_reach_mm": round(cg["reach"], 2), "clip_hook_on_flat_mm": args.clip_hook,
                              "clip_snap_deflection_mm": round(cg["defl"], 2),
                              # cantilever snap-fit root strain 1.5 t y / L^2 (PETG guideline <= 2 %)
                              "clip_strain_pct": round(150 * ct * cg["defl"] / L ** 2, 2)})
    return g


# ---------- planter ----------

def pot_outer_radius(args, geom, h):
    """Round shapes: outer radius of the pot at pot height h (array, clamped to the body)."""
    h = np.asarray(h, float)
    if args.shape == "cylinder":
        return geom.r_out(np.clip(h, 0.0, geom.H))
    # orb: the outer face is the mid ellipse offset by t/2 along its normal
    a, c, zc, t = geom.a, geom.c, geom.zc, geom.t
    ph = np.linspace(geom.phi_f - 0.2, math.pi / 2, 4000)
    nr, nz = c * np.sin(ph), -a * np.cos(ph)
    nn = np.hypot(nr, nz)
    rho = a * np.sin(ph) + t / 2 * nr / nn
    zo = zc - c * np.cos(ph) + t / 2 * nz / nn
    return np.interp(np.clip(h, 0.0, zc), zo, rho)


def planter_sdf(args, tg, geom, stats):
    """Closed planter (box, cylinder) or bowl (orb, up to the belly). Its inner wall follows
    the pot's INSERTION ENVELOPE plus --planter-gap: at each height the widest pot outline
    of everything below it (the pot slides in from above, so only lower parts pass a level).
    That one rule gives straight walls for the box, the taper for taper > 0, a straight
    cylinder for taper < 0 (the pot's floor is its widest part) and a bowl for the orb."""
    gap, wall, fl = args.planter_gap, args.planter_wall, args.planter_floor
    so = 6.0 if args.standoff is None else args.standoff
    zb = fl + so                                              # pot bottom
    if tg.shape == "orb":
        top_pot = geom.zc                                     # bowl up to the belly
    else:
        top_pot = args.height                                 # flush with the pot rim
    z_top = zb + top_pot
    ring = args.planter_eyelets > 0
    rt, Rt = args.eyelet_wire / 2, args.eyelet_hole / 2 + args.eyelet_wire / 2
    ez = z_top + 0.5 + args.eyelet_hole / 2
    if tg.shape == "box":
        env_at = lambda z: env_box
        r_top = None
    else:
        hs = np.linspace(0.0, top_pot, 2000)
        r_env = np.maximum.accumulate(pot_outer_radius(args, geom, hs))
        R_env = lambda z: np.interp(np.clip(z - zb, 0.0, top_pot), hs, r_env)
        dR = lambda z: (float(R_env(z + 0.2)) - float(R_env(z - 0.2))) / 0.4
        env_at = lambda z: rho2 - float(R_env(z))
        r_top = float(R_env(z_top))
    E = tg.extent + gap + wall + 3.0
    o = 0.137 * args.voxel
    g = Grid([-E + o, -E + o, -1.0 + o], [E, E, (ez + Rt + rt + 1.0) if ring else z_top + 1.0], args.voxel)
    X2, Y2 = np.meshgrid(g.ax[0], g.ax[1], indexing="ij")
    env_box = tg.outline(X2, Y2) if tg.shape == "box" else None
    rho2 = np.hypot(X2, Y2)
    th2 = np.arctan2(Y2, X2)
    feet = feet_fields(tg, X2, Y2, tg.foot_params)

    # centring ribs (box: side centres, round: 45 deg off the feet), 3 mm wide
    nrib = args.planter_ribs
    if tg.shape == "box":
        rp, rn = tg.band_at(tg.clip_targets(nrib)) if nrib else (np.zeros((0, 2)), np.zeros((0, 2)))
        ribs_s = [np.abs((X2 - p[0]) * -n[1] + (Y2 - p[1]) * n[0]) - 1.5 for p, n in zip(rp, rn)]
    else:
        ribs_s = [np.abs(rho2 * np.sin(th2 - (math.pi / 4 + 2 * math.pi * k / max(nrib, 1)))) - 1.5
                  + np.where(np.cos(th2 - (math.pi / 4 + 2 * math.pi * k / max(nrib, 1))) > 0, 0, 1e3)
                  for k in range(nrib)]
    clr = args.planter_rib_clearance

    for iz, z in enumerate(g.ax[2]):
        env = env_at(z)
        wh = wall * (math.sqrt(1 + dR(z) ** 2) if tg.shape != "box" else 1.0)  # normal thickness = wall
        wall2 = np.abs(env - gap - wh / 2) - wh / 2
        D = round_isect(wall2, np.maximum(z - z_top, -z), 1.0)
        env0 = env_at(0.0)
        D = smin(D, np.maximum(env0 - gap - wall, abs(z - fl / 2) - fl / 2), 1.0)
        for foot2, stop2 in feet:
            f = round_isect(foot2, np.maximum(z - zb, -z), 0.5)
            st = round_isect(stop2, np.maximum(z - (zb + FOOT_STOP_H), -z), 0.5)
            D = smin(D, np.minimum(f, st), 0.8)
        # rib: from the wall in to the envelope + clearance; 45 deg lead-in chamfer at the top
        cham = max(0.0, z - (z_top - 2.0 - (gap - clr)))
        for rs in ribs_s:
            rib = np.maximum.reduce([clr + cham - env, env - (gap + wh / 2), rs,
                                     np.full_like(env, fl - 0.5 - z), np.full_like(env, z - (z_top - 2.0))])
            D = smin(D, rib, 0.5)
        g.D[:, :, iz] = np.maximum(D, -z)

    if ring:                                                   # eyelets on the planter rim
        if tg.shape == "box":
            cp, cn = tg.band_at([tg.geom.side + tg.geom.r * math.pi / 4 + k * tg.geom.quad
                                 for k in range(4)][::max(1, 4 // args.planter_eyelets)])
            cen = cp + cn * (args.wall / 2 + gap + wall / 2)
        else:
            th = 2 * math.pi * np.arange(args.planter_eyelets) / args.planter_eyelets
            cn = np.column_stack([np.cos(th), np.sin(th)])
            cen = cn * (r_top + gap + wall / 2)
        rt_w = max(wall, 4.0)
        for c2, n in zip(cen, cn):
            c = np.array([c2[0], c2[1], ez])
            sl, X, Y, Z = g.box(c - Rt - rt - rt_w, c + Rt + rt + rt_w)
            a1 = (X - c[0]) * -n[1] + (Y - c[1]) * n[0]
            ar = (X - c[0]) * n[0] + (Y - c[1]) * n[1]
            ringd = round_isect(np.abs(np.hypot(a1, Z - c[2]) - Rt) - rt, np.abs(ar) - rt_w / 2, 0.8)
            g.D[sl] = smin(g.D[sl], ringd, 1.5)

    env0 = env_at(0.0)
    cav = (env0 < gap) & np.all([(f > 0) & (s_ > 0) for f, s_ in feet], axis=0)
    depth = so - 1.0
    stats["planter"] = {"type": "bowl" if tg.shape == "orb" else "closed",
                        "outer_extent_mm": round(2 * (tg.extent + gap + wall), 1) if tg.shape == "box"
                        else round(2 * (float(R_env(z_top)) + gap + wall), 1),
                        "height_mm": round(z_top, 1), "pot_bottom_z_mm": round(zb, 1),
                        "reservoir_depth_mm": depth,
                        "reservoir_ml": round(float(cav.sum()) * g.h ** 2 * depth / 1000, 0),
                        "ribs": nrib, "eyelets": args.planter_eyelets}
    return g


# CTX: The automatic proof that a part fits its pot; refuses a pot with another fingerprint.
# Depends: potlib.py::fingerprint (written into the pot report by build_pot / orb.run)
# Invariant: zb is the exact pot-bottom height; the rounded report value caused false failures.
def fit_checker(args, kind, fp, zb, clip_xy=()):
    """extra_checks for finish(): if the pot's STL and report are found (next to the part,
    or where --from-report points), verify the fingerprints and the real fit by mesh
    intersection, with the pot at the exact height zb (not the rounded report value).
    Seated: no overlap. Clip: lifted 1.5 mm, every clip must catch (intersection pieces are
    assigned to the nearest clip). Planter: the pot also slides in freely (10/30/60 mm
    above). Missing pot files: skipped, with a note in the report."""
    prefix = getattr(args, "pot_out", None) or args.out
    pot_stl, pot_rep = prefix + ".stl", prefix + "_report.json"

    def check(part, stats):
        if not (os.path.exists(pot_stl) and os.path.exists(pot_rep)):
            stats["fit_check"] = {"ok": None, "note": f"skipped: no pot STL/report at {prefix}.*"}
            return
        with open(pot_rep) as fh:
            pfp = json.load(fh).get("fingerprint")
        if pfp is None:
            stats["fit_check"] = {"ok": False, "note": f"{pot_rep} has no fingerprint (pot built "
                                                       f"before v0.7.2): rebuild the pot"}
            return
        bad = [k for k in fp if pfp.get(k) != fp[k]]
        if bad:
            stats["fit_check"] = {"ok": False, "note": f"pot does not match this part ({', '.join(bad)} "
                                                       f"fingerprint differs); use --from-report {pot_rep}"}
            return
        pot = trimesh.load(pot_stl)

        def overlap(dz):
            p = pot.copy(); p.apply_translation([0, 0, zb + dz])
            i = trimesh.boolean.intersection([p, part], engine="manifold")
            return (float(i.volume), i.split(only_watertight=False)) if len(i.faces) else (0.0, [])
        seated = overlap(0.02)[0]
        res = {"seated_overlap_mm3": round(seated, 2)}
        ok = seated <= 0.05
        if kind == "planter":
            ins = [overlap(dz)[0] for dz in (10.0, 30.0, 60.0)]
            res["insertion_overlap_mm3_at_10/30/60mm"] = [round(v, 2) for v in ins]
            ok = ok and max(ins) <= 0.05
        elif kind == "clip":
            cxy = np.asarray(clip_xy)
            hit = {int(np.argmin(np.linalg.norm(cxy - q.centroid[:2], axis=1))) for q in overlap(1.5)[1]}
            res.update({"lifted_1.5mm_clips_engaged": len(hit), "clips": len(cxy)})
            ok = ok and len(hit) == len(cxy)
        stats["fit_check"] = {"ok": ok, **res}
    return check


def report_fit(stats, out, part):
    """A part that fails the fit check is not left on disk as a printable STL."""
    fc = stats["fit_check"]
    if fc["ok"] is False:
        if os.path.exists(out + ".stl"):
            os.remove(out + ".stl")
        raise SystemExit(f"error: {part} does not fit the pot ({out}.stl removed, report kept): {fc}")
    if fc["ok"] is None and stats.get("tray", {}).get("clips"):
        print("WARNING: clip tray not checked against the pot (pot files not found); placement "
              "comes from the strut paths, the band-top height is calibrated, not measured")
    print("fit check:", fc)


def _with_height(a, args):
    """planter_sdf needs the pot's real height (a.height is parked at 1e6 for potlib)."""
    b = copy.copy(a)
    b.height = getattr(args, "height", 0.0)
    return b


def run(args, parts, geom, struts=None):
    t0 = time.time()
    tg = TrayGeometry(args, geom)
    if "tray" in parts:
        a = copy.copy(args)
        # height far above the part: potlib's rim/eyelet zones then never apply to a tray
        a.out = args.out + ("_tray" if args.tray == "saucer" else "_tray_clip")
        a.height, a.eyelets = 1e6, 0
        if a.standoff is None:
            a.standoff = 6.0 if a.tray == "saucer" else clip_standoff(a, tg)
        fp = fingerprint(args, struts if args.tray == "clip" else None)
        stats = {"params": vars(a), "part": "tray", "fingerprint": fp}
        cp = place_clips(tg, args, struts, stats) if args.tray == "clip" else ()
        g = tray_sdf(a, tg, stats, cp)
        print(f"tray sdf {time.time() - t0:.1f}s grid {g.n.tolist()}")
        pts, nrm = tg.band_at(cp) if len(cp) else (np.zeros((0, 2)), np.zeros((0, 2)))
        clip_xy = pts + nrm * (tg.band_w / 2 + args.clip_clearance)
        finish(g, a, stats, t0, extra_checks=fit_checker(args, "clip" if len(cp) else "saucer", fp,
                                                         a.tray_floor + a.standoff, clip_xy))
        report_fit(stats, a.out, "tray")
    if "planter" in parts:
        a = copy.copy(args)
        a.out, a.height, a.eyelets = args.out + "_planter", 1e6, 0
        fp = fingerprint(args)
        stats = {"params": vars(a), "part": "planter", "fingerprint": fp}
        g = planter_sdf(_with_height(a, args), tg, geom, stats)
        print(f"planter sdf {time.time() - t0:.1f}s grid {g.n.tolist()}")
        zb = a.planter_floor + (6.0 if a.standoff is None else a.standoff)
        finish(g, a, stats, t0, extra_checks=fit_checker(args, "planter", fp, zb))
        report_fit(stats, a.out, "planter")
