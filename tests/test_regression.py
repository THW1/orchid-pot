# SPDX-License-Identifier: MIT
# === CONTEXT ===
# Purpose:     Guard the box and cylinder outputs against unintended change, e.g. when the
#              experimental orb shape or shared helpers are added.
# Decisions:   Compares the STL byte hash and the full JSON report (minus `params`, which
#              holds the output path) for seed 7. The pipeline is deterministic, so any
#              difference is a real change. Goldens last regenerated for the intended v0.6.1
#              gusset change (later only report fields were added, STLs unchanged).
# Constraints: About 30 s per case (full mesh build).
# Related:     tests/golden/*_report.json, tests/golden/stl.sha256; regenerate them only
#              for an intended output change, and say so in the commit message.
# === END CONTEXT ===
import hashlib
import json
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
GOLDEN = ROOT / "tests" / "golden"
CASES = {"box": ["box"], "cylinder": ["cylinder"], "cylinder_t8": ["cylinder", "--taper", "8"]}
HASHES = dict(line.split()[::-1] for line in (GOLDEN / "stl.sha256").read_text().splitlines())


@pytest.mark.parametrize("name", CASES)
def test_output_unchanged(name, tmp_path):
    out = tmp_path / name
    subprocess.run([sys.executable, str(ROOT / "orchid_pot.py"), *CASES[name], "--out", str(out)],
                   check=True, capture_output=True)
    report = json.loads((tmp_path / f"{name}_report.json").read_text())
    golden = json.loads((GOLDEN / f"{name}_report.json").read_text())
    report.pop("params"), golden.pop("params")
    assert report == golden
    assert hashlib.sha256((tmp_path / f"{name}.stl").read_bytes()).hexdigest() == HASHES[name]


def test_orb_graph_smoke():
    """Experimental orb: no golden yet; the strut graph must need no node removal and keep
    all non-adjacent struts at least half a strut width apart."""
    r = subprocess.run([sys.executable, str(ROOT / "orchid_pot.py"), "orb", "--graph-only"],
                       check=True, capture_output=True, text=True).stdout
    stats = json.loads(r[r.index("{"):])
    assert stats["removed_unsupported_or_leaf_nodes"] == 0
    assert stats["strut_clearance_violations"] == 0
    assert stats["wall_metrics"]["inscribed_dia_mm_p50/p90/max"][2] < 35


@pytest.mark.parametrize("shape", ["box", "orb"])
def test_tray_saucer_checks(shape, tmp_path):
    """Accessory tray (no golden yet): one watertight body, no islands, holds water."""
    out = tmp_path / shape
    subprocess.run([sys.executable, str(ROOT / "orchid_pot.py"), shape, "--parts", "tray", "--out", str(out)],
                   check=True, capture_output=True)
    r = json.loads((tmp_path / f"{shape}_tray_report.json").read_text())
    assert r["watertight"] and r["bodies"] == 1 and r["unsupported_islands"] == 0
    assert r["tray"]["water_capacity_ml"] >= 80


@pytest.mark.parametrize("shape,clips", [("box", 4), ("orb", 3)])
def test_tray_clip_checks(shape, clips, tmp_path):
    """Clip tray (no golden yet): valid STL after the round trip, all clips placed, arm
    strain within the PETG snap-fit guideline."""
    out = tmp_path / shape
    subprocess.run([sys.executable, str(ROOT / "orchid_pot.py"), shape, "--parts", "tray", "--tray", "clip",
                    "--out", str(out)], check=True, capture_output=True)
    r = json.loads((tmp_path / f"{shape}_tray_clip_report.json").read_text())
    assert r["stl_roundtrip_watertight"] and r["bodies"] == 1 and r["unsupported_islands"] == 0
    assert r["tray"]["clips"] == clips and r["tray"]["clip_strain_pct"] <= 2.0


def test_tray_refuses_foreign_pot_and_from_report(tmp_path):
    """A clip tray for a pot built with another pattern is refused (fingerprint); with
    --from-report the tray takes the pot's parameters and passes the mesh fit check."""
    pot = tmp_path / "pot"
    run = lambda *a: subprocess.run([sys.executable, str(ROOT / "orchid_pot.py"), *a],
                                    capture_output=True, text=True)
    assert run("box", "--seed", "3", "--out", str(pot)).returncode == 0
    bad = run("box", "--parts", "tray", "--tray", "clip", "--out", str(pot))       # seed 7
    assert bad.returncode != 0 and "pattern fingerprint differs" in bad.stderr
    good = run("box", "--from-report", str(pot) + "_report.json", "--parts", "tray", "--tray", "clip")
    assert good.returncode == 0
    fc = json.loads((tmp_path / "pot_tray_clip_report.json").read_text())["fit_check"]
    assert fc["ok"] and fc["seated_overlap_mm3"] <= 0.05 and fc["lifted_clips_engaged"] == 4


def test_planter_box_fits_repo_pot(tmp_path):
    """Planter (no golden yet) around a copy of the repo's box pot: one watertight body, no
    islands; the pot sits without overlap and slides in freely (fit check)."""
    import shutil
    for ext in (".stl", "_report.json"):
        shutil.copy(ROOT / f"orchid_box{ext}", tmp_path / f"pot{ext}")
    subprocess.run([sys.executable, str(ROOT / "orchid_pot.py"), "box", "--parts", "planter",
                    "--out", str(tmp_path / "pot")], check=True, capture_output=True)
    r = json.loads((tmp_path / "pot_planter_report.json").read_text())
    assert r["stl_roundtrip_watertight"] and r["bodies"] == 1 and r["unsupported_islands"] == 0
    fc = r["fit_check"]
    assert fc["ok"] and max(fc["insertion_overlap_mm3_at_10/30/60mm"]) <= 0.05


def test_ignored_accessory_options_are_errors():
    """--clip-* without --tray clip (and similar) must fail loudly, not build a part without
    clips (seen in practice); a hook longer than the flat band top must be refused."""
    run = lambda *a: subprocess.run([sys.executable, str(ROOT / "orchid_pot.py"), *a],
                                    capture_output=True, text=True)
    r = run("cylinder", "--parts", "pot,tray", "--clip-width", "3", "--graph-only")
    assert r.returncode == 2 and "only apply to the clip tray" in r.stderr
    r = run("box", "--planter-gap", "6")
    assert r.returncode == 2 and "need planter in --parts" in r.stderr
    r = run("cylinder", "--wall", "4", "--parts", "tray", "--tray", "clip", "--clip-hook", "3")
    assert r.returncode != 0 and "exceeds the flat band top" in r.stderr


def test_clip_groove_pot_and_tray(tmp_path):
    """Snap groove (--clip-groove): the pot stays printable (no islands), the clip tray
    fits (seated: no contact; every clip catches) at <= 1.5 % arm strain, clips exactly even."""
    out = tmp_path / "g"
    subprocess.run([sys.executable, str(ROOT / "orchid_pot.py"), "box", "--clip-groove", "--parts", "pot,tray",
                    "--tray", "clip", "--out", str(out)], check=True, capture_output=True)
    pot = json.loads((tmp_path / "g_report.json").read_text())
    assert pot["params"]["band"] == 8.0 and pot["unsupported_islands"] == 0 and pot["stl_roundtrip_watertight"]
    t = json.loads((tmp_path / "g_tray_clip_report.json").read_text())
    assert t["fit_check"]["ok"] and t["fit_check"]["lifted_clips_engaged"] == 4
    assert t["tray"]["clip_mode"] == "groove" and t["tray"]["clip_strain_pct"] <= 1.5
    assert t["clip_spacing_min_pct_of_even"] == 100.0


def test_solid_wall_pot_with_saucer(tmp_path):
    """--wall-style solid: closed wall, solid floor with drain holes; prints without
    support; the saucer fits."""
    out = tmp_path / "s"
    subprocess.run([sys.executable, str(ROOT / "orchid_pot.py"), "box", "--wall-style", "solid",
                    "--parts", "pot,tray", "--out", str(out)], check=True, capture_output=True)
    pot = json.loads((tmp_path / "s_report.json").read_text())
    assert pot["params"]["wall"] == 2.4 and pot["params"]["drain_holes"] == 5
    assert pot["stl_roundtrip_watertight"] and pot["bodies"] == 1 and pot["unsupported_islands"] == 0
    assert json.loads((tmp_path / "s_tray_report.json").read_text())["fit_check"]["ok"]


def test_clip_hook_with_groove_is_an_error(tmp_path):
    """--clip-hook has no effect in groove mode (reach = groove depth - clearance): error,
    also when the groove comes from the pot's report via --from-report."""
    run = lambda *a: subprocess.run([sys.executable, str(ROOT / "orchid_pot.py"), *a],
                                    capture_output=True, text=True)
    r = run("box", "--clip-groove", "--parts", "tray", "--tray", "clip", "--clip-hook", "0.8")
    assert r.returncode == 2 and "only applies to lattice clips" in r.stderr
    rep = tmp_path / "p_report.json"
    rep.write_text(json.dumps({"params": {"shape": "box", "clip_groove": True, "out": str(tmp_path / "p")}}))
    r = run("box", "--from-report", str(rep), "--parts", "tray", "--tray", "clip", "--clip-hook", "0.8")
    assert r.returncode == 2 and "only applies to lattice clips" in r.stderr
