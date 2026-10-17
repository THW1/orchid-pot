# SPDX-License-Identifier: MIT
# === CONTEXT ===
# Purpose:     Guard the box and cylinder outputs against unintended change, e.g. when the
#              experimental orb shape or shared helpers are added.
# Decisions:   Compares the STL byte hash and the full JSON report (minus `params`, which
#              holds the output path) for seed 7. The pipeline is deterministic, so any
#              difference is a real change. The box hash equals the v0.4.0 release STL.
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
    assert fc["ok"] and fc["seated_overlap_mm3"] <= 0.05 and fc["lifted_1.5mm_hooks_engaged"] == 4
