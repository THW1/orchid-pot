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
