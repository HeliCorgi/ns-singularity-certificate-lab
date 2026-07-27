from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ns_certificate_lab.artifacts import load_candidate
from ns_certificate_lab.axis import check_axis_regularity
from ns_certificate_lab.convergence import assess_convergence, max_abs
from ns_certificate_lab.diagnostics import (
    load_diagnostics_json,
    write_diagnostics_json,
)
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.manufactured import manufactured_fields
from ns_certificate_lab.operators import divergence_physical, elliptic_defect


def test_detects_velocity_field_that_breaks_divergence_free_condition() -> None:
    grid = AxisymmetricGrid.uniform(nr=17, nz=32, r_max=1.0)
    r, _ = grid.mesh()
    broken_u_r = r.copy()
    broken_u_z = np.zeros(grid.shape)
    defect = divergence_physical(grid, broken_u_r, broken_u_z)
    assert max_abs(defect) > 1.9


def test_detects_sign_reversal_in_elliptic_equation() -> None:
    grid = AxisymmetricGrid.uniform(nr=33, nz=64, r_max=1.0)
    fields = manufactured_fields(grid)
    wrong_sign_omega1 = -fields.omega1
    defect = elliptic_defect(grid, fields.psi1, wrong_sign_omega1)
    assert max_abs(defect) > 10.0


def test_detects_axis_condition_violation() -> None:
    grid = AxisymmetricGrid.uniform(nr=17, nz=32, r_max=1.0)
    fields = manufactured_fields(grid)
    r, _ = grid.mesh()
    broken_psi1 = fields.psi1 + r
    report = check_axis_regularity(
        grid,
        even_fields={"broken_psi1": broken_psi1},
        tolerance=1e-6,
    )
    assert not report.passed
    assert report.metrics["broken_psi1.max_abs_dr_at_axis"] > 0.9


def test_rejects_invalid_candidate_file(scratch_dir: Path) -> None:
    path = scratch_dir / "candidate.npz"
    path.write_bytes(b"PK but deliberately not a valid candidate")
    with pytest.raises(ValueError):
        load_candidate(path)


def test_detects_intentionally_modified_diagnostic_data(
    scratch_dir: Path,
) -> None:
    path = scratch_dir / "diagnostics.json"
    write_diagnostics_json(
        path,
        records=[{"resolution": 33, "residual": 0.125}],
    )
    data = bytearray(path.read_bytes())
    index = data.index(b"0.125")
    data[index : index + 5] = b"0.001"
    path.write_bytes(data)
    with pytest.raises(ValueError, match="checksum"):
        load_diagnostics_json(path)


def test_detects_nonconvergent_resolution_sequence() -> None:
    report = assess_convergence(
        errors=[0.10, 0.11, 0.105],
        spacings=[0.1, 0.05, 0.025],
        minimum_order=1.5,
    )
    assert not report.passed
    assert not report.strictly_decreasing
