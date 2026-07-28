r"""Verification experiment for equation-audit entry **E-33**.

What E-33 claims
----------------
For the finite-cylinder problem of E-25,

.. math::

   -\mathcal L_5\psi_1
   =-\Big(\partial_{rr}+\frac3r\partial_r+\partial_{zz}\Big)\psi_1=\omega_1,
   \qquad 0\le r\le R,\quad z\ \text{periodic with period }L_z,
   \quad \psi_1(R,z)=0,

with ``omega_1`` compactly supported in ``r <= a``:

(a) per axial mode ``k`` the substitution ``psi = phi / r`` turns
    ``L_{5,k} = d_rr + (3/r) d_r - k^2`` into the modified Bessel equation of
    order one, so the homogeneous solutions are ``I_1(kr)/r`` (regular on the
    axis) and ``K_1(kr)/r`` (decaying);
(b) the wall-induced correction is
    ``psi^(R) - psi^inf = -A [K_1(kR)/I_1(kR)] I_1(kr)/r``, so at a fixed core
    radius its magnitude is exactly proportional to ``K_1(kR)/I_1(kR)`` and its
    radial shape is exactly proportional to ``I_1(kr)/r``;
(c) ``K_1/I_1 ~ pi e^{-2x} (1 + 3/(4x))`` for ``x >> 1`` and ``~ 2/x^2`` for
    ``x << 1``, with the crossover at ``kR ~ 1``;
(d) at ``k = 0`` exactly, the wall-dependent part of ``psi`` in the core is the
    **constant** ``-Q_inf / (2 R^2)`` with ``Q_inf = int_0^inf s^3 omega_1 ds``,
    so moving the wall changes ``psi`` over the whole core by exactly
    ``(Q_inf/2)(R'^-2 - R^-2)``, independent of ``r``.

How this experiment tests them
------------------------------
A compactly supported radial profile ``(1 - (r/a)^2)^p`` (exactly zero for
``r >= a``) is multiplied by a single axial Fourier mode, and the elliptic
problem is solved on nested cylinders that share their core grid nodes
bitwise.  Nothing in the measurement chain uses a Bessel function: the
predictions come from :mod:`ns_certificate_lab.bessel_reference`, an
independent power-series / quadrature oracle that imports nothing from this
package, and the measurements come from the finite-cylinder solvers.

Two structural facts make the comparison exact rather than approximate.

* **The reference wall need not be a good whole-space proxy.**  Writing
  ``rho(x) = K_1(x)/I_1(x)``, E-33(b) gives
  ``psi^(R) - psi^(R_ref) = -A [rho(kR) - rho(kR_ref)] I_1(kr)/r`` on
  ``r <= min(R, R_ref)``.  Every ratio compared below therefore has an *exact*
  closed-form prediction, with no truncation of the reference wall's own
  residual response.  The quality of the proxy is reported separately as
  ``reference_residual_relative`` so the reader can see it, not relied upon.
* **The ``k = 0`` constancy of E-33(d) is a discrete identity, not an
  approximation.**  The difference of two solutions on nested cylinders
  satisfies the *homogeneous* discrete equation on every row of the smaller
  cylinder, and for ``k = 0`` the axis row ``8(delta_1 - delta_0)/dr^2 = 0``
  forces ``delta_1 = delta_0``, after which the recursion propagates a
  constant.  Constancy is therefore expected at roundoff for every resolution,
  while the *magnitude* carries the ``O(dr^2)`` discretization error.  Both are
  measured separately below, and conflating them would hide a real fault.

Multi-solver cross-check
------------------------
A representative subset is repeated with all three finite-cylinder solvers.
This module is where the three-way comparison lives, because
``tests/test_poisson_cross_validation.py`` enforces that no test module imports
both solver A and solver B (its scan covers ``tests/test_*.py`` only).  The
expected agreement levels are the documented ones, not roundoff:

* A vs B differ by their radial truncation error, ``O(dr^2)``;
* A vs C differ by their axial symbol gap, ``O(dz^2)``: A applies ``k^2``
  exactly per Fourier mode while C applies ``(4/dz^2) sin^2(k dz/2)``.  At
  ``k = 0`` that gap vanishes identically, so A and C agree to roundoff there
  and the ``k = 0`` comparison is **not** evidence of their independence --
  ``docs/realspace_poisson.md`` records that C's radial stencil is a
  transcription of the same E-26 formulas that A implements.

What this experiment does not establish
---------------------------------------
Nothing about the continuum, nothing about singularity formation, and nothing
that upgrades E-25's declared Dirichlet wall into a decay condition inherited
from the whole-space problem.  It measures how strongly *this* discretization's
solution depends on where that artificial wall is placed, and checks that the
dependence follows the law E-33 derives.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable, Sequence

import numpy as np

from ns_certificate_lab import bessel_reference
from ns_certificate_lab._integrity import strict_json_loads
from ns_certificate_lab.finite_cylinder_poisson import solve_finite_cylinder_poisson
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.poisson import solve_streamfunction_poisson
from ns_certificate_lab.provenance import collect_runtime_provenance
from ns_certificate_lab.realspace_poisson import solve_realspace_poisson

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

EXPERIMENT_ID = "wall_truncation_scaling_v1"
EXPECTED_SCHEMA_VERSION = 1
MINIMUM_RESOLUTIONS = 3
MINIMUM_WALL_RADII = 3
EPS = float(np.finfo(np.float64).eps)

#: Shape metrics saturate at roundoff once the discrete and continuum shapes
#: agree to working precision; an observed order is only meaningful above this.
SHAPE_ORDER_FLOOR = 1.0e-12

#: The E-33(a) residual is a fourth-order difference of oracle values, so it is
#: a difference of nearly equal terms divided by ``h^2``.  Below this level it
#: measures the cancellation roundoff of that difference and not the truncation
#: error, and its "observed order" is noise.  Measured on this machine: the
#: ``I`` branch at ``k = 0.196`` sits at ``2e-11`` and gets *worse* under
#: refinement, which is the signature of exactly that saturation.
MODAL_ODE_ROUNDOFF_FLOOR = 1.0e-9

#: The vocabulary :func:`classify_regime` may return.  ``crossover`` is not a
#: failure: E-33(c) only names two limits, and the band between them is where
#: neither asymptote is quantitative.
REGIMES: tuple[str, ...] = ("exponential", "crossover", "algebraic")

TOP_LEVEL_KEYS: frozenset[str] = frozenset(
    {
        "schema_version",
        "experiment",
        "description",
        "interpretation",
        "source_profile",
        "core_radius",
        "axial_cases",
        "wall_radii",
        "reference_wall_radius",
        "axial_points",
        "radial_resolutions",
        "cross_check",
        "nonlinear_cross_link",
        "classification",
        "roundoff_floor_factor",
        "well_resolved_floor_margin",
        "acceptance",
    }
)

SOURCE_PROFILE_KEYS: frozenset[str] = frozenset({"support_radius", "exponent"})
AXIAL_CASE_KEYS: frozenset[str] = frozenset({"label", "z_period", "mode"})
CLASSIFICATION_KEYS: frozenset[str] = frozenset(
    {"exponential_min_kr", "algebraic_max_kr"}
)
CROSS_CHECK_KEYS: frozenset[str] = frozenset(
    {
        "points_per_unit_radius",
        "axial_points",
        "wall_radii",
        "axial_case_labels",
        "realspace_tolerance",
        "realspace_max_iterations",
    }
)
NONLINEAR_CROSS_LINK_KEYS: frozenset[str] = frozenset(
    {"summary_path", "member_group", "z_period", "mode"}
)
ACCEPTANCE_KEYS: frozenset[str] = frozenset(
    {
        "max_oracle_ratio_relative_error_finest",
        "max_oracle_ratio_relative_error_any",
        "min_oracle_ratio_observed_order",
        "max_oracle_ratio_observed_order",
        "max_zero_mode_closed_form_relative_error_finest",
        "min_zero_mode_closed_form_observed_order",
        "max_zero_mode_closed_form_observed_order",
        "max_zero_mode_constancy_relative_spread",
        "max_shape_one_minus_cosine",
        "min_shape_observed_order",
        "max_shape_observed_order",
        "max_exponential_slope_relative_deviation",
        "max_algebraic_slope_absolute_deviation",
        "max_measured_versus_oracle_slope_relative_error",
        "max_modal_ode_relative_residual",
        "min_modal_ode_observed_order",
        "max_wronskian_relative_defect",
        "max_published_value_relative_error",
        "max_cross_solver_relative_difference",
        "min_cross_solver_ab_difference_over_dr_squared",
        "max_cross_solver_ab_difference_over_dr_squared",
        "max_cross_solver_ratio_relative_deviation",
    }
)

LIMITATIONS: tuple[str, ...] = (
    "binary64 arithmetic without outward rounding; every tolerance below is a "
    "measured floating-point band, not a proof",
    "the elliptic problem solved here is E-25's finite cylinder with periodic z "
    "and a homogeneous Dirichlet wall; E-33 describes how its solution moves "
    "when that artificial wall moves, and says nothing about the whole-space "
    "problem beyond the R^-2 tail rate it quantifies",
    "the radial profile (1-(r/a)^2)^p is compactly supported and C^(p-1), not "
    "C^infinity; E-33 needs only compact support, but the profile is not the "
    "E-32 envelope used by the nonlinear wall-dependence run",
    "the axial direction carries a single grid-resolved Fourier mode, so "
    "solvers A and B are exact in z and only the radial discretization is "
    "refined; the reported observed orders are radial orders",
    "wall responses that E-33 predicts below the binary64 resolution of the "
    "solution itself are recorded and classified, but are excluded from the "
    "oracle-agreement acceptance because their measured value is roundoff",
)

KNOWN_GAPS: tuple[str, ...] = (
    "solver C shares E-26's radial stencil with solver A, so their agreement at "
    "k=0 is a transcription check and not independent evidence for the radial "
    "scheme; solver B's different radial stencil is what covers that",
    "the nonlinear cross-link compares ratios of two different quantities -- "
    "elliptic wall responses and nonlinear amplification separations -- and is "
    "recorded for interpretation only; it is deliberately not an acceptance "
    "gate",
    "no continuum error bound is produced; the convergence orders are observed "
    "orders over three grids",
)


# ------------------------------------------------------------------ helpers


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, data: object) -> None:
    text = (
        json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)
        + "\n"
    )
    path.write_text(text, encoding="utf-8")


def _prepare_output(output_dir: Path) -> None:
    if output_dir.exists():
        if not output_dir.is_dir():
            raise FileExistsError(f"output path is not a directory: {output_dir}")
        if any(output_dir.iterdir()):
            raise FileExistsError(
                f"refusing to overwrite nonempty output directory: {output_dir}"
            )
    else:
        output_dir.mkdir(parents=True)


def _inside_repository(path: Path) -> bool:
    try:
        path.relative_to(REPOSITORY_ROOT)
    except ValueError:
        return False
    return True


def _number(value: object, *, name: str, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    converted = float(value)
    if not math.isfinite(converted):
        raise ValueError(f"{name} must be a finite number")
    if positive and converted <= 0.0:
        raise ValueError(f"{name} must be strictly positive")
    return converted


def _integer(value: object, *, name: str, minimum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return int(value)


def _exact_keys(value: object, expected: frozenset[str], *, name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be an object")
    if set(value) != expected:
        missing = sorted(expected - set(value))
        unknown = sorted(set(value) - expected)
        raise ValueError(
            f"{name} has missing keys {missing} and unknown keys {unknown}"
        )
    return value


def _order(coarse: float, fine: float, coarse_h: float, fine_h: float) -> float:
    """Return the observed order implied by two errors on two spacings."""

    if not (coarse > 0.0 and fine > 0.0):
        raise ValueError("observed orders need strictly positive error levels")
    return math.log(coarse / fine) / math.log(coarse_h / fine_h)


def _max_abs(values: np.ndarray) -> float:
    return float(np.max(np.abs(values)))


def _finite(values: Iterable[float]) -> bool:
    return all(math.isfinite(float(value)) for value in values)


# ------------------------------------------------------- configuration


def radial_point_count(points_per_unit_radius: int, wall_radius: float) -> int:
    """Return ``nr`` such that ``dr = 1/points_per_unit_radius`` exactly.

    Every wall radius shares the same radial spacing, so the core nodes
    ``r_i = i/points_per_unit_radius`` are the *same* floating-point numbers on
    every member of the sweep and the core comparisons need no interpolation.
    The product must be an integer for that to hold; a fractional product is
    rejected rather than rounded.
    """

    product = points_per_unit_radius * wall_radius
    nearest = round(product)
    if abs(product - nearest) > 1.0e-9 * max(1.0, abs(product)):
        raise ValueError(
            f"points_per_unit_radius={points_per_unit_radius} does not divide "
            f"wall radius {wall_radius} into whole cells"
        )
    return int(nearest) + 1


def validate_config(config: object) -> dict[str, Any]:
    """Reject every missing, unknown, mislabeled or inconsistent config key."""

    config = _exact_keys(config, TOP_LEVEL_KEYS, name="wall truncation config")
    if _integer(config["schema_version"], name="schema_version") != (
        EXPECTED_SCHEMA_VERSION
    ):
        raise ValueError("unsupported wall truncation config schema")
    if config["experiment"] != EXPERIMENT_ID:
        raise ValueError(
            f"experiment must equal the audited canonical value {EXPERIMENT_ID!r}"
        )
    for key in ("description", "interpretation"):
        if not isinstance(config[key], str) or not config[key].strip():
            raise ValueError(f"{key} must be a non-empty string")

    profile = _exact_keys(
        config["source_profile"], SOURCE_PROFILE_KEYS, name="source_profile"
    )
    support_radius = _number(
        profile["support_radius"], name="source_profile.support_radius", positive=True
    )
    exponent = _integer(profile["exponent"], name="source_profile.exponent", minimum=1)

    core_radius = _number(config["core_radius"], name="core_radius", positive=True)

    radii = config["wall_radii"]
    if not isinstance(radii, list) or len(radii) < MINIMUM_WALL_RADII:
        raise ValueError(
            f"wall_radii must list at least {MINIMUM_WALL_RADII} radii"
        )
    wall_radii = [
        _number(item, name=f"wall_radii[{index}]", positive=True)
        for index, item in enumerate(radii)
    ]
    if any(
        wall_radii[index] >= wall_radii[index + 1]
        for index in range(len(wall_radii) - 1)
    ):
        raise ValueError("wall_radii must be strictly increasing")
    reference = _number(
        config["reference_wall_radius"], name="reference_wall_radius", positive=True
    )
    if reference != wall_radii[-1]:
        raise ValueError(
            "reference_wall_radius must be the largest entry of wall_radii"
        )
    if support_radius >= wall_radii[0]:
        # E-33's premise is that the source vanishes before the wall.  Without
        # this the modal decomposition into a fixed particular solution plus an
        # R-dependent regular homogeneous solution is simply false, and every
        # oracle comparison below would be measuring something else.
        raise ValueError(
            "source_profile.support_radius must be strictly smaller than the "
            "smallest wall radius, otherwise E-33's compact-support premise "
            "fails"
        )
    if core_radius >= support_radius:
        raise ValueError("core_radius must lie strictly inside the source support")

    cases = config["axial_cases"]
    if not isinstance(cases, list) or not cases:
        raise ValueError("axial_cases must be a non-empty list")
    labels: list[str] = []
    zero_mode_count = 0
    for index, item in enumerate(cases):
        case = _exact_keys(item, AXIAL_CASE_KEYS, name=f"axial_cases[{index}]")
        label = case["label"]
        if not isinstance(label, str) or not label.strip():
            raise ValueError(f"axial_cases[{index}].label must be a non-empty string")
        labels.append(label)
        _number(case["z_period"], name=f"axial_cases[{index}].z_period", positive=True)
        mode = _integer(case["mode"], name=f"axial_cases[{index}].mode", minimum=0)
        if mode == 0:
            zero_mode_count += 1
    if len(set(labels)) != len(labels):
        raise ValueError("axial_cases labels must be unique")
    if zero_mode_count != 1:
        raise ValueError(
            "exactly one axial case must have mode 0; it is the E-33(d) test"
        )

    axial_points = _integer(config["axial_points"], name="axial_points", minimum=5)

    resolutions = config["radial_resolutions"]
    if not isinstance(resolutions, list) or len(resolutions) < MINIMUM_RESOLUTIONS:
        raise ValueError(
            "radial_resolutions must contain at least "
            f"{MINIMUM_RESOLUTIONS} entries so every claim has an observed order"
        )
    points_per_unit = [
        _integer(item, name=f"radial_resolutions[{index}]", minimum=8)
        for index, item in enumerate(resolutions)
    ]
    if any(
        points_per_unit[index] >= points_per_unit[index + 1]
        for index in range(len(points_per_unit) - 1)
    ):
        raise ValueError("radial_resolutions must be strictly increasing")
    for ppur in points_per_unit:
        for radius in wall_radii:
            radial_point_count(ppur, radius)
        # The premise check needs at least one interior node on which the
        # source is exactly zero; otherwise "compactly supported inside the
        # wall" is not visible on the grid at all.
        interior_zero_rows = sum(
            1
            for index in range(radial_point_count(ppur, wall_radii[0]) - 1)
            if index / ppur >= support_radius
        )
        if interior_zero_rows < 1:
            raise ValueError(
                f"radial_resolutions entry {ppur} leaves no interior grid node "
                "between the source support and the smallest wall"
            )

    classification = _exact_keys(
        config["classification"], CLASSIFICATION_KEYS, name="classification"
    )
    exponential_min = _number(
        classification["exponential_min_kr"],
        name="classification.exponential_min_kr",
        positive=True,
    )
    algebraic_max = _number(
        classification["algebraic_max_kr"],
        name="classification.algebraic_max_kr",
        positive=True,
    )
    if not algebraic_max < 1.0 < exponential_min:
        raise ValueError(
            "the E-33(c) crossover sits at kR ~ 1, so algebraic_max_kr < 1 < "
            "exponential_min_kr is required for the two regimes to be disjoint "
            "and to bracket the crossover"
        )

    cross = _exact_keys(config["cross_check"], CROSS_CHECK_KEYS, name="cross_check")
    cross_ppur = _integer(
        cross["points_per_unit_radius"],
        name="cross_check.points_per_unit_radius",
        minimum=8,
    )
    _integer(cross["axial_points"], name="cross_check.axial_points", minimum=5)
    cross_radii_raw = cross["wall_radii"]
    if not isinstance(cross_radii_raw, list) or len(cross_radii_raw) < 2:
        raise ValueError("cross_check.wall_radii must list at least two radii")
    cross_radii = [
        _number(item, name=f"cross_check.wall_radii[{index}]", positive=True)
        for index, item in enumerate(cross_radii_raw)
    ]
    if any(
        cross_radii[index] >= cross_radii[index + 1]
        for index in range(len(cross_radii) - 1)
    ):
        raise ValueError("cross_check.wall_radii must be strictly increasing")
    if not set(cross_radii) <= set(wall_radii):
        raise ValueError("cross_check.wall_radii must be a subset of wall_radii")
    if cross_radii[-1] != reference:
        raise ValueError(
            "cross_check.wall_radii must end at the reference wall radius"
        )
    for radius in cross_radii:
        radial_point_count(cross_ppur, radius)
    cross_labels = cross["axial_case_labels"]
    if not isinstance(cross_labels, list) or not cross_labels:
        raise ValueError("cross_check.axial_case_labels must be a non-empty list")
    if not set(cross_labels) <= set(labels):
        raise ValueError("cross_check.axial_case_labels must name declared cases")
    _number(
        cross["realspace_tolerance"],
        name="cross_check.realspace_tolerance",
        positive=True,
    )
    _integer(
        cross["realspace_max_iterations"],
        name="cross_check.realspace_max_iterations",
        minimum=1,
    )

    link = _exact_keys(
        config["nonlinear_cross_link"],
        NONLINEAR_CROSS_LINK_KEYS,
        name="nonlinear_cross_link",
    )
    if not isinstance(link["summary_path"], str) or not link["summary_path"].strip():
        raise ValueError("nonlinear_cross_link.summary_path must be a string")
    if not isinstance(link["member_group"], str) or not link["member_group"].strip():
        raise ValueError("nonlinear_cross_link.member_group must be a string")
    _number(
        link["z_period"], name="nonlinear_cross_link.z_period", positive=True
    )
    _integer(link["mode"], name="nonlinear_cross_link.mode", minimum=1)

    _number(
        config["roundoff_floor_factor"], name="roundoff_floor_factor", positive=True
    )
    if (
        _number(
            config["well_resolved_floor_margin"],
            name="well_resolved_floor_margin",
            positive=True,
        )
        < 1.0
    ):
        raise ValueError(
            "well_resolved_floor_margin must be at least 1; it is the extra "
            "safety factor above the roundoff floor demanded before an "
            "observed order or a shape correlation is claimed"
        )

    acceptance = _exact_keys(
        config["acceptance"], ACCEPTANCE_KEYS, name="acceptance"
    )
    for key in sorted(ACCEPTANCE_KEYS):
        _number(acceptance[key], name=f"acceptance.{key}")
    return config


# ------------------------------------------------------------- the source


def radial_profile(
    radii: np.ndarray, *, support_radius: float, exponent: int
) -> np.ndarray:
    """Return ``(1 - (r/a)^2)^p`` inside ``r < a`` and exactly zero outside.

    The branch is on the grid coordinate, so the values at and beyond ``a`` are
    the floating-point constant ``0.0``, not a small residue.  E-33's premise is
    an exact statement about the support, and the guard below tests it as one.
    """

    values = np.zeros_like(radii, dtype=np.float64)
    inside = radii < support_radius
    scaled = radii[inside] / support_radius
    values[inside] = (1.0 - scaled * scaled) ** exponent
    return values


def source_moment(*, support_radius: float, exponent: int) -> float:
    r"""Return ``Q_inf = int_0^inf s^3 omega_1(s) ds`` in closed form.

    With ``omega_1 = (1 - (s/a)^2)^p`` on ``s < a``, substituting ``s = a u``
    and then ``v = u^2`` gives

    .. math::

       Q_\infty=a^4\int_0^1u^3(1-u^2)^p\,du
               =\frac{a^4}{2}\int_0^1v(1-v)^p\,dv
               =\frac{a^4}{2}B(2,p+1)
               =\frac{a^4}{2(p+1)(p+2)} .

    This is the only place a continuum integral enters the ``k = 0`` oracle, and
    it is exact, not quadrature.
    """

    return support_radius**4 / (2.0 * (exponent + 1) * (exponent + 2))


def build_source(
    grid: AxisymmetricGrid,
    *,
    support_radius: float,
    exponent: int,
    z_period: float,
    mode: int,
) -> np.ndarray:
    """Return the nodal right-hand side ``omega_1`` for one axial case."""

    profile = radial_profile(
        grid.r, support_radius=support_radius, exponent=exponent
    )
    if mode == 0:
        return np.repeat(profile[:, None], grid.nz, axis=1)
    axial = np.cos(2.0 * np.pi * mode * grid.z / z_period)
    return profile[:, None] * axial[None, :]


def assert_compact_support(
    grid: AxisymmetricGrid,
    omega: np.ndarray,
    *,
    support_radius: float,
    smallest_wall_radius: float,
) -> int:
    """Refuse to measure E-33 on a source that violates its premise.

    Returns the number of interior grid rows on which the source is exactly
    zero.  E-33(b) rests on splitting the solution into a wall-independent
    particular part and a regular homogeneous part whose coefficient is fixed
    by the trace at ``r = R``; that split requires the source to vanish before
    the wall.  A source that reaches the wall would still produce numbers here,
    and they would silently mean nothing, so this is an error and not a
    warning.
    """

    if support_radius >= smallest_wall_radius:
        raise ValueError(
            "E-33 requires the source to be compactly supported strictly "
            f"inside the smallest wall: support radius {support_radius!r} is "
            f"not smaller than wall radius {smallest_wall_radius!r}"
        )
    outside = grid.r >= support_radius
    interior_outside = outside[:-1]
    zero_rows = int(np.count_nonzero(interior_outside))
    if zero_rows < 1:
        raise ValueError(
            "no interior grid node lies between the source support and the "
            "wall, so compact support is not observable on this grid"
        )
    residue = _max_abs(omega[outside, :])
    if residue != 0.0:
        raise ValueError(
            "the source is not exactly zero outside its declared support "
            f"radius {support_radius!r}: max |omega| there is {residue!r}"
        )
    return zero_rows


# ----------------------------------------------------------- the solvers


def solve_a(grid: AxisymmetricGrid, omega: np.ndarray) -> np.ndarray:
    """Solver A: ``r^3``-flux finite volume with a Fourier axial direction."""

    return solve_streamfunction_poisson(
        grid, omega, 0.0, estimate_condition=False
    ).psi1


def solve_b(grid: AxisymmetricGrid, omega: np.ndarray) -> np.ndarray:
    """Solver B: direct non-divergence radial stencil, Fourier axial direction."""

    return solve_finite_cylinder_poisson(
        omega, grid, outer_boundary=0.0, condition_mode_indices=()
    ).psi


def solve_c(
    grid: AxisymmetricGrid,
    omega: np.ndarray,
    *,
    tol: float,
    max_iterations: int,
) -> np.ndarray:
    """Solver C: real-space centered differences, Jacobi-preconditioned CG."""

    return solve_realspace_poisson(
        grid, omega, 0.0, tol=tol, max_iterations=max_iterations
    ).psi1


# ------------------------------------------------------------ the oracle


def wall_response_amplitude(wavenumber: float, radius: float, reference: float) -> float:
    """Return ``|rho(kR) - rho(kR_ref)|`` with ``rho = K_1/I_1`` (E-33b).

    This is the exact ``R``-dependence of ``max |psi^(R) - psi^(ref)|`` over any
    fixed set of core nodes, because the two solutions differ only by the
    regular homogeneous mode ``I_1(kr)/r`` whose coefficient is
    ``-A rho(kR)``.
    """

    return abs(
        bessel_reference.k1_over_i1(wavenumber * radius)
        - bessel_reference.k1_over_i1(wavenumber * reference)
    )


def predicted_shape(wavenumber: float, radii: np.ndarray) -> np.ndarray:
    """Return the E-33(b) radial shape ``I_1(kr)/r`` on the given nodes."""

    return np.array(
        [
            wavenumber * bessel_reference.i1_over_argument(wavenumber * float(radius))
            for radius in radii
        ],
        dtype=np.float64,
    )


def zero_mode_closed_form(
    *, source_moment_value: float, radius: float, other_radius: float
) -> float:
    """Return the E-33(e) difference ``psi^(R) - psi^(R')`` in the core.

    E-33(d)/(e): ``psi^(R)(r) = psi^part(r) - Q_inf/(2 R^2)`` outside the
    support, so the difference of two wall radii is the ``r``-independent
    constant ``(Q_inf/2)(R'^-2 - R^-2)``.  The sign is fixed by the derivation
    -- shrinking the wall lowers ``psi`` in the core for a positive source --
    and a sign error here is exactly what the fault-injection test checks.
    """

    return 0.5 * source_moment_value * (other_radius**-2 - radius**-2)


def modal_ode_relative_residual(
    wavenumber: float, radius: float, *, step: float, branch: str
) -> float:
    r"""Return the relative residual of E-33(a) for one homogeneous branch.

    E-33(a) states that ``psi = I_1(kr)/r`` and ``psi = K_1(kr)/r`` both solve

    .. math::

       \mathcal L_{5,k}\psi=\psi''+\frac3r\psi'-k^2\psi=0 .

    The derivatives are taken with fourth-order centered differences of the
    oracle values, so the check uses **only** the Bessel module and never a
    Bessel derivative identity; the residual must fall at fourth order until
    roundoff, which is what ``tests/test_wall_truncation_scaling.py`` pins.
    """

    if branch == "I":
        def value(argument: float) -> float:
            return bessel_reference.bessel_i1(wavenumber * argument) / argument
    elif branch == "K":
        def value(argument: float) -> float:
            return bessel_reference.bessel_k1(wavenumber * argument) / argument
    else:
        raise ValueError("branch must be 'I' or 'K'")

    samples = [value(radius + offset * step) for offset in (-2, -1, 0, 1, 2)]
    second = (
        -samples[0] / 12.0
        + 4.0 * samples[1] / 3.0
        - 2.5 * samples[2]
        + 4.0 * samples[3] / 3.0
        - samples[4] / 12.0
    ) / (step * step)
    first = (
        samples[0] / 12.0
        - 2.0 * samples[1] / 3.0
        + 2.0 * samples[3] / 3.0
        - samples[4] / 12.0
    ) / step
    skew = 3.0 * first / radius
    zeroth = wavenumber * wavenumber * samples[2]
    residual = second + skew - zeroth
    scale = max(abs(second), abs(skew), abs(zeroth))
    return abs(residual) / scale


def classify_regime(
    *, smaller_kr: float, larger_kr: float, exponential_min_kr: float, algebraic_max_kr: float
) -> str:
    """Classify one radius pair against the E-33(c) crossover at ``kR ~ 1``."""

    if smaller_kr >= exponential_min_kr:
        return "exponential"
    if larger_kr <= algebraic_max_kr:
        return "algebraic"
    return "crossover"


# ------------------------------------------------------- the measurement


def _core_block(field: np.ndarray, core_row_count: int) -> np.ndarray:
    """Return the leading ``core_row_count`` radial rows of a field.

    The nested cylinders share their radial spacing, so row ``i`` is the same
    coordinate ``i * dr`` on every member and this slice is a nodewise
    restriction with no interpolation.  ``measure_axial_case`` verifies that
    bitwise before using it.
    """

    return field[:core_row_count, :]


def measure_axial_case(
    *,
    case: dict[str, Any],
    points_per_unit_radius: int,
    wall_radii: Sequence[float],
    reference_radius: float,
    axial_points: int,
    support_radius: float,
    exponent: int,
    core_radius: float,
    roundoff_floor_factor: float,
    well_resolved_floor_margin: float,
    classification: dict[str, Any],
    source_moment_value: float,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Solve one axial case on every wall radius and measure the response."""

    label = str(case["label"])
    z_period = float(case["z_period"])
    mode = int(case["mode"])
    wavenumber = 0.0 if mode == 0 else 2.0 * np.pi * mode / z_period

    solutions: dict[float, np.ndarray] = {}
    grids: dict[float, AxisymmetricGrid] = {}
    zero_rows: int | None = None
    for radius in wall_radii:
        grid = AxisymmetricGrid.uniform(
            nr=radial_point_count(points_per_unit_radius, radius),
            nz=axial_points,
            r_max=radius,
            z_min=0.0,
            z_max=z_period,
            periodic_z=True,
        )
        omega = build_source(
            grid,
            support_radius=support_radius,
            exponent=exponent,
            z_period=z_period,
            mode=mode,
        )
        observed_zero_rows = assert_compact_support(
            grid,
            omega,
            support_radius=support_radius,
            smallest_wall_radius=float(min(wall_radii)),
        )
        if zero_rows is None:
            # Record the count for the SMALLEST wall, which is the binding one:
            # every larger cylinder trivially has more source-free rows.
            zero_rows = observed_zero_rows
        grids[radius] = grid
        solutions[radius] = solve_a(grid, omega)

    smallest_grid = grids[wall_radii[0]]
    core_rows = int(np.count_nonzero(smallest_grid.r <= core_radius + 1.0e-12))
    if core_rows < 4:
        raise ValueError("the core must contain at least four radial nodes")
    core_r = smallest_grid.r[:core_rows]

    # The shared core nodes must be identical floating-point numbers, otherwise
    # the nodewise differences below would silently include an interpolation
    # error that has nothing to do with the wall.
    for radius in wall_radii:
        if not np.array_equal(grids[radius].r[:core_rows], core_r):
            raise ValueError(
                "core radial nodes are not bitwise identical across wall radii"
            )

    reference_field = _core_block(solutions[reference_radius], core_rows)
    reference_scale = _max_abs(reference_field)
    if reference_scale == 0.0:
        raise ValueError("the reference solution vanishes on the core")
    floor = roundoff_floor_factor * EPS * reference_scale

    shape = predicted_shape(wavenumber, core_r) if mode else None
    shape_norm = float(np.linalg.norm(shape)) if shape is not None else 0.0

    non_reference = [radius for radius in wall_radii if radius != reference_radius]

    radius_rows: list[dict[str, Any]] = []
    responses: dict[float, float] = {}
    stored: dict[str, np.ndarray] = {"core_r": core_r}
    calibration: float | None = None
    for radius in non_reference:
        difference = _core_block(solutions[radius], core_rows) - reference_field
        response = _max_abs(difference)
        responses[radius] = response
        flat_index = int(np.argmax(np.abs(difference)))
        argmax_row, argmax_column = np.unravel_index(flat_index, difference.shape)
        # The constancy signature of E-33(d) is a statement about the RADIAL
        # profile at fixed z.  Averaging over a full axial period would send the
        # mean of a cosine mode to zero and make the ratio meaningless, so the
        # spread and mean are taken along the radial line through the argmax.
        column = difference[:, argmax_column]
        mean = float(np.mean(column))
        spread = float(np.max(column) - np.min(column))
        constancy = spread / abs(mean) if abs(mean) > floor else None
        amplitude = (
            wall_response_amplitude(wavenumber, radius, reference_radius)
            if mode
            else abs(
                zero_mode_closed_form(
                    source_moment_value=source_moment_value,
                    radius=radius,
                    other_radius=reference_radius,
                )
            )
        )
        if calibration is None:
            calibration = response / amplitude if amplitude > 0.0 else 0.0
        predicted = float(calibration) * amplitude
        margin = predicted / floor if floor > 0.0 else math.inf
        comparable = predicted >= floor
        if margin >= 4.0:
            consistent = response > floor
        elif margin <= 0.25:
            consistent = response < 16.0 * floor
        else:
            consistent = True

        row: dict[str, Any] = {
            "axial_case": label,
            "points_per_unit_radius": points_per_unit_radius,
            "wall_radius": radius,
            "kr": wavenumber * radius,
            "nr": grids[radius].nr,
            "response_max_abs": response,
            "response_relative_to_reference_core_max": response / reference_scale,
            "argmax_r": float(core_r[argmax_row]),
            "argmax_z": float(smallest_grid.z[argmax_column]),
            "response_mean": mean,
            "response_spread": spread,
            "response_spread_over_abs_mean": constancy,
            "oracle_amplitude": amplitude,
            "oracle_predicted_response": predicted,
            "floating_point_floor": floor,
            "floor_margin": margin,
            "oracle_comparable": bool(comparable),
            # ``oracle_comparable`` only asks whether the response is above the
            # noise at all, which is enough to compare one number.  Claiming an
            # observed CONVERGENCE ORDER, or correlating a whole radial
            # profile whose smallest component is far below its maximum,
            # requires much more headroom: the discretization error being
            # measured must itself dominate the roundoff.  That is what this
            # second, stricter flag records.
            "well_resolved": bool(margin >= well_resolved_floor_margin),
            "floor_classification_consistent": bool(consistent),
        }
        if mode:
            column_norm = float(np.linalg.norm(column))
            if column_norm > 0.0 and shape is not None:
                cosine = float(np.dot(column, shape) / (column_norm * shape_norm))
            else:
                cosine = 0.0
            row["shape_cosine_similarity"] = cosine
            row["shape_one_minus_abs_cosine"] = 1.0 - abs(cosine)
            # E-33(b) fixes the sign as well as the shape: with a positive
            # source, shrinking the wall subtracts a positive multiple of
            # I_1(kr)/r from the core, so the correlation must be negative on
            # the axial slice where the source is positive.
            row["shape_sign_matches_e33b"] = bool(
                cosine < 0.0 if difference[0, argmax_column] < 0.0 else cosine > 0.0
            )
            row["reference_residual_relative"] = (
                bessel_reference.k1_over_i1(wavenumber * reference_radius)
                / bessel_reference.k1_over_i1(wavenumber * non_reference[0])
            )
        else:
            closed_form = zero_mode_closed_form(
                source_moment_value=source_moment_value,
                radius=radius,
                other_radius=reference_radius,
            )
            row["closed_form_difference"] = closed_form
            row["measured_difference_mean"] = mean
            row["closed_form_relative_error"] = abs(mean - closed_form) / abs(
                closed_form
            )
            if constancy is None:
                raise ArithmeticError(
                    "the k=0 wall response vanished, so E-33(d) constancy "
                    "cannot be measured"
                )
            row["constancy_relative_spread"] = constancy
        radius_rows.append(row)
        stored[f"response_R{radius:g}"] = difference

    # Adjacent-radius pairs, measured against the reference wall.
    pair_rows: list[dict[str, Any]] = []
    for index in range(len(non_reference) - 1):
        smaller = non_reference[index]
        larger = non_reference[index + 1]
        measured_small = responses[smaller]
        measured_large = responses[larger]
        oracle_small = radius_rows[index]["oracle_amplitude"]
        oracle_large = radius_rows[index + 1]["oracle_amplitude"]
        comparable = bool(
            radius_rows[index]["oracle_comparable"]
            and radius_rows[index + 1]["oracle_comparable"]
        )
        well_resolved = bool(
            radius_rows[index]["well_resolved"]
            and radius_rows[index + 1]["well_resolved"]
        )
        measured_ratio = (
            measured_small / measured_large if measured_large > 0.0 else None
        )
        oracle_ratio = oracle_small / oracle_large if oracle_large > 0.0 else None
        pair: dict[str, Any] = {
            "axial_case": label,
            "points_per_unit_radius": points_per_unit_radius,
            "smaller_radius": smaller,
            "larger_radius": larger,
            "smaller_kr": wavenumber * smaller,
            "larger_kr": wavenumber * larger,
            "measured_ratio": measured_ratio,
            "oracle_ratio": oracle_ratio,
            "ratio_relative_error": (
                abs(measured_ratio - oracle_ratio) / oracle_ratio
                if measured_ratio is not None and oracle_ratio
                else None
            ),
            "oracle_comparable": comparable,
            "well_resolved": well_resolved,
            "regime": (
                classify_regime(
                    smaller_kr=wavenumber * smaller,
                    larger_kr=wavenumber * larger,
                    exponential_min_kr=float(classification["exponential_min_kr"]),
                    algebraic_max_kr=float(classification["algebraic_max_kr"]),
                )
                if mode
                else "zero_mode"
            ),
            "slopes_measured": False,
        }
        if mode and comparable and measured_small > 0.0 and measured_large > 0.0:
            delta_radius = larger - smaller
            log_ratio = math.log(measured_large / measured_small)
            oracle_log_ratio = math.log(oracle_large / oracle_small)
            pure_small = bessel_reference.k1_over_i1(wavenumber * smaller)
            pure_large = bessel_reference.k1_over_i1(wavenumber * larger)
            pure_log_ratio = math.log(pure_large / pure_small)
            slope_radius = log_ratio / delta_radius
            slope_log_radius = log_ratio / math.log(larger / smaller)
            oracle_slope_radius = oracle_log_ratio / delta_radius
            oracle_slope_log_radius = oracle_log_ratio / math.log(larger / smaller)
            pair.update(
                {
                    "slopes_measured": True,
                    "two_k": 2.0 * wavenumber,
                    "slope_dlog_dradius": slope_radius,
                    "slope_dlog_dlogradius": slope_log_radius,
                    "oracle_slope_dlog_dradius": oracle_slope_radius,
                    "oracle_slope_dlog_dlogradius": oracle_slope_log_radius,
                    "pure_oracle_slope_dlog_dradius": pure_log_ratio / delta_radius,
                    "pure_oracle_slope_dlog_dlogradius": (
                        pure_log_ratio / math.log(larger / smaller)
                    ),
                    "exponential_slope_relative_deviation": abs(
                        slope_radius / (-2.0 * wavenumber) - 1.0
                    ),
                    "algebraic_slope_absolute_deviation": abs(
                        slope_log_radius + 2.0
                    ),
                    # A degenerate oracle -- one that predicts no decay at all
                    # between the two radii -- has no slope to compare against,
                    # so this is recorded as unavailable rather than as a
                    # division by zero.  The ratio agreement above still sees
                    # such a corruption, and the fault-injection test relies on
                    # exactly that.
                    "measured_versus_oracle_slope_relative_error": (
                        abs(slope_radius / oracle_slope_radius - 1.0)
                        if oracle_slope_radius != 0.0
                        else None
                    ),
                    "leading_asymptote_ratio": math.exp(
                        2.0 * wavenumber * delta_radius
                    ),
                    "two_term_asymptote_ratio": (
                        bessel_reference.k1_over_i1_large_argument_asymptote(
                            wavenumber * smaller, terms=2
                        )
                        / bessel_reference.k1_over_i1_large_argument_asymptote(
                            wavenumber * larger, terms=2
                        )
                    ),
                    "small_argument_asymptote_ratio": (
                        bessel_reference.k1_over_i1_small_argument_asymptote(
                            wavenumber * smaller
                        )
                        / bessel_reference.k1_over_i1_small_argument_asymptote(
                            wavenumber * larger
                        )
                    ),
                }
            )
        pair_rows.append(pair)

    # Reference-free adjacent differences: the quantity whose ratios the
    # nonlinear amplification separations are directly comparable to.
    direct_rows: list[dict[str, Any]] = []
    for index in range(len(wall_radii) - 1):
        smaller = wall_radii[index]
        larger = wall_radii[index + 1]
        difference = _core_block(solutions[smaller], core_rows) - _core_block(
            solutions[larger], core_rows
        )
        oracle = (
            wall_response_amplitude(wavenumber, smaller, larger)
            if mode
            else abs(
                zero_mode_closed_form(
                    source_moment_value=source_moment_value,
                    radius=smaller,
                    other_radius=larger,
                )
            )
        )
        direct_rows.append(
            {
                "axial_case": label,
                "points_per_unit_radius": points_per_unit_radius,
                "smaller_radius": smaller,
                "larger_radius": larger,
                "difference_max_abs": _max_abs(difference),
                "oracle_amplitude": oracle,
            }
        )

    summary = {
        "axial_case": label,
        "z_period": z_period,
        "mode": mode,
        "wavenumber": wavenumber,
        "points_per_unit_radius": points_per_unit_radius,
        "dr": 1.0 / points_per_unit_radius,
        "axial_points": axial_points,
        "core_row_count": core_rows,
        "core_radius_of_last_node": float(core_r[-1]),
        "interior_zero_source_rows": zero_rows,
        "reference_wall_radius": reference_radius,
        "reference_core_max_abs_psi": reference_scale,
        "floating_point_floor": floor,
        "radii": radius_rows,
        "pairs": pair_rows,
        "direct_pairs": direct_rows,
    }
    if mode:
        summary["reference_residual_relative"] = (
            bessel_reference.k1_over_i1(wavenumber * reference_radius)
            / bessel_reference.k1_over_i1(wavenumber * non_reference[0])
        )
    return summary, stored


# ------------------------------------------------------ convergence study


def observed_orders(
    values: Sequence[float], spacings: Sequence[float]
) -> list[float | None]:
    """Return observed orders between consecutive refinement levels."""

    orders: list[float | None] = []
    for index in range(len(values) - 1):
        coarse = float(values[index])
        fine = float(values[index + 1])
        if coarse > 0.0 and fine > 0.0:
            orders.append(
                _order(coarse, fine, float(spacings[index]), float(spacings[index + 1]))
            )
        else:
            orders.append(None)
    return orders


def convergence_report(
    case_rows: Sequence[dict[str, Any]]
) -> dict[str, Any]:
    """Assemble observed orders for every claimed quantity of one axial case."""

    spacings = [float(row["dr"]) for row in case_rows]
    mode = int(case_rows[0]["mode"])
    report: dict[str, Any] = {
        "axial_case": case_rows[0]["axial_case"],
        "mode": mode,
        "radial_spacings": spacings,
        "ratio_error_orders": [],
        "shape_orders": [],
        "zero_mode_closed_form_orders": [],
    }

    pair_count = len(case_rows[0]["pairs"])
    for index in range(pair_count):
        pairs = [row["pairs"][index] for row in case_rows]
        if not all(pair["well_resolved"] for pair in pairs):
            continue
        errors = [float(pair["ratio_relative_error"]) for pair in pairs]
        report["ratio_error_orders"].append(
            {
                "smaller_radius": pairs[0]["smaller_radius"],
                "larger_radius": pairs[0]["larger_radius"],
                "relative_errors": errors,
                "orders": observed_orders(errors, spacings),
            }
        )

    if mode:
        radius = case_rows[0]["radii"][0]["wall_radius"]
        shape_values = [
            float(row["radii"][0]["shape_one_minus_abs_cosine"]) for row in case_rows
        ]
        entry: dict[str, Any] = {
            "wall_radius": radius,
            "one_minus_abs_cosine": shape_values,
            "orders": None,
            "order_available": False,
            "order_floor": SHAPE_ORDER_FLOOR,
        }
        if min(shape_values) > 0.0 and shape_values[0] > SHAPE_ORDER_FLOOR:
            entry["orders"] = observed_orders(shape_values, spacings)
            entry["order_available"] = True
        report["shape_orders"].append(entry)
    else:
        for position in range(len(case_rows[0]["radii"])):
            errors = [
                float(row["radii"][position]["closed_form_relative_error"])
                for row in case_rows
            ]
            report["zero_mode_closed_form_orders"].append(
                {
                    "wall_radius": case_rows[0]["radii"][position]["wall_radius"],
                    "relative_errors": errors,
                    "orders": observed_orders(errors, spacings),
                }
            )
    return report


# ---------------------------------------------------- multi-solver check


def cross_solver_report(
    *,
    config: dict[str, Any],
    source_moment_value: float,
) -> dict[str, Any]:
    """Repeat a representative subset with solvers A, B and C."""

    cross = config["cross_check"]
    points_per_unit_radius = int(cross["points_per_unit_radius"])
    axial_points = int(cross["axial_points"])
    radii = [float(value) for value in cross["wall_radii"]]
    reference_radius = radii[-1]
    tol = float(cross["realspace_tolerance"])
    max_iterations = int(cross["realspace_max_iterations"])
    support_radius = float(config["source_profile"]["support_radius"])
    exponent = int(config["source_profile"]["exponent"])
    core_radius = float(config["core_radius"])
    dr = 1.0 / points_per_unit_radius

    cases = {str(case["label"]): case for case in config["axial_cases"]}
    entries: list[dict[str, Any]] = []
    for label in cross["axial_case_labels"]:
        case = cases[str(label)]
        z_period = float(case["z_period"])
        mode = int(case["mode"])
        wavenumber = 0.0 if mode == 0 else 2.0 * np.pi * mode / z_period
        fields: dict[str, dict[float, np.ndarray]] = {"A": {}, "B": {}, "C": {}}
        grids: dict[float, AxisymmetricGrid] = {}
        differences: list[dict[str, Any]] = []
        for radius in radii:
            grid = AxisymmetricGrid.uniform(
                nr=radial_point_count(points_per_unit_radius, radius),
                nz=axial_points,
                r_max=radius,
                z_min=0.0,
                z_max=z_period,
                periodic_z=True,
            )
            omega = build_source(
                grid,
                support_radius=support_radius,
                exponent=exponent,
                z_period=z_period,
                mode=mode,
            )
            assert_compact_support(
                grid,
                omega,
                support_radius=support_radius,
                smallest_wall_radius=float(min(radii)),
            )
            grids[radius] = grid
            psi_a = solve_a(grid, omega)
            psi_b = solve_b(grid, omega)
            psi_c = solve_c(grid, omega, tol=tol, max_iterations=max_iterations)
            fields["A"][radius] = psi_a
            fields["B"][radius] = psi_b
            fields["C"][radius] = psi_c
            scale = _max_abs(psi_a)
            gap_ab = _max_abs(psi_a - psi_b)
            differences.append(
                {
                    "axial_case": label,
                    "wall_radius": radius,
                    "solution_max_abs_a": scale,
                    "max_abs_difference_ab": gap_ab,
                    "max_abs_difference_ac": _max_abs(psi_a - psi_c),
                    "max_abs_difference_bc": _max_abs(psi_b - psi_c),
                    "relative_difference_ab": gap_ab / scale,
                    "relative_difference_ac": _max_abs(psi_a - psi_c) / scale,
                    "relative_difference_bc": _max_abs(psi_b - psi_c) / scale,
                    "difference_ab_over_dr_squared": gap_ab / (dr * dr),
                }
            )

        core_rows = int(
            np.count_nonzero(grids[radii[0]].r <= core_radius + 1.0e-12)
        )
        per_solver: dict[str, Any] = {}
        for name in ("A", "B", "C"):
            reference_block = _core_block(
                fields[name][reference_radius], core_rows
            )
            responses = []
            for radius in radii[:-1]:
                block = _core_block(fields[name][radius], core_rows)
                responses.append(_max_abs(block - reference_block))
            entry: dict[str, Any] = {
                "wall_radii": radii[:-1],
                "responses": responses,
                "ratios": [
                    responses[index] / responses[index + 1]
                    for index in range(len(responses) - 1)
                ],
            }
            if mode == 0:
                block_small = _core_block(fields[name][radii[0]], core_rows)
                block_large = _core_block(fields[name][radii[1]], core_rows)
                difference = block_small - block_large
                mean = float(np.mean(difference))
                closed_form = zero_mode_closed_form(
                    source_moment_value=source_moment_value,
                    radius=radii[0],
                    other_radius=radii[1],
                )
                entry["zero_mode_pair"] = [radii[0], radii[1]]
                entry["zero_mode_measured_difference"] = mean
                entry["zero_mode_closed_form"] = closed_form
                entry["zero_mode_relative_error"] = abs(mean - closed_form) / abs(
                    closed_form
                )
                entry["zero_mode_constancy_relative_spread"] = float(
                    np.max(difference) - np.min(difference)
                ) / abs(mean)
            per_solver[name] = entry

        if mode:
            oracle_ratios = []
            amplitudes = [
                wall_response_amplitude(wavenumber, radius, reference_radius)
                for radius in radii[:-1]
            ]
            for index in range(len(amplitudes) - 1):
                oracle_ratios.append(amplitudes[index] / amplitudes[index + 1])
        else:
            oracle_ratios = []

        ratio_deviation = 0.0
        for index in range(len(per_solver["A"]["ratios"])):
            values = [per_solver[name]["ratios"][index] for name in ("A", "B", "C")]
            spread = max(values) - min(values)
            centre = sum(values) / len(values)
            ratio_deviation = max(ratio_deviation, spread / abs(centre))

        entries.append(
            {
                "axial_case": label,
                "mode": mode,
                "wavenumber": wavenumber,
                "points_per_unit_radius": points_per_unit_radius,
                "axial_points": axial_points,
                "dr": dr,
                "dz": z_period / axial_points,
                "per_radius_differences": differences,
                "per_solver": per_solver,
                "oracle_ratios": oracle_ratios,
                "solver_ratio_relative_spread": ratio_deviation,
            }
        )

    all_relative = [
        max(
            float(row["relative_difference_ab"]),
            float(row["relative_difference_ac"]),
            float(row["relative_difference_bc"]),
        )
        for entry in entries
        for row in entry["per_radius_differences"]
    ]
    ab_scaled = [
        float(row["difference_ab_over_dr_squared"])
        for entry in entries
        for row in entry["per_radius_differences"]
    ]
    return {
        "interpretation": (
            "solvers A and B differ by their radial truncation error O(dr^2); "
            "solvers A and C differ by their axial symbol gap O(dz^2), which "
            "vanishes identically at k=0, so the k=0 agreement of A and C is a "
            "transcription check of the shared E-26 radial stencil and not "
            "evidence of independence (docs/realspace_poisson.md)"
        ),
        "entries": entries,
        "max_relative_difference": max(all_relative) if all_relative else 0.0,
        "min_ab_difference_over_dr_squared": min(ab_scaled) if ab_scaled else 0.0,
        "max_ab_difference_over_dr_squared": max(ab_scaled) if ab_scaled else 0.0,
        "max_solver_ratio_relative_spread": max(
            float(entry["solver_ratio_relative_spread"]) for entry in entries
        )
        if entries
        else 0.0,
    }


# ------------------------------------------------------- the Bessel oracle
# ------------------------------------------------------- self-check block


PUBLISHED_VALUES: tuple[tuple[str, float, float], ...] = (
    # (name, argument, published value).  Sources: Abramowitz & Stegun Table 9.8
    # tabulates e^{-x} I_1(x) = 0.2079104154 and e^{x} K_1(x) = 1.6361534863 at
    # x = 1; multiplying by e^{+1} and e^{-1} gives the values below, which are
    # the same ones DLMF 10.25/10.31 produce.
    ("I1", 1.0, 0.5651591039924850),
    ("K1", 1.0, 0.6019072301972346),
)


def oracle_selfcheck(*, wavenumbers: Sequence[float]) -> dict[str, Any]:
    """Record the module-A validation evidence inside the run's own summary."""

    published = []
    for name, argument, expected in PUBLISHED_VALUES:
        observed = (
            bessel_reference.bessel_i1(argument)
            if name == "I1"
            else bessel_reference.bessel_k1(argument)
        )
        published.append(
            {
                "function": name,
                "argument": argument,
                "published": expected,
                "observed": observed,
                "relative_error": abs(observed - expected) / abs(expected),
            }
        )

    wronskian = [
        {
            "argument": argument,
            "relative_defect": bessel_reference.wronskian_relative_defect(argument),
        }
        for argument in (1.0e-6, 1.0e-3, 0.1, 1.0, 5.0, 20.0, 60.0)
    ]
    step_halving = [
        {
            "argument": argument,
            "relative_change": bessel_reference.k_quadrature_step_halving_defect(
                argument
            ),
        }
        for argument in (1.0e-3, 0.05, 1.0, 20.0, 60.0)
    ]

    large = []
    for argument in (3.0, 5.0, 10.0, 20.0, 40.0, 60.0):
        exact = bessel_reference.k1_over_i1(argument)
        leading = bessel_reference.k1_over_i1_large_argument_asymptote(
            argument, terms=1
        )
        two_term = bessel_reference.k1_over_i1_large_argument_asymptote(
            argument, terms=2
        )
        three_term = bessel_reference.k1_over_i1_large_argument_asymptote(
            argument, terms=3
        )
        large.append(
            {
                "argument": argument,
                "exact": exact,
                "leading_relative_error": abs(exact / leading - 1.0),
                "two_term_relative_error": abs(exact / two_term - 1.0),
                "three_term_relative_error": abs(exact / three_term - 1.0),
                "scaled_two_term_remainder": argument
                * argument
                * (exact / two_term - 1.0),
                "log_slope_dlog_dx": (
                    math.log(
                        bessel_reference.k1_over_i1(argument + 1.0e-3)
                        / bessel_reference.k1_over_i1(argument)
                    )
                    / 1.0e-3
                ),
            }
        )

    small = []
    for argument in (1.0e-4, 1.0e-3, 1.0e-2, 0.1, 0.2):
        exact = bessel_reference.k1_over_i1(argument)
        asymptote = bessel_reference.k1_over_i1_small_argument_asymptote(argument)
        stepped = argument * 1.02
        small.append(
            {
                "argument": argument,
                "exact": exact,
                "relative_error_versus_two_over_x_squared": abs(exact / asymptote - 1.0),
                "log_slope_dlog_dlogx": (
                    math.log(bessel_reference.k1_over_i1(stepped) / exact)
                    / math.log(1.02)
                ),
            }
        )

    ode = []
    for wavenumber in wavenumbers:
        if wavenumber <= 0.0:
            continue
        for branch in ("I", "K"):
            residuals = []
            for step in (0.02, 0.01):
                residuals.append(
                    max(
                        modal_ode_relative_residual(
                            wavenumber, radius, step=step, branch=branch
                        )
                        for radius in (0.3, 0.5, 0.9, 1.5, 2.0)
                    )
                )
            meaningful = min(residuals) >= MODAL_ODE_ROUNDOFF_FLOOR
            ode.append(
                {
                    "wavenumber": wavenumber,
                    "branch": branch,
                    "steps": [0.02, 0.01],
                    "relative_residuals": residuals,
                    "order_meaningful": bool(meaningful),
                    "roundoff_floor": MODAL_ODE_ROUNDOFF_FLOOR,
                    "observed_order": (
                        _order(residuals[0], residuals[1], 0.02, 0.01)
                        if min(residuals) > 0.0
                        else None
                    ),
                }
            )

    return {
        "statement": (
            "E-33(a) is checked by the finite-difference residual of "
            "L_{5,k} psi for psi = I_1(kr)/r and psi = K_1(kr)/r; E-33(c) is "
            "checked directly against both asymptotic branches; the oracle "
            "itself is checked against published values and the DLMF 10.28.2 "
            "cross-product identity"
        ),
        "published_values": published,
        "max_published_relative_error": max(
            float(item["relative_error"]) for item in published
        ),
        "wronskian": wronskian,
        "max_wronskian_relative_defect": max(
            float(item["relative_defect"]) for item in wronskian
        ),
        "quadrature_step_halving": step_halving,
        "max_step_halving_relative_change": max(
            float(item["relative_change"]) for item in step_halving
        ),
        "large_argument_branch": large,
        "small_argument_branch": small,
        "modal_ode_residuals": ode,
        "max_modal_ode_relative_residual": max(
            float(item["relative_residuals"][-1]) for item in ode
        )
        if ode
        else 0.0,
        "min_modal_ode_observed_order": min(
            (
                float(item["observed_order"])
                for item in ode
                if item["order_meaningful"] and item["observed_order"] is not None
            ),
            default=0.0,
        ),
        "modal_ode_order_interpretation": (
            "orders are claimed only where both residuals exceed "
            f"{MODAL_ODE_ROUNDOFF_FLOOR!r}; below that the fourth-order "
            "difference of the oracle values is dominated by cancellation "
            "roundoff, which is why the I branch at the smallest wavenumbers "
            "shows no convergence at all"
        ),
    }


# -------------------------------------------------- nonlinear cross-link


def nonlinear_cross_link(
    *,
    config: dict[str, Any],
    case_rows: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    """Compare the nonlinear run's wall sensitivity with the elliptic law.

    ``outputs/wall_dependence_v1`` measured, for nested walls, the separation
    between core amplifications ``|A_R - A_R'|``.  If the whole wall sensitivity
    of the nonlinear run were the linear elliptic effect of E-33, then the ratio
    of consecutive separations would equal the ratio of consecutive elliptic
    *pairwise* wall differences, whose oracle is
    ``|rho(kR_i) - rho(kR_{i+1})| / |rho(kR_{i+1}) - rho(kR_{i+2})|``.

    This is recorded, never gated.  Agreement supports the reading that the
    nonlinear run's wall sensitivity is the linear elliptic mechanism; it does
    not prove it, because the two quantities are different functionals of
    different fields.
    """

    link = config["nonlinear_cross_link"]
    path = (REPOSITORY_ROOT / str(link["summary_path"])).resolve()
    z_period = float(link["z_period"])
    mode = int(link["mode"])
    wavenumber = 2.0 * np.pi * mode / z_period
    result: dict[str, Any] = {
        "summary_path": str(link["summary_path"]),
        "member_group": link["member_group"],
        "z_period": z_period,
        "mode": mode,
        "wavenumber": wavenumber,
        "available": False,
        "gated": False,
    }
    if not path.is_file():
        result["interpretation"] = (
            "the nonlinear wall-dependence summary is not present in this "
            "checkout, so the cross-link was not evaluated"
        )
        return result

    payload = strict_json_loads(
        path.read_text(encoding="utf-8"), label="wall dependence summary"
    )
    try:
        separations = payload["analysis"]["member_groups"][str(link["member_group"])][
            "amplification_separations"
        ]
    except (KeyError, TypeError):
        result["interpretation"] = (
            "the nonlinear summary does not contain the requested member group's "
            "amplification separations"
        )
        return result

    rows: list[dict[str, Any]] = []
    for entry in separations:
        rows.append(
            {
                "smaller_wall_radius": float(entry["smaller_wall_radius"]),
                "larger_wall_radius": float(entry["larger_wall_radius"]),
                "separation": float(entry["separation"]),
            }
        )
    comparisons: list[dict[str, Any]] = []
    for index in range(len(rows) - 1):
        first, second = rows[index], rows[index + 1]
        if first["separation"] <= 0.0 or second["separation"] <= 0.0:
            continue
        nonlinear_ratio = first["separation"] / second["separation"]
        oracle_first = wall_response_amplitude(
            wavenumber, first["smaller_wall_radius"], first["larger_wall_radius"]
        )
        oracle_second = wall_response_amplitude(
            wavenumber, second["smaller_wall_radius"], second["larger_wall_radius"]
        )
        oracle_ratio = oracle_first / oracle_second
        comparison = {
            "nonlinear_pairs": [
                [first["smaller_wall_radius"], first["larger_wall_radius"]],
                [second["smaller_wall_radius"], second["larger_wall_radius"]],
            ],
            "nonlinear_separation_ratio": nonlinear_ratio,
            "elliptic_oracle_ratio": oracle_ratio,
            "relative_difference": abs(nonlinear_ratio - oracle_ratio) / oracle_ratio,
        }
        measured = _measured_direct_ratio(
            case_rows,
            wavenumber=wavenumber,
            first_pair=(first["smaller_wall_radius"], first["larger_wall_radius"]),
            second_pair=(second["smaller_wall_radius"], second["larger_wall_radius"]),
        )
        if measured is not None:
            comparison["elliptic_measured_ratio"] = measured
            comparison["relative_difference_versus_measured"] = (
                abs(nonlinear_ratio - measured) / measured
            )
        comparisons.append(comparison)

    result["available"] = True
    result["nonlinear_separations"] = rows
    result["comparisons"] = comparisons
    result["max_relative_difference"] = (
        max(float(item["relative_difference"]) for item in comparisons)
        if comparisons
        else None
    )
    result["interpretation"] = (
        "The nonlinear run's consecutive amplification-separation ratios are "
        "compared with the ratios E-33(b) predicts for the purely linear "
        "elliptic wall response at the same wavenumber and the same wall radii. "
        "Close agreement means the wall sensitivity measured in the nonlinear "
        "run is quantitatively what the linear elliptic truncation alone would "
        "produce, i.e. the nonlinearity does not amplify or suppress the wall "
        "effect over that time window. This is an interpretive cross-link "
        "between two different functionals and is deliberately not an "
        "acceptance gate."
    )
    return result


def _measured_direct_ratio(
    case_rows: Sequence[dict[str, Any]],
    *,
    wavenumber: float,
    first_pair: tuple[float, float],
    second_pair: tuple[float, float],
) -> float | None:
    """Return the finest-grid measured ratio of two adjacent wall differences."""

    candidates = [
        row
        for row in case_rows
        if row["mode"] > 0 and math.isclose(row["wavenumber"], wavenumber, rel_tol=1e-12)
    ]
    if not candidates:
        return None
    finest = max(candidates, key=lambda row: row["points_per_unit_radius"])
    lookup = {
        (float(item["smaller_radius"]), float(item["larger_radius"])): float(
            item["difference_max_abs"]
        )
        for item in finest["direct_pairs"]
    }
    if first_pair not in lookup or second_pair not in lookup:
        return None
    if lookup[second_pair] <= 0.0:
        return None
    return lookup[first_pair] / lookup[second_pair]


# ---------------------------------------------------------------- driver


def evaluate(config: dict[str, Any]) -> dict[str, Any]:
    """Run every measurement and assemble the analysis (no file output)."""

    support_radius = float(config["source_profile"]["support_radius"])
    exponent = int(config["source_profile"]["exponent"])
    core_radius = float(config["core_radius"])
    wall_radii = [float(value) for value in config["wall_radii"]]
    reference_radius = float(config["reference_wall_radius"])
    axial_points = int(config["axial_points"])
    resolutions = [int(value) for value in config["radial_resolutions"]]
    roundoff_floor_factor = float(config["roundoff_floor_factor"])
    well_resolved_floor_margin = float(config["well_resolved_floor_margin"])
    classification = config["classification"]
    moment = source_moment(support_radius=support_radius, exponent=exponent)

    case_rows: list[dict[str, Any]] = []
    stored_arrays: dict[str, np.ndarray] = {}
    for case in config["axial_cases"]:
        for points_per_unit_radius in resolutions:
            row, stored = measure_axial_case(
                case=case,
                points_per_unit_radius=points_per_unit_radius,
                wall_radii=wall_radii,
                reference_radius=reference_radius,
                axial_points=axial_points,
                support_radius=support_radius,
                exponent=exponent,
                core_radius=core_radius,
                roundoff_floor_factor=roundoff_floor_factor,
                well_resolved_floor_margin=well_resolved_floor_margin,
                classification=classification,
                source_moment_value=moment,
            )
            case_rows.append(row)
            if points_per_unit_radius == resolutions[-1]:
                prefix = f"{case['label']}__"
                for name, array in stored.items():
                    stored_arrays[prefix + name] = array

    convergence = []
    for case in config["axial_cases"]:
        label = str(case["label"])
        rows = [row for row in case_rows if row["axial_case"] == label]
        convergence.append(convergence_report(rows))

    wavenumbers = sorted(
        {float(row["wavenumber"]) for row in case_rows if row["mode"] > 0}
    )
    selfcheck = oracle_selfcheck(wavenumbers=wavenumbers)
    cross = cross_solver_report(config=config, source_moment_value=moment)
    link = nonlinear_cross_link(config=config, case_rows=case_rows)

    crossover = [
        {
            "axial_case": row["axial_case"],
            "wavenumber": row["wavenumber"],
            "points_per_unit_radius": row["points_per_unit_radius"],
            "pairs": [
                {
                    "smaller_radius": pair["smaller_radius"],
                    "larger_radius": pair["larger_radius"],
                    "smaller_kr": pair["smaller_kr"],
                    "larger_kr": pair["larger_kr"],
                    "regime": pair["regime"],
                }
                for pair in row["pairs"]
            ],
        }
        for row in case_rows
        if row["mode"] > 0 and row["points_per_unit_radius"] == resolutions[-1]
    ]

    return {
        "source_moment_q_infinity": moment,
        "radial_resolutions": resolutions,
        "axial_case_measurements": case_rows,
        "convergence": convergence,
        "crossover_classification": crossover,
        "bessel_oracle_selfcheck": selfcheck,
        "cross_solver": cross,
        "nonlinear_cross_link": link,
        "stored_arrays": stored_arrays,
    }


def acceptance_report(
    analysis: dict[str, Any], config: dict[str, Any]
) -> dict[str, Any]:
    """Evaluate the preregistered acceptance checks against the measurements."""

    tolerances = config["acceptance"]
    resolutions = [int(value) for value in config["radial_resolutions"]]
    finest = resolutions[-1]
    rows = analysis["axial_case_measurements"]

    ratio_errors_any: list[float] = []
    ratio_errors_finest: list[float] = []
    shape_values: list[float] = []
    exponential_deviations: list[float] = []
    algebraic_deviations: list[float] = []
    slope_errors: list[float] = []
    zero_mode_finest: list[float] = []
    zero_mode_constancy: list[float] = []
    floor_consistent = True
    all_values_finite = True
    shape_signs_correct = True
    minimum_zero_source_rows = min(
        int(row["interior_zero_source_rows"]) for row in rows
    )

    for row in rows:
        for entry in row["radii"]:
            all_values_finite &= _finite(
                value
                for value in entry.values()
                if isinstance(value, (int, float)) and not isinstance(value, bool)
            )
            floor_consistent &= bool(entry["floor_classification_consistent"])
            if row["mode"]:
                # A response sitting near the roundoff floor has an arbitrary
                # shape: the max-norm of the difference is still meaningful,
                # because it is attained at the outer core node where the
                # predicted profile peaks, but the inner nodes -- which the
                # cosine similarity weighs equally -- are pure noise there.
                # Only well-resolved radii can carry the shape claim.
                if entry["well_resolved"]:
                    shape_values.append(float(entry["shape_one_minus_abs_cosine"]))
                    shape_signs_correct &= bool(entry["shape_sign_matches_e33b"])
            else:
                zero_mode_constancy.append(float(entry["constancy_relative_spread"]))
                if row["points_per_unit_radius"] == finest:
                    zero_mode_finest.append(
                        float(entry["closed_form_relative_error"])
                    )
        for pair in row["pairs"]:
            all_values_finite &= _finite(
                value
                for value in pair.values()
                if isinstance(value, (int, float)) and not isinstance(value, bool)
            )
            if not pair["oracle_comparable"]:
                continue
            error = float(pair["ratio_relative_error"])
            ratio_errors_any.append(error)
            if row["points_per_unit_radius"] == finest:
                ratio_errors_finest.append(error)
            if row["mode"] and pair["slopes_measured"]:
                if pair["measured_versus_oracle_slope_relative_error"] is not None:
                    slope_errors.append(
                        float(pair["measured_versus_oracle_slope_relative_error"])
                    )
                if pair["regime"] == "exponential":
                    exponential_deviations.append(
                        float(pair["exponential_slope_relative_deviation"])
                    )
                elif pair["regime"] == "algebraic":
                    algebraic_deviations.append(
                        float(pair["algebraic_slope_absolute_deviation"])
                    )

    ratio_orders = [
        float(order)
        for report in analysis["convergence"]
        for entry in report["ratio_error_orders"]
        for order in entry["orders"]
        if order is not None
    ]
    zero_mode_orders = [
        float(order)
        for report in analysis["convergence"]
        for entry in report["zero_mode_closed_form_orders"]
        for order in entry["orders"]
        if order is not None
    ]
    shape_orders = [
        float(order)
        for report in analysis["convergence"]
        for entry in report["shape_orders"]
        if entry["order_available"]
        for order in entry["orders"]
        if order is not None
    ]

    selfcheck = analysis["bessel_oracle_selfcheck"]
    cross = analysis["cross_solver"]

    checks = {
        "compact_support_premise_observable_on_every_grid": bool(
            minimum_zero_source_rows >= 1
        ),
        "all_reported_metrics_finite": bool(all_values_finite),
        "floor_classification_consistent": bool(floor_consistent),
        "shape_sign_matches_e33b": bool(shape_signs_correct),
        "oracle_ratio_agreement_finest": bool(
            ratio_errors_finest
            and max(ratio_errors_finest)
            <= float(tolerances["max_oracle_ratio_relative_error_finest"])
        ),
        "oracle_ratio_agreement_all_resolutions": bool(
            ratio_errors_any
            and max(ratio_errors_any)
            <= float(tolerances["max_oracle_ratio_relative_error_any"])
        ),
        "oracle_ratio_orders_in_band": bool(
            ratio_orders
            and min(ratio_orders) >= float(tolerances["min_oracle_ratio_observed_order"])
            and max(ratio_orders) <= float(tolerances["max_oracle_ratio_observed_order"])
        ),
        "zero_mode_closed_form_agreement_finest": bool(
            zero_mode_finest
            and max(zero_mode_finest)
            <= float(tolerances["max_zero_mode_closed_form_relative_error_finest"])
        ),
        "zero_mode_closed_form_orders_in_band": bool(
            zero_mode_orders
            and min(zero_mode_orders)
            >= float(tolerances["min_zero_mode_closed_form_observed_order"])
            and max(zero_mode_orders)
            <= float(tolerances["max_zero_mode_closed_form_observed_order"])
        ),
        "zero_mode_response_constant_in_radius": bool(
            zero_mode_constancy
            and max(zero_mode_constancy)
            <= float(tolerances["max_zero_mode_constancy_relative_spread"])
        ),
        "shape_matches_i1_over_r": bool(
            shape_values
            and max(shape_values) <= float(tolerances["max_shape_one_minus_cosine"])
        ),
        "shape_orders_in_band": bool(
            shape_orders
            and min(shape_orders) >= float(tolerances["min_shape_observed_order"])
            and max(shape_orders) <= float(tolerances["max_shape_observed_order"])
        ),
        "exponential_regime_slope_matches_two_k": bool(
            exponential_deviations
            and max(exponential_deviations)
            <= float(tolerances["max_exponential_slope_relative_deviation"])
        ),
        "algebraic_regime_slope_matches_minus_two": bool(
            algebraic_deviations
            and max(algebraic_deviations)
            <= float(tolerances["max_algebraic_slope_absolute_deviation"])
        ),
        "measured_slopes_match_exact_oracle": bool(
            slope_errors
            and max(slope_errors)
            <= float(tolerances["max_measured_versus_oracle_slope_relative_error"])
        ),
        "modal_ode_residual_small": bool(
            float(selfcheck["max_modal_ode_relative_residual"])
            <= float(tolerances["max_modal_ode_relative_residual"])
        ),
        "modal_ode_residual_converges": bool(
            float(selfcheck["min_modal_ode_observed_order"])
            >= float(tolerances["min_modal_ode_observed_order"])
        ),
        "bessel_wronskian_identity_holds": bool(
            float(selfcheck["max_wronskian_relative_defect"])
            <= float(tolerances["max_wronskian_relative_defect"])
        ),
        "bessel_published_values_reproduced": bool(
            float(selfcheck["max_published_relative_error"])
            <= float(tolerances["max_published_value_relative_error"])
        ),
        "cross_solver_agreement_within_band": bool(
            float(cross["max_relative_difference"])
            <= float(tolerances["max_cross_solver_relative_difference"])
        ),
        "cross_solver_ab_gap_is_genuine_truncation": bool(
            float(cross["min_ab_difference_over_dr_squared"])
            >= float(tolerances["min_cross_solver_ab_difference_over_dr_squared"])
            and float(cross["max_ab_difference_over_dr_squared"])
            <= float(tolerances["max_cross_solver_ab_difference_over_dr_squared"])
        ),
        "cross_solver_ratios_agree": bool(
            float(cross["max_solver_ratio_relative_spread"])
            <= float(tolerances["max_cross_solver_ratio_relative_deviation"])
        ),
        "nonlinear_cross_link_recorded": bool(
            analysis["nonlinear_cross_link"].get("available", False)
        ),
    }
    checks["all_passed"] = all(checks.values())
    return {
        "checks": checks,
        "observed": {
            "max_oracle_ratio_relative_error_finest": max(ratio_errors_finest)
            if ratio_errors_finest
            else None,
            "max_oracle_ratio_relative_error_any": max(ratio_errors_any)
            if ratio_errors_any
            else None,
            "oracle_ratio_observed_orders": ratio_orders,
            "max_zero_mode_closed_form_relative_error_finest": max(zero_mode_finest)
            if zero_mode_finest
            else None,
            "zero_mode_closed_form_observed_orders": zero_mode_orders,
            "max_zero_mode_constancy_relative_spread": max(zero_mode_constancy)
            if zero_mode_constancy
            else None,
            "max_shape_one_minus_cosine": max(shape_values) if shape_values else None,
            "shape_observed_orders": shape_orders,
            "max_exponential_slope_relative_deviation": max(exponential_deviations)
            if exponential_deviations
            else None,
            "max_algebraic_slope_absolute_deviation": max(algebraic_deviations)
            if algebraic_deviations
            else None,
            "max_measured_versus_oracle_slope_relative_error": max(slope_errors)
            if slope_errors
            else None,
            "exponential_pair_count": len(exponential_deviations),
            "algebraic_pair_count": len(algebraic_deviations),
            "minimum_interior_zero_source_rows": minimum_zero_source_rows,
            "comparable_pair_count": len(ratio_errors_any),
            "shape_sample_count": len(shape_values),
        },
    }


RESPONSE_CSV_FIELDS: tuple[str, ...] = (
    "axial_case",
    "points_per_unit_radius",
    "wall_radius",
    "kr",
    "nr",
    "response_max_abs",
    "response_relative_to_reference_core_max",
    "argmax_r",
    "argmax_z",
    "response_spread_over_abs_mean",
    "oracle_amplitude",
    "oracle_predicted_response",
    "floor_margin",
    "oracle_comparable",
)

PAIR_CSV_FIELDS: tuple[str, ...] = (
    "axial_case",
    "points_per_unit_radius",
    "smaller_radius",
    "larger_radius",
    "smaller_kr",
    "larger_kr",
    "regime",
    "measured_ratio",
    "oracle_ratio",
    "ratio_relative_error",
    "oracle_comparable",
)


def _write_csv(
    path: Path, fieldnames: Sequence[str], rows: Iterable[dict[str, Any]]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name) for name in fieldnames})


def run(config: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """Run the experiment and write the evidence bundle."""

    validate_config(config)
    _prepare_output(output_dir)
    provenance = collect_runtime_provenance()

    analysis = evaluate(config)
    stored_arrays = analysis.pop("stored_arrays")
    acceptance = acceptance_report(analysis, config)

    summary = {
        "schema_version": 1,
        "experiment_id": EXPERIMENT_ID,
        "audit_entry": "E-33",
        "statement": (
            "floating-point verification of the E-33 wall-truncation response "
            "law for the E-25 finite-cylinder elliptic problem; not a continuum "
            "bound and not singularity evidence"
        ),
        "equation": "-(d_rr + 3/r d_r + d_zz) psi1 = omega1, psi1(R,z)=0",
        "description": config["description"],
        "interpretation": config["interpretation"],
        "claims_tested": {
            "E-33a": (
                "modal reduction to the order-one modified Bessel equation; "
                "checked by the finite-difference residual of L_{5,k} applied "
                "to I_1(kr)/r and K_1(kr)/r"
            ),
            "E-33b": (
                "wall response proportional to K_1(kR)/I_1(kR) with radial "
                "shape I_1(kr)/r; checked by the oracle ratio agreement and the "
                "cosine similarity against the predicted shape"
            ),
            "E-33c": (
                "exponential decay e^{-2kR} for kR >> 1 and algebraic R^-2 for "
                "kR << 1 with the crossover at kR ~ 1; checked by the fitted "
                "local slopes and by the asymptotic branches of the oracle"
            ),
            "E-33d": (
                "the k=0 wall-dependent part is the constant -Q_inf/(2R^2); "
                "checked by the r-spread of the measured difference and by its "
                "agreement with the closed form"
            ),
        },
        "analysis": analysis,
        "acceptance": acceptance,
        "accepted": bool(acceptance["checks"]["all_passed"]),
        "limitations": list(LIMITATIONS),
        "known_gaps": list(KNOWN_GAPS),
        "reproducibility": {"runtime_provenance": provenance},
    }

    summary_path = output_dir / "summary.json"
    config_path = output_dir / "config.snapshot.json"
    response_path = output_dir / "wall_response.csv"
    pair_path = output_dir / "radius_pairs.csv"
    arrays_path = output_dir / "core_response.npz"

    _write_json(summary_path, summary)
    _write_json(config_path, config)
    _write_csv(
        response_path,
        RESPONSE_CSV_FIELDS,
        (
            entry
            for row in analysis["axial_case_measurements"]
            for entry in row["radii"]
        ),
    )
    _write_csv(
        pair_path,
        PAIR_CSV_FIELDS,
        (
            entry
            for row in analysis["axial_case_measurements"]
            for entry in row["pairs"]
        ),
    )
    np.savez(arrays_path, **stored_arrays)

    payloads = (summary_path, config_path, response_path, pair_path, arrays_path)
    manifest = {
        "schema_version": 1,
        "experiment_id": EXPERIMENT_ID,
        "files": {
            path.name: {"sha256": _sha256(path), "bytes": path.stat().st_size}
            for path in sorted(payloads, key=lambda item: item.name)
        },
    }
    manifest_path = output_dir / "manifest.json"
    _write_json(manifest_path, manifest)
    (output_dir / "manifest.json.sha256").write_text(
        _sha256(manifest_path) + "\n", encoding="ascii"
    )

    if not summary["accepted"]:
        failed = sorted(
            name for name, value in acceptance["checks"].items() if not value
        )
        raise RuntimeError(
            f"wall truncation scaling experiment failed acceptance: {failed}"
        )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Verify equation audit E-33.")
    parser.add_argument(
        "--config",
        type=Path,
        default=REPOSITORY_ROOT / "configs" / "wall_truncation_scaling.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "wall_truncation_scaling",
    )
    arguments = parser.parse_args(argv)
    config_path = arguments.config.resolve()
    output_dir = arguments.output_dir.resolve()
    if not _inside_repository(config_path) or not _inside_repository(output_dir):
        parser.error("config and output paths must remain inside this repository")
    config = strict_json_loads(
        config_path.read_text(encoding="utf-8"),
        label="wall truncation scaling config",
    )
    if not isinstance(config, dict):
        parser.error("config must be a JSON object")
    summary = run(config, output_dir)
    print(json.dumps(summary["acceptance"]["checks"], ensure_ascii=False, sort_keys=True))
    return 0 if summary["accepted"] else 2


if __name__ == "__main__":
    sys.exit(main())
