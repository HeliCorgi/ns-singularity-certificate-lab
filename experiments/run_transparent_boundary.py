r"""Acceptance evidence for milestone W-A: the exact modal transparent outer condition.

What is being tested
--------------------
``docs/whole_space_transition.md`` section 3 derives, from the E-33(b) modal
structure, the exact Dirichlet-to-Neumann condition at the artificial outer
radius ``R`` of the E-25 finite cylinder,

.. math::

   \partial_r\hat\psi_k(R)
   +\Big[\frac2R+k\,\frac{K_0(kR)}{K_1(kR)}\Big]\hat\psi_k(R)=0 ,
   \tag{W-1}

and section 4 (W-B) fixes six acceptance conditions *before* implementation.
This module measures exactly those six and nothing else:

1. ``k = 0`` exactness -- the transparent core solution is ``R``-independent up
   to discretization, and the Dirichlet offset ``(Q_inf/2)(R'^-2 - R^-2)`` of
   E-33(e) is gone;
2. ``k > 0`` exactness -- the same ``R``-independence across several ``L_z``,
   with the improvement factor over Dirichlet recorded, especially in the
   long-wavelength ``kR << 1`` regime where the Dirichlet error decays only
   like ``R^-2``;
3. manufactured second-order convergence **including the boundary row**;
4. agreement with a large-radius Dirichlet solve, split into its ``O(dr^2)``
   and ``O(R_big^-2)`` parts;
5. detection of three named faults in the bracket;
6. refusal of a source that is not compactly supported strictly inside ``R``.

Scope warning
-------------
The condition implemented here belongs to the **whole-space** problem.  The Hou
reproduction keeps the *physical* no-slip wall of E-27, where
``psi_1(t, 1, z) = 0`` is a statement about a real cylinder and not a truncation
artefact.  The two are never selected implicitly: the solver defaults to
Dirichlet, the transparent path must be named, and
``summary.json`` records ``outer_boundary_condition`` for every measurement.

Why the ``k = 0`` and ``k > 0`` exactness claims are not the same claim
----------------------------------------------------------------------
For ``k = 0`` the Dirichlet wall error is the ``r``-independent constant
``-Q_inf/(2R^2)`` (E-33(e)) and is known in closed form from the source moment
alone, so condition 1 can be checked against an exact number.  For ``k > 0``
the wall error is ``-A [K_1(kR)/I_1(kR)] I_1(kr)/r`` with an amplitude ``A``
that this experiment never needs: only *ratios* of wall responses are used, and
those are amplitude-free.

What this experiment does not establish
---------------------------------------
Nothing about the continuum, nothing about singularity formation.  (W-1) is
exact only under compact support of ``omega_1``; the discretized condition still
carries an ``O(dr^2)`` boundary error, which is what conditions 1-4 measure.  A
transparent outer condition makes the *elliptic truncation* honest; it says
nothing about axis resolution, about the nonlinear evolution, or about the Clay
problem.
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
from ns_certificate_lab.grid import AxisymmetricGrid
from ns_certificate_lab.poisson import solve_streamfunction_poisson
from ns_certificate_lab.provenance import collect_runtime_provenance
from ns_certificate_lab.transparent_boundary import (
    BRACKET_VARIANTS,
    DIRICHLET,
    TRANSPARENT,
    assert_compact_support,
    outer_bracket,
    solve_radial_mode,
    solve_streamfunction_poisson_outer,
)

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

EXPERIMENT_ID = "transparent_boundary_v1"
EXPECTED_SCHEMA_VERSION = 1
MINIMUM_RESOLUTIONS = 3
MINIMUM_WALL_RADII = 2

#: Tolerance used when selecting core nodes; the grids share ``dr`` exactly, so
#: the core node coordinates are bitwise identical and this only guards the
#: comparison ``r_i <= core_radius`` at the last included node.
CORE_TOLERANCE = 1.0e-12

#: Faults that must be detected (W-B condition 5), in the order they are
#: reported.  ``"exact"`` is deliberately excluded: it is the baseline.
INJECTED_FAULTS: tuple[str, ...] = (
    "sign_flipped",
    "no_curvature_term",
    "frozen_ratio",
)

TOP_LEVEL_KEYS: frozenset[str] = frozenset(
    {
        "schema_version",
        "experiment",
        "description",
        "interpretation",
        "physical_wall_warning",
        "source_profile",
        "core_radius",
        "axial_cases",
        "wall_radii",
        "axial_points",
        "radial_resolutions",
        "manufactured",
        "large_radius_dirichlet",
        "fault_injection",
        "support_leak",
        "long_wavelength_max_kr",
        "acceptance",
    }
)

SOURCE_PROFILE_KEYS: frozenset[str] = frozenset({"support_radius", "exponent"})
AXIAL_CASE_KEYS: frozenset[str] = frozenset({"label", "z_period", "mode"})
MANUFACTURED_KEYS: frozenset[str] = frozenset(
    {"wall_radius", "radial_resolutions", "axial_case_labels"}
)
LARGE_RADIUS_KEYS: frozenset[str] = frozenset(
    {"transparent_wall_radius", "dirichlet_wall_radii"}
)
FAULT_KEYS: frozenset[str] = frozenset(
    {"points_per_unit_radius", "wall_radii", "variants"}
)
SUPPORT_LEAK_KEYS: frozenset[str] = frozenset(
    {"points_per_unit_radius", "wall_radius", "leaking_support_radius"}
)
ACCEPTANCE_KEYS: frozenset[str] = frozenset(
    {
        "max_transparent_r_independence_over_dr_squared",
        "min_transparent_r_independence_observed_order",
        "max_transparent_r_independence_observed_order",
        "max_zero_mode_dirichlet_closed_form_relative_error",
        "min_dirichlet_over_transparent_improvement_factor",
        "min_long_wavelength_improvement_factor",
        "min_manufactured_observed_order",
        "max_manufactured_observed_order",
        "min_manufactured_boundary_row_observed_order",
        "max_manufactured_boundary_row_observed_order",
        "max_manufactured_error_over_dr_squared",
        "large_radius_discretization_constant",
        "max_large_radius_dirichlet_split_factor",
        "max_zero_mode_wall_corrected_residual_over_dr_squared",
        "min_zero_mode_wall_corrected_observed_order",
        "max_zero_mode_wall_corrected_observed_order",
        "min_fault_detection_factor",
        "max_outer_condition_defect_over_dr_squared",
        "max_k0_over_k1_published_relative_error",
        "max_k0_over_k1_derivative_identity_relative_defect",
        "max_wronskian_relative_defect",
        "max_k0_over_k1_small_argument_relative_error",
        "max_k0_over_k1_large_argument_relative_error",
        "max_zero_mode_bracket_relative_defect",
    }
)

LIMITATIONS: tuple[str, ...] = (
    "binary64 arithmetic without outward rounding; every tolerance below is a "
    "measured floating-point band, not a proof",
    "(W-1) is exact only while omega_1 vanishes for r >= R; the discretized "
    "condition still carries an O(dr^2) boundary error, so 'R-independent' "
    "here means 'R-dependent only at the discretization level', never zero",
    "the transparent condition is diagonal in the axial Fourier index and "
    "therefore nonlocal in z; it has no real-space z-stencil, so solver C "
    "cannot cross-validate it and the independent path used instead is the "
    "large-radius Dirichlet solve of condition 4",
    "the radial profile (1-(r/a)^2)^p used for conditions 1, 2, 4 and 6 is "
    "compactly supported and C^(p-1), not C^infinity",
    "the manufactured family of condition 3 is exactly the decaying branch "
    "A K_1(kr)/r outside its support, matched to an even sextic inside; it "
    "satisfies (W-1) identically at every R beyond the support, which is what "
    "makes the boundary row testable, but its source is only C^1 across the "
    "support edge",
    "observed orders are computed over the listed grids only; they are not "
    "continuum error bounds",
)

KNOWN_GAPS: tuple[str, ...] = (
    "the 'no_curvature_term' fault is detected as a refusal rather than as a "
    "wrong number: dropping 2/R makes the k=0 modal matrix a singular discrete "
    "Neumann operator and the unpivoted elimination stops; the k>0 modes are "
    "measured separately through the single-mode path so that the fault also "
    "has a numerical detection factor",
    "the 'frozen_ratio' fault (K_0/K_1 := 1) is invisible at k=0 by "
    "construction, because the bracket is 2/R there for both; it is detected "
    "only on the k>0 cases, and this experiment records that limitation as a "
    "measured bitwise identity rather than hiding it",
    "no continuum error bound is produced anywhere; the whole bundle is a "
    "finite-dimensional consistency measurement",
    "condition 4 compares against a Dirichlet solve at finite R_big, whose own "
    "wall error decays only algebraically for kR << 1; for the longest "
    "wavelength case the reference is therefore not converged and the split "
    "reports the wall term as dominant, which is the point of E-33 rather than "
    "a defect of the transparent condition",
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
    if not (coarse > 0.0 and fine > 0.0):
        raise ValueError("observed orders need strictly positive error levels")
    return math.log(coarse / fine) / math.log(coarse_h / fine_h)


def _max_abs(values: np.ndarray) -> float:
    return float(np.max(np.abs(values)))


def _finite(values: Iterable[Any]) -> bool:
    return all(math.isfinite(float(value)) for value in values)


def _observed_orders(errors: Sequence[float], spacings: Sequence[float]) -> list[float]:
    return [
        _order(errors[index], errors[index + 1], spacings[index], spacings[index + 1])
        for index in range(len(errors) - 1)
        if errors[index] > 0.0 and errors[index + 1] > 0.0
    ]


# ------------------------------------------------------- configuration


def radial_point_count(points_per_unit_radius: int, wall_radius: float) -> int:
    """Return ``nr`` such that ``dr = 1/points_per_unit_radius`` exactly.

    Every wall radius in a comparison shares the radial spacing, so the core
    nodes ``r_i = i/points_per_unit_radius`` are the *same* floating-point
    numbers on every member of the sweep and the core differences are exact
    nodewise comparisons with no interpolation.  A fractional product is
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

    config = _exact_keys(config, TOP_LEVEL_KEYS, name="transparent boundary config")
    if _integer(config["schema_version"], name="schema_version") != (
        EXPECTED_SCHEMA_VERSION
    ):
        raise ValueError("unsupported transparent boundary config schema")
    if config["experiment"] != EXPERIMENT_ID:
        raise ValueError(
            f"experiment must equal the audited canonical value {EXPERIMENT_ID!r}"
        )
    for key in ("description", "interpretation", "physical_wall_warning"):
        if not isinstance(config[key], str) or not config[key].strip():
            raise ValueError(f"{key} must be a non-empty string")
    if "E-27" not in config["physical_wall_warning"]:
        raise ValueError(
            "physical_wall_warning must name E-27 explicitly: the transparent "
            "condition is for the whole-space problem and must never be "
            "confused with the Hou physical no-slip wall"
        )

    profile = _exact_keys(
        config["source_profile"], SOURCE_PROFILE_KEYS, name="source_profile"
    )
    support_radius = _number(
        profile["support_radius"], name="source_profile.support_radius", positive=True
    )
    _integer(profile["exponent"], name="source_profile.exponent", minimum=1)

    core_radius = _number(config["core_radius"], name="core_radius", positive=True)
    if core_radius >= support_radius:
        raise ValueError("core_radius must lie strictly inside the source support")

    radii = config["wall_radii"]
    if not isinstance(radii, list) or len(radii) < MINIMUM_WALL_RADII:
        raise ValueError(f"wall_radii must list at least {MINIMUM_WALL_RADII} radii")
    wall_radii = [
        _number(item, name=f"wall_radii[{index}]", positive=True)
        for index, item in enumerate(radii)
    ]
    if any(
        wall_radii[index] >= wall_radii[index + 1]
        for index in range(len(wall_radii) - 1)
    ):
        raise ValueError("wall_radii must be strictly increasing")
    if support_radius >= wall_radii[0]:
        raise ValueError(
            "source_profile.support_radius must be strictly smaller than the "
            "smallest wall radius, otherwise (W-1) is not exact and every "
            "R-independence number below would be measuring something else"
        )

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
        if _integer(case["mode"], name=f"axial_cases[{index}].mode", minimum=0) == 0:
            zero_mode_count += 1
    if len(set(labels)) != len(labels):
        raise ValueError("axial_cases labels must be unique")
    if zero_mode_count != 1:
        raise ValueError(
            "exactly one axial case must have mode 0; it carries W-B condition 1"
        )
    if len(labels) < 2:
        raise ValueError(
            "at least one k>0 case is required for W-B condition 2"
        )

    axial_points = _integer(config["axial_points"], name="axial_points", minimum=5)
    if axial_points % 2 != 0:
        raise ValueError("axial_points must be even so the FFT has a Nyquist mode")

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

    manufactured = _exact_keys(
        config["manufactured"], MANUFACTURED_KEYS, name="manufactured"
    )
    manufactured_radius = _number(
        manufactured["wall_radius"], name="manufactured.wall_radius", positive=True
    )
    if manufactured_radius <= support_radius:
        raise ValueError(
            "manufactured.wall_radius must lie strictly outside the source "
            "support, otherwise the manufactured solution does not satisfy "
            "(W-1) at the wall"
        )
    manufactured_resolutions = manufactured["radial_resolutions"]
    if (
        not isinstance(manufactured_resolutions, list)
        or len(manufactured_resolutions) < MINIMUM_RESOLUTIONS
    ):
        raise ValueError(
            "manufactured.radial_resolutions needs at least "
            f"{MINIMUM_RESOLUTIONS} entries to produce observed orders"
        )
    manufactured_ppur = [
        _integer(item, name=f"manufactured.radial_resolutions[{index}]", minimum=8)
        for index, item in enumerate(manufactured_resolutions)
    ]
    if any(
        manufactured_ppur[index] >= manufactured_ppur[index + 1]
        for index in range(len(manufactured_ppur) - 1)
    ):
        raise ValueError("manufactured.radial_resolutions must be strictly increasing")
    for ppur in manufactured_ppur:
        radial_point_count(ppur, manufactured_radius)
    manufactured_labels = manufactured["axial_case_labels"]
    if not isinstance(manufactured_labels, list) or not manufactured_labels:
        raise ValueError("manufactured.axial_case_labels must be a non-empty list")
    if not set(manufactured_labels) <= set(labels):
        raise ValueError("manufactured.axial_case_labels must name declared cases")

    large = _exact_keys(
        config["large_radius_dirichlet"],
        LARGE_RADIUS_KEYS,
        name="large_radius_dirichlet",
    )
    transparent_radius = _number(
        large["transparent_wall_radius"],
        name="large_radius_dirichlet.transparent_wall_radius",
        positive=True,
    )
    if transparent_radius not in wall_radii:
        raise ValueError(
            "large_radius_dirichlet.transparent_wall_radius must be one of "
            "wall_radii so it shares their radial spacing"
        )
    reference_radii_raw = large["dirichlet_wall_radii"]
    if not isinstance(reference_radii_raw, list) or len(reference_radii_raw) < 2:
        raise ValueError(
            "large_radius_dirichlet.dirichlet_wall_radii must list at least two "
            "radii so the O(R_big^-2) part of the split can be separated"
        )
    reference_radii = [
        _number(
            item,
            name=f"large_radius_dirichlet.dirichlet_wall_radii[{index}]",
            positive=True,
        )
        for index, item in enumerate(reference_radii_raw)
    ]
    if any(
        reference_radii[index] >= reference_radii[index + 1]
        for index in range(len(reference_radii) - 1)
    ):
        raise ValueError(
            "large_radius_dirichlet.dirichlet_wall_radii must be strictly increasing"
        )
    if reference_radii[0] <= transparent_radius:
        raise ValueError(
            "every Dirichlet reference radius must exceed the transparent "
            "wall radius; otherwise it is not a large-radius reference"
        )
    for ppur in points_per_unit:
        for radius in reference_radii:
            radial_point_count(ppur, radius)

    faults = _exact_keys(config["fault_injection"], FAULT_KEYS, name="fault_injection")
    fault_ppur = _integer(
        faults["points_per_unit_radius"],
        name="fault_injection.points_per_unit_radius",
        minimum=8,
    )
    fault_radii_raw = faults["wall_radii"]
    if not isinstance(fault_radii_raw, list) or len(fault_radii_raw) != 2:
        raise ValueError(
            "fault_injection.wall_radii must list exactly the two radii whose "
            "core difference is the detection quantity"
        )
    fault_radii = [
        _number(item, name=f"fault_injection.wall_radii[{index}]", positive=True)
        for index, item in enumerate(fault_radii_raw)
    ]
    if fault_radii[0] >= fault_radii[1]:
        raise ValueError("fault_injection.wall_radii must be strictly increasing")
    if not set(fault_radii) <= set(wall_radii):
        raise ValueError("fault_injection.wall_radii must be a subset of wall_radii")
    for radius in fault_radii:
        radial_point_count(fault_ppur, radius)
    variants = faults["variants"]
    if not isinstance(variants, list) or list(variants) != list(INJECTED_FAULTS):
        raise ValueError(
            "fault_injection.variants must be exactly "
            f"{list(INJECTED_FAULTS)!r}; W-B condition 5 names all three and "
            "silently dropping one would make the check vacuous"
        )
    for variant in variants:
        if variant not in BRACKET_VARIANTS or variant == "exact":
            raise ValueError(f"unknown injected fault variant {variant!r}")

    leak = _exact_keys(config["support_leak"], SUPPORT_LEAK_KEYS, name="support_leak")
    leak_ppur = _integer(
        leak["points_per_unit_radius"],
        name="support_leak.points_per_unit_radius",
        minimum=8,
    )
    leak_radius = _number(
        leak["wall_radius"], name="support_leak.wall_radius", positive=True
    )
    radial_point_count(leak_ppur, leak_radius)
    leaking_support = _number(
        leak["leaking_support_radius"],
        name="support_leak.leaking_support_radius",
        positive=True,
    )
    if leaking_support <= leak_radius:
        raise ValueError(
            "support_leak.leaking_support_radius must reach the wall or beyond; "
            "otherwise the refusal being tested cannot be triggered"
        )

    long_wavelength = _number(
        config["long_wavelength_max_kr"], name="long_wavelength_max_kr", positive=True
    )
    if long_wavelength >= 1.0:
        raise ValueError(
            "long_wavelength_max_kr must be below the E-33(c) crossover kR ~ 1"
        )

    acceptance = _exact_keys(config["acceptance"], ACCEPTANCE_KEYS, name="acceptance")
    for key in sorted(ACCEPTANCE_KEYS):
        _number(acceptance[key], name=f"acceptance.{key}")
    _number(
        acceptance["large_radius_discretization_constant"],
        name="acceptance.large_radius_discretization_constant",
        positive=True,
    )
    if _number(
        acceptance["min_fault_detection_factor"],
        name="acceptance.min_fault_detection_factor",
    ) <= 1.0:
        raise ValueError(
            "acceptance.min_fault_detection_factor must exceed 1; a threshold "
            "of one or less would call any fault 'strongly detected'"
        )
    return config


# ------------------------------------------------------------- the source


def radial_profile(
    radii: np.ndarray, *, support_radius: float, exponent: int
) -> np.ndarray:
    """Return ``(1 - (r/a)^2)^p`` inside ``r < a`` and exactly zero outside."""

    values = np.zeros_like(radii, dtype=np.float64)
    inside = radii < support_radius
    scaled = radii[inside] / support_radius
    values[inside] = (1.0 - scaled * scaled) ** exponent
    return values


def source_moment(*, support_radius: float, exponent: int) -> float:
    r"""Return ``Q_inf = int_0^inf s^3 omega_1(s) ds = a^4/(2(p+1)(p+2))``.

    Derived in ``experiments/run_wall_truncation_scaling.py`` and reproduced
    here because the E-33(e) Dirichlet offset ``(Q_inf/2)(R'^-2 - R^-2)`` is the
    number condition 1 must show has disappeared.
    """

    return support_radius**4 / (2.0 * (exponent + 1) * (exponent + 2))


def wavenumber(*, z_period: float, mode: int) -> float:
    return 2.0 * math.pi * mode / z_period


def build_grid(
    *, points_per_unit_radius: int, wall_radius: float, axial_points: int,
    z_period: float,
) -> AxisymmetricGrid:
    return AxisymmetricGrid.uniform(
        nr=radial_point_count(points_per_unit_radius, wall_radius),
        nz=axial_points,
        r_max=wall_radius,
        z_min=0.0,
        z_max=z_period,
        periodic_z=True,
    )


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


def core_mask(grid: AxisymmetricGrid, core_radius: float) -> np.ndarray:
    return grid.r <= core_radius + CORE_TOLERANCE


# --------------------------------------------- the manufactured family


def decaying_branch_derivatives(k: float, radius: float) -> np.ndarray:
    r"""Return ``u, u', u'', u'''`` for the decaying branch ``u = K_1(kr)/r``.

    Differentiating with ``K_0' = -K_1`` and ``K_1'(x) = -K_0(x) - K_1(x)/x``
    (DLMF 10.29.2) gives

    .. math::

       u'   &= -\frac{kK_0}{r}-\frac{2K_1}{r^2}, \\
       u''  &= \frac{k^2K_1}{r}+\frac{3kK_0}{r^2}+\frac{6K_1}{r^3}, \\
       u''' &= -\frac{k^3K_0}{r}-\frac{5k^2K_1}{r^2}
               -\frac{12kK_0}{r^3}-\frac{24K_1}{r^4},

    with all Bessel functions evaluated at ``kr``.  The first line is exactly
    the derivative identity from which (W-1) is read off, so this function is
    also where a sign error in the derivation would surface; the tests check
    ``u'' + (3/r)u' - k^2 u = 0`` and check every line against high-order finite
    differences.  ``k = 0`` uses the elementary branch ``u = r^-2``.
    """

    if k == 0.0:
        return np.array(
            [
                radius**-2.0,
                -2.0 * radius**-3.0,
                6.0 * radius**-4.0,
                -24.0 * radius**-5.0,
            ],
            dtype=np.float64,
        )
    argument = k * radius
    k0 = bessel_reference.bessel_k0(argument)
    k1 = bessel_reference.bessel_k1(argument)
    return np.array(
        [
            k1 / radius,
            -k * k0 / radius - 2.0 * k1 / radius**2,
            k * k * k1 / radius + 3.0 * k * k0 / radius**2 + 6.0 * k1 / radius**3,
            -(k**3) * k0 / radius
            - 5.0 * k * k * k1 / radius**2
            - 12.0 * k * k0 / radius**3
            - 24.0 * k1 / radius**4,
        ],
        dtype=np.float64,
    )


def manufactured_polynomial_coefficients(k: float, support_radius: float) -> np.ndarray:
    r"""Return ``c_0..c_3`` of the even sextic matched to the decaying branch.

    The manufactured solution is

    .. math::

       \psi^*(r)=\begin{cases}
       \sum_{j=0}^{3}c_j r^{2j}, & r<a,\\[2pt]
       u(r)/u(a), & r\ge a,
       \end{cases}
       \qquad u(r)=K_1(kr)/r ,

    with ``c`` fixed by matching ``psi^*, psi^{*\prime}, psi^{*\prime\prime},
    psi^{*\prime\prime\prime}`` at ``r = a``.  Only even powers appear, so
    ``psi^*`` is even at the axis as E-16 requires.

    Two properties make this the right manufactured family for W-B condition 3.
    Outside ``a`` the solution *is* the decaying branch, so it satisfies (W-1)
    **identically** at every wall radius beyond ``a`` -- the boundary row is
    therefore being tested against an exact boundary condition, not against an
    approximation.  And ``-L5 psi^*`` vanishes identically outside ``a``
    (the branch is homogeneous) and is continuous at ``a`` (matching three
    derivatives forces ``psi^{*\prime\prime}+(3/a)psi^{*\prime}-k^2psi^*`` to
    agree with the branch's zero), so the source really is compactly supported.
    """

    derivatives = decaying_branch_derivatives(k, support_radius)
    normalized = derivatives / derivatives[0]
    matrix = np.zeros((4, 4), dtype=np.float64)
    for order in range(4):
        for column in range(4):
            power = 2 * column
            coefficient = 1.0
            for step in range(order):
                coefficient *= power - step
            exponent = power - order
            matrix[order, column] = (
                0.0 if exponent < 0 else coefficient * support_radius**exponent
            )
    return np.linalg.solve(matrix, normalized)


def manufactured_fields(
    k: float, radii: np.ndarray, *, support_radius: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return the radial factors ``(psi*, omega* = -L5 psi*)`` on ``radii``."""

    coefficients = manufactured_polynomial_coefficients(k, support_radius)
    psi = np.empty_like(radii, dtype=np.float64)
    omega = np.zeros_like(radii, dtype=np.float64)

    inside = radii < support_radius
    inner_radii = radii[inside]
    inner_psi = np.zeros_like(inner_radii)
    inner_laplacian = np.zeros_like(inner_radii)
    for power_index, coefficient in enumerate(coefficients):
        inner_psi += coefficient * inner_radii ** (2 * power_index)
        if power_index >= 1:
            # L5(r^{2j}) = 4 j (j+1) r^{2j-2}, verified in the tests.
            inner_laplacian += (
                coefficient
                * 4.0
                * power_index
                * (power_index + 1)
                * inner_radii ** (2 * power_index - 2)
            )
    psi[inside] = inner_psi
    omega[inside] = -(inner_laplacian - k * k * inner_psi)

    outer_radii = radii[~inside]
    scale = decaying_branch_derivatives(k, support_radius)[0]
    if k == 0.0:
        psi[~inside] = outer_radii**-2.0 / scale
    else:
        psi[~inside] = (
            np.array(
                [
                    bessel_reference.bessel_k1(k * float(value)) / float(value)
                    for value in outer_radii
                ],
                dtype=np.float64,
            )
            / scale
        )
    return psi, omega


def manufactured_case_fields(
    grid: AxisymmetricGrid,
    *,
    support_radius: float,
    z_period: float,
    mode: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the two-dimensional ``(psi*, omega*)`` for one axial case."""

    k = wavenumber(z_period=z_period, mode=mode)
    radial_psi, radial_omega = manufactured_fields(
        k, grid.r, support_radius=support_radius
    )
    if mode == 0:
        axial = np.ones(grid.nz, dtype=np.float64)
    else:
        axial = np.cos(2.0 * np.pi * mode * grid.z / z_period)
    return radial_psi[:, None] * axial[None, :], radial_omega[:, None] * axial[None, :]


# ----------------------------------------------------------- the solvers


def solve_transparent(
    grid: AxisymmetricGrid,
    omega: np.ndarray,
    *,
    support_radius: float,
    bracket_variant: str = "exact",
) -> Any:
    return solve_streamfunction_poisson_outer(
        grid,
        omega,
        boundary_condition=TRANSPARENT,
        support_radius=support_radius,
        bracket_variant=bracket_variant,
        estimate_condition=False,
    )


def solve_dirichlet(grid: AxisymmetricGrid, omega: np.ndarray) -> np.ndarray:
    return solve_streamfunction_poisson(
        grid, omega, 0.0, estimate_condition=False
    ).psi1


# --------------------------------------- W-B conditions 1 and 2: exactness


def measure_r_independence(config: dict[str, Any]) -> list[dict[str, Any]]:
    """Measure the core ``R``-dependence of both outer conditions."""

    profile = config["source_profile"]
    support_radius = float(profile["support_radius"])
    exponent = int(profile["exponent"])
    core_radius = float(config["core_radius"])
    wall_radii = [float(value) for value in config["wall_radii"]]
    axial_points = int(config["axial_points"])
    moment = source_moment(support_radius=support_radius, exponent=exponent)

    rows: list[dict[str, Any]] = []
    for case in config["axial_cases"]:
        label = str(case["label"])
        z_period = float(case["z_period"])
        mode = int(case["mode"])
        k = wavenumber(z_period=z_period, mode=mode)
        for ppur in [int(value) for value in config["radial_resolutions"]]:
            spacing = 1.0 / ppur
            transparent_core: dict[float, np.ndarray] = {}
            dirichlet_core: dict[float, np.ndarray] = {}
            defects: dict[float, float] = {}
            zero_source_rows = 0
            for radius in wall_radii:
                grid = build_grid(
                    points_per_unit_radius=ppur,
                    wall_radius=radius,
                    axial_points=axial_points,
                    z_period=z_period,
                )
                omega = build_source(
                    grid,
                    support_radius=support_radius,
                    exponent=exponent,
                    z_period=z_period,
                    mode=mode,
                )
                zero_source_rows = assert_compact_support(
                    grid, omega, support_radius=support_radius
                )
                transparent = solve_transparent(
                    grid, omega, support_radius=support_radius
                )
                mask = core_mask(grid, core_radius)
                transparent_core[radius] = transparent.psi1[mask]
                dirichlet_core[radius] = solve_dirichlet(grid, omega)[mask]
                defects[radius] = float(
                    transparent.metadata["outer_condition_defect_max_abs"]
                )
            base = wall_radii[0]
            for radius in wall_radii[1:]:
                transparent_difference = _max_abs(
                    transparent_core[radius] - transparent_core[base]
                )
                dirichlet_difference = _max_abs(
                    dirichlet_core[radius] - dirichlet_core[base]
                )
                # Sign convention of E-33(e) and of docs/equation_audit.md:
                # the difference is taken as psi^(R) - psi^(R') with R < R', so
                # the quoted zero-mode number is negative and is directly
                # comparable with the recorded -1.6969516537e-3 at R = 1 -> 2.
                signed_dirichlet = float(
                    np.mean(dirichlet_core[base] - dirichlet_core[radius])
                )
                signed_transparent = float(
                    np.mean(transparent_core[base] - transparent_core[radius])
                )
                closed_form = (
                    0.5 * moment * (radius**-2.0 - base**-2.0) if mode == 0 else None
                )
                rows.append(
                    {
                        "axial_case": label,
                        "z_period": z_period,
                        "mode": mode,
                        "wavenumber": k,
                        "points_per_unit_radius": ppur,
                        "dr": spacing,
                        "smaller_radius": base,
                        "larger_radius": radius,
                        "kr_smaller": k * base,
                        "interior_zero_source_rows": zero_source_rows,
                        "transparent_core_difference": transparent_difference,
                        "dirichlet_core_difference": dirichlet_difference,
                        "transparent_over_dr_squared": (
                            transparent_difference / spacing**2
                        ),
                        "improvement_factor": (
                            dirichlet_difference / transparent_difference
                            if transparent_difference > 0.0
                            else None
                        ),
                        "signed_dirichlet_core_difference": signed_dirichlet,
                        "signed_transparent_core_difference": signed_transparent,
                        "zero_mode_closed_form": closed_form,
                        "zero_mode_closed_form_relative_error": (
                            abs(signed_dirichlet - closed_form) / abs(closed_form)
                            if closed_form
                            else None
                        ),
                        "outer_condition_defect_max_abs": defects[base],
                        "outer_condition_defect_over_dr_squared": (
                            defects[base] / spacing**2
                        ),
                    }
                )
    return rows


def r_independence_convergence(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    """Return observed orders of the transparent ``R``-dependence."""

    grouped: dict[tuple[str, float, float], list[dict[str, Any]]] = {}
    for row in rows:
        key = (row["axial_case"], row["smaller_radius"], row["larger_radius"])
        grouped.setdefault(key, []).append(row)
    reports: list[dict[str, Any]] = []
    for (label, smaller, larger), entries in sorted(grouped.items()):
        ordered = sorted(entries, key=lambda item: item["points_per_unit_radius"])
        spacings = [float(item["dr"]) for item in ordered]
        errors = [float(item["transparent_core_difference"]) for item in ordered]
        reports.append(
            {
                "axial_case": label,
                "smaller_radius": smaller,
                "larger_radius": larger,
                "spacings": spacings,
                "transparent_core_differences": errors,
                "orders": _observed_orders(errors, spacings),
            }
        )
    return reports


# ---------------------------------- W-B condition 3: manufactured orders


def measure_manufactured(config: dict[str, Any]) -> list[dict[str, Any]]:
    """Measure the manufactured error including the outermost radial row."""

    support_radius = float(config["source_profile"]["support_radius"])
    settings = config["manufactured"]
    wall_radius = float(settings["wall_radius"])
    axial_points = int(config["axial_points"])
    selected = set(settings["axial_case_labels"])

    rows: list[dict[str, Any]] = []
    for case in config["axial_cases"]:
        label = str(case["label"])
        if label not in selected:
            continue
        z_period = float(case["z_period"])
        mode = int(case["mode"])
        errors: list[float] = []
        boundary_errors: list[float] = []
        spacings: list[float] = []
        case_rows: list[dict[str, Any]] = []
        for ppur in [int(value) for value in settings["radial_resolutions"]]:
            grid = build_grid(
                points_per_unit_radius=ppur,
                wall_radius=wall_radius,
                axial_points=axial_points,
                z_period=z_period,
            )
            exact, omega = manufactured_case_fields(
                grid, support_radius=support_radius, z_period=z_period, mode=mode
            )
            solution = solve_transparent(grid, omega, support_radius=support_radius)
            spacing = 1.0 / ppur
            all_rows_error = _max_abs(solution.psi1 - exact)
            boundary_error = _max_abs(solution.psi1[-1] - exact[-1])
            errors.append(all_rows_error)
            boundary_errors.append(boundary_error)
            spacings.append(spacing)
            case_rows.append(
                {
                    "axial_case": label,
                    "z_period": z_period,
                    "mode": mode,
                    "wavenumber": wavenumber(z_period=z_period, mode=mode),
                    "points_per_unit_radius": ppur,
                    "dr": spacing,
                    "wall_radius": wall_radius,
                    "max_error_all_rows": all_rows_error,
                    "max_error_boundary_row": boundary_error,
                    "max_error_over_dr_squared": all_rows_error / spacing**2,
                    "solution_scale": _max_abs(exact),
                    "outer_condition_defect_max_abs": float(
                        solution.metadata["outer_condition_defect_max_abs"]
                    ),
                    "observed_order_all_rows": None,
                    "observed_order_boundary_row": None,
                }
            )
        # Each order is attached to the *finer* grid of the pair that produced
        # it, so the coarsest row of every case reports ``None`` rather than a
        # borrowed number.
        for index in range(1, len(case_rows)):
            case_rows[index]["observed_order_all_rows"] = _order(
                errors[index - 1], errors[index], spacings[index - 1], spacings[index]
            )
            case_rows[index]["observed_order_boundary_row"] = _order(
                boundary_errors[index - 1],
                boundary_errors[index],
                spacings[index - 1],
                spacings[index],
            )
        rows.extend(case_rows)
    return rows


# ------------------------- W-B condition 4: large-radius Dirichlet split


def dirichlet_wall_factor(k: float, radius: float) -> float:
    """Return the exact ``R``-dependence of the Dirichlet wall error.

    E-33(b) gives ``K_1(kR)/I_1(kR)`` for ``k > 0``; E-33(e) gives the ``k = 0``
    analogue ``1/(2R^2)``, whose amplitude is the source moment ``Q_inf``.  Only
    *ratios* of this factor are used below, so the differing normalizations of
    the two branches never enter.
    """

    if k == 0.0:
        return 0.5 / (radius * radius)
    return bessel_reference.k1_over_i1(k * radius)


def measure_large_radius_agreement(config: dict[str, Any]) -> list[dict[str, Any]]:
    """Compare a small-``R`` transparent solve with large-``R`` Dirichlet solves."""

    profile = config["source_profile"]
    support_radius = float(profile["support_radius"])
    exponent = int(profile["exponent"])
    core_radius = float(config["core_radius"])
    axial_points = int(config["axial_points"])
    settings = config["large_radius_dirichlet"]
    transparent_radius = float(settings["transparent_wall_radius"])
    reference_radii = [float(value) for value in settings["dirichlet_wall_radii"]]
    moment = source_moment(support_radius=support_radius, exponent=exponent)

    rows: list[dict[str, Any]] = []
    for case in config["axial_cases"]:
        label = str(case["label"])
        z_period = float(case["z_period"])
        mode = int(case["mode"])
        k = wavenumber(z_period=z_period, mode=mode)
        for ppur in [int(value) for value in config["radial_resolutions"]]:
            spacing = 1.0 / ppur
            small_grid = build_grid(
                points_per_unit_radius=ppur,
                wall_radius=transparent_radius,
                axial_points=axial_points,
                z_period=z_period,
            )
            small_omega = build_source(
                small_grid,
                support_radius=support_radius,
                exponent=exponent,
                z_period=z_period,
                mode=mode,
            )
            transparent_core = solve_transparent(
                small_grid, small_omega, support_radius=support_radius
            ).psi1[core_mask(small_grid, core_radius)]

            differences: dict[float, float] = {}
            corrected: dict[float, float | None] = {}
            for radius in reference_radii:
                big_grid = build_grid(
                    points_per_unit_radius=ppur,
                    wall_radius=radius,
                    axial_points=axial_points,
                    z_period=z_period,
                )
                big_omega = build_source(
                    big_grid,
                    support_radius=support_radius,
                    exponent=exponent,
                    z_period=z_period,
                    mode=mode,
                )
                dirichlet_core = solve_dirichlet(big_grid, big_omega)[
                    core_mask(big_grid, core_radius)
                ]
                differences[radius] = _max_abs(transparent_core - dirichlet_core)
                if mode == 0:
                    # E-33(e): psi_D^(Rb) = psi_inf - Q_inf/(2 Rb^2) in the core,
                    # so removing that exact constant must leave only O(dr^2).
                    corrected[radius] = _max_abs(
                        transparent_core
                        - dirichlet_core
                        - moment / (2.0 * radius * radius)
                    )
                else:
                    corrected[radius] = None

            smallest_reference = reference_radii[0]
            largest_reference = reference_radii[-1]
            # Amplitude-free bound on the wall term left in the largest
            # reference: the whole difference at the smallest reference is an
            # upper bound for its wall term, and E-33(b) fixes the ratio of wall
            # terms exactly.
            wall_bound = differences[smallest_reference] * (
                dirichlet_wall_factor(k, largest_reference)
                / dirichlet_wall_factor(k, smallest_reference)
            )
            closed_form_wall_term = (
                moment / (2.0 * largest_reference * largest_reference)
                if mode == 0
                else None
            )
            rows.append(
                {
                    "axial_case": label,
                    "z_period": z_period,
                    "mode": mode,
                    "wavenumber": k,
                    "points_per_unit_radius": ppur,
                    "dr": spacing,
                    "transparent_wall_radius": transparent_radius,
                    "reference_radii": list(reference_radii),
                    "core_differences": [differences[r] for r in reference_radii],
                    "largest_reference_radius": largest_reference,
                    "difference_at_largest_reference": differences[largest_reference],
                    "oracle_wall_term_bound_at_largest_reference": wall_bound,
                    "closed_form_wall_term_at_largest_reference": closed_form_wall_term,
                    "wall_corrected_residual": corrected[largest_reference],
                    "wall_corrected_residual_over_dr_squared": (
                        None
                        if corrected[largest_reference] is None
                        else corrected[largest_reference] / spacing**2
                    ),
                    "difference_over_dr_squared": (
                        differences[largest_reference] / spacing**2
                    ),
                }
            )
    return rows


def large_radius_convergence(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(row["axial_case"], []).append(row)
    reports: list[dict[str, Any]] = []
    for label, entries in sorted(grouped.items()):
        ordered = sorted(entries, key=lambda item: item["points_per_unit_radius"])
        spacings = [float(item["dr"]) for item in ordered]
        totals = [float(item["difference_at_largest_reference"]) for item in ordered]
        corrected = [item["wall_corrected_residual"] for item in ordered]
        reports.append(
            {
                "axial_case": label,
                "spacings": spacings,
                "difference_orders": _observed_orders(totals, spacings),
                "wall_corrected_orders": (
                    _observed_orders([float(value) for value in corrected], spacings)
                    if all(value is not None for value in corrected)
                    else []
                ),
            }
        )
    return reports


# ---------------------------------------- W-B condition 5: fault injection


def measure_fault_injection(config: dict[str, Any]) -> list[dict[str, Any]]:
    """Inject each named bracket fault and record where it is detected."""

    profile = config["source_profile"]
    support_radius = float(profile["support_radius"])
    exponent = int(profile["exponent"])
    core_radius = float(config["core_radius"])
    axial_points = int(config["axial_points"])
    settings = config["fault_injection"]
    ppur = int(settings["points_per_unit_radius"])
    smaller, larger = (float(value) for value in settings["wall_radii"])

    rows: list[dict[str, Any]] = []
    for case in config["axial_cases"]:
        label = str(case["label"])
        z_period = float(case["z_period"])
        mode = int(case["mode"])
        k = wavenumber(z_period=z_period, mode=mode)

        grids = {}
        sources = {}
        for radius in (smaller, larger):
            grid = build_grid(
                points_per_unit_radius=ppur,
                wall_radius=radius,
                axial_points=axial_points,
                z_period=z_period,
            )
            grids[radius] = grid
            sources[radius] = build_source(
                grid,
                support_radius=support_radius,
                exponent=exponent,
                z_period=z_period,
                mode=mode,
            )

        def full_field_r_dependence(variant: str) -> float | None:
            cores = []
            for radius in (smaller, larger):
                grid = grids[radius]
                solution = solve_transparent(
                    grid,
                    sources[radius],
                    support_radius=support_radius,
                    bracket_variant=variant,
                )
                cores.append(solution.psi1[core_mask(grid, core_radius)])
            return _max_abs(cores[1] - cores[0])

        def single_mode_r_dependence(variant: str) -> float | None:
            """Same quantity restricted to this case's own axial mode.

            The full solve always contains the ``k = 0`` mode, even when the
            source has none, so a fault that makes the ``k = 0`` matrix singular
            stops the whole solve.  Restricting to the case's own wavenumber is
            what lets that fault still produce a *number* on the ``k > 0`` cases
            instead of only a refusal.
            """

            cores = []
            for radius in (smaller, larger):
                grid = grids[radius]
                # cos(2 pi m z / L_z) has Fourier coefficient nz/2 at mode m.
                radial = radial_profile(
                    grid.r, support_radius=support_radius, exponent=exponent
                ).astype(np.complex128)
                solved = solve_radial_mode(
                    grid, radial, k, bracket_variant=variant
                )
                cores.append(np.real(solved)[core_mask(grid, core_radius)])
            return _max_abs(cores[1] - cores[0])

        for variant in ["exact", *INJECTED_FAULTS]:
            entry: dict[str, Any] = {
                "axial_case": label,
                "mode": mode,
                "wavenumber": k,
                "kr_smaller": k * smaller,
                "points_per_unit_radius": ppur,
                "bracket_variant": variant,
                "smaller_radius": smaller,
                "larger_radius": larger,
                "bracket_value": outer_bracket(k, smaller, variant=variant),
                "exact_bracket_value": outer_bracket(k, smaller, variant="exact"),
            }
            try:
                entry["full_solve_r_dependence"] = full_field_r_dependence(variant)
                entry["full_solve_error"] = None
            except Exception as error:  # noqa: BLE001 - the refusal is the datum
                entry["full_solve_r_dependence"] = None
                entry["full_solve_error"] = f"{type(error).__name__}: {error}"
            try:
                entry["single_mode_r_dependence"] = single_mode_r_dependence(variant)
                entry["single_mode_error"] = None
            except Exception as error:  # noqa: BLE001
                entry["single_mode_r_dependence"] = None
                entry["single_mode_error"] = f"{type(error).__name__}: {error}"
            rows.append(entry)

    baselines = {
        (row["axial_case"], "full"): row["full_solve_r_dependence"]
        for row in rows
        if row["bracket_variant"] == "exact"
    }
    baselines.update(
        {
            (row["axial_case"], "single"): row["single_mode_r_dependence"]
            for row in rows
            if row["bracket_variant"] == "exact"
        }
    )
    for row in rows:
        if row["bracket_variant"] == "exact":
            row["full_detection_factor"] = None
            row["single_mode_detection_factor"] = None
            row["detected"] = None
            row["detected_where"] = None
            continue
        full_base = baselines[(row["axial_case"], "full")]
        single_base = baselines[(row["axial_case"], "single")]
        row["full_detection_factor"] = (
            row["full_solve_r_dependence"] / full_base
            if row["full_solve_r_dependence"] is not None
            and full_base
            else None
        )
        row["single_mode_detection_factor"] = (
            row["single_mode_r_dependence"] / single_base
            if row["single_mode_r_dependence"] is not None
            and single_base
            else None
        )
        row["bracket_is_bitwise_unchanged"] = bool(
            row["bracket_value"] == row["exact_bracket_value"]
        )
        where: list[str] = []
        if row["full_solve_error"] is not None:
            where.append("full solve refused (" + row["full_solve_error"] + ")")
        if row["single_mode_error"] is not None:
            where.append("single-mode solve refused (" + row["single_mode_error"] + ")")
        if (
            row["full_detection_factor"] is not None
            and row["full_detection_factor"] > 1.0
        ):
            where.append(
                "core R-dependence of the full solve grew by "
                f"{row['full_detection_factor']:.3g}x"
            )
        if (
            row["single_mode_detection_factor"] is not None
            and row["single_mode_detection_factor"] > 1.0
        ):
            where.append(
                "core R-dependence of this case's own mode grew by "
                f"{row['single_mode_detection_factor']:.3g}x"
            )
        row["detected_where"] = "; ".join(where) if where else None
        row["detected"] = bool(where)
    return rows


def fault_detection_summary(
    rows: Sequence[dict[str, Any]], *, minimum_factor: float
) -> dict[str, Any]:
    """Reduce the fault rows to one verdict per injected variant."""

    verdicts: dict[str, Any] = {}
    for variant in INJECTED_FAULTS:
        entries = [row for row in rows if row["bracket_variant"] == variant]
        factors = [
            float(value)
            for row in entries
            for value in (
                row["full_detection_factor"],
                row["single_mode_detection_factor"],
            )
            if value is not None
        ]
        refusals = [
            row["axial_case"]
            for row in entries
            if row["full_solve_error"] is not None
            or row["single_mode_error"] is not None
        ]
        unchanged = [
            row["axial_case"]
            for row in entries
            if row.get("bracket_is_bitwise_unchanged")
        ]
        # A case whose bracket the fault leaves bitwise unchanged cannot
        # possibly be detected, and pretending otherwise would be dishonest:
        # such cases are excluded from the detection requirement and named.
        detectable = [
            row for row in entries if not row.get("bracket_is_bitwise_unchanged")
        ]
        detected = [row for row in detectable if row["detected"]]
        strong = [
            row
            for row in detectable
            if (
                (
                    row["full_detection_factor"] is not None
                    and row["full_detection_factor"] >= minimum_factor
                )
                or (
                    row["single_mode_detection_factor"] is not None
                    and row["single_mode_detection_factor"] >= minimum_factor
                )
                or row["full_solve_error"] is not None
                or row["single_mode_error"] is not None
            )
        ]
        altered_factors = [
            float(value)
            for row in detectable
            for value in (
                row["full_detection_factor"],
                row["single_mode_detection_factor"],
            )
            if value is not None
        ]
        verdicts[variant] = {
            "cases_with_an_altered_bracket": [row["axial_case"] for row in detectable],
            "cases_whose_bracket_is_bitwise_unchanged": unchanged,
            "cases_detected": [row["axial_case"] for row in detected],
            "cases_detected_strongly": [row["axial_case"] for row in strong],
            "cases_refused_by_the_solver": sorted(set(refusals)),
            "min_detection_factor": min(factors) if factors else None,
            "max_detection_factor": max(factors) if factors else None,
            "min_detection_factor_over_altered_cases": (
                min(altered_factors) if altered_factors else None
            ),
            "all_altered_cases_detected": bool(
                detectable and len(detected) == len(detectable)
            ),
            "all_altered_cases_detected_strongly": bool(
                detectable and len(strong) == len(detectable)
            ),
        }
    return verdicts


# ------------------------------------ W-B condition 6: support-leak refusal


def measure_support_leak(config: dict[str, Any]) -> dict[str, Any]:
    """Two ways of violating compact support; both must be refused."""

    profile = config["source_profile"]
    exponent = int(profile["exponent"])
    declared_support = float(profile["support_radius"])
    settings = config["support_leak"]
    ppur = int(settings["points_per_unit_radius"])
    wall_radius = float(settings["wall_radius"])
    leaking_support = float(settings["leaking_support_radius"])
    axial_points = int(config["axial_points"])

    grid = build_grid(
        points_per_unit_radius=ppur,
        wall_radius=wall_radius,
        axial_points=axial_points,
        z_period=1.0,
    )
    leaking = build_source(
        grid,
        support_radius=leaking_support,
        exponent=exponent,
        z_period=1.0,
        mode=1,
    )
    residue_at_and_beyond_wall = _max_abs(leaking[grid.r >= wall_radius, :])
    residue_beyond_declared = _max_abs(leaking[grid.r >= declared_support, :])

    outcomes: dict[str, Any] = {
        "leaking_support_radius": leaking_support,
        "wall_radius": wall_radius,
        "declared_support_radius": declared_support,
        "max_abs_source_at_and_beyond_wall": residue_at_and_beyond_wall,
        "max_abs_source_beyond_declared_support": residue_beyond_declared,
    }

    # (i) The declared support reaches the wall: rejected on the declaration
    # alone, before any solve.
    try:
        solve_transparent(grid, leaking, support_radius=leaking_support)
    except ValueError as error:
        outcomes["declared_support_reaches_wall_refused"] = True
        outcomes["declared_support_reaches_wall_message"] = str(error)
    else:
        outcomes["declared_support_reaches_wall_refused"] = False
        outcomes["declared_support_reaches_wall_message"] = None

    # (ii) The declaration looks fine but the field itself leaks: rejected by
    # measuring the field, which is the case a warning would let through.
    try:
        solve_transparent(grid, leaking, support_radius=declared_support)
    except ValueError as error:
        outcomes["undeclared_residue_refused"] = True
        outcomes["undeclared_residue_message"] = str(error)
    else:
        outcomes["undeclared_residue_refused"] = False
        outcomes["undeclared_residue_message"] = None

    # (iii) The honest source on the same grid must still be accepted, so the
    # refusals above are not a solver that refuses everything.
    honest = build_source(
        grid,
        support_radius=declared_support,
        exponent=exponent,
        z_period=1.0,
        mode=1,
    )
    try:
        solve_transparent(grid, honest, support_radius=declared_support)
    except Exception as error:  # noqa: BLE001
        outcomes["compliant_source_accepted"] = False
        outcomes["compliant_source_message"] = f"{type(error).__name__}: {error}"
    else:
        outcomes["compliant_source_accepted"] = True
        outcomes["compliant_source_message"] = None

    outcomes["all_refusals_fired"] = bool(
        outcomes["declared_support_reaches_wall_refused"]
        and outcomes["undeclared_residue_refused"]
        and outcomes["compliant_source_accepted"]
    )
    return outcomes


# --------------------------------------------- oracle and bracket limits


def bessel_oracle_selfcheck() -> dict[str, Any]:
    """Check the ``K_0/K_1`` oracle inside the evidence bundle."""

    published_k0_at_one = 0.4210244382407083
    published_k1_at_one = 0.6019072301972346
    published_ratio = published_k0_at_one / published_k1_at_one
    measured_ratio = bessel_reference.k0_over_k1(1.0)

    arguments = [1.0e-3, 1.0e-2, 0.1, 0.5, 1.0, 2.0, 2.0 * math.pi, 20.0, 60.0]
    wronskian = max(
        bessel_reference.wronskian_relative_defect(value) for value in arguments
    )

    # K_1'(x) = -K_0(x) - K_1(x)/x, by a fourth-order central difference.  This
    # is the identity (W-1) is derived from.
    derivative_defects: list[float] = []
    for value in (0.25, 0.5, 1.0, 2.0, 4.0):
        step = 1.0e-3 * value
        derivative = (
            bessel_reference.bessel_k1(value - 2.0 * step)
            - 8.0 * bessel_reference.bessel_k1(value - step)
            + 8.0 * bessel_reference.bessel_k1(value + step)
            - bessel_reference.bessel_k1(value + 2.0 * step)
        ) / (12.0 * step)
        predicted = -bessel_reference.bessel_k0(value) - bessel_reference.bessel_k1(
            value
        ) / value
        derivative_defects.append(abs(derivative - predicted) / abs(predicted))

    small = [1.0e-6, 1.0e-5, 1.0e-4]
    small_errors = [
        abs(
            bessel_reference.k0_over_k1(value)
            - bessel_reference.k0_over_k1_small_argument_asymptote(value)
        )
        / bessel_reference.k0_over_k1(value)
        for value in small
    ]
    large = [20.0, 40.0, 60.0]
    large_errors = [
        abs(
            bessel_reference.k0_over_k1(value)
            - bessel_reference.k0_over_k1_large_argument_asymptote(value)
        )
        / bessel_reference.k0_over_k1(value)
        for value in large
    ]

    ladder = [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 30.0, 60.0]
    values = [bessel_reference.k0_over_k1(value) for value in ladder]

    return {
        "published_ratio_at_one": published_ratio,
        "measured_ratio_at_one": measured_ratio,
        "published_relative_error": abs(measured_ratio - published_ratio)
        / published_ratio,
        "max_wronskian_relative_defect": wronskian,
        "max_derivative_identity_relative_defect": max(derivative_defects),
        "small_argument_arguments": small,
        "max_small_argument_relative_error": max(small_errors),
        "large_argument_arguments": large,
        "max_large_argument_relative_error": max(large_errors),
        "monotone_increasing": bool(
            all(values[index] < values[index + 1] for index in range(len(values) - 1))
        ),
        "inside_unit_interval": bool(all(0.0 < value < 1.0 for value in values)),
    }


def bracket_limit_check(config: dict[str, Any]) -> dict[str, Any]:
    """Verify the two limits (W-1) is required to have."""

    radius = float(config["wall_radii"][0])
    zero_bracket = outer_bracket(0.0, radius)
    closed_form = 2.0 / radius
    # The k=0 decaying solution C/r^2 must satisfy the condition identically.
    amplitude = 0.7
    derivative = -2.0 * amplitude / radius**3
    residual = derivative + closed_form * (amplitude / radius**2)
    scale = abs(derivative)

    small = []
    for k in (1.0e-3, 1.0e-2, 1.0e-1):
        bracket = outer_bracket(k, radius)
        small.append(
            {
                "wavenumber": k,
                "kr": k * radius,
                "bracket": bracket,
                "curvature_only": closed_form,
                "excess_over_curvature": bracket - closed_form,
                "excess_relative_to_k": (bracket - closed_form) / k,
            }
        )
    large = []
    for k in (10.0, 50.0, 200.0):
        bracket = outer_bracket(k, radius)
        asymptote = closed_form + k
        large.append(
            {
                "wavenumber": k,
                "kr": k * radius,
                "bracket": bracket,
                "large_kr_asymptote": asymptote,
                "relative_gap": abs(bracket - asymptote) / asymptote,
                # K_0/K_1 = 1 - 1/(2x) + O(x^-2) with x = kR, so the bracket
                # falls short of 2/R + k by 1/(2R) to leading order.
                "predicted_relative_gap": (0.5 / radius) / asymptote,
            }
        )
    return {
        "radius": radius,
        "zero_mode_bracket": zero_bracket,
        "zero_mode_closed_form": closed_form,
        "zero_mode_bracket_is_bitwise_two_over_r": bool(zero_bracket == closed_form),
        "zero_mode_bracket_relative_defect": abs(zero_bracket - closed_form)
        / closed_form,
        "zero_mode_decaying_solution_residual": residual,
        "zero_mode_decaying_solution_relative_residual": abs(residual) / scale,
        "small_wavenumber_ladder": small,
        "large_wavenumber_ladder": large,
        # As k -> 0 the whole excess k K_0(kR)/K_1(kR) vanishes *faster than k*,
        # which is why beta_0 = 2/R is a limit and not just a special case.
        "small_limit_excess_over_k_vanishes": bool(
            all(
                small[index]["excess_relative_to_k"]
                > small[index - 1]["excess_relative_to_k"]
                for index in range(1, len(small))
            )
            and small[0]["excess_relative_to_k"] < 2.0e-2
        ),
        "large_limit_gap_decreases": bool(
            all(
                large[index]["relative_gap"] < large[index - 1]["relative_gap"]
                for index in range(1, len(large))
            )
        ),
    }


# ----------------------------------------------------------- orchestration


def evaluate(config: dict[str, Any]) -> dict[str, Any]:
    """Run every measurement and return the analysis block."""

    r_independence = measure_r_independence(config)
    manufactured = measure_manufactured(config)
    large_radius = measure_large_radius_agreement(config)
    faults = measure_fault_injection(config)
    minimum_factor = float(config["acceptance"]["min_fault_detection_factor"])
    return {
        "outer_boundary_conditions_compared": [DIRICHLET, TRANSPARENT],
        "r_independence": r_independence,
        "r_independence_convergence": r_independence_convergence(r_independence),
        "manufactured": manufactured,
        "large_radius_agreement": large_radius,
        "large_radius_convergence": large_radius_convergence(large_radius),
        "fault_injection": faults,
        "fault_verdicts": fault_detection_summary(faults, minimum_factor=minimum_factor),
        "support_leak": measure_support_leak(config),
        "bessel_oracle_selfcheck": bessel_oracle_selfcheck(),
        "bracket_limits": bracket_limit_check(config),
        "source_moment": source_moment(
            support_radius=float(config["source_profile"]["support_radius"]),
            exponent=int(config["source_profile"]["exponent"]),
        ),
    }


def acceptance_report(
    analysis: dict[str, Any], config: dict[str, Any]
) -> dict[str, Any]:
    """Evaluate the preregistered W-B acceptance conditions."""

    tolerances = config["acceptance"]
    resolutions = [int(value) for value in config["radial_resolutions"]]
    finest = resolutions[-1]
    long_wavelength_max_kr = float(config["long_wavelength_max_kr"])

    rows = analysis["r_independence"]
    all_finite = True
    for collection in (
        rows,
        analysis["manufactured"],
        analysis["large_radius_agreement"],
    ):
        for row in collection:
            all_finite &= _finite(
                value
                for value in row.values()
                if isinstance(value, (int, float)) and not isinstance(value, bool)
            )

    scaled_r_independence = [
        float(row["transparent_over_dr_squared"]) for row in rows
    ]
    improvement_finest = [
        float(row["improvement_factor"])
        for row in rows
        if row["points_per_unit_radius"] == finest
        and row["improvement_factor"] is not None
    ]
    long_wavelength_finest = [
        float(row["improvement_factor"])
        for row in rows
        if row["points_per_unit_radius"] == finest
        and row["improvement_factor"] is not None
        and row["mode"]
        and float(row["kr_smaller"]) <= long_wavelength_max_kr
    ]
    zero_mode_relative_errors = [
        float(row["zero_mode_closed_form_relative_error"])
        for row in rows
        if row["zero_mode_closed_form_relative_error"] is not None
        and row["points_per_unit_radius"] == finest
    ]
    defect_scaled = [
        float(row["outer_condition_defect_over_dr_squared"]) for row in rows
    ]
    r_independence_orders = [
        float(order)
        for report in analysis["r_independence_convergence"]
        for order in report["orders"]
    ]

    manufactured_orders = [
        float(row["observed_order_all_rows"])
        for row in analysis["manufactured"]
        if row["observed_order_all_rows"] is not None
    ]
    manufactured_boundary_orders = [
        float(row["observed_order_boundary_row"])
        for row in analysis["manufactured"]
        if row["observed_order_boundary_row"] is not None
    ]
    manufactured_scaled = [
        float(row["max_error_over_dr_squared"]) for row in analysis["manufactured"]
    ]

    discretization_budget = float(
        tolerances["large_radius_discretization_constant"]
    )
    split_factors: list[float] = []
    for row in analysis["large_radius_agreement"]:
        # The preregistered split: O(dr^2) + O(R_big^-2).  The second term is
        # bounded from the oracle without ever using the transparent solve, and
        # the first is a single preregistered constant times dr^2.
        predicted = (
            float(row["oracle_wall_term_bound_at_largest_reference"])
            + discretization_budget * float(row["dr"]) ** 2
        )
        split_factors.append(
            float(row["difference_at_largest_reference"]) / predicted
            if predicted > 0.0
            else math.inf
        )
    wall_corrected_scaled = [
        float(row["wall_corrected_residual_over_dr_squared"])
        for row in analysis["large_radius_agreement"]
        if row["wall_corrected_residual_over_dr_squared"] is not None
    ]
    wall_corrected_orders = [
        float(order)
        for report in analysis["large_radius_convergence"]
        for order in report["wall_corrected_orders"]
    ]

    verdicts = analysis["fault_verdicts"]
    selfcheck = analysis["bessel_oracle_selfcheck"]
    limits = analysis["bracket_limits"]
    leak = analysis["support_leak"]

    checks = {
        "all_reported_metrics_finite": bool(all_finite),
        "transparent_r_independence_at_discretization_level": bool(
            scaled_r_independence
            and max(scaled_r_independence)
            <= float(tolerances["max_transparent_r_independence_over_dr_squared"])
        ),
        "transparent_r_independence_orders_in_band": bool(
            r_independence_orders
            and min(r_independence_orders)
            >= float(tolerances["min_transparent_r_independence_observed_order"])
            and max(r_independence_orders)
            <= float(tolerances["max_transparent_r_independence_observed_order"])
        ),
        "dirichlet_zero_mode_offset_matches_e33e_closed_form": bool(
            zero_mode_relative_errors
            and max(zero_mode_relative_errors)
            <= float(tolerances["max_zero_mode_dirichlet_closed_form_relative_error"])
        ),
        "transparent_beats_dirichlet_everywhere": bool(
            improvement_finest
            and min(improvement_finest)
            >= float(tolerances["min_dirichlet_over_transparent_improvement_factor"])
        ),
        "long_wavelength_improvement_is_large": bool(
            long_wavelength_finest
            and min(long_wavelength_finest)
            >= float(tolerances["min_long_wavelength_improvement_factor"])
        ),
        "manufactured_orders_in_band": bool(
            manufactured_orders
            and min(manufactured_orders)
            >= float(tolerances["min_manufactured_observed_order"])
            and max(manufactured_orders)
            <= float(tolerances["max_manufactured_observed_order"])
        ),
        "manufactured_boundary_row_orders_in_band": bool(
            manufactured_boundary_orders
            and min(manufactured_boundary_orders)
            >= float(tolerances["min_manufactured_boundary_row_observed_order"])
            and max(manufactured_boundary_orders)
            <= float(tolerances["max_manufactured_boundary_row_observed_order"])
        ),
        "manufactured_error_scaled_bounded": bool(
            manufactured_scaled
            and max(manufactured_scaled)
            <= float(tolerances["max_manufactured_error_over_dr_squared"])
        ),
        "large_radius_agreement_within_predicted_split": bool(
            split_factors
            and max(split_factors)
            <= float(tolerances["max_large_radius_dirichlet_split_factor"])
        ),
        "zero_mode_wall_corrected_residual_is_second_order": bool(
            wall_corrected_scaled
            and max(wall_corrected_scaled)
            <= float(
                tolerances["max_zero_mode_wall_corrected_residual_over_dr_squared"]
            )
        ),
        "zero_mode_wall_corrected_orders_in_band": bool(
            wall_corrected_orders
            and min(wall_corrected_orders)
            >= float(tolerances["min_zero_mode_wall_corrected_observed_order"])
            and max(wall_corrected_orders)
            <= float(tolerances["max_zero_mode_wall_corrected_observed_order"])
        ),
        "every_injected_fault_is_detected": bool(
            all(
                verdicts[variant]["all_altered_cases_detected"]
                for variant in INJECTED_FAULTS
            )
        ),
        "every_injected_fault_is_detected_strongly": bool(
            all(
                verdicts[variant]["all_altered_cases_detected_strongly"]
                for variant in INJECTED_FAULTS
            )
        ),
        "support_leak_is_refused": bool(leak["all_refusals_fired"]),
        "outer_condition_defect_is_second_order": bool(
            defect_scaled
            and max(defect_scaled)
            <= float(tolerances["max_outer_condition_defect_over_dr_squared"])
        ),
        "k0_over_k1_matches_published_value": bool(
            float(selfcheck["published_relative_error"])
            <= float(tolerances["max_k0_over_k1_published_relative_error"])
        ),
        "k0_over_k1_satisfies_the_derivative_identity": bool(
            float(selfcheck["max_derivative_identity_relative_defect"])
            <= float(
                tolerances["max_k0_over_k1_derivative_identity_relative_defect"]
            )
        ),
        "bessel_wronskian_identity_holds": bool(
            float(selfcheck["max_wronskian_relative_defect"])
            <= float(tolerances["max_wronskian_relative_defect"])
        ),
        "k0_over_k1_asymptotic_branches_agree": bool(
            float(selfcheck["max_small_argument_relative_error"])
            <= float(tolerances["max_k0_over_k1_small_argument_relative_error"])
            and float(selfcheck["max_large_argument_relative_error"])
            <= float(tolerances["max_k0_over_k1_large_argument_relative_error"])
        ),
        "k0_over_k1_is_monotone_inside_the_unit_interval": bool(
            selfcheck["monotone_increasing"] and selfcheck["inside_unit_interval"]
        ),
        "zero_mode_bracket_is_exactly_two_over_r": bool(
            limits["zero_mode_bracket_is_bitwise_two_over_r"]
            and float(limits["zero_mode_bracket_relative_defect"])
            <= float(tolerances["max_zero_mode_bracket_relative_defect"])
        ),
        "zero_mode_decaying_solution_satisfies_the_condition": bool(
            float(limits["zero_mode_decaying_solution_relative_residual"])
            <= float(tolerances["max_zero_mode_bracket_relative_defect"])
        ),
        "small_kr_bracket_tends_to_two_over_r": bool(
            limits["small_limit_excess_over_k_vanishes"]
        ),
        "large_kr_bracket_tends_to_two_over_r_plus_k": bool(
            limits["large_limit_gap_decreases"]
            and float(limits["large_wavenumber_ladder"][-1]["relative_gap"]) < 1.0e-2
        ),
    }
    checks["all_passed"] = all(checks.values())
    return {
        "checks": checks,
        "observed": {
            "max_transparent_r_independence_over_dr_squared": max(
                scaled_r_independence
            )
            if scaled_r_independence
            else None,
            "transparent_r_independence_observed_orders": r_independence_orders,
            "max_zero_mode_dirichlet_closed_form_relative_error": max(
                zero_mode_relative_errors
            )
            if zero_mode_relative_errors
            else None,
            "min_improvement_factor_at_finest": min(improvement_finest)
            if improvement_finest
            else None,
            "max_improvement_factor_at_finest": max(improvement_finest)
            if improvement_finest
            else None,
            "min_long_wavelength_improvement_factor": min(long_wavelength_finest)
            if long_wavelength_finest
            else None,
            "long_wavelength_pair_count": len(long_wavelength_finest),
            "manufactured_observed_orders": manufactured_orders,
            "manufactured_boundary_row_observed_orders": manufactured_boundary_orders,
            "max_manufactured_error_over_dr_squared": max(manufactured_scaled)
            if manufactured_scaled
            else None,
            "max_large_radius_split_factor": max(split_factors)
            if split_factors
            else None,
            "max_zero_mode_wall_corrected_residual_over_dr_squared": max(
                wall_corrected_scaled
            )
            if wall_corrected_scaled
            else None,
            "zero_mode_wall_corrected_observed_orders": wall_corrected_orders,
            "max_outer_condition_defect_over_dr_squared": max(defect_scaled)
            if defect_scaled
            else None,
            "fault_detection_factor_range": {
                variant: [
                    verdicts[variant]["min_detection_factor"],
                    verdicts[variant]["max_detection_factor"],
                ]
                for variant in INJECTED_FAULTS
            },
        },
    }


R_INDEPENDENCE_CSV_FIELDS: tuple[str, ...] = (
    "axial_case",
    "mode",
    "wavenumber",
    "points_per_unit_radius",
    "dr",
    "smaller_radius",
    "larger_radius",
    "kr_smaller",
    "transparent_core_difference",
    "dirichlet_core_difference",
    "transparent_over_dr_squared",
    "improvement_factor",
    "signed_dirichlet_core_difference",
    "signed_transparent_core_difference",
    "zero_mode_closed_form",
    "zero_mode_closed_form_relative_error",
    "outer_condition_defect_max_abs",
)

MANUFACTURED_CSV_FIELDS: tuple[str, ...] = (
    "axial_case",
    "mode",
    "wavenumber",
    "points_per_unit_radius",
    "dr",
    "wall_radius",
    "max_error_all_rows",
    "max_error_boundary_row",
    "max_error_over_dr_squared",
    "observed_order_all_rows",
    "observed_order_boundary_row",
)

FAULT_CSV_FIELDS: tuple[str, ...] = (
    "axial_case",
    "mode",
    "kr_smaller",
    "bracket_variant",
    "bracket_value",
    "exact_bracket_value",
    "full_solve_r_dependence",
    "full_detection_factor",
    "full_solve_error",
    "single_mode_r_dependence",
    "single_mode_detection_factor",
    "single_mode_error",
    "detected",
    "detected_where",
)

LARGE_RADIUS_CSV_FIELDS: tuple[str, ...] = (
    "axial_case",
    "mode",
    "wavenumber",
    "points_per_unit_radius",
    "dr",
    "transparent_wall_radius",
    "largest_reference_radius",
    "difference_at_largest_reference",
    "oracle_wall_term_bound_at_largest_reference",
    "closed_form_wall_term_at_largest_reference",
    "wall_corrected_residual",
    "wall_corrected_residual_over_dr_squared",
    "difference_over_dr_squared",
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
    acceptance = acceptance_report(analysis, config)

    summary = {
        "schema_version": 1,
        "experiment_id": EXPERIMENT_ID,
        "audit_entry": "E-33",
        "design_document": "docs/whole_space_transition.md (W-1, W-A, W-B)",
        "outer_boundary_condition": TRANSPARENT,
        "outer_boundary_condition_baseline": DIRICHLET,
        "physical_wall_warning": config["physical_wall_warning"],
        "statement": (
            "floating-point verification that the exact modal transparent "
            "(Dirichlet-to-Neumann) outer condition W-1 removes the wall "
            "truncation error of the E-25 finite-cylinder elliptic problem "
            "down to its O(dr^2) boundary discretization error; not a "
            "continuum bound and not singularity evidence"
        ),
        "equation": (
            "-(d_rr + 3/r d_r + d_zz) psi1 = omega1 with "
            "d_r psi_hat_k(R) + [2/R + |k| K_0(|k|R)/K_1(|k|R)] psi_hat_k(R) = 0"
        ),
        "description": config["description"],
        "interpretation": config["interpretation"],
        "conditions_tested": {
            "W-B.1": (
                "k=0 exactness: the transparent core solution is R-independent "
                "up to O(dr^2) and the Dirichlet (Q_inf/2)(R'^-2 - R^-2) offset "
                "is gone; both numbers are reported side by side"
            ),
            "W-B.2": (
                "k>0 exactness across several L_z, with the Dirichlet/transparent "
                "improvement factor recorded, especially for kR << 1"
            ),
            "W-B.3": (
                "manufactured second-order convergence including the boundary "
                "row, against a solution that satisfies W-1 identically"
            ),
            "W-B.4": (
                "agreement with a large-radius Dirichlet solve, split into its "
                "O(dr^2) and O(R_big^-2) parts"
            ),
            "W-B.5": (
                "detection of a sign-flipped bracket, of a dropped 2/R term and "
                "of K_0/K_1 frozen to its kR >> 1 limit"
            ),
            "W-B.6": (
                "refusal of a source that is not compactly supported strictly "
                "inside R"
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
    r_independence_path = output_dir / "r_independence.csv"
    manufactured_path = output_dir / "manufactured_convergence.csv"
    fault_path = output_dir / "fault_injection.csv"
    large_radius_path = output_dir / "large_radius_agreement.csv"

    _write_json(summary_path, summary)
    _write_json(config_path, config)
    _write_csv(
        r_independence_path, R_INDEPENDENCE_CSV_FIELDS, analysis["r_independence"]
    )
    _write_csv(manufactured_path, MANUFACTURED_CSV_FIELDS, analysis["manufactured"])
    _write_csv(fault_path, FAULT_CSV_FIELDS, analysis["fault_injection"])
    _write_csv(
        large_radius_path,
        LARGE_RADIUS_CSV_FIELDS,
        analysis["large_radius_agreement"],
    )

    payloads = (
        summary_path,
        config_path,
        r_independence_path,
        manufactured_path,
        fault_path,
        large_radius_path,
    )
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
    for path in (*payloads, manifest_path):
        path.with_name(path.name + ".sha256").write_text(
            _sha256(path) + "\n", encoding="ascii"
        )

    if not summary["accepted"]:
        failed = sorted(
            name for name, value in acceptance["checks"].items() if not value
        )
        raise RuntimeError(
            f"transparent boundary experiment failed acceptance: {failed}"
        )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Verify the W-A transparent outer boundary condition."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=REPOSITORY_ROOT / "configs" / "transparent_boundary.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "transparent_boundary_v1",
    )
    arguments = parser.parse_args(argv)
    config_path = arguments.config.resolve()
    output_dir = arguments.output_dir.resolve()
    if not _inside_repository(config_path) or not _inside_repository(output_dir):
        parser.error("config and output paths must remain inside this repository")
    config = strict_json_loads(
        config_path.read_text(encoding="utf-8"),
        label="transparent boundary config",
    )
    if not isinstance(config, dict):
        parser.error("config must be a JSON object")
    summary = run(config, output_dir)
    print(
        json.dumps(
            summary["acceptance"]["checks"], ensure_ascii=False, sort_keys=True
        )
    )
    return 0 if summary["accepted"] else 2


if __name__ == "__main__":
    sys.exit(main())
