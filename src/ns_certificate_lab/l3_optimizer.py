r"""Search for initial data with a **positive** critical-norm generation rate.

The objective
-------------
:mod:`ns_certificate_lab.l3_generation` shows that

.. math::  J(u_0) = \underbrace{3\int p\,\nabla\!\cdot\!(|u|u)\,dx}_{P}
                    \;\underbrace{-\,3\nu\int|u|\bigl(|\nabla u|^2
                      + |\nabla|u||^2\bigr)dx}_{V\;\le\;0} ,

with the advective contribution exactly zero.  Two exact symmetries fix the
shape of the search.

**``P`` is odd and ``V`` is even under ``u\mapsto -u``.**  ``p`` is quadratic in
``u`` so it is unchanged, while ``\nabla\!\cdot\!(|u|u)`` changes sign; the
viscous integrand is a sum of even terms.  So for *any* datum with ``P\ne 0``,
one of the two signs gives ``P>0``.  Searching over shapes never has to search
over that sign — it is read off.

**Both terms are homogeneous of degree four at fixed Reynolds number.**  Under
``u\mapsto\lambda u`` with ``\mathrm{Re} = AL^2/\nu`` held fixed, ``\nu`` scales
like ``\lambda`` and ``P``, ``V`` both scale like ``\lambda^4``.  The overall
amplitude is therefore *not* a search direction either.

What is left is genuinely the shape, and the natural scale-free objective is the
**critical Reynolds number**

.. math::  \mathrm{Re}_{\rm crit}
    = \frac{A L}{\nu}\cdot\frac{|V|}{P}
    = A L\,\frac{3\int|u|(|\nabla u|^2+|\nabla|u||^2)dx}{P} ,

at which ``J`` changes sign: ``J>0`` exactly when ``\mathrm{Re} >
\mathrm{Re}_{\rm crit}``.  It is invariant under ``u\mapsto\lambda u`` and does
not mention ``\nu`` at all.  **Minimising it is the search.**

One power of ``L``, not two, because ``A`` here is the **physical velocity**
amplitude.  This is the same Reynolds number as the ``\mathrm{Re}=A_{u_1}L^2/\nu``
of :mod:`ns_certificate_lab.nondimensional`, whose amplitude is the swirl
*variable* ``u_1 = u^\theta/r``, since ``A = A_{u_1}L``.

A sharper statement of why the amplitude is not a search direction: the
Navier–Stokes scaling ``u\mapsto\lambda u(\lambda x)`` that preserves
``\|u\|_{L^3}`` sends ``A\mapsto\lambda A`` and ``L\mapsto L/\lambda``, so
``\mathrm{Re}`` is **invariant** and ``P`` and ``V`` both scale by
``\lambda^2``.  ``J>0`` cannot be manufactured by rescaling inside a critical
family; only the shape moves ``\mathrm{Re}_{\rm crit}``.

Why not maximise the vorticity
------------------------------
Because ``\max|\omega|`` is not a critical quantity: it is not invariant under
the scaling ``u_\lambda(x,t)=\lambda u(\lambda x,\lambda^2t)`` that leaves
Navier–Stokes invariant, so making it large is a statement about units rather
than about the solution.  ``\|u\|_{L^3}`` *is* invariant, which is why its
generation rate is the objective here.

The search space
----------------
A finite set of amplitudes over a **fixed** compactly supported generator basis.
Every point of the space is automatically ``C^\infty``, compactly supported,
divergence free, finite energy, finite ``L^3`` and axis regular, because those
properties hold for each basis element and are preserved by linear combination —
so the optimiser can never wander out of the admissible set and no feasibility
projection is needed.

Duplicates from pure rescaling are impossible for the same reason:
``\mathrm{Re}_{\rm crit}`` is scale invariant, so two amplitude vectors that
differ by a positive multiple give the same objective and are the same
candidate.  :func:`normalise_amplitudes` puts every vector on the unit sphere so
they are literally the same point.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
import numpy.typing as npt

from .grid import AxisymmetricGrid
from .l3_generation import (
    AxisymmetricPressureSolver,
    GenerationRate,
    l3_generation_rate,
)
from .mixed_initial_data import GeneratorComponent, MixedFamily

FloatArray = npt.NDArray[np.float64]

__all__ = [
    "CandidateScore",
    "OptimizationResult",
    "KILL_CONDITION_CRITICAL_REYNOLDS",
    "SEARCH_BASIS",
    "critical_reynolds",
    "require_clay_admissible",
    "evaluate_shape",
    "normalise_amplitudes",
    "optimise",
    "search_basis",
]


def require_clay_admissible(family) -> None:
    r"""Refuse to score or certify a family that is not a legitimate candidate.

    The Clay statements require ``C^\infty`` data.  A finite-``C^k`` surrogate
    (a spline basis, for instance) is a permitted *tool* — for checker
    development, HS-5 prototypes, or as an optimisation surrogate whose optimum
    is then mollified — but it is never itself a candidate, and nothing
    downstream of this guard may describe it as one.  A family opts out by
    setting ``clay_admissible = False``; absence of the attribute is treated as
    admissible so the existing ``C^\infty`` families need no change.
    """
    if not getattr(family, "clay_admissible", True):
        raise ValueError(
            f"family {getattr(family, 'name', '?')!r} is a finite-C^k surrogate "
            "(clay_admissible=False); it may not be scored or certified as a "
            "Clay candidate.  Mollify it into a C-infinity datum first and "
            "prove the transfer bound |J(q_eps) - J(q)| explicitly."
        )


def normalise_amplitudes(values: npt.ArrayLike) -> FloatArray:
    """Put an amplitude vector on the unit sphere.

    Scale is not a search direction, so two vectors differing by a positive
    multiple are the same candidate and must be the same point of the space.
    """
    array = np.asarray(values, dtype=np.float64).ravel()
    norm = float(np.linalg.norm(array))
    if norm <= 0.0:
        raise ValueError("the amplitude vector must be nonzero")
    return array / norm


def critical_reynolds(rate: GenerationRate, reference_length: float) -> float:
    r"""``\mathrm{Re}_{\rm crit}``, or ``inf`` when the pressure term is useless.

    ``J>0`` exactly when ``\mathrm{Re}>\mathrm{Re}_{\rm crit}``.  A non-positive
    ``P`` returns ``inf``: no viscosity makes such a shape work, and the caller
    should flip the sign of the datum instead.
    """
    if reference_length <= 0.0:
        raise ValueError("reference_length must be positive")
    if rate.pressure <= 0.0:
        return math.inf
    amplitude = rate.max_speed
    dissipation = abs(rate.viscous) / rate.viscosity
    return amplitude * reference_length * dissipation / rate.pressure


@dataclass(frozen=True)
class CandidateScore:
    """One evaluated shape."""

    amplitudes: tuple[float, ...]
    sign: float
    critical_reynolds: float
    pressure: float
    viscous_per_viscosity: float
    max_speed: float
    l3_cubed: float
    kinetic_energy: float
    relative_divergence: float
    poisson_residual: float
    transport_residual: float
    pressure_tail_bound: float
    outer_gap: float
    support_radius: float

    def as_dict(self) -> dict[str, object]:
        payload = {
            name: getattr(self, name)
            for name in self.__dataclass_fields__
            if name != "amplitudes"
        }
        payload["amplitudes"] = list(self.amplitudes)
        payload["critical_reynolds"] = (
            float(self.critical_reynolds)
            if math.isfinite(self.critical_reynolds)
            else None
        )
        return payload


def search_basis(
    *,
    swirl_widths: tuple[tuple[float, float, float], ...] = (
        (1.2, 1.5, 0.0), (0.8, 0.9, 0.0), (1.0, 0.8, 0.5), (0.7, 0.7, -0.5),
    ),
    stream_widths: tuple[tuple[float, float, float], ...] = (
        (1.2, 1.5, 0.0), (0.9, 0.8, 0.0), (0.6, 0.5, 0.0), (1.0, 0.9, 0.4),
    ),
    reference_length: float = 1.2,
) -> MixedFamily:
    r"""The fixed compactly supported basis the optimiser works over.

    Each entry is ``(radial support, axial support, axial centre)``.  Swirl
    components are odd in ``z`` about their own centre; stream components are
    odd as well, since the parity rule of
    :mod:`ns_certificate_lab.mixed_initial_data` kills the pressure term for an
    even stream generator on a symmetric domain.  Two entries of each group are
    centred away from ``z = 0``, which breaks the axial symmetry and lets the
    optimiser leave the symmetric subspace instead of being confined to it.
    """
    swirl = tuple(
        GeneratorComponent(
            amplitude=0.0, radial_support=radial, axial_support=axial,
            axial_center=centre, axial_concentration=0.5, odd_axial=True,
        )
        for radial, axial, centre in swirl_widths
    )
    stream = tuple(
        GeneratorComponent(
            amplitude=0.0, radial_support=radial, axial_support=axial,
            axial_center=centre, odd_axial=True,
        )
        for radial, axial, centre in stream_widths
    )
    return MixedFamily(
        name="basis", swirl=swirl, stream=stream,
        reference_length=reference_length,
    )


#: The preregistered basis: four swirl and four stream generators, eight
#: amplitudes in total.
SEARCH_BASIS = search_basis()

#: **Preregistered kill condition, recorded before the search was run.**  A
#: generic band-limited divergence-free field gives a pressure-to-viscous ratio
#: of order ``1e-3``, that is a critical Reynolds number of order ``1e3``.  If
#: shape optimisation over :data:`SEARCH_BASIS` cannot bring the best
#: critical Reynolds number below this value, the shape-factor lane is declared
#: dead: the scaling argument in the module docstring forbids rescuing it by
#: rescaling, and a larger basis would be a different, separately preregistered
#: experiment.
KILL_CONDITION_CRITICAL_REYNOLDS = 1.0e2


def evaluate_shape(
    amplitudes: npt.ArrayLike,
    *,
    basis: MixedFamily = SEARCH_BASIS,
    grid: AxisymmetricGrid,
    solver: AxisymmetricPressureSolver | None = None,
    probe_viscosity: float = 1.0e-3,
) -> CandidateScore:
    """Score one shape, taking whichever overall sign makes ``P`` positive.

    ``probe_viscosity`` only sets the scale of the reported ``V``; the objective
    divides it out, so the score does not depend on it.
    """
    require_clay_admissible(basis)
    values = normalise_amplitudes(amplitudes)
    pressure_solver = (
        AxisymmetricPressureSolver.build(grid) if solver is None else solver
    )
    best: CandidateScore | None = None
    for sign in (1.0, -1.0):
        family = basis.with_amplitudes(values * sign)
        field = family.field(grid)
        if float(np.max(field.speed)) <= 0.0:
            continue
        rate = l3_generation_rate(
            field, viscosity=probe_viscosity, solver=pressure_solver
        )
        score = CandidateScore(
            amplitudes=tuple(float(v) for v in values),
            sign=sign,
            critical_reynolds=critical_reynolds(rate, basis.reference_length),
            pressure=rate.pressure,
            viscous_per_viscosity=abs(rate.viscous) / rate.viscosity,
            max_speed=rate.max_speed,
            l3_cubed=rate.l3_cubed,
            kinetic_energy=rate.kinetic_energy,
            relative_divergence=rate.relative_divergence,
            poisson_residual=rate.poisson_residual,
            transport_residual=rate.transport_residual,
            pressure_tail_bound=rate.pressure_tail_bound,
            outer_gap=rate.outer_gap,
            support_radius=rate.support_radius,
        )
        if best is None or score.critical_reynolds < best.critical_reynolds:
            best = score
    if best is None:
        raise ValueError("the amplitude vector produces an identically zero field")
    return best


def _objective(score: CandidateScore) -> float:
    """``log Re_crit``, or a large finite penalty when ``P<=0``.

    The logarithm matters: ``Re_crit`` ranges over orders of magnitude, and a
    linear objective would spend the whole search flattening the largest value
    instead of improving the best one.
    """
    if not math.isfinite(score.critical_reynolds) or score.critical_reynolds <= 0.0:
        return 60.0
    return math.log(score.critical_reynolds)


@dataclass(frozen=True)
class OptimizationResult:
    """The outcome of a multi-start search."""

    best: CandidateScore
    starts: int
    evaluations: int
    history: tuple[float, ...]
    converged: bool

    def as_dict(self) -> dict[str, object]:
        return {
            "best": self.best.as_dict(),
            "starts": self.starts,
            "evaluations": self.evaluations,
            "history": list(self.history),
            "converged": self.converged,
        }


def optimise(
    *,
    grid: AxisymmetricGrid,
    basis: MixedFamily = SEARCH_BASIS,
    starts: int = 8,
    iterations: int = 40,
    seed: int = 20260729,
    step: float = 0.25,
    gradient_step: float = 1.0e-3,
    shrink: float = 0.6,
    minimum_step: float = 1.0e-3,
    probe_viscosity: float = 1.0e-3,
) -> OptimizationResult:
    r"""Multi-start projected-gradient descent on ``\log \mathrm{Re}_{\rm crit}``.

    The gradient is a central difference in each amplitude — a *direct* gradient
    rather than an adjoint.  With eight design variables that is sixteen extra
    evaluations per step, which is cheaper than deriving and validating an
    adjoint through a non-smooth objective: ``|u|`` is not differentiable on the
    zero set of ``u``, so an adjoint would need the same regularisation the
    identity does, and the saving would not pay for the extra failure mode.

    After each step the amplitude vector is projected back onto the unit sphere,
    which is exact rather than a penalty: scale is not a search direction.

    The starts are drawn from a symmetric Gaussian, so nothing biases the search
    toward the ``z``-symmetric subspace; the basis itself contains off-centre
    components so that subspace is not invariant either.
    """
    rng = np.random.default_rng(seed)
    solver = AxisymmetricPressureSolver.build(grid)
    dimension = len(basis.swirl) + len(basis.stream)
    evaluations = 0
    history: list[float] = []
    best: CandidateScore | None = None
    converged = False

    def score(vector: FloatArray) -> CandidateScore:
        nonlocal evaluations
        evaluations += 1
        return evaluate_shape(
            vector, basis=basis, grid=grid, solver=solver,
            probe_viscosity=probe_viscosity,
        )

    for start in range(starts):
        if start == 0:
            # One deterministic start from the shape Gate 8 measured by hand, so
            # a regression in the optimiser cannot hide behind randomness.
            current = np.zeros(dimension)
            current[0] = 1.0
            current[len(basis.swirl)] = 0.4
            current = normalise_amplitudes(current)
        else:
            current = normalise_amplitudes(rng.normal(size=dimension))
        current_score = score(current)
        current_value = _objective(current_score)
        local_step = step
        for _ in range(iterations):
            gradient = np.zeros(dimension)
            for index in range(dimension):
                probe = current.copy()
                probe[index] += gradient_step
                plus = _objective(score(normalise_amplitudes(probe)))
                probe = current.copy()
                probe[index] -= gradient_step
                minus = _objective(score(normalise_amplitudes(probe)))
                gradient[index] = (plus - minus) / (2.0 * gradient_step)
            magnitude = float(np.linalg.norm(gradient))
            if magnitude <= 0.0:
                break
            direction = gradient / magnitude
            trial = normalise_amplitudes(current - local_step * direction)
            trial_score = score(trial)
            trial_value = _objective(trial_score)
            if trial_value < current_value:
                current, current_score, current_value = trial, trial_score, trial_value
            else:
                local_step *= shrink
                if local_step < minimum_step:
                    converged = True
                    break
        history.append(current_value)
        if best is None or current_score.critical_reynolds < best.critical_reynolds:
            best = current_score

    assert best is not None
    return OptimizationResult(
        best=best, starts=starts, evaluations=evaluations,
        history=tuple(history), converged=converged,
    )
