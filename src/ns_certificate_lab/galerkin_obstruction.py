"""Track F obstruction: no finite-mode ansatz can be a forced Clay counterexample.

Context
-------
`START_NEW_SESSION_NAVIER_STOKES.md` §5 "Track F" proposes to reverse-engineer a
smooth force ``f`` for Clay statements (C)/(D): pick a velocity/pressure pair
``(u,p)`` that becomes singular at a finite time ``T``, define

    f := ∂_t u + (u·∇)u - ν Δu + ∇p,

and require ``f`` to stay smooth across ``T``.  §6 "優先候補A" then suggests
searching low-order divergence-free Fourier ansätze for such cancellations, and
its step 6 asks for a *no-go proof* if the search space is empty.

This module implements the machine-checkable half of that no-go proof.  The
mathematical statement and its full proof live in
``docs/research_notes/track_f_finite_mode_nogo.md``; the short version is:

    Let ``S ⊂ ℤ³`` be a finite symmetric mode set and let ``u(t)`` be a real
    divergence-free trigonometric polynomial with modes in ``S`` for every
    ``t ∈ [0,T)``.  If the Track-F residual satisfies ``∫₀ᵀ ‖f(t)‖_{L²} dt < ∞``
    then ``‖u(t)‖_{L²} ≤ ‖u(0)‖_{L²} + ∫₀ᵗ ‖f‖_{L²}`` for every ``t < T``,
    hence — all norms being equivalent on the finite-dimensional space ``V_S`` —
    every derivative of ``u`` stays bounded and ``u`` extends smoothly past
    ``T``.  No finite-mode ansatz is singular, so none can serve as the
    breakdown solution of Clay (C)/(D).

The single algebraic fact carrying the argument is the *energy neutrality* of
the Navier--Stokes nonlinearity restricted to ``V_S``,

    ⟨u, (u·∇)u⟩_{L²(𝕋³)} = 0                                            (T-1)

for every divergence-free trigonometric polynomial ``u``.  On paper (T-1) is one
integration by parts.  Here it is verified *exactly*, with integer arithmetic
and no floating point at all: the trilinear form is expanded into its complete
monomial list over ``ℤ[i]`` in coordinates adapted to the divergence-free
constraint, and every coefficient is asserted to vanish.  A fault-injection
entry point (``allow_longitudinal``) deliberately breaks the constraint so a
test can confirm the checker fires.

Scope caveats (see AGENTS.md).
* This module proves nothing about ``ℝ³`` solutions with infinitely many active
  modes, i.e. nothing about the Millennium problem itself.  It removes one
  explicitly delimited search class.
* The floating-point Galerkin integrator below is a *cross-check* of the exact
  certificate, not evidence for it.  The exact certificate is the theorem-side
  object; the integrator only demonstrates that the proved a priori bound is
  respected step by step and that injected faults break it.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from math import gcd, pi, sqrt
from typing import Callable, Iterable, Mapping, Sequence

import numpy as np

__all__ = [
    "BoundStreamReport",
    "DivergenceFreeModeSet",
    "GalerkinSystem",
    "NoGoCertificate",
    "TrilinearCertificate",
    "apriori_norm_bound",
    "build_galerkin_system",
    "build_mode_set",
    "derivative_amplification",
    "finite_mode_no_go_certificate",
    "sobolev_amplification",
    "stream_apriori_bound",
    "suggested_time_step",
    "transverse_integer_basis",
    "verify_trilinear_cancellation",
]

Mode = tuple[int, int, int]
#: Exact Gaussian integer ``re + i·im``.
Gauss = tuple[int, int]
#: Homogeneous linear form: variable index -> exact coefficient.
LinearForm = dict[int, Gauss]
#: Polynomial: sorted tuple of variable indices -> exact coefficient.
Polynomial = dict[tuple[int, ...], Gauss]

_ZERO: Gauss = (0, 0)
_UNIT_VECTORS: tuple[Mode, Mode, Mode] = ((1, 0, 0), (0, 1, 0), (0, 0, 1))


# --------------------------------------------------------------------------- #
# Exact Gaussian-integer polynomial arithmetic                                 #
# --------------------------------------------------------------------------- #


def _gauss_add(a: Gauss, b: Gauss) -> Gauss:
    return (a[0] + b[0], a[1] + b[1])


def _gauss_mul(a: Gauss, b: Gauss) -> Gauss:
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def _linear_scale(form: LinearForm, factor: int) -> LinearForm:
    if factor == 0:
        return {}
    return {v: (c[0] * factor, c[1] * factor) for v, c in form.items()}


def _linear_add(left: LinearForm, right: LinearForm) -> LinearForm:
    out = dict(left)
    for variable, coefficient in right.items():
        merged = _gauss_add(out.get(variable, _ZERO), coefficient)
        if merged == _ZERO:
            out.pop(variable, None)
        else:
            out[variable] = merged
    return out


def _linear_product(left: LinearForm, right: LinearForm) -> Polynomial:
    out: Polynomial = {}
    for v1, c1 in left.items():
        for v2, c2 in right.items():
            key = (v1, v2) if v1 <= v2 else (v2, v1)
            merged = _gauss_add(out.get(key, _ZERO), _gauss_mul(c1, c2))
            if merged == _ZERO:
                out.pop(key, None)
            else:
                out[key] = merged
    return out


def _polynomial_times_linear(poly: Polynomial, form: LinearForm) -> Polynomial:
    out: Polynomial = {}
    for monomial, c1 in poly.items():
        for variable, c2 in form.items():
            key = tuple(sorted(monomial + (variable,)))
            merged = _gauss_add(out.get(key, _ZERO), _gauss_mul(c1, c2))
            if merged == _ZERO:
                out.pop(key, None)
            else:
                out[key] = merged
    return out


def _polynomial_add_inplace(target: Polynomial, addend: Polynomial) -> None:
    for monomial, coefficient in addend.items():
        merged = _gauss_add(target.get(monomial, _ZERO), coefficient)
        if merged == _ZERO:
            target.pop(monomial, None)
        else:
            target[monomial] = merged


# --------------------------------------------------------------------------- #
# Integer geometry of the divergence-free constraint                            #
# --------------------------------------------------------------------------- #


def _cross(a: Mode, b: Mode) -> Mode:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: Mode, b: Mode) -> int:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _primitive(v: Mode) -> Mode:
    divisor = gcd(gcd(abs(v[0]), abs(v[1])), abs(v[2]))
    if divisor <= 1:
        return v
    return (v[0] // divisor, v[1] // divisor, v[2] // divisor)


def transverse_integer_basis(k: Mode) -> tuple[Mode, ...]:
    """Integer vectors spanning ``{v : k·v = 0}`` over ``ℚ``.

    For ``k = 0`` the constraint is empty and the three coordinate vectors are
    returned.  Otherwise two of the three cross products ``k × e_i`` are always
    independent, because ``v ↦ k × v`` has image exactly ``k^⊥``.  Working with
    an *integer* (not orthonormal) spanning set is what keeps
    :func:`verify_trilinear_cancellation` inside exact arithmetic.
    """
    if k == (0, 0, 0):
        return _UNIT_VECTORS
    candidates = [_primitive(_cross(k, e)) for e in _UNIT_VECTORS]
    for i in range(3):
        for j in range(i + 1, 3):
            if _cross(candidates[i], candidates[j]) != (0, 0, 0):
                return (candidates[i], candidates[j])
    raise ValueError(f"failed to build a transverse basis for mode {k!r}")


def _positive_representative(k: Mode) -> bool:
    """True when ``k`` is the chosen representative of the pair ``{k, -k}``."""
    for component in k:
        if component > 0:
            return True
        if component < 0:
            return False
    return True  # k == 0 is its own representative


@dataclass(frozen=True)
class DivergenceFreeModeSet:
    """A finite symmetric mode set ``S = -S ⊂ ℤ³`` with its transverse frames.

    ``modes`` is the full symmetric set; ``representatives`` is the chosen half
    (containing the zero mode when present).  ``transverse[k]`` holds integer
    vectors spanning ``k^⊥``; the same frame is used for ``k`` and ``-k`` so
    that the reality constraint ``a_{-k} = conj(a_k)`` is coordinate-wise
    conjugation.
    """

    modes: tuple[Mode, ...]
    representatives: tuple[Mode, ...]
    transverse: Mapping[Mode, tuple[Mode, ...]]

    @property
    def dimension(self) -> int:
        """Real dimension of ``V_S``."""
        total = 0
        for k in self.representatives:
            width = len(self.transverse[k])
            total += width if k == (0, 0, 0) else 2 * width
        return total

    @property
    def max_wavenumber(self) -> float:
        """``R_S = max_{k∈S} |k|`` (Euclidean, in lattice units)."""
        return max(sqrt(float(_dot(k, k))) for k in self.modes)

    def transverse_defect(self) -> int:
        """Exact count of violated constraints ``k · t = 0``.  Zero is required."""
        defect = 0
        for k in self.representatives:
            for t in self.transverse[k]:
                if _dot(k, t) != 0:
                    defect += 1
        return defect


def build_mode_set(
    seeds: Iterable[Sequence[int]],
    *,
    allow_longitudinal: bool = False,
) -> DivergenceFreeModeSet:
    """Build the symmetric closure of ``seeds`` with transverse frames.

    ``allow_longitudinal`` is a **fault-injection switch**: it appends the mode
    ``k`` itself to its own transverse frame, so the coordinates no longer
    parametrize divergence-free fields.  The exact trilinear certificate must
    then fail; ``tests/test_galerkin_obstruction.py`` asserts that it does.
    """
    collected: set[Mode] = set()
    for seed in seeds:
        k = tuple(int(component) for component in seed)
        if len(k) != 3:
            raise ValueError(f"mode {seed!r} is not three-dimensional")
        collected.add(k)  # type: ignore[arg-type]
        collected.add((-k[0], -k[1], -k[2]))
    if not collected:
        raise ValueError("at least one mode is required")

    modes = tuple(sorted(collected))
    representatives = tuple(k for k in modes if _positive_representative(k))
    transverse: dict[Mode, tuple[Mode, ...]] = {}
    for k in representatives:
        frame = transverse_integer_basis(k)
        if allow_longitudinal and k != (0, 0, 0):
            frame = frame + (k,)
        transverse[k] = frame
    return DivergenceFreeModeSet(
        modes=modes,
        representatives=representatives,
        transverse=transverse,
    )


# --------------------------------------------------------------------------- #
# (T-1): exact verification of the trilinear cancellation                       #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class TrilinearCertificate:
    """Result of the exact expansion of ``⟨u,(u·∇)u⟩`` over ``ℤ[i]``.

    The physical trilinear form equals ``2πi`` times the expanded polynomial,
    so the vanishing question is unaffected by the (irrational) prefactor and
    the whole check stays in exact integer arithmetic.
    """

    mode_count: int
    variable_count: int
    resonant_triples: int
    monomials_accumulated: int
    surviving_monomials: int
    max_abs_surviving_coefficient: int
    transverse_defect: int

    @property
    def cancels(self) -> bool:
        return self.surviving_monomials == 0 and self.transverse_defect == 0

    def as_dict(self) -> dict[str, object]:
        return {
            "mode_count": self.mode_count,
            "variable_count": self.variable_count,
            "resonant_triples": self.resonant_triples,
            "monomials_accumulated": self.monomials_accumulated,
            "surviving_monomials": self.surviving_monomials,
            "max_abs_surviving_coefficient": self.max_abs_surviving_coefficient,
            "transverse_defect": self.transverse_defect,
            "cancels": self.cancels,
        }


def _fourier_component_forms(
    mode_set: DivergenceFreeModeSet,
) -> tuple[dict[Mode, tuple[LinearForm, LinearForm, LinearForm]], int]:
    """Express every ``(a_k)_i`` as an exact linear form in the real unknowns.

    For a representative ``k ≠ 0`` the unknowns are ``x_{k,m}, y_{k,m}`` and
    ``a_k = Σ_m (x_{k,m} + i y_{k,m}) t_{k,m}``, ``a_{-k} = conj(a_k)``.  For
    ``k = 0`` reality forces ``a_0`` real, so only ``x_{0,m}`` appear.
    """
    forms: dict[Mode, tuple[LinearForm, LinearForm, LinearForm]] = {}
    next_variable = 0
    for k in mode_set.representatives:
        frame = mode_set.transverse[k]
        positive: list[LinearForm] = [{}, {}, {}]
        negative: list[LinearForm] = [{}, {}, {}]
        for t in frame:
            real_variable = next_variable
            next_variable += 1
            imaginary_variable: int | None = None
            if k != (0, 0, 0):
                imaginary_variable = next_variable
                next_variable += 1
            for i in range(3):
                if t[i] == 0:
                    continue
                positive[i] = _linear_add(positive[i], {real_variable: (t[i], 0)})
                negative[i] = _linear_add(negative[i], {real_variable: (t[i], 0)})
                if imaginary_variable is not None:
                    positive[i] = _linear_add(
                        positive[i], {imaginary_variable: (0, t[i])}
                    )
                    negative[i] = _linear_add(
                        negative[i], {imaginary_variable: (0, -t[i])}
                    )
        forms[k] = (positive[0], positive[1], positive[2])
        if k != (0, 0, 0):
            forms[(-k[0], -k[1], -k[2])] = (negative[0], negative[1], negative[2])
    return forms, next_variable


def verify_trilinear_cancellation(
    mode_set: DivergenceFreeModeSet,
) -> TrilinearCertificate:
    """Expand ``⟨u,(u·∇)u⟩_{L²(𝕋³)}`` exactly and count surviving monomials.

    With ``u = Σ_{k∈S} a_k e^{2πik·x}`` the Fourier orthogonality relations give

        ⟨u,(u·∇)u⟩ = 2πi · Σ_{k+l+m=0, k,l,m ∈ S} (a_l·m)(a_k·a_m),

    a cubic form in the real unknowns with Gaussian-integer coefficients once
    ``a`` is written in the integer transverse frames.  The theorem
    (``docs/research_notes/track_f_finite_mode_nogo.md``, Lemma 1) says every
    coefficient vanishes; this routine checks that claim monomial by monomial
    with no floating point involved.
    """
    forms, variable_count = _fourier_component_forms(mode_set)
    mode_index = set(mode_set.modes)

    total: Polynomial = {}
    resonant = 0
    accumulated = 0
    for k in mode_set.modes:
        for l in mode_set.modes:
            m = (-(k[0] + l[0]), -(k[1] + l[1]), -(k[2] + l[2]))
            if m not in mode_index:
                continue
            resonant += 1
            # (a_l · m): linear form.
            transport: LinearForm = {}
            for i in range(3):
                if m[i] != 0:
                    transport = _linear_add(transport, _linear_scale(forms[l][i], m[i]))
            if not transport:
                continue
            # (a_k · a_m): quadratic form.
            pairing: Polynomial = {}
            for i in range(3):
                _polynomial_add_inplace(
                    pairing, _linear_product(forms[k][i], forms[m][i])
                )
            if not pairing:
                continue
            term = _polynomial_times_linear(pairing, transport)
            accumulated += len(term)
            _polynomial_add_inplace(total, term)

    largest = 0
    for coefficient in total.values():
        largest = max(largest, abs(coefficient[0]), abs(coefficient[1]))
    return TrilinearCertificate(
        mode_count=len(mode_set.modes),
        variable_count=variable_count,
        resonant_triples=resonant,
        monomials_accumulated=accumulated,
        surviving_monomials=len(total),
        max_abs_surviving_coefficient=largest,
        transverse_defect=mode_set.transverse_defect(),
    )


# --------------------------------------------------------------------------- #
# The a priori bound and the norm-equivalence constants                         #
# --------------------------------------------------------------------------- #


def apriori_norm_bound(
    initial_l2_norm: float,
    force_l1_l2: float,
) -> float:
    """``‖u(t)‖_{L²} ≤ ‖u(0)‖_{L²} + ∫₀ᵗ ‖f(s)‖_{L²} ds`` (Theorem 1(i)).

    Neither the viscosity nor the mode set enters: dissipation only helps and
    the nonlinearity is energy-neutral by (T-1).
    """
    if initial_l2_norm < 0.0:
        raise ValueError("initial_l2_norm must be nonnegative")
    if force_l1_l2 < 0.0:
        raise ValueError("force_l1_l2 must be nonnegative")
    return initial_l2_norm + force_l1_l2


def sobolev_amplification(mode_set: DivergenceFreeModeSet, order: float) -> float:
    """Sharp constant in ``‖u‖_{H^s} ≤ C ‖u‖_{L²}`` on ``V_S``.

    ``C = max_{k∈S} (1 + 4π²|k|²)^{s/2}`` for ``s ≥ 0``.
    """
    if order < 0.0:
        raise ValueError("order must be nonnegative")
    return max(
        (1.0 + 4.0 * pi**2 * float(_dot(k, k))) ** (0.5 * order) for k in mode_set.modes
    )


def derivative_amplification(mode_set: DivergenceFreeModeSet, order: int) -> float:
    """Constant in ``‖∂^α u‖_{L^∞} ≤ C ‖u‖_{L²}`` for ``|α| = order`` on ``V_S``.

    ``|∂^α u(x)| ≤ Σ_k (2π)^{|α|}|k^α||a_k| ≤ (2πR_S)^{|α|} √|S| ‖u‖_{L²}`` by
    ``|k^α| ≤ |k|^{|α|}`` and Cauchy--Schwarz over the ``|S|`` retained modes.
    """
    if order < 0:
        raise ValueError("order must be nonnegative")
    return (2.0 * pi * mode_set.max_wavenumber) ** order * sqrt(len(mode_set.modes))


# --------------------------------------------------------------------------- #
# Floating-point Galerkin system (cross-check only, never the proof)            #
# --------------------------------------------------------------------------- #


def _orthonormal_frame(frame: Sequence[Mode]) -> np.ndarray:
    """Gram--Schmidt an integer frame into orthonormal rows (float64)."""
    rows: list[np.ndarray] = []
    for t in frame:
        v = np.asarray(t, dtype=np.float64)
        for previous in rows:
            v = v - float(previous @ v) * previous
        norm = float(np.linalg.norm(v))
        if norm <= 1.0e-12:
            raise ValueError("transverse frame is degenerate")
        rows.append(v / norm)
    return np.asarray(rows, dtype=np.float64)


@dataclass(frozen=True)
class GalerkinSystem:
    """The exact ``V_S``-projection of the momentum equation, in coordinates.

    Applying the ``L²``-orthogonal projection ``Π`` onto ``V_S`` to

        ∂_t u = f - (u·∇)u + νΔu - ∇p

    kills the pressure (gradients are orthogonal to ``V_S``), leaves ``∂_t u``
    and ``νΔu`` untouched (both stay in ``V_S``), and turns the nonlinearity
    into ``-Π[(u·∇)u]``.  The result is the closed ODE

        c' = Πf(t) + B(c,c) + A c,   A diagonal and negative semidefinite,

    which is *not* an approximation of the ansatz: it is the exact evolution of
    a trigonometric-polynomial ansatz whose residual is ``f``.  Modes of the
    residual outside ``S`` are carried entirely by ``f`` and never enter ``c``.

    The real basis is ``√2 ê_{k,m} cos(2πk·x)`` and ``√2 ê_{k,m} sin(2πk·x)``
    for representatives ``k ≠ 0``, and the constants ``e_m`` for ``k = 0``; it
    is ``L²(𝕋³)``-orthonormal, so ``‖u‖_{L²} = ‖c‖₂`` exactly.
    """

    mode_set: DivergenceFreeModeSet
    viscosity: float
    labels: tuple[tuple[Mode, int, str], ...]
    frames: Mapping[Mode, np.ndarray]
    diagonal: np.ndarray

    @property
    def dimension(self) -> int:
        return len(self.labels)

    def coefficients_to_fourier(self, c: np.ndarray) -> dict[Mode, np.ndarray]:
        """Map real coordinates to the complex Fourier coefficients ``a_k``."""
        c = np.asarray(c, dtype=np.float64)
        if c.shape != (self.dimension,):
            raise ValueError(f"expected shape {(self.dimension,)}, got {c.shape}")
        coefficients: dict[Mode, np.ndarray] = {
            k: np.zeros(3, dtype=np.complex128) for k in self.mode_set.modes
        }
        for index, (k, m, kind) in enumerate(self.labels):
            unit = self.frames[k][m]
            if k == (0, 0, 0):
                coefficients[k] = coefficients[k] + c[index] * unit
                continue
            negated = (-k[0], -k[1], -k[2])
            if kind == "cos":
                contribution = (c[index] / sqrt(2.0)) * unit
            else:
                contribution = (-1j * c[index] / sqrt(2.0)) * unit
            coefficients[k] = coefficients[k] + contribution
            coefficients[negated] = coefficients[negated] + np.conjugate(contribution)
        return coefficients

    def fourier_to_coefficients(self, field: Mapping[Mode, np.ndarray]) -> np.ndarray:
        """Project a real field given by its Fourier coefficients onto ``V_S``.

        Only the components along the transverse frames are read, so the Leray
        projection and the truncation to ``S`` happen simultaneously.
        """
        out = np.zeros(self.dimension, dtype=np.float64)
        for index, (k, m, kind) in enumerate(self.labels):
            unit = self.frames[k][m]
            value = complex(np.asarray(field[k], dtype=np.complex128) @ unit)
            if k == (0, 0, 0):
                out[index] = value.real
            elif kind == "cos":
                out[index] = sqrt(2.0) * value.real
            else:
                out[index] = -sqrt(2.0) * value.imag
        return out

    def advection_fourier(self, c: np.ndarray) -> dict[Mode, np.ndarray]:
        """Fourier coefficients of ``(u·∇)u`` restricted to the modes of ``S``."""
        a = self.coefficients_to_fourier(c)
        out: dict[Mode, np.ndarray] = {
            k: np.zeros(3, dtype=np.complex128) for k in self.mode_set.modes
        }
        index = set(self.mode_set.modes)
        for l in self.mode_set.modes:
            for m in self.mode_set.modes:
                k = (l[0] + m[0], l[1] + m[1], l[2] + m[2])
                if k not in index:
                    continue
                transport = complex(a[l] @ np.asarray(m, dtype=np.float64))
                out[k] = out[k] + (2.0j * pi) * transport * a[m]
        return out

    def nonlinear_term(self, c: np.ndarray) -> np.ndarray:
        """Coordinates of ``-Π[(u·∇)u]``."""
        return -self.fourier_to_coefficients(self.advection_fourier(c))

    def viscous_term(self, c: np.ndarray) -> np.ndarray:
        """Coordinates of ``νΔu`` (diagonal, ``-ν4π²|k|²``)."""
        return self.diagonal * np.asarray(c, dtype=np.float64)

    def energy_production(self, c: np.ndarray) -> float:
        """``⟨u, (u·∇)u⟩_{L²}`` evaluated in floating point.

        The exact certificate says this is identically zero; the value returned
        here is the roundoff shadow of that identity and is used only as a
        fault detector.
        """
        c = np.asarray(c, dtype=np.float64)
        return float(-c @ self.nonlinear_term(c))

    def rhs(
        self,
        time: float,
        c: np.ndarray,
        force: Callable[[float], np.ndarray] | None = None,
    ) -> np.ndarray:
        value = self.nonlinear_term(c) + self.viscous_term(c)
        if force is not None:
            value = value + np.asarray(force(time), dtype=np.float64)
        return value

    def velocity(self, c: np.ndarray, points: np.ndarray) -> np.ndarray:
        """Evaluate ``u`` at Cartesian points (shape ``(n,3)``); real by design."""
        points = np.asarray(points, dtype=np.float64)
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("points must have shape (n, 3)")
        a = self.coefficients_to_fourier(c)
        out = np.zeros((points.shape[0], 3), dtype=np.complex128)
        for k, coefficient in a.items():
            phase = np.exp(2.0j * pi * (points @ np.asarray(k, dtype=np.float64)))
            out = out + phase[:, None] * coefficient[None, :]
        return out


def build_galerkin_system(
    mode_set: DivergenceFreeModeSet,
    viscosity: float,
) -> GalerkinSystem:
    """Assemble the orthonormal-coordinate Galerkin system for ``V_S``."""
    if viscosity < 0.0:
        raise ValueError("viscosity must be nonnegative")
    labels: list[tuple[Mode, int, str]] = []
    frames: dict[Mode, np.ndarray] = {}
    diagonal: list[float] = []
    for k in mode_set.representatives:
        frame = _orthonormal_frame(mode_set.transverse[k])
        frames[k] = frame
        decay = -viscosity * 4.0 * pi**2 * float(_dot(k, k))
        for m in range(frame.shape[0]):
            if k == (0, 0, 0):
                labels.append((k, m, "cos"))
                diagonal.append(decay)
            else:
                labels.append((k, m, "cos"))
                diagonal.append(decay)
                labels.append((k, m, "sin"))
                diagonal.append(decay)
    return GalerkinSystem(
        mode_set=mode_set,
        viscosity=float(viscosity),
        labels=tuple(labels),
        frames=frames,
        diagonal=np.asarray(diagonal, dtype=np.float64),
    )


def suggested_time_step(
    system: GalerkinSystem,
    norm: float,
    *,
    safety: float = 0.2,
) -> float:
    """A resolved explicit step for the Galerkin ODE at state norm ``norm``.

    The fastest linear rate is ``ν4π²R_S²``; the nonlinear rate is bounded by
    ``2πR_S‖u‖_{L^∞} ≤ 2πR_S √|S| ‖c‖`` (see :func:`derivative_amplification`).
    RK4 needs ``dt`` times the total rate to stay well inside its stability
    region, so ``safety`` defaults to ``0.2``.

    This is a *numerical* convenience.  Theorem 1 needs no time step at all;
    an unresolved step corrupts only the cross-check, never the bound.
    """
    if norm < 0.0:
        raise ValueError("norm must be nonnegative")
    if safety <= 0.0:
        raise ValueError("safety must be positive")
    radius = system.mode_set.max_wavenumber
    nonlinear = 2.0 * pi * radius * sqrt(len(system.mode_set.modes)) * norm
    viscous = system.viscosity * 4.0 * pi**2 * radius**2
    rate = nonlinear + viscous
    if rate <= 0.0:
        return float("inf")
    return safety / rate


@dataclass(frozen=True)
class BoundStreamReport:
    """Every-step record of the proved bound during a Galerkin integration."""

    steps: int
    completed_steps: int
    final_time: float
    time_step: float
    initial_norm: float
    final_norm: float
    max_norm: float
    max_bound_excess: float
    max_bound_ratio: float
    max_relative_energy_production: float
    diverged: bool

    @property
    def bound_respected(self) -> bool:
        return not self.diverged and self.max_bound_excess <= 0.0

    def as_dict(self) -> dict[str, object]:
        """JSON-safe view.

        A diverged run has genuinely infinite entries.  The repository forbids
        ``Infinity``/``NaN`` tokens in stored JSON, so those entries are written
        as ``null``; ``diverged`` and ``bound_respected`` carry the verdict, so
        no information is lost.
        """

        def finite_or_none(value: float) -> float | None:
            return value if math.isfinite(value) else None

        return {
            "steps": self.steps,
            "completed_steps": self.completed_steps,
            "final_time": self.final_time,
            "time_step": self.time_step,
            "initial_norm": self.initial_norm,
            "final_norm": finite_or_none(self.final_norm),
            "max_norm": finite_or_none(self.max_norm),
            "max_bound_excess": finite_or_none(self.max_bound_excess),
            "max_bound_ratio": finite_or_none(self.max_bound_ratio),
            "max_relative_energy_production": finite_or_none(
                self.max_relative_energy_production
            ),
            "diverged": self.diverged,
            "bound_respected": self.bound_respected,
        }


def stream_apriori_bound(
    system: GalerkinSystem,
    initial: np.ndarray,
    *,
    final_time: float,
    steps: int,
    force: Callable[[float], np.ndarray] | None = None,
    force_l2_sup: float = 0.0,
    tolerance: float = 1.0e-8,
) -> BoundStreamReport:
    """Integrate the Galerkin ODE (RK4) and monitor the proved bound each step.

    The comparison is against ``‖c(0)‖ + force_l2_sup · t`` — an upper bound for
    ``‖c(0)‖ + ∫₀ᵗ‖f‖`` whenever ``force_l2_sup`` dominates ``‖Πf(s)‖`` — relaxed
    by the *relative* slack ``tolerance`` to absorb the integrator's own
    truncation error (RK4 conserves the quadratic invariant only to
    ``O(Δt⁴)`` per step, so an absolute slack would scale wrongly with the
    amplitude).

    A reported violation has exactly two possible causes and this routine does
    not distinguish them: the structural hypothesis of Theorem 1 was broken
    (fault injection), or the time step failed to resolve the dynamics.  Use
    :func:`suggested_time_step` to rule out the second, and read
    ``max_relative_energy_production`` — which is a per-step diagnostic
    independent of the step size — to confirm the first.
    """
    if steps <= 0:
        raise ValueError("steps must be positive")
    if final_time <= 0.0:
        raise ValueError("final_time must be positive")
    c = np.asarray(initial, dtype=np.float64).copy()
    if c.shape != (system.dimension,):
        raise ValueError(f"expected shape {(system.dimension,)}, got {c.shape}")

    dt = final_time / steps
    initial_norm = float(np.linalg.norm(c))
    max_norm = initial_norm
    max_excess = -float("inf")
    max_ratio = 0.0
    max_production = 0.0

    def evaluate(t: float, state: np.ndarray) -> np.ndarray:
        return system.rhs(t, state, force)

    completed = 0
    diverged = False
    with np.errstate(over="ignore", invalid="ignore"):
        for step in range(steps):
            t = step * dt
            k1 = evaluate(t, c)
            k2 = evaluate(t + 0.5 * dt, c + 0.5 * dt * k1)
            k3 = evaluate(t + 0.5 * dt, c + 0.5 * dt * k2)
            k4 = evaluate(t + dt, c + dt * k3)
            c = c + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            if not np.all(np.isfinite(c)):
                diverged = True
                max_excess = float("inf")
                max_ratio = float("inf")
                break
            completed = step + 1
            time = completed * dt
            norm = float(np.linalg.norm(c))
            max_norm = max(max_norm, norm)
            exact_bound = initial_norm + force_l2_sup * time
            bound = exact_bound + tolerance * max(exact_bound, 1.0)
            max_excess = max(max_excess, norm - bound)
            max_ratio = max(max_ratio, norm / bound)
            scale = max(norm, 1.0e-300)
            max_production = max(
                max_production, abs(system.energy_production(c)) / (scale**3)
            )

    return BoundStreamReport(
        steps=steps,
        completed_steps=completed,
        final_time=final_time,
        time_step=dt,
        initial_norm=initial_norm,
        final_norm=float(np.linalg.norm(c)) if not diverged else float("inf"),
        max_norm=max_norm,
        max_bound_excess=max_excess,
        max_bound_ratio=max_ratio,
        max_relative_energy_production=max_production,
        diverged=diverged,
    )


# --------------------------------------------------------------------------- #
# The packaged Track-F verdict                                                  #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class NoGoCertificate:
    """Verdict for one finite mode set, in the vocabulary of AGENTS.md §Review."""

    label: str
    modes: tuple[Mode, ...]
    dimension: int
    max_wavenumber: float
    trilinear: TrilinearCertificate
    sobolev_amplification_h1: float
    sobolev_amplification_h3: float
    gradient_amplification: float

    @property
    def verdict(self) -> str:
        if not self.trilinear.cancels:
            return "certificate_failed"
        return "rejected_as_clay_cd_candidate"

    def as_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "modes": [list(k) for k in self.modes],
            "dimension": self.dimension,
            "max_wavenumber": self.max_wavenumber,
            "trilinear": self.trilinear.as_dict(),
            "sobolev_amplification_h1": self.sobolev_amplification_h1,
            "sobolev_amplification_h3": self.sobolev_amplification_h3,
            "gradient_amplification": self.gradient_amplification,
            "verdict": self.verdict,
        }


def finite_mode_no_go_certificate(
    label: str,
    seeds: Iterable[Sequence[int]],
    *,
    allow_longitudinal: bool = False,
) -> NoGoCertificate:
    """Full Track-F rejection certificate for the ansatz class ``V_S``.

    A ``rejected_as_clay_cd_candidate`` verdict means: *provably* no velocity
    field of this class can be the breakdown solution of Clay (C)/(D), because
    Theorem 1 of ``docs/research_notes/track_f_finite_mode_nogo.md`` bounds
    every norm of the ansatz in terms of ``‖u(0)‖_{L²}`` and ``∫‖f‖_{L²}``.
    It does **not** mean that a search was run and found nothing.
    """
    mode_set = build_mode_set(seeds, allow_longitudinal=allow_longitudinal)
    trilinear = verify_trilinear_cancellation(mode_set)
    return NoGoCertificate(
        label=label,
        modes=mode_set.modes,
        dimension=mode_set.dimension,
        max_wavenumber=mode_set.max_wavenumber,
        trilinear=trilinear,
        sobolev_amplification_h1=sobolev_amplification(mode_set, 1.0),
        sobolev_amplification_h3=sobolev_amplification(mode_set, 3.0),
        gradient_amplification=derivative_amplification(mode_set, 1),
    )
