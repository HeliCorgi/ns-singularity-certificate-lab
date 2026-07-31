r"""Exact rational Fourier machinery on the periodic torus ``T^3`` (Track P).

Why this lane exists
--------------------
The Clay statements (B)/(D) are posed on the periodic torus, and the two
obstructions that block a rigorous a posteriori argument on ``\mathbb R^3`` both
disappear there.  The spectral gap ``|k| \ge 1`` restores the ``-\nu R`` damping
term of the control inequality (on ``\mathbb R^3`` the Laplacian's spectrum
reaches zero and the term vanishes — see
:doc:`../../docs/research_notes/a_posteriori_frameworks`).  And a Fourier–
Galerkin trajectory is an exact trigonometric polynomial, so its **continuous**
Navier–Stokes residual is itself a finite trigonometric polynomial, computable
exactly in rational arithmetic: the periodic version of the ``HS-5`` gap —
recovering a continuous residual from discrete data — closes *by construction*.
Nothing analogous is available on the whole space, where the approximate
solution lives on a grid; the two lanes must never be conflated.

The distinction that keeps this honest
--------------------------------------
The finite-mode no-go (``FiniteModeNoGo.lean``, Track F) excludes trajectories
that **remain** in a fixed finite band for all time.  A finite-band **initial
datum** does not satisfy that hypothesis: the true solution leaves the band
immediately — the Galerkin tail computed here is generically nonzero — so the
no-go excludes nothing about it.  Track P works precisely in the unexcluded
region: finite-band datum, infinite-band true solution, and the distance between
the true solution and the band-limited trajectory controlled by a certificate.
The Lean side of the distinction is ``TrackPFourier.lean``.

Representation
--------------
A scalar field is stored as a dictionary from a **canonical** wavevector ``k``
(first nonzero component positive) to the pair ``(A_k, B_k)`` of coefficients of
``A_k\cos(k\cdot x) + B_k\sin(k\cdot x)``.  Realness is automatic; the mean mode
``k = 0`` carries a cosine coefficient only.  A vector field is a triple of
scalars.  Every operation — product (by the product-to-sum identities), partial
derivative, Leray projection, Sobolev norm — is exact in ``Fraction``
arithmetic, and the same code runs unchanged on
:class:`~ns_certificate_lab.snapshot_certificate.Interval` coefficients, which
is how the slab enclosure of :mod:`ns_certificate_lab.torus_aposteriori`
evaluates the Galerkin field over a whole time slab.

Because the convolution is exact there is **no aliasing and no dealiasing
error**: the quantity a spectral code would call the dealiasing error appears
here instead as the exactly-computed Galerkin tail
``(I-P_G)P(u\cdot\nabla u)``, which is the differential residual of the
trajectory and the input to the control inequality.

Conventions
-----------
Measure normalised to ``(2\pi)^{-3}\,dx``, so Parseval reads
``\langle f,g\rangle = \sum_k \hat f_k\cdot\overline{\hat g_k}`` and a single
``\cos`` or ``\sin`` mode has ``L^2`` norm squared ``\tfrac12(A^2+B^2)``.
Homogeneous norms ``\|f\|_n^2 = \sum_k |k|^{2n}\tfrac12(A_k^2+B_k^2)`` on
mean-zero fields, monotone in ``n`` because ``|k|\ge 1``.  Energy
``= \tfrac12\|u\|_0^2`` and enstrophy ``= \|u\|_1^2`` in this normalisation;
the unnormalised torus values are ``(2\pi)^3`` times larger.

The initial data are finite trigonometric polynomials and therefore
``C^\infty`` — a fact certified in Lean (``TrackPFourier.lean``), not merely
asserted.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Callable, Iterable, Mapping

from .snapshot_certificate import DEFAULT_PRECISION_BITS, Interval

__all__ = [
    "TORUS_FAMILIES",
    "TrigScalar",
    "TrigVector",
    "adot_squared_upper",
    "advection",
    "euclidean_norm_upper",
    "family_P1",
    "family_P2",
    "family_P3",
    "galerkin_modes",
    "galerkin_rhs",
    "gradient_part",
    "leray",
    "residual_tail",
    "sup_derivative_bound",
    "taylor_coefficients",
]

Wavevector = tuple[int, int, int]


def _canonical(k: Wavevector) -> tuple[Wavevector, int]:
    """The canonical representative and the sine-flip sign.

    ``cos`` is even and ``sin`` is odd, so the term at ``-k`` equals the term at
    ``k`` with the sine coefficient negated.  Storing only the half-lattice with
    first nonzero component positive makes realness structural instead of a
    constraint to be checked.
    """
    if k == (0, 0, 0):
        return k, 1
    for component in k:
        if component > 0:
            return k, 1
        if component < 0:
            return (-k[0], -k[1], -k[2]), -1
    raise AssertionError("unreachable")


def _norm_sq(k: Wavevector) -> int:
    return k[0] * k[0] + k[1] * k[1] + k[2] * k[2]


class _FractionOps:
    """Coefficient ring: exact rationals."""

    zero = Fraction(0)

    @staticmethod
    def add(x, y):
        return x + y

    @staticmethod
    def mul(x, y):
        return x * y

    @staticmethod
    def neg(x):
        return -x

    @staticmethod
    def scale(x, q: Fraction):
        return x * q

    @staticmethod
    def of(q: Fraction):
        return q

    @staticmethod
    def is_zero(x) -> bool:
        return x == 0


class _IntervalOps:
    """Coefficient ring: exact rational intervals with outward rounding."""

    zero = Interval(Fraction(0), Fraction(0))

    def __init__(self, bits: int = DEFAULT_PRECISION_BITS) -> None:
        self.bits = bits

    def add(self, x: Interval, y: Interval) -> Interval:
        return (x + y).round_outward(self.bits)

    def mul(self, x: Interval, y: Interval) -> Interval:
        return (x * y).round_outward(self.bits)

    @staticmethod
    def neg(x: Interval) -> Interval:
        return -x

    def scale(self, x: Interval, q: Fraction) -> Interval:
        return x.scale(q).round_outward(self.bits)

    @staticmethod
    def of(q: Fraction) -> Interval:
        return Interval(q, q)

    @staticmethod
    def is_zero(x: Interval) -> bool:
        return x.lower == 0 and x.upper == 0


FRACTION_OPS = _FractionOps()


@dataclass
class TrigScalar:
    """One real scalar trigonometric polynomial with exact coefficients."""

    terms: dict[Wavevector, list]
    ops: object = FRACTION_OPS

    @classmethod
    def zero(cls, ops=FRACTION_OPS) -> "TrigScalar":
        return cls({}, ops)

    def _accumulate(self, k: Wavevector, cos_part, sin_part) -> None:
        key, sign = _canonical(k)
        if sign < 0:
            sin_part = self.ops.neg(sin_part)
        if key == (0, 0, 0):
            # sin(0) = 0: the sine part of the zero mode does not exist.  The
            # cosine part is the mean.
            sin_part = self.ops.zero
        entry = self.terms.get(key)
        if entry is None:
            self.terms[key] = [cos_part, sin_part]
        else:
            entry[0] = self.ops.add(entry[0], cos_part)
            entry[1] = self.ops.add(entry[1], sin_part)

    def cleaned(self) -> "TrigScalar":
        """Drop exactly-zero terms (meaningful for rational coefficients)."""
        kept = {
            k: pair
            for k, pair in self.terms.items()
            if not (self.ops.is_zero(pair[0]) and self.ops.is_zero(pair[1]))
        }
        return TrigScalar(kept, self.ops)

    def __add__(self, other: "TrigScalar") -> "TrigScalar":
        out = TrigScalar(dict(self.terms), self.ops)
        out.terms = {k: [pair[0], pair[1]] for k, pair in self.terms.items()}
        for k, (cos_part, sin_part) in other.terms.items():
            out._accumulate(k, cos_part, sin_part)
        return out

    def __neg__(self) -> "TrigScalar":
        return TrigScalar(
            {k: [self.ops.neg(a), self.ops.neg(b)] for k, (a, b) in self.terms.items()},
            self.ops,
        )

    def scale(self, q: Fraction) -> "TrigScalar":
        return TrigScalar(
            {k: [self.ops.scale(a, q), self.ops.scale(b, q)]
             for k, (a, b) in self.terms.items()},
            self.ops,
        )

    def derivative(self, direction: int) -> "TrigScalar":
        r"""``\partial_j``: ``(A, B) \mapsto (k_j B, -k_j A)`` at the same mode."""
        out = TrigScalar.zero(self.ops)
        for k, (a, b) in self.terms.items():
            factor = Fraction(k[direction])
            if factor == 0:
                continue
            out._accumulate(
                k, self.ops.scale(b, factor), self.ops.scale(a, -factor)
            )
        return out

    def __mul__(self, other: "TrigScalar") -> "TrigScalar":
        r"""Exact product by the product-to-sum identities.

        For modes ``k`` and ``l`` the product contributes at ``k+l`` and
        ``k-l``:

        .. math::
           A_{k+l} \mathrel{+}= \tfrac12(A_1A_2 - B_1B_2),\quad
           B_{k+l} \mathrel{+}= \tfrac12(A_1B_2 + B_1A_2),\\
           A_{k-l} \mathrel{+}= \tfrac12(A_1A_2 + B_1B_2),\quad
           B_{k-l} \mathrel{+}= \tfrac12(B_1A_2 - A_1B_2).

        Exact, hence alias-free: there is no grid and no folding of high modes
        onto low ones.
        """
        ops = self.ops
        half = Fraction(1, 2)
        out = TrigScalar.zero(ops)
        for k, (a1, b1) in self.terms.items():
            for l, (a2, b2) in other.terms.items():
                a1a2 = ops.mul(a1, a2)
                b1b2 = ops.mul(b1, b2)
                a1b2 = ops.mul(a1, b2)
                b1a2 = ops.mul(b1, a2)
                total = (k[0] + l[0], k[1] + l[1], k[2] + l[2])
                diff = (k[0] - l[0], k[1] - l[1], k[2] - l[2])
                out._accumulate(
                    total,
                    ops.scale(ops.add(a1a2, ops.neg(b1b2)), half),
                    ops.scale(ops.add(a1b2, b1a2), half),
                )
                out._accumulate(
                    diff,
                    ops.scale(ops.add(a1a2, b1b2), half),
                    ops.scale(ops.add(b1a2, ops.neg(a1b2)), half),
                )
        return out

    def evaluate(self, x: tuple[float, float, float]) -> float:
        """Float evaluation, for cross-checks only — never for certificates."""
        total = 0.0
        for k, (a, b) in self.terms.items():
            phase = k[0] * x[0] + k[1] * x[1] + k[2] * x[2]
            total += float(a) * math.cos(phase) + float(b) * math.sin(phase)
        return total


@dataclass
class TrigVector:
    """A real vector field as three scalar trigonometric polynomials."""

    components: tuple[TrigScalar, TrigScalar, TrigScalar]

    @classmethod
    def zero(cls, ops=FRACTION_OPS) -> "TrigVector":
        return cls((TrigScalar.zero(ops), TrigScalar.zero(ops), TrigScalar.zero(ops)))

    @classmethod
    def from_modes(
        cls,
        modes: Iterable[tuple[Wavevector, tuple, tuple]],
        ops=FRACTION_OPS,
    ) -> "TrigVector":
        r"""Build ``\sum_k a_k\cos(k\cdot x) + b_k\sin(k\cdot x)`` from triples.

        Each entry is ``(k, a, b)`` with ``a, b`` rational 3-vectors.  The mean
        mode is refused — every Track-P field is mean zero by construction.
        """
        scalars = [TrigScalar.zero(ops) for _ in range(3)]
        for k, a, b in modes:
            if k == (0, 0, 0):
                raise ValueError("Track P fields are mean-zero; k = 0 is not allowed")
            for i in range(3):
                scalars[i]._accumulate(k, ops.of(Fraction(a[i])), ops.of(Fraction(b[i])))
        return cls(tuple(scalars))

    @property
    def ops(self):
        return self.components[0].ops

    def __add__(self, other: "TrigVector") -> "TrigVector":
        return TrigVector(
            tuple(x + y for x, y in zip(self.components, other.components))
        )

    def __neg__(self) -> "TrigVector":
        return TrigVector(tuple(-x for x in self.components))

    def scale(self, q: Fraction) -> "TrigVector":
        return TrigVector(tuple(x.scale(q) for x in self.components))

    def cleaned(self) -> "TrigVector":
        return TrigVector(tuple(x.cleaned() for x in self.components))

    # -- structure ----------------------------------------------------------- #

    def coefficient_table(self) -> dict[Wavevector, tuple[list, list]]:
        """``k -> (a, b)`` with ``a, b`` 3-vectors of ring elements."""
        ops = self.ops
        keys: set[Wavevector] = set()
        for scalar in self.components:
            keys.update(scalar.terms)
        table = {}
        for k in keys:
            a = [
                self.components[i].terms.get(k, [ops.zero, ops.zero])[0]
                for i in range(3)
            ]
            b = [
                self.components[i].terms.get(k, [ops.zero, ops.zero])[1]
                for i in range(3)
            ]
            table[k] = (a, b)
        return table

    def divergence(self) -> TrigScalar:
        return (
            self.components[0].derivative(0)
            + self.components[1].derivative(1)
            + self.components[2].derivative(2)
        )

    def restrict(self, keep: Callable[[Wavevector], bool]) -> "TrigVector":
        return TrigVector(
            tuple(
                TrigScalar(
                    {k: [a, b] for k, (a, b) in scalar.terms.items() if keep(k)},
                    scalar.ops,
                )
                for scalar in self.components
            )
        )

    # -- norms ---------------------------------------------------------------- #

    def sobolev_sq(self, order: int):
        r"""``\|u\|_{\dot H^n}^2 = \sum_{k\ne 0}|k|^{2n}\tfrac12(|a_k|^2+|b_k|^2)``.

        The mean mode carries weight ``|k|^{2n} = 0`` for ``n \ge 1`` and is
        excluded for ``n = 0`` as well: Track-P fields are mean zero, and for
        interval coefficients the enclosure of the mean is slack around an exact
        zero, which must not leak into a norm bound.
        """
        ops = self.ops
        half = Fraction(1, 2)
        total = ops.zero
        for k, (a, b) in self.coefficient_table().items():
            if k == (0, 0, 0):
                continue
            weight = Fraction(_norm_sq(k)) ** order
            mode_sum = ops.zero
            for i in range(3):
                mode_sum = ops.add(mode_sum, ops.mul(a[i], a[i]))
                mode_sum = ops.add(mode_sum, ops.mul(b[i], b[i]))
            total = ops.add(total, ops.scale(mode_sum, weight * half))
        return total

    def energy(self):
        """``(1/2)||u||_0^2`` in the normalised measure."""
        return self.ops.scale(self.sobolev_sq(0), Fraction(1, 2))

    def enstrophy(self):
        """``||u||_1^2``, which equals ``||curl u||_0^2`` for divergence-free fields."""
        return self.sobolev_sq(1)

    def evaluate(self, x: tuple[float, float, float]) -> tuple[float, float, float]:
        return tuple(scalar.evaluate(x) for scalar in self.components)


# --------------------------------------------------------------------------- #
# operators                                                                    #
# --------------------------------------------------------------------------- #


def advection(u: TrigVector, v: TrigVector) -> TrigVector:
    r"""``(u\cdot\nabla)v``, exactly."""
    out = []
    for j in range(3):
        component = TrigScalar.zero(u.ops)
        for i in range(3):
            component = component + (u.components[i] * v.components[j].derivative(i))
        out.append(component)
    return TrigVector(tuple(out))


def leray(field: TrigVector) -> TrigVector:
    r"""The Leray projection, mode by mode: ``\hat a \mapsto \hat a - \frac{k(k\cdot\hat a)}{|k|^2}``.

    Exact and rational because ``|k|^2`` is an integer.  The mean mode is left
    untouched: for the fields this module produces it is exactly zero (the mean
    of ``(u\cdot\nabla)u`` vanishes for divergence-free ``u``), and the callers
    that need that fact assert it rather than assume it.
    """
    ops = field.ops
    scalars = [TrigScalar.zero(ops) for _ in range(3)]
    for k, (a, b) in field.coefficient_table().items():
        if k == (0, 0, 0):
            for i in range(3):
                scalars[i]._accumulate(k, a[i], b[i])
            continue
        k_sq = _norm_sq(k)
        dot_a = ops.zero
        dot_b = ops.zero
        for i in range(3):
            dot_a = ops.add(dot_a, ops.scale(a[i], Fraction(k[i])))
            dot_b = ops.add(dot_b, ops.scale(b[i], Fraction(k[i])))
        for i in range(3):
            ratio = Fraction(-k[i], k_sq)
            scalars[i]._accumulate(
                k,
                ops.add(a[i], ops.scale(dot_a, ratio)),
                ops.add(b[i], ops.scale(dot_b, ratio)),
            )
    return TrigVector(tuple(scalars))


def gradient_part(field: TrigVector) -> TrigVector:
    """``(I - P)`` applied to a field — the part a pressure would remove.

    Family P3 is designed to make this large; the test suite measures it, which
    is what "pressure-driven" means as a checkable property rather than a label.
    """
    return field + (-leray(field))


def galerkin_modes(cutoff_sq: int) -> list[Wavevector]:
    """Canonical wavevectors with ``1 <= |k|^2 <= cutoff_sq``."""
    bound = int(math.isqrt(cutoff_sq))
    modes = []
    for kx in range(-bound, bound + 1):
        for ky in range(-bound, bound + 1):
            for kz in range(-bound, bound + 1):
                k = (kx, ky, kz)
                if not 1 <= _norm_sq(k) <= cutoff_sq:
                    continue
                key, sign = _canonical(k)
                if sign > 0:
                    modes.append(key)
    return sorted(modes)


def galerkin_rhs(u: TrigVector, *, viscosity: Fraction, cutoff_sq: int) -> TrigVector:
    r"""``F(u) = \nu\Delta u - P_G P((u\cdot\nabla)u)`` — the Galerkin field.

    ``\Delta`` acts as ``-|k|^2`` per mode.  The result is supported in the
    Galerkin band ``1 \le |k|^2 \le`` ``cutoff_sq``.
    """
    ops = u.ops
    nonlinear = leray(advection(u, u)).restrict(
        lambda k: 1 <= _norm_sq(k) <= cutoff_sq
    )
    scalars = []
    for i in range(3):
        scalar = TrigScalar.zero(ops)
        for k, (a, b) in u.components[i].terms.items():
            factor = -viscosity * _norm_sq(k)
            scalar._accumulate(k, ops.scale(a, factor), ops.scale(b, factor))
        scalars.append(scalar)
    viscous = TrigVector(tuple(scalars)).restrict(
        lambda k: 1 <= _norm_sq(k) <= cutoff_sq
    )
    return viscous + (-nonlinear)


def residual_tail(u: TrigVector, *, cutoff_sq: int) -> TrigVector:
    r"""``e = (I - P_G)P((u\cdot\nabla)u)`` — the exact differential residual.

    This is the continuous Navier–Stokes residual of the Galerkin trajectory,
    computed in Fourier coefficients with no spatial discretisation anywhere:
    the periodic closure of ``HS-5``.  It is a finite trigonometric polynomial
    supported on ``cutoff_sq < |k|^2``, plus a mean mode that is exactly zero
    for divergence-free ``u`` (asserted by the tests, dropped by the norms).
    """
    return leray(advection(u, u)).restrict(lambda k: _norm_sq(k) > cutoff_sq)


def taylor_coefficients(
    datum: TrigVector, *, viscosity: Fraction, cutoff_sq: int, order: int
) -> list[TrigVector]:
    r"""Exact rational time-Taylor coefficients of the Galerkin trajectory.

    ``c_0`` is the datum and

    .. math::
       (m+1)\,c_{m+1} = \nu\Delta c_m
         - P_G P\Bigl(\sum_{i+j=m}(c_i\cdot\nabla)c_j\Bigr).

    Every coefficient is a finite trigonometric polynomial with rational
    coefficients — a deliverable in itself, and the seed for any future
    Taylor-with-remainder slab enclosure.
    """
    if order < 0:
        raise ValueError("order must be nonnegative")
    ops = datum.ops
    coefficients = [datum]
    keep = lambda k: 1 <= _norm_sq(k) <= cutoff_sq  # noqa: E731
    for m in range(order):
        convolution = TrigVector.zero(ops)
        for i in range(m + 1):
            convolution = convolution + advection(coefficients[i], coefficients[m - i])
        nonlinear = leray(convolution).restrict(keep)
        scalars = []
        for comp in range(3):
            scalar = TrigScalar.zero(ops)
            for k, (a, b) in coefficients[m].components[comp].terms.items():
                factor = -viscosity * _norm_sq(k)
                scalar._accumulate(k, ops.scale(a, factor), ops.scale(b, factor))
            scalars.append(scalar)
        viscous = TrigVector(tuple(scalars)).restrict(keep)
        step = (viscous + (-nonlinear)).scale(Fraction(1, m + 1))
        coefficients.append(step.cleaned())
    return coefficients


# --------------------------------------------------------------------------- #
# rigorous constants                                                           #
# --------------------------------------------------------------------------- #


def adot_squared_upper(lattice_cut: int = 20) -> Fraction:
    r"""A rigorous rational upper bound on ``\sum_{k\ne 0}|k|^{-4}``.

    Exact sum over ``|k|_\infty \le N`` plus the tail bound: the shell
    ``|k|_\infty = m`` has ``(2m+1)^3-(2m-1)^3 = 24m^2+2 \le 26m^2`` points,
    each with ``|k| \ge |k|_\infty = m``, so the tail is at most
    ``\sum_{m>N} 26 m^{-2} \le 26/N``.  This constant powers the embedding
    ``\|f\|_\infty \le \dot A\,\|f\|_{\dot H^2}`` by Cauchy–Schwarz, which is
    the only Sobolev-type input the whole ``H^4`` control inequality needs.
    """
    if lattice_cut < 1:
        raise ValueError("the lattice cut must be at least one")
    exact = Fraction(0)
    for kx in range(-lattice_cut, lattice_cut + 1):
        for ky in range(-lattice_cut, lattice_cut + 1):
            for kz in range(-lattice_cut, lattice_cut + 1):
                n = kx * kx + ky * ky + kz * kz
                if n:
                    exact += Fraction(1, n * n)
    return exact + Fraction(26, lattice_cut)


def euclidean_norm_upper(values, *, bits: int = DEFAULT_PRECISION_BITS) -> Fraction:
    """A rational upper bound on the Euclidean norm of a ring 3-vector."""
    from .l3_certificate import sqrt_interval

    total = Fraction(0)
    for value in values:
        magnitude = (
            value.magnitude if isinstance(value, Interval) else abs(Fraction(value))
        )
        total += magnitude * magnitude
    return sqrt_interval(Interval(total, total), bits=bits).upper


def sup_derivative_bound(
    field: TrigVector, order: int, *, bits: int = DEFAULT_PRECISION_BITS
) -> Fraction:
    r"""``M_j``: a rational bound on every ``j``-th derivative of the field.

    .. math::
       \max_{|\beta|=j}\ \sup_x |\partial^\beta u(x)|
         \;\le\; \sum_k |k|^j\,(|a_k| + |b_k|) ,

    with Euclidean ``|k|`` (which also dominates the Frobenius norm of the
    tensor ``k\otimes v`` that appears when a gradient is taken), Euclidean
    coefficient norms, and ``|k^\beta| \le |k|^{|\beta|}``.  Exact up to the
    outward-rounded square roots; valid verbatim for interval coefficients,
    which is how the slab bounds are produced.
    """
    from .l3_certificate import sqrt_interval

    total = Fraction(0)
    for k, (a, b) in field.coefficient_table().items():
        if k == (0, 0, 0):
            continue
        k_norm = sqrt_interval(
            Interval(Fraction(_norm_sq(k)), Fraction(_norm_sq(k))), bits=bits
        ).upper
        total += (k_norm**order) * (
            euclidean_norm_upper(a, bits=bits) + euclidean_norm_upper(b, bits=bits)
        )
    return total


# --------------------------------------------------------------------------- #
# the three preregistered families                                             #
# --------------------------------------------------------------------------- #


def _check_family(modes) -> None:
    for k, a, b in modes:
        dot_a = sum(Fraction(k[i]) * Fraction(a[i]) for i in range(3))
        dot_b = sum(Fraction(k[i]) * Fraction(b[i]) for i in range(3))
        if dot_a != 0 or dot_b != 0:
            raise ValueError(f"mode {k} is not divergence free: {dot_a}, {dot_b}")


def family_P1() -> TrigVector:
    r"""**P1** — a helical low-mode triad ``(1,0,0), (0,1,0), (1,1,0)``.

    Each mode's coefficient pair is orthogonal to its wavevector, so the field
    is exactly divergence free; the triad closes under convolution
    (``(1,0,0)+(0,1,0)=(1,1,0)``), so the nonlinearity acts inside and just
    outside the triad from the first instant.  The first two modes are Beltrami
    (helical) polarisations; the third breaks pure Beltrami structure so the
    nonlinear term does not degenerate to a gradient — a fact the tests
    measure rather than assume.
    """
    modes = [
        ((1, 0, 0), (0, Fraction(1, 2), 0), (0, 0, Fraction(1, 2))),
        ((0, 1, 0), (0, 0, Fraction(1, 2)), (Fraction(1, 2), 0, 0)),
        ((1, 1, 0), (0, 0, Fraction(1, 2)), (Fraction(1, 4), Fraction(-1, 4), 0)),
    ]
    _check_family(modes)
    return TrigVector.from_modes(modes)


def family_P2() -> TrigVector:
    r"""**P2** — two resonant triads sharing the mode ``(1,1,0)``.

    The triads ``{(1,0,0),(0,1,0),(1,1,0)}`` and ``{(1,1,0),(0,0,1),(1,1,1)}``
    are connected through the shared mode, which is the smallest structure in
    which energy can pass between triads rather than recirculate inside one.
    """
    modes = [
        ((1, 0, 0), (0, Fraction(1, 2), 0), (0, 0, Fraction(1, 2))),
        ((0, 1, 0), (0, 0, Fraction(1, 2)), (Fraction(1, 2), 0, 0)),
        ((1, 1, 0), (0, 0, Fraction(1, 2)), (Fraction(1, 4), Fraction(-1, 4), 0)),
        ((0, 0, 1), (Fraction(1, 3), 0, 0), (0, Fraction(1, 3), 0)),
        (
            (1, 1, 1),
            (Fraction(1, 6), Fraction(-1, 6), 0),
            (Fraction(1, 6), Fraction(1, 6), Fraction(-1, 3)),
        ),
    ]
    _check_family(modes)
    return TrigVector.from_modes(modes)


def family_P3() -> TrigVector:
    r"""**P3** — symmetry-broken and pressure-driven.

    No mirror or permutation symmetry survives the lopsided rational
    coefficients, and the datum is built so the gradient part
    ``(I-P)(u\cdot\nabla u)`` — the part the pressure removes — is large.
    "Pressure-driven" is checked by the tests as
    ``\|(I-P)(u\cdot\nabla u)\|_0 > \tfrac12\|P(u\cdot\nabla u)\|_0``, not
    asserted.
    """
    modes = [
        ((1, 0, 0), (0, Fraction(2, 3), Fraction(1, 5)), (0, Fraction(1, 7), Fraction(-1, 2))),
        ((0, 1, 0), (Fraction(1, 2), 0, Fraction(3, 7)), (Fraction(-1, 3), 0, Fraction(2, 9))),
        ((1, 1, 0), (Fraction(1, 4), Fraction(-1, 4), Fraction(1, 3)), (Fraction(1, 5), Fraction(-1, 5), Fraction(-1, 6))),
        ((1, -1, 1), (Fraction(1, 3), Fraction(1, 2), Fraction(1, 6)), (Fraction(1, 2), Fraction(1, 3), Fraction(-1, 6))),
    ]
    _check_family(modes)
    return TrigVector.from_modes(modes)


#: The preregistered periodic families, by name.
TORUS_FAMILIES: Mapping[str, Callable[[], TrigVector]] = {
    "P1": family_P1,
    "P2": family_P2,
    "P3": family_P3,
}
