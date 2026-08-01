# Mesoscopic coherent-cloud relay audit — 2026-08-01

**Label: SYMBOLIC CANDIDATE / NUMERICAL CANDIDATE / AUDIT REQUIRED**

This note refines the phase-coded Leray-cloud candidate by separating three
questions that an instantaneous best-response calculation conflates:

1. how a unit-\(L^2\) cloud's quadratic forcing scales;
2. how critical energy normalization changes that forcing; and
3. whether an initially empty child acquires a scale-independent fraction of
   the parent energy in one parabolic time.

The third test is the relay test. Growth of \(G_N\) alone is not sufficient.

## 1. Angled parent clouds and affine phase

Use

\[
p=(1,1,0),\qquad q=(1,0,1),\qquad c=p+q=(2,1,1).
\tag{1.1}
\]

For \(W_N=\lfloor N^\gamma\rfloor\), put \(d_N=W_N-1\) and define positive
parent representatives

\[
\Lambda_N^P=Np+[-d_N,d_N]^3,\qquad
\Lambda_N^Q=Nq+[-d_N,d_N]^3.
\tag{1.2}
\]

Negative representatives are their reflections. At each \(k\), project the
central polarization to the physical divergence-free plane:

\[
h_P(k)={P_ke_3\over|P_ke_3|},\qquad
h_Q(k)={P_ke_2\over|P_ke_2|},\qquad
P_k=I-{k\otimes k\over|k|^2}.
\tag{1.3}
\]

Choose positive-frequency phases

\[
\phi_P(k)=x_0\cdot k-\frac\pi2,\qquad
\phi_Q(k)=x_0\cdot k,
\tag{1.4}
\]

and impose \(\widehat u(-k)=\overline{\widehat u(k)}\). For every pair with
\(k+\ell=m\),

\[
\phi_P(k)+\phi_Q(\ell)=x_0\cdot m-\frac\pi2,
\tag{1.5}
\]

so all contributions to one sum mode have the same scalar phase. This is a
physical translation phase plus the sine/cosine phase of the central triad,
not an independent phase optimization.

The cross coefficient before the output Leray projection is

\[
i\sum_{k+\ell=m}\left[
(\ell\cdot\widehat u_P(k))\widehat u_Q(\ell)
+(k\cdot\widehat u_Q(\ell))\widehat u_P(k)\right].
\tag{1.6}
\]

Because the two carrier directions meet at an \(O(1)\) angle, the derivative
factor can be \(O(N)\). A single narrow cone interacting with itself instead
has the divergence-free depletion
\(\widehat u(k)\cdot\ell=\widehat u(k)\cdot(\ell-k)=O(W_N)\); the two cases
must not be mixed.

## 2. Three normalizations

Let \(\phi_N\) denote the whole two-parent shape with \(\|\phi_N\|_2=1\), and
let

\[
E_N={1\over2}\|u_N\|_2^2={c_E\over N},\qquad
u_N=s_N\phi_N,\qquad s_N^2={2c_E\over N}.
\tag{2.1}
\]

Define

\[
A_N^{\rm unit}
=\|P_C\mathbb P((\phi_N\cdot\nabla)\phi_N)\|_2,
\tag{2.2}
\]

\[
A_N^{\rm crit}
=\|P_C\mathbb P((u_N\cdot\nabla)u_N)\|_2
=s_N^2A_N^{\rm unit},
\tag{2.3}
\]

\[
G_N={A_N^{\rm crit}\over N^2\|u_N\|_2^2}
={A_N^{\rm unit}\over N^2}.
\tag{2.4}
\]

If \(M_N\asymp W_N^3\) and the angled affine cloud saturates coherent
capacity, the predicted exponents are

\[
M_N\asymp N^{3\gamma},\quad
A_N^{\rm unit}\asymp N^{1+3\gamma/2},\quad
A_N^{\rm crit}\asymp N^{3\gamma/2},\quad
G_N\asymp N^{3\gamma/2-1}.
\tag{2.5}
\]

The first \(A_N\) exponent is therefore a unit-shape exponent. Reporting it
for the critically normalized field would be a factor-\(N\) error.

## 3. Empty-child Duhamel identity

Let

\[
f_N=-P_C\mathbb P((u_N\cdot\nabla)u_N),\qquad v(0)=0,
\tag{3.1}
\]

and freeze the parent while evolving the child by heat flow. At
\(t=\tau N^{-2}\),

\[
\widehat v_k(t)=N^{-2}h_{\nu,\tau}(|k|/N)\widehat f_N(k),\qquad
h_{\nu,\tau}(r)={1-e^{-\nu\tau r^2}\over\nu r^2}.
\tag{3.2}
\]

Define the exact forcing-weighted heat factor

\[
H_N={\sum_{k\in C}h_{\nu,\tau}(|k|/N)^2|\widehat f_N(k)|^2
\over\sum_{k\in C}|\widehat f_N(k)|^2}.
\tag{3.3}
\]

Then

\[
\boxed{
D_N:={\|v(\tau N^{-2})\|_2^2\over\|u_N(0)\|_2^2}
={E_{\rm child}\over E_{\rm parent}}
={2c_E\over N}H_NG_N^2.}
\tag{3.4}
\]

Thus the two diagnostics requested in the review are

\[
R_N={\|v\|_2\over\|u_N\|_2}=\sqrt{D_N},\qquad
{E_{\rm child}\over E_{\rm parent}}=D_N.
\tag{3.4a}
\]

They are saved as separate columns, the verifier checks (R_N^2=D_N), and
the Duhamel identity residual is checked against the energy ratio only.
The legacy columns `divergence_relative` and `reality_relative` refer only to
the parent field.  The CSV repeats them under explicit `parent_*` names and
records separate `nonlinear_*` divergence and reality residuals from the
zero-padded full convolution.  These labels must not be interchanged.

There is also a phase-independent upper bound. Define the amplitude-effective
mode count

\[
M_N^{\rm eff}
:={\left(\sum_k|\widehat u_N(k)|\right)^2\over\|u_N\|_2^2}
\le M_N .
\]

If \(|k|\le\kappa N\) on the parent support, Bernstein and Parseval give

\[
\|P_C\mathbb P((u_N\cdot\nabla)u_N)\|_2
\le \|u_N\|_\infty\|\nabla u_N\|_2
\le \kappa N\sqrt{M_N^{\rm eff}}\|u_N\|_2^2 .
\tag{3.5}
\]

Heat contraction over \(t=\tau N^{-2}\) therefore yields the rigorous frozen
parent obstruction

\[
\boxed{
D_N\le
2\kappa^2\tau^2c_E\,{M_N^{\rm eff}\over N^3}
\le2\kappa^2\tau^2c_E\,{M_N\over N^3}.}
\tag{3.6}
\]

This bound uses no random-phase assumption and no assertion that coherent
capacity is attained. Thus a relay requires not just \(M_N\gtrsim N^3\), but
amplitude delocalisation \(M_N^{\rm eff}\gtrsim N^3\). In particular
\(M_N=o(N^3)\) is universally insufficient for a scale-independent
one-parabolic-time empty-child relay. If the optimistic capacity law is
saturated, its predicted value has the matching exponent

\[
D_N\asymp N^{3\gamma-3}.
\tag{3.7}
\]

Thus every \(\gamma<1\), including \(4/5\), is asymptotically rejected as a
critical empty-child relay even when \(G_N\) grows. The boundary family is

\[
W_N=\lfloor\eta N\rfloor,\qquad 0<\eta<1/3,
\tag{3.8}
\]

where \(D_N\) may have exponent zero. Its constant and persistence remain
unresolved.

## 4. Geometry and fitting rules

The finite boxes are treated as tagged convolution channels. We record
parent-parent overlap, old-parent/full-child overlap, core-channel overlap,
full-sumset overlap, actual \(W_N/N\), and raw/effective mode counts.  Every
available production row uses all coefficients of a local linear convolution:
an input side (L=2W_N-1) is zero-padded to
(K=2L-1=4W_N-3).  This is not a periodic FFT convolution, so it has no wrap
aliasing.  Overlapping carrier output boxes are summed on global integer
wave numbers before any norm or tag is evaluated.  Widths one and two, and a
deliberately overlapping output case, are independently equal to sparse
all-pairs convolution coefficient by coefficient.

The child core around \(Nc\) has half-width \(d_N\), as required by the
specification. The reachable convolution has half-width \(2d_N\); its outer
part is recorded as child spill, not child forcing. The same split is applied
to the difference, harmonic, and low channels.

The condition \(N>3d_N\) keeps the full reachable child away from the old
parents in this geometry. Full output sumsets can still overlap one another;
this is recorded separately and tagged ratios are not called orthogonal shell
budgets in that case. An overlapping row remains in the CSV with separate
core-fit and strict relay-fit eligibility flags.

The requested values \(N\le64,\gamma\ge0.70\) need not be asymptotic.  Core
quantities \(A_N,G_N,D_N\) may receive a diagnostic fit when the parent,
child core, and core channels are disjoint, even if outer full sumsets
overlap; such a fit is explicitly not relay-admissible.  Off-chain and
self/cross fits require disjoint full outputs (or an actual exact full
complement).  A fit with fewer than four points under its own eligibility
rule is reported as insufficient geometry, not as a measured exponent.

## 5. Channel diagnostics

For each valid row the exact local-FFT calculation reports the core and spill pieces
of child, difference, same-carrier harmonic, and low forcing. Derived ratios
include

\[
r_{\rm diff}={\|F_{\rm difference}\|_2\over\|F_{\rm child}\|_2},\qquad
r_{\rm low}={\|F_{\rm low}\|_2\over\|F_{\rm child}\|_2},
\tag{5.1}
\]

\[
r_{\rm core}^2={\|F\|_2^2-\|F_{C,\rm core}\|_2^2
\over\|F_{C,\rm core}\|_2^2},\quad
r_{\rm spill}^2={\|F_{C,\rm full}\|_2^2-\|F_{C,\rm core}\|_2^2
\over\|F_{C,\rm core}\|_2^2},
\tag{5.2}
\]

\[
r_{\rm out}^2={\|F\|_2^2-\|F_{C,\rm full}\|_2^2
\over\|F_{C,\rm core}\|_2^2},\qquad
r_{\rm core}^2=r_{\rm spill}^2+r_{\rm out}^2,
\tag{5.2a}
\]

and, when the corresponding geometric channels are disjoint,

\[
r_{\rm self/cross}=
{\sqrt{\|F_{\rm harmonic}\|_2^2+\|F_{\rm low}\|_2^2}
\over\sqrt{\|F_{\rm child}\|_2^2+\|F_{\rm difference}\|_2^2}}.
\tag{5.2b}
\]

The local-FFT backend also computes the self/cross interaction fields
algebraically before summing them, so its self/cross norm ratio remains
available when geometric tags overlap. Quadratic combinations of tag norms
are used as physical partitions only for disjoint recorded supports.

### 5.1 Registered scaling result

The full registered smoke run produced 42 power-width rows, of which 18 have
nonoverlapping parent boxes and exact measurements.  At the largest scale the
available rows are

\[
\begin{array}{c|c|r|c|c|c|c|c|c}
\gamma&W_N&M_N&G_N&R_N&D_N&r_{\rm core}&r_{\rm spill}&r_{\rm out}\\\hline
0.70&18&171500&0.525905&0.0862825&0.00744467&0.329127&0.248&0.216\\
0.75&22&318028&0.710241&0.116493&0.0135705&0.358963&0.244&0.263\\
0.80&27&595508&0.963601&0.157979&0.0249574&0.401087&0.241&0.321.
\end{array}
\tag{5.3}
\]

Here \(r_{\rm core}\) is the complement of the width-\(W_N\) child core,
\(r_{\rm spill}\) is the part inside the full width-\(2W_N\) child sumset but
outside that core, and \(r_{\rm out}\) is outside the full child sumset. They
satisfy \(r_{\rm core}^2=r_{\rm spill}^2+r_{\rm out}^2\). The nonvanishing
spill exposes a shape-broadening obstruction: the original prediction
\(r_{\rm core}=O(N^{\gamma-1})\) is rejected for the filled box. That
prediction is now compared only with \(r_{\rm out}\), while spill is audited
separately.

For \(\gamma=0.70\), the last four core-geometry-eligible points
\(N=16,32,48,64\) give diagnostic exponents
\(\beta_G=0.22315\), \(\beta_D=-0.55315\), and
\(\beta_{\rm out}=-0.19468\). The first two log-space \(R^2\) values are
0.889 and 0.925; the final registered bundle records the independently
recomputed \(r_{\rm out}\) fit quality. These are finite-range diagnostics,
not relay fits:
there are no strict relay-fit-eligible power-width scales.  For
\(\gamma\ge0.75\) fewer than four core-eligible points remain, so no exponent
is reported. In particular the single \(\gamma=4/5,N=64\) point cannot be
called evidence for the predicted \(1/5\) law, even though its \(G_N\) is near
one.

The fixed-relative controls at \(N=64\) give

\[
\begin{array}{c|c|c|c|c|c}
\eta&W_N&G_N&R_N&D_N&c_E\text{ required for }D_N=1/2\\\hline
0.10&6&0.0987762&0.0162139&0.000262889&1901.94\\
0.15&9&0.184325&0.0302541&0.000915309&546.264\\
0.20&12&0.285415&0.0468413&0.00219411&227.883.
\end{array}
\tag{5.4}
\]

The last column merely rescales the frozen quadratic response; it is not a
claim that a full Navier--Stokes evolution at that energy relays or remains in
the ansatz class.

The nominal \(\eta=0.20\) rows at \(N=16,32,48,64\) all realize
\(W_N/N=3/16\). Their energy ratios are respectively
\(0.00183536,0.00209627,0.00216543,0.00219411\). This is consistent with a
small nonzero continuum-response constant for that fixed profile, but four
binary64 Riemann sums neither prove convergence nor provide the two-stage
shape fixed point.

## 6. Small full-Galerkin cross-check

The registered small case uses \(N=4,W=2,\nu=1/40,c_E=1,\tau=1/4\), a cubic
Galerkin cutoff 10, a \(64^3\) padded FFT grid, and 16 RK4 steps. Products of
retained modes stay strictly below Nyquist. The binary64 result is

\[
D_{\rm frozen}=2.0407439931\times10^{-3},\qquad
D_{\rm full}=1.9761634559\times10^{-3},\qquad
{D_{\rm full}\over D_{\rm frozen}}=0.9683544151.
\tag{6.1}
\]

The total energy ratio is \(0.9737674808\), no positive stepwise energy defect
was observed, and final reality/divergence defects are \(6.3\times10^{-18}\)
and \(4.6\times10^{-17}\). Parent evolution therefore does not erase the
already small one-step response in this case, but the response is far below
the target child fraction \(1/2\).

## 7. Exact finite-carrier gate for a second stage

The one-triad child cannot be declared a relay merely because its first
forcing is coherent.  An exact rational search therefore enumerates a small
carrier alphabet before any two-stage phase optimization is interpreted.
Besides the known relay, the first same-shell pair found is

\[
 (0,1,-1),e_1\cos \quad\hbox{and}\quad
 (1,0,-1),e_2\cos
 \longrightarrow (1,1,-2),(1,1,1)\sin/3 .
\tag{7.1}
\]

The two children have nonzero Leray interactions at
((3,2,-1)) and ((1,0,3)), with isolated signed transfers (5/126) and
(7/90).  However, the diagonal parent pairs (A_1+B_2) and (A_2+B_1)
produce exact nonzero coefficients on the target child shell.  For the full
two-relay cross output,

\[
 \|B_{\rm cross}\|_2^2={2483\over1890},\qquad
 \|B_{\rm intended}\|_2^2={37\over315},\qquad
 {\|B_{\rm intended}\|_2^2\over\|B_{\rm cross}\|_2^2}
 ={222\over2483}\simeq0.0894 .
\tag{7.2}
\]

Moreover, the two aligned grandchildren interact exactly to zero.  The full
seven-mode populated field satisfies
(sum_k\Pi_k=0) in `Fraction` arithmetic.  The strict search fixes the known
relay and exhausts all 16 compatible orientations in the stated coordinate
polarization alphabet; none passes the diagonal cross-talk gate.  This is a
finite-scope negative certificate, not an algebraic no-go for larger carrier
sets.

Consequently the sequential one-stage-best-response procedure and the joint
two-stage objective do not yet have a common admissible carrier in this
search class.  Assigning a numerical value to (J_N) would optimize a gadget
already rejected by its exact constraint.  The honest comparison is therefore
`sequential: second stage algebraically present but leakage-dominated` versus
`joint: infeasible under the strict finite carrier gate`.  A larger carrier
alphabet is required before covariance rewards or phase variables are
meaningful.

### 7.1 Full two-stage carrier time evolution

To avoid calling an out-of-cutoff zero curve a second stage, a separate
Galerkin run starts the four exact parents of the rejected partial gadget and
uses a cube cutoff \(3N\).  At \(N=4\), grid \(64^3\), 16 RK4 steps and final
time \(2(\tau N^{-2})\), \(\tau=1/4\), every retained quadratic product lies
strictly below Nyquist.  The final ratios to initial parent energy are

\[
\begin{array}{c|c}
\text{channel}&E/E_{\rm parent}(0)\\\hline
\text{two intended first children}&5.74749\times10^{-4}\\
\text{two diagonal cross-talk children}&5.68886\times10^{-4}\\
\text{intended grandchild modes}&6.35542\times10^{-8}\\
\text{all remaining retained modes}&2.89955\times10^{-4}.
\end{array}
\tag{7.3}
\]

The grandchild signal exceeds the conservative binary64 floor and is stable
under time refinement, but it is **pathway-contaminated**: the cross-talk
children are almost as energetic as the intended children, and other
intermediate modes can enter the same grandchild.  It is therefore a resolved
two-stage shell observation and simultaneously a failure of causal relay
closure, not evidence for a clean cascade.

## 8. Falsification and survival rules

A sublinear-width family is rejected if geometry fails, \(D_N\) has no
scale-independent positive lower bound, low/off-chain forcing grows, the
evolving-parent response vanishes relative to the frozen response, or the
child density leaves the registered next-parent class. Equation (3.6) already
rejects every \(\gamma<1\) under the most favorable capacity scaling.

The fixed-relative family survives only that exponent obstruction. It must
still pass a constant-size child-filling threshold, two-stage density
recurrence, exact carrier cross-talk accounting, and an interval stage
budget. Nothing here constructs an invariant orbit.
