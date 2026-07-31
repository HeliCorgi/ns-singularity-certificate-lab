# CANDIDATE SOLUTION — Zeno critical packet relay

**Status: unverified solution candidate; not yet a PDE solution candidate**

**Label: FORMAL KINEMATIC SCALING SCENARIO / AUDIT REQUIRED**
**Created: 2026-08-01**

これは3次元 Navier–Stokes 特異点の証明ではない。有限 energy・有限総散逸と
臨界 \(L^3\) 発散を同時に満たす明示的な多尺度形、およびその形を生成し得る
pressure/triad relay の必要式を一つの反証可能な候補へまとめたものである。

## 1. 狙う Clay 公式命題

主標的は **Clay (D)**:

- 領域: \(\mathbb T^3=(\mathbb R/2\pi\mathbb Z)^3\)。
- 粘性: 固定 \(\nu>0\)。
- 初期値: 平均ゼロ、\(C^\infty\)、発散ゼロ。最初は有限 Fourier datum でもよい。
- 外力: 低波数有限台の \(C^\infty_{x,t}\) 発散ゼロ外力。
- 目標: ある有限 \(T\) まで古典解で、\(T\) を越えて同じ強解classへ延長不能。

高周波は外力で直接作らない。外力は coarse reservoir と低mode位相だけを制御し、
front は真の Leray 非線形項が作る。この条件は Clay smooth forcing の高周波急減衰と
整合する。外力を最終的に0へ落とせれば、同じ機構は (C) のより強い無外力版に
なるが、本候補では主張しない。

## 2. 完全な候補式

### 2.1 幾何学的 packet stack

\(W\in C_c^\infty(B(0,\rho);\mathbb R^3)\) を非零、発散ゼロとし、
\(\lambda_j=2^j\)。回転 \(R_j\in SO(3)\) と中心を

\[
x_j=x_*+C\lambda_j^{-1}e_1,\qquad
C>3\rho,\qquad x_j\to x_*
\tag{2.1}
\]

で選ぶ。静的な基準 packet は

\[
w_j(x)=\lambda_jR_j
W\!\left(\lambda_jR_j^T(x-x_j)\right).
\tag{2.2}
\]

各有限和は滑らか、発散ゼロで、support を互いに素にできる。

For the periodic target, the construction starts at a sufficiently large
\(j_0\) so all centers lie inside one injectivity chart.  Exact scale reuse on
the torus is not asserted: rotations must be lattice symmetries, or each stage
must absorb the periodic Green regular part as a scale-dependent perturbation.

### 2.2 Zeno schedule and smooth pre-seeds

one-cell rescaled timeを \(a_c>0\) とし、stage の終了時刻を

\[
t_j=T-{4a_c\over3}\lambda_j^{-2},\qquad
\Delta t_j=t_{j+1}-t_j=a_c\lambda_j^{-2}
\tag{2.3}
\]

とする。従って \(\sum_j\Delta t_j<\infty\)。Heaviside 的に非零 packet を
追加すれば時間不連続になるため、動的候補は全 future scale の滑らかな seed を
最初から含む収束級数としてのみ読む。local cascade cell \(V_j(y,\theta)\) を用いて

\[
u_{\rm ans}(x,t)
=u_{\rm res}(x,t)
+\sum_{j\ge j_0}
\lambda_jR_j
V_j\!\left(
\lambda_jR_j^T(x-x_j),
\theta_j(t)
\right),
\quad
\theta_j=\lambda_j^2(t-t_j).
\tag{2.4}
\]

と置く。\(u_{\rm res}\) は coarse energy reservoir。必要条件は
\(V_j(\cdot,\theta)\) が \(\theta\to-\infty\) で全必要normについて
super-algebraically small となり、\(\theta\ge a_c\) で wakeへ滑らかに接続すること。
従って各 \(t<T\) で和は無限だが収束し、front は exact support端でなく resolved
energy quantile である。連続 front index と観測 staircase を

\[
j_c(t)={1\over2\log2}\log\left({4a_c\over3(T-t)}\right),\qquad
J(t)=\lfloor j_c(t)\rfloor,\qquad N_c(t)=2^{j_c(t)}
\asymp(T-t)^{-1/2}
\tag{2.5}
\]

と置く。

### 2.3 relative-periodic cascade cell

必要な cell 条件は、親 packet と coarse background を初期状態に持つ full 3D
Navier–Stokes flow map \(\Phi\) が、rescaled時間 \(a_c\) 後に半scaleの子を
滑らかに作ること:

\[
\Phi_{a_c}(V_{\rm parent}+B_{\rm coarse})
=V_{\rm persist}
+2R_*V_{\rm child}(2R_*^T(\,\cdot-y_*))
+r_*.
\tag{2.6}
\]

ここで \(r_*\) は次段へ移すと幾何級数的にsummableでなければならない。
shell index と rescaled time の traveling-front 表現では

\[
b_{j+1}(s+2\log2)=b_j(s),\qquad
b_j(s)=q(j-J(s)),\qquad s=-\log(T-t),
\tag{2.7}
\]

となる。これは bounded periodic profile ではなく、front の後ろに packet wake を
一つずつ残す relative-periodic orbit である。torusでは (2.7) は leading cell
relationであり、実際には periodic Green regular partを含む \(V_j=V+o(1)\) の
stage-dependent perturbation theoremが必要になる。

### 2.4 smooth forcing cancellation 条件

外力を残差逆設計するだけでは不十分である。本候補は最終的に固定有限 low set
\(L\) だけを強制する強い版を採用する:

\[
P_{L^c}f\equiv0,\qquad
\sup_{t\ge0}\|\partial_t^qf(t)\|_{H^s}<\infty
\quad(\forall q,s\ge0).
\tag{2.8}
\]

従って high-scale cell は、単に小さい残差でなく

\[
P_{L^c}\left[
\partial_tu_{\rm ans}
+\mathbb P(u_{\rm ans}\cdot\nabla u_{\rm ans})
-\nu\Delta u_{\rm ans}
\right]
=0
\tag{2.9}
\]

を満たさなければならない。さらに low coefficients は \(t=T\) を越えて明示的に
\(C^\infty\) 延長される必要がある。各 leading block は pointwise \(O(N^3)\)
なので、(2.9) は高周波欠陥を外力へ隠すことを禁止し、事実上「各cellが
asymptotically unforced NS を解く」ことを要求する。段ごとの low-mode 位相変更は
時間微分を発散させ得るため、\(T\) 近傍では force を flat/fixed とし、phase locking
は cell 内の自律力学に委ねる。

## 3. shell law と finite-floor 訂正

周期 shell \(\lambda_j=2^j\), \(0\le j\le J\) に

\[
E_j=A(\lambda_j/N)^\beta,\qquad
N=\tau^{-\gamma},\quad A=\tau^\sigma,\quad \tau=T-t
\tag{3.1}
\]

を置く。最低波数が1なので低周波和は有限である。本候補の Fourier shell model は

\[
\boxed{\beta=-1,\quad\sigma=\gamma,\quad E_j=\lambda_j^{-1}}
\tag{3.2}
\]

を使う。すると

\[
\sum_{j=0}^JE_j<2,\qquad
\sum_{j=0}^J\lambda_j^2E_j\asymp N,\qquad
\sum_{j=0}^J\lambda_j^{1/2}E_j^{1/2}=J+1.
\tag{3.3}
\]

旧 shell classifier が \(\beta\le0\) を排除したのは、存在しない
\(\lambda<1\) shell まで和を延長したためである。完全な訂正は
`track_f_shell_constraints_finite_floor_erratum.md` にある。

ただし次節の compact packet は exact Fourier annulus に局在しない。不確定性原理に
より exact compact support と exact band limitation は両立しないから、(3.2) を
実際の Littlewood--Paley energyへ接続するには

\[
c\lambda_j^{-1}\le\|P_jU_J\|_2^2\le C\lambda_j^{-1}
\tag{3.4}
\]

と off-shell leakage/cross terms を証明する必要がある。現時点では (3.2) は shell
algebra、次節は同じnorm指数を持つ physical-space scaffoldであり、同一対象ではない。

## 4. 臨界 \(L^3\) の明示下限

静的 packet (2.2) は

\[
\|w_j\|_2^2=\lambda_j^{-1}\|W\|_2^2,\qquad
\|w_j\|_3^3=\|W\|_3^3,\qquad
\|\nabla w_j\|_2^2=\lambda_j\|\nabla W\|_2^2.
\tag{4.1}
\]

disjoint support なら

\[
\left\|\sum_{j=0}^Jw_j\right\|_3^3
=(J+1)\|W\|_3^3.
\tag{4.2}
\]

真の粘性 flow は直ちに空間 tail を持つため、最終証明では disjoint support を
使えない。代わりに互いに素な core balls \(B_j\) 上で

\[
\int_{B_j}|u_j(t)|^3dx\ge c_3>0,\qquad
\|u(t)-u_j(t)\|_{L^3(B_j)}
\le {1\over2}\|u_j(t)\|_{L^3(B_j)}
\tag{4.3}
\]

を一様に証明する必要がある。左辺は reservoir、cell defect、periodic/parabolic tailを
全て含む。(4.3) なら

\[
\|u(t)\|_3^3\ge{c_3\over8}(J(t)+1)
\asymp c\log{1\over T-t}\to\infty.
\tag{4.4}
\]

## 5. Moving Fourier front ODE

局在 packet に対する Bernstein estimate が与えるのは signed flux の上限だけである:

\[
|\Pi_N|\le C_{\rm NL}N^{5/2}E_N^{3/2}.
\tag{5.1}
\]

正符号は別に証明しなければならない。\(E_N=c_EN^{-1}\) とすると child energy は
\(E_{2N}=c_E/(2N)\)、shape係数 \(C_\nu>0\) を含む child粘性損失は
\(2\nu C_\nu c_EN\)。一段の wake/off-chain lossを \(C_LN\) 以下とし、真の
Leray cellについて

\[
q_*=c_{\rm flux}c_E^{3/2}-2\nu C_\nu c_E-C_L>0
\tag{5.2}
\]

という signed **lower** marginを仮定する。必要な interval budget は

\[
\int_{t_j}^{t_{j+1}}
 (\Pi_{j\to j+1}-D_{j+1}-L_{\rm wake}-L_{\rm off})\,dt
\ge E_{j+1}(t_{j+1}).
\tag{5.3}
\]

このmarginを child energyで割ると、連続 front envelope は

\[
\dot J=\kappa N^2,\qquad
\dot N=kN^3,\qquad
\kappa={2q_*\over c_E},\quad k=(\log2)\kappa.
\tag{5.4}
\]

従って

\[
N(t)=\left[N_0^{-2}-2k(t-t_0)\right]^{-1/2},
\tag{5.5}
\]

で finite-time front が得られる。離散 \(N=2^J\) 自体は微分不能なので、(5.4)
は shell activation times の連続包絡である。(2.3) と一致させるには
\(a_c=3/(8k)\)。この ODE の導出で未証明なのは (5.1) からは出ない正の signed
margin (5.2) と interval budget (5.3) である。

## 6. Pressure-Hessian relay

3D NS の gradient equation は

\[
D_tA=-A^2-\nabla^2p+\nu\Delta A+\nabla f,\qquad A=\nabla u.
\tag{6.1}
\]

periodic target の exact pressure は torus Green functionで

\[
p(x)=\int_{\mathbb T^3}\partial_i\partial_jG_{\mathbb T^3}(x-y)
u_i(y)u_j(y)\,dy,
\qquad G_{\mathbb T^3}(z)={1\over4\pi|z|}+H_{\mathbb T^3}(z)
\tag{6.2}
\]

と書く。localized template \(W\) の self-pressure Hessian は、親座標
\(x=x_j+\lambda_j^{-1}z\) で

\[
\nabla^2p_j(x)=\lambda_j^4\mathcal H_W(z)+O_{\mathbb T^3}(\lambda_j^{-1}),
\quad
\mathcal H_W(z)=\int_{\mathbb R^3}D^4{1\over4\pi|z-y|}:
[W(y)\otimes W(y)]\,dy.
\tag{6.3}
\]

遠方でのみ \(\mathcal H_W(z)=D^4(4\pi|z|)^{-1}:M(W)+R_W(z)\) と展開できる。
4階微分の leading tensor は

4階微分は

\[
\begin{aligned}
4\pi d^5\partial_{abij}|x|^{-1}
={}&105n_an_bn_in_j\\
&-15(\delta_{ab}n_in_j+\delta_{ai}n_bn_j+\delta_{aj}n_bn_i
+\delta_{bi}n_an_j+\delta_{bj}n_an_i+\delta_{ij}n_an_b)\\
&+3(\delta_{ab}\delta_{ij}+\delta_{ai}\delta_{bj}+\delta_{aj}\delta_{bi}).
\end{aligned}
\tag{6.4}
\]

rank-one moment は非零局在発散ゼロ場では実現不能である。代わりに明示的な
Schwartz vortex

\[
W_a(x)=(a\times x)e^{-|x|^2/2},\qquad a=e_1,\qquad
{M(W_a)\over\operatorname{tr}M(W_a)}={1\over2}\operatorname{diag}(0,1,1)
\tag{6.5}
\]

を使う。距離方向 \(n=e_1\) なら

\[
{4\pi d^5\nabla^2p\over\operatorname{tr}M}
=\operatorname{diag}(-12,6,6).
\tag{6.6}
\]

従って \(-\nabla^2p\) は \(e_1\) 方向へ normalized strength \(12\) の favourable
strain acceleration を持つ。しかし親半径と親子距離の比は全段で固定なので、
leading multipole 誤差は段とともに消えない。\(C/\rho\) を大きくすれば remainder は
減る一方、relay strength は \(C^{-5}\) で弱くなる。このtradeoffを exact
\(\mathcal H_W\) で通す必要がある。

\(M_j\asymp\lambda_j^{-1}\), \(d_j\asymp\lambda_j^{-1}\) なので

\[
\nabla^2p_{j\to j+1}\asymp\lambda_j^4,\qquad
\Delta t_j\nabla^2p_{j\to j+1}\asymp\lambda_j^2,
\tag{6.7}
\]

これは child strain scale と同じである。ただし他packet/reservoirとのcross pressure
も同じ次数になり得る。pressure は curl-free で新しい渦度や全 divergence-free shell
energyを直接作らないうえ、Leray triad の補成分であって独立のenergy源ではない。
従って (6.7) は exact projected cell evolutionの**診断的位相条件**としてのみ使い、
flux marginへ二重加算しない。

## 7. Scaling table

一般に \(N=\tau^{-\gamma}\)、本候補は moving-front law から \(\gamma=1/2\)。
ノルムは packet core が支配するとする。

| quantity | general \(\gamma\) | \(\gamma=1/2\) |
|---|---:|---:|
| energy \(\|u\|_2^2\) | \(O(1)\) | \(O(1)\) |
| enstrophy \(\|\omega\|_2^2\) | \(\tau^{-\gamma}\) | \(\tau^{-1/2}\) |
| global \(\|u\|_3^3\) | \(\asymp\log(1/\tau)\) | \(\asymp\log(1/\tau)\) |
| global \(\|u\|_3\) | \(\asymp\log(1/\tau)^{1/3}\) | same |
| \(\|\omega\|_\infty\) | \(\tau^{-2\gamma}\) | \(\tau^{-1}\) |
| dissipation rate \(\nu\|\nabla u\|_2^2\) | \(\nu\tau^{-\gamma}\) | \(\nu\tau^{-1/2}\) |
| remaining total dissipation | \(O(\tau^{1-\gamma})\) | \(O(\tau^{1/2})\) |
| \(\|(u\cdot\nabla)u\|_2\) | \(O(\tau^{-3\gamma/2})\) | \(O(\tau^{-3/4})\) |
| \(\|\nabla p\|_2\) | \(O(\tau^{-3\gamma/2})\) | \(O(\tau^{-3/4})\) |
| \(\|p\|_2\) | \(O(\tau^{-\gamma/2})\) | \(O(\tau^{-1/4})\) |
| velocity maximum | \(\tau^{-\gamma}\) | \(\tau^{-1/2}\) |
| physical stage time | \(N^{-2}\) when \(\gamma=1/2\) | \(\tau\) |
| Fourier bandwidth | \(\tau^{-\gamma}\) | \(\tau^{-1/2}\) |

BKM necessary integralが発散するscale条件は \(2\gamma\ge1\)、総散逸は有限 iff
\(\gamma<1\)。従って指数窓 \(1/2\le\gamma<1\) が残り、front ODE はその左端を
選ぶ。ただし \(\gamma=1/2\) は \(\sqrt\tau\|u\|_\infty\asymp1\) の Type-I 境界で
あり、一般3D ancient-limit/endpoint監査が別に必要である。

## 8. 特異feedback loop

\[
\begin{aligned}
&\text{low-mode reservoir/phase}
\xrightarrow{\text{true helical triads}}
\Pi_j>\nu\lambda_j^2E_j\\
&\xrightarrow{(5.4)}
\text{child packet at }2\lambda_j
\xrightarrow{\text{anisotropic }M_j}
-\nabla^2p_{j\to j+1}>0\\
&\xrightarrow{(6.1)}
\text{child strain/phase locking}
\xrightarrow{}
\text{next positive triad flux}.
\end{aligned}
\tag{8.1}
\]

各段の energy は \(\lambda_j^{-1}\) なので総和可能だが、各段の臨界 \(L^3\) mass
は一定なので段数とともに発散する。

## 9. 既存障害の回避と衝突点

| obstruction | audit |
|---|---|
| energy bounded | \(\sum2^{-j}<\infty\) |
| finite dissipation | \(\int\tau^{-1/2}dt<\infty\) |
| ESS \(L^3\) | (4.4) が閉じれば endpoint regularity仮定の外。発散から非延長を得る最後の論理はESSの逆ではなく smooth extensionの \(L^3\) 連続性 |
| fixed finite bandwidth | \(N(t)\to\infty\) |
| pure-swirl initial decrease | 一般3D helical packet、pure swirlでない |
| one-scale/DSS no-go | bounded fixed profileではないが、wakeだけで全Liouville仮定を外したとはまだ言えない。relative-front ancient limitを別監査 |
| Type-I endpoint | \(\sqrt\tau\|u\|_\infty\asymp1\)。軸対称Type-I no-goは直接適用外だが一般3D endpoint/ancient-limit監査が必要 |
| weak \(L^3\) / CKN | strong \(L^3\) は対数発散しても weak-\(L^3\) と scale-invariant local energy が有界になり得る。epsilon thresholdを別監査 |
| periodic image/wall | (D) は周期領域そのもの。packet間距離をboxより十分小さくし periodic Green regular part を別評価 |
| Galerkin global existence | 各cutoffはfrontを止めるだけ。cutoff時刻を特異時刻と読まない |
| smooth force high-frequency decay | high-mode forcingをexact 0、全time derivativesと \(T\) 後extensionを (2.8) で必須化 |
| front resolution | Fourier packet係数で直接表し、shell cutoffを増やして convergence |
| averaged-NS warning | energy cancellationだけでは足りない。本候補は平均化前の真のLeray triad係数でcellを作る義務を残す |

## 10. 数値・記号的支持

実装:

- `src/ns_certificate_lab/zeno_packet_relay.py`
- `experiments/run_zeno_packet_relay_pilot.py`
- `configs/zeno_packet_relay_pilot_v2.json`
- `outputs/zeno_packet_relay_pilot_v2/`

pilot の結果:

- \(J=8,16,24,32\): energy \(<2\)、enstrophy/\(N\to2\)、critical packet mass
  と modeled Besov shell sum はそれぞれ \(J+1\)（同一templateの同時exact性は
  主張しない）。
- 実現可能な rank-two Gaussian vortex momentで (6.6) の normalized stretch
  \(12\) を再現した。
- 20,000 random physical vortex orientations の最大 normalized stretch は
  \(11.9976793898\)。
- leading tensorの \(d^{-5}\) relative spread は0だが、これはhomogeneous closed
  formulaの回帰試験であり、exact pressureやremainderの検査ではない。
- front ODE の \(N\sqrt{T-t}\) relative spread は
  \(2.22\times10^{-16}\)。
- bundle verifier、tamper rejectionを含む新規 tests は22件合格。

初回 outputs/zeno_packet_relay_pilot_v1/ の rank-one moment は局在発散ゼロ場から
実現不能と判明したため、物理packet witnessとして棄却した。履歴は上書きせず保存した。
v2も kinematic/scaling/leading-tensor testにすぎず、packet-shell bridge、exact
pressure、cell (2.6)、または interval flux budget (5.3) を支持していない。

## 11. 最短の反証実験

一段だけを真の helical Fourier packet で試す。

1. 親/子 annulus に整数格子modeを置き、真の helical interaction coefficient を計算。
2. 親位相、回転、中心、coarse low modesを最適化。
3. 測定量を
   \[
   c_{\rm net}^{(j)}
   ={\Pi_{j\to j+1}-D_{j+1}-L_{\rm wake}-L_{\rm off}\over\lambda_j},
   \qquad
   \varepsilon_{\rm leak}^{(j)}
   ={\Pi_{\rm off-chain}\over\Pi_{j\to j+1}}
   \tag{11.1}
   \]
   とする。
4. \(j=2,3,4,5\) で \(c_{\rm net}^{(j)}\) が正の同一marginを持ち、
   \(\varepsilon_{\rm leak}^{(j)}\) が減るかを測る。

**Kill condition:** 最適化後も正marginが2段連続で得られない、phase coneが線形に
不安定、または off-chain leakage が主flux以上なら、この明示packet/mode familyを
`REJECTED` とする。他のpressure-relay全般を排除したとは読まない。

解像度は最初に各shell 30–300 modes、4 shells。探索はbinary64でよい。昇格時は
helical basis の平方根、triad係数、全off-chain sumを区間複素演算で囲う。

## 12. 未証明補題

1. 真の3D Leray係数で positive invariant phase cone が存在する。
2. physical packet energy と Littlewood--Paley shell law の bridge (3.4)。
3. interval flux budget (5.3) が全段で一様正。
4. off-chain leakage と coarse backreaction がsummable。
5. exact torus pressure (6.3) で remainder/strength tradeoff後も favourable strain。
6. smooth future seedsの無限和と局所 \(L^3\) lower bound (4.3)。
7. low-mode smooth force/finite Fourier datum から最初のcellへ入る。
8. relative-front orbitが full PDE の古典解として全 \(t<T\) で存在。
9. forceが (2.8) を満たし \(t=T\) を越えて \(C^\infty\)、かつlocal energy/pressureを制御。
10. classical uniquenessと \(L^3\) 連続性によるClay (D)への忠実な接続。

## 13. 最終証明鎖

1. packet/triad graphと低mode forcingを明示する。
2. 一段cellを Newton–Kantorovich/radii polynomial で検証する。
3. phase cone、positive flux、leakage contractionを区間演算で閉じる。
4. torus Green regular partを含むstage-dependent perturbation theoremを全段へ適用する。
5. Zeno時刻 (2.3) と全 \(t<T\) の古典解を構成する。
6. energy・総散逸・局所energy・forcing smoothnessを証明する。
7. packet lower boundsを足し \(\|u(t)\|_3^3\to\infty\) を証明する。
8. global smooth alternativeがあれば classical uniquenessでこの解と一致し、\(T\)
   近傍の \(L^3\) 連続有界性が (4.4) と矛盾することを示す。
9. \(\widetilde u(y,s)=2\pi u(2\pi y,(2\pi)^2s)\),
   \(\widetilde f(y,s)=(2\pi)^3f(2\pi y,(2\pi)^2s)\) で period 1 のClay公式 (D)へ写す。
10. 独立実装、区間証明、最終論理接続を監査する。

## 14. 現時点の判定

この候補の新しい点は、(i) repo の shell 有限和誤りを修正して現れた
\(E_j\sim\lambda_j^{-1}\) の等臨界wake、(ii) front後方へscale-critical packetを
残すことで bounded periodic profile no-go を外すこと、(iii) pressure quadrupoleを
次段strainのrelayとして使うことである。

最大の穴は明確である: **真の Navier–Stokes triad で (2.6), (5.3) を満たす
stage-dependent cascade cell が存在するか。** ここが否定されればこのpacket familyは
終了する。
