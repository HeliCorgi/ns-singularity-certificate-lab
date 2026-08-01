# Track F shell 制約の finite-floor 訂正と対数臨界階段

作成: 2026-08-01

状態: **SYMBOLIC CANDIDATE / AUDIT REQUIRED**
対象: `docs/research_notes/track_f_shell_constraints.md` §3–§5 の低周波和

## 1. 訂正点

周期領域 \(\mathbb T^3\) の平均ゼロ場では最低非零波数は固定され、二進 shell は
\(\lambda_j=b^j\), \(j=0,\ldots,J\), \(N=b^J\) である。したがって

\[
E_j=A(\lambda_j/N)^\beta,\qquad
N=\tau^{-\gamma},\quad A=\tau^\sigma,\quad \tau=T-t
\tag{1.1}
\]

に対する低周波側の添字 \(m=J-j\) は \(0\le m\le J\) であり、有限時刻
\(t<T\) ごとに有限である。元ノートはこれを \(m\ge0\) と無限まで延長し、
\(\beta\le0\) をその場で排除した。この延長は \(\beta>0\) なら漸近的に無害だが、
\(\beta\le0\) では主項を変える。

この訂正は既存の固定有限帯域 no-go を弱めない。各 \(t<T\) で帯域は有限でも、
\(J(t)\to\infty\) なので軌道全体は一つの固定有限次元空間に留まらない。

## 2. 三つの正しい有限和

### 2.1 Energy

\[
\sum_{j=0}^J E_j
=A\sum_{m=0}^Jb^{-\beta m}
\asymp
\begin{cases}
A,&\beta>0,\\
A\log N,&\beta=0,\\
AN^{-\beta},&\beta<0.
\end{cases}
\tag{2.1}
\]

よって energy 有界条件は

\[
\begin{cases}
\sigma\ge0,&\beta>0,\\
\sigma>0,&\beta=0,\\
\sigma+\gamma\beta\ge0,&\beta<0.
\end{cases}
\tag{2.2}
\]

である。特に \(\beta\le0\) は一律には排除されない。

### 2.2 Enstrophy と総散逸

\[
\sum_{j=0}^J\lambda_j^2E_j
=AN^2\sum_{m=0}^Jb^{-(2+\beta)m}
\asymp
\begin{cases}
AN^2,&\beta>-2,\\
AN^2\log N,&\beta=-2,\\
AN^{-\beta},&\beta<-2.
\end{cases}
\tag{2.3}
\]

従って \(\nu\int_0^T\sum\lambda_j^2E_jdt<\infty\) には

\[
\begin{cases}
\sigma-2\gamma>-1,&\beta\ge-2,\\
\sigma+\gamma\beta>-1,&\beta<-2
\end{cases}
\tag{2.4}
\]

が必要である。\(\beta=-2\) の対数因子は等号を救わない。

### 2.3 臨界 Besov 上界

Bernstein と三角不等式による既存の上界を有限和のまま計算すると

\[
\begin{aligned}
\|u\|_3
&\lesssim \sum_{j=0}^J\lambda_j^{1/2}E_j^{1/2}\\
&=A^{1/2}N^{1/2}
\sum_{m=0}^Jb^{-(1+\beta)m/2}\\
&\asymp
\begin{cases}
A^{1/2}N^{1/2},&\beta>-1,\\
A^{1/2}N^{1/2}\log N,&\beta=-1,\\
A^{1/2}N^{-\beta/2},&\beta<-1.
\end{cases}
\end{aligned}
\tag{2.5}
\]

この右辺が一様有界なら ESS 端点定理で正則となる。従ってこの**上界による
排除を避ける**条件は

\[
\begin{cases}
\sigma<\gamma,&\beta>-1,\\
\sigma\le\gamma,&\beta=-1,\\
\sigma+\gamma\beta<0,&\beta<-1.
\end{cases}
\tag{2.6}
\]

である。ただし右辺の発散だけでは \(\|u\|_3\) の発散を証明しない。

## 3. 修正後の実現可能領域

(2.2), (2.4), (2.6), \(\gamma>0\) を合わせる。

| spectral slope | kinematically surviving exponents |
|---|---|
| \(\beta>0\) | \(0<\gamma<1\), \(\max(0,2\gamma-1)<\sigma<\gamma\), energy由来の \(\sigma=0\) endpoint は \(\gamma<1/2\) で可 |
| \(\beta=0\) | \(0<\gamma<1\), \(\max(0,2\gamma-1)<\sigma<\gamma\), ただし \(\sigma=0\) は対数energyで不可 |
| \(-1<\beta<0\) | \(0<\gamma<1\), \(\max(-\gamma\beta,2\gamma-1)<\sigma<\gamma\); \(-\gamma\beta\) 側だけが最大なら等号可 |
| \(\beta=-1\) | **\(0<\gamma<1,\ \sigma=\gamma\)** の一本の対数臨界境界 |
| \(\beta<-1\) | energy 条件と臨界上界回避条件が矛盾するため空 |

従って \(\gamma<1\) は残るが、\(\beta>0\) は必要ではない。最も鋭い新branchが
\(\beta=-1,\sigma=\gamma\) である。

## 4. 同じnorm指数を持つ physical-space staircase

非零の発散ゼロ template \(W\in C_c^\infty(B(0,\rho);\mathbb R^3)\) を取り、
\(\lambda_j=2^j\) とする。中心 \(x_j\to x_*\) を

\[
x_j=x_*+C\lambda_j^{-1}e_1,\qquad C>3\rho
\tag{4.1}
\]

と選べば supports を互いに素にできる。packet

\[
w_j(x)=\lambda_j W(\lambda_j(x-x_j)),\qquad
U_J=\sum_{j=0}^Jw_j
\tag{4.2}
\]

は各有限 \(J\) で滑らか、発散ゼロ、compact support であり、disjointness により

\[
\begin{aligned}
\|U_J\|_2^2
&=\|W\|_2^2\sum_{j=0}^J\lambda_j^{-1}<2\|W\|_2^2,\\
\|U_J\|_3^3
&=(J+1)\|W\|_3^3,\\
\|\nabla U_J\|_2^2
&=\|\nabla W\|_2^2\sum_{j=0}^J\lambda_j\asymp N,\\
\|\omega_J\|_\infty&\asymp N^2,\\
\|(U_J\cdot\nabla)U_J\|_2&=O(N^{3/2}).
\end{aligned}
\tag{4.3}
\]

ここで \(N=2^J\)。\(N(\tau)=\tau^{-\gamma}\) なら

\[
J\sim {\gamma\over\log2}\log{1\over\tau},\qquad
\|U_J\|_3^3\asymp\log{1\over\tau}.
\tag{4.4}
\]

総散逸は \(\int_0\tau^{-\gamma}d\tau<\infty\) iff \(\gamma<1\)、BKM必要発散は
\(\int_0\tau^{-2\gamma}d\tau=\infty\) iff \(\gamma\ge1/2\)。従って

\[
\boxed{\frac12\le\gamma<1}
\tag{4.5}
\]

は energy、総散逸、ESS、BKM の**指数条件だけ**とは整合する。

式 (4.2) は時刻ごとの kinematic field であり、NS 軌道ではない。さらに compact
support と exact Fourier band limitation は両立しないため、\(\|w_j\|_2^2\) は
packet全energyであって §2 の exact shell energyではない。両者の接続には

\[
c\lambda_j^{-1}\le\|P_jU_J\|_2^2\le C\lambda_j^{-1}
\tag{4.6}
\]

および off-shell leakage/cross-term boundが必要である。滑らかな初期値から packetを
順に増幅する力学と合わせ、これは未証明義務である。

## 5. Moving front と pressure relay

packet energy は \(E_N=c_E/N\)。Bernstein estimate が与えるのは signed flux の
上限

\[
|\Pi_N|\le C_{\rm NL}N^{5/2}E_N^{3/2}
=C_{\rm NL}c_E^{3/2}N
\tag{5.1}
\]

であり、正符号ではない。child粘性、wake、off-chain lossを引いた interval budgetが
\(q_*N\), \(q_*>0\) という別仮定を満たすときだけ、child energy \(c_E/(2N)\) で
割って連続front envelope

\[
\dot J=\kappa N^2,\qquad
\dot N=(\log2)\kappa N^3.
\tag{5.2}
\]

従って

\[
N(t)=\left[N_0^{-2}-2(\log2)\kappa(t-t_0)\right]^{-1/2},
\tag{5.3}
\]

すなわち \(\gamma=1/2\)。離散 \(N=2^J\) は微分不能なので、これは activation
times の連続包絡である。front の有限時刻発散を**仮定した signed flux から導いた
有効 ODE**であって、PDE flux 下界ではない。

非局所 pressure が次段の strain を補助できる符号も leading order では存在する。
transmitter の \(M_{ij}=\int u_i u_jdx\) に対し

\[
p(x)=M_{ij}\partial_i\partial_j(4\pi|x|)^{-1}+O(|x|^{-4}),
\quad
\nabla^2p=M_{ij}\nabla^4_{ij}(4\pi|x|)^{-1}+O(|x|^{-6}).
\tag{5.4}
\]

rank-one \(M\) は非零局在発散ゼロ場では実現不能である。明示的な
\(W_a=(a\times x)e^{-|x|^2/2}\), \(a=e_1\) は
\(M/\operatorname{tr}M=\operatorname{diag}(0,1/2,1/2)\) を持ち、距離方向
\(n=e_1\) なら

\[
{4\pi d^5\nabla^2p\over\operatorname{tr}M}
=\operatorname{diag}(-12,6,6).
\tag{5.5}
\]

従って forced gradient equation
\(D_t\nabla u=-(\nabla u)^2-\nabla^2p+\nu\Delta\nabla u+\nabla f\)
の \(-\nabla^2p\) は favourable direction を持つ。source packetで
\(M\asymp\lambda^{-1}\)、距離 \(d\asymp\lambda^{-1}\) なら
\(\nabla^2p\asymp\lambda^4\)。stage 時間 \(\Delta t\asymp\lambda^{-2}\) に
strain \(O(\lambda^2)\) を変え得るため、scale 次元は閉じる。ただし support/separation
比は固定で multipole remainder は段とともに減らず、pressureはLeray triadと独立の
energy源ではない。

## 6. Pilot と判定

実装:

- `src/ns_certificate_lab/zeno_packet_relay.py`
- `experiments/run_zeno_packet_relay_pilot.py`
- `configs/zeno_packet_relay_pilot_v2.json`
- `outputs/zeno_packet_relay_pilot_v2/`

binary64 pilot は \(J=8,16,24,32\) で energy \(<2\)、critical shell sum
\(=J+1\)、enstrophy/\(N\to2\) を確認した。実現可能な rank-two Gaussian vortex
momentで (5.5) の12、20,000 orientation の最大 normalized stretch
\(11.9976793898\) を再現した。
front ODE の \(N\sqrt{T-t}\) は相対 spread \(2.3\times10^{-16}\)。

v1のrank-one tensorは物理packet witnessとして棄却済みである。v2の \(d^{-5}\)
checkもleading homogeneous formulaの回帰で、exact pressure/remainderではない。
これらは中心仮定をまだ試していない。最短 kill test は、真の divergence-free
wave-packet basis について一段の NS triad map を構成し、

\[
\inf_n{\Pi_{n\to n+1}-\nu\lambda_{n+1}^2E_{n+1}\over\lambda_n}>0
\tag{6.1}
\]

と off-chain leakage の幾何級数上界を同時に満たせるか調べること。満たせなければ
pressure-relay branch を `REJECTED` とする。

## 7. 証明へ進む場合の鎖

1. compact divergence-free packet template と disjoint accumulating centers を固定。
2. 一段の full 3D NS cascade cell（pressureを含む）を構成。
3. 正の scale-normalized child flux と viscous margin を証明。
4. off-chain modes、coarse backreaction、packet overlap を summable に抑える。
5. smooth finite-band datum または smooth low-mode force から最初の cell へ接続。
6. relative-periodic shift orbit を全段で Newton–Kantorovich 検証。
7. Zeno 時間 \(\sum\Delta t_j<\infty\)、energy、局所energy、総散逸を閉じる。
8. disjoint-packet lower boundから \(\|u(t)\|_3^3\to\infty\) を閉じる。
9. \(t<T\) の古典解存在・一意性と \(T\) での延長不能を接続。
10. Clay (C) または (D) の初期値・外力条件へ最終接続。
