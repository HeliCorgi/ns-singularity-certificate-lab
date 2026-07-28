# 方程式監査

## 1. 監査規則

本書の状態は、指定された次の5語だけを用いる。

| 状態 | 意味 |
|---|---|
| **導出済み** | 本書で定義から独立に代数導出し、符号・係数・次元を検算した。 |
| **一次資料で確認済み** | 独立導出に加え、査読論文、著者プレプリント、または Clay 公式記述の同じ規約の式と照合した。 |
| **未確認** | 導出または一次資料照合が未完了。実装の前提にしてはならない。 |
| **不整合** | 規約または資料間で一致せず、原因を切り分け中。実装の前提にしてはならない。 |
| **誤り** | 反例または代数検算により誤りと確定した。実装してはならない。 |

この監査で実装の規範となる式には **未確認・不整合・誤りはない**。資料中の表記上の不整合は §11 に隔離した。

符号規約は

\[
\boldsymbol e_\theta=(-\sin\theta,\cos\theta,0),\qquad
\boldsymbol\omega=\nabla\times\boldsymbol u,\qquad
\omega^\theta=\partial_zu^r-\partial_ru^z,
\tag{A0}
\]

\[
u^r=-\partial_z\psi^\theta,\qquad
u^z=r^{-1}\partial_r(r\psi^\theta)
\tag{A1}
\]

で固定する。別の \(\psi^\theta\) 符号規約の式を混ぜない。

## 2. 監査要約

| ID | 対象 | 状態 | 主な根拠 |
|---|---|---|---|
| E-01 | 3D Navier–Stokes と \(\nabla\cdot u=0\) | **一次資料で確認済み** | [F, eqs. (1)–(3)] |
| E-02 | 軸対称の3D物理発散 | **一次資料で確認済み** | 直接導出、[HL09, eq. (2.5)]、[LW, eq. (2.3)] |
| E-03 | 円柱3成分運動方程式 | **一次資料で確認済み** | 基底微分から導出、[HLW, eq. (2.1a)]、[Z23, eq. (1.12)] |
| E-04 | 渦度3成分 | **一次資料で確認済み** | curl の直接計算、[HL, eq. (12)] |
| E-05 | \(u^\theta\) 発展式 | **一次資料で確認済み** | E-03、[HL, eq. (13)]、[HLW, eq. (2.1a)] |
| E-06 | \(\omega^\theta\) 発展式 | **一次資料で確認済み** | 交差微分、[HL, eq. (14)]、[HLW, eq. (2.1b)] |
| E-07 | 流れ関数から速度回復 | **一次資料で確認済み** | curl の直接計算、[HL, eq. (16)]、[HLW, eq. (2.1d)] |
| E-08 | \(\psi^\theta\) 楕円式 | **一次資料で確認済み** | curl–curl の直接計算、[HL, eq. (15)] |
| E-09 | \(u_1,\omega_1,\psi_1\) 定義 | **一次資料で確認済み** | [HLW, eq. (2.2)] |
| E-10 | \((\Delta_0-r^{-2})(rf)=r\mathcal L_5f\) | **導出済み** | 積微分を展開 |
| E-11 | \(u_1\) 発展式 | **一次資料で確認済み** | E-05/E-10、[HLW, eq. (2.3a)] |
| E-12 | \(\omega_1\) 発展式 | **一次資料で確認済み** | E-06/E-10、[HLW, eq. (2.3b)] |
| E-13 | \(\psi_1\) 楕円式 | **一次資料で確認済み** | E-08/E-10、[HLW, eq. (2.3c)] |
| E-14 | \(u^r,u^z\) 回復式 | **一次資料で確認済み** | E-07、[HLW, eq. (2.3d)] |
| E-15 | 回復速度の3D発散ゼロ | **導出済み** | 項別相殺 |
| E-16 | 軸の極条件・偶奇性 | **一次資料で確認済み** | [LW, Cor. 1, Lemma 2]、[HL, eqs. (18)–(21)] |
| E-17 | \(\mathcal L_5\) の軸上極限 | **導出済み** | 偶 Taylor 展開 |
| E-18 | Cartesian 3D速度・渦度復元 | **導出済み** | 円柱基底へ代入 |
| E-19 | 圧力 Poisson 式 | **導出済み** | E-01 の発散 |
| E-20 | 物理エネルギーと測度 | **一次資料で確認済み** | [F, bounded energy]、[HLW, eq. (1.4a)] |
| E-21 | Navier–Stokes スケーリング | **一次資料で確認済み** | 代入、[HLW, eqs. (1.3), (2.5)] |
| E-22 | 変数の物理次元 | **導出済み** | 次元解析 |
| E-23 | 循環 \(\Gamma\) 方程式 | **導出済み** | E-05 へ \(\Gamma=ru^\theta\) を代入 |
| E-24 | 閉じた系から3D原始変数系への同値性 | **一次資料で確認済み** | 正則性・極条件込みで [LW] |
| E-25 | 有限円柱Poisson対照問題 | **導出済み** | E-13へ周期 \(z\)・外側Dirichletを明示 |
| E-26 | \(r^3\)-flux離散式・軸係数8 | **導出済み** | control-volume積分とFourier変換 |
| E-27 | Hou 円柱の壁条件(\(r=1\)) | **一次資料で確認済み** | 独立導出、[Hou21, eqs. (2.4)–(2.5)] |
| E-28 | Hou 円柱の対称性・半周期領域 | **一次資料で確認済み** | [Hou21, §2]、[LW] |
| E-29 | Hou 初期値と導出ノルム | **一次資料で確認済み** | [Hou21, eq. (2.2)]、ノルム値は導出済み |
| E-30 | Hou 数値プロトコル(二段階粘性等) | **一次資料で確認済み** | [Hou21, §1.4, §3, App. A]、[HH21, App. A] |
| E-31 | 壁渦度条件の2次離散式 | **導出済み** | 壁での Taylor 展開 |
| E-32 | 壁依存性試験用 \(C^\infty\) envelope 初期値族 | **導出済み** | E-29 に compact support カットオフを乗算、数値検算済み |

## 3. 元の3次元方程式

### E-01: 原始変数系

\[
\partial_tu_i+\sum_{j=1}^3u_j\partial_ju_i
=-\partial_ip+\nu\Delta u_i,\qquad
\sum_{i=1}^3\partial_i u_i=0,\qquad \nu>0.
\tag{E-01}
\]

**状態: 一次資料で確認済み。**
Clay 公式記述 [F, eqs. (1)–(3)] と一致する。圧力は密度で割っている。探索は外力 \(f=0\) を採用する。

次元は各運動項が \(LT^{-2}\)、発散が \(T^{-1}\) で一致する。

### E-02: 軸対称の物理的な発散

\[
\nabla\cdot\boldsymbol u
={1\over r}\partial_r(ru^r)+\partial_z u^z
=\partial_ru^r+{u^r\over r}+\partial_zu^z=0.
\tag{E-02}
\]

**状態: 一次資料で確認済み。**
3次元円柱座標の発散から \(\partial_\theta u^\theta=0\) を入れて導出。[HL09, eq. (2.5)] と [LW, eq. (2.3)] に一致。

監査警告: \(\partial_ru^r+3u^r/r+\partial_zu^z=0\) はこの問題の発散条件ではなく **誤り**。

## 4. 円柱成分と渦度

### E-03: 3成分運動方程式

\[
\begin{aligned}
\partial_tu^r+u^ru^r_r+u^zu^r_z-{(u^\theta)^2\over r}
&=-p_r+\nu(\Delta_0-r^{-2})u^r,\\
\partial_tu^\theta+u^ru^\theta_r+u^zu^\theta_z+{u^ru^\theta\over r}
&=\nu(\Delta_0-r^{-2})u^\theta,\\
\partial_tu^z+u^ru^z_r+u^zu^z_z
&=-p_z+\nu\Delta_0u^z,
\end{aligned}
\tag{E-03}
\]

\[
\Delta_0=\partial_r^2+r^{-1}\partial_r+\partial_z^2.
\]

**状態: 一次資料で確認済み。**
\(\partial_\theta e_r=e_\theta,\partial_\theta e_\theta=-e_r\) を用いて移流項を独立導出した。第2式は [HLW, eq. (2.1a)]、全成分は [Z23, eq. (1.12)] と一致する。遠心項は radial 式の左辺で負、swirl の幾何移流項は左辺で正。

### E-04: curl の成分

\[
\omega^r=-u^\theta_z,\qquad
\omega^\theta=u^r_z-u^z_r,\qquad
\omega^z={1\over r}(ru^\theta)_r.
\tag{E-04}
\]

**状態: 一次資料で確認済み。**
円柱 curl determinant から導出し [HL, eq. (12)] と照合。

### E-05: swirl 速度式

\[
\partial_tu^\theta+u^ru^\theta_r+u^zu^\theta_z
=\nu(\Delta_0-r^{-2})u^\theta-{u^ru^\theta\over r}.
\tag{E-05}
\]

**状態: 一次資料で確認済み。**
E-03 第2式の移項。[HL, eq. (13)]、[HLW, eq. (2.1a)] と一致。

誤り検出条件:

- 最後の符号を正にすると E-03 と不一致。
- 粘性から \(-r^{-2}u^\theta\) を落とすと vector Laplacian でなくなる。

### E-06: swirl 渦度式

\[
\partial_t\omega^\theta+u^r\omega^\theta_r+u^z\omega^\theta_z
=\nu(\Delta_0-r^{-2})\omega^\theta
+{1\over r}\partial_z[(u^\theta)^2]
+{u^r\omega^\theta\over r}.
\tag{E-06}
\]

**状態: 一次資料で確認済み。**
radial 式の \(\partial_z\) から axial 式の \(\partial_r\) を引いて独立導出。遠心項 \(-(u^\theta)^2/r\) を右へ移すため forcing は正。[HL, eq. (14)]、[HLW, eq. (2.1b)] と一致。

誤り検出条件: \(-r^{-1}\partial_z[(u^\theta)^2]\) は遠心項の移項と不一致。

## 5. 流れ関数、速度回復、楕円式

### E-07: 速度回復

\[
u^r=-\psi^\theta_z,\qquad
u^z={1\over r}(r\psi^\theta)_r.
\tag{E-07}
\]

**状態: 一次資料で確認済み。**
\(\boldsymbol u_{\rm mer}=\nabla\times(\psi^\theta e_\theta)\) の定義から導出。[HL, eq. (16)]、[HLW, eq. (2.1d)] と一致。E-02 は恒等的にゼロ。

### E-08: 楕円式

\[
-(\Delta_0-r^{-2})\psi^\theta=\omega^\theta.
\tag{E-08}
\]

**状態: 一次資料で確認済み。**
E-04 の \(\omega^\theta=u^r_z-u^z_r\) に E-07 を代入して導出。[HL, eq. (15)]、[HLW, eq. (2.1c)] と一致。

独立符号検査:

\[
u^r_z-u^z_r
=-\psi^\theta_{zz}
-\partial_r[r^{-1}(r\psi^\theta)_r]
=-(\Delta_0-r^{-2})\psi^\theta.
\]

したがって \(+(\Delta_0-r^{-2})\psi^\theta=\omega^\theta\) は、本書の E-07 と組み合わせる限り **誤り**。

## 6. 正規化変数

### E-09: 定義

\[
u_1={u^\theta\over r},\qquad
\omega_1={\omega^\theta\over r},\qquad
\psi_1={\psi^\theta\over r}.
\tag{E-09}
\]

**状態: 一次資料で確認済み。**
[HLW, eq. (2.2)] と一致。

### E-10: 作用素変換

\[
(\Delta_0-r^{-2})(rf)
=r\left(f_{rr}+{3\over r}f_r+f_{zz}\right)
=r\mathcal L_5f.
\tag{E-10}
\]

**状態: 導出済み。**
\((rf)_{rr}=2f_r+rf_{rr}\) を展開すると \(f/r\) が相殺し \(3f_r\) が残る。全項の次元は \([f]/L\)。

### E-11: \(u_1\) 発展式

\[
\partial_tu_1+u^ru_{1,r}+u^zu_{1,z}
=2\psi_{1,z}u_1+\nu\mathcal L_5u_1.
\tag{E-11}
\]

**状態: 一次資料で確認済み。**
E-05 へ \(u^\theta=ru_1\)、E-10、\(u^r/r=-\psi_{1,z}\) を代入。幾何項が左右から2個生じるため係数は \(2\)。[HLW, eq. (2.3a)] と一致。

項別次元:

\[
[u_{1,t}]=[\psi_{1,z}u_1]=[\nu\mathcal L_5u_1]=T^{-2}.
\]

### E-12: \(\omega_1\) 発展式

\[
\partial_t\omega_1+u^r\omega_{1,r}+u^z\omega_{1,z}
=\partial_z(u_1^2)+\nu\mathcal L_5\omega_1.
\tag{E-12}
\]

**状態: 一次資料で確認済み。**
E-06 へ \(\omega^\theta=r\omega_1,u^\theta=ru_1\) を代入。material derivative から出る \(u^r\omega_1\) と stretching \(u^r\omega_1\) が相殺。[HLW, eq. (2.3b)] と一致。

項別次元は \(L^{-1}T^{-2}\)。

### E-13: \(\psi_1\) 楕円式

\[
-\mathcal L_5\psi_1=\omega_1.
\tag{E-13}
\]

**状態: 一次資料で確認済み。**
E-08 と E-10 を \(r\) で割って導出。[HLW, eq. (2.3c)] と一致。両辺の次元は \(L^{-1}T^{-1}\)。

### E-14: 正規化流れ関数からの速度回復

\[
u^r=-r\psi_{1,z},\qquad
u^z=2\psi_1+r\psi_{1,r}.
\tag{E-14}
\]

**状態: 一次資料で確認済み。**
E-07 へ \(\psi^\theta=r\psi_1\) を代入。[HLW, eq. (2.3d)] と一致。

### E-15: 回復速度の3D発散

\[
\partial_ru^r+{u^r\over r}+\partial_zu^z
=(-\psi_{1,z}-r\psi_{1,rz})-\psi_{1,z}
+(2\psi_{1,z}+r\psi_{1,rz})=0.
\tag{E-15}
\]

**状態: 導出済み。**
ここで使うのは E-02 の3次元発散であり、\(\mathcal L_5\) ではない。

## 7. 軸条件

### E-16: 偶奇・極条件

符号付き \(r\) 延長で

\[
u^r,u^\theta,\omega^\theta,\psi^\theta\text{ は奇},\qquad
u^z,p,u_1,\omega_1,\psi_1\text{ は偶}.
\tag{E-16a}
\]

\[
u^r(0,z)=u^\theta(0,z)=\omega^\theta(0,z)=\psi^\theta(0,z)=0,
\tag{E-16b}
\]

\[
\partial_r^{2k+1}u_1(0,z)
=\partial_r^{2k+1}\omega_1(0,z)
=\partial_r^{2k+1}\psi_1(0,z)=0.
\tag{E-16c}
\]

**状態: 一次資料で確認済み。**
[LW, Corollary 1 and Lemma 2] は smooth axisymmetric vector field の axial 成分に奇数階、radial/swirl 成分と角流れ関数に偶数階の消失条件を与える（[LW] の軸方向変数 \(x\) を本書の \(z\) と読み替える）。[HL, eqs. (18)–(21)] も \(u^\theta,\omega^\theta,\psi^\theta\) の奇 Taylor 展開を明記する。

監査注意:

- 値だけの \(u^\theta=\omega^\theta=\psi^\theta=0\) は高階正則性の十分条件ではない。
- \(u_1(0,z)=0\) は一般には要求されない。要求される最低条件は \(u_{1,r}(0,z)=0\)。
- \(z\) 方向の偶奇は別途課した対称性であり、軸対称性からは出ない。

### E-17: 軸上の \(\mathcal L_5\)

\[
(\mathcal L_5f)(0,z)=4f_{rr}(0,z)+f_{zz}(0,z)
\quad\text{for even smooth }f.
\tag{E-17}
\]

**状態: 導出済み。**
\(f=f_0+\frac12f_{rr}(0)r^2+O(r^4)\) を代入。軸で \(3f_r/r\) を浮動小数点除算してはならない。

## 8. 3次元への復元

### E-18: Cartesian 速度と全渦度

\[
u_x=-x\psi_{1,z}-yu_1,\qquad
u_y=-y\psi_{1,z}+xu_1,\qquad
u_z=2\psi_1+r\psi_{1,r},
\tag{E-18a}
\]

\[
\omega^r=-ru_{1,z},\qquad
\omega^\theta=r\omega_1,\qquad
\omega^z=2u_1+ru_{1,r}.
\tag{E-18b}
\]

**状態: 導出済み。**
E-14、\(u^\theta=ru_1\)、円柱基底、および E-04 から得る。E-16 により軸まで滑らかに延長できる。

楕円式を使う独立 curl 検査:

\[
\partial_zu^r-\partial_ru^z
=-r\mathcal L_5\psi_1=r\omega_1.
\]

### E-19: 圧力回復

\[
-\Delta p=\partial_i\partial_j(u_i u_j).
\tag{E-19}
\]

**状態: 導出済み。**
E-01 の発散を取り E-01 の発散ゼロを使用。全空間では減衰条件と加法定数の規格化が必要。

### E-24: 原始変数系との同値性

十分な関数空間、無限遠条件、E-16 の極条件を満たす vorticity–stream 解は E-18 と E-19 を介して3次元原始変数解へ戻る。

**状態: 一次資料で確認済み。**
[LW] は軸対称 solenoidal vector field の vorticity–stream 形式と primitive 形式の同値性を、極条件を含めて扱う。極条件を外した半平面上の形式解については、この同値性を主張してはならない。

### E-25: 有限円柱Poisson対照問題

独立楕円solverの数値対照では、E-13そのものに対して

\[
0\le r\le R,\qquad z\in[0,L_z)\ \text{periodic},\qquad
\psi_1(R,z)=g(z),\qquad \partial_r\psi_1(0,z)=0
\tag{E-25}
\]

を課す。\(g\) は明示的な滑らかなDirichlet traceである。

**状態: 導出済み。**
PDEの符号はE-13、軸条件はE-16から得る。周期条件と外側Dirichlet条件は
manufactured controlのために宣言した有限領域境界であり、元の
\(\mathbb R^3\) 問題から導かれた無限遠条件ではない。従ってE-25の解を
全空間解へ同一視すること、または外側境界誤差が制御済みとすることは
**誤り**。

### E-26: 独立Poisson solverの離散式

\[
\mathcal L_{5,r}\psi
=r^{-3}\partial_r(r^3\partial_r\psi)
\]

について、\(r_i=i\Delta r\) と
\[
V_i=\int_{r_{i-1/2}}^{r_{i+1/2}}r^3\,dr
=\frac{r_{i+1/2}^4-r_{i-1/2}^4}{4}
\]
を用い、
\[
(\mathcal L_{5,r}^{\,h}\psi)_i=
\frac{
r_{i+1/2}^3(\psi_{i+1}-\psi_i)/\Delta r
-r_{i-1/2}^3(\psi_i-\psi_{i-1})/\Delta r
}{V_i}.
\tag{E-26a}
\]

軸cellでは \(r_{-1/2}=0\)、\(r_{1/2}=\Delta r/2\) なので
\[
(\mathcal L_{5,r}^{\,h}\psi)_0
=\frac{8(\psi_1-\psi_0)}{\Delta r^2}.
\tag{E-26b}
\]

\(z\) のFourier波数 \(k_m=2\pi m/L_z\) ごとに、未知radial行の
\(-\mathcal L_5\) matrixは
\[
-a_i^-\widehat\psi_{i-1}
+(a_i^-+a_i^++k_m^2)\widehat\psi_i
-a_i^+\widehat\psi_{i+1}
=\widehat\omega_i,
\tag{E-26c}
\]
\[
a_i^\pm=\frac{r_{i\pm1/2}^3}{\Delta r\,V_i}.
\]
最後の未知行では既知の外側値を右辺へ
\(+a_i^+\widehat g_m\) として移す。

**状態: 導出済み。**
E-26aはfluxのcontrol-volume積分、E-26bはその軸極限、E-26cは
\(-\partial_{zz}\mapsto+k_m^2\) から得た。実装はRHSの節点値をcell平均の
近似として使用するため、滑らかな偶関数に対する二次整合性は
manufactured refinementで数値検査できるが、厳密cell積分ではない。
また、この非対称な座標基底matrixの通常の条件数は、\(r^3dr\) 重み付き
coercivity定数ではない。

## 8b. Hou 有限円柱設定(E-27–E-31)

以下は [Hou21](arXiv:2107.06509v2、LaTeX 原文で監査)の有限円柱計算を
再現するための境界条件・対称性・初期値・数値プロトコルの監査である。
詳細な一次資料監査は `docs/hou_setup_audit.md` に記録した。
[Hou21, eqs. (2.1a)–(2.1d)] は E-11–E-14 と符号込みで完全に一致することを
原文照合で確認した(既存の [HLW, eqs. (2.3a–d)] 照合に追加)。

### E-27: 壁条件(\(r=1\)、no-slip/no-flow)

\[
\psi_1(t,1,z)=0,\qquad
u_1(t,1,z)=0,\qquad
\omega_1(t,1,z)=-\psi_{1,rr}(t,1,z).
\tag{E-27}
\]

導出: no-flow \(u^r(1,z)=-1\cdot\psi_{1,z}(1,z)=0\) は
\(\psi_1(1,z)=0\)(z 全体で定数 0)から従う [Hou21, eq. (2.4)]。
no-slip \(u^\theta(1,z)=1\cdot u_1(1,z)=0\) から \(u_1(1,z)=0\)。
no-slip \(u^z(1,z)=2\psi_1(1,z)+\psi_{1,r}(1,z)=0\) と \(\psi_1(1,z)=0\) から
\(\psi_{1,r}(1,z)=0\)。このとき E-13 を壁上で評価すると、壁に沿って
\(\psi_1\equiv0\) なので \(\psi_{1,zz}(1,z)=0\)、\(\psi_{1,r}(1,z)=0\) より
\(-\psi_{1,rr}(1,z)=\omega_1(1,z)\)、すなわち Thom 型渦度境界条件
[Hou21, eq. (2.5)] を得る。

**状態: 一次資料で確認済み。**

監査警告: Poisson solve には Dirichlet 条件 \(\psi_1(1,z)=0\) **のみ**を課す。
\(\psi_{1,r}(1,z)=0\) は渦度境界条件の生成にのみ使う。両方を E-13 の
境界条件として課すと2階楕円問題が過剰決定になり **誤り**。
[Hou21, §2] は "We will enforce the no-slip boundary condition for
\(\omega_1\) as a vorticity boundary condition by discretizing
\(\omega_1(t,1,z)=-\psi_{1,rr}(t,1,z)\) and imposing
\(\psi_{1,r}(t,1,z)=0\)" と明記する。

### E-28: 対称性と半周期領域

[Hou21, §2] の設定: \(z\) 周期 1、\(u_1,\omega_1,\psi_1\) は \(r\) について偶
(E-16 と一致)、\(z\) について奇。周期性と奇対称性から計算領域は半周期
\(D_1=\{(r,z):0\le r\le1,\ 0\le z\le1/2\}\) に縮約でき、境界
\(z=0,1/2\) では \(u_1=\omega_1=\psi_1=0\)、\(u^z=0\)。

**状態: 一次資料で確認済み。**

監査警告: \(z\) 奇対称性は初期値 E-29 の性質が力学で保存されるもので
あり(\(\omega_1\) の奇性は \(\partial_z(u_1^2)\) を通じて動的に誘導される
[Hou21, §2])、軸対称性からは従わない(E-16 の注意と同じ)。フル周期
\(z\in[0,1)\) で計算する実装では、奇対称性の保存を課すのではなく診断として
監視し、破れを故障として扱う。

### E-29: Hou 初期値と導出ノルム

\[
u_1(0,r,z)=\frac{12000\,(1-r^2)^{18}\,\sin(2\pi z)}{1+12.5\,\sin^2(\pi z)},
\qquad \omega_1(0,r,z)=0.
\tag{E-29}
\]

[Hou21, eq. (2.2)] を原文照合で確認(係数 12000、指数 18、分母
\(1+12.5\sin^2(\pi z)\) を含む)。\(\omega_1(0)=0\) と E-27 の同次境界条件
から \(\psi_1(0)\equiv0\)、従って \(u^r(0)=u^z(0)=0\)(導出済み。論文の
"The other two velocity components are set to zero initially" と整合)。

導出値(論文には**記載がない**。E-04/E-18b から独立導出):
\(\omega(0)\) は \(\omega^\theta(0)=0\)、\(\omega^r=-r\,\partial_zu_1\)、
\(\omega^z=2u_1+r\,u_{1,r}\) のみ。最大値は \(z=0\)、\(r=1/\sqrt{37}\) の
\(\omega^r\) 成分で

\[
\|\omega(0)\|_\infty
=24000\pi\cdot37^{-1/2}\cdot(36/37)^{18}
\approx 7569.62,
\qquad
\|u_1(0)\|_\infty\approx3265.9863\ \text{at}\ (r,z)\approx(0,0.0845843).
\tag{E-29b}
\]

**状態: 一次資料で確認済み(E-29)、導出済み(E-29b)。**
E-29b は再現実行の増幅率換算(論文は比のみ記載)に必須であり、実装時に
数値最大化で再検証する。

### E-30: Hou 数値プロトコル

一次資料に記録された再現必須事実:

1. **二段階粘性**: \(\nu=5\times10^{-4}\) を \(t\in[0,t_0]\)、
   \(t_0=0.00227375\)、その後 \(\nu=5\times10^{-3}\) [Hou21, §3]。
   \(\nu=5\times10^{-3}\) 一定ではない。
2. 空間 2 次有限差分(適応写像座標)、時間 2 次陽的 Runge–Kutta
   (Heun、Butcher 表 \(c=(0,1)\), \(a_{21}=1\), \(b=(1/2,1/2)\))
   [Hou21, §1.4; HH21, App. A]。
3. 適応時間刻み: CFL 定数 0.1 の対流・粘性制約の最小
   [HH21, App. A]。
4. \(t_0\) 以後の粗格子 run は 1536\(^2\) 状態からの restart であり、
   独立 run ではない [Hou21, §3.3.2]。
5. 早期再現ターゲット(1536\(^2\)):
   \(t=T_1=0.002191729\)(45,000 steps)で
   \(\|\omega\|_\infty/\|\omega(0)\|_\infty\approx20.5235\)、
   \(t=T_2=0.002261605\)(60,000 steps)で \(\approx139.5777\)、
   \(t_0\) で \(498.42\) [Hou21, App. A.2, §3.2.2]。
   ごく早期には \(\|u_1\|_\infty\) は**減少**する [Hou21, §2]。

**状態: 一次資料で確認済み(プロトコルの記録として)。**
これは Hou の手法記述の監査であって、本リポジトリ実装の正当性や
"potentially singular" 主張の検証ではない。Poisson solver
(B-spline Galerkin)、filter、適応 mesh 写像は [HH21] に委ねられており、
一様固定格子の本実装はそれらを複製しない(差異は再現報告に明記する)。

### E-31: 壁渦度条件の2次離散式

一様格子 \(r_i=i\Delta r\)、\(r_{n}=1\)(壁)、\(h=\Delta r\) とする。
E-27 の \(\psi_1(1,z)=0\)、\(\psi_{1,r}(1,z)=0\) を用いた壁での Taylor 展開

\[
\psi_{n-1}=\tfrac{h^2}{2}\psi_{1,rr}(1,z)-\tfrac{h^3}{6}\psi_{1,rrr}(1,z)+O(h^4),
\qquad
\psi_{n-2}=2h^2\psi_{1,rr}(1,z)-\tfrac{4h^3}{3}\psi_{1,rrr}(1,z)+O(h^4)
\]

から \(\psi_{1,rrr}\) を消去して

\[
\omega_1(t,1,z)
=-\psi_{1,rr}(1,z)
=-\frac{8\psi_{n-1}-\psi_{n-2}}{2h^2}+O(h^2).
\tag{E-31}
\]

1次の Thom 式は \(\omega_1(1,z)=-2\psi_{n-1}/h^2+O(h)\)。

**状態: 導出済み。**
[Hou21] は離散式そのものを印刷しておらず([HH21] は wall 行で
\((\omega_1)_{n,j}=-(\psi_{1,rr})_{n,j}\) と外挿
\(v_{n+1,j}=3v_{n,j}-3v_{n-1,j}+v_{n-2,j}\) を記載)、E-31 の具体的
離散化は本リポジトリの選択である。manufactured 収束テストで
2次を確認しなければ実装受入としない。

### E-32: 壁依存性試験用の \(C^\infty\) compact-support 初期値族

壁半径 \(R_{\mathrm{wall}}\) を変える実験では、E-29 の式を \(r>1\) へ
形式的に延長してはならない(\((1-r^2)^{18}\) は \(r>1\) で増大する別の関数に
なる)。単純なゼロ延長も \(r=1\) で \(C^{17}\) 止まりである。代わりに

\[
\tilde u_1(0,r,z)=u_1^{E29}(0,r,z)\,\chi_c(r^2),
\qquad
\tilde\omega_1(0,r,z)=0,
\tag{E-32a}
\]

\[
\chi_c(\rho)=
\begin{cases}
1,&\rho\le\rho_1,\\[2pt]
\dfrac{\theta\!\left(\frac{\rho_2-\rho}{\rho_2-\rho_1}\right)}
      {\theta\!\left(\frac{\rho_2-\rho}{\rho_2-\rho_1}\right)
      +\theta\!\left(\frac{\rho-\rho_1}{\rho_2-\rho_1}\right)},
&\rho_1<\rho<\rho_2,\\[10pt]
0,&\rho\ge\rho_2,
\end{cases}
\qquad
\theta(s)=\begin{cases}e^{-1/s},&s>0\\0,&s\le0\end{cases}
\tag{E-32b}
\]

を用いる。既定パラメータは \(\rho_1=0.81\)(\(r=0.9\))、
\(\rho_2=0.9025\)(\(r=0.95\))。

**性質(すべて検算済み)。**

1. **\(C^\infty\)**: \(\theta\) は標準的な平坦関数で \(C^\infty\)、
   分母は \(\rho\in[\rho_1,\rho_2]\) で正(両引数のいずれかが正)なので
   \(\chi_c\in C^\infty(\mathbb R)\)。\(\chi_c(\rho_1)=1\)、
   \(\chi_c(\rho_2)=0\)(端点で \(\theta(0)=0\) から直接)。
2. **軸互換・偶対称**: \(\chi_c\) を \(\rho=r^2\) の関数として定義するため
   \(\tilde u_1\) は \(r\) について偶。E-16 の軸条件は E-29 と同じく満たす。
3. **core 同一性**: \(r\le0.9\) で \(\chi_c\equiv1\)(分岐による厳密値)。
   浮動小数点でも \(1.0\) 倍は厳密なので、**core は E-29 と bit 一致**する。
   全壁半径で core 初期値は同一である。
4. **偏差上界**: \(0\le\chi_c\le1\) より
   \(\sup_{r,z}|\tilde u_1-u_1^{E29}|
   \le 12000\,(1-\rho_1)^{18}\max_z|g(z)|\)、
   ここで \(g(z)=\sin(2\pi z)/(1+12.5\sin^2(\pi z))\)。
   数値検算: \(\max|g|=0.272165526975908\)(\(z=0.0845842\))、
   \(12000\max|g|=3265.986323710896\)(E-29b と一致)、
   \((1-0.81)^{18}=1.0412735\times10^{-13}\) より
   \[
   \sup|\tilde u_1-u_1^{E29}|\le3.4008\times10^{-10},
   \qquad\text{相対}\ 1.0413\times10^{-13}.
   \tag{E-32c}
   \]
   実格子(\(n_r=193,385\))での実測 sup 偏差は \(2.578\times10^{-12}\)。
5. **compact support**: \(r\ge0.95\) で厳密に 0。したがって
   \(R_{\mathrm{wall}}\ge1\) のすべての壁で、初期値は壁より手前で消える。
6. **有限エネルギー**: support が有界、\(z\) 周期有限長なので E-20 測度で
   有限。
7. **高階微分の摂動**: 遷移帯は E-29 の振幅が \(O(10^{-13})\) 相対の領域
   にあるため、離散高階差分への影響は測定限界以下である。実測(\(z\) を
   \(g\) の最大点に固定、遷移帯 \(0.85<r<1.0\)):
   4 階中心差分の帯内最大は E-29 単体と envelope で**同一**
   (\(n_r=193\): 21.739、\(n_r=385\): 26.362)。全域最大
   \(\approx1.19\times10^7\) は軸近傍で達成され、envelope とは無関係。

**状態: 導出済み。**
\(\chi_c\) は標準的な \(C^\infty\) partition-of-unity 構成であり、E-29 の
既知性質と組み合わせて上記 1–7 を導いた。これは E-29 の再現用初期値では
なく、**壁依存性を測るための別の初期値族**である。この族での成長が E-29
の成長を再現する保証はなく、core 同一性(性質 3)と偏差上界(E-32c)が
その比較を意味あるものにする根拠である。

監査警告: \(\rho_1,\rho_2\) を core 半径へ近づけると性質 3・4 の保証が
失われる。パラメータを変更する場合は E-32c を再計算し、
`docs/wall_dependence_prereg.md` の受入検査を再実行しなければならない。

## 9. エネルギー、次元、スケーリング

### E-20: 物理エネルギー

\[
E(t)={1\over2}\int_{\mathbb R^3}|u|^2dx
=\pi\int_{\mathbb R}\int_0^\infty
\big[(u^r)^2+(u^\theta)^2+(u^z)^2\big]r\,dr\,dz,
\tag{E-20a}
\]

\[
E(t)+\nu\int_0^t\int_{\mathbb R^3}|\nabla u|^2dx\,ds=E(0).
\tag{E-20b}
\]

**状態: 一次資料で確認済み。**
有限エネルギー要件は [F]、エネルギー等式は [HLW, eq. (1.4a)]。円柱積分は \(d^3x=r\,d\theta\,dr\,dz\) から導出。

監査警告: \(r^3dr\,dz\) は形式的5次元測度であり E-20 の代用にはならない。

### E-21: スケーリング

\[
u^{(\lambda)}(x,t)=\lambda u(\lambda x,\lambda^2t),\qquad
p^{(\lambda)}(x,t)=\lambda^2p(\lambda x,\lambda^2t),
\tag{E-21a}
\]

\[
u_1^{(\lambda)}=\lambda^2u_1(\lambda r,\lambda z,\lambda^2t),\quad
\omega_1^{(\lambda)}=\lambda^3\omega_1(\lambda r,\lambda z,\lambda^2t),\quad
\psi_1^{(\lambda)}=\lambda\psi_1(\lambda r,\lambda z,\lambda^2t).
\tag{E-21b}
\]

**状態: 一次資料で確認済み。**
E-01 と E-11–E-14 へ代入して導出。[HLW, eqs. (1.3), (2.5)] は同じ変換を \(\tau=\lambda^{-2}\) で記載。

### E-22: 次元

\[
[u]=LT^{-1},\ [p]=L^2T^{-2},\ [\nu]=L^2T^{-1},\
[\omega]=T^{-1},\ [\psi^\theta]=L^2T^{-1},
\tag{E-22a}
\]

\[
[u_1]=T^{-1},\quad
[\omega_1]=L^{-1}T^{-1},\quad
[\psi_1]=LT^{-1}.
\tag{E-22b}
\]

**状態: 導出済み。**
E-07 と E-09 から得る。E-11–E-13 の全項で一致を確認した。

### E-23: 循環

\[
\Gamma=ru^\theta=r^2u_1,
\tag{E-23a}
\]

\[
\partial_t\Gamma+u^r\Gamma_r+u^z\Gamma_z
=\nu\left(\partial_r^2-{1\over r}\partial_r+\partial_z^2\right)\Gamma.
\tag{E-23b}
\]

**状態: 導出済み。**
E-05 に \(u^\theta=\Gamma/r\) を代入。軸正則性から \(\Gamma=O(r^2)\)。候補の独立診断として有用。

## 10. 独立な符号・係数チェック一覧

| チェック | 計算 | 期待結果 |
|---|---|---|
| 流れ関数の符号 | \(u^r_z-u^z_r\) に E-07 | \(-(\Delta_0-r^{-2})\psi^\theta\) |
| swirl forcing の符号 | radial momentum を \(z\) 微分 | \(+r^{-1}\partial_z[(u^\theta)^2]\) |
| \(u_1\) stretching 係数 | \(D_t(ru_1)\) と E-05 | \(-2u^r u_1/r=+2\psi_{1,z}u_1\) |
| \(3/r\) 粘性係数 | \((\Delta_0-r^{-2})(rf)\) 展開 | \(r(f_{rr}+3f_r/r+f_{zz})\) |
| \(\omega_1\) stretching 相殺 | \(D_t(r\omega_1)\) と E-06 | \(u^r\omega_1\) が左右で相殺 |
| 3D発散 | E-14 を E-02 へ代入 | 恒等的に0 |
| 軸作用素 | 偶 Taylor 展開 | \(4f_{rr}+f_{zz}\) |
| 次元 | E-11–E-13 の各項 | 方程式内で一致 |

## 11. 資料上の注意点と採用しない式

### 11.1 Hou–Lei (2009) の表示上の \(u_1\) 粘性係数

[HL09, printed eq. (2.9)] のオンライン抽出では、\(u_1\) の粘性項だけが

\[
\nu(\partial_r^2+r^{-1}\partial_r+\partial_z^2)u_1
\]

と表示される版がある。しかし同論文が直後に提示する convection-free model の \(u_1\) 粘性項は \(\nu\mathcal L_5u_1\) であり、E-10 の直接恒等式、Hou–Li の元式、後続の査読論文 [HLW, eq. (2.3a)] はすべて \(3/r\) を与える。

上の \(1/r\) 表示の状態: **誤り**。
採用式 E-11 の状態: **一次資料で確認済み**。

### 11.2 5次元作用素の誤読

\[
\mathcal L_5=\partial_r^2+3r^{-1}\partial_r+\partial_z^2
\]

を5次元 radial scalar Laplacian として利用する代数的対応自体は正しい。一方、

\[
\partial_ru^r+3u^r/r+\partial_zu^z=0
\]

を物理非圧縮条件とすることの状態: **誤り**。
物理条件は E-02。

### 11.3 1D Hou–Li model との区別

[HL, eqs. (22)–(24)] は軸近傍展開と追加の異方性仮定から radial coupling を落とした1D model であり、一般の3次元軸対称 Navier–Stokes と同じではない。本実装の規範は E-11–E-14 の2変数 \((r,z)\) 系である。1D model を完全な3D方程式として用いることの状態: **誤り**。

### 11.4 Nowakowski–Zajączkowski (2023) の \(\omega^\theta\) 式

[Z23, printed eq. (1.9)] は、同論文自身の
\(\omega^\theta=u^r_z-u^z_r\) という規約にもかかわらず、swirl forcing と
\(u^r\omega^\theta/r\) を左辺で正符号に表示している。この表示をそのまま右辺形へ移すと E-06 と逆符号になる。一方、同論文の変換後の式 [Z23, eq. (1.15)] は
\[
\omega_{1,t}+\boldsymbol v\cdot\nabla\omega_1
-\nu\mathcal L_5\omega_1=2u_1u_{1,z}+F_1
\]
であり E-12 と一致する。E-06 の符号は独立な交差微分と [HL]・[HLW] で確定している。

[Z23, eq. (1.9)] の二つの正符号表示の状態: **誤り**。
[Z23, eqs. (1.12), (1.15)] のうち本書が照合に用いた成分式・変換後式の状態: **一次資料で確認済み**。

## 12. 実装受入条件

実装は次を満たさなければ監査合格としない。

1. E-11 と E-12 の粘性に同じ \(\mathcal L_5\) を用いる。
2. E-13 の負符号を保存する。
3. 速度回復は E-14、物理発散検査は E-02 を用いる。
4. 軸では E-16 と E-17 を用い、\(1/r\) を直接評価しない。
5. E-13 残差に加え、E-18 から独立に E-04 の \(\omega^\theta\) を再計算する。
6. 物理エネルギーは E-20 の \(r\,dr\,dz\) 測度で計算する。
7. 候補保存物から E-18 の3次元 Cartesian 場を再構成できる。
8. 符号反転した楕円式、壊した E-02、壊した E-16 を必ずテストで拒否する。
9. 既存の円柱差分を共有しない一様 Cartesian \(x,y,z\) 実装で、E-01 の
   3成分発散、full curl、vector Laplacian、primitive momentum residualを
   再計算し、局所的な故障もRMSと最大誤差の双方で拒否する。
10. 有限円柱Poisson solveはE-25–E-26を明示し、E-13の全体符号、軸係数8、
    非零外側trace、複数Fourier mode、代数残差と独立物理空間残差を分離して
    少なくとも3解像度で検査する。

**現在の実装状態。** 条件 1–4, 6, 8 は production finite-difference
経路と故障注入で実装・テスト済み。条件 7 は保存候補をchecksum検証後に
読み、任意の \(\theta\) 標本上で E-18 の
\((x,y,z,u_x,u_y,u_z)\) を生成するreaderとして実装した。条件 5 は、その
Cartesian成分を円柱半径方向へ再投影し、production \(L_5\) を呼ばない独立な
4次 \(r,z\) stencilで E-04 を評価して \(r\omega_1\) と比較する。滑らかな
manufactured fieldで defect は解像度とともに4次で減少し、誤った
\(\omega_1\) 符号を拒否する。

条件 9 は `cartesian_validation.py` の一様 \(x,y,z\) 格子と、このmoduleだけが
所有する2次差分で実装した。ここで計算する

\[
\nabla\cdot u,\qquad \nabla\times u,\qquad
u_t+(u\cdot\nabla)u+\nabla p-\nu\Delta u-f
\]

は既存 `operators.py`、`pde.py`、円柱差分結果を呼ばない。周期解析場を
\(12^3,24^3,48^3\) と細分した数値試験では、各微分、各運動量項、全残差の
誤差が約2次で減少した。非周期片側closureも二次多項式の閉形式微分で
全境界を含めて検査済みである。

artifact-level試験では、candidateを保存しchecksum検証付きで再読込した後、
`cartesian_candidate_adapter.py` 固有の2次 \(r,z\) stencilとbilinear補間から
E-18a, E-18bを一様Cartesian格子へ写す。checkerはCartesian配列だけを受け取り、
物理3次元発散とfull curlを再計算する。adapterの出力は別途、閉形式
manufactured oracleへ直接照合し、円柱radial符号、成分写像、渦度符号、
発散汚染、局所一点、周期seamの故障を拒否する。

条件10は、既存 `operators.py` をimportしない `poisson.py` の
\(r^3\)-flux/Fourier/Thomas経路で実装した。test側の直接
非発散形stencil、閉形式manufactured pair、独立dense matrix solveを用いて、
符号、軸係数、非零境界、FFT mode、三重対角解を交差検査した。

従って、ここに列挙した条件 1–10 の初期受入状態は **導出済み・数値テスト
済み** である。ただし、これは有限格子・binary64・許容差付きの独立監査で
あり、E-01 の連続体解、全空間の境界条件、離散化誤差上界、または特異点を
証明しない。将来の未知候補には、候補固有の圧力・時間微分・境界・精度検査を
改めて適用しなければならない。

## 13. 一次資料

- **[F]** Charles L. Fefferman, “Existence and Smoothness of the Navier–Stokes Equation,” Clay Mathematics Institute.
  <https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf>
- **[HL]** Thomas Y. Hou and Congming Li, “Dynamic Stability of the 3D Axi-symmetric Navier–Stokes Equations with Swirl,” *Communications on Pure and Applied Mathematics* 61 (2008), 661–697.
  DOI: <https://doi.org/10.1002/cpa.20212>
  arXiv: <https://arxiv.org/abs/math/0608295>
- **[HL09]** Zhen Lei and Thomas Y. Hou, “On the Stabilizing Effect of Convection in Three-Dimensional Incompressible Flows,” *Communications on Pure and Applied Mathematics* 62 (2009), 501–564.
  DOI: <https://doi.org/10.1002/cpa.20254>
  著者 PDF: <https://users.cms.caltech.edu/~hou/papers/3DModelNS-final.pdf>
- **[HLW]** Thomas Y. Hou, Pengfei Liu, and Fei Wang, “Global regularity for a family of 3D models of the axi-symmetric Navier–Stokes equations,” *Nonlinearity* 31 (2018), 1940–1954.
  DOI: <https://doi.org/10.1088/1361-6544/aaaa0b>
  arXiv: <https://arxiv.org/abs/1708.07536>
- **[LW]** Jian-Guo Liu and Wei-Cheng Wang, “Characterization and Regularity for Axisymmetric Solenoidal Vector Fields with Application to Navier–Stokes Equation,” *SIAM Journal on Mathematical Analysis* 41 (2009), 1825–1850.
  DOI: <https://doi.org/10.1137/080739744>
  PDF: <https://archive.ymsc.tsinghua.edu.cn/pacm_download/200/8347-Liu_Wang_SIMA_2009.pdf>
- **[Z23]** Bernard Nowakowski and Wojciech M. Zajączkowski, “Global Regular Axially-Symmetric Solutions to the Navier–Stokes Equations with Small Swirl,” *Journal of Mathematical Fluid Mechanics* 25 (2023), article 73.
  DOI: <https://doi.org/10.1007/s00021-023-00793-9>
- **[Hou21]** Thomas Y. Hou, “Potentially Singular Behavior of the 3D Navier–Stokes Equations,” *Foundations of Computational Mathematics* (2022), DOI 10.1007/s10208-022-09578-4(= arXiv:2107.06509v2)。
  arXiv: <https://arxiv.org/abs/2107.06509>
  監査は v1(2021-07-14)と v2(2022-05-26)の LaTeX 原文
  (`arxiv.org/e-print/2107.06509`)、および**出版版 PDF**
  (SHA-256 `b8ad5ed5...e765acd`、`docs/hou_setup_audit.md` §12)に対して
  行った。出版版で計算設定の変更はなく、方程式・境界条件・初期値・粘性
  プロトコルの引用は arXiv v2 の式番号のまま有効(付録 (A.1)→(4.1) のみ
  改番)。v1 と v2 で結論の解釈が反転している(同 §9)。
- **[HH21]** Thomas Y. Hou and De Huang, “Potential Singularity of the 3D Euler Equations in the Interior Domain,” arXiv:2102.06663。
  arXiv: <https://arxiv.org/abs/2102.06663>
  [Hou21, §2] が数値手法の詳細(B-spline Galerkin Poisson solver、
  CFL 定数、mesh 密度関数、filter)を本論文 Appendix A/B に委ねている。
