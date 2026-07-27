# 研究仕様: 3次元非圧縮 Navier–Stokes 特異点候補と検証証明書

## 1. 目的と主張の境界

本プロジェクトは、3次元非圧縮 Navier–Stokes 方程式の有限時間特異点**候補**を再現可能に探索し、その候補を将来の区間演算・コンピューター支援証明へ渡せる明示的データにするための研究基盤である。

次は本プロジェクトの成果とはみなさない。

- 大きな渦度、高い勾配、または急速なノルム増大だけを有限時間特異点と呼ぶこと。
- 離散残差またはニューラルネットの訓練残差が小さいことだけを PDE 解の存在証明と呼ぶこと。
- 有限領域、有限解像度、浮動小数点計算の発散を数学的発散と同一視すること。
- 変換後の形式的な5次元スカラー Laplacian を、5次元流体または5次元の非圧縮条件とみなすこと。

厳密な存在・誤差評価・発散評価がすべて完成するまでは、出力の状態は常に「数値的候補」または「未確認」である。

## 2. 数学的対象

### 2.1 物理領域と方程式

主対象は、外力なし、密度で割った圧力を用いる全空間問題

\[
\partial_t\boldsymbol u+(\boldsymbol u\cdot\nabla)\boldsymbol u
=-\nabla p+\nu\Delta\boldsymbol u,\qquad
\nabla\cdot\boldsymbol u=0
\quad\text{on }\mathbb R^3\times(0,T),
\tag{S-NS}
\]

\[
\boldsymbol u(\boldsymbol x,0)=\boldsymbol u_0(\boldsymbol x),\qquad \nu>0
\tag{S-IC}
\]

である。数学上の領域は \(\mathbb R^3\) であり、数値計算の有限長方形領域は単なる打ち切り近似である。

初期データは

\[
\boldsymbol u_0\in C_c^\infty(\mathbb R^3;\mathbb R^3),\qquad
\nabla\cdot\boldsymbol u_0=0,
\tag{S-DATA}
\]

を既定とする。より広い Schwartz 級を許す実験は可能だが、データ級をメタデータに明記しなければならない。探索対象では軸対称かつ旋回あり、すなわち \(u_0^\theta\not\equiv0\) とする。

### 2.2 解と有限時間特異点

この仕様でいう滑らかな解は、各 \(T'<T\) について

\[
\boldsymbol u,p\in C^\infty(\mathbb R^3\times[0,T'])
\]

で (S-NS)–(S-IC) を各点で満たし、各時刻で有限エネルギーを持ち、空間無限遠で圧力・流れ関数を一意に回復できるだけの減衰を持つ解である。同値な局所強解の枠組みとして、整数 \(m>5/2\) に対する

\[
\boldsymbol u\in C([0,T);H^m)\cap
C^1([0,T);H^{m-2})
\]

を用いてよい。

最大存在時刻 \(T_*\) が有限で、その強解を同じクラスで \(T_*\) より先へ延長できないときにのみ有限時間特異点という。数値時系列から推定した \(T_*\)、有限個の点での増大、または回帰したべき指数は、この定義を満たす証明ではない。

## 3. 座標、符号、変数

\[
r=\sqrt{x^2+y^2},\quad
\boldsymbol e_r=(\cos\theta,\sin\theta,0),\quad
\boldsymbol e_\theta=(-\sin\theta,\cos\theta,0),\quad
\boldsymbol e_z=(0,0,1).
\]

軸対称速度場を

\[
\boldsymbol u
=u^r(r,z,t)\boldsymbol e_r
+u^\theta(r,z,t)\boldsymbol e_\theta
+u^z(r,z,t)\boldsymbol e_z,\qquad \partial_\theta(\cdot)=0
\tag{S-AXI}
\]

と書く。旋回ありとは \(u^\theta\not\equiv0\) である。渦度の符号規約は

\[
\boldsymbol\omega=\nabla\times\boldsymbol u,\qquad
\omega^\theta=\partial_z u^r-\partial_r u^z
\tag{S-CURL}
\]

で固定する。

角方向ベクトルポテンシャルを

\[
\boldsymbol A=\psi^\theta\boldsymbol e_\theta,\qquad
(u^r\boldsymbol e_r+u^z\boldsymbol e_z)=\nabla\times\boldsymbol A
\]

と定義する。この規約により

\[
u^r=-\partial_z\psi^\theta,\qquad
u^z={1\over r}\partial_r(r\psi^\theta),\qquad
-\left(\partial_r^2+{1\over r}\partial_r+\partial_z^2-{1\over r^2}\right)
\psi^\theta=\omega^\theta.
\tag{S-BS}
\]

逆符号の流れ関数規約を同じ候補ファイル内で混用してはならない。

## 4. 正規化変数と規範系

\[
u_1={u^\theta\over r},\qquad
\omega_1={\omega^\theta\over r},\qquad
\psi_1={\psi^\theta\over r},
\tag{S-VAR}
\]

および

\[
\mathcal L_5
:=\partial_r^2+{3\over r}\partial_r+\partial_z^2
\tag{S-L5}
\]

を用いる。実装が従うべき閉じた系は

\[
\partial_tu_1+u^r\partial_ru_1+u^z\partial_zu_1
=2(\partial_z\psi_1)u_1+\nu\mathcal L_5u_1,
\tag{S-U1}
\]

\[
\partial_t\omega_1+u^r\partial_r\omega_1+u^z\partial_z\omega_1
=\partial_z(u_1^2)+\nu\mathcal L_5\omega_1,
\tag{S-W1}
\]

\[
-\mathcal L_5\psi_1=\omega_1,
\tag{S-P1}
\]

\[
u^r=-r\partial_z\psi_1,\qquad
u^z=2\psi_1+r\partial_r\psi_1,\qquad
u^\theta=ru_1.
\tag{S-REC}
\]

である。これらの符号と \(3/r\) 係数は [数学的定式化](docs/mathematical_formulation.md) で3次元式から導出し、[方程式監査](docs/equation_audit.md) で一次資料と照合している。

\(\mathcal L_5\) は、\(r\) を \(\mathbb R^4\) の半径と形式的に読めば5次元の軸対称**スカラー** Laplacian と同じ式になる。しかし物理速度は3次元であり、非圧縮条件は常に

\[
\partial_ru^r+{u^r\over r}+\partial_zu^z=0.
\tag{S-DIV3}
\]

物理体積要素も \(2\pi r\,dr\,dz\) であって \(r^3dr\,dz\) ではない。

## 5. 軸上の正則性

\(r\) を符号付き変数として延長したとき、Cartesian 空間で滑らかな軸対称場には次を要求する。

- \(u^r,u^\theta,\omega^\theta,\psi^\theta\) は \(r\) の奇関数。
- \(u^z,p,u_1,\omega_1,\psi_1\) は \(r\) の偶関数。
- したがって \(u^r=u^\theta=\omega^\theta=\psi^\theta=0\) at \(r=0\)。
- \(\partial_r^{2k+1}u_1=\partial_r^{2k+1}\omega_1
  =\partial_r^{2k+1}\psi_1=0\) at \(r=0\)（存在するすべての次数）。

特に Neumann 条件 \(\partial_r u_1=\partial_r\omega_1=\partial_r\psi_1=0\) は最低次数の必要条件にすぎない。\(z\) 方向の偶奇性は軸対称性からは生じない。

偶関数 \(f\) に対する軸上の作用素は特異な式を直接評価せず、

\[
(\mathcal L_5f)(0,z)=4\,\partial_r^2f(0,z)+\partial_z^2f(0,z)
\tag{S-AXOP}
\]

という正則極限で評価する。

## 6. 元の3次元場への復元

\(r>0\) では (S-REC) と円柱基底から復元する。Cartesian 成分では見かけの \(1/r\) を消去して

\[
u_x=-x\,\partial_z\psi_1-yu_1,\qquad
u_y=-y\,\partial_z\psi_1+xu_1,\qquad
u_z=2\psi_1+r\partial_r\psi_1.
\tag{S-CART}
\]

渦度は

\[
\omega^r=-r\partial_zu_1,\qquad
\omega^\theta=r\omega_1,\qquad
\omega^z=2u_1+r\partial_ru_1.
\tag{S-VORTREC}
\]

軸正則性、無限遠条件、楕円関係を満たす \((u_1,\omega_1,\psi_1)\) は、(S-CART) により滑らかな3次元速度へ戻せる。圧力は規格化定数を除き

\[
-\Delta p=\partial_i\partial_j(u_i u_j)
\tag{S-PRESS}
\]

から回復する。閉じた軸対称系から原始変数系への同値性を主張するときは、極条件と関数空間を省略してはならない。

## 7. 次元とスケーリング

基本次元を長さ \(L\)、時間 \(T\) とする。

| 量 | 次元 |
|---|---:|
| \(r,z\) | \(L\) |
| \(t\) | \(T\) |
| \(u^r,u^\theta,u^z\) | \(L\,T^{-1}\) |
| \(p\)（密度で規格化） | \(L^2T^{-2}\) |
| \(\nu\) | \(L^2T^{-1}\) |
| \(\boldsymbol\omega\) | \(T^{-1}\) |
| \(\psi^\theta\) | \(L^2T^{-1}\) |
| \(u_1\) | \(T^{-1}\) |
| \(\omega_1\) | \(L^{-1}T^{-1}\) |
| \(\psi_1\) | \(LT^{-1}\) |

Navier–Stokes の拡大縮小は

\[
\boldsymbol u^{(\lambda)}(\boldsymbol x,t)
=\lambda\boldsymbol u(\lambda\boldsymbol x,\lambda^2t),\qquad
p^{(\lambda)}(\boldsymbol x,t)
=\lambda^2p(\lambda\boldsymbol x,\lambda^2t),
\]

\[
u_1^{(\lambda)}=\lambda^2u_1(\lambda r,\lambda z,\lambda^2t),\quad
\omega_1^{(\lambda)}=\lambda^3\omega_1(\lambda r,\lambda z,\lambda^2t),\quad
\psi_1^{(\lambda)}=\lambda\psi_1(\lambda r,\lambda z,\lambda^2t).
\tag{S-SCALE}
\]

実験が無次元量を保存する場合は、基準長 \(L_0\)、基準速度 \(U_0\)、基準時間 \(L_0/U_0\)、Reynolds 数 \(Re=U_0L_0/\nu\) を必ずメタデータに保存する。

## 8. 候補データに対する必須不変条件

候補を「監査可能」とする最低条件は次である。

1. \(u_1,\omega_1,\psi_1\) を明示的配列または明示的基底係数として保存する。
2. 座標、領域、単位、\(\nu\)、時刻、基底規約、乱数 seed、コード版を保存する。
3. (S-DIV3)、(S-P1)、(S-U1)、(S-W1) の残差を総和だけでなく項別に保存する。
4. 軸の偶奇・極条件を次数別に検査する。
5. 物理測度 \(2\pi r\,dr\,dz\) で有限エネルギーを検査する。
6. 解像度、時間刻み、領域サイズ、精度を独立に変えた収束列を残す。
7. 保存後の再読み込みで係数と診断値が一致することを検査する。
8. 元の3次元 Cartesian 場を (S-CART) から再構成し、独立な Cartesian 微分でも発散と PDE 残差を検査できる形式にする。

これらは必要条件であって特異点の証明ではない。

## 9. 数値候補と証明の区別

| 段階 | 許される結論 |
|---|---|
| 単一解像度の増大 | 高勾配イベントの観測 |
| 複数解像度・領域・精度で収束 | 数値候補の再現 |
| 明示係数と小さな検証残差 | 近似軌道または近似プロファイル |
| 区間演算で作用素・尾部・打切り誤差を包含 | コンピューター支援証明の一部 |
| 初期データから有限物理時刻のノルム発散まで全義務を閉じる | はじめて反例または特異点の証明候補として審査可能 |

## 10. 規範的資料

1. Charles L. Fefferman, “Existence and Smoothness of the Navier–Stokes Equation,” Clay Mathematics Institute official problem description.
   <https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf>
2. Thomas Y. Hou and Congming Li, “Dynamic Stability of the 3D Axi-symmetric Navier–Stokes Equations with Swirl,” *Communications on Pure and Applied Mathematics* 61 (2008), 661–697. DOI: <https://doi.org/10.1002/cpa.20212>; arXiv: <https://arxiv.org/abs/math/0608295>.
3. Thomas Y. Hou, Pengfei Liu, and Fei Wang, “Global regularity for a family of 3D models of the axi-symmetric Navier–Stokes equations,” *Nonlinearity* 31 (2018), 1940–1954. DOI: <https://doi.org/10.1088/1361-6544/aaaa0b>; arXiv: <https://arxiv.org/abs/1708.07536>.
4. Jian-Guo Liu and Wei-Cheng Wang, “Characterization and Regularity for Axisymmetric Solenoidal Vector Fields with Application to Navier–Stokes Equation,” *SIAM Journal on Mathematical Analysis* 41 (2009), 1825–1850. DOI: <https://doi.org/10.1137/080739744>; author-hosted PDF: <https://archive.ymsc.tsinghua.edu.cn/pacm_download/200/8347-Liu_Wang_SIMA_2009.pdf>.
