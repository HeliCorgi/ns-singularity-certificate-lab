# 数学的定式化と独立導出

## 1. 対象と論理上の位置づけ

対象は外力なし、粘性 \(\nu>0\) の3次元非圧縮 Navier–Stokes 初期値問題

\[
\partial_t\boldsymbol u+(\boldsymbol u\cdot\nabla)\boldsymbol u
=-\nabla p+\nu\Delta\boldsymbol u,\qquad
\nabla\cdot\boldsymbol u=0,\qquad
\boldsymbol u|_{t=0}=\boldsymbol u_0
\tag{1.1}
\]

on \(\mathbb R^3\times(0,T)\) である。圧力 \(p\) は密度で割った運動学的圧力である。既定の初期データは

\[
\boldsymbol u_0\in C_c^\infty(\mathbb R^3;\mathbb R^3),\qquad
\nabla\cdot\boldsymbol u_0=0
\tag{1.2}
\]

で、軸対称かつ \(u_0^\theta\not\equiv0\) とする。Clay の公式問題記述は全空間で滑らかかつ急減少する発散ゼロ初期データ、\(\nu>0\)、外力ゼロ、滑らかな有限エネルギー解を明記している [F]。コンパクト台はその急減少条件を満たす、ここでの再現可能な既定選択である。

各 \(T'<T\) 上で \(\boldsymbol u,p\in C^\infty(\mathbb R^3\times[0,T'])\) かつ有限エネルギーであるものを滑らかな解と呼ぶ。強解の定式化では整数 \(m>5/2\) に対する

\[
\boldsymbol u\in C([0,T);H^m)\cap C^1([0,T);H^{m-2})
\tag{1.3}
\]

を用いてよい。最大強解の存在時刻 \(T_*<\infty\) で同じクラスの延長が不可能な場合を有限時間特異点と定義する。数値値の増大、回帰された「発散時刻」、小さい離散残差は、この延長不可能性を証明しない。

## 2. 円柱座標の規約

\[
r=\sqrt{x^2+y^2},\quad
\boldsymbol e_r=(\cos\theta,\sin\theta,0),\quad
\boldsymbol e_\theta=(-\sin\theta,\cos\theta,0),\quad
\boldsymbol e_z=(0,0,1).
\tag{2.1}
\]

したがって

\[
\partial_\theta\boldsymbol e_r=\boldsymbol e_\theta,\qquad
\partial_\theta\boldsymbol e_\theta=-\boldsymbol e_r.
\tag{2.2}
\]

軸対称・旋回ありの速度を

\[
\boldsymbol u
=u^r(r,z,t)\boldsymbol e_r
+u^\theta(r,z,t)\boldsymbol e_\theta
+u^z(r,z,t)\boldsymbol e_z,\qquad
\partial_\theta u^r=\partial_\theta u^\theta=\partial_\theta u^z=0
\tag{2.3}
\]

とする。「旋回あり」は \(u^\theta\not\equiv0\) を意味し、\(\theta\) 依存性があるという意味ではない。

軸対称スカラーに対する3次元 scalar Laplacian を

\[
\Delta_0:=\partial_r^2+{1\over r}\partial_r+\partial_z^2
\tag{2.4}
\]

と書く。

## 3. 物理的な非圧縮条件

円柱座標の3次元発散は

\[
\nabla\cdot\boldsymbol u
={1\over r}\partial_r(ru^r)+{1\over r}\partial_\theta u^\theta+\partial_z u^z.
\]

軸対称性から

\[
\boxed{\partial_r u^r+{u^r\over r}+\partial_z u^z=0}
\quad\Longleftrightarrow\quad
\partial_r(ru^r)+\partial_z(ru^z)=0.
\tag{3.1}
\]

これが物理的な3次元非圧縮条件である。後に現れる \(\partial_r^2+3r^{-1}\partial_r+\partial_z^2\) から5次元の発散条件を作ってはならない。

## 4. 3成分の運動方程式

軸対称係数に対しても基底は \(\theta\) で変化する。方向微分

\[
\boldsymbol u\cdot\nabla
=u^r\partial_r+{u^\theta\over r}\partial_\theta+u^z\partial_z
\]

を (2.2) とともに \(\boldsymbol u\) に作用させると

\[
\begin{aligned}
[(\boldsymbol u\cdot\nabla)\boldsymbol u]^r
&=u^r\partial_ru^r+u^z\partial_zu^r-{(u^\theta)^2\over r},\\
[(\boldsymbol u\cdot\nabla)\boldsymbol u]^\theta
&=u^r\partial_ru^\theta+u^z\partial_zu^\theta+{u^ru^\theta\over r},\\
[(\boldsymbol u\cdot\nabla)\boldsymbol u]^z
&=u^r\partial_ru^z+u^z\partial_zu^z.
\end{aligned}
\tag{4.1}
\]

同様に vector Laplacian は

\[
(\Delta\boldsymbol u)^r=(\Delta_0-r^{-2})u^r,\quad
(\Delta\boldsymbol u)^\theta=(\Delta_0-r^{-2})u^\theta,\quad
(\Delta\boldsymbol u)^z=\Delta_0u^z.
\tag{4.2}
\]

よって

\[
\partial_tu^r+u^r\partial_ru^r+u^z\partial_zu^r-{(u^\theta)^2\over r}
=-\partial_rp+\nu(\Delta_0-r^{-2})u^r,
\tag{4.3}
\]

\[
\partial_tu^\theta+u^r\partial_ru^\theta+u^z\partial_zu^\theta
+{u^ru^\theta\over r}
=\nu(\Delta_0-r^{-2})u^\theta,
\tag{4.4}
\]

\[
\partial_tu^z+u^r\partial_ru^z+u^z\partial_zu^z
=-\partial_zp+\nu\Delta_0u^z.
\tag{4.5}
\]

特に (4.4) を右辺に幾何項を置く形で書けば

\[
\partial_tu^\theta+u^r\partial_ru^\theta+u^z\partial_zu^\theta
=\nu(\Delta_0-r^{-2})u^\theta-{u^ru^\theta\over r}.
\tag{4.6}
\]

この符号は Hou–Li [HL, eq. (13)] および Hou–Liu–Wang [HLW, eq. (2.1a)] と一致する。

## 5. 渦度

符号規約を

\[
\boldsymbol\omega=\nabla\times\boldsymbol u
\tag{5.1}
\]

と固定すると、軸対称場では

\[
\boxed{
\omega^r=-\partial_zu^\theta,\qquad
\omega^\theta=\partial_zu^r-\partial_ru^z,\qquad
\omega^z={1\over r}\partial_r(ru^\theta)
}
\tag{5.2}
\]

である。

### 5.1 \(\omega^\theta\) 方程式の符号検算

(4.3) を \(\partial_z\)、(4.5) を \(\partial_r\) で微分し、前者から後者を引く。圧力の混合微分は消える。遠心項は (4.3) の左辺で \(-(u^\theta)^2/r\) なので、右辺へ移してから \(\partial_z\) を取ると

\[
+\partial_z\!\left({(u^\theta)^2\over r}\right)
=+{1\over r}\partial_z[(u^\theta)^2]
\]

となる。残りの移流微分を (3.1) と
\(\omega^\theta=\partial_zu^r-\partial_ru^z\) で整理すると

\[
\boxed{
\partial_t\omega^\theta+u^r\partial_r\omega^\theta
+u^z\partial_z\omega^\theta
=\nu(\Delta_0-r^{-2})\omega^\theta
+{1\over r}\partial_z[(u^\theta)^2]
+{u^r\omega^\theta\over r}.
}
\tag{5.3}
\]

したがって swirl forcing の符号は正であり、stretching 項も右辺で正である。これは [HL, eq. (14)]、[HLW, eq. (2.1b)] と一致する。

## 6. 流れ関数と楕円式

角方向ベクトルポテンシャル

\[
\boldsymbol A=\psi^\theta(r,z,t)\boldsymbol e_\theta
\tag{6.1}
\]

を取り、meridional velocity を

\[
u^r\boldsymbol e_r+u^z\boldsymbol e_z=\nabla\times\boldsymbol A
\tag{6.2}
\]

と**定義**する。円柱座標の curl を直接計算すると

\[
\boxed{
u^r=-\partial_z\psi^\theta,\qquad
u^z={1\over r}\partial_r(r\psi^\theta).
}
\tag{6.3}
\]

(6.3) を (3.1) に代入すれば

\[
\partial_r(-r\partial_z\psi^\theta)
+\partial_z[\partial_r(r\psi^\theta)]=0
\]

となり、3次元発散ゼロは恒等的に満たされる。

次に (5.2) の \(\theta\) 成分へ (6.3) を代入する。

\[
\begin{aligned}
\omega^\theta
&=\partial_z(-\partial_z\psi^\theta)
-\partial_r\!\left[{1\over r}\partial_r(r\psi^\theta)\right]\\
&=-\left(\partial_r^2+{1\over r}\partial_r+\partial_z^2-{1\over r^2}\right)\psi^\theta.
\end{aligned}
\]

よって符号を含む楕円式は

\[
\boxed{-(\Delta_0-r^{-2})\psi^\theta=\omega^\theta.}
\tag{6.4}
\]

この符号は、流れ関数の定義 (6.2) を逆にすれば同時に反転する。したがって (6.3) と (6.4) は一組として監査する。

## 7. \(u_1,\omega_1,\psi_1\) 系の導出

### 7.1 定義と作用素恒等式

\[
u_1={u^\theta\over r},\qquad
\omega_1={\omega^\theta\over r},\qquad
\psi_1={\psi^\theta\over r},
\tag{7.1}
\]

\[
\mathcal L_5:=\partial_r^2+{3\over r}\partial_r+\partial_z^2.
\tag{7.2}
\]

任意の滑らかな \(f(r,z)\) について

\[
\boxed{(\Delta_0-r^{-2})(rf)=r\mathcal L_5f.}
\tag{7.3}
\]

実際、

\[
\partial_r^2(rf)=2f_r+rf_{rr},\quad
{1\over r}\partial_r(rf)={f\over r}+f_r,\quad
-{rf\over r^2}=-{f\over r},
\]

なので \(f/r\) が消え、\(3f_r+rf_{rr}+rf_{zz}\) が残る。

### 7.2 速度回復

(6.3) と \(\psi^\theta=r\psi_1\) から

\[
\boxed{
u^r=-r\partial_z\psi_1,\qquad
u^z={1\over r}\partial_r(r^2\psi_1)
=2\psi_1+r\partial_r\psi_1.
}
\tag{7.4}
\]

したがって

\[
{u^r\over r}=-\partial_z\psi_1.
\tag{7.5}
\]

### 7.3 \(u_1\) 方程式

\(D_t=\partial_t+u^r\partial_r+u^z\partial_z\) とする。\(u^\theta=ru_1\) なので

\[
D_t(ru_1)=rD_tu_1+u^ru_1.
\tag{7.6}
\]

(4.6)、(7.3) より

\[
rD_tu_1+u^ru_1
=\nu r\mathcal L_5u_1-u^ru_1.
\]

\(r>0\) で割り、(7.5) を使うと

\[
\boxed{
\partial_tu_1+u^r\partial_ru_1+u^z\partial_zu_1
=2(\partial_z\psi_1)u_1+\nu\mathcal L_5u_1.
}
\tag{7.7}
\]

粘性項の係数は \(1/r\) ではなく **\(3/r\)** である。

### 7.4 \(\omega_1\) 方程式

\(\omega^\theta=r\omega_1\) から

\[
D_t(r\omega_1)=rD_t\omega_1+u^r\omega_1.
\tag{7.8}
\]

(5.3) の右辺は

\[
\nu r\mathcal L_5\omega_1
+{1\over r}\partial_z(r^2u_1^2)
+u^r\omega_1
=\nu r\mathcal L_5\omega_1+r\partial_z(u_1^2)+u^r\omega_1.
\]

(7.8) の \(u^r\omega_1\) と打ち消し合うので

\[
\boxed{
\partial_t\omega_1+u^r\partial_r\omega_1+u^z\partial_z\omega_1
=\partial_z(u_1^2)+\nu\mathcal L_5\omega_1.
}
\tag{7.9}
\]

ここで \(\partial_z(u_1^2)=2u_1\partial_zu_1\) である。

### 7.5 \(\psi_1\) 楕円式

(6.4)、(7.3)、\(\psi^\theta=r\psi_1\)、\(\omega^\theta=r\omega_1\) から

\[
\boxed{-\mathcal L_5\psi_1=\omega_1.}
\tag{7.10}
\]

(7.7)、(7.9)、(7.10)、(7.4) は [HLW, eqs. (2.3a)–(2.3d)] と完全に一致する。

## 8. \(\mathcal L_5\) と3次元物理の区別

\(y\in\mathbb R^4\)、\(r=|y|\) として \(f=f(|y|,z)\) を形式的に延長すれば

\[
\Delta_{\mathbb R^5}f
=\partial_r^2f+{3\over r}\partial_rf+\partial_z^2f
=\mathcal L_5f.
\tag{8.1}
\]

これはスカラー作用素の同型にすぎない。次は変更されない。

- 物理速度は \(\mathbb R^3\) のベクトル。
- 非圧縮条件は (3.1)。
- 物理体積要素は \(d^3x=2\pi r\,dr\,dz\)。
- 物理エネルギーは \(\frac12\int_{\mathbb R^3}|\boldsymbol u|^2dx\)。

\(\mathbb R^5\) の放射測度 \(r^3dr\,dz\) や5次元ベクトル発散を3次元診断へ代入すると、別の問題になる。

## 9. 軸 \(r=0\) の正則性と偶奇性

軸は物理境界ではなく、円柱座標の特異集合である。Cartesian で滑らかな軸対称速度は局所的に

\[
u^r=r\,a(r^2,z,t),\qquad
u^\theta=r\,b(r^2,z,t),\qquad
u^z=c(r^2,z,t)
\tag{9.1}
\]

と書ける。Liu–Wang [LW, Corollary 1 and Lemma 2] は、軸方向成分の奇数階 radial derivative と、radial/swirl 成分および角流れ関数の偶数階 radial derivative が軸で消えるという極条件を与え、これが原始変数の滑らかさに本質的であることを示している。[LW] は軸方向座標を \(x\) と書くが、本書ではそれを \(z\) と読み替える。

本規約に直すと、符号付き \(r\) への延長で

| 奇関数 | 偶関数 |
|---|---|
| \(u^r,u^\theta,\omega^\theta,\psi^\theta\) | \(u^z,p,u_1,\omega_1,\psi_1\) |

となる。したがって

\[
u^r=u^\theta=\omega^\theta=\psi^\theta=0\quad(r=0),
\tag{9.2}
\]

\[
\partial_r^{2k+1}u_1
=\partial_r^{2k+1}\omega_1
=\partial_r^{2k+1}\psi_1=0\quad(r=0)
\tag{9.3}
\]

（存在する全次数）である。とくに \(k=0\) は homogeneous Neumann 条件であるが、高階の Cartesian 正則性には高階極条件も必要である。軸対称性だけから \(z\) の偶奇は従わない。

偶関数 \(f=f_0+\frac12f_{rr}(0)r^2+O(r^4)\) なら

\[
\lim_{r\downarrow0}\left(f_{rr}+{3\over r}f_r+f_{zz}\right)
=4f_{rr}(0,z)+f_{zz}(0,z).
\tag{9.4}
\]

数値実装では \(3f_r/r\) を軸で直接除算せず、この極限または偶拡張 stencil を使う。

## 10. 元の3次元速度・渦度への復元

\(r>0\) では

\[
\boldsymbol u
=(-r\psi_{1,z})\boldsymbol e_r
+ru_1\boldsymbol e_\theta
+(2\psi_1+r\psi_{1,r})\boldsymbol e_z.
\tag{10.1}
\]

Cartesian 形式なら

\[
\boxed{
u_x=-x\psi_{1,z}-yu_1,\qquad
u_y=-y\psi_{1,z}+xu_1,\qquad
u_z=2\psi_1+r\psi_{1,r}.
}
\tag{10.2}
\]

(9.3) があれば (10.2) は軸まで滑らかに延長される。渦度の全成分も

\[
\boxed{
\omega^r=-r u_{1,z},\qquad
\omega^\theta=r\omega_1,\qquad
\omega^z=2u_1+r u_{1,r}
}
\tag{10.3}
\]

と回復できる。

### 10.1 発散の直接検算

(7.4) を3次元物理発散へ入れると

\[
\begin{aligned}
\partial_ru^r+{u^r\over r}+\partial_zu^z
&=(-\psi_{1,z}-r\psi_{1,rz})-\psi_{1,z}
 +(2\psi_{1,z}+r\psi_{1,rz})\\
&=0.
\end{aligned}
\tag{10.4}
\]

この検算は5次元作用素を一切使わない。

### 10.2 楕円関係の直接検算

(10.1) から

\[
\partial_zu^r-\partial_ru^z
=-r\left(\psi_{1,rr}+{3\over r}\psi_{1,r}+\psi_{1,zz}\right)
=r\omega_1=\omega^\theta.
\tag{10.5}
\]

よって \(-\mathcal L_5\psi_1=\omega_1\) の符号は元の curl と一致する。

### 10.3 圧力と原始変数系

十分な正則性・減衰と極条件の下で、閉じた vorticity–stream 系は原始変数系と同値である [LW]。圧力は加法定数を除き

\[
-\Delta p=\sum_{i,j=1}^{3}\partial_i\partial_j(u_i u_j)
\tag{10.6}
\]

から回復する。したがって候補証明書は、\((u_1,\omega_1,\psi_1)\) だけでなく、(10.2)、(10.6) を用いた3次元 Cartesian 残差も独立に検査できなければならない。

## 11. 次元

\([x]=L,[t]=T\) とすると

\[
[u]=LT^{-1},\quad[p]=L^2T^{-2},\quad[\nu]=L^2T^{-1},\quad
[\omega]=T^{-1},\quad[\psi^\theta]=L^2T^{-1},
\]

\[
[u_1]=T^{-1},\qquad
[\omega_1]=L^{-1}T^{-1},\qquad
[\psi_1]=LT^{-1}.
\tag{11.1}
\]

例えば (7.9) の全項は \(L^{-1}T^{-2}\)、(7.10) の両辺は \(L^{-1}T^{-1}\) で一致する。

基準長 \(L_0\)、基準速度 \(U_0\) で

\[
\tilde x={x\over L_0},\quad
\tilde t={tU_0\over L_0},\quad
\tilde u={u\over U_0},\quad
\tilde p={p\over U_0^2}
\]

とすれば、無次元粘性は \(Re^{-1}\)、\(Re=U_0L_0/\nu\) である。出力には dimensional / nondimensional の別と \(L_0,U_0\) を保存する。

## 12. Navier–Stokes スケーリング

\(\lambda>0\) に対し

\[
\boldsymbol u^{(\lambda)}(\boldsymbol x,t)
=\lambda\boldsymbol u(\lambda\boldsymbol x,\lambda^2t),\qquad
p^{(\lambda)}(\boldsymbol x,t)
=\lambda^2p(\lambda\boldsymbol x,\lambda^2t)
\tag{12.1}
\]

は同じ \(\nu\) の解である。これから

\[
u_1^{(\lambda)}(r,z,t)
=\lambda^2u_1(\lambda r,\lambda z,\lambda^2t),
\tag{12.2}
\]

\[
\omega_1^{(\lambda)}(r,z,t)
=\lambda^3\omega_1(\lambda r,\lambda z,\lambda^2t),
\tag{12.3}
\]

\[
\psi_1^{(\lambda)}(r,z,t)
=\lambda\psi_1(\lambda r,\lambda z,\lambda^2t).
\tag{12.4}
\]

[HLW, eq. (2.5)] は (12.2)–(12.3) を \(\tau=\lambda^{-2}\) の表記で確認している。循環

\[
\Gamma=ru^\theta=r^2u_1
\tag{12.5}
\]

は振幅不変に

\[
\Gamma^{(\lambda)}(r,z,t)=\Gamma(\lambda r,\lambda z,\lambda^2t)
\]

と変換される。また

\[
\partial_t\Gamma+u^r\partial_r\Gamma+u^z\partial_z\Gamma
=\nu\left(\partial_r^2-{1\over r}\partial_r+\partial_z^2\right)\Gamma.
\tag{12.6}
\]

これは独立診断に使える。

## 13. エネルギーと有限エネルギー条件

物理エネルギーは

\[
E(t)={1\over2}\int_{\mathbb R^3}|\boldsymbol u|^2\,d^3x
=\pi\int_{\mathbb R}\int_0^\infty
\left[(u^r)^2+(u^\theta)^2+(u^z)^2\right]r\,dr\,dz.
\tag{13.1}
\]

滑らかで十分減衰する解は

\[
E(t)+\nu\int_0^t\int_{\mathbb R^3}|\nabla\boldsymbol u|^2\,dx\,ds=E(0)
\tag{13.2}
\]

を満たす。\(u_1,\psi_1\) が局所的に滑らかでも、無限遠で十分減衰しなければ (13.1) は発散し得る。したがって局所プロファイルの正則性だけでは有限エネルギー初期データへの接続にならない。

## 14. 数値候補と数学的証明

数値候補が示せるのは、有限個の離散点・有限精度・有限領域で特定の近似軌道が観測されたという事実である。本当の有限時間特異点には少なくとも次が別途必要である。

1. 滑らかな有限エネルギー初期データからの厳密解の存在。
2. 離散候補近傍に厳密な軌道またはプロファイルが存在すること。
3. 空間・時間離散化、領域打ち切り、スペクトル尾部の包含評価。
4. 候補軌道へ実際に入ることと不安定方向の制御。
5. 再スケーリング時間ではなく有限の物理時刻に対応すること。
6. 座標変換の見かけではない、3次元の適切なノルムの発散。

この文書は方程式と変換を定義・監査するものであり、上記を証明していない。

## 15. 一次資料

- **[F]** Charles L. Fefferman, “Existence and Smoothness of the Navier–Stokes Equation,” Clay Mathematics Institute official problem description.
  <https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf>
- **[HL]** Thomas Y. Hou and Congming Li, “Dynamic Stability of the 3D Axi-symmetric Navier–Stokes Equations with Swirl,” *Communications on Pure and Applied Mathematics* 61 (2008), 661–697.
  DOI: <https://doi.org/10.1002/cpa.20212>
  arXiv: <https://arxiv.org/abs/math/0608295>
  著者 PDF: <https://users.cms.caltech.edu/~hou/papers/2007CPAMHouLiNSEaxis-sym.pdf>
- **[HLW]** Thomas Y. Hou, Pengfei Liu, and Fei Wang, “Global regularity for a family of 3D models of the axi-symmetric Navier–Stokes equations,” *Nonlinearity* 31 (2018), 1940–1954.
  DOI: <https://doi.org/10.1088/1361-6544/aaaa0b>
  arXiv: <https://arxiv.org/abs/1708.07536>
- **[LW]** Jian-Guo Liu and Wei-Cheng Wang, “Characterization and Regularity for Axisymmetric Solenoidal Vector Fields with Application to Navier–Stokes Equation,” *SIAM Journal on Mathematical Analysis* 41 (2009), 1825–1850.
  DOI: <https://doi.org/10.1137/080739744>
  PDF: <https://archive.ymsc.tsinghua.edu.cn/pacm_download/200/8347-Liu_Wang_SIMA_2009.pdf>
- **[LWN]** Jian-Guo Liu and Wei-Cheng Wang, “Convergence Analysis of the Energy and Helicity Preserving Scheme for Axisymmetric Flows,” *SIAM Journal on Numerical Analysis* 44 (2006), 2456–2480.
  DOI: <https://doi.org/10.1137/050639314>
