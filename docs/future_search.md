# Future search design: Type II and dynamic rescaling

## 1. 目的

この文書は、基礎実装と manufactured tests が完成した後に行う候補探索の設計であり、今回この探索を実行するものではない。目的は画像や極端な最大値ではなく、次の変換可能な成果物を得ることである。

1. 元の三次元 Navier–Stokes 解へ戻せる候補軌道。
2. ニューラルネットの重みではなく、明示的な基底係数と scale functions。
3. 各 PDE 項、軸条件、境界/tail、有限エネルギー、継続判定量を再計算できる checkpoint。
4. 区間演算による Newton–Kantorovich/radii polynomial 型検証へ渡せる有限次元方程式と tail bound。

既知の障害 [known_obstructions.md](known_obstructions.md) により、有限局所エネルギーを持つ厳密な一尺度後方自己相似解と、軸対称全空間の Type I 上界は主探索対象にしない。主系統は Type II、時間依存 profile、異方的集中、周期・準定常軌道である。

## 2. 探索を開始するための前提ゲート

次がすべて満たされるまで未知候補の最適化を開始しない。

- `docs/equation_audit.md` で使う式が「導出済み」または「一次資料で確認済み」。
- 物理三次元の発散ゼロと、形式的五次元スカラー作用素を区別したテストが通る。
- manufactured solution で空間・時間の収束次数、速度回復、楕円関係、軸 parity が確認済み。
- 非特異基準実験で空間・時間・領域収束が確認済み。
- candidate schema、config、seed、diagnostics、hash、全 run manifest が実装済み。
- [threat_model.md](threat_model.md) の TM-01 から TM-15 に対応する検査を実行可能。

## 3. 探索する力学的対象

候補を一つの固定 profile に限定せず、再スケーリングされた力学系の次の不変対象として探す。

### 3.1 定常固定点

rescaled fields と modulation rates が \(\tau\to\infty\) で定数へ収束する。これは連続自己相似または漸近自己相似を表す。ただし、元の物理解が有限局所エネルギーを持つ厳密な Leray 一尺度形なら既知非存在定理に入るため、固定点を見つけても integrability、異方性、Type II 性を直ちに監査する。

### 3.2 周期軌道

profile と modulation rates が rescaled time \(\tau\) で周期的になる。物理空間では離散自己相似または log-periodic modulation に対応する。周期軌道は backward self-similar 非存在定理を自動的に回避するわけではなく、discrete self-similar 解に対する Liouville 型制限もあるため、[Seregin2024] の仮定をフィルターにする。

### 3.3 準定常・緩慢ドリフト軌道

profile が低次元 manifold 上をゆっくり移動し、scale ratios や viscosity coefficient が \(\tau\) とともに drift する。有限の観測窓で「固定点」と誤認しやすいので、drift rate を truncation error と比較し、より長い \(\tau\) 区間・高解像度で追跡する。

### 3.4 heteroclinic / connecting orbit

滑らかな有限エネルギー初期データから、候補固定点・周期軌道・準定常 manifold の安定方向へ入る軌道を探す。profile の存在だけでは反例にならないため、最終的にはこの接続が核心になる。

## 4. 二尺度動的再スケーリング

### 4.1 監査済み物理変数

標準の軸対称変数

\[
u_1=\frac{u^\theta}{r},\qquad
\omega_1=\frac{\omega^\theta}{r},\qquad
\psi_1=\frac{\psi^\theta}{r}
\]

を想定する。以下は

\[
-\left(\partial_{rr}+\frac3r\partial_r+\partial_{zz}\right)\psi_1=\omega_1,
\quad
u^r=-r\partial_z\psi_1,\quad
u^z=2\psi_1+r\partial_r\psi_1
\tag{4.1}
\]

という符号規約に対応する。実装はこの文書から式をコピーせず、`equation_audit.md` で確定した規約から変換式を機械的・記号的に再導出する。監査の符号規約が異なる場合は、以下も同時に変更する。

### 4.2 軸中心の異方スケール

動的な半径方向スケールと動的な軸方向スケールを別々に持つ異方的な二尺度集中を、再スケーリング時間 \(\tau\) 上の力学系として扱う。

動的な半径方向スケール \(L_r(\tau)>0\)、軸方向スケール
\(L_z(\tau)>0\)、軸方向中心 \(z_c(\tau)\)、stream-function amplitude
\(C(\tau)>0\) を置き、

\[
\xi=\frac{r}{L_r(\tau)},\qquad
\eta=\frac{z-z_c(\tau)}{L_z(\tau)}
\]

\[
\psi_1=C\Psi,\qquad
u_1=\frac{C}{L_z}U,\qquad
\omega_1=\frac{C}{L_z^2}\Omega,\qquad
\frac{d\tau}{dt}=\frac{C}{L_z}.
\tag{4.2}
\]

この amplitude の組は、楕円関係、\(2u_1\partial_z\psi_1\)、\(\partial_z(u_1^2)\) の係数を正規化する。異方比と有効粘性を

\[
\delta=\frac{L_z}{L_r},\qquad
\mu(\tau)=\frac{\nu}{C L_z}
\tag{4.3}
\]

とすると、楕円関係は正確に

\[
-\left[
\delta^2\left(\partial_{\xi\xi}+\frac3\xi\partial_\xi\right)
+\partial_{\eta\eta}
\right]\Psi=\Omega.
\tag{4.4}
\]

回復される rescaled meridional velocity は

\[
V^\xi=-\xi\Psi_\eta,\qquad
V^\eta=2\Psi+\xi\Psi_\xi,
\tag{4.5}
\]

であり、物理速度は

\[
u^r=\frac{C L_r}{L_z}V^\xi,\qquad
u^z=C V^\eta,\qquad
u^\theta=\frac{C L_r}{L_z}\,\xi U.
\tag{4.6}
\]

従って、rescaled amplitude が有界でも物理速度の Type I/II 判定は
\(C\)、\(L_r/L_z\)、\(t(\tau)\) を戻して行わなければならない。

modulation rates を

\[
c_r=-\frac{d}{d\tau}\log L_r,\quad
c_z=-\frac{d}{d\tau}\log L_z,\quad
c_C=\frac{d}{d\tau}\log C,\quad
s=\frac{1}{L_z}\frac{dz_c}{d\tau}
\tag{4.7}
\]

と定義すると、監査済みの標準 \(u_1,\omega_1,\psi_1\) 系から次の rescaled system を得る。

\[
\begin{aligned}
U_\tau
&+(c_r\xi+V^\xi)U_\xi
+(c_z\eta-s+V^\eta)U_\eta\\
&=-(c_C+c_z)U+2U\Psi_\eta
+\mu\left[
\delta^2\left(U_{\xi\xi}+\frac3\xi U_\xi\right)+U_{\eta\eta}
\right],
\end{aligned}
\tag{4.8}
\]

\[
\begin{aligned}
\Omega_\tau
&+(c_r\xi+V^\xi)\Omega_\xi
+(c_z\eta-s+V^\eta)\Omega_\eta\\
&=-(c_C+2c_z)\Omega+\partial_\eta(U^2)
+\mu\left[
\delta^2\left(\Omega_{\xi\xi}+\frac3\xi\Omega_\xi\right)+\Omega_{\eta\eta}
\right],
\end{aligned}
\tag{4.9}
\]

と (4.4)。scale の復元は

\[
\frac{dL_r}{d\tau}=-c_rL_r,\quad
\frac{dL_z}{d\tau}=-c_zL_z,\quad
\frac{dC}{d\tau}=c_CC,\quad
\frac{dz_c}{d\tau}=sL_z,\quad
\frac{dt}{d\tau}=\frac{L_z}{C}.
\tag{4.10}
\]

式 (4.8)–(4.10) は将来実装前に symbolic differentiation と元の物理 residual への往復で再監査する。小さい rescaled residual のみを採用条件にせず、復元した物理場で residual を独立計算する。

### 4.3 Type I と Type II の識別

等方一尺度 \(L_r\asymp L_z=L\)、\(C\asymp L^{-1}\) なら
\(\mu=\nu/(CL)\) は \(O(1)\)、物理速度は \(O(L^{-1})\) である。
\(L\asymp\sqrt{T-t}\) は Type I の自然レートに対応する。

候補は、復元した物理場について

\[
\sqrt{T-t(\tau)}\,
\|u(t(\tau))\|_{L^\infty(\mathbb R^3)}
\tag{4.11}
\]

が有界という [CSTY2009] の排除領域に入らないことを示す必要がある。Type II を「\(c_r\) が大きい」「\(\mu\) が小さい」など一つの rescaled 指標だけで定義しない。

また、物理時刻が有限であるための必要条件は

\[
T-t(\tau_0)
=\int_{\tau_0}^{\infty}\frac{L_z(\sigma)}{C(\sigma)}\,d\sigma<\infty.
\tag{4.12}
\]

この積分の数値的有限性は証明ではなく、最終的には scale rate の厳密な上下界が必要である。

### 4.4 軸へ移動する ring と radial two-scale

集中核が \(r=0\) に中心を持たず、ring 半径 \(R(\tau)>0\) と厚さ
\(\ell_r(\tau)\) が別レートで縮む可能性を別 branch とする。

\[
\xi=\frac{r-R(\tau)}{\ell_r(\tau)},\qquad
\eta=\frac{z-z_c(\tau)}{\ell_z(\tau)},\qquad
\rho(\tau)=\frac{R(\tau)}{\ell_r(\tau)}.
\tag{4.13}
\]

このとき軸は moving boundary \(\xi=-\rho\) であり、radial scalar operator は

\[
\partial_{rr}+\frac3r\partial_r+\partial_{zz}
=\frac1{\ell_z^2}
\left[
\delta^2\left(\partial_{\xi\xi}
+\frac{3}{\rho+\xi}\partial_\xi\right)
+\partial_{\eta\eta}
\right],
\qquad \delta=\frac{\ell_z}{\ell_r}.
\tag{4.14}
\]

\(R'\) による drift を含む変換を元の方程式から再導出する。二つの重要な regime は次である。

- \(\rho=O(1)\): ring の軸までの距離と厚さが同程度。
- \(\rho\to\infty\) かつ \(R\to0\): ring は軸へ近づくが、厚さはさらに速く縮む真の二尺度集中。

軸対称特異点は軸上に限られるため、\(R(\tau)\not\to0\) の branch は有限時刻特異点候補として棄却する。moving core 座標で軸が遠方に見えても、物理軸 parity と有限エネルギーを捨ててはならない。

## 5. modulation gauge

scale を場から一意に決めないと、profile drift と座標 drift を区別できない。点wise maximum は格子ノイズに敏感なので、主 gauge は滑らかな重み付き moment とし、maximum gauge は補助にする。

例として非負 monitor

\[
Q=U^2+\kappa_\Omega\Omega^2
\]

を compact weight \(w(\xi,\eta)\) とともに用い、

\[
\int \eta\,wQ\,\xi\,d\xi d\eta=0,\qquad
\frac{\int \eta^2wQ\,\xi\,d\xi d\eta}{\int wQ\,\xi\,d\xi d\eta}=1,
\]

\[
\frac{\int \xi^2wQ\,\xi\,d\xi d\eta}{\int wQ\,\xi\,d\xi d\eta}=1,\qquad
\int wU^2\,\xi\,d\xi d\eta=1
\tag{5.1}
\]

のような条件で \(z_c,L_z,L_r,C\) を決める。実際の gauge は Jacobian が非特異かを数値・解析の両方で確認する。

ring branch では radial centroid を \(R\)、radial variance を \(\ell_r\) とする。gauge を変えても復元した物理場が一致することを必須テストにする。

## 6. 明示的な候補表現

### 6.1 空間基底

軸中心 branch では、軸で滑らかな \(u_1,\omega_1,\psi_1\) は \(r\) の偶関数として拡張できる。従って次のいずれかを使う。

- \(r^2\) の mapped Chebyshev/Jacobi 基底。
- 偶次数 Chebyshev 基底。
- 軸正則性を組み込んだ compactly supported B-spline/FEM 基底。

\(z\in\mathbb R\) には rational Chebyshev、複数 mapped Chebyshev patch、または「有限核 + 解析 tail」の分割を使う。周期境界を全空間問題の代用にしない。追加の \(z\)-反射対称性は初期値が本当に持つ場合だけ課す。

候補ファイルは最低限次を持つ。

- \(U,\Omega,\Psi\) の基底名、map、次数、係数。
- \(L_r,L_z,C,z_c\) または \(R,\ell_r,\ell_z,C,z_c\) の明示係数。
- fixed/periodic/quasi-stationary の種別と \(\tau\) 区間・周期。
- viscosity、単位、符号規約、軸 parity、tail model。
- residual norm、oversampling grid、truncation、hash、生成 commit。

### 6.2 rescaled time 基底

- fixed point: 時間係数なし。
- periodic orbit: Fourier series in \(\tau\) と位相条件。
- finite connecting segment: piecewise Chebyshev in \(\tau\) と continuity constraints。
- quasi-stationary orbit: slow parameter を追加し、二変数係数または短い Chebyshev segment の列として保存。

scale functions を unrestricted spline のままにせず、最終候補では有理係数多項式、Chebyshev 係数、または区間演算可能な ODE と初期区間へ変換する。

### 6.3 tail

有限領域のゼロ境界だけでは証明書にならない。各変数について、

- algebraic/exponential decay exponent の候補、
- weighted \(\ell^1\) または Sobolev norm の係数 tail、
- 外部領域の楕円 Green contribution、
- 物理三次元エネルギーの tail

を明示する。未知 tail をゼロと置いた profile は探索 seed には使えても証明書化候補にはしない。

## 7. AI / neural network の役割

AI は探索器または低次元 manifold の提案器としてのみ使い、最終表現にはしない。

### 7.1 許される用途

- 多数の初期値・scale gauge から有望 branch を順位付けする。
- rescaled profile の初期 guess を生成する。
- slow manifold の座標、周期候補、unstable direction の近似を学習する。
- continuation の predictor を作る。

### 7.2 hard constraints

可能な限り parameterization 自体に次を組み込む。

- 軸 parity。
- stream function からの速度回復による三次元発散ゼロ。
- \(\Gamma=ru^\theta=O(r^2)\)。
- gauge 条件。
- 初期・遠方条件。

PDE、楕円関係、局所エネルギーを loss の小ささだけで満たしたことにしない。

### 7.3 明示基底への蒸留

1. 学習出力を network collocation 点と独立な oversampled 点で評価する。
2. parity-compatible basis へ weighted projection する。
3. 係数 tail を見て次数を増やす。
4. ネットワークを捨て、係数だけを Newton–Krylov/trust-region で離散方程式へ補正する。
5. 自動微分を使わない独立 residual で検証する。
6. 保存成果物を coefficients + maps + scales + metadata とし、重みは任意の探索履歴に格下げする。

蒸留後に候補が消える場合は失敗 run として残す。[Rahaman2019] のスペクトルバイアスと [Krishnapriyan2021] の optimizer failure を [threat_model.md](threat_model.md) の受入れ試験へ反映する。

## 8. 数値探索パイプライン

### Stage A: 物理座標での seed 探索

- 滑らか・有限エネルギー・軸適合の明示初期データ族だけを使う。
- 低～中解像度 parameter sweep では blow-up claim を行わず、集中、伸長/散逸比、scale separation を順位付けする。
- 全 run を manifest に残す。

### Stage B: scale extraction

- 集中幅、center、amplitude を複数 gauge で抽出する。
- (4.11) の Type II 指標、(4.12) の物理時間、ring 比 \(\rho\)、有効粘性 \(\mu\) を追跡する。
- gauge 変更後も物理場が一致する branch だけを残す。

### Stage C: rescaled evolution

- 物理 solver と rescaled solver を重複区間で走らせ、往復誤差を測る。
- fixed/periodic/quasi-stationary のいずれに近いかを、profile difference、autocorrelation、Floquet/linearized spectrum で分類する。
- 小 residual だけで分類しない。

### Stage D: coefficient solve

- 抽出 profile を明示基底へ射影する。
- fixed point なら nonlinear algebraic system、periodic orbit なら space-time system、connecting orbit なら multi-segment boundary-value systemを解く。
- gauge/phase conditions を含め、Jacobian の null direction を除く。

### Stage E: convergence and obstruction screening

- 空間次数、時間次数、領域、tail、precision を系統的に増やす。
- [known_obstructions.md](known_obstructions.md) の O-SS～O-V を全て評価する。
- [threat_model.md](threat_model.md) の停止規則に一つでも触れれば branch を昇格しない。

### Stage F: certificate preparation

- approximate inverse、linearized operator bound、nonlinear Lipschitz bound、tail bound を作る。
- interval arithmetic で residual と inverse defect を外向き丸め評価する。
- fixed/periodic profile の局所存在だけでなく、物理初期データからの connecting orbit と不安定方向を別証明義務として残す。

## 9. 不安定方向と滑らかな初期データからの接続

### 9.1 線形化

fixed point \(\bar X=(\bar U,\bar\Omega,\bar\Psi,\bar c)\) では、gauge を含む rescaled operator の線形化

\[
\partial_\tau h=\mathcal L h
\]

を作る。周期軌道では一周期 monodromy operator の Floquet multiplier を求める。

必須の区別:

- gauge/対称性による中立方向。
- 真の不安定方向。
- truncation で生じた spurious eigenvalue。
- essential/tail spectrum。

固有値は次数・領域を変え、adjoint mode と residual を確認する。unstable dimension は候補へ入る初期値集合の codimension を決める。

### 9.2 edge tracking / multiple shooting

1. 明示初期データ族 \(u_0(a_1,\ldots,a_k)\) を選ぶ。
2. 候補 manifold の unstable coordinate の符号を event function にし、regularizing side と concentrating side の境界を bracket する。
3. bisection/continuation で edge trajectory を追う。
4. 物理座標と rescaled 座標の両方で接続を確認する。
5. 最終的には multiple shooting と interval enclosure で、有限時間区間の orbit segment を囲う。

初期データをネットワークの潜在変数だけで表さず、compactly supported/rational spectral coefficients として保存する。

### 9.3 nonlinear stability

候補が一つの精密調整された軌道でも反例にはなりうるが、接続の証明には stable/center-stable manifold が必要になる可能性が高い。少なくとも次を分離する。

- profile の局所存在。
- rescaled flow での局所安定性。
- 物理有限エネルギー空間での摂動制御。
- tail と外部領域からの feedback。
- 初期値族が stable manifold と交わること。

## 10. 主要診断

各候補 branch で次を同じ物理時刻に揃えて保存する。

- \(L_r,L_z,C,z_c\)、ring branch では \(R,\ell_r,\ell_z,\rho\)。
- \(c_r,c_z,c_C,s,\delta,\mu\)。
- 物理 \(\|u\|_\infty\)、\(\|\omega\|_\infty\)、エネルギー、散逸。
- \(\sqrt{T-t}\|u\|_\infty\) と BKM/Serrin/渦度積分量。
- \(t(\tau)\) の tail integral estimate。
- rescaled profile drift \(\|X(\tau+\Delta)-X(\tau)\|\)。
- fixed point residual、periodic return residual、phase/gauge residual。
- physical PDE residual と rescaled PDE residual。
- 係数 spectrum、tail bound、境界 flux、軸 parity。
- linearized eigenvalues/Floquet multipliers と truncation sensitivity。

## 11. 候補の段階ラベル

| レベル | 意味 |
|---|---|
| S0 | 単一 run で集中が見えた。研究上の候補ではない |
| S1 | 空間・時間・領域系列が収束し、脅威モデルの基本試験を通った |
| S2 | 動的再スケールで fixed/periodic/quasi-stationary 構造が再現した |
| S3 | 明示基底係数へ蒸留し、独立実装で residual と診断を再現した |
| S4 | profile/orbit の局所存在を区間演算で検証した |
| S5 | nonlinear stability、初期値からの接続、有限物理時刻、物理ノルム発散を検証した |
| S6 | 元の三次元問題への完全なコンピューター支援証明と独立監査が完成した |

S0–S3 は数値結果であり特異点の証明ではない。S4 も profile の局所存在だけなら反例の証明ではない。必要な依存関係は [proof_obligations.md](proof_obligations.md) に従う。

## 12. 先行研究との位置づけ

[Hou2023] は、滑らかな有限エネルギー初期値から軸対称 Navier–Stokes の potentially singular behavior を高解像度で調べた重要な数値先行例である。しかし同論文も数値的証拠として述べており、証明ではない。本設計では、同様の大振幅観測を出発点にしても、最終成果を明示係数、independent residual、tail bounds、interval proof obligations へ変換する点を必須にする。

[Seregin2024] と [Seregin2026] は、Type II blow-up から Euler scaling の ancient limit を取り出す際の必要条件・Liouville 障害を与える。いずれも探索空間を狭める有用なフィルターだが、プレプリントであること、仮定が局所混合ノルムや追加可積分性を含むことを明示し、要旨だけで候補を排除しない。

## 13. 最初に実装すべき最小の一手

未知候補探索の最初の追加実装は、長時間・高価な最適化ではない。次の小さな機能を順に作る。

1. 既存の非特異基準 run から \(L_r,L_z,C,z_c\) を moment gauge で抽出する pure diagnostic。
2. (4.2) で profile を作り、(4.10) で物理場へ完全に戻す round-trip test。
3. 復元場で divergence、楕円関係、エネルギー、physical PDE residual が元 checkpoint と一致するテスト。
4. \(t(\tau)\)、Type I 指標、profile drift を JSON/CSV に出す。

この round-trip が manufactured field と非特異基準流で収束するまで、(4.8)–(4.9) の rescaled time stepper、AI探索、periodic orbit solve へ進まない。
