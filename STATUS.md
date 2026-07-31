# Project status

最終更新: 2026-07-31 第 10 便(branch `fable5-mainline`)
状態: **第 10 便で Track P の単発スラブを連結し、条件付き certified existence
interval を実際に延長した。** 各スラブは厳密有理・厳密発散ゼロの再中心化点から
開始(Taylor 終端包絡 → dyadic 丸め → 厳密 Leray 射影; 捨てた幅は
δ_{n+1} = R_n(t_{n+1}) + transfer としてスカラー半径に課金)するため、区間 box が
スラブ境界を越えて伝播せず **wrapping は構造的に不発生**(前登録の Lohner/QR
導入条件は実測で不発火 — transfer は丸め床 ~1e-8 に張り付き、tube 半径 1e-3
に対して無視可能)。スイープ 12 連結(P1/P2/P3 × ν∈{1/4,1/10,1/40,1/100}、h₀=1/2048)は スラブ数 12〜15、certified horizon 2.01e-03〜4.88e-03(単発スラブ長の 4.1〜10.0 倍)。長尺(P1、ν=1/10、h=1/8192、倍増なし)は 41 スラブ連結、T=4.90e-03(1/2048 基準で 10.0 倍)。停止分類は ['control_linear_coefficient']、checker 全数合格=True。
停止は前登録分類法で必ず分類され、**一貫して「control ODE の粗い線形係数
9(K₁+K₂) が縛り」= 解の性質ではなく方法の限界**(Riccati 天井
T* ≈ (1/a)log(a²/(bε)) ≈ 5.7e-3 @P1)。**証明区間の終了は特異点の主張ではない**
(checker が文言を強制)。Lean は連結の有限不等式骨格 9 定理
(`TrackPChain.lean`: 2 スラブ合成・n スラブ帰納法・転送三角不等式・離散
Grönwall・Lagrange 終端剰余)+ 二次 ODE の Picard–Lindelöf 5 定理
(`GalerkinPicard.lean`: 明示 Lipschitz 定数・存在区間半幅 ε=1/(L+1)・一意性 —
EXT-P1 の Galerkin 半分の有限次元核)を追加、`lake build` 8670 jobs、監査
110/110 が古典 3 公理のみ。EXT-P1 は**完全な紙上証明(未監査)**と 14 行依存表を
`ext_dependencies.md` に整備(一様 H⁴ 界は本 repo の可換子代数の再利用:
d/dt Y ≤ −2νY + 270Ȧ Y^{3/2}、T* = 1/(270Ȧ‖u₀‖_Ḣ⁴))— **payload の
proved:false は監査完了まで不変**。EXT-P2 の各時刻 Dini 節は未達で外部のまま
(正直に記録)。R³ 側 J>0 は勾配形式 −3∫|u|u·∇p_h で再構成:
**圧力×流束の非相関化は成功**(P 下界が各解像度で 3〜4 倍改善、依存過大評価は
負 = 搾れる相関はもう無い)が、**閉じない** — 拘束幅は ν=1e-3 では粘性上界の
過大評価(16〜90 倍)、ν→0 極限では速度因子の cell 包絡。しかも float 再計算で
J>0 となるのは ν < ~1e-4 のみ(第 9 便の参照値 −2.5e-3/+1.1e-5 はコードと
不整合と判明 — 本便の測定が上書き)。最良 P 下界 −5.13e-6 vs 真値 +3.25e-7、
必要精細化 ~40 倍(~10⁴ 倍コスト)。**候補昇格ゼロ、δ_J なし、短時間継続
ドラフトは発火条件未達のため書かず。** 新離散圧力仮定 P1G は payload に
proved:false で記録され checker が強制。
**特異点証明ではなく、Clay 問題も解決していない。**

前便まで: 第 9 便で Track P 新設(12/12 スラブ)+ Gaussian–Hermite 基底、
第 8 便で L3 生成恒等式・純粋旋回 no-go、第 7 便で Picard 領域離脱、
第 3〜6 便で有限帯域 no-go・Gate 群・区間証明書層。

## 2026-07-31 第 10 便の結果(fable5-mainline)

### 無条件に証明したこと(Python 厳密有理 + Lean)

- **スラブ連結機構**(`torus_chain.py`): Galerkin 軌道の Taylor 終端包絡
  (Lagrange 剰余は Picard box 上の係数漸化式評価)、dyadic 丸め + 厳密 Leray
  射影による再中心化(射影はモード毎直交なので Ḣⁿ を縮小 — 誤差を増やさない)、
  δ 漸化式、全リンク再計算 checker(datum 完全一致・係数組立・δ 漸化式・文言
  契約・EXT 記録を全数再検証)、前登録適応則(失敗時半減 → Taylor 次数昇格 →
  相対半径 cap)と停止分類法。
- **Lean 14 定理追加(計 110、全て古典 3 公理のみ)**: 上記の有限不等式骨格
  (`TrackPChain.lean` 9)+ 二次 ODE Picard–Lindelöf(`GalerkinPicard.lean` 5)。
- **J>0 勾配形式の離散包含**(`gaussian_gradient_certificate.py`): ∇p_h を
  節点値の厳密有理線形結合として扱い、p と ∇p_h に独立 hull を与えない。
  離散 box 積分の包含は無条件; 真の圧力への接続は新記録仮定 P1G
  (proved:false、checker 強制)の条件付き。

### EXT 条件付きで証明したこと

- スイープ 12 連結(P1/P2/P3 × ν∈{1/4,1/10,1/40,1/100}、h₀=1/2048)は スラブ数 12〜15、certified horizon 2.01e-03〜4.88e-03(単発スラブ長の 4.1〜10.0 倍)。長尺(P1、ν=1/10、h=1/8192、倍増なし)は 41 スラブ連結、T=4.90e-03(1/2048 基準で 10.0 倍)。停止分類は ['control_linear_coefficient']、checker 全数合格=True。各連結は (i) Picard box(無条件)、(ii) スラブ上有効な全定数
  (無条件)、(iii) control ODE 管 δ_n スタート(無条件)、(iv) EXT-P1/P2/P3
  (忠実記録・Lean 公理化なし・スラブ境界の解の同一性は EXT-P1 の一意性節)の
  条件付きで「真の周期強解が [0,T] 全体に存在し ‖u−u_a‖_Ḣ⁴ ≤ R(t)」。
- 物理量追跡: 各再中心化点で energy/enstrophy/Ḣ³/Ḣ⁴/shell energy/**厳密 shell
  flux** を厳密有理で、L³/渦度 sup を Lipschitz 転送半径付き
  (|Q(u)−Q(u_a)| ≤ L_Q R、L=1/Ȧ/(‖u_a‖₀+R/2) 等の表を証明書に同梱)で保存。
  追跡量はすべて減衰 — 「certified growth」と呼べるものは無い。

### 棄却・不成立(誠実記録)

- J>0 区間証明は勾配形式でも不成立(上記)。候補昇格ゼロ。
- Lohner/QR frame は**導入しない**: 前登録 3 条件が全て不発火(設計上 box が
  境界を越えないため wrapping 増幅が存在しない)。導入すれば削除対象だった。
- 第 9 便の J 参照値(margin −2.5e-3、float P +1.1e-5)はコードと不整合と判明。
  本便の測定(P ≈ +3.3e-7、float ν_crit ≈ 1.02e-4)で上書きし、原因調査は
  次便の宿題として記録。

### 未証明のまま(従来どおり + 新規)

- EXT-P1/P2/P3(EXT-P1 は紙上証明あり・未監査)、P1G(新規)、HS-5 全空間版、
  NT-N1、P1、H3。n=3 閉鎖・鋭い Kato 定数(= 地平 T* を延ばす唯一のレバー)は
  未着手。


## 2026-07-30 第 9 便の結果(fable5-mainline)

### 数学的に証明したこと

- **H⁴ control 不等式の自前導出**(`docs/research_notes/track_p_periodic.md`)。
  平均ゼロ場のスペクトルギャップ `|k| ≥ 1` による減衰項 `−νR` の回復、
  多項係数 81(√81 = 9)、multi-Vandermonde による二項係数、
  埋め込み `‖f‖_∞ ≤ Ȧ‖f‖_Ḣ²` の格子定数
  `Ȧ² ≤ Σ_{|k|_∞≤20}|k|⁻⁴ + 26/20 ≤ 17.33`(尾部は殻計数 `24m²+2 ≤ 26m²`)。
  n=4 で閉じる理由(|β|=1,2 は自因子に、|β|=3,4 は相手因子に sup を置く)と
  n=3 では閉じない理由(L⁶×L³ の Gagliardo–Nirenberg 定数 = MP の G₃ = 0.438
  が必要)を明記。
- **残差の厳密性**: Galerkin 軌道の連続 PDE 残差はちょうど Galerkin tail であり
  有限三角多項式(HS-5 の周期版の構成的閉鎖)。dealiasing 誤差は存在しない
  (convolution が厳密)— その役割は厳密計算される tail が担う。
- **固定帯域軌道 vs 有限帯域初期値の区別**と、後者が前者を含意しない反例。
- **H³ 変換**: `‖w‖_H³ ≤ √8 ‖w‖_Ḣ⁴`((1+|k|²)³ ≤ 8|k|⁸)。

### Lean で証明したこと(21 定理追加、全 96 定理が古典 3 公理のみ)

- `TrackPFourier.lean`(14): Leray 乗数の有限代数 5 定理(k=0 でも真)、
  単一モード slot-divergence 定理(cos/sin 両方)、有限三角多項式の
  `ContDiff ℝ ⊤`、`FixedBandTrajectory → FiniteBandDatum` と**逆の反例**
  `exists_finiteBandDatum_not_fixedBandTrajectory`(u t = (1,t))、既存 no-go
  との合成 `FixedBandwidthCandidate.fixedBand_scope`、Ḣ 梯子単調性 2 定理、
  control ODE 層との合成 `trackP_slab_error_le`。
- `GaussianTransfer.lean`(7): 多項式×Gaussian の微分閉包(witness
  `p′−2αXp` 明示)、J 連続性の有限不等式 3 定理(spline→smooth transfer の
  部品)、`torus_control_bound`。
- `lake build` 8668 jobs 成功。`sorry`/`admit`/project 公理ゼロ。
  EXT-P1/2/3 は Lean に**入れていない**(全 Lean 定理は無条件に真)。

### 区間証明書で証明したこと

- **Track P スラブ証明書 12/12**(下記)。independent checker は係数組立
  (K₁,K₂→linear、135Ȧ→quadratic)・√8 変換・埋め込まれた control 証明書を
  再検証し、EXT の「proved/公理化」偽装・係数改竄・免責削除を拒否する。
- **Gaussian–Hermite 区間包絡**: divergence 包絡は 1/8 セル・bits=48 で
  ±0.38(勾配スケール ~1)— flat bump の ±145(スケール 0.05)から構造的に
  解消。粘性積分の厳密下界は真値の 7 倍以内(flat bump は 47 桁外)。
- **`J>0` は依然として区間で閉じない**: margin = −2.5e-3(真の P = +1.1e-5)。
  残る障害は離散圧力 × 流束の積包絡の幅で、次の一手として `−3∫|u|u·∇p_h`
  形(∇p_h は節点値の厳密線形結合)への書換えを記録。**規則どおり候補は
  昇格させていない。**

### 周期軌道について厳密に保証したこと

`outputs/track_p_slab_v1/`: 3 族 × スラブ長 h ∈ {1/2048, 1/1024, 1/512,
1/256} の 12 証明書すべてで、(i) Galerkin 軌道の Picard 包含(厳密)、
(ii) スラブ上有効な全定数(厳密)、(iii) control ODE 管(厳密)、
(iv) EXT-P1/2/3 条件付きで「真の周期強解が存在し ‖u−u_a‖_Ḣ⁴ ≤ R(t)」。
相対距離: P1 6.9e-4〜7.8e-4、P2 4.4e-4〜5.1e-4、P3 7.7e-4〜9.0e-4。
**これは軌道近傍の正則性の証明であり、特異点証明の反対物である。**

### R³ 候補について保証したこと

- Gaussian–Hermite 族は Cartesian C^∞(r² の関数)・厳密発散ゼロ(恒等式)・
  急減衰・有限 energy/L³・全微分閉形式(Lean 定理)・解析的 tail 上界
  (Mills 比)・有理パラメータ化を満たす(`gaussian_hermite.py`、テスト済み)。
- J 最適化(同一目的 Re_crit、同一グリッド規約): 316 (41×81, 4start×12iter) → 検証 331 (97×193) / 333 (129×257); flat-bump 最良 1410 の 4.2 倍改善だが kill 条件 1e2 は未達。
  flat-bump 最適値(1.41e3 @129×257)に対し **4.2 倍の改善**(333 @129×257、
  解像度安定)。ただし flat-bump 実験で登録した kill 条件(Re_crit < 1e2)には
  依然届かず、**バーは動かさない**: 形状因子レーンの棄却は Gaussian 基底でも
  維持する(バーの再登録は、するなら結果を見る前に行うべきだった)。

### 未証明の関数解析定数と残る義務

- EXT-P1/P2/P3(古典; 忠実記録のみ)。HS-5 の全空間版、NT-N1、P1(圧力離散
  化)、H3(半離散→連続)は従来どおり未解決。スラブ連結(Lohner frame)は
  未着手。n=3 での control 閉鎖は MP 型 Kato 定数の自前証明が必要。


## 2026-07-29 第 8 便の結果(fable5-mainline)

### 数学的に証明したこと

- **L3 生成恒等式**(`docs/research_notes/l3_generation_rate.md`)。仮説は
  `|u| + |∇u| ≤ C_J <x>^{-4}` を時間について局所一様に、`p = R_i R_j(u_i u_j)`。
  **Schwartz 類は誤り**: Navier-Stokes は Schwartz 性を保存しない
  (Brandolese の局在化障害。純粋旋回自身が `M = diag(a,a,0)` で反例)。
  正則化 `s_ε`、DCT による項別極限、`ε → 0` の一様収束による `F ∈ C^1`、
  境界項 `O(ε R^{-5})` まで明示。**この恒等式は新しくない**: ESS の背後にある
  古典的 `L^p` エネルギー恒等式の `p=3` の場合である。
- **移流項の厳密消滅**: `|u| u·(u·∇)u = (1/3) u·∇(|u|^3)` は純粋な発散。
  全ての発散ゼロ場、全ての `ε > 0` で厳密にゼロ。極限も小ささも不要。
- **純粋旋回 no-go**(軸条件 `u^θ = r g(r^2,z)` 込み)。等号条件は連結性の議論
  (`∫|u||∇u|^2 = 0 ⇒ ∇u = 0` という一行は non-sequitur)。
- **軸パリティ選択則**とその証明。
- **スケーリング**: `P ~ A^4 L^2`、`|V| ~ ν A^3 L`、比は `AL/ν`。臨界スケーリング
  `A=λ, L=1/λ` で比は不変 ⇒ 臨界族の内部で再スケールしても `J>0` は作れない。

### Lean で証明したこと

- `formal/NSSingularity/L3Generation.lean`(18 定理): `d/dt||f||^3` の連鎖律
  (零点でも有効)、移流相殺 `||f|| <f,f'> = (1/3) deriv(||f||^3)`、正則化速度の
  `C^∞` 性と挟み込み、Kato 分割不等式(区間証明書が `|u|` で割らずに済む根拠)、
  `0 < J ⇒ 0 < P`、`Re_crit` の同値、**純粋旋回の Cartesian 発散ゼロ**、
  **純粋旋回が任意の軸対称スカラーの勾配と各点で直交**、等号条件の位相的部分。
- `formal/NSSingularity/ControlODE.lean`: Chaplygin-Dini 比較、**HS-6(変係数
  Grönwall)を積分形まで**、完全 Riccati 比較、二次爆発時刻。
- `lake build` 成功(8666 jobs)、`#print axioms` は**全 75 定理**が
  `[propext, Classical.choice, Quot.sound]` のみ。`sorry`/`admit`/project 固有
  axiom はゼロ。

### 区間証明書で証明したこと

- **スラブ証明書の H1 を定理化**: 節点値から導関数を推測する代わりに、補間子を
  **定義**し、その bicubic Hermite 係数を Bernstein 基底へ変換する。テンソル
  Bernstein 基底は非負かつ単位分割なので、補間子はセル上で係数の凸結合であり、
  16 係数の凸包が厳密な包絡になる。**膨張因子ゼロ、仮説ゼロ**。
  divergence 上界は 2.13e-4 → 4.30e-5 と、むしろ狭くなった。
- **H2 は撤回**。`|y - H| ≤ Δ^4 M_4/384` は `H` が**一本の**軌道を両端で補間する
  ときのみ有効だが、リポジトリの終端状態は RK4 の出力であって始端を通る軌道の
  `t_{n+1}` での値ではない。`M_4` を証明しても救えない。枠組み自体が不健全だった。
- 残る仮説は **H3(半離散 → 連続)** のみで `proved: false`。checker は定理を
  格下げする payload も仮説を格上げする payload も拒否する。
- **`J > 0` の厳密証明書は失敗**。構造と checker は正しいが、包絡が 60 桁緩い
  (粘性積分の上界 1.44e53 対 真値 4.4e-2)。原因は台の縁での dependency problem:
  `χ''` が `(1-σ)^{-4}` を持ち、区間評価は `χ` の超多項式減衰との相殺を見られない。
  厳密にゼロの divergence が `±1.45e2` に囲まれることが診断指標。
  **規則どおり候補は昇格させていない。**

### 数値的に観測したこと

- 混合族の離散 divergence は `3e-16` 相対(解析微分による構成)。
- `P` は `u -> -u` で符号反転、`V` は全桁一致で不変(表示桁すべて)。
- 最適化: 1.27e4(手設計)→ 1.13e3(粗格子)→ 検証 2.04e3 → 1.30e3(中格子)
  → 検証 1.37e3、1.41e3。**粗格子の最適解は格子に過適合する**。
- 領域独立性: 4 種の box で `Re_crit` が 2093/2042/2059/2045。
- 移流残差は `P` の 0.2-0.9%(恒等式の検証として十分小さい)。

### 正の L3 生成候補と棄却

- **候補**: 最適形状は `P > 0` かつ `Re > 約 1.4e3` で `J > 0`。`Re_crit` の半分で
  負、2 倍で正であることを直接確認。
- **棄却**: 事前登録キル条件 `Re_crit < 1e2` は満たされない(1.41e3)。形状最適化は
  手設計比 9 倍しか稼げず、汎用場に対してはほぼ何も稼げなかった。
  **形状因子レーンは事前登録の規則により棄却。**
- **厳密性による非昇格**: 区間で `J > 0` を示せないので昇格させない。

### HS-5 と a posteriori 框架

- `docs/research_notes/hs5_function_space.md`: 主空間を**物理 3 次元の整数階
  `H^n(R^3)`, n = 3 or 4** に固定。5D lift の `H^s` 代数性から 3D PDE 安定性を
  主張してはならない理由を列挙(5D lift は `div_5 = 2u^r/r ≠ 0` で圧縮性、
  5D 圧力は存在しない、`u_1` は階数 n だが `ψ_1` は n+1、`r` 倍は `s ≥ 2` で
  `H^s(R^5)` 乗数でなく障害は軸であって無限遠ではない)。
  **norm transfer は厳密等距(定数 1)** であり、不等式ではない。
- `docs/research_notes/a_posteriori_frameworks.md`: **Morosi-Pizzocchero は
  移植しない**。control ODE の `-ν R_n` 項はトーラスの `|k| ≥ 1` にのみ由来し、
  `R^3` では `-Δ` のスペクトルが `[0,∞)` なので消える。これが最も重い障害で、
  smallness からの大域存在も減衰系も全て失われる。再利用可能なのは比較補題、
  項別導出、BKM 基準、`n > d/2+1` の閾値のみ。


## 2026-07-29 第 7 便の結果(fable5-mainline)

### 主結果

- **無次元化**(`src/ns_certificate_lab/nondimensional.py`): 監査済み系への
  `r=Lρ, z=Lζ, τ=At, u₁=AU, ω₁=(A/L)W, ψ₁=ALΨ` 代入を項別に実行し、残る
  パラメータが `Re = AL²/ν` のみであることを確定。振幅・長さ・粘性を独立に
  探索することは以後 `deduplicate_settings` が禁止する。
- **32 点スイープの再分類**: 重複ペアは 0(無駄な重複はなかった)だが、
  異なる無次元形状は **6 種のみ**、到達 `τ` は **0.00163 – 0.0233**。
  目標 `τ = 1` に対して **43 倍不足**していた。
- **Picard 梯子**(`picard_continuation.py`): level 0(純拡散)/1/2 と完全解を
  同一 RK4・同一 accepted step で同時積分。軌道間補間なし、snapshot 差分なし。
  accepted step ごとに solver 内部の正確な RHS を保存。
- **τ 継続**: 前登録チェックポイント τ = 0.025/0.05/0.1/0.2/0.4/0.7/1.0 に
  adaptive CFL + diffusion stability で到達。18 run すべて τ=1 完走。
- **Re 継続**: Re = 10/25/50/100/200/400 × 族 S/A/H。乖離は Re とともに単調増加
  (S: 0.0094 → 0.2947)。
- **前線解像度の修正**: 初回は前線点数 2 項目のみ不合格(6.5/7.0 対 10)。
  **閾値は変えず**、145×289 格子を追加して同じ 10 点バーを満たした。修正理由と
  「閾値不変」は config の `amendments` に記録し `summary.json` へ複製。

### 結果は否定的である(隠していない)

18 run すべてで:

- 臨界 `L³` ノルムは**単調減少**(比 0.037 – 0.919)。Re→∞ で 1 に近づくが
  増分は縮小しており(S: +0.093, +0.054)、1 未満で飽和する方向。
- dyadic shell 数は 16 run で**増加**(集中ではなく拡散)。減少は族 A の
  Re≥200 のみ(−0.03, −0.06)。
- 幅比はほぼ 1。最小は A_Re10 の 0.82 で、これは最粘性設定の拡散由来。

**昇格候補ゼロ。** 最良の A_Re400 でも前登録 6 条件中 1 条件のみ充足
(Δshell = −0.0612 ≤ −0.05)、`L³` 成長 0.8953(要 ≥1.05)と幅比 0.9925
(要 ≤0.90)で不合格。棄却理由は全 run 分が `summary.json` に数値つきで残る。

### Lean

- **F-7c を閉じた**(`formal/NSSingularity/TimeDependentGalerkin.lean`)。
  第 4 便の「mathlib の局所存在定理は自励系専用」という判定は、その*定理*に
  ついては正しいが API 全体については誤りだった。pin 済み mathlib の
  `IsPicardLindelof` は `f : ℝ → E → E` の時間依存場に対して述べられている。
  2 経路を比較し、直接経路が仮定・行数ともに厳密に少ないため採用
  (自励化経路は `E × ℝ` の instance 義務と 70 行の還元定理を要する)。
  実質は `B x x` の局所 Lipschitz 評価 `2‖B‖(‖x₀‖+a)` 1 点。
- 全 46 定理の `#print axioms` が `[propext, Classical.choice, Quot.sound]` のみ。
  `sorry`/`admit`/project 固有 axiom はゼロ。公理化していない。

### 証明書と数学

- **時空スラブ証明書**(`slab_certificate.py`): `[t_n, t_{n+1}]` 上で **cell
  内部・全時刻**を包含。Hermite 基底の値域を厳密有理数で使用
  (`h₀₀+h₀₁≡1` かつ両者 `[0,1]` なので値部分は端点の凸結合、接線基底の極値は
  `±4/27`)。flagship で 706 cell、Poisson 残差 1.57e-15、trapezoid 局所欠損
  2.35e-6、Simpson 8.25e-10。独立 checker は改竄 4 種を拒否し、壊れた入力には
  例外ではなく判定を返す。
  **仮説 H1(cell Lipschitz)と H2(Hermite 剰余)は未証明**であり、payload に
  `proved: false` として明記、checker はこれを真と主張する payload を拒否する。
- **`H^s` 誤差伝播**(`docs/research_notes/hs_error_propagation.md`): `s > 5/2`
  への移行を導出。圧力項の消去(§3.1)と粘性項の符号(§3.4)は証明済み、
  積差恒等式は Lean 既証。未証明部は **HS-1…HS-6** として文と未証明定数つきで
  列挙。**最大の欠落は HS-5**(離散残差 → `‖R‖_{H^s}`)であり、これが無い限り
  証明書は PDE についての言明にならない。`L^∞` 証明書を無条件 PDE 証明として
  提示していない旨を §7 に明記した。


## 2026-07-28 セッションの追加結果(fable5-mainline)

### 照合と統合

- ユーザ展開の Poisson バンドルは未統合(未追跡)だったため、全 6 ファイルの
  byte-identity を検証して正規パスへ統合した。パッチの二重適用はしていない。
  パッケージング残骸は `archive/poisson_gate_packaging/` に provenance 込みで
  保存。バンドル同梱の旧 snapshot 証拠は
  `outputs/poisson_gate_v1_bundle_snapshot/` として保存した。
- 統合直後の全テスト: **119 passed**(Python 3.11.9, Windows)。
- Poisson ゲート新規実行 `outputs/poisson_gate_fable5`: 全 7 受入検査合格、
  観測次数 1.9569/1.9782、manifest と全 payload の SHA-256 検証合格。
  バンドル記録と有効数字 10 桁一致(ビット単位では環境差)。

### Poisson 2 実装の相互監査と相互検証テスト

- `poisson.py`(\(r^3\)-flux 有限体積)と `finite_cylinder_poisson.py`
  (非発散形直接差分)の規約は完全一致(負作用素、軸係数 8、外側 Dirichlet
  の意味論)。数学的中核に欠陥なし。
- **独立性は部分的**: radial 離散化と Thomas 解法は真に独立、
  z 方向 Fourier 処理(`numpy.fft.fftfreq` の波数配列がビット同一)と
  `AxisymmetricGrid` は共有の単一障害点。この限定を文書と PLAN に明記した。
- `tests/test_poisson_cross_validation.py`: CV-1(同一入力での
  \(O(\Delta r^2)\) 実測一致 \(D\approx0.115\Delta r^2\)、独立性が崩れると
  発火する下限クローズ付き)、CV-2(radially exact 場での丸め一致
  \(\le4\times10^{-15}\) — FFT 規約・Nyquist・符号・軸係数のピン留め)、
  CV-3(対故障注入)、import 独立性ガード。
- 監査で検出した周辺欠陥を修正: B の複素入力の暗黙切り捨て(D2)、
  `run_poisson_gate` の config 未検証・provenance 欠如・解像度数規則
  (D3/D6/D7)と実験テスト不在(D4)、B の行列が row 1 で M-matrix で
  ないことの明文化+構造ピン留めテスト(D1)。
- この時点の全テスト: **146 passed**。

### Hou (arXiv:2107.06509) 一次資料監査

- v1・v2 の LaTeX 原文と、数値手法の委譲先 Hou–Huang (arXiv:2102.06663) を
  取得して監査した(`docs/hou_setup_audit.md`)。式 (2.1a–d) は E-11–E-14 と
  符号込みで一致。
- 新規監査エントリ E-27–E-31(壁条件 \(\psi_1=0\)・\(u_1=0\)・Thom 型
  \(\omega_1=-\psi_{1,rr}\)、半周期奇対称、初期値式 (2.2)、二段階粘性
  \(5\times10^{-4}\to5\times10^{-3}\) at \(t_0=0.00227375\)、壁渦度の
  2 次離散式)。
- 導出値 \(\|\omega(0)\|_\infty=24000\pi\cdot37^{-1/2}(36/37)^{18}
  \approx7569.62\)、\(\|u_1(0)\|_\infty\approx3265.99\)(論文は比のみ記載の
  ため再現換算に必須。数値最大化で独立検証済み)。
- v1(非爆発の主張)→ v2(potentially singular)の**結論反転**を記録。
  計算設定は両版で同一であり、変わったのは解釈のみ。判定量は後期の
  \(R/Z\) と \(\int\|\omega\|_{L^2}^4\)。
- 論文の誤植 4 箇所、取得不能事項(絶対値時系列、filter 保持の曖昧さ等)を
  記録し、FoCM 出版版の入手をユーザへ依頼した。
- **FoCM 出版版照合(ユーザ提供 PDF、同日追記)**: 計算設定は出版版で
  一切変更なし(本文数値トークン全 219 種の機械照合+項目別逐語照合)。
  実質的変更は「vacuum region」→「no-spinning region」の用語 1 件のみ。
  疑義 4 箇所を含む誤り計 9 件が出版版にも残存し、正誤表なし。
  当方の判読(\(t_3=0.0022868453\) 等)が最終根拠として確定。
  図の解像度地図(最強の爆発判定図の多くは 1024² 実行、Fig. 12 下段は
  Euler 計算で NS ターゲット不可)と、arXiv v2 PDF のベクタ軸目盛から
  回収した絶対値アンカー \(\|u\|_{L^3}\approx46.84\)–\(46.86\) を
  `docs/hou_setup_audit.md` §12 に記録。PDF はハッシュのみ記録し
  コミットしない。

### 非線形 production ソルバ(`nonlinear_cylinder.py`)

- 設計は `docs/nonlinear_solver_design.md` に確定(前提はすべて監査済み式)。
- 実装: Heun/RK2、段ごとの拘束順序(u1 壁ピン → \(\psi_1(1,z)=0\) の楕円
  solve → E-31 Thom 壁渦度 → E-14 速度回復)、壁行は代数拘束として発展
  させない、適応 CFL(0.1、軸係数 4 の拡散余裕)+`max_time_step`、
  二段階粘性 schedule、schema v2 checkpoint と restart。
- テスト 39 件+実験テスト 28 件: manufactured 空間次数(u1: 1.990/1.998、
  ψ1: 2.005/2.003、ω1: 1.845/1.902 — sup 誤差は E-31 壁行に乗り 2 へ下から
  接近)、時間次数 ≈2.0、E-31 単体 1.950/1.975、零場不動点、小振幅の線形
  拡散極限、z 奇対称保存(丸めレベル 1.2e-15 相対)、軸 parity defect の
  閉形式 \(O(\Delta r^3)\) 一致、循環・エネルギー単調性、故障注入 5 種
  (検出比 18.7〜4865)、restart 忠実性。
- **既知の前登録済み注意**: 全振幅 12000 では離散循環最大原理が
  \(O(10^{-4})\) 相対で破れる(中心差分移流+陽的 RK2 の離散化 artifact、
  細分で 2 次超で消失)。受入閾値 1e-3 は前登録値であり、超過時は
  acceptance 失敗として正直に報告される。
- この時点の全テスト: **193 passed**。

### 早期 Hou 実行(`outputs/hou_early_time_v1`)

E-29 初期値(振幅 12000)、\(\nu=5\times10^{-4}\)(E-30 第 1 段のみ、
\(t_0\) 前なので切替なし)、フル周期 \(z\)、\(t\in[0,T_1=0.002191729]\)、
CFL 0.1、`max_time_step` 1e-6。全 8 受入検査合格。manifest+全 50 payload
の SHA-256 検証合格。

| 格子 | steps | \(\min\Delta t\) | 増幅率 \(\|\omega\|_\infty/\|\omega(0)\|_\infty\) at \(T_1\) | \(\max\|u_1\|\) at \(T_1\) | argmax \(u_1\) at \(T_1\) |
|---|---:|---:|---:|---:|---|
| 65×128 | 2192 | 7.29e-7 | 6.1148 | 7605.1 | (0.0469, 0.0156) |
| 129×256 | 2192 | 7.29e-7 | 12.6957 | 14742.9 | (0.0391, 0.0117) |
| 193×384 | 2205 | 2.76e-7 | 15.6280 | 18718.1 | (0.0312, 0.0104) |

観察(すべて **numerical observation** の語彙水準):

- 増幅率は解像度で単調増加し、Hou の 1536² 適応格子公表値 20.5235 へ
  **下から接近**する。ただし隣接差(6.58, 2.93)から見かけの収束次数は
  1 未満であり、**収束していない**。外挿で 20.52 への一致を主張することは
  できない。〔改版 2026-07-28(P0-D): 旧文の「grid-scale での飽和はない」
  という表現を撤回する。増幅が plateau しないことは grid-scale saturation
  の不在を意味しない。実際ピークの radial FWHM は 65/129/193 格子で約
  4/5/6 点、ピークは軸から 3/5/6 セルであり、構造は grid-scale に近い。〕
- ごく早期の \(\|u_1\|_\infty\) 減少([Hou21, §2] の定性ターゲット)を全
  解像度で確認: 3265.6 → 最小 2011〜2025 → その後成長し \(T_1\) で初期値の
  4.7〜5.7 倍。
- \(u_1\) の argmax は軸付近・\(z\) 小の領域へ移動(原点方向への伝播と
  定性的に整合)。ただし \(T_1\) での front 位置 \(r\approx0.031\)(193 格子で
  軸から約 6 セル)であり、解像度は限界的。
- エネルギー増加 0.0(全解像度)。循環最大原理の破れは
  7.6e-4 → 2.2e-4 → 4.5e-5 と細分で減少し、前登録閾値 1e-3 以内。
- z 奇対称 defect 比 ≤ 2.0e-9(課さずに監視、保存された)。
- 独立 solver B との \(\psi_1\) cross-check 相対差:
  1.19e-2 → 7.78e-3 → 4.61e-3。急峻化する front 上では見かけの次数が
  2 を下回る(記録のみ、gate ではない)。
- E-02 発散残差最大は 1072〜1764 で単調減少せず、front の解像度不足を
  反映する。初期ノルムの E-29b 一致は 9.95e-3 / 1.98e-3 / 1.07e-3(相対)。

**ラダー延長(`outputs/hou_early_time_v2_hires`、前登録
`configs/hou_early_time_hires.json`)**: 同一プロトコルで 129×256 を再実行し
257×512 へ延長した。全 8 受入検査合格。

- **再現性**: 129×256 の増幅率は前回実行と **bit 単位で同一**
  (12.6956952437)。決定論的再現を確認。
- 257×512: 増幅率 **17.2588** at \(T_1\)(2298 steps、\(\min\Delta t\)
  2.85e-7)、\(\max|u_1|(T_1)=20306.5\)、argmax (0.0312, 0.0098)。
- 4 点ラダー 6.11 → 12.70 → 15.63 → 17.26 は公表 1536² 値 20.5235 へ
  単調接近を続ける。隣接差 6.58, 2.93, 1.63 の見かけの収束次数は依然
  1 未満であり**未収束**。外挿による一致主張は行わない。front 位置
  \(r\approx0.031\) は 257 格子で約 8 セルであり解像度限界が支配的。
  〔改版 2026-07-28(P0-D): 旧文の「grid-scale 飽和はない」という表現を
  撤回。plateau の不在は saturation の不在を意味しない。〕
- 初期ノルム E-29b との相対誤差は 1.98e-3 → 4.99e-4(比 3.97、
  きれいな 2 次接近)。solver B との \(\psi_1\) cross-check 相対差は
  7.78e-3 → 2.95e-3。
- 実効 advective CFL の設定値 0.1 に対する微小超過(+0.16%)を 257 でも
  記録(原因調査と受入条件の明文化は本セッションの別タスク)。

### 拘束誤差の相対化と CFL 受入条件(タスク 1+5)

絶対残差だけでは場の増幅と数値破綻を区別できないため、TM-09 に従い
全拘束診断を相対化した(記録のみ、既存 gate は不変)。

- **相対発散**: E-02 残差最大 / 打消し項和最大
  \(\max(|\partial_ru^r|+|u^r/r|+|\partial_zu^z|)\)。分母・argmax 位置・
  各点比も保存。実測(振幅 12000、短時間 run):
  33×64 で 2.92e-2 → 65×128 で **8.04e-3**(≈O(h²) 改善)。
- **相対軸パリティ**: 軸片側微分 defect / \(\max|\partial_r\) 場\(|\)。
  実測: u1 7.48e-3 → **9.39e-4**、ω1 → 2.75e-3(≈O(h³) 改善)。
- いずれも**離散化オーダーで減少しており、増幅由来の破綻ではない**。
- **CFL 超過の原因確定**: \(\Delta t\) は step 開始時点の速度最大で決定し、
  実効 CFL は step 後の状態で評価するため、超過分は step 内成長率に
  厳密に一致する。v1 実測 +0.227%(193×384、拘束が効き始めた step 2085
  以降)。HH21 の <1%/step 指針内。受入条件を
  \(\text{CFL}\le C(1+\varepsilon)\)、\(\varepsilon=0.05\)(config で変更可)
  として明文化し、\(\varepsilon=0\) では v1 値が不合格になることを
  テストで固定した。

### 保存 snapshot の独立 Cartesian 検査(タスク 2、`outputs/hou_snapshot_cartesian_audit_v1`)

円柱演算子を import しない(AST テストで固定)独立経路で、checkpoint を
一様 Cartesian 箱 \([-0.7,0.7]^2\times[0,1)\) へ復元し検査した。
全 7 hard gate 合格、manifest+payload SHA-256 検証済み。

- t=0: 相対発散 RMS 1.5e-4〜2.5e-4、curl RMS ≈4.6e-4(離散化レベル)。
- \(T_1\): RMS はゲート合格(div 6.1e-4〜9.2e-4、curl 4.9e-3〜5.9e-3)、
  監査格子細分で ≈2 次減少。ただし **pointwise 最大は front 近傍で悪い**
  (curl defect 最大 = 勾配スケールの 0.41、方位一致最大/信号 = 0.62、
  nr193)。監査格子(dx=0.022)が nr193 の source 格子より 4.3 倍粗い
  ことと整合するが、このレベルでの snapshot 欠陥を排除するものではない
  (summary に明記、gate は RMS のみ)。
- 故障注入: u_y 符号反転(×52)、成分入替(×52)、E-18b 符号、
  軸条件違反(奇 r kink、×45; 軸正則な整合対照は不変)を検出。
  発散検査は E-18a 恒等性より ψ1 編集に構造的に盲目であることも記録。
- **primitive NS 残差**: checkpoint 対+圧力から組み立てる関数を実装し
  manufactured 場で時間次数 2.004/2.001 を検証。実 checkpoint への適用は
  「圧力未保存・snapshot 間隔 ≈500 step」のため未実施
  (`primitive_navier_stokes_residual_gap` として全 summary に明記)。

### Poisson 第三経路 solver C(`realspace_poisson.py`)

A/B が共有する axial Fourier(`fftfreq` bit 同一)+Thomas の単一障害点を
破る実空間経路: E-26a の \(r^3\)-flux 行(独立転写)+周期 2 次 z 差分
(`np.roll`)+\(V_i\) 重み付き SPD の Jacobi 前処理 CG。`numpy.fft`
不使用をソース文字列レベルでテスト固定。

- C-vs-A 差 = **0.0901 Δz²**(純 axial 離散化ギャップ; 下限 90 倍・
  上限 11 倍のマージンで帯域拘束。経路統合が起きると下限が発火)。
- 重み付き対称性 1.9e-16、軸行の基底ベクトルプローブは厳密に 8/Δr²。
- 発見 2 点を文書化: (1) \(V_0=\Delta r^4/64\) のため重み付き CG 残差は
  軸セルをほぼ無視 → 非重み付き代数残差も併record・gate;
  (2) 大域誤差ゲートは軸係数 12/Δr² への破壊に盲目(誤差変化 0.4%)
  → 直接プローブが必須。
- 範囲: C は A と radial 数学を共有するため A-vs-C は axial 経路のみを
  検証。radial 独立性は B の役割。grid class と binary64 は 3 経路共有。

### 時間刻み収束(タスク 3、`outputs/hou_time_refinement_v1`)

同一格子 65×128・同一終了時刻 \(T_1\)、固定
\(\Delta t=6\times10^{-7},3\times10^{-7},1.5\times10^{-7}\)
(step 数 3653/7306/14612、実効刻みは厳密に半減)。
全 9 受入検査合格、manifest+全 14 payload SHA-256 検証済み。

| 量 | 観測時間次数 | step-doubling 差(粗→中/中→細) |
|---|---:|---|
| 増幅率 | **1.998** | 2.67e-5 / 6.69e-6 |
| energy | **2.034** | 1.23e-5 / 3.01e-6 |
| \(\max|\omega_1|\) | **1.999** | 0.197 / 0.0493(値 5.29e5 に対し相対 ~1e-7) |
| \(\max|u_1|\) | ≈2.0 | 7.55e-3 / 1.89e-3 |

- argmax 位置(u1・Cartesian 渦度)は 3 水準で**格子点単位で同一**。
- **時間/空間分離**: 拘束系指標(相対発散 4.245e-2、相対軸パリティ
  u1 0.706、循環違反 7.664e-4、solver-B cross-check 1.193e-2)は
  3 水準で相対スプレッド <2e-5 の**dt 非依存** — すべて空間離散化
  支配であることを直接確認。増幅率の時間誤差(最細 ~6.7e-6)対
  空間ラダー差(1.63〜6.58)の比は \(10^{5{-}6}\) で
  `dominant_error_source = "spatial"`。
- 固定 dt の実効 CFL 最大 0.0149(適応拘束の余裕内)。
- **注記**: \(T_1\) での相対軸パリティ 0.706(u1、65×128)は軸近傍
  front(\(r\approx0.047\)= 3 セル)の空間的未解像を相対指標が正しく
  可視化したもの(dt 非依存が示すとおり時間積分の破綻ではない)。
  Cartesian 監査の pointwise 所見と整合し、65×128 は \(T_1\) 近傍で
  軸近傍が深刻に解像不足であることを定量化する。

### 圧力非依存の primitive 残差(ギャップ解消、`outputs/hou_primitive_residual_v1`)

Cartesian 監査で「圧力未保存・snapshot 間隔 ≈500 step」として明示していた
ギャップを、**圧力を発明せずに**閉じた。運動量残差
\(R=u_t+(u\cdot\nabla)u-\nu\Delta u\) は厳密解で \(-\nabla p\) に等しいので
\(\operatorname{curl}R=0\)、同値に渦度輸送残差 \(S\) が消える。両者は
圧力を含まない。

- 実験は 7 状態(offset \(0,\pm1,\pm2,\pm4\) step)を保存するため、中心差分
  \(u_t\) の Richardson **次数が測定できる**(誤差推定ではない):
  実測 **1.99997**。
- 評価は保存 artifact を `load_candidate` で読み直す経路のみ(非循環性)、
  stencil は `cartesian_validation` 固有、円柱 module 不 import を AST 固定。
- 129×256 実測(内部領域、分母は**同一領域**の項和最大 — 甘くない正規化):
  相対 RMS \(\operatorname{curl}R=5.835\times10^{-4}\)、
  \(S=9.346\times10^{-4}\)、\(\nabla\cdot u=2.815\times10^{-4}\)。
  監査格子細分の次数 1.885 / 1.835 / 2.279。
- 全 10 受入検査+5 record-only 合格、manifest+26 payload SHA-256 検証済み。
- **正直に記録した所見**: (1) \(\nu=5\times10^{-4}\) では粘性項が項和の
  \(2.9\times10^{-4}\) しかなく、実データ上で粘性符号反転は**検出できない**
  (その旨を主張するテストを置き、決定的な注入は manufactured 場で実施);
  (2) 減衰 Taylor–Green は移流が厳密勾配のため非線形項に対して退化 →
  別の非自明 solenoidal 場を追加; (3) 2 形式は実データで相対 RMS 2.4% 相違、
  次数 1.28(バンドル中最弱); (4) 65×128 は次数 0.79/1.13 で bilinear 復元
  律速のため 129×256 を出荷; (5) 既存 `hou_early_time` checkpoint は
  この観点では**未監査のまま**。
- 圧力回復は scoped・record-only(転置を随伴恒等式 \(3\times10^{-15}\) で
  検証、零空間を射影除去。\(x,y\) 端の閉じ方により診断であって検証済み
  圧力ではない)。

### Lean 4 段階 1 — F-3(`formal/NSSingularity/VelocityRecovery.lean`)

E-14 ⇒ E-15(速度回復から物理発散ゼロ)を機械検証した。

- `divergence_of_recovered_velocity_eq_zero`(\(r\neq0\))と
  `divergence_of_recovered_velocity_eq_zero'`(\(u^r/r\) を連続延長に
  置換、軸込み)、`mixed_partial_comm`(Schwarz)。
- `sorry`・`admit`・新規 axiom なし。`#print axioms` は
  `[propext, Classical.choice, Quot.sound]` のみ(独立に再確認)。
- 独立検証: `lake build` 成功(8658 jobs)、軸の向き(`partialR` が
  \(r\)、`partialZ` が \(z\))を別途 Lean で確認 — 軸取り違えで別命題に
  なる罠を排除。
- 形式化により、E-15 の相殺が**混合偏微分の一致のみ**に依存することが
  構造的に露出した(滑らかさ仮定の唯一の使用箇所)。
- **非スコープ**: 任意の \(C^2\) スカラーに対する座標表現の恒等式であり、
  E-13 も運動方程式も使わない。E-18/E-24(3D 場との対応)は未形式化で、
  `ClayStatement.lean` の `DivergenceFree` へは未接続。

### 壁依存性(E-32 + 実装完了、実行中)

- **E-32**(`docs/equation_audit.md`): \(C^\infty\) compact-support
  envelope 初期値族。core(\(r\le0.9\))は E-29 と **bit 一致**、
  sup 偏差 \(\le3.4008\times10^{-10}\)(実測 \(2.577\times10^{-12}\))、
  \(r\ge0.95\) で厳密 0、遷移帯の 4 階差分は平の E-29 と bit 一致。
- 実装(`wall_dependence.py` + 実験 + 60 テスト)。`nonlinear_cylinder` が
  \(r_{\max}\) について既に一般であることを検証(4 半径で壁行の挙動が正しく、
  初期エネルギーは compact support ゆえ全半径同一 4008.106)。
- **前登録の完全性**: 閾値は一つも変更していない。受入検査 3 件の文言が
  一意でなかったため、実装した読み方を
  `docs/wall_dependence_prereg.md` §8 に**改版として明記**し、
  前登録されていなかった実装判断(整合性許容幅 0.25 等)も別枠で列挙した。

### 壁依存性の実行結果(`outputs/wall_dependence_v1`)

6 member(主群 \(\Delta r=1/192\) の \(R_{\mathrm{wall}}=1,1.5,2,3\)、
粗群 \(\Delta r=1/128\) の \(R=1,2\))、\(\nu=5\times10^{-4}\)、
\(t\in[0,T_1]\)。全 12 受入検査合格、manifest+97 payload SHA-256 検証済み。

| \(R_{\mathrm{wall}}\) | \(n_r\) | core 増幅率 \(A(T_1)\) | \(R=3\) との差 |
|---:|---:|---:|---:|
| 1.0 | 193 | 15.627954940635 | 2.4915e-3 |
| 1.5 | 289 | 15.630441984776 | 4.4673e-6 |
| 2.0 | 385 | 15.630446443875 | 8.1784e-9 |
| 3.0 | 577 | 15.630446452053 | — |

- 隣接対の分離 \(S\): **1.594e-4 → 2.853e-7 → 5.232e-10**、いずれも前登録
  閾値 0.05 を大きく下回る。argmax 変位は**厳密に 0**(全半径で同一格子点)。
- **§4 literal の判定は `wall_effect_small`**。ただし保守的な §2 読み
  (`classification_with_section_2_hold`)は **`undecided`**(下記の
  解像度整合性が不安定なため)。
- 楕円非局所寄与(core 上の \(\max|\Delta\psi_1|\))も同じ比で減衰:
  相対 2.968e-4 → 5.321e-7 → 9.760e-10。argmax は **\(r=0.8958\)**、
  すなわち core 外縁 \(r\le0.9\) の際 — 壁の影響が core 内で最大になる
  場所として整合的。
- **交差検証**: \(R=1\) の envelope member の増幅率
  15.627954940635405 は、平の E-29 初期値による `hou_early_time_v1` の
  193×384 の値と**完全一致**。E-32 性質 3・4(core bit 一致・ノルム不変)が
  力学レベルでも成り立つことの確認。

**機構の独立検証(本オーケストレータによる)。** 上表の差は
\(\Delta R=0.5\) ごとに **557.7 倍 / 546.2 倍**で減衰する。これは
\(z\)-Fourier モード \(k\) における \(\mathcal L_5\) の径方向解
\(r^{-1}I_1(kr),\ r^{-1}K_1(kr)\) に対する Dirichlet 壁の像応答
\(K_1(kR)/I_1(kR)\sim\pi e^{-2kR}(1+\tfrac3{4kR})\) が予言する比である。
\(z\) 周期 1 かつ \(\omega_1\) が \(z\) 奇なので \(k=0\) は存在せず最低
モードは \(k=2\pi\)。予言値は **555.2 / 545.5**(先頭項のみなら 535.5)で、
実測との差は **0.45% / 0.12%**。壁効果の正体は最低軸方向モードの
指数的像応答であると同定した。

**この結果が意味しないこと(最重要)。**

1. **全空間の壁独立性を示していない。** 指数減衰は \(z\) 周期 1 が
   \(k\ge2\pi>0\) を強制することの帰結である。\(\mathbb R^3\) では軸方向
   波数は \(k\to0\) まで連続で、楕円減衰は指数的ではなく代数的になる。
   本実験は「軸周期を固定したまま半径を後退させた」測定であり、
   前登録 §6 が明記するとおり壁を無限遠へ送る極限ではない。
2. **\(k=0\) モードの不在は Hou 初期値の奇対称性に依存する。** 一般の
   データでは代数減衰する単極子成分が残り、この指数的鈍感さは失われる。
3. **`wall_effect_small` は Clay 候補の証拠ではない。** summary の
   `interpretation` にも同文を記録した。増幅率自体は依然として解像度
   未収束(6.11→12.70→15.63→17.26)であり、壁を退けても
   軸近傍の未解像(相対軸パリティ 0.706 at 65×128)は解消しない。
4. **解像度整合性は不安定。** 共有半径対 \((1,2)\) の \(S\) は
   \(\Delta r=1/192\) で 1.594e-4、\(\Delta r=1/128\) で 2.504e-4 と
   **57.1% 変化**し、(前登録外の)許容幅 0.25 を超える。両値とも閾値
   0.05 を遥かに下回るため「既に無視できる量の 3 桁目が未収束」という
   状況だが、保守的読みでは判定を保留する。
5. 早期区間 \([0,T_1]\) のみ。粘性切替後・中成長段の壁依存性は未測定。

### E-33: 壁打切り応答の閉形式(`outputs/wall_truncation_scaling_v1`)

壁依存性の**結果を説明する**ために導出・監査・検証した。
\(\psi=\varphi/r\) の代入で \(\mathcal L_5\) は各 \(z\)-Fourier モードで
厳密に 1 位の変形 Bessel 方程式になり、compact support の \(\omega_1\) に
対する壁誘起誤差は core 上で厳密に \(K_1(kR)/I_1(kR)\) に比例する。

- **二つの漸近域**: \(kR\gg1\) で指数 \(\pi e^{-2kR}\)、\(kR\ll1\) で
  **代数 \(R^{-2}\)**、交差は \(kR\sim1\)。
- **\(k=0\) の厳密閉形式**: 壁依存部分は core 上で定数
  \(-Q_\infty/(2R^2)\)、\(Q_\infty=\int_0^\infty s^3\omega_1\,ds\)。
- 独立 oracle `bessel_reference.py`(scipy 非依存、級数 \(I_\nu\)+積分表示
  \(K_\nu\)、Wronskian 恒等式を \(6.7\times10^{-16}\) で検証、求積は
  spectral)に対し、全 23 受入検査合格。測定と oracle の一致は最細格子で
  \(4.2\times10^{-4}\)(次数 1.997–2.002)、形状の \(I_1(kr)/r\) 相関ずれ
  \(5.5\times10^{-8}\)(次数 3.98–4.02)、\(k=0\) 閉形式は次数 2.0000。
- **非線形との連結**: 壁依存性実行の分離比 557.746 / 545.231 は楕円 oracle
  と \(1.2\times10^{-4}\) / \(2.8\times10^{-4}\) で一致。**非線形力学の壁
  感受性は線形楕円 Bessel 機構そのもの**である。
- 記録した整合性訂正: 初期の探索値は参照壁の規約を明記していなかった。
  検証は参照 \(R=3\) の規約で厳密に再現し、最大半径参照の正しい値
  (5.939/5.306、3.242/2.711、4.5006)を E-33 に記載した。
- 丸め床の扱いを first-class 化(応答が \(10^{-19}\) を切る半径は oracle
  から決定論的に分類し、測定値で比を取らない)。

**プログラム上の帰結。** Hou 設定は \(kR=2\pi\approx6.28\) と指数域の奥に
あるため壁が効かない。しかし \(\mathbb R^3\) では \(k\to0\) まで連続で、
長波長成分の打切り誤差は \(R^{-2}\) でしか減らない
(\(\varepsilon=10^{-6}\) に core 半径の \(10^3\) 倍)。**大半径 Dirichlet
円柱は全空間近似として代数的にしか収束せず実用的でない。**
同じ modal 構造から**厳密な透過条件**
\(\partial_r\hat\psi_k(R)+[2/R+kK_0(kR)/K_1(kR)]\hat\psi_k(R)=0\) が導け、
compact support の下で打切り誤差が厳密に消える。設計と受入条件を
`docs/whole_space_transition.md` に実装前固定した。

### W-A: 透過境界条件ソルバ(`outputs/transparent_boundary_v1`)

実装 `transparent_boundary.py`、全 25 受入検査合格(3.5 秒)、
W-B の受入条件 6 件すべて合格。既定は全経路で Dirichlet のまま。

- **\(k=0\)**: core の \(R=1\to2\) 差は透過で \(9.07\text{e-}7\to
  2.27\text{e-}7\to5.66\text{e-}8\)(次数 2.0002)。Dirichlet は
  \(-1.697\text{e-}3\) で不変。改善率 1872 → 7489 → **29959**。
  E-33e の offset が消え、2 次で 0 へ収束する量に置き換わった。
- **\(k>0\)**: 改善率は \(L_z=1\) で ×11013、\(L_z=8\) で ×31787、
  長波長 \(L_z=32\)(代数域)で ×30510、最大 ×34231。全 24 の
  \(R\) 非依存性の観測次数は [1.948, 2.001]。
- manufactured 収束は全行 1.9975–2.0038、**境界行のみ 2.0099–2.0150**。
  境界行の局所打切りは \(O(\Delta r)\) だが M-matrix 性により全体 2 次
  (仮定ではなく測定)。
- 故障注入: 符号反転 3628–4152 倍、\(2/R\) 脱落はソルバ拒否、
  \(K_0/K_1\!:=\!1\) は 22.5–186 倍。support 漏れは 2 経路で拒否。
- 正直に記録した限界: `frozen_ratio` は \(k=0\) で bit 単位不可視、
  条件 4 の \(L_z=32\) は比較対象の大半径 Dirichlet 自体が未収束
  (E-33 の代数 \(R^{-2}\) そのもの)、\(K_\nu\) 求積は \(x>10^3\) で
  文書化した漸近形へ切替。

限界: 一様固定格子は Hou の適応最小格子幅 \(O(10^{-8})\) に遠く及ばず、
これは公表増幅率の再現主張ではない。FABLE5_HANDOFF §7.2 の受入基準
(3 空間解像度+3 時間刻み+peak 位置・振幅傾向の一致等)のうち、
時間刻み系列と published-diagnostic 比較の大部分は未実施である。
中成長段(§7.1 stage 2)へ進む前に、より高解像度または適応格子・
半周期 sine 実装の設計判断を行う(「次に行うべき最小の一手」参照)。

一様 \(x,y,z\) 格子上の独立な3成分発散、full curl、vector Laplacian、
primitive PDE項別残差を追加した。保存候補のchecksum検証付き再読込から、
既存円柱差分を呼ばないadapterを経てCartesian物理場を検査するend-to-end
経路も実装した。非特異対照では空間格子を固定し、時間刻みだけを
\(\Delta t,\Delta t/2,\Delta t/4\) とする収束系列を保存した。

このリポジトリは有限時間特異点を発見・証明しておらず、Navier–Stokes
ミレニアム問題を解決していない。

## 2026-07-28 第 2 セッション: FABLE5_NEXT_TASK_AUDIT の P0/P1 ゲート

指示書 `FABLE5_NEXT_TASK_AUDIT.md`(リポジトリルート、コミット済み)の
記載順に実施した。開始時の確認: 作業ツリーはリモートと一致(2cb8e48)、
全 565 テスト合格、`lake build` 成功(8659 jobs)、`formal` 内に
sorry/admit/新規 axiom ゼロ。

### P0-A: von Neumann 安定性監査(`von_neumann.py`、11 tests)

- Heun+中心差分は純移流で厳密増幅(\(|G(i\alpha)|^2=1+\alpha^4/4\))。
  凍結係数 advection–diffusion の全波数 scan、予測子段の別評価、
  snapshot 監査 API、独立参照 propagator(シンボル経路と実配列経路の
  相互一致 3.7e-16)を実装した。
- **出荷済み運転点の判定は「stability-unverified」**: 出荷 v1 の記録
  (min dt=2.76e-7、max advective CFL=0.10023)と自己整合な読みで
  Heun worst \(|G|\) は 1.0000035(radial 支配)/ 1.0(axial 支配)/
  1.000152(両方向同時)— tolerance 1e-12 では不合格。全波数 pass に
  必要な dt は約 2.4e-8(出荷値の 1/11.5)。粘性 5e-4 の寄与は
  増幅率 −6.1e-5 に留まる。詳細は `docs/numerical_stability_audit.md`。
- 過去の Heun 実行はすべて stability-unverified に再分類。**Heun 単独の
  増幅を候補判定に使うことを禁止**(決定規則として文書化)。

### P0-A: 交差検証積分器(SSPRK3 / RK4)

`take_step` に SSPRK3・古典 RK4 を追加(空間離散化・拘束順序・楕円
solve は Heun と完全共有)。毎段 `constrain_state` 射影込みの実測時間
次数: **Heun 1.97/2.00、SSPRK3 3.00/3.00、RK4 3.95/3.98** — 射影は
観測次数を落とさなかった。ゼロ場不動点・小 dt 相互一致・誤差の大小
関係(rk4 < ssprk3 < heun)をテストで固定。

### P0-B/P0-C: 全 step streaming gate(`test_integrator_gates.py`、20 tests)

- `IntegrationResult` に `step_stream`(全 accepted step × 28 量)と
  `gate_summary`(streaming 極値)を追加。エネルギー増分、循環 defect、
  奇対称相対、軸 parity 相対、相対発散、壁拘束、**solver A が実際に
  解いた線形系の代数残差**(全 step 相対 <1e-12)、pre/predictor/post
  CFL、粘性安定数、dt を縛った拘束の名前、エネルギー収支 defect を
  全 step 保存。出力間引きは gate に作用しない。
- 中間段 CFL 超過での step 棄却(`stage_cfl_limit`、dt 半減再試行)を
  実装、受入 step の段 CFL が閾値以下であることをテストで固定。
- **義務の合成違反テスト**: 記録行間だけの強制パルス注入→抽出で、
  間引き history は単調減衰のみを示すが streaming は丸め床の
  3×10¹² 倍の増加を捕捉する。

### P1-C: エネルギー収支と viscosity_sign fault

- E-27 壁は swirl のみ no-slip(\(u^z(1,z)\ne0\) の滑り壁)なので、
  正しい恒等式は \(dE/dt=-\nu\int|\omega|^2dV-\nu\oint u^z\omega^\theta dS\)。
  壁項込み/なしの両 defect と swirl エネルギーの項別仕事率
  (移流・stretching・粘性)を全 step 記録。
- **整合性の記録**: 初回実装は \(\int|\omega|^2dV\)(\(2\pi\) 測度)を
  E-20 の \(\pi\) 正規化 enstrophy と取り違え、相対 defect が理論値
  どおり厳密に 0.5 へ飽和(fault 時 1.5)。これを因子 2 の欠陥として
  特定・修正した。修正後: 滑らか control で 5.9e-2 → 1.6e-2 → 4.2e-3
  (空間時間同時細分で収束)。
- 新 fault `viscosity_sign`: Hou 運転点(ν=5e-4)では項和の 3e-4 で
  不可視(既知)だが、拡散支配 control(ν=2e-2)で相対 defect
  clean 2.21e-2 vs 反転 **2.000**(理論値 2、比 90 倍)、エネルギー
  単調性反転により確実に棄却される。

### P0-D: core-width / points-per-scale(`core_width.py`、15 tests)

- radial/axial FWHM、10–90% front thickness、subgrid 二次ピーク、
  勾配長スケール、高周波 tail(z rfft / r DCT-II)、共通格子への
  Catmull-Rom 補間比較、manufactured tanh front 研究を実装。
- **前登録閾値 `PREREGISTERED_MIN_POINTS_PER_FRONT = 7`**(tanh front
  で 10–90 幅の相対誤差 ≤2% となる最小整数点数。テストが研究を再計算し
  定数との整合を強制するため事後調整不可)。
- 出荷済み 193×384 の T₁ snapshot: points_per_fwhm_r 6.92、
  **points_per_front 4.36**、ピークは軸から 6 セル、勾配スケールは
  0.49dr / 0.98dz(1 セル未満)→ `fit_precondition` は**不合格**
  (3 理由)。「grid-scale 飽和なし」という旧表現は撤回した(本文中に
  改版注記)。**収束 fit は現データでは禁止**が機械化された。

### P1-A: blind 外挿(`extrapolation.py`、10 tests)

- 3 点厳密解・最小二乗 power law・固定次数 Richardson・全部分列感度。
  署名は \((h, A)\) のみで外部 anchor を構造的に受け付けない
  (`20.5235`・`20.52` の不在をテストが assert)。
- 実ラダー (6.11, 12.70, 15.63, 17.26): 指示書の見積(全4点 27.38/0.54、
  前 3 点 28.85/0.49、後 3 点 24.60/0.70)を独立実装+brute force で
  すべて再現。A_inf 散らばり相対 **0.488**(前登録閾値 0.05 の 9.8 倍)、
  p 幅 0.203 → **判定 not_in_asymptotic_range: 極限値は一切引用不可**。
- 記録した偶然: 次数 1 固定 Richardson は 20.528(公表 20.5235 の
  0.02% 以内)を与えるが、これは相互に矛盾する 7 外挿値の一つに
  すぎず(次数 2 固定は 16.844)、確認としては読めない。テスト
  docstring に非結果として明記。

### Poisson 第三経路の SPOF 破壊(P0 §5、`test_realspace_poisson.py` 追加)

solver A/B が共有する Fourier 機構の故障 3 種(モード正規化 slip、
周期 seam ずれ、Nyquist 混入)を solver A 解の事後改変として合成し、
実空間第三経路 C が **>10×**(実測 48× 等)で全て検出することを固定。
A/B の一致を「continuum 精度」「完全独立」と呼ばない limitation は従来
どおり。

### P0-E: 語彙修正と Gate 4 仕様

- `docs/whole_space_transition.md` §0 を新設: W-A と壁依存性の全結果を
  「**periodic-z radial-wall sensitivity observation**」と固定し、
  「whole-space validation」「R³ wall independence」の表現を禁止。
- 同 §7 に真の全空間移行 gate(非周期 z、z 方向 C∞ compact 台、
  free-space 楕円経路、R_max/Z_max 独立拡大、低波数 stress、
  Cartesian 体積測度での有限エネルギー直接検査、有限円柱解との
  同一視禁止)を**未実装の仕様**として定義。
- STATUS/設計文書の「grid-scale 飽和なし」表現を撤回(改版注記付き)。

### Gate 1 実行結果(`outputs/integrator_comparison_v1`、全 8 受入合格)

65×128、E-29 datum(振幅 12000)、\(\nu=5\times10^{-4}\)、\(T_1\) まで、
固定 dt ∈ {6e-7, 3e-7} × {heun, ssprk3, rk4}(空間離散化・拘束順序・
楕円 solve は完全共有)。前登録許容: 増幅相対差 1e-3、対差の dt 縮小、
argmax 1 セル以内。

| dt | 対 | 増幅相対差 | 場 L∞ 相対差 | argmax 差 |
|---|---|---:|---:|---|
| 6e-7 | heun vs ssprk3 | 5.823e-6 | 4.39e-6 | 0 セル |
| 6e-7 | heun vs rk4 | 5.828e-6 | 4.39e-6 | 0 セル |
| 6e-7 | ssprk3 vs rk4 | **5.09e-9** | 4.40e-9 | 0 セル |
| 3e-7 | heun vs ssprk3 | 1.458e-6 | 1.10e-6 | 0 セル |
| 3e-7 | heun vs rk4 | 1.459e-6 | 1.10e-6 | 0 セル |
| 3e-7 | ssprk3 vs rk4 | **6.35e-10** | 5.50e-10 | 0 セル |

- heun と高次法の差は dt 半減で正確に 1/4(heun の \(O(dt^2)\) 誤差
  そのもの)。ssprk3 と rk4 は互いに 5e-9 で一致。**von Neumann worst-case
  上界(累積 ~0.8%)は実測では発現せず、時間積分スキーム依存は
  ~6 ppm** — 現在の増幅値に対する時間離散化リスクの実測上界。
  (これは凍結係数監査の「stability-unverified」分類を置き換えるもの
  ではなく、補完する経験的上界である。固定 dt=6e-7 の段 CFL は最大
  0.0149 で、適応 run の 0.1 よりはるかに安定側にある点にも注意。)
- heun dt=6e-7 の A_grid=6.1147053247 は `hou_time_refinement_v1` の
  粗 dt 値と**ビット単位一致**(実験間決定論再現)。
- **P1-B 両正規化の併記**: A_grid(離散初期最大で正規化)6.1147 に対し
  A_common(連続参照 \(24000\pi/\sqrt{37}(36/37)^{18}=7569.6227\) で
  正規化)**6.0539** — 65×128 の離散初期最大は連続値より約 1% 低い。
  分母の格子依存だけで増幅率が 1% 動くことが定量化された。
- **正直な記録**: 全 step streaming の相対エネルギー収支 defect 最大は
  **0.936**(全 6 run とも)。滑らか control では 4.2e-3 まで収束する
  計装なので、これは E-29 front(65×128 で FWHM 約 4 点)の空間
  未解像がエネルギー恒等式を閉じさせないという**解像度の言明**である。
- エネルギーは**全 accepted step で単調非増加**(max step increase < 0。
  旧「snapshot 対でゼロ」より強い)。循環 defect 7.67e-4(前登録閾値
  1e-3 以内)、Poisson 代数残差は全 step 相対 5e-15 以下。

### 出荷済み証拠への von Neumann 監査適用(`outputs/von_neumann_audit_v1`)

出荷済み 3 bundle の全診断行(1470 行、stride 25、step-0 除外 8)を
`audit_snapshot` で監査(721² scan × 2940 回、51.8 秒)。

| run | 行数 | 不合格行 | worst Heun max\|G\| | 判定 |
|---|---:|---:|---|---|
| v1 65×128 | 88 | 0 | 1.0(厳密) | verified-at-recorded-rows |
| v1 129×256 | 88 | 0 | 1.0(厳密) | verified-at-recorded-rows |
| **v1 193×384** | 89 | **4** | **1.0000312** | **stability-unverified** |
| v2 129×256 | 88 | 0 | 1.0(厳密) | verified-at-recorded-rows |
| v2 257×512 | 92 | 0 | 1.0(厳密) | verified-at-recorded-rows |
| refinement dt=6e-7/3e-7/1.5e-7 | 147/293/585 | 0 | 1.0(厳密) | verified-at-recorded-rows |

- 不合格 4 行は 193×384 の \(T_1\) 直前(t≈2.12–2.19e-3、CFL_z≈0.1002)。
  strided 外挿値 1.0016 は「bound ではない」と明示ラベル付き。
- Euler 予測子段は全 run で 1 を超える(最大 1.0176)— 記録のみ。
  完成 step の Heun のみを gate する。
- 制約: stride 25 の記録行のみの被覆。行間の 24 step は未監査
  (将来 run は `step_stream` で全 step 被覆)。判定語彙は
  「stability-unverified であって不安定ではない」を全箇所で維持。

### 既存 snapshot の core-width / P1-B 再正規化(`outputs/core_width_audit_v1`)

全 4 解像度 × 5 時刻 × 2 場(u1、\(|\omega|\))の points-per-scale 監査
(入力 manifest 検証済み、v1/v2 の共有 129×256 は byte 一致を確認)。

- **\(T_1\) の fit 前提は全解像度・両場で不合格**: points_per_front は
  u1 で 2.57/3.36/4.36/5.43、\(|\omega|\) で 2.55/3.29/4.43/5.70
  (閾値 7)。**T₁ 増幅ラダーの収束 fit 禁止が機械的に確定**
  (`convergence_fit_precondition_satisfied_at_final_snapshot = False`)。
- 計算された正直な例外: t=5e-4 の \(|\omega|\) ラダーのみ全解像度で
  前提を満たす(front 7.06/14.0/21.0/28.0 点)。早期時刻の front は
  まだ広いという整合的な結果。
- **P1-B 表(離散初期最大 vs 連続参照 7569.6226982)**:

  | nr | 離散初期max | a/b | A_grid | A_common |
  |---|---:|---:|---:|---:|
  | 65 | 7494.31 | 0.9901 | 6.1148 | **6.0539** |
  | 129 | 7554.61 | 0.9980 | 12.6957 | **12.6705** |
  | 193 | 7561.48 | 0.9989 | 15.6280 | **15.6112** |
  | 257 | 7565.84 | 0.9995 | 17.2588 | **17.2502** |

  分母の格子依存だけで最大 1% 動く。以後の主比較は絶対値と
  A_common を用い、A_grid は補助値とする(P1-B)。
- 連続初期最大位置 \(r^*=1/\sqrt{37}\) を数値最大化で 1.2e-16 まで再現。
  初期最大位置の格子誤差は最大 0.48 セル。
- 隣接解像度の共通格子差(\(T_1\)、\(|\omega|\) L∞): 6.27e4 → 3.38e4 →
  1.97e4 と値は減少するが、**微分 L∞ は 4.5e6/5.9e6/5.0e6 と減少しない**
  — 微分レベルの収束は現ラダーに存在しない(正直に記録)。

### Lean 監査(P0 §7)

`formal/AxiomAudit.lean` を追加し `lake env lean AxiomAudit.lean` を実行:
**9 定理すべて `[propext, Classical.choice, Quot.sound]` のみ**に依存
(記録は `docs/formalization_map.md`)。sorry/admit/新規 axiom ゼロ、
toolchain/mathlib は v4.32.1 固定を再確認。「8659 jobs = 8659 定理」
という読みの禁止を明文化。

## 2026-07-29: 外部セッション成果物の統合(臨界 \(L^3\)・Type-II・自由空間 Poisson)

外部セッション(ChatGPT)が作成した 3 バンドル(ZIP+patch+文書+参照出力)を
`fable5-mainline` へ統合した。**統合作業のみで、新規計算・設計変更・
数学的内容の改変は行っていない。** 3 つの patch はいずれも現行ツリーへ
競合なく適用でき、ZIP 内実装と patch は同一内容(ハッシュ照合済み)。

統合物: `critical_l3.py` / `scaling_constraints.py` / `scaling_fit.py` /
`free_space_poisson.py` / `wall_sensitivity.py`、対応する 5 テストファイル、
5 実験スクリプト、6 文書、参照出力。統合後の全テスト
**747 passed, 1 skipped**(統合前 711 passed。skip は scipy backend の
照合テストで、scipy はこのリポジトリの依存に含まれない — 下記のとおり
scipy-free オラクルによる同等の照合を別途追加した)。

### 数学的に証明されたこと(紙上の導出。Lean 未形式化)

- **臨界 \(L^3\) 障害**(`docs/research_notes/critical_l3_obstruction.md`):
  等方的動的再スケーリング \(u(x,t)=L^{-1}U((x-x_*)/L,s)\) の下で
  \(\|u(t)\|_{L^3}=\|U(s(t))\|_{L^3}\) は**厳密な恒等式**(変数変換)。
  よって \(\sup_s\|U(s)\|_{L^3}<\infty\) なら
  \(u\in L^\infty_tL^3_x\) となり、Escauriaza–Seregin–Šverák の端点正則性
  定理により \(T\) は特異時刻ではない。
  **帰結(除外定理)**: 一様に \(L^3\) 有界な一スケール再スケーリング候補
  (定常 profile、周期軌道、準周期軌道、有界な Type-I core)は
  \(\mathbb R^3\) の有限時間爆発を与えない。これは**爆発の構成ではなく
  探索クラスの除外**である。
- **異方的版**: \(\|u(t)\|_{L^3}^3=A^3L_r^2L_z\|U(s)\|_{L^3}^3\)。
  したがって \(U\) が \(L^3\) 有界なら、特異点の必要条件は
  \(A^3L_r^2L_z\to\infty\)。標準等方放物型 \(A=L^{-1},L_r=L_z=L\) では
  この積は恒等的に 1 なので**除外される**。
- **有限円柱の壁補正の閉形式**(`docs/free_space_l5_poisson.md`):
  軸方向モード \(k\) ごとに
  \(\psi_R-\psi_\infty=-\psi_\infty(R)\,[I_1(kr)/r]/[I_1(kR)/R]\)(\(k>0\))、
  \(k=0\) では定数 \(-\psi_\infty(R)\) で、compact support 源に対し
  \(\psi_\infty(R)=(2R^2)^{-1}\int\rho^3f\,d\rho\)。既存 E-33 と整合し、
  ゼロモードの \(R^{-2}\) 代数尾部を独立に再導出している。

### 数値的に観測されたこと(すべて有限周期円柱上の浮動小数点観測)

- **早期 Hou 窓は臨界位相で陰性**
  (`docs/hou_early_critical_l3_result.md`、
  `outputs/imported_chatgpt_results/critical_l3_old_uploaded_snapshot/`):
  193×384 で \(t=0\to T_1\) に \(\max|u|\) は 327.8→679.1 と増加する一方、
  表現領域の \(L^3\) ノルムは 118.57→107.94(**8.96% 減少**)、臨界密度の
  RMS 幅は \(L_r\) 0.2085→0.2584、\(L_z\) 0.1399→0.1664 と**拡大**、
  実効 shell 数は 1.915→1.814 と微減、外側 radial 臨界質量は 3.4e-8。
  すなわち critical mass の増大・スケールの増殖・外側尾部蓄積の
  いずれも観測されない。\(Q=A^3L_r^2L_z\) の 16 倍増は幅の収縮ではなく
  振幅増加が支配しており、Type-II 濃縮の証拠ではない。
- **候補特異時刻の走査は全滅**(`docs/early_hou_scaling_gate_result.md`):
  \(T\in[T_{\rm last}+10^{-6},T_{\rm last}+2\times10^{-3}]\) と 3 点以上の
  全 suffix 窓を走査して、3 解像度いずれも gate 通過 fit **0 件**。
  最良 fit でも \(\beta_r,\beta_z<0\)(幅が拡大)、\(\alpha\le0.139<1/2\)
  で Serrin 正則域内。
- **自由空間 radial Green solver の manufactured gate**
  (`docs/free_space_l5_poisson.md`): 5 次元 Gauss 試験の内部相対 \(L^2\)
  誤差は z padding 1/2/4 で 1.189e-2 → 4.643e-3 → 4.195e-3。最初の倍化で
  2.56 倍改善し、その後は radial 求積と source 打切りが支配。
  radial 壁は解析的に除去されるが、**\(z\) は zero padding 後も周期のまま**
  なので、これは period-image 診断であって自由空間 \(z\) の厳密証明ではない。
- **低波数壁 gate**(`docs/low_frequency_wall_obstruction.md`):
  軸周期を伸ばして \(kR\) を下げると、\(r\le2\) の有限壁相対誤差は
  9.06e-4(\(kR=25.1\))→ 2.25e-4(6.28)→ 6.35e-3(3.14)→ 3.29e-2(1.57)、
  ゼロモードで 9.13e-2。**短周期での壁感度の小ささは Fourier スペクトル
  ギャップで説明される**という、当方の P0-E の結論(periodic-z radial-wall
  sensitivity observation)を独立経路で裏づける。
- 統合時に追加した唯一のテスト: `free_space_poisson` の cephes Bessel 多項式
  は scipy 依存のため常に skip されていたので、当リポジトリの scipy-free
  オラクル `bessel_reference` との照合を追加した(相対 3e-7 で一致、
  ビット一致しないことも assert)。取り込んだ数学には手を触れていない。

### 条件付きのスケーリング分類(解ではなく探索領域)

`docs/type_ii_scaling_constraints.md`。有限エネルギー、有限散逸、臨界ノルム
増大、単一の最薄微分スケール、先頭項の非相殺を**仮定**した局所冪則 core に
対し、生き残る 3 族(\(B=2\beta_r+\beta_z\)、\(\gamma=\max(\beta_r,\beta_z)\)):

1. Euler 型 Type-II: \(2/5\le\gamma<1/2\)、\(\alpha=1-\gamma\)、
   \(2-2\gamma\le B\le3\gamma\)
2. 異方的放物型: \(\alpha=\gamma=1/2\)、\(1<B<3/2\)
3. 準定常粘性–慣性型: \(1/2<\gamma<1\)、\(\alpha=\gamma\)、
   \(4\gamma-1<B<3\gamma\)

走査結果は 20181 格子点中 438 点が条件付き実行可能
(`outputs/imported_chatgpt_results/type_ii_scaling/`)。
**これらは条件付き漸近スケーリング点であって Navier–Stokes の解ではない。**
除外されていない escape route(先頭相殺、成分別スケール、多重 core、
shell cascade、外場が担う \(L^3\) 増大、対数補正、非冪則)も明記されている。

### 未解決の証明義務・限界

- 上記の除外定理は ESS 端点定理を**引用**している。Lean 最終経路へ載せる
  なら忠実な形式化か明示的に監査された定理インタフェースが必要で、
  未証明の project 固有 axiom として挿入してはならない。
- ノートが提案する Lean 識別子 `F-4`〜`F-7` は、当リポジトリで既に
  `F-4`(証明書の有限次元不等式)を使用済みのため**採番が衝突**する。
  形式化着手時に `docs/formalization_map.md` 側で振り直す(未実施)。
- 参照出力は**古いリポジトリ snapshot** に対する外部計算であり、
  `257x512` を含まない。現行 `fable5-mainline` での再計算は未実施
  (`outputs/imported_chatgpt_results/README.md` に provenance 限界を明記)。
- 自由空間 solver はまだ非線形時間発展に接続されていない。Gate 4
  (非周期 \(z\)、\(R_{\max}\)/\(Z_{\max}\) 独立拡大、四方向収束表)は
  依然として未実装である。
- 臨界 \(L^3\) 診断は full-step streaming ではなく checkpoint 後処理のまま。

### Clay 問題について

**これらの統合物は Navier–Stokes ミレニアム問題を解決していない。**
臨界 \(L^3\) 障害は特異点の**構成ではなく候補クラスの除外**であり、
早期 Hou 窓の結果は**陰性の診断**である。後期 Hou 発展の正則性を
証明したわけでも、特異点を発見したわけでもない。過去の結果は
書き換えていない。

## 2026-07-29(第 2 便): PDE 項別釣合い診断の統合

外部セッションの追加バンドル(`term-balance-progress`)を統合した。
**統合作業のみ。新規計算・設計変更・数学的内容の改変はない。**
patch は現行ツリーへ競合なく適用でき、ZIP 内 7 ファイルすべてが
リポジトリへ反映されたことをハッシュで確認した。統合後の全テスト
**755 passed, 1 skipped**(統合前 747 passed)。

統合物: [term_balance.py](src/ns_certificate_lab/term_balance.py)、
[test_term_balance.py](tests/test_term_balance.py)(8 件)、
[analyze_term_balance.py](experiments/analyze_term_balance.py)、
ノート [term_balance_progress.md](docs/research_notes/term_balance_progress.md)、
参照出力 `outputs/imported_chatgpt_results/term_balance_old_snapshot/`。

### 入力の同一性(統合時に照合した新事実)

参照出力の CSV は入力 checkpoint ごとに `checkpoint_sha256` を記録している。
統合時に、**臨界 \(L^3\) 表と項別釣合い表の全 15 行(65×128 / 129×256 /
193×384 の 5 時刻)を現行リポジトリの `outputs/hou_early_time_v1` /
`v2_hires` の実バイトと照合し、両表とも 15/15 一致**した(mismatch 0)。
外部環境で走った診断だが、入力は現在コミットされている checkpoint と
ビット単位で同一である。したがって 2026-07-29 第 1 便の README にあった
「入力の同一性は未確認」という記述は**解消**した(該当 README を更新済み)。
ただし `257x512` はどの表にも含まれず、解析コードのバージョン同一性も
保証されないため、**現行コードでの再計算は依然として未実施**である。

### 数値的に観測されたこと(すべて有限周期円柱上の浮動小数点観測)

- **早期区間の中盤は非粘性項が支配的**: \(t=0.001\)、193×384 の臨界コア内で
  粘性項 / (移流+生成) は \(u_1\) 式 **2.85e-4**、\(\omega_1\) 式 **3.85e-4**。
  正規化残差は 0.0622 / 0.0668、移流–生成の相殺率は 0.352 / 0.243 で、
  分類は両式とも `time_inviscid`(時間項 ≈ −(移流−生成))。3 解像度で
  ほぼ一致する。これは既知の「\(\nu=5\times10^{-4}\) では粘性項が項和の
  \(O(10^{-4})\) しかなく符号反転すら実データで識別できない」という
  本リポジトリの P1-C 所見と独立に整合する。
- **ただし前登録ゲートは全件不合格**: 正規化残差 ≤ 0.10 かつ時間微分感度
  ≤ 0.20 かつ未解像分類でないこと、を要求すると各解像度で **0/8 合格**。
  主因は checkpoint 間隔の粗さで、時間微分感度は \(t=0.0005\) で既に
  0.31–0.58、\(t=0.0015\) 付近で 0.69–0.96 に達する。最終保存時刻の
  193×384 コア残差は \(u_1\) 0.318、\(\omega_1\) 0.444。
  **これは PDE が破れているという判定ではなく、保存時刻が疎すぎて
  時間微分の補間誤差と真の釣合い変化を区別できないという判定である。**

### 条件付きの位置づけ

早期 Hou 窓は「非粘性的な集中機構を示すが、Type-II・異方的・準定常粘性型の
いずれにも昇格できない」。第 1 便の結果(\(L^3\) が 9% 減少、臨界幅は拡大、
Type-II 指数 fit 全件不合格、shell 数は増えない)と合わせて、
**この早期区間は全空間特異点候補ではない**という陰性の診断が積み上がった。

### 未解決の証明義務・次に必要なこと

項別釣合いを数学的理由で判定するには、最終時刻の延長ではなく
**時間微分の解像**が要る。ノートが挙げる 4 経路の一致確認
(solver 内部 RHS / 隣接状態からの独立微分 / 別積分器 / 時間刻み半減系列)は
本リポジトリの既存資産で大部分が満たせる: `step_stream` の全 step 計装、
`take_step` の 3 積分器、Gate 1 の dt 系列。ただし
**時間微分を保存する経路は未実装**であり、これは
[STATUS 末尾の「次に行うべき最小の一手」項目 2](STATUS.md) と同じ作業に属する。

### Clay 問題について

**この統合物は Navier–Stokes ミレニアム問題を解決していない。**
項別釣合いは離散モデル上の**項の大きさの観測**であって、Euler 型特異点の
証明でも Navier–Stokes 特異点の証明でもない。過去の結果は書き換えていない。

## 2026-07-29(第 3 便): Track F 有限モード ansatz の除外定理と Lean F-6

`START_NEW_SESSION_NAVIER_STOKES.md` §6 の選択基準(全空間非周期 solver が
大規模未完成、外力あり反例の十分条件が未形式化、低次 ansatz の記号探索が
短時間で実装可能)はいずれも現 HEAD で成立していたので **優先候補A
(Track F)** を選んだ。ただし step 1–5 の**探索は実行していない**。
探索空間が空であることを step 6 として直接証明したためである。

### 数学的に証明したこと

新規ノート
[track_f_finite_mode_nogo.md](docs/research_notes/track_f_finite_mode_nogo.md)。

- **Lemma 1(T-1)**: 任意の実数値・発散ゼロな三角多項式 `u` に対し
  `⟨u,(u·∇)u⟩_{L²(𝕋³)} = 0`。Fourier 展開版の証明は
  `Σ_{k+l+m=0}(a_l·m)(a_k·a_m)` で `k↔m` を入れ替えて足すと
  `Σ(a_k·a_m)(a_l·(k+m)) = −Σ(a_k·a_m)(a_l·l) = 0`(発散ゼロ条件)。
- **Theorem 1**: 有限対称モード集合 `S ⊂ ℤ³` に対し、`u:[0,T)→V_S` が `C¹`、
  残差 `f = ∂_t u+(u·∇)u−νΔu+∇p` が `L¹((0,T);L²)` なら
  (i) `‖u(t)‖_{L²} ≤ ‖u(0)‖_{L²}+∫₀ᵗ‖f‖`、
  (ii) `sup_{t<T}‖u(t)‖_{H^s} ≤ (1+4π²R_S²)^{s/2}(‖u(0)‖+M)` および
  `sup‖∂^α u‖_{L^∞} ≤ (2πR_S)^{|α|}√|S|(‖u(0)‖+M)`、
  (iii) `Πf` が `T` を含む閉区間で連続なら `u` は `[0,T+δ)` へ滑らかに延長。
  **`ν` も `S` も (i) には現れない**(散逸は助けるだけ、非線形項は Lemma 1 で無害)。
- **Corollary 1**: 有限モード ansatz は Clay (D) 反例の破綻解になれない。
- **Corollary 2(帯域幅)**: 台が固定有限集合に留まる限り破綻しない。対偶として
  **Track-F 反例の速度場は `t→T⁻` で Fourier 帯域幅が非有界でなければならない**。
- **Corollary 3**: `ℝ³`(Clay (C))でも、固定有限次元の急減衰発散ゼロ空間に
  留まる ansatz は同じ理由で除外される。
- **必要条件 (N-1)(N-2)**(有限次元性を使わない): 任意の Track-F 反例は
  `sup_{t<T}‖u‖_{L²} < ∞` かつ `ν∫₀ᵀ‖∇u‖²dt < ∞`。すなわち
  **Clay 条件 (7) は外力を付けても自動的に満たされ、破綻は劣臨界ノルムでは
  起こり得ない**。「外力を使えば楽になる」という当初の想定は否定された。
- 新規性の主張はしない。Theorem 1 (i) は Galerkin 近似の大域可解性
  (Leray 1934 / Hopf 1951)の言い換えである。寄与は
  (a) それが Track F の探索空間をちょうど空にするという接続、
  (b) Lemma 1 の厳密算術による機械検証、(c) Theorem 1 (i) の Lean 化。

### Lean 4 で証明したこと(F-6)

`formal/NSSingularity/GalerkinNoBlowup.lean`(新規)。実 inner product 空間上で
`EnergyNeutral B := ∀x, ⟪x,B x x⟫ = 0`、`Dissipative A := ∀x, ⟪x,A x⟫ ≤ 0` と
定義し、5 定理を証明した。

| 定理 | 主張 |
|---|---|
| `norm_le_of_energy_inequality` | `⟪u,u'⟫ ≤ ‖u‖F` ⟹ `‖u b‖ ≤ ‖u 0‖+∫₀ᵇF` |
| `inner_galerkin_le` | `⟪x, g+B x x+A x⟫ ≤ ‖x‖‖g‖`(PDE 構造が入る唯一の箇所) |
| `galerkin_norm_le` | F-6 本体 |
| `galerkin_norm_le_of_mem` | `[0,T]` 上の一様上界 |
| `galerkin_not_tendsto_atTop` | `t→T⁻` で `‖u t‖ → +∞` は起こらない |

有限次元性は仮定していない。`F ≥ 0` も仮定していない(`‖g t‖ ≤ F t` から従う)。
証明は `t ↦ √(⟪u,u⟫+ε) − ∫₀ᵗF` の単調減少性による(`ε>0` の正則化は必須:
`‖·‖` は原点で微分不能で軌道は零点を通ってよい)。

`lake build` は **8660 jobs で成功**。`sorry`・`admit`・新規 `axiom` はゼロ。
`lake env lean AxiomAudit.lean` は全 14 定理について
`[propext, Classical.choice, Quot.sound]` のみを報告した。

**Lean 化していないこと**: (a) NS の移流項が実際に `EnergyNeutral` である
こと(= Lemma 1)、(b) 常微分方程式の延長(= Theorem 1(iii)、新規義務 `F-7`)、
(c) 有限次元空間上のノルム同値、(d) `ClayStatement.lean` との接続。

### 厳密算術で検証したこと(浮動小数点なし)

`src/ns_certificate_lab/galerkin_obstruction.py`(新規、`verify_trilinear_cancellation`)。
各 `k ∈ S` の `k^⊥` を **整数**基底(`k × e_i` の独立な 2 本)で座標付けし、
`Σ_{k+l+m=0}(a_l·m)(a_k·a_m)` を `ℤ[i]` 係数の 3 次多項式へ完全展開して
**全単項式係数がゼロ**であることを確認する。実験 `outputs/track_f_finite_mode_scan_v1`
での実測(10 族、全件合格):

| 族 | モード数 | 実次元 | 共鳴三つ組 | 展開単項式 | 残存単項式 | 判定 |
|---|---:|---:|---:|---:|---:|---|
| single_mode | 2 | 4 | 0 | 0 | 0 | rejected |
| planar_triad | 6 | 12 | 12 | 160 | 0 | rejected |
| oblique_triad | 8 | 16 | 12 | 576 | 0 | rejected |
| taylor_green | 8 | 16 | 0 | 0 | 0 | rejected |
| with_zero_mode | 7 | 15 | 31 | 192 | 0 | rejected |
| anisotropic_ladder | 10 | 20 | 24 | 608 | 0 | rejected |
| ball_one (\|k\|²≤1) | 6 | 12 | 0 | 0 | 0 | rejected |
| ball_two (\|k\|²≤2) | 18 | 36 | 120 | 3264 | 0 | rejected |
| ball_three (\|k\|²≤3) | 26 | 52 | 264 | 7872 | 0 | rejected |
| ball_four (\|k\|²≤4) | 32 | 64 | 426 | 11136 | 0 | rejected |

**故障注入**: 発散ゼロ拘束を外して `k` 自身を横断基底へ加えると、共鳴を持つ
全族で残存単項式が生じる(planar_triad 32、ball_two 352、ball_three 992、
ball_four 1506)。検出は `k·t = 0` 監査とは独立に、単項式検査だけで成立する。

### 数値的に観測したこと(浮動小数点クロスチェック。証明ではない)

`stream_apriori_bound`(RK4、全 step 監視)。`ball_two`(36 次元):

- **非粘性・無外力・振幅 200**(爆発を狙う設定): 2000 step 全てで上界を遵守。
  終端ノルム 199.99999995(相対 2.5e-10 の減少)、
  相対エネルギー生成 max **3.13e-16**、`bound_respected = true`。
- **粘性 1e-3・定方向外力 |f|=4**: 3000 step 全てで上界を遵守し、
  **`max_bound_ratio = 0.99995`** — 証明した上界は空虚ではなくほぼ sharp。
- **故障注入(縦方向成分を許す)**: 同じ設定で 200 step 目に発散
  (max ノルム 1.17e36)、相対エネルギー生成 **1.634**。
  クリーン実行との検出比 **5.2e15**。
- Parseval(`‖u‖_{L²} = ‖c‖₂`)、実数性、`k·a_k = 0` はいずれも 1e-12 以下で成立。

### 変更・追加したファイル

- 新規: `src/ns_certificate_lab/galerkin_obstruction.py`、
  `tests/test_galerkin_obstruction.py`(52 件)、
  `experiments/run_track_f_finite_mode_scan.py`、
  `configs/track_f_finite_mode_scan.json`、
  `outputs/track_f_finite_mode_scan_v1/`(summary/CSV/config snapshot/manifest+SHA-256)、
  `docs/research_notes/track_f_finite_mode_nogo.md`、
  `docs/final_target.md`、`formal/NSSingularity/GalerkinNoBlowup.lean`。
- 更新: `formal/NSSingularity.lean`、`formal/AxiomAudit.lean`、`formal/README.md`、
  `docs/formalization_map.md`(F-6/F-7 追加、公理監査更新)、
  `docs/known_obstructions.md`(§8.5 と O-FD/O-FE フィルター)、
  `docs/research_notes/README.md`(**Lean 識別子の採番衝突を解消**)、
  `PLAN.md`(Phase 2.9)、`README.md`(再現コマンド 10)、本書。

### Lean 識別子の採番衝突を解消

外部ノート `critical_l3_obstruction.md` §9 が提案していた `F-4`〜`F-7` と、
本リポジトリ既存の `F-4`(証明書不等式)/`F-5`(Clay 命題定義)の衝突を、
ノート本文を改変せずに `docs/final_target.md` §4 の登録簿で確定した:
ノートの `F-4`→`F-8`、`F-5`→`F-9`、`F-6`→`F-10`、`F-7`→`F-11`。
本セッションの 2 件は `F-6`(Galerkin 上界、形式化済み)と
`F-7`(ODE 延長、未着手)。以後 `docs/final_target.md` §4 が唯一の権威である。

### テスト結果

```text
PYTHONPATH=src .venv/Scripts/python.exe -m pytest -q
  -> 755 passed, 1 skipped(セッション開始時、HEAD fa6c2a5 で再確認)
  -> 807 passed, 1 skipped(本便追加後)
```

skip は `scipy.special` 依存の照合テスト 1 件(scipy は本リポジトリの依存ではない)。

### 未解決の証明義務・限界

- Corollary 2 の対偶は「帯域幅が非有界でなければならない」と言うだけで、
  そのような ansatz の**存在は何も示さない**。残る Track F の探索空間は
  Track U の動的再スケーリングと同じ困難を持つ。
- (N-3)(Ladyzhenskaya–Prodi–Serrin による臨界ノルム発散の必要性)は
  **引用**であり、本セッションで証明も形式化もしていない。
- Theorem 1(iii) の常微分方程式延長は紙上のみ(`F-7`)。
- Corollary 3(全空間版)は紙上のみ。部分積分の境界項評価は自明だが未形式化。
- 厳密算術の証明書は列挙した有限個のモード集合を覆うだけである。
  すべての `S` に対する主張は Lemma 1 の 2 行の証明であって、この scan ではない。

### Clay 問題について

**本便は Navier–Stokes ミレニアム問題を解決していない。** 有限モード除外は
特異点の**構成ではなく候補クラスの除外**であり、しかも証明した内容は
Galerkin 近似の大域可解性という古典的事実の言い換えである。
Track F が容易な近道ではないことが確定した、というのが実質的な成果である。
過去の結果は書き換えていない。

## 2026-07-29(第 4 便): Track U Gate 4 合格、Lean F-7a/F-7b/F-12/F-13、帯域幅発散候補の指数領域

主作業は **Track U の Gate 4**(非周期 `z` と自由空間楕円回復を備えた線形
全空間ゲート)、副作業は**固定有限モード no-go の Lean 完全化**、
第 3 の作業は **Track F の次の探索空間の必要条件導出**。

### 数学的に証明したこと

**(1) Gate 4 の厳密参照解(新規)。** `L₅ = ∂_rr + (3/r)∂_r + ∂_zz` は
5 次元 Laplacian の軸対称版なので、5 次元半径 `R = √(r²+z²)` のみに依存する
場に対し `L₅ψ = ψ'' + 4ψ'/R`。したがってコンパクト台の 5 次元動径 bump
`ω = c(1−(R/a)²)^p` に対し Newton の公式

```text
ψ_∞(R) = ∫_R^∞ s^{-4} m(s) ds,   m(s) = ∫_0^s t⁴ω(t)dt,   ψ_∞ = M/(3R³) (R ≥ a)
```

がすべて初等的な有限和になる。`L₅` の `z` 平行移動不変性から、異なる `z₀` に
置いた bump の**重ね合わせも厳密**である。求積も離散化も含まないので、
これは真に独立な参照解である。実測でも、参照解の PDE 残差は格子細分に対し
2 次で 0 へ収束する(9.04e-3 → 2.31e-3 → 5.80e-4、次数 1.968/1.993)。

**(2) a posteriori tail bound(新規、U-X5)。** `G(X) = 1/(8π²|X|³)` を `ℝ⁵` の
`−Δ` の基本解とし、源が `|Y| ≤ a` に台を持つとき

```text
|ψ_∞(X) − M₅G(X)| ≤ 3a‖ω‖_{L¹(dV₅)} / (8π²(|X|−a)⁴).
```

`|∇G(Z)| = 3/(8π²|Z|⁴)` と `|X−Y| ≥ |X|−a` による。5 次元 Laplacian の
**最大値原理**により、この値は単極子境界条件を課した box の
**連続レベルの領域打切り誤差**を内部の全点で抑える。
**離散化誤差は抑えない**(半径行 `i=1` は M-matrix でなく離散最大値原理がない)。

**(3) 周期 `z` のゼロ軸モードの厳密な過大評価率(新規、U-X4)。**
周期 `L` の `k=0` 軸モードは 4 次元動径 Poisson 方程式に従い、その全空間核は
`G₀(r,ρ) = 1/(2max(r,ρ)²)`。源の外で

```text
ψ_{k=0}(r) = M₅/(4π²L r²)     （代数的 r⁻² 尾部）
```

であり、真の自由空間場 `M₅/(8π²R³)` との比は **厳密に `2R/L`**。すなわち
周期箱は遠方場を、半径に**線形に比例して**過大評価する。実測は
`R=4,8,16,32`(`L=12`)で 0.667 / 1.333 / 2.667 / 5.333 と `2R/L` に
相対 1e-9 で一致した。**これが「周期箱を大きくして全空間に近づける」設計が
成立しない理由である。** 一方、非周期 `z` の Dirichlet 作用素の固有値は
`λ_m = −4sin²(mπ/2(N+1))/Δz²` で常に狭義負、**ゼロモードは存在しない**
(実測: `Z=4,8,16` で `|λ|_min` が `(π/2Z)²` と相対 1e-3 以内で一致、
ゼロモード検出ゼロ)。

**(4) 帯域幅発散候補の必要条件と実現可能領域(新規)。** 新規ノート
[track_f_shell_constraints.md](docs/research_notes/track_f_shell_constraints.md)。
二進シェル `E_j = A(t)(λ_j/N(t))^β`、`N = (T−t)^{−γ}`、`A = (T−t)^σ` に対し

| 必要条件 | 不等式 |
|---|---|
| エネルギー有界 (N-1) | `β > 0` かつ `σ ≥ 0` |
| 総散逸有限 (N-2) | `σ − 2γ > −1` |
| 臨界ノルム発散(ESS を**引用**) | `σ < γ` |
| 帯域幅発散(固定有限モード no-go の対偶) | `γ > 0` |

> **Proposition 3(仮定付き): 現在の shell ansatz と非退化仮定の下で `γ ≥ 1` は排除される。** `σ > 2γ−1` と `σ < γ` が
> 両立するには `2γ−1 < γ`、すなわち `γ < 1`。冪則・単一コア・対数補正なし・
> `β` 時間非依存という非退化仮定の下でのみ成り立つ。

実現可能領域は `0 < γ < 1`、`max(0,2γ−1) ≤ σ < γ`、`β > 0` の三角形のみ
(下端が到達可能なのは `γ < 1/2` のときだけ)。古典的放物型自己相似
`γ = 1/2, σ = 0` は散逸積分が対数発散するため**境界上で排除**される。

**(5) 滑らかな外力は高波数で無力(新規、F-N4)。** Clay の外力は
`𝕋³×[0,∞)` 上 `C^∞` なので、コンパクト時間区間上で
`‖f_j(t)‖_{L²} ≤ C_m λ_j^{−m}`(∀m)。シェル収支に入れると外力の寄与は
任意の多項式より速く消える。正確には: **滑らかな外力は高周波を直接供給できないが、低周波を制御して非線形カスケードを間接駆動する可能性は残る。**
排除されるのは直接注入経路 `⟨u_j,f_j⟩` のみで、低周波整形 → 非線形 triad →
高シェルという間接経路は排除されていない。
なお `Π_j` の粗い上界(`|Π_J| ≲ N^{5/2}A·E^{1/2}` 対 必要量 `≳ νN²A`)は
`(γ,σ)` に追加制約を与えないので、**追加制約を捏造していない**。

### Lean 4 で証明したこと(新規 10 定理、`FiniteModeNoGo.lean`)

| ID | 定理 | 主張 |
|---|---|---|
| **F-12** | `advectionForm_eq_zero` | `k_i·a_i = 0` の下で共鳴 3 次形式 `Σ_{p+q+s=0}(a_q·k_s)(a_p·a_s) = 0`。`p↔s` の対合を `Finset.sum_nbij'` で使う |
| **F-13** | `weighted_sq_sum_le` | `Σ w_i c_i² ≤ W Σ c_i²` — `H^s` 対 `L²` の定数 |
| **F-13** | `sq_sum_abs_le_card_mul_sum_sq`, `sum_abs_le_sqrt_card_mul_sqrt_sum_sq` | `(Σ|c|)² ≤ |S|Σc²` — `L^∞` 対 `L²` の定数 |
| **F-7a** | `exists_tendsto_nhdsWithin_of_norm_deriv_le` | `[0,T)` 上で導関数が連続かつ有界なら `t→T⁻` で極限が存在 |
| **F-7a** | `intervalIntegrable_of_continuousOn_bounded` | 有界連続導関数の区間可積分性 |
| **F-7b** | `contDiff_galerkinField`, `exists_local_galerkin_solution` | 自励 Galerkin 場は `C¹` で、任意点を通る両側局所解を持つ |
| — | `not_isBreakdownCandidate_of_galerkin` | **F-6 の状態上界 → 速度上界 → F-7a** を連結し、固定有限モード候補が破綻候補になれないことを結論 |
| — | `galerkin_bounded_and_reaches_endpoint` | 有界性と端点到達を梱包した形 |

`lake build` は **8661 jobs で成功**。`sorry`・`admit`・新規 `axiom` はゼロ。
`lake env lean AxiomAudit.lean` は**全 24 定理**について
`[propext, Classical.choice, Quot.sound]` のみを報告。

**F-7c(時間依存外力での局所延長)は未着手**。当初の見立て「F-7 は既存定理の
単純な適用」は**誤りだった**。実際に mathlib API を確認した結果、
`ContDiffAt.exists_forall_mem_closedBall_exists_eq_forall_mem_Ioo_hasDerivAt₀`
は **`f : E → E` の自励系専用**であり、`g(t)` が時間依存だと直接適用できない。
必要な継続補題は (i) `E × ℝ` 上への自励化(第 2 成分が `s(t)=t+c` である
ことの証明を含む)、または (ii) 場 `(t,x) ↦ g t + B x x + A x` に対する
`IsPicardLindelof` の直接構成(uncurry の連続性、球上で `x` について一様な
Lipschitz 定数、ノルム上界の 3 点)。詳細は
[final_target.md](docs/final_target.md) §4.1。**公理化していない。**

### 数値的に観測したこと(すべて binary64。証明ではない)

`outputs/whole_space_gate4_v1`、**前登録受入検査 20 件すべて合格**。

**Gate 4 収束表(格子細分、`R_max=Z_max=8`、内部 `R≤3`)**

| 境界条件 | `n_r` | 相対 `L²` 誤差 | 観測次数 | 離散残差 |
|---|---:|---:|---:|---:|
| monopole | 65 | 1.7123e-02 | — | 1.81e-14 |
| monopole | 129 | 4.2126e-03 | **2.023** | 9.41e-14 |
| monopole | 257 | 1.0488e-03 | **2.006** | 5.70e-13 |
| zero | 65 | 2.4049e-02 | — | |
| zero | 129 | 1.2295e-02 | 0.968 | |
| zero | 257 | 9.8522e-03 | **0.320** | |

零トレースは**収束が飽和する**(打切り誤差が即座に支配)。これは欠陥ではなく
所見であり、単極子トレースを使う定量的理由そのものである。

**領域拡大表(同一刻み `Δ=0.125`、共通内部差分で打切り成分を分離)**

| 境界条件 | 拡大方向 | 観測 rate |
|---|---|---|
| zero | radial | 3.11, 3.13(理論 `D⁻³`) |
| zero | axial | 4.22, 4.36 |
| monopole | radial | 3.85, 3.74(理論 `D⁻⁴`) |
| monopole | axial | 5.45, 4.95 |

**内部誤差を参照解と直接比べると離散化誤差の床で飽和する**ため、領域拡大の
効果は同一刻みでの共通内部差分でしか見えない。この事実自体を記録する。

**tail 誤差表(単極子トレース)**

| 拡大 | `R_max` | `Z_max` | 境界データ誤差 | 内部打切り差 | a posteriori 上界 | 鋭さ |
|---|---:|---:|---:|---:|---:|---:|
| radial | 5.00 | 20.0 | 1.812e-06 | 1.614e-06 | 1.635e-04 | 0.011 |
| radial | 7.50 | 20.0 | 3.533e-07 | 3.390e-07 | 1.713e-05 | 0.021 |
| radial | 11.25 | 20.0 | 6.949e-08 | 7.454e-08 | 2.331e-06 | 0.030 |
| axial | 20.0 | 5.00 | 8.529e-06 | 1.625e-06 | 1.635e-04 | 0.052 |
| axial | 20.0 | 7.50 | 1.518e-06 | 1.781e-07 | 1.713e-05 | 0.089 |
| axial | 20.0 | 11.25 | 2.791e-07 | 2.390e-08 | 2.331e-06 | 0.120 |

上界は全行で境界データ誤差と内部打切り差の**両方を支配**する。
鋭さ 0.01〜0.12 なので上界は 8〜90 倍保守的である(勾配評価が粗いため)。

**周期像誤差と radial 壁誤差の分離(連続レベル、`R_max=8`、内部 `R≤2.5`)**

| 半周期 | 周期像誤差(厳密) | 壁誤差上界 | 非周期打切り上界 | 実測 有限円柱 | 実測 非周期単極子 |
|---:|---:|---:|---:|---:|---:|
| 6.0 | 7.115e-06 | 1.073e-05 | 5.670e-05 | 2.406e-03 | 2.913e-03 |
| 12.0 | 7.119e-07 | 7.969e-06 | 1.231e-05 | 2.405e-03 | 2.913e-03 |

**失敗した gate として明記する**: 4 経路の**実測**誤差(~2.4e-3)は
本リポジトリが払える刻み幅では**すべて離散化誤差に支配され**、
周期像の寄与(~7e-6)を実測では分離できなかった。したがって分離は
厳密参照解から連続レベルで行った。非周期 `z` では周期像成分は
**近似的にゼロではなく厳密にゼロ**である(像が存在しない)。

**独立 Cartesian 検査**(円柱演算子を一切呼ばない経路)

| 点数 | `div u` 最大 | 次数 | `curl` 誤差 | 次数 | 相対 |
|---:|---:|---:|---:|---:|---:|
| 33 | 6.332e-03 | — | 2.089e-02 | — | 1.051e-01 |
| 65 | 1.775e-03 | 1.84 | 6.330e-03 | 1.72 | 3.141e-02 |
| 129 | 4.503e-04 | **1.98** | 1.604e-03 | **1.98** | 7.802e-03 |

検査した恒等式は `∇·u = 0`(E-15)と **`∇×u = ω₁·(−y, x, 0)`**。後者は
`ω^θ = rω₁` と `ê_θ = (−y,x,0)/r` から従い、**楕円解と物理 curl を
solver とコードを共有しない経路で結ぶ**。軸方向 curl 成分は速度の 1e-3 以下。

**故障注入(4 種、すべて検出)**

| 注入 | 軸誤差(65/129/257) | 観測次数 | 検出 |
|---|---|---|---|
| healthy | 7.39e-04 / 1.85e-04 / 4.64e-05 | 1.99, 2.00 | — |
| 軸係数 8→4 | 1.92e-03 / 4.80e-04 / 1.20e-04 | 2.00, 2.00 | **振幅比 2.6 倍** |
| 半径 drift 3→1 | 6.90e-02 / 7.00e-02 / 7.03e-02 | −0.02, −0.01 | **収束停止** |
| 外側結合の脱落 | — | — | 誤差 5 倍超 |
| 源ノルム 1000 倍過小申告 | — | — | tail bound が破れる |

**重要な限界を記録する**: 軸係数を 8 から 4 に壊しても**観測次数は 2 のまま**
である。欠陥行の制御体積は 5 次元重み付けで `O(dr⁴)` しかないため寄与は
依然収束する。**収束次数だけを見るゲートはこの故障を見逃す**ので、
振幅の比較を併用しなければならない。

### 変更・追加したファイル

新規: [whole_space_gate.py](src/ns_certificate_lab/whole_space_gate.py)、
[shell_constraints.py](src/ns_certificate_lab/shell_constraints.py)、
[test_whole_space_gate.py](tests/test_whole_space_gate.py)(34 件)、
[test_shell_constraints.py](tests/test_shell_constraints.py)(21 件)、
[run_whole_space_gate4.py](experiments/run_whole_space_gate4.py)、
[whole_space_gate4.json](configs/whole_space_gate4.json)、
`outputs/whole_space_gate4_v1/`(summary、4 CSV、config snapshot、manifest+SHA-256)、
[track_f_shell_constraints.md](docs/research_notes/track_f_shell_constraints.md)、
[FiniteModeNoGo.lean](formal/NSSingularity/FiniteModeNoGo.lean)。

更新: `formal/NSSingularity.lean`、`formal/AxiomAudit.lean`、`formal/README.md`、
`docs/final_target.md`(F-7a/b/c・F-12・F-13 追加、Gate 4 合格反映、§4.1 新設)、
`docs/formalization_map.md`、`docs/whole_space_transition.md`(Gate 4 合格記録)、
`docs/known_obstructions.md`、`docs/research_notes/track_f_finite_mode_nogo.md`、
`PLAN.md`(Phase 2.85 新設)、`README.md`(再現コマンド 11)、
`.github/workflows/tests.yml`、本書。

### 表現の修正(第 3 便の記述を訂正)

1. 「Track F は空」→ **「固定有限次元・固定帯域 Track F は除外」**。
   帯域幅発散族は除外されておらず、その必要条件を第 4 便で導出した。
2. 「劣臨界ノルムでは破綻しない」→ **「`L²` エネルギーと時間積分された `H¹`
   散逸は、滑らかな外力の下で有限時間内に自動的に有界であり、特異点の
   直接的な発散指標にはできない」**。
3. 「F-7 は既存定理の単純な適用」→ 実際に mathlib API を確認し、
   F-7a/F-7b は閉じたが **F-7c には自励化または `IsPicardLindelof` の
   直接構成が必要**であることを明記した。

### テスト結果

```text
PYTHONPATH=src .venv/Scripts/python.exe -m pytest -q
  -> 807 passed, 1 skipped(第 3 便時点)
  -> 862 passed, 1 skipped(本便追加後)
```

`lake build`: **Build completed successfully (8661 jobs)**。
`lake env lean AxiomAudit.lean`: 全 24 定理が古典公理のみ。
`git ls-files formal | xargs grep -InE '\bsorry\b|\badmit\b|^[[:space:]]*axiom '`:
文書コメント中の言及のみ。

### 未解決の証明義務・限界

- **Gate 4 は線形ゲートである。** 非線形時間発展への結合は許可されたが、
  得られる発展については何も主張しない。次の律速は「一様全空間格子で
  Hou 型 front を解像できるか」という解像度設計であり、
  この判断の前に非線形全空間実行を開始しない。
- a posteriori 上界は**連続レベル**の打切り誤差のみを抑える。離散化誤差は
  測定しただけであり、区間演算ではない。
- 製造源は `C^{p−1}`(`p=6` で `C⁵`)であって `C^∞` ではない。`C^∞` の bump は
  初等的な自由空間ポテンシャルを持たず、参照解が独立でなくなる。
- 周期像誤差と壁誤差の**実測**分離は失敗した(離散化誤差に支配)。
- F-12 は **Fourier 表現の代数的恒等式**であり、`∫_{𝕋³}u·(u·∇)u = 0` の
  形式化ではない。両者を繋ぐトーラス関数空間は mathlib にない。
- F-7c は未着手。`ClayStatement.lean` との接続もない。
- シェル制約の §4.3 は **ESS 端点定理の引用**に依存する。冪則 ansatz の外
  (対数補正、複数コア、時間依存 `β`)は覆っていない。`ℝ³` 版も未導出。

### Clay 問題について

**本便は Navier–Stokes ミレニアム問題を解決していない。** Gate 4 は
**線形楕円ソルバの検証**であり、特異点にも正則性にも触れていない。
Lean の追加分は固定有限モードという限定クラスの除外を機械検証したもので、
Clay 命題とは接続されていない。シェル制約は**必要条件**であって存在ではなく、
候補は依然として一つも構成されていない。過去の結果は書き換えていない。

## 2026-07-29(第 5 便): 微分 tail 上界、自由空間速度回復 API、小振幅全空間非線形 run

主作業は**線形 Gate 4 を「微分まで保証された全空間非線形入口」へ昇格**させること。
副作業は Lean(F-14/F-15/F-16/F-7c 還元)と、Track F の外力役割分離。

### 表現の修正(第 4 便の記述を訂正)

1. 「滑らかな外力は Track U に対する優位性を与えない」→
   **「滑らかな外力は高周波を直接供給できないが、低周波を制御して非線形
   カスケードを間接駆動する可能性は残る」**。減衰評価が閉じるのは 5 経路のうち
   直接注入 1 本だけである。§「cascade 模型」で模型検査した。
2. 「`γ≥1` は全面排除」→ **「現在の shell ansatz と非退化仮定の下で `γ≥1` は排除」**。
   非退化仮定(単一冪 `β`、単一コア、対数補正なし、`β` 時間非依存、
   `≍` の隠れ因子が発散しない)を明示した。Lean では `ShellAdmissible` 構造体
   としてフィールドに固定した(F-16)。
3. 「楕円側の障害は消えた」→ **「自由空間ポテンシャル**値**の線形 Gate を
   通過した。速度回復に必要な微分 tail と非線形結合は未検証」**。
   本便でその未検証部分を扱った。

### 数学的に証明したこと

**(1) Green 核の微分の斉次性定数(新規、U-X6 の基礎)。** `G₅ = c₅|X|^{−3}` は
`(−3)` 斉次なので `‖D^mG₅(Z)‖ = A_m|Z|^{−3−m}` が**恒等式**。閉形式の微分
テンソルから `A_m/c₅ = 1, 3, 12, 150, 1620, 21420`(`m=0..5`)。
`A_1=3`、`A_2=12` は**厳密**(Hessian `−3|Z|^{−5}I + 15|Z|^{−7}Z⊗Z` の固有値は
`12` と `−3`)。`m≥3` は三角不等式による上界で、`A_3` の sharp 値 60 に対し
2.5 倍保守的であることまで分かっている。**「ある定数」ではない。**

**(2) 微分 tail 上界(新規、U-X6)。** 源が `|Y| ≤ a` に台を持ち `d = |X|−a > 0` のとき

```text
|D^kψ_∞ − D^k(M₅G₅)| ≤ A_{k+1} I₁ / d^{4+k}                (k = 0,1,2,3)
|D^kψ_∞ − D^k(M₅G₅ − P·∇G₅)| ≤ ½ A_{k+2} I₂ / d^{5+k}
```

`I₁ = ∫|Y||ω|dV₅`、`I₂ = ∫|Y|²|ω|dV₅`。証明は積分記号下の `k` 回微分と、
線分 `X−θY` 上の平均値定理のみ。**最大値原理は一切使っていない。**

**(3) 内部への微分伝播(新規、U-X7)。** 切断誤差 `e` は `ℝ⁵` 上調和なので、
古典的**内部楕円評価** `|D^kh| ≤ (nk/ρ)^k sup|h|`(`n=5`)から
`|D^ke| ≤ (5k/ρ)^kε₀`。`k=0` の場合だけが最大値原理である。
速度への翻訳は `|δu| ≤ 2ε₀ + R_max(5/ρ)ε₀`。

**(4) 軸対称簡約核と求積の一致(新規)。** `∫_{S³}f(x̂·ŷ)dσ = |S²|∫_{−1}^1
f(t)√(1−t²)dt` により簡約核は 1 次元積分になり、重み `√(1−t²)` は
**第 2 種 Gauss–Chebyshev** が厳密に一致する。これが独立評価経路の基礎。

**(5) 外力の役割分離(新規、F-N4 の訂正)。** Clay 外力が使えない経路は
**高周波への直接注入 1 本だけ**である。低周波注入、非線形 triad 経由の
shell flux、位相・整列制御、粘性減衰を超える cascade 維持の 4 経路は
**排除されていない**。

### Lean 4 で証明したこと(新規 10 定理、`GreenAndCascade.lean`)

| ID | 定理 | 主張 |
|---|---|---|
| **F-14** | `greenProfile_radial_laplace_eq_zero` ほか 3 件 | `R^{−3}` は `f''+4f'/R` に消される(5 次元動径 Laplace) |
| **F-15** | `flux_newtonSlope`, `hasDerivAt_flux` | Newton の flux 恒等式 `R⁴ψ'(R) = −m(R)` とその微分 |
| **F-16** | `ShellAdmissible` 構造体 + 3 定理 | shell 指数の 4 条件を**構造体で明示**し、`γ<1`、`σ ∈ Ico(max 0 (2γ−1)) γ`、`γ≥1` の排除 |
| — | `breakdown_time_set_empty` | 固定有限帯域候補の**破綻時刻の集合は空**(Clay への限定的接続。1 点の否定から集合の空性へ強化) |
| — | `galerkin_solution_of_autonomised` | **F-7c 還元**: 自励化場の局所流があれば時間依存系にも局所解がある。第 2 成分が `s(t)=t` であることを導関数ゼロから示す |

`lake build` **8662 jobs 成功**。`sorry`・`admit`・新規 `axiom` ゼロ。
`lake env lean AxiomAudit.lean` は**全 34 定理**について
`[propext, Classical.choice, Quot.sound]` のみ。

**F-7c は閉じていない。** 残るのは自励化した場
`F(x,s) = (g s + B x x + A x, 1)` の局所流の構成であり、**公理化していない**。

### 数値的に観測したこと(すべて binary64。証明ではない)

`outputs/whole_space_gate5_v1`、**前登録受入検査 25 件すべて合格**。

**微分 tail 上界 対 実測誤差**(閉形式参照解、球面 `|X|=R` 上 121 点)

| 打切り | `R` | `d` | 値 | 勾配 | Hessian |
|---|---:|---:|---:|---:|---:|
| monopole | 2.5 | 0.804 | 0.0087 | 0.0032 | 0.0005 |
| monopole | 4.0 | 2.304 | 0.0712 | 0.0449 | 0.0113 |
| monopole | 6.0 | 4.304 | 0.1507 | 0.1151 | 0.0351 |
| dipole | 2.5 | 0.804 | 0.0020 | 0.0003 | 0.0001 |
| dipole | 4.0 | 2.304 | 0.0322 | 0.0077 | 0.0025 |
| dipole | 6.0 | 4.304 | 0.0921 | 0.0269 | 0.0109 |

(値は実測誤差/上界。全 18 セルで `<1`。)距離を 2 倍にすると上界は
値 `2⁻⁴`、勾配 `2⁻⁵`、Hessian `2⁻⁶` で減少(実測 4.000000/5.000000/6.000000)。
**源ノルムを 1000 倍過小申告すると上界が破れる**ことも確認した。

**独立 Green 経路**(核を解析微分してから求積、ソルバ不使用)

| 格子 | `ψ` | `∂_r` | `∂_z` | `∂_{rr}` | `∂_{zz}` | `∂_{rz}` |
|---|---:|---:|---:|---:|---:|---:|
| 161×321 | 2.51e-07 | 2.26e-07 | 2.78e-07 | 2.16e-07 | 3.21e-07 | 3.74e-07 |
| 321×641 | 1.57e-08 | 1.41e-08 | 1.73e-08 | 1.35e-08 | 2.00e-08 | 2.33e-08 |

観測次数は 6 成分すべてで **4.00**。

**速度回復 API の空間収束**(`R=Z=8`、内部 `R≤3`、最大誤差と観測次数)

| `n_r` | `ψ₁` | `∂_rψ₁` | `∂_zψ₁` | `u^r` | `u^z` |
|---:|---|---|---|---|---|
| 65 | 7.39e-04 | 5.71e-04 | 9.16e-04 | 2.36e-04 | 1.48e-03 |
| 129 | 1.85e-04 (1.99) | 1.40e-04 (2.03) | 2.26e-04 (2.02) | 6.16e-05 (1.94) | 3.71e-04 (1.99) |
| 257 | 4.64e-05 (2.00) | 3.49e-05 (2.01) | 5.61e-05 (2.01) | 1.54e-05 (2.00) | 9.28e-05 (2.00) |

軸正則性は**厳密**(`max|u^r|_{r=0}| = 0`、`max|ψ_{1,r}|_{r=0}| = 0`)。
領域拡大では勾配の共通内部差分が `7.51e-07 → 9.59e-08 → 1.44e-08`(radial)、
`1.11e-06 → 1.04e-07 → 1.11e-08`(axial)と減少し、いずれも微分 tail 上界
(`1.31e-03 → 1.45e-04 → 2.04e-05`)に支配される。
故障注入: 軸係数 `8→4` で勾配誤差 `1.40e-04 → 2.47e-03`(17.6 倍)、
半径 drift `3→1` で `4.51e-02`(321 倍)。両方検出。

**小振幅全空間非線形 run**(純粋旋回、振幅 0.05、`ν=5e-3`、`T=0.1`)

| 系列 | 設定 | 最終 `max|ω₁|` | rel div | curl defect |
|---|---|---|---|---|
| joint 細分 | 33×65 | 4.57017781e-04 | 2.47e-02 | 7.16e-02 |
| joint 細分 | 65×129 | 4.89858127e-04 | 7.08e-03 | 1.80e-02 |
| joint 細分 | 129×257 | 4.94903624e-04 | 1.83e-03 | 4.59e-03 |
| `dr` のみ | 33/65/129 × 129 | — | 2.47e-02→7.08e-03→1.83e-03 | 3.70e-02→1.80e-02→**1.40e-02(飽和)** |
| `dz` のみ | 65 × 65/129/257 | — | **7.06e-03→7.08e-03→7.08e-03(飽和)** | 5.31e-02→1.80e-02→9.42e-03 |

joint 細分の観測次数は **1.989 / 1.974**。`dr` のみ細分すると divergence が
収束し curl が飽和、`dz` のみ細分すると逆になる — **2 つの誤差源が分離できている**。
積分器 3 種(Heun/SSPRK3/RK4)の相対差は **2.40e-07**。
エネルギーは全 run で単調減少(増加量の最大 0)、Poisson 残差 `<1e-12`、
外側帯の質量比 `~1e-87`、閾値で落とした源質量比 `<1e-10`、
非線形 RHS の境界感度 `~1e-25`、低周波成分比 0.81。
独立 Cartesian 検査(円柱演算子を呼ばない経路)は
`div` 相対 3.16e-03、`curl` 相対 2.86e-02。

**cascade 模型**(`N=10`、外力はシェル 0–1 のみ、測定帯は `j≥4`)

| ケース | 帯域幅最大 | 到達最高シェル | 高シェル振幅最大 | 直接注入 |
|---|---:|---:|---:|---:|
| 無外力・粘性 | 0.0012 | 2 | 5.42e-28 | 0 |
| **低モード外力・粘性** | **1.2132** | **6** | **1.31e-02** | **0** |
| 低モード外力・粘性 10 倍 | 0.5752 | 4 | 4.15e-07 | 0 |
| 低モード外力・非粘性 | 4.8076 | 9 | 8.70e-02 | 0 |
| 無外力・非粘性 | 0.0016 | 2 | 3.69e-26 | 0 |

**直接注入が厳密にゼロのまま、低周波のみの外力が高シェル振幅を 26 桁動かす。**
粘性を 10 倍にすると `1.31e-02 → 4.15e-07` に落ちるので、勝負は
「外力の減衰」ではなく「triad flux 対 `νk²`」である。
非線形転送のエネルギー欠損は全ケース `<3e-18`。

### 失敗した gate(隠さずに記録)

**小振幅ではこの 3 つの部分ゲートが情報を持たなかった。**

1. **時間細分**(`dt`, `dt/2`, `dt/4`): 最終 `max|ω₁|` が
   `4.89858127e-04` で印字桁まで一致。RK4 の時間誤差がこの振幅では
   空間誤差より桁違いに小さく、**時間次数を測定できない**。
2. **領域拡大**(`R_max ∈ {3,4,6}`、`Z_max ∈ {3,4,6}`): 同じく完全一致。
   場が外側境界に到達しない(outer band fraction `~1e-49`)ため、
   **領域独立性の確認にはなるが外側境界条件の検証にはならない**。
3. **境界次数**(monopole 対 zero): 相対差 **厳密に 0**。同じ理由で
   **多重極の次数比較がこの振幅では空虚**。

いずれも「正しく無関係」であることの確認ではあるが、**外側境界条件そのものは
未検証のまま**である。場が境界に届く振幅で 3 つを再実行するまで、
強振幅 Hou 型全空間候補へ進んではならない。

さらに、微分 tail 上界の内部評価 `(5k/ρ)^k` は既定の `interior_radius` が
`boundary_radius − max(dr,dz)` なので、**格子を細かくすると評価点が境界に
近づいて上界が悪化する**(`2.21e-04 → 4.50e-04 → 9.01e-04`)。これは実装既定値の
性質であり、比較時は `interior_radius` を固定しなければならない。

### 変更・追加したファイル

新規: [free_space_recovery.py](src/ns_certificate_lab/free_space_recovery.py)、
[whole_space_evolution.py](src/ns_certificate_lab/whole_space_evolution.py)、
[cascade_toy.py](src/ns_certificate_lab/cascade_toy.py)、
[test_free_space_recovery.py](tests/test_free_space_recovery.py)(23 件)、
[test_whole_space_evolution.py](tests/test_whole_space_evolution.py)(22 件)、
[run_whole_space_gate5.py](experiments/run_whole_space_gate5.py)、
[whole_space_gate5.json](configs/whole_space_gate5.json)、
`outputs/whole_space_gate5_v1/`、
[green_derivative_tail_bounds.md](docs/research_notes/green_derivative_tail_bounds.md)、
[cascade_toy_model.md](docs/research_notes/cascade_toy_model.md)、
[GreenAndCascade.lean](formal/NSSingularity/GreenAndCascade.lean)。

更新: [whole_space_gate.py](src/ns_certificate_lab/whole_space_gate.py)(厳密 2 階微分)、
`formal/NSSingularity.lean`、`formal/AxiomAudit.lean`、`formal/README.md`、
`docs/final_target.md`、`docs/formalization_map.md`、
`docs/whole_space_transition.md`、`docs/research_notes/README.md`、
`docs/research_notes/track_f_shell_constraints.md`(表現修正)、
`docs/proof_obligations.md`、`PLAN.md`(Phase 2.86)、`README.md`(再現コマンド 12)、
`.github/workflows/tests.yml`、本書。

### テスト結果

```text
PYTHONPATH=src .venv/Scripts/python.exe -m pytest -q
  -> 862 passed, 1 skipped(第 4 便時点)
  -> 907 passed, 1 skipped(本便追加後)
```

`lake build`: **Build completed successfully (8662 jobs)**。
`lake env lean AxiomAudit.lean`: 全 34 定理が古典公理のみ。

### 未解決の証明義務・限界

- **微分 tail 上界は連続レベルの命題である。** 離散解の微分誤差は測定した
  だけで、抑えていない。半径行 `i=1` は M-matrix でないので離散内部評価もない。
- 製造源は `C^{p−1}`(`p=6` で `C⁵`)であって `C^∞` ではない。
- Green 求積路は**源の台の外でのみ**有効(内部で簡約核が特異)。
- 台半径は数値的閾値に依存し、落とした質量比を毎回報告している。
  **落とした質量が無視できない実行では上界は成立しない。**
- 小振幅 run は外側境界条件を検証していない(上記「失敗した gate」)。
- F-7c は未着手(還元のみ)。`ClayStatement.lean` との接続もない。
- cascade 模型は Navier–Stokes ではない。幾何も位相もなく、
  非粘性版が有限時間爆発することが知られているので**爆発の証拠としては無価値**。
  判定できるのは設計可否だけである。

### Clay 問題について

**本便は Navier–Stokes ミレニアム問題を解決していない。** Gate 5 は
**線形楕円回復の微分保証と、小振幅での非線形結合の健全性検査**であり、
特異点にも正則性にも触れていない。Lean の追加分は Green 核の動径恒等式と
指数算術と ODE 還元であり、Clay 命題とは接続されていない。cascade 模型は
外力の役割についての設計可否判定であって PDE の定理ではない。
候補は依然として一つも構成されていない。過去の結果は書き換えていない。

## 2026-07-29(第 6 便): 非線形 tail 伝播、明示的初期値族、Gate 6、区間証明書

主作業は Gate 5 の小振幅健全性試験を**中振幅 Gate 6** へ昇格させ、その後
明示的候補族の**振幅・形状継続**を実行すること。副作業は Lean 証明書層と
区間演算の開始。

### 数学的に証明したこと

**(1) 非線形領域打切り誤差の伝播(新規、U-X8)。**
`ε₀ ≥ |δψ₁|`、`ε₁ ≥ |∇δψ₁|`、`ε₂ ≥ |D²δψ₁|` から、回復が線型であることから

```text
|δu^r| ≤ R_max·ε₁,        |δu^z| ≤ 2ε₀ + R_max·ε₁
```

積差恒等式 `ab − ãb̃ = (a−ã)b + ã(b−b̃)` から移流項について

```text
|δ(u^r f_r + u^z f_z)| ≤ |δu^r|·‖f_r‖ + |δu^z|·‖f_z‖ + ‖ũ^r‖·‖δf_r‖ + ‖ũ^z‖·‖δf_z‖
```

(末尾 2 項は状態誤差に比例するので Lipschitz 定数へ集約)。旋回源
`2u₁ψ_{1,z}` は `2‖ũ₁‖ε₁ + 2‖ψ_{1,z}‖|δu₁|`、伸長源 `∂_z(u₁²)` は
`(u₁−ũ₁)(u₁+ũ₁)` の分解から得られ、**tail は直接には入らない**(`u₁` は
状態変数でありポテンシャル微分ではない)ので定数項はゼロ、寄与は
Lipschitz 定数のみ。最後に `L^∞` 最大値原理(発散ゼロ移流 + `L₅` 拡散)から

```text
d/dt‖e‖_∞ ≤ D + Λ‖e‖_∞  ⟹  ‖e(t)‖_∞ ≤ (‖e(0)‖ + D t) e^{Λ t}.
```

**これが要求された「証明書形式」である。** 全定数は
`src/ns_certificate_lab/tail_propagation.py` に明示実装した。
**限界を明記する**: 状態誤差 `e_ω` から速度誤差への変換は自由空間解作用素の
ノルムを要し、Biot–Savart は `L^∞→L^∞` で有界ではないので、その定数は
**入力**として持ち回り、厳密供給は未解決義務である。

**(2) 明示的全空間初期値族(新規)。** `src/ns_certificate_lab/initial_data.py`:

```text
u₁(0,r,z) = A·χ(r²/R₀²)·χ(z²/Z₀²)·(z/Z₀)/(1+c(z/Z₀)²),   ω₁(0)=0,
χ(s) = exp(−1/(1−s))  (0≤s<1),  0  (s≥1).
```

**radial 因子は `r²` の関数**である。物理旋回は `u^θ = ru₁`、Cartesian では
`u₁·(−y,x,0)` なので、`u₁` が `(r²,z)` の滑らかな関数であるときにのみ原点で
`C^∞` になる。`χ(r/R₀)` と書くと軸上で Lipschitz どまりになる。
紙上とテストの両方で確認: Cartesian 滑らかさ(軸を跨ぐ 4 階差分 `<1e-9`)、
**厳密な発散ゼロ**(純旋回なので恒等的に 0、Cartesian ステンシルでは
2 次収束を確認)、コンパクト台、有限エネルギー、有限 `L³`、軸正則性。
粘性は**固定正定数**で、Hou の二段階粘性は候補計算に使わない。

**(3) 対称性による多重極の退化(新規、U-X9)。** コンパクト台の初期値では
`∫ω₁dV₅ = 0`(実測 `1.5e-17` 相対)。したがって **zero トレースと monopole
トレースは同一**であり、軸方向四重極も消えるので **dipole と quadrupole も同一**
(退化/非退化の比 `5.2e-11`)。有効な比較は 1 段だけである。
**これが第 5 便で「monopole 対 zero が空虚だった」ことの真因**であり、
場が境界に届かなかったこと以上の構造的理由である。

### Lean 4 で証明したこと(新規 10 定理、`CertificateLayer.lean`)

| ID | 定理 | 主張 |
|---|---|---|
| **F-17** | `velocity_radial_error_le`, `velocity_axial_error_le` | ポテンシャル誤差 → 速度誤差 |
| **F-18** | `product_difference`, `product_error_le`, `advection_error_le` | 積差恒等式と移流項誤差 |
| **F-19** | `gronwallBound_le_simple`, `norm_le_simple_gronwall` | mathlib の Grönwall を `(δ+εt)e^{Kt}` へ(`K` で割らない) |
| — | `FixedBandwidthCandidate` 構造体 + 2 定理 | 固定有限帯域候補の**全仮定を構造体化**し、破綻時刻集合の空性と各時刻での到達性 |

`lake build` **8663 jobs 成功**。`sorry`・`admit`・新規 `axiom` ゼロ。
`lake env lean AxiomAudit.lean` は**全 43 定理**が
`[propext, Classical.choice, Quot.sound]` のみ。
解析的 Green 積分を一度に形式化しようとはしていない — 証明書が供給できるのは
有限個の非負上界だけなので、この層はその**合成**だけを担う。

### 区間演算で証明したこと(新規、PO-13 着手)

`src/ns_certificate_lab/snapshot_certificate.py`。**単一 snapshot** について、
すべて**厳密有理数演算**(binary64 値は二進有理数なので入力は無近似)で:

| 量 | 結果 |
|---|---|
| Poisson 残差上界 | **1.626e-19**(25×49 格子、382 内部節点) |
| 発散上界 | 3.479e-08 |
| エネルギー区間 | 浮動小数点値を包含 |
| `L³` ノルム区間 | 立方根の整数二分法による包含 |
| 移流項区間 | `[−3.605e-08, 3.605e-08]` |
| 独立 checker | **10 検査すべて合格** |

生成器と checker は分離されており、checker は保存された有理数データだけから
不等式を再検証する(元の浮動小数点配列を一切読まない)。改竄したペイロードが
拒否されることもテストで確認した。
**限界**: 包含しているのは**離散量**であって連続解ではない。離散化誤差そのものは
包含していないので **PO-05 は未解決のまま**である。

### 数値的に観測したこと

`outputs/whole_space_gate6_v1`。

**多重極階層**(閉形式参照解、`R=5`): monopole 誤差 > dipole > quadrupole、
全段で対応する tail 上界が支配。

**Gate 6 校正表**(振幅 10、`ω₁` コア差の相対値)

| box | `n_r` | zero↔mono | mono↔dip | dip↔quad | Richardson | 比 |
|---|---:|---:|---:|---:|---:|---:|
| 2.50×2.50 | 193 | 1.66e-24 | 1.99e-07 | 8.5e-27 | 9.78e-04 | 2.03e-04 |
| 1.60×1.90 | 193 | 2.22e-18 | 6.44e-07 | 4.5e-24 | 4.47e-04 | 1.44e-03 |
| 1.35×1.65 | 193 | 7.07e-25 | 1.375e-06 | 7.1e-17 | 1.74e-04 | **7.93e-03** |

境界差は**解像度に依存しない**(3 桁一致、解像度間広がり `9.97e-2`)。

**一因子分離**(振幅 10、`1.6×1.9` box): joint 細分差 `1.17e-2 → 1.66e-3`、
`dr` のみ `3.19e-7 → 8.02e-8`、`dz` のみ `1.17e-2 → 1.66e-3`、
時間細分差 `3.13e-15 → 1.28e-15`、積分器差 `7.11e-13`、
`R_max` 広がり `1.01e-8`、`Z_max` 広がり `4.23e-4`。

**領域拡大**: トリガ発火、`ω_max` は完全一致、エネルギー相対差 `2.3e-9`、
`L³` 相対差 `4e-10`、補間欠損 `2.66e-15`(格子座標の丸め)。

**振幅・形状継続**(32 組合せ、上位 3 を 3 解像度へ昇格):
すべての候補で `L³` は**減少**(`0.9977`)、臨界幅はほぼ不変、
実効 shell 数の変化 `+0.0073`。**上位 3 候補すべて明示的に棄却**
(理由: `global_l3_increases` 不成立、`shell_count_increases` 不成立、
`tail_bound_below_signal` 不成立)。

**非線形 tail 伝播**: Grönwall 上界 `3.489e-03` < 信号 `8.027e-03`(比 0.43)。

### 合格した gate

`spatial_refinement_converges`、`time_refinement_converges`、
`integrators_agree_within_1e_3`、`radial_domain_independent_within_1e_3`、
`axial_domain_independent_within_1e_3`、`expansion_trigger_fires`、
`expansion_preserves_invariants`、`snapshot_certificate_verifies`、
`tail_propagation_bound_below_signal`、`every_candidate_resolved`(全 11 件)。

### 失敗した gate(2 件。閾値は結果を見た後に変更していない)

**1. `boundary_difference_exceeds_richardson`。** 前登録基準は「境界条件差 ≥
Richardson 推定の 8 倍」。最良値は **`7.93e-3`**、要求の **1/1000**。
比が 8 に達するには 2 次収束のまま更に 5〜6 段の細分が要り、実行不能。
**解釈**: この設定で律速なのは離散化であって外側境界ではない。良い境界条件を
使う限り「境界誤差が支配する領域」は存在せず、基準はその領域を要求している。

**2. `continuation_left_the_quadratic_regime`。** `max|ω₁|` は
`A = 2→5→10→20` で厳密に `A²` 比例(相対残差 `2.9e-6`, `1.0e-5`, `4.1e-5`)。
第一 Picard 反復が支配しており、**順位付けは非線形挙動を一切見ていない**。
どの候補も実力で昇格することはあり得なかった。離れるのに必要なのは振幅ではなく
**時間**であり、現在 `max|ω₁|·T ≈ 7e-4 ≪ 1`。

**さらに**: 境界条件差が測れるほど box を絞ると、源が評価球に到達して
**多重極 tail 上界がそもそも存在しない**。`recover_free_space_velocity` は
`tail_bound_available = False` を返して数値を捏造しない。これは 2 つの
Gate-6 基準の間の実在する緊張である。

### 作成した候補

**明示的初期値族 1 つ**(`SwirlFamily`)を固定し、`(A, R₀, Z₀, c)` の
32 点で継続を実行した。**特異点候補は 1 つも生成されていない。**

### 棄却した候補

継続の上位 3 候補 `A2_R1.2_Z1.5_c2`、`A5_R1.2_Z1.5_c2`、`A10_R1.2_Z1.5_c2` を
3 解像度で検証し、**すべて棄却**した。棄却理由は各候補について
`global_l3_increases`(`L³` は増加せず減少)、`shell_count_increases`(shell 数は
実質不変)、`tail_bound_below_signal`(pilot box が tail 上界の存在条件を満たさない)。
残る 29 点は複合スコアが下位で昇格対象にならなかった。

### 変更・追加したファイル

新規: [tail_propagation.py](src/ns_certificate_lab/tail_propagation.py)、
[initial_data.py](src/ns_certificate_lab/initial_data.py)、
[domain_expansion.py](src/ns_certificate_lab/domain_expansion.py)、
[snapshot_certificate.py](src/ns_certificate_lab/snapshot_certificate.py)、
[test_gate6_modules.py](tests/test_gate6_modules.py)(20 件)、
[run_whole_space_gate6.py](experiments/run_whole_space_gate6.py)、
[whole_space_gate6.json](configs/whole_space_gate6.json)、
`outputs/whole_space_gate6_v1/`、
[CertificateLayer.lean](formal/NSSingularity/CertificateLayer.lean)。

更新: [free_space_recovery.py](src/ns_certificate_lab/free_space_recovery.py)
(`A_6`、四重極、`tail_bound_available`)、
[whole_space_gate.py](src/ns_certificate_lab/whole_space_gate.py)(境界 4 種)、
`formal/NSSingularity.lean`、`formal/AxiomAudit.lean`、`formal/README.md`、
`docs/final_target.md`、`docs/formalization_map.md`、
`docs/whole_space_transition.md`、`PLAN.md`(Phase 2.87)、
`README.md`(再現コマンド 13)、`.github/workflows/tests.yml`、本書。

### テスト結果

```text
PYTHONPATH=src .venv/Scripts/python.exe -m pytest -q
  -> 907 passed, 1 skipped(第 5 便時点)
  -> 927 passed, 1 skipped(本便追加後)
```

`lake build`: **Build completed successfully (8663 jobs)**。
`lake env lean AxiomAudit.lean`: 全 43 定理が古典公理のみ。

### 未解決の証明義務・限界

- **非線形 tail 伝播の Grönwall は解作用素ノルムを入力として要求する。**
  Biot–Savart は `L^∞→L^∞` で有界でないので、既定値 0 は「tail 寄与のみを
  抑える」意味であり、状態誤差経由の寄与は抑えていない。
- 区間証明書は**離散量**の包含であり、離散化誤差は未包含。PO-05 は未解決。
- 区間証明書は単一 snapshot のみ。時間発展の証明書は未着手。
- 振幅継続が二次応答領域を離れていないので、複合ゲートは非線形挙動を
  判定していない。長時間積分の安定性と計算量が次の設計判断である。
- 外側境界条件は「律速ではない」ことが分かっただけで、
  境界差が離散化誤差を超える設定での検証は依然できていない。
- F-7c、`ClayStatement.lean` への橋、`L^∞` 最大値原理の形式化は未着手。

### Clay 問題について

**本便は Navier–Stokes ミレニアム問題を解決していない。** Gate 6 は
入口の計装と校正であり、2 つの前登録基準に**不合格**である。
振幅継続は実行したが、二次応答領域を離れておらず、**候補は 1 つも昇格せず
上位 3 件はすべて明示的に棄却**した。区間証明書は離散量の包含にとどまる。
過去の結果は書き換えていない。

## 数学的に確認できたこと

- 外力なしの3次元非圧縮 Navier–Stokes 方程式から、指定した curl と
  stream-function 規約の下で軸対称・旋回ありの成分式を導出した。
- \(u_1=u^\theta/r\)、\(\omega_1=\omega^\theta/r\)、
  \(\psi_1=\psi^\theta/r\) の閉じた系、粘性の \(3/r\) 係数、source の符号、
  速度回復、楕円式を代数導出し一次資料と照合した。
- 物理3次元発散
  \(\partial_ru^r+u^r/r+\partial_zu^z=0\) と、形式的な5次元スカラー作用素
  \(\mathcal L_5\) を分離した。
- 軸上の偶奇性・極条件、軸作用素の極限、物理次元、Navier–Stokes
  スケーリング、Cartesian復元式、物理エネルギー測度を監査した。
- 後方一尺度自己相似、Type I、BKM型、Serrin型、軸対称旋回なし、
  有限/局所エネルギー、渦度可積分性に関する既知障害を、一次資料を優先して
  仮定と適用範囲付きで整理した。本文まで確認したものと書誌・要旨のみを
  確認したものは `REFERENCES.md` で区別している。
- 旧Poisson試作の \(-\mathcal L_5\) 符号、内部 \(3/r\) 係数、軸行の係数8、
  周期 \(z\)、外側Dirichlet identity行、manufactured pairは静的に再導出して
  整合を確認した。ただし旧solver自体の正しさや安定性は認証していない。

これらは方程式・既知定理の確認であり、特異点の存在確認ではない。

## 数値的に確認できたこと

### Manufactured solution

独立に書いた解析微分を参照値とし、\(N_r=17,33,65\) で全診断の誤差が
減少した。

| 診断 | 誤差 \(17\to33\to65\) | 隣接観測次数 |
|---|---|---|
| 速度回復 RMS | 6.223e-3 → 1.404e-3 → 3.267e-4 | 2.148, 2.103 |
| 物理発散 RMS | 1.128e-2 → 2.364e-3 → 5.097e-4 | 2.255, 2.214 |
| 楕円 defect RMS | 3.963e-2 → 9.131e-3 → 2.169e-3 | 2.118, 2.074 |
| Cartesian復元後の独立curl defect RMS | 1.065e-5 → 6.778e-7 → 4.272e-8 | 3.973, 3.988 |
| \(u_1\) forced residual RMS | 1.428e-2 → 3.532e-3 → 8.768e-4 | 2.015, 2.010 |
| \(\omega_1\) forced residual RMS | 4.690e-2 → 1.158e-2 → 2.862e-3 | 2.018, 2.017 |

軸 parity は合格し、保存した配列・設定・seed・診断はchecksum検証後に
同一値で再読込できた。candidate/run-config v2は単位、無次元化、物理時刻、
粘性、基底規約、Python/NumPy/platform、Git状態、実行入力のsource
fingerprintを記録し、同一runの3成果物が同じpre-write provenanceを持つことも
確認した。

### 一様Cartesian独立検証

`cartesian_validation.py` はNumPy以外の数値演算実装を共有せず、一様
\((x,y,z)\) 格子上で3成分を直接差分する。周期解析場の
\(12^3,24^3,48^3\) refinementでは次の隣接観測次数を得た。

| 診断 | 観測次数 |
|---|---|
| divergence | 1.925, 1.981 |
| full curl | 1.944, 1.986 |
| vector Laplacian | 1.964, 1.991 |
| advection | 1.946, 1.986 |
| pressure gradient | 1.985, 1.996 |
| viscous term | 1.964, 1.991 |
| forced primitive defect \(R_0-f\) | 1.951, 1.988 |

保存candidate \(65\times128\) を再読込して一様
\(25\times25\times64\) 格子へ復元した検査では、

- divergence RMS/max: `2.613327e-3 / 6.328976e-3`
- full-curl defect RMS/max: `9.369950e-3 / 2.975832e-2`

となり、RMSと最大誤差の両gateを通った。adapter出力をtest側の閉形式
Cartesian oracleへ直接比較すると、source grid
\(33\times64\to65\times128\) で速度誤差は
`6.810473e-4 → 1.701222e-4`（次数 `2.001`）、全渦度誤差は
`2.501129e-3 → 6.240860e-4`（次数 `2.003`）へ減少した。

円柱radial符号、成分写像、保存 `omega1` の符号、発散汚染に加え、RMSだけ
なら埋もれる一点curl故障と、周期 \(z\) seamの故障も拒否した。非周期
one-sided closureは全境界を含む二次多項式exact testに合格した。

### 非特異基準実験

滑らかな旋回のみのガウス拡散解を、production operatorを呼ばない独立な
Crank–Nicolson法で計算した。

- \(N_r=33,65,129\) の相対 \(L^2\) 誤差:
  `8.804028753e-4`, `2.209989848e-4`, `5.531034051e-5`
- 観測次数: `1.994124`, `1.998419`
- 最細格子のエネルギー/単位 \(z\) 長:
  `0.3926991427 → 0.3562013981`
- 最大の相対エネルギー増加: `0.0`
- 半径 \(R=2,3\) の同次外側境界を比較した \(r\le1\) の最大差:
  `8.94057e-11`
- peak physical vorticity は持続増大せず、blow-up fit は開始条件で拒否。
- 全7 acceptance checks: 合格。
- \(L_0=U_0=1\)、\(L_0/U_0=1\)、\(Re=20\)、単位、有限周期
  \(L_z=2\pi\) を設定snapshotとsummaryへ記録。

これは既知に滑らかな負の対照であり、一般解の正則性を証明しない。

### 固定空間格子の時間収束

非特異ガウス旋回拡散対照を \(N_r=513\)、\(R=5\)、\(T=1\) に固定し、
\(\Delta t=0.5,0.25,0.125\)（step数 `2,4,8`）だけを変更した。

| \(\Delta t\) | 重み付き相対 \(L^2\) 誤差 | 最大絶対誤差 | 最終energy/\(z\)長 | 最終最大渦度 |
|---:|---:|---:|---:|---:|
| 0.5 | `8.5860759e-4` | `1.6269454e-3` | `0.2725504862` | `1.3856349981` |
| 0.25 | `2.0489425e-4` | `3.8828669e-4` | `0.2726696354` | `1.3881123155` |
| 0.125 | `4.3822753e-5` | `8.4939737e-5` | `0.2726992741` | `1.3887190094` |

解析解誤差の観測次数は `2.067119, 2.225128`、同一格子上のstep-doubling
差は `6.537258e-4, 1.611782e-4`、その観測次数は `2.020029` だった。
全runのenergyは初期値 `0.3926990819` から減少し、履歴中の最大渦度は
初期値 `2.0` を超えなかった。

補助的な同次境界 \(R=3,4\) 比較で \(r\le1.5\) の全時刻最大差/最終差は、
各刻みについて `2.155502e-8`, `5.382548e-9`, `2.031933e-9` だった。
これは主計算 \(R=5\) の打切り誤差を直接評価または証明する値ではない。
全11 acceptance checksは合格した。

### 自動テスト(歴史的記録: 2026-07-27 時点のマイルストーン)

**注意: 以下の `69 passed` は 2026-07-27 のマイルストーン時点の歴史的
記録である。現在の正式なテスト数は本書冒頭「2026-07-28 セッションの
追加結果」を参照(統合直後 119 → CV 追加後 146 → 非線形ソルバ追加後
193、以降のセッション内追加はセッション末尾の記録が正)。**

2026-07-27 時点の結果:

```text
69 passed in 3.21s
```

要求された故障注入はすべて検出した。

- 発散ゼロを壊した速度
- 符号反転した楕円関係
- 軸条件を破った場
- 不正/改変候補archive
- 改変診断データ
- 非収束の解像度系列
- 既存証拠ディレクトリの上書き試行
- 再署名された非canonical `float32` candidate
- 再署名された不正provenance、JSON/CSVのNaN・Infinity、JSON重複key
- 独立4次監査に不足する4点radial grid
- Cartesian一点故障をRMSだけで見逃す判定
- periodic \(z\) seamの故障
- 保存済み `omega1` の符号反転
- 非正のCartesian RMS/max許容差
- 不正な時間刻み系列とtime-convergence成果物の上書き
- time-convergence manifest欠落・hash不整合を許す実装回帰

## 仮説にすぎないこと

- Type II、異方的二尺度、周期軌道、準定常軌道、connecting orbit が
  有望な探索空間であるという設計判断。
- 動的再スケーリングで安定した候補力学が見つかる可能性。
- AIが探索初期値や低次元構造の提案に役立つ可能性。

これらは未実装の研究案であり、候補の存在を示していない。

## 未確認・未解決

- 元の3次元 Navier–Stokes 方程式に有限時間特異点が存在するか。
- 非自明な候補profile、軌道、または滑らかな初期データからの接続。
- Hou の公表増幅率(1536² 適応格子)と本実装(一様固定格子)の定量一致。
  一様格子は Hou の最小格子幅 \(O(10^{-8})\) に遠く及ばない。
- 適応 mesh、半周期 sine 対称実装、非 Fourier の独立 z 方向経路。
- 全空間の領域打切り誤差、楕円Green tail、スペクトル尾部の厳密評価。
- 圧力回復、射影、または原始変数時間発展を別実装すること。
- 候補近傍の非線形安定性、不安定方向、finite physical time、物理ノルム発散。
- 区間演算、validated inverse、radii polynomial、形式証明
  (`docs/formalization_map.md` に Lean 化対応表を開始済み)。
- 古典的な軸対称旋回なし定理の原ロシア語関数空間を、現在のSobolev記法へ
  完全に逐語対応させる作業。
- 二段階粘性を含む Hou プロトコルの \(t_0\) 以後(第 2 粘性段)の再現、
  および壁依存性実験(壁半径を広げた系列)。

## 実装したもの

- install可能なNumPy中心のPython packageとpytest設定。
- uniform axisymmetric grid、2次の \(r,z\) 微分、軸極限。
- stream functionからの \(u^r,u^z\) 回復。
- 物理3D発散、楕円defect、PDE項別残差。
- 保存候補からのCartesian速度復元、独立4次stencilによる
  \(\omega^\theta=\partial_zu^r-\partial_ru^z\) 検査。
- 一様Cartesian gridと、既存円柱operatorを呼ばない独立2次gradient、
  3成分divergence、full curl、vector Laplacian。
- primitive residualの `time`, `advection`, `pressure_gradient`,
  `viscous=-nu*Laplacian(u)`, unforced `total`, forced defect `total-f` の
  成分別配列。
- 保存candidateのverified load後に、固有の \(r,z\) stencilとbilinear補間で
  E-18a速度・E-18b全渦度を一様Cartesian格子へ写すadapter。
- 軸正則性の必要条件検査。
- canonical `<f8` explicit NPZ candidate、v2 schema manifest、
  array/archive SHA-256、必須単位・無次元化・物理解釈。
- config/seed、JSON/CSV diagnostics、checksum、pre-write runtime/source
  provenance。v1は由来欠如を明示するread-only互換。
- 独立解析式を持つ manufactured solution。
- 独立Crank–Nicolson非特異対照、CSV/JSON/NPZ/SVG、manifest。
- 同一空間格子の \(\Delta t,\Delta t/2,\Delta t/4\) 時間収束、解析誤差、
  step-doubling、energy、最大渦度、補助境界感度を保存する実験。
- 旧ZIPの指定6ファイルだけをread-only監査した
  `docs/legacy_reuse_review.md`。旧コードは移植していない。
- 正常系、round-trip、tamper、故障注入、自動上書き防止テスト。
- Python 3.10/3.12でテストと3実験を再生するGitHub Actions workflow。

## 実行した主なコマンド

2026-07-28 セッション(venv `.venv`、Python 3.11.9、Windows 11):

```text
PYTHONPATH=src PYTHONDONTWRITEBYTECODE=1 .venv/Scripts/python.exe -m pytest -q
  -> 119 passed(統合直後)/ 146 passed(CV 追加後)/ 193 passed(非線形ソルバ追加後)
PYTHONPATH=src .venv/Scripts/python.exe -m experiments.run_poisson_gate --config configs/poisson_gate.json --output-dir outputs/poisson_gate_fable5
PYTHONPATH=src .venv/Scripts/python.exe -m experiments.run_hou_early_time --config configs/hou_early_time.json --output-dir outputs/hou_early_time_v1
```

それ以前(2026-07-27 まで):

```text
$env:PYTHONPATH='src'; $env:PYTHONDONTWRITEBYTECODE='1'; python -m pytest -q
$env:PYTHONPATH='src'; $env:PYTHONDONTWRITEBYTECODE='1'; python experiments/run_manufactured.py --config configs/manufactured.json --output-dir outputs/manufactured_v5
$env:PYTHONPATH='src'; $env:PYTHONDONTWRITEBYTECODE='1'; python experiments/run_baseline.py --config configs/baseline.json --output-dir outputs/baseline_v5
$env:PYTHONPATH='src'; $env:PYTHONDONTWRITEBYTECODE='1'; python -m experiments.run_time_convergence --config configs/baseline_time_convergence.json --output-dir outputs/time_convergence_v1
python -m compileall -q src experiments tests
git diff --check
```

加えて、Pythonから candidate/config/JSON/CSV をchecksum付きで再読込し、
baseline manifestの全7 payload、time-convergence manifestの全5 payloadに
ついてhash・byte lengthとmanifest sidecarを検査した。両NPZの有限性と
`allow_pickle=False` 読込、両SVGのXML妥当性、全summary acceptanceも検証した。
3つの現行runのsource fingerprintは
`bac2077ff5e7333e6c0201d4ebcc319c2fb93604fdd2c8c721bad6e01f8333cf`
で一致した。

時間収束scriptを最初にfile pathとして直接起動した試行は、sibling
`experiments` namespaceを解決できず実験開始前に終了した。証拠directoryは
作られなかった。再現コマンドとCIを `python -m experiments.run_time_convergence`
へ統一して再実行し、上記成果物を生成した。

## 出力ファイル

現行成果物 `outputs/manufactured_v5/`:

- `diagnostics.json`, `diagnostics.csv` と各SHA-256 sidecar
- `run_config.json` とSHA-256 sidecar
- `manufactured_candidate.npz`
- `manufactured_candidate.manifest.json` とSHA-256 sidecar

現行成果物 `outputs/baseline_v5/`:

- `summary.json`, `convergence.csv`, `diagnostics.csv`
- `profiles.npz`
- `energy.svg`, `profiles.svg`
- `config.snapshot.json`, `manifest.json` とmanifest SHA-256 sidecar

現行成果物 `outputs/time_convergence_v1/`:

- `summary.json`
- `time_convergence.csv`, `time_diagnostics.csv`
- `final_profiles.npz`
- `config.snapshot.json`
- `manifest.json` とmanifest SHA-256 sidecar

2026-07-28 追加成果物:

- `outputs/poisson_gate_fable5/` — 統合ツリーでの Poisson ゲート新規実行
  (全 7 受入合格、manifest+payload SHA-256 検証済み)
- `outputs/poisson_gate_v1_bundle_snapshot/` — バンドル同梱の旧 snapshot
  証拠(改変なし保存)
- `outputs/hou_early_time_v1/` — 早期 Hou 実行(summary、diagnostics.csv、
  snapshots.csv、trajectories.npz、3 解像度×5 checkpoint、manifest+
  全 50 payload SHA-256 検証済み)
- `archive/poisson_gate_packaging/` — 統合済みバンドルのパッケージング残骸
  (provenance 用 README 付き)

`outputs/manufactured/`, `outputs/baseline/`, `outputs/*_v2/`,
`outputs/*_v3/`, `outputs/*_v4/` は、途中段階の結果を隠さないため保存している。最初の
`manufactured/` はlegacy v1でprovenanceを持たず、現行成果物には用いない。

## 既知の問題

- package operatorの外側 radial stencil は診断用で、全空間境界条件や
  楕円solveを提供しない。
- 旧Poisson試作は指定ファイルの静的監査だけを行った。旧Grid、boundary
  helper、旧 `l5` は指定外のため未確認で、solverは移植していない。
- axis checkは有限個の必要条件を検査するだけで、Cartesian滑らかさの証明ではない。
- manufactured fieldは強制付きであり、未知の無外力解ではない。
- 一様Cartesian checkerは有限box上の2次binary64差分で、candidate adapterは
  bilinear補間を使う。観測収束と故障検出はあるが、補間・離散化・領域打切りの
  厳密上界はない。圧力は入力であり、独立pressure solveはない。
- baselineは旋回のみの特殊な滑らかな対照で、非線形meridional dynamicsを
  試していない。有限長の \(z\)-周期円柱では有限エネルギーだが、同じ
  \(z\)-不変場を \(\mathbb R^3\) 全体へ延長すると総エネルギーは無限なので、
  主対象の全空間有限エネルギー対照ではない。
- 時間収束の主系列は \(R=5\) で固定したが、保存した境界感度は別の
  \(R=3,4\) 補助比較である。主領域の打切り誤差boundではない。
- SHA-256は改変検出であって、数値の正しさ・作者・実数包含を証明しない。
- 現行v5/v1成果物のprovenanceは `git_head=1af29b9...` と
  `git_dirty=true` を記録する。今回の未commit変更をsource fingerprintで
  固定しているが、署名や信頼時刻証明はない。
- `float64` の誤差は測定しただけで、外向き丸め区間にはなっていない。
- 非収束の故障注入は判定器へ与える合成誤差系列である。将来のproduction
  solverには、意図的に壊した時間発展を通すend-to-end拒否試験も必要。
- 文献の適用はまだ存在しない将来候補に対しては行えない。
- (2026-07-28 追加)出荷済み run の von Neumann 監査は stride 25 の
  記録行のみを被覆する。行間の step は未監査であり、将来 run の
  `step_stream` 全 step 被覆でのみ閉じる。
- (2026-07-28 追加)193×384 run は記録行監査でも stability-unverified
  (T₁ 直前 4 行、max|G| ≤ 1.0000312)。65×128 の Gate 1 相互比較は
  時間スキーム依存 ~6 ppm を示したが、193×384 の適応 CFL 0.1 運転点
  そのものでの相互比較は未実施。
- (2026-07-28 追加)65×128 の E-29 run では全 step 相対エネルギー収支
  defect が 0.936 に達する(滑らか control では 4.2e-3 へ収束)。現行
  解像度では離散エネルギー恒等式が front 上で閉じない。
- (2026-07-28 追加)隣接解像度の共通格子差は値では減少するが微分 L∞
  では減少しない。T₁ の全解像度が fit 前提(front ≥ 7 点)を満たさず、
  増幅ラダーの外挿は blind 判定でも not_in_asymptotic_range。
- (2026-07-28 追加)`hou_early_time` / `wall_dependence` 実験本体は
  まだ旧計装(snapshot 系 gate)のまま。次回実行前に `gate_summary`
  読み出しへの移行が必要(新規 run のみ。既存証拠は不変)。
- (2026-07-28 追加)GitHub が既定 branch に dependabot 警告 4 件
  (moderate 2、low 2)を報告している。数値結果には影響しないが未対処。
- (2026-07-29 第 3 便追加)`galerkin_obstruction.stream_apriori_bound` の
  上界違反は「故障」と「時間刻み未解像」を区別しない。区別には
  刻みに依存しない診断 `max_relative_energy_production` を読む必要がある
  (クリーン 3e-16 対 故障 1.6)。この設計上の限界は docstring に明記した。
- (2026-07-29 第 3 便追加)厳密算術の相殺証明書は列挙した有限個の `S` を
  覆うだけである。全 `S` に対する主張は Lemma 1 の紙上証明であって、
  scan の合格ではない。

## 次に行うべき最小の一手

〔改版 2026-07-29 第 3 便。最終目標までの単一依存グラフは
`docs/final_target.md`(そこに Clay A–D、Track U/F の最終定理、全証明義務、
Lean 識別子 F-1〜F-11 の確定登録簿がある)。
FABLE5_NEXT_TASK_AUDIT の Gate 順序(Gate 1–4 が通るまで中後期成長・
blow-up fit・AI 候補探索へ進まない)を引き続き最上位の拘束とする。
Gate 1 は合格、Gate 2/3 は既存証拠+新監査で部分的、Gate 4 は未実装。

**第 3 便で閉じた項目**: Track F の有限モード ansatz 族は除外定理として
閉鎖した(`docs/research_notes/track_f_finite_mode_nogo.md`、Lean `F-6`)。
`START_NEW_SESSION_NAVIER_STOKES.md` §6「優先候補A」の step 1–5
(低次 Fourier ansatz の記号探索)は**証明により空なので実行しない**。
今後この探索を再登録してはならない。〕

0. 〔第 6 便で更新〕**弱非線形領域を離れることが Track U の新しい律速である。**
   Gate 6 の振幅継続は `max|ω₁| ∝ A²` を相対 `5e-5` で満たし、二次応答領域を
   一度も離れていない。必要なのは振幅ではなく**時間**(現在
   `max|ω₁|·T ≈ 7e-4 ≪ 1`)であり、長時間積分の安定性と計算量が次の設計判断。
   外側境界は律速ではないことが測定で分かった(境界差は離散化誤差の `8e-3` 倍)。
   区間証明書は単一 snapshot から時間発展へ広げる。
   〔第 5 便の記述、参考〕**中振幅での再ゲート**
   Gate 5 は小振幅で通ったが、時間細分・領域拡大・境界次数の 3 部分ゲートは
   この振幅では情報を持たなかった(答えが桁まで一致)。場が外側境界に到達する
   振幅で 3 つを再実行するまで、**外側境界条件は未検証**であり、
   強振幅 Hou 型全空間候補へ進んではならない。
   その後 `core_width.fit_precondition`(front ≥ 7 点)を全空間 box で評価し、
   適応 mesh か半周期 sine 実装かを決める。
   Lean 側は **F-7c**(自励化場の局所流の構成)を閉じる。
   Track F の帯域幅発散族については、`Π_j` の**鋭い**上界を導く
   (現在の粗い評価は `(γ,σ)` に追加制約を与えず、cascade 模型と PDE を繋げない)。

1. **Gate 4 の実装(最優先)**: 非周期 \(z\) の有限 box、\(z\) 方向も
   \(C^\infty\) compact な初期値族、free-space 楕円経路(W-1 の \(z\)
   非周期版/Green 積分/Hankel)、\(R_{\max}\)/\(Z_{\max}\) 独立拡大、
   低波数 stress test(`docs/whole_space_transition.md` §7)。その入口
   として \(L_z\in\{1,2,4\}\) 族(既存実装の config 変更で測定可能)で
   指数→代数遷移を実測する。
2. **実験本体の新計装への移行**: `run_hou_early_time` /
   `run_wall_dependence` の受入検査を `gate_summary`(全 step streaming)
   読み出しへ切替え、`stage_cfl_limit` の使用を判断する。以後の新規 run
   はすべて全 step 被覆+3 積分器のうち 2 つ以上での交差確認を要件とする
   (Heun 単独増幅の候補判定使用は禁止済み)。
3. **軸近傍解像度の設計判断**: T₁ の fit 前提不合格(front ≤ 5.7 点)、
   相対軸パリティ 0.706、微分レベルの共通格子差の非減少はいずれも
   front 未解像を指す。適応 mesh または半周期 sine 実装の設計判断を、
   中成長段の前に行う。
4. **既存 checkpoint の primitive 監査の完結**: snapshot 前後 1 step の
   追加保存 option(既存実験へ)。
5. **Lean 段階 1 の継続**: F-1(再スケーリング恒等式)と F-4(証明書
   不等式)。F-2/F-3 と同じ方針(定義明示、非スコープ明記、
   `#print axioms` 記録、`AxiomAudit.lean` へ追記)。

これらが通っても、長時間探索、AI最適化、特異点fitへ自動的には進まない。
動的再スケーリング探索の前に、全空間tailと候補用離散化の証明可能な設計を
再評価する(PO-05〜PO-07)。
