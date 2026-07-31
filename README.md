# ns-singularity-certificate-lab

3次元非圧縮 Navier–Stokes 方程式の有限時間特異点**候補**を、将来の
区間演算・コンピューター支援証明へ接続できる形で研究するための監査可能な
基盤です。

> **重要:** このリポジトリはミレニアム懸賞問題を解決していません。
> 特異点を発見・証明しておらず、収録した実験は既知に滑らかな人工場と
> 非特異な減衰対照だけです。大きな数値、発散らしい回帰、小さい残差は
> 数学的証明ではありません。

<!-- MCR:BEGIN (machine-checkable results; guarded by tests/test_readme_claims.py) -->
## Machine-checkable results / 第三者が機械検証できる成果

第三者が短時間で「何を無条件に検証できるか・何が外部定理に依存するか・
どのコマンドを実行すればよいか・どの程度の信頼基盤が必要か」を判断する
ための節です。詳細は成果ごとに安定 ID を振った
**[Verified Results Registry](docs/verified_results.md)** にあります。

**用語規則(本 README と登録簿で厳守):**

- **Lean-verified** — Lean kernel が定理を検査し、project 固有の未証明
  axiom がない場合のみ。
- **certificate-verified** — 独立 Python checker が証明書を**完全再計算**
  する場合。Lean 証明ではない。
- **conditional PDE certificate** — EXT 等の外部仮定が残る場合。
- **numerically observed** — 浮動小数点計算のみ。証明ではない。
- **candidate** — Clay 公式命題への橋が未検証の場合。

「formally verified / fully verified / machine-checked PDE theorem」は
依存鎖が実際に閉じた場合以外では使いません。

### A. Lean 4 で無条件に検証できる結果(Lean-verified)

現在の HEAD で `sorry` / `admit` / project 固有 axiom なしに証明済み。
公理監査は全定理が `[propext, Classical.choice, Quot.sound]`(= Lean/mathlib
の通常基盤; **無公理証明の主張ではない**)のみを報告します
(`cond_to_uncond` のみ公理非依存)。全 124 定理の一覧は
[formal/AxiomAudit.lean](formal/AxiomAudit.lean) と
[登録簿](docs/verified_results.md)。主要例:

| Result | Lean theorem | File | Verified claim | Not claimed |
|---|---|---|---|---|
| Galerkin energy/norm bound | `galerkin_norm_le`, `norm_le_of_energy_inequality` | [GalerkinNoBlowup.lean](formal/NSSingularity/GalerkinNoBlowup.lean) | エネルギー不等式下の有限次元ノルム有界 | PDE のエネルギー不等式の成立 |
| 固定有限帯域軌道の no-go | `FixedBandwidthCandidate.breakdown_times_empty`, `.reaches_every_time`, `.fixedBand_scope` | [FiniteModeNoGo.lean](formal/NSSingularity/FiniteModeNoGo.lean), [TrackPFourier.lean](formal/NSSingularity/TrackPFourier.lean) | 固定帯域に**留まる軌道**は破綻時刻を持たない | 有限帯域**初期値**についての主張(反例 `exists_finiteBandDatum_not_fixedBandTrajectory` が区別を証明) |
| Fourier 三線形相殺の代数 | `advectionForm_eq_zero`, `inner_leray_eq_zero` ほか Leray 代数 5 定理 | [FiniteModeNoGo.lean](formal/NSSingularity/FiniteModeNoGo.lean), [TrackPFourier.lean](formal/NSSingularity/TrackPFourier.lean) | 発散ゼロ場での advection 形式の消滅・Leray 乗数の直交性 | 無限次元 Kato–Ponce 可換子評価(監査済み紙上証明; C 節) |
| 有限次元 Picard–Lindelöf | `galerkin_local_solution`, `quadratic_ode_local_solution`, `quadratic_ode_unique` | [GalerkinNoBlowup.lean](formal/NSSingularity/GalerkinNoBlowup.lean), [GalerkinPicard.lean](formal/NSSingularity/GalerkinPicard.lean) | 二次 ODE `u' = Au + B(u,u)` の局所存在・一意性(明示 Lipschitz 定数・存在区間半幅) | EXT-P1(PDE 命題)そのもの |
| control ODE comparison | `riccati_comparison`, `gronwall_variable_coefficient(_integral)`, `roughEnclosure_*` 3 定理 | [ControlODE.lean](formal/NSSingularity/ControlODE.lean) | スカラー Riccati 比較と rough enclosure の PL 存在・一意性 | — |
| 積分形比較(EXT-P2-INT のスカラー半分) | `integral_comparison`, `integral_riccati_comparison` | [ChainAnalysis.lean](formal/NSSingularity/ChainAnalysis.lean) | 2 パラメータ積分不等式を満たす連続 φ は ODE 解に支配される(Dini 微分不要) | PDE 差分場がその積分不等式を満たすこと(監査済み紙上証明) |
| two-slab composition | `two_slab_composition` | [TrackPChain.lean](formal/NSSingularity/TrackPChain.lean) | 抽象 tube の 2 スラブ合成(piecewise 中心・半径) | tube の PDE 的供給 |
| finite chain composition | `chain_composition`, `chain_composition_union` | [TrackPChain.lean](formal/NSSingularity/TrackPChain.lean) | `δ_out = δ_end + transfer` に忠実な n スラブ帰納法 | 〃 |
| transfer triangle inequality | `transfer_triangle` | [TrackPChain.lean](formal/NSSingularity/TrackPChain.lean) | 再中心化予算の 3 項三角不等式 | — |
| Taylor endpoint remainder | `taylor_endpoint_remainder_bound` | [TrackPChain.lean](formal/NSSingularity/TrackPChain.lean) | Lagrange 剰余 `≤ M h^{m+1}/(m+1)!` | 係数上界 `M` の数値(certificate 層) |
| Kato 証明書の有限代数 | `cube_diff_bound`, `am_gm_split`, `shifted_ratio_bound`, `inv_pow_tail_bound`, `g3_of_a4` | [KatoConstant.lean](formal/NSSingularity/KatoConstant.lean) | `G₃ ≤ 12√A₄` 組立の全代数ステップと格子 tail | 可換子評価本体(C 節) |
| EXT-P3 貼り合わせ論理 | `glued_continuous`, `exists_continuousOn_Icc_extension` | [ChainAnalysis.lean](formal/NSSingularity/ChainAnalysis.lean) | 連続貼り合わせ・完備空間での端点延長 | 延長関数が方程式を満たすこと(監査済み紙上証明) |
| Green/scaling/certificate 層 | `greenProfile_radial_laplace_eq_zero`, `physicalTime_lt_blowupTime`, `velocity_radial_error_le` ほか | [GreenAndCascade.lean](formal/NSSingularity/GreenAndCascade.lean), [FiniteTime.lean](formal/NSSingularity/FiniteTime.lean), [CertificateLayer.lean](formal/NSSingularity/CertificateLayer.lean) | 各層の有限不等式(登録簿参照) | — |
| 純粋旋回の等号ケース | `pure_swirl_equality_case`, `swirl_cartesianDiv_eq_zero` | [L3Generation.lean](formal/NSSingularity/L3Generation.lean) | L³ 生成恒等式の代数部分 | 積分論の無限次元部分 |

### B. 独立 checker で検証できる証明書(certificate-verified)

checker は payload を**完全再計算**します(builder と状態を共有しない)。
信頼基盤: Python 処理系と checker 実装(checker 論理の一部は上記 Lean 層に
鏡写し)。全証明書に改竄拒否テストがあります。

| 証明書 | 保存場所 | 実行コマンド | checker が再計算するもの | 未証明外部仮定 |
|---|---|---|---|---|
| Track P 単発スラブ | [outputs/track_p_slab_v1/](outputs/track_p_slab_v1/) | `python -m experiments.run_track_p_slab --config configs/track_p_slab.json --output-dir <dir>` | datum 検査・Picard box・全定数・control 管 | EXT-P1/P2/P3(v1 条件付き文言) |
| Track P chain(H⁴, v1) | [outputs/track_p_chain_v1/](outputs/track_p_chain_v1/) | `python experiments/run_track_p_chain.py` | 全リンク: box・定数組立・管・Taylor 終端・再中心化点の厳密一致・δ 漸化式・文言契約 | 〃 |
| Track P chain(n=3 Kato) | [outputs/track_p_chain_h3_v1/](outputs/track_p_chain_h3_v1/) | `python experiments/run_track_p_chain_h3.py` | 〃 + G₃ 証明書再検証 + `C_kato`/`C_shift` 再計算 | 監査済み紙上証明(C 節) |
| 監査済み再発行(v2) | [outputs/track_p_chain_reissued_v2/](outputs/track_p_chain_reissued_v2/) | `python experiments/reissue_chains.py` | 〃 + 閉鎖メタデータ整合・新旧文言の混在拒否 | G-DINI(未消費・open) |
| Kato 定数証明書 | [outputs/track_p_chain_h3_v1/kato_certificate.json](outputs/track_p_chain_h3_v1/kato_certificate.json) | `python -m pytest tests/test_kato_constant.py` | A₄/A₆ を自前格子ループで再計算・√ と 12 倍率・単調性 | 可換子評価本体(C 節) |
| R³ スペクトル圧力証明書 | builder/checker: `gaussian_spectral_pressure.py` | `python -m pytest tests/test_gaussian_spectral_pressure.py` | 厳密閉形式 ∇p(Δp+σ≡0 の厳密有理自己検証)と J 下界包含 | **なし(hypotheses 空)** — ただし結果は負の下界(candidate ではない) |
| snapshot / 時空スラブ(Gate 6/7) | outputs/(README 13–14 節) | 各 replay コマンド(下記) | cell 内包含・改竄拒否 | 各 payload の hypotheses 欄に明記 |

### C. 条件付き成果(conditional PDE certificate)

- **v1 証明書(第 9–10 便)**: 有限次元軌道、残差、control ODE、スラブ
  連結は機械検査済み。**真の周期 Navier–Stokes 解の存在と tube 包含は
  EXT-P1/P2/P3 に条件付き。**
- **第 11 便の監査後(再発行 v2)**: EXT-P1★/EXT-P2-INT+Lemma C/EXT-P3★/
  系 P3-3 は 3 名の敵対的監査+修理+再監査 2 名(拒否権つき)を経た
  **監査済み紙上証明**として閉鎖([監査文書](docs/research_notes/ext_p1_p2_p3_audit.md))。
  payload の `proved: true` は**この意味であり Lean 形式化ではない**
  (`lean_formalised: false` 固定、公理化は禁止のまま)。旧 Dini 節は
  G-DINI として open(どこからも未消費、checker 強制)。
- **Kato 不等式 `G₃ ≤ 12√A₄`**: 格子和は certificate-verified、代数
  ステップは Lean-verified、無限次元可換子評価は監査済み紙上証明
  ([導出ノート](docs/research_notes/kato_h3_constants.md))。
- HS-5 全空間版・NT-N1・H3 は open(条件付きのまま)。

<!-- PROMOTION:BEGIN (guard: requires all EXTERNAL_THEOREMS_AUDITED proved:true) -->
**昇格成果(第 11 便、監査閉鎖後)** — For the explicitly listed
finite-Fourier initial data, viscosity, and time interval, the repository
provides a machine-checkable certificate — modulo the audited classical
theorems named above — that a unique periodic strong Navier–Stokes solution
exists and remains within the stated Sobolev-radius tube around the
certified Galerkin approximation. / 明示された有限 Fourier 初期値・粘性・
時間区間について、**監査済み古典定理を法として**、一意な周期強解が存在し
証明書化された Galerkin 近似の Sobolev 半径 tube 内に留まることを機械検証
できる。

代表例(全 27 本は登録簿 VR-COND 系列): 初期値 = P1(helical triad、
厳密有理 Fourier 係数)、粘性 ν = 1/10、certified interval
[0, 5/256]、Sobolev 次数 Ḣ³(H³ は √8 倍)、誤差半径 ≤ 0.0309
([reissued_h3chain_strict_same_step_P1_nu_1over10.json](outputs/track_p_chain_reissued_v2/reissued_h3chain_strict_same_step_P1_nu_1over10.json))。
使用定理: EXT-P1★/EXT-P2-INT+Lemma C/系 P3-3(監査済み紙上証明)+
`chain_composition`/`integral_riccati_comparison` ほか(Lean-verified)。
役割分担: Lean = 合成・比較の骨格、checker = 全リンクの厳密有理再計算、
監査文書 = 古典解析。

**必ず併記する範囲**: 特定初期値・特定短時間区間の結果であり、大域正則性
ではなく、特異点構成ではなく、Clay 問題の解決ではない。certified horizon
は**証明手法の到達範囲**(粗い定数の Riccati 天井)であって解の性質では
ない。
<!-- PROMOTION:END -->

### D. このリポジトリが証明していないこと

| 範囲外 | 状態 |
|---|---|
| すべての滑らかな NS 初期値の大域正則性 | 主張していない |
| 有限時間特異点(存在・構成) | 主張していない |
| Clay 公式命題 (A)–(D) のいずれか | 主張していない |
| 数値軌道の長時間正確性一般 | tube 証明の範囲外は numerically observed のみ |
| 条件付き EXT 群を仮定なしで使えること | 監査済み紙上証明を法とする(trust model 参照) |
| certificate chain の終了時刻が PDE 特異時刻であること | 終了は前登録分類法で方法の限界として分類される |

## Reproduce the verified results / 再現手順

**Quick verification(〜30 分 + mathlib 初回ビルド):**

```bash
git clone https://github.com/HeliCorgi/ns-singularity-certificate-lab.git
cd ns-singularity-certificate-lab
git checkout fable5-mainline

# Lean: build + axiom audit
cd formal
lake build                      # 期待: "Build completed successfully"
lake env lean AxiomAudit.lean   # 期待: 各定理につき 1 行、標準 3 公理のみ
cd ..

# Python: 依存導入 + 全テスト
python -m pip install -e ".[dev]"
python -m pytest                # 期待: 全件 pass(scipy 不在時 1 skip)
```

テスト件数・Lean job 数は開発とともに変動します(固定値は書かない方針;
CI が毎 push で再検証)。

**Full verification(証明書リプレイ、計数時間):** README 13–14 節の各
コマンド(Track P スラブ数分、chain 各 1–2 時間、Gate 4–7 各数分〜20 分)。
それぞれ出力ディレクトリに payload + `summary.json` + sha256 `manifest.json`
を生成し、独立 checker の verdict を表示します。

## Trust model

| Layer | Trusted components |
|---|---|
| Lean theorem checking | Lean kernel、pinned mathlib([formal/lean-toolchain](formal/lean-toolchain), [formal/lake-manifest.json](formal/lake-manifest.json))、`#print axioms` が報告する標準 3 公理(propext / Classical.choice / Quot.sound — Lean/mathlib の通常基盤であり、無公理証明の主張ではない) |
| Rational certificate checking | Python 処理系と checker 実装(checker 論理のうち Lean に鏡写しされた部分は上記に還元) |
| Floating-point exploration | ハードウェア/ランタイム; **証明ではない** |
| Audited classical theorems | EXT-P1★/EXT-P2-INT/EXT-P3★/P3-3 と Kato 可換子評価: [監査文書](docs/research_notes/ext_p1_p2_p3_audit.md)の敵対的監査プロセスを信頼(Lean 形式化されるまでこの区分に留まる; 公理化は禁止) |

Lean 監査の内訳: `lake env lean AxiomAudit.lean`(全定理の `#print axioms`)、
`sorry`/`admit`/`axiom` の grep(CI)、pinned `lean-toolchain` と
`lake-manifest.json`。

陳腐化防止: [tests/test_readme_claims.py](tests/test_readme_claims.py) が
本節の定理名実在・成果物パス実在・コマンド対象実在・条件付き/無条件の
区分整合・EXT 状態と昇格文の整合を毎 CI で検査します。
<!-- MCR:END -->

## 現在できること

- 元の3次元方程式から軸対称・旋回ありの
  \((u_1,\omega_1,\psi_1)\) 系を導出・監査する。
- 二次有限差分で微分、速度回復、物理3次元発散、楕円関係、PDE項別残差、
  軸正則性を検査する。
- 保存候補からCartesian速度成分を復元し、production作用素と別の4次
  \(r,z\) stencilでazimuthal curlを照合する。
- 既存の円柱差分を呼ばない一様 \((x,y,z)\) 格子実装で、3成分divergence、
  full curl、vector Laplacian、元のprimitive Navier--Stokes方程式の
  3成分項別残差を検査する。
- candidateをいったん保存・再読込し、専用の独立stencilと補間だけで
  Cartesian物理場へ復元して、divergenceとfull curlをend-to-endで照合する。
- 実装と独立な解析微分を持つ manufactured solution で収束次数を測る。
- 明示的canonical `<f8` 配列を candidate v2 NPZ、単位・無次元化・物理時刻・
  粘性・基底規約・runtime/source provenance付きmanifest、設定、seed、
  SHA-256とともに保存・再読込する。
- 改変した診断、壊れた発散、楕円符号反転、軸条件破壊などを自動テストで
  拒否する。
- 独立なCrank–Nicolson実装で、非特異な旋回拡散対照を実行する。
- 空間格子を固定し、\(\Delta t,\Delta t/2,\Delta t/4\) の解析解誤差、
  step-doubling次数、エネルギー、最大渦度、境界感度を保存する。
- 有限円柱 \(-\mathcal L_5\psi_1=\omega_1\) を **2 つの独立実装**
  (`poisson.py` の \(r^3\)-flux 有限体積と
  `finite_cylinder_poisson.py` の非発散形直接差分)で解き、相互検証
  テスト(CV-1/2/3)で \(O(\Delta r^2)\) 一致・規約一致・故障検出を
  検査する。両者の z 方向 Fourier 処理と grid は共有であり、その範囲の
  独立性はない(明文化済み)。
- Hou (arXiv:2107.06509) の壁付き有限円柱で、完全非線形
  \((u_1,\omega_1,\psi_1)\) 系を Heun/RK2・Thom 型壁渦度条件(E-31)・
  二段階粘性(E-30)付きで時間発展させる(`nonlinear_cylinder.py`)。
  強制 manufactured・時間細分・対称性保存・循環最大原理・故障注入 5 種・
  restart のテストを持つ。
- E-29 初期値の早期 Hou 実行を複数解像度で行い、独立楕円 solver B との
  cross-check、E-02 発散残差、奇対称 defect、増幅率軌跡を保存する
  (`run_hou_early_time.py`)。

候補探索、動的再スケーリング、区間演算、厳密な打切り誤差評価、全空間
(\(\mathbb R^3\))楕円処理、圧力回復の独立実装はまだありません。

## 数学的対象

主対象は無外力・全空間 \(\mathbb R^3\) の

\[
\partial_tu+(u\cdot\nabla)u=-\nabla p+\nu\Delta u,\qquad
\nabla\cdot u=0,\qquad \nu>0
\]

で、滑らか・発散ゼロ・有限エネルギーの軸対称初期データを考えます。規約、
導出、次元、スケーリング、軸の偶奇性は [SPEC.md](SPEC.md) と
[数学的定式化](docs/mathematical_formulation.md) に、各式の状態は
[方程式監査](docs/equation_audit.md) にあります。

変換変数は

\[
u_1=u^\theta/r,\qquad \omega_1=\omega^\theta/r,\qquad
\psi_1=\psi^\theta/r
\]

です。形式的作用素

\[
\mathcal L_5=\partial_{rr}+3r^{-1}\partial_r+\partial_{zz}
\]

はスカラー式の略記であり、物理的な5次元流体を意味しません。非圧縮条件と
体積測度は常に3次元のものを使います。

## 新しい環境での再現

### 1. 取得とインストール

```console
git clone https://github.com/HeliCorgi/ns-singularity-certificate-lab.git
cd ns-singularity-certificate-lab
python -m venv .venv
```

PowerShell / Windows:

```console
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

POSIX:

```console
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

依存は実行時の NumPy と、テスト用pytestだけです。GPUは不要です。
以降のコマンドは、この仮想環境を有効にした同じshellで実行します。

### 2. 全テスト

```console
python -m pytest
```

期待結果は全件成功です。正常系だけでなく、意図的な故障を検出できた場合も
テスト成功として数えます。

### 3. manufactured solution 監査

既存証拠を上書きしないよう、新しい出力ディレクトリを指定します。

```console
python experiments/run_manufactured.py --config configs/manufactured.json --output-dir outputs/manufactured_replay
```

速度回復、物理発散、楕円関係、両PDE残差の誤差と観測収束次数をJSON/CSVへ
保存し、最細格子の明示配列を再読込可能な候補形式で保存します。これは
**強制付き滑らかな人工場**の整合性試験です。

### 4. 非特異基準実験

```console
python experiments/run_baseline.py --config configs/baseline.json --output-dir outputs/baseline_replay
```

滑らかな旋回のみのガウス場を独立なCrank–Nicolson法で減衰させ、解析解への
収束、エネルギー、外側境界感度、有限勾配を発散と誤認しない診断を保存します。
これは有限長の周期円柱上の対照であり、全空間有限エネルギー解の対照では
ありません。
スクリプトは非空の証拠ディレクトリを上書きしません。再実行には別名を使って
ください。

### 5. 固定空間格子での時間収束

```console
python -m experiments.run_time_convergence --config configs/baseline_time_convergence.json --output-dir outputs/time_convergence_replay
```

同一の \(513\) 点半径格子で時間刻みだけを
\(0.5,0.25,0.125\) と変えます。解析解への重み付き相対 \(L^2\) 誤差と
隣接次数に加え、共通の空間誤差を概ね相殺するstep-doubling差、各刻みの
エネルギー、最大物理渦度、有限領域境界感度をJSON/CSV/NPZへ保存します。
これは滑らかな負の対照の時間離散化試験であり、特異点の証拠ではありません。

### 6. 独立 Poisson ゲート

```console
python -m experiments.run_poisson_gate --config configs/poisson_gate.json --output-dir outputs/poisson_gate_replay
python -m experiments.run_poisson_manufactured --config configs/poisson_manufactured.json --output-dir outputs/poisson_manufactured_replay
```

### 7. 早期 Hou 実行(数十分〜数時間)

```console
python -m experiments.run_hou_early_time --config configs/hou_early_time.json --output-dir outputs/hou_early_time_replay
```

E-29 監査済み初期値・二段階粘性で \(t=T_1=0.002191729\) まで 3 解像度を
実行します。これは一様固定格子上の解像度制限つき数値観察であり、Hou の
適応格子計算の再現主張ではありません。

### 8. 保存 snapshot の独立 Cartesian 検査

```console
python experiments/run_hou_snapshot_cartesian_audit.py --config configs/hou_snapshot_cartesian_audit.json --output-dir outputs/hou_snapshot_cartesian_audit_replay
```

円柱演算子を使わない独立経路で、保存 checkpoint の発散・full curl・
渦度一致を相対化指標付きで検査します。

### 9. 時間刻み収束(強非線形 Hou 実行)

```console
python -m experiments.run_hou_time_refinement --config configs/hou_time_refinement.json --output-dir outputs/hou_time_refinement_replay
```

同一空間格子・同一終了時刻で固定 \(\Delta t,\Delta t/2,\Delta t/4\) を
比較し、時間誤差と空間誤差を分離します。

### 10. Track F 有限モード除外証明書(数秒)

```console
python -m experiments.run_track_f_finite_mode_scan --config configs/track_f_finite_mode_scan.json --output-dir outputs/track_f_finite_mode_scan_replay
```

滑らかな外力を使う Clay (C)/(D) 反例の「有限モード ansatz」族について、
三線型形式の相殺 \(\langle u,(u\cdot
abla)u
angle=0\) を**厳密整数演算**で
検証し(浮動小数点を一切使わない)、除外判定を出力します。これは探索の
陰性結果ではなく**除外定理**です
([docs/research_notes/track_f_finite_mode_nogo.md](docs/research_notes/track_f_finite_mode_nogo.md))。

### 11. 全空間 Gate 4(線形楕円ゲート、約 1 分)

```console
python -m experiments.run_whole_space_gate4 --config configs/whole_space_gate4.json --output-dir outputs/whole_space_gate4_replay
```

非周期 \(z\) の有限 box 上で \(-\mathcal L_5\psi_1=\omega_1\) を解き、
**閉形式の厳密な自由空間参照解**に対して格子細分・領域拡大・尾部上界・
周期像分離・独立 Cartesian 検査を測定します。軸方向は FFT を使わない
密な離散サイン変換で、既存ソルバと規約を共有しません。これは**線形**
ゲートであり、非線形発展について何も主張しません
([docs/whole_space_transition.md](docs/whole_space_transition.md))。

### 12. 全空間 Gate 5(微分 tail・速度回復・小振幅非線形、数分)

```console
python -m experiments.run_whole_space_gate5 --config configs/whole_space_gate5.json --output-dir outputs/whole_space_gate5_replay
```

Green 核の解析微分から導いた**微分 tail 上界**を閉形式参照解に対して検査し、
自由空間速度回復 API の空間・領域収束、軸正則性、独立 Cartesian 検査、
故障注入を測定し、滑らか・コンパクト台・発散ゼロの**小振幅純粋旋回**初期値から
非周期 \(z\) の全空間非線形短時間発展を回します。最後に、低周波のみの
滑らかな外力が非線形 triad 経由で高シェルを駆動しうるかを有限 cascade 模型で
判定します
([docs/research_notes/green_derivative_tail_bounds.md](docs/research_notes/green_derivative_tail_bounds.md)、
[docs/research_notes/cascade_toy_model.md](docs/research_notes/cascade_toy_model.md))。

### 13. 全空間 Gate 6(中振幅校正・振幅継続・区間証明書、十数分)

```console
python -m experiments.run_whole_space_gate6 --config configs/whole_space_gate6.json --output-dir outputs/whole_space_gate6_replay
```

境界条件 4 種(zero / monopole / dipole / quadrupole)の core 差を Richardson
離散化誤差と比較して校正し、`dr`/`dz`/joint/`dt`/積分器/`Rmax`/`Zmax` を一因子ずつ
分離し、明示的な初期値族について振幅・形状継続を実行して複合ゲートで順位付けし、
動的領域拡大と**厳密有理数区間演算による snapshot 証明書**を生成・独立検査します。
**2 つの前登録基準は不合格として記録されます**
([docs/whole_space_transition.md](docs/whole_space_transition.md))。

### 13.5 Track P 周期スラブ証明書(数分)

```console
python -m experiments.run_track_p_slab --config configs/track_p_slab.json --output-dir outputs/track-p-replay
```

周期 T³ 上の有理 Fourier 初期値 3 族(P1/P2/P3)について、厳密有理数演算で
Galerkin 軌道の Picard 包含・厳密な連続 NS 残差(= Galerkin tail)・H⁴ control
不等式・control ODE 管を組み立て、「真の周期強解がスラブ全体に存在し
‖u−u_a‖_Ḣ⁴ ≤ R(t)」の条件付き証明書(古典外部定理 EXT-P1/P2/P3 は忠実記録、
Lean 公理化なし)を 12 スラブ分生成・独立検査します
([docs/research_notes/track_p_periodic.md](docs/research_notes/track_p_periodic.md))。
**これは特異点証明ではなく、軌道近傍の正則性の証明です。**

### 13.6 Track P スラブ連結(certified horizon、約 1 時間)

```console
python experiments/run_track_p_chain.py
```

第 9 便の単発スラブをスカラー H⁴ 誤差半径で連結します: 各スラブは**厳密有理
再中心化点**から開始し(区間 box はスラブ境界を越えて伝播しない = wrapping の
入る場所が構造的にない)、Taylor 終端包絡+dyadic 丸め+厳密 Leray 射影で捨てた
幅は δ_{n+1} = R_n(h) + transfer としてスカラー半径に課金されます。P1/P2/P3 ×
ν ∈ {1/4, 1/10, 1/40, 1/100} の 12 連結 + 長尺 1 本(h = 1/8192、48 スラブ予算)
を前登録 config([configs/track_p_chain.json](configs/track_p_chain.json))で実行し、
独立 checker が全リンクを再計算します。停止は前登録分類法で必ず分類され、
**証明区間の終了は特異点の主張ではありません**
([docs/research_notes/track_p_chain.md](docs/research_notes/track_p_chain.md))。

### 13.7 Track P チェーン n=3(Kato 定数、約 2 時間)+ 監査済み再発行

```console
python experiments/run_track_p_chain_h3.py
python experiments/reissue_chains.py
```

第 11 便: 正規化完全一致で自前導出した `G₃ ≤ 12√A₄`
([docs/research_notes/kato_h3_constants.md](docs/research_notes/kato_h3_constants.md)、
独立 checker 付き証明書)と厳密帯域和 `C_kato`/`C_shift` による n=3 control
不等式でチェーンを再実行し、certified horizon を旧 `9(K₁+K₂)` 比で実測
約 11〜13 倍に延長。EXT-P1/P2/P3 は 3 名の敵対的監査+組立修理+再監査 2 名
(拒否権つき)を経て**監査済み紙上証明として閉鎖**
([docs/research_notes/ext_p1_p2_p3_audit.md](docs/research_notes/ext_p1_p2_p3_audit.md))、
全チェーンを閉鎖メタデータつきで再発行(`proved:true` = 監査済み紙上証明の
意味であり Lean 形式化ではない; Lean 公理化は不変で禁止; checker が新旧の
混在を拒否)。**特異点主張ではない。**

### 14. 全空間 Gate 7(Picard 領域からの離脱・τ/Re 継続・時空スラブ証明書、数分)

```console
python -m experiments.run_tau_continuation --config configs/tau_continuation_gate7.json --output-dir outputs/tau_continuation_gate7_replay
```

第 6 便の 32 点スイープを無次元座標 `(Re, aspect, c, τ)` へ再分類し(到達 `τ` は
最大 0.0233 だったことが判明)、Picard 梯子(level 0/1/2 + 完全解を同時積分)で
第一 Picard 反復からの乖離を**実測**し、前登録 τ = {0.025 … 1.0} と
Re = {10 … 400} × 族 S/A/H の 18 run を実行し、乖離ゲート 9 項目と昇格 2 基準で
判定し、`[t_n, t_{n+1}]` の**時空スラブ証明書**(cell 内部・全時刻を包含、
厳密有理数、独立 checker + 改竄拒否)を生成します。
**乖離ゲートは全項目合格、昇格候補はゼロ**
([docs/research_notes/tau_continuation_gate7.md](docs/research_notes/tau_continuation_gate7.md))。

より厳密な再現プロトコルは [docs/reproducibility.md](docs/reproducibility.md)
を参照してください。

## 現行の小規模結果

Manufactured solution の \(17/33/65\) radial refinement で得た隣接観測次数:

| 診断 | 観測次数 |
|---|---:|
| 速度回復 | 2.148, 2.103 |
| 物理3D発散 | 2.255, 2.214 |
| 楕円関係 | 2.118, 2.074 |
| 円柱samplingからCartesian復元後の独立azimuthal curl | 3.973, 3.988 |
| \(u_1\) 強制付き残差 | 2.015, 2.010 |
| \(\omega_1\) 強制付き残差 | 2.018, 2.017 |

これとは別に、一様Cartesian格子上の解析的manufactured fieldでは、
3成分divergence、full curl、vector Laplacian、移流、圧力勾配、粘性項、
primitive残差がすべて約2次で収束しました。保存・再読込した軸対称候補の
end-to-end検査では、divergence RMS/maxが
\(2.613327\times10^{-3}/6.328976\times10^{-3}\)、full-curl defect
RMS/maxが \(9.369950\times10^{-3}/2.975832\times10^{-2}\) でした。
これらは有限格子上の許容差付き数値検査であり、恒等式の厳密証明では
ありません。

非特異基準実験の \(33/65/129\) refinement では、相対 \(L^2\) 誤差

\[
8.8040\times10^{-4},\quad2.2100\times10^{-4},\quad5.5310\times10^{-5}
\]

と観測次数 \(1.994,1.998\) を得ました。エネルギーは減少し、内側領域の
境界半径感度は \(8.94\times10^{-11}\)、発散時刻fitは成長条件を満たさない
ため実行されませんでした。詳細は [STATUS.md](STATUS.md) と
[`outputs/baseline_v5/summary.json`](outputs/baseline_v5/summary.json) に
あります。現行manufactured成果物は
[`outputs/manufactured_v5/diagnostics.json`](outputs/manufactured_v5/diagnostics.json)
です。

固定空間格子の時間収束では、\(\Delta t=0.5,0.25,0.125\) に対する
重み付き相対 \(L^2\) 誤差がそれぞれ
\(8.586076\times10^{-4},2.048943\times10^{-4},4.382275\times10^{-5}\)、
解析解に対する観測次数が \(2.0671,2.2251\)、step-doubling次数が
\(2.0200\) でした。各runでエネルギーは減少し、最大物理渦度は初期値
\(2.0\) を超えず、内側領域の境界感度は
\(2.16\times10^{-8},5.38\times10^{-9},2.03\times10^{-9}\) でした。
この境界感度は補助的な \(R=3,4\) 比較であり、主計算 \(R=5\) の
打切り誤差を直接評価するものではありません。
機械可読な詳細は
[`outputs/time_convergence_v1/summary.json`](outputs/time_convergence_v1/summary.json)
にあります。以前の出力は、途中結果を隠さないためそのまま保存しています。

## リポジトリ構成

| パス | 内容 |
|---|---|
| `SPEC.md` | 数学的対象、変数、解・特異点の定義 |
| `LEAN4_VERIFICATION_POLICY.md` | Lean 4 / mathlib4 による最終検証の必須規約 |
| `docs/equation_audit.md` | 符号・係数・境界・同値性の式別監査(E-01–E-31) |
| `docs/hou_setup_audit.md` | Hou (arXiv:2107.06509) v1/v2 LaTeX 原文の一次資料監査 |
| `docs/nonlinear_solver_design.md` | production 非線形ソルバの設計と受入条件 |
| `docs/formalization_map.md` | 各マイルストーンの Lean 化対応表 |
| `docs/legacy_reuse_review.md` | 旧試作の限定的なread-only監査と非移植判断 |
| `docs/known_obstructions.md` | 既知の非存在・正則性・継続定理 |
| `docs/threat_model.md` | 偽特異点の原因、検出試験、停止規則 |
| `docs/future_search.md` | Type II・動的再スケーリング探索設計 |
| `docs/proof_obligations.md` | 数値候補から反例までの証明義務 |
| `src/ns_certificate_lab/` | 小さなNumPy数値・保存・診断基盤 |
| `tests/` | manufactured、round-trip、故障注入、相互検証 |
| `experiments/` | 安価な監査・非特異対照・早期 Hou 実行 |
| `configs/` | 固定された実験入力 |
| `outputs/` | 機械可読診断・checkpoint・グラフ |
| `archive/` | 統合済みバンドルのパッケージング残骸(provenance 用) |
| `certificates/` | 将来の明示候補証明書用（現在候補なし） |

## 研究上のゲート

既知障害は [docs/known_obstructions.md](docs/known_obstructions.md)、誤検出対策は
[docs/threat_model.md](docs/threat_model.md)、証明までの依存関係は
[docs/proof_obligations.md](docs/proof_obligations.md) を参照してください。
ニューラルネットを将来使う場合も、保持する候補は明示基底係数へ変換し、
ネットワークを使わない独立残差評価に合格させます。

プロジェクトの現在地と未解決事項は [STATUS.md](STATUS.md) にのみ集約し、
数値的確認と数学的証明を混同しません。
