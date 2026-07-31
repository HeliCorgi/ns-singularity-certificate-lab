# Final target — Clay 命題までの単一依存グラフ

作成: 2026-07-29(branch `fable5-mainline`)
最終更新: 2026-07-29 第 6 便(Gate 6: 中振幅校正と振幅継続。2 基準不合格。非線形 tail 伝播、明示的初期値族、動的領域拡大、区間 snapshot 証明書、F-17/F-18/F-19)
根拠: `START_NEW_SESSION_NAVIER_STOKES.md` §6 Step 2。

本書は、Clay 公式命題 (A)〜(D)、Track U(外力なし全空間反例)の最終定理、
Track F(滑らかな外力を使う (C)/(D) 反例)の最終定理、および全証明義務を
**一つの依存グラフ**にまとめ、各義務を次の 5 状態のいずれかに分類する。

| 記号 | 意味 |
|---|---|
| **M** | 数学的に閉じた(紙上の証明が完了。Lean 未形式化でもよい) |
| **L** | Lean 4 で閉じた(`lake build` 成功、`sorry`/`admit`/新規 axiom なし、`#print axioms` 記録済み) |
| **I** | 区間演算(validated numerics)が必要 |
| **N** | 数値的観測のみ(浮動小数点。証明ではない) |
| **O** | 未着手 |

**重要**: 状態 **N** はどれも証明ではない。**M** と **L** だけが証明であり、
現時点で Clay 命題そのものに到達している項目は一つもない。

---

## 1. 最終命題(Lean 定義は `formal/NSSingularity/ClayStatement.lean`)

| ID | 命題 | Lean 定義名 | 状態 |
|---|---|---|---|
| CLAY-A | `ℝ³` 全空間・外力なしの大域正則性 | `ClayWholeSpaceRegularity` | **O**(定義のみ) |
| CLAY-B | 周期領域・外力なしの大域正則性 | `ClayPeriodicRegularity` | **O**(定義のみ) |
| CLAY-C | `ℝ³` 全空間・滑らかな外力ありの破綻 | `ClayWholeSpaceBreakdown` | **O**(定義のみ) |
| CLAY-D | 周期領域・滑らかな外力ありの破綻 | `ClayPeriodicBreakdown` | **O**(定義のみ) |
| TARGET-U | **外力ゼロ**の全空間破綻((C) より真に強い) | `UnforcedWholeSpaceBreakdown` | **O**(定義のみ) |

定義は Fefferman の公式問題文に対応づけて監査済み(`formal/README.md`)。
**定義を書いたことは証明ではない。**

## 2. Track U — 外力なし全空間反例(TARGET-U → CLAY-C)

```text
TARGET-U
 ├─ U-1  明示的な滑らか・発散ゼロ・急減衰初期値の許容性        [PO-02]
 ├─ U-2  縮約系(軸対称+swirl)と 3D NS の同値性               [PO-01]
 ├─ U-3  候補 profile / 再スケーリング軌道の厳密存在           [PO-04, PO-13]
 │    ├─ U-3a 離散化誤差の厳密上界                             [PO-05]
 │    ├─ U-3b 領域打切り誤差の厳密上界(全空間 tail)          [PO-06]
 │    ├─ U-3c スペクトル尾部の厳密上界                         [PO-07]
 │    └─ U-3d Newton–Kantorovich / radii polynomial            [F-4]
 ├─ U-4  候補近傍の非線形安定性                                [PO-08]
 ├─ U-5  滑らかな初期値から候補軌道への到達                    [PO-03, PO-09]
 ├─ U-6  有限物理時刻 T < ∞                                    [F-2 ✔L]
 ├─ U-7  物理ノルムの発散(座標由来でない)                    [PO-11, PO-12]
 └─ U-8  最終論理接続                                          [PO-15]
```

| ID | 状態 | 現況 |
|---|---|---|
| U-1 | **O** | 候補が存在しないため未着手 |
| U-2 | **M**(部分) | `docs/equation_audit.md` E-11〜E-31 で導出照合済み。Lean は F-3 のみ **L** |
| U-3 | **O** | 候補自体が未発見。線形 Gate 4(値)と Gate 5(微分 tail・速度回復・**小振幅**非線形 run)は合格。強振幅の全空間発展は未着手 |
| U-3a/b/c | **I**(着手) | 第 6 便で**単一 snapshot の区間証明書**を実装(厳密有理数演算、独立 checker つき)。離散量の包含のみで、離散化誤差そのものは未包含 |
| U-3d | **O** | `F-4` として形式化予定、未着手 |
| U-4, U-5 | **O** | — |
| U-6 | **L** | `NSSingularity.tendsto_physicalTime` ほか。ただし可積分性仮定は **I** で供給する義務が残る |
| U-7 | **O** | blow-up criterion の形式化は未設計 |
| U-8 | **O** | — |

**Track U の現在の律速**(第 4 便で更新): **線形** Gate 4(非周期 `z`、
自由空間楕円経路、`R_max`/`Z_max` 独立拡大、低波数 stress test、
a posteriori tail bound、独立 Cartesian 検査)は
`src/ns_certificate_lab/whole_space_gate.py` として実装され、前登録受入検査
20 件すべてに合格した(`outputs/whole_space_gate4_v1`)。したがって
**非線形時間発展への結合は許可された**。ただし線形ゲートの合格は
非線形発展について何も主張しない。次の律速は、一様全空間格子で
Hou 型 front を解像できるかという解像度設計そのものである。

### Track U で既に閉じた除外定理(候補の構成ではない)

| ID | 内容 | 状態 |
|---|---|---|
| U-X1 | 臨界 `L³` 障害: 一様 `L³` 有界な一尺度再スケーリング軌道は有限時間爆発を与えない(ESS 端点定理を**引用**) | **M**(引用依存) |
| U-X2 | 異方的版 `‖u‖³_{L³} = A³L_r²L_z‖U‖³_{L³}`;標準等方放物型では積が恒等的に 1 なので除外 | **M** |
| U-X3 | 有限円柱壁補正の閉形式(ゼロモードの `R⁻²` 代数尾部) | **M** |
| U-X4 | 周期 `z` のゼロ軸モードは自由場を `2R/L` 倍だけ過大評価する(厳密な比)。非周期 `z` の Dirichlet 作用素にはゼロモードが存在しない | **M** |
| U-X5 | 単極子境界条件の連続レベル打切り誤差 `≤ 3a‖ω‖_{L¹(dV₅)}/(8π²(R−a)⁴)`(5 次元 Laplacian の最大値原理) | **M** |
| U-X6 | **微分 tail 上界** `|D^k(ψ−ψ_multipole)| ≤ A_{k+1}I₁/d^{4+k}`(単極子)、`≤ ½A_{k+2}I₂/d^{5+k}`(双極子)。`A_m` は `D^mG₅` の斉次性定数。Green 核の解析微分と Taylor 剰余のみ、最大値原理不使用 | **M**、`A_1,A_2` は Lean F-14 に接続 |
| U-X7 | box 内部への微分伝播 `|D^k e| ≤ (5k/ρ)^k ε₀`(調和関数の内部楕円評価) | **M**(連続レベル) |
| U-X8 | **非線形 tail 伝播**: `ε₀,ε₁,ε₂ → |δu^r|,|δu^z| → 移流項誤差 → 短時間 Grönwall`。全定数明示、積差恒等式のみ使用 | **M**、合成層は Lean F-17/F-18/F-19 で **L** |
| U-X9 | コンパクト台の初期値では `ω₁` の 5 次元単極子と軸方向四重極が**対称性で消える**ので、zero と monopole のトレースは一致し、dipole と quadrupole も一致する | **M**(実測 `5.2e-11` の比で確認) |

## 3. Track F — 滑らかな外力を逆設計する反例(→ CLAY-C / CLAY-D)

```text
CLAY-C / CLAY-D  (via Track F)
 ├─ F-α  特異 ansatz u,p の設計(t<T で滑らか、T で延長不能)
 │    └─ F-α1  有限モード / 有限次元 ansatz 族        ✘ 除外済み(本節)
 ├─ F-β  残差 f = ∂_t u + (u·∇)u − νΔu + ∇p が全時間滑らか
 ├─ F-γ  初期値が滑らか・発散ゼロ(+ (C) では急減衰)
 ├─ F-δ  局所一意性により同じデータの別の大域滑らか解を排除   [PO-03]
 └─ F-ε  Clay (C)/(D) への最終接続                            [PO-15]
```

| ID | 状態 | 現況 |
|---|---|---|
| F-α1(**固定有限次元・固定帯域**のみ) | **M + L** | **除外定理として閉じた**。`docs/research_notes/track_f_finite_mode_nogo.md` Theorem 1 と Corollary 1–3。Lean 側は `F-6`(a priori 上界)、`F-12`(三線型相殺の Fourier 版)、`F-13`(有限 Fourier 空間のノルム同値)、`F-7a`(端点到達)、`F-7b`(自励系の局所延長)がすべて **L**。時間依存外力での延長 `F-7c` も第 7 便で **L**(`TimeDependentGalerkin.lean`)|
| F-α(残り = 帯域幅発散族) | **M(必要条件のみ)** | `docs/research_notes/track_f_shell_constraints.md`。**現在の shell ansatz と非退化仮定の下で**実現可能領域は `0<γ<1`、`max(0,2γ−1) ≤ σ < γ`、`β>0` の三角形に限られ、同じ仮定の下で `γ ≥ 1` は排除される。候補は一つも構成していない |
| F-β | **O** | — |
| F-γ | **O** | — |
| F-δ | **O** | 局所適切性(Kato 型)は mathlib になく、自前形式化か証明書側での回避が必要 |
| F-ε | **O** | — |

### Track F で確立した必要条件(任意の (C)/(D) 反例が満たすべきもの)

| ID | 内容 | 状態 |
|---|---|---|
| F-N1 | エネルギー有界: `sup_{t<T}‖u(t)‖_{L²} ≤ ‖u(0)‖_{L²} + ∫₀ᵀ‖f‖_{L²}`。**外力を許しても Clay 条件 (7) は自動的に満たされる** | **M** |
| F-N2 | 総散逸有限: `ν∫₀ᵀ‖∇u‖²_{L²} dt < ∞`。破綻は `L²_tH¹_x` にも見えない | **M** |
| F-N3 | 破綻には Ladyzhenskaya–Prodi–Serrin 臨界ノルムの発散が必要 | **引用**(未証明・未形式化) |
| F-N4 | 滑らかな外力の Fourier 係数は任意の多項式より速く減衰するので、**高シェルへの直接注入 `⟨u_j,f_j⟩` は不可視**。低周波制御による間接駆動は**排除されていない** | **M**(直接注入の排除のみ) |
| F-N5 | 有限 cascade 模型で、低周波のみの外力が直接注入ゼロのまま高シェル振幅を 26 桁駆動することを確認。粘性が競合相手 | **N**(模型。PDE ではない) |

**帰結**: 「外力を使えば反例が作りやすい」という当初の想定は F-N1/F-N2 により
否定された。`L²` エネルギーと時間積分された `H¹` 散逸は、滑らかな外力の下で有限時間内に自動的に有界であり、特異点の直接的な発散指標にはできない。
外力を許しても、破綻の判定は臨界ノルムでしか行えない。

## 4. Lean 4 形式化の登録簿(採番の確定)

外部セッション由来のノートが `F-4`〜`F-7` を独自に提案しており、本リポジトリの
既存採番と衝突していた(`docs/research_notes/README.md` §3 で記録済み)。
**本書をもって次のとおり確定する。** 以後この表が唯一の権威である。

| ID | 内容 | 状態 | ファイル |
|---|---|---|---|
| F-1 | 再スケーリングの代数的恒等式(E-21b) | **O** | — |
| F-2 | 有限物理時間条件 | **L** | `formal/NSSingularity/FiniteTime.lean` |
| F-3 | 速度回復式と発散ゼロ(E-14/E-15) | **L** | `formal/NSSingularity/VelocityRecovery.lean` |
| F-4 | 候補証明書の有限次元不等式(radii polynomial) | **O** | — |
| F-5 | Clay 命題の定義と依存 DAG | **O**(定義のみ。証明はない) | `formal/NSSingularity/ClayStatement.lean` |
| **F-6** | **Galerkin エネルギー上界(Track F 有限モード除外の中核)** | **L** | `formal/NSSingularity/GalerkinNoBlowup.lean` |
| **F-7a** | **有界な軌道が有限時刻端点で極限を持つ** | **L** | `formal/NSSingularity/FiniteModeNoGo.lean` |
| **F-7b** | **自励 Galerkin 系の端点からの局所延長(Picard–Lindelöf)** | **L** | 同上 |
| **F-7c** | **時間依存外力での局所延長** | **L** | `TimeDependentGalerkin.lean`。mathlib の `IsPicardLindelof` が時間依存なので自励化は不要だった。放棄した自励化経路の還元は `GreenAndCascade.lean` に記録として残置 |
| **F-14** | **5 次元 Green 動径プロファイルの調和性** | **L** | `formal/NSSingularity/GreenAndCascade.lean` |
| **F-15** | **Newton の flux 恒等式 `R⁴ψ' = −m`** | **L** | 同上 |
| **F-16** | **shell 指数領域(仮定を構造体で明示)** | **L** | 同上 |
| **F-17** | **ポテンシャル誤差 → 速度誤差** | **L** | `formal/NSSingularity/CertificateLayer.lean` |
| **F-18** | **積差恒等式と移流項誤差** | **L** | 同上 |
| **F-19** | **短時間 Grönwall(証明書向け形）** | **L** | 同上 |
| F-8 | 等方的 `L³` スケーリング恒等式(旧ノート `F-4`) | **O** | — |
| F-9 | 異方的 `L³` スケーリング恒等式(旧ノート `F-5`) | **O** | — |
| F-10 | 有界スケーリング積 + 有界 profile ⇒ 物理 `L³` 有界(旧ノート `F-6`) | **O** | — |
| F-11 | 端点正則性定理から候補クラス除外への橋(旧ノート `F-7`) | **O** | — |
| **F-12** | **三線型相殺 `⟨u,(u·∇)u⟩=0` の Fourier 表現版** | **L** | `formal/NSSingularity/FiniteModeNoGo.lean` |
| **F-13** | **有限 Fourier 空間のノルム同値定数(重み付き和と Cauchy–Schwarz)** | **L** | 同上 |

`F-11` は Escauriaza–Seregin–Šverák の忠実な形式化または明示的に監査された
定理インタフェースを要する。**未証明の project 固有 axiom として挿入しては
ならない**(`LEAN4_VERIFICATION_POLICY.md`)。

### 4.1 F-7c — **閉じた**(直接経路)

「F-7 は既存定理の単純な適用」という当初の見立ては誤りだったが、「mathlib の
局所存在定理は自励系専用」という第 4 便の判定**も誤りだった**。実測:

- **F-7a** は `intervalIntegral.integral_eq_sub_of_hasDerivAt` +
  `intervalIntegral.continuousOn_primitive_interval` +
  `Integrable.of_bound` で閉じた。
- **F-7b** は `ContDiffAt.exists_forall_mem_closedBall_exists_eq_forall_mem_Ioo_hasDerivAt₀`
  が直接使えた。これは確かに `f : E → E` の自励系専用である。
- **F-7c(閉じた、`formal/NSSingularity/TimeDependentGalerkin.lean`)**:
  自励系専用なのは*この定理*であって API 全体ではない。pin されている
  mathlib の `IsPicardLindelof` は最初から `f : ℝ → E → E` の**時間依存**場に
  対して述べられており、`IsPicardLindelof.exists_eq_forall_mem_Icc_hasDerivWithinAt₀`
  が解を直接与える。したがって自励化は不要だった。

  **2 経路の比較(セッション指示による)**:

  | | 自励化経路 | 直接経路(採用) |
  |---|---|---|
  | `IsPicardLindelof` を作る空間 | `E × ℝ` | `E` |
  | 追加の instance 義務 | `CompletePace (E × ℝ)`、Prod ノルムでの Lipschitz 評価 | なし |
  | 追加で必要な定理 | `galerkin_solution_of_autonomised`(70 行) | なし |
  | 第 2 成分の恒等性の議論 | 必要 | 不要 |

  直接経路が仮定・行数ともに厳密に少ないため、こちらを採用した。

  実質的な内容は `B x x` の**局所** Lipschitz 評価 1 点である:
  `B x x - B y y = B x (x-y) + B (x-y) y` なので Lipschitz 定数は
  `2‖B‖(‖x₀‖+a)` となり、球半径 `a` が必ず入る。これを大域定数にすることは
  できない。`galerkin_isPicardLindelof` がこの 4 条件を与え、
  `galerkin_local_solution` が両側局所解を出す。`#print axioms` は 3 定理とも
  `[propext, Classical.choice, Quot.sound]` のみ。

  `galerkin_solution_of_autonomised`(`GreenAndCascade.lean`)は放棄した経路の
  記録として残すが、依存するものは無くなった。

  これで固定有限モード除外は Lean だけで完結する(残る仮定は `advectionForm`
  と実 PDE の同一視のみ)。

## 4.5 Track P — 周期 T³ の厳密 a posteriori レーン(第 9 便で新設)

Clay 公式命題 (B)/(D) は周期領域の主張であり、全空間で証明を阻んでいた
2 つの障害が周期側では消える:

1. **スペクトルギャップ**: 平均ゼロ場では `|k| ≥ 1` なので control 不等式の
   減衰項 `−νR` が存在する(全空間では `−Δ` のスペクトルが 0 に達し消滅、
   `a_posteriori_frameworks.md` 参照)。
2. **HS-5 の周期版は構成的に閉じる**: Fourier–Galerkin 軌道は厳密な三角多項式
   なので、その**連続** Navier–Stokes 残差 `(I−P_G)P(u·∇u)` は有限三角多項式で
   あり、有理数演算で厳密に計算できる。空間補間も離散→連続の復元も不要。
   全空間レーンの HS-5 ギャップはそのまま残る。両レーンを混同してはならない。

**除外領域の区別(誤読防止)**: 有限モード no-go(Track F)が除外するのは
**全時刻で固定有限帯域に留まる軌道**である。**有限帯域の初期値**はこの仮定を
満たさない — 真の解は即座に帯域を離れる(Galerkin tail が非ゼロ)— ので、
初期値については何も除外されていない。Track P はこの未除外領域で作業する:
有限帯域初期値、無限帯域の真の解、両者の距離を証明書で制御する。
Lean 側は `TrackPFourier.lean` の `FixedBandTrajectory` / `FiniteBandDatum`。

**第 9 便の成果**: 有理 Fourier 初期値 3 族(P1 helical triad、P2 連結 2 triad、
P3 対称性破れ・圧力駆動)について、`H⁴` control 不等式
`d⁺R/dt ≤ (−ν+9(K₁+K₂))R + 135Ȧ R² + ‖e‖_Ḣ⁴`(全定数自前導出、
`docs/research_notes/track_p_periodic.md`)による**12/12 スラブ証明書**が成立:
古典外部定理 EXT-P1/P2/P3(忠実記録、Lean 公理化なし)の下で、真の周期強解が
スラブ全体に存在し `‖u−u_a‖_Ḣ⁴ ≤ R(t)`(相対 5e-4〜9e-4)。
**これは軌道近傍の正則性の証明であり、特異点証明の反対物である。**

**第 10 便の成果(スラブ連結)**: 単発スラブをスカラー `H⁴` 誤差半径で連結する
レーンを実装(`torus_chain.py`、`docs/research_notes/track_p_chain.md`)。
各スラブは**厳密有理・厳密発散ゼロの再中心化点**から開始し、区間 box は
スラブ境界を越えて伝播しない(wrapping の入る場所が構造的にない — 前登録の
Lohner/QR 導入条件は実測で不発火)。δ 漸化式
`δ_{n+1} = R_n(t_{n+1}) + transfer` と連結の有限不等式骨格は Lean
(`TrackPChain.lean`)で証明。3 族 × 4 粘性 + 長尺 1 本の連結が
`outputs/track_p_chain_v1/` にあり、全て独立 checker が全リンクを再計算して
検証。到達地平は Riccati 天井 `T* ≈ (1/a)log(a²/(bε))` に律速され、停止分類は
一貫して `control_linear_coefficient`(粗い自前定数 `9(K₁+K₂)` が縛り —
解の性質ではない)。**証明区間の終了は特異点の主張ではない**(checker 強制)。
EXT-P1 については完全な紙上証明(未監査)と依存表が
`docs/research_notes/ext_dependencies.md` に整備され、有限次元核
(二次 ODE の Picard–Lindelöf)は `GalerkinPicard.lean` で無条件に証明済み。
EXT-P1/P2/P3 は引き続き全 payload で `proved: false`。

## 5. 証明義務の状態一覧(`docs/proof_obligations.md` の PO と対応)

| PO | 内容 | 状態 | 備考 |
|---|---|---|---|
| PO-01 | 元の 3D PDE との同値性 | **M**(部分) | 紙上照合済み。Lean は F-3 のみ |
| PO-02 | 滑らかな有限エネルギー初期データ | **O** | 候補未定 |
| PO-03 | 局所解の存在・一意性 | **O** | mathlib に NS の定義自体がない |
| PO-04 | 候補軌道 / profile の存在 | **O** | Gate 4 未実装が律速 |
| PO-05 | 離散化誤差の制御 | **I** | 区間証明書は**離散量の包含**であって離散化誤差の包含ではない。依然未着手 |
| PO-06 | 領域打切り誤差の制御 | **I** | 壁補正の閉形式は **M**。第 4 便で連続レベルの a posteriori 上界(U-X5)が **M** になったが、離散化誤差の上界は依然 **N** |
| PO-07 | スペクトル尾部の評価 | **I** | — |
| PO-08 | 非線形安定性 | **O** | — |
| PO-09 | 初期データから候補軌道への進入 | **O** | — |
| PO-10 | 有限物理時刻 | **L**(条件付き) | F-2。可積分性の供給は **I** |
| PO-11 | 適切なノルムの発散 | **O** | F-N1/F-N2 により `L²`/`L²_tH¹` は判定量として使えないと判明 |
| PO-12 | 座標変換由来の見かけの発散でないこと | **O** | — |
| PO-13 | 区間演算による検証 | **I**(着手) | 第 6 便で単一 snapshot の証明書を厳密有理数演算で生成し、**独立 checker** で検査(10 検査合格)。時間発展の証明書は未着手 |
| PO-14 | 独立実装による再現 | **N**(部分) | Poisson は 4 経路(第 4 便で非 FFT の軸方向 DST 経路が加わった)、時間発展は 1 実装のみ |
| PO-15 | 形式証明・最終定理 | **O** | F-2/F-3/F-6 が部品 |

## 6. 数値的観測のみの項目(証明経路には乗らない)

- 早期 Hou 窓の増幅・`L³`・幅・shell 数・Type-II fit・項別釣合い(すべて **N**)。
- Poisson ゲート、積分器相互比較、von Neumann 監査、core-width、blind 外挿(**N**)。
- 壁依存性・透過境界・低波数壁ゲート(**N**、ただし壁補正の閉形式は **M**)。
- 本セッションの `stream_apriori_bound` クロスチェック(**N**。証明は Lean 側)。

## 7. 現時点で Clay 命題までに残る距離(率直な要約、第 4 便更新)

1. **候補が一つも存在しない。** Track U は線形 Gate 4 を通過したが
   非線形全空間発展は未着手。Track F は固定有限次元・固定帯域族が除外され、
   残るのは帯域幅発散族(必要条件のみ導出済み)。
2. 区間演算の器(**I** の 5 項目)は設計のみで未実装。第 4 便の
   a posteriori tail bound は連続レベルの解析的上界であって区間演算ではない。
3. Lean で閉じているのは F-2 / F-3 / F-6 / F-7a / F-7b / F-12 / F-13 /
   F-14 / F-15 / F-16 と F-7c 還元、および Clay 命題の定義のみ。
   固定有限帯域の鎖(F-6→F-7a→`breakdown_time_set_empty`)だけが
   内部で接続されており、`ClayStatement.lean` へは依然未接続。
4. したがって **(A)〜(D) のいずれについても、証明の骨格すら存在しない**。

## 8. 本書から導かれる次の最小の一手(第 4 便更新)

1. **弱非線形領域を離れること**(Track U、最優先): Gate 6 の振幅継続は
   `max|ω₁| ∝ A²` を相対 `5e-5` で満たし、**二次応答領域を一度も離れて
   いない**。したがって候補の順位付けは非線形挙動を見ていない。
   必要なのは振幅ではなく**時間**(現在 `max|ω₁|·T ≈ 7e-4 ≪ 1`)であり、
   長時間積分の安定性と計算量が次の設計判断である。
   境界条件の検証は別問題で、Gate 6 の測定によれば境界差は離散化誤差の
   `8e-3` 倍しかない。すなわち**外側境界は律速ではない**。
   既存の `core_width.fit_precondition`(front ≥ 7 点)を全空間 box で
   評価し、適応 mesh か半周期 sine 実装かの設計判断を行う。
   **この判断の前に非線形全空間実行を開始しない。**
2. ~~**F-7c の形式化**~~ — 第 7 便で完了(§4.1)。
3. **F-1 の形式化**: 再スケーリング恒等式。F-8/F-9 の前提でもある。
4. Track F の帯域幅発散族については、`Π_j` の**鋭い**上界を導く。
   現在の粗い評価は `(γ,σ)` に追加制約を与えない。
