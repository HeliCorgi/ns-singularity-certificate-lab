# NSSingularity — Lean 4 formalization (stage 0)

`LEAN4_VERIFICATION_POLICY.md` の段階 0: Clay 公式問題文
(C. Fefferman, *Existence and Smoothness of the Navier–Stokes Equation*,
Clay Mathematics Institute)に対応する命題を Lean 4 上に**定義のみ**で
固定する。証明・`sorry`・公理は含まない。

- 環境: `lean-toolchain` = `leanprover/lean4:v4.32.1`、
  mathlib4 `v4.32.1`(`lakefile.toml` で git rev 固定)
- ビルド: `lake exe cache get && lake build`

## 自然言語対応表(独立査読用)

| Fefferman 公式記述 | Lean 定義 (`NSSingularity/ClayStatement.lean`) | 忠実性の注記 |
|---|---|---|
| 方程式 (1) 運動量方程式(\(t>0\) で成立) | `SatisfiesMomentum` | 成分別・座標方向の Fréchet 微分で記述 |
| 方程式 (2) \(\mathrm{div}\,u=0\) | `SatisfiesIncompressibility` / `DivergenceFree` | — |
| 条件 (3) 初期値 | `IsGlobalSmoothSolution` 内の `u 0 x = u₀ x` | — |
| 条件 (4) 初期値の急減衰 \(|\partial^\alpha u^0|\le C_{\alpha K}(1+|x|)^{-K}\) | `RapidlyDecaying` | 多重指数 \(\alpha\) の代わりに全微分次数 \(n\) の `iteratedFDeriv` ノルムで量化(多重指数形を優越) |
| 条件 (5) 外力の空間時間減衰 | `AdmissibleForce` | 同上、\((1+|x|+t)^{-K}\) |
| 条件 (6) \(p,u\in C^\infty(\mathbb R^3\times[0,\infty))\) | `SmoothOnHalfTime` / `SmoothScalarOnHalfTime` | \([0,\infty)\times\mathbb R^3\) 上の joint `ContDiffOn` |
| 条件 (7) 有界エネルギー | `BoundedEnergy` | Bochner 積分・volume 測度 |
| (8) 空間周期性 | `SpacePeriodic` | 単位格子 \(\mathbb Z^3\)、各座標方向の平行移動不変 |
| (9) 周期外力の有界性 | `AdmissiblePeriodicForce` | 全微分の一様有界 |
| (10) 周期初期値 | `AdmissiblePeriodicData` | 減衰要求なし |
| **命題 (A)** 全空間大域正則性 | `ClayWholeSpaceRegularity` | \(\forall\nu>0\)、外力 0 |
| **命題 (B)** 周期大域正則性 | `ClayPeriodicRegularity` | 同上 |
| **命題 (C)** 全空間 breakdown | `ClayWholeSpaceBreakdown` | 外力あり許容 |
| **命題 (D)** 周期 breakdown | `ClayPeriodicBreakdown` | — |
| 本リポジトリの優先目標(C より強い) | `UnforcedWholeSpaceBreakdown` | 外力恒等的に 0 |

## 既知の忠実性論点(査読対象)

1. **微分可能性の前提**: `fderiv` は微分不能点で junk 値を返すが、
   すべての使用箇所は同一命題内の滑らかさ仮定で保護される。
   breakdown 命題では「滑らかな解が存在しない」ことを主張するため、
   非滑らかな candidate ペア \((u,p)\) は `IsGlobalSmoothSolution` の
   滑らかさ節で自動的に排除される。
2. **急減衰の量化**: Fefferman の多重指数 \(\alpha\) ごとの評価を、
   `iteratedFDeriv` の作用素ノルム(全次数 \(n\))で置いた。作用素ノルムは
   各方向微分を優越するので、この形は公式条件を**含意する**。逆向きの
   同値性は成り立つ(有限次元なので定数倍の違い)が、段階 0 では
   「この Lean 命題を証明すれば公式命題が従う」方向のみを要求する。
3. **時間微分の扱い**: 運動量方程式の \(\partial_tu\) は各点の
   `deriv`(1 変数微分)で表す。`SmoothOnHalfTime` の joint 滑らかさが
   その存在を保証する。
4. **圧力の規約**: 公式文どおり \(p\) に減衰・正規化を課さない。
5. **周期設定のエネルギー**: (B)/(D) では公式文に従い有界エネルギー条項を
   置かない(周期+滑らかさから自動)。

これらの論点は独立した数学者による対応表査読(policy §「必須の監査」)の
対象である。

## 構成(現状)

```text
formal/
├── lakefile.toml        # mathlib v4.32.1 固定
├── lean-toolchain       # leanprover/lean4:v4.32.1
├── NSSingularity.lean   # ルート import
└── NSSingularity/
    ├── ClayStatement.lean      # 段階 0(定義のみ、sorry なし)
    ├── VelocityRecovery.lean   # 段階 1 / F-3(E-14/E-15、証明あり、sorry なし)
    └── FiniteTime.lean         # 段階 1 / F-2(有限物理時間、証明あり、sorry なし)
```

以降の段階(1: 有限次元恒等式 F-1〜F-4、2: 解析的な橋、3: 数値証明書、
4: 最終定理)は `docs/formalization_map.md` を参照。

## 段階 1 — F-3: 速度回復式と発散ゼロ(`NSSingularity/VelocityRecovery.lean`)

### Lean で証明したこと

`docs/equation_audit.md` の E-14 / E-15 を機械検証した。子午面を素直な
積 `ℝ × ℝ`(`p.1 = r`、`p.2 = z`)でモデル化し、偏微分は
`ClayStatement.lean` と同じ規約(座標方向に評価した Fréchet 微分)で

```text
partialR f p = fderiv ℝ f p (1,0)      partialZ f p = fderiv ℝ f p (0,1)
uR ψ p = -(p.1 * partialZ ψ p)         uZ ψ p = 2 * ψ p + p.1 * partialR ψ p
```

と定義する(E-14)。主定理は

| 定理 | 主張 |
|---|---|
| `divergence_of_recovered_velocity_eq_zero` | `ContDiff ℝ 2 ψ`、`p.1 ≠ 0` のとき `partialR (uR ψ) p + uR ψ p / p.1 + partialZ (uZ ψ) p = 0`(E-15 そのもの) |
| `divergence_of_recovered_velocity_eq_zero'` | `uʳ/r` をその連続延長 `uROverR ψ = -partialZ ψ` に置き換えた形。`r ≠ 0` 不要で軸上 `r = 0` でも成立 |
| `mixed_partial_comm` | `partialR (partialZ ψ) p = partialZ (partialR ψ) p`。**滑らかさが入る唯一の箇所** |
| `partialR_uR` / `partialZ_uZ` | E-15 の第 1・第 3 括弧(積の微分則のみ) |
| `uROverR_eq_div` | `p.1 ≠ 0` で `uROverR ψ p = uR ψ p / p.1` |

`mixed_partial_comm` は mathlib の
`ContDiffAt.isSymmSndFDerivAt : ContDiffAt 𝕜 n f x → minSmoothness 𝕜 2 ≤ n →
IsSymmSndFDerivAt 𝕜 f x` を使う。E-15 の相殺は**混合偏微分の一致だけ**に
依存しており、Lean 化によりその依存が構造的に露出している。

`sorry`・`admit`・新規 `axiom` はない。
`#print axioms divergence_of_recovered_velocity_eq_zero` は
`[propext, Classical.choice, Quot.sound]`(mathlib 標準の古典公理のみ)。

### Lean で証明していないこと(重要な限界)

1. これは**選んだ座標表現に関する恒等式**であって、Navier–Stokes 方程式の
   解についての主張ではない。`ψ` は任意の `C²` スカラー関数でよく、
   E-13(`-𝓛₅ψ₁ = ω₁`)も運動方程式も一切使っていない。
2. `uR` / `uZ` が実際に3次元軸対称ベクトル場の円柱成分であること、すなわち
   E-18 の Cartesian 復元と E-24 の同値性は**未形式化**である。したがって
   本定理から「`ℝ³` 上の発散ゼロベクトル場が得られる」とは言えない。
   `ClayStatement.lean` の `DivergenceFree` との接続も未着手。
3. 軸 `r = 0` 上の扱い: `divergence_of_recovered_velocity_eq_zero'` は
   `uʳ/r` を連続延長 `-∂_zψ₁` に**定義として置き換えた**主張であり、
   極限が実際にその値へ収束することの証明(E-16 の偶奇性を要する
   軸正則性の議論)は含まない。Lean の `/` は全域(`x / 0 = 0`)なので、
   素の商を使う主定理には `p.1 ≠ 0` が必要である。
4. `fderiv` は微分不能点で junk 値を返すが、本ファイルのすべての命題は
   同一命題内の `ContDiff ℝ 2` 仮定で保護されている
   (`ClayStatement.lean` と同じ論点)。

## 段階 1 — F-2: 有限物理時間条件(`NSSingularity/FiniteTime.lean`)

### Lean で証明したこと

動的再スケーリング軌道では物理時刻 `t` と再スケーリング時刻 `s` が
`dt/ds = L(s)^2` で結ばれる。`L : ℝ → ℝ` をスケール関数として

```text
scaleRate L σ           = L σ ^ 2                       -- 再スケーリング率 dt/ds
physicalTime t₀ s₀ L s  = t₀ + ∫ σ in s₀..s, scaleRate L σ   -- 区間積分
blowupTime  t₀ s₀ L     = t₀ + ∫ σ in Set.Ioi s₀, scaleRate L σ
```

と定義する。`physicalTime` を**区間積分**、`blowupTime` を **`Ioi` 上の集合
積分**で書くのは意図的で、前者は微積分学の基本定理
(`intervalIntegral.integral_hasDerivAt_right`)がその形で述べられているため、
後者は広義積分の収束定理
(`MeasureTheory.intervalIntegral_tendsto_integral_Ioi`)がその形だからである。
両者を接続するのがまさにその収束定理であり、解析的内容が入る唯一の箇所である。

| 定理 | 主張 |
|---|---|
| `tendsto_physicalTime` | `IntegrableOn (scaleRate L) (Set.Ici s₀)` のとき `Tendsto (physicalTime t₀ s₀ L) atTop (𝓝 (blowupTime t₀ s₀ L))`(**F-2 本体**) |
| `physicalTime_le_blowupTime` | 同仮定のもと `∀ s, physicalTime t₀ s₀ L s ≤ blowupTime t₀ s₀ L` |
| `physicalTime_lt_blowupTime` | さらに `∀ σ, 0 < L σ` なら `∀ s, physicalTime … s < blowupTime …`(**到達しない**) |
| `exists_finite_blowupTime` | 上記 3 つの梱包形: `∃ T, Tendsto … (𝓝 T) ∧ (∀ s, … ≤ T) ∧ MonotoneOn … (Ici s₀)` |
| `physicalTime_monotoneOn` | 局所仮定 `∀ b, s₀ ≤ b → IntervalIntegrable (scaleRate L) volume s₀ b` のもと `MonotoneOn (physicalTime t₀ s₀ L) (Set.Ici s₀)` |
| `physicalTime_strictMonoOn` | 加えて `∀ σ, 0 < L σ` なら `StrictMonoOn … (Set.Ici s₀)`。**正値性が使われる唯一の実質的箇所** |
| `tendsto_physicalTime_atTop` | 局所可積分だが `¬ IntegrableOn (scaleRate L) (Set.Ioi s₀)` なら `Tendsto (physicalTime t₀ s₀ L) atTop atTop`(**逆向き**) |
| `hasDerivAt_physicalTime` | `Continuous L` のとき `HasDerivAt (physicalTime t₀ s₀ L) (scaleRate L s) s`。定義が確かに `dt/ds = L²` を解いていることの証明 |
| 補助 | `intervalIntegrable_scaleRate`、`le_physicalTime_of_le`、`physicalTime_le_of_le_base`、`physicalTime_le_physicalTime`、`physicalTime_base`、`scaleRate_nonneg`、`scaleRate_pos` |

被積分関数を `L σ ^ 2` と書くことで非負性が `sq_nonneg` から無償で出るため、
**単調性・有限極限・上界 `t(s) ≤ T` のいずれにも `L` の正値性は不要**である。
正値性は狭義単調性と `t(s) < T` にのみ使う。

`tendsto_physicalTime` と `tendsto_physicalTime_atTop` を合わせると、
可積分性は「有限爆発時刻」の十分条件であるだけでなく**分水嶺そのもの**である
ことが Lean 上で確定する。

`sorry`・`admit`・新規 `axiom` はない。
本ファイルの全 15 定理について `#print axioms` は
`[propext, Classical.choice, Quot.sound]`(mathlib 標準の古典公理のみ)。

主に使用した mathlib 補題:
`MeasureTheory.intervalIntegral_tendsto_integral_Ioi`、
`MeasureTheory.integrableOn_Ioi_of_intervalIntegral_norm_bounded`、
`MeasureTheory.setIntegral_pos_iff_support_of_nonneg_ae`、
`intervalIntegral.integral_mono_interval`、
`intervalIntegral.integral_add_adjacent_intervals`、
`intervalIntegral.integral_hasDerivAt_right`。

### Lean で証明していないこと(重要な限界)

1. これは**与えられたスケール関数 `L` についての命題**である。`L` を構成しない。
   任意の可測な `L : ℝ → ℝ` に対して成り立つ主張であり、爆発の存在とは無関係。
2. Navier–Stokes 方程式との接続は**一切ない**。`L` が実際の解の自己相似スケール
   であること、`physicalTime` が実際の解の物理時刻であることは形式化されていない。
   `ClayStatement.lean` の定義とも未接続。
3. **仮定 `IntegrableOn (fun σ => L σ ^ 2) (Set.Ici s₀)` こそが本質**であり、
   本ファイルはそれを一切検証しない。将来の区間証明書が厳密な数値上界として
   供給すべき対象がまさにこの仮定である
   (`docs/formalization_map.md` の M5 / F-4 参照)。
   なお `scaleRate` は可測性を仮定していないが、`IntegrableOn` が
   `AEStronglyMeasurable` を含むので可測性はこの仮定に内包されている。
4. `t(s) < T` は「物理時刻が `T` に達しない」ことを言うだけで、`T` において
   解が特異になること(ノルム発散)は含まない。それは別の blow-up criterion の
   形式化(未着手)を要する。
5. `hasDerivAt_physicalTime` は `Continuous L` を仮定する。主定理群は連続性を
   使わないので、この補題は定義の妥当性確認のためだけに存在する。
