# Formalization map — Lean 4 検証への対応表

`LEAN4_VERIFICATION_POLICY.md` の要求に基づき、各マイルストーンについて
「最終的に Lean 化する命題 / その仮定 / 現在は数値的にしか確認できていない部分 /
区間証明書へ変換すべき部分 / mathlib で利用可能な基盤 / 新たに形式化が必要な
解析定理 / 最終 Clay 命題までの依存関係」を記録する。

状態語彙は `FABLE5_HANDOFF.md` §12 の 8 段階
(derivation checked / implementation tested / numerical observation / candidate /
validated numerical object / connected physical solution /
finite-time breakdown theorem / Clay-complete result)を用いる。

最終更新: 2026-07-29(branch `fable5-mainline`、第 6 便)

**Lean 識別子 `F-1`〜`F-11` の確定登録簿は `docs/final_target.md` §4 にある。**
外部ノート由来の採番衝突はそこで解消済みで、本書の節見出しはその登録簿に従う。

---

## M0 — Clay 命題の固定(段階0、未着手)

- **Lean 化する命題:** Clay 公式問題文の (A)(B)(C)(D) の 4 命題、および本リポジトリの
  優先目標である「外力なし・滑らかな発散ゼロ有限エネルギー初期値からの
  \(\mathbb R^3\) 有限時間爆発」((C) より強い形)。
- **仮定:** Fefferman の公式問題文(claymath.org PDF)の定義
  (滑らかさ、急減衰、有界エネルギー、外力条件)。
- **数値のみの部分:** なし(これは仕様の固定であり計算を含まない)。
- **区間証明書:** 不要。
- **mathlib 基盤:** `Mathlib.Analysis.InnerProductSpace.EuclideanDist`、
  `Mathlib.MeasureTheory`(\(L^2\) エネルギー)、`ContDiff ℝ ⊤`(滑らかさ)、
  Schwartz 空間 `Mathlib.Analysis.Distribution.SchwartzSpace`(急減衰)。
- **新規形式化が必要:** Navier–Stokes 方程式そのものの述語
  (`NavierStokes.IsSolution`)、発散ゼロ、圧力の存在的取り扱い。
- **Clay までの依存:** 最終定理 `FinalTheorem.lean` はこの命題を結論とする。
  すべての後続マイルストーンがここへ接続される。

## M1 — 軸対称+swirl 成分系の導出(derivation checked)

- **Lean 化する命題:** 3 次元非圧縮 NS(外力なし)の滑らかな軸対称+swirl 解と、
  変換変数 \(u_1=u^\theta/r,\ \omega_1=\omega^\theta/r,\ \psi_1=\psi^\theta/r\) が満たす
  閉じた系
  \[
  \partial_tu_1+u^r\partial_ru_1+u^z\partial_zu_1=2u_1\partial_z\psi_1+\nu\mathcal L_5u_1,\quad
  \partial_t\omega_1+u^r\partial_r\omega_1+u^z\partial_z\omega_1=\partial_z(u_1^2)+\nu\mathcal L_5\omega_1,\quad
  -\mathcal L_5\psi_1=\omega_1
  \]
  (\(\mathcal L_5=\partial_{rr}+\tfrac3r\partial_r+\partial_{zz}\)、
  \(u^r=-r\partial_z\psi_1\)、\(u^z=2\psi_1+r\partial_r\psi_1\))との同値性。
- **仮定:** 解の滑らかさ、軸対称性、軸上の正則性(偶奇性)、円柱座標変換の定義。
- **数値のみの部分:** なし(代数的導出は `docs/equation_audit.md` で照合済み)。
  ただし Lean 上では未検証。
- **区間証明書:** 不要(恒等式)。
- **mathlib 基盤:** `fderiv`/`deriv` の合成則、円柱座標は新規定義が必要
  (mathlib に極座標の一部はあるが 3D 円柱座標の体系的理論はない)。
- **新規形式化が必要:** 円柱⇔Cartesian 成分変換、\(u^r=-r\partial_z\psi_1,\
  u^z=2\psi_1+r\partial_r\psi_1\) から物理発散
  \(\partial_ru^r+u^r/r+\partial_zu^z=0\) が従う恒等式(**F-3 として形式化済み**、
  `formal/NSSingularity/VelocityRecovery.lean`)、
  \(\mathcal L_5\) と軸対称 Laplacian の関係(未着手)。
- **Clay までの依存:** M1 は「縮約系で構成した解が元の 3D 方程式を満たす」ことの
  橋であり、最終定理の再構成段階(段階4の項目 4)が直接依存する。

## M2 — 数値プリミティブと独立 Cartesian 検証(implementation tested)

- **Lean 化する命題:** なし(有限格子の浮動小数点診断は最終証明経路に乗らない)。
- **役割:** 候補探索の信頼性向上のみ。Lean 経路では、これらの数値検査は
  段階3の区間証明書に置き換えられる。
- **区間証明書:** 将来、候補の残差・作用素ノルムを外向き丸め区間で再計算する
  独立実装が必要(現行 binary64 実装は流用しない)。
- **Clay までの依存:** なし(発見支援のみ)。

## M3 — 有限円柱 Poisson ゲート(implementation tested、2 独立実装)

- **Lean 化する命題:** 将来の候補検証で用いる楕円ソルバの誤差評価そのものは
  Lean 化しない。Lean 化するのは「与えられた有限データ(スペクトル係数等)が
  \(-\mathcal L_5\psi_1=\omega_1\) を残差 \(\le\varepsilon\) で満たす」という
  区間証明書検査(段階3の `IntervalChecker.lean` / `OperatorBounds.lean`)。
- **仮定:** 軸極限 \(\lim_{r\to0}\mathcal L_5\psi_1 = 4\partial_{rr}\psi_1(0)\)
  (滑らかな偶関数に対して)。この極限自体は段階1で形式化可能な小さな解析命題。
- **数値のみの部分:** 観測 2 次収束、条件数推定、cross-stencil defect。
  いずれも binary64 であり厳密性なし。
- **区間証明書へ変換すべき部分:** 楕円逆作用素のノルム上界
  (validated inverse)。現行の 2 実装(`poisson.py`,
  `finite_cylinder_poisson.py`)は発見用であり、証明書検査には第三の
  区間演算実装が必要。
- **mathlib 基盤:** 有理数・区間演算の基礎は `Mathlib.Data.Rat`、
  `Mathlib.Analysis.SpecialFunctions`。区間演算パッケージは実質新規。
- **Clay までの依存:** 候補存在証明(段階3→4)の楕円成分。

## M4 — 非線形有限円柱ソルバと Hou 再現(計画中 → numerical observation)

- **Lean 化する命題:** なし(数値再現は証拠であって定理ではない)。
- **役割:** 候補発見の前提。Hou 機構の壁依存性判定が探索方針を決める。
- **数値のみの部分:** すべて。再現成功でも「numerical observation」止まり。
- **区間証明書:** この段階では作らない。候補が確定した後、動的再スケーリング
  軌道に対して段階3形式の証明書を設計する。
- **Clay までの依存:** 直接依存はないが、候補の物理的出所として
  connected-solution 段階(段階4の項目 1–2)の探索を方向付ける。

## M5 — 動的再スケーリング軌道の候補(未着手)

- **Lean 化する命題(将来):** 「明示的有限データ \(\hat U\) の近傍に、再スケーリング
  方程式の厳密解 \(U\) が存在する」(Newton–Kantorovich / radii polynomial 不等式)。
- **仮定:** 再スケーリング方程式が物理 PDE から正しく導出されていること
  (M1 型の同値性命題に還元)。ゲージ条件の非退化性。
- **区間証明書:** 残差上界、線形化逆作用素上界、非線形 Lipschitz 定数、
  収縮不等式の 4 点。`CertificateFormat.lean` の設計対象。
- **新規形式化が必要:** 関数空間(重み付き Sobolev / 解析的 Banach 空間)、
  スペクトル尾部評価、Newton–Kantorovich 定理(mathlib に一般形はない)。
- **Clay までの依存:** 段階4 最終定理の項目 3(軌道の厳密存在)。

## 最終 Clay 命題までの依存関係(全体)

```text
M0 (Clay命題固定)
 └─ FinalTheorem ← 段階4
     ├─ 初期値の許容性     ← M1 の座標恒等式 + 明示的初期値の検査
     ├─ 軌道の厳密存在     ← M5 の区間証明書 + Newton–Kantorovich
     ├─ 物理解への再構成   ← M1 の同値性定理
     ├─ 有限物理時間       ← F-2(形式化済み)+ L(s) 可積分性(区間証明書)
     ├─ t<T での滑らかさ   ← connected-solution 定理(未設計)
     └─ ノルム発散         ← blow-up criterion(新規形式化)
```

## 段階 0–1 の具体的形式化対象(2026-07-28 追加)

数値探索を止めない範囲で先行形式化する小対象。環境は
`formal/`(Lean 4 v4.32.1 + mathlib4 v4.32.1、`lean-toolchain` で固定)。

### F-1 再スケーリングの代数的恒等式(E-21b)

- **命題:** \(u^{(\lambda)}(x,t)=\lambda u(\lambda x,\lambda^2t)\) が E-01 を
  満たすことと \(u\) が満たすことの同値、および
  \(u_1^{(\lambda)}=\lambda^2u_1(\lambda r,\lambda z,\lambda^2t)\) 等の
  変換則(E-21b)。
- **仮定:** 場の微分可能性のみ。純代数+連鎖律。
- **mathlib 基盤:** `fderiv` の合成則、`deriv_comp`。新規理論不要。
- **数値のみの部分:** なし。

### F-2 有限物理時間条件 — **形式化済み**(2026-07-28)

- **状態:** Lean 4 で証明完了。`formal/NSSingularity/FiniteTime.lean`
  (`lake build` 成功、`sorry`・`admit`・新規 `axiom` なし)。
- **命題:** \(L:\mathbb R\to\mathbb R\) に対し
  \(t(s)=t_0+\int_{s_0}^sL(\sigma)^2d\sigma\) は、\(L^2\) が \([s_0,\infty)\)
  上可積分なら \(s\to\infty\) で有限極限
  \(T=t_0+\int_{(s_0,\infty)}L^2\) を持ち、\([s_0,\infty)\) 上単調で
  \(t(s)\le T\)(さらに \(L>0\) なら \(t(s)<T\))。逆に \(L^2\) が
  \((s_0,\infty)\) 上可積分でなければ \(t(s)\to+\infty\)。
  すなわち可積分性は有限爆発時刻の**分水嶺そのもの**である。
- **Lean 定理名:**
  - `NSSingularity.tendsto_physicalTime`
    — `(hint : IntegrableOn (scaleRate L) (Set.Ici s₀)) (t₀ : ℝ) :
      Tendsto (physicalTime t₀ s₀ L) atTop (𝓝 (blowupTime t₀ s₀ L))`(**F-2 本体**)
  - `NSSingularity.physicalTime_le_blowupTime`
    — 同仮定で `∀ s, physicalTime t₀ s₀ L s ≤ blowupTime t₀ s₀ L`
  - `NSSingularity.physicalTime_lt_blowupTime`
    — `(hpos : ∀ σ, 0 < L σ)` を追加した狭義版(`T` に到達しない)
  - `NSSingularity.exists_finite_blowupTime`
    — 梱包形 `∃ T, Tendsto … (𝓝 T) ∧ (∀ s, … ≤ T) ∧ MonotoneOn … (Set.Ici s₀)`
  - `NSSingularity.physicalTime_monotoneOn` / `physicalTime_strictMonoOn`
    — 局所仮定 `∀ b, s₀ ≤ b → IntervalIntegrable (scaleRate L) volume s₀ b`
      のもとでの単調性 / `L>0` 追加時の狭義単調性
  - `NSSingularity.tendsto_physicalTime_atTop`
    — **逆向き**: 局所可積分かつ `¬ IntegrableOn (scaleRate L) (Set.Ioi s₀)`
      なら `Tendsto (physicalTime t₀ s₀ L) atTop atTop`
  - `NSSingularity.hasDerivAt_physicalTime`
    — `Continuous L` のとき `HasDerivAt (physicalTime t₀ s₀ L) (scaleRate L s) s`。
      定義が実際に \(dt/ds=L^2\)、\(t(s_0)=t_0\) を解くことの確認
  - 補助: `intervalIntegrable_scaleRate`、`physicalTime_le_physicalTime`、
    `le_physicalTime_of_le`、`physicalTime_le_of_le_base`、`physicalTime_base`、
    `scaleRate_nonneg`、`scaleRate_pos`
- **定義:** `scaleRate L σ = L σ ^ 2`、
  `physicalTime t₀ s₀ L s = t₀ + ∫ σ in s₀..s, scaleRate L σ`(**区間積分**)、
  `blowupTime t₀ s₀ L = t₀ + ∫ σ in Set.Ioi s₀, scaleRate L σ`(**`Ioi` 上の
  集合積分**)。前者は FTC が、後者は広義積分収束定理が mathlib でその形を
  取るための選択で、両者を橋渡しするのが解析的内容の入る唯一の箇所。
  被積分を \(L^2\) と平方で書いたため非負性は `sq_nonneg` から無償であり、
  単調性・有限極限・`t(s) ≤ T` に \(L\) の正値性は**不要**。
- **依存公理:** 本ファイルの全 15 定理について `#print axioms` は
  `[propext, Classical.choice, Quot.sound]`。mathlib 標準の古典公理のみ。
- **mathlib 基盤(実際に使用):**
  `MeasureTheory.intervalIntegral_tendsto_integral_Ioi`
  (`Mathlib/MeasureTheory/Integral/IntegralEqImproper.lean`)、
  `MeasureTheory.integrableOn_Ioi_of_intervalIntegral_norm_bounded`(逆向き)、
  `MeasureTheory.setIntegral_pos_iff_support_of_nonneg_ae`(狭義単調性)、
  `intervalIntegral.integral_mono_interval`、
  `intervalIntegral.integral_add_adjacent_intervals`、
  `intervalIntegral.integral_nonneg`、`intervalIntegral.integral_symm`、
  `intervalIntegral.integral_hasDerivAt_right`、
  `MeasureTheory.IntegrableOn.intervalIntegrable`、
  `MeasureTheory.intervalIntegrable_iff_integrableOn_Ioc_of_le`、
  `ge_of_tendsto`、`Real.volume_Ioc`。
- **形式化していないこと:** これは**与えられた `L` についての命題**であり、
  `L` を構成しない。Navier–Stokes の解との接続は皆無で、`ClayStatement.lean`
  とも未接続。仮定 \(\int_{s_0}^\infty L^2<\infty\) は本ファイルでは検証されず、
  将来の区間証明書が厳密上界として供給すべき対象そのものである。
  また \(t(s)<T\) は物理時刻が \(T\) に達しないことのみを言い、\(T\) での
  ノルム発散(blow-up criterion)は含まない。
- **区間証明書へ:** \(\int L^2\) の数値上界を厳密上界へ変換する部分
  (= `IntegrableOn (fun σ => L σ ^ 2) (Set.Ici s₀)` の供給)。

### F-3 速度回復式と発散ゼロ(E-14/E-15)— **形式化済み**(2026-07-28)

- **状態:** Lean 4 で証明完了。`formal/NSSingularity/VelocityRecovery.lean`
  (`lake build` 成功、`sorry`・`admit`・新規 `axiom` なし)。
- **命題:** \(C^2\) な \(\psi_1\) に対し \(u^r=-r\psi_{1,z}\)、
  \(u^z=2\psi_1+r\psi_{1,r}\) は \(r\neq0\) で
  \(\partial_ru^r+u^r/r+\partial_zu^z=0\) を満たし、\(u^r/r\) を連続延長
  \(-\psi_{1,z}\) に置き換えれば軸 \(r=0\) を含む全点で 0。
- **Lean 定理名:**
  - `NSSingularity.divergence_of_recovered_velocity_eq_zero`
    — `(ψ : ℝ × ℝ → ℝ) (hψ : ContDiff ℝ 2 ψ) (p : ℝ × ℝ) (hr : p.1 ≠ 0) :
      partialR (uR ψ) p + uR ψ p / p.1 + partialZ (uZ ψ) p = 0`(E-15 本体)
  - `NSSingularity.divergence_of_recovered_velocity_eq_zero'`
    — 軸込みの形(`uROverR ψ p = -(partialZ ψ p)` を使用、`r ≠ 0` 不要)
  - `NSSingularity.mixed_partial_comm`
    — `partialR (partialZ ψ) p = partialZ (partialR ψ) p`。
      滑らかさ仮定が入る**唯一**の箇所(Schwarz/Clairaut)
  - 補助: `partialR_uR`、`partialZ_uZ`、`uROverR_eq_div`、
    `fderiv_apply_const`、`differentiable_fderiv_of_contDiff_two`
- **定義:** `partialR f p = fderiv ℝ f p (1,0)`、
  `partialZ f p = fderiv ℝ f p (0,1)`(`ClayStatement.lean` の
  `partialDeriv` と同規約)、`uR ψ p = -(p.1 * partialZ ψ p)`、
  `uZ ψ p = 2 * ψ p + p.1 * partialR ψ p`(= E-14)。
  子午面は `ℝ × ℝ`(`p.1 = r`、`p.2 = z`)。
- **依存公理:** `#print axioms divergence_of_recovered_velocity_eq_zero` は
  `[propext, Classical.choice, Quot.sound]`。mathlib 標準の古典公理のみで、
  それ以外の公理には依存しない。
- **mathlib 基盤(実際に使用):**
  `ContDiffAt.isSymmSndFDerivAt : ContDiffAt 𝕜 n f x → minSmoothness 𝕜 2 ≤ n →
  IsSymmSndFDerivAt 𝕜 f x`(`Mathlib/Analysis/Calculus/FDeriv/Symmetric.lean`)、
  `HasFDerivAt.clm_apply`、`HasFDerivAt.mul`、`ContDiff.fderiv_right`。
- **形式化していないこと:** これは選んだ座標表現に関する恒等式であり、
  \(u^r,u^z\) が実際に3次元軸対称ベクトル場の円柱成分であること
  (E-18 の Cartesian 復元、E-24 の同値性)、および軸極限が E-16 の
  偶奇性から従うことは未形式化。`ClayStatement.lean` の `DivergenceFree`
  との接続も未着手。詳細は `formal/README.md` の当該節。
- **意義:** M1 の最初の Lean 補題。E-15 の紙上導出の機械検証であり、
  相殺が混合偏微分の一致のみに依存することが構造的に露出した。

### F-4 候補証明書の有限次元不等式

- **命題:** radii-polynomial / Newton–Kantorovich 型
  \(Y+Zr+\tfrac K2r^2\le r\)(\(Y,Z,K,r\in\mathbb Q_{\ge0}\)、
  \(Z<1\))の検査器と、その成立が縮小写像条件を含意する抽象補題。
- **仮定:** 抽象 Banach 空間での標準仮定(段階 3 で実体化)。
- **mathlib 基盤:** `Rat` 算術、`norm_num`;縮小写像は
  `ContractingWith` が既存。
- **区間証明書へ:** \(Y,Z,K\) の数値を二進有理数上界として出力する
  発見コード側の変換器(独立実装、未着手)。

### F-6 Galerkin エネルギー上界 — **形式化済み**(2026-07-29)

- **状態:** Lean 4 で証明完了。`formal/NSSingularity/GalerkinNoBlowup.lean`
  (`lake build` 成功、`sorry`・`admit`・新規 `axiom` なし)。
- **役割:** Track F(滑らかな外力による Clay (C)/(D) 反例)の**有限モード
  ansatz 族を閉じる除外定理**の中核。数学的全体像は
  `docs/research_notes/track_f_finite_mode_nogo.md`。
- **命題:** 実 inner product 空間 `E` 上の可微分曲線 `u` が
  `⟪u t, u' t⟫ ≤ ‖u t‖ F t` を満たせば `‖u b‖ ≤ ‖u 0‖ + ∫₀ᵇ F`。
  とくに `u' = g + B(u,u) + A u` で `⟪x,B x x⟫ = 0`(エネルギー中立)、
  `⟪x,A x⟫ ≤ 0`(散逸)、`‖g t‖ ≤ F t` なら、`‖u‖` は `[0,T]` 上
  `‖u 0‖ + ∫₀ᵀ F` で一様に抑えられ、`t → T⁻` で `+∞` へ発散しない。
- **Lean 定理名:**
  - `NSSingularity.norm_le_of_energy_inequality` — 解析的中核(Grönwall)
  - `NSSingularity.inner_galerkin_le` — PDE 構造(`EnergyNeutral`+`Dissipative`)
    が入る唯一の箇所
  - `NSSingularity.galerkin_norm_le` — **F-6 本体**
  - `NSSingularity.galerkin_norm_le_of_mem` — `[0,T]` 上の一様上界
  - `NSSingularity.galerkin_not_tendsto_atTop` — 有限時間爆発の不可能性
- **定義:** `EnergyNeutral B := ∀ x, ⟪x, B x x⟫ = 0`、
  `Dissipative A := ∀ x, ⟪x, A x⟫ ≤ 0`。
- **仮定:** 有限次元性は**不要**(状態空間は任意の実 inner product 空間)。
  `F` の非負性も仮定しない(`‖g t‖ ≤ F t` から従う)。
- **依存公理:** 本ファイルの全 5 定理について `#print axioms` は
  `[propext, Classical.choice, Quot.sound]`。
- **mathlib 基盤(実際に使用):** `HasDerivAt.inner`、`HasDerivAt.sqrt`、
  `antitoneOn_of_deriv_nonpos`、
  `intervalIntegral.integral_hasDerivAt_right`、
  `intervalIntegral.continuousOn_primitive_interval`、
  `intervalIntegral.integral_mono_interval`、`real_inner_le_norm`、
  `ContinuousOn.stronglyMeasurableAtFilter`。
- **形式化していないこと:** (a) Navier–Stokes の移流項が実際に
  `EnergyNeutral` であること(= Lemma 1。紙上 2 行、個別モード集合について
  `galerkin_obstruction.py` が厳密整数演算で機械検証)、(b) 有界性から
  `T` を越える滑らかな延長を導く常微分方程式の議論(= `F-7`)、
  (c) 有限次元空間上のノルム同値、(d) `ClayStatement.lean` との接続。
- **区間証明書へ:** 不要(この命題自体は解析的で、数値入力を持たない)。

### F-7 Galerkin ODE の延長(未着手)

- **命題:** `c' = g(t) + B(c,c) + A c`(`g` 連続、`B` 双線型、`A` 線型)の解が
  `[0,T)` 上有界なら、最大解は `T` を越えて延長される。`g ∈ C^∞` なら
  延長も `C^∞`。
- **役割:** `F-6` の有界性を「特異でない」へ変える最後の一歩。
  `track_f_finite_mode_nogo.md` Theorem 1(iii)。
- **mathlib 基盤:** `Mathlib.Analysis.ODE.PicardLindelof`、
  `ODE_solution_unique` 系。多項式右辺は局所 Lipschitz なので追加理論は不要。
- **数値のみの部分:** なし(純粋な常微分方程式論)。

### F-7a / F-7b Galerkin ODE の延長 — **部分的に形式化済み**(2026-07-29 第 4 便)

- **状態:** `formal/NSSingularity/FiniteModeNoGo.lean`。`lake build` 成功、
  `sorry`・`admit`・新規 `axiom` なし。
- **F-7a(端点到達、形式化済み)**: `[0,T)` 上で微分可能・導関数が連続かつ
  有界な曲線は `t → T⁻` で極限を持つ。
  `NSSingularity.exists_tendsto_nhdsWithin_of_norm_deriv_le`。
  補助: `intervalIntegrable_of_continuousOn_bounded`。
  使用した mathlib: `intervalIntegral.integral_eq_sub_of_hasDerivAt`、
  `intervalIntegral.continuousOn_primitive_interval`、
  `Integrable.of_bound`、`intervalIntegrable_iff_integrableOn_Ioo_of_le`。
- **F-7b(自励系の局所延長、形式化済み)**:
  `NSSingularity.exists_local_galerkin_solution`、
  `NSSingularity.contDiff_galerkinField`。使用した mathlib:
  `ContDiffAt.exists_forall_mem_closedBall_exists_eq_forall_mem_Ioo_hasDerivAt₀`、
  `ContDiff.clm_apply`、`ContinuousLinearMap.contDiff`。
- **F-7c(時間依存外力、形式化済み)**: 上記 API は確かに `f : E → E` の自励系
  専用だが、それは*その定理*の性質であって mathlib の ODE API 全体の性質では
  なかった。`IsPicardLindelof` 自体は `f : ℝ → E → E` の時間依存場に対して
  述べられている。下の「F-7c」節を参照。**公理化していない。**
- **論理接続**: `NSSingularity.not_isBreakdownCandidate_of_galerkin` が
  F-6 の状態上界 → 速度上界 → F-7a を連結し、
  「固定有限モード候補は破綻候補になれない」を結論する。
  `IsBreakdownCandidate u T := 0 < T ∧ ¬∃L, Tendsto u (𝓝[<]T) (𝓝 L)`。
- **接続していないこと**: `ClayStatement.lean` との接続。それには
  `V_S` と係数空間の Fourier 同型、`⟨u,(u·∇)u⟩` と `advectionForm` の同一視
  (解析側は mathlib にトーラス関数空間がなく未着手)、および
  Navier–Stokes の局所一意性理論(mathlib になし)が要る。

### F-12 三線型相殺の Fourier 表現 — **形式化済み**(2026-07-29 第 4 便)

- **命題:** `k_i · a_i = 0` を満たす振幅族に対し、共鳴 3 次形式
  `Σ_{p+q+s=0} (a_q·k_s)(a_p·a_s)` は恒等的にゼロ。
- **Lean 定理名:** `NSSingularity.advectionForm_eq_zero`。
  定義 `dotAmp`(双線型内積)、`resonantTriples`、`advectionForm`。
- **証明:** `p ↔ s` の入れ替えが共鳴集合の対合であることを
  `Finset.sum_nbij'` で使い、`2·form = Σ (a_p·a_s)(a_q·(k_p+k_s))` を
  `k_p+k_s = −k_q` と発散ゼロ条件で消す。
- **形式化していないこと:** これは **Fourier 表現での代数的恒等式**であり、
  `∫_{𝕋³} u·(u·∇)u = 0` という多様体上の積分の形式化ではない。
  両者を繋ぐには mathlib にないトーラス関数空間と Fourier 同型が要る。
- 対応する厳密整数演算の機械検証は
  `src/ns_certificate_lab/galerkin_obstruction.py` の
  `verify_trilinear_cancellation`。

### F-13 有限 Fourier 空間のノルム同値 — **形式化済み**(2026-07-29 第 4 便)

- **命題:** (a) 重み付き和 `Σ w_i c_i² ≤ W Σ c_i²`(`w_i ≤ W`)、
  (b) Cauchy–Schwarz `(Σ|c_i|)² ≤ |S| Σ c_i²` とその平方根形。
- **Lean 定理名:** `NSSingularity.weighted_sq_sum_le`、
  `NSSingularity.sq_sum_abs_le_card_mul_sum_sq`、
  `NSSingularity.sum_abs_le_sqrt_card_mul_sqrt_sum_sq`。
- これらは note の `‖u‖_{H^s} ≤ (1+4π²R_S²)^{s/2}‖u‖` と
  `‖∂^α u‖_∞ ≤ (2πR_S)^{|α|}√|S|‖u‖` が使う定数そのものである。
  「有限次元空間ではノルムが同値」という抽象命題ではなく、
  実際に使う明示定数の形で形式化した。
- 使用した mathlib: `Finset.sq_sum_le_card_mul_sum_sq`(Chebyshev)、
  `Real.le_sqrt_of_sq_le`、`Real.sqrt_mul`。

### F-14 / F-15 5 次元 Green 核の動径恒等式 — **形式化済み**(2026-07-29 第 5 便)

- **状態:** `formal/NSSingularity/GreenAndCascade.lean`。`lake build` 成功、
  `sorry`・`admit`・新規 `axiom` なし。
- **F-14:** `greenProfile R = R^(-3:ℤ)` が `f'' + 4f'/R = 0` を満たす
  (`greenProfile_radial_laplace_eq_zero`、微分値は
  `hasDerivAt_greenProfile` / `hasDerivAt_greenProfileDeriv`)。
  使用した mathlib: `hasDerivAt_zpow`、`zpow_neg`、`zpow_natCast`。
- **F-15:** `flux_newtonSlope`(`R⁴ψ'(R) = −m(R)`)と `hasDerivAt_flux`。
- **形式化していないこと:** `Δ₅G₅ = −δ₀` の分布論的形式化。mathlib に展開がなく、
  本リポジトリのどの上界も Dirac 側を使わない。これらは**動径プロファイルの
  1 次元恒等式**である。
- **数値側との対応:** `src/ns_certificate_lab/free_space_recovery.py` の
  微分 tail 上界と、`whole_space_gate.py` の閉形式参照解の出発点。

### F-16 shell 指数領域(仮定明示)— **形式化済み**(2026-07-29 第 5 便)

- **状態:** 同ファイル。`ShellAdmissible γ σ β` を**構造体**として定義し、
  4 条件(帯域幅発散・スペクトル可和・エネルギー有界・散逸可積分・臨界ノルム発散)を
  フィールドとして明示。どれも黙って落とせない。
- **定理:** `ShellAdmissible.bandwidth_lt_one`(`γ < 1`)、
  `ShellAdmissible.sigma_mem`(`σ ∈ Ico (max 0 (2γ−1)) γ`)、
  `not_shellAdmissible_of_one_le`(`γ ≥ 1` なら admissible な点はない)。
- **形式化していないこと:** 4 つの不等式を PDE から導く部分。それは紙上であり
  端点正則性定理を**引用**する(`F-11`)。本定理は指数の算術である。

### F-7c 時間依存外力での局所延長 — **形式化済み**(2026-07-29 第 7 便)

- **ファイル:** `formal/NSSingularity/TimeDependentGalerkin.lean`。
- **定理:**
  - `galerkin_isPicardLindelof` — 射影 Galerkin 場
    `f t x = g t + B x x + A x` に対する mathlib の `IsPicardLindelof` の
    4 条件。実質的な内容は `B x x` の**局所** Lipschitz 評価であり、
    `B x x - B y y = B x (x-y) + B (x-y) y` から定数 `2‖B‖(‖x₀‖+a)` を得る。
    球半径 `a` が必ず入り、大域定数にはできない。
  - `galerkin_local_solution` — 両側局所解の存在。
  - `galerkin_local_solution_of_continuous` — 外力の上界を連続性から取り出す版。
- **使用した mathlib:** `IsPicardLindelof`(**時間依存**構造体)、
  `IsPicardLindelof.exists_eq_forall_mem_Icc_hasDerivWithinAt₀`、
  `ContinuousLinearMap.le_opNorm₂`、`LipschitzOnWith.of_dist_le_mul`。
- **2 経路の比較**(セッション指示による):

  | | 自励化経路 | 直接経路(採用) |
  |---|---|---|
  | 構成空間 | `E × ℝ` | `E` |
  | 追加 instance 義務 | `CompleteSpace (E × ℝ)`、Prod ノルム評価 | なし |
  | 追加定理 | `galerkin_solution_of_autonomised`(70 行) | なし |
  | 第 2 成分の恒等性 | 必要 | 不要 |

  直接経路が仮定・行数ともに厳密に少ない。

### F-7c 還元(放棄した自励化経路の記録) — 第 5 便

- **定理:** `galerkin_solution_of_autonomised`。自励化した場
  `F(x,s) = (g s + B x x + A x, 1)` の局所流 `α` が `(L,T)` を通るなら、
  時間依存 Galerkin 系は `L` を通る局所解を持つ。
  証明は第 2 成分が `s' = 1`、`s(T) = T` から `s(t) = t` であることによる
  (`Convex.norm_image_sub_le_of_norm_hasDerivWithin_le` で導関数ゼロ ⇒ 定数)。
- **状態:** 定理として正しく、公理を含まないので残置するが、F-7c は上の直接
  経路で閉じており、これに依存するものは無い。

### Clay への限定的接続 — **形式化済み**(2026-07-29 第 5 便)

- **定理:** `breakdown_time_set_empty`。固定有限帯域 Galerkin 軌道について
  `{T | IsBreakdownCandidate u T} = ∅`。1 点の否定から**時刻集合の空性**へ
  強めた形。
- **接続していないこと:** `ClayStatement.lean`。必要なのは Fourier 同型、
  `⟨u,(u·∇)u⟩` と `advectionForm` の同一視の解析側、および
  Navier–Stokes の局所一意性理論(いずれも mathlib になし)。

### F-17 / F-18 / F-19 証明書合成層 — **形式化済み**(2026-07-29 第 6 便)

- **状態:** `formal/NSSingularity/CertificateLayer.lean`。`lake build` 成功、
  `sorry`・`admit`・新規 `axiom` なし。
- **F-17:** ポテンシャル誤差 → 速度誤差。`velocity_radial_error_le`、
  `velocity_axial_error_le`。回復が線型なので積の規則は不要。
- **F-18:** 積差恒等式 `ab − a'b' = (a−a')b + a'(b−b')`(`product_difference`)と
  そこから `product_error_le`、`advection_error_le`。
- **F-19:** 短時間 Grönwall。mathlib の
  `norm_le_gronwallBound_of_norm_deriv_right_le` と `gronwallBound_of_K_ne_0` を
  使い、証明書が検査しやすい `(δ+εt)e^{Kt}` 形へ落とす
  (`K` で割らないので微小 Lipschitz 定数でも検査可能)。
- **Clay 制限の梱包:** `FixedBandwidthCandidate` 構造体と
  `breakdown_times_empty` / `reaches_every_time`。
- **形式化していないこと:** 上界の**計算**、`L^∞` 最大値原理、
  `ClayStatement.lean` への橋。

### F-5 最終 Clay 反例命題までの依存関係

段階 0 として `formal/NSSingularity/ClayStatement.lean` に (A)–(D) と
無外力強化版反例命題を **定義のみ**(証明なし・`sorry` なし)で固定する。
依存 DAG:

```text
F-5 (Clay命題定義)
 ← F-3 (発散ゼロ恒等式) … 初期値の許容性検査
 ← F-1 (スケーリング)   … 再スケーリング軌道→物理解
 ← F-2 (有限物理時間)   … T<∞
 ← F-4 (証明書不等式)   … 軌道存在(段階3)
 → FinalTheorem(段階4、未着手)
```

最終証明経路の `sorry`・`admit`・未証明 axiom 禁止は
`LEAN4_VERIFICATION_POLICY.md` のとおり維持する。段階 0 は定義のみで
この制約に抵触しない。

## 公理監査の記録(2026-07-28、P0 Lean gate)

`formal/AxiomAudit.lean`(ライブラリ root からは import されない監査専用
ファイル)を `lake env lean AxiomAudit.lean` で実行した記録:

```text
'NSSingularity.mixed_partial_comm' depends on axioms: [propext, Classical.choice, Quot.sound]
'NSSingularity.divergence_of_recovered_velocity_eq_zero'' depends on axioms: [propext, Classical.choice, Quot.sound]
'NSSingularity.divergence_of_recovered_velocity_eq_zero' depends on axioms: [propext, Classical.choice, Quot.sound]
'NSSingularity.hasDerivAt_physicalTime' depends on axioms: [propext, Classical.choice, Quot.sound]
'NSSingularity.physicalTime_strictMonoOn' depends on axioms: [propext, Classical.choice, Quot.sound]
'NSSingularity.tendsto_physicalTime' depends on axioms: [propext, Classical.choice, Quot.sound]
'NSSingularity.physicalTime_lt_blowupTime' depends on axioms: [propext, Classical.choice, Quot.sound]
'NSSingularity.exists_finite_blowupTime' depends on axioms: [propext, Classical.choice, Quot.sound]
'NSSingularity.tendsto_physicalTime_atTop' depends on axioms: [propext, Classical.choice, Quot.sound]
```

2026-07-29 第 7 便で F-7c の 3 定理をさらに追記し、
`lake env lean AxiomAudit.lean` は**全 46 行**について同じ古典公理
`[propext, Classical.choice, Quot.sound]` のみを報告した(第 6 便は 43 行)。
新規 3 行:

```text
'NSSingularity.galerkin_isPicardLindelof' depends on axioms: [propext, Classical.choice, Quot.sound]
'NSSingularity.galerkin_local_solution' depends on axioms: [propext, Classical.choice, Quot.sound]
'NSSingularity.galerkin_local_solution_of_continuous' depends on axioms: [propext, Classical.choice, Quot.sound]
```

第 3 便で追記した
F-6 の 5 定理は次のとおり:

```text
'NSSingularity.norm_le_of_energy_inequality' depends on axioms: [propext, Classical.choice, Quot.sound]
'NSSingularity.inner_galerkin_le' depends on axioms: [propext, Classical.choice, Quot.sound]
'NSSingularity.galerkin_norm_le' depends on axioms: [propext, Classical.choice, Quot.sound]
'NSSingularity.galerkin_norm_le_of_mem' depends on axioms: [propext, Classical.choice, Quot.sound]
'NSSingularity.galerkin_not_tendsto_atTop' depends on axioms: [propext, Classical.choice, Quot.sound]
```

同日、`git ls-files formal | xargs grep -InE '\bsorry\b|\badmit\b|^[[:space:]]*axiom '`
は文書コメント中の言及のみを返し、コード中の `sorry`/`admit`/新規 `axiom` は
ゼロ。`lake build` は 2026-07-28 に 8659 jobs、F-6 追加後は 8660 jobs、
F-7a/F-7b/F-12/F-13 追加後は 8661 jobs、F-14/F-15/F-16 追加後は 8662 jobs、F-17/F-18/F-19 追加後は 8663 jobs、
F-7c 追加後は 8664 jobs で成功 —
これは**コンパイルジョブ数**であり、「8664 個の定理を証明した」ことを
意味しない。証明済みの命題は本書に列挙されたもののみである。
`lean-toolchain` は leanprover/lean4:v4.32.1、mathlib は v4.32.1 タグに固定。

## 未解決の形式化上の論点

1. mathlib には NS 方程式の定義も局所適切性もない。局所存在
   (Kato 型 mild solution)を自前で形式化するか、証明書側で回避する設計かは
   M5 以降に決定する。
2. 円柱座標での軸正則性(\(u^\theta/r\) の滑らかさ)を Cartesian 滑らかさへ
   翻訳する補題群は、早期に段階1で着手する価値がある。
3. 浮動小数点→厳密データ変換(二進有理数化)の器は `CertificateFormat.lean`
   設計時に固定する。発見コード側の出力仕様(NPZ v2 schema)からの変換器は
   独立実装とする。


### Track P / Gaussian transfer 層 — **形式化済み**(2026-07-30 第 9 便)

- **ファイル:** `formal/NSSingularity/TrackPFourier.lean`(14 定理)、
  `formal/NSSingularity/GaussianTransfer.lean`(7 定理)。
- **TrackPFourier:** Leray 乗数の有限代数(直交性・冪等性・自己共役性・縮小性、
  `k = 0` でも Lean の `0/0 = 0` 規約で恒等写像となり全主張が成立)、
  単一モードの slot-divergence 定理(`k·a = 0` ⇒ 発散ゼロ)、有限三角多項式の
  `ContDiff ℝ ⊤`、**固定帯域軌道と有限帯域初期値の区別**
  (`FixedBandTrajectory → FiniteBandDatum` は自明、逆は反例
  `u t = (1, t)` で棄却: `exists_finiteBandDatum_not_fixedBandTrajectory`)、
  重み付き和の Ḣ 梯子単調性、`trackP_slab_error_le`(control ODE 層との合成)。
- **GaussianTransfer:** 多項式×Gaussian の微分閉包(witness 多項式
  `p' − 2αXp` を明示)、J 連続性の有限不等式
  (`|a³−b³| ≤ 3max²|a−b|`、`|‖u‖³−‖v‖³| ≤ 3(‖u‖+‖v‖)²‖u−v‖`、
  `‖‖u‖•u−‖v‖•v‖ ≤ (‖u‖+‖v‖)‖u−v‖`)、`torus_control_bound`
  (Riccati 比較の Track-P 形への特殊化)。
- **形式化していないこと(モジュール docstring に記録):** H⁴ エネルギー評価
  そのもの(mathlib にトーラス Sobolev/Fourier 等距がない; F6 では仮定として
  消費)、Galerkin tail の作用素恒等式、Ȧ 格子和(Python 層の有限有理計算)、
  EXT-P1/2/3(古典外部定理; **Lean 公理としては決して挿入しない** — 全 Lean
  定理は記述どおり無条件に真)。
- `lake build` 8668 jobs 成功、`#print axioms` 全 96 定理が
  `[propext, Classical.choice, Quot.sound]` のみ。

### Track P chain / Galerkin Picard 層 — **形式化済み**(2026-07-31 第 10 便)

- **ファイル:** `formal/NSSingularity/TrackPChain.lean`(9 定理)、
  `formal/NSSingularity/GalerkinPicard.lean`(5 定理)。
- **TrackPChain:** スラブ連結の有限不等式骨格 — `two_slab_composition`
  (2 スラブ合成、piecewise 中心・半径は文字通り `if t ≤ t₁ then … else …`)、
  `transfer_triangle`(再中心化予算の 3 項三角不等式: tube 終端 + Taylor 剰余 +
  丸め・Leray 射影)、`ChainLink`/`LinkCertified`/`LinkComposable` +
  `chain_composition`(リスト帰納法による n スラブ合成; `LinkComposable` は
  Python の `delta_out = delta_end + transfer` の Lean 転写)、
  `chain_composition_union`(スラブ和集合上の被覆形)、`discrete_gronwall`
  (`x_{n+1} ≤ A x_n + B ⇒ x_n ≤ Aⁿx₀ + B Σ Aⁱ`)、
  `piecewise_radius_le_max`/`le_foldr_max`/`chain_radius_le_foldr_max`
  (連結半径 ≤ スラブ半径の最大)、`taylor_endpoint_remainder_bound`
  (mathlib の Lagrange 剰余定理の特殊化: `|f(t₀+h) − Taylor_m| ≤ M h^{m+1}/(m+1)!`)。
- **GalerkinPicard:** 二次 ODE `u' = A u + B u u` の有限次元局所存在・一意性 —
  `quadratic_field_lipschitzOnWith`(閉球上の明示 Lipschitz 定数
  `‖A‖ + 2‖B‖(‖x₀‖+r)`)、`quadratic_ode_local_solution`(存在区間半幅
  `ε = 1/(L+1)` を明示、mathlib Picard–Lindelöf 経由)、
  `quadratic_ode_unique`。**これは EXT-P1 の Galerkin 半分の有限次元核**であり、
  Python の Picard box テストはその具体的インスタンス。EXT-P1 自体(PDE 命題)は
  未証明のままで、主張しない(docstring に明記)。
- 真の解の per-slab tube(`LinkCertified`)は仮定として入る: それを
  Navier–Stokes に供給するのは解析層 + EXT-P1/2/3(忠実記録、公理化なし)。
- `lake build` 8670 jobs 成功、`#print axioms` 全 110 定理が
  `[propext, Classical.choice, Quot.sound]` のみ。

### Kato 定数 / Chain 解析層 — **形式化済み**(2026-07-31 第 11 便)

- **ファイル:** `formal/NSSingularity/KatoConstant.lean`(7 定理)、
  `formal/NSSingularity/ChainAnalysis.lean`(7 定理)。
- **KatoConstant:** `G₃ ≤ 12√A₄` 証明書の有限代数 — `cube_diff_bound`
  (`|x³−y³| ≤ 3|x−y|max(x,y)²`)、`am_gm_split`(`p²j² ≤ (pj³+p³j)/2`)、
  `shifted_ratio_bound`(`‖m‖≥1 ⇒ ‖m+j‖³ ≤ (1+‖j‖)³‖m‖³`)、
  格子 tail の telescoping(`Σ_{m>N} m⁻⁴ ≤ 1/(3N³)`、`Finset.Icc` 上の
  有限和として)、checker の組立単調性 `g3_assembly_mono`/`g3_of_a4`。
- **ChainAnalysis:** **積分形比較定理**(`integral_comparison`: 連続 φ が
  2 パラメータ積分不等式を満たせば ODE 解 R が上界 — Dini 微分なし、
  mathlib の liminf-slope 境界補題経由で strict-supersolution fencing)、
  その Riccati 具体化 `integral_riccati_comparison`(`f r = ar+br²+e`、
  Lipschitz 定数 `|a|+b(C_φ+C_R)` を明示)— **EXT-P2-INT + Lemma C の
  スカラー半分の Lean 化**。EXT-P3 用の貼り合わせ論理
  (`glued_continuous`、一様連続度による Cauchy 性、`extendFrom` による
  端点延長)、`cond_to_uncond`(条件付き証明書の仮定放電の命題論理;
  公理ゼロ)。
- **MesoscopicDuhamelNoGo:** empty-child 上界の有限代数核
  \((V/U)^2\le2\kappa^2\tau^2c_EM_{\rm eff}/N^3\)、
  \(M_{\rm eff}\le M\) 版、F-13 の有限 Cauchy--Schwarz を再利用した
  \(M_{\rm eff}\le|\operatorname{support}|\) 版。Fourier/Leray/PDE から
  bilinear・Duhamel 仮定を導く bridge は形式化していない。
- **形式化していないこと:** 無限次元 Kato–Ponce 可換子評価そのもの
  (紙上・監査済み、`kato_h3_constants.md` §4; 公理化はしない)、
  延長関数が方程式を満たすことの同定(EXT-P3 の解析半分)。
- `lake build` 8673 jobs 成功、監査 129 定理 = 128 が古典 3 公理のみ +
  `cond_to_uncond` は公理非依存。
