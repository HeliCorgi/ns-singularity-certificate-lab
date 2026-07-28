# Formalization map — Lean 4 検証への対応表

`LEAN4_VERIFICATION_POLICY.md` の要求に基づき、各マイルストーンについて
「最終的に Lean 化する命題 / その仮定 / 現在は数値的にしか確認できていない部分 /
区間証明書へ変換すべき部分 / mathlib で利用可能な基盤 / 新たに形式化が必要な
解析定理 / 最終 Clay 命題までの依存関係」を記録する。

状態語彙は `FABLE5_HANDOFF.md` §12 の 8 段階
(derivation checked / implementation tested / numerical observation / candidate /
validated numerical object / connected physical solution /
finite-time breakdown theorem / Clay-complete result)を用いる。

最終更新: 2026-07-28(branch `fable5-mainline`)

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
     ├─ 有限物理時間       ← L(s) 可積分性(区間証明書)
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

### F-2 有限物理時間条件

- **命題:** \(L:[s_0,\infty)\to(0,\infty)\) 可測に対し
  \(t(s)=\int_{s_0}^sL(\sigma)^2d\sigma\) が
  \(\int_{s_0}^\infty L^2<\infty\) のとき有限極限 \(T\) を持ち、
  \(s\to\infty\) が物理時刻 \(t\uparrow T<\infty\) に対応する。
- **仮定:** \(L\) の可測性・可積分性(候補証明書から供給)。
- **mathlib 基盤:** `MeasureTheory.integral`、単調収束。ほぼ既存。
- **区間証明書へ:** \(\int L^2\) の数値上界を厳密上界へ変換する部分。

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

## 未解決の形式化上の論点

1. mathlib には NS 方程式の定義も局所適切性もない。局所存在
   (Kato 型 mild solution)を自前で形式化するか、証明書側で回避する設計かは
   M5 以降に決定する。
2. 円柱座標での軸正則性(\(u^\theta/r\) の滑らかさ)を Cartesian 滑らかさへ
   翻訳する補題群は、早期に段階1で着手する価値がある。
3. 浮動小数点→厳密データ変換(二進有理数化)の器は `CertificateFormat.lean`
   設計時に固定する。発見コード側の出力仕様(NPZ v2 schema)からの変換器は
   独立実装とする。
