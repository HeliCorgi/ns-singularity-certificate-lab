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
    └── ClayStatement.lean   # 段階 0(定義のみ、sorry なし)
```

以降の段階(1: 有限次元恒等式 F-1〜F-4、2: 解析的な橋、3: 数値証明書、
4: 最終定理)は `docs/formalization_map.md` を参照。
