# Fable 5 追加タスク指示書
## Navier–Stokes Clay ゴールへ向けた数値監査・全空間移行・厳密化ゲート

対象リポジトリ:

- Repository: `https://github.com/HeliCorgi/ns-singularity-certificate-lab`
- Branch: `fable5-mainline`
- このブランチを本流として扱うこと。

---

## 0. 最終目的

このプロジェクトの目的は、コードや図を増やすことではない。

Clay Mathematics Institute の3次元非圧縮Navier–Stokes問題について、次のいずれかへ到達することである。

1. `R^3` 上の所定の滑らか・発散ゼロ・急減衰初期値について、滑らかな解が全時間存在することを証明する。
2. 3次元周期領域について同様の大域正則性を証明する。
3. Clayの条件を満たす滑らかな初期値と必要なら滑らかな外力について、大域滑らかな解が存在しない反例を構成する。
4. より強い結果として、外力なし・滑らか・有限エネルギー初期値から有限時間特異点を構成する。

現在の主攻撃経路は4である。ただし、存在を前提にしないこと。否定的結果や機構の棄却も正式な成果として保存すること。

最終的な反例には、少なくとも次が必要である。

- 元の3次元方程式との厳密な同値性
- 明示的な滑らか・発散ゼロ・有限エネルギー初期値
- `t < T` での滑らかな物理解
- 有限の物理時刻 `T`
- `t ↑ T` での適切な物理ノルムの発散
- 離散化、領域打切り、丸め、スペクトル尾部の厳密な上界
- 区間演算または同等のvalidated numerics
- 最終的なLean 4検証

数値増幅、小さい残差、Houの図との類似、`lake build`成功だけを解決と呼ばないこと。

---

# 1. 現在の成果の語彙水準

現在の結果は **numerical observation** である。

報告されている早期Hou系列:

- amplification: `6.11 -> 12.70 -> 15.63 -> 17.26`
- Houの1536×1536適応格子公表値: 約`20.5235`
- 増幅は単調だが未収束
- `u1`最大値の初期減少と、その後の成長を観測
- `argmax(u1)`が原点方向へ移動
- energy increaseは記録上ゼロ
- circulation maximum-principle defectは細分で低下
- `z`奇対称は強制せず保存
- 2つのPoisson solverはradial方向では独立だが、`z` Fourier処理とgridを共有
- `nu=5e-4`のHou実データでは粘性項が小さく、粘性符号反転を実データだけで識別できない
- 周期`z`が最小波数`k >= 2 pi`を強制するため、現在の小さいradial-wall effectは全空間独立性を意味しない

これらを特異点、収束済み再現、Clayへの反例と呼ばないこと。

---

# 2. 今回の監査で判明した最重要リスク

以下は、次の計算を始める前に検査・修正すること。

## P0-A: Heun/RK2＋中心差分移流の安定性が未証明

現在の既知のproduction設計は、空間一階微分に中心差分、時間積分にexplicit Heun/RK2を使用している。

定数係数純移流

```text
q_t + c q_z = 0
```

を中心差分すると、Fourier modeの固有値は純虚数になる。Heunの増幅因子は

```text
G(z) = 1 + z + z^2/2
```

であり、`z = i alpha`なら

```text
|G(i alpha)|^2 = 1 + alpha^4/4 > 1
```

となる。

したがって、無粘性の中心差分移流に対してHeunは、非零modeを厳密には増幅する。実際の計算では粘性が高波数を安定化する可能性があるが、それは実測snapshotごとに確認する必要がある。物理項の大きさの比だけでは、数値安定性を保証しない。

### 必須タスク

1. frozen-coefficient advection–diffusionのvon Neumann解析をコード化する。
2. 実際の各snapshotから
   - `max |u^r|`
   - `max |u^z|`
   - `dr`
   - `dz`
   - `dt`
   - `nu`
   を取り、全離散波数についてHeun増幅率の最大値を評価する。
3. predictor stageとcorrector stageの両方を評価する。
4. `max |G| > 1 + tolerance`なら、そのrunを「安定性未合格」とする。
5. 比較実装として少なくとも一つ追加する。
   - classical RK4
   - SSPRK3
   - diffusion implicitのIMEX-RK
   のいずれか。
6. 同じ空間離散化・同じdt系列でHeunと比較し、増幅、ピーク位置、core width、energy balanceが一致するか確認する。
7. Heunだけで得た増幅を今後の候補判定に使わない。

テスト例:

- 単一Fourier modeの純移流
- 移流拡散modeの解析増幅因子
- 複数modeの長時間位相・振幅誤差
- 故意に`nu=0`へした場合にHeun＋中心差分の増幅を検出するテスト

---

## P0-B: CFL診断が「実際にstepを選んだ状態」を測っているか確認する

以前の実装では、`dt`はstep前の速度から選び、記録CFLはstep後の速度から計算していた。この場合、表示されたCFLはstep選択時のCFLではない。

Heunではpredictor stateの速度がさらに大きくなる可能性がある。

### 必須タスク

各stepで次を別々に保存する。

- `cfl_pre_state`
- `cfl_predictor_stage`
- `cfl_post_state`
- `viscous_stability_number_pre`
- `viscous_stability_number_predictor`
- 実際にdtを制限した項
- target/snapshotへ着地するためのstep clippingの有無

`cfl_predictor_stage`が閾値を超えたらstepをrejectして縮小し、再実行する設計を検討する。

CFL acceptanceをdiagnostic snapshotだけでなく、全stepのstreaming maximumで判定する。

---

## P0-C: acceptance-critical quantitiesを間引いて判定しない

以前の実装では`diagnostic_stride=25`であり、energy、circulation、CFL、divergence、parityなどの最大値を記録snapshotだけから計算していた可能性がある。

これでは、記録間に発生した違反を見逃す。

### 必須タスク

出力ファイルは間引いてよいが、次のgate値は全accepted stepで計算し、streaming maximum/minimumを保存する。

- energy growth
- circulation maximum-principle defect
- odd symmetry defect
- axis regularity/parity defect
- relative divergence
- wall constraint defect
- pre/predictor/post CFL
- finite-value check
- Poisson algebraic residual
- energy-balance defect

テストで、間引き点の間だけ違反する合成trajectoryを与え、gateが確実に失敗することを確認する。

---

## P0-D: 「grid-scale saturationなし」という表現を使用しない

最後にアップロードされた旧snapshotでは、`u1`正ピークのradial FWHMは概ね次の点数しかなかった。

- 65×128: 約4点
- 129×256: 約5点
- 193×384: 約6点

ピーク位置も軸から約3、5、6 radial cellだった。

最新版で再計算する必要があるが、少なくとも旧結果はgrid-scaleに近い。増幅がplateauしないことは「grid-scale saturationがない」ことを意味しない。

### 必須タスク

全snapshot、全解像度で次を保存する。

- radial/axial FWHM
- 10–90% front thickness
- peakから軸までのcell数
- narrowest physical scale / `dr`
- narrowest physical scale / `dz`
- local gradient scale `|f|/|grad f|`
- Fourierまたは適切な基底でのhigh-frequency tail
- coarse/fineを共通物理gridへ補間した局所`L∞`, weighted `L2`, derivative error
- subgrid quadratic peak locationとpeak value
- HouのMEMに対応する局所mesh effectiveness指標

定量比較の事前条件を設定する。例として、narrowest scaleが十分な点数で表現され、two-grid差とparity defectが同時に低下するまで、amplificationの収束fitを行わない。

点数の閾値は文献・manufactured front testから決め、結果を見て後付けしない。

---

## P0-E: 現在のwall testは`R^3`移行試験ではない

周期`z`では非零Fourier modeが`|k| >= 2 pi/Lz`を持つ。現在`Lz=1`なら低波数gapが大きく、各modeのradial elliptic tailは急速に減衰しうる。

したがって、radial wallを遠ざけて差が小さかったことは、

- 現在の周期問題でradial wall感度が小さい

ことを示すだけであり、

- 非周期`z`
- `R^3`
- 連続低波数スペクトル
- free-space Green tail

に対する壁独立性を示さない。

### 必須タスク

`docs/whole_space_transition.md`のW-A結果を次の語彙へ修正する。

> periodic-z radial-wall sensitivity observation

「whole-space validation」「R3 wall independence」と表現しない。

次に、真の全空間移行gateを実装する。

1. `z`非周期の有限boxを追加する。
2. 初期データを`z`方向にも`C∞`コンパクト台にする。
3. radial envelopeはCartesian smoothnessのため、軸近傍で`r`ではなく`r^2`の滑らかな関数として設計する。
4. `-L5 psi1 = omega1`についてfree-space境界条件を実装する。
   - 5次元Poisson対応を用いたGreen積分
   - Hankel/Fourier変換
   - multipole/artificial boundary condition
   のいずれか。
5. finite boxの外側境界誤差を、domain doublingだけでなく解析tail boundまたは独立Green solverで評価する。
6. `Rmax`と`Zmax`を独立に増やす。
7. 低波数stress testを追加する。
8. periodic-image effectとradial-wall effectを分離する。
9. `R^3`有限エネルギーを、Cartesian体積測度で直接検査する。
10. 全空間候補を有限円柱no-slip解と同一視しない。

このgateが通るまで、現在のHou機構をClayの`R^3`候補と呼ばない。

---

# 3. 増幅率の扱いを修正する

## P1-A: Houの20.5235を収束先としてfitしない

提示された4点

```text
6.11, 12.70, 15.63, 17.26
```

を

```text
A(h) = A_inf + C h^p
```

へ盲目的にfitすると、使用点によって推定が大きく動く。

概算:

- 全4点: `A_inf ≈ 27.38`, `p ≈ 0.54`
- 最初の3点: `A_inf ≈ 28.85`, `p ≈ 0.49`
- 最後の3点: `A_inf ≈ 24.60`, `p ≈ 0.70`

これは真の極限推定ではなく、系列がまだ漸近領域にないことを示す診断である。

Houの`20.5235`は別の適応格子・別のmesh mapを使った外部anchorであり、現在のuniform-grid系列の既知の極限値ではない。

### 必須タスク

- `20.5235`を固定したfitを禁止する。
- referenceを見ないblind extrapolationを複数modelで行う。
- subset、解像度、fit windowを変えた感度を保存する。
- core resolution gateが通るまで`A_inf`を科学的結論にしない。
- 単調接近は記録してよいが、一致または収束と呼ばない。

---

## P1-B: amplificationの分母をgrid依存のまま比較しない

現在のamplificationが各grid上の離散初期最大渦度で正規化されている場合、分母自体が解像度依存である。

既知のE-29初期値について、continuum referenceは概ね

```text
max |u1(0)| ≈ 3265.9863237
max |omega(0)| ≈ 7569.6226982
```

である。

### 必須タスク

各解像度で次を同時に報告する。

- discrete initial max vorticity
- common continuum-reference initial max vorticity
- absolute final max vorticity
- grid-normalized amplification
- common-reference amplification
- 初期最大点の位置誤差

主要な収束比較にはabsolute final valueとcommon-reference amplificationを使用する。grid-normalized ratioは補助値として残す。

---

# 4. 粘性符号とenergy balance

## P1-C: energy non-increasingだけでは粘性項の正しさを検証できない

数値散逸や境界再構成でもenergyは減少しうる。実データで粘性項が小さい場合、粘性符号反転が観測量に現れないこともある。

### 必須タスク

1. `viscosity_sign` faultを追加する。
2. diffusion-dominated manufactured/control problemを追加し、符号反転を必ず拒否する。
3. 無外力no-slip問題について、離散energy balanceを記録する。

連続系の目標形:

```text
dE/dt = -nu ∫ |omega|^2 dx
```

離散系では、少なくとも

```text
(E_{n+1} - E_n)/dt + nu * dissipation_measure
```

のdefectを保存し、空間・時間細分で減少することを確認する。

4. advection、stretching、diffusion、wall reconstructionがenergyへ与える寄与を可能な範囲で項別に分解する。
5. 「energy increase 0.0」を「energy identityが正しい」と読み替えない。

---

# 5. Poisson solverの独立性

2 solverはradial方向の離散化が異なるが、`z` grid、FFT、波数規約を共有している。

### 必須タスク

第三経路を追加する。

候補:

- `z`実空間有限差分＋block tridiagonal solve
- sine/cosine basis on half-period
- dense small-grid reference matrix
- arbitrary precision direct solve

最低限small manufactured gridsで、FFT規約、Nyquist、normalization、periodic seamの単一障害点を破ること。

solver A/Bの一致をcontinuum精度や完全独立性と表現しない。

---

# 6. Hou早期計算を再実行する前の順序

次の順序を守る。

## Gate 1: 数値安定性

- Heun＋中心差分von Neumann監査
- pre/predictor/post CFL
- 全step acceptance monitoring
- RK4/SSPRK3/IMEX比較

## Gate 2: 時間収束

同じ空間gridで少なくとも

```text
dt, dt/2, dt/4
```

またはadaptive toleranceの3段階を比較する。

比較対象:

- absolute max vorticity
- common-reference amplification
- `max |u1|`
- peak location
- radial/axial core width
- energy-balance defect
- circulation defect
- relative divergence
- relative parity
- Poisson cross-check
- high-frequency tail

## Gate 3: 空間収束

- 65×128
- 129×256
- 193×384
- 257×512
- 計算可能なら追加解像度

単なるscalar値だけでなく、common physical coordinates上のprofile差を比較する。

## Gate 4: 全空間移行

- nonperiodic `z`
- compact support in `r,z`
- free-space elliptic solve
- independent `Rmax`, `Zmax` enlargement
- rigorousまたは少なくともa posteriori tail estimate

Gate 1–4が通るまで、中後期成長、blow-up time fit、AI候補探索へ進まない。

---

# 7. Lean 4監査

`F-2`, `F-3`が証明済みとされていても、以下を必ず実行する。

```bash
grep -RInE '\bsorry\b|\badmit\b|axiom ' formal
lake build
```

最終定理と中間定理について:

```lean
#print axioms theoremName
```

を記録する。

### 必須確認

- `lean-toolchain`固定
- mathlib commit固定
- `sorry`, `admit`なし
- 核心を仮定するproject-specific `axiom`なし
- Clay自然言語命題とLean命題の対応表
- pressure、smoothness、decay、finite energy、divergence-free、global timeの定義が忠実
- finite-dimensional identityの証明とPDE存在証明を混同しない
- `8659 jobs`を`8659個の数学定理を証明した`と表現しない

早期に形式化してよいもの:

- rescaling algebra
- finite physical time condition
- exact initial-data smoothness/divergence-free/finite-energy
- certificate schema
- finite rational inequalities

候補発見前にPDE全体の形式化だけへ移行しない。

---

# 8. 次の科学的マイルストーン

上記gateが合格した後、次へ進む。

## Milestone A: `R^3`初期値族

次を満たす明示初期値を生成する。

- axisymmetric with swirl
- Cartesian `C∞`
- divergence-free
- compact supportまたはClay条件を満たす急減衰
- finite energy
- axisで`u_theta = O(r)`
- radial dependenceは軸近傍でsmooth in `r^2`
- exact formulaとrational/interval parameter representationを持つ

この性質を数値サンプリングだけでなく、Lean 4または解析で証明する。

## Milestone B: dynamic rescaling candidate

探索対象:

- rescaled orbit `U(s)`
- scale `L(s)`
- core position
- anisotropic scales
- periodic/quasiperiodic/slowly drifting orbit
- unstable eigendirections

AIは候補生成にのみ使う。候補は明示的なspectral coefficientsへexportする。

## Milestone C: rigorous certificate

- interval residual
- truncation tail
- inverse/operator norm
- nonlinear Lipschitz bound
- Newton–Kantorovichまたはradii polynomial
- orbit existence
- connection from smooth initial data
- finite physical time
- norm blow-up
- Lean 4 final bridge

---

# 9. 停止条件

次のいずれかが起きたら、結果を隠さず停止して報告する。

- Heun＋中心差分の安定性gateが不合格
- RK4/SSPRK3/IMEXと増幅が一致しない
- time refinementで結果が安定しない
- core widthがgrid点数へlockする
- relative parity/divergenceが細分で減らない
- common-grid profile差が減らない
- wallを遠ざけると機構が消える
- `z` periodを増やすと機構が大きく変わる
- free-space solverでHou型成長が消える
- high-frequency tailが減衰しない
- energy-balance defectが収束しない
- Lean最終経路に未証明の核心的公理が必要になる
- 既知の非存在・正則性定理が探索class全体を排除する

否定的結果も正式な成果物としてcommitする。

---

# 10. 作業完了時の必須成果物

作業後、次を更新する。

- `STATUS.md`
- `PLAN.md`
- `docs/equation_audit.md`
- `docs/numerical_stability_audit.md`
- `docs/whole_space_transition.md`
- `docs/formalization_map.md`
- `docs/proof_obligations.md`
- `docs/threat_model.md`

機械可読出力:

- full-step streaming gate summary
- pre/predictor/post CFL
- von Neumann stability scan
- integrator cross-comparison
- time convergence
- core-width and points-per-scale table
- common-grid profile differences
- absolute and common-reference vorticity amplification
- energy-balance defect
- domain/period sensitivity
- Poisson third-path cross-check

全テストを実行する。

```bash
python -m pytest
cd formal
lake build
```

最後の報告では、次を明確に分ける。

1. 数学的に証明したこと
2. 数値的に観測したこと
3. gateに合格したこと
4. gateに失敗したこと
5. 未確認事項
6. Clay最終命題まで残る証明義務
7. 次の最小の一手

質問して停止せず、合理的な範囲で進めること。ただし、不明点を推測で「合格」にしないこと。
