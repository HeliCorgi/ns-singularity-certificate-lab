# 再現手順(完全版)

本書は「新しい環境でこのリポジトリの成果を最初から再生する」ための唯一の
完全手順書です。README には短いクイック検証ブロックと索引だけを置き、
各リプレイの完全なコマンド・所要時間・出力・留保はすべてここに集約します。

> **重要:** このリポジトリはミレニアム懸賞問題を解決しておらず、その進展も
> 主張しません。本書のどの手順も特異点の主張ではありません。証明書チェーンの
> 停止時刻は**手法の限界**であって解の性質ではなく、大きな数値・発散らしい
> 回帰・小さい残差はいずれも数学的証明ではありません。

**本書で使う語彙(README・[登録簿](verified_results.md)と同一。厳守):**

- **Lean-verified** — Lean kernel が定理を検査し、project 固有の未証明 axiom が
  ない場合のみ。
- **certificate-verified** — 独立 Python checker が payload を**完全再計算**する
  場合。Lean 証明ではない。
- **hybrid** — Lean-verified な有限論理 + certificate-verified な有理計算 +
  AUDITED-BUT-NOT-FORMALISED な古典解析、の組み合わせ。
- **conditional** — 名前のついた未証明の外部仮定が残る場合。
- **numerically observed** — 浮動小数点のみ。証明ではない。

Python checker の結果を Lean-verified と呼ぶことは、いかなる文脈でもしません。

急ぐ査読者は §6「最初に走らせる順序」から読んでください。

---

## 1. Python 環境の構築

```console
git clone https://github.com/HeliCorgi/ns-singularity-certificate-lab.git
cd ns-singularity-certificate-lab
git checkout main
python -m venv .venv
```

PowerShell / Windows:

```console
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

POSIX(bash / zsh):

```console
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

以降のコマンドはすべて、この仮想環境を有効にしたまま**リポジトリ root** から
実行します。

### 1.1 依存関係の実態

[pyproject.toml](../pyproject.toml) に書かれていることがすべてです。

| 項目 | 実際の内容 |
|---|---|
| `requires-python` | `>=3.10` |
| 実行時依存 | `numpy>=1.24,<3` のみ |
| `[test]` / `[dev]` extra | どちらも `pytest>=8,<10` のみ |
| GPU | 不要 |
| ネットワーク | Python 側は実行時に不要(Lean の `lake exe cache get` のみ必要) |

- **SciPy は依存に入っていません。** 未インストールでも全実験・全テストが動きます。
  SciPy を使うテストは 1 件だけで、
  [tests/test_free_space_poisson.py](../tests/test_free_space_poisson.py) の
  `test_cephes_scaled_bessel_product_matches_scipy_when_available` が
  `pytest.importorskip("scipy.special")` で skip されます。同じ検査を
  リポジトリ内蔵の `I`/`K` オラクルに対して行う
  `test_cephes_scaled_bessel_product_matches_the_in_repo_oracle` は SciPy なしで
  常に実行されるので、**SciPy 不在でも検査自体は失われません**。
  したがって `python -m pytest` の期待結果は「SciPy ありなら skip 0、SciPy なしなら
  skip 1」です。
- `[tool.pytest.ini_options]` は `filterwarnings = ["error"]` を設定しています。
  警告は失敗になります。NumPy のバージョン差で新しい DeprecationWarning が出ると
  テストが赤くなり得るので、赤が出たら**まず警告か本当の失敗かを見分けて**ください。

### 1.2 全テスト

```console
python -m pytest
```

期待結果は全件成功(SciPy 不在なら 1 skip)です。正常系だけでなく、意図的な
故障注入を**検出できた**場合もテスト成功として数えます。テスト件数は開発とともに
変動するため、本書は固定値を書きません。

CI([.github/workflows/tests.yml](../.github/workflows/tests.yml))は Python
3.10 と 3.12 の 2 本で `python -m pip install -e ".[dev]"` → `python -m pytest`
→ 後述のリプレイ 9 本を毎 push 実行します。

---

## 2. Lean 4 環境と公理監査

### 2.1 必要なもの

- **elan**(Lean のツールチェーン管理)。`formal/lean-toolchain` に書かれた
  `leanprover/lean4:v4.32.1` を elan が自動取得します。
- **lake**(elan が同梱)。
- ネットワーク(mathlib の olean キャッシュ取得のため)。

固定されているもの:

| ファイル | 固定内容 |
|---|---|
| [formal/lean-toolchain](../formal/lean-toolchain) | `leanprover/lean4:v4.32.1` |
| [formal/lakefile.toml](../formal/lakefile.toml) | mathlib4 を `rev = "v4.32.1"` で要求 |
| [formal/lake-manifest.json](../formal/lake-manifest.json) | mathlib の git rev を `520045ab14e26149ee970e2e617ca04b09bde5d6` に固定(依存パッケージの rev も同様) |

### 2.2 ビルド

```console
cd formal
lake exe cache get
lake build
cd ..
```

`lake exe cache get` はビルド済み mathlib olean を取得します。**これを省くと
mathlib をソースからビルドすることになり、数時間かかります。** `lake build` の
期待出力は `Build completed successfully` です。

### 2.3 公理監査

```console
cd formal
lake env lean AxiomAudit.lean
cd ..
```

このコマンドは本書執筆時点の HEAD で**実際に実行して確認済み**です。
[formal/AxiomAudit.lean](../formal/AxiomAudit.lean) はライブラリ root
(`NSSingularity.lean`)から import されていないため、**`lake build` だけでは
監査は走りません**。上のコマンドで明示的に elaborate する必要があります。

期待出力は 1 定理 1 行の `#print axioms` 報告です。実測(HEAD `fd1e1c5`):

- 全 **124** 定理。
- うち **123** 定理が `[propext, Classical.choice, Quot.sound]` のみを報告
  (出力上は 121 行が 1 行完結、2 定理は行が折り返されて 3 行で 1 件になります)。
- **1** 定理 `NSSingularity.cond_to_uncond` は
  `does not depend on any axioms`。

受入規則([LEAN4_VERIFICATION_POLICY.md](../LEAN4_VERIFICATION_POLICY.md))は
「各行が高々この標準 3 公理のみを報告し、project 固有 axiom・`sorry`(`sorryAx`)・
`admit` を決して報告しない」ことです。**標準 3 公理の報告は「無公理証明」の主張
ではありません。** これは Lean/mathlib が乗っている通常の古典基盤です。

### 2.4 Lean 側の CI 状況(正確に)

`.github/workflows/` にあるワークフローは
[tests.yml](../.github/workflows/tests.yml) 1 本だけで、**Lean のビルドジョブは
ありません**。CI が毎 push で機械検査している Lean 関連事項は、
[tests/test_readme_claims.py](../tests/test_readme_claims.py) が Python テストとして
行う次の 2 つです。

- `formal/**/*.lean` 全体に `sorry` / `admit` / 行頭 `axiom` 宣言が無いこと
  (`test_no_sorry_admit_axiom_in_lean_sources`)。
- README が引用する Lean 定理名が `formal/` に実在すること。

したがって **kernel チェック自体(`lake build` と公理監査)は、査読者がローカルで
実行して確認する必要があります。** 本書 §2.2–2.3 がその手順です。

---

## 3. リプレイ実験一覧

所要時間の列は 2 種類あります。

- **実測** — 本書執筆時に実際に計測した値(計測環境: Windows 11 / AMD64 /
  20 論理 CPU / CPython 3.11.9 / NumPy 2.4.6)。他機では数倍動きます。
- **見積り** — 本書では計測していない値。README の記載を引き継いだ概算です。

| 実験 | コマンド(リポジトリ root から) | 時間 | 出力先(既存参照物) | 何を示すか |
|---|---|---|---|---|
| 全テスト | `python -m pytest` | 実測 約 11 分(別の重い実験と並走した計測) | — | 正常系 + 故障注入の検出 |
| 3. manufactured 監査 | `python experiments/run_manufactured.py --config configs/manufactured.json --output-dir outputs/manufactured_replay` | 実測 1 秒未満 | `outputs/manufactured_v5` | 強制付き滑らか人工場での収束次数(numerically observed) |
| 4. 非特異基準 | `python experiments/run_baseline.py --config configs/baseline.json --output-dir outputs/baseline_replay` | 実測 1 秒未満 | `outputs/baseline_v5` | 有限周期円柱上の減衰対照(numerically observed) |
| 5. 時間収束 | `python -m experiments.run_time_convergence --config configs/baseline_time_convergence.json --output-dir outputs/time_convergence_replay` | 実測 1 秒未満 | `outputs/time_convergence_v1` | 固定空間格子での時間離散化次数(numerically observed) |
| 6a. Poisson ゲート | `python -m experiments.run_poisson_gate --config configs/poisson_gate.json --output-dir outputs/poisson_gate_replay` | 実測 1 秒未満 | `outputs/poisson_gate_fable5` | 有限円柱 $-\mathcal L_5\psi_1=\omega_1$(非発散形実装)の人工解ゲート |
| 6b. Poisson manufactured | `python -m experiments.run_poisson_manufactured --config configs/poisson_manufactured.json --output-dir outputs/poisson_manufactured_replay` | 実測 1 秒未満 | (参照物なし) | 独立 Poisson 実装の人工解収束 |
| 7. 早期 Hou 実行 | `python -m experiments.run_hou_early_time --config configs/hou_early_time.json --output-dir outputs/hou_early_time_replay` | 見積り 数十分〜数時間 | `outputs/hou_early_time_v1` | E-29 初期値の解像度制限つき数値観察 |
| 8. snapshot Cartesian 監査 | `python experiments/run_hou_snapshot_cartesian_audit.py --config configs/hou_snapshot_cartesian_audit.json --output-dir outputs/hou_snapshot_cartesian_audit_replay` | 見積り 数分〜十数分 | `outputs/hou_snapshot_cartesian_audit_v1` | 保存 checkpoint の独立 Cartesian 検査 |
| 9. Hou 時間刻み細分 | `python -m experiments.run_hou_time_refinement --config configs/hou_time_refinement.json --output-dir outputs/hou_time_refinement_replay` | 見積り 数十分〜数時間 | `outputs/hou_time_refinement_v1` | 強非線形実行の時間誤差と空間誤差の分離 |
| 10. Track F 有限モード除外 | `python -m experiments.run_track_f_finite_mode_scan --config configs/track_f_finite_mode_scan.json --output-dir outputs/track_f_finite_mode_scan_replay` | 実測 約 13 秒 | `outputs/track_f_finite_mode_scan_v1` | 厳密整数演算による**除外定理**の証明書 |
| 11. 全空間 Gate 4 | `python -m experiments.run_whole_space_gate4 --config configs/whole_space_gate4.json --output-dir outputs/whole_space_gate4_replay` | 実測 約 3 秒 | `outputs/whole_space_gate4_v1` | 非周期 $z$ の**線形**楕円ゲート |
| 12. 全空間 Gate 5 | `python -m experiments.run_whole_space_gate5 --config configs/whole_space_gate5.json --output-dir outputs/whole_space_gate5_replay` | 実測 約 4 分 | `outputs/whole_space_gate5_v1` | 微分 tail 上界・速度回復・小振幅非線形・cascade 模型 |
| 13. 全空間 Gate 6 | `python -m experiments.run_whole_space_gate6 --config configs/whole_space_gate6.json --output-dir outputs/whole_space_gate6_replay` | 実測 約 7 分 | `outputs/whole_space_gate6_v1` | 中振幅校正・振幅継続・区間 snapshot 証明書(**前登録 2 基準は不合格**) |
| 13.5 Track P 単発スラブ | `python -m experiments.run_track_p_slab --config configs/track_p_slab.json --output-dir outputs/track-p-replay` | 実測 約 11 秒 | `outputs/track_p_slab_v1` | 周期 $\mathbb T^3$ 上 12 スラブの条件付き証明書 |
| 13.6 Track P チェーン(H⁴) | `python experiments/run_track_p_chain.py` | 見積り 約 2 時間 | `outputs/track_p_chain_v1`(**固定・上書き**) | スラブ連結による certified horizon |
| 13.7a Track P チェーン(n=3 Kato) | `python experiments/run_track_p_chain_h3.py` | 見積り 約 2.5 時間 | `outputs/track_p_chain_h3_v1`(**固定・上書き**) | Kato 定数 `G₃ ≤ 12√A₄` を使った再実行 |
| 13.7b 監査済み再発行 | `python experiments/reissue_chains.py` | 見積り 約 2.3 時間 | `outputs/track_p_chain_reissued_v2`(**固定・上書き**) | 監査閉鎖メタデータ付きの再発行 + 全再検査 |
| 14. 全空間 Gate 7 | `python -m experiments.run_tau_continuation --config configs/tau_continuation_gate7.json --output-dir outputs/tau_continuation_gate7_replay` | 実測 約 27 分(一部テストと並走) | `outputs/tau_continuation_gate7_v1` | Picard 領域からの離脱・τ/Re 継続・時空スラブ証明書 |

**表の読み方に関する注意**

- 13.6 / 13.7a / 13.7b の 3 本には CLI 引数がありません。config と出力先が
  スクリプト内に固定されており、**リポジトリに commit 済みの証拠ディレクトリを
  上書きします**(§5.2)。
- 13.6 / 13.7a / 13.7b の「見積り」は、参照出力の `summary.json` に記録された
  `build_seconds` + `verify_seconds` の合計から算出した値です
  (H⁴ チェーン計 6470 秒、n=3 チェーン計 8998 秒、再発行の再検査計 8321 秒)。
  記録した機械は本書の計測環境とは限りません。
- 上表のスクリプト・config は全件、HEAD の作業ツリー上で実在を確認済みです。

---

## 4. 各実験の詳細

### 4.1 manufactured solution 監査

```console
python experiments/run_manufactured.py --config configs/manufactured.json --output-dir outputs/manufactured_replay
```

速度回復、物理3次元発散、楕円関係、両 PDE 残差の誤差と**観測収束次数**を
JSON/CSV に保存し、最細格子の明示配列を再読込可能な候補形式で保存します。
実装と独立な解析微分を持つ人工場を使うため、これは
**強制付き滑らかな人工場の整合性試験**であって、候補探索でも特異点の証拠でも
ありません。既存証拠を上書きしないよう、新しい出力ディレクトリ名を指定します。

config の細分は半径 17/33/65 の 3 段です。参照値は
[`outputs/manufactured_v5/diagnostics.json`](../outputs/manufactured_v5/diagnostics.json)
にあり、観測次数の要約は README と [STATUS.md](../STATUS.md) に載っています。

### 4.2 非特異基準実験(負の対照)

```console
python experiments/run_baseline.py --config configs/baseline.json --output-dir outputs/baseline_replay
```

滑らかな旋回のみのガウス場を独立な Crank–Nicolson 実装で減衰させ、解析解への
収束、エネルギー、外側境界感度、そして「有限勾配を発散と誤認しない」故障回避
診断を保存します。

**これは有限長の周期円柱上の対照であり、全空間有限エネルギー解の対照では
ありません。** この区別は
[docs/whole_space_transition.md](whole_space_transition.md) の語彙規約に従う
必須の留保です。周期 $z$ は非零 Fourier モードに $|k|\ge 2\pi/L_z$ を強制する
ため、壁を遠ざけた際の差の小ささは主としてこの低波数 gap の帰結であり、
$\mathbb R^3$ の壁独立性を意味しません。

スクリプトは非空の証拠ディレクトリを上書きしません。再実行には別名を使います。

### 4.3 固定空間格子での時間収束

```console
python -m experiments.run_time_convergence --config configs/baseline_time_convergence.json --output-dir outputs/time_convergence_replay
```

同一の 513 点半径格子で時間刻みだけを $0.5, 0.25, 0.125$ と変えます。解析解に
対する重み付き相対 $L^2$ 誤差と隣接次数に加え、共通の空間誤差を概ね相殺する
step-doubling 差、各刻みのエネルギー、最大物理渦度、有限領域の境界感度を
JSON/CSV/NPZ に保存します。**滑らかな負の対照の時間離散化試験であり、特異点の
証拠ではありません。** 記録される境界感度は補助的な $R=3,4$ 比較であり、
主計算 $R=5$ の打切り誤差を直接評価するものではありません。

### 4.4 独立 Poisson ゲート

```console
python -m experiments.run_poisson_gate --config configs/poisson_gate.json --output-dir outputs/poisson_gate_replay
python -m experiments.run_poisson_manufactured --config configs/poisson_manufactured.json --output-dir outputs/poisson_manufactured_replay
```

有限円柱の $-\mathcal L_5\psi_1=\omega_1$ には **2 つの独立実装**があります。
上の 2 コマンドはそれぞれ別の実装の人工解ゲートです。

| コマンド | 使う実装 | 測るもの |
|---|---|---|
| `run_poisson_gate` | `finite_cylinder_poisson.py`(非発散形の直接差分) | 観測次数、cross-stencil defect の次数、外側 Dirichlet 境界誤差、Fourier 行列残差、選択モードの条件数 |
| `run_poisson_manufactured` | `poisson.py`($r^3$-flux 有限体積) | 観測次数、独立残差次数、最細格子 RMS、離散残差・虚部漏れ・境界 defect |

2 実装の**相互検証**そのものは実験ではなくテスト側にあり、
`tests/test_poisson_cross_validation.py` の CV-1($O(\Delta r^2)$ でのみ一致)、
CV-2(両ステンシルが厳密に積分する場では丸め誤差まで一致)、CV-3(符号反転・
外側トレース欠落・軸上限定の故障をそれぞれ検出)として `python -m pytest` で
走ります。

**独立性の範囲を明記します:** 両実装は $z$ 方向 Fourier 処理と格子を共有して
おり、その範囲については独立ではありません。また `run_poisson_gate` の
`summary.json` は自身の限界を
「binary64 で外向き丸めなし」「周期 $z$ の有限円柱」「全空間 Green tail や
領域打切りの評価はない」「条件数は選択した Fourier モードの密行列推定のみ」
と記録します。詳細は
[docs/finite_cylinder_poisson.md](finite_cylinder_poisson.md)。

なお形式作用素 $\mathcal L_5=\partial_{rr}+3r^{-1}\partial_r+\partial_{zz}$ は
スカラー式の略記であり、物理的な 5 次元流体を意味しません。非圧縮条件と体積
測度は常に 3 次元のものです。

### 4.5 早期 Hou 実行

```console
python -m experiments.run_hou_early_time --config configs/hou_early_time.json --output-dir outputs/hou_early_time_replay
```

E-29 の監査済み初期値と E-30 の二段階粘性で、Hou (arXiv:2107.06509) の壁付き
有限円柱を $t=T_1=0.002191729$ まで 3 解像度
(65×128、129×256、193×384)で走らせ、独立楕円 solver B との cross-check、
E-02 発散残差、奇対称 defect、増幅率軌跡を保存します。

**これは一様固定格子上の解像度制限つき数値観察(numerically observed)であり、
Hou の適応格子計算の再現主張ではなく、公表された増幅値の再現でもなく、
特異点候補でもなく、証明でもありません。** この文言は config 自身の
`interpretation` フィールドに前登録されています。

### 4.6 保存 snapshot の独立 Cartesian 検査

```console
python experiments/run_hou_snapshot_cartesian_audit.py --config configs/hou_snapshot_cartesian_audit.json --output-dir outputs/hou_snapshot_cartesian_audit_replay
```

円柱演算子を一切呼ばない独立経路(一様 $(x,y,z)$ 格子)で、保存済み
checkpoint の発散・full curl・$\omega^\theta = r\,\omega_1$ 一致を、絶対値と
相対化指標の両方で検査します。監査 box は half-width 0.7(角が円柱内に残る
最大の丸い値: $0.7\sqrt2=0.98995<r_{\max}=1$)、$z$ は周期全域。監査打切りと
snapshot 欠陥を分離するため、粗細 2 解像度で走ります。

入力は `outputs/hou_early_time_v1/checkpoints/*.npz` で、これらは
**リポジトリに commit 済み**です。したがって §4.5 を先に走らせなくても、
clone 直後にこの監査だけを単独実行できます。

### 4.7 時間刻み収束(強非線形 Hou 実行)

```console
python -m experiments.run_hou_time_refinement --config configs/hou_time_refinement.json --output-dir outputs/hou_time_refinement_replay
```

同一空間格子(65×128)・同一終了時刻 $T_1$ で、固定
$\Delta t,\ \Delta t/2,\ \Delta t/4$ を比較し、時間誤差と空間誤差を分離します。
独立 solver B の楕円 cross-check 付き。**これは出荷している Heun 積分器の時間
離散化誤差の測定であり、増幅値の再現でも特異点候補でも証明でもありません。**

### 4.8 Track F 有限モード除外証明書

```console
python -m experiments.run_track_f_finite_mode_scan --config configs/track_f_finite_mode_scan.json --output-dir outputs/track_f_finite_mode_scan_replay
```

滑らかな外力を使う Clay (C)/(D) 反例の「有限モード ansatz」族について、
三線型形式の相殺

$$\langle u,(u\cdot\nabla)u\rangle = 0$$

を**厳密整数演算**で検証し(浮動小数点を一切使わない)、除外判定を出力します。
Fourier 表現では、$k_i\cdot a_i=0$ の下で共鳴 3 次形式
$\sum_{p+q+s=0}(a_q\cdot k_s)(a_p\cdot a_s)$ が恒等的にゼロになる、という代数
恒等式です(Lean 側の対応定理は `advectionForm_eq_zero`)。

**これは探索の陰性結果ではなく除外定理です。** 詳細は
[docs/research_notes/track_f_finite_mode_nogo.md](research_notes/track_f_finite_mode_nogo.md)。
ただし除外できるのは明示的に限定された ansatz クラスだけで、Clay 命題そのものに
ついては何も主張しません。

### 4.9 全空間 Gate 4(線形楕円ゲート)

```console
python -m experiments.run_whole_space_gate4 --config configs/whole_space_gate4.json --output-dir outputs/whole_space_gate4_replay
```

非周期 $z$ の有限 box 上で $-\mathcal L_5\psi_1=\omega_1$ を解き、**閉形式の
厳密な自由空間参照解**に対して格子細分・領域拡大・尾部上界・周期像分離・独立
Cartesian 検査を測定します。軸方向は FFT を使わない密な離散サイン変換で、既存
ソルバと規約を共有しません。

**これは線形ゲートであり、非線形発展については何も主張しません。**
([docs/whole_space_transition.md](whole_space_transition.md))

### 4.10 全空間 Gate 5(微分 tail・速度回復・小振幅非線形)

```console
python -m experiments.run_whole_space_gate5 --config configs/whole_space_gate5.json --output-dir outputs/whole_space_gate5_replay
```

Green 核の解析微分から導いた**微分 tail 上界**を閉形式参照解に対して検査し、
自由空間速度回復 API の空間・領域収束、軸正則性、独立 Cartesian 検査、故障注入を
測定し、滑らか・コンパクト台・発散ゼロの**小振幅純粋旋回**初期値から非周期 $z$ の
全空間非線形短時間発展を回します。最後に、低周波のみの滑らかな外力が非線形
triad 経由で高シェルを駆動しうるかを有限 cascade 模型で判定します。

詳細:
[green_derivative_tail_bounds.md](research_notes/green_derivative_tail_bounds.md)、
[cascade_toy_model.md](research_notes/cascade_toy_model.md)。

### 4.11 全空間 Gate 6(中振幅校正・振幅継続・区間証明書)

```console
python -m experiments.run_whole_space_gate6 --config configs/whole_space_gate6.json --output-dir outputs/whole_space_gate6_replay
```

境界条件 4 種(zero / monopole / dipole / quadrupole)の core 差を Richardson
離散化誤差と比較して校正し、`dr` / `dz` / joint / `dt` / 積分器 / `Rmax` / `Zmax`
を一因子ずつ分離し、明示的な初期値族について振幅・形状継続を実行して複合ゲートで
順位付けし、動的領域拡大と**厳密有理数区間演算による snapshot 証明書**を生成・
独立検査します。

**前登録基準のうち 2 件は不合格として記録されます。この不合格は消さずに残します。**
参照出力 `outputs/whole_space_gate6_v1/summary.json` の `gate6` ブロックで
`false` を報告するのは

- `boundary_difference_exceeds_richardson`
- `continuation_left_the_quadratic_regime`

の 2 件で、集約フラグ `all_passed` も `false` です。本書執筆時のリプレイは
この `gate6` ブロックを参照出力と**完全一致**で再現しました。
昇格候補はゼロです([whole_space_transition.md](whole_space_transition.md))。

同 `summary.json` の `limitations` は、この不合格について
「校正基準『core 境界差 ≥ 8× Richardson』は**満たされておらず、データを見た後で
閾値を再調整していない**」と明記し、理由を `calibration_finding` に記録します。
他に記録されている限界は「snapshot 証明書以外は binary64」「snapshot 証明書は
離散量のみを包含し、離散化誤差自体は包含しないので PO-05 は open のまま」などです。
**この節の不合格記録を削除・緩和してはいけません。**

### 4.12 Track P 周期スラブ証明書

```console
python -m experiments.run_track_p_slab --config configs/track_p_slab.json --output-dir outputs/track-p-replay
```

周期 $\mathbb T^3$ 上の有理 Fourier 初期値 3 族(P1/P2/P3)について、厳密有理数
演算で Galerkin 軌道の Picard 包含・厳密な連続 NS 残差(= Galerkin tail)・
H⁴ control 不等式・control ODE 管を組み立て、
「真の周期強解がスラブ全体に存在し $\|u-u_a\|_{\dot H^4}\le R(t)$」という
**条件付き**証明書を 12 スラブ分(3 族 × ステップ 1/256, 1/512, 1/1024, 1/2048)
生成し、独立 checker が再検査します。

古典外部定理 EXT-P1/P2/P3 は payload に忠実記録され、**Lean 公理化はしていません**。
**これは特異点証明ではなく、軌道近傍の正則性の証明です。**
([track_p_periodic.md](research_notes/track_p_periodic.md))

### 4.13 Track P スラブ連結(certified horizon)

```console
python experiments/run_track_p_chain.py
```

第 9 便の単発スラブをスカラー H⁴ 誤差半径で連結します。各スラブは**厳密有理
再中心化点**から開始するため、区間 box はスラブ境界を越えて伝播せず、wrapping の
入る場所が構造的にありません。Taylor 終端包絡 + dyadic 丸め + 厳密 Leray 射影で
捨てた幅は $\delta_{n+1}=R_n(h)+\text{transfer}$ としてスカラー半径に課金されます。

P1/P2/P3 × $\nu\in\{1/4,1/10,1/40,1/100\}$ の 12 連結 + 長尺 1 本
(`h = 1/8192`、48 スラブ予算)を前登録 config
([configs/track_p_chain.json](../configs/track_p_chain.json))で実行し、独立
checker が全リンクを再計算します。

**停止は前登録分類法で必ず分類され、証明区間の終了は特異点の主張ではありません。**
([track_p_chain.md](research_notes/track_p_chain.md))

### 4.14 Track P チェーン n=3(Kato 定数)と監査済み再発行

```console
python experiments/run_track_p_chain_h3.py
python experiments/reissue_chains.py
```

第 11 便。正規化完全一致で自前導出した `G₃ ≤ 12√A₄`
([kato_h3_constants.md](research_notes/kato_h3_constants.md)、独立 checker 付き
証明書)と厳密帯域和 `C_kato` / `C_shift` による n=3 control 不等式でチェーンを
再実行し、certified horizon を旧 `9(K₁+K₂)` 比で実測約 11〜13 倍に延長しました
(参照出力 `outputs/track_p_chain_h3_v1/summary.json` の
`comparison_vs_h4_baseline` に記録された 12 組の比は 11.09〜12.57)。

`run_track_p_chain_h3.py` は `outputs/track_p_chain_v1/summary.json` を H⁴ 基準線
として読むので、**§4.13 を先に走らせるか、commit 済みの参照出力を残しておく必要が
あります。**

`reissue_chains.py` は `outputs/track_p_chain_v1` と
`outputs/track_p_chain_h3_v1` の全チェーン payload(計 27 本)を、EXT-P1★ /
EXT-P2-INT + Lemma C / EXT-P3★ / 系 P3-3 の**監査済み紙上証明**閉鎖メタデータ
付きで再発行し、再発行後の payload を独立 checker で**全リンク再計算**します。
量的内容は一切変わりません。

**語彙の厳守事項:** payload の `proved: true` は「監査済み紙上証明」の意味であり、
**Lean 形式化ではありません**(`lean_formalised: false` 固定、公理化は禁止のまま。
checker が新旧文言の混在を拒否します)。旧 Dini 節は G-DINI として open のまま、
どこからも未消費です。**特異点主張ではありません。**
([ext_p1_p2_p3_audit.md](research_notes/ext_p1_p2_p3_audit.md))

### 4.15 全空間 Gate 7(Picard 領域からの離脱・τ/Re 継続)

```console
python -m experiments.run_tau_continuation --config configs/tau_continuation_gate7.json --output-dir outputs/tau_continuation_gate7_replay
```

第 6 便の 32 点スイープを無次元座標 `(Re, aspect, c, τ)` へ再分類し(到達 `τ` は
最大 0.0233 だったことが判明)、Picard 梯子(level 0/1/2 + 完全解を同時積分)で
第一 Picard 反復からの乖離を**実測**し、前登録 $\tau=\{0.025,\dots,1.0\}$ と
$\mathrm{Re}=\{10,\dots,400\}$ × 族 S/A/H の 18 run を実行し、乖離ゲート 9 項目と
昇格 2 基準で判定し、$[t_n,t_{n+1}]$ の**時空スラブ証明書**(cell 内部・全時刻を
包含、厳密有理数、独立 checker + 改竄拒否)を生成します。参照出力の `runs` 配列は
この前登録 18 run に解像度確認用の追加 run を足した 23 行になります。

**乖離ゲートは全項目合格、昇格候補はゼロです。**
([tau_continuation_gate7.md](research_notes/tau_continuation_gate7.md))
本書執筆時のリプレイは `departure_gate` ブロックを参照出力と**完全一致**で再現し、
`promoted` は空、時空スラブ証明書は独立 checker で `verified: true` でした。

2 つの再現上の注意:

- 再分類は `outputs/whole_space_gate6_v1/continuation.csv`(commit 済み)を読みます。
  無い場合は `legacy_sweep_reclassification.available` が `false` になり、実験自体は
  続行します。
- config の `amendments` ブロックに前登録の**修正記録**が 1 件あります: 初回パスは
  前登録どおり 2 件(front 解像度)で不合格になり、閾値は変えずに 145×289 の
  ultra 格子を追加して評価し直した、と記録されています。閾値を動かしていないこと・
  修正が結果を見る前に記録されたことは config 自身が主張しており、査読者は
  `outputs/tau_continuation_gate7_v1/preregistration.json` と
  `configs/tau_continuation_gate7.json` を突き合わせて確認できます。

---

## 5. 決定性・来歴・上書き規約

### 5.1 出力に付く来歴

- **config snapshot** — ほぼ全実験が `config.snapshot.json` を出力先に書き、
  実際に使った入力を固定します。
- **前登録** — Track P スラブと Gate 7 は `preregistration.json` を別ファイルで
  出力し、`recorded_before_any_run` が真でない config は**実行前に拒否**されます
  (`experiments/run_track_p_slab.py` の検証)。Track P チェーンは
  `configs/track_p_chain.json` の `preregistered` / `registered_at` /
  `acceptance` / `lohner_triggers_preregistered` に閾値を前登録しています。
  Gate 6 は config に前登録ブロックを持たず、判定基準はドライバ側に固定され、
  出力の `limitations` に「データを見た後で閾値を再調整していない」と明記されます。
- **runtime provenance v2** — candidate / run-config / diagnostics の v2 writer は
  Python・NumPy・OS・Git 状態・安定な source fingerprint、および BLAS/OpenMP スレッド
  環境変数と NumPy ビルド構成の要約を自動記録します
  ([release_process.md](release_process.md))。

**sha256 の残し方は 3 通りあります。**

| 形式 | 対象 |
|---|---|
| `manifest.json`(ファイル名 → sha256 / バイト数)+ `manifest.json.sha256` | `run_baseline`, `run_time_convergence`, `run_poisson_gate`, `run_poisson_manufactured`, `run_hou_early_time`, `run_hou_snapshot_cartesian_audit`, `run_hou_time_refinement`, `run_track_f_finite_mode_scan`, `run_whole_space_gate4/5/6`, `run_tau_continuation`, `run_track_p_slab` |
| `manifest.json` のみ(manifest 自身のダイジェストなし) | `run_track_p_chain`, `run_track_p_chain_h3`, `reissue_chains` |
| 実行単位の manifest ではなく**ファイル毎の `.sha256` サイドカー**(`diagnostics.json.sha256` など) | `run_manufactured` |

`run_manufactured` の出力にある `manufactured_candidate.manifest.json` は実行単位の
manifest ではなく、保存した candidate v2 NPZ 自身の単位・無次元化・物理時刻・粘性・
基底規約・来歴を記録した**成果物 manifest** です(これにも `.sha256` が付きます)。

**バイト一致は約束されません。** NPZ の圧縮メタデータや浮動小数点ライブラリの差で、
プラットフォームをまたぐとハッシュは変わり得ます。数値は**許容差付きで比較**し、
sha256 は「その 1 回の実行を認証する」目的にだけ使ってください。

### 5.2 上書きガードの実態(スクリプトごとに違います)

| 挙動 | 対象スクリプト |
|---|---|
| **非空**の出力ディレクトリを拒否(`FileExistsError`) | `run_manufactured.py`, `run_baseline.py`, `run_time_convergence.py`, `run_poisson_manufactured.py`, `run_hou_early_time.py`, `run_hou_snapshot_cartesian_audit.py`, `run_hou_time_refinement.py` |
| ディレクトリの**存在自体**を拒否 | `run_poisson_gate.py`, `run_track_f_finite_mode_scan.py`, `run_whole_space_gate4.py`, `run_whole_space_gate5.py`, `run_whole_space_gate6.py` |
| **ガードなし(既存を上書きします)** | `run_tau_continuation.py`, `run_track_p_slab.py`, `run_track_p_chain.py`, `run_track_p_chain_h3.py`, `reissue_chains.py` |

最後の行が重要です。とくに `run_track_p_chain.py` / `run_track_p_chain_h3.py` /
`reissue_chains.py` は **CLI 引数を持たず、出力先が
`outputs/track_p_chain_v1` / `outputs/track_p_chain_h3_v1` /
`outputs/track_p_chain_reissued_v2` に固定**されています。これらは commit 済みの
証拠ディレクトリなので、実行すると**その場で上書きされます**。走らせる前に
`git status` が綺麗であることを確認し、走らせた後の差分は
「同じ payload が再生成されたか」として読んでください。

### 5.3 出力パスの制約

一部のスクリプトは、config と出力先が**リポジトリ内**にあることを要求し、外に
向けると `config and output paths must remain inside this repository` で停止します:
`run_manufactured.py`, `run_baseline.py`, `run_time_convergence.py`,
`run_poisson_manufactured.py`, `run_hou_early_time.py`,
`run_hou_snapshot_cartesian_audit.py`, `run_hou_time_refinement.py`。
`/tmp` などリポジトリ外にリプレイ出力を置くことはできません。

### 5.4 再現性の 3 段階(混同しないこと)

1. **成果物の完全性** — バイトが記録済み sha256 と一致する。
2. **環境リプレイ** — 同じソース・パッケージ・設定で同じ数値が出る。
3. **独立再現** — 別に書かれた実装(できれば別言語)で同じ結論が出る。

本リポジトリの Python Cartesian checker は意図的に分離した経路ですが、
**第 3 段階ではありません**。数値的に近いリプレイは独立実装ではなく、数学的証明
でもありません([release_process.md](release_process.md)、
[proof_obligations.md](proof_obligations.md))。

---

## 6. 最初に走らせる順序(新しい査読者向け)

**第 1 段階(テスト 10 分強 + mathlib キャッシュ取得とビルド。ここまでで
「何が無条件に検証できるか」が分かる)**

1. `python -m pip install -e ".[dev]"` → `python -m pytest`
   — 故障注入を含む全テスト。
2. `cd formal && lake exe cache get && lake build`
   — Lean 側の kernel チェック。
3. `cd formal && lake env lean AxiomAudit.lean`
   — 124 定理の公理報告。project 固有 axiom / `sorry` / `admit` が 1 つも出ないこと。
   ここまでが **Lean-verified** の範囲です。

**第 2 段階(数分。certificate-verified の実物を見る)**

4. §4.8 Track F 除外証明書(約 13 秒、厳密整数演算)。
5. §4.12 Track P 単発スラブ(約 11 秒、厳密有理演算 + 独立 checker)。
6. `python -m pytest tests/test_kato_constant.py tests/test_gaussian_spectral_pressure.py`
   — Kato 定数証明書と $\mathbb R^3$ スペクトル圧力証明書の再計算。

**第 3 段階(合計 40 分程度。数値ゲートの実態と、記録された不合格を見る)**

7. §4.9 Gate 4(約 3 秒)→ §4.10 Gate 5(約 4 分)。
8. §4.11 Gate 6(約 7 分)— **2 件の前登録不合格がそのまま再現すること**を確認。
9. §4.15 Gate 7(約 27 分)— 乖離ゲート全合格・昇格ゼロ。

**第 4 段階(時間がある場合のみ)**

10. §4.1〜4.4 の小規模対照(いずれも 1 秒未満)。
11. §4.13 / §4.14 のチェーン再実行(各 2 時間規模、**commit 済み出力を上書き**)。
12. §4.5 / §4.7 の Hou 実行(数十分〜数時間)。

最初の 3 段階を終えた時点で、README の A 節(Lean-verified)・B 節
(certificate-verified)・C1 節(conditional)・C2 節(hybrid)の区分が、実際に
手元で再現された出力によって裏付けられます。C2/C3 に残る EXT 群と Kato 可換子評価は
**監査済み紙上証明**であり、機械検証ではありません
([ext_p1_p2_p3_audit.md](research_notes/ext_p1_p2_p3_audit.md))。

---

## 参照

- 成果ごとの安定 ID と区分: [docs/verified_results.md](verified_results.md)
- 誤検出の脅威モデル: [docs/threat_model.md](threat_model.md)
- 証明までの依存関係: [docs/proof_obligations.md](proof_obligations.md)
- リリース工程と来歴スキーマ: [docs/release_process.md](release_process.md)
- 既知障害: [docs/known_obstructions.md](known_obstructions.md)
