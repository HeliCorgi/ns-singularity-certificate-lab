# ns-singularity-certificate-lab

3次元非圧縮 Navier–Stokes 方程式の有限時間特異点**候補**を、将来の
区間演算・コンピューター支援証明へ接続できる形で研究するための監査可能な
基盤です。

> **重要:** このリポジトリはミレニアム懸賞問題を解決していません。
> 特異点を発見・証明しておらず、収録した実験は既知に滑らかな人工場と
> 非特異な減衰対照だけです。大きな数値、発散らしい回帰、小さい残差は
> 数学的証明ではありません。

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
