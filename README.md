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
- 実装と独立な解析微分を持つ manufactured solution で収束次数を測る。
- 明示的canonical `<f8` 配列を candidate v2 NPZ、単位・無次元化・物理時刻・
  粘性・基底規約・runtime/source provenance付きmanifest、設定、seed、
  SHA-256とともに保存・再読込する。
- 改変した診断、壊れた発散、楕円符号反転、軸条件破壊などを自動テストで
  拒否する。
- 独立なCrank–Nicolson実装で、非特異な旋回拡散対照を実行する。

非線形の本番時間発展、候補探索、動的再スケーリング、区間演算、厳密な
打切り誤差評価はまだ実装していません。

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

より厳密な再現プロトコルは [docs/reproducibility.md](docs/reproducibility.md)
を参照してください。

## 初回結果

Manufactured solution の \(17/33/65\) radial refinement で得た隣接観測次数:

| 診断 | 観測次数 |
|---|---:|
| 速度回復 | 2.148, 2.103 |
| 物理3D発散 | 2.255, 2.214 |
| 楕円関係 | 2.118, 2.074 |
| Cartesian復元後の独立curl | 3.973, 3.988 |
| \(u_1\) 強制付き残差 | 2.015, 2.010 |
| \(\omega_1\) 強制付き残差 | 2.018, 2.017 |

非特異基準実験の \(33/65/129\) refinement では、相対 \(L^2\) 誤差

\[
8.8040\times10^{-4},\quad2.2100\times10^{-4},\quad5.5310\times10^{-5}
\]

と観測次数 \(1.994,1.998\) を得ました。エネルギーは減少し、内側領域の
境界半径感度は \(8.94\times10^{-11}\)、発散時刻fitは成長条件を満たさない
ため実行されませんでした。詳細は [STATUS.md](STATUS.md) と
[`outputs/baseline_v4/summary.json`](outputs/baseline_v4/summary.json) に
あります。現行manufactured成果物は
[`outputs/manufactured_v4/diagnostics.json`](outputs/manufactured_v4/diagnostics.json)
です。以前の出力は、途中結果を隠さないためそのまま保存しています。

## リポジトリ構成

| パス | 内容 |
|---|---|
| `SPEC.md` | 数学的対象、変数、解・特異点の定義 |
| `docs/equation_audit.md` | 符号・係数・境界・同値性の式別監査 |
| `docs/known_obstructions.md` | 既知の非存在・正則性・継続定理 |
| `docs/threat_model.md` | 偽特異点の原因、検出試験、停止規則 |
| `docs/future_search.md` | Type II・動的再スケーリング探索設計 |
| `docs/proof_obligations.md` | 数値候補から反例までの証明義務 |
| `src/ns_certificate_lab/` | 小さなNumPy数値・保存・診断基盤 |
| `tests/` | manufactured、round-trip、故障注入 |
| `experiments/` | 安価な監査・非特異対照 |
| `configs/` | 固定された実験入力 |
| `outputs/` | 初回の機械可読診断とグラフ |
| `certificates/` | 将来の明示候補証明書用（現在候補なし） |

## 研究上のゲート

既知障害は [docs/known_obstructions.md](docs/known_obstructions.md)、誤検出対策は
[docs/threat_model.md](docs/threat_model.md)、証明までの依存関係は
[docs/proof_obligations.md](docs/proof_obligations.md) を参照してください。
ニューラルネットを将来使う場合も、保持する候補は明示基底係数へ変換し、
ネットワークを使わない独立残差評価に合格させます。

プロジェクトの現在地と未解決事項は [STATUS.md](STATUS.md) にのみ集約し、
数値的確認と数学的証明を混同しません。
