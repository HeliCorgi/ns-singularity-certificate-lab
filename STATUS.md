# Project status

最終更新: 2026-07-27
状態: **数学式監査、初期数値監査、非特異対照を完了。未知候補探索は未開始。**

保存候補からのCartesian成分復元と、それを再投影して行う独立
\(r,z\)-stencilのazimuthal curl検査までは完了した。一様 \(x,y,z\) 格子上の
独立な3成分発散・full curl・primitive PDE residualは未実装なので、将来候補
に必要な上位の独立実装ゲートはまだ合格としていない。

このリポジトリは有限時間特異点を発見・証明しておらず、Navier–Stokes
ミレニアム問題を解決していない。

## 数学的に確認できたこと

- 外力なしの3次元非圧縮 Navier–Stokes 方程式から、指定した curl と
  stream-function 規約の下で軸対称・旋回ありの成分式を導出した。
- \(u_1=u^\theta/r\)、\(\omega_1=\omega^\theta/r\)、
  \(\psi_1=\psi^\theta/r\) の閉じた系、粘性の \(3/r\) 係数、source の符号、
  速度回復、楕円式を代数導出し一次資料と照合した。
- 物理3次元発散
  \(\partial_ru^r+u^r/r+\partial_zu^z=0\) と、形式的な5次元スカラー作用素
  \(\mathcal L_5\) を分離した。
- 軸上の偶奇性・極条件、軸作用素の極限、物理次元、Navier–Stokes
  スケーリング、Cartesian復元式、物理エネルギー測度を監査した。
- 後方一尺度自己相似、Type I、BKM型、Serrin型、軸対称旋回なし、
  有限/局所エネルギー、渦度可積分性に関する既知障害を、一次資料を優先して
  仮定と適用範囲付きで整理した。本文まで確認したものと書誌・要旨のみを
  確認したものは `REFERENCES.md` で区別している。

これらは方程式・既知定理の確認であり、特異点の存在確認ではない。

## 数値的に確認できたこと

### Manufactured solution

独立に書いた解析微分を参照値とし、\(N_r=17,33,65\) で全診断の誤差が
減少した。

| 診断 | 誤差 \(17\to33\to65\) | 隣接観測次数 |
|---|---|---|
| 速度回復 RMS | 6.223e-3 → 1.404e-3 → 3.267e-4 | 2.148, 2.103 |
| 物理発散 RMS | 1.128e-2 → 2.364e-3 → 5.097e-4 | 2.255, 2.214 |
| 楕円 defect RMS | 3.963e-2 → 9.131e-3 → 2.169e-3 | 2.118, 2.074 |
| Cartesian復元後の独立curl defect RMS | 1.065e-5 → 6.778e-7 → 4.272e-8 | 3.973, 3.988 |
| \(u_1\) forced residual RMS | 1.428e-2 → 3.532e-3 → 8.768e-4 | 2.015, 2.010 |
| \(\omega_1\) forced residual RMS | 4.690e-2 → 1.158e-2 → 2.862e-3 | 2.018, 2.017 |

軸 parity は合格し、保存した配列・設定・seed・診断はchecksum検証後に
同一値で再読込できた。candidate/run-config v2は単位、無次元化、物理時刻、
粘性、基底規約、Python/NumPy/platform、Git状態、実行入力のsource
fingerprintを記録し、同一runの3成果物が同じpre-write provenanceを持つことも
確認した。

### 非特異基準実験

滑らかな旋回のみのガウス拡散解を、production operatorを呼ばない独立な
Crank–Nicolson法で計算した。

- \(N_r=33,65,129\) の相対 \(L^2\) 誤差:
  `8.804028753e-4`, `2.209989848e-4`, `5.531034051e-5`
- 観測次数: `1.994124`, `1.998419`
- 最細格子のエネルギー/単位 \(z\) 長:
  `0.3926991427 → 0.3562013981`
- 最大の相対エネルギー増加: `0.0`
- 半径 \(R=2,3\) の同次外側境界を比較した \(r\le1\) の最大差:
  `8.94057e-11`
- peak physical vorticity は持続増大せず、blow-up fit は開始条件で拒否。
- 全7 acceptance checks: 合格。
- \(L_0=U_0=1\)、\(L_0/U_0=1\)、\(Re=20\)、単位、有限周期
  \(L_z=2\pi\) を設定snapshotとsummaryへ記録。

これは既知に滑らかな負の対照であり、一般解の正則性を証明しない。

### 自動テスト

最終結果:

```text
42 passed in 2.01s
```

要求された故障注入はすべて検出した。

- 発散ゼロを壊した速度
- 符号反転した楕円関係
- 軸条件を破った場
- 不正/改変候補archive
- 改変診断データ
- 非収束の解像度系列
- 既存証拠ディレクトリの上書き試行
- 再署名された非canonical `float32` candidate
- 再署名された不正provenance、JSON/CSVのNaN・Infinity、JSON重複key
- 独立4次監査に不足する4点radial grid

## 仮説にすぎないこと

- Type II、異方的二尺度、周期軌道、準定常軌道、connecting orbit が
  有望な探索空間であるという設計判断。
- 動的再スケーリングで安定した候補力学が見つかる可能性。
- AIが探索初期値や低次元構造の提案に役立つ可能性。

これらは未実装の研究案であり、候補の存在を示していない。

## 未確認・未解決

- 元の3次元 Navier–Stokes 方程式に有限時間特異点が存在するか。
- 非自明な候補profile、軌道、または滑らかな初期データからの接続。
- 非線形production solver、時間刻みだけを独立に変えた収束系列。
- 全空間の領域打切り誤差、楕円Green tail、スペクトル尾部の厳密評価。
- 一様Cartesian \(x,y,z\) 格子で別実装した3成分発散・full curl・
  primitive PDE residual。
- 候補近傍の非線形安定性、不安定方向、finite physical time、物理ノルム発散。
- 区間演算、validated inverse、radii polynomial、形式証明。
- 古典的な軸対称旋回なし定理の原ロシア語関数空間を、現在のSobolev記法へ
  完全に逐語対応させる作業。

## 実装したもの

- install可能なNumPy中心のPython packageとpytest設定。
- uniform axisymmetric grid、2次の \(r,z\) 微分、軸極限。
- stream functionからの \(u^r,u^z\) 回復。
- 物理3D発散、楕円defect、PDE項別残差。
- 保存候補からのCartesian速度復元、独立4次stencilによる
  \(\omega^\theta=\partial_zu^r-\partial_ru^z\) 検査。
- 軸正則性の必要条件検査。
- canonical `<f8` explicit NPZ candidate、v2 schema manifest、
  array/archive SHA-256、必須単位・無次元化・物理解釈。
- config/seed、JSON/CSV diagnostics、checksum、pre-write runtime/source
  provenance。v1は由来欠如を明示するread-only互換。
- 独立解析式を持つ manufactured solution。
- 独立Crank–Nicolson非特異対照、CSV/JSON/NPZ/SVG、manifest。
- 正常系、round-trip、tamper、故障注入、自動上書き防止テスト。
- Python 3.10/3.12でテストと両実験を再生するGitHub Actions workflow。

## 実行した主なコマンド

```text
$env:PYTHONPATH='src'; $env:PYTHONDONTWRITEBYTECODE='1'; python -m pytest -q
$env:PYTHONPATH='src'; python experiments/run_manufactured.py --config configs/manufactured.json --output-dir outputs/manufactured_v4
$env:PYTHONPATH='src'; python experiments/run_baseline.py --config configs/baseline.json --output-dir outputs/baseline_v4
```

加えて、Pythonから candidate/config/JSON/CSV をchecksum付きで再読込し、
baseline manifestの全7 hash・byte length、そのmanifest sidecar、NPZ有限性、
両SVGのXML妥当性を検証した。両最終runのsource fingerprintは
`ef2978fb8720090f2e4db341978e7a0c204c647ececd89b2dd2e235ad5015fd0`
で一致した。

## 出力ファイル

現行成果物 `outputs/manufactured_v4/`:

- `diagnostics.json`, `diagnostics.csv` と各SHA-256 sidecar
- `run_config.json` とSHA-256 sidecar
- `manufactured_candidate.npz`
- `manufactured_candidate.manifest.json` とSHA-256 sidecar

現行成果物 `outputs/baseline_v4/`:

- `summary.json`, `convergence.csv`, `diagnostics.csv`
- `profiles.npz`
- `energy.svg`, `profiles.svg`
- `config.snapshot.json`, `manifest.json` とmanifest SHA-256 sidecar

`outputs/manufactured/`, `outputs/baseline/`, `outputs/*_v2/`,
`outputs/*_v3/` は、途中段階の結果を隠さないため保存している。最初の
`manufactured/` はlegacy v1でprovenanceを持たず、現行成果物には用いない。

## 既知の問題

- package operatorの外側 radial stencil は診断用で、全空間境界条件や
  楕円solveを提供しない。
- axis checkは有限個の必要条件を検査するだけで、Cartesian滑らかさの証明ではない。
- manufactured fieldは強制付きであり、未知の無外力解ではない。
- baselineは旋回のみの特殊な滑らかな対照で、非線形meridional dynamicsを
  試していない。有限長の \(z\)-周期円柱では有限エネルギーだが、同じ
  \(z\)-不変場を \(\mathbb R^3\) 全体へ延長すると総エネルギーは無限なので、
  主対象の全空間有限エネルギー対照ではない。
- SHA-256は改変検出であって、数値の正しさ・作者・実数包含を証明しない。
- 現行v4成果物は初回Git commitの直前に生成したため、そのprovenanceの
  `git_head` は `unborn/not-a-git-checkout`、`git_dirty` はtrueである。
  その後、成果物を含むリポジトリ全体を初回commitへ記録した。実行入力内容は
  上記source fingerprintで固定しているが、署名や信頼時刻証明はない。
- `float64` の誤差は測定しただけで、外向き丸め区間にはなっていない。
- 非収束の故障注入は判定器へ与える合成誤差系列である。将来のproduction
  solverには、意図的に壊した時間発展を通すend-to-end拒否試験も必要。
- 文献の適用はまだ存在しない将来候補に対しては行えない。

## 次に行うべき最小の一手

未知候補探索へ進む前に、現在の円柱差分コードを共有しない一様
Cartesian \(x,y,z\) または係数空間実装で3成分 divergence、full curl、
primitive PDE residualを再計算する。その後、非特異対照で時間刻みだけを
\(\Delta t,\Delta t/2,\Delta t/4\) と変える収束試験を追加する。

この二つが通るまで、長時間探索、AI最適化、特異点fitは開始しない。
