# Project status

最終更新: 2026-07-27
状態: **数学式監査、独立Cartesian監査、空間・時間収束を持つ非特異対照を完了。未知候補探索は未開始。**

一様 \(x,y,z\) 格子上の独立な3成分発散、full curl、vector Laplacian、
primitive PDE項別残差を追加した。保存候補のchecksum検証付き再読込から、
既存円柱差分を呼ばないadapterを経てCartesian物理場を検査するend-to-end
経路も実装した。非特異対照では空間格子を固定し、時間刻みだけを
\(\Delta t,\Delta t/2,\Delta t/4\) とする収束系列を保存した。

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
- 旧Poisson試作の \(-\mathcal L_5\) 符号、内部 \(3/r\) 係数、軸行の係数8、
  周期 \(z\)、外側Dirichlet identity行、manufactured pairは静的に再導出して
  整合を確認した。ただし旧solver自体の正しさや安定性は認証していない。

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

### 一様Cartesian独立検証

`cartesian_validation.py` はNumPy以外の数値演算実装を共有せず、一様
\((x,y,z)\) 格子上で3成分を直接差分する。周期解析場の
\(12^3,24^3,48^3\) refinementでは次の隣接観測次数を得た。

| 診断 | 観測次数 |
|---|---|
| divergence | 1.925, 1.981 |
| full curl | 1.944, 1.986 |
| vector Laplacian | 1.964, 1.991 |
| advection | 1.946, 1.986 |
| pressure gradient | 1.985, 1.996 |
| viscous term | 1.964, 1.991 |
| forced primitive defect \(R_0-f\) | 1.951, 1.988 |

保存candidate \(65\times128\) を再読込して一様
\(25\times25\times64\) 格子へ復元した検査では、

- divergence RMS/max: `2.613327e-3 / 6.328976e-3`
- full-curl defect RMS/max: `9.369950e-3 / 2.975832e-2`

となり、RMSと最大誤差の両gateを通った。adapter出力をtest側の閉形式
Cartesian oracleへ直接比較すると、source grid
\(33\times64\to65\times128\) で速度誤差は
`6.810473e-4 → 1.701222e-4`（次数 `2.001`）、全渦度誤差は
`2.501129e-3 → 6.240860e-4`（次数 `2.003`）へ減少した。

円柱radial符号、成分写像、保存 `omega1` の符号、発散汚染に加え、RMSだけ
なら埋もれる一点curl故障と、周期 \(z\) seamの故障も拒否した。非周期
one-sided closureは全境界を含む二次多項式exact testに合格した。

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

### 固定空間格子の時間収束

非特異ガウス旋回拡散対照を \(N_r=513\)、\(R=5\)、\(T=1\) に固定し、
\(\Delta t=0.5,0.25,0.125\)（step数 `2,4,8`）だけを変更した。

| \(\Delta t\) | 重み付き相対 \(L^2\) 誤差 | 最大絶対誤差 | 最終energy/\(z\)長 | 最終最大渦度 |
|---:|---:|---:|---:|---:|
| 0.5 | `8.5860759e-4` | `1.6269454e-3` | `0.2725504862` | `1.3856349981` |
| 0.25 | `2.0489425e-4` | `3.8828669e-4` | `0.2726696354` | `1.3881123155` |
| 0.125 | `4.3822753e-5` | `8.4939737e-5` | `0.2726992741` | `1.3887190094` |

解析解誤差の観測次数は `2.067119, 2.225128`、同一格子上のstep-doubling
差は `6.537258e-4, 1.611782e-4`、その観測次数は `2.020029` だった。
全runのenergyは初期値 `0.3926990819` から減少し、履歴中の最大渦度は
初期値 `2.0` を超えなかった。

補助的な同次境界 \(R=3,4\) 比較で \(r\le1.5\) の全時刻最大差/最終差は、
各刻みについて `2.155502e-8`, `5.382548e-9`, `2.031933e-9` だった。
これは主計算 \(R=5\) の打切り誤差を直接評価または証明する値ではない。
全11 acceptance checksは合格した。

### 自動テスト

最終結果:

```text
69 passed in 3.21s
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
- Cartesian一点故障をRMSだけで見逃す判定
- periodic \(z\) seamの故障
- 保存済み `omega1` の符号反転
- 非正のCartesian RMS/max許容差
- 不正な時間刻み系列とtime-convergence成果物の上書き
- time-convergence manifest欠落・hash不整合を許す実装回帰

## 仮説にすぎないこと

- Type II、異方的二尺度、周期軌道、準定常軌道、connecting orbit が
  有望な探索空間であるという設計判断。
- 動的再スケーリングで安定した候補力学が見つかる可能性。
- AIが探索初期値や低次元構造の提案に役立つ可能性。

これらは未実装の研究案であり、候補の存在を示していない。

## 未確認・未解決

- 元の3次元 Navier–Stokes 方程式に有限時間特異点が存在するか。
- 非自明な候補profile、軌道、または滑らかな初期データからの接続。
- 非線形production solverと、その候補固有の独立な空間・時間収束系列。
- 全空間の領域打切り誤差、楕円Green tail、スペクトル尾部の厳密評価。
- 圧力回復、射影、または原始変数時間発展を別実装すること。
- 有限円柱上の独立 \(-\mathcal L_5\) Poisson solver、その条件数・境界誤差・
  全空間tail。旧試作からは移植していない。
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
- 一様Cartesian gridと、既存円柱operatorを呼ばない独立2次gradient、
  3成分divergence、full curl、vector Laplacian。
- primitive residualの `time`, `advection`, `pressure_gradient`,
  `viscous=-nu*Laplacian(u)`, unforced `total`, forced defect `total-f` の
  成分別配列。
- 保存candidateのverified load後に、固有の \(r,z\) stencilとbilinear補間で
  E-18a速度・E-18b全渦度を一様Cartesian格子へ写すadapter。
- 軸正則性の必要条件検査。
- canonical `<f8` explicit NPZ candidate、v2 schema manifest、
  array/archive SHA-256、必須単位・無次元化・物理解釈。
- config/seed、JSON/CSV diagnostics、checksum、pre-write runtime/source
  provenance。v1は由来欠如を明示するread-only互換。
- 独立解析式を持つ manufactured solution。
- 独立Crank–Nicolson非特異対照、CSV/JSON/NPZ/SVG、manifest。
- 同一空間格子の \(\Delta t,\Delta t/2,\Delta t/4\) 時間収束、解析誤差、
  step-doubling、energy、最大渦度、補助境界感度を保存する実験。
- 旧ZIPの指定6ファイルだけをread-only監査した
  `docs/legacy_reuse_review.md`。旧コードは移植していない。
- 正常系、round-trip、tamper、故障注入、自動上書き防止テスト。
- Python 3.10/3.12でテストと3実験を再生するGitHub Actions workflow。

## 実行した主なコマンド

```text
$env:PYTHONPATH='src'; $env:PYTHONDONTWRITEBYTECODE='1'; python -m pytest -q
$env:PYTHONPATH='src'; $env:PYTHONDONTWRITEBYTECODE='1'; python experiments/run_manufactured.py --config configs/manufactured.json --output-dir outputs/manufactured_v5
$env:PYTHONPATH='src'; $env:PYTHONDONTWRITEBYTECODE='1'; python experiments/run_baseline.py --config configs/baseline.json --output-dir outputs/baseline_v5
$env:PYTHONPATH='src'; $env:PYTHONDONTWRITEBYTECODE='1'; python -m experiments.run_time_convergence --config configs/baseline_time_convergence.json --output-dir outputs/time_convergence_v1
python -m compileall -q src experiments tests
git diff --check
```

加えて、Pythonから candidate/config/JSON/CSV をchecksum付きで再読込し、
baseline manifestの全7 payload、time-convergence manifestの全5 payloadに
ついてhash・byte lengthとmanifest sidecarを検査した。両NPZの有限性と
`allow_pickle=False` 読込、両SVGのXML妥当性、全summary acceptanceも検証した。
3つの現行runのsource fingerprintは
`bac2077ff5e7333e6c0201d4ebcc319c2fb93604fdd2c8c721bad6e01f8333cf`
で一致した。

時間収束scriptを最初にfile pathとして直接起動した試行は、sibling
`experiments` namespaceを解決できず実験開始前に終了した。証拠directoryは
作られなかった。再現コマンドとCIを `python -m experiments.run_time_convergence`
へ統一して再実行し、上記成果物を生成した。

## 出力ファイル

現行成果物 `outputs/manufactured_v5/`:

- `diagnostics.json`, `diagnostics.csv` と各SHA-256 sidecar
- `run_config.json` とSHA-256 sidecar
- `manufactured_candidate.npz`
- `manufactured_candidate.manifest.json` とSHA-256 sidecar

現行成果物 `outputs/baseline_v5/`:

- `summary.json`, `convergence.csv`, `diagnostics.csv`
- `profiles.npz`
- `energy.svg`, `profiles.svg`
- `config.snapshot.json`, `manifest.json` とmanifest SHA-256 sidecar

現行成果物 `outputs/time_convergence_v1/`:

- `summary.json`
- `time_convergence.csv`, `time_diagnostics.csv`
- `final_profiles.npz`
- `config.snapshot.json`
- `manifest.json` とmanifest SHA-256 sidecar

`outputs/manufactured/`, `outputs/baseline/`, `outputs/*_v2/`,
`outputs/*_v3/`, `outputs/*_v4/` は、途中段階の結果を隠さないため保存している。最初の
`manufactured/` はlegacy v1でprovenanceを持たず、現行成果物には用いない。

## 既知の問題

- package operatorの外側 radial stencil は診断用で、全空間境界条件や
  楕円solveを提供しない。
- 旧Poisson試作は指定ファイルの静的監査だけを行った。旧Grid、boundary
  helper、旧 `l5` は指定外のため未確認で、solverは移植していない。
- axis checkは有限個の必要条件を検査するだけで、Cartesian滑らかさの証明ではない。
- manufactured fieldは強制付きであり、未知の無外力解ではない。
- 一様Cartesian checkerは有限box上の2次binary64差分で、candidate adapterは
  bilinear補間を使う。観測収束と故障検出はあるが、補間・離散化・領域打切りの
  厳密上界はない。圧力は入力であり、独立pressure solveはない。
- baselineは旋回のみの特殊な滑らかな対照で、非線形meridional dynamicsを
  試していない。有限長の \(z\)-周期円柱では有限エネルギーだが、同じ
  \(z\)-不変場を \(\mathbb R^3\) 全体へ延長すると総エネルギーは無限なので、
  主対象の全空間有限エネルギー対照ではない。
- 時間収束の主系列は \(R=5\) で固定したが、保存した境界感度は別の
  \(R=3,4\) 補助比較である。主領域の打切り誤差boundではない。
- SHA-256は改変検出であって、数値の正しさ・作者・実数包含を証明しない。
- 現行v5/v1成果物のprovenanceは `git_head=1af29b9...` と
  `git_dirty=true` を記録する。今回の未commit変更をsource fingerprintで
  固定しているが、署名や信頼時刻証明はない。
- `float64` の誤差は測定しただけで、外向き丸め区間にはなっていない。
- 非収束の故障注入は判定器へ与える合成誤差系列である。将来のproduction
  solverには、意図的に壊した時間発展を通すend-to-end拒否試験も必要。
- 文献の適用はまだ存在しない将来候補に対しては行えない。

## 次に行うべき最小の一手

未知候補探索へ進む前に、有限円柱上の
\(-\mathcal L_5\psi_1=\omega_1\) を解く独立Poisson solverを、旧コードの
コピーではなく本流規約から新規実装する。

最初の小さな受入系列は、周期 \(z\)・明示的外側Dirichlet条件に限定し、
軸行の係数8、全体符号、境界行を直接testする。非零解析境界を持つ
manufactured解を少なくとも3解像度で検査し、matrix residualと解析解誤差を
別経路で保存する。条件数、領域半径感度、重み付き安定性は未確認として
分離する。

この楕円gateが通っても、長時間探索、AI最適化、特異点fitへ自動的には
進まない。次の研究判断では、全空間tailと候補用production time integratorの
証明可能な離散化設計を先に再評価する。
