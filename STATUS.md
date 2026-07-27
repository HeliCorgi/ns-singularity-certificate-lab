# Project status

最終更新: 2026-07-28(branch `fable5-mainline`)
状態: **Poisson ゲート統合・2 実装相互検証、Hou 設定の一次資料監査、壁付き非線形 production ソルバ、早期 Hou 実行までを完了。未知候補探索は未開始。**

## 2026-07-28 セッションの追加結果(fable5-mainline)

### 照合と統合

- ユーザ展開の Poisson バンドルは未統合(未追跡)だったため、全 6 ファイルの
  byte-identity を検証して正規パスへ統合した。パッチの二重適用はしていない。
  パッケージング残骸は `archive/poisson_gate_packaging/` に provenance 込みで
  保存。バンドル同梱の旧 snapshot 証拠は
  `outputs/poisson_gate_v1_bundle_snapshot/` として保存した。
- 統合直後の全テスト: **119 passed**(Python 3.11.9, Windows)。
- Poisson ゲート新規実行 `outputs/poisson_gate_fable5`: 全 7 受入検査合格、
  観測次数 1.9569/1.9782、manifest と全 payload の SHA-256 検証合格。
  バンドル記録と有効数字 10 桁一致(ビット単位では環境差)。

### Poisson 2 実装の相互監査と相互検証テスト

- `poisson.py`(\(r^3\)-flux 有限体積)と `finite_cylinder_poisson.py`
  (非発散形直接差分)の規約は完全一致(負作用素、軸係数 8、外側 Dirichlet
  の意味論)。数学的中核に欠陥なし。
- **独立性は部分的**: radial 離散化と Thomas 解法は真に独立、
  z 方向 Fourier 処理(`numpy.fft.fftfreq` の波数配列がビット同一)と
  `AxisymmetricGrid` は共有の単一障害点。この限定を文書と PLAN に明記した。
- `tests/test_poisson_cross_validation.py`: CV-1(同一入力での
  \(O(\Delta r^2)\) 実測一致 \(D\approx0.115\Delta r^2\)、独立性が崩れると
  発火する下限クローズ付き)、CV-2(radially exact 場での丸め一致
  \(\le4\times10^{-15}\) — FFT 規約・Nyquist・符号・軸係数のピン留め)、
  CV-3(対故障注入)、import 独立性ガード。
- 監査で検出した周辺欠陥を修正: B の複素入力の暗黙切り捨て(D2)、
  `run_poisson_gate` の config 未検証・provenance 欠如・解像度数規則
  (D3/D6/D7)と実験テスト不在(D4)、B の行列が row 1 で M-matrix で
  ないことの明文化+構造ピン留めテスト(D1)。
- この時点の全テスト: **146 passed**。

### Hou (arXiv:2107.06509) 一次資料監査

- v1・v2 の LaTeX 原文と、数値手法の委譲先 Hou–Huang (arXiv:2102.06663) を
  取得して監査した(`docs/hou_setup_audit.md`)。式 (2.1a–d) は E-11–E-14 と
  符号込みで一致。
- 新規監査エントリ E-27–E-31(壁条件 \(\psi_1=0\)・\(u_1=0\)・Thom 型
  \(\omega_1=-\psi_{1,rr}\)、半周期奇対称、初期値式 (2.2)、二段階粘性
  \(5\times10^{-4}\to5\times10^{-3}\) at \(t_0=0.00227375\)、壁渦度の
  2 次離散式)。
- 導出値 \(\|\omega(0)\|_\infty=24000\pi\cdot37^{-1/2}(36/37)^{18}
  \approx7569.62\)、\(\|u_1(0)\|_\infty\approx3265.99\)(論文は比のみ記載の
  ため再現換算に必須。数値最大化で独立検証済み)。
- v1(非爆発の主張)→ v2(potentially singular)の**結論反転**を記録。
  計算設定は両版で同一であり、変わったのは解釈のみ。判定量は後期の
  \(R/Z\) と \(\int\|\omega\|_{L^2}^4\)。
- 論文の誤植 4 箇所、取得不能事項(絶対値時系列、filter 保持の曖昧さ等)を
  記録し、FoCM 出版版の入手をユーザへ依頼した。
- **FoCM 出版版照合(ユーザ提供 PDF、同日追記)**: 計算設定は出版版で
  一切変更なし(本文数値トークン全 219 種の機械照合+項目別逐語照合)。
  実質的変更は「vacuum region」→「no-spinning region」の用語 1 件のみ。
  疑義 4 箇所を含む誤り計 9 件が出版版にも残存し、正誤表なし。
  当方の判読(\(t_3=0.0022868453\) 等)が最終根拠として確定。
  図の解像度地図(最強の爆発判定図の多くは 1024² 実行、Fig. 12 下段は
  Euler 計算で NS ターゲット不可)と、arXiv v2 PDF のベクタ軸目盛から
  回収した絶対値アンカー \(\|u\|_{L^3}\approx46.84\)–\(46.86\) を
  `docs/hou_setup_audit.md` §12 に記録。PDF はハッシュのみ記録し
  コミットしない。

### 非線形 production ソルバ(`nonlinear_cylinder.py`)

- 設計は `docs/nonlinear_solver_design.md` に確定(前提はすべて監査済み式)。
- 実装: Heun/RK2、段ごとの拘束順序(u1 壁ピン → \(\psi_1(1,z)=0\) の楕円
  solve → E-31 Thom 壁渦度 → E-14 速度回復)、壁行は代数拘束として発展
  させない、適応 CFL(0.1、軸係数 4 の拡散余裕)+`max_time_step`、
  二段階粘性 schedule、schema v2 checkpoint と restart。
- テスト 39 件+実験テスト 28 件: manufactured 空間次数(u1: 1.990/1.998、
  ψ1: 2.005/2.003、ω1: 1.845/1.902 — sup 誤差は E-31 壁行に乗り 2 へ下から
  接近)、時間次数 ≈2.0、E-31 単体 1.950/1.975、零場不動点、小振幅の線形
  拡散極限、z 奇対称保存(丸めレベル 1.2e-15 相対)、軸 parity defect の
  閉形式 \(O(\Delta r^3)\) 一致、循環・エネルギー単調性、故障注入 5 種
  (検出比 18.7〜4865)、restart 忠実性。
- **既知の前登録済み注意**: 全振幅 12000 では離散循環最大原理が
  \(O(10^{-4})\) 相対で破れる(中心差分移流+陽的 RK2 の離散化 artifact、
  細分で 2 次超で消失)。受入閾値 1e-3 は前登録値であり、超過時は
  acceptance 失敗として正直に報告される。
- この時点の全テスト: **193 passed**。

### 早期 Hou 実行(`outputs/hou_early_time_v1`)

E-29 初期値(振幅 12000)、\(\nu=5\times10^{-4}\)(E-30 第 1 段のみ、
\(t_0\) 前なので切替なし)、フル周期 \(z\)、\(t\in[0,T_1=0.002191729]\)、
CFL 0.1、`max_time_step` 1e-6。全 8 受入検査合格。manifest+全 50 payload
の SHA-256 検証合格。

| 格子 | steps | \(\min\Delta t\) | 増幅率 \(\|\omega\|_\infty/\|\omega(0)\|_\infty\) at \(T_1\) | \(\max\|u_1\|\) at \(T_1\) | argmax \(u_1\) at \(T_1\) |
|---|---:|---:|---:|---:|---|
| 65×128 | 2192 | 7.29e-7 | 6.1148 | 7605.1 | (0.0469, 0.0156) |
| 129×256 | 2192 | 7.29e-7 | 12.6957 | 14742.9 | (0.0391, 0.0117) |
| 193×384 | 2205 | 2.76e-7 | 15.6280 | 18718.1 | (0.0312, 0.0104) |

観察(すべて **numerical observation** の語彙水準):

- 増幅率は解像度で単調増加し、Hou の 1536² 適応格子公表値 20.5235 へ
  **下から接近**する。grid-scale での飽和はない。ただし隣接差
  (6.58, 2.93)から見かけの収束次数は 1 未満であり、**収束していない**。
  外挿で 20.52 への一致を主張することはできない。
- ごく早期の \(\|u_1\|_\infty\) 減少([Hou21, §2] の定性ターゲット)を全
  解像度で確認: 3265.6 → 最小 2011〜2025 → その後成長し \(T_1\) で初期値の
  4.7〜5.7 倍。
- \(u_1\) の argmax は軸付近・\(z\) 小の領域へ移動(原点方向への伝播と
  定性的に整合)。ただし \(T_1\) での front 位置 \(r\approx0.031\)(193 格子で
  軸から約 6 セル)であり、解像度は限界的。
- エネルギー増加 0.0(全解像度)。循環最大原理の破れは
  7.6e-4 → 2.2e-4 → 4.5e-5 と細分で減少し、前登録閾値 1e-3 以内。
- z 奇対称 defect 比 ≤ 2.0e-9(課さずに監視、保存された)。
- 独立 solver B との \(\psi_1\) cross-check 相対差:
  1.19e-2 → 7.78e-3 → 4.61e-3。急峻化する front 上では見かけの次数が
  2 を下回る(記録のみ、gate ではない)。
- E-02 発散残差最大は 1072〜1764 で単調減少せず、front の解像度不足を
  反映する。初期ノルムの E-29b 一致は 9.95e-3 / 1.98e-3 / 1.07e-3(相対)。

限界: 一様固定格子は Hou の適応最小格子幅 \(O(10^{-8})\) に遠く及ばず、
これは公表増幅率の再現主張ではない。FABLE5_HANDOFF §7.2 の受入基準
(3 空間解像度+3 時間刻み+peak 位置・振幅傾向の一致等)のうち、
時間刻み系列と published-diagnostic 比較の大部分は未実施である。
中成長段(§7.1 stage 2)へ進む前に、より高解像度または適応格子・
半周期 sine 実装の設計判断を行う(「次に行うべき最小の一手」参照)。

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

### 自動テスト(歴史的記録: 2026-07-27 時点のマイルストーン)

**注意: 以下の `69 passed` は 2026-07-27 のマイルストーン時点の歴史的
記録である。現在の正式なテスト数は本書冒頭「2026-07-28 セッションの
追加結果」を参照(統合直後 119 → CV 追加後 146 → 非線形ソルバ追加後
193、以降のセッション内追加はセッション末尾の記録が正)。**

2026-07-27 時点の結果:

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
- Hou の公表増幅率(1536² 適応格子)と本実装(一様固定格子)の定量一致。
  一様格子は Hou の最小格子幅 \(O(10^{-8})\) に遠く及ばない。
- 適応 mesh、半周期 sine 対称実装、非 Fourier の独立 z 方向経路。
- 全空間の領域打切り誤差、楕円Green tail、スペクトル尾部の厳密評価。
- 圧力回復、射影、または原始変数時間発展を別実装すること。
- 候補近傍の非線形安定性、不安定方向、finite physical time、物理ノルム発散。
- 区間演算、validated inverse、radii polynomial、形式証明
  (`docs/formalization_map.md` に Lean 化対応表を開始済み)。
- 古典的な軸対称旋回なし定理の原ロシア語関数空間を、現在のSobolev記法へ
  完全に逐語対応させる作業。
- 二段階粘性を含む Hou プロトコルの \(t_0\) 以後(第 2 粘性段)の再現、
  および壁依存性実験(壁半径を広げた系列)。

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

2026-07-28 セッション(venv `.venv`、Python 3.11.9、Windows 11):

```text
PYTHONPATH=src PYTHONDONTWRITEBYTECODE=1 .venv/Scripts/python.exe -m pytest -q
  -> 119 passed(統合直後)/ 146 passed(CV 追加後)/ 193 passed(非線形ソルバ追加後)
PYTHONPATH=src .venv/Scripts/python.exe -m experiments.run_poisson_gate --config configs/poisson_gate.json --output-dir outputs/poisson_gate_fable5
PYTHONPATH=src .venv/Scripts/python.exe -m experiments.run_hou_early_time --config configs/hou_early_time.json --output-dir outputs/hou_early_time_v1
```

それ以前(2026-07-27 まで):

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

2026-07-28 追加成果物:

- `outputs/poisson_gate_fable5/` — 統合ツリーでの Poisson ゲート新規実行
  (全 7 受入合格、manifest+payload SHA-256 検証済み)
- `outputs/poisson_gate_v1_bundle_snapshot/` — バンドル同梱の旧 snapshot
  証拠(改変なし保存)
- `outputs/hou_early_time_v1/` — 早期 Hou 実行(summary、diagnostics.csv、
  snapshots.csv、trajectories.npz、3 解像度×5 checkpoint、manifest+
  全 50 payload SHA-256 検証済み)
- `archive/poisson_gate_packaging/` — 統合済みバンドルのパッケージング残骸
  (provenance 用 README 付き)

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

早期 Hou 実行の複数解像度結果を評価した上で、次の順に進める。

1. **中成長段の再現判断**: 早期実行で解像度整合な成長傾向が得られた場合に
   限り、\(t_0=0.00227375\) の粘性切替を跨ぐ第 2 段と、より高解像度
   (例 257×512)を別実験として前登録・実行する。得られない場合は
   一様格子の限界として記録し、適応 mesh または半周期 sine 対称実装を
   先に設計する。
2. **壁依存性実験の前登録**: FABLE5_HANDOFF §8.1 の入れ子壁半径系列
   (core を固定解像度で保ちながら \(r_{\max}\) を拡大)の設計と受入基準を
   `docs/` に前登録する。Hou 機構が壁依存かどうかは Clay 目標への
   適用可能性を左右する分岐点である。
3. **非 Fourier の独立 z 経路**: Poisson 相互検証で残った共有単一障害点
   (axial Fourier、grid)を潰す第三経路(実空間差分または sine 基底)を
   検証側に追加する。
4. **formalization_map の更新**: 非線形マイルストーンで確定した命題
   (E-27–E-31 の恒等式群、E-31 の収束次数)を段階 1 形式化候補として
   追記する。

これらが通っても、長時間探索、AI最適化、特異点fitへ自動的には進まない。
動的再スケーリング探索の前に、全空間tailと候補用離散化の証明可能な設計を
再評価する(PO-05〜PO-07)。
