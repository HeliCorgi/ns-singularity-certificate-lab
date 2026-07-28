# Project status

最終更新: 2026-07-28 第 2 セッション(branch `fable5-mainline`)
状態: **FABLE5_NEXT_TASK_AUDIT の P0 ゲート(von Neumann 安定性監査、全 step streaming acceptance、core-width fit 前提、blind 外挿、エネルギー収支、Gate 1 積分器相互比較)を実装・実行済み。出荷済み Heun 実行は stability-unverified に再分類(時間スキーム依存の実測上界は ~6 ppm)。T₁ 増幅ラダーの収束 fit は前提不合格により機械的に禁止。Gate 4(真の全空間移行)は仕様のみで未実装。未知候補探索は未開始。**

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
  **下から接近**する。ただし隣接差(6.58, 2.93)から見かけの収束次数は
  1 未満であり、**収束していない**。外挿で 20.52 への一致を主張することは
  できない。〔改版 2026-07-28(P0-D): 旧文の「grid-scale での飽和はない」
  という表現を撤回する。増幅が plateau しないことは grid-scale saturation
  の不在を意味しない。実際ピークの radial FWHM は 65/129/193 格子で約
  4/5/6 点、ピークは軸から 3/5/6 セルであり、構造は grid-scale に近い。〕
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

**ラダー延長(`outputs/hou_early_time_v2_hires`、前登録
`configs/hou_early_time_hires.json`)**: 同一プロトコルで 129×256 を再実行し
257×512 へ延長した。全 8 受入検査合格。

- **再現性**: 129×256 の増幅率は前回実行と **bit 単位で同一**
  (12.6956952437)。決定論的再現を確認。
- 257×512: 増幅率 **17.2588** at \(T_1\)(2298 steps、\(\min\Delta t\)
  2.85e-7)、\(\max|u_1|(T_1)=20306.5\)、argmax (0.0312, 0.0098)。
- 4 点ラダー 6.11 → 12.70 → 15.63 → 17.26 は公表 1536² 値 20.5235 へ
  単調接近を続ける。隣接差 6.58, 2.93, 1.63 の見かけの収束次数は依然
  1 未満であり**未収束**。外挿による一致主張は行わない。front 位置
  \(r\approx0.031\) は 257 格子で約 8 セルであり解像度限界が支配的。
  〔改版 2026-07-28(P0-D): 旧文の「grid-scale 飽和はない」という表現を
  撤回。plateau の不在は saturation の不在を意味しない。〕
- 初期ノルム E-29b との相対誤差は 1.98e-3 → 4.99e-4(比 3.97、
  きれいな 2 次接近)。solver B との \(\psi_1\) cross-check 相対差は
  7.78e-3 → 2.95e-3。
- 実効 advective CFL の設定値 0.1 に対する微小超過(+0.16%)を 257 でも
  記録(原因調査と受入条件の明文化は本セッションの別タスク)。

### 拘束誤差の相対化と CFL 受入条件(タスク 1+5)

絶対残差だけでは場の増幅と数値破綻を区別できないため、TM-09 に従い
全拘束診断を相対化した(記録のみ、既存 gate は不変)。

- **相対発散**: E-02 残差最大 / 打消し項和最大
  \(\max(|\partial_ru^r|+|u^r/r|+|\partial_zu^z|)\)。分母・argmax 位置・
  各点比も保存。実測(振幅 12000、短時間 run):
  33×64 で 2.92e-2 → 65×128 で **8.04e-3**(≈O(h²) 改善)。
- **相対軸パリティ**: 軸片側微分 defect / \(\max|\partial_r\) 場\(|\)。
  実測: u1 7.48e-3 → **9.39e-4**、ω1 → 2.75e-3(≈O(h³) 改善)。
- いずれも**離散化オーダーで減少しており、増幅由来の破綻ではない**。
- **CFL 超過の原因確定**: \(\Delta t\) は step 開始時点の速度最大で決定し、
  実効 CFL は step 後の状態で評価するため、超過分は step 内成長率に
  厳密に一致する。v1 実測 +0.227%(193×384、拘束が効き始めた step 2085
  以降)。HH21 の <1%/step 指針内。受入条件を
  \(\text{CFL}\le C(1+\varepsilon)\)、\(\varepsilon=0.05\)(config で変更可)
  として明文化し、\(\varepsilon=0\) では v1 値が不合格になることを
  テストで固定した。

### 保存 snapshot の独立 Cartesian 検査(タスク 2、`outputs/hou_snapshot_cartesian_audit_v1`)

円柱演算子を import しない(AST テストで固定)独立経路で、checkpoint を
一様 Cartesian 箱 \([-0.7,0.7]^2\times[0,1)\) へ復元し検査した。
全 7 hard gate 合格、manifest+payload SHA-256 検証済み。

- t=0: 相対発散 RMS 1.5e-4〜2.5e-4、curl RMS ≈4.6e-4(離散化レベル)。
- \(T_1\): RMS はゲート合格(div 6.1e-4〜9.2e-4、curl 4.9e-3〜5.9e-3)、
  監査格子細分で ≈2 次減少。ただし **pointwise 最大は front 近傍で悪い**
  (curl defect 最大 = 勾配スケールの 0.41、方位一致最大/信号 = 0.62、
  nr193)。監査格子(dx=0.022)が nr193 の source 格子より 4.3 倍粗い
  ことと整合するが、このレベルでの snapshot 欠陥を排除するものではない
  (summary に明記、gate は RMS のみ)。
- 故障注入: u_y 符号反転(×52)、成分入替(×52)、E-18b 符号、
  軸条件違反(奇 r kink、×45; 軸正則な整合対照は不変)を検出。
  発散検査は E-18a 恒等性より ψ1 編集に構造的に盲目であることも記録。
- **primitive NS 残差**: checkpoint 対+圧力から組み立てる関数を実装し
  manufactured 場で時間次数 2.004/2.001 を検証。実 checkpoint への適用は
  「圧力未保存・snapshot 間隔 ≈500 step」のため未実施
  (`primitive_navier_stokes_residual_gap` として全 summary に明記)。

### Poisson 第三経路 solver C(`realspace_poisson.py`)

A/B が共有する axial Fourier(`fftfreq` bit 同一)+Thomas の単一障害点を
破る実空間経路: E-26a の \(r^3\)-flux 行(独立転写)+周期 2 次 z 差分
(`np.roll`)+\(V_i\) 重み付き SPD の Jacobi 前処理 CG。`numpy.fft`
不使用をソース文字列レベルでテスト固定。

- C-vs-A 差 = **0.0901 Δz²**(純 axial 離散化ギャップ; 下限 90 倍・
  上限 11 倍のマージンで帯域拘束。経路統合が起きると下限が発火)。
- 重み付き対称性 1.9e-16、軸行の基底ベクトルプローブは厳密に 8/Δr²。
- 発見 2 点を文書化: (1) \(V_0=\Delta r^4/64\) のため重み付き CG 残差は
  軸セルをほぼ無視 → 非重み付き代数残差も併record・gate;
  (2) 大域誤差ゲートは軸係数 12/Δr² への破壊に盲目(誤差変化 0.4%)
  → 直接プローブが必須。
- 範囲: C は A と radial 数学を共有するため A-vs-C は axial 経路のみを
  検証。radial 独立性は B の役割。grid class と binary64 は 3 経路共有。

### 時間刻み収束(タスク 3、`outputs/hou_time_refinement_v1`)

同一格子 65×128・同一終了時刻 \(T_1\)、固定
\(\Delta t=6\times10^{-7},3\times10^{-7},1.5\times10^{-7}\)
(step 数 3653/7306/14612、実効刻みは厳密に半減)。
全 9 受入検査合格、manifest+全 14 payload SHA-256 検証済み。

| 量 | 観測時間次数 | step-doubling 差(粗→中/中→細) |
|---|---:|---|
| 増幅率 | **1.998** | 2.67e-5 / 6.69e-6 |
| energy | **2.034** | 1.23e-5 / 3.01e-6 |
| \(\max|\omega_1|\) | **1.999** | 0.197 / 0.0493(値 5.29e5 に対し相対 ~1e-7) |
| \(\max|u_1|\) | ≈2.0 | 7.55e-3 / 1.89e-3 |

- argmax 位置(u1・Cartesian 渦度)は 3 水準で**格子点単位で同一**。
- **時間/空間分離**: 拘束系指標(相対発散 4.245e-2、相対軸パリティ
  u1 0.706、循環違反 7.664e-4、solver-B cross-check 1.193e-2)は
  3 水準で相対スプレッド <2e-5 の**dt 非依存** — すべて空間離散化
  支配であることを直接確認。増幅率の時間誤差(最細 ~6.7e-6)対
  空間ラダー差(1.63〜6.58)の比は \(10^{5{-}6}\) で
  `dominant_error_source = "spatial"`。
- 固定 dt の実効 CFL 最大 0.0149(適応拘束の余裕内)。
- **注記**: \(T_1\) での相対軸パリティ 0.706(u1、65×128)は軸近傍
  front(\(r\approx0.047\)= 3 セル)の空間的未解像を相対指標が正しく
  可視化したもの(dt 非依存が示すとおり時間積分の破綻ではない)。
  Cartesian 監査の pointwise 所見と整合し、65×128 は \(T_1\) 近傍で
  軸近傍が深刻に解像不足であることを定量化する。

### 圧力非依存の primitive 残差(ギャップ解消、`outputs/hou_primitive_residual_v1`)

Cartesian 監査で「圧力未保存・snapshot 間隔 ≈500 step」として明示していた
ギャップを、**圧力を発明せずに**閉じた。運動量残差
\(R=u_t+(u\cdot\nabla)u-\nu\Delta u\) は厳密解で \(-\nabla p\) に等しいので
\(\operatorname{curl}R=0\)、同値に渦度輸送残差 \(S\) が消える。両者は
圧力を含まない。

- 実験は 7 状態(offset \(0,\pm1,\pm2,\pm4\) step)を保存するため、中心差分
  \(u_t\) の Richardson **次数が測定できる**(誤差推定ではない):
  実測 **1.99997**。
- 評価は保存 artifact を `load_candidate` で読み直す経路のみ(非循環性)、
  stencil は `cartesian_validation` 固有、円柱 module 不 import を AST 固定。
- 129×256 実測(内部領域、分母は**同一領域**の項和最大 — 甘くない正規化):
  相対 RMS \(\operatorname{curl}R=5.835\times10^{-4}\)、
  \(S=9.346\times10^{-4}\)、\(\nabla\cdot u=2.815\times10^{-4}\)。
  監査格子細分の次数 1.885 / 1.835 / 2.279。
- 全 10 受入検査+5 record-only 合格、manifest+26 payload SHA-256 検証済み。
- **正直に記録した所見**: (1) \(\nu=5\times10^{-4}\) では粘性項が項和の
  \(2.9\times10^{-4}\) しかなく、実データ上で粘性符号反転は**検出できない**
  (その旨を主張するテストを置き、決定的な注入は manufactured 場で実施);
  (2) 減衰 Taylor–Green は移流が厳密勾配のため非線形項に対して退化 →
  別の非自明 solenoidal 場を追加; (3) 2 形式は実データで相対 RMS 2.4% 相違、
  次数 1.28(バンドル中最弱); (4) 65×128 は次数 0.79/1.13 で bilinear 復元
  律速のため 129×256 を出荷; (5) 既存 `hou_early_time` checkpoint は
  この観点では**未監査のまま**。
- 圧力回復は scoped・record-only(転置を随伴恒等式 \(3\times10^{-15}\) で
  検証、零空間を射影除去。\(x,y\) 端の閉じ方により診断であって検証済み
  圧力ではない)。

### Lean 4 段階 1 — F-3(`formal/NSSingularity/VelocityRecovery.lean`)

E-14 ⇒ E-15(速度回復から物理発散ゼロ)を機械検証した。

- `divergence_of_recovered_velocity_eq_zero`(\(r\neq0\))と
  `divergence_of_recovered_velocity_eq_zero'`(\(u^r/r\) を連続延長に
  置換、軸込み)、`mixed_partial_comm`(Schwarz)。
- `sorry`・`admit`・新規 axiom なし。`#print axioms` は
  `[propext, Classical.choice, Quot.sound]` のみ(独立に再確認)。
- 独立検証: `lake build` 成功(8658 jobs)、軸の向き(`partialR` が
  \(r\)、`partialZ` が \(z\))を別途 Lean で確認 — 軸取り違えで別命題に
  なる罠を排除。
- 形式化により、E-15 の相殺が**混合偏微分の一致のみ**に依存することが
  構造的に露出した(滑らかさ仮定の唯一の使用箇所)。
- **非スコープ**: 任意の \(C^2\) スカラーに対する座標表現の恒等式であり、
  E-13 も運動方程式も使わない。E-18/E-24(3D 場との対応)は未形式化で、
  `ClayStatement.lean` の `DivergenceFree` へは未接続。

### 壁依存性(E-32 + 実装完了、実行中)

- **E-32**(`docs/equation_audit.md`): \(C^\infty\) compact-support
  envelope 初期値族。core(\(r\le0.9\))は E-29 と **bit 一致**、
  sup 偏差 \(\le3.4008\times10^{-10}\)(実測 \(2.577\times10^{-12}\))、
  \(r\ge0.95\) で厳密 0、遷移帯の 4 階差分は平の E-29 と bit 一致。
- 実装(`wall_dependence.py` + 実験 + 60 テスト)。`nonlinear_cylinder` が
  \(r_{\max}\) について既に一般であることを検証(4 半径で壁行の挙動が正しく、
  初期エネルギーは compact support ゆえ全半径同一 4008.106)。
- **前登録の完全性**: 閾値は一つも変更していない。受入検査 3 件の文言が
  一意でなかったため、実装した読み方を
  `docs/wall_dependence_prereg.md` §8 に**改版として明記**し、
  前登録されていなかった実装判断(整合性許容幅 0.25 等)も別枠で列挙した。

### 壁依存性の実行結果(`outputs/wall_dependence_v1`)

6 member(主群 \(\Delta r=1/192\) の \(R_{\mathrm{wall}}=1,1.5,2,3\)、
粗群 \(\Delta r=1/128\) の \(R=1,2\))、\(\nu=5\times10^{-4}\)、
\(t\in[0,T_1]\)。全 12 受入検査合格、manifest+97 payload SHA-256 検証済み。

| \(R_{\mathrm{wall}}\) | \(n_r\) | core 増幅率 \(A(T_1)\) | \(R=3\) との差 |
|---:|---:|---:|---:|
| 1.0 | 193 | 15.627954940635 | 2.4915e-3 |
| 1.5 | 289 | 15.630441984776 | 4.4673e-6 |
| 2.0 | 385 | 15.630446443875 | 8.1784e-9 |
| 3.0 | 577 | 15.630446452053 | — |

- 隣接対の分離 \(S\): **1.594e-4 → 2.853e-7 → 5.232e-10**、いずれも前登録
  閾値 0.05 を大きく下回る。argmax 変位は**厳密に 0**(全半径で同一格子点)。
- **§4 literal の判定は `wall_effect_small`**。ただし保守的な §2 読み
  (`classification_with_section_2_hold`)は **`undecided`**(下記の
  解像度整合性が不安定なため)。
- 楕円非局所寄与(core 上の \(\max|\Delta\psi_1|\))も同じ比で減衰:
  相対 2.968e-4 → 5.321e-7 → 9.760e-10。argmax は **\(r=0.8958\)**、
  すなわち core 外縁 \(r\le0.9\) の際 — 壁の影響が core 内で最大になる
  場所として整合的。
- **交差検証**: \(R=1\) の envelope member の増幅率
  15.627954940635405 は、平の E-29 初期値による `hou_early_time_v1` の
  193×384 の値と**完全一致**。E-32 性質 3・4(core bit 一致・ノルム不変)が
  力学レベルでも成り立つことの確認。

**機構の独立検証(本オーケストレータによる)。** 上表の差は
\(\Delta R=0.5\) ごとに **557.7 倍 / 546.2 倍**で減衰する。これは
\(z\)-Fourier モード \(k\) における \(\mathcal L_5\) の径方向解
\(r^{-1}I_1(kr),\ r^{-1}K_1(kr)\) に対する Dirichlet 壁の像応答
\(K_1(kR)/I_1(kR)\sim\pi e^{-2kR}(1+\tfrac3{4kR})\) が予言する比である。
\(z\) 周期 1 かつ \(\omega_1\) が \(z\) 奇なので \(k=0\) は存在せず最低
モードは \(k=2\pi\)。予言値は **555.2 / 545.5**(先頭項のみなら 535.5)で、
実測との差は **0.45% / 0.12%**。壁効果の正体は最低軸方向モードの
指数的像応答であると同定した。

**この結果が意味しないこと(最重要)。**

1. **全空間の壁独立性を示していない。** 指数減衰は \(z\) 周期 1 が
   \(k\ge2\pi>0\) を強制することの帰結である。\(\mathbb R^3\) では軸方向
   波数は \(k\to0\) まで連続で、楕円減衰は指数的ではなく代数的になる。
   本実験は「軸周期を固定したまま半径を後退させた」測定であり、
   前登録 §6 が明記するとおり壁を無限遠へ送る極限ではない。
2. **\(k=0\) モードの不在は Hou 初期値の奇対称性に依存する。** 一般の
   データでは代数減衰する単極子成分が残り、この指数的鈍感さは失われる。
3. **`wall_effect_small` は Clay 候補の証拠ではない。** summary の
   `interpretation` にも同文を記録した。増幅率自体は依然として解像度
   未収束(6.11→12.70→15.63→17.26)であり、壁を退けても
   軸近傍の未解像(相対軸パリティ 0.706 at 65×128)は解消しない。
4. **解像度整合性は不安定。** 共有半径対 \((1,2)\) の \(S\) は
   \(\Delta r=1/192\) で 1.594e-4、\(\Delta r=1/128\) で 2.504e-4 と
   **57.1% 変化**し、(前登録外の)許容幅 0.25 を超える。両値とも閾値
   0.05 を遥かに下回るため「既に無視できる量の 3 桁目が未収束」という
   状況だが、保守的読みでは判定を保留する。
5. 早期区間 \([0,T_1]\) のみ。粘性切替後・中成長段の壁依存性は未測定。

### E-33: 壁打切り応答の閉形式(`outputs/wall_truncation_scaling_v1`)

壁依存性の**結果を説明する**ために導出・監査・検証した。
\(\psi=\varphi/r\) の代入で \(\mathcal L_5\) は各 \(z\)-Fourier モードで
厳密に 1 位の変形 Bessel 方程式になり、compact support の \(\omega_1\) に
対する壁誘起誤差は core 上で厳密に \(K_1(kR)/I_1(kR)\) に比例する。

- **二つの漸近域**: \(kR\gg1\) で指数 \(\pi e^{-2kR}\)、\(kR\ll1\) で
  **代数 \(R^{-2}\)**、交差は \(kR\sim1\)。
- **\(k=0\) の厳密閉形式**: 壁依存部分は core 上で定数
  \(-Q_\infty/(2R^2)\)、\(Q_\infty=\int_0^\infty s^3\omega_1\,ds\)。
- 独立 oracle `bessel_reference.py`(scipy 非依存、級数 \(I_\nu\)+積分表示
  \(K_\nu\)、Wronskian 恒等式を \(6.7\times10^{-16}\) で検証、求積は
  spectral)に対し、全 23 受入検査合格。測定と oracle の一致は最細格子で
  \(4.2\times10^{-4}\)(次数 1.997–2.002)、形状の \(I_1(kr)/r\) 相関ずれ
  \(5.5\times10^{-8}\)(次数 3.98–4.02)、\(k=0\) 閉形式は次数 2.0000。
- **非線形との連結**: 壁依存性実行の分離比 557.746 / 545.231 は楕円 oracle
  と \(1.2\times10^{-4}\) / \(2.8\times10^{-4}\) で一致。**非線形力学の壁
  感受性は線形楕円 Bessel 機構そのもの**である。
- 記録した整合性訂正: 初期の探索値は参照壁の規約を明記していなかった。
  検証は参照 \(R=3\) の規約で厳密に再現し、最大半径参照の正しい値
  (5.939/5.306、3.242/2.711、4.5006)を E-33 に記載した。
- 丸め床の扱いを first-class 化(応答が \(10^{-19}\) を切る半径は oracle
  から決定論的に分類し、測定値で比を取らない)。

**プログラム上の帰結。** Hou 設定は \(kR=2\pi\approx6.28\) と指数域の奥に
あるため壁が効かない。しかし \(\mathbb R^3\) では \(k\to0\) まで連続で、
長波長成分の打切り誤差は \(R^{-2}\) でしか減らない
(\(\varepsilon=10^{-6}\) に core 半径の \(10^3\) 倍)。**大半径 Dirichlet
円柱は全空間近似として代数的にしか収束せず実用的でない。**
同じ modal 構造から**厳密な透過条件**
\(\partial_r\hat\psi_k(R)+[2/R+kK_0(kR)/K_1(kR)]\hat\psi_k(R)=0\) が導け、
compact support の下で打切り誤差が厳密に消える。設計と受入条件を
`docs/whole_space_transition.md` に実装前固定した。

### W-A: 透過境界条件ソルバ(`outputs/transparent_boundary_v1`)

実装 `transparent_boundary.py`、全 25 受入検査合格(3.5 秒)、
W-B の受入条件 6 件すべて合格。既定は全経路で Dirichlet のまま。

- **\(k=0\)**: core の \(R=1\to2\) 差は透過で \(9.07\text{e-}7\to
  2.27\text{e-}7\to5.66\text{e-}8\)(次数 2.0002)。Dirichlet は
  \(-1.697\text{e-}3\) で不変。改善率 1872 → 7489 → **29959**。
  E-33e の offset が消え、2 次で 0 へ収束する量に置き換わった。
- **\(k>0\)**: 改善率は \(L_z=1\) で ×11013、\(L_z=8\) で ×31787、
  長波長 \(L_z=32\)(代数域)で ×30510、最大 ×34231。全 24 の
  \(R\) 非依存性の観測次数は [1.948, 2.001]。
- manufactured 収束は全行 1.9975–2.0038、**境界行のみ 2.0099–2.0150**。
  境界行の局所打切りは \(O(\Delta r)\) だが M-matrix 性により全体 2 次
  (仮定ではなく測定)。
- 故障注入: 符号反転 3628–4152 倍、\(2/R\) 脱落はソルバ拒否、
  \(K_0/K_1\!:=\!1\) は 22.5–186 倍。support 漏れは 2 経路で拒否。
- 正直に記録した限界: `frozen_ratio` は \(k=0\) で bit 単位不可視、
  条件 4 の \(L_z=32\) は比較対象の大半径 Dirichlet 自体が未収束
  (E-33 の代数 \(R^{-2}\) そのもの)、\(K_\nu\) 求積は \(x>10^3\) で
  文書化した漸近形へ切替。

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

## 2026-07-28 第 2 セッション: FABLE5_NEXT_TASK_AUDIT の P0/P1 ゲート

指示書 `FABLE5_NEXT_TASK_AUDIT.md`(リポジトリルート、コミット済み)の
記載順に実施した。開始時の確認: 作業ツリーはリモートと一致(2cb8e48)、
全 565 テスト合格、`lake build` 成功(8659 jobs)、`formal` 内に
sorry/admit/新規 axiom ゼロ。

### P0-A: von Neumann 安定性監査(`von_neumann.py`、11 tests)

- Heun+中心差分は純移流で厳密増幅(\(|G(i\alpha)|^2=1+\alpha^4/4\))。
  凍結係数 advection–diffusion の全波数 scan、予測子段の別評価、
  snapshot 監査 API、独立参照 propagator(シンボル経路と実配列経路の
  相互一致 3.7e-16)を実装した。
- **出荷済み運転点の判定は「stability-unverified」**: 出荷 v1 の記録
  (min dt=2.76e-7、max advective CFL=0.10023)と自己整合な読みで
  Heun worst \(|G|\) は 1.0000035(radial 支配)/ 1.0(axial 支配)/
  1.000152(両方向同時)— tolerance 1e-12 では不合格。全波数 pass に
  必要な dt は約 2.4e-8(出荷値の 1/11.5)。粘性 5e-4 の寄与は
  増幅率 −6.1e-5 に留まる。詳細は `docs/numerical_stability_audit.md`。
- 過去の Heun 実行はすべて stability-unverified に再分類。**Heun 単独の
  増幅を候補判定に使うことを禁止**(決定規則として文書化)。

### P0-A: 交差検証積分器(SSPRK3 / RK4)

`take_step` に SSPRK3・古典 RK4 を追加(空間離散化・拘束順序・楕円
solve は Heun と完全共有)。毎段 `constrain_state` 射影込みの実測時間
次数: **Heun 1.97/2.00、SSPRK3 3.00/3.00、RK4 3.95/3.98** — 射影は
観測次数を落とさなかった。ゼロ場不動点・小 dt 相互一致・誤差の大小
関係(rk4 < ssprk3 < heun)をテストで固定。

### P0-B/P0-C: 全 step streaming gate(`test_integrator_gates.py`、20 tests)

- `IntegrationResult` に `step_stream`(全 accepted step × 28 量)と
  `gate_summary`(streaming 極値)を追加。エネルギー増分、循環 defect、
  奇対称相対、軸 parity 相対、相対発散、壁拘束、**solver A が実際に
  解いた線形系の代数残差**(全 step 相対 <1e-12)、pre/predictor/post
  CFL、粘性安定数、dt を縛った拘束の名前、エネルギー収支 defect を
  全 step 保存。出力間引きは gate に作用しない。
- 中間段 CFL 超過での step 棄却(`stage_cfl_limit`、dt 半減再試行)を
  実装、受入 step の段 CFL が閾値以下であることをテストで固定。
- **義務の合成違反テスト**: 記録行間だけの強制パルス注入→抽出で、
  間引き history は単調減衰のみを示すが streaming は丸め床の
  3×10¹² 倍の増加を捕捉する。

### P1-C: エネルギー収支と viscosity_sign fault

- E-27 壁は swirl のみ no-slip(\(u^z(1,z)\ne0\) の滑り壁)なので、
  正しい恒等式は \(dE/dt=-\nu\int|\omega|^2dV-\nu\oint u^z\omega^\theta dS\)。
  壁項込み/なしの両 defect と swirl エネルギーの項別仕事率
  (移流・stretching・粘性)を全 step 記録。
- **整合性の記録**: 初回実装は \(\int|\omega|^2dV\)(\(2\pi\) 測度)を
  E-20 の \(\pi\) 正規化 enstrophy と取り違え、相対 defect が理論値
  どおり厳密に 0.5 へ飽和(fault 時 1.5)。これを因子 2 の欠陥として
  特定・修正した。修正後: 滑らか control で 5.9e-2 → 1.6e-2 → 4.2e-3
  (空間時間同時細分で収束)。
- 新 fault `viscosity_sign`: Hou 運転点(ν=5e-4)では項和の 3e-4 で
  不可視(既知)だが、拡散支配 control(ν=2e-2)で相対 defect
  clean 2.21e-2 vs 反転 **2.000**(理論値 2、比 90 倍)、エネルギー
  単調性反転により確実に棄却される。

### P0-D: core-width / points-per-scale(`core_width.py`、15 tests)

- radial/axial FWHM、10–90% front thickness、subgrid 二次ピーク、
  勾配長スケール、高周波 tail(z rfft / r DCT-II)、共通格子への
  Catmull-Rom 補間比較、manufactured tanh front 研究を実装。
- **前登録閾値 `PREREGISTERED_MIN_POINTS_PER_FRONT = 7`**(tanh front
  で 10–90 幅の相対誤差 ≤2% となる最小整数点数。テストが研究を再計算し
  定数との整合を強制するため事後調整不可)。
- 出荷済み 193×384 の T₁ snapshot: points_per_fwhm_r 6.92、
  **points_per_front 4.36**、ピークは軸から 6 セル、勾配スケールは
  0.49dr / 0.98dz(1 セル未満)→ `fit_precondition` は**不合格**
  (3 理由)。「grid-scale 飽和なし」という旧表現は撤回した(本文中に
  改版注記)。**収束 fit は現データでは禁止**が機械化された。

### P1-A: blind 外挿(`extrapolation.py`、10 tests)

- 3 点厳密解・最小二乗 power law・固定次数 Richardson・全部分列感度。
  署名は \((h, A)\) のみで外部 anchor を構造的に受け付けない
  (`20.5235`・`20.52` の不在をテストが assert)。
- 実ラダー (6.11, 12.70, 15.63, 17.26): 指示書の見積(全4点 27.38/0.54、
  前 3 点 28.85/0.49、後 3 点 24.60/0.70)を独立実装+brute force で
  すべて再現。A_inf 散らばり相対 **0.488**(前登録閾値 0.05 の 9.8 倍)、
  p 幅 0.203 → **判定 not_in_asymptotic_range: 極限値は一切引用不可**。
- 記録した偶然: 次数 1 固定 Richardson は 20.528(公表 20.5235 の
  0.02% 以内)を与えるが、これは相互に矛盾する 7 外挿値の一つに
  すぎず(次数 2 固定は 16.844)、確認としては読めない。テスト
  docstring に非結果として明記。

### Poisson 第三経路の SPOF 破壊(P0 §5、`test_realspace_poisson.py` 追加)

solver A/B が共有する Fourier 機構の故障 3 種(モード正規化 slip、
周期 seam ずれ、Nyquist 混入)を solver A 解の事後改変として合成し、
実空間第三経路 C が **>10×**(実測 48× 等)で全て検出することを固定。
A/B の一致を「continuum 精度」「完全独立」と呼ばない limitation は従来
どおり。

### P0-E: 語彙修正と Gate 4 仕様

- `docs/whole_space_transition.md` §0 を新設: W-A と壁依存性の全結果を
  「**periodic-z radial-wall sensitivity observation**」と固定し、
  「whole-space validation」「R³ wall independence」の表現を禁止。
- 同 §7 に真の全空間移行 gate(非周期 z、z 方向 C∞ compact 台、
  free-space 楕円経路、R_max/Z_max 独立拡大、低波数 stress、
  Cartesian 体積測度での有限エネルギー直接検査、有限円柱解との
  同一視禁止)を**未実装の仕様**として定義。
- STATUS/設計文書の「grid-scale 飽和なし」表現を撤回(改版注記付き)。

### Gate 1 実行結果(`outputs/integrator_comparison_v1`、全 8 受入合格)

65×128、E-29 datum(振幅 12000)、\(\nu=5\times10^{-4}\)、\(T_1\) まで、
固定 dt ∈ {6e-7, 3e-7} × {heun, ssprk3, rk4}(空間離散化・拘束順序・
楕円 solve は完全共有)。前登録許容: 増幅相対差 1e-3、対差の dt 縮小、
argmax 1 セル以内。

| dt | 対 | 増幅相対差 | 場 L∞ 相対差 | argmax 差 |
|---|---|---:|---:|---|
| 6e-7 | heun vs ssprk3 | 5.823e-6 | 4.39e-6 | 0 セル |
| 6e-7 | heun vs rk4 | 5.828e-6 | 4.39e-6 | 0 セル |
| 6e-7 | ssprk3 vs rk4 | **5.09e-9** | 4.40e-9 | 0 セル |
| 3e-7 | heun vs ssprk3 | 1.458e-6 | 1.10e-6 | 0 セル |
| 3e-7 | heun vs rk4 | 1.459e-6 | 1.10e-6 | 0 セル |
| 3e-7 | ssprk3 vs rk4 | **6.35e-10** | 5.50e-10 | 0 セル |

- heun と高次法の差は dt 半減で正確に 1/4(heun の \(O(dt^2)\) 誤差
  そのもの)。ssprk3 と rk4 は互いに 5e-9 で一致。**von Neumann worst-case
  上界(累積 ~0.8%)は実測では発現せず、時間積分スキーム依存は
  ~6 ppm** — 現在の増幅値に対する時間離散化リスクの実測上界。
  (これは凍結係数監査の「stability-unverified」分類を置き換えるもの
  ではなく、補完する経験的上界である。固定 dt=6e-7 の段 CFL は最大
  0.0149 で、適応 run の 0.1 よりはるかに安定側にある点にも注意。)
- heun dt=6e-7 の A_grid=6.1147053247 は `hou_time_refinement_v1` の
  粗 dt 値と**ビット単位一致**(実験間決定論再現)。
- **P1-B 両正規化の併記**: A_grid(離散初期最大で正規化)6.1147 に対し
  A_common(連続参照 \(24000\pi/\sqrt{37}(36/37)^{18}=7569.6227\) で
  正規化)**6.0539** — 65×128 の離散初期最大は連続値より約 1% 低い。
  分母の格子依存だけで増幅率が 1% 動くことが定量化された。
- **正直な記録**: 全 step streaming の相対エネルギー収支 defect 最大は
  **0.936**(全 6 run とも)。滑らか control では 4.2e-3 まで収束する
  計装なので、これは E-29 front(65×128 で FWHM 約 4 点)の空間
  未解像がエネルギー恒等式を閉じさせないという**解像度の言明**である。
- エネルギーは**全 accepted step で単調非増加**(max step increase < 0。
  旧「snapshot 対でゼロ」より強い)。循環 defect 7.67e-4(前登録閾値
  1e-3 以内)、Poisson 代数残差は全 step 相対 5e-15 以下。

### 出荷済み証拠への von Neumann 監査適用(`outputs/von_neumann_audit_v1`)

出荷済み 3 bundle の全診断行(1470 行、stride 25、step-0 除外 8)を
`audit_snapshot` で監査(721² scan × 2940 回、51.8 秒)。

| run | 行数 | 不合格行 | worst Heun max\|G\| | 判定 |
|---|---:|---:|---|---|
| v1 65×128 | 88 | 0 | 1.0(厳密) | verified-at-recorded-rows |
| v1 129×256 | 88 | 0 | 1.0(厳密) | verified-at-recorded-rows |
| **v1 193×384** | 89 | **4** | **1.0000312** | **stability-unverified** |
| v2 129×256 | 88 | 0 | 1.0(厳密) | verified-at-recorded-rows |
| v2 257×512 | 92 | 0 | 1.0(厳密) | verified-at-recorded-rows |
| refinement dt=6e-7/3e-7/1.5e-7 | 147/293/585 | 0 | 1.0(厳密) | verified-at-recorded-rows |

- 不合格 4 行は 193×384 の \(T_1\) 直前(t≈2.12–2.19e-3、CFL_z≈0.1002)。
  strided 外挿値 1.0016 は「bound ではない」と明示ラベル付き。
- Euler 予測子段は全 run で 1 を超える(最大 1.0176)— 記録のみ。
  完成 step の Heun のみを gate する。
- 制約: stride 25 の記録行のみの被覆。行間の 24 step は未監査
  (将来 run は `step_stream` で全 step 被覆)。判定語彙は
  「stability-unverified であって不安定ではない」を全箇所で維持。

### 既存 snapshot の core-width / P1-B 再正規化(`outputs/core_width_audit_v1`)

全 4 解像度 × 5 時刻 × 2 場(u1、\(|\omega|\))の points-per-scale 監査
(入力 manifest 検証済み、v1/v2 の共有 129×256 は byte 一致を確認)。

- **\(T_1\) の fit 前提は全解像度・両場で不合格**: points_per_front は
  u1 で 2.57/3.36/4.36/5.43、\(|\omega|\) で 2.55/3.29/4.43/5.70
  (閾値 7)。**T₁ 増幅ラダーの収束 fit 禁止が機械的に確定**
  (`convergence_fit_precondition_satisfied_at_final_snapshot = False`)。
- 計算された正直な例外: t=5e-4 の \(|\omega|\) ラダーのみ全解像度で
  前提を満たす(front 7.06/14.0/21.0/28.0 点)。早期時刻の front は
  まだ広いという整合的な結果。
- **P1-B 表(離散初期最大 vs 連続参照 7569.6226982)**:

  | nr | 離散初期max | a/b | A_grid | A_common |
  |---|---:|---:|---:|---:|
  | 65 | 7494.31 | 0.9901 | 6.1148 | **6.0539** |
  | 129 | 7554.61 | 0.9980 | 12.6957 | **12.6705** |
  | 193 | 7561.48 | 0.9989 | 15.6280 | **15.6112** |
  | 257 | 7565.84 | 0.9995 | 17.2588 | **17.2502** |

  分母の格子依存だけで最大 1% 動く。以後の主比較は絶対値と
  A_common を用い、A_grid は補助値とする(P1-B)。
- 連続初期最大位置 \(r^*=1/\sqrt{37}\) を数値最大化で 1.2e-16 まで再現。
  初期最大位置の格子誤差は最大 0.48 セル。
- 隣接解像度の共通格子差(\(T_1\)、\(|\omega|\) L∞): 6.27e4 → 3.38e4 →
  1.97e4 と値は減少するが、**微分 L∞ は 4.5e6/5.9e6/5.0e6 と減少しない**
  — 微分レベルの収束は現ラダーに存在しない(正直に記録)。

### Lean 監査(P0 §7)

`formal/AxiomAudit.lean` を追加し `lake env lean AxiomAudit.lean` を実行:
**9 定理すべて `[propext, Classical.choice, Quot.sound]` のみ**に依存
(記録は `docs/formalization_map.md`)。sorry/admit/新規 axiom ゼロ、
toolchain/mathlib は v4.32.1 固定を再確認。「8659 jobs = 8659 定理」
という読みの禁止を明文化。

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
- (2026-07-28 追加)出荷済み run の von Neumann 監査は stride 25 の
  記録行のみを被覆する。行間の step は未監査であり、将来 run の
  `step_stream` 全 step 被覆でのみ閉じる。
- (2026-07-28 追加)193×384 run は記録行監査でも stability-unverified
  (T₁ 直前 4 行、max|G| ≤ 1.0000312)。65×128 の Gate 1 相互比較は
  時間スキーム依存 ~6 ppm を示したが、193×384 の適応 CFL 0.1 運転点
  そのものでの相互比較は未実施。
- (2026-07-28 追加)65×128 の E-29 run では全 step 相対エネルギー収支
  defect が 0.936 に達する(滑らか control では 4.2e-3 へ収束)。現行
  解像度では離散エネルギー恒等式が front 上で閉じない。
- (2026-07-28 追加)隣接解像度の共通格子差は値では減少するが微分 L∞
  では減少しない。T₁ の全解像度が fit 前提(front ≥ 7 点)を満たさず、
  増幅ラダーの外挿は blind 判定でも not_in_asymptotic_range。
- (2026-07-28 追加)`hou_early_time` / `wall_dependence` 実験本体は
  まだ旧計装(snapshot 系 gate)のまま。次回実行前に `gate_summary`
  読み出しへの移行が必要(新規 run のみ。既存証拠は不変)。
- (2026-07-28 追加)GitHub が既定 branch に dependabot 警告 4 件
  (moderate 2、low 2)を報告している。数値結果には影響しないが未対処。

## 次に行うべき最小の一手

〔改版 2026-07-28 第 2 セッション。FABLE5_NEXT_TASK_AUDIT の Gate 順序
(Gate 1–4 が通るまで中後期成長・blow-up fit・AI 候補探索へ進まない)を
最上位の拘束とする。Gate 1 は合格、Gate 2/3 は既存証拠+新監査で部分的、
Gate 4 は未実装。〕

1. **Gate 4 の実装(最優先)**: 非周期 \(z\) の有限 box、\(z\) 方向も
   \(C^\infty\) compact な初期値族、free-space 楕円経路(W-1 の \(z\)
   非周期版/Green 積分/Hankel)、\(R_{\max}\)/\(Z_{\max}\) 独立拡大、
   低波数 stress test(`docs/whole_space_transition.md` §7)。その入口
   として \(L_z\in\{1,2,4\}\) 族(既存実装の config 変更で測定可能)で
   指数→代数遷移を実測する。
2. **実験本体の新計装への移行**: `run_hou_early_time` /
   `run_wall_dependence` の受入検査を `gate_summary`(全 step streaming)
   読み出しへ切替え、`stage_cfl_limit` の使用を判断する。以後の新規 run
   はすべて全 step 被覆+3 積分器のうち 2 つ以上での交差確認を要件とする
   (Heun 単独増幅の候補判定使用は禁止済み)。
3. **軸近傍解像度の設計判断**: T₁ の fit 前提不合格(front ≤ 5.7 点)、
   相対軸パリティ 0.706、微分レベルの共通格子差の非減少はいずれも
   front 未解像を指す。適応 mesh または半周期 sine 実装の設計判断を、
   中成長段の前に行う。
4. **既存 checkpoint の primitive 監査の完結**: snapshot 前後 1 step の
   追加保存 option(既存実験へ)。
5. **Lean 段階 1 の継続**: F-1(再スケーリング恒等式)と F-4(証明書
   不等式)。F-2/F-3 と同じ方針(定義明示、非スコープ明記、
   `#print axioms` 記録、`AxiomAudit.lean` へ追記)。

これらが通っても、長時間探索、AI最適化、特異点fitへ自動的には進まない。
動的再スケーリング探索の前に、全空間tailと候補用離散化の証明可能な設計を
再評価する(PO-05〜PO-07)。
