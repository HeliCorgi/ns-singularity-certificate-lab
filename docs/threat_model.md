# Threat model for false singularities

## 1. 守るべき主張

このプロジェクトが数値段階で許される主張は、「指定した離散化・領域・精度で、既知の障害と直ちに矛盾しない候補挙動を観測した」までである。次は数値結果だけから主張しない。

- 有限時間特異点が存在する。
- Clay Millennium Problem を解決した。
- 小さい collocation residual が連続 PDE の解の存在を示す。
- フィットした \(T\) が真の発散時刻である。

保護対象は、候補場、実験設定、乱数 seed、全診断時系列、失敗実験、コード版、依存関係、再現手順である。脅威は「誤って特異点らしく見えること」と「都合の悪い反証データが記録から落ちること」の両方を含む。

## 2. 判定の原則

### 2.1 三段階

1. **単一実行の健全性**: NaN、境界、軸条件、solver residual、保存データを検査する。
2. **精度系列の整合性**: 空間、時間、領域、精度、アルゴリズムを独立に変える。
3. **独立検証**: 同じ微分行列や同じ residual 関数を共有しない別実装で再計算する。

段階 1 を通っても候補ではない。段階 2 を通ったものを「数値候補」、段階 3 を通ったものだけを「証明書化候補」と呼ぶ。

### 2.2 事前登録する量

実験開始前に config へ次を固定する。

- 比較する解像度列 \((N_r,N_z,N_t)\)。
- 時間刻み列と CFL 上限。
- 領域列 \((R_{\max},Z_{\max})\)。
- 候補受入れに使う主要量と補助量。
- 比較ノルム、共通物理領域、時刻補間法。
- 収束と不合格の閾値。
- blow-up fit に使う時刻窓の候補集合。
- seed、optimizer、停止条件、全 run を保存する規則。

閾値を結果を見た後で動かした場合は、変更理由と旧結果を残し、新しい実験系列として扱う。

## 3. 必須の誤検出シナリオ

### TM-01 空間解像度不足

**発生機構。** 集中幅が格子幅に近づくと、補間ピーク、差分勾配、渦度最大値が格子依存で増大する。adaptive mapping が物理空間の不足を隠す場合もある。

**具体的検出テスト。**

1. 同じ物理初期値を少なくとも三つの入れ子解像度
   \(N,\,\lceil3N/2\rceil,\,2N\) で計算する。
2. 高解像度解を低解像度の**共通物理格子**へ制限し、
   \(L^2\)、重み付き軸対称 \(L^2(2\pi r\,dr\,dz)\)、\(L^\infty\)、勾配ノルムを比較する。
3. manufactured solution で既知の空間次数 \(p\) を測り、実流でも
   \[
   p_{\rm obs}=\log_2(E_N/E_{2N})
   \]
   が事前指定下限を保つか調べる。
4. 最小集中幅 \(\ell_r,\ell_z\) あたりの有効点数と mapping Jacobian の最小・最大を保存する。
5. スペクトル法では最後の 10–20% の mode energy、係数包絡線、解析性幅推定を保存する。

**不合格ゲート。** ピーク値または推定 \(T\) が解像度とともに単調に移動し、場の差が減少しない、集中幅が事前指定点数未満、または tail が丸め誤差に達する前に打切りへ衝突した場合。

### TM-02 時間刻み不足

**発生機構。** 急成長を跨ぐ大きな \(\Delta t\)、adaptive step の拒否漏れ、時間補間の overshoot が偽ピークや偽の有限 \(T\) を作る。

**具体的検出テスト。**

1. 同一空間格子で \(\Delta t,\Delta t/2,\Delta t/4\) を比較する。
2. 固定刻み法と、埋め込み誤差推定を持つ adaptive 法を比較する。
3. 全 accepted/rejected step、局所誤差推定、CFL、拡散数
   \(\nu\Delta t(\Delta r^{-2}+\Delta z^{-2})\) を保存する。
4. 時刻補間を使わない output と dense output を比較し、最大値の位置・値を照合する。
5. 既知の滑らかな基準解で時間収束次数を独立に確認する。

**現在の前段証拠。** 滑らかな旋回拡散対照については、同一の513点空間格子で
\(\Delta t=0.5,0.25,0.125\) を比較し、解析解誤差とstep-doubling差が約2次で
減少することを確認した。これは上記1と5の実装試験に限られ、adaptive step、
非線形production solver、未知候補の時間誤差を検査したものではない。

**不合格ゲート。** \(\Delta t\) 半減で主要診断の差が所定比で減らない、step rejection が集中直前に連続する、または推定 \(T\) が時間刻みに比例して移動する場合。

### TM-03 エイリアシング

**発生機構。** 二次非線形項の unresolved convolution が低波数へ折り返し、偽のエネルギー移送、渦度増幅、または小 residual を生む。

**具体的検出テスト。**

1. Fourier 方向では \(3/2\)-padding または同等の exact dealiasing を既定とし、未処理版と比較する。
2. Chebyshev/Jacobi 方向では overintegration、padding、または係数空間 convolution を比較する。
3. 非線形項を物理空間積と係数空間 convolution の二経路で計算し、十分低い mode で一致を確認する。
4. inviscid/nonlinear substep の離散エネルギー収支
   \(\langle u,(u\cdot\nabla)u\rangle\approx0\) を検査する。
5. 位相シフト dealiasing または 2 倍 oversampling を spot check として使う。

**不合格ゲート。** dealiasing の有無で成長率・スペクトル flux・推定 \(T\) が許容差を超えて変わる、または非線形エネルギー収支誤差が離散化誤差より大きい場合。

### TM-04 スペクトルブロッキング

**発生機構。** 散逸が cutoff 前で不足すると、最高 mode 付近にエネルギーが堆積し、Gibbs 振動や “tyger” 状構造が物理的集中に見える。

**具体的検出テスト。**

1. shell/mode ごとのエネルギーと enstrophy、最高 10% mode の割合を時系列保存する。
2. tail の対数係数包絡線が単調減少するか、cutoff 前で平坦化・上昇しないか検査する。
3. 解像度を上げたとき、pile-up の波数が cutoff と共に移動するか、物理波数に固定されるかを見る。
4. ごく弱い spectral viscosity/filter を複数強度で加え、低 mode と候補時刻の感度を測る。filter 付き結果だけを採用しない。
5. 実空間の高周波振動と局所解析性幅を照合する。

**不合格ゲート。** ピーク成長が tail pile-up と同時に始まる、振動の波長が常に格子幅に比例する、または弱い filter で候補が消える場合。

### TM-05 有限領域境界の反射・汚染

**発生機構。** \(\mathbb R^3\) を有限円柱で切ると、Dirichlet/Neumann/人工境界が圧力・楕円速度回復を瞬時に内部へ伝え、放物拡散も境界誤差を運ぶ。

**具体的検出テスト。**

1. 同じ内部初期値で \((R_{\max},Z_{\max})\)、\(1.5\) 倍、\(2\) 倍の領域を比較する。
2. 候補核から境界までの距離を集中幅で割った
   \(d_{\partial\Omega}/\max(\ell_r,\ell_z)\) を保存する。
3. 境界上の速度、渦度、圧力、法線エネルギー flux、楕円 Green tail の推定を保存する。
4. 少なくとも二種類の遠方条件（解析 tail matching、同次条件など）で共通内部領域を比較する。
5. 初期データを cutoff する半径も変え、cutoff 層の渦度が核へ影響しないか確認する。

**不合格ゲート。** 領域拡大または境界条件変更で核の成長率・位置・推定 \(T\) が収束しない、境界 flux が内部散逸と同程度、または楕円解の差が核で支配的な場合。

### TM-06 不正確な軸境界条件

**発生機構。** \(1/r\) 項を安易に評価したり、偶奇性を破ると、最初の格子点に \(O(1/r)\) の偽ソースが生じる。軸上の非零 \(u^\theta\) は見かけ上の \(u^\theta/r\) 発散を作る。

**具体的検出テスト。**

1. 各変数の軸 Taylor 展開と偶奇性を明示し、ghost/基底が構造を厳密に満たすかテストする。
2. \(u^\theta(0,z)=u^r(0,z)=0\)、\(\partial_r u^z(0,z)=0\)、\(\Gamma=ru^\theta=O(r^2)\) を離散的に検査する。
3. \(u_1,\omega_1,\psi_1\) の \(\partial_r=0\) と、Cartesian 再構成の滑らかさを検査する。
4. 軸を横切る Cartesian manufactured field から円柱成分を作り、軸近傍の収束次数を測る。
5. 軸点を極限公式で扱う実装と、偶拡張した独立実装を比較する。
6. PDE residual を軸からの格子層別に集計し、第一層だけが支配していないか見る。

**不合格ゲート。** parity defect が収束しない、軸近傍 residual が全体を支配する、または Cartesian 再構成が方向依存になる場合。

### TM-07 stiffness

**発生機構。** 小さい粘性長、adaptive mapping、拡散作用素の大きな固有値、非線形増幅が、陽解法の安定領域や nonlinear solve の許容誤差を超える。数値不安定が物理成長に見える。

**具体的検出テスト。**

1. explicit、IMEX、完全 implicit のうち少なくとも二方式で短区間を重複計算する。
2. 線形化作用素の最大固有値または Gershgorin/FFT bound から安定刻みを見積もる。
3. implicit solve の反復 residual と update norm を各 step 保存する。
4. 時間刻みを安定限界の一定割合以下へ下げ、増幅率が変わらないか確認する。
5. 既知の stiff manufactured problem で order reduction を測る。

**不合格ゲート。** 主要成長が安定限界を越えた時刻から始まる、solver residual が PDE residual より大きい、または方式変更で候補が消える場合。

### TM-08 浮動小数点のオーバーフロー (overflow) / underflow / 非有限値

**発生機構。** 大振幅、微小スケール、座標 Jacobian、指数的再スケーリングが表現範囲を超える。`inf` の直前値を物理的発散と誤認しうる。

**具体的検出テスト。**

1. すべての state、derivative、Jacobian、diagnostic に `isfinite` assertion を置く。
2. 最大指数、最小正規数、subnormal 数、再スケーリング係数の log を保存する。
3. 値そのものと \(\log |f|\) を併用し、積・商は必要に応じて log-domain で計算する。
4. 候補末端の checkpoint を `float64`、`longdouble`（実装が真に拡張精度を提供する場合）、任意精度評価で再診断する。
5. 非有限値を含む run は明示的失敗として保存し、成功 run に分類しない。

**不合格ゲート。** 非有限値、subnormal の大量発生、指数余裕の不足、または高精度再評価で主要桁が一致しない場合。

### TM-09 桁落ち・相殺

**発生機構。** 大きな移流、伸長、粘性項の差として小 residual を計算すると、有効桁が失われる。小 residual が「方程式をよく満たす」ことを意味しない。

**具体的検出テスト。**

1. 各 PDE 項を別々に保存し、
   \[
   \kappa_{\rm res}
   =\frac{\sum_j|F_j|}{|\sum_jF_j|+\epsilon_{\rm scale}}
   \]
   を点wise・ノルムで保存する。
2. compensated/pairwise summation と通常加算を比較する。
3. checkpoint を任意精度で読み、微分は独立な高精度係数評価で再計算する。
4. 有次元値と無次元化値の両方で residual を計算する。
5. 差分による時間微分と、PDE右辺からの時間微分を独立に比較する。

**不合格ゲート。** residual が丸め誤差推定以下なのに \(\kappa_{\rm res}\) が精度の逆数に近い、加算順序で符号が変わる、または高精度評価で residual が大幅に増える場合。

### TM-10 発散時刻 \(T\) の過剰フィッティング

**発生機構。** 短い時系列に \(A(T-t)^{-\alpha}\) を当てると、滑らかな急増や指数成長にも有限 \(T\) を割り当てられる。\(T,\alpha,A\) の相関が強い。

**具体的検出テスト。**

1. fit window の開始・終了を格子状に変え、\(T,\alpha\) の profile likelihood と covariance を保存する。
2. 末尾 20–30% を holdout とし、fit に使わず予測誤差を評価する。
3. power law、exponential、double exponential、有限値への飽和、log-corrected law を情報量基準と holdout error で比較する。
4. \(\|\omega\|_\infty\)、\(\|\nabla u\|_\infty\)、集中幅、再スケール時間、複数 Serrin 量から得た \(T\) を比較する。
5. 解像度・時間刻み・領域ごとに別々に fit し、最後に \(T\) の収束を検査する。
6. 推定 \(T\) を固定して fit し直した残差と、自由 \(T\) の改善量を報告する。

**不合格ゲート。** \(T\) が window 終端に追随する、holdout を予測できない、別診断の信頼区間が重ならない、または代替モデルが同等以上に説明する場合。

### TM-11 一時的な高勾配

**発生機構。** vortex stretching と粘性の競合で大きな一過性増幅が起きても、その後に飽和・拡散しうる。短時間だけを見ると発散に見える。

**具体的検出テスト。**

1. 安定に追跡できる限り、ピーク後も複数の局所 turnover time / diffusive time まで計算する。
2. 振幅だけでなく集中幅、スペクトル tail、散逸率、BKM 累積量、Serrin 累積量を同時追跡する。
3. 動的再スケール時間 \(\tau\) で profile が固定・周期・準定常になるか、元へ戻るかを見る。
4. 初期値の小摂動でピーク時刻・増幅率が連続に変わるか調べる。
5. 同じ診断を既知の非特異な高勾配基準流へ適用し、分類器の偽陽性率を測る。

**不合格ゲート。** 幅が下げ止まる、ピークが解像されたまま減衰へ転じる、または singular fit が観測窓延長で崩れる場合。

### TM-12 ニューラルネットのスペクトルバイアス

**発生機構。** ネットワークは低周波を先に学習しやすく、薄い層・高周波 tail・軸境界層を落とした滑らかな偽プロファイルでも平均 residual を小さくしうる。[Rahaman2019]

**具体的検出テスト。**

1. 学習後の場を oversampled 格子上で Fourier/Chebyshev/Jacobi 係数へ射影し、mode ごとの真の residual を測る。
2. collocation 点とは独立なランダム点、低 discrepancy 点、軸・ピーク・tail を重点化した adversarial 点で検証する。
3. 高周波 manufactured solution を同じ architecture/loss で学習し、回収可能帯域を測る。
4. ネットワーク出力を明示基底へ蒸留した後、係数を直接 Newton/least-squares で再最適化し、ネットワークなしで residual を再評価する。
5. architecture、activation、Fourier feature 帯域を変え、物理 profile の一致を確認する。

**不合格ゲート。** collocation residual だけが小さい、oversampling で residual が増える、基底射影の tail が未収束、または蒸留後に候補が維持できない場合。

### TM-13 オプティマイザー (optimizer) が作る偽解

**発生機構。** 複合 loss の重み、局所極小、勾配消失、soft boundary/PDE constraints が、物理解でない低 loss 点を作る。[Krishnapriyan2021]

**具体的検出テスト。**

1. 複数 seed・初期化・optimizer で multi-start し、失敗を含む全 run を保存する。
2. loss の総和だけでなく、各 PDE 項、境界、初期値、正規化、発散、楕円関係を別々に報告する。
3. loss weight を各方向に 10 倍変え、同じ物理候補へ収束するか見る。
4. 自動微分 residual を、明示係数微分または finite-difference spot check で独立検証する。
5. 得られた係数を Newton–Krylov / trust-region 法へ渡し、真の離散方程式の根へ改善するか確認する。
6. zero/trivial solution、境界だけ満たす解、時間に依らない偽解を adversarial unit test にする。

**不合格ゲート。** loss weight 変更で解が質的に変わる、PDE の一項だけ大きい、勾配は小さいが equation residual が大きい、または独立 root solver が近傍の根を見つけられない場合。

### TM-14 同じコードによる循環的な検証

**発生機構。** solver、保存、diagnostic、test が同じ誤った微分行列・符号・速度回復式を共有すると、完全に一致しても誤りを検出できない。

**具体的検出テスト。**

1. 重要恒等式ごとに独立経路を持つ。
   - 発散: 円柱公式と Cartesian 再構成。
   - 渦度: 解析式と Cartesian curl。
   - 楕円関係: 係数微分と実空間差分。
   - 時間 residual: stepper RHS と checkpoint 差分。
2. manufactured solution の期待値を solver の同じ関数から生成しない。
3. 少なくとも一つは別離散化（有限差分対スペクトル）、別 library、可能なら別言語で再現する。
4. 方程式の符号反転、軸 parity 破壊、診断 JSON 改変を行う negative test が必ず失敗することを確認する。
5. 保存ファイルに schema version、shape、dtype、単位、hash を持たせる。

**不合格ゲート。** 検証コードが production derivative/RHS をそのまま呼ぶだけ、negative test が通る、または別離散化で主要量が一致しない場合。

### TM-15 選択バイアス・都合のよい実験だけの採用

**発生機構。** seed、初期値、fit window、解像度、境界条件のうち「発散らしい」ものだけを残すと、偶然・数値不安定・多重比較を発見と誤認する。

**具体的検出テスト。**

1. run 開始時に一意 ID と config hash を発行し、成功・失敗・中断を同じ manifest に追記する。
2. 出力を上書きせず、除外理由を machine-readable に残す。
3. parameter sweep の全点を表または欠測理由付きで公開する。
4. primary endpoint、停止規則、fit window 集合を事前登録する。
5. 探索用 data と確認用 data/seed を分離し、確認系列を一度だけ評価する。
6. 多数の候補から選んだ場合は、選択数と ranking rule を報告し、best run だけでなく分布を示す。
7. git commit、dirty state、Python/依存版、CPU/GPU、乱数生成器状態を保存する。

**不合格ゲート。** 欠落 run の理由が説明できない、出力が上書きされている、結果を見た後に primary metric が変更された、または確認用系列で再現しない場合。

## 4. 追加の高リスク脅威

### TM-16 adaptive mapping / 再メッシュ誤差

- mapping の単調性と Jacobian 正値を各 step 検査する。
- 再メッシュ前後の質量的積分、エネルギー、\(\Gamma\) 最大値、PDE residual を比較する。
- 異なる monitor function と固定格子の短時間比較を行う。
- mapping 座標で滑らかでも、物理座標微分を必ず再計算する。

### TM-17 圧力・楕円 solve の不十分な収束

- Poisson/stream-function residual と境界 residual を別保存する。
- tolerance を 10 倍、100 倍厳しくし、速度回復と PDE residual の感度を調べる。
- null space/gauge を明示し、compatibility condition を検査する。
- direct solve と iterative solve、または異なる preconditioner を比較する。

### TM-18 データ破損・診断の後編集

- candidate と diagnostics に暗号学的 hash を付け、manifest に記録する。
- 診断は candidate checkpoint から再生成可能にする。
- 意図的に一値を変更した JSON/CSV が hash または再計算比較で拒否される negative test を持つ。
- plot は診断原本から生成し、plot 内の値を手入力しない。

### TM-19 単位・正規化・座標変換の混同

- 各配列に物理変数名、単位/無次元化、座標、scale factor を metadata として付ける。
- rescaled field から物理 \(u,\omega\) に戻した後、3D divergence、エネルギー、継続判定量を再計算する。
- 形式的な五次元作用素の体積要素と、物理三次元の \(2\pi r\,dr\,dz\) を混同しないテストを置く。

### TM-20 時間積分スキーム由来の偽増幅(2026-07-28 追加)

- Heun+中心差分は純移流で厳密に増幅する(\(|G(i\alpha)|^2=1+\alpha^4/4\))。
  粘性項の相対的大きさから安定性を推論しない。
- 各 run の運転点(全 step の \(\max|u^r|,\max|u^z|,dt,\nu\))で凍結係数
  von Neumann scan(`von_neumann.py`)を実行し、`max|G|>1+tol` なら
  「stability-unverified」と分類する(不安定の証明ではない)。
- 増幅・ピーク位置・core width・エネルギー収支を、虚軸安定な比較積分器
  (SSPRK3/RK4、同一空間離散化・同一 dt 系列)と突き合わせる
  (`run_integrator_comparison.py`)。**Heun 単独の増幅を候補判定に
  使わない。**

### TM-21 診断間引きによる違反見逃し(2026-07-28 追加)

- acceptance-critical 量(エネルギー増分、循環 defect、parity、発散、
  CFL 三点、Poisson 代数残差、エネルギー収支 defect)は全 accepted step で
  計算し、streaming 極値(`gate_summary`)で判定する。出力間引き
  (`diagnostic_stride`)は gate に作用しない。
- 記録行の間だけ違反する合成 trajectory を gate が捕捉することをテストで
  固定する(`test_gate_catches_violation_between_history_rows`)。
- CFL は選択時(pre)・中間段(predictor)・終了時(post)を別々に保存し、
  中間段超過での step 棄却(`stage_cfl_limit`)を利用可能にする。

### TM-22 未解像ピークでの収束外挿(2026-07-28 追加)

- 増幅が plateau しないことを「grid-scale saturation なし」と読まない。
- ピークの points-per-scale(radial/axial FWHM、10–90% front、軸距離
  セル数、勾配スケール、高周波 tail)を全 snapshot で保存する
  (`core_width.py`)。
- 収束 fit は manufactured front で導出した前登録閾値
  (`PREREGISTERED_MIN_POINTS_PER_FRONT = 7`)を満たすまで**禁止**
  (`fit_precondition`)。外部公表値(20.5235 等)を fit の anchor に
  しない(`extrapolation.py` は構造的に anchor を受け付けない)。

## 5. 候補昇格の停止規則

次のいずれかがあれば、候補を昇格せず「偽陽性または未解像」と記録する。

- 空間・時間・領域のいずれかで収束系列がない。
- 軸条件、発散ゼロ、楕円関係、局所エネルギー収支の一つが許容差を超える。
- tail が未解像、または blocking/aliasing の指標が作動する。
- blow-up fit が window、診断量、解像度に安定でない。
- checkpoint から診断を再生成できない。
- independent calculation がない。
- 既知の正則性定理の仮定を満たすのに「特異」と分類している。

逆に、すべてを通っても得られるのは「高品質な数値候補」である。連続 PDE の解の存在、誤差の厳密上界、有限物理時刻、ノルム発散、初期値から候補軌道への接続は、[proof_obligations.md](proof_obligations.md) の別義務である。

## 6. 各 run が保存すべき最小監査記録

- config とその hash、seed、git commit/dirty state。
- 格子、mapping、\(\Delta t\)、accepted/rejected step。
- 候補場 checkpoint と schema/hash。
- PDE 各項、発散、楕円、軸 parity、境界 residual。
- エネルギー、散逸、境界 flux、\(\Gamma\)、BKM/Serrin/渦度診断。
- mode/shell spectrum、tail 指標、集中幅と境界距離。
- fit の全 window、全代替モデル、holdout error。
- run status: success / failed / interrupted / excluded と理由。
- 再生成コマンドとソフトウェア環境。
