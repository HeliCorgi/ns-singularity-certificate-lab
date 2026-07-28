# 非線形有限円柱ソルバ設計(production)

状態: 設計(実装前)。実装前提はすべて `docs/equation_audit.md` の
監査済み式(E-11–E-17、E-20、E-23、E-25–E-31)であり、
未確認・不整合・誤りの式は使用しない。

目的: Hou 幾何(壁付き有限円柱)での完全非線形軸対称+swirl 系
\((u_1,\omega_1,\psi_1)\) の時間発展を、固定一様格子・2 次精度で計算し、
複数解像度の早期 Hou 再現(E-30 プロトコル)を可能にする。

これは候補生成側の実装である。検証側(独立 Poisson solver B、
`pde.py` 残差、Cartesian adapter、4 次 stencil 監査)とは楕円 solve・
微分 stencil を共有しない経路で snapshot を再検査する。

## 1. 幾何と状態

- 領域: \(r\in[0,1]\)、\(z\in[0,1)\) **フル周期**(周期 1)。
  Hou の半周期領域 \(D_1\)(E-28)は使わず、奇対称性は**課さずに監視**する
  (破れの成長は数値誤差の診断になる。sine 級数による半周期実装は
  将来の独立検証経路として保留)。
- 格子: `AxisymmetricGrid.uniform(nr, nz, r_max=1.0, z_min=0.0, z_max=1.0,
  periodic_z=True)`。配列は \((n_r,n_z)\) の float64。
- 発展変数: \(u_1,\omega_1\)。\(\psi_1\) は各段で E-13 を解いて得る診断量。

## 2. 半離散系

監査済み演算子(`operators.py` の 2 次中心差分、軸は
`even_at_axis=True` と E-17)を \(D_r,D_z,L_5^h\) と書く。

\[
\partial_tu_1=-u^rD_ru_1-u^zD_zu_1+2u_1D_z\psi_1+\nu L_5^hu_1+F_{u_1},
\]
\[
\partial_t\omega_1=-u^rD_r\omega_1-u^zD_z\omega_1+D_z(u_1^2)+\nu L_5^h\omega_1+F_{\omega_1}.
\]

\(F_{u_1},F_{\omega_1}\) は manufactured 検証専用の強制項で、production run
では 0。各 RHS 評価で

1. \(-L_5\psi_1=\omega_1\) を solver A(`poisson.py`、外側 Dirichlet
   \(g\equiv0\)= E-27)で解く。
2. \(u^r=-r\,D_z\psi_1\)、\(u^z=2\psi_1+r\,D_r\psi_1\)(E-14)。

## 3. 境界処理

- **壁 \(r=1\)(index \(n_r-1\))**: 各段の状態更新後に
  \(u_1[n_r{-}1,:]=0\)(E-27)、
  \(\omega_1[n_r{-}1,:]=-(8\psi_1[n_r{-}2,:]-\psi_1[n_r{-}3,:])/(2\Delta r^2)\)
  (E-31、2 次 Thom 型)。壁行は PDE で発展させない
  (`run_nonlinear_control.py` の外側 trace 上書きと同型)。
- **Poisson**: Dirichlet \(\psi_1(1,z)=0\) のみ(E-27 監査警告)。
- **軸 \(r=0\)**: 追加の行上書きなし。演算子の偶対称処理(E-16/E-17)が
  正則性を保つ。\(u^r(0,z)=0\) は \(u^r=-r\psi_{1,z}\) から恒等的に成立。
- **z**: 周期(演算子内蔵)。奇対称 defect
  \(\max|f(r,z)+f(r,1-z)|\)(格子対称点)を診断として保存。

## 4. 時間積分と刻み

- Heun / 陽的 RK2(E-30 と同型; Butcher \(c=(0,1)\), \(a_{21}=1\),
  \(b=(1/2,1/2)\))。
- 適応刻み:
  \[
  \Delta t=\min\!\Big(
  C\,\frac{\Delta r}{\max|u^r|+\varepsilon},\
  C\,\frac{\Delta z}{\max|u^z|+\varepsilon},\
  C\,\frac{\min(\Delta r^2,\Delta z^2)}{4\nu}\Big),\qquad C=0.1
  \]
  (係数 4 は E-17 の軸極限 \(4\partial_{rr}\) の安定性余裕)。
  manufactured 収束テスト用に固定 \(\Delta t\) モードを持つ。
- 二段階粘性: 設定 `viscosity_schedule = [[0.0, 5e-4], [0.00227375, 5e-3]]`
  形式(E-30)。切替は step 境界で行い、切替時刻を診断へ記録する。

## 5. 診断(step/snapshot ごと)

エネルギー \(E(t)\)(E-20 の \(r\,dr\,dz\) 測度)、
\(\|u_1\|_\infty\)、\(\|\omega_1\|_\infty\)、Cartesian 渦度最大
\(\|\omega\|_\infty\)(E-18b)、\(u_1\) の argmax \((R(t),Z(t))\)(E-30)、
増幅率 \(\|\omega\|_\infty/\|\omega(0)\|_\infty\)、
enstrophy、物理発散残差(E-02)、独立 solver B による楕円 defect、
軸 parity defect、z 奇対称 defect、循環最大 \(\|r^2u_1\|_\infty\)
(E-23: 非増加が厳密な必要条件 — 増加は即 stop condition)、
CFL 実効値、\(\Delta t\) 履歴、壁残差。

無強制 run では \(E(t)\) 非増加(粘性散逸+壁 no-slip で流入なし)を
acceptance に含める。

## 6. checkpoint と証拠

- candidate schema v2(`artifacts.py`)で \(u_1,\omega_1,\psi_1\) を
  指定時刻と最終時刻に保存(単位・無次元化・\(\nu\)・物理時刻・provenance
  込み)。
- 実験は非空出力ディレクトリを拒否。config snapshot、summary.json
  (acceptance、limitations、再現コマンド)、CSV 時系列、manifest+SHA-256。
- 再開(restart)テスト: checkpoint から再開した run が連続 run と
  許容差内で一致すること。

## 7. テスト(実装受入)

1. **E-31 単体収束**: 壁条件を満たす解析場
   \(\psi=(1-r^2)^2\phi(r,z)\) で Thom 式が \(-\psi_{rr}(1,z)\) へ 2 次収束。
2. **強制 manufactured(全結合非線形+壁)**: \(\psi_1^*=(1-r^2)^2\,
   a(t)\cos(2\pi z)\cdot\rho(r)\)、\(u_1^*=(1-r^2)^2\,b(t)\sin(2\pi z)
   \cdot\sigma(r)\) 型(壁条件を厳密に満たし \(r\) 偶)。強制項は閉形式解析
   微分から独立導出。3 解像度で全変数の観測次数 \(\ge1.85\)。
3. **固定格子時間収束**: \(\Delta t,\Delta t/2,\Delta t/4\) で次数 ≈2。
4. **零場不動点**: \(u_1=\omega_1=0\) が厳密に保存。
5. **小振幅極限**: \(u_1\) のみ小振幅 → 線形 swirl 拡散
   (既存 Crank–Nicolson 対照)と一致(振幅 2 乗のオーダーで乖離)。
6. **パリティ・対称性保存**: 軸偶対称・z 奇対称の defect が丸め+離散化
   オーダーに留まる。
7. **循環最大原理**: 無強制 run で \(\|r^2u_1\|_\infty\) 非増加。
8. **故障注入**(`test_detects_*`、清浄 run の 10 倍超で検出):
   Thom 符号反転、stretching 項符号反転、軸係数 8→4、壁 \(u_1\) 拘束解除、
   Poisson 符号反転、粘性切替無視。
9. **restart 忠実性**。
10. **独立残差**: 保存 snapshot に対し solver B の楕円 defect、
    `pde.py` 残差、E-02 発散、Cartesian adapter 検査。

いずれかの収束が失敗した場合、Hou 初期値の run は実行しない
(handoff §6.3 のブロック規定)。

## 8. 早期 Hou 再現実験(`run_hou_early_time`)

- 初期値 E-29、\(\nu\) スケジュール E-30(第 1 段 \(5\times10^{-4}\))。
- 解像度 3 段(例: \(129\times256\)、\(193\times384\)、\(257\times512\))、
  フル周期 z。
- 目標時刻: まず \(t\in[0,T_1=0.002191729]\)。資源が許せば \(t_0\) まで。
- 記録: 増幅率軌跡、\((R(t),Z(t))\)、早期の \(\|u_1\|_\infty\) 減少
  (E-30 の定性ターゲット)、解像度間の一致区間、1536² 論文値
  (増幅 20.5235 at \(T_1\))との差。〔改版 2026-07-28(P0-D): 旧項目
  「grid-scale 飽和の有無」を撤回し、`core_width` の points-per-scale
  指標(radial/axial FWHM、10–90% front thickness、ピーク軸距離セル数、
  高周波 tail)へ置き換える。増幅の plateau 有無は saturation の
  判定材料にしない。〕
- 一様固定格子は Hou の適応格子(最小 \(10^{-8}\))に遠く及ばないため、
  ここでの主張上限は「解像度依存性を持つ数値観察」であり、
  一致は期待値ではなく測定対象。乖離も証拠として保存する。

## 9. stop conditions(このマイルストーン固有)

- 循環最大 \(\|r^2u_1\|_\infty\) の増加(許容差超)。
- 無強制 run のエネルギー増加。
- 奇対称 defect が場の振幅に対して \(O(1)\) へ成長。
- 解像度で growth が系統的に変わり続け、傾向が同定できない。
- manufactured 収束の失敗。

発生時は run を停止し、証拠を保存して STATUS に記録する。
