# 数値安定性監査(P0-A/P0-B/P0-C/P1-C)

状態: 2026-07-28(branch `fable5-mainline`)。本書は FABLE5_NEXT_TASK_AUDIT.md
§2 の P0 ゲート群に対する実装と実測の記録である。ここに書かれた数値は
すべて数値観察であり、安定性の証明でも不安定性の証明でもない。

## 1. P0-A: Heun+中心差分の von Neumann 監査

### 1.1 数学的事実

定数係数純移流 \(q_t + c q_z = 0\) の中心差分は純虚数シンボル
\(\lambda = -i\,c\sin\theta/\Delta z\) を持ち、Heun の安定性多項式
\(G(z)=1+z+z^2/2\) は

\[
|G(i\alpha)|^2 = 1 + \alpha^4/4 > 1 \quad (\alpha\neq0)
\]

を満たす。**無粘性の中心差分移流に対して Heun は厳密に増幅する。**
粘性がすべての離散波数を安定化するかは運転点ごとの検査事項であり、
「粘性項が項和に占める割合が小さい/大きい」からは何も従わない。

実装: `src/ns_certificate_lab/von_neumann.py`(凍結係数
advection–diffusion シンボル、Heun/Euler 予測子/SSPRK3/RK4 の安定性
多項式、全離散波数 scan、snapshot 監査 API、検証用の独立 1D/2D 参照
propagator)。テスト 11 件(シンボル vs 実配列 propagate の一致
3.7e-16、RK4 1000 step の半離散厳密解一致 2.0e-11、\(\nu=0\) での Heun
増幅検出、他)。

### 1.2 出荷済み運転点の実測(65×128 相当の悲観点と自己整合点)

`outputs/hou_early_time_v1` の記録(`min dt = 2.7587e-7`、
`max advective CFL = 0.10023`)と整合する運転点での凍結係数 worst case:

| 読み方 | Heun full step max\|G\| | 判定(tol 1e-12) |
|---|---|---|
| radial 方向が CFL を運ぶ | 1.000003460085721 | **FAIL** |
| axial 方向が CFL を運ぶ | 1.000000000000000 | PASS |
| 両方向とも記録 CFL | 1.000152277605197 | **FAIL** |
| 悲観的仮定(max\|u^r\|=2000, max\|u^z\|=7000 同時) | 1.0625773 | FAIL(予測子段 1.3110) |

- 判定は 5〜6 桁目で割れる。**分類は「stability-unverified」**であり、
  「不安定」でも「安定」でもない(凍結係数モデルは可変係数の交換子、
  軸行 E-17、壁行 E-27/E-31、段ごとの楕円 solve、0 階生成項を含まない)。
- 同運転点で Heun が全波数 pass する dt は二分法で \(\le2.40\times10^{-8}\)
  (出荷 dt の約 1/11.5)。
- \(\nu=5\times10^{-4}\) の寄与は増幅率を 6.07e-5 だけ下げるに留まる
  (悲観点、advective 増幅 0.129 に対して)。
- scan は 721×721 標本の下界(\(\theta=\pi,0\) は厳密に含む)。181 標本との
  差は相対 5.6e-6。
- \(3\nu/r\) 凍結項は移流虚部と**加算**する向き(保守的読み)を採用。

### 1.3 帰結(決定規則)

1. **Heun 単独で得た増幅を候補判定に使わない**(P0-A item 7)。以後の
   候補判定は積分器相互比較(§1.4)またはその後継を引用しなければ
   ならない。
2. 過去の Heun 実行(`hou_early_time_v1`/`v2_hires`、
   `hou_time_refinement_v1`、`wall_dependence_v1`)は本監査により
   「stability-unverified」と再分類する。無効ではないが、増幅値の
   引用には §1.4 の相互比較の裏づけを添える。
3. 凍結係数 worst case の累積上界(1 step あたり 3.5e-6、2200 step で
   \(\lesssim0.8\%\))は、相互比較の実測差と突き合わせる。

### 1.4 交差検証積分器(SSPRK3 / RK4)

`nonlinear_cylinder.take_step` に SSPRK3(虚軸安定区間 \(|\alpha|\le\sqrt3\))
と古典 RK4(\(|\alpha|\le2\sqrt2\))を追加した。空間離散化・拘束順序・
楕円 solve は Heun と完全に共有する。射影(毎段 `constrain_state`)込みの
実測時間次数: **Heun 1.97/2.00、SSPRK3 3.00/3.00、RK4 3.95/3.98**
(振幅 3000、17×32、\(t=3.2\times10^{-4}\)、同法細 dt 参照)。射影は
観測次数を落とさなかった。

Gate 1 実験 `experiments/run_integrator_comparison.py`
(config `configs/integrator_comparison.json`、前登録許容: 増幅相対差
1e-3、対差は dt 半減で縮小、argmax 1 セル以内): 65×128、E-29 datum、
dt ∈ {6e-7, 3e-7} 固定、3 積分器。結果は `outputs/integrator_comparison_v1`
と STATUS を参照。

## 2. P0-B: CFL 三点測定と段棄却

従来は「step 開始時の状態で dt を選び、step 終了状態で CFL を記録」して
いた(既知・文書化済みの意味論)。本改修で全 accepted step について

- `cfl_pre_state`(選択された dt × 開始状態の方向比)
- `cfl_predictor_stage`(中間段状態での最大。Heun は Euler 予測子、
  SSPRK3 は 2 段、RK4 は 3 段の最大)
- `cfl_post_state`(終了状態)
- `viscous_stability_number`(\(4\nu\,dt/\min(\Delta r^2,\Delta z^2)\))
- dt を制限した拘束の名前(`advective_r`/`advective_z`/`viscous`/
  `max_time_step`/`target_clip`/`fixed`)と棄却回数

を `step_stream` に保存する。`stage_cfl_limit` を与えると、中間段 CFL が
閾値を超えた step は**棄却して dt を半減し再実行**する(最大 25 回、
adaptive モード限定)。受入 step の段 CFL が閾値以下であることはテストで
固定した。

## 3. P0-C: 全 step streaming gate

acceptance-critical 量(エネルギー増分、循環最大原理 defect、奇対称、
軸 parity 相対、相対発散、壁拘束、Poisson 代数残差、CFL 三点、有限値、
エネルギー収支 defect)は **全 accepted step で計算**し、`step_stream`
(全量)と `gate_summary`(streaming 極値)に保存する。出力の間引き
(`diagnostic_stride`)は history 行にのみ作用し、gate は間引けない。

義務づけられた合成 trajectory テスト
(`test_gate_catches_violation_between_history_rows`): 記録行の間だけ
強制パルスでエネルギーを注入→抽出すると、間引き history は単調減衰しか
見せないが、streaming は丸め床の \(3\times10^{12}\) 倍の増加を捕捉する。

Poisson 代数残差は solver A が**実際に解いた線形系**への再適用残差
(丸め量)であり、PDE 離散化誤差ではない。全 step で相対 1e-12 以下を
gate する。

## 4. P1-C: エネルギー収支と粘性符号

### 4.1 恒等式(壁項込み)

E-27 壁は swirl のみ no-slip(\(u_1=0\))で \(u^z(R,z)=R\,\psi_{1,r}\neq0\)
の滑り壁である。したがって正しい連続恒等式は

\[
\frac{dE}{dt} = -\nu\int|\omega|^2\,dV - \nu\oint_{r=R} u^z\,\omega^\theta\,dS
\]

であり、壁項を落とすと defect が物理境界項で汚染される。離散 defect は
台形時間平均で全 step 記録し、壁項なしの読みも併記する
(`energy_balance_defect` / `energy_balance_defect_no_wall`)。

**整合性の記録(隠さない)**: 初回実装は \(\int|\omega|^2dV=2\pi\int
|\omega|^2 r\,dr\,dz\) を E-20 の \(\pi\) 正規化 measure と取り違え、
相対 defect が厳密に 0.5 へ飽和した(符号反転 fault では 1.5)。飽和値が
理論予測と一致したことで因子 2 の欠陥として特定・修正した。修正後の
clean 実測: 滑らか control で相対 5.9e-2 → 1.6e-2 → 4.2e-3(空間時間
同時細分、収束)。

### 4.2 viscosity_sign fault

\(\nu=5\times10^{-4}\) の Hou 実データでは粘性項は項和の \(3\times10^{-4}\)
程度で、符号反転は実データからは識別できない(既知)。新 fault
`viscosity_sign` は拡散支配 control(\(\nu=2\times10^{-2}\)、滑らか datum)で

- 相対 defect: clean **2.21e-2** vs 反転 **2.000**(理論値 2、比 90 倍)
- エネルギー: clean 単調減衰 vs 反転 単調成長

により確実に棄却される。「energy increase 0.0」を「energy identity が
正しい」と読み替えないこと — 恒等式の検査は defect の収束で行う。

### 4.3 項別分解

swirl エネルギー \(E_s=\pi\int(ru_1)^2 r\,dr\,dz\) への移流・stretching・
粘性の仕事率を全 step 記録する(`swirl_power_*`)。meridional 側は
\(\omega_1\mapsto\psi_1\) の楕円診断を通じて駆動されるため分解せず、
全体収支のみ gate する(限界として明記)。

## 5. 本監査が主張しないこと

- von Neumann pass は凍結係数モデルの pass であり、非線形・可変係数系の
  安定性証明ではない。
- 積分器間一致は時間離散化リスクの上界であり、空間離散化誤差(現在
  支配的)には何も言わない。
- エネルギー収支 defect の収束は使用した恒等式と離散化の整合性を示す
  だけで、解の正しさを示さない。
