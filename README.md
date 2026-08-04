# Spectral-front identities for periodic Navier–Stokes<br>and an obstruction to uniform pointwise Osgood closure

> ## Not a Clay solution
>
> **This repository does not solve the Clay Navier–Stokes problem, in either
> direction.** No global regularity for all smooth data, no blow-up, no partial
> solution.
>
> **Proven** (modulo three classical inputs — Kato-type local theory,
> subcritical Serrin regularity, mean-zero Sobolev embedding): a
> bandwidth–dissipation dichotomy on the torus, obtained from exact modal
> identities, together with an exact representation of the action it turns on;
> and, separately, an obstruction showing that the dichotomy cannot be closed
> by any uniform pointwise Osgood majorant. **Rejected:** several blow-up
> mechanisms, each kept on record with the equation or theorem that killed it.
>
> Review basis: exact-rational machine certificates and adversarial referee
> passes run by **AI agents — this is not human peer review.**

**Five-minute path:** the theorem below → the status table → the paper
([PDF](docs/paper_lambda_dichotomy/paper_draft.pdf) ·
[statement](docs/paper_lambda_dichotomy/theorem_statement.md) ·
[proof](docs/paper_lambda_dichotomy/complete_proof.md) ·
[gap audit](docs/paper_lambda_dichotomy/dependency_and_gap_audit.md)).
Questions and corrections are welcome — please open an issue.
日本語の詳細記録は本ページ後半にあります。

---

## Main theorem in one screen

Unforced Navier–Stokes on $\mathbb{T}^3$, viscosity $\nu>0$, maximal strong
solution $u$ on $[0,T_{\max})$ from $u_0\in H^m_\sigma$, $m>5/2$, $u_0\neq0$.
Write $H_r=\sum_k |k|^{2r}|\hat u_k|^2$, $D=H_1$, bandwidth $N_0^2=H_1/H_0$,
$N_1^2=H_2/H_1$, $z=\log N_0^2\ge 0$, and

$$\mathcal{N}=-\mathbb{P}(u\cdot\nabla u),\qquad K=\frac{\lVert\mathcal{N}\rVert_2^2}{\lVert\nabla u\rVert_2^4}.$$

**Theorem.**

- **(a)** $\Lambda(t)=\log N_0^2(t)-\dfrac{1}{2\nu}\displaystyle\int_0^t K D\,ds$ is non-increasing on $[0,T_{\max})$, with an explicit nonnegative two-part defect.
- **(b)** Either $\int_0^{T_{\max}}K D\,dt<\infty$ — and then $z$ is bounded, $u\in L^\infty(0,T_{\max};H^1)$ and $T_{\max}=\infty$ — or the integral diverges. In particular $T_{\max}<\infty$ forces divergence. *(No claim that divergence forces blow-up.)*
- **(c)** Exact identity: $\displaystyle\int_0^{T'} K D\,dt=\int_0^{T'}\frac{\lVert\partial_t u\rVert_2^2}{D}\,dt+\nu^2\int_0^{T'}N_1^2\,dt+\nu\log\frac{D(T')}{D(0)}$, whence by AM–GM finiteness of the $\dot H^1$-bandwidth action implies globality. One-sided only: an explicit decaying solution kills the converse.
- **(d)** $K D\le\lVert u\rVert_{L^\infty}^2$ and $K D\le C_S^2\lVert\nabla u\rVert_{L^3}^2$, so the criterion in (b) is implied by the Serrin $(\infty,2)$ and the critical $L^3$ gradient actions. It is **not** claimed to be new or best: the paper's bibliography records that Cheskidov–Shvydkoy already give a wavenumber criterion provably weaker than *every* Ladyzhenskaya–Prodi–Serrin condition, that $N_0$ and $N_1$ are the classical Taylor-microscale inverse lengths, and that $K$ measures the known "nonlinear depletion".
- **(e)** Any majorant $K D\le\Phi(z)D+R$ with $\Phi$ nondecreasing, $\int^\infty ds/\Phi=\infty$ and $\int_0^{T_{\max}}R\,dt<\infty$ puts the solution in the global case of (b).

**Theorem (obstruction).** Let $\Phi$ be nondecreasing with
$K(u)\le\Phi(\log N_0^2(u))$ for **every** real zero-mean divergence-free
trigonometric field $u$. Then $\Phi(s)\ge c\,e^{s}$ for all large $s$, so
$\int^\infty ds/\Phi<\infty$: **no Osgood-admissible $\Phi$ exists**, and the
remainder-free route into (e) is closed. The proof exhibits one explicit
family — the coherent critical-spectrum field
$\hat u_N(k)=\chi(|k|/N)\,P_k v_0/|k|^2$ — for which
$\lVert\mathbb{P}(u_N\cdot\nabla u_N)\rVert_2^2\ge c_0N^3$ holds for **every**
admissible cutoff profile $\chi$ and **every** nonzero seed vector $v_0$, with
$c_0=c_0(\chi,v_0)$ and the threshold $N\ge N_*(\chi,v_0)$; the constants are
not uniform in $(\chi,v_0)$, and they are non-effective. Solution-adapted
remainders $R(t)$ are **not** excluded.

## Proven / conditional / rejected

| Status | Content |
|---|---|
| **Proven** | (a)–(e) above; the exact action identity; the obstruction theorem, including the capacity lower bound for the smoothly truncated family that drives it; the exact spectral laws of that family. External inputs: Kato-type local theory, subcritical Serrin, mean-zero Sobolev — nothing else. Constants in (a)–(e) are explicit; the obstruction's constant is **not** effective. |
| **Conditional** | Every PDE tube certificate in section C1/C2 below, on its named external theorems. |
| **Rejected, kept on record** | Steady self-similar front (Tsai 1998, and Chae–Wolf at the weak-$L^3$ level); the exact Leray cycle gate; continuum-to-lattice shadowing; the action–bandwidth *equivalence* (withdrawn after refereeing); the uniform deficit floor; a proven no-go for the constant-vector-split route to the capacity bound. See [VERDICTS.md](docs/research_notes/verification_sprint_v1/VERDICTS.md). |
| **Still open, no longer used** | The capacity bound for the *sharply* truncated family, the repository's original Hypothesis L\*. The obstruction no longer depends on it, because its hypothesis quantifies over all fields and the smooth family suffices. |

## Verify it yourself — one command

```bash
python scripts/verify.py
```

It installs the two dependencies, runs the whole suite, and then prints the
scope of what a green run does and does not establish — in particular that it
does **not** verify Hypothesis L\*, the paper's infinite-dimensional analysis,
or the Lean development. Lean is deliberately outside this command; the
optional second command is in the reproduction section below.

## Where things are

| You want | Go to |
|---|---|
| the paper as one file | [paper_draft.pdf](docs/paper_lambda_dichotomy/paper_draft.pdf) (built from source by CI) |
| the theorem, its proof, its gaps | [docs/paper_lambda_dichotomy/](docs/paper_lambda_dichotomy/theorem_statement.md) |
| what was tried and killed | [VERDICTS.md](docs/research_notes/verification_sprint_v1/VERDICTS.md) |
| the AI referee reports | [A](docs/paper_lambda_dichotomy/referee_report_A.md), [B](docs/paper_lambda_dichotomy/referee_report_B.md) |
| every verified result, with assumptions | [docs/verified_results.md](docs/verified_results.md) |
| current state and open obligations | [STATUS.md](STATUS.md) |

## The position and value of this result

*A calibration note: what is new here, what is not, and what this work is
worth. 日本語版はこの節の後半にあります。*

**In one sentence.** For the Navier–Stokes equations on a periodic box, this
repository proves exact identities that meter how fast a solution's energy can
migrate to finer scales, and proves that one natural way of turning those
identities into a regularity proof — a uniform pointwise Osgood bound — cannot
work.

**For the non-specialist.** The Navier–Stokes equations describe how fluids
move. Whether their three-dimensional solutions can spontaneously develop
infinitely fine structure in finite time is a famous open problem, and this
repository does not settle it, in either direction. What it contains instead
is exact bookkeeping: an identity tracking a single number that measures how
fine the flow's structure currently is; a proof that a certain integral
controls all possible growth of that number; and a proof that a tempting
shortcut — bounding that integral by a slowly growing function of the number
itself, which a classical lemma would then convert into a full regularity
proof — is impossible, because explicit velocity fields violate every such
bound. Knowing precisely why a natural route fails, with the counterexamples
written down, is a modest but real piece of mathematical information.

**What is clarified, precisely.** Three things. (i) An exact modal identity
for the spectral bandwidth $N_0^2=H_1/H_0$, whose Cauchy–Schwarz and
square-completion defects are exactly computable — the identity and its
defect decomposition are verified in exact rational arithmetic on finite
fields. (ii) A dichotomy — finiteness of $\int K D\,dt$ forces global
existence — together with the exact representation
$\int KD=\int\lVert\partial_tu\rVert_2^2/D+\nu^2\int H_2/H_1+\nu\log(D/D_0)$,
which isolates the supercritical content of the criterion in the
$\dot H^1$-bandwidth action and a harmless logarithm. (iii) An obstruction:
no nondecreasing $\Phi$ with $\int^\infty ds/\Phi=\infty$ satisfies
$K\le\Phi(\log N_0^2)$ over all divergence-free trigonometric fields,
witnessed by an explicit spectral family with capacity of order $N^3$.

**What is not new.** The building blocks are classical: $N_0$ and $N_1$ are
the Taylor-microscale inverse lengths (Foias–Guillopé–Temam and later work),
$K$ quantifies the long-studied depletion of the nonlinearity, and
logarithmically refined regularity criteria go back to Kozono–Ogawa–Taniuchi
and Planchon. Cheskidov–Shvydkoy proved a wavenumber criterion strictly
weaker than every Ladyzhenskaya–Prodi–Serrin condition; this work does not
claim to match or beat it — the criterion here is only shown to be *implied
by* two classical actions. The contribution is the exact packaging:
identities with exactly computable defects, the action representation, and a
certified counterexample family. A closed door mapped precisely, not a new
door opened.

**What it can be used for.** Structure: the action representation shows
exactly where the supercriticality of the problem sits. Negative knowledge:
the pointwise-uniform Osgood route is now provably dead, which spares future
effort and sharpens what any successful attempt must do differently
(solution-adapted remainders $R(t)$ are untouched by the obstruction). Raw
material: an exact-rational certificate layer with independent checkers, a
coherent spectral family with closed-form laws, and a ledger of rejected
blow-up mechanisms, each recorded with the equation or theorem that killed
it.

**Value this work does not have.** The Clay problem is unresolved in both
directions and nothing here bears on it beyond the statements above. There
is no practical or engineering application. Review was adversarial and by AI
agents, with reports committed unedited — it is not human peer review, and
during it several claims (a quantifier, a periodisation identity, a
positivity argument) were broken and corrected: evidence that the process
catches errors, and equally evidence that the material contained them. The
mechanized part covers only finite-sum algebra; the monotone quantity itself
requires a time integration that is not formalized. The obstruction's
constants are non-effective, and the capacity bound for the sharply
truncated family — the original hypothesis — remains open, though nothing
now depends on it.

**Method, briefly.** Every load-bearing finite computation is an
exact-rational certificate with an independent checker and a
tamper-rejection battery; the paper-level proofs went through adversarial
review by AI agents whose reports are in the repository; corrections are
published as tagged releases rather than silent edits.

### この成果の位置づけと価値（日本語）

**一言で。** 周期箱上の Navier–Stokes 方程式について、解のエネルギーが
より細かいスケールへ移る速さを厳密に測る恒等式群を証明し、それを正則性
証明へ変換する自然な方法の一つ(帯域の関数による一様な各点 Osgood 上界)
が不可能であることを証明した。

**専門外の方へ。** Navier–Stokes 方程式は流体の運動を記述します。その
3次元の解が有限時間で無限に細かい構造を作り得るかは有名な未解決問題で、
本リポジトリは**どちらの方向にもこれを解決していません**。ここにあるのは
厳密な帳簿づけです: 流れの構造の細かさを測る一つの数を追跡する恒等式、
その数の増大をある積分が完全に支配することの証明、そして魅力的な近道 —
その積分を「細かさ自身のゆっくり増える関数」で抑え、古典的な補題で正則性
証明に変換する方法 — が不可能であることの証明(あらゆるそうした上界を
破る速度場を明示的に構成)。自然な経路がなぜ失敗するかを反例つきで正確に
知ることは、小さいが実在する数学的情報です。

**数学的に明確になったこと。** (i) スペクトル帯域 $N_0^2=H_1/H_0$ の厳密
なモード恒等式。その Cauchy–Schwarz 欠損と平方完成欠損は厳密に計算可能
で、有限場上で有理演算により検証済み。 (ii) 二分法($\int KD\,dt<\infty$
なら解は大域的)と、厳密な作用表現
$\int KD=\int\lVert\partial_tu\rVert_2^2/D+\nu^2\int H_2/H_1+\nu\log(D/D_0)$。
後者は判定条件の超臨界性が $\dot H^1$ 帯域作用に局在し、残りは無害な対数
だけであることを示す。 (iii) 障害: $\int^\infty ds/\Phi=\infty$ を満たす
非減少 $\Phi$ で $K\le\Phi(\log N_0^2)$ を全ての発散ゼロ三角多項式場で
満たすものは存在しない。証人は容量が $N^3$ 次の明示的スペクトル族。

**新しくないもの。** 部品は古典的です: $N_0,N_1$ は Taylor 微視スケール
の逆長さ(Foias–Guillopé–Temam ほか)、$K$ は長く研究されてきた非線形性
の枯渇の定量化、対数的に精密化された正則性判定は Kozono–Ogawa–Taniuchi
や Planchon に遡ります。Cheskidov–Shvydkoy は「あらゆる
Ladyzhenskaya–Prodi–Serrin 条件より真に弱い」波数判定を既に証明しており、
本研究はそれに並ぶとも勝るとも主張しません — ここでの判定は古典的な
2つの作用に*含意される*ことしか示していません。貢献は厳密なパッケージング
です: 欠損が厳密計算できる恒等式、作用表現、認証つき反例族。新しい扉を
開いたのではなく、閉じている扉を正確に地図にしたものです。

**何の役に立つか。** 構造の理解: 作用表現は問題の超臨界性がどこに座って
いるかを正確に示します。負の知識: 各点一様な Osgood 経路は証明つきで
死んだので、後続の試みが何を変えなければならないかが明確になります
(解に適応した剰余 $R(t)$ の経路はこの障害の対象外のまま残ります)。
素材: 独立検査器つきの厳密有理証明書層、閉形式の法則を持つスペクトル族、
そして棄却された爆発機構の台帳(各項目に、それを殺した式または定理を
記録)。

**持たない価値。** Clay 問題は両方向とも未解決で、本成果は上記の言明を
超えてそれに関与しません。実用的・工学的応用はありません。検証は AI
エージェントによる敵対的レビューで、報告書は未編集のまま収録されています
が、人間の査読ではありません — その過程で複数の主張(量化子、周期化
恒等式、正値性の議論)が破られ修正されました。これは過程が誤りを捕まえる
証拠であると同時に、素材に誤りが含まれていた証拠でもあります。機械検証は
有限和の代数のみで、単調量そのものには形式化されていない時間積分が必要
です。障害の定数は非有効で、鋭いカットオフ族の容量下界(元の仮説)は
未解決のまま残っています(ただし現在は何もそれに依存しません)。

**方法論。** 負荷のかかる有限計算はすべて、独立検査器と改竄拒否テスト
つきの厳密有理証明書です。紙上証明は AI エージェントの敵対的レビューを
経ており、その報告書はリポジトリに収録されています。訂正は黙った編集では
なく、タグ付き Release として公開されます。

---

## このリポジトリについて（日本語）

3次元非圧縮 Navier–Stokes 方程式について、**数値候補探索・候補除外・有限次元
形式証明・厳密有理証明書・computer-assisted PDE verification** を同一の監査
可能な基盤で扱うリポジトリです。

現在、Clay 公式命題を満たす特異点候補も、大域正則性の証明もありません。
本リポジトリは **数値観測 / Lean 証明 / 証明書検査 / 監査済み紙上解析** を
常に区別して記録します。各成果の正確な仮定・再現コマンド・信頼基盤・
適用範囲は下記および [docs/verified_results.md](docs/verified_results.md)
にあります。

---

## Current representative verified result (hybrid) / 現在の代表的検証結果

本結果の区分は下記語彙表の **hybrid**(§C2)です。**Lean だけの結果でも、
条件付き証明書でもありません。**

**明示された有限 Fourier 初期値・粘性・短時間区間について、Galerkin 軌道・
連続 PDE 残差・control ODE・スラブ連結を機械検査し、監査済み(= **AI エージェントによる
独立監査パス**。人間の査読ではありません)の古典的 PDE 解析と組み合わせる
ことで、一意な周期強解が指定された Sobolev tube 内に存在することを検証
できます。**

| 項目 | 値 |
|---|---|
| 領域 | 周期 $\mathbb{T}^3 = (\mathbb{R}/2\pi\mathbb{Z})^3$、正規化測度 $(2\pi)^{-3}dx$ |
| 初期値 | 族 P1 — 厳密有理係数の有限 Fourier 場、モード $(1,0,0),(0,1,0),(1,1,0)$、平均ゼロ・厳密発散ゼロ(有限三角多項式なので $C^\infty$) |
| 粘性 | $\nu = 1/10$(固定・正) |
| Galerkin 帯域 | $1 \le \lvert k\rvert^2 \le 4$ |
| certified interval | $[0,\ 5/256]$($\approx 1.953\times10^{-2}$、$h = 1/2048$ のスラブ 40 本) |
| Sobolev 次数 | $\dot H^3$(参考: $H^3 \le \sqrt{8}\,\dot H^3$) |
| tube 半径 | $\lVert u - u_a\rVert_{\dot H^3} \le 0.030903$ — 初期値ノルム $\lVert u_0\rVert_{\dot H^3} = \sqrt2$ の約 **2.19%**($H^3$ 換算 $\le 0.087407$) |
| 停止理由 | `slab_budget_exhausted`(前登録の 40 スラブ予算に到達。証明の破綻ではない) |
| 成果物 | [`outputs/track_p_chain_reissued_v2/reissued_h3chain_strict_same_step_P1_nu_1over10.json`](outputs/track_p_chain_reissued_v2/reissued_h3chain_strict_same_step_P1_nu_1over10.json) |

三層の役割分担:

- **Lean が検査する部分** — スラブ合成と連結の有限論理(`two_slab_composition`、
  `chain_composition`)、Taylor 終端剰余、積分形比較(`integral_comparison`、
  `integral_riccati_comparison`)、Kato 定数組立の有限代数、有限次元 ODE の
  局所存在・一意性。
- **Python checker が検査する部分** — 各リンクの厳密有理再計算: 初期値の
  発散ゼロ・帯域所属、Galerkin 軌道の Picard 包含、control 不等式の定数、
  厳密な連続 NS 残差(Galerkin tail)、control ODE 管、再中心化点の一致、
  $\delta$ 漸化式、文言契約(チェーン checker は各リンクの定数を実際に
  再計算します。単発スラブ checker の守備範囲は §B の表を参照)。
- **監査済み紙上定理に依存する部分** — 周期 $H^4$ 局所存在・一意性
  (EXT-P1★)、積分形エネルギー不等式(EXT-P2-INT + 比較補題)、継続原理
  (EXT-P3★ と系 P3-3)。**これらは Lean 形式化されていません。**

**三層の接続は機械化されていません。** Lean 定理はいずれも抽象ノルム空間上の
言明であり、証明書の数値が Lean へ入力されることはありません。対応付けは
人手です(登録簿 VR-L-001 / VR-L-014)。

同じ再発行セットには、同一の hybrid 信頼基盤で $T \approx 5.86\times10^{-2}$
(P1、$\nu = 1/100$)に達するチェーンもあります。本表が strict same-step 実行を
選ぶのは、刻みを固定した旧定数との厳密比較のためであり、tube が最も締まって
いる(相対 2.19%)ためです。

この結果は特定初期値・特定短時間区間についてのものであり、大域正則性でも
特異点構成でもありません。certified interval は証明手法の到達範囲です。

厳密な再現は、成果物を生成したコミットを指定してください:

```bash
git checkout 17d41df
```

---

<!-- MCR:BEGIN (machine-checkable results; guarded by tests/test_readme_claims.py) -->
## Machine-checkable results / 第三者が機械検証できる成果

**用語規則(本 README と登録簿で厳守):**

| 用語 | 意味 |
|---|---|
| **Lean-verified** | Lean kernel が検査。project 固有の未証明 axiom なし |
| **certificate-verified** | 独立 Python checker が証明書を再計算・検査(範囲は §B の行ごと)。**Lean 証明ではない** |
| **hybrid** | Lean-verified 有限論理 + certificate-verified 有理計算 + **監査済みだが未形式化**の古典解析 |
| **conditional** | 明示された未証明の外部仮定が残る |
| **numerically observed** | 浮動小数点計算のみ。証明ではない |

<!-- NOTCLAIMED:BEGIN (phrases listed here are explicitly disclaimed, not asserted) -->
「fully Lean-verified PDE theorem」「Lean だけで周期 Navier–Stokes 解の存在を
証明」「Clay 問題への部分解」「大域正則性証明」「特異点証明」は、いずれも
**本リポジトリの成果ではなく、使用しません**。
<!-- NOTCLAIMED:END -->

以下は各区分の代表例です。**全定理・全証明書の完全な一覧は
[docs/verified_results.md](docs/verified_results.md)**(安定 ID 付き登録簿)に
あります。

### A. Lean-verified results

Lean kernel のみで検査され、`sorry` / `admit` / project 固有 axiom を含みません。
公理監査の全項目は [formal/AxiomAudit.lean](formal/AxiomAudit.lean) を参照して
ください(件数は開発とともに変動するため本文には書きません)。

| Result | Lean theorem | File | Exact verified claim | Explicit limitation |
|---|---|---|---|---|
| 2 スラブ合成 | `two_slab_composition` | [TrackPChain.lean](formal/NSSingularity/TrackPChain.lean) | 抽象 tube・転送上界・予算不等式から piecewise tube が従う | tube を PDE に対して供給することは含まない |
| n スラブ連結 | `chain_composition`, `chain_composition_union` | [TrackPChain.lean](formal/NSSingularity/TrackPChain.lean) | $\delta_{n+1} \ge R_n + \text{transfer}$ に忠実なリスト帰納法(証明書層の等式はその特例) | 同上 |
| 転送三角不等式 | `transfer_triangle` | [TrackPChain.lean](formal/NSSingularity/TrackPChain.lean) | 再中心化予算の 3 項分解 | 各項の数値は証明書層 |
| Taylor 終端剰余 | `taylor_endpoint_remainder_bound` | [TrackPChain.lean](formal/NSSingularity/TrackPChain.lean) | $\lvert f(t_0+h) - \text{Taylor}_m\rvert \le M h^{m+1}/(m+1)!$ | 上界 $M$ の値は証明書層 |
| 積分形比較 | `integral_comparison`, `integral_riccati_comparison` | [ChainAnalysis.lean](formal/NSSingularity/ChainAnalysis.lean) | 2 パラメータ積分不等式を満たす連続関数は ODE 解に支配される(**Dini 微分不要**) | PDE 差分場がその積分不等式を満たすことは含まない。Lean 版は $W$ レベルの不等式を仮定し、EXT-P2-INT が供給するのは $W^2$ レベル — 両者を繋ぐ $\eta$-retreat(比較補題本体)は紙上のみ |
| 有限次元 Picard–Lindelöf | `quadratic_ode_local_solution`, `quadratic_ode_unique`, `galerkin_local_solution` | [GalerkinPicard.lean](formal/NSSingularity/GalerkinPicard.lean), [TimeDependentGalerkin.lean](formal/NSSingularity/TimeDependentGalerkin.lean) | 二次場 $u' = Au + B(u,u)$ の局所存在・一意性(Lipschitz 定数と存在区間を明示) | EXT-P1(PDE 命題)そのものではない |
| Galerkin ノルム有界 | `galerkin_norm_le`, `norm_le_of_energy_inequality` | [GalerkinNoBlowup.lean](formal/NSSingularity/GalerkinNoBlowup.lean) | エネルギー不等式の下で有限次元ノルムが有界 | PDE のエネルギー不等式の成立は仮定 |
| 固定有限帯域の no-go | `FixedBandwidthCandidate.breakdown_times_empty`, `.reaches_every_time` | [CertificateLayer.lean](formal/NSSingularity/CertificateLayer.lean) | 固定帯域に**留まる軌道**は破綻時刻を持たない | 抽象ノルム空間上の**係数軌道**についての言明であり NS 解についてではない(発展方程式は structure の field = 仮定)。有限帯域**初期値**も除外しない(下記反例) |
| 帯域の区別 | `exists_finiteBandDatum_not_fixedBandTrajectory`, `FixedBandwidthCandidate.fixedBand_scope` | [TrackPFourier.lean](formal/NSSingularity/TrackPFourier.lean) | 有限帯域初期値は固定帯域軌道を含意しない(反例つき) | — |
| 三線形相殺・Leray 代数 | `advectionForm_eq_zero`, `inner_leray_eq_zero`, `leray_leray`, `norm_leray_le` | [FiniteModeNoGo.lean](formal/NSSingularity/FiniteModeNoGo.lean), [TrackPFourier.lean](formal/NSSingularity/TrackPFourier.lean) | 発散ゼロ場での advection 形式の消滅、Leray 乗数の直交性・冪等性・縮小性 | 無限次元 Kato–Ponce 可換子評価は含まない(C2/C3) |
| Kato 定数の有限代数 | `cube_diff_bound`, `am_gm_split`, `shifted_ratio_bound`, `inv_pow_tail_bound`, `g3_of_a4` | [KatoConstant.lean](formal/NSSingularity/KatoConstant.lean) | $12\sqrt{A_4}$ 組立に使う個別の代数補題と格子 tail の telescoping | **不等式 $G_3 \le 12\sqrt{A_4}$ 自体は Lean にない**(`g3_of_a4` は $x\mapsto 12\sqrt{x}$ の単調性のみ)。可換子評価本体は紙上(C2/C3) |
| control ODE 比較 | `riccati_comparison`, `gronwall_variable_coefficient`, `roughEnclosure_solution_unique` | [ControlODE.lean](formal/NSSingularity/ControlODE.lean) | スカラー Riccati 比較、rough enclosure の存在・一意性 | 微分不等式の供給は外部 |
| 継続の貼り合わせ論理 | `glued_continuous`, `exists_continuousOn_Icc_extension` | [ChainAnalysis.lean](formal/NSSingularity/ChainAnalysis.lean) | 連続貼り合わせ、完備空間での端点延長 | 延長関数が方程式を満たすことは含まない |
| 仮定放電の形 | `cond_to_uncond` | [ChainAnalysis.lean](formal/NSSingularity/ChainAnalysis.lean) | 条件付き結論と仮定から結論が従う(**公理非依存**) | 仮定の実例は Lean 側に一切ない |
| $L^3$ 生成・純粋旋回 no-go | `pure_swirl_equality_case`, `viscous_contribution_nonpos`, `transport_eq_one_third_deriv` | [L3Generation.lean](formal/NSSingularity/L3Generation.lean) | 生成恒等式の点ごとの代数、粘性寄与の符号、等号ケース | $\mathbb{R}^3$ 上の積分論(部分積分)は含まない |
| Green / スケーリング / 証明書層 | `greenProfile_radial_laplace_eq_zero`, `physicalTime_lt_blowupTime`, `velocity_radial_error_le` | [GreenAndCascade.lean](formal/NSSingularity/GreenAndCascade.lean), [FiniteTime.lean](formal/NSSingularity/FiniteTime.lean), [CertificateLayer.lean](formal/NSSingularity/CertificateLayer.lean) | 各層の有限不等式・恒等式 | 詳細な限界は登録簿の各項目 |
| 主定理の front 恒等式 (I.1)–(I.4) | `log_bandwidth_derivative_identity`, `covariance_sq_le_variance_mul_action`, `log_bandwidth_derivative_le`, `spectral_front_defect_decomposition` | [SpectralFrontIdentities.lean](formal/NSSingularity/SpectralFrontIdentities.lean) | 有限和の実代数のみ: 中心化分散恒等式、モードごと Cauchy–Schwarz からの共分散上界、平方完成による $G/(2\nu H_1)$ 上界、非負2部分への欠損分解 | **時間積分を含まない。** $\Lambda$ の単調性(主定理の看板)は微分不等式の時間積分を要し、それは形式化していない。$\dot H_0,\dot H_1$ は導関数ではなく ledger 関係式で拘束された実数 |
| Lemma K の有限核 | `sq_sum_abs_le_sum_inv_mul_sum_mul_sq`, `lemmaK_bound` | [SpectralFrontIdentities.lean](formal/NSSingularity/SpectralFrontIdentities.lean) | 重み付き Cauchy–Schwarz $(\sum\lvert c\rvert)^2\le(\sum x^{-1})(\sum xc^2)$ と、それを front 波数の上界へ包む有限代数 | sup ノルム評価と双線形評価は**名前つき仮説**として露出(公理化していない)。候補文書側の命題であり、論文の補題ではない |
| 帯対称性補題 | `sum_offDiagonal_12_eq_zero`, `sum_sq_fst_eq_sum_sq_snd`, `three_mul_sum_sq_fst_eq_sum_normSq`, `band_symmetry` | [BandSymmetry.lean](formal/NSSingularity/BandSymmetry.lean) | 符号反転の対合で非対角和が消え、座標置換で対角和が一致し、$3\sum k_1^2w=\sum\lvert k\rvert^2w$ | 有限対称集合上の和についての言明。格子点の漸近評価は一切含まない |

### B. Certificate-verified results

独立 Python checker が payload を検査します(builder と状態を共有しません)。
**再計算の範囲は行ごとに異なります** — 下表の「What is recomputed」欄が
各 checker の実際の守備範囲であり、それ以外は「記録値どうしの整合性検査」です。
**これらは Lean-verified ではありません。**

**注意:** 下記 3 本のチェーン駆動(`run_track_p_chain.py`、
`run_track_p_chain_h3.py`、`reissue_chains.py`)は CLI 引数を持たず、
**commit 済みの出力ディレクトリを上書きします**。実行前に `git status` が
綺麗であることを確認してください。

| Certificate | Artifact path | Replay / check command | What is recomputed | Trusted implementation | External analytic dependency | Tamper-rejection |
|---|---|---|---|---|---|---|
| Track P 単発スラブ | [`outputs/track_p_slab_v1/`](outputs/track_p_slab_v1/) | `python -m experiments.run_track_p_slab --config configs/track_p_slab.json --output-dir <dir>` | 初期値検査、Picard box、control 不等式の**定数組立**と control 管(**$M_j$ の区間畳み込みと残差上界は builder 側で、checker は再計算しない** — モジュール docstring に明記) | `torus_aposteriori.verify_torus_slab_certificate` | EXT-P1/P2/P3(v1 は条件付き) | あり |
| Track P chain($H^4$) | [`outputs/track_p_chain_v1/`](outputs/track_p_chain_v1/) | `python experiments/run_track_p_chain.py` | 全リンク: box・定数組立・管・Taylor 終端・再中心化点の厳密一致・$\delta$ 漸化式・文言契約 | `torus_chain.verify_chain_certificate` | 同上 | あり |
| Track P chain($n=3$ Kato) | [`outputs/track_p_chain_h3_v1/`](outputs/track_p_chain_h3_v1/) | `python experiments/run_track_p_chain_h3.py` | 上記 + $G_3$ 証明書の再検証 + $C_{\text{kato}}$/$C_{\text{shift}}$ の再計算 | 同上 | 監査済み紙上解析(C2) | あり |
| 監査済み再発行 | [`outputs/track_p_chain_reissued_v2/`](outputs/track_p_chain_reissued_v2/) | `python experiments/reissue_chains.py` | 上記 + 閉鎖メタデータ整合 + 新旧文言の混在拒否 | 同上 | 同上 | あり |
| Kato 定数証明書 | [`outputs/track_p_chain_h3_v1/kato_certificate.json`](outputs/track_p_chain_h3_v1/kato_certificate.json) | `python -m pytest tests/test_kato_constant.py` | $A_4$/$A_6$ を独自の格子ループで再計算、$\sqrt{\ }$ と倍率、単調性 | `kato_constant.verify_kato_certificate` | 可換子評価本体は紙上(C2) | あり |
| $\mathbb{R}^3$ スペクトル圧力 | builder/checker: [`gaussian_spectral_pressure.py`](src/ns_certificate_lab/gaussian_spectral_pressure.py) | `python -m pytest tests/test_gaussian_spectral_pressure.py` | 厳密閉形式 $\nabla p$($\Delta p + \sigma \equiv 0$ の厳密有理自己検証)と $J$ 下界の包含 | 同モジュールの verifier | **なし**(仮定ブロックは空) | あり |
| Track F 有限モード除外 | [`outputs/track_f_finite_mode_scan_v1/`](outputs/track_f_finite_mode_scan_v1/) | `python -m experiments.run_track_f_finite_mode_scan --config configs/track_f_finite_mode_scan.json --output-dir <dir>` | 三線形相殺 $\langle u,(u\cdot\nabla)u\rangle = 0$ を**厳密整数演算**で(浮動小数点なし) | 同 experiment の checker | — | あり |
| snapshot / 時空スラブ | [`outputs/whole_space_gate6_v1/`](outputs/whole_space_gate6_v1/)、[`outputs/tau_continuation_gate7_v1/`](outputs/tau_continuation_gate7_v1/) | [docs/reproducibility.md](docs/reproducibility.md) 参照 | cell 内包含、区間演算の再計算 | `slab_certificate` ほか | 各 payload の hypotheses 欄に明記 | あり(ただし **Gate 7 の保存済み payload は現行 checker のスキーマでは reject される** — 検証には replay での再生成が必要。登録簿 VR-C-008) |

### C1. Conditional certificates

第 9–10 便に発行した v1 証明書は、監査前の文言のまま保存しています。

checker が強制する文言は次の 2 文です(`torus_chain.ALLOWED_WORDING` 逐語):

> conditional PDE certificate assuming EXT-P1/P2/P3
>
> the finite-dimensional Galerkin enclosure and the scalar control ODE are
> verified unconditionally

すなわち **finite-dimensional trajectory・residual・control ODE・chain
accounting は機械検査済みだが、PDE 存在・一意性・tube 包含は明示された EXT
仮定に条件付き**です(payload が無条件と呼ぶのは前者 2 項目のみ。残差と
$\delta$ 漸化式の再計算は checker の実装として行われます)。

対象: [`outputs/track_p_slab_v1/`](outputs/track_p_slab_v1/)、
[`outputs/track_p_chain_v1/`](outputs/track_p_chain_v1/)。これらは履歴上の
成果物として意図的に据え置いており、checker は現在も条件付き文言を強制します。

### C2. Hybrid verified results

<!-- PROMOTION:BEGIN (guard: requires all EXTERNAL_THEOREMS_AUDITED proved:true) -->
次の三層を組み合わせた結果です。**conditional でも Lean-only でもありません。**

1. Lean-verified な有限論理・代数(A 節)
2. certificate-verified な厳密有理計算(B 節)
3. **監査済みだが Lean 形式化されていない**古典解析。**精査の厚みは一様では
   ありません**:
   - **EXT-P1★** — 3 パスの敵対的監査 + 修理 + 再監査 2 名(最も厚い)
   - **EXT-P2-INT + 比較補題 / EXT-P3★ / 系 P3-3** — 第 11 便の新規執筆で、
     再監査 2 名のみ(監査文書 §4 は執筆時点で
     "closed at the paper level, new this turn" と記録し、§6 の再監査で
     closed へ更新)

checker が強制する結論文言は `torus_chain.AUDITED_KIND` 逐語です:

> unconditional PDE certificate modulo audited classical theorems:
> EXT-P1/P2-INT/P3 closed by audited paper proofs
> (docs/research_notes/ext_p1_p2_p3_audit.md); Lean formalisation of the
> infinite-dimensional analysis remains open and is never axiomatised

言い換えれば **machine-checkable certificate modulo explicitly named audited
classical analysis** — 明示された監査済み古典解析を法とする機械検査可能証明書

> **`proved:true` は監査済み紙上証明として閉鎖されたことを意味し、
> Lean 形式化済みを意味しません。**

対象: [`outputs/track_p_chain_reissued_v2/`](outputs/track_p_chain_reissued_v2/)
の再発行チェーン(冒頭の代表的検証結果はこの一つ)。各 payload は閉鎖の
根拠文書と `lean_formalised: false` を必ず携行し、checker は「閉鎖メタデータ
なしの `proved:true`」「条件付き文言との混在」を拒否します。

監査の実施主体は **独立したエージェント監査パス**(AI エージェントによる
複数回の独立監査)であり、**人間の査読ではありません**。監査回数は証明の
根拠ではありません。根拠となるのは、証明文書
[docs/research_notes/ext_p1_p2_p3_audit.md](docs/research_notes/ext_p1_p2_p3_audit.md)
に記録された完全な命題・修正履歴・未解決依存(下記 C3)です。

将来の機械可読形式として、外部定理レコードを次のフィールドへ移行する方針です
(現行 payload の互換性は維持します):

```text
paper_proof_status: audited
lean_formalised: false
external_theorem_status: audited_not_formalised
```
<!-- PROMOTION:END -->

### C3. Open analytic dependencies

現在も未解決の依存です。

| 依存 | 状態 |
|---|---|
| **G-DINI** — EXT-P2 の各時刻 Dini 微分節 | open。ただし**どこからも消費されていない**(積分形 EXT-P2-INT が代替、checker 強制) |
| **HS-5 全空間版** — 離散残差 → 連続 PDE 残差の橋 | open(周期版は Track P で構成的に閉鎖済み) |
| **NT-N1** | open |
| **H3** — 半離散 → 連続の空間補間 | open |
| 無限次元 Kato 型可換子評価の Lean 形式化 | open(紙上は監査済み、有限代数のみ Lean 化) |
| EXT-P1/P2/P3 の PDE 部分の Lean 形式化 | open(**Lean 公理としての挿入は禁止**) |
| **P1 / P1G** — 旧 $J$ 証明書の離散圧力・圧力勾配仮定 | open(Gaussian–Hermite スペクトル経路では消滅したが、旧証明書では checker 強制のまま) |
| $\mathbb{R}^3$ の $J>0$ 候補 | 無条件下界は負のまま。閉鎖には桁違いの計算量が必要 |

### C4. Numerically observed only

証明ではなく浮動小数点による観測です。**いかなる結論の根拠にもしません。**

| 観測 | 内容 | 登録簿 |
|---|---|---|
| 全空間 Gate 4–7 の掃引 | $\tau$/Re 継続、離脱ゲート、振幅継続の数値挙動 | VR-N-001 |
| $L^3$ 減衰の観測 | 純粋旋回族で $\tau>0$ でも減衰が続くこと(定理が保証するのは初期時刻の符号のみ) | VR-N-002 |
| 形状最適化の float 値 | $\mathrm{Re}_{\text{crit}}$、$J$ の float 評価 | VR-N-003 |
| Hou 早期実行 | 解像度制限つき数値観察(公表値の再現主張ではない) | docs/reproducibility.md |

### D. Not proved by this repository

| 範囲外 | 状態 |
|---|---|
| すべての滑らかな初期値に対する大域正則性 | 主張しない |
| 有限時間特異点(存在・構成) | 主張しない |
| Clay 公式命題 (A)–(D) のいずれか | 主張しない |
| certified interval 外での数値軌道の正確性 | numerically observed のみ |
| chain 終了時刻が特異時刻であること | 終了は前登録分類法で**方法の限界**として分類される |
| hybrid result が Lean-only の PDE 定理であること | 該当なし(C2 の三層構成) |
<!-- MCR:END -->

---

## Reproduce / 再現手順

### Quick verification

```bash
git clone https://github.com/HeliCorgi/ns-singularity-certificate-lab.git
cd ns-singularity-certificate-lab
git checkout main

python -m venv .venv
source .venv/bin/activate             # POSIX

python -m pip install --upgrade pip
python -m pip install -e ".[dev]"

python -m pytest                      # 全テスト(証明書 checker と改竄拒否を含む)
python -m pytest tests/test_readme_claims.py   # 本 README の整合性監査

cd formal
lake exe cache get
lake build                            # 期待: Build completed successfully
# 主要定理の #print axioms。標準 3 公理を超える報告・sorryAx・admit が皆無で
# あることを確認(行は一部折り返し、`cond_to_uncond` は「公理非依存」と出る)
lake env lean AxiomAudit.lean
cd ..
```

PowerShell では仮想環境の有効化のみ次に置き換えてください:

```powershell
.\.venv\Scripts\Activate.ps1
```

テスト件数や Lean の job 数は開発とともに変動するため、本 README には固定値を
書きません。**CI([.github/workflows/tests.yml](.github/workflows/tests.yml))が
再実行するのは Python 側(全テストと証明書リプレイ)のみで、
Lean の build と公理監査は CI に含まれていません** — 上記のコマンドで
手元で実行してください。CI は `main` とタグへの push、および全 pull request で
走ります(作業ブランチへの push で二重に起動しないようにするためで、
検証範囲は PR 側で確保されます)。

### Full certificate replay

各実験の完全なコマンド・所要時間・出力・限界の説明は
**[docs/reproducibility.md](docs/reproducibility.md)** にまとめてあります。
概要のみ:

| 系統 | 内容 | 目安 |
|---|---|---|
| manufactured / baseline / 時間収束 | 収束次数と非特異対照 | 各 1 秒未満(実測) |
| Poisson ゲート | 独立 2 実装の相互検証 | 各 1 秒未満(実測) |
| Hou 早期実行・snapshot 監査 | 解像度制限つき数値観察 | 数十分〜 |
| Track F 有限モード除外 | 厳密整数演算の除外証明書 | 数秒 |
| 全空間 Gate 4–7 | 線形/中振幅/離脱ゲートと時空スラブ証明書 | Gate 4 は数秒、Gate 5–6 は数分、Gate 7 は約 30 分 |
| Track P スラブ / chain / $n=3$ Kato / 再発行 | 周期 PDE verification 本体 | スラブ 11 秒(実測)、チェーン 3 本は各 約 2〜2.5 時間 |

### ブランチと固定コミット

| 用途 | 参照 |
|---|---|
| 最新の安定版 | `main`(リポジトリ既定ブランチ) |
| 開発版 | `fable5-mainline`(現在も使用中。`main` へ随時 merge) |
| 本 README 掲載成果の厳密再現 | コミット `17d41df`(代表結果の成果物を生成) |

```bash
git checkout main        # 最新の安定版
git checkout 17d41df     # README 掲載結果の厳密な再現
```

---

## Trust model

| Layer | Trusted components |
|---|---|
| **Lean-verified theorem** | Lean kernel、pinned toolchain と mathlib([formal/lean-toolchain](formal/lean-toolchain)、[formal/lake-manifest.json](formal/lake-manifest.json))、`#print axioms` が表示する標準公理 |
| **Certificate verification** | Python 処理系、独立 checker 実装、入力証明書(payload) |
| **Hybrid PDE result** | 上記 2 層 **+ 明示された監査済み古典解析**([監査文書](docs/research_notes/ext_p1_p2_p3_audit.md)) |
| **Floating-point exploration** | ランタイムとハードウェア。**証明ではない** |
| **External references** | 一次文献の定理内容と、本リポジトリの規約への対応表([kato_h3_constants.md](docs/research_notes/kato_h3_constants.md) §1 ほか) |

公理監査が `[propext, Classical.choice, Quot.sound]` のみを報告することは、
Lean / mathlib の**通常の基盤**であることを意味します。**「無公理」ではありません。**
Lean 側の監査は `lake env lean AxiomAudit.lean`(主要定理の `#print axioms`。
補助補題は個別には列挙されず、主要定理経由で推移的に覆われます)、
`sorry` / `admit` / `axiom` の grep(CI)、および pinned toolchain によります。

陳腐化防止として [tests/test_readme_claims.py](tests/test_readme_claims.py) が
本節の定理名・成果物パス・区分整合・EXT 状態・代表結果の数値・ブランチ記述を
毎 CI で検査します。

---

## 現在できること

### Formal mathematics / Lean

- スラブ合成・連結の有限論理、積分形比較、貼り合わせ・端点延長。
- 有限次元 Galerkin 系の局所存在・一意性(Picard–Lindelöf)。
- 固定有限帯域軌道の破綻不能性と、有限帯域初期値との区別(反例つき)。
- Kato 定数組立、Leray 代数、$L^3$ 生成恒等式の点ごとの代数。
- 全定理の公理監査 → [formal/AxiomAudit.lean](formal/AxiomAudit.lean)、
  規約 → [LEAN4_VERIFICATION_POLICY.md](LEAN4_VERIFICATION_POLICY.md)。

### Rational and interval certificates

- 厳密有理数・有理区間演算による証明書の生成と、**独立 checker による再計算・検査**
  (守備範囲は証明書ごとに §B の表で明示)。
- 改竄拒否テスト(係数改竄、仮定の偽装、免責削除、区間の反転など)。
- Kato 定数 $G_3 \le 12\sqrt{A_4}$ の格子和証明書([導出](docs/research_notes/kato_h3_constants.md))。
- $\mathbb{R}^3$ Gaussian–Hermite 基底の**厳密閉形式スペクトル圧力**(仮定なしの離散包含)。
- sha256 マニフェスト付きの成果物と、前登録された config。

### Periodic PDE verification (Track P)

- 有理 Fourier 初期値 3 族と、厳密有理 Fourier 機構(積和公式・per-mode Leray・
  $H^n$ ノルム・time-Taylor 係数。FFT なし、エイリアシングなし)。
- Galerkin 軌道の**連続** NS 残差を有限三角多項式として厳密計算(周期版 HS-5)。
- 自前導出の control 不等式($H^4$ 粗定数版と $n=3$ Kato 定数版)。
- スラブ連結による certified interval の延長([設計](docs/research_notes/track_p_chain.md))。

### Whole-space numerical gates

- Gate 4–7: 線形楕円ゲート、微分 tail と速度回復、中振幅校正と振幅継続、
  Picard 領域からの離脱と $\tau$/Re 継続、時空スラブ証明書。
- **前登録基準の不合格は不合格として記録**されます([whole_space_transition.md](docs/whole_space_transition.md))。

### Candidate search and rejection

- 有限モード ansatz 族の**除外定理**(厳密整数演算、[Track F](docs/research_notes/track_f_finite_mode_nogo.md))。
- 純粋旋回の $L^3$ 生成 no-go と、混合 meridional+swirl 族の探索。
- 前登録された kill 条件による棄却(閾値は結果を見る前に固定)。

### Reproducibility and audit infrastructure

- 全実験が config 固定・seed 固定・sha256 付き。**CI が再実行するのは
  全テストと 9 本のリプレイのみ**で、Track P chain / $n=3$ / 再発行(各 約 2 時間)と
  Hou 実行は CI 対象外です。テストとリプレイは並列の2ジョブに分かれています
  (最長工程は Gate 7 の約 12 分半)。
- 数値観測と証明を混同しないための STATUS / 証明義務台帳 / 登録簿の分離。
- README 自体の整合性を CI が検査([tests/test_readme_claims.py](tests/test_readme_claims.py))。

---

## 数学的対象

主対象は無外力の Navier–Stokes 方程式

$$
\partial_t u + (u\cdot\nabla)u = -\nabla p + \nu\Delta u,\qquad
\nabla\cdot u = 0,\qquad \nu > 0
$$

です。Track P は周期領域 $\mathbb{T}^3$(Clay 公式命題 (B)/(D) の設定)を、
全空間ゲートは $\mathbb{R}^3$ の軸対称・旋回ありの設定を扱います。規約・導出・
次元・スケーリング・軸の偶奇性は [SPEC.md](SPEC.md) と
[数学的定式化](docs/mathematical_formulation.md)、各式の状態は
[方程式監査](docs/equation_audit.md) にあります。

軸対称レーンの変換変数は

$$u_1 = u^\theta/r,\qquad \omega_1 = \omega^\theta/r,\qquad \psi_1 = \psi^\theta/r$$

で、形式的作用素

$$\mathcal{L}_5 = \partial_{rr} + 3r^{-1}\partial_r + \partial_{zz}$$

はスカラー式の略記であり、物理的な 5 次元流体を意味しません。非圧縮条件と
体積測度は常に 3 次元のものを使います。

---

## リポジトリ構成

| パス | 内容 |
|---|---|
| [`docs/verified_results.md`](docs/verified_results.md) | **検証成果の登録簿**(安定 ID、仮定、再現コマンド、限界) |
| [`docs/reproducibility.md`](docs/reproducibility.md) | **全実験の再現手順**(環境構築、コマンド、所要時間、限界) |
| [`STATUS.md`](STATUS.md) | 現在地と未解決事項の唯一の集約先 |
| [`SPEC.md`](SPEC.md) | 数学的対象、変数、解・特異点の定義 |
| [`LEAN4_VERIFICATION_POLICY.md`](LEAN4_VERIFICATION_POLICY.md) | Lean 4 / mathlib4 による最終検証の必須規約 |
| `formal/` | Lean 4 形式化と公理監査 |
| `docs/research_notes/` | 導出・監査文書(EXT 監査、Kato 定数、Track P、chain 設計ほか) |
| `docs/proof_obligations.md` | 数値候補から反例までの証明義務 |
| `docs/known_obstructions.md` / `docs/threat_model.md` | 既知の非存在定理 / 偽特異点の検出試験 |
| `docs/equation_audit.md` | 符号・係数・境界・同値性の式別監査(E-01–E-33) |
| `src/ns_certificate_lab/` | 数値・証明書・checker の実装 |
| `tests/` | manufactured、round-trip、故障注入、改竄拒否、README 監査 |
| `experiments/` / `configs/` / `outputs/` | 実験ドライバ / 前登録入力 / 成果物 |
| `archive/` / `certificates/` | 履歴 provenance / 将来の明示候補証明書用(現在候補なし) |

---

## 研究上のゲート

既知障害は [docs/known_obstructions.md](docs/known_obstructions.md)、誤検出対策は
[docs/threat_model.md](docs/threat_model.md)、証明までの依存関係は
[docs/proof_obligations.md](docs/proof_obligations.md) にあります。
ニューラルネットを将来使う場合も、保持する候補は明示基底係数へ変換し、
ネットワークを使わない独立残差評価に合格させます。

プロジェクトの現在地と未解決事項は [STATUS.md](STATUS.md) にのみ集約し、
数値的確認と数学的証明を混同しません。

---

## License / ライセンス

[Apache License 2.0](LICENSE)([NOTICE](NOTICE) を含む)。

この選択は任意ではなく既存の宣言と整合させたものです: 本リポジトリの Lean
ファイルは以前から Apache 2.0 を宣言して「the file LICENSE」を参照していました。
また Lean 4 と mathlib4 自体が Apache 2.0 であり、エコシステムの標準です。
特許許諾条項と NOTICE 機構を持つ点も、第三者による検証・引用・再利用を前提と
する本リポジトリに適合します。

対象範囲はソースコード(Python / Lean)、証明書成果物、研究ノートを含む全成果物
です。依存する mathlib4(Apache 2.0)と NumPy(BSD 3-Clause)は再配布物に含まれず、
それぞれのライセンスに従います。

再配布・引用の際は [NOTICE](NOTICE) の scope notice を必ず引き継いでください。
