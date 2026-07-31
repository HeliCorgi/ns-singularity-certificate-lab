# Navier–Stokes Millennium Project — New Session Bootstrap

## 0. この文書の使い方

この文書だけを新しいClaude Code / Fable 5セッションへ渡す。

ユーザーは追加の会話履歴、ZIP、patch、個別ファイルを貼らない。必要な情報は、以下のGitリポジトリとこの文書から取得すること。

- Repository: `https://github.com/HeliCorgi/ns-singularity-certificate-lab`
- Working branch: `fable5-mainline`

最初にリポジトリを開くかcloneし、必ず`fable5-mainline`をcheckoutすること。

```bash
git clone -b fable5-mainline https://github.com/HeliCorgi/ns-singularity-certificate-lab.git
cd ns-singularity-certificate-lab
git status
git branch --show-current
git log --oneline -10
```

既知の最終確認commitは`83f5301`だが、これは参考情報にすぎない。GitHub上の現在の`origin/fable5-mainline`を真実の基準にし、HEADを必ず確認すること。

新しいセッションでは、過去チャットの内容を推測しない。不明点はリポジトリ内のコード、テスト、文書、出力から再構成すること。

---

# 1. 最終目的

このプロジェクトの目的は、Clay Mathematics Instituteの3次元非圧縮Navier–Stokes問題について、公式ステートメント(A)〜(D)のいずれか一つを完全に証明することである。

最終成果は、次のいずれかでなければならない。

1. `R^3`上で、任意の滑らか・発散ゼロ・十分速く減衰する初期値から、外力なしの滑らかな有限エネルギー解が全時間存在することを証明する。
2. 3次元周期領域で、任意の滑らか・発散ゼロ初期値から、外力なしの滑らかな解が全時間存在することを証明する。
3. `R^3`上で、Clay条件を満たす滑らかな初期値と滑らかな外力を明示し、大域滑らかな有限エネルギー解が存在しないことを証明する。
4. 3次元周期領域で、滑らかな初期値と滑らかな外力を明示し、大域滑らかな周期解が存在しないことを証明する。

現在の主攻撃経路は、公式要件より強い次の反例である。

> `R^3`上、外力なし、滑らか・発散ゼロ・有限エネルギーの初期値から有限時間特異点を構成する。

ただし、賞金獲得が目的であるため、滑らかな外力を使える公式(C)/(D)も並行して検討する。

---

# 2. 重要な認識規則

以下を絶対に混同しないこと。

- 数値的増幅は証明ではない。
- テスト成功はPDE定理ではない。
- `lake build`成功はClay問題の証明ではない。
- Hou論文の図との類似は反例ではない。
- 小さい浮動小数点残差は厳密解の存在証明ではない。
- Lean 4で証明された有限次元恒等式は、PDEの大域存在・特異点証明ではない。
- 有限円柱・周期`z`の計算は、そのまま`R^3`の候補ではない。
- 一様に`L^3`有界な標準一尺度再スケーリング軌道は、全空間有限時間特異点の最終候補にならない。
- 否定的結果や候補棄却も正式な研究成果として保存する。

未知の特異点、Clay問題の解決、Hou公表値への収束を主張してはならない。

---

# 3. 最初に読むファイル

次の順序で読むこと。ファイル名や配置が変わっている場合は、意味的に対応する最新ファイルを探す。

1. `STATUS.md`
2. `README.md`
3. `AGENTS.md`
4. `PLAN.md`
5. `docs/equation_audit.md`
6. `docs/hou_setup_audit.md`
7. `docs/whole_space_transition.md`
8. `docs/numerical_stability_audit.md`
9. `docs/formalization_map.md`
10. `docs/proof_obligations.md`
11. `docs/threat_model.md`
12. `formal/README.md`
13. `docs/research_notes/`
14. `outputs/imported_chatgpt_results/`

その後、次の実装を確認する。

- 軸対称旋回方程式
- 有限円柱Poisson solver群
- 独立Cartesian検査
- Hou早期実行
- streaming gate
- 積分器比較
- 臨界`L^3`診断
- shell分解
- Type-II / 異方的スケーリング制約
- scaling fit gate
- 自由空間radial Green Poisson solver
- 低周波・ゼロモード壁補正
- Lean 4形式化
- 各テスト

---

# 4. 現在までに確認されている研究状況

リポジトリHEADを優先するが、既知の状態は次の通り。

## 数値基盤

- 数百件規模のPythonテストが存在する。
- 最後に報告された統合直後の結果は`747 passed, 1 skipped`。
- Leanでは`lake build`成功が報告されている。
- `sorry`、`admit`、新規の核心的`axiom`はないと報告されている。
- ただし現在のHEADで必ず再検証すること。

## Hou早期計算

- Houの有限円柱・周期`z`・軸対称旋回計算の早期区間を部分再現している。
- 渦度増幅は解像度とともに増加したが、未収束。
- frontは十分な格子点数で解像されていない。
- グローバル物理`L^3`は早期区間で減少した。
- 臨界密度の幅は縮小せず、早期区間では拡大した。
- Type-II / 異方的スケーリングfitは早期snapshotでは不合格。
- 多尺度shell増殖も早期snapshotでは観測されていない。
- したがって早期Hou現象は反例候補へ昇格していない。

## 臨界`L^3`障害

標準一尺度再スケーリング

```text
u(x,t) = L(t)^(-1) U((x-x*(t))/L(t), s)
```

では、

```text
||u(t)||_L3 = ||U(s)||_L3
```

が成り立つ。

一様に`L^3`有界な全空間軌道は、端点`L∞_t L3_x`正則性により有限時間特異点を作れない。

したがって全空間反例には、少なくとも次のいずれかが必要。

- Type-II速度増幅
- 異方的集中
- 再スケーリング後`L^3`の増大
- 多尺度カスケード
- 外側tailによるグローバル臨界ノルム増大
- 単一グローバル再スケーリングでは表せない機構

## 全空間移行

- 周期`z`の有限円柱壁感度が小さいことは、`R^3`の壁独立性を意味しない。
- 短周期では最小非零波数が大きく、radial elliptic tailが指数的に小さくなる。
- 低波数・ゼロモードでは壁補正が大きくなりうる。
- 自由空間radial Green Poisson solverと壁補正式の実装が存在するはずなので確認すること。
- 非周期`z`を含む完全な全空間非線形時間発展は、まだ核心的未完了項目である。

## Lean 4

既存の形式証明は、有限物理時間条件、速度回復と発散ゼロなどの中間結果が中心。

Lean 4は次に使用する。

- Clay命題の忠実な形式化
- スケーリング指数条件
- 候補クラスの排除
- 特異項相殺条件
- 初期値の発散ゼロ・滑らかさ・有限エネルギー
- 有限次元証明書の検査
- 最終論理接続

Lean 4だけで大規模PDE候補探索を行わない。

---

# 5. 研究を二本立てにする

## Track U — 外力なし全空間反例

目標:

> 明示的な滑らか・発散ゼロ・有限エネルギー初期値から、有限時間で延長不能な`R^3` Navier–Stokes解を構成する。

必要な段階:

1. 非周期`z`と自由空間楕円回復を持つ全空間軸対称solver。
2. Cartesianで滑らかなコンパクト台または急減衰初期値。
3. 空間・時間・領域・積分器の独立収束。
4. グローバル`L^3`、shell別`L^3`、Type-II指数、PDE項釣合い。
5. 候補となる非定常再スケーリング軌道。
6. 区間残差、tail bound、逆作用素ノルム、非線形Lipschitz bound。
7. Newton–Kantorovich、radii polynomial、縮小写像等による真のPDE軌道の存在。
8. 滑らかな初期値から候補軌道へ入る接続。
9. 有限物理時間とノルム発散。
10. Lean 4による最終論理鎖。

## Track F — 滑らかな外力を逆設計する公式(C)/(D)反例

目標:

> 特異ansatz`u,p`を設計し、
>
> `f = ∂t u + (u·∇)u - νΔu + ∇p`
>
> が特異時刻を越えて滑らかになるように主要特異項を相殺する。

最低条件:

- `u`は`t<T`で滑らか。
- `u`は有限時刻`T`で滑らかに延長不能。
- 初期値は滑らか・発散ゼロ。
- `f`は全時間滑らか。
- `R^3`ならClayの減衰条件、周期領域なら周期性。
- 局所滑らか解の一意性を使い、同じデータの別の大域滑らか解を排除。
- 最終的にClayの(C)または(D)へ接続。

探索方法:

- 低次Fourier ansatz
- 自己相似・異方的ansatz
- 時間依存有限モードansatz
- 記号計算による特異次数相殺
- Python/SymPy等で有限次元候補探索
- 有理化された代数条件をLean 4で検証

注意:

単に任意の特異`u`を選んで残差を`f`と定義するだけでは不十分。`f`自身が全時間滑らかでなければClay反例にならない。

---

# 6. 新しいセッションで最初に行う作業

最初のセッションは、計画だけで終わらせない。ただし、長時間の未知シミュレーションを無制限に開始しない。

## Step 1 — 現在状態の確定

```bash
git status
git branch --show-current
git log --oneline -10
python -m pytest
cd formal
lake build
cd ..
```

さらに以下を監査する。

```bash
grep -RInE '\bsorry\b|\badmit\b|^[[:space:]]*axiom ' formal
```

テスト数、skip理由、Lean依存公理を記録する。

## Step 2 — 最終定理の固定

`docs/final_target.md`を作成または更新し、次を一つの依存グラフにする。

- Clay A〜D
- Track Uの最終定理
- Track Fの最終定理
- 全証明義務
- 数学的に閉じた義務
- Leanで閉じた義務
- 数値的観測だけの義務
- 区間証明が必要な義務
- 未着手義務

## Step 3 — 研究上の次の一手

リポジトリHEADを確認したうえで、次のどちらかを選ぶ。

### 優先候補A: Track Fの有限次元相殺探索

選ぶ条件:

- 全空間非周期solverがまだ大規模未完成。
- 外力あり反例の十分条件がまだ形式化されていない。
- 低次ansatzの記号探索を短時間で実装可能。

実施内容:

1. divergence-free周期Fourier ansatzを定義。
2. 有限時刻で一部係数が発散する時間依存を置く。
3. 残差`f`の全特異次数を記号的に展開。
4. `f`が滑らかになるための代数相殺条件を生成。
5. 低次モードで解が存在するか探索。
6. 存在しなければ、そのansatz族のno-go条件を証明・Lean化。
7. 存在すれば、初期値、`u`、`p`、`f`を完全に明示し、局所一意性への接続義務を列挙。

### 優先候補B: Track Uの非周期`z`自由空間入口

選ぶ条件:

- 既存free-space radial solverが安定している。
- 非周期`z`の線形楕円検証を小さい範囲で完成可能。

実施内容:

1. 非周期`z`のcompact-support manufactured source。
2. 自由空間`-L5`解の独立参照。
3. `Rmax,Zmax`独立拡大。
4. 低波数・ゼロモード・tail bound。
5. Cartesian復元後のdivergence/curl検査。
6. 線形gate合格後のみ非線形時間発展へ接続。

最初のセッションでは、AかBのどちらかについて、コード、テスト、文書、具体的な数学結果を一つ以上完成させること。

---

# 7. 実装・証明の品質規則

## 数値

- 浮動小数点結果を証明と呼ばない。
- 候補生成と検証を同一実装だけで完結させない。
- 空間、時間、領域、積分器を独立に変える。
- 受入値は全step streamingで監視する。
- 結果を見てからfit窓や閾値を変更しない。
- 不合格結果を隠さない。

## Lean 4

最終証明経路では禁止:

- `sorry`
- `admit`
- 核心を仮定するproject-specific `axiom`
- 浮動小数点値をそのまま厳密実数として使用
- Pythonの`assert`を証明とみなす
- 外部solverの成功表示だけを証明とする

各重要定理について`#print axioms`を記録する。

## Git

- 作業前後に`git status`を確認。
- `.venv`、cache、巨大な不要中間出力をcommitしない。
- 結果を再現するconfig、seed、hashを保存。
- 作業後にcommitし、`origin/fable5-mainline`へpushする。
- pushできなければlocal commit hashと実行すべきpushコマンドを報告。

---

# 8. 完了報告の形式

毎回、次を明確に分ける。

1. 数学的に証明したこと
2. Lean 4で証明したこと
3. 区間演算で証明したこと
4. 数値的に観測したこと
5. 棄却した候補
6. 新しく作った候補
7. 変更ファイル
8. テスト結果
9. Lean buildと`#print axioms`
10. Clay最終定理まで残る証明義務
11. 次の最小の一手
12. commit hashとpush結果

「進展した」「面白い」などの曖昧な表現だけで終了しない。

---

# 9. セッション開始命令

この文書を読み終えたら、質問して停止せず、次を実行する。

1. `fable5-mainline`の現在状態を確定する。
2. 指定されたファイルを読む。
3. 全テストとLean buildを確認する。
4. `docs/final_target.md`を整備する。
5. Track FまたはTrack Uから、現HEADに対して最も短く具体的に前進できる作業を一つ選ぶ。
6. コード、テスト、文書、数学結果を完成させる。
7. 全テストを再実行する。
8. commitしてpushする。
9. 上記の完了報告形式で報告する。

目的は診断数を増やすことではない。

> Clay公式(A)〜(D)のいずれかの完全証明へ向かう、候補構成または証明義務を毎回少なくとも一つ前進させること。
