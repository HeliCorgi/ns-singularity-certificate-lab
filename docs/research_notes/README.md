# Research notes

本ディレクトリには 2 種類のノートがある。

**(a) 外部由来(2026-07-29 統合)** — 外部セッション(ChatGPT)が作成し
`fable5-mainline` へ統合した進捗ノート。**内容は統合時に一切書き換えて
いない**(ファイル名のみ `docs/` の小文字 snake_case 規約に合わせた)。
実装・テスト・実験スクリプトは通常の `src/`、`tests/`、`experiments/`、
`docs/` へ統合済みで、本ディレクトリはその経緯と外部由来の主張を原文の
まま保存する場所である。

**(b) 本リポジトリ内で作成** — 当セッションで証明・実装した結果のノート。
外部ノートと混同しないよう表を分けてある。

| 現ファイル名 | 元ファイル名 | 内容 |
|---|---|---|
| `critical_l3_obstruction.md` | `CRITICAL_L3_OBSTRUCTION.md` | 臨界 \(L^3\) 障害(ESS 端点正則性定理)と探索目標の改訂 |
| `type_ii_scaling_progress.md` | `TYPE_II_SCALING_PROGRESS.md` | Type-II・異方的スケーリングの条件付き族と早期 Hou への適用 |
| `whole_space_poisson_progress.md` | `WHOLE_SPACE_POISSON_PROGRESS.md` | 自由空間 radial Green solver と低波数壁補正 |
| `term_balance_progress.md` | `TERM_BALANCE_PROGRESS.md` | 変換 PDE の項別釣合い診断と早期 Hou への適用(2026-07-29 第 2 便) |

## (b) 本リポジトリ内で作成したノート

| ファイル名 | 内容 |
|---|---|
| `track_f_finite_mode_nogo.md` | Track F **固定有限帯域** ansatz の除外定理(2026-07-29 第 3 便)。Lemma 1 の厳密算術検証、Theorem 1 の a priori 上界(Lean `F-6`)、任意の (C)/(D) 反例が満たすべき必要条件 |
| `track_f_shell_constraints.md` | 帯域幅発散候補の必要条件と実現可能指数領域(第 4 便)。Lean `F-16` |
| `green_derivative_tail_bounds.md` | 5 次元 Green 核の微分 tail 上界と 2 つの独立評価経路(第 5 便)。Lean `F-14`/`F-15` |
| `cascade_toy_model.md` | 外力の役割分離と有限 cascade 模型(第 5 便)。低周波のみの外力が間接駆動しうることの模型検査 |

## 読む際の注意

1. ノート内の数値は外部環境での計算である。ただし統合時の照合で、
   `critical_l3_snapshots.csv` と `term_balance.csv` が記録する入力
   checkpoint の SHA-256 は**全 15 行が現行 `fable5-mainline` の
   checkpoint と一致**した(`outputs/imported_chatgpt_results/README.md`)。
   すなわち入力バイトは同一である。一方で `257x512` はどの表にも含まれず、
   解析コードのバージョン同一性も保証されないため、**現行コードでの
   再計算は依然として未実施**である。
2. ノート内の「テスト N passed」等の記述は外部環境での実行結果である。
   統合後の本リポジトリでの実測は `STATUS.md` に記録する。
3. `critical_l3_obstruction.md` §9 は Lean 形式化の識別子として `F-4`〜`F-7`
   を提案しているが、本リポジトリの `docs/formalization_map.md` では
   `F-4`(証明書の有限次元不等式)と `F-5`(Clay 命題定義)が先に割り当て
   済みであった。**この衝突は 2026-07-29 第 3 便で解消した**:
   ノート本文は変更せず、`docs/final_target.md` §4 の登録簿で
   ノートの `F-4`→`F-8`、`F-5`→`F-9`、`F-6`→`F-10`、`F-7`→`F-11` と
   確定した。`F-6`/`F-7` は Track F の Galerkin 除外定理に割り当てた。
   以後は `docs/final_target.md` §4 が唯一の権威である。
4. ノート内で参照される参照出力のパスは、統合時に
   `outputs/imported_chatgpt_results/` 配下へ移した。特に
   `term_balance_progress.md` の「Files」節が挙げる
   `outputs/term_balance_old_snapshot/...` は、実際には
   `outputs/imported_chatgpt_results/term_balance_old_snapshot/...` にある
   (ノート本文は改変していない)。移動後も manifest の SHA-256 は一致する。
5. これらのノートは Clay 問題を解決していない。§3 の no-go 定理は
   「一様に \(L^3\) 有界な一スケール再スケーリング候補は \(\mathbb R^3\) の
   有限時間爆発を与えない」という**除外**であり、爆発の構成ではない。
