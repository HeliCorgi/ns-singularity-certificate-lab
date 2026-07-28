# Research notes (imported, 2026-07-29)

外部セッション(ChatGPT)が作成し、2026-07-29 に `fable5-mainline` へ統合した
進捗ノート。**内容は統合時に一切書き換えていない**(ファイル名のみ `docs/`
の小文字 snake_case 規約に合わせた)。実装・テスト・実験スクリプトは通常の
`src/`、`tests/`、`experiments/`、`docs/` へ統合済みで、本ディレクトリは
その経緯と外部由来の主張を原文のまま保存する場所である。

| 現ファイル名 | 元ファイル名 | 内容 |
|---|---|---|
| `critical_l3_obstruction.md` | `CRITICAL_L3_OBSTRUCTION.md` | 臨界 \(L^3\) 障害(ESS 端点正則性定理)と探索目標の改訂 |
| `type_ii_scaling_progress.md` | `TYPE_II_SCALING_PROGRESS.md` | Type-II・異方的スケーリングの条件付き族と早期 Hou への適用 |
| `whole_space_poisson_progress.md` | `WHOLE_SPACE_POISSON_PROGRESS.md` | 自由空間 radial Green solver と低波数壁補正 |

## 読む際の注意

1. ノート内の数値は**古いリポジトリ snapshot** に対する外部計算である
   (`outputs/imported_chatgpt_results/README.md` の provenance 限界を参照)。
   `257x512` は含まれず、現行 `fable5-mainline` での再計算は未実施。
2. ノート内の「テスト N passed」等の記述は外部環境での実行結果である。
   統合後の本リポジトリでの実測は `STATUS.md` に記録する。
3. `critical_l3_obstruction.md` §9 は Lean 形式化の識別子として `F-4`〜`F-7`
   を提案しているが、本リポジトリの `docs/formalization_map.md` では
   `F-4` は既に**証明書の有限次元不等式**に割り当て済みである。統合時に
   ノート本文は変更していないので、**識別子の衝突は未解決**であり、
   実際に形式化へ着手する際に `formalization_map.md` 側で採番し直す。
4. これらのノートは Clay 問題を解決していない。§3 の no-go 定理は
   「一様に \(L^3\) 有界な一スケール再スケーリング候補は \(\mathbb R^3\) の
   有限時間爆発を与えない」という**除外**であり、爆発の構成ではない。
