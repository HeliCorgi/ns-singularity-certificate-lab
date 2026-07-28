# Imported reference outputs (external session, 2026-07-29)

これらのファイルは**このリポジトリの実験ハーネスが生成したものではない**。
外部セッション(ChatGPT)が、当時アップロードされていた**より古い**
リポジトリ snapshot の checkpoint に対して計算した参照出力である。
統合の経緯は `docs/research_notes/` の 3 つのノートを参照。

## provenance の限界(重要)

- `outputs/` の他のディレクトリと異なり、ここには本リポジトリの
  `collect_runtime_provenance()` による runtime provenance ブロック、
  config snapshot、SHA-256 sidecar 一式が**ない**。
  `critical_l3_old_uploaded_snapshot/manifest.json` は外部セッションが
  書いた 2 ファイル分の manifest であり、統合時に両ファイルの SHA-256 が
  一致することは確認した。それ以外の完全性保証はない。
- 入力 checkpoint は**古い snapshot**であり、現行 `fable5-mainline` の
  `outputs/hou_early_time_v1` / `v2_hires` とビット単位の同一性は
  確認していない。特に `257x512` は含まれない。
- したがってこれらは**再実行で置き換えるべき参照値**であり、
  受入検査の根拠にはしない。現行コードで再計算する経路は
  `experiments/analyze_critical_l3.py`、
  `experiments/fit_critical_scaling.py`、
  `experiments/scan_scaling_constraints.py`、
  `experiments/run_free_space_poisson_gate.py`、
  `experiments/run_low_frequency_wall_gate.py`。

## 内容

| ディレクトリ | 内容 | 生成元スクリプト |
|---|---|---|
| `critical_l3_old_uploaded_snapshot/` | 3 解像度 × 5 時刻の臨界 \(L^3\)・shell 分解診断(`critical_l3_snapshots.csv`、`critical_l3_summary.json`、外部 manifest) | `experiments/analyze_critical_l3.py` |
| `type_ii_scaling/` | 条件付き指数制約の走査結果(`feasible_scalings.csv` 全 20181 格子点中 438 点、`summary.json`)と早期 Hou snapshot への fit 結果(`scaling_fit_old_snapshot.json`) | `experiments/scan_scaling_constraints.py`、`experiments/fit_critical_scaling.py` |
| `whole_space_poisson/` | 自由空間 radial Green solver の manufactured gate と低波数壁補正 gate | `experiments/run_free_space_poisson_gate.py`、`experiments/run_low_frequency_wall_gate.py` |

## 主張しないこと

これらは有限周期円柱上の浮動小数点観測であり、\(L^3(\mathbb R^3)\) の
上下界ではない。特異点の証拠でも正則性の証拠でもなく、Clay 問題に対する
いかなる主張の根拠でもない。
