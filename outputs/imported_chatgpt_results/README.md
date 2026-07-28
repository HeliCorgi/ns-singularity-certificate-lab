# Imported reference outputs (external session, 2026-07-29)

これらのファイルは**このリポジトリの実験ハーネスが生成したものではない**。
外部セッション(ChatGPT)が、当時アップロードされていたリポジトリ snapshot の
checkpoint に対して計算した参照出力である。統合の経緯は
`docs/research_notes/` の 4 つのノートを参照。

## 入力の同一性(2026-07-29 に照合済み)

`critical_l3_snapshots.csv` と `term_balance.csv` は入力 checkpoint の
`checkpoint_sha256` を各行に記録している。統合時に、記録された全 15 行
(65×128 / 129×256 / 193×384 の 5 時刻)を現行 `fable5-mainline` の
`outputs/hou_early_time_v1` / `hou_early_time_v2_hires` の実バイトと照合し、
**両表とも 15/15 が一致**した(mismatch 0、missing 0)。

したがってこれらの診断は、外部環境で走ったにもかかわらず、**現在
コミットされている checkpoint とビット単位で同じ入力**に対する計算である。
`checkpoint` 列が外部実行環境の絶対パス(`/tmp/ns-old/...`、`/tmp/nswork/...`)
を記録している点は表示上の差異にすぎない。

**ただし `257x512` は 4 解像度目として存在するが、どの表にも含まれていない。**
臨界 \(L^3\)・項別釣合いのラダーは 3 解像度どまりである。

## provenance の限界(重要)

- `outputs/` の他のディレクトリと異なり、ここには本リポジトリの
  `collect_runtime_provenance()` による runtime provenance ブロック、
  config snapshot、SHA-256 sidecar 一式が**ない**。
  各 `manifest.json` は外部セッションが書いた 2 ファイル分の manifest で
  あり、統合時に記載 SHA-256 が実ファイルと一致することは確認した。
  それ以外の完全性保証はない。
- 入力が同一でも、**解析コードのバージョンが現行 `src/` と同一である保証は
  ない**(外部セッションは会話内の ZIP を基盤にしている)。数値を引用する
  前に現行コードで再計算すること。再計算経路:
  `experiments/analyze_critical_l3.py`、
  `experiments/fit_critical_scaling.py`、
  `experiments/scan_scaling_constraints.py`、
  `experiments/run_free_space_poisson_gate.py`、
  `experiments/run_low_frequency_wall_gate.py`、
  `experiments/analyze_term_balance.py`。
- したがってこれらは**再実行で置き換えるべき参照値**であり、
  受入検査の根拠にはしない。

## 内容

| ディレクトリ | 内容 | 生成元スクリプト |
|---|---|---|
| `critical_l3_old_uploaded_snapshot/` | 3 解像度 × 5 時刻の臨界 \(L^3\)・shell 分解診断(`critical_l3_snapshots.csv`、`critical_l3_summary.json`、外部 manifest) | `experiments/analyze_critical_l3.py` |
| `type_ii_scaling/` | 条件付き指数制約の走査結果(`feasible_scalings.csv` 全 20181 格子点中 438 点、`summary.json`)と早期 Hou snapshot への fit 結果(`scaling_fit_old_snapshot.json`) | `experiments/scan_scaling_constraints.py`、`experiments/fit_critical_scaling.py` |
| `whole_space_poisson/` | 自由空間 radial Green solver の manufactured gate と低波数壁補正 gate | `experiments/run_free_space_poisson_gate.py`、`experiments/run_low_frequency_wall_gate.py` |
| `term_balance_old_snapshot/` | 15 checkpoint(3 解像度 × 5 時刻)の項別釣合い診断(`term_balance.csv`、`term_balance.json`、外部 manifest)。統合時に manifest の SHA-256 2 件が一致することを確認済み | `experiments/analyze_term_balance.py` |

## 主張しないこと

これらは有限周期円柱上の浮動小数点観測であり、\(L^3(\mathbb R^3)\) の
上下界ではない。特異点の証拠でも正則性の証拠でもなく、Clay 問題に対する
いかなる主張の根拠でもない。
