# 旧試作リポジトリ限定再利用監査

## 1. 目的と隔離範囲

監査日: 2026-07-27

参照元は
`C:\Users\corgi\Desktop\navier-stokes-singularity-search.zip` である。
ZIP全体の展開、合流、コピーは行わず、archiveから次の6エントリだけを
read-onlyで読んだ。

1. `src/ns_singularity/solvers/poisson.py`
2. `tests/integration/test_poisson.py`
3. `src/ns_singularity/diagnostics/operator_validation.py` のPoisson検査部分
4. `src/ns_singularity/equations/dynamic_rescaling.py`
5. `tests/unit/test_scaling.py`
6. `tests/property/test_finite_inputs.py`

旧`.git`、`.venv`、cache、outputs、設定、manifest、provenance、旧演算子、
旧速度回復、旧動的再スケーリング文書は読込・移植対象にしていない。
本リポジトリの [equation_audit.md](equation_audit.md) の規約と現在の
artifact形式を常に優先する。

この文書は静的な参考実装監査である。旧環境を復元して旧test suiteを実行した
ものではなく、旧コードの正しさを認証しない。

## 2. 旧Poisson試作品の対象

旧実装は一様な有限円柱

\[
0\le r\le R,\qquad z\ \text{periodic}
\]

上で

\[
-\mathcal L_5\psi_1=\omega_1,\qquad
\mathcal L_5=\partial_{rr}+\frac3r\partial_r+\partial_{zz}
\]

を解く疎行列試作品である。外側 \(r=R\) はcaller指定Dirichlet値、軸は
偶正則性を使う。非周期 \(z\) は明示的に `NotImplementedError` とする。
これは全空間Poisson問題でも、領域打切り誤差を評価するsolverでもない。

## 3. 行列式の静的監査

| 項目 | 旧実装 | 監査 |
|---|---|---|
| 方程式の符号 | matrixは \(-\mathcal L_5\)、RHSは内部で \(\omega_1\) | **導出済み: 整合** |
| 周期 \(z\) 行 | 隣接係数 \(-1/\Delta z^2\)、対角へ \(+2/\Delta z^2\) | **導出済み: 整合** |
| 内部radial行 | 外側隣接 \(-[1/\Delta r^2+3/(2r\Delta r)]\)、内側隣接 \(-[1/\Delta r^2-3/(2r\Delta r)]\)、対角 \(2/\Delta r^2\) | **導出済み: 整合** |
| 軸行 | 対角 \(8/\Delta r^2+2/\Delta z^2\)、\(r=\Delta r\) へ \(-8/\Delta r^2\) | **導出済み: 整合** |
| 外側境界行 | \(r=R\) の行をidentityに置換 | **導出済み: Dirichlet行として整合** |
| RHS境界 | 最外行をcallerの境界値へ置換 | **導出済み: identity行と整合** |
| flatten順序 | `radial_index * nz + axial_index` | **確認済み: 行列assembly内では一貫** |

軸行は、偶関数に対する

\[
(\mathcal L_5f)(0,z)=4f_{rr}(0,z)+f_{zz}(0,z),\qquad
f_{rr}(0,z)\approx\frac{2(f_1-f_0)}{\Delta r^2}
\]

から

\[
-(\mathcal L_5f)_0
\approx
\frac{8f_0-8f_1}{\Delta r^2}
\frac{2f_{0,j}-f_{0,j-1}-f_{0,j+1}}{\Delta z^2}
\]

を得るため、係数8と符号は本リポジトリのE-13、E-17に一致する。

## 4. Manufactured solution監査

旧試験は

\[
\psi_1=(1-r^2)^2\cos z
\]

を使う。\(q=1-2r^2+r^4\) とおけば

\[
q_{rr}+\frac3r q_r=-16+24r^2,\qquad
\partial_{zz}(q\cos z)=-q\cos z,
\]

したがって

\[
-\mathcal L_5\psi_1
=(16-24r^2+q)\cos z
\]

であり、旧testの \(\omega_1\) は符号・係数とも整合する。さらに
\(\psi_1(R=1,z)=0\) なのでdefaultの同次外側Dirichlet条件とも整合する。

旧検査は次を含む。

- \(9\to17\to33\) 型のcoupled \(r,z\) refinement。ただしPoisson解誤差に
  実際に使うのは最後の \(17\to33\) だけで、Poissonの観測次数は1個;
- Poisson解のRMS誤差と観測次数;
- discrete residualの \(L^\infty\);
- 非零Dirichlet traceの直接検査;
- 非周期 \(z\) を推測せず拒否する検査。

ただし `poisson_residual` は旧 `equations.operators.l5` を呼び、matrix assembly
と同じ有限差分規約を共有する。小さいdiscrete residualはlinear solveの
自己整合性として有用だが、符号・軸係数の独立検証ではない。解析解への
solution errorとrefinementだけが、この循環から独立した主要検査である。
また、軸行の係数8をmatrixから直接assertするtestはなく、軸専用検査は有限値
確認だけである。非零境界testも最終行の一致を確認するだけで、非零境界を持つ
内部解の精度・収束は測らない。

内部のpointwise central stencilは通常のEuclidean内積では非対称である。
特に \(r_1=\Delta r\) では
\(1/\Delta r^2-3/(2r_1\Delta r)<0\) となり、\(-\mathcal L_5\) matrixの一つの
off-diagonalが正になる。これは偶smooth fieldに対する局所整合性を直ちに
否定しないが、M-matrix性、\(r^3dr\) 重みでのcoercivity、条件数、安定性は
旧検査からは分からない。

## 5. 再利用できる設計上の要点

将来、本リポジトリへ有限円柱用の独立楕円solverを追加する際には、コードを
コピーせず次だけを設計要件として再利用できる。

1. PDE行と外側Dirichlet行を明確に分離する。
2. 軸行を \(4f_{rr}+f_{zz}\) の極限から再導出し、係数8をtestで固定する。
3. solver名とRHS規約に \(-\mathcal L_5\) の負符号を明示する。
4. periodic \(z\) wrapと非対応境界条件の明示拒否を行う。
5. 非零の解析的Dirichlet traceを試験し、同次境界だけに依存しない。
6. 解析manufactured solutionの解誤差を少なくとも3解像度で測る。
7. PDE residualからDirichlet行を分離して報告する。
8. matrix、RHS、reshape順序をmanifestまたはsolver metadataへ記録する。
9. \(r^{-3}\partial_r(r^3\partial_r f)\) のflux/finite-volume離散化とも比較し、
   重み付き安定性を監査する。
10. \(r,z\) を別々にrefineし、複数のaxial Fourier modeと非零解析境界を使う。

## 6. 移植しないものと理由

- `poisson.py` 本体: SciPy、旧grid、旧boundary helper、旧`l5`へ結合しており、
  現在のNumPy最小依存と独立監査経路へそのまま適合しない。
- 旧 `poisson_residual`: matrixと旧円柱operatorが規約を共有し、独立性が弱い。
  `include_outer_boundary=True` もidentity境界行の残差ではなく外端での
  differential residualを返すため、境界診断として意味が曖昧である。
- 旧operator/速度回復一式: 現在の監査済み規約を上書きし、同じ誤りの相殺を
  招くため禁止する。
- `operator_validation.py` 全体: 現在のconvergence、diagnostics、provenance
  schemaと重複し、旧manifestへ再結合するため移植しない。
- 旧dynamic-rescaling文書・実装: 読んだPythonはscaling exponentの整合検査
  だけで、動的再スケーリング方程式や時間発展solverではない。本リポジトリの
  E-21と [future_search.md](future_search.md) の方が対象と義務を明確にする。
- Hypothesis property test: 「有限定数へ \(\mathcal L_5\) を作用させると有限な
  0になる」という着想は有用だが、この1件のために依存を増やさない。必要なら
  現行pytestのparameterized testとして独立に書き直す。

## 7. 静的監査だけでは未確認の点

- 旧archiveの依存環境を起動していないため、旧testsの実行結果は **未確認**。
- 指定外である旧 `AxisymmetricGrid`、`outer_dirichlet_values`、旧`l5` の
  実装は読んでいない。そのためgrid endpoint規約、boundary shape変換、
  residual実装の詳細は **未確認**。
- 疎行列の条件数、roundoff増幅、factorization再利用、収束率の別実装再現は
  **未確認**。
- \(r^3dr\) 重みに関する離散coercivity、一意可解性、M-matrix性は
  **未確認**。
- 非零Dirichlet値を持つ内部解の精度と収束は **未確認**。
- 有限円柱から全空間への領域打切り誤差、外側境界感度、Green tailは
  **未確認**。
- 区間演算、外向き丸め、solver誤差の厳密上界はない。

## 8. 現在の採否

**結論: 旧Poisson solverは移植していない。設計監査後、現行規約だけから
独立solverを新規実装した。**

先にCartesian原始変数監査と固定格子時間収束を完了した後、
`src/ns_certificate_lab/poisson.py` にNumPyだけの別実装を追加した。旧試作の
pointwise sparse matrix、SciPy、旧grid、旧boundary helper、旧`l5`はコピーも
importもしていない。新実装は
\(r^{-3}\partial_r(r^3\partial_r\psi_1)\) のcontrol-volume flux、
\(z\) のFFT、各modeのローカルThomas solveを用いる。

上記設計要点のうち、PDE/境界行の分離、軸係数8、全体符号、周期wrap、
非零解析Dirichlet trace、3解像度・複数modeのmanufactured convergence、
代数残差と独立物理空間残差の分離、metadata、異なるflux離散化は
実装・自動テスト済みである。

一方、\(r,z\) の完全に分離したrefinement、領域半径感度、全空間tail、
重み付きcoercivity、区間演算は **未確認** のままである。従って新solverも
有限円柱上の浮動小数点試作品であり、旧実装または連続体Poisson問題の
正しさを認証するものではない。
