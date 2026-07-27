# References and source policy

この文書は、本プロジェクトで用いる定理・問題設定・数値的警告の出典台帳である。数学的主張には、可能な限り原論文、著者公開版、査読版、または Clay Mathematics Institute の公式問題文を用いる。検索結果の要約、ブログ、百科事典は根拠にしない。

最終照合日: 2026-07-27

## 状態ラベル

- **本文確認済み**: 論文本文または著者公開版で、ここで使う仮定と結論を確認した。
- **書誌・要旨確認済み**: 出版者または著者ページで書誌情報と要旨を確認したが、定理の全証明・全補助仮定までは再照合していない。
- **注意**: 適用範囲や版について、本文中で明示的に限定する必要がある。

## 公式問題設定・弱解・部分正則性

### [CMI-NS] Charles L. Fefferman, “Existence and Smoothness of the Navier–Stokes Equation”

- 種別: Clay Mathematics Institute 公式問題記述
- URL: <https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf>
- 用途: \(\mathbb R^3\) と \(\mathbb R^3/\mathbb Z^3\) における問題の正確な定式化、滑らかで発散ゼロの初期データ、有限エネルギー・急減少条件。
- 状態: **本文確認済み**

### [Leray1934] Jean Leray, “Sur le mouvement d’un liquide visqueux emplissant l’espace”

- Acta Mathematica 63 (1934), 193–248.
- DOI: <https://doi.org/10.1007/BF02547354>
- 公開版: <https://warwick.ac.uk/fac/sci/maths/people/staff/james_robinson/lf/leray.pdf>
- 用途: 全空間での有限エネルギー弱解、エネルギー不等式、自己相似変数、正則時刻と潜在的特異時刻。
- 状態: **本文確認済み**

### [CKN1982] Luis Caffarelli, Robert Kohn, Louis Nirenberg, “Partial Regularity of Suitable Weak Solutions of the Navier–Stokes Equations”

- Communications on Pure and Applied Mathematics 35 (1982), 771–831.
- DOI: <https://doi.org/10.1002/cpa.3160350604>
- 出版者PDF: <https://onlinelibrary.wiley.com/doi/pdf/10.1002/cpa.3160350604>
- 用途: suitable weak solution、局所エネルギー不等式、特異集合の一次元放物型 Hausdorff 測度がゼロであること。
- 状態: **本文・書誌確認済み**

## 軸対称方程式・極条件の監査

### [HouLi2008] Thomas Y. Hou, Congming Li, “Dynamic Stability of the 3D Axi-symmetric Navier–Stokes Equations with Swirl”

- Communications on Pure and Applied Mathematics 61 (2008), 661–697.
- DOI: <https://doi.org/10.1002/cpa.20212>
- 著者プレプリント: <https://arxiv.org/abs/math/0608295>
- 用途: 円柱成分、\(u_1,\omega_1,\psi_1\) 系、軸 Taylor 展開、流れ関数からの速度回復。
- 状態: **本文確認済み**（詳細は `docs/equation_audit.md`）

### [LiuWang2009] Jian-Guo Liu, Wei-Cheng Wang, “Characterization and Regularity for Axisymmetric Solenoidal Vector Fields with Application to Navier–Stokes Equation”

- SIAM Journal on Mathematical Analysis 41 (2009), 1825–1850.
- DOI: <https://doi.org/10.1137/080739744>
- 公開版: <https://archive.ymsc.tsinghua.edu.cn/pacm_download/200/8347-Liu_Wang_SIMA_2009.pdf>
- 用途: 滑らかな軸対称 solenoidal vector field の極条件、偶奇性、vorticity–stream 形式と三次元 primitive 形式の同値性。
- 状態: **本文確認済み**

### [HouLiuWang2018] Thomas Y. Hou, Pengfei Liu, Fei Wang, “Global Regularity for a Family of 3D Models of the Axi-symmetric Navier–Stokes Equations”

- Nonlinearity 31 (2018), 1940–1954.
- DOI: <https://doi.org/10.1088/1361-6544/aaaa0b>
- プレプリント: <https://arxiv.org/abs/1708.07536>
- 用途: 監査する正規化軸対称系と Navier–Stokes スケーリングの独立照合。
- 状態: **本文確認済み**

## 継続・正則性判定

### [Serrin1962] James Serrin, “On the Interior Regularity of Weak Solutions of the Navier–Stokes Equations”

- Archive for Rational Mechanics and Analysis 9 (1962), 187–195.
- DOI: <https://doi.org/10.1007/BF00253344>
- 用途: Ladyzhenskaya–Prodi–Serrin 型の速度の混合 \(L_t^qL_x^p\) 正則性条件。
- 状態: **書誌・定理内容確認済み**
- 注意: 現在通常使われる等号ケースと全大域版は Prodi、Ladyzhenskaya、後続研究を合わせた定式化である。端点 \((p,q)=(3,\infty)\) は [ESS2003] を使う。

### [Prodi1959] Giovanni Prodi, “Un teorema di unicità per le equazioni di Navier–Stokes”

- Annali di Matematica Pura ed Applicata 48 (1959), 173–182.
- DOI: <https://doi.org/10.1007/BF02410664>
- 用途: Prodi–Serrin 条件の古典的起源、弱解の一意性・正則性。
- 状態: **書誌確認済み**

### [ESS2003] Luis Escauriaza, Gregory Seregin, Vladimír Šverák, “\(L_{3,\infty}\)-Solutions of the Navier–Stokes Equations and Backward Uniqueness”

- Russian Mathematical Surveys 58(2) (2003), 211–250.
- DOI: <https://doi.org/10.1070/RM2003v058n02ABEH000609>
- 公開版: <https://www.mathnet.ru/eng/rm609>
- 用途: \(u\in L^\infty(0,T;L^3(\mathbb R^3))\) という速度の臨界端点条件が正則性を与えること。
- 状態: **本文・要旨確認済み**
- 注意: 論文名の \(L_{3,\infty}\) は空間 \(L^3\)、時間 \(L^\infty\) の混合ノルムを表す。Lorentz 空間 \(L^{3,\infty}_x\) ではない。

### [BKM1984] J. Thomas Beale, Tosio Kato, Andrew Majda, “Remarks on the Breakdown of Smooth Solutions for the 3-D Euler Equations”

- Communications in Mathematical Physics 94 (1984), 61–66.
- DOI: <https://doi.org/10.1007/BF01212349>
- 用途: \(\int_0^T\|\omega(t)\|_\infty\,dt\) による古典的な渦度継続判定の原型。
- 状態: **本文・書誌確認済み**
- **注意**: 原論文の定理は Euler 方程式に対するものである。本プロジェクトで使う Navier–Stokes の継続判定の直接出典は [KatoPonce1988] とし、「BKM原論文が Navier–Stokes を直接証明した」とは記述しない。

### [KatoPonce1988] Tosio Kato, Gustavo Ponce, “Commutator Estimates and the Euler and Navier–Stokes Equations”

- Communications on Pure and Applied Mathematics 41(7) (1988), 891–907.
- DOI: <https://doi.org/10.1002/cpa.3160410704>
- 出版者PDF: <https://onlinelibrary.wiley.com/doi/pdf/10.1002/cpa.3160410704>
- 用途: \(1<p<\infty\), \(s>n/p+1\) の Bessel-potential 空間
  \(L_s^p(\mathbb R^n)=(I-\Delta)^{-s/2}L^p(\mathbb R^n)\) の古典解に対する局所適切性と、Euler および Navier–Stokes の
  \(\int_0^T\|\omega(t)\|_{L^\infty}\,dt\) 型継続判定。
- 状態: **書誌・定理内容確認済み**
- 注意: 本プロジェクトでは \(n=3\), \(\nu>0\) の Navier–Stokes 版を使う。
  \(p=2\) なら \(L_s^2=H^s\) で、閾値は \(s>5/2\) となる。関数空間の
  閾値 \(s>n/p+1\) を落として、任意の弱解に対する判定のように引用しない。

### [BdV1995] Hugo Beirão da Veiga, “A New Regularity Class for the Navier–Stokes Equations in \(\mathbb R^n\)”

- Chinese Annals of Mathematics, Series B 16 (1995), 407–412.
- 出版誌アーカイブ・要旨・PDF: <https://camath.fudan.edu.cn/cambcn/ch/reader/view_abstract.aspx?file_no=16B401&flag=1>
- 用途: \(\nabla u\)、したがって \(1<p<\infty\) で Riesz 変換を介した渦度の臨界混合ノルムによる正則性条件。
- 状態: **要旨・定理内容確認済み**
- 注意: 端点 \(L_t^1L_x^\infty\) は別扱いとし、本文では安全な範囲 \(3/2<p<\infty\) を用いる。

## 後方自己相似解と Type I の障害

### [NRS1996] Jiří Nečas, Michael Růžička, Vladimír Šverák, “On Leray’s Self-Similar Solutions of the Navier–Stokes Equations”

- Acta Mathematica 176 (1996), 283–294.
- DOI: <https://doi.org/10.1007/BF02551584>
- 公開版: <https://archive.ymsc.tsinghua.edu.cn/pacm_download/117/6533-11511_2006_Article_BF02551584.pdf>
- 用途: Leray の後方自己相似プロファイル方程式の弱解 \(U\in L^3(\mathbb R^3)\) は \(U=0\) であること。
- 状態: **本文確認済み**

### [Tsai1998] Tai-Peng Tsai, “On Leray’s Self-Similar Solutions of the Navier–Stokes Equations Satisfying Local Energy Estimates”

- Archive for Rational Mechanics and Analysis 143 (1998), 29–51.
- DOI: <https://doi.org/10.1007/s002050050099>
- 著者公開版: <https://personal.math.ubc.ca/~ttsai/publications/leray.pdf>
- 書誌: <https://cir.nii.ac.jp/crid/1360292621610149248>
- 用途: 後方自己相似弱解が有限局所エネルギー条件を満たす場合の非存在、ならびに \(U\in L^q\), \(3<q\leq\infty\) の分類。
- 状態: **本文・定理内容確認済み**
- 正誤表: Tai-Peng Tsai, Erratum, Archive for Rational Mechanics and Analysis 147 (1999), 363. <https://personal.math.ubc.ca/~ttsai/publications/leray-erratum.pdf>。比較関数の積分に \(s^{-2}\) を補う訂正で、主要結論は変更されない。

### [CSTY2008] Chiun-Chuan Chen, Robert M. Strain, Tai-Peng Tsai, Horng-Tzer Yau, “Lower Bound on the Blow-Up Rate of the Axisymmetric Navier–Stokes Equations”

- International Mathematics Research Notices 2008, article rnn016.
- DOI: <https://doi.org/10.1093/imrn/rnn016>
- プレプリント: <https://arxiv.org/abs/math/0701796>
- 用途: 軸対称強解について、\(|v(x,t)|\le C_*(r^2-t)^{-1/2}\) 型の尺度不変上界が時刻 \(0\) の正則性を与えること。
- 状態: **本文・要旨確認済み**

### [CSTY2009] Chiun-Chuan Chen, Robert M. Strain, Tai-Peng Tsai, Horng-Tzer Yau, “Lower Bounds on the Blow-Up Rate of the Axisymmetric Navier–Stokes Equations II”

- Communications in Partial Differential Equations 34(3) (2009), 203–232.
- DOI: <https://doi.org/10.1080/03605300902793956>
- 著者プレプリント: <https://arxiv.org/abs/0709.4230>
- 用途: 軸対称強解について、\(|v|\le C_*|t|^{-1/2}\) または \(|v|\le C_*r^{-1+\varepsilon}|t|^{-\varepsilon/2}\) が時刻 \(0\) の正則性を与えること。
- 状態: **本文・要旨確認済み**

## 軸対称流

### [Ladyzhenskaya1968] O. A. Ladyzhenskaya, “On Unique Solvability ‘in the Large’ of the Three-Dimensional Cauchy Problem for the Navier–Stokes Equations with Axial Symmetry”

- Zapiski Nauchnykh Seminarov LOMI 7 (1968), 155–177.
- 公開版: <https://www.mathnet.ru/eng/znsl2240>
- 用途: 旋回なし軸対称 Navier–Stokes 流の古典的大域正則性。
- 状態: **書誌確認済み**
- 注意: 原文はロシア語であり、現代的 Sobolev 空間での完全な言い換えは本プロジェクトでは未照合。適用時には [UY1968] または現代的再証明も併記する。

### [UY1968] M. R. Ukhovskii, V. I. Yudovich, “Axially Symmetric Flows of Ideal and Viscous Fluids Filling the Whole Space”

- Journal of Applied Mathematics and Mechanics 32(1) (1968), 52–62.
- DOI: <https://doi.org/10.1016/0021-8928(68)90147-0>
- 用途: 全空間の旋回なし軸対称流の大域可解性と \(\omega^\theta/r\) の構造。
- 状態: **書誌・要旨確認済み**

### [KNSS2009] Gabriel Koch, Nikolai Nadirashvili, Gregory Seregin, Vladimír Šverák, “Liouville Theorems for the Navier–Stokes Equations and Applications”

- Acta Mathematica 203 (2009), 83–105.
- DOI: <https://doi.org/10.1007/s11511-009-0039-6>
- 著者公開版: <https://www-users.cse.umn.edu/~sverak/publications/liouville.pdf>
- 用途: bounded ancient mild/weak solutions、旋回なし軸対称 ancient solution の Liouville 定理、軸対称流の尺度不変点wise上界からの正則性、Type I シナリオ。
- 状態: **本文確認済み**

## 数値探索と Type II（証明ではない）

### [Hou2023] Thomas Y. Hou, “Potentially Singular Behavior of the 3D Navier–Stokes Equations”

- Foundations of Computational Mathematics 23 (2023), 2251–2299.
- DOI: <https://doi.org/10.1007/s10208-022-09578-4>
- 著者公開版: <https://users.cms.caltech.edu/~hou/papers/FoCM-Navier-Stokes-2022.pdf>
- プレプリント: <https://arxiv.org/abs/2107.06509>
- 用途: 軸対称・旋回ありの高解像度候補探索、適応格子、複数の継続判定量、局所 \(L^3\) 診断の先行例。
- 状態: **本文・要旨確認済み**
- **注意**: 論文自身が “potentially singular” と述べる数値的証拠であり、特異点の証明ではない。

### [Seregin2024] Gregory Seregin, “A Note on Potential Type II Blowups of Axisymmetric Solutions to the Navier–Stokes Equations”

- プレプリント: <https://arxiv.org/abs/2402.13229>
- 用途: Type II の Euler スケーリング極限、軸対称性、可積分性・自己相似性に対する必要条件。
- 状態: **本文確認済み**
- 注意: プレプリントとして扱い、査読済み定理と混同しない。

### [Seregin2026] Gregory Seregin, “On Potential Type II Blowups for the Navier–Stokes Equations”

- プレプリント: <https://arxiv.org/abs/2606.29468>
- 用途: 局所 Type II シナリオ、Euler スケーリングと ancient limit に対する最新の必要条件を探索フィルターとして検討する。
- 状態: **本文確認済み**
- 注意: 2026年6月公開のプレプリント。探索空間を排除する最終根拠にはせず、仮定を原文どおり個別確認する。

## AI候補探索の失敗モード

### [Rahaman2019] Nasim Rahaman et al., “On the Spectral Bias of Neural Networks”

- Proceedings of the 36th ICML, PMLR 97 (2019), 5301–5310.
- 公開版: <https://proceedings.mlr.press/v97/rahaman19a.html>
- 用途: ニューラルネットが低周波成分を先に学習しやすいスペクトルバイアス。
- 状態: **本文・要旨確認済み**

### [Krishnapriyan2021] Aditi S. Krishnapriyan et al., “Characterizing Possible Failure Modes in Physics-Informed Neural Networks”

- Advances in Neural Information Processing Systems 34 (2021), 26548–26560.
- 公開版: <https://proceedings.neurips.cc/paper_files/paper/2021/hash/df438e5206f31600e6ae4af72f2725f1-Abstract.html>
- 用途: PDE損失が悪条件化した最適化地形を作り、表現能力とは別に偽収束を起こしうること。
- 状態: **本文・要旨確認済み**

## 引用上の規則

1. 「数値的証拠」「候補」「必要条件」「十分条件」「証明」を区別する。
2. 定理を適用するときは、領域、解のクラス、時間区間、積分指数、端点、空間無限遠条件を省略しない。
3. [BKM1984] を Navier–Stokes の原論文と呼ばず、Navier–Stokes の
   継続判定には [KatoPonce1988] を直接引用する。
4. [ESS2003] の \(L_{3,\infty}\) を Lorentz \(L^{3,\infty}\) と読み違えない。
5. [KNSS2009] の「形式的な5次元ラプラシアン」はスカラー作用素の表現であり、3次元速度場の物理的発散条件を5次元発散条件に置換しない。
6. プレプリント [Seregin2024], [Seregin2026] は査読状況を明記する。
7. 参照先が訂正された場合は、版・正誤表・取得日をこの台帳に追記する。
