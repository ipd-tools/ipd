
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/ipd-tools/ipd/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ipd-tools/ipd/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

## <img src="man/figures/ipd.png" align="right" height="200" style="float:right; height:200px;"/>

### 概要

# ipd: Inference on Predicted Data

ipdは、人工知能（AI）や機械学習（ML）によって予測されたデータを用いた統計的推論（Inference
on Predicted Data: IPD）を支援する、オープンソースのRパッケージです。

現代の多くの応用分野では、ラベル付きデータ（結果が観測されているデータ）が限られており、代わりにAIモデルが予測したラベルを使用して解析を行う場面が増えています。ipdパッケージは、こうした「予測された結果（outcomes）」を使って、共変量（説明変数）と結果の関係を推定する際に生じるバイアスや不確実性を適切に考慮した推論を可能にします。

本パッケージでは、最近提案された複数の手法（PostPI、PPI、PSPAなど）を、共通のインターフェースを通じて簡単に使えるように統合しています。また、print、summary、tidy、glance、augmentといったモデル出力の整形や確認に便利なメソッドも提供しています。

### 背景

AIや機械学習モデルによる予測結果は、近年さまざまな分野で活用が進んでいます。特に、医療や社会科学、自然科学などでは、ラベル付きデータの収集が高コストまたは困難であるため、AIが予測した結果をそのまま下流の統計解析に利用するケースが増えています。

しかし、このような予測値を本来の「観測値」の代わりに使うことで、以下のような統計的な問題が生じる可能性があります：

予測値と真の値のズレ：予測された結果（たとえば疾患の有無や治療効果）は、実際に観測されたデータとは異なる性質を持ちます。このズレを考慮しないと、推定に偏りが生じます。
不確実性の過小評価：予測には必ず誤差が含まれますが、それを無視して推論を行うと、標準誤差が小さく見積もられ、信頼区間や仮説検定が正しく機能しません。
AIモデル自体のバイアス：予測に使われるAI/MLモデルが訓練されたデータや手法によっては、特定のサンプルに偏った結果を出すことがあり、そのバイアスが下流の解析結果に影響します。
こうした問題に対処するため、近年さまざまな方法論が提案されてきました。たとえば：

PostPI（予測後推論；Wang et al.,
2020）：予測された結果と実際の結果の差を補正する方法 PPI /
PPI++（予測駆動推論；Angelopoulos et al.,
2023）：予測モデルを再利用して不確実性を伝搬させる手法
PSPA（予測後適応推論；Miao et al.,
2023）：データの構造と予測精度を利用して適応的に推論を行う方法
これらの手法は共通して、観測されていない結果をAIで補完したうえで、信頼できる統計的結論を導くことを目的としています。

<figure>
<img src="man/figures/ipd_overview.png"
alt="Overview of data and setup for IPD" />
<figcaption aria-hidden="true">Overview of data and setup for
IPD</figcaption>
</figure>

ラベルなしデータの特徴量と予測結果を活用して、ラベル付きデータを補完することで、IPD手法を用いて補正された推定値と標準誤差を得ることができます。

こうした最先端の手法に関心のある研究者や実務者が簡単に利用できるよう、私たちはそれらをIPDという枠組みのもとで実装したRパッケージipdを開発しました。本READMEでは、本パッケージの概要、インストール方法、基本的な使用例、および詳細なドキュメントへのリンクを紹介しています。使用例では、データの生成、モデルの適合、そしてパッケージが提供するカスタムメソッドの使い方を説明しています。

\##インストール方法

GitHubからipdの開発版をインストールするには、devtoolsパッケージを使用してください。

``` r
    #-- devtoolsがインストールされていない場合はインストール

    install.packages("devtools")   

    #-- GitHubからipdパッケージをインストール

    devtools::install_github("ipd-tools/ipd")
```

\##使用方法

ipdパッケージに含まれるメソッドの使い方を示すシンプルな例を紹介します。

\###データの生成
simdat関数を使えば、各種回帰モデルに対応した合成データを生成できます。データセットのサイズ、効果量、残差分散、そしてモデルの種類を指定することで、自由にデータを作成できます。現在サポートされているモデルの種類は、「mean（平均）」「quantile（分位数）」「ols（線形回帰）」「logistic（ロジスティック回帰）」「poisson（ポアソン回帰）」です。

simdat関数は、次の3種類のサブセットを含むdata.frameを生成します：

予測モデルの学習に使用される追加の観測値を含む独立した「training（訓練）」セット
観測された結果と予測された結果を含む「labeled（ラベル付き）」セット
同じく観測された結果と予測値、および特徴量を含む「unlabeled（ラベルなし）」セット

``` r
#– ipdライブラリを読み込む

library(ipd)

#– 線形回帰用のサンプルデータを生成

set.seed(123)

n <- c(10000, 500, 1000)

dat <- simdat(n = n, effect = 1, sigma_Y = 4, model = 'ols')

#– 訓練、ラベル付き、ラベルなしサブセットの最初の6行を表示

options(digits=2)

head(dat[dat$set_label == 'training',])
#>       X1    X2    X3     X4     Y  f set_label
#> 1 -0.560 -0.56  0.82 -0.356 -0.15 NA  training
#> 2 -0.230  0.13 -1.54  0.040 -4.49 NA  training
#> 3  1.559  1.82 -0.59  1.152 -1.08 NA  training
#> 4  0.071  0.16 -0.18  1.485 -3.67 NA  training
#> 5  0.129 -0.72 -0.71  0.634  2.19 NA  training
#> 6  1.715  0.58 -0.54 -0.037 -1.42 NA  training

head(dat[dat$set_label == 'labeled',])
#>          X1      X2    X3    X4     Y     f set_label
#> 10001  2.37 -1.8984  0.20 -0.17  1.40  3.24   labeled
#> 10002 -0.17  1.7428  0.26 -2.05  3.56  1.03   labeled
#> 10003  0.93 -1.0947  0.76  1.25 -3.66  2.37   labeled
#> 10004 -0.57  0.1757  0.32  0.65 -0.56  0.58   labeled
#> 10005  0.23  2.0620 -1.35  1.46 -0.82 -0.15   labeled
#> 10006  1.13 -0.0028  0.23 -0.24  7.30  2.16   labeled

head(dat[dat$set_label == 'unlabeled',])
#>          X1     X2    X3    X4    Y     f set_label
#> 10501  0.99 -3.280 -0.39  0.97  8.4  1.25 unlabeled
#> 10502 -0.66  0.142 -1.36 -0.22 -7.2 -1.08 unlabeled
#> 10503  0.58 -1.368 -1.73  0.15  5.6 -0.31 unlabeled
#> 10504 -0.14 -0.728  0.26 -0.23 -4.2  0.91 unlabeled
#> 10505 -0.17 -0.068 -1.10  0.58  2.2 -0.39 unlabeled
#> 10506  0.58  0.514 -0.69  0.97 -1.2  0.76 unlabeled
```

`simdat`関数は、ラベル付きとラベルなしの両方のデータセットに対して、観測された結果と観測されていない結果を提供しますが、実際にはラベルなしセットには観測された結果はありません。これらの変数間の関係を視覚化できます：

<img src="man/figures/README-plot-1.png" width="100%" />

グラフからは次のような傾向が読み取れます：

-プロットA：予測値は実測値よりも共変量との相関が高くなる傾向があります。

-プロットB：予測値は実測値と異なる分布を持つ。

### モデルの適合

データを分析するための2つの非IPDアプローチと`ipd`パッケージに含まれる手法を比較します。以下の表に比較概要表を示し、その後に各手法の具体的な呼び出し方法を示します：

    #>                            Estimate Std.Error
    #> ナイーブ                       0.98      0.03
    #> クラシック                     1.10      0.19
    #> PostPI\n(ブートストラップ)     1.16      0.18
    #> PostPI (解析的)                1.15      0.18
    #> PPI++                          1.12      0.19
    #> PSPA                           1.12      0.19

IPD手法では、推定値や標準誤差は概ね一致しますが、ナイーブ手法では推定値が異なり、標準誤差も過小評価される傾向があります。

#### 0.1 予測された結果を使用した「ナイーブ」回帰

``` r
    #--- ナイーブ回帰を適合

    lm(f ~ X1, data = dat[dat$set_label == "unlabeled",]) |> 
      
      summary()
#> 
#> Call:
#> lm(formula = f ~ X1, data = dat[dat$set_label == "unlabeled", 
#>     ])
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.5426 -0.6138 -0.0153  0.6345  2.8907 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)   0.8391     0.0297    28.3   <2e-16 ***
#> X1            0.9848     0.0296    33.3   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.94 on 998 degrees of freedom
#> Multiple R-squared:  0.527,  Adjusted R-squared:  0.526 
#> F-statistic: 1.11e+03 on 1 and 998 DF,  p-value: <2e-16
```

#### 0.2 ラベル付きデータのみを使用した「クラシック」回帰

``` r
#— クラシック回帰を適合

lm(Y ~ X1, data = dat[dat$set_label == 'labeled',]) |>

summary()
#> 
#> Call:
#> lm(formula = Y ~ X1, data = dat[dat$set_label == "labeled", ])
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -15.262  -2.828  -0.094   2.821  11.685 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)    0.908      0.187    4.86  1.6e-06 ***
#> X1             1.097      0.192    5.71  1.9e-08 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 4.2 on 498 degrees of freedom
#> Multiple R-squared:  0.0614, Adjusted R-squared:  0.0596 
#> F-statistic: 32.6 on 1 and 498 DF,  p-value: 1.95e-08
```

提供されるラッパー関数`ipd()`を使用して、様々なIPD手法をデータに適合させ、要約を取得できます：

#### 1.1 PostPIブートストラップ補正（Wang et al., 2020）

``` r
    #-- 式を指定

    formula <- Y - f ~ X1

    #-- PostPIブートストラップ補正を適合

    nboot <- 200

    ipd::ipd(formula, 
             
      method = "postpi_boot", model = "ols", data = dat, label = "set_label", 
      
      nboot = nboot) |> 
      
      summary()
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Method: postpi_boot 
#> Model: ols 
#> Intercept: Yes 
#> 
#> Coefficients:
#>             Estimate Std.Error Lower.CI Upper.CI
#> (Intercept)    0.866     0.183    0.507     1.22
#> X1             1.164     0.183    0.806     1.52
```

#### 1.2 PostPI解析的補正（Wang et al., 2020）

``` r
#– PostPI解析的補正を適合

ipd::ipd(formula,

method = 'postpi_analytic', model = 'ols', data = dat, label =
'set_label') |>

summary()
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Method: postpi_analytic 
#> Model: ols 
#> Intercept: Yes 
#> 
#> Coefficients:
#>             Estimate Std.Error Lower.CI Upper.CI
#> (Intercept)    0.865     0.183    0.505     1.22
#> X1             1.145     0.182    0.788     1.50
```

#### 2. 予測駆動推論（PPI; Angelopoulos et al., 2023）

``` r
    #-- PPI補正を適合

    ipd::ipd(formula, 
             
      method = "ppi", model = "ols", data = dat, label = "set_label") |> 
      
      summary()
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Method: ppi 
#> Model: ols 
#> Intercept: Yes 
#> 
#> Coefficients:
#>             Estimate Std.Error Lower.CI Upper.CI
#> (Intercept)    0.871     0.182    0.514     1.23
#> X1             1.122     0.195    0.740     1.50
```

#### 3. PPI++（Angelopoulos et al., 2023）

``` r
#– PPI++補正を適合

ipd::ipd(formula,

method = 'ppi_plusplus', model = 'ols', data = dat, label = 'set_label')|>

summary()
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Method: ppi_plusplus 
#> Model: ols 
#> Intercept: Yes 
#> 
#> Coefficients:
#>             Estimate Std.Error Lower.CI Upper.CI
#> (Intercept)    0.881     0.182    0.524     1.24
#> X1             1.116     0.187    0.750     1.48
```

#### 4. 予測後適応推論（PSPA; Miao et al., 2023）

``` r
    #-- PSPA補正を適合

    ipd::ipd(formula, 
             
      method = "pspa", model = "ols", data = dat, label = "set_label") |> 
      
      summary()
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Method: pspa 
#> Model: ols 
#> Intercept: Yes 
#> 
#> Coefficients:
#>             Estimate Std.Error Lower.CI Upper.CI
#> (Intercept)    0.881     0.182    0.524     1.24
#> X1             1.109     0.187    0.743     1.47
```

### 出力と整形

モデルの検証を容易にするため、print、summary、tidy、glance、augmentといったカスタムメソッドも利用できます：

``` r
#– PostPIブートストラップ補正を適合

nboot <- 200

fit_postpi <- ipd::ipd(formula,

method = 'postpi_boot', model = 'ols', data = dat, label = 'set_label',

nboot = nboot)

#– モデルを出力

print(fit_postpi)
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Coefficients:
#> (Intercept)          X1 
#>        0.86        1.15

#– モデルを要約

summ_fit_postpi <- summary(fit_postpi)

#– モデル要約を出力

print(summ_fit_postpi)
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Method: postpi_boot 
#> Model: ols 
#> Intercept: Yes 
#> 
#> Coefficients:
#>             Estimate Std.Error Lower.CI Upper.CI
#> (Intercept)    0.860     0.183    0.502     1.22
#> X1             1.148     0.182    0.790     1.50

#– モデル出力を整形

tidy(fit_postpi)
#>                    term estimate std.error conf.low conf.high
#> (Intercept) (Intercept)     0.86      0.18     0.50       1.2
#> X1                   X1     1.15      0.18     0.79       1.5

#– モデルの1行要約を取得

glance(fit_postpi)
#>        method model include_intercept nobs_labeled nobs_unlabeled       call
#> 1 postpi_boot   ols              TRUE          500           1000 Y - f ~ X1

#– 元のデータに適合値と残差を追加

augmented_df <- augment(fit_postpi)

head(augmented_df)
#>          X1     X2    X3    X4    Y     f set_label .fitted .resid
#> 10501  0.99 -3.280 -0.39  0.97  8.4  1.25 unlabeled   1.992    6.5
#> 10502 -0.66  0.142 -1.36 -0.22 -7.2 -1.08 unlabeled   0.099   -7.3
#> 10503  0.58 -1.368 -1.73  0.15  5.6 -0.31 unlabeled   1.522    4.1
#> 10504 -0.14 -0.728  0.26 -0.23 -4.2  0.91 unlabeled   0.702   -4.9
#> 10505 -0.17 -0.068 -1.10  0.58  2.2 -0.39 unlabeled   0.667    1.5
#> 10506  0.58  0.514 -0.69  0.97 -1.2  0.76 unlabeled   1.521   -2.7
```

## ビネット

さらなる詳細として、パッケージのビネットでより多くのユースケースと例を提供しています：

``` r
    vignette("ipd")
```

## フィードバック

質問、コメント、その他のフィードバックについては、開発者（<ssalerno@fredhutch.org>）にお問い合わせください。

## 貢献

貢献を歓迎します！[GitHub](https://github.com/ipd-tools/ipd)で問題を開くか、プルリクエストを提出してください。現在、以下の手法/モデルの組み合わせが実装されています：

| 手法                                                            | 平均推定           | 分位数推定         | 線形回帰           | ロジスティック回帰 | ポアソン回帰       | 多クラス回帰 |
|-----------------------------------------------------------------|--------------------|--------------------|--------------------|--------------------|--------------------|--------------|
| [PostPI](https://www.pnas.org/doi/full/10.1073/pnas.2001238117) | :x:                | :x:                | :white_check_mark: | :white_check_mark: | :x:                | :x:          |
| [PPI](https://www.science.org/doi/10.1126/science.adi6000)      | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | :x:                | :x:          |
| [PPI++](https://arxiv.org/abs/2311.01453)                       | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | :x:                | :x:          |
| [PSPA](https://arxiv.org/abs/2311.14220)                        | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: | :x:          |
| [PSPS](https://arxiv.org/abs/2405.20039)                        | :x:                | :x:                | :x:                | :x:                | :x:                | :x:          |
| [PDC](https://arxiv.org/abs/2312.06478)                         | :x:                | :x:                | :x:                | :x:                | :x:                | :x:          |
| [Cross-PPI](https://www.pnas.org/doi/10.1073/pnas.2322083121)   | :x:                | :x:                | :x:                | :x:                | :x:                | :x:          |
| [PPBoot](https://arxiv.org/abs/2405.18379)                      | :x:                | :x:                | :x:                | :x:                | :x:                | :x:          |
| [DSL](https://naokiegami.com/paper/dsl.pdf)                     | :x:                | :x:                | :x:                | :x:                | :x:                | :x:          |

## ライセンス

このパッケージはMITライセンスの下でライセンスされています。
