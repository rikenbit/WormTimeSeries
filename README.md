# WormTimeSeries
- Since the neural activity patterns showed little correlation when viewed at the same time, we investigated whether there was any correlation with neural activity shifted by time τ.
- There are 7 analysis.  
1. Time series neural activity plot (each cell)  
2. Correlogram (each cell, each c.elegans)  
3. Time series dimensionality reduction and clustering (each c.elegans)  
4. Spectrogram (each cell)  
5. Time-delay embedding (each cell)  
6. Eigenvalue distribution of the history matrix of all cells (each c.elegans)  
7. Recurrence plot (each cell)


Installation
======
```bash
git clone https://github.com/rikenbit/WormTimeSeries/
cd WormTimeSeries
snakemake -j 1
```
# Prerequisites
- Snakemake (5.32.0 or higher)

# Summary
## 1. Time series neural activity plot  
- DAG  
![DAG](output/WTS1/dag.png)
- Plot(exsample)  
- Pattern  
- Artifact

## 2. Correlogram  
## 3. Time series dimensionality reduction and clustering  
## 4. Spectrogram  
## 5. Time-delay embedding  
## 6. Eigenvalue distribution of the history matrix of all cells  
## 7. Recurrence plot  

# How to reproduce this workflow
## 1. Time series neural activity plot  
- Files to prepare in advance (data/)
	- 15個体の神経活動データ
		- cleandata/xxx_ratio.csv (xxx：1,…,15)
	- 15個体の塩刺激データ
		- stimulation_timing.csv

- Function each files (src/)
	- WTS1_load.R
		神経活動データ，塩刺激データ，全個体全細胞のスプレッドシートを作成.Snakefileの「WTS1 dataload」のセクション使用
	- WTS1_1cell.R
		「1. Time series neural activity plot」のggplotを行う.Snakefileの「WTS1 1cell plot」のセクション使用を使用することで，全個体全細胞のggplotを実行

- Memo
	- snakemakeは「###################################################」で囲われているエリアごとに実行．正しいやり方とは思えないので，そのうちRuleファイルを分けるなど要対応


- Conda recipe (envs/myenv_WTS1.yaml)
	- command
		- conda install -c r r-tidyverse -y
		- conda install -c r r-rcpproll -y

- sessionInfo(R env)  

R version 4.0.2 (2020-06-22)  
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)  

Matrix products: default  
BLAS/LAPACK: /home/yamaken/miniconda3/envs/WTS1/lib/libopenblasp-r0.3.10.so  

locale:
 [1] LC_CTYPE=ja_JP.UTF-8       LC_NUMERIC=C  
 [3] LC_TIME=ja_JP.UTF-8        LC_COLLATE=ja_JP.UTF-8  
 [5] LC_MONETARY=ja_JP.UTF-8    LC_MESSAGES=ja_JP.UTF-8   
 [7] LC_PAPER=ja_JP.UTF-8       LC_NAME=C  
 [9] LC_ADDRESS=C               LC_TELEPHONE=C  
[11] LC_MEASUREMENT=ja_JP.UTF-8 LC_IDENTIFICATION=C  

attached base packages:  
[1] stats     graphics  grDevices utils     datasets  methods   base  

other attached packages:
 [1] RcppRoll_0.3.0  forcats_0.5.0   stringr_1.4.0   dplyr_1.0.2  
 [5] purrr_0.3.4     readr_1.4.0     tidyr_1.1.2     tibble_3.0.4  
 [9] ggplot2_3.3.3   tidyverse_1.3.0

loaded via a namespace (and not attached):  
 [1] Rcpp_1.0.5        pillar_1.4.7      compiler_4.0.2    cellranger_1.1.0  
 [5] dbplyr_2.0.0      tools_4.0.2       digest_0.6.27     jsonlite_1.7.2  
 [9] lubridate_1.7.9.2 lifecycle_0.2.0   gtable_0.3.0      pkgconfig_2.0.3  
[13] rlang_0.4.10      reprex_0.3.0      cli_2.2.0         rstudioapi_0.13  
[17] DBI_1.1.0         haven_2.3.1       withr_2.3.0       xml2_1.3.2  
[21] httr_1.4.2        generics_0.1.0    vctrs_0.3.6       fs_1.5.0  
[25] hms_1.0.0         grid_4.0.2        tidyselect_1.1.0  glue_1.4.2  
[29] R6_2.5.0          fansi_0.4.1       readxl_1.3.1      farver_2.0.3  
[33] modelr_0.1.8      magrittr_2.0.1    ps_1.5.0          backports_1.2.1  
[37] scales_1.1.1      ellipsis_0.3.1    assertthat_0.2.1  rvest_0.3.6  
[41] colorspace_2.0-0  labeling_0.4.2    stringi_1.4.6     munsell_0.5.0  
[45] broom_0.7.3       crayon_1.3.4  

# Reference
1. [rikenbit/WormSCE](https://github.com/rikenbit/WormSCE)

# License
Copyright (c) 2020 Kentaro Yamamoto and RIKEN Bioinformatics Research Unit Released under the [Artistic License 2.0](http://www.perlfoundation.org/artistic_license_2_0).

# Authors
- Kentaro Yamamoto
- Koki Tsuyuzaki
- Itoshi Nikaido