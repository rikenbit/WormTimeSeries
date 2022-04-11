#library
##################################################
library(tidyverse)
##################################################
random.sampling <- function(A, n){
  # 各行(細胞)x各列のランダムサンプリング
  out <- unlist(
    # 各行(各細胞)だけ、割り当てる細胞数分実行
    lapply(seq(n), 
           function(x){
             # MARGIN = 2,各列（各TimeFrame）に対して関数適用。
             apply(A, 2, 
                   # replace=FALSEなので重複なしでランダムサンプリング。
                   # 緑（数字名の細胞達）の「ある1TimeFrame」の「ある細胞」のデータを取り出す
                   function(xx){sample(xx, 1)}
             )
           }
    )
  )
  # 数値ベクトルの形になっているので、行列の形式に当てはめる
  dim(out) <- c(n, length(out)/n)
  out
}