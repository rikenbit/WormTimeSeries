#library
##################################################
library(tidyverse)
library(openxlsx)
##################################################
# 同じラベルのメンバーと同じクラスタに入った回数
# 一回、接続行列に変換
.H <- function(vec){
  uniq_vec <- unique(vec)
  out <- matrix(0, nrow=length(vec), ncol=length(uniq_vec))
  rownames(out) <- names(vec)
  colnames(out) <- uniq_vec
  for(i in seq_len(length(vec))){
    idx <- which(vec[i] == uniq_vec)
    out[i, idx] <- 1
  }
  out
}

# 自分以外のメンバー数 正規化用
.other_member <- function(H_label){
  out <- rep(0, length=nrow(H_label))
  for(i in seq_len(nrow(H_label))){
    idx1 <- which(H_label[i, ] == 1)
    out[i] <- sum(H_label[,idx1]) - 1
  }
  out
}

.df_count = function(x) {
  #### 各個体の行動ラベル回数データフレーム作成####
  cls_n <- C[[x]]
  
  df_cls <- as.data.frame(cls_n)
  df_cls <- data.frame(cell_type = rownames(df_cls),
                       Cluster = df_cls$cls_n,
                       stringsAsFactors = FALSE
  )
  df_cls_label <- merge(df_eval_label, 
                        df_cls, 
                        by.x = "CellType", 
                        by.y = "cell_type", 
                        all.y = TRUE)
  # NA to others
  df_cls_label$Classes <- ifelse(is.na(df_cls_label$Classes), "others", df_cls_label$Classes)
  
  label <- df_cls_label$Classes
  cluster <- df_cls_label$Cluster
  H_label <- .H(label)
  H_cluster <- .H(cluster)
  
  other_member <- .other_member(H_label)
  
  # カウント
  count1 <- rep(0, length=length(label))
  for(i in seq_len(nrow(H_label))){
    idx1 <- which(H_label[i, ] == 1)
    idx2 <- which(H_cluster[i, ] == 1)
    count1[i] <- sum(H_label[,idx1] * H_cluster[,idx2]) - 1
  }
  # 正規化
  count1 <- count1 / other_member
  # 0/0のNaN対策
  count1 <- ifelse(is.nan(count1), 0, count1)
  # 1個体分 count1 
  df_cls_label %>% 
    mutate(Count=count1) -> df_count
  df_count <- df_count[grep("^[0-9]",df_count$CellType, invert=TRUE),]
	return(df_count)
}

.df_count_rbind = function(x) {
  df_count_list[[x]] %>% 
    mutate(SampleNumber=as.character(c(rep(sample_sort_num[x])))) -> return_object
  return(return_object)
}