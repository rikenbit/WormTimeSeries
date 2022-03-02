#library
##################################################
library(openxlsx)
library(tidyverse)
library(factoextra)
library(ggpubr)
library(patchwork)
library(cluster)
##################################################

#### myfviz_silhouette####
# https://stackoverflow.com/questions/61820019/how-can-i-change-the-color-to-a-variable-other-than-cluster-number-in-fviz-silho
myfviz_silhouette <- function (sil.obj, var.col, label = FALSE, print.summary = TRUE, ...) {
  if (inherits(sil.obj, c("eclust", "hcut", "pam", "clara", 
                          "fanny"))) {
    df <- as.data.frame(sil.obj$silinfo$widths, stringsAsFactors = TRUE)
  }
  else if (inherits(sil.obj, "silhouette")) 
    df <- as.data.frame(sil.obj[, 1:3], stringsAsFactors = TRUE)
  else stop("Don't support an oject of class ", class(sil.obj))
  df <- df[order(df$cluster, -df$sil_width), ]
  if (!is.null(rownames(df))) 
    df$name <- factor(rownames(df), levels = rownames(df))
  else df$name <- as.factor(1:nrow(df))
  df$cluster <- as.factor(df$cluster)
  # fix var_col -> cluster_col
  df$cluster_col <- var.col
  mapping <- aes_string(x = "name", y = "sil_width", color = "cluster_col", 
                        fill = "cluster_col")
  p <- ggplot(df, mapping) + geom_bar(stat = "identity") + 
    labs(y = "Silhouette width", x = "", title = "", color ="Cluster",fill = "Cluster") + ggplot2::ylim(c(NA, 1)) + geom_hline(yintercept = mean(df$sil_width), linetype = "dashed", color = "red")
  p <- ggpubr::ggpar(p, ...)
  if (!label) 
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  else if (label) 
    p <- p + theme(axis.text.x = element_text(angle = 45))
  ave <- tapply(df$sil_width, df$cluster, mean)
  n <- tapply(df$cluster, df$cluster, length)
  sil.sum <- data.frame(cluster = names(ave), size = n, ave.sil.width = round(ave, 
                                                                              2), stringsAsFactors = TRUE)
  if (print.summary) 
    print(sil.sum)
  p
}
# myfviz_silhouette <- function (sil.obj, var.col, label = FALSE, print.summary = TRUE, ...) {
#   if (inherits(sil.obj, c("eclust", "hcut", "pam", "clara", 
#                           "fanny"))) {
#     df <- as.data.frame(sil.obj$silinfo$widths, stringsAsFactors = TRUE)
#   }
#   else if (inherits(sil.obj, "silhouette")) 
#     df <- as.data.frame(sil.obj[, 1:3], stringsAsFactors = TRUE)
#   else stop("Don't support an oject of class ", class(sil.obj))
#   df <- df[order(df$cluster, -df$sil_width), ]
#   if (!is.null(rownames(df))) 
#     df$name <- factor(rownames(df), levels = rownames(df))
#   else df$name <- as.factor(1:nrow(df))
#   df$cluster <- as.factor(df$cluster)
#   # fix var_col -> cluster_col
#   df$cluster_col <- var.col
#   mapping <- aes_string(x = "name", y = "sil_width", color = "cluster_col", 
#                         fill = "cluster_col")
#   p <- ggplot(df, mapping) + geom_bar(stat = "identity") + 
#     labs(y = "Silhouette width Si", x = "", title = paste0("Clusters silhouette plot ", 
#                                                            "\n Average silhouette width: ", round(mean(df$sil_width), 
#                                                                                                   2))) + ggplot2::ylim(c(NA, 1)) + geom_hline(yintercept = mean(df$sil_width), 
#                                                                                                                                               linetype = "dashed", color = "red")
#   p <- ggpubr::ggpar(p, ...)
#   if (!label) 
#     p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
#   else if (label) 
#     p <- p + theme(axis.text.x = element_text(angle = 45))
#   ave <- tapply(df$sil_width, df$cluster, mean)
#   n <- tapply(df$cluster, df$cluster, length)
#   sil.sum <- data.frame(cluster = names(ave), size = n, ave.sil.width = round(ave, 
#                                                                               2), stringsAsFactors = TRUE)
#   if (print.summary) 
#     print(sil.sum)
#   p
# }