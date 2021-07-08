source("src/functions_WTS3_DimReduc.R")

#### args setting####
#### test args####
# select animal number 個体番号の指定
args_sample <- c("1")
# input 距離データ
args_dist <- c("output/WTS3/normalize_1/all/SBD/SampleNumber_1/SBD.RData")
# input SBD yshift df
args_yshift_value <- c("output/WTS3/normalize_1/all/SBD/SampleNumber_1/yshift_value.RData")
# input label_table
args_label <- c("output/WTS3/normalize_1/all/SBD/ARI/SampleNumber_1/label_table.RData")
# input cutree_table
args_cutree <- c("output/WTS3/normalize_1/all/SBD/ARI/SampleNumber_1/cutree_table.RData")
# 次元圧縮手法
args_DimReduc <- c("tsne")
# クラスタリング評価
args_eval <- c("ARI")
# output plot
args_output <- c("output/WTS3/normalize_1/all/SBD/ARI/DimReduc/SampleNumber_1.png")

### load SBD####
load(args_dist)
#### Dimensionality Reduction####
df_cord <- switch(args_DimReduc,
          "tsne" = .wts_tsne(d),
          "umap" = .wts_umap(d),
          stop("Only can use tsne,")
          )
#### load yshift####
load(args_yshift_value)
#### load label####
load(args_label)
#### prepare ggplot table####
df_label_cord <- merge(df_label, 
                       df_cord, 
                       by.x = "cell_type", 
                       by.y = "cell_type", 
                       all.x = TRUE)

df_label_cord_yshift <- merge(df_label_cord, 
                              yshift_value_table, 
                              by.x = "cell_type", 
                              by.y = "cell_type", 
                              all.x = TRUE)

#### ggplot NeuronType x label_acf####
df_label_cord_yshift$NeuronType <- as.character(df_label_cord_yshift$NeuronType)
df_label_cord_yshift$label_acf <- as.character(df_label_cord_yshift$label_acf)
colar_col <- "NeuronType"
shape_col <- "label_acf"
gg_label_acf <- .gg_label(colar_col, shape_col)
#### ggplot NeuronType x label_acf x yshift####
colar_col <- "NeuronType"
shape_col <- "label_acf"
gg_label_acf_yshift <- .gg_label_yshift(colar_col, shape_col)
  
#### ggplot clustering number x label_cls####
df_label_cord_yshift$label_cls <- as.character(df_label_cord_yshift$label_cls)
df_label_cord_yshift$cls <- as.character(df_label_cord_yshift$cls)
colar_col <- "cls"
shape_col <- "label_cls"
gg_label_cls <- .gg_label(colar_col, shape_col)
#### ggplot clustering number x label_cls x yshift####
colar_col <- "cls"
shape_col <- "label_cls"
gg_label_cls_yshift <- .gg_label_yshift(colar_col, shape_col)
#### ggpubr table of eval_value####
load(args_cutree)
df_cutree %>%
    dplyr::select(set_cutree, eval_value) %>% 
        dplyr::summarise_all(list(round), digits=5) %>% 
            ggtexttable(rows = NULL, theme = ttheme(base_size = 20)) -> gg_cutree # ggpubr
gg_cutree <- table_cell_bg(gg_cutree, 
                           row = which(df_cutree$select_cutree==1) + 1,
                           column = 1:2, 
                           linewidth = 5,
                           fill="darkolivegreen1", 
                           color = "darkolivegreen4")
#### patchwork & save####
eval(parse(text=paste0("plot_title <- c('SampleNumber_",args_sample,"_evalMethod_",args_eval,"')")))
gg <- gg_label_acf + gg_label_acf_yshift + plot_spacer() + gg_label_cls + gg_label_cls_yshift + gg_cutree
gg <- gg + 
    plot_annotation(
        title = plot_title,
        caption = 'made with patchwork',
        theme = theme(plot.title = element_text(size = 60, hjust = 0.5))
    )
gg <- gg + plot_layout(ncol = 3,widths = c(2, 2, 1))
ggsave(filename = args_output, 
       plot = gg, 
       dpi = 100, 
       width = 25.0, 
       height = 20.0,
       limitsize = FALSE
       )