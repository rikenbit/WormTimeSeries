source("src/functions_WTS4_Membership_vis.R")

#### args setting####
#### test args####
args_input <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_6/SampleNumber_1.RData")
args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_vis/k_Number_6/SampleNumber_1.png")
#### load yshift_mem_df####
load(args_input)
#### filter col cell_cell,mem,shift####
yshift_mem_df %>% 
    dplyr::select(cell_cell, member, yshift) -> ggdata
#### ggplot####
gg_box <- ggplot(ggdata, aes(x=factor(member), y=yshift)) +
    guides(colour="none") +
    # geom_boxplot(alpha =0.8) +
    #外れ値のプロットを省く https://stats.biopapyrus.jp/r/ggplot/geom-boxplot.html
    geom_boxplot(outlier.shape = NA, alpha =0.8) +
    # geom_jitterを使うとかぶらない http://sakananoiroiro.seesaa.net/article/454739891.html
    geom_point(position = position_jitter(width=0.05), size = 3.0 ,alpha = 0.7) +
    # geom_jitter(size = 0.8,alpha = 0.7) +
    xlab("Membership") + 
    ylab("yshift") + 
    scale_y_continuous(limits = c(-6000, 6000)) +
    theme(legend.position = "none") +
    theme(text = element_text(size = 60)) 

gg_vio <- ggplot(ggdata, aes(x=factor(member), y=yshift)) +
    guides(colour="none") +
    # geom_boxplot(alpha =0.8) +
    #外れ値のプロットを省く https://stats.biopapyrus.jp/r/ggplot/geom-boxplot.html
    # geom_boxplot(outlier.shape = NA, alpha =0.8) +
    geom_violin() +
    # geom_jitterを使うとかぶらない http://sakananoiroiro.seesaa.net/article/454739891.html
    # geom_point(position = position_jitter(width=0.05), size = 3.0 ,alpha = 0.7) +
    # geom_jitter(size = 0.8,alpha = 0.7) +
    xlab("Membership") +
    ylab("yshift") +
    scale_y_continuous(limits = c(-6000, 6000)) +
    theme(legend.position = "none") +
    theme(text = element_text(size = 60))
#### t-Test####
# ggdata
T_Result <- t.test(yshift ~ member, data = ggdata)
if(T_Result$p.value==0){
    T_pvalue <- c("p-value < 2.2e-16")
}else{
    T_pvalue <- signif(T_Result$p.value, digits = 3)
}

#### u-Test####
U_Result <- wilcox.test(yshift ~ member, data = ggdata)
if(U_Result$p.value==0){
    U_pvalue <- c("p-value < 2.2e-16")
}else{
    U_pvalue <- signif(U_Result$p.value, digits = 3)
}

#### F-Test####
F_Result <- var.test(yshift ~ member, data = ggdata)
if(F_Result$p.value==0){
    F_pvalue <- c("p-value < 2.2e-16")
}else{
    F_pvalue <- signif(F_Result$p.value, digits = 3)
}
#### Test table####
data.frame(T_pvalue = T_pvalue,
           U_pvalue = U_pvalue,
           F_pvalue = F_pvalue,
           stringsAsFactors = FALSE,
           row.names = NULL) -> t_table
t_table %>% 
    ggtexttable(rows = NULL, 
                theme = ttheme(base_size = 60, 
                               base_style ="mBlue")
                ) -> gg_table
#### graph title####
args_output |> 
    str_remove("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_vis/") |> 
    str_remove(".png") -> plot_title

#### patchwork####
# gg <- gg_box +
#     gg_table +
#     plot_layout(ncol = 1, heights = c(3, 1)) +
#     plot_annotation(title = plot_title,
#                     caption = 'made with patchwork',
#                     theme = theme(plot.title = element_text(size = 60, hjust = 0.5))
#     )
gg <- gg_box +
    gg_vio -
    gg_table +
    plot_layout(ncol = 1, heights = c(3, 1)) +
    plot_annotation(title = plot_title,
                    caption = 'made with patchwork',
                    theme = theme(plot.title = element_text(size = 60, hjust = 0.5))
    )
#### ggsave####
ggsave(filename = args_output, 
       plot = gg,
       dpi = 50, 
       width = 20.0, 
       height = 20.0,
       limitsize = FALSE
)