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
gg <- ggplot(ggdata, aes(x=factor(member), y=yshift)) +
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
    
#### t-Test####
#### u-Test####
#### patchwork####
#### patchworkでタイトル追加####
#### ggsave####
ggsave(filename = args_output, 
       plot = gg,
       dpi = 50, 
       width = 20.0, 
       height = 20.0,
       limitsize = FALSE
)