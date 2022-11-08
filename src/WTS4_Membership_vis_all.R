source("src/functions_WTS4_Membership_vis_all.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input <- args[1]
args_output <- args[2]
#### test args####
# args_input <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_6/DFs.RData")
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_vis_all/k_Number_6.png")

#### load####
load(args_input)
#### yshift list to df####
dfr_yshift <- purrr::map_dfr(names(DFs), .DFs_yshift)

#### ggplot####
cord_x <- c("Membership")
cord_y <- c("yshift")

gg_box <- ggplot(dfr_yshift, 
                 aes(x=factor(member), 
                     y=yshift)
                 ) +
    geom_boxplot(outlier.shape = NA, alpha =0.8) +
    geom_point(position = position_jitter(width=0.1), size = 3.0 ,alpha = 0.7) +
    xlab("Membership") + 
    ylab("yshift") + 
    scale_y_continuous(limits = c(-6000, 6000)) +
    theme(text = element_text(size = 120)) 

gg_vio <- ggplot(dfr_yshift, 
                 aes(x=factor(member), 
                     y=yshift)
                 ) +
    geom_violin() +
    xlab("Membership") +
    ylab("yshift") +
    scale_y_continuous(limits = c(-6000, 6000)) +
    theme(legend.position = "none") +
    theme(text = element_text(size = 120))
#### t-Test####
# ggdata
T_Result <- t.test(yshift ~ member, data = dfr_yshift)
if(T_Result$p.value==0){
    T_pvalue <- c("p-value < 2.2e-16")
}else{
    T_pvalue <- signif(T_Result$p.value, digits = 3)
}

#### u-Test####
U_Result <- wilcox.test(yshift ~ member, data = dfr_yshift)
if(U_Result$p.value==0){
    U_pvalue <- c("p-value < 2.2e-16")
}else{
    U_pvalue <- signif(U_Result$p.value, digits = 3)
}

#### F-Test####
F_Result <- var.test(yshift ~ member, data = dfr_yshift)
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
                theme = ttheme(base_size = 120, base_style ="mBlue")
    ) -> gg_table
#### graph title####
# args_output |> 
#     str_remove("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_vis_all/") |> 
#     str_remove(".png") -> plot_title
plot_title <- c("Membership_vis_all")
#### patchwork####
gg <- gg_box +
    gg_vio -
    gg_table +
    plot_layout(ncol = 1, heights = c(3, 1)) +
    plot_annotation(title = plot_title,
                    caption = 'made with patchwork',
                    theme = theme(plot.title = element_text(size = 120, hjust = 0.5))
    )
#### ggsave####
ggsave(filename = args_output, 
       plot = gg,
       dpi = 50, 
       width = 50.0, 
       height = 50.0,
       limitsize = FALSE)