source("src/functions_WTS4_Membership_vis_all.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input <- args[1]
args_output <- args[2]
args_output_csv <- args[3]
#### test args####
# args_input <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_6/DFs.RData")
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_vis_all/k_Number_6.png")
# args_output_csv <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_vis_all/k_Number_6.csv")

#### load####
load(args_input)

#### yshift list to df####
dfr_yshift <- purrr::map_dfr(names(DFs), .DFs_yshift)

# #### ggplot####
# cord_x <- c("Membership")
# cord_y <- c("yshift")

#### ggplot density all fill_Member####
# gg_dens <- ggplot(dfr_yshift,
#                   aes(x = yshift,
#                       fill=factor(member))
# ) +
#     geom_density(alpha=0.2) +
#     labs(x = cord_x,
#          y = cord_y,
#          fill = "Membership") +
#     ylim(c(0, 0.05)) +
#     xlim(c(-1000, 1000))

# #### ggplot density all fill_Member facet_wrap####
# gg_dens <- ggplot(dfr_yshift,
#                   aes(x = yshift)
# ) +
#     geom_density(fill="blue", alpha=0.2) +
#     labs(x = cord_x,
#          y = cord_y,
#          fill = "Animal_all") +
#     ylim(c(0, 0.05)) +
#     xlim(c(-1000, 1000)) +
#     facet_wrap(member ~ ., nrow=2)

#### ggplot density group facet_wrap####
cord_x <- c("yshift")
cord_y <- c("density")

gg_dens <- ggplot(dfr_yshift,
             aes(x = yshift,
                 fill=factor(animal, levels = names(DFs))
                 )
             ) +
    geom_density(alpha=0.2) +
    # geom_density(fill="blue", alpha=0.2) +
    # labs(x = cord_x,
    #      y = cord_y,
    #      fill = "Membership") +
    labs(x = cord_x,
         y = cord_y,
         fill = "Animal_No.") +
    # ylim(c(0, 0.02)) +
    xlim(c(-1000, 1000)) +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 60),
          legend.key.size = unit(1, 'cm')
          ) +
    facet_wrap(member ~ ., nrow=2)

#### ggplot density group####
# gg_dens <- ggplot(dfr_yshift[dfr_yshift$member==1,], 
#                   aes(x = yshift,
#                       fill=factor(animal, levels = names(DFs))
#                       )
#                   ) + 
#     geom_density(alpha=0.2) +
#     labs(x = cord_x,
#          y = cord_y,
#          fill = "Animal_No.") +
#     ylim(c(0, 1.00)) +
#     xlim(c(-1000, 1000)) +
#     theme(plot.title = element_text(hjust = 0.5),
#           text = element_text(size = 60),
#           legend.key.size = unit(1, 'cm')
#     ) 


#### t-Test####
# ggdata
T_Result <- t.test(yshift ~ member, data = dfr_yshift)
if(T_Result$p.value==0){
    T_pvalue <- c("p-value < 2.2e-16")
    T_qvalue <- c("q-value < 6.6e-16")
}else{
    T_pvalue <- signif(T_Result$p.value, digits = 3)
    T_qvalue <- p.adjust(T_pvalue, "BH", n = 3)
}

#### u-Test####
U_Result <- wilcox.test(yshift ~ member, data = dfr_yshift)
if(U_Result$p.value==0){
    U_pvalue <- c("p-value < 2.2e-16")
    U_qvalue <- c("q-value < 6.6e-16")
}else{
    U_pvalue <- signif(U_Result$p.value, digits = 3)
    U_qvalue <- p.adjust(U_pvalue, "BH", n = 3)
}

#### F-Test####
F_Result <- var.test(yshift ~ member, data = dfr_yshift)
if(F_Result$p.value==0){
    F_pvalue <- c("p-value < 2.2e-16")
    F_qvalue <- c("q-value < 6.6e-16")
}else{
    F_pvalue <- signif(F_Result$p.value, digits = 3)
    F_qvalue <- p.adjust(F_pvalue, "BH", n = 3)
}
#### SD####
dfr_yshift[dfr_yshift$member=="0",]$yshift |> 
    sd() |> 
    round() ->SD_0 
dfr_yshift[dfr_yshift$member=="1",]$yshift |> 
    sd() |> 
    round() ->SD_1
#### Test table####
data.frame(T_pvalue = T_pvalue,
           U_pvalue = U_pvalue,
           F_pvalue = F_pvalue,
           SD_Mem0 = SD_0,
           SD_Mem1 = SD_1,
           stringsAsFactors = FALSE,
           row.names = NULL) -> t_table
t_table %>% 
    ggtexttable(rows = NULL, 
                theme = ttheme(base_size = 50, 
                               base_style ="mBlue",
                               padding = unit(c(10, 10), "mm")
                               )
    ) -> gg_table
#### graph title####
plot_title <- c("Animal_all")
#### patchwork####
# gg <- gg_box -
#     # gg_box +
#     # gg_vio -
#     gg_table +
#     plot_layout(ncol = 1, heights = c(5, 1)) +
#     plot_annotation(title = plot_title,
#                     caption = 'made with patchwork',
#                     theme = theme(plot.title = element_text(size = 60, hjust = 0.5))
#     )

gg <- gg_dens -
    gg_table +
    plot_layout(ncol = 1, heights = c(2, 1)) +
    plot_annotation(title = plot_title,
                    caption = 'made with patchwork',
                    theme = theme(plot.title = element_text(size = 60, hjust = 0.5))
    )
#### ggsave####
ggsave(filename = args_output, 
       plot = gg,
       dpi = 50, 
       width = 30.0, 
       height = 30.0,
       limitsize = FALSE)
# ggsave(filename = args_output, 
#        plot = gg_dens,
#        dpi = 50, 
#        width = 30.0, 
#        height = 30.0,
#        limitsize = FALSE)

#### ggsave p q value####
purrr::map_dfr(names(DFs), .DFs_test) -> Test_table
column_to_rownames(Test_table, "Animal") -> Test_table 
Test_table_LOG <- Test_table
Test_table_LOG[,1:6] <- -log10(Test_table_LOG[,1:6]) 
write.csv(Test_table_LOG, 
          args_output_csv, 
          row.names=TRUE)