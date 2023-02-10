
# data <- list(ACF = prot1, Clustering = prot2, Proximity = prot3)
data_5 <- list(ACF = prot1,
               Clustering = prot2,
               Proximity = prot3,
               Zaslaver_et_al = prot4,
               Kaufman_et_al = prot5
               )

#### ggvenn package####
install.packages("ggvenn")
library(ggvenn)

prot4_or_prot5 <- union(prot4,prot5)

data_4 <- list(ACF = prot1,
               Clustering = prot2,
               Proximity = prot3,
               Zaslaver_Kaufman = prot4_or_prot5
               )
gg_F <- ggvenn(data_4,
               show_elements = FALSE,
               text_size = 5,
               stroke_size = 0.5,
               set_name_size = 3,
               label_sep = "\n")
# + theme(plot.margin= unit(c(2, 2, 2, 2), "lines"))
gg_T <- ggvenn(data_4,
               show_elements = TRUE,
               text_size = 5,
               stroke_size = 0.5,
               set_name_size = 3,
               label_sep = "\n")


# ggsave
ggsave(filename = "output/WTS3/normalize_1/stimAfter/SBD_abs/ARI/venn5/WTS3_venn4_count.png",
       plot = gg_F,
       width = 10.0,
       height = 8.0,
       limitsize = FALSE
       )
ggsave(filename = "output/WTS3/normalize_1/stimAfter/SBD_abs/ARI/venn5/WTS3_venn4_element.png",
       plot = gg_T,
       width = 10.0,
       height = 8.0,
       limitsize = FALSE
       )