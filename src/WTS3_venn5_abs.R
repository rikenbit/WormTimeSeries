# prot1 <- c("AIBR","ASEL","ASER","ASGR","ASHL","ASHR","ASIR","ASJL","ASKR","AUAL","AWAL","AWAR","AWBL","AWBR","AWCL","AWCR","BAGL","BAGR","CEPVL","CEPVR","OLQVR","RIAR","RIR","SIBDR","URXL","URXR")
prot1 <- c("ASER","BAGL","BAGR","AWCR","AWCL","AWBL","ASEL","AWBR","ASIR","ASKR","ASHL","OLQVR","AWAL","ASGR","ASHR","SIBDR","RIAR","AUAL","AWAR","URXL","URXR","CEPVR","AIBR","CEPVL")

# prot2 <- c("AIMR","ASEL","ASER","ASHL","ASHR","ASIL","ASKR","AUAL","AUAR","AWBL","AWCR","BAGL","BAGR","OLLL","OLQVR","RIVR")
prot2 <- c("ASER","BAGR","BAGL","ASEL","AWBL","AWCR","ASHR","ASHL","AWAL","RMED","AWCL","AUAL","AUAR","I3","AWBR")

# prot3 <- c("ASHL","ASIL","ASKR","AUAL","AUAR","AWAL","BAGL","BAGR","CEPVR","OLLL","OLQVR","RIVR","SIBDR","URXL","URXR")
prot3 <- c("BAGL","BAGR","ASKR","OLQVR","ASHL","AWAL","SIBDR","AUAL","URXL","URXR","CEPVR","AUAR","OLLL")

# Zaslaver et al. 2015
prot4 <- c("ASEL","ASER","AWCL","AWCR","AWBL","AWBR","AFDL","AFDR","ASHL","ASHR","ASJL","ASJR")
# Kaufman et al. 2005
prot5 <- c("ASEL","ASER","ADFL","ADFR","ASHL","ASHR","ASIL","ASIR","ASJL","ASJR","ASKL","ASKR","ADLL","ADLR")

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