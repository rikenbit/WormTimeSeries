source("src/functions_WTS3_venn.R")

prot1 <- c("AIBR","ASEL","ASER","ASGR","ASHL","ASHR","ASIR","ASJL","ASKR","AUAL","AWAL","AWAR","AWBL","AWBR","AWCL","AWCR","BAGL","BAGR","CEPVL","CEPVR","OLQVR","RIAR","RIR","SIBDR","URXL","URXR")
prot2 <- c("AIMR","ASEL","ASER","ASHL","ASHR","ASIL","ASKR","AUAL","AUAR","AWBL","AWCR","BAGL","BAGR","OLLL","OLQVR","RIVR")
prot3 <- c("ASHL","ASIL","ASKR","AUAL","AUAR","AWAL","BAGL","BAGR","CEPVR","OLLL","OLQVR","RIVR","SIBDR","URXL","URXR")
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

#### venn package####
venn(5, ilab=TRUE, zcolor = "style")
venn(7, ilab=TRUE, zcolor = "style")
fiveCellVenn()

venn(data, ilab=TRUE, zcolor = "style")
venn("A", snames = "A, B, C")
venn("A, E", snames = "A, B, C, D, E", zcolor = "red, blue")
venn("A + B~C", snames = "A, B, C, D")
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

data_2 <- list(Zaslaver_et_al_2015 = prot4, 
               Kaufman_et_al_2005 = prot5
               )
gg_article_F <- ggvenn(data_2, 
                       show_elements = FALSE,
                       text_size = 5,
                       stroke_size = 0.5,
                       set_name_size = 6,
                       label_sep = "\n")

gg_article_T <- ggvenn(data_2, 
                       show_elements = TRUE, 
                       text_size = 5, 
                       stroke_size = 0.5,
                       set_name_size = 6,
                       label_sep = "\n")

# ggsave
ggsave(filename = "output/WTS3/normalize_1/stimAfter/SBD/ARI/venn5/WTS3_venn4_count.png", 
       plot = gg_F,
       width = 10.0, 
       height = 8.0,
       limitsize = FALSE
       )
ggsave(filename = "output/WTS3/normalize_1/stimAfter/SBD/ARI/venn5/WTS3_venn4_element.png", 
       plot = gg_T,
       width = 10.0, 
       height = 8.0,
       limitsize = FALSE
       )

ggsave(filename = "output/WTS3/normalize_1/stimAfter/SBD_abs/ARI/venn5/WTS3_venn2_count.png", 
       plot = gg_article_F,
       width = 10.0, 
       height = 8.0,
       limitsize = FALSE
       )
ggsave(filename = "output/WTS3/normalize_1/stimAfter/SBD_abs/ARI/venn5/WTS3_venn2_element.png", 
       plot = gg_article_T,
       width = 10.0, 
       height = 8.0,
       limitsize = FALSE
       )

#### VennDiagram####
library("VennDiagram")
# move to new plotting page
grid.newpage()
# create pairwise Venn diagram
draw.quintuple.venn(prot1,prot2,prot3,prot4,prot5)

library(VennDiagram)

# グループを準備
# a <- read.csv('Sample1.csv')
# data1 <- as.vector(a$id)
# b <- read.csv('Sample2.csv')
# data2 <- as.vector(b$id)
# c <- read.csv('Sample3.csv')
# data3 <- as.vector(c$id)
# d <- read.csv('Sample4.csv')
# data4 <- as.vector(d$id)
# e <- read.csv('Sample5.csv')
# data5 <- as.vector(e$id)

# データをリスト型に変換
data_VennD <- list("ACF" = prot1, 
                   "Clustering" = prot2, 
                   "Proximity" = prot3, 
                   "Zaslaver_et_al" = prot4, 
                   "Kaufman_et_al" = prot5)

# Plot
# pdf("Benn.pdf")
venn.diagram(data_VennD, 
             filename="output/WTS3/normalize_1/stimAfter/SBD/ARI/venn5/WTS3_VennD.tiff", 
             margin = 0.5,
             fill=c(2,3,4,5,6), 
             alpha=0.4, 
             lty=1)





