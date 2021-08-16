source("src/functions_WTS3_venn.R")

prot1 <- c("AIBR","ASEL","ASER","ASGR","ASHL","ASHR","ASIR","ASJL","ASKR","AUAL","AWAL","AWAR","AWBL","AWBR","AWCL","AWCR","BAGL","BAGR","CEPVL","CEPVR","OLQVR","RIAR","RIR","SIBDR","URXL","URXR")
prot2 <- c("AIMR","ASEL","ASER","ASHL","ASHR","ASIL","ASKR","AUAL","AUAR","AWBL","AWCR","BAGL","BAGR","OLLL","OLQVR","RIVR")
prot3 <- c("ASHL","ASIL","ASKR","AUAL","AUAR","AWAL","BAGL","BAGR","CEPVR","OLLL","OLQVR","RIVR","SIBDR","URXL","URXR")
data <- list(ACF = prot1, Clustering = prot2, Proximity = prot3)

# Venn ACF VS クラスタ VS 空間的近さ
##############################################################
venn.diagram(
	data,
	fill = c(3, 2, 1),           # background color
	alpha = c(0.5, 0.5, 0.5),   # transparency
	lty = 1,            # border line type
	margin = 0.1 ,         #ベン図画像の余白の指定
	cat.cex=c(1, 1, 1),      #カテゴリー名のフォントサイズの指定
	cat.dist=c(0.1, 0.1, 0.1), #カテゴリー名の枠線からの距離の指定
	filename = "output/WTS3/normalize_1/stimAfter/SBD/ARI/venn/WTS3_venn.tiff"       # file name
	)
##############################################################

intersect(prot3,intersect(prot1,prot2))
intersect(prot3,intersect(prot1,prot2))
w <- venn.diagram(
	data,
	fill = c(3, 2, 1),           # background color
	alpha = c(0.5, 0.5, 0.5),   # transparency
	lty = 1,            # border line type
	margin = 0.1 ,         #ベン図画像の余白の指定
	cat.cex=c(1, 1, 1),      #カテゴリー名のフォントサイズの指定
	cat.dist=c(0.1, 0.1, 0.1), #カテゴリー名の枠線からの距離の指定
	filename = NULL       # file name
	)
grid.newpage()
grid.draw(w)

w[[7]]$label  <- paste(setdiff(prot1,union(prot2,prot3)),collapse="\n")
# # Intersecção A, B e C ou AnBnC
# inters <- intersect(aa,intersect(bb,cc))
# # A inter B - inters
# w[[8]]$label <- paste(setdiff(intersect(aa,bb), inters),collapse = "")
# # B-(A+C)
# w[[9]]$label <- paste(setdiff(bb,union(aa,cc)),collapse = "")
# # A inter C - inters
# w[[10]]$label <- paste(setdiff(intersect(aa,cc), inters),collapse = "")
# # Intersecção A, B e C: posso usar inters ou ww[[1,5]]
# w[[11]]$label <- paste(ww[[1,5]], collapse = "")
# # B inter C - inters
# w[[12]]$label <- paste(setdiff(intersect(cc,bb), inters),collapse = "")
# # C-(A+B)
# w[[13]]$label <- paste(setdiff(cc,union(aa,bb)),collapse = "")

# # plot  
# grid.newpage()
# grid.draw(w)

# png("venn1.png", 300, 300)
# venn(data)

