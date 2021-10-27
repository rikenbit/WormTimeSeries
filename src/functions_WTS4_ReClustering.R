#library
##################################################
library()
##################################################
.func = function(x) {
	return(return_object)
}

CSPA <- function(Hs){
	# Merge
	mergedHs <- as.matrix(unlist(Hs))
	dim(mergedHs) <- c(nrow(Hs[[1]]), length(mergedHs) / nrow(Hs[[1]]))
	mergedHs %*% t(mergedHs) / length(Hs)
}