#library
##################################################
# MC-MI-HOOI
library("rTensor")
library("einsum")
##################################################

##################################################
# Approach 1: Consensus Clustering
##################################################
CSPA <- function(Hs){
	# Merge
	mergedHs <- as.matrix(unlist(Hs))
	dim(mergedHs) <- c(nrow(Hs[[1]]), length(mergedHs) / nrow(Hs[[1]]))
	mergedHs %*% t(mergedHs) / length(Hs)
}

##################################################
# Approarch 2: Tensor Decompotision → Clustering
##################################################
# Orthogonal Individual Differences Scaling (O-INDSCAL)
OINDSCAL <- function(S, k, num.iter=30, thr=1E-10, verbose=FALSE){
    # Initialization
    initS <- as.matrix(unlist(S))
    dim(initS) <- c(nrow(S[[1]]), length(initS) / nrow(S[[1]]))
    X <- svd(initS)$u[, seq(k)]
  	D <- .calcD(S, X)
  	G <- .calcG(S, X, D)
  	# Iteration
  	iter <- 1
  	RecError <- c()
  	RelChange <- c()
  	RecError[iter] <- 10^2
  	RelChange[iter] <- 10^2
  	while(RelChange[iter] >= thr && iter <= num.iter){
  		X <- .calcX(G)
	  	D <- .calcD(S, X)
	  	G <- .calcG(S, X, D)
		iter <- iter + 1
		S_bar <- .recMatrices(X, D)
		RecError[iter] <- .recError(S, S_bar)
		RelChange[iter] <- abs(RecError[iter-1] - RecError[iter]) / RecError[iter]
        if(verbose){
             cat(paste0(iter - 1, " / ", num.iter,
                " |Previous Error - Error| / Error = ", RelChange[iter], "\n"))
        }
  	}
  	# Output
  	list(X=X, D=D, G=G, RecError=RecError, RelChange=RelChange)
}

.recMatrices <- function(X, D){
	lapply(seq_along(D), function(s){
		X %*% D[[s]] %*% t(X)
	})	
}

.recError <- function(S, S_bar){
	out <- lapply(seq_along(S), function(s){
		S[[s]] - S_bar[[s]]
	})
	sqrt(sum(unlist(out)^2))
}

.calcX <- function(G){
	res <- svd(G)
	res$u %*% t(res$v)
}

.calcD <- function(S, X){
	lapply(S, function(Ss){
		out <- matrix(0, nrow=ncol(X), ncol=ncol(X))
		diag(out) <- diag(t(X) %*% Ss %*% X)
		out[which(out < 0)] <- 0
		out
	})
}

.calcG <- function(S, X, D){
	out <- lapply(seq_along(S), function(s){
		S[[s]] %*% X %*% D[[s]]
	})
	Reduce("+", out)
}

##################################################
# Approarch 3: MC-MI-HOOI with Graph Laplacian
##################################################
MCMIHOOI <- function(A, k, num.iter=30, thr=1E-20, verbose=FALSE){
	# Initialization
	A1 <- rs_unfold(as.tensor(A), m=1)@data
	U <- svd(A1)$u[, seq(k)]
  	# Iteration
  	iter <- 1
  	RecError <- c()
  	RelChange <- c()
  	RecError[iter] <- 10^2
  	RelChange[iter] <- 10^2
  	while(RelChange[iter] >= thr && iter <= num.iter){
		# Step1
		A3 <- rs_unfold(as.tensor(A), m=3)@data
		W <- as.matrix(svd(A3 %*% kronecker(U, U))$u[,1])
		# Step2
		S_ <- array(0, dim=dim(A)[1:2])
		for(i in seq_len(dim(A)[3])){
			S_ <- S_ + W[i] * A[,,i]
		}
		# Step3
		U <- svd(S_)$u[,seq(k)]
		# Update
		iter <- iter + 1
		G <- einsum('ijk,il,jm,kn->lmn', A, U, U, W)
		A_bar <- einsum('lmn,il,jm,kn->ijk', G, U, U, W)
		RecError[iter] <- .recError(A, A_bar)
		RelChange[iter] <- abs(RecError[iter-1] - RecError[iter]) / RecError[iter]
        if(verbose){
             cat(paste0(iter - 1, " / ", num.iter,
                " |Previous Error - Error| / Error = ", RelChange[iter], "\n"))
        }
	}
	# Output
	list(U=U, W=W, G=G)
}