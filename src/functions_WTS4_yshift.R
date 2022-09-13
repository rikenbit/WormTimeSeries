#library
##################################################
library(tidyverse)
library(openxlsx)
library(dtwclust)
library(usedist)
##################################################
# setting Random number generator
RNGkind(kind = "Mersenne-Twister")

.ReadData_stimAfter = function(x,y) {
    #### load stim timing extra####
    stimtimng_sheet <- read.xlsx(y,
                                 sheet = "Sheet1",
                                 rowNames = FALSE,
                                 colNames =TRUE)
    stimtimng_sheet %>% 
        dplyr::select(sample_number = 1, 
                      stim_first = 7,
        ) %>% 
        filter(sample_number == args_sample) %>% 
        .$stim_first %>% 
        trunc() -> stimtimng #切り捨て
    # ceiling() -> stimtimng #切り上げ
    return_object <- x[stimtimng:nrow(x),]
    return(return_object)
}

.ReadData_all = function(x) {
    return_object <- x
    return(return_object)
}


.mSBD <- function(x, y, znorm = FALSE, error.check = TRUE,
                  return.shifted = TRUE) {
    nx <- length(x)
    ny <- length(y)
    if (nx > ny) {
        swap <- x
        x <- y
        y <- swap
    } else {
        swap <- NULL
    }
    if (znorm) {
        CCseq <- NCCc(zscore(x, error.check = FALSE),
                      zscore(y, error.check = FALSE),
                      error.check = FALSE
        )
    } else {
        CCseq <- NCCc(x, y, error.check = FALSE)
    }
    m <- max(abs(CCseq))
    if (!return.shifted) {
        return(1 - m) # nocov
    }
    shift <- which.max(abs(CCseq)) - max(nx, ny)
    if (is.null(swap)) {
        if (shift < 0L) {
            yshift <- y[(-shift + 1L):ny]
        } else {
            yshift <- c(rep(0, shift), y)
        }
    } else {
        if (shift < 0L) {
            yshift <- c(rep(0, -shift), x)
        } else {
            yshift <- x[(shift + 1L):ny]
        }
    }
    nys <- length(yshift)
    if (nys < nx) {
        yshift <- c(yshift, rep(0, nx - nys))
    } else {
        yshift <- yshift[1L:nx]
    }
    # return
    list(dist = 1 - m, yshift = yshift, shift_value = shift)
}