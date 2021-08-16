#library
##################################################
# library(dtwclust)
# RNGkind(kind = "Mersenne-Twister")
library(tidyverse)
library(openxlsx)
##################################################
# setting Random number generator
# RNGkind(kind = "Mersenne-Twister")
# dtwclust_rngkind <- "Mersenne-Twister"

# https://github.com/asardaes/dtwclust/blob/f4a5978ce051585d3c989e658c84d6b4b0fb77d4/R/DISTANCES-sbd.R#L96
#### SBD_abs####
#' Shape-based distance
#'
#' Distance based on coefficient-normalized cross-correlation as proposed by Paparrizos and Gravano
#' (2015) for the k-Shape clustering algorithm.
#'
#' @export
#'
#' @param x,y Univariate time series.
#' @param znorm Logical. Should each series be z-normalized before calculating the distance?
#' @template error-check
#' @param return.shifted Logical. Should the shifted version of `y` be returned? See details.
#'
#' @details
#'
#' This distance works best if the series are *z-normalized*. If not, at least they should have
#' appropriate amplitudes, since the values of the signals **do** affect the outcome.
#'
#' If `x` and `y` do **not** have the same length, it would be best if the longer sequence is
#' provided in `y`, because it will be shifted to match `x`. After matching, the series may have to
#' be truncated or extended and padded with zeros if needed.
#'
#' The output values lie between 0 and 2, with 0 indicating perfect similarity.
#'
#' @return For `return.shifted = FALSE`, the numeric distance value, otherwise a list with:
#'
#' - `dist`: The shape-based distance between `x` and `y`.
#' - `yshift`: A shifted version of `y` so that it optimally matches `x` (based on [NCCc()]).
#'
#' @template proxy
#' @template symmetric
#' @section Proxy version:
#'
#'   In some situations, e.g. for relatively small distance matrices, the overhead introduced by the
#'   logic that computes only half the distance matrix can be bigger than just calculating the whole
#'   matrix.
#'
#' @note
#'
#' If you wish to calculate the distance between several time series, it would be better to use the
#' version registered with the `proxy` package, since it includes some small optimizations. See the
#' examples.
#'
#' This distance is calculated with help of the Fast Fourier Transform, so it can be sensitive to
#' numerical precision. Thus, this function (and the functions that depend on it) might return
#' different values in 32 bit installations compared to 64 bit ones.
#'
#' @references
#'
#' Paparrizos J and Gravano L (2015). ``k-Shape: Efficient and Accurate Clustering of Time Series.''
#' In *Proceedings of the 2015 ACM SIGMOD International Conference on Management of Data*, series
#' SIGMOD '15, pp. 1855-1870. ISBN 978-1-4503-2758-9, \url{http://doi.org/10.1145/2723372.2737793}.
#'
#' @seealso
#'
#' [NCCc()], [shape_extraction()]
#'
#' @examples
#'
#' # load data
#' data(uciCT)
#'
#' # distance between series of different lengths
#' sbd <- SBD(CharTraj[[1]], CharTraj[[100]], znorm = TRUE)$dist
#'
#' # cross-distance matrix for series subset (notice the two-list input)
#' sbD <- proxy::dist(CharTraj[1:10], CharTraj[1:10], method = "SBD", znorm = TRUE)
#'
SBD <- function(x, y, znorm = FALSE, error.check = TRUE, return.shifted = TRUE) {
    if (is_multivariate(list(x, y))) stop("SBD does not support multivariate series.")
    if (error.check) {
        check_consistency(x, "ts")
        check_consistency(y, "ts")
        x <- as.numeric(x)
        y <- as.numeric(y)
    }
    nx <- length(x)
    ny <- length(y)
    
    if (nx > ny) {
        # The order in which I provide the arguments to NCCc affects 'shift'
        swap <- x
        x <- y
        y <- swap
    }
    else {
        swap <- NULL
    }
    
    if (znorm)
        CCseq <- NCCc(zscore(x, error.check = FALSE),
                      zscore(y, error.check = FALSE),
                      error.check = FALSE)
    else
        CCseq <- NCCc(x, y, error.check = FALSE)
    #### test code####
    CCseq <- CCseq[CCseq < 0.1]
    ##########
    # m <- max(CCseq)
    
    # improvement
    m <- max(abs(CCseq))
    if (!return.shifted) return(1 - m) # nocov
    # shift <- which.max(CCseq) - max(nx, ny)
    shift <- which.max(abs(CCseq)) - max(nx, ny)
    
    if (is.null(swap)) {
        if (shift < 0L)
            yshift <- y[(-shift + 1L):ny]
        else
            yshift <- c( rep(0, shift), y )
    }
    else {
        # Remember, if I swapped them, then I have to shift what is now saved in 'x'
        if (shift < 0L)
            yshift <- c( rep(0, -shift), x )
        else
            yshift <- x[(shift + 1L):ny]
    }
    
    nys <- length(yshift)
    if (nys < nx)
        yshift <- c( yshift, rep(0, nx-nys) )
    else
        yshift <- yshift[1L:nx]
    # return
    list(dist = 1 - m, yshift = yshift)
}

#' @rdname SBD
#' @export
#'
sbd <- SBD

# ==================================================================================================
# Wrapper for proxy::dist
# ==================================================================================================

#' @importFrom stats fft
#' @importFrom stats nextn
#'
sbd_proxy <- function(x, y = NULL, znorm = FALSE, ..., error.check = TRUE, pairwise = FALSE) {
    x <- tslist(x)
    
    if (error.check) check_consistency(x, "vltslist")
    if (znorm) x <- zscore(x, ..., error.check = FALSE) # nocov
    
    if (is.null(y)) {
        symmetric <- TRUE
        y <- x
        
        # Precompute FFTs, padding with zeros as necessary, which will be compensated later
        L <- max(lengths(x)) * 2L - 1L
        fftlen <- stats::nextn(L, 2L)
        fftx <- lapply(x, function(u) { stats::fft(c(u, rep(0, fftlen - length(u)))) })
        ffty <- lapply(fftx, Conj)
    }
    else {
        symmetric <- FALSE
        y <- tslist(y)
        if (error.check) check_consistency(y, "vltslist")
        if (znorm) y <- zscore(y, ..., error.check = FALSE) # nocov
        
        # Precompute FFTs, padding with zeros as necessary, which will be compensated later
        L <- max(lengths(x)) + max(lengths(y)) - 1L
        fftlen <- stats::nextn(L, 2L)
        fftx <- lapply(x, function(u) { stats::fft(c(u, rep(0, fftlen - length(u)))) })
        ffty <- lapply(y, function(v) { Conj(stats::fft(c(v, rep(0, fftlen - length(v))))) })
    }
    
    if (is_multivariate(c(x,y))) stop("SBD does not support multivariate series.") # nocov
    fill_type <- mat_type <- dim_names <- NULL # avoid warning about undefined globals
    eval(prepare_expr) # UTILS-expressions.R
    
    # calculate distance matrix
    distance <- "SBD" # read in C++, can't be temporary!
    distargs <- list(
        fftlen = fftlen,
        fftx = fftx,
        ffty = ffty
    )
    num_threads <- get_nthreads()
    .Call(C_distmat_loop,
          D, x, y, distance, distargs, fill_type, mat_type, num_threads,
          PACKAGE = "dtwclust")
    
    # adjust D's attributes
    if (pairwise) {
        dim(D) <- NULL
        class(D) <- "pairdist"
    }
    else {
        dimnames(D) <- dim_names
        class(D) <- "crossdist"
    }
    attr(D, "method") <- "SBD"
    # return
    D
}


########

# https://github.com/asardaes/dtwclust/blob/f4a5978ce051585d3c989e658c84d6b4b0fb77d4/R/CLUSTERING-tsclust.R
#### tsclust####
# ==================================================================================================
# Main function
# ==================================================================================================

#' Time series clustering
#'
#' This is the main function to perform time series clustering. See the details and the examples for
#' more information, as well as the included package vignettes (which can be found by typing
#' `browseVignettes("dtwclust")`). A convenience wrapper is available in [compare_clusterings()],
#' and a shiny app in [interactive_clustering()].
#'
#' @export
#' @importFrom dtw symmetric1
#' @importFrom dtw symmetric2
#' @importFrom methods new
#' @importFrom methods slot
#' @importFrom parallel splitIndices
#' @importFrom proxy pr_DB
#' @importFrom stats as.dist
#' @importFrom stats as.hclust
#' @importFrom stats cutree
#' @importFrom stats hclust
#'
#' @param series A list of series, a numeric matrix or a data frame. Matrices and data frames are
#'   coerced to a list row-wise (see [tslist()]).
#' @param type What type of clustering method to use: `"partitional"`, `"hierarchical"`, `"tadpole"`
#'   or `"fuzzy"`.
#' @param k Number of desired clusters. It can be a numeric vector with different values.
#' @param ... Arguments to pass to preprocessing, centroid **and** distance functions (added to
#'   `args`). Also passed to `method` from [hierarchical_control()] if it happens to be a function,
#'   and to [stats::hclust()] if it contains the `members` parameter.
#' @param preproc Function to preprocess data. Defaults to [zscore()] *only* if `centroid` `=`
#'   `"shape"`, but will be replaced by a custom function if provided.
#' @param distance A registered distance from [proxy::dist()]. Ignored for `type` `=` `"tadpole"`.
#' @param centroid Either a supported string, or an appropriate function to calculate centroids when
#'   using partitional/hierarchical/tadpole methods. See Centroids section.
#' @param control An appropriate list of controls. See [tsclust-controls].
#' @param args An appropriate list of arguments for preprocessing, distance and centroid functions.
#'   See [tsclust_args()] and the examples.
#' @param seed Random seed for reproducibility.
#' @param trace Logical flag. If `TRUE`, more output regarding the progress is printed to screen.
#' @template error-check
#'
#' @details
#'
#' Partitional and fuzzy clustering procedures use a custom implementation. Hierarchical clustering
#' is done with [stats::hclust()] by default. TADPole clustering uses the [TADPole()] function.
#' Specifying `type` = `"partitional"`, `preproc` = `zscore`, `distance` = `"sbd"` and `centroid` =
#' `"shape"` is equivalent to the k-Shape algorithm (Paparrizos and Gravano 2015).
#'
#' The `series` may be provided as a matrix, a data frame or a list. Matrices and data frames are
#' coerced to a list, both row-wise. Only lists can have series with different lengths or multiple
#' dimensions. Most of the optimizations require series to have the same length, so consider
#' reinterpolating them to save some time (see Ratanamahatana and Keogh 2004; [reinterpolate()]). No
#' missing values are allowed.
#'
#' In the case of multivariate time series, they should be provided as a list of matrices, where
#' time spans the rows of each matrix and the variables span the columns (see [CharTrajMV] for an
#' example). All included centroid functions should work with the aforementioned format, although
#' `shape` is *not* recommended. Note that the `plot` method will simply append all dimensions
#' (columns) one after the other.
#'
#' @return
#'
#' An object with an appropriate class from [TSClusters-class].
#'
#' If `control$nrep > 1` and a partitional procedure is used, `length(method)` `> 1` and
#' hierarchical procedures are used, or `length(k)` `>` `1`, a list of objects is returned.
#'
#' @section Centroid Calculation:
#'
#'   In the case of partitional/fuzzy algorithms, a suitable function should calculate the cluster
#'   centroids at every iteration. In this case, the centroids may also be time series. Fuzzy
#'   clustering uses the standard fuzzy c-means centroid by default.
#'
#'   In either case, a custom function can be provided. If one is provided, it will receive the
#'   following parameters with the shown names (examples for partitional clustering are shown in
#'   parentheses):
#'
#'   - `x`: The *whole* data list (`list(ts1, ts2, ts3)`)
#'   - `cl_id`: An integer vector with length equal to the number of series in `data`, indicating
#'   which cluster a series belongs to (`c(1L, 2L, 2L)`)
#'   - `k`: The desired number of total clusters (`2L`)
#'   - `cent`: The current centroids in order, in a list (`list(centroid1, centroid2)`)
#'   - `cl_old`: The membership vector of the *previous* iteration (`c(1L, 1L, 2L)`)
#'   - The elements of `...` that match its formal arguments
#'
#'   In case of fuzzy clustering, the membership vectors (2nd and 5th elements above) are matrices
#'   with number of rows equal to amount of elements in the data, and number of columns equal to the
#'   number of desired clusters. Each row must sum to 1.
#'
#'   The other option is to provide a character string for the custom implementations. The following
#'   options are available:
#'
#'   - "mean": The average along each dimension. In other words, the average of all \eqn{x^j_i}
#'   among the \eqn{j} series that belong to the same cluster for all time points \eqn{t_i}.
#'   - "median": The median along each dimension. Similar to mean.
#'   - "shape": Shape averaging. By default, all series are z-normalized in this case, since the
#'   resulting centroids will also have this normalization. See [shape_extraction()] for more
#'   details.
#'   - "dba": DTW Barycenter Averaging. See [DBA()] for more details.
#'   - "sdtw_cent": Soft-DTW centroids, See [sdtw_cent()] for more details.
#'   - "pam": Partition around medoids (PAM). This basically means that the cluster centroids are
#'   always one of the time series in the data. In this case, the distance matrix can be
#'   pre-computed once using all time series in the data and then re-used at each iteration. It
#'   usually saves overhead overall for small datasets (see [tsclust-controls]).
#'   - "fcm": Fuzzy c-means. Only supported for fuzzy clustering and used by default in that case.
#'   - "fcmdd": Fuzzy c-medoids. Only supported for fuzzy clustering. It **always** precomputes/uses
#'   the whole cross-distance matrix.
#'
#'   The `dba`, `shape` and `sdtw_cent` implementations check for parallelization. Note that only
#'   `shape`, `dba`, `sdtw_cent`, `pam` and `fcmdd` support series of different length. Also note
#'   that for `shape`, `dba` and `sdtw_cent`, this support has a caveat: the final centroids' length
#'   will depend on the length of those series that were randomly chosen at the beginning of the
#'   clustering algorithm. For example, if the series in the dataset have a length of either 10 or
#'   15, 2 clusters are desired, and the initial choice selects two series with length of 10, the
#'   final centroids will have this same length.
#'
#'   As special cases, if hierarchical or tadpole clustering is used, you can provide a centroid
#'   function that takes a list of series as first input. It will also receive the contents of
#'   `args$cent` that match its formal arguments, and should return a single centroid series. These
#'   centroids are returned in the `centroids` slot. By default, the medoid of each cluster is
#'   extracted (similar to what [pam_cent()] does).
#'
#'   In the following cases, the `centroids` list will have an attribute `series_id` with an integer
#'   vector indicating which `series` were chosen as centroids:
#'
#'   - Partitional clustering using "pam" centroid.
#'   - Fuzzy clustering using "fcmdd" centroid.
#'   - Hierarchical clustering with the default centroid extraction.
#'   - TADPole clustering with the default centroid extraction.
#'
#' @section Distance Measures:
#'
#'   The distance measure to be used with partitional, hierarchical and fuzzy clustering can be
#'   modified with the `distance` parameter. The supported option is to provide a string, which must
#'   represent a compatible distance registered with `proxy`'s [proxy::dist()]. Registration is done
#'   via [proxy::pr_DB()], and extra parameters can be provided in `args$dist` (see the examples).
#'
#'   Note that you are free to create your own distance functions and register them. Optionally, you
#'   can use one of the following custom implementations (all registered with `proxy`):
#'
#'   - `"dtw"`: DTW, optionally with a Sakoe-Chiba/Slanted-band constraint. Done with [dtw::dtw()].
#'   - `"dtw2"`: DTW with L2 norm and optionally a Sakoe-Chiba/Slanted-band constraint. See
#'   [dtw2()].
#'   - `"dtw_basic"`: A custom version of DTW with less functionality, but faster. See
#'   [dtw_basic()].
#'   - `"dtw_lb"`: DTW with L1 or L2 norm and a Sakoe-Chiba constraint. Some computations are
#'   avoided by first estimating the distance matrix with Lemire's lower bound and then
#'   iteratively refining with DTW. See [dtw_lb()]. Not suitable for `pam.precompute` = `TRUE` nor
#'   hierarchical clustering.
#'   - `"lbk"`: Keogh's lower bound for DTW with either L1 or L2 norm for the Sakoe-Chiba
#'   constraint. See [lb_keogh()].
#'   - `"lbi"`: Lemire's lower bound for DTW with either L1 or L2 norm for the Sakoe-Chiba
#'   constraint. See [lb_improved()].
#'   - `"sbd"`: Shape-based distance. See [sbd()].
#'   - `"gak"`: Global alignment kernels. See [gak()].
#'   - `"sdtw"`: Soft-DTW. See [sdtw()].
#'
#'   Out of the aforementioned, only the distances based on DTW lower bounds *don't* support series
#'   of different length. The lower bounds are probably unsuitable for direct clustering unless
#'   series are very easily distinguishable.
#'
#'   If you know that the distance function is symmetric, and you use a hierarchical algorithm, or a
#'   partitional algorithm with PAM centroids, or fuzzy c-medoids, some time can be saved by
#'   calculating only half the distance matrix. Therefore, consider setting the symmetric control
#'   parameter to `TRUE` if this is the case (see [tsclust-controls]).
#'
#' @section Preprocessing:
#'
#'   It is strongly advised to use z-normalization in case of `centroid = "shape"`, because the
#'   resulting series have this normalization (see [shape_extraction()]). Therefore, [zscore()] is
#'   the default in this case. The user can, however, specify a custom function that performs any
#'   transformation on the data, but the user must make sure that the format stays consistent, i.e.
#'   a list of time series.
#'
#'   Setting to `NULL` means no preprocessing (except for `centroid = "shape"`). A provided function
#'   will receive the data as first argument, followed by the contents of `args$preproc` that match
#'   its formal arguments.
#'
#'   It is convenient to provide this function if you're planning on using the [stats::predict()]
#'   generic (see also [TSClusters-methods]).
#'
#' @section Repetitions:
#'
#'   Due to their stochastic nature, partitional clustering is usually repeated several times with
#'   different random seeds to allow for different starting points. This function uses
#'   [parallel::nextRNGStream()] to obtain different seed streams for each repetition, utilizing the
#'   `seed` parameter (if provided) to initialize it. If more than one repetition is made, the
#'   streams are returned in an attribute called `rng`.
#'
#'   Multiple values of `k` can also be provided to get different partitions using any `type` of
#'   clustering.
#'
#'   Repetitions are greatly optimized when PAM centroids are used and the whole distance matrix is
#'   precomputed, since said matrix is reused for every repetition.
#'
#' @template parallel
#'
#' @section Parallel Computing:
#'
#'   Multi-processing is used in partitional and fuzzy clustering for multiple values of `k` and/or
#'   `nrep` (in [partitional_control()]). See [TADPole()] to know how it uses parallelization. For
#'   cross-distance matrix calculations, the parallelization strategy depends on whether the
#'   distance is included with \pkg{dtwclust} or not, see the caveats in [tsclustFamily-class].
#'
#'   If you register a parallel backend and special packages must be loaded, provide their names in
#'   the `packages` element of `control`. Note that "dtwclust" is always loaded in each parallel
#'   worker, so that doesn't need to be included. Alternatively, you may want to pre-load
#'   \pkg{dtwclust} in each worker with [parallel::clusterEvalQ()].
#'
#' @note
#'
#' The lower bounds are defined only for time series of equal length. They are **not** symmetric,
#' and `DTW` is not symmetric in general.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @references
#'
#' Please refer to the package vignette references (which can be loaded by typing
#' `vignette("dtwclust")`).
#'
#' @seealso
#'
#' [TSClusters-class], [TSClusters-methods], [tsclustFamily-class], [tsclust-controls],
#' [compare_clusterings()], [interactive_clustering()], [ssdtwclust()].
#'
#' @example man-examples/tsclust-examples.R
#'
tsclust <- function(series = NULL, type = "partitional", k = 2L, ...,
                    preproc = NULL, distance = "dtw_basic",
                    centroid = ifelse(type == "fuzzy", "fcm", "pam"),
                    control = do.call(paste0(type, "_control"), list()),
                    args = tsclust_args(),
                    seed = NULL, trace = FALSE, error.check = TRUE)
{
    # ==============================================================================================
    # Start
    # ==============================================================================================

    tic <- proc.time()
    handle_rngkind() # UTILS-rng.R
    if (!is.null(seed)) {
        if (length(seed) == 1L)
            set.seed(seed)
        else if (length(seed) == 7L)
            assign(".Random.seed", seed, .GlobalEnv)
        else
            stop("Invalid seed provided") # nocov
    }
    type <- match.arg(type, c("partitional", "hierarchical", "tadpole", "fuzzy"))
    series <- tslist(series, error.check) # coerce to list if necessary
    if (any(k < 2L)) stop("At least two clusters must be defined") # nocov start
    if (any(k >= length(series))) stop("Cannot have more clusters than series in the dataset")
    if (!is.list(control)) stop("Invalid control argument") # nocov end
    MYCALL <- match.call(expand.dots = TRUE)
    dots <- list(...)
    args <- adjust_args(args, dots) # UTILS-utils.R

    # ----------------------------------------------------------------------------------------------
    # Preprocess
    # ----------------------------------------------------------------------------------------------

    if (!is.null(preproc) && is.function(preproc)) {
        series <- quoted_call(preproc, series, dots = subset_dots(args$preproc, preproc))
        preproc_char <- as.character(substitute(preproc))[1L]
    }
    else if (type == "partitional" && is.character(centroid) && centroid == "shape") {
        preproc <- zscore
        preproc_char <- "zscore"
        series <- quoted_call(zscore, series, dots = args$preproc)
    }
    else if (is.null(preproc)) {
        preproc <- function(x, ...) { x } # nocov
        environment(preproc) <- .GlobalEnv
        preproc_char <- "none"
    }
    else stop("Invalid preprocessing")

    if (error.check) check_consistency(series, "vltslist")

    # ----------------------------------------------------------------------------------------------
    # Further options
    # ----------------------------------------------------------------------------------------------

    # after preprocessing!
    distance_missing <- missing(distance)
    diff_lengths <- different_lengths(series)
    check_consistency(distance, "dist", trace = trace, diff_lengths = diff_lengths, silent = FALSE)
    distance <- tolower(distance)
    cent_missing <- missing(centroid)
    cent_char <- check_consistency(centroid, "cent", clus_type = type,
                                   diff_lengths = diff_lengths, cent_missing = cent_missing,
                                   cent_char = as.character(substitute(centroid))[1L])

    if (type != "tadpole") {
        # symmetric versions of dtw that I know of
        # unconstrained and with symmetric1/symmetric2 is always symmetric, regardless of lengths
        # constrained and same lengths with symmetric1/symmetric2 is also symmetric
        symmetric_pattern <- is.null(args$dist$step.pattern) ||
            identical(args$dist$step.pattern, dtw::symmetric1) ||
            identical(args$dist$step.pattern, dtw::symmetric2)

        if (distance %in% c("dtw", "dtw2", "dtw_basic"))
            control$symmetric <- symmetric_pattern && (is.null(args$dist$window.size) || !diff_lengths)
        else if (distance %in% c("lbk", "lbi"))
            control$symmetric <- FALSE
        else if (distance %in% c("sbd", "gak", "sdtw"))
            control$symmetric <- TRUE

        if (distance == "dtw_lb" && isTRUE(args$dist$nn.margin != 1L)) {
            warning("Using dtw_lb in tsclust() always uses row-wise nearest neighbors.")
            args$dist$nn.margin <- 1L
        }
    }

    RET <- switch(
        type,
        partitional =, fuzzy = {
            # ======================================================================================
            # Partitional or fuzzy
            # ======================================================================================

            if (!inherits(control, "PtCtrl") && !inherits(control, "FzCtrl"))
                stop("Invalid control provided") # nocov
            nrep <- if (is.null(control$nrep)) 1L else control$nrep
            if (!is.character(centroid) || !(cent_char %in% c("pam", "fcmdd")))
                control$distmat <- NULL

            # --------------------------------------------------------------------------------------
            # Family creation, see initialization in S4-tsclustFamily.R
            # --------------------------------------------------------------------------------------

            family <- methods::new("tsclustFamily",
                                   dist = distance,
                                   allcent = centroid,
                                   preproc = preproc,
                                   control = control,
                                   fuzzy = isTRUE(type == "fuzzy"))

            if (!all(c("x", "cl_id", "k", "cent", "cl_old") %in% names(formals(family@allcent))))
                stop("The provided centroid function must have at least the following ",
                     "arguments with the shown names:\n\t",
                     paste(c("x", "cl_id", "k", "cent", "cl_old"), collapse = ", "))

            # --------------------------------------------------------------------------------------
            # PAM precompute?
            # --------------------------------------------------------------------------------------

            # precompute distance matrix?
            if (cent_char %in% c("pam", "fcmdd")) {
                dm <- pam_distmat(series, control, distance, cent_char, family, args, trace)
                distmat <- dm$distmat
                distmat_provided <- dm$distmat_provided
                # Redefine new distmat if appropriate
                control$distmat <- distmat
                environment(family@allcent)$control$distmat <- distmat
                if (!(distance == "dtw_lb" && !isTRUE(control$pam.precompute)))
                    environment(family@dist)$control$distmat <- distmat
            }
            else {
                distmat <- NULL
                distmat_provided <- FALSE
            }

            # --------------------------------------------------------------------------------------
            # Cluster
            # --------------------------------------------------------------------------------------

            .rng_ <- dots$.rng_
            if (length(k) == 1L && nrep == 1L) {
                # UTILS-rng.R
                if (is.null(.rng_)) .rng_ <- rng_seq(1L, seed = seed, simplify = TRUE)
                assign(".Random.seed", .rng_, .GlobalEnv)
                # just one repetition,
                # done like this so dist/cent functions can take advantage of parallelization
                pc_list <- list(quoted_call(
                    pfclust,
                    x = series,
                    k = k,
                    family = family,
                    control = control,
                    fuzzy = isTRUE(type == "fuzzy"),
                    cent = cent_char,
                    trace = trace,
                    args = args
                ))
            }
            else {
                # I need to re-register any custom distances in each parallel worker
                dist_entry <- proxy::pr_DB$get_entry(distance)
                export <- c("pfclust", "check_consistency", "quoted_call")
                if (is.null(.rng_))
                    .rng_ <- rng_seq(length(k) * nrep, seed = seed, simplify = FALSE) # UTILS-rng.R
                # if %do% is used, the outer loop replaces values in this envir
                rng0 <- lapply(parallel::splitIndices(length(.rng_), length(k)),
                               function(i) { .rng_[i] })
                k0 <- k
                # sequential allows the sparse matrix to be updated iteratively
                if (inherits(control$distmat, "SparseDistmat")) {
                    `%this_op%` <- `%do%`
                    packages <- setdiff(control$packages, "dtwclust")
                }
                else {
                    `%this_op%` <- `%op%`
                    packages <- control$packages
                }
                # CHECK complains about non-initialization and globals
                rng <- list()
                i <- integer()
                # cluster
                pc_list <- foreach(k = k0, rng = rng0,
                                   .combine = c, .multicombine = TRUE,
                                   .packages = packages,
                                   .export = export) %:%
                    foreach(i = 1L:nrep,
                            .combine = c, .multicombine = TRUE,
                            .packages = packages,
                            .export = export) %this_op%
                            {
                                if (trace) message("Repetition ", i, " for k = ", k)
                                assign(".Random.seed", rng[[i]], .GlobalEnv)

                                if (!check_consistency(dist_entry$names[1L], "dist"))
                                    do.call(proxy::pr_DB$set_entry, dist_entry, TRUE) # nocov

                                # return
                                list(quoted_call(
                                    pfclust,
                                    x = series,
                                    k = k,
                                    family = family,
                                    control = control,
                                    fuzzy = isTRUE(type == "fuzzy"),
                                    cent = cent_char,
                                    trace = trace,
                                    args = args
                                ))
                            }
            }

            # --------------------------------------------------------------------------------------
            # Prepare results
            # --------------------------------------------------------------------------------------

            # Replace distmat with NULL so that, if the distance function is called again,
            # it won't subset it
            environment(family@dist)$control$distmat <- NULL
            # If distmat was provided, let it be shown in the results
            if (distmat_provided) {
                dist_method <- attr(distmat, "method")
                distance <- if (is.null(dist_method)) "unknown" else dist_method
            }
            if (inherits(distmat, "Distmat")) distmat <- distmat$distmat
            if (inherits(distmat, "uninitializedField")) distmat <- NULL
            # Create objects
            RET <- lapply(pc_list, function(pc) {
                if (type == "partitional") {
                    methods::new("PartitionalTSClusters",
                                 call = MYCALL,
                                 family = family,
                                 control = control,
                                 datalist = series,

                                 type = type,
                                 distance = distance,
                                 centroid = cent_char,
                                 preproc = preproc_char,

                                 k = pc$k,
                                 cluster = pc$cluster,
                                 centroids = pc$centroids,
                                 distmat = distmat,

                                 dots = dots,
                                 args = args,

                                 iter = pc$iter,
                                 converged = pc$converged,
                                 clusinfo = pc$clusinfo,
                                 cldist = pc$cldist,

                                 override.family = FALSE)
                }
                else {
                    methods::new("FuzzyTSClusters",
                                 call = MYCALL,
                                 family = family,
                                 control = control,
                                 datalist = series,

                                 type = type,
                                 distance = distance,
                                 centroid = cent_char,
                                 preproc = preproc_char,

                                 k = pc$k,
                                 cluster = pc$cluster,
                                 centroids = pc$centroids,
                                 distmat = distmat,

                                 dots = dots,
                                 args = args,

                                 iter = pc$iter,
                                 converged = pc$converged,
                                 fcluster = pc$fcluster,

                                 override.family = FALSE)
                }
            })
            if (any(!sapply(RET, methods::slot, name = "converged")))
                warning("At least one clustering did not converge within the allowed iterations.")
            # return partitional/fuzzy
            RET
        },
        hierarchical = {
            # ======================================================================================
            # Hierarchical
            # ======================================================================================

            if (!inherits(control, "HcCtrl")) stop("Invalid control provided") # nocov
            method <- control$method
            distmat <- control$distmat
            if (!is.function(centroid)) centroid <- NA
            if (distance == "dtw_lb")
                warning("Using dtw_lb with hierarchical clustering is not advised.") # nocov

            # --------------------------------------------------------------------------------------
            # Calculate distance matrix
            # --------------------------------------------------------------------------------------

            # Take advantage of the function I defined for the partitional methods
            # Which can do calculations in parallel if appropriate
            distfun <- ddist2(distance = distance, control = control)

            if (!is.null(distmat)) {
                if (nrow(distmat) != length(series) || ncol(distmat) != length(series))
                    stop("Dimensions of provided cross-distance matrix don't correspond to ",
                         "length of provided data")
                if (trace) cat("\nDistance matrix provided...\n")

                if (is.null(attr(distmat, "method")))
                    stop("Provided distance matrix does not include the 'method' attribute")
                else
                    distance <- attr(distmat, "method")
            }
            else {
                if (trace) cat("\nCalculating distance matrix...\n")
                distmat <- quoted_call(distfun, x = series, centroids = NULL, dots = args$dist)
            }

            # --------------------------------------------------------------------------------------
            # Cluster
            # --------------------------------------------------------------------------------------

            if (trace) cat("Performing hierarchical clustering...\n")
            if (!base::isSymmetric(base::as.matrix(distmat)))
                warning("Distance matrix is not symmetric, ",
                        "and hierarchical clustering assumes it is ",
                        "(it ignores the upper triangular).")
            if (is.character(method)) {
                # Using hclust
                hc <- lapply(method, function(method) {
                    stats::hclust(stats::as.dist(distmat), method, members = dots$members)
                })
            }
            else {
                # Using provided function
                hc <- list(quoted_call(method, stats::as.dist(distmat), dots = subset_dots(dots, method)))
                method <- attr(method, "name")
            }

            # --------------------------------------------------------------------------------------
            # Prepare results
            # --------------------------------------------------------------------------------------

            if (trace) cat("Extracting centroids...\n\n")
            RET <- lapply(k, function(k) {
                lapply(hc, function(hc) {
                    # cutree and corresponding centroids
                    cluster <- stats::cutree(stats::as.hclust(hc), k)
                    if (is.function(centroid)) {
                        allcent <- function(...) { list(centroid(...)) }
                        environment(allcent) <- new.env(parent = .GlobalEnv)
                        assign("centroid", centroid, environment(allcent))
                        centroids <- lapply(1L:k, function(kcent) {
                            quoted_call(centroid,
                                        series[cluster == kcent],
                                        dots = subset_dots(args$cent, centroid))
                        })
                    }
                    else {
                        allcent <- function(...) {} # dummy
                        centroid_ids <- sapply(1L:k, function(kcent) {
                            id_k <- cluster == kcent
                            d_sub <- distmat[id_k, id_k, drop = FALSE]
                            id_centroid <- which.min(apply(d_sub, 1L, sum))
                            which(id_k)[id_centroid]
                        })
                        centroids <- series[centroid_ids]
                        attr(centroids, "series_id") <- unname(centroid_ids)
                    }

                    methods::new("HierarchicalTSClusters",
                                 stats::as.hclust(hc),
                                 call = MYCALL,
                                 family = methods::new("tsclustFamily",
                                                       dist = distfun,
                                                       allcent = allcent,
                                                       preproc = preproc),
                                 control = control,
                                 datalist = series,

                                 type = type,
                                 distance = distance,
                                 centroid = cent_char,
                                 preproc = preproc_char,

                                 k = as.integer(k),
                                 cluster = cluster,
                                 centroids = centroids,
                                 distmat = distmat,

                                 dots = dots,
                                 args = args,

                                 method = if (!is.null(hc$method)) hc$method else method,

                                 override.family = !is.function(centroid))
                })
            })
            RET <- unlist(RET, recursive = FALSE)
            # return hierarchical
            RET
        },
        tadpole = {
            # ======================================================================================
            # TADPole
            # ======================================================================================

            if (!inherits(control, "TpCtrl")) stop("Invalid control provided") # nocov start
            if (!distance_missing)
                warning("The distance argument is ignored for TADPole.") # nocov end

            # --------------------------------------------------------------------------------------
            # Parameters
            # --------------------------------------------------------------------------------------

            # for predict and cvi
            distfun <- ddist2("dtw_basic", control = control)
            # for family@dist
            args$dist$window.size <- control$window.size
            args$dist$norm <- "L2"
            args$dist$window.type <- "sakoechiba"

            # --------------------------------------------------------------------------------------
            # Cluster
            # --------------------------------------------------------------------------------------

            if (trace) cat("\n\tEntering TADPole...\n\n")
            R <- TADPole(series,
                         k = k,
                         dc = control$dc,
                         window.size = control$window.size,
                         lb = control$lb,
                         trace = trace)
            if (length(k) == 1L && length(control$dc) == 1L) R <- list(R)

            # --------------------------------------------------------------------------------------
            # Prepare results
            # --------------------------------------------------------------------------------------

            # seeds (UTILS-rng.R)
            .rng_ <- dots$.rng_
            if (is.null(.rng_))
                .rng_ <- rng_seq(length(k) * length(control$dc), seed = seed, simplify = FALSE)

            RET <- Map(R, .rng_, f = function(R, rng) {
                assign(".Random.seed", rng, .GlobalEnv)
                k <- length(R$centroids)
                if (is.function(centroid)) {
                    allcent <- function(...) { list(centroid(...)) }
                    environment(allcent) <- new.env(parent = .GlobalEnv)
                    assign("centroid", centroid, environment(allcent))
                    centroids <- lapply(1L:k, function(kcent) {
                        quoted_call(centroid,
                                    series[R$cl == kcent],
                                    dots = subset_dots(args$cent, centroid))
                    })
                }
                else {
                    allcent <- function(...) {}
                    centroids <- series[R$centroids]
                    attr(centroids, "series_id") <- R$centroids
                }

                obj <- methods::new("PartitionalTSClusters",
                                    call = MYCALL,
                                    family = methods::new("tsclustFamily",
                                                          dist = distfun,
                                                          allcent = allcent,
                                                          preproc = preproc),
                                    control = control,
                                    datalist = series,

                                    type = type,
                                    distance = "dtw_lb",
                                    centroid = cent_char,
                                    preproc = preproc_char,

                                    k = as.integer(k),
                                    cluster = R$cl,
                                    centroids = centroids,
                                    distmat = NULL,

                                    dots = dots,
                                    args = args,

                                    override.family = !is.function(centroid))
                obj@distance <- "LB+DTW2"
                obj
            })
            # return tadpole
            RET
        }
    )

    # ==============================================================================================
    # Finish
    # ==============================================================================================

    toc <- proc.time() - tic
    RET <- lapply(RET, function(ret) {
        ret@proctime <- toc
        ret@seed <- as.integer(seed)
        ret
    })
    if (length(RET) == 1L)
        RET <- RET[[1L]]
    else if (type %in% c("partitional", "fuzzy"))
        attr(RET, "rng") <- unlist(rng0, recursive = FALSE, use.names = FALSE)
    if (trace) cat("\tElapsed time is", toc["elapsed"], "seconds.\n\n")
    RET
}
########

# handle_rngkind() でエラー: 関数 "handle_rngkind" を見つけることができませんでした 
# https://github.com/asardaes/dtwclust/blob/f4a5978ce051585d3c989e658c84d6b4b0fb77d4/R/UTILS-rng.R
########
#' @importFrom parallel nextRNGStream
#'
rng_seq <- function(n, seed = NULL, simplify = FALSE) {
    stopifnot(RNGkind()[1L] == dtwclust_rngkind, n > 0L)

    if (!is.null(seed)) {
        prev_seed <- get(".Random.seed", .GlobalEnv)
        if (length(seed) == 1L) {
            set.seed(as.integer(seed))
        }
        else if (length(seed) == 7L) {
            assign(".Random.seed", seed, .GlobalEnv)
        }
        else {
            stop("Invalid seed provided") # nocov
        }
        on.exit(assign(".Random.seed", prev_seed, .GlobalEnv))
    }

    first_seed <- get(".Random.seed", .GlobalEnv)
    if (n == 1L) {
        if (!isTRUE(simplify)) first_seed <- list(first_seed)
        return(first_seed) # return 1
    }

    seed_seq <- vector("list", n)
    seed_seq[[1L]] <- first_seed
    for (i in 2L:n) seed_seq[[i]] <- parallel::nextRNGStream(seed_seq[[i-1L]])
    # return
    seed_seq
}

# see https://stackoverflow.com/a/20998531/5793905 and the link there
handle_rngkind <- function() {
    previous_rngkind <- RNGkind(dtwclust_rngkind)[1L]
    if (previous_rngkind != dtwclust_rngkind) {
        # evaluate on.exit on the caller's environment, otherwise it would execute immediately
        do.call(on.exit, envir = parent.frame(), args = list(substitute({
            RNGkind(previous_rngkind)
        })))
    }
    invisible()
}
########

# tslist(series, error.check) でエラー:関数 "tslist" を見つけることができませんでした
# https://github.com/asardaes/dtwclust/blob/f4a5978ce051585d3c989e658c84d6b4b0fb77d4/R/UTILS-tslist.R
########
#' Coerce matrices or data frames to a list of time series
#'
#' Change a matrix or data frame to a list of univariate time series
#'
#' @export
#'
#' @param series A matrix or data frame where each row is a time series.
#' @param simplify Coerce all series in the resulting list to either matrix (multivariate) or
#'   numeric (univariate).
#'
#' @details
#'
#' Almost all functions in \pkg{dtwclust} work internally with lists of time series. If you want to
#' avoid constant coercion, create a list of time series once by calling this function.
#'
#' For matrices and data frames, each **row** is considered as one time series. A list input is
#' simply passed through.
#'
#' @return
#'
#' A list of time series.
#'
#' @note
#'
#' The function assumes that matrix-like objects can be first coerced via [base::as.matrix()], so
#' that the result can be indexed with `series[i, ]`.
#'
#' No consistency checks are performed by this function.
#'
tslist <- function(series, simplify = FALSE) {
    if (is.matrix(series) || is.data.frame(series)) {
        rnms <- rownames(series)
        mat <- unname(base::as.matrix(series))
        series <- vector("list", nrow(mat))
        if (!is.null(rnms)) setnames_inplace(series, rnms)
        for (i in 1L:nrow(mat)) {
            series[[i]] <- mat[i,]
        }
    }
    else if (is.numeric(series))
        series <- list(series) # nocov
    else if (!is.list(series))
        stop("Unsupported data type.")
    # coerce to simple types that are known to work
    if (simplify) {
        if (is_multivariate(series))
            series <- lapply(series, base::as.matrix)
        else
            series <- lapply(series, base::as.numeric)
    }
    # return
    series
}
########

# is_multivariate(series) でエラー:関数 "is_multivariate" を見つけることができませんでした
# https://github.com/asardaes/dtwclust/blob/bc6008c26b46ff7105b59135dd9ad004a99f7d87/R/UTILS-utils.R
########
# ==================================================================================================
# Miscellaneous
# ==================================================================================================

check_consistency <- function(obj, case, ..., clus_type,
                              diff_lengths = FALSE, cent_missing, cent_char,
                              trace = FALSE, silent = TRUE)
{
    case <- match.arg(case, c("ts", "tslist", "vltslist", "window", "dist", "cent"))
    if (case == "ts") {
        if (!is.numeric(obj)) stop("The series must be numeric")
        if (length(obj) < 1L) stop("The series must have at least one point")
        if (anyNA(obj)) stop("There are missing values in the series")
    }
    else if (case %in% c("tslist", "vltslist")) {
        if (!is.list(obj)) stop("Oops, data should already be a list by this point...") # nocov
        if (length(obj) < 1L) stop("Data is empty")
        if (case == "tslist" && different_lengths(obj)) stop("All series must have the same length")
        sapply(obj, check_consistency, case = "ts", ...)
    }
    else if (case == "window") {
        if (is.null(obj)) stop("Please provide the 'window.size' parameter")
        if (any(obj < 0L)) stop("Window size must be non-negative")
        return(as.integer(obj))
    }
    else if (case == "dist") {
        if (!is.character(obj) || !pr_DB$entry_exists(obj)) {
            if (silent)
                return(FALSE)
            else
                stop("Please provide the name of a valid distance function registered with the ",
                     "'proxy' package.")
        }
        if (diff_lengths) {
            obj <- tolower(obj)
            if ((obj %in% distances_known) && !(obj %in% distances_difflength))
                stop("Only the following distances are supported for series with different length:\n\t",
                     paste(distances_difflength, collapse = "\t"))
            else if (!(obj %in% distances_known) && trace)
                message("Series have different lengths. ", # nocov start
                        "Please confirm that the provided distance function supports this.") # nocov end
        }
        # valid registered distance
        return(TRUE)
    }
    else if (case == "cent") {
        cent_char <- switch(
            clus_type,
            partitional = {
                if (is.character(obj)) {
                    cent_char <- match.arg(obj, centroids_nonfuzzy)
                    if (diff_lengths &&
                        cent_char %in% centroids_included &&
                        !(cent_char %in% centroids_difflength))
                        stop("Only the following centroids are supported for ",
                             "series with different lengths:\n\t",
                             paste(centroids_difflength, collapse = "\t"))
                }
                else {
                    force(cent_char)
                }
                # return partitional switch
                cent_char
            },
            fuzzy = {
                if (is.character(obj)) {
                    cent_char <- match.arg(obj, centroids_fuzzy)
                    if (diff_lengths && cent_char == "fcm")
                        stop("Fuzzy c-means does not support series with different length.")
                }
                else {
                    force(cent_char)
                }
                # return fuzzy switch
                cent_char
            },
            hierarchical =, tadpole = {
                if (is.function(obj))
                    force(cent_char)
                else if (!cent_missing)
                    warning("The 'centroid' argument was provided but it wasn't a function, ",
                            "so it was ignored.",
                            call. = FALSE, immediate. = TRUE)
                else
                    cent_char <- paste("PAM",
                                       switch(clus_type,
                                              hierarchical = "(Hierarchical)",
                                              tadpole = "(TADPole)"))
                # return hierarchical/tadpole switch
                cent_char
            }
        )
        # cent case
        return(cent_char)
    }
    invisible(NULL)
}

# Check if list of series have different length
different_lengths <- function(x) { any(diff(lengths(x)) != 0L) }

# Enlist parameters for do.calls
enlist <- function(..., dots = NULL) { c(list(...), dots) }

# Check if a function has the ellipsis in its formals
has_dots <- function(foo) { is.function(foo) && !is.null(formals(foo)$`...`) }

# Subset dots for do.calls of functions without ellipsis
subset_dots <- function(dots = list(), foo) {
    if (has_dots(foo))
        dots
    else if (length(dots) > 0L)
        dots[intersect(names(dots), names(formals(foo)))]
    else
        list()
}

# Adjust args for tsclust() and TSClusters-class
adjust_args <- function(args, dots) {
    lapply(args, function(arg) {
        arg <- c(arg, dots)
        arg[!duplicated(names(arg))]
    })
}

# Like dynGet() I assume, but that one is supposed to be experimental...
get_from_callers <- function(obj_name, mode = "any") {
    ret <- get0(obj_name, mode = mode, inherits = TRUE)
    if (!is.null(ret)) return(ret)
    for (env in sys.frames()) {
        ret <- get0(obj_name, env, mode = mode, inherits = FALSE)
        if (!is.null(ret)) return(ret)
    }
    stop("Could not find object '", obj_name, "' of mode '", mode, "'") # nocov
}

# do.call but always quoted
quoted_call <- function(fun, ..., dots = NULL) {
    do.call(fun, enlist(..., dots = dots), quote = TRUE)
}

# ==================================================================================================
# Helper C/C++ functions
# ==================================================================================================

# Create combinations of all possible pairs
call_pairs <- function(n = 2L) {
    n <- as.integer(n)
    if (n < 2L) stop("At least two elements are needed to create pairs.") # nocov
    .Call(C_pairs, n, PACKAGE = "dtwclust")
}

# Modify names in place
setnames_inplace <- function(vec, names) {
    if (!is.vector(vec) || !is.vector(names)) stop("Both 'vec' and 'names' must be vectors.")
    if (length(vec) != length(names)) stop("Length mismatch when changing names in place.")
    if (!is.character(names)) stop("Trying to set names in place with non-character names.") # nocov
    invisible(.Call(C_setnames_inplace, vec, names, PACKAGE = "dtwclust"))
}

# ==================================================================================================
# Parallel helper functions
# ==================================================================================================

# Custom binary operator for %dopar% to avoid unnecessary warnings and adjust available threads
`%op%` <- function(obj, ex) {
    # check to see if the workers have specified how many threads to use
    num_workers <- foreach::getDoParWorkers()
    backend_name <- foreach::getDoParName()
    if (num_workers > 1L) {
        reset_workers <- foreach::foreach(
            i = 1L:num_workers,
            .combine = c,
            .multicombine = TRUE,
            .inorder = FALSE,
            .packages = "dtwclust"
        ) %dopar% {
            reset <- TRUE
            if (nzchar(Sys.getenv("RCPP_PARALLEL_NUM_THREADS")))
                reset <- FALSE # nocov
            else
                RcppParallel::setThreadOptions(1L)
            reset
        }
    }
    # check the RNGkind of the workers
    if (backend_name != "doSEQ") {
        rng_kind <- foreach::foreach(
            i = 1L:num_workers,
            .combine = c,
            .multicombine = TRUE,
            .inorder = FALSE
        ) %dopar% {
            RNGkind("L'Ecuyer-CMRG")[1L]
        }
    }
    # do not load dtwclust in sequential cases
    if (backend_name == "doSEQ" || (backend_name == "doFuture" && num_workers == 1L)) {
        obj$packages <- setdiff(obj$packages, "dtwclust")
    }
    # evaluate expression
    withCallingHandlers({
        ret <- tryCatch(eval.parent(substitute(obj %dopar% ex)), error = function(e) { e })
    },
    warning = function(w) {
        if (grepl("package:dtwclust", w$message, ignore.case = TRUE))
            invokeRestart("muffleWarning") # nocov
    })
    # reset parallel workers if needed
    if (num_workers > 1L && any(reset_workers)) {
        foreach::foreach(
            i = 1L:num_workers,
            .combine = c,
            .multicombine = TRUE,
            .inorder = FALSE,
            .packages = "dtwclust"
        ) %dopar% {
            RcppParallel::setThreadOptions("auto")
        }
    }
    # reset RNGkind if needed
    if (backend_name != "doSEQ" && any(rng_kind != dtwclust_rngkind)) {
        foreach::foreach(
            rng_kind = rng_kind,
            .combine = c,
            .multicombine = TRUE,
            .inorder = FALSE
        ) %dopar% {
            RNGkind(rng_kind)
        }
    }
    if (inherits(ret, "error") && obj$errorHandling != "pass") stop(ret)
    # return
    ret
}

# Split a given object into chunks for parallel workers
#' @importFrom parallel splitIndices
#'
split_parallel <- function(obj, margin = NULL) {
    num_workers <- foreach::getDoParWorkers()
    if (num_workers == 1L) return(list(obj))

    num_tasks <- if (is.null(margin)) length(obj) else dim(obj)[margin]
    if (!is.integer(num_tasks))
        stop("Invalid attempt to split an object into parallel tasks") # nocov

    num_tasks <- parallel::splitIndices(num_tasks, num_workers)
    num_tasks <- num_tasks[lengths(num_tasks, use.names = FALSE) > 0L]

    if (is.null(margin))
        ret <- lapply(num_tasks, function(id) obj[id])
    else
        ret <- switch(EXPR = margin,
                      lapply(num_tasks, function(id) obj[id, , drop = FALSE]),
                      lapply(num_tasks, function(id) obj[ , id, drop = FALSE]))
    # return
    ret
}

# This only works if it's used after split_parallel()
validate_pairwise <- function(x, y) {
    if (!identical(lengths(x, use.names = FALSE), lengths(y, use.names = FALSE)))
        stop("Pairwise distances require the same amount of series in 'x' and 'y'.") # nocov

    invisible(NULL)
}

# Number of configured/available threads according to RcppParallel
#' @importFrom RcppParallel defaultNumThreads
#'
get_nthreads <- function() {
    as.integer(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", RcppParallel::defaultNumThreads()))
}

# ==================================================================================================
# Helper distance-related
# ==================================================================================================

# allocate distance matrix for custom proxy loops
allocate_distmat <- function(x_len, y_len, pairwise, symmetric) {
    if (pairwise)
        D <- matrix(0, x_len, 1L)
    else if (symmetric)
        D <- matrix(0, x_len, x_len)
    else
        D <- matrix(0, x_len, y_len)
    # return
    D
}

# Euclidean norm
l2norm <- function(x) { sqrt(sum(x * x)) }

# Get the clusinfo slot for TSClusters
compute_clusinfo <- function(k, cluster, cldist) {
    cluster <- factor(cluster, levels = 1L:k)
    # return
    data.frame(
        size = as.integer(table(cluster)),
        av_dist = as.numeric(tapply(cldist[, 1L], list(cluster), mean, default = 0))
    )
}

# ==================================================================================================
# Multivariate helpers
# ==================================================================================================

is_multivariate <- function(x) {
    if (length(x) == 0L) stop("Empty list of series received.") # nocov
    ncols <- sapply(x, NCOL)
    if (any(diff(ncols) != 0L)) stop("Inconsistent dimensions across series.")
    any(ncols > 1L)
}

reshape_multivariate <- function(series, cent) {
    ncols <- ncol(series[[1L]])
    series <- lapply(1L:ncols, function(idc) {
        lapply(series, function(s) { s[ , idc, drop = TRUE] })
    })
    cent <- lapply(1L:ncols, function(idc) {
        if (is.null(cent))
            NULL
        else
            cent[ , idc, drop = TRUE]
    })
    list(series = series, cent = cent)
}
########

# partitional_control() でエラー:関数 "partitional_control" を見つけることができませんでした
# https://github.com/asardaes/dtwclust/blob/bc6008c26b46ff7105b59135dd9ad004a99f7d87/R/CLUSTERING-tsclust-controls.R
########
#' Control parameters for clusterings with [tsclust()]
#'
#' Control parameters for fine-grained control.
#'
#' @name tsclust-controls
#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @param pam.precompute Logical flag. Precompute the whole distance matrix once and reuse it on
#'   each iteration if using PAM centroids. Otherwise calculate distances at every iteration. See
#'   details.
#' @param iter.max Integer. Maximum number of allowed iterations for partitional/fuzzy clustering.
#' @param nrep Integer. How many times to repeat clustering with different starting points (i.e.,
#'   different random seeds).
#' @param symmetric Logical flag. Is the distance function symmetric? In other words, is `dist(x,y)`
#'   == `dist(y,x)`? If `TRUE`, only half the distance matrix needs to be computed. Automatically
#'   detected and overridden for the distances included in \pkg{dtwclust}.
#' @param packages Character vector with the names of any packages required for custom `proxy`
#'   functions. Relevant for parallel computation, although since the distance entries are
#'   re-registered in each parallel worker if needed, this is probably useless, but just in case.
#' @param distmat If available, the cross-distance matrix can be provided here. Only relevant for
#'   partitional with PAM centroids, fuzzy with FCMdd centroids, or hierarchical clustering.
#' @param pam.sparse Attempt to use a sparse matrix for PAM centroids. See details.
#' @param version Which version of partitional/fuzzy clustering to use. See details.
#'
#' @details
#'
#' The functions essentially return their function arguments in a classed list, although some checks
#' are performed.
#'
#' Regarding parameter \code{version}: the first version of partitional/fuzzy clustering implemented
#' in the package always performed an extra iteration, which is unnecessary. Use version 1 to mimic
#' this previous behavior.
#'
#' @section Partitional:
#'
#'   When `pam.precompute = FALSE`, using `pam.sparse = TRUE` defines a sparse matrix (refer to
#'   [Matrix::sparseMatrix()]) and updates it every iteration (except for `"dtw_lb"` distance). For
#'   most cases, precomputing the whole distance matrix is still probably faster. See the timing
#'   experiments in `browseVignettes("dtwclust")`.
#'
#'   Parallel computations for PAM centroids have the following considerations:
#'
#'   - If `pam.precompute` is `TRUE`, both distance matrix calculations and repetitions are done in
#'   parallel, regardless of `pam.sparse`.
#'   - If `pam.precompute` is `FALSE` and `pam.sparse` is `TRUE`, repetitions are done sequentially,
#'   so that the distance calculations can be done in parallel and the sparse matrix updated
#'   iteratively.
#'   - If both `pam.precompute` and `pam.sparse` are `FALSE`, repetitions are done in parallel, and
#'   each repetition performs distance calculations sequentially, but the distance matrix cannot be
#'   updated iteratively.
#'
partitional_control <- function(pam.precompute = TRUE,
                                iter.max = 100L,
                                nrep = 1L,
                                symmetric = FALSE,
                                packages = character(0L),
                                distmat = NULL,
                                pam.sparse = FALSE,
                                version = 2L)
{
    if (any(iter.max <= 0L)) stop("Maximum iterations must be positive") # nocov start
    if (any(nrep < 1L)) stop("Number of repetitions must be at least one") # nocov end
    # return
    structure(
        list(pam.precompute = as.logical(pam.precompute),
             pam.sparse = as.logical(pam.sparse),
             iter.max = as.integer(iter.max),
             nrep = as.integer(nrep)[1L],
             symmetric = as.logical(symmetric)[1L],
             packages = unique(c("dtwclust", as.character(packages))),
             distmat = distmat,
             version = as.integer(version)),
        "class" = c(control_classes[["partitional"]])
    )
}

#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @param method Character vector with one or more linkage methods to use in hierarchical procedures
#'   (see [stats::hclust()]), the character `"all"` to use all of the available ones, or a function
#'   that performs hierarchical clustering based on distance matrices (e.g. [cluster::diana()]). See
#'   details.
#'
#' @section Hierarchical:
#'
#'   There are some limitations when using a custom hierarchical function in `method`: it will
#'   receive the lower triangular of the distance matrix as first argument (see [stats::as.dist()])
#'   and the result should support the [stats::as.hclust()] generic. This functionality was added
#'   with the \pkg{cluster} package in mind, since its functions follow this convention, but other
#'   functions could be used if they are adapted to work similarly.
#'
hierarchical_control <- function(method = "average",
                                 symmetric = FALSE,
                                 packages = character(0L),
                                 distmat = NULL)
{
    if (is.character(method)) {
        method <- match.arg(method,
                            c("ward.D", "ward.D2", "single", "complete",
                              "average", "mcquitty", "median", "centroid",
                              "all"),
                            several.ok = TRUE)
        if ("all" %in% method)
            method <- c("ward.D", "ward.D2", "single", "complete",
                        "average", "mcquitty", "median", "centroid")
    }
    else if (!is.function(method))
        stop("Argument 'method' must be either a supported character or a function.") # nocov
    else
        attr(method, "name") <- as.character(substitute(method))[1L]
    # return
    structure(
        list(method = method,
             symmetric = as.logical(symmetric),
             packages = unique(c("dtwclust", as.character(packages))),
             distmat = distmat),
        "class" = c(control_classes[["hierarchical"]])
    )
}

#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @param fuzziness Numeric. Exponent used for fuzzy clustering. Commonly termed `m` in the
#'   literature.
#' @param delta Numeric. Convergence criterion for fuzzy clustering.
#'
fuzzy_control <- function(fuzziness = 2,
                          iter.max = 100L,
                          delta = 1e-3,
                          packages = character(0L),
                          symmetric = FALSE,
                          version = 2L,
                          distmat = NULL)
{
    if (any(fuzziness <= 1)) stop("Fuzziness exponent should be greater than one") # nocov start
    if (any(iter.max <= 0L)) stop("Maximum iterations must be positive")
    if (any(delta < 0)) stop("Delta should be positive") # nocov end
    # return
    structure(
        list(fuzziness = fuzziness,
             iter.max = as.integer(iter.max),
             delta = delta,
             symmetric = as.logical(symmetric)[1L],
             packages = unique(c("dtwclust", as.character(packages))),
             version = as.integer(version),
             distmat = distmat),
        "class" = c(control_classes[["fuzzy"]])
    )
}

#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @param dc The cutoff distance for the TADPole algorithm.
#' @param window.size The window.size specifically for the TADPole algorithm.
#' @param lb The lower bound to use with TADPole. Either `"lbk"` or `"lbi"`.
#'
#' @section TADPole:
#'
#'   When using TADPole, the `dist` argument list includes the `window.size` and specifies `norm =
#'   "L2"`.
#'
tadpole_control <- function(dc,
                            window.size,
                            lb = "lbk")
{
    if (any(dc <= 0)) stop("Cutoff distance 'dc' must be positive") # nocov
    window.size <- check_consistency(window.size, "window")
    lb <- match.arg(lb, c("lbk", "lbi"), several.ok = TRUE)
    # return
    structure(
        list(dc = dc,
             window.size = window.size,
             lb = lb),
        "class" = c(control_classes[["tadpole"]])
    )
}

#' @rdname tsclust-controls
#' @aliases tsclust-controls
#' @export
#'
#' @param preproc A list of arguments for a preprocessing function to be used in [tsclust()].
#' @param dist A list of arguments for a distance function to be used in [tsclust()].
#' @param cent A list of arguments for a centroid function to be used in [tsclust()].
#'
tsclust_args <- function(preproc = list(), dist = list(), cent = list()) {
    structure(
        list(preproc = preproc, dist = dist, cent = cent),
        "class" = c(control_classes[["args"]])
    )
}
########

# structure(list(pam.precompute = as.logical(pam.precompute), pam.sparse = as.logical(pam.sparse),  でエラー:オブジェクト 'control_classes' がありません
# https://github.com/asardaes/dtwclust/blob/bc6008c26b46ff7105b59135dd9ad004a99f7d87/R/UTILS-globals-internal.R
########
# ==================================================================================================
# Internal global variables
# ==================================================================================================

supported_clusterings <- c("partitional", "hierarchical", "fuzzy", "tadpole")
dtwclust_rngkind <- "L'Ecuyer-CMRG"

distances_known <- c("dtw", "dtw2", "dtw_lb", "lbk", "lbi", "sbd", "dtw_basic", "gak", "sdtw")
distances_included <- c("dtw_lb", "lb_keogh", "lb_improved", "sbd", "dtw_basic", "gak", "sdtw")
distances_difflength <- c("dtw", "dtw2", "sbd", "dtw_basic", "gak", "sdtw")
distances_multivariate <- c("dtw", "dtw2", "dtw_basic", "gak", "sdtw")

centroids_included <- c("mean", "median", "shape", "dba", "pam", "fcm", "fcmdd", "sdtw_cent")
centroids_fuzzy <- c("fcm", "fcmdd")
centroids_nonfuzzy <- setdiff(centroids_included, centroids_fuzzy)
centroids_difflength <- c("dba", "pam", "shape", "fcmdd", "sdtw_cent")

control_classes <- c(partitional = "PtCtrl",
                     hierarchical = "HcCtrl",
                     fuzzy = "FzCtrl",
                     tadpole = "TpCtrl",
                     args = "TscArgs")
########

# check_consistency(distance, "dist", trace = trace, diff_lengths = diff_lengths,  でエラー:オブジェクト 'pr_DB' がありません
# https://github.com/asardaes/dtwclust/blob/bc6008c26b46ff7105b59135dd9ad004a99f7d87/R/pkg.R
########
#' Time series clustering along with optimizations for the Dynamic Time Warping distance
#'
#' Time series clustering with a wide variety of strategies and a series of optimizations specific
#' to the Dynamic Time Warping (DTW) distance and its corresponding lower bounds (LBs).
#'
#' @name dtwclust-package
#'
#' @details
#'
#' Many of the algorithms implemented in this package are specifically tailored to DTW, hence its
#' name. However, the main clustering function is flexible so that one can test many different
#' clustering approaches, using either the time series directly, or by applying suitable
#' transformations and then clustering in the resulting space. Other implementations included in the
#' package provide some alternatives to DTW.
#'
#' DTW is a dynamic programming algorithm that tries to find the optimum warping path between two
#' series. Over the years, several variations have appeared in order to make the procedure faster or
#' more efficient. Please refer to the included references for more information, especially Giorgino
#' (2009), which is a good practical introduction.
#'
#' Most optimizations require equal dimensionality, which means time series should have equal
#' length. DTW itself does not require this, but it is relatively expensive to compute. Other
#' distance definitions may be used, or series could be reinterpolated to a matching length
#' (Ratanamahatana and Keogh 2004).
#'
#' The main clustering function and entry point for this package is [tsclust()], with a convenience
#' wrapper for multiple tests in [compare_clusterings()], and a shiny app in
#' [interactive_clustering()]. There is another less-general-purpose shiny app in [ssdtwclust()].
#'
#' Please note the random number generator is set to L'Ecuyer-CMRG when \pkg{dtwclust} is attached
#' in an attempt to preserve reproducibility. You are free to change this afterwards if you wish
#' (see [base::RNGkind()]), but \pkg{dtwclust} will always use L'Ecuyer-CMRG internally.
#'
#' For more information, please read the included package vignettes, which can be accessed by typing
#' `browseVignettes("dtwclust")`.
#'
#' @note
#'
#' This software package was developed independently of any organization or institution that is or
#' has been associated with the author.
#'
#' This package can be used without attaching it with [base::library()] with some caveats:
#'
#' - The \pkg{methods} [package][methods::methods-package] must be attached. `R` usually does this
#'   automatically, but [utils::Rscript()] only does so in R versions 3.5.0 and above.
#' - If you want to use the \pkg{proxy} version of [dtw::dtw()] (e.g. for clustering), you have to
#'   attach the \pkg{dtw} package manually.
#'
#' Be careful with reproducibility, `R`'s random number generator is only changed session-wide if
#' \pkg{dtwclust} is attached.
#'
#' @author Alexis Sarda-Espinosa
#'
#' @references
#'
#' Please refer to the package's vignette's references.
#'
#' @seealso
#'
#' [tsclust()], [compare_clusterings()], [interactive_clustering()], [ssdtwclust()], [dtw_basic()],
#' [proxy::dist()].
#'
#' @useDynLib dtwclust, .registration = TRUE
#' @import foreach
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
"_PACKAGE"

# PREFUN for some of my proxy distances so that they support 'pairwise' directly
proxy_prefun <- function(x, y, pairwise, params, reg_entry) {
    params$pairwise <- pairwise
    list(x = x, y = y, pairwise = pairwise, p = params, reg_entry = reg_entry)
}

.onLoad <- function(lib, pkg) {
    # Register DTW2
    if (!check_consistency("DTW2", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = dtw2_proxy, names=c("DTW2", "dtw2"),
                               loop = TRUE, type = "metric", distance = TRUE,
                               description = "DTW with L2 norm",
                               PACKAGE = "dtwclust")
    # Register DTW_BASIC
    if (!check_consistency("DTW_BASIC", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = dtw_basic_proxy, names=c("DTW_BASIC", "dtw_basic"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Basic and maybe faster DTW distance",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)
    # Register LB_Keogh
    if (!check_consistency("LB_Keogh", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = lb_keogh_proxy, names=c("LBK", "LB_Keogh", "lbk"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Keogh's DTW lower bound for the Sakoe-Chiba band",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)
    # Register LB_Improved
    if (!check_consistency("LB_Improved", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = lb_improved_proxy, names=c("LBI", "LB_Improved", "lbi"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Lemire's improved DTW lower bound for the Sakoe-Chiba band",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)
    # Register SBD
    if (!check_consistency("SBD", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = sbd_proxy, names=c("SBD", "sbd"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Paparrizos and Gravanos' shape-based distance for time series",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun,
                               convert = function(d) { 2 - d })
    # Register DTW_LB
    if (!check_consistency("DTW_LB", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = dtw_lb, names=c("DTW_LB", "dtw_lb"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "DTW distance aided with Lemire's lower bound",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)
    # Register GAK
    if (!check_consistency("GAK", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = gak_proxy, names=c("GAK", "gak"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Fast (triangular) global alignment kernel distance",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun,
                               convert = function(d) { 1 - d })
    # Register uGAK
    if (!check_consistency("uGAK", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = gak_simil, names=c("uGAK", "ugak"),
                               loop = FALSE, type = "metric", distance = FALSE,
                               description = "Fast (triangular) global alignment kernel similarity",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)
    # Register soft-DTW
    if (!check_consistency("sdtw", "dist", silent = TRUE))
        proxy::pr_DB$set_entry(FUN = sdtw_proxy, names=c("sdtw", "SDTW", "soft-DTW"),
                               loop = FALSE, type = "metric", distance = TRUE,
                               description = "Soft-DTW",
                               PACKAGE = "dtwclust", PREFUN = proxy_prefun)

    # avoids default message if no backend exists
    if (is.null(foreach::getDoParName())) foreach::registerDoSEQ()
}

#' @importFrom utils packageVersion
#'
.onAttach <- function(lib, pkg) {
    RNGkind(dtwclust_rngkind)

    packageStartupMessage("dtwclust:\n",
                          "Setting random number generator to L'Ecuyer-CMRG (see RNGkind()).\n",
                          'To read the included vignettes type: browseVignettes("dtwclust").\n',
                          'See news(package = "dtwclust") after package updates.')

    if (grepl("\\.9000$", utils::packageVersion("dtwclust")))
        packageStartupMessage("This is a developer version of dtwclust.")
}

.onUnload <- function(libpath) {
    # Unegister distances
    if (check_consistency("DTW2", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("DTW2")
    if (check_consistency("DTW_BASIC", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("DTW_BASIC")
    if (check_consistency("LB_Keogh", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("LB_Keogh")
    if (check_consistency("LB_Improved", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("LB_Improved")
    if (check_consistency("SBD", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("SBD")
    if (check_consistency("DTW_LB", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("DTW_LB")
    if (check_consistency("GAK", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("GAK")
    if (check_consistency("uGAK", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("uGAK")
    if (check_consistency("sdtw", "dist", silent = TRUE)) proxy::pr_DB$delete_entry("sdtw")
    library.dynam.unload("dtwclust", libpath)
}

release_questions <- function() {
    c(
        "Changed .Rbuildignore to exclude test rds files?",
        "Built the binary with --compact-vignettes?",
        "Set vignette's cache to FALSE?"
    )
}
########

# pr_DB
# https://github.com/asardaes/dtwclust/blob/bc6008c26b46ff7105b59135dd9ad004a99f7d87/R/CLUSTERING-ddist2.R
########
# ==================================================================================================
# Helpers
# ==================================================================================================

# Use existing distmat if available
use_distmat <- function(distmat, x, centroids) {
    if (!inherits(distmat, "Distmat"))
        stop("Invalid distance matrix in control.") # nocov
    # internal class, sparse or full
    i <- 1L:length(x)
    j <- if (is.null(centroids)) i else distmat$id_cent
    distmat[i, j, drop = FALSE]
}

# Extra distance parameters in case of parallel computation
# They can be for the function or for proxy::dist
#' @importFrom proxy dist
#'
get_dots <- function(dist_entry, x, centroids, ...) {
    dots <- list(...)

    # Added defaults
    if (is.null(dots$window.size)) {
        dots$window.type <- "none"
    }
    else if (is.null(dots$window.type)) {
        dots$window.type <- "slantedband"
    }

    dots$error.check <- FALSE

    # dtw uses L2 by default, but in dtwclust I want dtw to use L1 by default
    # Important for multivariate series
    if (tolower(dist_entry$names[1L]) == "dtw" && is.null(dots$dist.method) && is_multivariate(c(x, centroids))) {
        dots$dist.method <- "L1" # nocov
    }

    # If the function doesn't have '...', remove invalid arguments from 'dots'
    valid_args <- names(dots)
    if (is.function(dist_entry$FUN)) {
        if (!has_dots(dist_entry$FUN)) {
            valid_args <- union(names(formals(proxy::dist)), names(formals(dist_entry$FUN)))
        }
    }
    else {
        valid_args <- names(formals(proxy::dist))
    }

    dots[intersect(names(dots), valid_args)]
}

# Function to split indices for the symmetric, parallel, proxy case
#' @importFrom parallel splitIndices
#'
split_parallel_symmetric <- function(n, num_workers, adjust = 0L) {
    if (num_workers <= 2L || n <= 4L) {
        mid_point <- as.integer(n / 2)
        # indices for upper part of the lower triangular
        ul_trimat <- 1L:mid_point + adjust
        # indices for lower part of the lower triangular
        ll_trimat <- (mid_point + 1L):n + adjust
        # put triangular parts together for load balance
        trimat <- list(ul = ul_trimat, ll = ll_trimat)
        attr(trimat, "trimat") <- TRUE
        trimat <- list(trimat)
        mid_point <- mid_point + adjust
        attr(ul_trimat, "rows") <- ll_trimat
        mat <- list(ul_trimat)
        ids <- c(trimat, mat)
    }
    else {
        mid_point <- as.integer(n / 2)
        # recursion
        rec1 <- split_parallel_symmetric(mid_point, as.integer(num_workers / 4), adjust)
        rec2 <- split_parallel_symmetric(n - mid_point, as.integer(num_workers / 4), mid_point + adjust)
        endpoints <- parallel::splitIndices(mid_point, max(length(rec1) + length(rec2), num_workers))
        endpoints <- endpoints[lengths(endpoints) > 0L]
        mat <- lapply(endpoints, function(ids) {
            ids <- ids + adjust
            attr(ids, "rows") <- (mid_point + 1L):n + adjust
            ids
        })
        ids <- c(rec1, rec2, mat)
    }
    chunk_sizes <- unlist(lapply(ids, function(x) {
        if (is.null(attr(x, "trimat"))) length(x) else median(lengths(x))
    }))
    # return
    ids[sort(chunk_sizes, index.return = TRUE)$ix]
}

# calculate only half the distance matrix in parallel
#' @importFrom bigmemory attach.big.matrix
#' @importFrom proxy dist
#'
parallel_symmetric <- function(d_desc, ids, x, distance, dots) {
    dd <- bigmemory::attach.big.matrix(d_desc)
    if (isTRUE(attr(ids, "trimat"))) {
        # assign upper part of lower triangular
        ul <- ids$ul
        if (length(ul) > 1L) {
            dd[ul,ul] <- base::as.matrix(quoted_call(
                proxy::dist,
                x = x[ul],
                y = NULL,
                method = distance,
                dots = dots
            ))
        }

        # assign lower part of lower triangular
        ll <- ids$ll
        if (length(ll) > 1L) {
            dd[ll,ll] <- base::as.matrix(quoted_call(
                proxy::dist,
                x = x[ll],
                y = NULL,
                method = distance,
                dots = dots
            ))
        }
    }
    else {
        # assign matrix chunks
        rows <- attr(ids, "rows")
        mat_chunk <- base::as.matrix(quoted_call(
            proxy::dist,
            x = x[rows],
            y = x[ids],
            method = distance,
            dots = dots
        ))

        dd[rows,ids] <- mat_chunk
        dd[ids,rows] <- t(mat_chunk)
    }
}

# ==================================================================================================
# Return a custom distance function that calls registered functions of proxy
# ==================================================================================================

#' @importFrom bigmemory big.matrix
#' @importFrom bigmemory describe
#' @importFrom proxy dist
#' @importFrom proxy pr_DB
#'
ddist2 <- function(distance, control) {
    # I need to re-register any custom distances in each parallel worker
    dist_entry <- proxy::pr_DB$get_entry(distance)
    symmetric <- isTRUE(control$symmetric)

    # variables/functions from the parent environments that should be exported
    export <- c("check_consistency", "quoted_call", "parallel_symmetric", "distance", "dist_entry")

    ret <- function(result, ...) {
        ret <- structure(result, method = toupper(distance), ...)
        if (!is.null(attr(ret, "call"))) {
            attr(ret, "call") <- NULL
        }
        ret
    }

    # Closures capture the values of the objects from the environment where they're created
    distfun <- function(x, centroids = NULL, ...) {
        x <- tslist(x)
        if (!is.null(centroids)) centroids <- tslist(centroids)

        if (length(x) == 1L && is.null(centroids)) {
            return(ret(base::matrix(0, 1L, 1L),
                       class = "crossdist",
                       dimnames = list(names(x), names(x))))
        }

        if (!is.null(control$distmat)) {
            return(ret(use_distmat(control$distmat, x, centroids)))
        }

        dots <- get_dots(dist_entry, x, centroids, ...)

        if (!dist_entry$loop) {
            # CUSTOM LOOP, LET THEM HANDLE OPTIMIZATIONS
            dm <- base::as.matrix(quoted_call(
                proxy::dist, x = x, y = centroids, method = distance, dots = dots
            ))

            if (isTRUE(dots$pairwise)) {
                dim(dm) <- NULL
                return(ret(dm, class = "pairdist"))
            }
            else {
                return(ret(dm, class = "crossdist"))
            }
        }

        if (is.null(centroids) && symmetric && !isTRUE(dots$pairwise)) {
            if (foreach::getDoParWorkers() > 1L) {
                # WHOLE SYMMETRIC DISTMAT IN PARALLEL
                # Only half of it is computed
                # proxy can do this if y = NULL, but not in parallel
                len <- length(x)

                # undo bigmemory's seed change, backwards reproducibility
                seed <- get0(".Random.seed", .GlobalEnv, mode = "integer")
                d <- bigmemory::big.matrix(len, len, "double", 0)
                d_desc <- bigmemory::describe(d)
                assign(".Random.seed", seed, .GlobalEnv)

                ids <- integer() # "initialize", so CHECK doesn't complain about globals
                foreach(
                    ids = split_parallel_symmetric(len, foreach::getDoParWorkers()),
                    .combine = c,
                    .multicombine = TRUE,
                    .noexport = c("d"),
                    .packages = c(control$packages, "bigmemory"),
                    .export = export
                ) %op% {
                    if (!check_consistency(dist_entry$names[1L], "dist")) {
                        do.call(proxy::pr_DB$set_entry, dist_entry, TRUE) # nocov
                    }

                    parallel_symmetric(d_desc, ids, x, distance, dots)
                    NULL
                }

                # coerce to normal matrix
                return(ret(d[,], class = "crossdist", dimnames = list(names(x), names(x))))
            }
            else {
                # WHOLE SYMMETRIC DISTMAT WITH CUSTOM LOOP OR SEQUENTIAL proxy LOOP
                dm <- base::as.matrix(quoted_call(
                    proxy::dist, x = x, y = NULL, method = distance, dots = dots
                ))

                return(ret(dm, class = "crossdist"))
            }
        }

        # WHOLE DISTMAT OR SUBDISTMAT OR NOT SYMMETRIC
        if (is.null(centroids)) centroids <- x
        dim_names <- list(names(x), names(centroids)) # x and centroids may change in parallel!
        x <- split_parallel(x)

        if (isTRUE(dots$pairwise)) {
            centroids <- split_parallel(centroids)
            validate_pairwise(x, centroids)
            combine <- c
        }
        else {
            centroids <- lapply(1L:foreach::getDoParWorkers(), function(dummy) { centroids })
            if (length(centroids) > length(x)) centroids <- centroids[1L:length(x)] # nocov
            combine <- rbind
        }

        d <- foreach(
            x = x, centroids = centroids,
            .combine = combine,
            .multicombine = TRUE,
            .packages = control$packages,
            .export = export
        ) %op% {
            if (!check_consistency(dist_entry$names[1L], "dist")) {
                do.call(proxy::pr_DB$set_entry, dist_entry, TRUE)
            }

            quoted_call(proxy::dist, x = x, y = centroids, method = distance, dots = dots)
        }

        if (isTRUE(dots$pairwise)) {
            attr(d, "class") <- "pairdist"
        }
        else {
            attr(d, "class") <- "crossdist"
            attr(d, "dimnames") <- dim_names
        }
        # return
        ret(d)
    }

    # return closure
    distfun
}
########

# pr_DB
# https://github.com/asardaes/dtwclust/blob/bc6008c26b46ff7105b59135dd9ad004a99f7d87/R/CLUSTERING-compare-clusterings.R
########
#' Helper function for preprocessing/distance/centroid configurations
#'
#' Create preprocessing, distance and centroid configurations for [compare_clusterings_configs()].
#'
#' @export
#' @importFrom dplyr bind_rows
#'
#' @param type Which type of function is being targeted by this configuration.
#' @param ... Any number of named lists with functions and arguments that will be shared by all
#'   clusterings. See details.
#' @param partitional A named list of lists with functions and arguments for partitional
#'   clusterings.
#' @param hierarchical A named list of lists with functions and arguments for hierarchical
#'   clusterings.
#' @param fuzzy A named list of lists with functions and arguments for fuzzy clusterings.
#' @param tadpole A named list of lists with functions and arguments for TADPole clusterings.
#' @param share.config A character vector specifying which clusterings should include the shared
#'   lists (the ones specified in `...`). It must be any combination of (possibly abbreviated):
#'   partitional, hierarchical, fuzzy, tadpole.
#'
#' @details
#'
#' The named lists are interpreted in the following way: the name of the list will be considered to
#' be a function name, and the elements of the list will be the possible parameters for the
#' function. Each function must have at least an empty list. The parameters may be vectors that
#' specify different values to be tested.
#'
#' For preprocessing, the special name `none` signifies no preprocessing.
#'
#' For centroids, the special name `default` leaves the centroid unspecified.
#'
#' Please see the examples in [compare_clusterings()] to see how this is used.
#'
#' @return
#'
#' A list for each clustering, each of which includes a data frame with the computed configurations.
#'
pdc_configs <- function(type = c("preproc", "distance", "centroid"), ...,
                        partitional = NULL, hierarchical = NULL, fuzzy = NULL, tadpole = NULL,
                        share.config = c("p", "h", "f", "t"))
{
    type <- match.arg(type)
    shared <- list(...)
    specific <- list(partitional = partitional,
                     hierarchical = hierarchical,
                     fuzzy = fuzzy,
                     tadpole = tadpole)
    specific <- specific[!sapply(specific, is.null)]
    share_missing <- missing(share.config)
    share.config <- match.arg(share.config, supported_clusterings, TRUE)
    if (type == "distance") {
        if (!is.null(specific$tadpole) || (!share_missing && "tadpole" %in% share.config))
            warning("TADPole ignores distance configurations.")
        specific$tadpole <- NULL
        share.config <- setdiff(share.config, "tadpole")
    }

    # ==============================================================================================
    # Shared configs
    # ==============================================================================================

    if (length(shared) > 0L && length(share.config) > 0L) {
        # careful, singular and plural below
        shared_cfg <- Map(shared, names(shared), f = function(shared_args, fun) {
            cfg <- quoted_call(expand.grid, foo = fun, stringsAsFactors = FALSE, dots = shared_args)
            names(cfg)[1L] <- type
            cfg
        })
        shared_cfg <- dplyr::bind_rows(shared_cfg)
        shared_cfgs <- lapply(share.config, function(dummy) { shared_cfg })
        names(shared_cfgs) <- share.config
        shared_cfgs <- shared_cfgs[setdiff(share.config, names(specific))]
    }
    else {
        shared_cfg <- NULL
        shared_cfgs <- list()
    }

    # ==============================================================================================
    # Specific configs
    # ==============================================================================================

    if (length(specific) > 0L) {
        cfgs <- Map(specific, names(specific), f = function(config, clus_type) {
            config_names <- names(config)
            if (!is.list(config) || is.null(config_names))
                stop("All parameters must be named lists.") # nocov
            cfg <- Map(config, config_names, f = function(config_args, fun) {
                cfg <- quoted_call(expand.grid, foo = fun, stringsAsFactors = FALSE, dots = config_args)
                names(cfg)[1L] <- type
                cfg
            })
            cfg <- dplyr::bind_rows(cfg)
            if (clus_type %in% share.config)
                cfg <- dplyr::bind_rows(shared_cfg, cfg) # singular shared
            cfg
        })
        cfgs <- c(cfgs, shared_cfgs) # plural shared
    }
    else {
        cfgs <- shared_cfgs
    }
    # return
    cfgs
}

#' Create clustering configurations.
#'
#' Create configurations for [compare_clusterings()]
#'
#' @export
#' @importFrom dplyr bind_cols
#'
#' @param k A numeric vector with one or more elements specifying the number of clusters to test.
#' @param types Clustering types. It must be any combination of (possibly abbreviated): partitional,
#'   hierarchical, fuzzy, tadpole.
#' @param controls A named list of [tsclust-controls]. `NULL` means defaults. See details.
#' @param preprocs Preprocessing configurations. See details.
#' @param distances Distance configurations. See details.
#' @param centroids Centroid configurations. See details.
#' @param no.expand A character vector indicating parameters that should *not* be expanded between
#'   [pdc_configs()] configurations. See examples.
#'
#' @details
#'
#' Preprocessing, distance and centroid configurations are specified with the helper function
#' [pdc_configs()], refer to the examples in [compare_clusterings()] to see how this is used.
#'
#' The controls list may be specified with the usual [tsclust-controls] functions. The names of the
#' list must correspond to "partitional", "hierarchical", "fuzzy" or "tadpole" clustering. Again,
#' please refer to the examples in [compare_clusterings()].
#'
#' @return
#'
#' A list for each clustering type, each of which includes a data frame with the computed and merged
#' configurations. Each data frame has an extra attribute `num.configs` specifying the number of
#' configurations.
#'
#' @examples
#'
#' # compare this with leaving no.expand empty
#' compare_clusterings_configs(
#'     distances = pdc_configs("d", dtw_basic = list(window.size = 1L:2L, norm = c("L1", "L2"))),
#'     centroids = pdc_configs("c", dba = list(window.size = 1L:2L, norm = c("L1", "L2"))),
#'     no.expand = c("window.size", "norm")
#' )
#'
compare_clusterings_configs <- function(types = c("p", "h", "f"), k = 2L, controls = NULL,
                                        preprocs = pdc_configs("preproc", none = list()),
                                        distances = pdc_configs("distance", dtw_basic = list()),
                                        centroids = pdc_configs("centroid", default = list()),
                                        no.expand = character(0L))
{
    # ==============================================================================================
    # Start
    # ==============================================================================================

    types <- match.arg(types, supported_clusterings, TRUE)

    # ----------------------------------------------------------------------------------------------
    # Check controls specification
    # ----------------------------------------------------------------------------------------------

    if (is.null(controls)) {
        controls <- lapply(types, function(type) { do.call(paste0(type, "_control"), list()) })
        names(controls) <- types
    }
    else if (!is.list(controls) || is.null(names(controls))) {
        stop("The 'controls' argument must be NULL or a named list")
    }
    else if (!all(types %in% names(controls))) {
        stop("The names of the 'controls' argument do not correspond to the provided 'types'")
    }
    else {
        controls <- controls[intersect(names(controls), types)]
    }

    # ----------------------------------------------------------------------------------------------
    # Check preprocessings specification
    # ----------------------------------------------------------------------------------------------

    if (missing(preprocs))
        force(preprocs)
    else if (!is.list(preprocs) || (length(preprocs) > 0L && is.null(names(preprocs))))
        stop("The 'preprocs' argument must be a list with named elements")
    else if (!all(types %in% names(preprocs)))
        stop("The names of the 'preprocs' argument do not correspond to the provided 'types'")

    preprocs <- preprocs[intersect(names(preprocs), types)]

    # ----------------------------------------------------------------------------------------------
    # Check distance specification
    # ----------------------------------------------------------------------------------------------

    if (missing(distances))
        force(distances)
    else if (!is.list(distances) || (length(distances) > 0L && is.null(names(distances))))
        stop("The 'distances' argument must be a list with named elements")
    else if (!all(setdiff(types, "tadpole") %in% names(distances)))
        stop("The names of the 'distances' argument do not correspond to the provided 'types'")

    distances <- distances[intersect(names(distances), types)]

    # ----------------------------------------------------------------------------------------------
    # Check centroids specification
    # ----------------------------------------------------------------------------------------------

    if (missing(centroids))
        force(centroids)
    else if (!is.list(centroids) || (length(centroids) > 0L && is.null(names(centroids))))
        stop("The 'centroids' argument must be a list with named elements")
    else if (!all(types %in% names(centroids)))
        stop("The names of the 'centroids' argument do not correspond to the provided 'types'")

    centroids <- centroids[intersect(names(centroids), types)]

    # ==============================================================================================
    # Create configs
    # ==============================================================================================

    # return here
    Map(types, controls[types], preprocs[types], distances[types], centroids[types],
        f = function(type, control, preproc, distance, centroid) {
            # --------------------------------------------------------------------------------------
            # Control configs
            # --------------------------------------------------------------------------------------

            if (class(control) != control_classes[type]) stop("Invalid ", type, " control") # nocov

            # if it's within a list, it's to prevent expansion
            cfg <- switch(
                type,
                partitional = {
                    quoted_call(
                        expand.grid,
                        k = list(k),
                        pam.precompute = control$pam.precompute,
                        iter.max = control$iter.max,
                        nrep = control$nrep,
                        symmetric = control$symmetric,
                        version = control$version,
                        stringsAsFactors = FALSE
                    )
                },
                hierarchical = {
                    quoted_call(
                        expand.grid,
                        k = list(k),
                        method = list(control$method),
                        symmetric = control$symmetric,
                        stringsAsFactors = FALSE
                    )
                },
                fuzzy = {
                    quoted_call(
                        expand.grid,
                        k = list(k),
                        fuzziness = control$fuzziness,
                        iter.max = control$iter.max,
                        delta = control$delta,
                        symmetric = control$symmetric,
                        version = control$version,
                        stringsAsFactors = FALSE
                    )
                },
                tadpole = {
                    quoted_call(
                        expand.grid,
                        k = list(k),
                        dc = list(control$dc),
                        window.size = control$window.size,
                        lb = control$lb,
                        stringsAsFactors = FALSE
                    )
                }
            )

            # --------------------------------------------------------------------------------------
            # Merge configs
            # --------------------------------------------------------------------------------------

            need_adjustment <- character(0L)

            # preproc
            if (!is.null(preproc) && nrow(preproc)) {
                nms <- names(preproc)
                if (any(nms %in% no.expand)) need_adjustment <- c(need_adjustment, "preproc")
                nms_args <- nms != "preproc" & !(nms %in% no.expand)
                if (any(nms_args)) names(preproc)[nms_args] <- paste0(nms[nms_args], "_preproc")
                cfg <- base::merge(cfg, preproc, all = TRUE)
            }
            # distance
            if (type != "tadpole" && !is.null(distance) && nrow(distance)) {
                nms <- names(distance)
                if (any(nms %in% no.expand)) need_adjustment <- c(need_adjustment, "distance")
                nms_args <- nms != "distance" & !(nms %in% no.expand)
                if (any(nms_args)) names(distance)[nms_args] <- paste0(nms[nms_args], "_distance")
                cfg <- base::merge(cfg, distance, all = TRUE)
            }
            # centroid
            if (!is.null(centroid) && nrow(centroid)) {
                nms <- names(centroid)
                if (any(nms %in% no.expand)) need_adjustment <- c(need_adjustment, "centroid")
                nms_args <- nms != "centroid" & !(nms %in% no.expand)
                if (any(nms_args)) names(centroid)[nms_args] <- paste0(nms[nms_args], "_centroid")
                cfg <- base::merge(cfg, centroid, all = TRUE)
            }
            # special case: tadpole
            tadpole_controls <- names(formals(tadpole_control))
            if (type == "tadpole" && any(no.expand %in% tadpole_controls)) {
                need_adjustment <- c(need_adjustment, "tadpole")
                tadpole <- cfg[, tadpole_controls, drop = FALSE]
            }

            # adjust no.expand columns
            if (length(need_adjustment) > 0L) {
                adjust_cols <- cfg[, no.expand, drop = FALSE]
                cfg <- cfg[, -which(names(cfg) %in% no.expand), drop = FALSE]
                adjusted_cols <- lapply(need_adjustment, function(suffix) {
                    pdc_cfg <- get_from_callers(suffix, mode = "list")
                    cols <- intersect(names(pdc_cfg), no.expand)
                    adjusted_cols <- adjust_cols[, cols, drop = FALSE]
                    if (suffix != "tadpole")
                        names(adjusted_cols) <- paste0(names(adjusted_cols), "_", suffix)
                    adjusted_cols
                })
                cfg <- dplyr::bind_cols(cfg, adjusted_cols)
            }

            # for info
            attr(cfg, "num.configs") <- switch(
                type,
                partitional = length(k) * sum(cfg$nrep),
                hierarchical = length(k) * length(cfg$method[[1L]]) * nrow(cfg),
                fuzzy = length(k) * nrow(cfg),
                tadpole = length(k) * length(cfg$dc[[1L]]) * nrow(cfg)
            )
            # return Map
            cfg
        })
}

#' Compare different clustering configurations
#'
#' Compare many different clustering algorithms with support for parallelization.
#'
#' @export
#' @importFrom dplyr bind_rows
#' @importFrom dplyr inner_join
#' @importFrom proxy pr_DB
#'
#' @param series A list of series, a numeric matrix or a data frame. Matrices and data frames are
#'   coerced to a list row-wise (see [tslist()]).
#' @param types Clustering types. It must be any combination of (possibly abbreviated):
#'   "partitional", "hierarchical", "fuzzy", "tadpole."
#' @param configs The list of data frames with the desired configurations to run. See
#'   [pdc_configs()] and [compare_clusterings_configs()].
#' @param seed Seed for random reproducibility.
#' @param trace Logical indicating that more output should be printed to screen.
#' @param ... Further arguments for [tsclust()], `score.clus` or `pick.clus`.
#' @param score.clus A function that gets the list of results (and `...`) and scores each one. It
#'   may also be a named list of functions, one for each type of clustering. See Scoring section.
#' @param pick.clus A function to pick the best result. See Picking section.
#' @param shuffle.configs Randomly shuffle the order of configs, which can be useful to balance load
#'   when using parallel computation.
#' @param return.objects Logical indicating whether the objects returned by [tsclust()] should be
#'   given in the result.
#' @param packages A character vector with the names of any packages needed for any functions used
#'   (distance, centroid, preprocessing, etc.). The name "dtwclust" is added automatically. Relevant
#'   for parallel computation.
#' @param .errorhandling This will be passed to [foreach::foreach()]. See Parallel section below.
#'
#' @details
#'
#' This function calls [tsclust()] with different configurations and evaluates the results with the
#' provided functions. Parallel support is included. See the examples.
#'
#' Parameters specified in `configs` whose values are `NA` will be ignored automatically.
#'
#' The scoring and picking functions are for convenience, if they are not specified, the `scores`
#' and `pick` elements of the result will be `NULL`.
#'
#' See [repeat_clustering()] for when `return.objects = FALSE`.
#'
#' @return
#'
#' A list with:
#'
#' - `results`: A list of data frames with the flattened configs and the corresponding scores
#' returned by `score.clus`.
#' - `scores`: The scores given by `score.clus`.
#' - `pick`: The object returned by `pick.clus`.
#' - `proc_time`: The measured execution time, using [base::proc.time()].
#' - `seeds`: A list of lists with the random seeds computed for each configuration.
#'
#' The cluster objects are also returned if `return.objects` `=` `TRUE`.
#'
#' @section Parallel computation:
#'
#'   The configurations for each clustering type can be evaluated in parallel (multi-processing)
#'   with the \pkg{foreach} package. A parallel backend can be registered, e.g., with
#'   \pkg{doParallel}.
#'
#'   If the `.errorhandling` parameter is changed to "pass" and a custom `score.clus` function is
#'   used, said function should be able to deal with possible error objects.
#'
#'   If it is changed to "remove", it might not be possible to attach the scores to the results data
#'   frame, or it may be inconsistent. Additionally, if `return.objects` is `TRUE`, the names given
#'   to the objects might also be inconsistent.
#'
#'   Parallelization can incur a lot of deep copies of data when returning the cluster objects,
#'   since each one will contain a copy of `datalist`. If you want to avoid this, consider
#'   specifying `score.clus` and setting `return.objects` to `FALSE`, and then using
#'   [repeat_clustering()].
#'
#' @section Scoring:
#'
#'   The clustering results are organized in a *list of lists* in the following way (where only
#'   applicable `types` exist; first-level list names in bold):
#'
#'   - **partitional** - list with
#'     + Clustering results from first partitional config
#'     + etc.
#'   - **hierarchical** - list with
#'     + Clustering results from first hierarchical config
#'     + etc.
#'   - **fuzzy** - list with
#'     + Clustering results from first fuzzy config
#'     + etc.
#'   - **tadpole** - list with
#'     + Clustering results from first tadpole config
#'     + etc.
#'
#'   If `score.clus` is a function, it will be applied to the available partitional, hierarchical,
#'   fuzzy and/or tadpole results via:
#'
#'   ```
#'   scores <- lapply(list_of_lists, score.clus, ...)
#'   ```
#'
#'   Otherwise, `score.clus` should be a list of functions with the same names as the list above, so
#'   that `score.clus$partitional` is used to score `list_of_lists$partitional` and so on (via
#'   [base::Map()]).
#'
#'   Therefore, the scores returned shall always be a list of lists with first-level names as above.
#'
#' @section Picking:
#'
#'   If `return.objects` is `TRUE`, the results' data frames and the list of [TSClusters-class]
#'   objects are given to `pick.clus` as first and second arguments respectively, followed by `...`.
#'   Otherwise, `pick.clus` will receive only the data frames and the contents of `...` (since the
#'   objects will not be returned by the preceding step).
#'
#' @section Limitations:
#'
#'   Note that the configurations returned by the helper functions assign special names to
#'   preprocessing/distance/centroid arguments, and these names are used internally to recognize
#'   them.
#'
#'   If some of these arguments are more complex (e.g. matrices) and should *not* be expanded,
#'   consider passing them directly via the ellipsis (`...`) instead of using [pdc_configs()]. This
#'   assumes that said arguments can be passed to all functions without affecting their results.
#'
#'   The distance matrices (if calculated) are not re-used across configurations. Given the way the
#'   configurations are created, this shouldn't matter, because clusterings with arguments that can
#'   use the same distance matrix are already grouped together by [compare_clusterings_configs()]
#'   and [pdc_configs()].
#'
#' @author Alexis Sarda-Espinosa
#'
#' @seealso
#'
#' [compare_clusterings_configs()], [tsclust()]
#'
#' @example man-examples/comparison-examples.R
#'
compare_clusterings <- function(series = NULL, types = c("p", "h", "f", "t"),
                                configs = compare_clusterings_configs(types),
                                seed = NULL, trace = FALSE, ...,
                                score.clus = function(...) stop("No scoring"),
                                pick.clus = function(...) stop("No picking"),
                                shuffle.configs = FALSE, return.objects = FALSE,
                                packages = character(0L), .errorhandling = "stop")
{
    # ==============================================================================================
    # Start
    # ==============================================================================================

    tic <- proc.time()
    handle_rngkind() # UTILS-rng.R
    set.seed(seed)
    score_missing <- missing(score.clus)
    pick_missing <- missing(pick.clus)

    # nocov start
    if (is.null(series))
        stop("No series provided.")

    if (!return.objects && score_missing)
        stop("Returning no objects and specifying no scoring function would return no useful results.")

    types <- match.arg(types, supported_clusterings, TRUE)
    .errorhandling <- match.arg(.errorhandling, c("stop", "remove", "pass"))

    # coerce to list if necessary
    if (is.data.frame(series) || !is.list(series))
        series <- tslist(series, TRUE)
    check_consistency(series, "vltslist")

    if (!is.function(score.clus) && !(is.list(score.clus) && all(sapply(score.clus, is.function))))
        stop("Invalid evaluation function(s)")
    else if (is.list(score.clus)) {
        if (!all(types %in% names(score.clus)))
            stop("The names of the 'score.clus' argument do not correspond to the provided 'types'")

        score.clus <- score.clus[types]
    }

    if (!is.function(pick.clus))
        stop("Invalid pick function") # nocov end

    # ----------------------------------------------------------------------------------------------
    # Misc parameters
    # ----------------------------------------------------------------------------------------------

    packages <- unique(c("dtwclust", packages))
    dots <- list(...)
    configs <- configs[types]
    if (any(sapply(configs, is.null)))
        stop("The configuration for one of the chosen clustering types is missing.") # nocov
    if (shuffle.configs) {
        configs <- lapply(configs, function(config) {
            config[sample(nrow(config)), , drop = FALSE]
        })
    }

    # ----------------------------------------------------------------------------------------------
    # Obtain random seeds
    # ----------------------------------------------------------------------------------------------

    num_seeds <- cumsum(sapply(configs, nrow))
    seeds <- rng_seq(num_seeds[length(num_seeds)], seed = seed, simplify = FALSE) # UTILS-rng.R
    seeds <- Map(c(1L, num_seeds[-length(num_seeds)] + 1L), num_seeds,
                 f = function(first, last) { seeds[first:last] })
    setnames_inplace(seeds, names(configs)) # UTILS-utils.R

    # ==============================================================================================
    # Preprocessings
    # ==============================================================================================

    if (trace) message("=================================== Preprocessing ",
                       "series ===================================\n")

    processed_series <- Map(configs, types, f = function(config, type) {
        preproc_cols <- grepl("_?preproc$", names(config))
        preproc_df <- unique(config[, preproc_cols, drop = FALSE])
        preproc_args <- grepl("_preproc$", names(preproc_df))

        if (trace) {
            message("-------------- Applying ", type, " preprocessings: --------------")
            print(preproc_df)
        }

        config$.preproc_id_ <- seq_len(nrow(config))
        lapply(seq_len(nrow(preproc_df)), function(i) {
            preproc_char <- preproc_df$preproc[i]
            if (preproc_char != "none") {
                # find all configs that have this preproc to assign them as attribute at the end
                df <- dplyr::inner_join(config,
                                        preproc_df[i, , drop = FALSE],
                                        by = names(config)[preproc_cols])

                preproc_fun <- get_from_callers(preproc_char, "function")

                if (any(preproc_args)) {
                    this_config <- preproc_df[i, preproc_args, drop = FALSE]
                    names(this_config) <- sub("_preproc$", "", names(this_config))
                    preproc_args <- as.list(this_config)
                    preproc_args <- preproc_args[!sapply(preproc_args, is.na)]
                }
                else
                    preproc_args <- list()

                ret <- quoted_call(preproc_fun, series, dots = preproc_args)
                attr(ret, "config_ids") <- df$.preproc_id_
            }
            else {
                ret <- series
            }
            ret
        })
    })

    # UTILS-utils.R
    setnames_inplace(processed_series, names(configs))

    # ==============================================================================================
    # Clusterings
    # ==============================================================================================

    if (trace) cat("\n")
    objs_by_type <- Map(configs, names(configs), seeds, f = function(config, type, seeds) {
        if (trace) message("=================================== Performing ",
                           type,
                           " clusterings ===================================\n")
        series <- processed_series[[type]]

        # ------------------------------------------------------------------------------------------
        # distance entries to re-register in parallel workers
        # ------------------------------------------------------------------------------------------

        if (type != "tadpole") {
            dist_names <- unique(config$distance)
            dist_entries <- lapply(dist_names, function(dist) { proxy::pr_DB$get_entry(dist) })
            setnames_inplace(dist_entries, dist_names)
        }

        # ------------------------------------------------------------------------------------------
        # export any necessary preprocessing and centroid functions
        # ------------------------------------------------------------------------------------------

        custom_preprocs <- setdiff(unique(config$preproc), "none")
        custom_centroids <- setdiff(unique(config$centroid), c("default", centroids_included))

        for (custom_preproc in custom_preprocs)
            assign(custom_preproc, get_from_callers(custom_preproc, "function"))

        for (custom_centroid in custom_centroids)
            assign(custom_centroid, get_from_callers(custom_centroid, "function"))

        export <- c("trace", "score.clus", "return.objects",
                    "dots",
                    "centroids_included",
                    "check_consistency", "quoted_call", "enlist", "subset_dots", "get_from_callers",
                    "setnames_inplace",
                    custom_preprocs, custom_centroids)

        # ------------------------------------------------------------------------------------------
        # perform clusterings
        # ------------------------------------------------------------------------------------------

        force(seeds)
        i <- nrow(config)
        objs <- foreach::foreach(
            i = seq_len(i),
            .combine = c,
            .multicombine = TRUE,
            .packages = packages,
            .export = export,
            .errorhandling = .errorhandling
        ) %op% {
            cfg <- config[i, , drop = FALSE]
            seed <- seeds[[i]]
            if (trace) {
                message("-------------- Using configuration: --------------")
                print(cfg)
            }

            # ----------------------------------------------------------------------------------
            # obtain args from configuration
            # ----------------------------------------------------------------------------------

            args <- lapply(c("preproc", "distance", "centroid"), function(func) {
                col_ids <- grepl(paste0("_", func, "$"), names(cfg))
                if (cfg[[func]] != "none" && any(col_ids)) {
                    this_args <- as.list(cfg[, col_ids, drop = FALSE])
                    names(this_args) <- sub(paste0("_", func, "$"),
                                            "",
                                            names(this_args))
                    # return
                    this_args[!sapply(this_args, is.na)]
                }
                else list()
            })

            setnames_inplace(args, c("preproc", "dist", "cent"))
            args <- do.call(tsclust_args, args, TRUE)

            # ----------------------------------------------------------------------------------
            # controls for this configuration
            # ----------------------------------------------------------------------------------

            control_fun <- match.fun(paste0(type, "_control"))
            control_args <- subset_dots(as.list(cfg), control_fun)
            control_args <- lapply(control_args, unlist, recursive = FALSE)
            control <- do.call(control_fun, control_args, TRUE)

            # ----------------------------------------------------------------------------------
            # get processed series
            # ----------------------------------------------------------------------------------

            preproc_char <- cfg$preproc

            config_ids <- lapply(series, attr, which = "config_ids")
            if (preproc_char == "none")
                this_series <- which(sapply(config_ids, is.null))
            else
                this_series <- which(sapply(config_ids, function(cfg_id) { i %in% cfg_id }))

            if (length(this_series) > 1L) # nocov start
                stop("Could not find unique processed series for ", type,
                     " clustering and config row=", i)
            else # nocov end
                this_series <- series[[this_series]]

            # ----------------------------------------------------------------------------------
            # distance entry to re-register in parallel worker
            # ----------------------------------------------------------------------------------

            if (type != "tadpole") {
                distance <- cfg$distance
                dist_entry <- dist_entries[[distance]]
                if (!check_consistency(dist_entry$names[1L], "dist"))
                    do.call(proxy::pr_DB$set_entry, dist_entry, TRUE) # nocov
            }
            else distance <- NULL # dummy

            # ----------------------------------------------------------------------------------
            # centroid for this configuration
            # ----------------------------------------------------------------------------------

            centroid_char <- cfg$centroid

            # ----------------------------------------------------------------------------------
            # call tsclust
            # ----------------------------------------------------------------------------------

            this_args <- enlist(series = this_series,
                                type = type,
                                k = unlist(cfg$k),
                                distance = distance,
                                seed = seed,
                                trace = trace,
                                args = args,
                                control = control,
                                error.check = FALSE,
                                dots = dots)

            if (type == "tadpole")
                this_args$distance <- NULL

            if (centroid_char == "default") {
                # do not specify centroid
                tsc <- do.call(tsclust, this_args, TRUE)
            }
            else if (type %in% c("partitional", "fuzzy") && centroid_char %in% centroids_included) {
                # with included centroid
                tsc <- quoted_call(tsclust, centroid = centroid_char, dots = this_args)
            }
            else {
                # with centroid function
                tsc <- quoted_call(tsclust,
                                   centroid = get_from_callers(centroid_char, "function"),
                                   dots = this_args)
            }

            if (inherits(tsc, "TSClusters"))
                tsc <- list(tsc)

            ret <- lapply(tsc, function(tsc) {
                tsc@preproc <- preproc_char
                if (preproc_char != "none")
                    tsc@family@preproc <- get_from_callers(preproc_char, "function")
                if (centroid_char != "default")
                    tsc@centroid <- centroid_char
                tsc
            })

            # ----------------------------------------------------------------------------------
            # evaluate
            # ----------------------------------------------------------------------------------

            if (!return.objects) {
                if (!is.function(score.clus)) score.clus <- score.clus[[type]]
                ret <- list(quoted_call(score.clus, ret, dots = dots))
            }
            # return config result from foreach()
            ret
        }

        class(objs) <- NULL
        if (.errorhandling == "pass") {
            failed_cfgs <- sapply(objs, function(obj) { !inherits(obj, "TSClusters") })
            if (any(failed_cfgs)) {
                warning("At least one of the ", type, " configurations resulted in an error.")

                # a simple error is a list with 2 elements: message and call, so I need to re-pack
                # each pair of elements in a single element of the objs list
                which_failed <- which(failed_cfgs)[seq(from = 1L, by = 2L, length.out = sum(failed_cfgs) / 2)]
                for (failed_cfg_id in which_failed) {
                    objs[[failed_cfg_id]] <- structure(objs[failed_cfg_id:(failed_cfg_id + 1L)],
                                                       class = c("simpleError", "error", "condition"))
                }
                names(objs)[which_failed] <- paste0("failure_", seq_along(which_failed))
                objs[which_failed + 1L] <- NULL
            }
        }

        objs
    })

    # ==============================================================================================
    # Evaluations
    # ==============================================================================================

    if (return.objects) {
        if (is.function(score.clus))
            scores <- try(lapply(objs_by_type, score.clus, ...), silent = TRUE)
        else
            scores <- try(mapply(objs_by_type, score.clus[names(objs_by_type)],
                                 SIMPLIFY = FALSE,
                                 MoreArgs = dots,
                                 FUN = function(objs, score_fun, ...) { score_fun(objs, ...) }),
                          silent = TRUE)

        if (inherits(scores, "try-error")) {
            if (!score_missing) warning("The score.clus function(s) did not execute successfully:\n",
                                        attr(scores, "condition")$message)
            scores <- NULL
        }
    }
    else {
        scores <- lapply(objs_by_type, function(objs) {
            failed_cfgs <- sapply(objs, function(obj) { inherits(obj, "error") })

            if (any(failed_cfgs))
                passed_objs <- objs[!failed_cfgs]
            else
                passed_objs <- objs

            if (length(passed_objs) == 0L) return(objs)

            if (any(sapply(passed_objs, function(score) { is.null(dim(score)) })))
                unlist(passed_objs, recursive = FALSE)
            else
                dplyr::bind_rows(lapply(passed_objs, base::as.data.frame))
        })
    }

    # ==============================================================================================
    # Data frame with results
    # ==============================================================================================

    # create initial IDs
    i_cfg <- 1L
    config_ids <- lapply(sapply(configs, nrow), function(nr) {
        ids <- seq(from = i_cfg, by = 1L, length.out = nr)
        i_cfg <<- i_cfg + nr
        ids
    })

    # change to config names and assign to seeds (before flattening)
    config_ids <- Map(config_ids, seeds, f = function(ids, seed) {
        nms <- paste0("config", ids)
        try(setnames_inplace(seed, nms), silent = TRUE)
        nms
    })

    # flatten
    configs_out <- Map(configs, config_ids, types, f = function(config, ids, type) {
        config <- data.frame(config_id = ids, config, stringsAsFactors = FALSE)
        k <- unlist(config$k[1L])
        dfs <- switch(
            type,
            partitional = {
                lapply(seq_len(nrow(config)), function(i) {
                    this_config <- config[i, , drop = FALSE]
                    rep <- 1L:this_config$nrep
                    this_config <- this_config[setdiff(names(this_config), c("k", "nrep"))]
                    df <- expand.grid(rep = rep, k = k)
                    make_unique_ids(df, this_config) # see EOF
                })
            },
            hierarchical = {
                lapply(seq_len(nrow(config)), function(i) {
                    this_config <- config[i, , drop = FALSE]
                    method <- unlist(this_config$method)
                    this_config <- this_config[setdiff(names(this_config), c("k", "method"))]
                    df <- expand.grid(k = k, method = method, stringsAsFactors = FALSE)
                    make_unique_ids(df, this_config) # see EOF
                })
            },
            fuzzy = {
                lapply(seq_len(nrow(config)), function(i) {
                    this_config <- config[i, , drop = FALSE]
                    this_config <- this_config[setdiff(names(this_config), c("k"))]
                    df <- expand.grid(k = k)
                    make_unique_ids(df, this_config) # see EOF
                })
            },
            tadpole = {
                lapply(seq_len(nrow(config)), function(i) {
                    this_config <- config[i, , drop = FALSE]
                    dc <- unlist(this_config$dc)
                    this_config <- this_config[setdiff(names(this_config), c("k", "dc"))]
                    df <- expand.grid(k = k, dc = dc)
                    make_unique_ids(df, this_config) # see EOF
                })
            }
        )
        # return Map
        dplyr::bind_rows(dfs)
    })

    # ----------------------------------------------------------------------------------------------
    # Add scores and pick
    # ----------------------------------------------------------------------------------------------

    # in case ordering is required below
    if (shuffle.configs)
        configs_cols <- lapply(configs_out, function(config) {
            setdiff(colnames(config), c("config_id", "rep", "k", "method", "dc",  "window.size", "lb"))
        })

    if (!is.null(scores)) {
        results <- try(Map(configs_out, scores,
                           f = function(config, score) {
                               cbind(config, base::as.data.frame(score))
                           }),
                       silent = TRUE)

        if (inherits(results, "try-error")) {
            warning("The scores could not be appended to the results data frame:\n",
                    attr(results, "condition")$message)
            results <- configs_out
            pick <- NULL
        }
        else {
            if (return.objects) {
                pick <- try(pick.clus(results, objs_by_type, ...), silent = TRUE)
                if (inherits(pick, "try-error")) {
                    if (!pick_missing) warning("The pick.clus function did not execute successfully:\n",
                                               attr(pick, "condition")$message)
                    pick <- NULL
                }
            }
            else {
                pick <- try(pick.clus(results, ...), silent = TRUE)
                if (inherits(pick, "try-error")) {
                    if (!pick_missing) warning("The pick.clus function did not execute successfully:\n",
                                               attr(pick, "condition")$message)
                    pick <- NULL
                }
            }
        }
    }
    else {
        results <- configs_out
        pick <- NULL
    }

    # ==============================================================================================
    # List with all results
    # ==============================================================================================

    results <- list(results = results,
                    scores = scores,
                    pick = pick,
                    proc_time = proc.time() - tic,
                    seeds = seeds)

    if (return.objects) {
        setnames_res <- try(Map(objs_by_type, results$results,
                                f = function(objs, res) { setnames_inplace(objs, res$config_id) }),
                            silent = TRUE)
        if (inherits(setnames_res, "try-error"))
            warning("Could not assign names to returned objects:\n",
                    attr(setnames_res, "condition")$message)
        results <- c(results, objects = objs_by_type)
    }

    if (shuffle.configs)
        results$results <- Map(results$results, configs_cols[names(results$results)],
                               f = function(result, cols) {
                                   order_args <- as.list(result[cols])
                                   names(order_args) <- NULL
                                   result[do.call(base::order, order_args, TRUE), , drop = FALSE]
                               })
    # return results
    results
}

# ==================================================================================================
# compare_clusterings helpers
# ==================================================================================================

make_unique_ids <- function(df, this_config) {
    rownames(this_config) <- NULL
    this_config <- cbind(this_config[, 1L, drop = FALSE], df, this_config[, -1L, drop = FALSE])
    nr <- nrow(this_config)
    if (nr > 1L) this_config$config_id <- paste0(this_config$config_id, "_", 1L:nr)
    this_config
}
########

# pr_DB
# https://github.com/asardaes/dtwclust/blob/bc6008c26b46ff7105b59135dd9ad004a99f7d87/R/S4-tsclustFamily.R
########
#' Class definition for `tsclustFamily`
#'
#' Formal S4 class with a family of functions used in [tsclust()].
#'
#' @exportClass tsclustFamily
#' @importFrom methods setClass
#'
#' @details
#'
#' The custom implementations also handle parallelization.
#'
#' Since the distance function makes use of \pkg{proxy}, it also supports any extra [proxy::dist()]
#' parameters in `...`.
#'
#' The prototype includes the `cluster` function for partitional methods, as well as a pass-through
#' `preproc` function. The initializer expects a control from [tsclust-controls]. See more below.
#'
#' @slot dist The function to calculate the distance matrices.
#' @slot allcent The function to calculate centroids on each iteration.
#' @slot cluster The function used to assign a series to a cluster.
#' @slot preproc The function used to preprocess the data (relevant for [stats::predict()]).
#'
#' @section Distance function:
#'
#'   The family's dist() function works like [proxy::dist()] but supports parallelization and
#'   optimized symmetric calculations. If you like, you can use the function more or less directly,
#'   but provide a control argument when creating the family (see examples). However, bear in mind
#'   the following considerations.
#'
#'   - The second argument is called `centroids` (inconsistent with [proxy::dist()]).
#'   - If `control$distmat` is *not* `NULL`, the function will try to subset it.
#'   - If `control$symmetric` is `TRUE`, `centroids` is `NULL`, *and* there is no argument
#'   `pairwise` that is `TRUE`, only half the distance matrix will be computed.
#'     + If the distance was registered in [proxy::pr_DB] with `loop = TRUE` and more than one
#'     parallel worker is detected, the computation will be in parallel (using multi-processing with
#'     [foreach::foreach()]), otherwise it will be sequential with [proxy::dist()].
#'   - The function always returns a `crossdist` matrix.
#'
#'   Note that all distances implemented as part of \pkg{dtwclust} have custom proxy loops that use
#'   multi-threading independently of \pkg{foreach}, so see their respective documentation to see
#'   what optimizations apply to each one.
#'
#'   For distances *not* included in \pkg{dtwclust}, the symmetric, parallel case mentioned above
#'   makes chunks for parallel workers, but they are not perfectly balanced, so some workers might
#'   finish before the others.
#'
#' @section Centroid function:
#'
#'   The default partitional allcent() function is a closure with the implementations of the
#'   included centroids. The ones for [DBA()], [shape_extraction()] and [sdtw_cent()] can use
#'   multi-process parallelization with [foreach::foreach()]. Its formal arguments are described in
#'   the Centroid Calculation section from [tsclust()].
#'
#' @note
#'
#' This class is meant to group together the relevant functions, but they are **not** linked with
#' each other automatically. In other words, neither `dist` nor `allcent` apply `preproc`. They
#' essentially don't know of each other's existence.
#'
#' @seealso
#'
#' [dtw_basic()], [dtw_lb()], [gak()], [lb_improved()], [lb_keogh()], [sbd()], [sdtw()].
#'
#' @examples
#'
#' \dontrun{
#' data(uciCT)
#' # See "GAK" documentation
#' fam <- new("tsclustFamily", dist = "gak")
#'
#' # This is done with symmetric optimizations, regardless of control$symmetric
#' crossdist <- fam@@dist(CharTraj, window.size = 18L)
#'
#' # This is done without symmetric optimizations, regardless of control$symmetric
#' crossdist <- fam@@dist(CharTraj, CharTraj, window.size = 18L)
#'
#' # For non-dtwclust distances, symmetric optimizations only apply
#' # with an appropriate control AND a single data argument:
#' fam <- new("tsclustFamily", dist = "dtw",
#'            control = partitional_control(symmetric = TRUE))
#' fam@@dist(CharTraj[1L:5L])
#'
#' # If you want the fuzzy family, use fuzzy = TRUE
#' ffam <- new("tsclustFamily", control = fuzzy_control(), fuzzy = TRUE)
#' }
#'
tsclustFamily <- methods::setClass(
    "tsclustFamily",
    slots = c(dist = "function",
              allcent = "function",
              cluster = "function",
              preproc = "function"),
    prototype = prototype(preproc = function(x, ...) { x },
                          cluster = function(distmat = NULL, ...) {
                              max.col(-distmat, "first")
                          })
)

# ==================================================================================================
# Membership update for fuzzy clustering
# ==================================================================================================

f_cluster <- function(distmat, m) {
    cprime <- apply(distmat, 1L, function(dist_row) { sum( (1 / dist_row) ^ (2 / (m - 1)) ) })
    u <- 1 / apply(distmat, 2L, function(dist_col) { cprime * dist_col ^ (2 / (m - 1)) })
    if (is.null(dim(u))) u <- rbind(u) # for predict generic
    u[is.nan(u)] <- 1 # in case fcmdd is used
    u
}

# ==================================================================================================
# Custom initialize
# ==================================================================================================

#' @importFrom methods callNextMethod
#' @importFrom methods initialize
#' @importFrom methods setMethod
#'
setMethod("initialize", "tsclustFamily",
          function(.Object, dist, allcent, ..., control = list(), fuzzy = FALSE) {
              dots <- list(...)
              dots$.Object <- .Object
              if (!missing(dist)) {
                  if (is.character(dist))
                      dots$dist <- ddist2(dist, control)
                  else
                      dots$dist <- dist
              }
              if (fuzzy) {
                  dots$cluster <- f_cluster
                  if (!missing(allcent) && is.character(allcent))
                      allcent <- match.arg(allcent, c("fcm", "fcmdd"))
              }
              if (!missing(allcent)) {
                  if (is.character(allcent)) {
                      if (allcent %in% c("pam", "fcmdd")) {
                          if (!is.null(control$distmat) && !inherits(control$distmat, "Distmat"))
                              control$distmat <- Distmat$new( # see S4-Distmat.R
                                  distmat = base::as.matrix(control$distmat)
                              )
                      }
                      dots$allcent <- all_cent2(allcent, control)
                  }
                  else if (is.function(allcent))
                      dots$allcent <- allcent
                  else
                      stop("Centroid definition must be either a function or a character")
              }
              do.call(methods::callNextMethod, dots, TRUE)
          })
########

# pr_DB
# https://github.com/asardaes/dtwclust/blob/master/R/DISTANCES-dtw-basic.R
########
#' Basic DTW distance
#'
#' This is a custom implementation of the DTW algorithm without all the functionality included in
#' [dtw::dtw()]. Because of that, it should be faster, while still supporting the most common
#' options.
#'
#' @export
#' @importFrom dtw symmetric1
#' @importFrom dtw symmetric2
#'
#' @param x,y Time series. Multivariate series must have time spanning the rows and variables
#'   spanning the columns.
#' @param window.size Size for slanted band window. `NULL` means no constraint.
#' @param norm Norm for the LCM calculation, "L1" for Manhattan or "L2" for (squared) Euclidean. See
#'   notes.
#' @param step.pattern Step pattern for DTW. Only `symmetric1` or `symmetric2` supported here. Note
#'   that these are *not* characters. See [dtw::stepPattern].
#' @param backtrack Also compute the warping path between series? See details.
#' @param normalize Should the distance be normalized? Only supported for `symmetric2`.
#' @param sqrt.dist Only relevant for `norm = "L2"`, see notes.
#' @param ... Currently ignored.
#' @template error-check
#'
#' @details
#'
#' If `backtrack` is `TRUE`, the mapping of indices between series is returned in a list.
#'
#' @template window
#'
#' @return The DTW distance. For `backtrack` `=` `TRUE`, a list with:
#'
#'   - `distance`: The DTW distance.
#'   - `index1`: `x` indices for the matched elements in the warping path.
#'   - `index2`: `y` indices for the matched elements in the warping path.
#'
#' @template proxy
#' @template symmetric
#' @section Proxy version:
#'
#'   In order for symmetry to apply here, the following must be true: no window constraint is used
#'   (`window.size` is `NULL`) or, if one is used, all series have the same length.
#'
#' @note
#'
#' The elements of the local cost matrix are calculated by using either Manhattan or squared
#' Euclidean distance. This is determined by the `norm` parameter. When the squared Euclidean
#' version is used, the square root of the resulting DTW distance is calculated at the end (as
#' defined in Ratanamahatana and Keogh 2004; Lemire 2009; see vignette references). This can be
#' avoided by passing `FALSE` in `sqrt.dist`.
#'
#' The DTW algorithm (and the functions that depend on it) might return different values in 32 bit
#' installations compared to 64 bit ones.
#'
#' An infinite distance value indicates that the constraints could not be fulfilled, probably due to
#' a too small `window.size` or a very large length difference between the series.
#'
#' @example man-examples/multivariate-dtw.R
#'
dtw_basic <- function(x, y, window.size = NULL, norm = "L1",
                      step.pattern = dtw::symmetric2, backtrack = FALSE,
                      normalize = FALSE, sqrt.dist = TRUE, ..., error.check = TRUE)
{
    if (error.check) {
        check_consistency(x, "ts")
        check_consistency(y, "ts")
    }

    if (is.null(window.size))
        window.size <- -1L
    else
        window.size <- check_consistency(window.size, "window")

    if (NCOL(x) != NCOL(y)) stop("Multivariate series must have the same number of variables.")

    if (identical(step.pattern, dtw::symmetric1))
        step.pattern <- 1
    else if (identical(step.pattern, dtw::symmetric2))
        step.pattern <- 2
    else
        stop("step.pattern must be either symmetric1 or symmetric2 (without quotes)")

    norm <- match.arg(norm, c("L1", "L2"))
    norm <- switch(norm, "L1" = 1, "L2" = 2)
    backtrack <- isTRUE(backtrack)
    normalize <- isTRUE(normalize)
    sqrt.dist <- isTRUE(sqrt.dist)
    if (normalize && step.pattern == 1) stop("Unable to normalize with chosen step pattern.")

    if (backtrack)
        gcm <- matrix(0, NROW(x) + 1L, NROW(y) + 1L)
    else
        gcm <- matrix(0, 2L, NROW(y) + 1L)

    d <- .Call(C_dtw_basic, x, y, window.size,
               NROW(x), NROW(y), NCOL(x),
               norm, step.pattern, backtrack, normalize, sqrt.dist,
               gcm, PACKAGE = "dtwclust")

    if (backtrack) {
        d$index1 <- d$index1[d$path:1L]
        d$index2 <- d$index2[d$path:1L]
        d$path <- NULL
    }
    # return
    d
}

# ==================================================================================================
# Wrapper for proxy::dist
# ==================================================================================================

#' @importFrom dtw symmetric1
#' @importFrom dtw symmetric2
#'
dtw_basic_proxy <- function(x, y = NULL, window.size = NULL, norm = "L1",
                            step.pattern = dtw::symmetric2,
                            normalize = FALSE, sqrt.dist = TRUE, ...,
                            error.check = TRUE, pairwise = FALSE)
{
    x <- tslist(x)
    if (error.check) check_consistency(x, "vltslist")
    if (is.null(y)) {
        y <- x
        symmetric <- is.null(window.size) || !different_lengths(x)
    }
    else {
        y <- tslist(y)
        if (error.check) check_consistency(y, "vltslist")
        symmetric <- FALSE
    }

    fill_type <- mat_type <- dim_names <- NULL # avoid warning about undefined globals
    eval(prepare_expr) # UTILS-expressions.R

    # adjust parameters for this distance
    if (is.null(window.size))
        window.size <- -1L
    else
        window.size <- check_consistency(window.size, "window")

    if (identical(step.pattern, dtw::symmetric1))
        step.pattern <- 1
    else if (identical(step.pattern, dtw::symmetric2))
        step.pattern <- 2
    else
        stop("step.pattern must be either symmetric1 or symmetric2 (without quotes)")

    normalize <- isTRUE(normalize)
    sqrt.dist <- isTRUE(sqrt.dist)

    if (normalize && step.pattern == 1) stop("Unable to normalize with chosen step pattern.")

    norm <- match.arg(norm, c("L1", "L2"))
    norm <- switch(norm, "L1" = 1, "L2" = 2)

    mv <- is_multivariate(c(x, y))
    backtrack <- FALSE

    # calculate distance matrix
    distance <- "DTW_BASIC" # read in C++, can't be temporary!
    distargs <- list(
        window.size = window.size,
        norm = norm,
        step.pattern = step.pattern,
        backtrack = backtrack,
        normalize = normalize,
        sqrt.dist = sqrt.dist
    )
    num_threads <- get_nthreads()
    .Call(C_distmat_loop,
          D, x, y, distance, distargs, fill_type, mat_type, num_threads,
          PACKAGE = "dtwclust")

    # adjust D's attributes
    if (pairwise) {
        dim(D) <- NULL
        class(D) <- "pairdist"
    }
    else {
        dimnames(D) <- dim_names
        class(D) <- "crossdist"
    }
    attr(D, "method") <- "DTW_BASIC"
    # return
    D
}
########

# pr_DB
# https://github.com/asardaes/dtwclust/blob/master/R/S4-Distmat.R
########
# ==================================================================================================
# Distmat RC and methods to transparently handle PAM centroids
# ==================================================================================================

#' Distance matrix
#'
#' Reference class that is used internally for cross-distance matrices.
#'
#' @importFrom methods setRefClass
#'
#' @field distmat A distance matrix.
#' @field series Time series list.
#' @field distfun The distance function to calculate the distance matrix.
#' @field dist_args Arguments for the distance function.
#' @field id_cent Indices of the centroids (if any).
#'
#' @keywords internal
#'
Distmat <- methods::setRefClass("Distmat",
                       fields = list(
                           distmat = "ANY",
                           series = "list",
                           distfun = "function",
                           dist_args = "list",
                           id_cent = "integer"
                       ),
                       methods = list(
                           initialize = function(..., distmat, series, distance, control, error.check = TRUE) {
                               "Initialization based on needed parameters"

                               if (missing(distmat)) {
                                   if (tolower(distance) == "dtw_lb") distance <- "dtw_basic"
                                   if (error.check) {
                                       check_consistency(series, "vltslist")
                                       check_consistency(distance,
                                                         "dist",
                                                         trace = FALSE,
                                                         diff_lengths = different_lengths(series),
                                                         silent = FALSE)
                                       if (class(control) != "PtCtrl")
                                           stop("Invalid control provided.") # nocov
                                   }
                                   # need another dist closure, otherwise it would be recursive
                                   control$distmat <- NULL
                                   initFields(...,
                                              series = series,
                                              distfun = ddist2(distance, control))
                               }
                               else
                                   initFields(..., distmat = distmat)
                               # return
                               invisible(NULL)
                           }
                       )
)

#' Generics for `Distmat`
#'
#' Generics with methods for [Distmat-class].
#'
#' @name Distmat-generics
#' @rdname Distmat-generics
#' @keywords internal
#' @importFrom methods setMethod
#'
NULL

#' @rdname Distmat-generics
#' @aliases [,Distmat,ANY,ANY,ANY
#'
#' @param x A [Distmat-class] object.
#' @param i Row indices.
#' @param j Column indices.
#' @param ... Ignored.
#' @param drop Logical to drop dimensions after subsetting.
#'
#' @details
#'
#' Accessing matrix elements with `[]` first calculates the values if necessary.
#'
setMethod(`[`, "Distmat", function(x, i, j, ..., drop = TRUE) {
    if (inherits(x$distmat, "uninitializedField")) {
        if (inherits(x$distfun, "uninitializedField"))
            stop("Invalid internal Distmat instance.") # nocov

        centroids <- if (identical(i,j)) NULL else x$series[j]
        dm <- quoted_call(x$distfun, x = x$series[i], centroids = centroids, dots = x$dist_args)
    }
    else {
        dm <- x$distmat[i, j, drop = drop]
        if (identical(dim(dm), dim(x$distmat))) attributes(dm) <- attributes(x$distmat)
    }
    # return
    dm
})

dim.Distmat <- function(x) { dim(x$distmat) } # nocov
########

#
#
########
########

.sbd_y = function(x) {
    shift_2 <- ReadData.list[[x]] %>% as.numeric()
    sbd <- dtwclust::SBD(shift_1,
                         shift_2, 
                         znorm = FALSE, 
                         error.check = TRUE, 
                         return.shifted = TRUE)
    return(sbd$yshift)
}

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