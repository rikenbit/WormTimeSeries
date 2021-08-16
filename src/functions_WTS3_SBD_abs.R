#library
##################################################
# library(dtwclust)
# RNGkind(kind = "Mersenne-Twister")
library(tidyverse)
library(openxlsx)
##################################################

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

# setting Random number generator
RNGkind(kind = "Mersenne-Twister")
dtwclust_rngkind <- "Mersenne-Twister"
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

# 
# 
########
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