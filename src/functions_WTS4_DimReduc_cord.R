# library
##################################################
library(Rtsne)
library(uwot)
##################################################
set.seed(1234)

#### TSNE#### 
.wts_tsne = function(x) {
    dist_data <- x
    set.seed(1234)
    tSNE <- Rtsne(dist_data, 
                  is_distance = TRUE, 
                  dims = 2, 
                  perplexity = 15, 
                  verbose = TRUE, 
                  max_iter = 1000)
    return_object <- data.frame(cord_1 = tSNE$Y[,1],
                                cord_2 = tSNE$Y[,2],
                                cell_type = attr(dist_data, "Labels")
                                )
    return(return_object)
}
#### UMAP####
.wts_umap = function(x) {
    dist_data <- x
    # attr(dist_data, "Labels") %>% length() -> lab_length
    length(attr(dist_data, "Labels")) -> lab_length
    set.seed(1234)
    umap_d <- uwot::umap(dist_data,
                         metric = "precomputed", 
                         nn_method = uwot:::dist_nn(dist_data, k = lab_length),
                         n_neighbors = 15,
                         n_components = 2
                         )
    return_object <- data.frame(cord_1 = umap_d[,1],
                                cord_2 = umap_d[,2],
                                cell_type = attr(dist_data, "Labels")
                                )
    return(return_object)
}