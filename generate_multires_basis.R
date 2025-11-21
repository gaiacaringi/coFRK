build_BAUs_basis <- function(
    spatial_domain, 
    nres ,
    cellsize ,
    grid_BAUs = NULL
) {
  library(sp)
  library(FRK)
  library(Matrix)
  library(data.table)

  # BAUs
  BAUs <- auto_BAUs(manifold = plane(), data = spatial_domain,cellsize=cellsize)
  BAUs$fs <- 1
  #BAUs <- grid_BAUs
  nres <- nres
  # Basis functions
  G_all <- auto_basis(
    manifold = plane(),
    data = spatial_domain,
    regular = 1,
    nres = nres,
    #type     = "bisquare",
    #subsamp  = 5000,                 
    prune    = 0,                  
    #max_basis = 400,
    #verbose = 1L
  )
  tbl <- as.data.frame(G_all@df)
  
  #Computing dimensions of the rectangluar grid of the basis centroids per resolutional level
  tbl_splitted <- split(tbl, tbl$res)
  get_dims <- function(d, digits = 6) {
    nb   <- nrow(d)
    cnt  <- table(round(as.numeric(d$loc2), digits))  # counts per identical loc2
    n_y  <- max(as.integer(cnt))
    n_x  <- as.integer(nb / n_y)
    c(n_x = n_x, n_y = n_y, n_basis = nb)
  }
  dims_mat <- t(sapply(tbl_splitted, get_dims))
  dims_basis_lattice <- data.frame(res = as.integer(rownames(dims_mat)), dims_mat, row.names = NULL)
  
  
  # Evaluate basis functions at BAU centroids
  BAU_centroids <- coordinates(BAUs)
  n_BAUs <- nrow(BAUs)
  n_bases <- length(G_all@fn)
  
  S0 <- Matrix(0, nrow = n_BAUs, ncol = n_bases, sparse = TRUE)
  for (j in seq_len(n_bases)) {
    S0[, j] <- G_all@fn[[j]](BAU_centroids)
  }
  
  # Grouping per resolution level
  resolution_levels <- G_all@df$res
  unique_levels <- sort(unique(resolution_levels))
  S0_list <- vector("list", length(unique_levels))
  names(S0_list) <- paste0("level_", unique_levels)
  
  for (i in seq_along(unique_levels)) {
    lvl <- unique_levels[i]
    col_idx <- which(resolution_levels == lvl)
    S0_list[[i]] <- S0[, col_idx, drop = FALSE]
  }
  
  return(list(
    BAUs = BAUs,
    S0 = S0,
    S0_list = S0_list,
    basis_object = G_all,
    dims_basis_lattice = dims_basis_lattice
  ))
  
}

