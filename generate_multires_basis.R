generate_multires_basis <- function(
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
    #subsamp  = 5000,                  # limit points used during placement
    prune    = 0,                    # drop basis with little support from data
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
  
  # Raggruppamento corretto per livelli di risoluzione
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

generate_multires_basis <- function(
    spatial_domain,
    scales = c(7000,3000),
    cellsize = 2000,
    crs_proj = CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs")
) {
  library(sp)
  library(FRK)
  library(Matrix)
  
  # 1. Proiezione (lo faccio fuori)
  #spatial_domain <- spTransform(spatial_domain, crs_proj)
  
  # 2. BAUs
  BAUs <- auto_BAUs(manifold = plane(), data = spatial_domain, cellsize = cellsize)
  proj4string(BAUs) <- crs_proj
  BAUs$fs <- 1
  
  # 3. Costruzione delle basi con nres > 1
  G_all <- auto_basis(
    manifold = plane(),
    data = list(spatial_domain),
    nres = length(scales),
    type = "bisquare",
    scale = scales,
    regular = 1
  )
  
  # 4. Valutazione manuale su centroidi
  BAU_centroids <- coordinates(BAUs)
  n_BAUs <- nrow(BAUs)
  n_bases <- length(G_all@fn)
  
  S0 <- Matrix(0, nrow = n_BAUs, ncol = n_bases, sparse = TRUE)
  for (j in seq_len(n_bases)) {
    S0[, j] <- G_all@fn[[j]](BAU_centroids)
  }
  
  # 5. Calcolo split automatico (assumendo regular griglia → divisione equa)
  nres <- length(scales)
  nb_per_level <- as.integer(n_bases / nres)
  S0_list <- vector("list", nres)
  col_start <- 1
  
  for (i in seq_len(nres)) {
    col_end <- if (i < nres) col_start + nb_per_level - 1 else n_bases
    S0_list[[i]] <- S0[, col_start:col_end, drop = FALSE]
    col_start <- col_end + 1
  }
  
  return(list(
    BAUs = BAUs,
    S0 = S0,
    S0_list = S0_list,
    basis_object = G_all
  ))
}




#Funzione fatta usando LatticeKrig -> non funziona questo approccio
generate_multires_basis2 <- function(
    spatial_domain, 
    delta_level1 = 10,      # Distanza nodi della griglia più grezza     
    nlevel = 2,              # Numero di livelli di risoluzione
    overlap = 2.5,           # Fattore di overlap (quanto le basi si sovrappongono)
    buffer = 0.01,            # Buffer extra fuori dal dominio
    cellsize = 0.05
) {
  library(sp)
  library(LatticeKrig)
  library(Matrix)
  
  # BAUs
  BAUs <- auto_BAUs(manifold = plane(), data = spatial_domain, cellsize = 2)
  proj4string(BAUs) <- crs_utm
  
  #BAUs <- auto_BAUs(manifold = plane(), data = spatial_domain, cellsize = 0.1)
  BAUs$fs <- 1
  
  BAU_centroids <- SpatialPoints(coordinates(BAUs), proj4string = CRS(proj4string(BAUs)))
  
  bbox <- bbox(BAU_centroids)
  x_seq <- seq(bbox[1,1] - buffer, bbox[1,2] + buffer, by = cellsize)
  y_seq <- seq(bbox[2,1] - buffer, bbox[2,2] + buffer, by = cellsize)
  grid_pts <- expand.grid(x = x_seq, y = y_seq)
  coordinates(grid_pts) <- ~x + y
  proj4string(grid_pts) <- proj4string(BAU_centroids)
  
  
  # Struttura multiresoluzionale con LKrig (Kleiber)
  nx <- ceiling((max(x_seq) - min(x_seq)) / delta_level1)
  ny <- ceiling((max(y_seq) - min(y_seq)) / delta_level1)
  NC <- c(nx, ny)
  
  LKinfo <- LKrigSetup(
    x = coordinates(grid_pts),
    NC = NC,
    nlevel = nlevel,
    delta = rep(delta_level1, nlevel),
    a.wght = as.list(rep(4.5, nlevel)),     
    alpha = rep(1, nlevel),       
    normalize = TRUE,
    overlap = overlap,
    edge = FALSE
  )
  
  # Valuto le funzioni base sui centroidi delle BAUs
  S0 <- LKrig.basis(coordinates(BAU_centroids), LKinfo)
  
  return(list(
    BAUs = BAUs,
    S0 = S0,
    LKinfo = LKinfo
  ))
  
}



