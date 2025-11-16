library(sp)
library(Matrix)

# Function that maps new data to BAUs
# Inputs: 1) data_sp -> SpatialPoints/ SpatialDataFrame 
#         2) BAUs -> SpatialPolygons
# Output: - if points -> SptialPointsDataFrame with a BAU_name columns
#         - if polygons -> same polygons back
        
map_data_to_BAUs <- function(data_sp, BAUs) {
  if (is.null(row.names(BAUs))) {
    rownames(BAUs) <- as.character(seq_len(length(BAUs)))
  }
  
  if (inherits(data_sp, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    hit <- sp::over(data_sp, BAUs)                       # data.frame with BAU attributes
    BAU_name <- rownames(BAUs)[sp::over(data_sp, BAUs)]  # map each point to BAU rowname
    keep <- !is.na(BAU_name)

    if (inherits(data_sp, "SpatialPointsDataFrame")) {
      df <- data_sp@data[keep, , drop = FALSE]
      df$BAU_name <- BAU_name[keep]
    } else {
      df <- data.frame(BAU_name = BAU_name[keep], stringsAsFactors = FALSE)
    }
    
    return(sp::SpatialPointsDataFrame(coords      = sp::coordinates(data_sp)[keep, , drop = FALSE],
                                      data        = df,
                                      proj4string = data_sp@proj4string))
  }

  if (inherits(data_sp, c("SpatialPolygons", "SpatialPolygonsDataFrame"))) {
    return(data_sp)
  }
  
  stop("map_data_to_BAUs: only SpatialPoints* or SpatialPolygons* are supported.")
}

# Function that creates the triples (i,j,x) describing which prediciton regions (rows) overlap which BAUs (columns)
# Inputs: - if points -> data: SpatialPointsDataFrame with column BAU_name 
#         - if polygons -> a SpatialPolygons* object with prediction regions
# Outputs: A list with three vectors:
# - i_idx: row indices (prediction regions).
# - j_idx: column indices (BAUs).
# - x_idx: weights (currently just 1; row normalization happens later).

buildC <- function(data, BAUs) {

  if (inherits(data, "SpatialPointsDataFrame") && "BAU_name" %in% names(data)) {
    BAU_index <- data.frame(row.names = row.names(BAUs),
                            col = seq_len(length(BAUs)))
    i_idx <- seq_len(length(data))
    j_idx <- as.integer(BAU_index[data$BAU_name, 1])
    return(list(i_idx = i_idx, j_idx = j_idx, x_idx = rep(1, length(i_idx))))
  }

  if (inherits(data, c("SpatialPolygons", "SpatialPolygonsDataFrame"))) {
    centroids <- sp::SpatialPoints(sp::coordinates(BAUs),
                                   proj4string = BAUs@proj4string)
    
    i_idx <- integer(0)
    j_idx <- integer(0)
    x_idx <- numeric(0)
    
    for (i in seq_len(length(data))) {
      poly_i <- sp::SpatialPolygons(
        list(data@polygons[[i]]),
        1L,
        proj4string = sp::CRS(sp::proj4string(data))
      )
      inside <- which(sp::over(centroids, poly_i) == 1)   # BAU centr. dentro poligono
      
      if (length(inside) > 0) {
        i_idx <- c(i_idx, rep(i, length(inside)))
        j_idx <- c(j_idx, inside)
        x_idx <- c(x_idx, rep(1, length(inside)))  # peso=1 (normalizzato dopo)
      }
    }
    return(list(i_idx = i_idx, j_idx = j_idx, x_idx = x_idx))
  }
  

  stop("buildC: only SpatialPointsDataFrame(with BAU_name) or SpatialPolygons* are supported.")
}


make_CP <- function(newdata, BAUs, normalize = TRUE) {
  newdata2 <- map_data_to_BAUs(newdata, BAUs)  # adds BAU_name for points; passthrough for polygons
  C_idx <- buildC(newdata2, BAUs)
  
  CP <- Matrix::sparseMatrix(i = C_idx$i_idx,
                             j = C_idx$j_idx,
                             x = C_idx$x_idx,
                             dims = c(length(newdata), length(BAUs)))
  
  if (normalize) {
    rs <- Matrix::rowSums(CP)
    nonempty <- rs > 0
    CP[nonempty, ] <- CP[nonempty, ] / rs[nonempty]
  }
  
  CP
}

buildC <- function(data, BAUs) {
  baus_sf <- sf::st_as_sf(BAUs) |> sf::st_make_valid()
  data_sf <- sf::st_as_sf(data) |> sf::st_make_valid()
  sf::sf_use_s2(FALSE)  # projected CRS recommended
  
  hits <- sf::st_intersects(data_sf, baus_sf)  # spatial index used
  i_idx <- rep(seq_len(nrow(data_sf)), lengths(hits))
  j_idx <- unlist(hits)
  list(i_idx = i_idx, j_idx = j_idx, x_idx = rep(1, length(i_idx)))
}

