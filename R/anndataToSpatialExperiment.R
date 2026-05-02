#' @title anndataToSpatialExperiment
#' @description
#' Reads in an AnnData H5AD object from a provided file path, and converts to a SpatialExperiment/SingleCellExperiment Object. Compatible with single-cell and spatial objects, including with and without H & E and immunofluorescent images.
#' @param adata_path String file path to AnnData H5AD Object.
#' @param image Logical Import H & E/IHC Image. FALSE by default.
#' @return SingleCellExperiment or SpatialExperiment Object

#' @export
anndataToSpatialExperiment <- function(adata_path, image = FALSE){
  ## Read in H5AD file
  file <- rhdf5::H5Fopen(adata_path)
  on.exit(rhdf5::h5closeAll())

  # SingleCellExperiment Object from counts matrix and metadata
  spe <- SingleCellExperiment::SingleCellExperiment(list(counts=h5_createCountsMatrix(file)), colData = h5_createMetadata(file), rowData = h5_createRowData(file))

  # Extra layers

  layer_contents <- rhdf5::h5ls(file, recursive = 2)
  layer_contents <- layer_contents[layer_contents$group == '/layers','name']
  if (length(layer_contents) > 0){
    for (layer in layer_contents){
      SummarizedExperiment::assay(spe, layer) <- h5_createAssayMatrix(file, layer)
    }
  }


  ## Dimension Reductions (if any)
  red_list <- list()

  reduc_contents <- rhdf5::h5ls(file, recursive = 2)
  reduc_contents <- reduc_contents[reduc_contents$group == '/obsm','name']
  for (reduc_name in reduc_contents){
    if (grepl("X_", reduc_name)){
      red_list[[substring(reduc_name, 3)]] <- h5_createDimReduc(file, reduc_name, colnames(spe))
    }
  }

  if (length(red_list) > 0){
    SingleCellExperiment::reducedDims(spe) <- red_list
  }

  # Spatial Information
  if ("spatial" %in% reduc_contents){ # Check if Spatial Data exists
    spe_orig <- spe
    spe <- SpatialExperiment::SpatialExperiment(SummarizedExperiment::assays(spe_orig), colData = SingleCellExperiment::colData(spe_orig), rowData = SingleCellExperiment::rowData(spe_orig), reducedDims = SingleCellExperiment::reducedDims(spe_orig)) # convert to SpatialExperiment Object
    rm(spe_orig)
    SpatialExperiment::spatialCoords(spe) <- h5_createCentroids_SPE(file, colnames(spe))

    # Images (if any)
    if (image == TRUE){
      uns_contents <- rhdf5::h5ls(file, recursive = 2)
      uns_contents <- uns_contents[uns_contents$group == '/uns', 'name']
    if ("spatial" %in% uns_contents){
      img_df_total <- S4Vectors::DataFrame()
      for (img_name in names(file$'/uns/spatial')){
        if ("lowres" %in% names(file$'/uns/spatial'[[img_name]][['images']])){
          if ("scalefactors" %in% names(file$'/uns/spatial'[[img_name]])){
            img_raw <- aperm(file$'/uns/spatial'[[img_name]]$images$lowres, c(3,2,1))
            storage.mode(img_raw) <- "double"
            if (max(img_raw) != 1){
              img_raw <- img_raw/max(img_raw)
            }
            img <- SpatialExperiment::SpatialImage(x = grDevices::as.raster(img_raw))
            imgdf <- S4Vectors::DataFrame(sample_id = 'sample01', image_id = img_name, data = I(list(img)), scaleFactor = h5_createScalefactors_SPE(file, img_name, 'lowres'))
            img_df_total <- rbind(img_df_total, imgdf)
            }
        }
        else if ("hires" %in% names(file$'/uns/spatial'[[img_name]][['images']])){
          if ('scalefactors' %in% names(file$'/uns/spatial'[[img_name]])){
            img_raw <- aperm(file$'/uns/spatial'[[img_name]]$images$hires, c(3,2,1))
            storage.mode(img_raw) <- "double"
            if (max(img_raw) != 1){
              img_raw <- img_raw/max(img_raw)
            }
            img <- SpatialExperiment::SpatialImage(x = grDevices::as.raster(img_raw))
            imgdf <- S4Vectors::DataFrame(sample_id = 'sample01', image_id = img_name, data = I(list(img)), scaleFactor = h5_createScalefactors_SPE(file, img_name, 'hires'))
            img_df_total <- rbind(img_df_total, imgdf)
            }
        }

      }
      if (nrow(img_df_total) > 0){
        SpatialExperiment::imgData(spe) <- img_df_total
      }
    }
    }
  }

  rhdf5::h5closeAll()
  return(spe)
}

