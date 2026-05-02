#' @title spatialExperimentToAnndata
#' @description
#' Provided a SingleCellExperiment/SpatialExperiment Object, converts to an AnnData object, which is exported as a h5ad file to a provided filepath. Compatible with single-cell and spatial objects, including with and without H & E and immunofluorescent images.
#' @param spe SingleCellExperiment/SpatialExperiment Object
#' @param export_path String file path to export file to.
#' @param raw_counts_name String Name of assay with unnormalised counts data, to be used as raw counts matrix 'X' by AnnData. "counts" by default.
#' @param image Logical Import H & E/IHC Image. FALSE by default.
#' @importFrom grDevices col2rgb
#' @return None

spatialExperimentToAnndata <- function(spe, export_path, raw_counts_name = 'counts', image = FALSE){
  # Make sure export path ends in '.h5ad'
  if (!endsWith(export_path,".h5ad")){
    export_path <- paste(export_path,".h5ad",sep = "")
  }

  # Add counts matrix X, cell (obs) and gene (var) metadata
  anndata <- anndataR::AnnData(X = Matrix::t(spe@assays@data[[raw_counts_name]]), obs = as.data.frame(spe@colData), var = as.data.frame(SingleCellExperiment::rowData(spe)))

  # Add extra assay layers
  layers <- spe@assays@data
  # Remove counts matrix given as 'X'
  layers[[raw_counts_name]] <- NULL
  if (length(layers) > 0){
    for (layer in names(layers)){
      # transpose matrices
      layers[[layer]] <- Matrix::t(layers[[layer]])
    }
    anndata$layers <- as.list(layers)
  }

  # Dimension Reductions
  reductions <- SingleCellExperiment::reducedDimNames(spe)
  if (length(reductions) > 0){
    obsm <- list()
    for (reduc in reductions){
      mat <- SingleCellExperiment::reducedDim(spe, reduc)
      colnames(mat) <- NULL
      obsm[[paste("X", reduc, sep = "_")]] <- mat
    }
  }
  else{
    obsm <- list()
  }

  # Append spatial data to obsm
  if (inherits(spe,"SpatialExperiment")){
    spatial_coords <- SpatialExperiment::spatialCoords(spe)
    colnames(spatial_coords) <- NULL
    obsm[["spatial"]] <- spatial_coords
  }

  # Add obsm to anndata obj
  if (length(obsm) > 0){
    anndata$obsm <- obsm
  }

  # Image data (if required)
  if (image == TRUE){
    if (length(SpatialExperiment::imgData(spe)) > 0){
      spatial <- list()
      for (i in seq_len(nrow(SpatialExperiment::imgData(spe)))){
        img_raw <- as.array(SpatialExperiment::imgRaster(SpatialExperiment::imgData(spe)[[i,'data']]))
        img <- aperm(array(t(col2rgb(img_raw)), dim = c(ncol(img_raw),nrow(img_raw),3))/max(col2rgb(img_raw)),c(2,1,3))
        if (all(SpatialExperiment::imgData(spe)[,'sample_id'] == 'sample01')){
          img_name <- SpatialExperiment::imgData(spe)[[i,'image_id']]
        }
        else{
          img_name <- paste(SpatialExperiment::imgData(spe)[[i,'sample_id']],SpatialExperiment::imgData(spe)[[i,'image_id']],sep = "_")
        }


        spatial[[img_name]] <- list()
        spatial[[img_name]][['images']] <- list()
        spatial[[img_name]][['images']][['lowres']] <- img
        sf_orig <-
        spatial[[img_name]][['scalefactors']] <- list(spot_diameter_fullres = 1, fiducial_diameter_fullres = 1, tissue_hires_scalef = 1, tissue_lowres_scalef = SpatialExperiment::imgData(spe)[[i,'scaleFactor']])
      }
      anndata$uns <- list(spatial = spatial)
    }
  }
  anndataR::write_h5ad(anndata, export_path)
}
