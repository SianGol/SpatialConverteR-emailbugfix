#' @title seuratToAnndata
#' @description
#' Provided a Seurat object, converts to an Anndata object, which is exported as a h5ad file to a provided filepath.Compatible with single-cell and spatial objects, including with and without H & E and immunofluorescent images.
#' @param seurat_obj Seurat Object.
#' @param export_path String file path to export file to.
#' @param raw_counts_names Vector c(assay_name, layer_name) for assay layer to be used as raw counts matrix 'X' by AnnData. c('RNA', 'counts') by default.
#' @param image Logical Import H & E/IHC Image. FALSE by default.
#' @return None
#'
#' @export
seuratToAnndata <- function(seurat_obj, export_path, raw_counts_names=c('RNA','counts'), image = FALSE){
  # Make sure export path ends in '.h5ad'
  if (!endsWith(export_path,".h5ad")){
    export_path <- paste(export_path,".h5ad",sep = "")
  }

  # Extract counts matrix
  counts_matrices <- seurat_assayMatrix(seurat_obj)

  # Extract matrix for X
  raw_counts_names <- paste(raw_counts_names, collapse = "_")
  X <- counts_matrices[[raw_counts_names]]
  counts_matrices[[raw_counts_names]] <- NULL

  # Create anndata object
  anndata <- anndataR::AnnData(X = Matrix::t(X), obs = seurat_obj[[]])

  # If any extra assays/layers, load as layers IF size is same.
  if (length(counts_matrices) > 0){
    for (matrix in names(counts_matrices)){
      if (any(dim(counts_matrices[[matrix]]) != dim(X))){
        warning(sprintf('Dimensions of all array matrices must be the same as raw matrix. Matrix %s will be skipped.',matrix))
        counts_matrices[[matrix]] <- NULL
      }
      else{
        # Transpose matrices
        counts_matrices[[matrix]] <- Matrix::t(counts_matrices[[matrix]])
      }
    }

    if (length(counts_matrices) > 0){
      anndata$layers <- counts_matrices
    }
  }

  # Dimension Reductions for obsm
  if (length(SeuratObject::Reductions(seurat_obj)) > 0){
    obsm <- list()
    for (reduc in SeuratObject::Reductions(seurat_obj)){
      mat <- SeuratObject::Embeddings(seurat_obj, reduc)
      colnames(mat) <- NULL
      obsm[[paste("X", reduc, sep = "_")]] <- mat
    }
  }
  else{
    obsm <- list()
  }

  # Append spatial data to obsm
  if (length(SeuratObject::Images(seurat_obj)) > 0){
    spatial_coords <- seurat_createSpatialCoords(seurat_obj)
    colnames(spatial_coords) <- NULL
    obsm[["spatial"]] <- spatial_coords
  }

  if (length(obsm) > 0){
    anndata$obsm <- obsm
  }

  # Image data (if required)
  if (image == TRUE){
    spatial <- list()
    for (img_name in names(seurat_obj@images)){
      if ('image' %in% methods::slotNames(seurat_obj@images[[img_name]])){
        if ('scale.factors' %in% methods::slotNames(seurat_obj@images[[img_name]])){
          img <- seurat_obj@images[[img_name]]@image
          if (max(img) != 1){
            img <- img/max(img)
          }
          spatial[[img_name]] <- list()
          spatial[[img_name]][['images']] <- list()
          spatial[[img_name]][['images']][['lowres']] <- img
          spatial[[img_name]][['scalefactors']] <- seurat_createScalefactors_ad(seurat_obj, img_name)
        }
      }
    }
    anndata$uns <- list(spatial = spatial)
  }

  anndataR::write_h5ad(anndata, export_path)
}
