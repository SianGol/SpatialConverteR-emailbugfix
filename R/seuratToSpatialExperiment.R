#' @title seuratToSpatialExperiment
#' @description
#' Provided a Seurat object, converts to a SpatialExperiment/SingleCellExperiment Object. Compatible with single-cell and spatial objects, including with and without H & E and immunofluorescent images.
#' @param seurat_obj Seurat Object.
#' @param image Logical Import H & E/IHC Image. FALSE by default.
#' @return SingleCellExperiment or SpatialExperiment Object
#'
#' @export
seuratToSpatialExperiment <- function(seurat_obj, image = FALSE){

  assay_list <- seurat_assayMatrix(seurat_obj)


  # Create SingleCellExperiment object, with counts matrix, and metadata
  spe <- SingleCellExperiment::SingleCellExperiment(assay_list, colData = seurat_obj[[]])

  # Dimension Reductions
  if (length(SeuratObject::Reductions(seurat_obj)) > 0){
    red_list <- list()
    for (reduc in SeuratObject::Reductions(seurat_obj)){
      red_list[[reduc]] <- SeuratObject::Embeddings(seurat_obj, reduc)
    }
    SingleCellExperiment::reducedDims(spe) <- red_list
  }

  # Spatial Information
  if (length(SeuratObject::Images(seurat_obj)) > 0){ #check for spatial info
    spe_orig <- spe
    spe <- SpatialExperiment::SpatialExperiment(SummarizedExperiment::assays(spe_orig), colData = SingleCellExperiment::colData(spe_orig), rowData = SingleCellExperiment::rowData(spe_orig), reducedDims = SingleCellExperiment::reducedDims(spe_orig)) # convert to SpatialExperiment Object
    rm(spe_orig)
    SpatialExperiment::spatialCoords(spe) <- seurat_createSpatialCoords(seurat_obj)

    # Image Import (if any)
    if (image == TRUE){
      img_df_total <- S4Vectors::DataFrame()
      for (img_name in names(seurat_obj@images)){
        if ('image' %in% methods::slotNames(seurat_obj@images[[img_name]])){
          if ('scale.factors' %in% methods::slotNames(seurat_obj@images[[img_name]])){
            img_raw <- aperm(seurat_obj@images[[img_name]]@image, c(2,1,3))
            if (max(img_raw) != 1){
              img_raw <- img_raw/max(img_raw)
            }
            img <- SpatialExperiment::SpatialImage(x = grDevices::as.raster(img_raw))
            imgdf <- S4Vectors::DataFrame(sample_id = 'sample01', image_id = img_name, data = I(list(img)), scaleFactor = seurat_createScalefactors_SPE(seurat_obj, img_name))
            img_df_total <- rbind(img_df_total, imgdf)
          }

        }
      }
      if (nrow(img_df_total) > 0){
        SpatialExperiment::imgData(spe) <- img_df_total
        for (i in seq_len(nrow(SpatialExperiment::imgData(spe)))) {
          spe <- SpatialExperiment::rotateImg(spe,
                           sample_id = SpatialExperiment::imgData(spe)$sample_id[i],
                           image_id = SpatialExperiment::imgData(spe)$image_id[i],
                           degrees = 90)
          spe <- SpatialExperiment::mirrorImg(spe, sample_id = SpatialExperiment::imgData(spe)$sample_id[i],
                                                 image_id = SpatialExperiment::imgData(spe)$image_id[i],
                                                 axis = 'v')
        }
      }
    }
  }


  return(spe)
}

