#' @title spatialExperimentToSeurat
#' @description
#' Provided a SingleCellExperiment/SpatialExperiment Object, converts it to a Seurat Object. Compatible with single-cell and spatial objects, including with and without H & E and immunofluorescent images.
#' @param spe SingleCellExperiment/SpatialExperiment Object
#' @param raw_counts_name String Name of assay with unnormalised counts data, to be provided as assay counts. "counts" by default.
#' @param image Logical Import H & E/IHC Image. FALSE by default.
#' @importFrom grDevices col2rgb
#' @return Seurat Object.
#'
#' @export
spatialExperimentToSeurat <- function(spe, raw_counts_name = "counts", image = FALSE){
  # Create Seurat Object with counts assay
  seurat_obj <- SeuratObject::CreateSeuratObject(spe@assays@data[[raw_counts_name]])

  # Add metadata
  seurat_obj <- SeuratObject::AddMetaData(seurat_obj, as.data.frame(spe@colData))

  # Extra layers
  layers <- spe@assays@data
  layers[[raw_counts_name]] <- NULL
  if (length(layers) > 0){
    for (layer in names(layers)){
      seurat_obj[['RNA']][layer] <- spe@assays@data[[layer]]
    }
  }

  # Dimension Reductions
  reductions <- SingleCellExperiment::reducedDimNames(spe)
  if (length(reductions) > 0){
    for (reduc in reductions){
      seurat_obj[[reduc]] <- SeuratObject::CreateDimReducObject(embeddings = SingleCellExperiment::reducedDim(spe, reduc), key = paste(reduc,"_", sep = ""), assay = SeuratObject::DefaultAssay(seurat_obj))
    }
  }

  # Spatial Coordinates (if any)
  if (inherits(spe,"SpatialExperiment")){
    centroids <- SeuratObject::CreateCentroids(SpatialExperiment::spatialCoords(spe))

    # Images (if required)
    if (image == TRUE){
      # Check if imgData has any info
      if (length(SpatialExperiment::imgData(spe)) > 0){
        for (i in seq_len(nrow(SpatialExperiment::imgData(spe)))){
          img_raw <- as.array(SpatialExperiment::imgRaster(SpatialExperiment::imgData(spe)[[i,'data']]))
          img <- array(t(col2rgb(img_raw)), dim = c(nrow(img_raw),ncol(img_raw),3))/max(col2rgb(img_raw))
          spatial <-  methods::new("VisiumV2", image = img, scale.factors = Seurat::scalefactors(spot = 1, fiducial = 1, hires = 1, lowres = SpatialExperiment::imgData(spe)[[1,'scaleFactor']]), boundaries = list(centroids = centroids), assay = SeuratObject::DefaultAssay(seurat_obj), coords_x_orientation = "vertical")
          if (all(SpatialExperiment::imgData(spe)[,'sample_id'] == 'sample01')){
            seurat_obj[[SpatialExperiment::imgData(spe)[[i,'image_id']]]] <- spatial
          }
          else{
            seurat_obj[[paste(SpatialExperiment::imgData(spe)[[i,'sample_id']],SpatialExperiment::imgData(spe)[[i,'image_id']],sep = "_")]] <- spatial
          }
        }
      }
      else{
        spatial <- SeuratObject::CreateFOV(centroids, assay = SeuratObject::DefaultAssay(seurat_obj))
        seurat_obj[["spatial"]] <- spatial
      }

    }
    else{
      spatial <- SeuratObject::CreateFOV(centroids, assay = SeuratObject::DefaultAssay(seurat_obj))
      seurat_obj[["spatial"]] <- spatial
    }
  }
  seurat_obj <- Seurat::UpdateSeuratObject(seurat_obj)
  return(seurat_obj)
}
