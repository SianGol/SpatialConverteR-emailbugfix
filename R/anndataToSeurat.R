#' @title anndataToSeurat
#' @description
#' Reads in an AnnData H5AD object from a provided file path, and converts to a Seurat Object. Compatible with single-cell and spatial objects, including with and without H & E and immunofluorescent images.
#' @param adata_path String file path to AnnData H5AD Object.
#' @param image Logical Import H & E/IHC Image. FALSE by default.
#' @importClassesFrom Seurat VisiumV2
#' @return Seurat Object

#' @export
anndataToSeurat <- function(adata_path, image = FALSE){
  ## Read in H5AD file
  file <- rhdf5::H5Fopen(adata_path)
  on.exit(rhdf5::h5closeAll())

  ## Counts Matrix into Seurat Object
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = h5_createCountsMatrix(file))

  ## Load metadata
  seurat_obj <- SeuratObject::AddMetaData(seurat_obj, h5_createMetadata(file))

  ## Extra layers
  layer_contents <- rhdf5::h5ls(file, recursive = 2)
  layer_contents <- layer_contents[layer_contents$group == "/layers","name"]
  if (length(layer_contents) > 0){
    for (layer in layer_contents){
      seurat_obj[['RNA']][layer] <- h5_createAssayMatrix(file, layer)
    }
  }

  ## Dimension Reductions (if any)
  reduc_contents <- rhdf5::h5ls(file, recursive = 2)
  reduc_contents <- reduc_contents[reduc_contents$group == "/obsm","name"]
  for (reduc_name in reduc_contents){
  if (grepl("X_",reduc_name)){
    seurat_obj[[substring(reduc_name, 3)]] <- SeuratObject::CreateDimReducObject(embeddings = h5_createDimReduc(file,reduc_name, SeuratObject::Cells(seurat_obj)), key = paste(substring(reduc_name, 3), '_', sep = ""), assay = SeuratObject::DefaultAssay(seurat_obj))
  }
  }

  # Spatial Coordinates (if any)
  if ("spatial" %in% reduc_contents){
    centroids <- h5_createCentroids_Seurat(file, SeuratObject::Cells(seurat_obj))

    # Images (if any)
    if (image == TRUE){
      uns_contents <- rhdf5::h5ls(file, recursive = 2)
      uns_contents <- uns_contents[uns_contents$group == "/uns", "name"]
      if ("spatial" %in% uns_contents){
        for (img_name in names(file$'uns/spatial')){
        if ("lowres" %in% names(file$'uns/spatial'[[img_name]]$images)){
          if ("scalefactors" %in% names(file$'/uns/spatial'[[img_name]])){
            img <- aperm(file$'/uns/spatial'[[img_name]]$images$lowres, c(2,3,1))
            storage.mode(img) <- "double"
            if (max(img) != 1){
              img <- img/max(img)
            }
            #
            spatial = methods::new("VisiumV2", image = img, scale.factors = h5_createScalefactors_Seurat(file, img_name), boundaries = list('centroids' = centroids), assay = SeuratObject::DefaultAssay(seurat_obj), coords_x_orientation = "vertical")
          }
          else{
            spatial <- SeuratObject::CreateFOV(centroids, assay = SeuratObject::DefaultAssay(seurat_obj))
          }
          }
        else if ("hires" %in% names(file$'uns/spatial'[[img_name]]$images)){
          warning("Low resolution image missing. Using high resolution image instead.")
          if ('scalefactors' %in% names(file$'/uns/spatial'[[img_name]])){
            img <- aperm(file$'/uns/spatial'[[img_name]]$images$hires, c(2,3,1))
            storage.mode(img) <- "numeric"
            if (max(img) != 1){
              img <- img/max(img)
            }
            spatial = methods::new("VisiumV2", image = img, scale.factors = h5_createScalefactors_Seurat(file, img_name), boundaries = list(centroids = centroids), assay = SeuratObject::DefaultAssay(seurat_obj), coords_x_orientation = "vertical")
          }
          else{
            spatial <- SeuratObject::CreateFOV(centroids, assay = SeuratObject::DefaultAssay(seurat_obj))
          }
      }
        seurat_obj[[img_name]] <- spatial
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

  rhdf5::h5closeAll()
  return(seurat_obj)
}
