#' @title flipImageSeurat
#' @description
#' When importing immunofluorescent images from ScanPy into Seurat, the images may require flipping, which can be done using this function.
#' If needed for SpatialExperiment objects, use existing SpatialExperiment::mirrorObject function.
#' @param seurat_obj Seurat Object with immunofluorescent images.
#' @return Seurat Object with flipped images.

#' @export
flipImageSeurat <- function(seurat_obj){

  # Check if SpatialImage slot exists in Seurat Object
  if (length(SeuratObject::Images(seurat_obj)) == 0){
    stop("Error: Seurat Object does not contain any SpatialImage objects.")
  }

  # Check if each SpatialImage obj has an "image" slot (location of image array)
  for (image in SeuratObject::Images(seurat_obj)){
    if ("image" %in% methods::slotNames(seurat_obj@images[[image]])){
      seurat_obj@images[[image]]@image <- flipFunc(seurat_obj@images[[image]]@image)
    }
    else{
      warning(paste(image, "has no image data."))
    }
  }
  return(seurat_obj)
}
