#' Write in return properties for all functions
#' Internal helper to extract counts matrix from h5ad
#'
#' @param h5ad_file Loaded H5File using
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom methods as
#' @noRd
h5_createCountsMatrix <- function(h5ad_file){
  # Pull out counts matrix
  data <- h5ad_file$'X/data'

  n_rows <- length(h5ad_file$'obs/_index')
  n_cols <- length(h5ad_file$'var/_index')

  X <- methods::new("dgRMatrix",
                    x = as.double(data),
                    j = as.vector(h5ad_file$"X/indices"),
                    p = as.vector(h5ad_file$"X/indptr"),
                    Dim = c(n_rows, n_cols))
  X <- as(Matrix::t(X), 'CsparseMatrix')
  rownames(X) <- as.vector(h5ad_file$'var/_index')
  colnames(X) <- as.vector(h5ad_file$'obs/_index')

  return(X)
}

#' Internal helper to extract other layer matrices from h5ad
#' k
#' @param h5ad_file Loaded H5File using rhdf5
#' @param assay_name Name of assay to extract
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom methods as
#' @noRd
h5_createAssayMatrix <- function(h5ad_file, assay_name){
  # Pull out counts matrix
  data <- h5ad_file$layers[[assay_name]]$data

  n_rows <- length(h5ad_file$'obs/_index')
  n_cols <- length(h5ad_file$'var/_index')

  X <- methods::new("dgRMatrix",
                    x = as.double(data),
                    j = as.vector(h5ad_file$layers[[assay_name]]$indices),
                    p = as.vector(h5ad_file$layers[[assay_name]]$indptr),
                    Dim = c(n_rows, n_cols))
  X <- as(Matrix::t(X), 'CsparseMatrix')
  rownames(X) <- as.vector(h5ad_file$'var/_index')
  colnames(X) <- as.vector(h5ad_file$'obs/_index')

  return(X)
}

#' Internal helper to extract metadata from h5ad
#'
#' @param h5ad_file Loaded H5File using rhdf5
#' @noRd
h5_createMetadata <- function(h5ad_file){
  metadata_contents <- rhdf5::h5ls(h5ad_file, recursive = 2)
  metadata_items <- metadata_contents[metadata_contents$group == "/obs", "name"]

  metadata <- lapply(metadata_items, function(nm){

    item_type <- metadata_contents[metadata_contents$group == '/obs' & metadata_contents$name == nm, 'otype']
    if (!(any(item_type == "H5I_GROUP"))){
      return(as.vector(h5ad_file$obs[[nm]]))
    }
    else{
      cats  <- as.vector(h5ad_file$obs[[nm]][["categories"]])
      codes <- as.vector(h5ad_file$obs[[nm]][["codes"]])

      # fix string encoding
      if (is.list(cats)) {cats <- unlist(cats)}
      if (is.character(cats)) {Encoding(cats) <- "UTF-8"}

      # convert codes → levels
      return(factor(cats[codes + 1], levels = cats))
    }
  })

  metadata <- as.data.frame(metadata)
  colnames(metadata) <- metadata_items
  rownames(metadata) <- metadata$'_index'
  metadata$'_index' <- NULL

  return(metadata)
}

#' Internal helper to extract feature metadata from h5ad
#'
#' @param h5ad_file Loaded H5File using rhdf5
#' @noRd
h5_createRowData <- function(h5ad_file){
  var_contents <- rhdf5::h5ls(h5ad_file, recursive = 2)
  var_items <- var_contents[var_contents$group == "/var",'name']

  vardata <- lapply(var_items, function(nm){

    item_type <- var_contents[var_contents$group == '/var' & var_contents$name == nm, 'otype']
    if (!(any(item_type == "H5I_GROUP"))){
      return(as.vector(h5ad_file$var[[nm]]))
    }
    else{
      cats  <- as.vector(h5ad_file$var[[nm]][["categories"]])
      codes <- as.vector(h5ad_file$var[[nm]][["codes"]])

      # fix string encoding
      if (is.list(cats)) {cats <- unlist(cats)}
      if (is.character(cats)) {Encoding(cats) <- "UTF-8"}

      # convert codes → levels
      return(factor(cats[codes + 1], levels = cats))
    }
  })
  vardata <- as.data.frame(vardata)
  colnames(vardata) <- var_items
  rownames(vardata) <- vardata$'_index'
  vardata$'_index' <- NULL

  return(vardata)
}

#' Helper function to extract dimension reductions from H5AD
#'
#' @param h5ad_file Loaded H5File using rhdf5
#' @param dim_name String. Dimension Reduction to extract (PCA/UMAP)
#' @param cell_names Cells from SeuratObject::Cells()
#' @noRd
h5_createDimReduc <- function(h5ad_file, dim_name, cell_names){
  X_dim <- t(h5ad_file$'/obsm'[[dim_name]])
  colnames(X_dim) <- sapply(1:ncol(X_dim), function(x){paste(paste(unlist(strsplit(dim_name, ''))[3:nchar(dim_name)],collapse = ''), x, sep = "_")})
  rownames(X_dim) <- cell_names

  return(X_dim)
}

#' Helper function to extract centroids from H5AD for Seurat Object
#'
#' @param h5ad_file Loaded H5File using rhdf5
#' @param cell_names Cells from SeuratObject::Cells()
#' @noRd
h5_createCentroids_Seurat <- function(h5ad_file, cell_names){
  spatial <- t(h5ad_file$obsm[['spatial']])
  colnames(spatial) <- sapply(1:ncol(spatial), function(x){paste("spatial", x, sep = "_")})
  rownames(spatial) <- cell_names
  centroids <- SeuratObject::CreateCentroids(spatial[,c(1,2)])

  return(centroids)
}

#' Helper function to extract scale factors from H5AD for Seurat Object
#'
#' @param h5ad_file Loaded H5File using rhdf5
#' @param img_name String. Name of image to extract scale factors for.
#' @noRd
h5_createScalefactors_Seurat <- function(h5ad_file, img_name){
  if ("spot_diameter_fullres" %in% names(h5ad_file$'/uns/spatial'[[img_name]][["scalefactors"]])){
    spot_h5 <- as.numeric(h5ad_file$'/uns/spatial'[[img_name]]$scalefactors$spot_diameter_fullres)
  }
  else{
    spot_h5 <- 1
  }
  if ("fiducial_diameter_fullres" %in% names(h5ad_file$'/uns/spatial'[[img_name]][["scalefactors"]])){
    fiducial_h5 <- as.numeric(h5ad_file$"/uns/spatial"[[img_name]]$scalefactors$fiducial_diameter_fullres)
  }
  else{
    fiducial_h5 <- 1
  }
  if ("tissue_hires_scalef" %in% names(h5ad_file$'/uns/spatial'[[img_name]][["scalefactors"]])){
    hires_h5 <- as.numeric(h5ad_file$"/uns/spatial"[[img_name]]$scalefactors$tissue_hires_scalef)
  }
  else{
    hires_h5 <- 1
  }
  if ("tissue_lowres_scalef" %in% names(h5ad_file$'/uns/spatial'[[img_name]][["scalefactors"]])){
    lowres_h5 <- as.numeric(h5ad_file$"/uns/spatial"[[img_name]]$scalefactors$tissue_lowres_scalef)
  }
  else{
    lowres_h5 <- 1
  }
  sf <- Seurat::scalefactors(spot = spot_h5, fiducial = fiducial_h5, hires = hires_h5, lowres = lowres_h5)
  return(sf)
}
#' Helper function to extract centroids from H5AD for SpatialExperiment Object
#'
#' @param h5ad_file Loaded H5File using rhdf5
#' @param cell_names Cells from colnames(SpatialExperiment Obj)
#' @noRd
h5_createCentroids_SPE <- function(h5ad_file, cell_names){
  spatial <- t(h5ad_file$'/obsm/spatial')
  rownames(spatial) <- cell_names
  centroids <- spatial[,c(1,2)]
  colnames(centroids) <- c('x','y')
  return(centroids)
}

#' Helper function to extract scale factors from H5AD for SpatialExperiment Object
#'
#' @param h5ad_file Loaded H5File using rhdf5
#' @param img_name String. Name of image to extract scale factors for.
#' @param res Which resolution (lowres, hires) to extract scale factors for.
#' @noRd
h5_createScalefactors_SPE <- function(h5ad_file, img_name, res){
  if (res == 'hires'){
    if ("tissue_hires_scalef" %in% names(h5ad_file$'/uns/spatial'[[img_name]][["scalefactors"]])){
      sf <- as.vector(h5ad_file$"/uns/spatial"[[img_name]]$scalefactors$tissue_hires_scalef)
    }
    else{
      sf <- 1
    }
  }
  else if (res == 'lowres'){
    if ("tissue_lowres_scalef" %in% names(h5ad_file$'/uns/spatial'[[img_name]][["scalefactors"]])){
      sf <- as.vector(h5ad_file$"/uns/spatial"[[img_name]]$scalefactors$tissue_lowres_scalef)
    }
    else{
      sf <- 1
    }
  }
  return(sf)
}

#' Helper function for image for Seurat IHC images
#' @param img image as 3D array
#' @noRd
flipFunc <- function(img){
  return(img[,ncol(img):1,])
}

#' Helper function to extract arrays/layers from seurat object
#' @param seurat_obj Seurat Object
#' @noRd
seurat_assayMatrix <- function(seurat_obj){

  # Extract count matrices as list
  assay_names <- SeuratObject::Assays(seurat_obj)
  layer_names <- lapply(assay_names, function(x){SeuratObject::Layers(seurat_obj, assay = x)})
  names(layer_names) <- assay_names
  assay_list <- vector("list")
  for (i in assay_names){
    for (j in layer_names[[i]]){
      assay_list[[paste(i, j, sep = "_")]] <- SeuratObject::LayerData(seurat_obj, assay = i, layer = j)
    }
  }

  # Check if assays are the same size/have the same features and genes
  filter_req = FALSE
  if (!all(sapply(assay_list, function(x) identical(dim(x), dim(assay_list[[1]]))))){
    stop("Assays are not of the same size, as required by SingleCellExperiment")
  }

  all_rows_match <- all(sapply(assay_list, function(x) identical(rownames(x), rownames(assay_list[[1]]))))
  all_cols_match <- all(sapply(assay_list, function(x) identical(colnames(x), colnames(assay_list[[1]]))))

  if (!(all_rows_match & all_cols_match)){
    stop('Assays do not contain the same feature/cell names, as required by SingleCellExperiment.')
  }

  return(assay_list)
}

#' Helper function to extract spatial coords from Seurat as matrix
#' @param seurat_obj Seurat Object
#' @noRd
seurat_createSpatialCoords <- function(seurat_obj){
  spatial <- as.matrix(Seurat::GetTissueCoordinates(seurat_obj))
  spatial <- spatial[,-3]
  mode(spatial) <- "numeric"
  return(spatial)
}

#' Helper function to extract scale factors from Seurat Object for SpatialExperiment Object
#'
#' @param seurat_obj Seurat Object
#' @param img_name String. Name of image to extract scale factors for.
#' @noRd
seurat_createScalefactors_SPE <- function(seurat_obj, img_name){
    if ("lowres" %in% names(seurat_obj@images[[img_name]]@scale.factors)){
      sf <- seurat_obj@images[[img_name]]@scale.factors$lowres
    }
    else{
      sf <- 1
    }
  return(sf)
}

#' Helper function to extract scale factors from Seurat Object for AnnData Object
#'
#' @param seurat_obj Seurat Object
#' @param img_name String. Name of image to extract scale factors for.
#' @noRd
seurat_createScalefactors_ad <- function(seurat_obj, img_name){
  sf_orig <- seurat_obj@images[[img_name]]@scale.factors
  sf <- list()
  if ("spot" %in% names(sf_orig)){
    sf[["spot_diameter_fullres"]] <- sf_orig$spot
  }
  else{
    sf[["spot_diameter_fullres"]] <- 1
  }

  if ("fiducial" %in% names(sf_orig)){
    sf[["fiducial_diameter_fullres"]] <- sf_orig$fiducial
  }
  else{
    sf[["fiducial_diameter_fullres"]] <- 1
  }

  if ("hires" %in% names(sf_orig)){
    sf[["tissue_hires_scalef"]] <- sf_orig$hires
  }
  else{
    sf[["tissue_hires_scalef"]] <- 1
  }

  if("lowres" %in% names(sf_orig)){
    sf[["tissue_lowres_scalef"]] <- sf_orig$lowres
  }
  else{
    sf[["tissue_lowres_scalef"]] <- 1
  }

  return(sf)
}

