# this function creates all the files for data visualization tool

#' Prepares Data from Seurat object and launches SpatialView web application locally.
#'
#' @param seuratObj Seurat object.
#' @param dataPaths String array of paths for the raw samples.
#' Number of elements in the array should same as number of samples used in the Seurat object.
#' @param exportPath Directory path for exporting the files.
#' If downloadRepo = TRUE (default), then the SpatialView application will be created at this place.
#' If downloadRepo = FALSE, then download SpatialView from GitHub manually and set exportPath to SpatialView/data.
#' @param  projectName Name of the project.
#' If downloadRepo= TRUE (default), the projectName is used to create the peoject directory at exportPath, else ignores.
#' @param clusterColumn Cluster column name in the Seurat object metadata.
#' @param clusterNames Array of names for clusters matching to clusters in 'clusterColumn'.
#' @param clusterGenes List of arrays for genes specific to clusters, if interested to see in the visualization.
#' @param clusterColors Array of hex colors for clusters, else default colors will be used.
#' @param dimPlot Reduced dimension to be used for cluster visualization.
#' @param slot Slot name in seuratObject for extracting data. Default 'data'.
#' @param assayName Assay name in seuratObject. If not provided then default assay in seuratObject will be used.
#' @param multiSamplePattern Separator used for multiple samples in seuratObject. Default '_'.
#' @param exprRound The expressions are rounded for effectively saving space. Default rounding will be three decimal places.
#' @param spatialSubDir Sub-directory name where spatial files are located for 10x.
#' @param sampleInfo A dataframe with metadata information for the samples. Rows are for samples.
#' @param downloadRepo
#' If TRUE (default), downloads SpatialView from the GitHub repository,
#' and runs the SpatialView application on a local computer.
#'
#' When TRUE, the updated SpatialView files are downloaded from GitHub (spatialviewRepo URL).
#' To download the files, SpatialViewR uses the 'wget' utility. The 'wget' package is pre-installed on most Linux distributions.
#' For Mac, see \link{https://stackoverflow.com/questions/33886917/how-to-install-wget-in-macos}
#'
#' For Windows, see   \link{https://gnuwin32.sourceforge.net/packages/wget.htm}
#'
#' When downloadRepo set to FALSE, SpatialView files may be downloaded manually from the GitHub repository
#' \link{https://github.com/kendziorski-lab/spatialview/archive/refs/tags/spatialview-latest.zip}
#' and exportPath to be set to <PATH TO SPATIALVIEW DIR>/data/
#'
#'For details please refer to the \href{https://raw.githubusercontent.com/kendziorski-lab/kendziorski-lab.github.io/main/projects/spatialview/user-guide.pdf}{user-guide}
#'
#' @param spatialviewRepo Downloads SpatialView from GitHub repository (ignored if downloadRepo = FALSE).
#' @param spatialviewVersion Downloads the specified version of SpatialView (ignored if downloadRepo = FALSE).
#' @param port Port to be used for running SpatialView (ignored if downloadRepo = FALSE).
#' @param launchApp If TRUE (default), then launches the web application (ignored if downloadRepo = FALSE).
#' @param export_sparse If TRUE (default, recommended), sparse expression matrix will be exported, else dense csv format will be used.
#' @param data_file_name_expressions Name of the exported expression matrix in csv format (ignored if export_sparse = TRUE).
#' @param data_file_name_expressions_sparse Name of the exported expression matrix in sparse format (ignored if export_sparse = FALSE).
#' @param data_file_name_genes Name of the exported file containing names of the genes (ignored if export_sparse = FALSE).
#' @param data_file_name_barcodes Name of the exported file containing barcodes (ignored if export_sparse = FALSE).
#' @param configList A list of named strings for changing the default configuration.
#'
#' Names should be the names of keys and values should be the corresponding values.
#'
#' This works if 'downloadRepo' is set to TRUE, else ignored.
#' Configuration changes are having low precedence as compared to other parameters.
#'
#' @param verbose TRUE/FALSE whether to print comments.
#'
#' @return None
#' @export
#'
prepare10xVisium_from_seurat <- function(seuratObj,
                                   dataPaths,
                                   exportPath,
                                   projectName = "spatial",
                                   clusterColumn = "seurat_clusters",
                                   clusterNames = NULL,
                                   clusterGenes = NULL,
                                   clusterColors = NULL,
                                   dimPlot = "umap",
                                   slot = "data",
                                   assayName = NULL,
                                   multiSamplePattern = "_",
                                   exprRound = 3,
                                   spatialSubDir = "spatial",
                                   sampleInfo = NULL,
                                   downloadRepo = TRUE,
                                   spatialviewRepo = "https://github.com/kendziorski-lab/spatialview/archive/refs/tags/",
                                   spatialviewVersion = "spatialview-latest",
                                   port = 8878,
                                   launchApp = TRUE,
                                   export_sparse = TRUE,
                                   data_file_name_expressions = "expression_matrix.csv",
                                   data_file_name_expressions_sparse = "expression_matrix_sparse.txt",
                                   data_file_name_genes = "genes.csv",
                                   data_file_name_barcodes = "barcodes.csv",
                                   configList = NA,
                                   verbose = FALSE) {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package Seurat must be installed to use this function.",
         call. = FALSE)
  }

  orign.exportPath <- exportPath
  seurat.idents <- unique(seuratObj$orig.ident)
  sel.ident <- NA
  if (length(seurat.idents) > 1 &
      length(seurat.idents) != length(dataPaths)) {
    stop(
      "Multiple indents found in the Seurat object. Length of dataPaths is not same as number of indents.\n"
    )
  }

  # check if file path exists
  if (!dir.exists(exportPath))
    stop(paste(exportPath, "Does not exist!"))

  stopifnot(is.logical(downloadRepo))
  stopifnot(is.logical(launchApp))
  stopifnot(is.logical(export_sparse))
  stopifnot(is.logical(verbose))

  if (downloadRepo) {
    tryCatch(
      download.file2(
        url = paste0(spatialviewRepo, spatialviewVersion, '.zip'),
        destfile = file.path(exportPath, paste0(projectName, ".zip")),
        method = "wget",
        cacheOK = FALSE,
        quiet = verbose
      ),
      finally = function(e)
        (stop(e)
        ))


    # unzip the .zip file
    dir.create(file.path(exportPath, "spatialview_temp_"), recursive = TRUE)
    utils::unzip(
      zipfile = file.path(exportPath, paste0(projectName, ".zip")),
      exdir = file.path(exportPath, "spatialview_temp_"),
      overwrite = TRUE
    )
    #remove the zip file once it is un-zipped
    file.remove(file.path(exportPath, paste0(projectName, ".zip")))

    spatialview_temp_dir_files <-
       list.dirs(file.path(exportPath, "spatialview_temp_"), recursive = FALSE)
    unlink(file.path(exportPath, projectName), recursive = TRUE)
    file.rename(
      file.path(spatialview_temp_dir_files),
      file.path(exportPath, projectName)
    )

    file.remove(file.path(exportPath, "spatialview_temp_"), recursive = TRUE)
    config.path <-
      file.path(exportPath, projectName, "config", "app_config.json")
    dataHTML.path <-
      file.path(exportPath, projectName, "config", "data_location.html")
    exportPath <- file.path(exportPath, projectName, "data")

    #Updating the configuration file
    spatialview.config <-
      rjson::fromJSON(file = config.path)

    config.attrs <- names(spatialview.config)
    if (!is.na(configList)) {
      if (!is.list(configList) | is.null(names(configList))) {
        config.change.attrs <- names(configList)
        for (i in seq_along(config.change.attrs)) {
          conf.att <- config.change.attrs[i]
          if (conf.att %in% config.attrs) {
            spatialview.config[conf.att] <- configList[[i]]
          }
        }
      } else{
          stop("configList should be a named list")
        }
    }

    spatialview.config["data_file_name_expressions"] <-
      data_file_name_expressions
    spatialview.config["data_file_name_expressions_sparse"] <-
      data_file_name_expressions_sparse
    spatialview.config["data_file_name_genes"] <- data_file_name_genes
    spatialview.config["data_file_name_barcodes"] <-
      data_file_name_barcodes
    spatialview.config["data_cluster_column"] <- clusterColumn
    spatialview.config["dim_plot"] <- dimPlot

    write(rjson::toJSON(spatialview.config), file = config.path)

    data.info <- data.frame(Samples = seurat.idents)
    if (!is.null(sampleInfo)) {
      if (nrow(sampleInfo) == length(seurat.idents)) {
        data.info <- cbind(data.info, sampleInfo)
      }else{
        stop('number of rows in sampleInfo is not matching with number of samples')
      }
    }else{
      sampleInfo <- data.info
    }

    print(xtable::xtable(data.info),
          type = "html",
          file = dataHTML.path)
  }


  for (i in seq_along(seurat.idents)) {
    sel.ident <- seurat.idents[i]
    data.path <- dataPaths[i]
    seuratObj.temp <-
      subset(seuratObj, subset = orig.ident == sel.ident)

    ## step 0
    # create the export directory
    dir_name <- sel.ident
    # remove space in names
    dir_name <-
      stringr::str_replace_all(string = dir_name,
                               pattern = " ",
                               replacement = "")
    #first char should not not a number
    if (!is.na(as.numeric(substring(dir_name, 1, 1)))) {
      dir_name <- paste0('X', dir_name)
    }
    #A name should be at least two characters long.
    if (nchar(dir_name) == 1) {
      dir_name <- paste0('X', dir_name)
    }

    if (verbose &
        dir_name != sel.ident)
      message(paste(dir_name, 'is created for', sel.ident))

    exportPath_sample <- file.path(exportPath, dir_name)
    dir.create(exportPath_sample, showWarnings = FALSE, recursive = TRUE)

    spatialDir.files <- list.dirs(file.path(data.path, spatialSubDir))
    # Step 1
    # checking scalefactors_json.json, this file is a must have one else error
    if (!file.exists(file.path(data.path, spatialSubDir, "scalefactors_json.json"))) {
      if (file.exists(file.path(data.path, spatialSubDir, "scalefactors_json.json.gz"))) {
        require(R.utils, quietly = TRUE)
        R.utils::gunzip(file.path(
          data.path,
          spatialSubDir,
          "scalefactors_json.json.gz"
        ),
        remove = FALSE)
      } else {
        stop(paste0(
          "scalefactors_json.json file not found at ",
          file.path(data.path, spatialSubDir, "scalefactors_json.json")
        ))
      }
    }
    file.copy(
      from = file.path(data.path, spatialSubDir, "scalefactors_json.json"),
      to = exportPath_sample,
      overwrite = TRUE,
      recursive = FALSE,
      copy.mode = TRUE
    )

    # checking tissue_positions.csv or tissue_positions_list.csv
    # this file is a must have one else error
    #ref: https://www.10xgenomics.com/support/software/space-ranger/analysis/outputs/spatial-outputs


    if (file.exists(file.path(data.path, spatialSubDir, "tissue_positions.csv"))) {
      file.copy(
        from = file.path(data.path, spatialSubDir, "tissue_positions.csv"),
        to = data.path,
        overwrite = TRUE,
        recursive = FALSE,
        copy.mode = TRUE
      )
    }else if (file.exists(file.path(data.path, spatialSubDir, "tissue_positions_list.csv"))) {
      file.copy(
        from = file.path(data.path, spatialSubDir, "tissue_positions_list.csv"),
        to = data.path,
        overwrite = TRUE,
        recursive = FALSE,
        copy.mode = TRUE
      )
    } else if (file.exists(file.path(data.path, spatialSubDir, "tissue_positions.csv.gz"))) {
      if (!requireNamespace("R.utils", quietly = TRUE)) {
        stop("Package R.utils must be installed to process your data.",
             call. = FALSE)
      }

      R.utils::gunzip(
        file.path(
          data.path,
          spatialSubDir,
          "tissue_positions.csv.gz"
        ),
        remove = FALSE
      )

      file.copy(
        from = file.path(data.path, spatialSubDir, "tissue_positions.csv"),
        to = data.path,
        overwrite = TRUE,
        recursive = FALSE,
        copy.mode = TRUE
      )

    } else if (file.exists(file.path(data.path, spatialSubDir, "tissue_positions_list.csv.gz"))) {
      if (!requireNamespace("R.utils", quietly = TRUE)) {
        stop("Package R.utils must be installed to process your data.",
             call. = FALSE)
      }
      R.utils::gunzip(
        file.path(
          data.path,
          spatialSubDir,
          "tissue_positions_list.csv.gz"
        ),
        remove = FALSE
      )
      file.copy(
        from = file.path(data.path, spatialSubDir, "tissue_positions_list.csv"),
        to = data.path,
        overwrite = TRUE,
        recursive = FALSE,
        copy.mode = TRUE
      )
    }else {
      stop(paste0(
        "tissue_positions.csv or tissue_positions_list.csv file not found in ",
        file.path(data.path, spatialSubDir)
      ))
    }



    # Step 2
    # Check either tissue_lowres_image.png or tissue_hires_image.png image present
    if (!(file.exists(
      file.path(data.path, spatialSubDir, "tissue_lowres_image.png")
    ) |
    file.exists(
      file.path(data.path, spatialSubDir, "tissue_hires_image.png")
    ))) {
      stop(paste0(
        "Expected either ",
        file.path(data.path, spatialSubDir, "tissue_lowres_image.png"),
        " or ",
        file.path(data.path, spatialSubDir, "tissue_hires_image.png"),
        ", but not found"
      ))
    }

    # Step 2.a
    # Checking if tissue_lowres_image.png present or not, if not then trying to use tissue_hires_image.png
    if ((!file.exists(
      file.path(data.path, spatialSubDir, "tissue_lowres_image.png")
    )) &
    file.exists(file.path(data.path, "spatial", "tissue_hires_image.png"))) {
      # copy tissue_hires_image.png as tissue_lowres_image.png
      file.copy(
        file.path(data.path, spatialSubDir, "tissue_hires_image.png"),
        file.path(data.path, spatialSubDir, "tissue_lowres_image.png"),
        overwrite = FALSE
      )


      # updating the scalefactors_json.json file
      scalefactors <-
        rjson::fromJSON(file = file.path(data.path, spatialSubDir, "scalefactors_json.json"))

      scalefactors$tissue_lowres_scalef <-
        scalefactors$tissue_hires_scalef
      write(
        rjson::toJSON(scalefactors),
        file.path(
          exportPath_sample,
          spatialSubDir,
          "scalefactors_json.json"
        )
      )
    }

    # Step 2.b
    # If tissue_hires_image.png is not present, then trying to use tissue_lowres_image.png
    if ((!file.exists(
      file.path(data.path, spatialSubDir, "tissue_hires_image.png")
    )) &
    file.exists(file.path(data.path, spatialSubDir, "tissue_lowres_image.png"))) {
      # copy tissue_lowres_image.png as tissue_hires_image.png
      file.copy(
        file.path(data.path, spatialSubDir, "tissue_lowres_image.png"),
        file.path(data.path, spatialSubDir, "tissue_hires_image.png"),
        overwrite = FALSE
      )


      # updating the scalefactors_json.json file
      scalefactors$tissue_hires_scalef <-
        scalefactors$tissue_lowres_scalef
      write(
        toJSON(scalefactors),
        file.path(
          exportPath_sample,
          spatialSubDir,
          "scalefactors_json.json"
        )
      )
    }

    file.copy(
      from = file.path(data.path, spatialSubDir, "tissue_hires_image.png"),
      to = file.path(exportPath_sample, "tissue_hires_image.png"),
      overwrite = TRUE,
      copy.mode = TRUE
    )

    # Step 3
    # creating sampleInfo
    if (!is.null(sampleInfo) & is.data.frame(sampleInfo)) {
      if (nrow(sampleInfo) >= i) {
        write.csv(
          sampleInfo[i,],
          file.path(exportPath_sample, "sample_info.csv"),
          row.names = FALSE,
          quote = TRUE
        )
      }
    }

    # Step 4
    # creating metadata
    metadata <- seuratObj.temp@meta.data
    if (clusterColumn %in% names(metadata)) {
      metadata$cluster <- metadata[, clusterColumn]
    } else {
      metadata$cluster <- 1
    }

    metadata$barcode <-
      stringr::str_split(rownames(metadata),
                         pattern = multiSamplePattern,
                         simplify = TRUE)[, 1]

    if (!is.null(dimPlot) & stringr::str_to_upper(dimPlot) %in% stringr::str_to_upper(Reductions(seuratObj.temp))) {
      existing.reductions <- Reductions(seuratObj.temp)

      reductions.df <- NULL
      for (r in existing.reductions) {
        if (stringr::str_to_upper(dimPlot) == stringr::str_to_upper(r)) {
          reductions.df <- Embeddings(seuratObj.temp, reduction = r)
          break
        }
      }

      if (dim(reductions.df)[2] > 2) warning("Only first two reduced dimentions used")
      if (dim(reductions.df)[2] < 2) stop("At least 2 reduced dimentions required.")

      metadata[paste0(dimPlot, "_1")] <- reductions.df[,1]
      metadata[paste0(dimPlot, "_2")] <- reductions.df[,2]

    }else{
      warning(paste(dimPlot, "is not present in the Seurat object. You may choose from available reductions", Reductions(seuratObj.temp)))
      metadata[,paste0(dimPlot, "_1")] <- 0
      metadata[,paste0(dimPlot, "_2")] <- 0
    }

    if (file.exists(file.path(data.path, spatialSubDir, "tissue_positions.csv"))) {
      spot_info <-
        read.csv(file.path(data.path, spatialSubDir, "tissue_positions.csv"),
                 header = TRUE)
    }else{
      spot_info <-
        read.csv(file.path(data.path, spatialSubDir, "tissue_positions_list.csv"),
                 header = FALSE)
      colnames(spot_info) <- c(
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres"
      )
    }

    metadata <- merge(x = metadata, y = spot_info)

    write.csv(
      metadata,
      file.path(exportPath_sample, "metadata.csv"),
      quote = TRUE,
      row.names = FALSE
    )


    # Step 5
    # normalized_counts <- seuratObj.temp@assays$Spatial@data
    normalized_counts <-
      Seurat::GetAssayData(seuratObj.temp, slot = slot, assay = assayName)
    colnames(normalized_counts) <-
      stringr::str_split(colnames(normalized_counts),
                         pattern = multiSamplePattern,
                         simplify = TRUE)[, 1]


    # creating expression_matrix data
    if (export_sparse) {
      normalized_counts <- round(normalized_counts, exprRound)
      Matrix::writeMM(
        normalized_counts,
        file = file.path(exportPath_sample, data_file_name_expressions_sparse)
      )

      # genes
      write.csv(
        rownames(normalized_counts),
        file.path(exportPath_sample, data_file_name_genes),
        row.names = FALSE,
        quote = FALSE
      )
      # barcodes
      write.csv(
        colnames(normalized_counts),
        file.path(exportPath_sample, data_file_name_barcodes),
        row.names = FALSE,
        quote = FALSE
      )

    } else{
      normalized_counts <-
        as.data.frame(t(as.matrix(normalized_counts)))
      normalized_counts <- round(normalized_counts, exprRound)

      normalized_counts$barcode <- rownames(normalized_counts)
      write.csv(
        normalized_counts,
        file.path(exportPath_sample, data_file_name_expressions),
        row.names = FALSE,
        quote = TRUE
      )
    }



    # Step 6
    # creating cluster info
    unique_clusters <- 1
    if (clusterColumn %in% names(seuratObj[[]])) {
      unique_clusters <- unique(seuratObj[[]][, clusterColumn])
      if (is.factor(unique_clusters)) {
        unique_clusters <- levels(unique_clusters)
      } else{
        unique_clusters <- sort(unique_clusters)
      }
    }
    if (is.null(clusterColors)) {
      a <- ggplot2::scale_color_hue(h.start = 0)
      clusterColors <-
        a$palette(length(unique_clusters)) # number of colors
    }


    if (is.null(clusterNames)) {
      # clusterNames <- rep("undefined", length(unique_clusters))
      clusterNames <- unique_clusters
    }
    if (is.null(clusterGenes)) {
      clusterGenes <- rep("undefined", length(unique_clusters))
    }

    clusterInfo.df <- data.frame(
      cluster = unique_clusters,
      color = clusterColors,
      name = clusterNames,
      genes = unlist(lapply(clusterGenes, toString))
    )
    write.csv(
      clusterInfo.df,
      file.path(exportPath_sample, "cluster_info.csv"),
      row.names = FALSE,
      quote = TRUE
    )
    remove(seuratObj.temp)
    if (verbose) {
      message(paste("All required files for", sel.ident , "are found/prepared."))
    }
  }
  if (downloadRepo & launchApp) {
    start.httpserver(file.path(orign.exportPath, projectName),
                     port = port,
                     verbose = verbose)
  }

}


# this function creates all the files for data visualization tool from SpatialExperiment

#' Prepares Data from SpatialExperiment object and launches SpatialView web application locally.
#'
#' @param speObj SpatialExperiment object.
#' @param dataPaths String array of paths for the raw samples.
#' Number of elements in the array should same as number of samples used in the SpatialExperiment object.
#' @param exportPath Directory path for exporting the files.
#' If downloadRepo = TRUE, then the SpatialView application will be created at this place.
#' If downloadRepo = FALSE, then download SpatialView from GitHub manually and set exportPath to SpatialView/data.
#' @param  projectName Name of the project.
#' If downloadRepo= TRUE, the projectName is used to create the peoject directory at exportPath, else ignores.
#' @param clusterColumn Cluster column name in the SpatialExperiment object metadata.
#' @param clusterNames Array of names for clusters matching to clusters in 'clusterColumn'.
#' @param clusterGenes List of arrays for genes specific to clusters, if interested to see in the visualization.
#' @param clusterColors Array of hex colors for clusters, else default colors will be used.
#' @param dimPlot Reduced dimension to be used for cluster visualization.
#' @param slot Slot name in speObject for extracting data. Default 'data'.
#' @param assayName Assay name in speObject. If not provided then default Assay in speObject will be used.
#' @param multiSamplePattern Separator used for multiple samples in speObject. Default '_'.
#' @param exprRound The expressions are rounded for effectively saving space. Default rounding will be three decimal places.
#' @param spatialSubDir Sub-directory name where spatial files are located for 10x.
#' @param sampleInfo A dataframe with metadata information for the samples. Rows are for samples.
#' @param downloadRepo
#' If TRUE (default), downloads SpatialView from the GitHub repository,
#' and runs the SpatialView application on a local computer.
#'
#' When TRUE, the updated SpatialView files are downloaded from GitHub (spatialviewRepo URL).
#' To download the files, SpatialViewR uses the 'wget' utility. The 'wget' package is pre-installed on most Linux distributions.
#' For Mac, see \link{https://stackoverflow.com/questions/33886917/how-to-install-wget-in-macos}
#'
#' For Windows, see   \link{https://gnuwin32.sourceforge.net/packages/wget.htm}
#'
#' When downloadRepo set to FALSE, SpatialView files may be downloaded manually from the GitHub repository
#' \link{https://github.com/kendziorski-lab/spatialview/archive/refs/tags/spatialview-latest.zip}
#' and exportPath to be set to <PATH TO SPATIALVIEW DIR>/data/
#'
#'For details please refer to the \href{https://raw.githubusercontent.com/kendziorski-lab/kendziorski-lab.github.io/main/projects/spatialview/user-guide.pdf}{user-guide}
#'
#' @param spatialviewRepo  Downloads the SpatialView from GitHub repository (ignored if downloadRepo = FALSE).
#' @param spatialviewVersion Downloads the specified version of SpatialView (ignored if downloadRepo = FALSE).
#'                  The latest  version is maintained in the 'spatialview-latest'.
#'
#' @param port Port to be used for running SpatialView (ignored if downloadRepo = FALSE).
#' @param launchApp If TRUE (default), then launches the web application (ignored if downloadRepo = FALSE).
#' @param export_sparse If TRUE (default, recommended), sparse expression matrix will be exported, else dense csv format will be used.
#' @param data_file_name_expressions Name of the exported expression matrix in csv format (ignored if export_sparse = TRUE).
#' @param data_file_name_expressions_sparse Name of the exported expression matrix in sparse format (ignored if export_sparse = FALSE).
#' @param data_file_name_genes Name of the exported file containing names of the genes (ignored if export_sparse = FALSE).
#' @param data_file_name_barcodes Name of the exported file containing barcodes (ignored if export_sparse = FALSE).
#' @param configList A list of named strings for changing the configuration.
#'
#' Names should be the names of keys and values should be the corresponding values.
#' This works if 'downloadRepo' is set to TRUE, else ignored.
#' Configuration changes are having low precedence as compared to other parameters.
#'
#' @param verbose TRUE/FALSE whether to print comments.
#'
#' @return None
#' @export
#'
#' @examples
prepare10xVisium_from_SpatialExperiment <- function(speObj,
                                              dataPaths,
                                              exportPath,
                                              projectName = "spatial",
                                              clusterColumn = "label",
                                              clusterNames = NULL,
                                              clusterGenes = NULL,
                                              clusterColors = NULL,
                                              dimPlot = "umap",
                                              slot = "data",
                                              assayName = "logcounts",
                                              multiSamplePattern = "_",
                                              exprRound = 3,
                                              spatialSubDir = "spatial",
                                              sampleInfo = NULL,
                                              downloadRepo = TRUE,
                                              spatialviewRepo = "https://github.com/kendziorski-lab/spatialview/archive/refs/tags/",
                                              spatialviewVersion = "spatialview-latest",
                                              port = 8878,
                                              launchApp = TRUE,
                                              export_sparse = TRUE,
                                              data_file_name_expressions = "expression_matrix.csv",
                                              data_file_name_expressions_sparse = "expression_matrix_sparse.txt",
                                              data_file_name_genes = "genes.csv",
                                              data_file_name_barcodes = "barcodes.csv",
                                              configList = NA,
                                              verbose = FALSE) {

  if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
    stop("Package SpatialExperiment must be installed to use this function.",
         call. = FALSE)
  }

  orign.exportPath <- exportPath
  spe.idents <- unique(speObj$sample_id)
  sel.ident <- NA
  if (length(spe.idents) > 1 &
      length(spe.idents) != length(dataPaths)) {
    stop(
      "Multiple sample_id found in the SpatialExperiment object. Length of dataPaths is not same as number of sample_id.\n"
    )
  }

  # check if file path exists
  if (!dir.exists(exportPath))
    stop(paste(exportPath, "Does not exist!"))

  if (downloadRepo) {
    # if(file.exists(file.path(exportPath, paste0( ".zip")))){
    #   file.remove(file.path(exportPath, paste0(projectName, ".zip")))
    # }
    # TODO:: remove the url and provide the public url
    tryCatch(
      download.file2(
        url = paste0(spatialviewRepo, spatialviewVersion, '.zip'),
        destfile = file.path(exportPath, paste0(projectName, ".zip")),
        method = "wget",
        cacheOK = FALSE,
        quiet = verbose
      ),
      finally = function(e)
        (stop(e)
        ))


    # unzip the .zip file
    dir.create(file.path(exportPath, "spatialview_temp_"), recursive = TRUE)
    utils::unzip(
      zipfile = file.path(exportPath, paste0(projectName, ".zip")),
      exdir = file.path(exportPath, "spatialview_temp_"),
      overwrite = TRUE
    )
    file.remove(file.path(exportPath, paste0(projectName, ".zip")))

    spatialview_temp_dir <-
      list.files(file.path(exportPath, "spatialview_temp_"))
    unlink(file.path(exportPath, projectName), recursive = TRUE)
    file.rename(
      file.path(exportPath, "spatialview_temp_", spatialview_temp_dir),
      file.path(exportPath, projectName)
    )

    file.remove(file.path(exportPath, "spatialview_temp_"))
    config.path <-
      file.path(exportPath, projectName, "config", "app_config.json")
    dataHTML.path <-
      file.path(exportPath, projectName, "config", "data_location.html")
    exportPath <- file.path(exportPath, projectName, "data")

    #Updating the configuration file
    spatialview.config <-
      rjson::fromJSON(file = config.path)

    config.attrs <- names(spatialview.config)
    if (!is.na(configList)) {
      if (!is.list(configList) | is.null(names(configList))) {
        config.change.attrs <- names(configList)
        for (i in seq_along(config.change.attrs)) {
          conf.att <- config.change.attrs[i]
          if (conf.att %in% config.attrs) {
            spatialview.config[conf.att] <- configList[[i]]
          }
        }
      } else{
        stop("configList should be a named list")
      }
    }

    spatialview.config["data_file_name_expressions"] <-
      data_file_name_expressions
    spatialview.config["data_file_name_expressions_sparse"] <-
      data_file_name_expressions_sparse
    spatialview.config["data_file_name_genes"] <- data_file_name_genes
    spatialview.config["data_file_name_barcodes"] <-
      data_file_name_barcodes
    spatialview.config["data_cluster_column"] <- clusterColumn
    spatialview.config["dim_plot"] <- dimPlot


    write(rjson::toJSON(spatialview.config), file = config.path)

    data.info <- data.frame(Samples = spe.idents)
    if (!is.null(sampleInfo)) {
      if (nrow(sampleInfo) == length(spe.idents)) {
        data.info <- cbind(data.info, sampleInfo)
      }else{
        stop('number of rows in sampleInfo is not matching with number of samples')
      }
    }else{
      sampleInfo <- data.info
    }

    print(xtable::xtable(data.info),
          type = "html",
          file = dataHTML.path)
  }


  for (i in seq_along(spe.idents)) {
    sel.ident <- spe.idents[i]
    data.path <- dataPaths[i]
    speObj.temp <-
      speObj[, speObj$sample_id == sel.ident]

    ## step 0
    # create the export directory
    dir_name <- sel.ident
    # remove space in names
    dir_name <-
      stringr::str_replace_all(string = dir_name,
                               pattern = " ",
                               replacement = "")
    #first char should not not a number
    if (!is.na(as.numeric(substring(dir_name, 1, 1)))) {
      dir_name <- paste0('X', dir_name)
    }
    #A name should be at least two characters long.
    if (nchar(dir_name) == 1) {
      dir_name <- paste0('X', dir_name)
    }

    if (verbose &
        dir_name != sel.ident)
      message(paste(dir_name, 'is created for', sel.ident))

    exportPath_sample <- file.path(exportPath, dir_name)
    dir.create(exportPath_sample, showWarnings = FALSE, recursive = TRUE)

    spatialDir.files <- list.dirs(file.path(data.path, spatialSubDir))
    # Step 1
    # checking scalefactors_json.json, this file is a must have one else error
    if (!file.exists(file.path(data.path, spatialSubDir, "scalefactors_json.json"))) {
      if (file.exists(file.path(data.path, spatialSubDir, "scalefactors_json.json.gz"))) {
        require(R.utils, quietly = TRUE)
        R.utils::gunzip(file.path(
          data.path,
          spatialSubDir,
          "scalefactors_json.json.gz"
        ),
        remove = FALSE)
      } else {
        stop(paste0(
          "scalefactors_json.json file not found at ",
          file.path(data.path, spatialSubDir, "scalefactors_json.json")
        ))
      }
    }
    file.copy(
      from = file.path(data.path, spatialSubDir, "scalefactors_json.json"),
      to = exportPath_sample,
      overwrite = TRUE,
      recursive = FALSE,
      copy.mode = TRUE
    )

    # checking tissue_positions.csv or tissue_positions_list.csv
    # this file is a must have one else error
    #ref: https://www.10xgenomics.com/support/software/space-ranger/analysis/outputs/spatial-outputs


    if (file.exists(file.path(data.path, spatialSubDir, "tissue_positions.csv"))) {
      file.copy(
        from = file.path(data.path, spatialSubDir, "tissue_positions.csv"),
        to = data.path,
        overwrite = TRUE,
        recursive = FALSE,
        copy.mode = TRUE
      )
    }else if (file.exists(file.path(data.path, spatialSubDir, "tissue_positions_list.csv"))) {
      file.copy(
        from = file.path(data.path, spatialSubDir, "tissue_positions_list.csv"),
        to = data.path,
        overwrite = TRUE,
        recursive = FALSE,
        copy.mode = TRUE
      )
    } else if (file.exists(file.path(data.path, spatialSubDir, "tissue_positions.csv.gz"))) {
      if (!requireNamespace("R.utils", quietly = TRUE)) {
        stop("Package R.utils must be installed to process your data.",
             call. = FALSE)
      }

      R.utils::gunzip(
        file.path(
          data.path,
          spatialSubDir,
          "tissue_positions.csv.gz"
        ),
        remove = FALSE
      )

      file.copy(
        from = file.path(data.path, spatialSubDir, "tissue_positions.csv"),
        to = data.path,
        overwrite = TRUE,
        recursive = FALSE,
        copy.mode = TRUE
      )

    } else if (file.exists(file.path(data.path, spatialSubDir, "tissue_positions_list.csv.gz"))) {
      if (!requireNamespace("R.utils", quietly = TRUE)) {
        stop("Package R.utils must be installed to process your data.",
             call. = FALSE)
      }
      R.utils::gunzip(
        file.path(
          data.path,
          spatialSubDir,
          "tissue_positions_list.csv.gz"
        ),
        remove = FALSE
      )
      file.copy(
        from = file.path(data.path, spatialSubDir, "tissue_positions_list.csv"),
        to = data.path,
        overwrite = TRUE,
        recursive = FALSE,
        copy.mode = TRUE
      )
    }else {
      stop(paste0(
        "tissue_positions.csv or tissue_positions_list.csv file not found in ",
        file.path(data.path, spatialSubDir)
      ))
    }

    # Step 2
    # Check either tissue_lowres_image.png or tissue_hires_image.png image present
    if (!(file.exists(
      file.path(data.path, spatialSubDir, "tissue_lowres_image.png")
    ) |
    file.exists(
      file.path(data.path, spatialSubDir, "tissue_hires_image.png")
    ))) {
      stop(paste0(
        "Expected either ",
        file.path(data.path, spatialSubDir, "tissue_lowres_image.png"),
        " or ",
        file.path(data.path, spatialSubDir, "tissue_hires_image.png"),
        ", but not found"
      ))
    }

    # Step 2.a
    # Checking if tissue_lowres_image.png present or not, if not then trying to use tissue_hires_image.png
    if ((!file.exists(
      file.path(data.path, spatialSubDir, "tissue_lowres_image.png")
    )) &
    file.exists(file.path(data.path, "spatial", "tissue_hires_image.png"))) {
      # copy tissue_hires_image.png as tissue_lowres_image.png
      file.copy(
        file.path(data.path, spatialSubDir, "tissue_hires_image.png"),
        file.path(data.path, spatialSubDir, "tissue_lowres_image.png"),
        overwrite = FALSE
      )


      # updating the scalefactors_json.json file
      scalefactors <-
        rjson::fromJSON(file = file.path(data.path, spatialSubDir, "scalefactors_json.json"))

      scalefactors$tissue_lowres_scalef <-
        scalefactors$tissue_hires_scalef
      write(
        rjson::toJSON(scalefactors),
        file.path(
          exportPath_sample,
          spatialSubDir,
          "scalefactors_json.json"
        )
      )
    }

    # Step 2.b
    # If tissue_hires_image.png is not present, then trying to use tissue_lowres_image.png
    if ((!file.exists(
      file.path(data.path, spatialSubDir, "tissue_hires_image.png")
    )) &
    file.exists(file.path(data.path, spatialSubDir, "tissue_lowres_image.png"))) {
      # copy tissue_lowres_image.png as tissue_hires_image.png
      file.copy(
        file.path(data.path, spatialSubDir, "tissue_lowres_image.png"),
        file.path(data.path, spatialSubDir, "tissue_hires_image.png"),
        overwrite = FALSE
      )


      # updating the scalefactors_json.json file
      scalefactors$tissue_hires_scalef <-
        scalefactors$tissue_lowres_scalef
      write(
        toJSON(scalefactors),
        file.path(
          exportPath_sample,
          spatialSubDir,
          "scalefactors_json.json"
        )
      )
    }

    file.copy(
      from = file.path(data.path, spatialSubDir, "tissue_hires_image.png"),
      to = file.path(exportPath_sample, "tissue_hires_image.png"),
      overwrite = TRUE,
      copy.mode = TRUE
    )

    # Step 3
    # creating sampleInfo
    if (!is.null(sampleInfo) & is.data.frame(sampleInfo)) {
      if (nrow(sampleInfo) >= i) {
        write.csv(
          sampleInfo[i,],
          file.path(exportPath_sample, "sample_info.csv"),
          row.names = FALSE,
          quote = TRUE
        )
      }
    }

    # Step 4
    # creating metadata
    metadata <- speObj.temp@colData
    if (clusterColumn %in% names(metadata)) {
      metadata$cluster <- metadata[, clusterColumn]
    } else {
      metadata$cluster <- 1
    }

    metadata$barcode <-
      stringr::str_split(rownames(metadata),
                         pattern = multiSamplePattern,
                         simplify = TRUE)[, 1]

    if (!is.null(dimPlot) & stringr::str_to_upper(dimPlot) %in% stringr::str_to_upper(reducedDimNames(seuratObj.temp))) {
      existing.reductions <- reducedDimNames(seuratObj.temp)

      reductions.df <- NULL
      for (r in existing.reductions) {
        if (stringr::str_to_upper(dimPlot) == stringr::str_to_upper(r)) {
          reductions.df <- reducedDim(seuratObj.temp, type = r)
          break
        }
      }

      if (dim(reductions.df)[2] > 2) warning("Only first two reduced dimentions used")
      if (dim(reductions.df)[2] < 2) stop("At least 2 reduced dimentions required.")

      metadata[paste0(dimPlot, "_1")] <- reductions.df[,1]
      metadata[paste0(dimPlot, "_2")] <- reductions.df[,2]

    }else{
      warning(paste(dimPlot, "is not present in the Seurat object. You may choose from available reductions", reducedDimNames(seuratObj.temp)))
      metadata[,paste0(dimPlot, "_1")] <- 0
      metadata[,paste0(dimPlot, "_2")] <- 0
    }

    if (file.exists(file.path(data.path, spatialSubDir, "tissue_positions.csv"))) {
      spot_info <-
        read.csv(file.path(data.path, spatialSubDir, "tissue_positions.csv"),
                 header = TRUE)
    }else{
      spot_info <-
        read.csv(file.path(data.path, spatialSubDir, "tissue_positions_list.csv"),
                 header = FALSE)
      colnames(spot_info) <- c(
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres"
      )
    }

    metadata <- base::merge(x = metadata, y = spot_info)

    write.csv(
      metadata,
      file.path(exportPath_sample, "metadata.csv"),
      quote = TRUE,
      row.names = FALSE
    )


    # Step 5
    if (is.null(assayName) | assayName == '') {
      assayName <- SummarizedExperiment::assayNames(speObj.temp)[1]
    }
    normalized_counts <-
      as(SummarizedExperiment::assay(speObj.temp, i = assayName), "dgCMatrix")
    colnames(normalized_counts) <-
      stringr::str_split(colnames(normalized_counts),
                         pattern = multiSamplePattern,
                         simplify = TRUE)[, 1]


    # creating expression_matrix data
    if (export_sparse) {
      normalized_counts <- round(normalized_counts, exprRound)
      Matrix::writeMM(
        normalized_counts,
        file = file.path(exportPath_sample, data_file_name_expressions_sparse)
      )

      # genes
      write.csv(
        rownames(normalized_counts),
        file.path(exportPath_sample, data_file_name_genes),
        row.names = FALSE,
        quote = FALSE
      )
      # barcodes
      write.csv(
        colnames(normalized_counts),
        file.path(exportPath_sample, data_file_name_barcodes),
        row.names = FALSE,
        quote = FALSE
      )

    } else{
      normalized_counts <-
        as.data.frame(t(as.matrix(normalized_counts)))
      normalized_counts <- round(normalized_counts, exprRound)

      normalized_counts$barcode <- rownames(normalized_counts)
      write.csv(
        normalized_counts,
        file.path(exportPath_sample, data_file_name_expressions),
        row.names = FALSE,
        quote = TRUE
      )
    }



    # Step 6
    # creating cluster info
    unique_clusters <- 1
    if (clusterColumn %in% names(colData(speObj))) {
      unique_clusters <- unique(colData(speObj)[, clusterColumn])
      if (is.factor(unique_clusters)) {
        unique_clusters <- levels(unique_clusters)
      } else{
        unique_clusters <- sort(unique_clusters)
      }
    }
    if (is.null(clusterColors)) {
      a <- ggplot2::scale_color_hue(h.start = 0)
      clusterColors <-
        a$palette(length(unique_clusters)) # number of colors
    }


    if (is.null(clusterNames)) {
      # clusterNames <- rep("undefined", length(unique_clusters))
      clusterNames <- unique_clusters
    }
    if (is.null(clusterGenes)) {
      clusterGenes <- rep("undefined", length(unique_clusters))
    }

    clusterInfo.df <- data.frame(
      cluster = unique_clusters,
      color = clusterColors,
      name = clusterNames,
      genes = unlist(lapply(clusterGenes, toString))
    )
    write.csv(
      clusterInfo.df,
      file.path(exportPath_sample, "cluster_info.csv"),
      row.names = FALSE,
      quote = TRUE
    )
    remove(speObj.temp)
    if (verbose) {
      message(paste("All required files for", sel.ident , "are found/prepared."))
    }
  }
  if (downloadRepo & launchApp) {
    start.httpserver(file.path(orign.exportPath, projectName),
                     port = port,
                     verbose = verbose)
  }

}




start.httpserver <- function(path, port = NULL, verbose = FALSE) {
  #if port is not provided, then 8000 port is used.
  http.port <- 8000
  if (!is.null(port))
    http.port <- port
  tryCatch({
    v = NULL
    v = system("python3 --version", intern = TRUE)
    if (is.null(v)) {
      stop("Python not found. Please install python 3.X.X https://www.python.org/")
    } else{
      http.command <- "python3 -m http.server"

      http.command <-
        paste(http.command, http.port, "--directory", path)
      a <- system(http.command , wait = FALSE)
      browseURL(paste0("http://localhost:", http.port))
      if (verbose) message(paste0("Application is served at http://localhost:",http.port))
    }
  }
  , finally = function(e) {
    warning(e)
    stop("Python not found. Please install python https://www.python.org/")
  })
}

download.file2 <-
  function(url,
           destfile,
           method,
           quiet = FALSE,
           mode = "w",
           cacheOK = TRUE,
           extra = getOption("download.file.extra"),
           headers = NULL,
           ...) {
    if (quiet)
      extra <- c(extra, "--quiet")
    system(paste(
      "wget",
      paste(extra, collapse = " "),
      shQuote(url),
      "-O",
      shQuote(path.expand(destfile))
    ))

  }
