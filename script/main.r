#
#' @title Load Reference Atlas
#' @description Load or download the reference map for dataset projection. By the default it downloads the reference
#' atlas of tumor-infiltrating lymphocytes (TILs)
#' @param ref Reference atlas as a Seurat object (by default downloads a mouse reference TIL atlas)
#' To use a custom reference atlas, provide a `.rds` or a URL to a `.rds` object, storing a Seurat object.
#' @examples 
#' load.reference.map()
# ---- main ----
load.reference.map <- function(ref='referenceTIL'){
  if(identical(ref,'referenceTIL')){
    cat('Loading Default Reference Atlas...\n')
    refFileName <- 'reference/ref_TILAtlas_mouse_v1.rds'
    refUrl <- 'https://ndownloader.figshare.com/files/23136746'
    if(file.exists(refFileName)){
      print('Loading dataset: ref_TILAtlas_mouse_v1.rds')
    }else{
      print(paste0(refFileName, ' not found, downloading reference TIL map from the sever...'))
      download.file(refUrl,refFileName)
    }
    ref <- readRDS(refFileName)
    print(paste0('Loaded Reference map: ',ref@misc$ProjecTME))
  }else if(identical(ref,'referenceLCMVatlas')){
    cat('Loading Reference LCMV Atlas...\n')
    refFileName <- 'reference/ref_LCMV_Atlas_mouse_v1.rds'
    refUrl <- 'https://ndownloader.figshare.com/files/23166794'
    if(file.exists(refFileName)){
      print('Loading dataset: ref_LCMV_Atlas_mouse_v1.rds')
    }else{
      print(paste0(refFileName, ' not found, downloading reference TIL map from the sever...'))
      download.file(refUrl,refFileName)
    }
    ref <- readRDS(refFileName)
    print(paste0('Loaded Reference map: ',ref@misc$ProjecTME))
  }else if(identical(ref,'referenceLCMVCD4')){
    cat('Loading Reference LCMV CD4...\n')
    refFileName <- 'reference/ref_LCMV_CD4_mouse_release_v1.rds'
    refUrl <- 'https://ndownloader.figshare.com/files/31057081'
    if(file.exists(refFileName)){
      print('Loading dataset: ref_LCMV_CD4_mouse_release_v1.rds')
    }else{
      print(paste0(refFileName, ' not found, downloading reference TIL map from the sever...'))
      download.file(refUrl,refFileName)
    }
    ref <- readRDS(refFileName)
    print(paste0('Loaded Reference map: ',ref@misc$ProjecTME))
  }else{
    if(identical(ref,'CustomReference')){
      print('Loading Custom Reference Atlas...')
      refFileName <- 'reference/custom_reference.rds'
      ref <- readRDS(refFileName)
      print(paste0('Loaded Reference map ',ref@misc$ProjecTME))
    }else if(grepl('^[ftp|http]', ref, perl = T)){
      refUrl <- ref
      refFileName <- paste0('reference/custom_reference.rds')
      cat(paste0('Trying to download custom reference from ',refUrl,'...'))
      download.file(refUrl, refFileName)
      ref <- readRDS(refFileName)
      print(paste0('Loaded Reference map: ',ref@misc$ProjecTME))
    }else{
      stop('Provide ref is not a valid reference or a valid URL.')
    }
  }
  return(ref)
}


#' @title Read to memory a query expression matrix
#' @description Load a query expression matrix to be projected onto the reference atlas.
#' Several formats (10x, hdf5, raw and log counts) are supported.
#' @param filename Path to expression matrix file or folder
#' @param type Expression matrix format (10x, hdf5, raw, raw.log2)
#' @param project.name Title for the project
#' @param min.cells Only keep genes represented in at least min.cells number of cells
#' @param min.features Only keep cells expressing at least min.features genes
#' @param gene.column.10x For 10x format - which column of genes.tsv or features.tsv to use for gene names
#' @param raw.rownames For raw matrix format - A vector of row names, or a single number giving the column
#' of the table which contains the row names.
#' @param raw.sep For raw matrix format - Separator for raw expression matrix
#' @param raw.header For raw matrix format - Use headers in expression matrix
#' @param use.readmtx Use ReadMtx function to read in 10x files with custom names
#' @return A Seurat object populated with raw counts and normalized counts for single-cell expression
#' fname <- 'data/Bortezomib_6hr_expt1'
#' querydata <- read.sc.query(fname, type='10x')
read.sc.query <- function(filename,type=c('10x','hdf5','raw','raw.log2'),project.name='Query',
                          min.cells = 3, min.features = 50, gene.column.10x=2,raw.rownames=1,
                          raw.sep=c('auto',' ','\t',','), raw.header = T, use.readmtx = F){
  if(is.null(filename)){
    stop('Please provide a query dataset in one of the supported formats: 10x, hdf5, raw, raw.log2')
  }
  if(tolower(type)=='10x'){
    fn <- list.files(filename)
    mtx.file <- grep('matrix.mtx', fn, value = T)[1]
    fea.file <- grep('features.tsv|genes.tsv',fn,value = T)[1]
    bar.file <- grep('barcodes.tsv', fn, value = T)[1]
    if(is.na(mtx.file)) stop('Cannot find matrix file')
    if(is.na(fea.file)) stop('Cannot find gens file')
    if(is.na(bar.file)) stop('Cannot find barcode file')
    mtx.file <- paste0(filename,'/',mtx.file)
    fea.file <- paste0(filename, '/', fea.file)
    bar.file <- paste0(filename, '/', bar.file)
    if(use.readmtx){
      query.exp <- readMTX(mtx = mtx.file,cells = bar.file,features = fea.file,feature.column = gene.column.10x)
    }else{
      query.exp <- Read10X(filename, gene.column = gene.column.10x)
    }
  }else if(identical(tolower(type),'hdf5')){
    query.exp <- Read10X_h5(filename)
  }else if(tolower(type) == 'raw' | tolower(type) == 'raw.log2'){
    raw.sep <- raw.sep[1]
    if(raw.sep == 'auto'){
      raw.sep <- guess_raw_separator(filename)
      if(is.null(raw.sep)){
        stop('Could not guess separator for raw matrix format. Try sepcifying manually with raw.sep parameter.\n')
      }
    }
    p <- regexpr('\\.([[:alnum:]]+)$', filename)
    extension <- ifelse(p > -1L, substring(filename, p + 1L), '')
    if(extension=='gz'){
      query.exp <- read.table(gzfile(filename),row.names = raw.rownames, sep = raw.sep, header = raw.header)
    }else{
      query.exp <- read.table(filename,row.names = raw.rownames,sep = raw.sep,header = raw.header)
    }
    query.exp[is.na(query.exp)] <- 0
    if(identical(tolower(type),'raw.log2')){ query.exp <- 2^(query.exp)-1 }
    # Also try to determine whether genes are on rows or columns
    load('data/Hs2Mm.convert.table.RData')
    gnames <- c(Hs2Mm.convert.table$Gene.MM, Hs2Mm.convert.table$Gene.stable.ID.HS,Hs2Mm.convert.table$Gene.HS)
    gr <- length(intersect(rownames(query.exp), gnames))
    gc <- length(intersect(colnames(query.exp), gnames))
    gmax <- max(gr, gc)
    if(gmax==0) stop('Could not find gene names in matrix. Check matrix format.\n')
    if(gc > gr) query.exp <- t(query.exp)
  }else{
    stop('Please provide a query dataset in one of the supported formats: 10x, hdf5, raw, raw.log2')
  }
  query.seurat <- CreateSeuratObject(counts=query.exp, project=project.name, min.cells=min.cells, min.features=min.features)
  query.seurat <- NormalizeData(query.seurat)
  return(query.seurat)
}


#' @title Project a query scRNA-seq dataset onto a reference atlas
#' @description This function allows projecting ('query') single-cell RNA-seq datasets onto a 
#' reference map (i.e. a curated and annotated scRNA-seq dataset).
#' To project multiple datasets, submit a list of Seurat objects with the query parameter.
#' The projection consists of 3 steps:
#' \itemize{
#'   \item{pre-processing: optional steps which might include pre-filtering of cells by markers
#'   using 'scGate', data normalization, and ortholog conversion.}
#'   \item{batch-effect correction: uses built-in STACAS algorithm to detect and correct for batch
#'   effects (this step assumes that at least a fraction of the cells in the query are in the same
#'   state than cells in the reference)}
#'   \item{embedding of corrected query data in the reduced-dimensionality spaces (PCA and UMAP)
#'   of the reference map.}
#' }
#' @param query Query data, either as single Seurat object or as a list of Seurat object
#' @param ref Reference Atlas - if NULL, downloads default TIL reference atlas
#' @param query.assay Which assay slot to use for the query (default to DefaultAssay(query))
#' @param direct.projection If TRUE, apply PCA transformation directly without alignment
#' @param STACAS.anchor.coverage Focus on few robust anchors (low STACAS.anchor.coverage) or on
#' a large amount of anchors (high STACAS.anchor.coverage). Must be number between 0 and 1.
#' @param STACAS.correction.scale Slope of sigmoid function used to determine strength of batch effect correction.
#' @param STACAS.k.anchor Integer. For alignment, how many neighbors (k) to use when picking anchors.
#' @param STACAS.k.wight Number of neighbors to consider when weighting anchors. Default is 'max',
#' which disables local anchor weighting.
#' @param skip.normalize By default, log-normalize the count data. If you have already normalized
#' your data, you can skip normalization.
#' @param filter.cells Pre-filter cells using 'scGate'. Only set to FALSE if the dataset has bee
#' previously subset to cell types represented in the reference.
#' @param scGate.model scGate model used to filter target cell type from query data
#' (if NULL, use the model stored in the `ref@misc$scGate`)
#' @param ortholog_table Dataframe for conversion between ortholog genes. (Default `Hs2Mm.covert.table`)
#' @param fast.umap.predict Fast approximation for UMAP projection. Uses coordinates of nearest neighbors
#' in PCA space to assign UMAP coordinates (credits to Changsheng Li for the implementation).
#' @param ncores Number of cores for parallel execution (requires `BiocParallel`)
#' @return An augmented Seurat object with projected UMAP coordinates on the reference map
#' @examples 
#' load('data/query_example_seurat.RData')
#' ref <- load.reference.map()
#' make.projection(query_example_seurat, ref = ref)
make.projection <- function(query,ref=NULL,filter.cells=TRUE,query.assay=NULL,direct.projection=FALSE,
                            STACAS.anchor.coverage=0.7, STACAS.correction.scale=100,
                            STACAS.k.anchor=5, STACAS.k.weight='max', skip.normalize=FALSE,
                            fast.umap.predict=FALSE, ortholog_table=NULL, scGate_model=NULL,
                            ncores = 1){
  if(is.null(ref)) ref <- load.reference.map()
  projected.list <- list()
  if(is.null(ortholog_table)){
    load('data/Hs2Mm.convert.table.RData')
    ortholog_table <- Hs2Mm.convert.table
  }
  if(!is.list(query)){
    query.list <- list(query = query)
  }else{
    query.list <- query
    if(is.null(names(query.list))){
      names(query.list) <- paste0('query', c(1:length(query.list)))
    }
  }
  rm(query);gc()
  # Parallelize
  param <- BiocParallel::MulticoreParam(workers = ncores)
  projected.list <- BiocParallel::bplapply(
    X = 1:length(query.list), BPPARAM = param, FUN = function(i){
      projection.helper(
        query = query.list[[i]], ref = ref, filter.cells = filter.cells, query.assay = query.assay,
        direct.projection = direct.projection, fast.umap.predict = fast.umap.predict, k.anchor = STACAS.k.anchor,
        k.weight = STACAS.k.weight, anchor.coverage = STACAS.anchor.coverage, correction.scale = STACAS.correction.scale,
        remove.thr = 0, alpha = 0.5, ncores = ncores, ortholog_table = ortholog_table,
        skip.normalize = skip.normalize, id = names(query.list)[i], scGate_model = scGate_model
      )
    }
  )
  names(projected.list) <- names(query.list)
  # De-list of single object was submitted
  if(length(projected.list)==1) projected.list <- projected.list[[1]]
  return(projected.list)
}


#' @title Predict cell states of a projected dataset
#' @description This function uses a nearest-neighbor algorithm to predict a feature (e.g. the cell state)
#' of the query cells. Distances between cells in the reference map and cells in the query are calculated
#' in a reduced space (PCA or UMAP) and the feature is assigned to query cells based on a consensus of its
#' nearest neighbors in the reference object.
#' @param ref Reference Atlas
#' @param query Seurat object with query data
#' @param reduction The dimensionality reduction used to calculate pairwise distances. One of 'pca' or 'umap'
#' @param ndim How many dimensions in the reduced space to be used for distance calculations
#' @param k Number of neighbors to assign the cell type
#' @param labels.col The metadata field of the reference to annotate the clusters (default: functional.cluster)
#' @return The query object submitted as parameter, with two additional metadata slots for predicted state
#' ans its confidence score
#' @examples 
#' load('data/query_example_seurat.RData')
#' ref <- load.reference.map()
#' q <- make.projection(query_example_seurat, ref = ref)
#' q <- cellstate.predict(ref, query = q)
#' table(q$functional.cluster)
cellstate.predict <- function(ref,query,reduction='pca',ndim=NULL,k=20,labels.col='functional.cluster'){
  if(is.null(ndim)){
    if(!is.null(ref@misc$umap_object$data)){
      ndim <- ncol(ref@misc$umap_object$data)
    }else{ stop('Please specify ndim parameter.')}
  }
  pca.dim <- dim(ref@misc$umap_object$data)[2]
  tdim <- dim(ref@reductions[[reduction]]@cell.embeddings)[2]
  if(ndim > tdim){
    warning(sprintf('Number of dimensions ndim=%i is larger than the dimensions in reduction %s 
                    - Using only first %i dimensions', ndim, reduction, tdim))
    ndim <- tdim
  }
  labels <- ref[[labels.col]][,1]
  ref.space <- ref@reductions[[reduction]]@cell.embeddings[,1:ndim]
  query.space <- query@reductions[[reduction]]@cell.embeddings[,1:ndim]
  pred.type <- rep('Unknown', dim(query.space)[1])
  pred.conf <- numeric(dim(query.space)[1])
  # Use NN-search wrapper in Seurat
  nn.method <- 'rann'
  nn.ranked <- Seurat:::NNHelper(data = ref.space, query = query.space, k = k, method = nn.method)
  for(r in 1:dim(nn.ranked@nn.idx)[1]){
    top.k <- nn.ranked@nn.idx[r,]
    scores <- sort(table(labels[top.k]), decreasing = T)/k
    pred.type[r] <- names(scores)[1]
    pred.conf[r] <- scores[1]
  }
  pred <- data.frame(cbind(rownames(query.space),pred.type,pred.conf),stringsAsFactors = F)
  rownames(pred) <- rownames(query.space)
  colnames(pred) <- c('id', 'pred.state', 'confidence')
  pred <- transform(pred, confidence = as.numeric(confidence))
  query <- AddMetaData(query, metadata = pred$pred.state, col.name = labels.col)
  query <- AddMetaData(query, metadata = pred$confidence, col.name = paste0(labels.col,'.conf'))
  message(paste0('Creating slots ',labels.col, ' and ',labels.col, '.conf in query object'))
  return(query)
}


#' @title Merge Seurat objects, including reductions (e.g. PCA, UMAP, ICA)
#' @description Given two Seurat objects, merge counts and data as well as dim reduction (PCA, UMAP, ICA, etc)
#' @param x First object to merge
#' @param y Second object to merge
#' @param ... More parameters to \link{merge} function
#' @return A merged Seurat object
#' @examples 
#' seurat.merged <- merge.Seurat.embeddings(obj.1, obj.2)
#' # To merge multiple object stored in a list
#' seurat.merged <- Reduce(merge.Seurat.embeddings, obj.list)
merge.Seurat.embeddings <- function(x,y,...){
  # first regular Seurat merge, inheriting parameters
  m <- merge(x, y, ...)
  # preserve reductions (PCA, UMAP, ...)
  reds <- intersect(names(x@reductions), names(y@reductions))
  for(r in reds){
    message(sprintf('Merging %s embeddings...',r))
    m@reductions[[r]] <- x@reductions[[r]]
    if(dim(y@reductions[[r]]@cell.embeddings)[1] > 0){
      m@reductions[[r]]@cell.embeddings <- rbind(m@reductions[[r]]@cell.embeddings, y@reductions[[r]]@cell.embeddings)
    }
    if(dim(y@reductions[[r]]@feature.loadings)[1] > 0){
      m@reductions[[r]]@feature.loadings <- rbind(m@reductions[[r]]@feature.loadings,y@reductions[[r]]@feature.loadings)
    }
    if(dim(y@reductions[[r]]@feature.loadings.projected)[1] > 0){
      m@reductions[[r]]@feature.loadings.projected <- rbind(m@reductions[[r]]@feature.loadings.projected,y@reductions[[r]]@feature.loadings.projected)
    }
  }
  return(m)
}


#' @title Project a query scRNA-seq dataset onto a reference atlas
#' @description This function allows projecting ('query') single-cell RNA-seq datasets onto a reference
#' map (i.e., a curated and annotated scRNA-seq dataset). To project multiple datasets, submit a list
#' of Seurat objects with the query parameter.
#' The projection consists of 3 steps:
#' \itemize{
#'   \item{pre-processing: optional steps which might include pre-filtering of cells by markers,
#'   data normalization, and ortholog conversion.}
#'   \item{batch-effect correction: use built-in STACAS algorithm to detect and correct for batch effects
#'   (this step assumes that at least a fraction of the cells in the qeury are in the same state than
#'   cells in the reference)}
#'   \item{embedding of corrected query data in the reduced-dimensionality space (PCA and UMAP) of the reference map.}
#' }
#' This function acts as a wrapper for \link{make.projection} and \link{cellstate.predict}
#' @param query Query data, either as single Seurat object or as a list of Seurat object
#' @param ref Reference Atlas
#' @param filter.cells Pre-filter cells. Only set to FALSE if the dataset has been previously subset to
#' cell types in the reference.
#' @param split.by Grouping variable to split the query object (e.g., if the object contains multiple samples)
#' @param reduction The dimensionality reduction used to assign cell type labels.
#' @param ndim The number of dimensions used for cell type classification
#' @param k Number of neighbors for cell type classification
#' @param labels.col The metadata field of the reference to annotate the clusters
#' @param ... Additional parameters to \link{make.projection}
#' @return An augmented Seurat object with projected UMAP coordinates on the reference map and cell classificaitons
#' @examples 
#' load('data/query_example_seurat.RData')
#' ref <- load.reference.map()
#' projected <- Run.ProjecTME(query_example_seurat, ref = ref, fast.umap.predict = T)
#' plot.projection(ref = ref, query = projected)
#' circMarkCelltype(seu = projected,reduction = 'umap',annot = 'functional.cluster',sub = c('CD8_Tex'),palette = rainbow(6))
#' 
#' ## Visualization
#' if(!require(dittoSeq)) BiocManager::install('dittoSeq')
#' library(dittoSeq)
#' # violin plot
#' dittoPlot(projected, 'Cd14', group.by = 'functional.cluster')
#' 
#' # jittered boxplot
#' dittoPlot(projected, 'Cd14', group.by = 'functional.cluster', plots = c('boxplot','jitter'))
#' 
#' # jittered ridgeplot
#' dittoPlot(projected, 'Cd14', group.by = 'functional.cluster', plots = c('ridgeplot','jitter'))
#' 
#' # Hex Dimplot
#' dittoDimHex(projected, 'Cd14')
#' 
#' # DimPlot
#' dittoDimPlot(projected, 'functional.cluster', size = 3)
#' 
#' # DimPlot (genes)
#' dittoDimPlot(projected, 'Cd14', size = 3)
#' 
#' # Cell proportion barplot
#' dittoBarPlot(projected, 'functional.cluster', group.by = 'Sample')
#' dittoBarPlot(projected, 'functional.cluster', group.by = 'Sample', scale = 'count')
#' 
#' # DimHeatmap
#' htmap.colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdYlBu')))(50)
#' dittoHeatmap(projected, genes = dittoSeq::getGenes(projected)[1:20],heatmap.colors = htmap.colors)
#' # binary heatmap
#' dittoHeatmap(projected,genes = dittoSeq::getGenes(projected)[1:20],
#' annot.by = c('Sample','SampleLab','nFeature_RNA','functional.cluster'),
#' scaled.to.max = T, treeheight_row = 10)
#' 
Run.ProjecTME <- function(query,ref=NULL,filter.cells=T,split.by=NULL,reduction='pca',ndim=NULL,
                           k = 20, labels.col='functional.cluster', ...){
  if(!is.null(split.by)){
    if(!split.by %in% colnames(query[[]])){
      stop(sprintf('No grouping variable %s available in metadata', split.by))
    }
    query <- SplitObject(query, split.by = split.by)
  }
  # Run projection
  query <- make.projection(query = query, ref = ref, filter.cells = filter.cells, ...)
  if(!is.list(query)) query <- list(query = query)
  # Cell type classification
  query <- lapply(query, function(x){
    cellstate.predict(ref=ref, query=x, reduction=reduction, ndim=ndim, k=k, labels.col=labels.col)
  })
  # Merge embeddings
  if(length(query)==1 || !is.null(split.by)) query <- Reduce(merge.Seurat.embeddings, query)
  return(query)
}


#' @title plot.projection
#' @description Show UMAP projection of query on reference map.
#' Plot the UMAP representation of the reference map, together with the projected coordinates of a query dataset.
#' @param ref Reference Atlas
#' @param query Seurat object with query data
#' @param labels.col The metadata field to annotate the clusters (default: functional.cluster)
#' @param cols Custom color palette for clusters
#' @param linesize Contour line thickness for projected query
#' @param pointsize Point size for cells in projected query
#' @param ref.alpha Transparency parameter for reference cells
#' @param ref.size Adjusted point size for reference cells
#' @return UMAP plot of reference map with projected query set in the same space
#' @examples 
#' load('data/query_example_seurat.RData')
#' ref <- load.reference.map()
#' q <- Run.ProjecTME(query_example_seurat, ref = ref, fast.umap.predict = T)
#' plot.projection(ref = ref, query = q)
plot.projection <- function(ref,query=NULL,labels.col='functional.cluster',cols=NULL,linesize=1,
                            pointsize=1, ref.alpha=0.3, ref.size=NULL){
  labels <- ref[[labels.col]][,1]
  states_all <- levels(factor(labels))
  nstates <- length(states_all)
  if(!is.null(cols)){
    if(nstates > length(cols)){
      warning('Not enough colors provided. Making an automatic palette')
      palette <- rainbow(n = nstates)
    }else{ palette <- cols }
  }else{
    if(!is.null(ref@misc$atlas.palette)){
      palette <- ref@misc$atlas.palette
    }else{
      palette <- c("#edbe2a","#A58AFF","#53B400","#F8766D","#00B6EB","#d1cfcc","#FF0000","#87f6a5","#e812dd")
      if(nstates > length(palette)){
        palette <- rainbow(n = nstates)
      }
    }
  }
  # apply transparency to ref cells
  cols_use <- scales::alpha(palette, ref.alpha)
  if(is.null(query)){
    p <- DimPlot(ref,reduction = 'umap',label = F,group.by = labels.col,repel = T,
                 pt.size = ref.size,cols = cols_use) + ggtitle('Reference map') + 
      theme(aspect.ratio = 1)
  }else{
    p <- DimPlot(ref, reduction = 'umap', label = F, group.by = labels.col, 
                 repel = T, pt.size = ref.size, cols = cols_use) +
      geom_point(data.frame(query@reductions$umap@cell.embeddings),
                 mapping = aes(x = UMAP_1, y = UMAP_2),
                 alpha = 0.6, size = pointsize, shape = 17, color = 'gray10') +
      geom_density_2d(data = data.frame(query@reductions$umap@cell.embeddings),
                      mapping = aes(x = UMAP_1, y = UMAP_2),
                      color='black',n=200,h=2,size=linesize) +
      ggtitle('Projection of query on reference map') + theme(aspect.ratio = 1)
  }
  return(p)
}


#' @title plot.statepred.composition
#' @description Summarize the predicted cell states of an object. Makes an barplot of the
#' frequency of cell states in a query object.
#' @param ref Seurat object with the reference object
#' @param query Seurat object with query data
#' @param labels.col The metadata field used to annotate the clusters (default: functional.cluster)
#' @param metric One of 'Count' or 'Percent', 'Count' plots the absolute number of cells,
#' 'Percent' the fraction on the total number of cells.
#' @param cols Custom color palette for clusters
#' @return Barplot of predicted state composition
#' @examples 
#' load('data/query_example_seurat.RData')
#' ref <- load.reference.map()
#' q <- Run.ProjecTME(query_example_seurat, ref = ref, fast.umap.predict = T)
#' plot.statepred.composition(ref = ref, query = q)
plot.statepred.composition <- function(ref,query,labels.col='functional.cluster',
                                       cols=NULL,metric=c('Count','Percent')){
  metric <- tolower(metric[1])
  labels <- ref[[labels.col]][,1]
  states_all <- levels(factor(labels))
  nstates <- length(states_all)
  if(!is.null(cols)){
    if(nstates <= length(cols)){
      palette <- cols[1:nstates]
    }else{
      warning('Not enough colors provided. Making an automatic palette')
      palette <- rainbow(nstates)
    }
  }else{
    if(!is.null(ref@misc$atlas.palette)){
      palette <- ref@misc$atlas.palette
    }else{
      palette <- c("#edbe2a","#A58AFF","#53B400","#F8766D","#00B6EB","#d1cfcc","#FF0000","#87f6a5","#e812dd")
      if(nstates <= length(palette)){
        palette <- palette[1:nstates]
      }else{
        palette <- rainbow(nstates)
      }
    }
  }
  names(palette) <- states_all
  cols_use <- palette[states_all]
  tb <- table(factor(query[[labels.col]][,1],levels = states_all))
  if(metric=='percent'){
    tb <- tb*100/sum(tb)
    tb.m <- reshape2::melt(tb)
    colnames(tb.m) <- c('Cell_state', 'Perc_cells')
    p <- ggplot(tb.m, aes(x = Cell_state, y = Perc_cells, fill = Cell_state)) +
      geom_bar(stat = 'identity') + theme_bw() +
      scale_fill_manual(values = cols_use) +
      theme(axis.text.x = element_blank(), legend.position = 'left')
  }else if(metric=='count'){
    tb.m <- reshape2::melt(tb)
    colnames(tb.m) <- c('Cell_state', 'Ncells')
    p <- ggplot(tb.m, aes(x = Cell_state, y = Ncells, fill = Cell_state)) +
      geom_bar(stat = 'identity') + theme_bw() +
      scale_fill_manual(values = cols_use) +
      theme(axis.text.x = element_blank(), legend.position = 'left')
  }else{
    stop('Unknown metric specified (Must be either "count" or "percent")')
  }
  return(p)
}


#' @title plot.states.radar
#' @description Show expression level of key genes. Make a radar plot of the expression level of a
#' set of genes. It can be useful to compare the gene expression profile of different cell states
#' in the reference atlas vs. a projected set.
#' @param ref Reference Atlas
#' @param query Query data, either as a Seurat object or as a list of Seurat objects
#' @param labels.col The metadata field used to annotate the clusters
#' @param genes4radar Which genes to use for plotting
#' @param meta4radar Which metadata columns (numeric) to use for plotting. If not NULL,
#' `genes4radar` are ignored.
#' @param min.cells Only display cell states with a minimum number of cells
#' @param cols Custom color palette for samples in radar plot
#' @param ref.assay The assay to pull the reference expression data
#' @param query.assay The assay to pull the query expression data
#' @param return.as.list Return plots in a list, instead of combining them in a single plot
#' @return Radar plot of gene expression of key genes by cell subtype
#' @examples 
#' load('data/query_example_seurat.RData')
#' ref <- load.reference.map()
#' plot.states.radar(ref)
#' q <- Run.ProjecTME(query_example_seurat, ref = ref, fast.umap.predict = T)
#' plot.states.radar(ref,q)
plot.states.radar = function(ref, query=NULL,labels.col="functional.cluster",ref.assay='RNA', 
                             query.assay='RNA',genes4radar=NULL,meta4radar=NULL,min.cells=50, 
                             cols=NULL,return=FALSE, return.as.list=FALSE){
  #Make sure query is a list
  if(!is.null(query) & !is.list(query)) query <- list(Query=query)
  #Check assays exist
  if(!ref.assay %in% Assays(ref)){
    stop(sprintf("Assay %s not found in reference object. Please check ref.assay parameter", ref.assay))
  }
  if(is.null(genes4radar)){
    genes4radar <- c("Foxp3","Cd4","Cd8a","Tcf7","Ccr7","Gzmb","Gzmk","Pdcd1","Havcr2","Tox","Mki67")
  }else{
    genes4radar <- genes4radar
  }
  #Whether to use gene expression or metadata
  if(!is.null(meta4radar)){
    refmat <- t(ref[[]])
    feat.use <- intersect(meta4radar, row.names(refmat))
    if(length(feat.use)==0){
      stop("None of the provided meta columns were found - check option 'meta4radar'")
    }
    feat.missing <- setdiff(meta4radar, feat.use)
    if(length(feat.missing)>0){
      to.print <- paste(feat.missing, sep=",", collapse = ",")
      warning(sprintf("Some metadata columns were not found:\n%s", to.print))
    }
  }else{
    refmat <- ref@assays[[ref.assay]]@data
    #Check gene names/feature names
    feat.use <- intersect(genes4radar, row.names(refmat))
    #If overlap is zero, first check whether wrong species was used (upper case to human)
    if(length(feat.use)==0){
      genes4radar <- toupper(genes4radar)
      feat.use <- intersect(genes4radar, row.names(refmat))
      if (length(feat.use)==0) {
        stop("None of the provided genes were found - check option 'genes4radar'")
      }
    }
    feat.missing <- setdiff(genes4radar, feat.use)
    if(length(feat.missing)>0){
      to.print <- paste(feat.missing, sep=",", collapse = ",")
      warning(sprintf("Some gene symbols were not found:\n%s", to.print))
    }
  }
  order <- match(feat.use, row.names(refmat))
  rr <- as.matrix(refmat[order,])
  rr <- matrix(as.numeric(rr), ncol=ncol(rr))
  #Set colors
  ncolors <- 1+length(query)
  if(ncolors==1){
    radar.colors <- "black"
  }else{
    if(is.null(cols)){
      radar.colors <- c("black", scales::hue_pal()(ncolors-1))
    }else{
      cols <- c("black", cols)
      if(ncolors <= length(cols)){
        radar.colors <- cols[1:ncolors]
      }else{  
        warning("Not enough colors provided. Making an automatic palette")
        radar.colors <- c("black", scales::hue_pal()(ncolors-1))
      }
    }
  }  
  names(radar.colors) <- c("Reference", names(query))
  labels <- ref[[labels.col]][,1]
  states_all <- levels(factor(labels))
  nstates <- length(states_all)
  if(!is.null(query)){
    if(is.null(names(query))){
      for(i in 1:length(query)){
        names(query)[[i]] <- paste0("Query",i)
      }
    }
    labels.q <- list()
    qq <- list()
    for(i in 1:length(query)){
      if(!labels.col %in% colnames(query[[i]]@meta.data)){
        message1 <- sprintf("Could not find %s column in query object metadata.",labels.col)
        message2 <- "Did you run cellstate.predict() on this object to predict cell states?"
        stop(paste(message1, message2, sep="\n"))
      }
      labels.q[[i]] <- query[[i]][[labels.col]][,1]
      
      if(!is.null(meta4radar)){
        qmat <- t(query[[i]][[]])
      }else{
        if(!query.assay %in% Assays(query[[i]])){
          stop(sprintf("Assay %s not found in query object. Please check ref.assay parameter", query.assay))
        }
        qmat <- query[[i]]@assays[[query.assay]]@data
      }
      order <- match(feat.use, row.names(qmat))
      qq[[i]] <- as.matrix(qmat[order,])
      qq[[i]] <- matrix(as.numeric(qq[[i]]), ncol=ncol(qq[[i]]))
    }
  }
  
  #Get raw expression means, to normalize by gene
  m <- matrix(, nrow = length(states_all), ncol = length(feat.use))
  rownames(m) <- states_all
  colnames(m) <- feat.use
  for(i in 1:length(states_all)){
    s <- states_all[i]
    m[i,] <- apply(rr[, labels == s], MARGIN=1, function(x){mean(x, na.rm=T)})
  }
  normfacs <- apply(m, MARGIN=2, function(x) {max(c(1,x), na.rm=T)})
  
  pll <- list()
  for(j in 1:length(states_all)){
    s <- states_all[j]
    this.mean <- apply(rr[, labels == s], MARGIN=1, function(x){mean(x, na.rm=T)})
    this.mean <- this.mean/normfacs
    this.df <- data.frame(t(rbind(names(this.mean), this.mean, "Reference")))
    colnames(this.df) <- c("Gene","Expression","Dataset")
    this.df$Expression <- as.numeric(as.character(this.df$Expression))
    i <- 1
    while(i <= length(query)){
      ll <- labels.q[[i]]
      m <- as.matrix(qq[[i]][, !is.na(ll) & ll == s])
      if(ncol(m) >= min.cells){
        q.mean <- apply(m, MARGIN=1, function(x){mean(x, na.rm=T)})
        q.mean <- q.mean/normfacs
        q.df <- data.frame(t(rbind(names(q.mean), q.mean, names(query)[[i]])))
        colnames(q.df) <- c("Gene","Expression","Dataset")
        q.df$Expression <- as.numeric(as.character(q.df$Expression))
        this.df <- rbind(this.df, q.df)
      }
      i=i+1
    }
    ymin <- min(c(-0.1, min(this.df$Expression)))
    ymax <- max(c(1, max(this.df$Expression)))
    levs <- unique(this.df$Dataset)
    this.df$Dataset <- factor(this.df$Dataset, levels=levs)
    this.df$Gene <- factor(this.df$Gene, levels=feat.use)
    pll[[j]] <- ggplot(data=this.df,  aes(x=Gene, y=Expression, group= Dataset, 
                                          colour=Dataset, fill=Dataset)) +
      geom_point(size=2) + geom_polygon(size = 0.75, alpha= 0.1) +
      ylim(ymin, ymax) + ggtitle(s)  + scale_x_discrete() +
      scale_fill_manual(values= radar.colors) +
      scale_colour_manual(values= radar.colors) +
      theme_light() + theme(axis.text.x=element_blank()) +
      annotate(geom="text", x=seq(1,length(feat.use)), y=ymax-0.05*ymax, label=feat.use, size=3) +
      coord_polar()
  }
  #Return plots
  if(return.as.list){
    return(pll)
  }else{
    g <- patchwork::wrap_plots(pll) + patchwork::plot_annotation(paste0("Radar plots for ", labels.col))
    return(g)
  }
}


#' @title find.discriminant.dimensions
#' @description Find discriminant dimensions. Searches PCA or ICA dimensions where the query set deviates
#' the most from a control set or from the reference map. It can be useful to suggest novel cell states
#' that escape from the main axes of diversity of the UMAP.
#' @param ref Seurat object with reference atlas
#' @param query Seurat object with query data
#' @param query.control Optionally, you can compare your query with a control sample, instead of the reference
#' @param query.assay The data slot to be used for enrichment analysis
#' @param state Perform discriminant analysis on this cell state. Can be either:
#' \itemize{
#'   \item{"largest" - Perform analysis on the cell state most represented in the query set(s)}
#'   \item{"all" - Performs analysis on the complete dataset, using all cells}
#'   \item{A specific cell state, one of the states in metadata field labels.col}
#' }
#' @param labels.col The metadata field used to annotate the clusters (default: `functional.cluster`)
#' @param reduction Which dimensionality reduction to use (either ICA or PCA)
#' @param test Which test to perform between the dataset distribution in each ICA/PCA dimension.
#' One of `ks` (Kolmogorov-Smirnov) or `t.test` (T-test)
#' @param ndim How many dimensions to consider in the reduced ICA/PCA space
#' @param print.n The number of top dimensions to return to STDOUT
#' @param verbose Print results to STDOUT
#' @return A dataframe, where rows are ICA/PCA dimensions. ICA/PCA are ranked by statistical significance
#' when comparing their distribution between query and control (or query vs. reference map)
#' @examples 
#' load('data/query_example_seurat.RData')
#' ref <- load.reference.map()
#' q <- Run.ProjecTME(query_example_seurat, ref = ref, fast.umap.predict = T)
#' find.discriminant.dimensions(ref, query = q)
#' find.discriminant.dimensions(ref, query = q, query.control = control.set)
find.discriminant.dimensions <- function(ref,query,query.control=NULL,query.assay='RNA',state='largest',
                                         labels.col='functional.cluster',reduction='ICA',
                                         test=c('ks','t.test'),ndim=50,print.n=3,verbose=T){
  reduction <- tolower(reduction)
  test <- test[1]
  if(is.null(ref)) stop('Please provide the reference object (ref)')
  if(is.null(query)) stop('Please provide a query object (query)')
  # determine cell state for analysis
  if(is.null(state) | state == 'largest'){
    if(!is.null(query.control)){
      ss <- table(rbind(query[[labels.col]], query.control[[labels.col]]))
    }else{
      ss <- table(query[[labels.col]])
    }
    state <- names(sort(ss, decreasing = T))[1]
    message(paste0('Performing discriminant analysis using the most abundant cell state in the query - ',state))
  }else if(state == 'all'){
    message('Performing discriminant analysis using all cells')
  }else{
    if(!state %in% query[[labels.col]][,1]){
      stop(sprintf('State %s not found in query metadata column %s', state, labels.col))
    }
    if(!is.null(query.control) & !state %in% query.control[[labels.col]][,1]){
      stop(sprintf('State %s not found in query.control metadata column %s', state, labels.col))
    }
    message(paste0('Performing discriminant analysis with user-specified state - ', state))
  }
  # subset query data on specific state
  if(state=='all'){
    ref.cells <- seq(1, dim(ref)[2])
    query.cells <- seq(1, dim(query)[2])
    if(!is.null(query.control)) query.c.cells <- seq(1,dim(query.control)[2])
  }else{
    ref.cells <- which(ref[[labels.col]] == state)
    query.cells <- which(query[[labels.col]] == state)
    if(is.null(query.control)) query.c.cells <- which(query.control[[labels.col]] == state)
  }
  if(reduction=='ica'){
    if(is.null(ref@reductions[[reduction]])){
      message('Reduction ICA not found. Calculating ICA for reference object')
      ref <- run.ica(ref, ndim = ndim)
    }
    ref_dimRed <- ref@misc$ica
    perturb_dimRed <- apply.ica.obj(query = query, query.assay = query.assay,ica.obj = ref_dimRed)
    if(!is.null(query.control)){
      control_dimRed <- apply.ica.obj(query = query.control,query.assay = query.assay,ica.obj = ref_dimRed)
    }
  }else if(reduction=='pca'){
    ref_dimRed <- ref@misc$pca_object
    perturb_dimRed <- apply.pca.obj.2(query = query, query.assay = query.assay, pca.obj = ref_dimRed)
    if(is.null(query.control)){
      control_dimRed <- apply.pca.obj.2(query = query.control,query.assay = query.assay,pca.obj = ref_dimRed)
    }
  }else{
    stop(paste0('Unrecognized reduction slot: ',reduction))
  }
  ndim <- min(ndim, length(colnames(perturb_dimRed)))
  df <- data.frame(matrix(ncol = 3, nrow = ndim))
  colnames(df) <- c('stat', 'stat_abs', 'p_val')
  rownames(df) <- colnames(perturb_dimRed)[1:ndim]
  if(!is.null(query.control)){
    message('Query and control datasets was provided. Determining discriminant components of Query vs. Control...')
  }else{
    message('Single query dataset was provided. Determining discriminant components of Query vs. Reference...')
  }
  for(pc in 1:ndim){
    d1 <- perturb_dimRed[query.cells,pc]
    if(!is.null(query.control)){
      d2 <- control_dimRed[query.c.cells,pc]
    }else{
      d2 <- ref@reductions[[reduction]]@cell.embeddings[ref.cells,pc]
    }
    ttest <- t.test(d1, d2)
    if(test=='ks'){
      this.test <- ks.test(d1, d2, alternative='two.sided')
      # KS test statistic has no sign
      this.test.signed <- this.test$statistic * sign(ttest$statistic)
    }else{
      this.test <- ttest
      this.test.signed <- this.test$statistic
    }
    df[pc,'stat'] <- this.test.signed
    df[pc,'stat_abs'] <- abs(this.test.signed)
    df[pc,'p_val'] <- this.test$p.value * ndim
  }
  df <- df[with(df, order(stat_abs, decreasing = T)),]
  buffer <- ''
  for(i in 1:print.n){
    topPC <- rownames(df)[i]
    pc.index <- match(topPC, colnames(perturb_dimRed))
    feats <- ref@reductions[[reduction]]@feature.loadings[,pc.index]
    topgenes <- names(head(sort(abs(feats), decreasing = T), 10))
    topgenes.sign <- feats[topgenes]
    if(df[i,'stat'] > 0){
      topgenes.p <- topgenes.sign[topgenes.sign >= 0]
      topgenes.n <- topgenes.sign[topgenes.sign < 0]
    }else{
      topgenes.p <- topgenes.sign[topgenes.sign < 0]
      topgenes.n <- topgenes.sign[topgenes.sign >= 0]
    }
    pval2print <- ifelse(df[i,'p_val'] < 0.001, sprintf('%.1e',df[i,'p_val']),sprintf('%.4f',df[i,'p_val']))
    buffer <- paste0(buffer,sprintf('-----\nTop %s component %i: %s Stat %.3f p.val %s\n',
                                    reduction, i, topPC, df[i,'stat'], pval2print))
    buffer <- paste0(buffer, 'Driver genes for this component:\n')
    buffer <- paste0(buffer, paste(c('Higher in query +++ ',names(topgenes.p)),collapse = ' '), '\n')
    buffer <- paste0(buffer, paste(c('Lower in query --- ',names(topgenes.n)), collapse = ' '), '\n')
  }
  if(verbose) message(buffer)
  return(df)
}


#' @title plot.discriminant.3d
#' @description 3D plot of reference map with extra discriminant dimension
#' Add an extra dimension to the reference map (it can be suggested by `find.discriminant.dimensions`),
#' to explore addiitonal axes of variability in a query dataset compared to the reference map.
#' @param ref Seurat object with reference atlas
#' @param query Seurat object with query data
#' @param query.control Optionally, you can compare your query with a control sample, instead of the reference
#' @param query.assay The data slot to be used for enrichment analysis
#' @param labels.col The metadata field used to annotate the clusters
#' @param query.state Only plot the query cells from this specific state
#' @param extra.dim The additional dimension to be added on the z-axis of the plot. Can be either:
#' \itemize{
#'   \item{An ICA or PCA dimension (e.g. 'ICA_10')}
#'   \item{Any numeric metadata field associated to the cells (e.g. 'cycling.score')}
#' }
#' @return A three dimensional plot with UMAP_1 and UMAP_2 on the x and y axis respectively, and the
#' specified `extra.dim` on the z-score.
#' @examples 
#' load('data/query_example_seurat.RData')
#' ref <- load.reference.map()
#' q <- Run.ProjecTME(query_example_seurat, ref = ref, fast.umap.predict = T)
#' plot.discriminant.3d(ref, query = q, extra.dim = 'ICA_19')
#' plot.discriminant.3d(ref, query = query.set, query.control=control.set, extra.dim = 'ICA_2')
plot.discriminant.3d <- function(ref,query,query.control=NULL,query.assay='RNA',
                                 labels.col='functional.cluster',
                                 extra.dim='ICA_1', query.state=NULL){
  reduction <- NULL
  message(paste0('Generating UMAP with 3D dimension on ', extra.dim))
  if(grepl('^ica_\\d+', tolower(extra.dim), perl = T)){
    reduction <- 'ica'
  }else if(grepl('^pca_\\d+', tolower(extra.dim), perl = T)){
    reduction <- 'pca'
  }
  ref.3d <- ref
  ref.3d$queryGroup <- 'Reference'
  # Only show cells of a specific state for the query
  if(!is.null(query.state)){
    query.cells <- colnames(query)[which(query[[labels.col]] == query.state)]
    query <- subset(query, cells = query.cells)
    if(!is.null(query.control)){
      query.c.cells <- colnames(query.control)[which(query.control[[labels.col]] == query.state)]
      query.control <- subset(query.control, cells = query.c.cells)
    }
  }
  # Prepare embeddings
  if(!is.null(query.control)){
    q.labs <- c(rep('Control', dim(query.control)[2]), rep('Query', dim(query)[2]))
    q.umaps <- rbind(query.control@reductions$umap@cell.embeddings, query@reductions$umap@cell.embeddings)
    query <- merge(query.control, query)
    query$queryGroup <- q.labs
    query[['umap']] <- CreateDimReducObject(embeddings = q.umaps, key = 'UMAP_', assay = query.assay)
  }else{
    query$queryGroup <- 'Query'
  }
  query@meta.data[,labels.col] <- 'Query'
  metacols <- intersect(names(ref.3d@meta.data), names(query@meta.data))
  ref.3d@meta.data <- rbind(ref.3d@meta.data[,metacols], query@meta.data[,metacols])
  ref.3d@reductions$umap@cell.embeddings <- rbind(ref.3d@reductions$umap@cell.embeddings,query@reductions$umap@cell.embeddings)
  ref.3d@meta.data <- cbind(ref.3d@meta.data, ref.3d@reductions$umap@cell.embeddings)
  # Calculate ICA or PCA embeddings for query
  if(!is.null(reduction)){
    if(reduction=='ica'){
      projected <- apply.ica.obj(query = query, query.assay = query.assay, ica.obj = ref.3d@misc$ica)
      ref.3d@reductions[[reduction]]@cell.embeddings <- rbind(ref.3d@reductions[[reduction]]@cell.embeddings,projected)
    }else if(reduction=='pca'){
      projected <- apply.pca.obj.2(query = query,query.assay = query.assay,pca.obj = ref.3d@misc$pca_object)
      ref.3d@reductions[[reduction]]@cell.embeddings <- rbind(ref.3d@reductions[[reduction]]@cell.embeddings,projected)
    }
  }
  # Add metadata column
  if(is.null(reduction)){
    if(extra.dim %in% colnames(ref.3d@meta.data)){
      ref.3d@meta.data <- cbind(ref.3d@meta.data, Discriminant = ref.3d@meta.data[,extra.dim])
    }else{
      stop(sprintf('extra.dim %s not present in meta.data',extra.dim))
    }
  }else{
    ref.3d@meta.data <- cbind(ref.3d@meta.data,Discriminant=ref.3d@reductions[[reduction]]@cell.embeddings[,extra.dim])
  }
  if(is.null(query.control)){
    cols <- c(Reference = 'gray50', Query = 'red')
  }else{
    cols <- c(Reference = 'gray50', Control = 'green', Query = 'red')
  }
  plotting.data <- ref.3d@meta.data[,c('UMAP_1','UMAP_2','Discriminant',labels.col,'queryGroup')]
  plotting.data$size <- ifelse(plotting.data$queryGroup == 'Reference', 0.3, 6)
  g <- plotly::plot_ly(data = plotting.data, x = ~UMAP_1, y = ~UMAP_2, z = ~Discriminant,
                       color = ~queryGroup, type = 'scatter3d', mode = 'markers',
                       text = ~queryGroup, hoverinfo = 'text', alpha = 0.6, alpha_stroke = 0.6,
                       size = ~size, colors = cols) %>%
    plotly::layout(
      title = sprintf('Projection of query on reference map + dimension %s',extra.dim),
      scene = list(xaxis=list(title='UMAP_1'),yaxis=list(title='UMAP_2'),zaxis=list(title=extra.dim))
    )
  print(g)
  return(g)
}


#' @title find.discriminant.genes
#' @description Find discriminant genes
#' Based on `FindMarkers`, it performs differential expression analysis between a projected query and a
#' control (either the reference map or a control sample) for a given cell type. Useful to detect whether
#' specific cell states over/under-express genes between conditions or with respect to the reference.
#' @param ref Seurat object with reference atlas
#' @param query Seurat object with query data
#' @param query.control Optionally, you can compare your query with a control sample, instead of the reference
#' @param query.assay The data slot to be used for enrichment analysis
#' @param state Perform discriminant analysis on this cell state. Can be either:
#' \itemize{
#'   \item{"largest" - Performs analysis on the cell state most represented in the query set(s)}
#'   \item{"all" - Performs analysis on the complete dataset, using all cells}
#'   \item{A specific cell state, one of the states in metadata field labels.col}
#' }
#' @param labels.col The metadata field used to annotate the cluster (default: functional.cluster)
#' @param test Type of test for DE analysis.
#' @param min.cells Minimum number of cells in the cell type to proceed with anlaysis
#' @param genes.use What subset of genes to consider for DE analysis
#' \itemize{
#'   \item{"variable" - Only consider variable genes of the reference}
#'   \item{"all" - Use inersection of all genes in query and control}
#'   \item{A custom list of genes}
#' }
#' @param ... Adding parameters for `FindMarkers`
#' @return A dataframe with a ranked list of genes as rows, and statistics as columns (e.g. log fold-change, p-values)
#' @examples 
#' # Discriminant genes between query and reference in cell type "Tex"
#' markers <- find.discriminant.genes(ref, query = query.set, state = 'Tex')
#' 
#' # Discriminant genes between query and control sample in most represented cell type
#' markers <- find.discriminant.genes(ref, query = query.set, query.control = control.set)
#' 
#' # Pass results to EnhancedVolcano for visual results
#' library(EnhancedVolcano)
#' EnhancedVolcano(markers, lab = rownames(markers), x = 'avg_log2FC', y = 'p_val')
find.discriminant.genes <- function(ref,query,query.control=NULL,query.assay='RNA',state='largest',
                                    labels.col='functional.cluster', test='wilcox', min.cells=10,
                                    genes.use=c('variable','all'), ...){
  if(is.null(ref)) stop('Please provide the reference object (ref)')
  if(is.null(query)) stop('Please provide a query object (query)')
  # Determine cell state for analysis
  if(is.null(state) | state=='largest'){
    if(!is.null(query.control)){
      ss <- table(rbind(query[[labels.col]],query.control[[labels.col]]))
    }else{
      ss <- table(query[[labels.col]])
    }
    state <- names(sort(ss, decreasing = T))[1]
    message(paste0('Performing differential expression analysis using the most abundant cell state in the query - ',state))
  }else if(state=='all'){
    message('Performing differential expression analysis using all cells')
  }else{
    if(!state %in% query[[labels.col]][,1]){
      stop(sprintf('State %s not found in query metadata column %s', state, labels.col))
    }
    if(!is.null(query.control) & !state %in% query.control[[labels.col]][,1]){
      stop(sprintf('State %s not found in query.control metadata column %s', state, labels.col))
    }
    message(paste0('Performing differential expression analysis with user-specified state - ',state))
  }
  # Subset query data on specific state
  if(state=='all'){
    s1.cells <- colnames(query)
    if(!is.null(query.control)){
      s2.cells <- colnames(query.control)
    }else{
      s2.cells <- colnames(ref)
    }
  }else{
    s1.cells <- colnames(query)[which(query[[labels.col]]==state)]
    if(!is.null(query.control)){
      s2.cells <- colnames(query.control)[which(query.control[[labels.col]]==state)]
    }else{
      s2.cells <- colnames(ref)[which(ref[[labels.col]]==state)]
    }
  }
  # check we have enough cells
  if(length(s1.cells) < min.cells){
    stop(sprintf('Too few cells for state %s in query. Exit.',state))
  }
  if(!is.null(query.control) & length(s2.cells) < min.cells){
    stop(sprintf('Too few cells for state %s in query control. Exit.',state))
  }
  # Subset on subtype
  DefaultAssay(query) <- query.assay
  s1 <- subset(query, cells = s1.cells)
  s1$Group <- 'Query'
  if(!is.null(query.control)){
    DefaultAssay(query.control) <- query.assay
    s2 <- subset(query.control, cells = s2.cells)
  }else{
    DefaultAssay(ref) <- query.assay
    s2 <- subset(ref, cells = s2.cells)
  }
  s2$Group <- 'Control'
  s.m <- merge(s1, s2)
  Idents(s.m) <- 'Group'
  # use all genes or only variable genes from the reference
  if(genes.use[1]=='all'){
    which.genes <- NULL
  }else if(genes.use[1] == 'variable'){
    which.genes <- intersect(ref@assays$integrated@var.features, rownames(s.m))
  }else{
    which.genes <- intersect(genes.use, rownames(s.m))
  }
  markers <- FindMarkers(s.m,slot='data',ident.1 = 'Query', ident.2 = 'Control', only.pos = F,
                         test.use = test, assay = query.assay, features = which.genes, ...)
  return(markers)
}


#' @title make a ProjecTME reference
#' @description Converts a Seurat objec to a ProjecTME reference atlas. You can preserve your 
#' low-dimensionality embeddings (e.g. UMAP) in the reference atlas by setting `recalculate.umap=F`,
#' or recalculate the UMAP using one of the two methods \link{umap::umap} or \link{uwot::umap}.
#' Recalculation allows exploting the 'predict' funcioanlities of these methods for embedding of
#' new points; skipping recalculation will make the projection use an approximation for UMAP 
#' embedding of the query.
#' @param ref Seurat object with reference atlas
#' @param assay The data slot where to pull the expression data
#' @param atlas.name An optional name for your reference
#' @param annotation.column The metadata column with the cluster annotations for this atlas
#' @param recalculate.umap If TRUE, run `umap` or `uwot` algorithm to generate embeddings.
#' Otherwise use the embeddings stored in the `dimred` slot.
#' @param umap.method Which method to use for calculating the umap reduction
#' @param metric Distance metric to use to find nearest neighbors for UMAP
#' @param min_dist Effective minimum distance between UMAP embedded points
#' @param n_neighbors Size of local neighborhood for UMAP
#' @param ndim Number of PCA dimensions
#' @param dimred Use the pre-calculated embeddings stored at `Embeddings(ref,dimred)`
#' @param nfeatures Number of variable features (only calculated if not already present)
#' @param color.palette A (named) vector of colors for the reference plotting functions.
#' One color for each cell type in `functional.cluster`
#' @param scGate.model.human A human scGate model to purify the cell types represented in
#' the map. For example, if the map contains CD4 T cell subtype, specify an scGate model
#' for CD4 T cells.
#' @param scGate.model.mouse A mouse scGate model to purify the cell types represented in the map.
#' @param store.markers Whether to store the top differentially expressed genes in `ref@misc$gene.panel`
#' @param n.markers Store the top `n.markers` for each subtype given by differential expression analysis
#' @param seed Random seed
#' @return A reference atlas compatible with ProjecTME
#' @examples 
#' custom_reference <- make.reference(myref, recalculate.umap = T)
#' colon <- make.reference(colon,atlas.name='crc',annotation.column='majorCluster',recalculate.umap=T)
make.reference <- function(ref,assay=NULL,atlas.name='custom_reference',annotation.column='functional.cluster',
                           recalculate.umap=F,umap.method=c('umap','uwot'),metric='cosine',min_dist=0.3,
                           n_neighbors=30,ndim=20,dimred='umap',nfeatures=1000,color.palette=NULL,
                           scGate.model.human=NULL,scGate.model.mouse=NULL,store.markers=F,
                           n.markers=10,seed=123){
  if(is.null(assay)) assay <- DefaultAssay(ref)
  if(is.null(ref@assays[[assay]])){
    stop(sprintf('Assay %s not found in reference object. Select a different assay',assay))
  }
  DefaultAssay(ref) <- assay
  if(length(ref@reductions)==0){
    print('Perform Dimensionaly Reduction in raw data ...')
    ref <- FindVariableFeatures(ref, assay = assay, nfeatures = nfeatures, verbose=F)
    varfeat <- ref@assays[[assay]]@var.features
    ref <- ScaleData(ref, verbose = F)
    ref <- RunPCA(ref, features = varfeat, verbose = F)
    ref <- FindNeighbors(ref, dims = 1:ndim, reduction = 'pca', verbose = F)
    ref <- RunUMAP(ref, dims = 1:ndim, reduction = 'pca', verbose = F)
  }else{
    varfeat <- ref@assays[[assay]]@var.features
  }
  # Recompute PCA embeddings using prcomp
  print('Recompute PCA embeddings using prcomp ...')
  ref <- prcomp.seurat(ref, ndim = ndim, assay = assay)
  # Recalculate UMAP, or use an already-present dimensionality reduction
  print('Recalculate UMAP, or use an already-present dimensionality reduction ...')
  if(!recalculate.umap){
    if(dimred %in% names(ref@reductions)){
      ref.pca <- ref@misc$pca_object
      cell.order <- rownames(ref.pca$x)
      low.emb <- ref@reductions[[dimred]]@cell.embeddings[cell.order,]
      colnames(low.emb) <- c('UMAP_1', 'UMAP_2')
      # save these embeddings
      ref@misc$umap_object <- list()
      ref@misc$umap_object$data <- ref.pca$x
      ref@misc$umap_object$layout <- low.emb
    }else{
      stop(sprintf('Dimred %s not found in reference object. Select a different dimensionality reduction, 
                   or recalculate.umap=TRUE to compute UMAP coordinates', dimred))
    }
  }else{
    umap.method <- umap.method[1]
    ref.pca <- ref@misc$pca_object
    if(umap.method == 'umap'){
      # generate UMAP embeddings
      print('generate UMAP embeddings ...')
      ref.umap <- run.umap.2(ref.pca,ndim=ndim,seed=seed,n.neighbors=n_neighbors,min.dist=min_dist,metric=metric)
      # save UMAP object
      print('save UMAP object ...')
      ref@misc$umap_object <- ref.umap
      ref@reductions$umap@cell.embeddings <- ref.umap$layout
    }else if(umap.method=='uwot'){
      warning('There are known issues with saving a loading uwot models. If you plan to save your reference as an .rds file,
              please use umap.method = "umap"')
      print('generate UMAP embeddings ...')
      # generate UMAP embeddings
      ref.umap <- run.umap.uwot(ref.pca,ndim=ndim,seed=seed,n.neighbors=n_neighbors,min.dist=min_dist,metric=metric)
      ref.umap$data <- ref.pca$x
      # save UMAP object
      print('save UMAP object ...')
      ref@misc$umap_object <- ref.umap
      ref@reductions$umap@cell.embeddings <- ref.umap$embedding
    }else{
      stop('Unsupported UMAP method.')
    }
  }
  if(!annotation.column=='functional.cluster'){
    ref$functional.cluster <- ref@meta.data[,annotation.column]
  }
  ref$functional.cluster <- factor(ref$functional.cluster)
  lvls <- levels(ref$functional.cluster)
  # Reference name
  ref@misc$ProjecTME <- atlas.name
  # Add color palette
  if(!is.null(color.palette)){
    if(length(color.palette) != length(lvls)) stop('Length of color palette must match number of cell types')
    if(!is.null(names(color.palette))){
      d <- setdiff(names(color.palette), levels(ref$functional.cluster))
      if(length(d)>0){
        stop('Names of color palette do not match annotation levels')
      }
    }else{ names(color.palette) <- lvls }
  }else{
    color.palette <- rainbow(n = length(lvls))
    names(color.palette) <- lvls
  }
  ref@misc$atlas.palette <- color.palette
  # Add scGate models
  ref@misc$scGate <- list()
  if(!is.null(scGate.model.human)) ref@misc$scGate$human <- scGate.model.human
  if(!is.null(scGate.model.mouse)) ref@misc$scGate$mouse <- scGate.model.mouse
  # Add DEGs
  Idents(ref) <- 'functional.cluster'
  if(store.markers){
    DefaultAssay(ref) <- 'RNA'
    print('Find markers for each cluster ...')
    cluster.markers <- FindAllMarkers(
      ref,only.pos = T,min.pct = 0.1,min.diff.pct = 0.1,logfc.threshold = 0.25,
      max.cells.per.ident = 500,test.use = 'wilcox',base = exp(1),verbose = T
    )
    all <- cluster.markers %>% group_by(cluster) %>% top_n(n = n.markers, wt = abs(avg_log2FC))
    markers <- list()
    for(i in levels(ref@active.ident)) markers[[i]] <- subset(all, cluster == i)$gene
    ref@misc$gene.panel <- markers
  }
  # Store in integrated assay, to be understood by ProjecTME
  names(ref@assays)[names(ref@assays)==assay] <- 'integrated'
  DefaultAssay(ref) <- 'integrated'
  return(ref)
}


#' @title Recalculate low dimensional embeddings after projection
#' @description Given a reference object and a (list of) projected objects, recalculate low-dim embeddings
#' accounting for the projected cells.
#' @param ref Reference map
#' @param projected A projected object (or list of projected objects) generated using \link{make.projection}
#' @param ref.assay Assay for reference object
#' @param proj.assay Assay for projected object(s)
#' @param ndim Number of dimensions for recalculating dimensionality reductions
#' @param n.neighbors Number of neighbors for UMAP embedding
#' @param metric Distance metric to use to find nearest neighbors for UMAP
#' @param recalc.pca Whether to recalculate the PCA embeddings with the combined reference and projected data
#' @param umap.method Which method should be used to calculate UMAP embeddings
#' @param resol Resolution for unsupervised clustering
#' @param k.param Number of nearest neighbors for clustering
#' @param seed Random seed for reproducibility
#' @return A combined reference object of reference and projected object(s), with new low dimensional embeddings
#' @examples 
#' load('data/query_example_seurat.RData')
#' ref <- load.reference.map()
#' q <- Run.ProjecTME(query_example_seurat, ref = ref, fast.umap.predict = T)
#' combined <- recalculate.embeddings(ref, q, ndim = 10)
recalculate.embeddings <- function(ref,projected,ref.assay='integrated',proj.assay='integrated',ndim=NULL,
                                   n.neighbors=20,min.dist=0.3,recalc.pca=F,resol=0.4,k.param=15,
                                   metric='cosine',umap.method=c('umap','uwot'),seed=123){
  if(is.null(ref) | is.null(projected)){
    stop('Please provide a reference and a projected object (or list of projected objects)')
  }
  if(is.list(projected) | is.vector(projected)){
    projected <- Reduce(f = merge.Seurat.embeddings, x = projected)
  }
  projected$ref_or_query <- 'query'
  ref$ref_or_query <- 'ref'
  umap.method <- umap.method[1]
  projected <- RenameCells(object = projected, add.cell.id = 'Q')
  merged <- merge.Seurat.embeddings(ref,projected)
  DefaultAssay(merged) <- proj.assay
  DefaultAssay(ref) <- ref.assay
  if(is.null(ndim)){
    ndim <- ncol(ref@misc$umap_object$data)
    if(is.null(ndim)) stop('Please provide number of dimensions for dimensionality reduction (ndim)')
  }
  VariableFeatures(merged) <- VariableFeatures(ref)
  merged@misc <- ref@misc
  merged@misc$pca_object$x <- merged@reductions$pca@cell.embeddings
  if(recalc.pca){
    message('Recalculating PCA embeddings...')
    merged <- prcomp.seurat(merged, assay = proj.assay, ndim = ndim)
  }
  ref.pca <- merged@misc$pca_object
  if(umap.method=='uwot'){
    message('Recalculating UMAP embeddings using uwot...')
    warning('There are known issues with saving a loading uwot models. 
    If you plan to save your reference as an .rds file,please use umap.method = "umap"')
    ref.umap <- run.umap.uwot(ref.pca, ndim = ndim, n.neighbors = n.neighbors, min.dist = min.dist,
                              seed = seed, metric = metric)
    ref.umap$data <- ref.pca$x
    # save UMAP object
    merged@misc$umap_object <- ref.umap
    merged@reductions$umap@cell.embeddings <- ref.umap$embeddings
  }else if(umap.method=='umap'){
    message('Recalculating UMAP embeddings using "umap" package...')
    ref.umap <- run.umap.2(ref.pca,ndim = ndim,n.neighbors = n.neighbors,min.dist = min.dist,seed = seed,metric = metric)
    # save UMAP object
    merged@misc$umap_object <- ref.umap
    merged@reductions$umap@cell.embeddings <- ref.umap$layout
  }else{
    stop('UMAP method not supported.')
  }
  # Did any new clusters arise after adding projected data?
  DefaultAssay(merged) <- ref.assay
  merged <- FindNeighbors(merged, reduction = 'pca', dims = 1:ndim, k.param = k.param, verbose = F)
  merged <- FindClusters(merged, resolution = resol, verbose = F)
  tab <- table(merged$seurat_clusters, merged$ref_or_query)
  glob.freq <- table(merged$ref_or_query)['query']/ncol(merged)
  freq <- apply(tab, 1, function(x){x/sum(x)})['query',] - glob.freq
  freq[freq < 0] <- 0
  merged$newclusters <- freq[merged$seurat_clusters]
  Idents(merged) <- 'funtional.cluster'
  return(merged)
}


#' @title Calculate Silhouette coefficient
#' @description Given a projected object and its reference, calculate silhouette coefficient for 
#' query cells with respect to reference cells with the same cell labels.
#' @param ref Reference object
#' @param query Query object. If not specified, the silhouette coefficient of only the reference will be calculated
#' @param reduction Which dimensionality reduction to use for euclidean distance calculation
#' @param ndim Number of dimensions in the dimred to use for distance calculation. If NULL, use all dimensions.
#' @param label_col Metadata column with cell type annotations. Must be present both in reference and query
#' @param normalize.scores Whether to normalize silhouette scores by the average cell type silhouettes of the reference
#' @param min.cells Only report silhouette scores for cell type with at least this number of cells.
#' @return A dataframe with average silhouette coefficient for each cell type
#' @examples 
#' load('data/query_example_seurat.RData')
#' ref <- load.reference.map()
#' q <- Run.ProjecTME(query_example_seurat, ref = ref, fast.umap.predict = T)
#' sil <- compute.silhouette(ref, query = q)
compute.silhouette <- function(ref,query=NULL,reduction='pca',ndim=NULL,label_col='functional.cluster',
                               normalize.scores=F, min.cells = 20){
  y <- Reductions(ref, slot = reduction)
  if(is.null(ndim)) ndim <- ncol(y@cell.embeddings)
  y <- y@cell.embeddings[, 1:ndim]
  lvls <- levels(as.factor(ref@meta.data[, label_col]))
  labs.y <- ref@meta.data[,label_col]
  if(is.null(query)){
    x <- y; labs.x <- labs.y
  }else{
    subtypes <- table(query@meta.data[, label_col])
    subtypes <- names(subtypes[subtypes > min.cells])
    cid <- Cells(query)[query@meta.data[,label_col] %in% subtypes]
    query <- subset(query, cells = cid)
    x <- Reductions(query, slot = reduction)
    x <- x@cell.embeddings[, 1:ndim]
    labs.x <- factor(query@meta.data[,label_col], levels = lvls)
  }
  dists <- pracma::distmat(x,y)
  sil <- silhouette_2sets(dists, labs.x, labs.y)
  # summarize by cluster
  means <- aggregate(sil$sil_width, list(sil$cluster), mean)
  colnames(means) <- c('Cluster', 'Silhouette')
  if(normalize.scores){
    # normalize by silhouette of the reference
    dist.ref <- pracma::distmat(y,y)
    sil.ref <- silhouette_2sets(dist.ref, labs.y, labs.y)
    means.ref <- aggregate(sil.ref$sil_width, list(sil.ref$cluster), mean)
    colnames(means.ref) <- c('Cluster', 'Silhouette')
    for(i in 1:nrow(means)){
      subset <- means[i, 'Cluster']
      j <- which(means.ref$Cluster == subset)
      means[i,'Silhouette.norm'] <- means[i,'Silhouette']/means.ref[j,'Silhouette']
      if(means[i,'Silhouette.norm'] > 1) means[i,'Silhouette.norm'] = 1
      if(means[i,'Silhouette.norm'] < 0) means[i, 'Silhouette.norm'] = 0
    }
  }
  return(means)
}


#' @title Annotate query dataset using a reference object
#' @description Apply label transfer to annotate a query dataset with the cell types of a reference
#' object. Compared to \link{Run.ProjeTILs}, only cell labels are returned. The low-dim embeddings
#' of the query object (PCA, UMAP) are not modified.
#' @param query Query data stored in a Seurat object
#' @param ref Reference Atlas.
#' @param filter.cells Pre-filter cells using `filterCells`. Only set to FALSE if the dataset has been
#' previously subset to cell types represented in the reference.
#' @param splity.by Grouping variable to split the query object (e.g. if the object contains multiple samples)
#' @param reduction The dimensionality reduction used to assign cell type labels.
#' @param ndim The number of dimensions used for cell type classification
#' @param k Number of neighbors for cell type classification
#' @param labels.col The metadata field with label annotations of the reference, which will be transfered
#' to the query dataset.
#' @param ... Additional parameters to \link{make.projection}
#' @return The query object with an additional metadata column containing predicted cell labels.
#' If cells are filtered prior to projection, they will be labeled as 'NA'.
#' @examples 
#' load('data/query_example_seurat.RData')
#' ref <- load.reference.map()
#' q <- ProjecTME.classifier(query_example_seurat, ref = ref)
#' table(q$functional.cluster, useNA = 'ifany')
#' 
ProjecTME.classifier <- function(query, ref=NULL,filter.cells=F,split.by=NULL,reduction='pca',
                                   ndim=NULL, k=20, labels.col='functional.cluster', ...){
  fast.umap.predict <- TRUE
  # only needed if we want to predict labels based on UMAP neighbors
  if(reduction=='umap') fast.umap.predict <- FALSE
  if(is.list(query)) stop('Query must be a single Seurat object')
  if(!is.null(split.by)){
    if(!split.by %in% colnames(query[[]])){
      stop(sprintf('No grouping variable %s available in metadata', split.by))
    }
    q <- SplitObject(query, split.by = split.by)
  }else{
    q <- query
  }
  q <- make.projection(query = q, ref = ref, filter.cells = filter.cells,
                       fast.umap.predict = fast.umap.predict, ...)
  if(!is.list(q)) q <- list(query = q)
  # cell type classification
  q <- lapply(q, function(x){
    cellstate.predict(ref=ref, query=x, reduction=reduction, ndim=ndim, k=k, labels.col=labels.col)
  })
  # merge embeddings
  q <- Reduce(merge.Seurat.embeddings, q)
  # Transfer labels to original query object
  labs <- q@meta.data[, labels.col]
  names(labs) <- rownames(q@meta.data)
  query@meta.data[,labels.col] <- NA
  query@meta.data[names(labs), labels.col] <- labs
  return(query)
}


#' @title compute.cluster.purity
#' @description Calculate the ROGUE value of each putative cluster for each sample
#' @param query Query data stored in a Seurat object
#' @param labels A vector of cell cluster labels for all cells.
#' @param samples A vector of samples (e.g. patients) to which each cell belongs.
#' @param platform The platform ('UMI' or 'full-length') used for generating the tested dataset.
#' @param k The scaling factor for calculating ROGUE. The default value of k is set to 45 and
#' 500 for droplet-based ('UMI') and 'full-length' based datasets, respectively. When specifying
#' a custom k value, the 'platform' argument is redundant.
#' @param min.cell.n Only clusters containing at least this many cells will receive ROGUE values.
#' @param remove.outlier.n Rmove this many outliers cells when calculating ROGUE
#' @param span The parameter α which controls the dergee of smoothing.
#' @param r A small fixed value to avoid log(0) of mean gene expression levels. The default value of
#' r is set to 1, but can also be set to other values such as 0.1 and 0.01.
#' @param filter Logical, whether to filter out low-abundance genes and low-quality cells.
#' @param min.cells if parameter filter is `TRUE`, include genes detected in at least this many cells.
#' @param min.genes if parameter filter is `TRUE`, include cells where at least this many genes are detected.
#' @param mt.method The multiple testing method used in p.adjsut.
#' @return A dataframe where rows represent samples, cols represent clusters, and values represent corresponding ROGUEs.
#' @examples 
#' load('data/query_example_seurat.RData')
#' ref <- load.reference.map()
#' q <- Run.ProjecTME(query_example_seurat, ref = ref, fast.umap.predict = T)
#' res <- compute.cluster.purity(query = q, labels = q$functional.cluster,
#' samples = q$SampleLab, platform = 'UMI')
#' 
#' dc <- readr::read_rds('data/DC.rds.gz')
#' res2 <- compute.cluster.purity(query = dc, labels = dc$ct, 
#' samples=dc$Patient, platform = 'UMI', span = 0.6)
compute.cluster.purity <- function(query,labels,samples,platform = NULL, k = NULL,min.cell.n = 10,
                                   remove.outlier.n = 2,span = 0.5,r = 1,filter = F,min.cells = 10,
                                   min.genes = 10, mt.method = 'fdr'){
  expr <- as.matrix(query@assays$RNA@counts)
  clusters <- unique(as.character(labels))
  meta <- dplyr::tibble(CellID = 1:ncol(expr), ct = labels, sample = samples)
  sample.rogue <- function(meta,cluster,samples,min.cell.n,filter,platform,k,remove.outlier.n,
                           span, r, min.cells, min.genes, mt.method){
    tmp <- meta %>% dplyr::filter(ct == cluster)
    s <- unique(as.character(samples))
    rogue <- c()
    for(j in 1:length(s)){
      idx1 <- tmp %>% dplyr::filter(sample == s[j]) %>% dplyr::pull(CellID)
      if(length(idx1) >= min.cell.n){
        tmp.matr <- expr[,idx1]
        if(isTRUE(filter)){
          print('Filtering out low-abundance genes and low-quality cells...')
          gene_count <- colSums(tmp.matr > 0, na.rm = T)
          cell_count <- rowSums(tmp.matr > 0, na.rm = T)
          lq1 <- cell_count < min.cells
          lq2 <- gene_count < min.genes
          tmp.matr <- tmp.matr[!lq1, !lq2]
        }else{
          tmp.matr <- tmp.matr
        }
        print('Identify highly informative genes using S-E model...')
        tmp.res <- Entropy(tmp.matr, r = r)
        tmp.res <- entropy_fit(tmp.res,span = span, mt.method = mt.method)
        # tmp.res <- tmp.res %>% dplyr::mutate(p.adj = p.value)
        
        print('S-E plot...')
        tmp.res.plot <- tmp.res %>% dplyr::mutate(sig = ifelse(p.adj <= 0.05, 1, 0))
        p1 <- tmp.res.plot %>% ggplot(aes(mean.expr, entropy)) +
          geom_point(aes(colour = factor(sig)), size = 1) +
          geom_line(aes(mean.expr, fit), lwd = 0.7) +
          scale_color_manual(values = c("#1E90FF", "red")) +
          theme_bw() +
          theme(legend.position = 'top', axis.title = element_text(size=15,color='black'),
                axis.text = element_text(size = 15, colour = 'black'),
                legend.title = element_text(size = 0), 
                legend.text = element_text(size = 0)) +
          labs(x = 'log(mean expression)', y = 'expression entropy')
        ggsave(paste0(cluster,'.SE.plot.pdf'),p1,width = 5,height = 5)
        rm(tmp.res.plot);gc()
        
        print('Remove outlier cells when calculating ROGUE...')
        sig.gene <- tmp.res %>% dplyr::filter(p.adj < 0.05) %>% dplyr::pull(Gene)
        tmp.matr <- tmp.matr[sig.gene,]
        mean.v <- c(); entr.v <- c()
        for(ii in 1:length(sig.gene)){
          x <- sort(as.numeric(tmp.matr[ii,]), decreasing = T)
          x <- x[-c(1:remove.outlier.n)]
          mean.v[ii] <- log(mean(x) + r)
          entr.v[ii] <- mean(log(x + 1))
        }
        mean.cut <- min(tmp.res$mean.expr)
        tmp.res$mean.expr[1:length(sig.gene)] <- mean.v
        tmp.res$entropy[1:length(sig.gene)] <- entr.v
        tmp.res <- tmp.res %>% dplyr::select(-p.adj) %>% dplyr::filter(mean.expr > mean.cut)
        tmp.res <- entropy_fit(tmp.res, span = span, mt.method = mt.method)
        
        print('calculating ROGUE...')
        rogue[j] <- calculateRogue(tmp.res, platform = platform, k = k)
      }else{
        rogue[j] <- NA
      }
    }
    return(rogue)
  }
  res <- list()
  for(i in 1:length(clusters)){
    cat(paste0('Processing for cluster: ',clusters[i],'\n'))
    res[[i]] <- sample.rogue(
      meta = meta, cluster = clusters[i],samples = samples,min.cell.n = min.cell.n,filter = filter,
      platform = platform, k = k, remove.outlier.n = remove.outlier.n, span = span, r = r,
      min.cells = min.cells, min.genes = min.genes, mt.method = mt.method
    )
  }
  res.tib <- Reduce(rbind,res) %>% as.matrix() %>% t() %>% as.data.frame()
  colnames(res.tib) <- clusters
  rownames(res.tib) <- unique(as.character(samples))
  
  cat('Visualize ROGUE values on a boxplot...')
  p <- res.tib %>% tidyr::gather(key = clusters, value = ROGUE) %>%
    ggplot(aes(clusters, ROGUE)) +
    geom_boxplot(color = '#FF3E96', outlier.shape = NA) +
    geom_point(color = '#FF3E96', size = 1.5) +
    theme_bw() +
    theme(axis.text = element_text(size = 12, colour = 'black'),
          axis.title = element_text(size = 13, colour = 'black')) +
    labs(x = 'Clusters', y = 'ROGUE')
  ggsave('ROGUE.pdf',p,width = 6,height = 4)
  
  return(res.tib)
}


#' @title Accurate and fast cell marker gene identification
#' @description Marker gene identification for cell groups in a given dataset
#' @param object Seurat object
#' @param groups Groups to perform comparation. Default: All.
#' @param assay Assay to use in marker gene identification. Default: `RNA` assay
#' @param slot Slot to pull data from. Default: `data` slot
#' @param mu The penalty fcator to penalize gene expression in cells not belonging to the cluster of interest. Default: 1
#' @param remove_lowly_expressed Logical to determine whether to remove lowly expressed genes. Default: TRUE
#' @param expressed_pct Numeric value to set percentage of expressed genes. Default: 0.1
#' @param n_top_genes Number of top ranked genes returned in the result. Default: 100
#' @return A list containing two dataframes for ranked marker genes' names and scores, respectively.
#' @examples 
#' library(Seurat)
#' # Case 1:
#' data("pbmc_small")
#' # all groups
#' markers <- ProjecTME.find.markers(pbmc_small,groups='all',assay='RNA',slot='data',mu=1,n_top_genes=100)
#' # group: 0, 2
#' markers2 <- ProjecTME.find.markers(pbmc_small,groups=c('0','2'),assay='RNA',slot='data',mu=1,n_top_genes=100)
#' head(markers$names)
#' head(markers$scores)
#' top_list <- c()
#' for(group in colnames(markers$names)){
#'   top_i <- markers$names[group][1:5,1]
#'   top_list <- c(top_list,top_i)
#' }
#' DotPlot(pbmc_small, assay = 'RNA', features = unique(top_list)) + RotatedAxis()
#' DoHeatmap(pbmc_small, assay = 'RNA', features = top_list)
#' 
#' # Case 2:
#' if(!require(dittoSeq)) BiocManager::install('dittoSeq')
#' load('data/query_example_seurat.RData')
#' ref <- load.reference.map()
#' projected <- Run.ProjecTME(query_example_seurat, ref = ref, fast.umap.predict = T)
#' Idents(projected) <- 'functional.cluster'
#' markers <- ProjecTME.find.markers(projected,groups=c('CD8_Tex','CD8_Tpex'),assay='RNA',slot='data',mu=1,n_top_genes=100)
#' markers.name <- unique(markers$names$CD8_Tex, markers$names$CD8_Tpex)
#' top_list <- c()
#' for(group in colnames(markers$names)){
#'   top_i <- markers$names[group][1:10,1]
#'   top_list <- c(top_list,top_i)
#' }
#' htmap.colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdYlBu')))(50)
#' dittoHeatmap(projected,genes = markers.name,annot.by = c('functional.cluster','Sample','SampleLab','nFeature_RNA'),
#'   scale.to.max = F, treeheight_row = 10, heatmap.colors = htmap.colors,
#'   show_rownames = F, highlight.features = top_list)
ProjecTME.find.markers <- function(object,groups='all',assay='RNA',slot='data',mu=1,
                                    remove_lowly_expressed = T, expressed_pct = 0.1,
                                    n_top_genes = 100){
  if(!require(proxyC)) install.packages('proxyC')
  # Obtain the cellxgene data
  genexcell <- Seurat::GetAssayData(object = object[[assay]], slot = slot)
  if(length(groups) == 1){
    if(groups == 'all'){
      group_info <- Seurat::Idents(object = object)
    }else{
      stop('Cannot perform marker gene identification on a single cluster.')
    }
  }else{
    object <- subset(object, idents = groups)
    group_info <- Seurat::Idents(object = object)
  }
  # unique groups
  groups_order <- sort(unique(group_info))
  n_cluster <- length(groups_order)
  if(n_cluster==1) stop('Cannot perform marker gene identification on a single cluster.')
  n_cell <- ncol(genexcell)
  n_gene <- nrow(genexcell)
  gene_name <- rownames(genexcell)
  # If specifying too many genes to return
  if(n_top_genes > n_gene) n_top_genes <- n_gene
  cluster_mat <- matrix(0, nrow = n_cluster, ncol = n_cell)
  order_i <- 1
  # Set gene lambda and gene omega
  for(group in groups_order){
    idx <- group_info == group
    cluster_mat[order_i, idx] <- 1
    order_i <- order_i + 1
  }
  cluster_mat_sparse <- as(cluster_mat, 'dgCMatrix')
  # Calculate the cosine similarity
  cosine_sim <- proxyC::simil(genexcell, cluster_mat_sparse, method = 'cosine', drop0 = T)
  pos_nonzero <- cosine_sim != 0
  pos_nonzero <- which(as.matrix(pos_nonzero), arr.ind = T)
  # Second stage
  genexlambda <- cosine_sim * cosine_sim
  e_power2_sum <- Matrix::rowSums(genexlambda)
  if(mu==1){
    genexlambda[pos_nonzero] <- genexlambda[pos_nonzero]/(replicate(ncol(genexlambda),e_power2_sum)[as.matrix(pos_nonzero)])
  }else{
    genexlambda[pos_nonzero] <- genexlambda[pos_nonzero]/(((1-mu)*genexlambda[pos_nonzero]+mu*(replicate(ncol(genexlambda),e_power2_sum)[as.matrix(pos_nonzero)])))
  }
  genexlambda <- genexlambda * cosine_sim
  rank_stats_names <- data.frame(
    matrix(matrix(), n_top_genes, length(groups_order), dimnames = list(seq(1,n_top_genes),groups_order)),
    stringsAsFactors = F
  )
  rank_stats_scores <- data.frame(
    matrix(matrix(),n_top_genes,length(groups_order),dimnames = list(seq(1,n_top_genes),groups_order)),
    stringsAsFactors = F
  )
  order_i <- 1
  # Set gene lambda and gene omega
  for(group in groups_order){
    idx <- group_info == group
    scores <- genexlambda[,order_i]
    # Mask these genes expressed in less than given percentage of cells in the cluster of interest
    if(remove_lowly_expressed){
      # https://stackoverflow.com/questions/51560456/r-package-matrix-get-number-of-non-zero-entries-per-rows-columns-of-a-sparse
      n_cells_expressed <- tabulate(genexcell[,idx]@i + 1)
      n_cells <- sum(idx)
      scores[n_cells_expressed < n_cells * expressed_pct] <- -1
    }
    dd <- data.frame(x = data.table::copy(scores), indice = seq(1, length(scores)))
    data.table::setDT(dd)
    data.table::setorder(dd, -x)
    global_indices <- dd$indice[1:n_top_genes]
    rank_stats_names[,order_i] <- gene_name[global_indices]
    rank_stats_scores[,order_i] <- scores[global_indices]
    # save the group names
    order_i <- order_i + 1
  }
  colnames(rank_stats_names) <- groups_order
  colnames(rank_stats_scores) <- groups_order
  ranks_stats <- list(names = rank_stats_names, scores = rank_stats_scores)
  return(ranks_stats)
}


#' @title Functional enrichment analysis for DE markers
#' @param markers Output of `ProjecTME.find.markers`
#' @param top_genes Number of top ranked genes returned in the result. Default: 100
#' @param project Suffix of the output
#' @examples 
#' library(Seurat)
#' data("pbmc_small")
#' # group: 0, 2
#' markers <- ProjecTME.find.markers(pbmc_small,groups=c('0','2'),assay='RNA',slot='data',mu=1,n_top_genes=100)
#' ProjecTME.enrichment(markers = markers, top_genes = 100, project = 'pbmc')
ProjecTME.enrichment <- function(markers,species = 'human',top_genes = 100, project){
  pacman::p_load(clusterProfiler,org.Hs.eg.db,org.Mm.eg.db,ReactomePA,ggplot2,dplyr)
  symbols_list <- as.list(as.data.frame(apply(markers$names,2,head,top_genes)))
  print('Convert gene symbols to ENTREZID...')
  gcSample <- lapply(symbols_list, function(y){
    if(species=='human'){
      as.character(na.omit(select(org.Hs.eg.db,keys=y,columns='ENTREZID',keytype='SYMBOL')[,2]))
    }else if(species=='mouse'){
      as.character(na.omit(select(org.Mm.eg.db,keys=y,columns='ENTREZID',keytype='SYMBOL')[,2]))
    }else{
      stop('Currently only "human" and "mouse" are supported!')
    }
  })
  print('KEGG enrichment analysis...')
  if(species=='human'){
    xx <- compareCluster(gcSample,fun = 'enrichKEGG',organism = 'hsa', pvalueCutoff = 0.05)
  }else if(species=='mouse'){
    xx <- compareCluster(gcSample,fun = 'enrichKEGG',organism = 'mmu', pvalueCutoff = 0.05)
  }else{
    stop('Currently only "human" and "mouse" are supported!')
  }
  p1 <- dotplot(xx) + scale_y_discrete(labels=function(x) stringr::str_wrap(x,width = 50))
  ggsave(paste0(project,'_kegg.pdf'),p1,width = 10,height = 8)
  
  print('ReactomePA enrichment analysis...')
  if(species=='human'){
    xx <- compareCluster(gcSample,fun = 'enrichPathway',organism = 'human', pvalueCutoff = 0.05)
  }else if(species=='mouse'){
    xx <- compareCluster(gcSample,fun = 'enrichPathway',organism = 'mouse', pvalueCutoff = 0.05)
  }else{
    stop('Currently only "human" and "mouse" are supported!')
  }
  p2 <- dotplot(xx) + scale_y_discrete(labels=function(x) stringr::str_wrap(x,width = 50))
  ggsave(paste0(project,'_ReactomePA.pdf'),p2,width = 10,height = 8)
  
  print('GO enrichment analysis for: BP...')
  if(species=='human'){
    formula_res <- compareCluster(gcSample, fun = 'enrichGO', OrgDb = 'org.Hs.eg.db', ont = 'BP',
                                  pAdjustMethod = 'BH', qvalueCutoff = 0.05)
  }else if(species=='mouse'){
    formula_res <- compareCluster(gcSample, fun = 'enrichGO', OrgDb = 'org.Mm.eg.db', ont = 'BP',
                                  pAdjustMethod = 'BH', qvalueCutoff = 0.05)
  }else{
    stop('Currently only "human" and "mouse" are supported!')
  }
  cat(paste0('test and merge terms that are close to each other to remove result redundancy'))
  lin_ego <- simplify(formula_res, cutoff = 0.5, by = 'p.adjust', select_fun = min)
  pdf(paste0(project,'_GO_BP_cluster_simplified.pdf'),width = 15,height = 8)
  print(dotplot(lin_ego,showCategory=5) + scale_y_discrete(labels=function(x) stringr::str_wrap(x,width = 50)))
  dev.off()
  write.csv(lin_ego@compareClusterResult, file = paste0(project,'_GO_BP_cluster_simplified.csv'))
  
  print('GO enrichment analysis for: CC...')
  if(species=='human'){
    formula_res <- compareCluster(gcSample, fun = 'enrichGO', OrgDb = 'org.Hs.eg.db', ont = 'CC',
                                  pAdjustMethod = 'BH', qvalueCutoff = 0.05)
  }else if(species=='mouse'){
    formula_res <- compareCluster(gcSample, fun = 'enrichGO', OrgDb = 'org.Mm.eg.db', ont = 'CC',
                                  pAdjustMethod = 'BH', qvalueCutoff = 0.05)
  }else{
    stop('Currently only "human" and "mouse" are supported!')
  }
  cat(paste0('test and merge terms that are close to each other to remove result redundancy'))
  lin_ego <- simplify(formula_res, cutoff = 0.5, by = 'p.adjust', select_fun = min)
  pdf(paste0(project,'_GO_CC_cluster_simplified.pdf'),width = 15,height = 8)
  print(dotplot(lin_ego,showCategory=5) + scale_y_discrete(labels=function(x) stringr::str_wrap(x,width = 50)))
  dev.off()
  write.csv(lin_ego@compareClusterResult, file = paste0(project,'_GO_CC_cluster_simplified.csv'))
  
  print('GO enrichment analysis for: MF...')
  if(species=='human'){
    formula_res <- compareCluster(gcSample, fun = 'enrichGO', OrgDb = 'org.Hs.eg.db', ont = 'MF',
                                  pAdjustMethod = 'BH', qvalueCutoff = 0.05)
  }else if(species=='mouse'){
    formula_res <- compareCluster(gcSample, fun = 'enrichGO', OrgDb = 'org.Mm.eg.db', ont = 'MF',
                                  pAdjustMethod = 'BH', qvalueCutoff = 0.05)
  }else{
    stop('Currently only "human" and "mouse" are supported!')
  }
  cat(paste0('test and merge terms that are close to each other to remove result redundancy'))
  lin_ego <- simplify(formula_res, cutoff = 0.5, by = 'p.adjust', select_fun = min)
  pdf(paste0(project,'_GO_MF_cluster_simplified.pdf'),width = 15,height = 8)
  print(dotplot(lin_ego,showCategory=5) + scale_y_discrete(labels=function(x) stringr::str_wrap(x,width = 50)))
  dev.off()
  write.csv(lin_ego@compareClusterResult, file = paste0(project,'_GO_MF_cluster_simplified.csv'))
}


#' @title calculate Pathway Activity Score
#'
#' @param object Seurat object
#' @param method calculation method
#' @param gmt_file pathways in GMT format
#' @param species species, human or mouse, default by human
#' @param pathway abbreviaiton for pathway database
#' @param filter whether filtering for genes expressed in less than 5 percent cells
#' @param normalize normalization method
#' @examples 
#' library(Seurat)
#' data('pbmc_small')
#' seu <- ProjecTME.calculate_PAS(object = pbmc_small, method = 'AUCell', normalize='log',species='human',pathway='kegg')
#' markers2 <- ProjecTME.find.markers(seu,groups=c('0','2'),assay='PAS',slot='data',mu=1,n_top_genes=100)
#' top_list <- c()
#' for(group in colnames(markers$names)){
#'   top_i <- markers$names[group][1:5,1]
#'   top_list <- c(top_list,top_i)
#' }
#' DotPlot(pbmc_small, assay = 'PAS', features = unique(top_list)) + RotatedAxis()
#' DoHeatmap(pbmc_small, assay = 'PAS', features = top_list, slot = 'data')
ProjecTME.calculate_PAS <- function(object,method=c('AUCell','VISION','GSVA','ssGSEA','plage','zscore'),
                                     species='human',pathway='none',gmt_file='none',
                                     normalize='sctransform',n_cores=3,rand_seed=123){
  if(gmt_file=='none'){
    gSets.path <- paste0('reference/gmtFiles/',species,'/',pathway,'.gmt')
  }else{
    if(species != 'none' || pathway != 'none'){
      warning(" 'gmt_file' is already present. Ignore 'species' and 'pathways'.")
    }
    gSets.path <- gmt_file
  }
  if(method %in% c('pagoda2','VISION')){
    warning("'log' transform will counter error for 'pagoda2' or 'VISION'.
            force the 'normalize' to 'none' or you could modify your code
            to call other normalization function.\n")
    normalize <- 'none'
  }
  score <- cal_all_tools(
    object = object, gSets_path = gSets.path, method = method,
    normalize = normalize, n_cores = n_cores, rand_seed = rand_seed
  )
  PAS <- CreateAssayObject(score[[method]])
  object[['PAS']] <- PAS
  DefaultAssay(object) <- 'PAS'
  return(object)
}


#' @title ProjecTME.pseudotime
#' @description Calculating pseudotime trajectory and pathway analysis
#' @param seu Seurat object obtained from `Run.ProjecTME()`
#' @param n_var Top variable genes, default: 1000.
#' @param ti_method One or mor methods. Must be one of:
#' 1) an object or list of ti_... objects (eg. `dymethods::ti_comp1()`).
#' 2) a character vector containing the names of methods to execute (eg. "scorpius").
#' 3) a character vector containing dockerhub repositories (eg. `dynverse/paga`), or
#' 4) a dynguidelines data frame.
#' @param root Logical value to determine whether infer trajectory biologically. Default: FALSE.
#' @param root.markers Only run when `root = T`.
#' @param milestone_labelling Logical value to perform milestone labelling analysis.
#' @param milestone_marker_list Only run when `milestone_labelling = T`
#' @param add_dr Logical value to determine whether to add dimensionality reduction analysis
#' peformed by `dyndimred::dimred_mds`. Default: TRUE.
#' @param features Features of interest to visualize. Default: NULL.
#' @param interest_features Interating features to draw heatmap. If NULL, plot top 50 features.
#' Default: NULL.
#' @param milestone_celltype Only run when `milestone_labelling = T`.
#' @param cut_point Intereting branch cut point. Default: 2.
#' @param gmt File path for pathway in gmt format.
#' @param highlight.pathways Selecting interesting pathway for highlighting.
#' @param show.row.names Logical value to determine whether to draw all pathway names in heatmap.
#' Default: FALSE.
#' @param base.point.size Default: 2.5.
#' @param highlight.point.size Default: 3.
#' 
#' @examples 
#' load('data/query_example_seurat.RData')
#' ref <- load.reference.map()
#' projected <- Run.ProjecTME(query_example_seurat, ref = ref, fast.umap.predict = T)
#' ProjecTME.pseudotime(
#'  seu = projected, n_var = 1000, ti_method = 'slingshot', root = F, root.markers = NULL,
#'  milestone_labelling = F, milestone_marker_list = NULL, add_dr = T, features = 'Gzmf',
#'  interest_features = NULL, cut_point = 2, gmt = 'reference/gmtFiles/human/hallmarker.gmt',
#'  highlight.pathways = 'mtorc', show.row.names = F, base.point.size = 2.5,
#'  highlight.point.size = 3
#' )
ProjecTME.pseudotime <- function(seu, n_var = 1000, ti_method, root = F, root.markers = NULL,
                                 milestone_labelling = F, milestone_marker_list = NULL,
                                 add_dr = T, features = NULL, interest_features = NULL,
                                 milestone_celltype, cut_point = 2, gmt, highlight.pathways,
                                 show.row.names = F,
                                 base.point.size = 2.5, highlight.point.size = 3){
  if(!require(dyno)) devtools::install_github("dynverse/dyno")
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  df <- as.matrix(seu@assays$RNA@data)
  var_genes <- names(sort(apply(df, 1, var),decreasing = T))[1:n_var]
  count <- Matrix::t(as(as.matrix(seu@assays$RNA@counts[var_genes,]),'sparseMatrix'))
  expression <- Matrix::t(as(as.matrix(seu@assays$RNA@data[var_genes,]),'sparseMatrix'))
  dataset <- wrap_expression(expression = expression, counts = count)
  print('Building trajectory model...')
  model <- infer_trajectory(dataset, method = ti_method, verbose = T)
  if(root == T){
    print('Interpreting the trajectory biologically...')
    if(!is.null(root.markers)){
      model <- model %>% add_root_using_expression(root.markers, dataset$expression)
    }else{
      stop('Expect features as root.markers, e.g., "root.markers = \'Vim\'"')
    }
  }
  if(milestone_labelling == T){
    print('Milestone labelling...')
    if(!is.null(milestone_marker_list)){
      model <- label_milestones_markers(model,markers = milestone_marker_list, dataset$expression)
    }else{
      stop(
        'Expect gene lists as milestone_marker_list, 
        e.g., "milestone_marker_list = list(MEF = c(\'Vim\'), Myocyte = c(\'Myl1\'), Neuron = c(\'Stmn3\'))"'
      )
    }
  }
  
  print('Generating pseudotime trajectory plot...')
  p1 <- plot_dimred(
    model, 'pseudotime', pseudotime = calculate_pseudotime(model),
    hex_cells = F, plot_trajectory = T, size_cells = 1, alpha_cells = 0.8
  ) + theme(aspect.ratio = 1)
  
  print('Generating cell cluster plot...')
  p2 <- plot_dimred(model, grouping = seu$functional.cluster, hex_cells = F,
                    plot_trajectory = T, size_cells = 1, alpha_cells = 0.8) + theme(aspect.ratio = 1)
  print('Generating milestone grouping plot...')
  p3 <- plot_dimred(model, grouping = group_onto_nearest_milestones(model), hex_cells = F,
                    plot_trajectory = T, size_cells = 1, alpha_cells = 0.8) + theme(aspect.ratio = 1)
  p <- cowplot::plot_grid(plotlist = list(p1,p2,p3),nrow = 1)
  ggsave(filename = 'pseudo.pdf',p,width = 12,height = 4)
  
  if(add_dr == TRUE){
    model <- model %>% add_dimred(dyndimred::dimred_mds, expression_source = dataset$expression)
    p1 <- plot_dimred(model, 'pseudotime', pseudotime = calculate_pseudotime(model),
                      hex_cells = F, plot_trajectory = T, size_cells = 1, alpha_cells = 0.8) + 
      theme(aspect.ratio = 1)
    p2 <- plot_dimred(model, expression_source = dataset$expression,grouping = projected$functional.cluster)
    p3 <- plot_dimred(model, grouping = group_onto_nearest_milestones(model), hex_cells = F,
                      plot_trajectory = T, size_cells = 1, alpha_cells = 0.8) + theme(aspect.ratio = 1)
    p <- cowplot::plot_grid(plotlist = list(p1,p2,p3),nrow = 1)
    ggsave(filename = 'pseudo.mds.pdf',p,width = 12,height = 4)
  }
  print('Saving model...')
  saveRDS(model, file = 'pseudo.model.rds')
  
  if(!is.null(features)){
    p.gene <- plot_dimred(model, expression_source = dataset$expression, feature_oi = features)
    p.g.genes <- plot_dimred(
      model, expression_source = dataset$expression, color_cells = 'feature',
      feature_oi = features, color_density = 'grouping',
      grouping = seu$functional.cluster, label_milestones = milestone_labelling
    )
    p <- cowplot::plot_grid(plotlist = list(p.gene, p.g.genes), nrow = 1)
    ggsave('features_of_interest.pdf',p,width = 8,height = 4)
  }else{ next }
  
  print('Predicting and visualising genes of interest...')
  # At default, the overall most important genes are calculated when plotting a heatmap.
  pdf('heatmap_of_predictive_genes.pdf',width = 12,height = 8)
  if(is.null(interest_features)){
    message('At default, visualize top 50 most important genes...')
    p1 <- plot_heatmap(model,expression_source = dataset$expression,
                       grouping = seu$functional.cluster,features_oi = 50)
  }else{
    p1 <- plot_heatmap(model, expression_source = dataset$expression,
                       grouping = seu$functional.cluster, features_oi = interest_features)
  }
  print(p1)
  dev.off()
  
  print('Lineage/branch markers...')
  print('Extract features specific for a branch ...')
  branch_feature_importance <- calculate_branch_feature_importance(model,expression_source = dataset$expression)
  if(milestone_labelling == T){
    if(!is.null(milestone_celltype)){
      features <- branch_feature_importance %>%
        filter(to == which(model$milestone_labelling == milestone_celltype)) %>%
        top_n(50, importance) %>% pull(feature_id)
    }else{
      stop('Please provide celltype markers: e.g. 
           milestone_marker_list <- list(MEF = c(\'Vim\'), Myocyte = c(\'Myl1\'), Neuron = c(\'Stmn3\'))')
    }
  }else{
    features <- branch_feature_importance %>%
      filter(to == which(model$milestone_ids == cut_point)) %>%
      top_n(50, importance) %>% pull(feature_id)
  }
  pdf('heatmap_of_predictive_genes.2.pdf',width = 12,height = 8)
  p2 <- plot_heatmap(model,expression_source = dataset$expression,features_oi = features)
  print(p2)
  dev.off()
  
  print('Genes important at bifurcation points...')
  branching_milestone <- model$milestone_network %>% group_by(from) %>% pull(from)
  branching_milestone <- branching_milestone[cut_point]
  branch_feature_importance <- calculate_branching_point_feature_importance(
    model, expression_source = dataset$expression,
    milestones_oi = branching_milestone
  )
  branching_point_features <- branch_feature_importance %>% top_n(20,importance) %>% pull(feature_id)
  pdf('branch.point.features.pdf',p,width = 12,height = 8)
  p <- plot_heatmap(model,expression_source = dataset$expression,features_oi = branching_point_features)
  print(p)
  dev.off()
  
  space <- dyndimred::dimred_mds(dataset$expression)
  plist <- map(branching_point_features[1:12], function(feature_oi){
    plot_dimred(model, dimred = space, expression_source = dataset$expression, 
                feature_oi = feature_oi, label_milestones = F) +
      theme(legend.position = 'none') + ggtitle(feature_oi)
  })
  p <- patchwork::wrap_plots(plist)
  ggsave('features.plot.pdf',p,width = 9,height = 6)
  
  print('Comparing pathway over pseudotime...')
  # extract the cells based on this grouping
  mile_group <- data.frame(group_onto_nearest_milestones(model)) %>% 
    magrittr::set_colnames('milestone') %>% rownames_to_column('cell')
  seu$milestone <- as.character(mile_group$milestone)
  # Extract expression data for the populations we're comparing
  pseudo <- list()
  for(i in unique(seu$milestone)){
    pseudo[[i]] <- seurat_extract(seu, meta1 = 'milestone', value_meta1 = i)
  }
  print('Define pathways and run comparison...')
  if(is.null(gmt)){
    stop('Please provide file path for gmt...')
  }else{
    gsets <- getGMT(gmt)
  }
  pathway <- list()
  for(i in 1:length(gsets)){
    pathway[[i]] <- data.frame(Pathway = rep(names(gsets[i]),length(gsets[i])), Genes = gsets[[i]])
  }
  
  metabolism <- compare_pathways(samples = pseudo, pathways = pathway)
  saveRDS(metabolism, file = 'metabolism.rds')
  print('Generating pathway heatmap...')
  pdf('global_summary_pathway_heatmap.pdf',width = 4,height = 6)
  generate_heatmap(out = metabolism, highlight.pathways = highlight.pathways,show.row.names = show.row.names)
  dev.off()
  print('Highlight pathway by rank...')
  p <- plot_rank(out = metabolism, pathway = highlight.pathways, base.point.size = base.point.size,highlight.point.size = highlight.point.size)
  ggsave('plot.rank.pdf',p,width = 4,height = 5)
}


#' @param seu Seurat object generated from `Run.ProjecTME()`.
#' @param signature_collection Signature collection prepared for generating geneset, 
#' either PSc (Drug Perturbation Signature Collection), or SSc (Drug Sensitivity Signature Collection).
#' @param include.pathways Logical, return \code{beyondcell}'s pre-computed signatures for functional pathways.
#' @param expr.thres Minimum fraction of signature genes that must be expressed in a cell to compute its BCS.
#' Cells with a number of expressed genes below this fraction will have a \code{NaN} BCS.
#' @param pc Number of principal components to use (which can be estimated from an elbow plot). Default: 10.
#' @param k.neighbors (\code{\link[Seurat]{FindNeighbors}}' \code{k.param}). Defines \emph{k} for the k-Nearest
#' Neighbour algorithm.
#' @param res (\code{\link[Seurat]{FindClusters}}' \code{resolution}) Value of the resolution parameter, 
#' use a value above/below 1.0 if you want to obtain a larger/smaller number of communities. 
#' Can be a single number of a numeric vector.
#' @param method (\code{\link[Seurat]{RunUMAP}}'s \code{umap.method}) UMAP implementation to run. Can be:
#' \itemize{
#'   \item{\code{uwot}}: Runs the UMAP via the \code{\link[uwot]{umap}} function of the \code{uwot} package.
#'   \item{\code{umap-learn}}: Runs the UMAP via the \code{Seurat} wrapper of the python \code{umap-learn} package.
#' }
#' @param UMAP UMAP reduction to plot. Either \code{beyondcell}, computed using \code{\link[beyondcell]{bcUMAP}},
#' or \code{Seurat}, obtained using \code{Seurat}'s functions.
#' @param isRegressOut Logical value determine whether perform regression.
#' @param vars.to.regress Character stored in metadata to regress out.
#' @param idents Name of the metadata column of interest. I \code{idents = NULL}, the function computes the ranks
#' using all cells. If \code{idents != NULL}, the signatures' ranks are computed for each level in \code{idents}.
#' @param drugs Name of drugs to perform `FindDrugs()`.
#' @param signatures Name of signature to perform `bcSignatures()`
#' @param genes Name of signature to perform `bcSignatures()`
#' @param subcluster Character vector with the \code{idents}' level(s) of interest. 
#' @param topgenes Number of top drugs per quadrant to be labelled.
#' 
#' @examples 
#' if(!file.exists('data/JerbyArnon.rds')){
#'   # Make reference
#'   if(!file.exists('reference/ref_TICAtlas_human_downsampled.rds')){
#'     tica <- readRDS('reference/TICAtlas_downsampled.rds')
#'     ref <- make.reference(tica,atlas.name = 'TICA',annotation.column = 'lv1_annot',recalculate.umap = F)
#'     saveRDS(ref, file = 'reference/ref_TICAtlas_human_downsampled.rds')
#'   }else{
#'     ref <- readRDS('reference/ref_TICAtlas_human_downsampled.rds')
#'   }
#'   test <- readRDS('data/testing.datasets.2k.rds')
#'   JerbyArnon <- test$JerbyArnon
#'   projected <- Run.ProjecTME(JerbyArnon, ref = ref, fast.umap.predict = T, filter.cells = F)
#'   unique(projected$functional.cluster)
#'   old.idents <- unique(projected$functional.cluster)
#'   new.idents <- c('Bcells','CD8.pre.exhausted','Monocytes','Bcells.proliferative','CD4.effector.memory',
#'                   'Plasma.Bcells','CD8.terminally.exhausted','cDC','CD4.recently.activated','NK',
#'                   'Tcells.proliferative','Macrophages.SPP1','CD4.naive.memory','Mast.cells',
#'                   'Macro.and.mono.prolif.','Tcells.naive','TAMs.C1QC','CD8.cytotoxic','pDC',
#'                   'TAMs.proinflamatory','Tcells.regulatory',
#'                   'T.helper.cells','mDC','CD8.effector.memory','CD4.transitional.memory')
#'   projected$functional.cluster <- plyr::mapvalues(projected$functional.cluster,from = old.idents,to = new.idents)
#'   saveRDS(projected, file = 'data/JerbyArnon.rds')
#' }else{
#'   projected <- readRDS('data/JerbyArnon.rds')
#' }
#' ProjecTME.FindDrugs(seu = projected, signature_collection = 'PSc', vars.to.regress = 'nFeature_RNA',
#' idents = 'functional.cluster', drugs = 'BORTEZOMIB', signatures = 'sig_18868', genes = 'PSMA5',
#' subcluster = 'Bcells', topgenes = 5)
ProjecTME.FindDrugs <- function(seu, signature_collection, include.pathways = T,
                                expr.thres = 0.1, pc = 10, k.neighbors = 20, res = 0.2,
                                method = 'uwot', UMAP = 'beyondcell', isRegressOut = TRUE,
                                vars.to.regress = 'nFeature_RNA', idents = 'functional.cluster',
                                drugs, signatures,genes,subcluster, topgenes = 5){
  
  if(!file.exists('beyondcell.obj.rds')){
    cat('Loading in-bult datasets, containing: "PSc", "SSc", "DSS", "drugInfo" and \"pathways\"')
    if(signature_collection=='PSc'){
      load('data/beyondcell/PSc.RData')
    }else if(signature_collection=='SSc'){
      load('data/beyondcell/SSc.RData')
    }else{
      stop('Currently only support PSc and SSc!!!')
    }
    cat(paste0('\nGenerate geneset object with: ',signature_collection,'.\n'))
    sig_collec <- get(load(paste0('data/beyondcell/',signature_collection, '.RData')))
    gs <- GenerateGenesets(sig_collec, include.pathways = include.pathways)
    
    DefaultAssay(seu) <- 'RNA'
    message(paste0('--- Compute the BCS (Beyondcell Score) ---\n'))
    bc <- bcScore(sc = seu, gs = gs, expr.thres = expr.thres)
    
    message(paste0('--- Compute the TCs (Therapeutic Clusters) ---\n'))
    bc <- bcUMAP(bc, pc = pc, k.neighbors = k.neighbors, res = res, method = method)
    if(isRegressOut==TRUE){
      print('-- Regress out unwanted sources of variation --')
      bc <- bcRegressOut(bc, vars.to.regress = vars.to.regress)
      cat('Recompute the UMAP...')
      bc <- bcUMAP(bc, pc = pc, k.neighbors = k.neighbors, res = res, method = method)
    }
    
    message(paste0('--- Compute ranks ---\n'))
    cat(paste0('Obtain general statistics\n'))
    bc <- bcRanks(bc)
    cat(paste0('Obtain ',idents,'-based statistics\n'))
    bc <- bcRanks(bc, idents = idents)
    cat(paste0('Obtain unextended therapeutic bc_clusters_res.',res,'-based statistics'))
    bc <- bcRanks(bc, idents = paste0('bc_clusters_res.',res), extended = F)
    saveRDS(bc, file = 'beyondcell.obj.rds')
  }else{
    bc <- readRDS('beyondcell.obj.rds')
  }
  
  message(paste0('--- Visualization ---\n'))
  p1 <- bcSignatures(bc, UMAP = UMAP, signatures = list(values = signatures),pt.size = 1.5)
  ggsave(filename = paste0('signature_',signatures,'.pdf'), p1[[1]], height = 6,width = 6.5)
  
  drug_ids <- FindDrugs(bc, drugs)$IDs
  drug <- bcSignatures(bc, UMAP = UMAP, signatures = list(values = drug_ids), pt.size = 1.5)
  pdf(file = paste0('signature_', drug_ids, '.pdf'), height = 5,width = 15)
  p2 <- cowplot::plot_grid(plotlist = drug, ncol = 3)
  print(p2)
  dev.off()
  p3 <- bcSignatures(bc, UMAP = UMAP, genes = list(values = genes))
  ggsave(filename = paste0('signature_',genes,'.pdf'), p3[[1]], height = 6,width = 6.5)
  p4 <- bc4Squares(bc, idents = idents, lvl = subcluster, top = topgenes)
  ggsave(filename = paste0('squares_',subcluster, '.pdf'), p4[[1]], height = 6,width = 9)
  p5 <- bcHistogram(bc, signatures = signatures, idents = NULL)
  ggsave(filename = paste0('histogram_',signatures,'.pdf'), p5[[1]], height = 6,width = 6.5)
  pdf(file = paste0('histogram_',signatures,'_celltype.pdf'), width = 5,height = 3 * length(unique(bc@meta.data[[idents]])))
  p6 <- bcHistogram(bc, signatures = signatures, idents = idents)
  print(p6)
  dev.off()
}
