#' @title readMTX
# ---- utils ----
readMTX <- function(mtx,cells,features,cell.column=1,feature.column=2,skip.cell=0,
                    skip.feature=0,unique.features = T, strip.suffix = F){
  all.files <- list('expression matrix'=mtx, 'barcode list'=cells, 'feature list'=features)
  for(i in seq_along(along.with = all.files)){
    all.files[[i]] <- normalizePath(all.files[[i]], mustWork = F)
  }
  cell.barcodes <- read.table(all.files[['barcode list']],header = F,sep = '\t',row.names = NULL,skip = skip.cell)
  feature.names <- read.table(all.files[['feature list']],header = F,sep = '\t',row.names = NULL,skip = skip.feature)
  bcols <- ncol(cell.barcodes)
  if(bcols < cell.column){
    stop('cell.column was set to ',cell.column,' but ',cells, ' only has ', bcols, ' columns.',
         'Try setting the cell.column argument to a value <= to ', bcols,'.')
  }
  cell.names <- cell.barcodes[, cell.column]
  if(all(grepl('\\-1$',cell.names)) & strip.suffix){
    cell.names <- stringr::str_split(cell.names,'-',simplify = T)[,1]
  }
  # read features
  fcols <- ncol(feature.names)
  if(fcols < feature.column){
    stop('feature.column was set to ',feature.column,' but ',features, ' only has ',fcols,' column(s).',
         'Try setting the feature.column argument to a value <= to',fcols,'.')
  }
  if(any(is.na(feature.names[,feature.column]))){
    na.features <- which(is.na(feature.names[,feature.column]))
    replace.column <- ifelse(feature.column==2, 1, 2)
    
    if(replace.column > fcols){
      stop('Some feature names are NA in column ',feature.column,'. Try sepcifying a different column.',call. = F)
    }else{
      warning('Some feature names are NA in column ',feature.column,'. Replacing NA names with ID from column ',replace.column,'.',call. = F)
    }
    feature.names[na.features,feature.column] <- feature.names[na.features,replace.column]
  }
  feature.names <- feature.names[,feature.column]
  if(unique.features) feature.names <- make.unique(feature.names)
  data <- Matrix::readMM(all.files[['expression matrix']])
  if(length(cell.names) != ncol(data)){
    stop('Matrix has ',ncol(data),' columns but found ',length(cell.names),' barcodes. ',
         ifelse(length(cell.names) > ncol(data), 'Trying increasing `skip.cell`. ', ''), call. = F)
  }
  if(length(feature.names) != nrow(data)){
    stop('Matrix has ',ncol(data), ' rows but found ',length(feature.names),' features. ',
         ifelse(length(feature.names)>nrow(data), 'Trying increasing `skip.feature`.', ''),call. = F)
  }
  colnames(data) <- cell.names
  rownames(data) <- feature.names
  data <- as(data, Class = 'dgCMatrix')
  return(data)
}

guess_raw_separator <- function(f,sep = c(' ','\t',',')){
  lines <- readLines(f, n = 10)
  if(length(lines) ==0) return(NULL)
  spl <- lapply(sep, grep, lines)
  counts <- unlist(lapply(spl, length))
  if(max(counts)==0) return(NULL)
  sep.idx <- which(counts == max(counts))[1]
  return(sep[sep.idx])
}

# Automatically determine species and gene ID column
get.species <- function(genes,table=Hs2Mm.convert.table){
  g.mm <- length(intersect(genes, table$Gene.MM))
  g.hs1 <- length(intersect(genes, table$Gene.stable.ID.HS))
  g.hs2 <- length(intersect(genes, table$Gene.HS))
  gg <- c(g.mm, g.hs1, g.hs2)
  if(max(gg) == g.mm){
    species <- 'mouse'
    col.id <- 'Gene.MM'
  }else{
    species <- 'human'
    col.id <- ifelse(g.hs1 > g.hs2, 'Gene.stable.ID.HS', 'Gene.HS')
  }
  res <- list('species' = species, 'col.id' = col.id)
  return(res)
}

# Internal function to filter cells
filterCells <- function(query.object, species="mouse", gating.model=NULL){
  load('data/cell.cycle.obj.RData')
  query.object <- suppressWarnings(
    scGate(data=query.object, model = gating.model, verbose=FALSE, assay=DefaultAssay(query.object),
           additional.signatures = cell.cycle.obj[[species]])
  )
  ncells <- ncol(query.object)
  ncells.keep <- sum(query.object$is.pure == 'Pure')
  if(ncells.keep > 0){
    query.object <- subset(query.object, subset=is.pure=='Pure') 
  }else{
    query.object <- NULL
  }
  message <- sprintf("%i out of %i ( %i%% ) non-pure cells removed. Use filter.cells=FALSE to avoid pre-filtering",
                     ncells - ncells.keep, ncells, round(100*(ncells-ncells.keep)/ncells))
  print(message)
  if(ncells.keep == 0){stop("Stopping. All cells were removed by cell filter!")}
  #Parse metadata columns
  query.object$cycling.score <- query.object$cycling_UCell
  query.object$cycling.score.G1_S <- query.object$cycling_G1.S_UCell
  query.object$cycling.score.G2_M <- query.object$cycling_G2.M_UCell
  to_remove <- grep("is.pure", colnames(query.object@meta.data))
  to_remove <- c(to_remove, grep("_UCell$", colnames(query.object@meta.data), perl=T))
  query.object@meta.data <- query.object@meta.data[,-to_remove]
  return(query.object)
}

# Internal function for mouse-human ortholog conversion
convert.orthologs <- function(obj,table,from="Gene.HS",to="Gene.MM",query.assay="RNA",slot="counts"){
  exp.mat <- slot(obj@assays[[query.assay]], name=slot)
  genes.select <- rownames(exp.mat) %in% table[[from]]
  if(length(genes.select) < 100){
    message("Warning: fewer than 100 genes with orthologs were found. Check your matrix format and gene names")
  }
  if(length(genes.select) > 0){
    exp.mat <- exp.mat[rownames(exp.mat) %in% table[[from]], ]
  }else{
    stop(paste0("Error: No genes found in column ", from))
  }
  #Convert
  ortho.genes <- table[[to]][match(row.names(exp.mat), table[[from]])]
  #Update matrix gene names
  row.names(exp.mat) <- ortho.genes
  slot(obj@assays[[query.assay]], name=slot) <- exp.mat
  if(slot=="data"){  #keep same size of data matrices
    slot(obj@assays[[query.assay]], name="counts") <- exp.mat
  }
  return(obj)
}

# Internal function to randomly split an object into subsets
randomSplit <- function(obj, n = 2, seed = 44, verbose = F){
  set.seed(seed)
  lgt <- dim(obj)[2]
  ind <- sample.int(n,lgt,replace = T)
  cell.list <- split(colnames(obj), ind)
  seurat.list <- list()
  if(verbose==T) message(sprintf('Splitting object into %i random subsets', n))
  for(h in 1:n) seurat.list[[h]] <- subset(obj, cells = cell.list[[h]])
  return(seurat.list)
}

# Find integration anchors using reciprocal PCA
FindIntegrationAnchors_local <- function(object.list = NULL,assay = NULL,anchor.coverage = 1,
                                         correction.scale = 100,alpha=0.5,anchor.features = 2000,
                                         sct.clip.range = NULL,l2.norm = TRUE,dims = 1:30, k.anchor = 5,
                                         k.filter = NA,k.score = 30,remove.thr = 0,max.features = 200,
                                         nn.method = "annoy",n.trees = 50,eps = 0,verbose = TRUE){
  normalization.method <- "LogNormalize"
  reference <- NULL
  reduction <- "pca"
  object.ncells <- sapply(object.list, function(x) dim(x = x)[2])
  if(any(object.ncells <= max(dims))){
    bad.obs <- which(object.ncells <= max(dims))
    stop("Max dimension too large: objects ", paste(bad.obs, collapse = ", "),
         " contain fewer than ", max(dims), " cells. \n Please specify a",
         " maximum dimensions that is less than the number of cells in any ",
         "object (", min(object.ncells), ").")
  }
  if(!is.null(assay)){
    if(length(assay) != length(object.list)){
      stop("If specifying the assay, please specify one assay per object in the object.list")
    }
    object.list <- sapply(1:length(object.list), function(x){
      DefaultAssay(object = object.list[[x]]) <- assay[x]
      return(object.list[[x]])
    })
  }else{
    assay <- sapply(X = object.list, FUN = DefaultAssay)
  }
  object.list <- CheckDuplicateCellNames_local(object.list = object.list)
  slot <- "data"
  nn.reduction <- reduction
  internal.neighbors <- list()
  if(verbose){message("Computing within dataset neighborhoods")}
  k.neighbor <- max(k.anchor, k.score)
  internal.neighbors <- lapply(1:length(object.list), function(x){
    Seurat:::NNHelper(
      data = Embeddings(object = object.list[[x]][[nn.reduction]])[, dims],
      k = k.neighbor + 1,method = nn.method,n.trees = n.trees,eps = eps
    )
  })
  # determine the proper offsets for indexing anchors
  objects.ncell <- sapply(object.list, ncol)
  offsets <- as.vector(cumsum(c(0, objects.ncell)))[1:length(object.list)]
  if (verbose) {message("Finding all pairwise anchors")}
  i <- 1; j <- 2
  object.1 <- DietSeurat(object = object.list[[i]], assays = assay[i], features = anchor.features,
                         counts = FALSE, scale.data = TRUE, dimreducs = reduction)
  object.2 <- DietSeurat(object = object.list[[j]],assays = assay[j],features = anchor.features,
                         counts = FALSE, scale.data = TRUE, dimreducs = reduction)
  # suppress key duplication warning
  suppressWarnings(object.1[["ToIntegrate"]] <- object.1[[assay[i]]])
  DefaultAssay(object.1) <- "ToIntegrate"
  if(reduction %in% Reductions(object = object.1)){
    slot(object = object.1[[reduction]], name = "assay.used") <- "ToIntegrate"
  }
  object.1 <- DietSeurat(object = object.1,assays = "ToIntegrate",counts = FALSE,scale.data = TRUE,dimreducs = reduction)
  suppressWarnings(object.2[["ToIntegrate"]] <- object.2[[assay[j]]])
  DefaultAssay(object.2) <- "ToIntegrate"
  if(reduction %in% Reductions(object = object.2)){
    slot(object = object.2[[reduction]], name = "assay.used") <- "ToIntegrate"
  }
  object.2 <- DietSeurat(object = object.2,assays = "ToIntegrate",counts = FALSE,scale.data = TRUE,dimreducs = reduction)
  #Reciprocal PCA
  common.features <- intersect(rownames(Loadings(object.1[["pca"]])),rownames(Loadings(object.2[["pca"]])))
  common.features <- intersect(common.features, anchor.features)
  object.pair <- merge(object.1, object.2, merge.data = TRUE)
  projected.embeddings.1<- t(GetAssayData(object = object.1, slot = "scale.data")[common.features, ]) %*%
    Loadings(object = object.2[["pca"]])[common.features, ]
  object.pair[['projectedpca.1']] <- CreateDimReducObject(
    embeddings = rbind(projected.embeddings.1, Embeddings(object = object.2[["pca"]])),
    assay = DefaultAssay(object = object.1),
    key = "projectedpca1_"
  )
  projected.embeddings.2 <- t(GetAssayData(object = object.2, slot = "scale.data")[common.features, ]) %*%
    Loadings(object = object.1[["pca"]])[common.features, ]
  object.pair[['projectedpca.2']] <- CreateDimReducObject(
    embeddings = rbind(projected.embeddings.2, Embeddings(object = object.1[["pca"]])),
    assay = DefaultAssay(object = object.2),key = "projectedpca2_")
  object.pair[["pca"]] <- CreateDimReducObject(
    embeddings = rbind(Embeddings(object = object.1[["pca"]]),Embeddings(object = object.2[["pca"]])),
    assay = DefaultAssay(object = object.1), key = "pca_")
  reduction <- "projectedpca.1"; reduction.2 <- "projectedpca.2"
  if (l2.norm){
    slot(object = object.pair[["projectedpca.1"]], name = "cell.embeddings") <- Sweep_local(
      x = Embeddings(object = object.pair[["projectedpca.1"]]), MARGIN = 2,
      STATS = apply(X = Embeddings(object = object.pair[["projectedpca.1"]]),MARGIN = 2,FUN = sd),
      FUN = "/"
    )
    slot(object = object.pair[["projectedpca.2"]], name = "cell.embeddings") <- Sweep_local(
      x = Embeddings(object = object.pair[["projectedpca.2"]]), MARGIN = 2,
      STATS = apply(X = Embeddings(object = object.pair[["projectedpca.2"]]), MARGIN = 2, FUN = sd),
      FUN = "/"
    )
    object.pair <- L2Dim(object = object.pair, reduction = "projectedpca.1")
    object.pair <- L2Dim(object = object.pair, reduction = "projectedpca.2")
    reduction <- paste0(reduction, ".l2")
    reduction.2 <- paste0(reduction.2, ".l2")
  }
  internal.neighbors <- internal.neighbors[c(i, j)]
  anchors <- FindAnchors_local(
    object.pair = object.pair,assay = c("ToIntegrate", "ToIntegrate"),slot = slot,
    cells1 = colnames(x = object.1),cells2 = colnames(x = object.2),
    internal.neighbors = internal.neighbors,reduction = reduction,reduction.2 = reduction.2,
    nn.reduction = nn.reduction,dims = dims,k.anchor = k.anchor,k.filter = k.filter,
    k.score = k.score,max.features = max.features,nn.method = nn.method,
    n.trees = n.trees,eps = eps,verbose = verbose)
  anchors[, 1] <- anchors[, 1] + offsets[i]
  anchors[, 2] <- anchors[, 2] + offsets[j]
  #Average distances
  anchors <- as.data.frame(anchors)
  anchors$dist.mean <- apply(anchors[,c("dist1.2","dist2.1")], MARGIN=1, mean)
  message(sprintf("  SD on anchor distances: %.3f",sd(anchors$dist.mean)))
  if(anchor.coverage < 1){
    #Combine anchor distance with anchor score
    sigmoid_center <- unname(quantile(anchors$dist.mean, probs = anchor.coverage, na.rm = T))
    distance_factors <-  sigmoid(x = anchors$dist.mean, center = sigmoid_center, scale = correction.scale)
    #Multiply distance factors by score
    anchors$score <- anchors$score * distance_factors
    ##Remove distant anchors
    anchors <- anchors[distance_factors > remove.thr,] 
  }
  nanchors <- nrow(anchors)
  ##Include reciprocal anchors
  anchors <- rbind(anchors[, c("cell1","cell2","score","dist.mean")],anchors[, c("cell2","cell1","score","dist.mean")])  
  anchors <- AddDatasetID_local(anchor.df = anchors, offsets = offsets, obj.lengths = objects.ncell)
  command <- LogSeuratCommand(object = object.list[[1]], return.command = TRUE)
  anchor.set <- new(Class = "IntegrationAnchorSet", object.list = object.list,
                    reference.objects = seq_along(object.list),
                    anchors = anchors, offsets = offsets,
                    anchor.features = anchor.features,command = command)
  return(anchor.set)
}

# Ensure no duplicate cell names
CheckDuplicateCellNames_local <- function(object.list,verbose=T,stop=F){
  cell.names <- unlist(lapply(object.list,colnames))
  if(any(duplicated(cell.names))){
    if(stop) stop('Duplicate cell names present across objects provided.')
    if(verbose) warning('Some cell names are duplicated across object provided. Renaming to enforce unique cell names.')
    object.list <- lapply(1:length(object.list),function(x){
      return(RenameCells(object.list[[x]], new.names = paste0(Cells(object.list[[x]]),'_',x)))
    })
  }
  return(object.list)
}

Sweep_local <- function(x,MARGIN,STATS,FUN='-',check.margin=T,...){
  if(any(grepl('X',names(formals(fun = sweep))))){
    return(sweep(x = x, MARGIN = MARGIN, STATS = STATS, FUN = FUN, check.margin = check.margin, ...))
  }else{
    return(sweep(x = x, MARGIN = MARGIN, STATS = STATS, FUN = FUN, check.margin = check.margin, ...))
  }
}

# Find anchors between a pair of objects
FindAnchors_local <- function(object.pair,assay,slot,cells1,cells2,internal.neighbors,reduction,
                              reduction.2 = character(),nn.reduction = reduction,dims = 1:10,k.anchor = 5,
                              k.filter = NA,k.score = 30,max.features = 200,nn.method = "annoy",n.trees = 50,
                              nn.idx1 = NULL, nn.idx2 = NULL, eps = 0, verbose = TRUE){
  # compute local neighborhoods, use max of k.anchor and k.score if also scoring to avoid recomputing neighborhoods
  k.neighbor <- k.anchor
  if (!is.na(k.score)) { k.neighbor <- max(k.anchor, k.score) }
  object.pair <- FindNN_local(object = object.pair,cells1 = cells1,cells2 = cells2,internal.neighbors = internal.neighbors,
                              dims = dims,reduction = reduction,reduction.2 = reduction.2,nn.reduction = nn.reduction,
                              k = k.neighbor,nn.method = nn.method,n.trees = n.trees,nn.idx1 = nn.idx1,nn.idx2 = nn.idx2,
                              eps = eps, verbose = verbose)
  object.pair <- FindAnchorPairs_local(object = object.pair,integration.name = "integrated",k.anchor = k.anchor,verbose = verbose)
  if (!is.na(k.score)){
    object.pair = ScoreAnchors_local(object = object.pair,assay = DefaultAssay(object = object.pair),
                                     integration.name = "integrated",verbose = verbose,k.score = k.score)
  }
  ###Return distances
  anc.tab <- object.pair@tools$integrated@anchors
  d1.2 <- numeric(length = dim(anc.tab)[1])
  d2.1 <- numeric(length = dim(anc.tab)[1])
  for (r in 1:dim(anc.tab)[1]) {
    c1 <- anc.tab[r,"cell1"]
    c2 <- anc.tab[r,"cell2"]
    d1.2[r] <- object.pair@tools$integrated@neighbors$nnab@nn.dist[c1, which(object.pair@tools$integrated@neighbors$nnab@nn.idx[c1,] == c2)]
    d2.1[r] <- object.pair@tools$integrated@neighbors$nnba@nn.dist[c2, which(object.pair@tools$integrated@neighbors$nnba@nn.idx[c2,] == c1)]
  }
  object.pair@tools$integrated@anchors <- cbind(object.pair@tools$integrated@anchors, dist1.2=d1.2)
  object.pair@tools$integrated@anchors <- cbind(object.pair@tools$integrated@anchors, dist2.1=d2.1)
  anchors <- GetIntegrationData(object = object.pair,integration.name = 'integrated',slot = 'anchors')
  return(anchors)
}

#Find anchor pairs
FindAnchorPairs_local <- function(object,integration.name = 'integrated',k.anchor = 5,verbose = TRUE){
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = 'neighbors')
  max.nn <- c(ncol( neighbors$nnab), ncol(neighbors$nnba))
  if(any(k.anchor > max.nn)){
    message(paste0('warning: requested k.anchor = ', k.anchor, ', only ', min(max.nn), ' in dataset'))
    k.anchor <- min(max.nn)
  }
  if(verbose){message("Finding anchors")}
  # convert cell name to neighbor index
  nn.cells1 <- neighbors$cells1
  nn.cells2 <- neighbors$cells2
  cell1.index <-  suppressWarnings(which(colnames(object) == nn.cells1, arr.ind = TRUE))
  ncell <- 1:nrow(neighbors$nnab)
  ncell <- ncell[ncell %in% cell1.index]
  anchors <- list()
  # pre allocate vector
  anchors$cell1 <- rep(0, length(ncell) * 5)
  anchors$cell2 <- anchors$cell1
  anchors$score <- anchors$cell1 + 1
  idx <- 0
  indices.ab <- Indices(object = neighbors$nnab)
  indices.ba <- Indices(object = neighbors$nnba)
  for(cell in ncell){
    neighbors.ab <- indices.ab[cell, 1:k.anchor]
    mutual.neighbors <- which(indices.ba[neighbors.ab, 1:k.anchor, drop = FALSE] == cell,arr.ind = TRUE)[, 1]
    for (i in neighbors.ab[mutual.neighbors]){
      idx <- idx + 1
      anchors$cell1[idx] <- cell
      anchors$cell2[idx] <- i
      anchors$score[idx] <- 1
    }
  }
  anchors$cell1 <- anchors$cell1[1:idx]
  anchors$cell2 <- anchors$cell2[1:idx]
  anchors$score <- anchors$score[1:idx]
  anchors <- t(do.call(rbind, anchors))
  anchors <- as.matrix(anchors)
  object <- SetIntegrationData(object=object,integration.name=integration.name,slot='anchors',new.data=anchors)
  if(verbose){message(paste0("\tFound ", nrow(anchors), " anchors"))}
  return(object)
}

# Find nearest neighbors
FindNN_local <- function(object,cells1 = NULL,cells2 = NULL,internal.neighbors,grouping.var = NULL,dims = 1:10,
                         reduction = "cca.l2",reduction.2 = character(),nn.dims = dims,nn.reduction = reduction,
                         k = 300,nn.method = "annoy",n.trees = 50,nn.idx1 = NULL,nn.idx2 = NULL,eps = 0,
                         integration.name = 'integrated',verbose = TRUE){
  if(xor(is.null(cells1),is.null(cells2))) stop('cells1 and cells2 must both be specified')
  if(!is.null(cells1) && !is.null(cells2) && !is.null(grouping.var)){
    stop("Specify EITHER grouping.var or cells1/2.")
  }
  if(is.null(cells1) && is.null(cells2) && is.null(grouping.var)){
    stop("Please set either cells1/2 or grouping.var")
  }
  if(!is.null(grouping.var)){
    if(nrow(unique(object[[grouping.var]])) != 2){
      stop("Number of groups in grouping.var not equal to 2.")
    }
    groups <- names(sort(table(object[[grouping.var]]), decreasing = TRUE))
    cells1 <- colnames(object)[object[[grouping.var]] == groups[[1]]]
    cells2 <- colnames(object)[object[[grouping.var]] == groups[[2]]]
  }
  if(verbose) message("Finding neighborhoods")
  dim.data.self <- Embeddings(object = object[[nn.reduction]])[, nn.dims]
  if(!is.null(internal.neighbors[[1]])){
    nnaa <- internal.neighbors[[1]]
  }else{
    dims.cells1.self <- dim.data.self[cells1, ]
    nnaa <- Seurat:::NNHelper(data=dims.cells1.self,k=k+1,method=nn.method,n.trees=n.trees,eps=eps,index=nn.idx1)
  }
  if(!is.null(internal.neighbors[[2]])){
    nnbb <- internal.neighbors[[2]]
  }else{
    dims.cells2.self <- dim.data.self[cells2, ]
    nnbb <- Seurat:::NNHelper(data=dims.cells2.self,k=k+1,method=nn.method,n.trees=n.trees,eps=eps,index=nn.idx1)
  }
  if(length(reduction.2) > 0){
    nnab <- Seurat:::NNHelper(
      data = Embeddings(object = object[[reduction.2]])[cells2, ],
      query = Embeddings(object = object[[reduction.2]])[cells1, ],
      k = k,method = nn.method,n.trees = n.trees,eps = eps,index = nn.idx2)
    nnba <- Seurat:::NNHelper(
      data = Embeddings(object = object[[reduction]])[cells1, ],
      query = Embeddings(object = object[[reduction]])[cells2, ],
      k = k,method = nn.method,n.trees = n.trees,eps = eps,index = nn.idx1)
  }else{
    dim.data.opposite <- Embeddings(object = object[[reduction]])[ ,dims]
    dims.cells1.opposite <- dim.data.opposite[cells1, ]
    dims.cells2.opposite <- dim.data.opposite[cells2, ]
    nnab <- Seurat:::NNHelper(
      data = dims.cells2.opposite,query = dims.cells1.opposite,k = k,
      method = nn.method,n.trees = n.trees,eps = eps,index = nn.idx2)
    nnba <- Seurat:::NNHelper(
      data = dims.cells1.opposite,query = dims.cells2.opposite,k = k,
      method = nn.method,n.trees = n.trees,eps = eps,index = nn.idx1)
  }
  object <- SetIntegrationData(
    object = object,integration.name = integration.name,slot = 'neighbors',
    new.data = list('nnaa'=nnaa, 'nnab'=nnab, 'nnba'=nnba, 'nnbb'=nnbb, 'cells1'=cells1, 'cells2'=cells2)
  )
  return(object)
}

# Find anchor pairs
ScoreAnchors_local <- function(object,assay=NULL,integration.name='integrated',verbose=TRUE,k.score=30){
  if (is.null(assay)) assay <- DefaultAssay(object)
  anchor.df <- as.data.frame(GetIntegrationData(object=object,integration.name=integration.name,slot='anchors'))
  neighbors <- GetIntegrationData(object=object,integration.name=integration.name,slot="neighbors")
  offset <- length(neighbors$cells1)
  indices.aa <- Indices(object = neighbors$nnaa)
  indices.bb <- Indices(object = neighbors$nnbb)
  indices.ab <- Indices(object = neighbors$nnab)
  indices.ba <- Indices(object = neighbors$nnba)
  nbrsetA <- function(x) c(indices.aa[x, 1:k.score], indices.ab[x, 1:k.score] + offset)
  nbrsetB <- function(x) c(indices.ba[x, 1:k.score], indices.bb[x, 1:k.score] + offset)
  # score = number of shared neighbors
  anchor.new <- data.frame(
    'cell1' = anchor.df[, 1], 'cell2' = anchor.df[, 2],
    'score' = mapply(function(x,y){length(intersect(nbrsetA(x),nbrsetB(y)))},anchor.df[,1],anchor.df[,2])
  )
  # normalize the score
  max.score <- quantile(anchor.new$score, 0.9)
  min.score <- quantile(anchor.new$score, 0.01)
  anchor.new$score <- anchor.new$score - min.score
  anchor.new$score <- anchor.new$score / (max.score - min.score)
  anchor.new$score[anchor.new$score > 1] <-  1
  anchor.new$score[anchor.new$score < 0] <- 0
  anchor.new <- as.matrix(anchor.new)
  object <- SetIntegrationData(object=object,integration.name=integration.name,slot='anchors',new.data=anchor.new)
  return(object)
}

sigmoid <- function(x,scale,center){
  sigm <- 1/(1 + exp(scale*(x-center)))
  return(sigm)
}

# Add dataset ID
AddDatasetID_local <- function(anchor.df,offsets,obj.lengths){
  ndataset <- length(offsets)
  row.offset <- rep.int(x = offsets, times = obj.lengths)
  dataset <- rep.int(x = 1:ndataset, times = obj.lengths)
  anchor.df <- data.frame(
    'cell1' = anchor.df[,'cell1'] - row.offset[anchor.df[,'cell1']],
    'cell2' = anchor.df[,'cell2'] - row.offset[anchor.df[,'cell2']],
    'score' = anchor.df[,'score'], 'dataset1' = dataset[anchor.df[,'cell1']],
    'dataset2' = dataset[anchor.df[,'cell2']], 'dist.mean' = anchor.df[,'dist.mean']
  )
  return(anchor.df)
}

# Calculate top features across a set of dimensions
TopDimFeatures_local <- function(object,reduction,dims=1:10,features.per.dim=100,max.features=200,projected=F){
  dim.reduction <- object[[reduction]]
  max.features <- max(length(dims)*2, max.features)
  num.features <- sapply(1:features.per.dim,function(y){
    length(unique(as.vector(sapply(dims, function(x){
      unlist(TopFeatures(object = dim.reduction, dim = x, nfeatures = y, balanced = T, projected = projected))
    }))))
  })
  max.per.pc <- which.max(num.features[num.features < max.features])
  features <- unique(as.vector(sapply(dims, function(x){
    unlist(TopFeatures(object = dim.reduction, dim = x, nfeatures = max.per.pc, balanced = T, projected = projected))
  })))
  features <- unique(features)
  return(features)
}

# Helper for projecting individual data sets
projection.helper <- function(query, ref=NULL, filter.cells=TRUE, query.assay=NULL, 
                              direct.projection=FALSE, fast.umap.predict=FALSE, ortholog_table=NULL,
                              k.weight=100, k.anchor=5, skip.normalize=FALSE, id="query1",
                              anchor.coverage=1, correction.scale=100, alpha=0.5, remove.thr=0,
                              scGate_model=NULL, ncores=ncores){
  retry.direct <- FALSE
  do.orthology <- FALSE
  # Reference
  DefaultAssay(ref) <- "integrated"
  ref.var.features <- ref@assays$integrated@var.features
  # If query.assay not specified, use the default
  if(is.null(query.assay)){
    query.assay <- DefaultAssay(query)
  }else{
    DefaultAssay(query) <- query.assay
  }
  print(paste0("Using assay ",query.assay," for ",id))
  if(!is.null(ref@misc$umap_object$data)){ 
    pca.dim=dim(ref@misc$umap_object$data)[2] #use the number of PCs used to build the reference
  }else{ pca.dim=10 }
  species.ref <- get.species(genes = rownames(ref), table = ortholog_table)
  species.query <- get.species(genes = rownames(query), table = ortholog_table)
  if(species.ref$species != species.query$species){ do.orthology <- TRUE }
  if(filter.cells){
    message("Pre-filtering cells with scGate...")
    if(is.null(scGate_model)){  #read filter model from atlas
      if(!is.null(ref@misc$scGate[[species.query$species]])){
        scGate_model <- ref@misc$scGate[[species.query$species]]
      }else{   #if no model was specified, and no model was found in the atlas, use a default filter
        message("No scGate model specified: using default filter for T cells")
        models <- readRDS('data/models.rds')
        scGate_model <- models[[species.query$species]]$generic$Tcell  
      }
    }
    query <- filterCells(query, species=species.query$species, gating.model=scGate_model)
  }
  if(is.null(query)){
    message(sprintf("Warning! Skipping %s - all cells were removed by cell filter", id))   #Update text
    return(NULL)
  }
  
  # Check if slots are populated, and normalize data.
  if(skip.normalize){
    slot <- "data"
    exp.mat <-  slot(query@assays[[query.assay]], name=slot)
    if(dim(exp.mat)[1]==0){
      stop("Data slot not found in your Seurat object. Please normalize the data")
    }
    if(do.orthology){
      print("Transforming expression matrix into space of orthologs") 
      query <- convert.orthologs(query,table=ortholog_table,query.assay=query.assay,slot=slot,from=species.query$col.id,to=species.ref$col.id)
    }        
  }else{
    slot <- "counts"
    exp.mat <-  slot(query@assays[[query.assay]], name=slot)
    if(dim(exp.mat)[1]==0){
      stop("Counts slot not found in your Seurat object. If you already normalized your data, re-run with option skip.normalize=TRUE")
    }
    if(do.orthology){
      print("Transforming expression matrix into space of orthologs") 
      query <- convert.orthologs(query,table=ortholog_table,query.assay=query.assay,slot=slot,from=species.query$col.id,to=species.ref$col.id)
    }
    query@assays[[query.assay]]@data <- query@assays[[query.assay]]@counts
    query <- NormalizeData(query) 
  }
  rm(exp.mat)
  query <- RenameCells(query, add.cell.id = "Q")
  query.metadata <- query@meta.data   #back-up metadata (and re-add it after projection)
  genes4integration <- intersect(ref.var.features, row.names(query))
  if(length(genes4integration)/length(ref.var.features)<0.5){ stop("Too many genes missing. Check input object format") }
  if (length(genes4integration)/length(ref.var.features)<0.8) {
    print("Warning! more than 20% of variable genes not found in the query")
  }
  if(direct.projection){
    projected <- query
    print("DIRECTLY projecting query onto Reference PCA space")
    query.pca.proj <-apply.pca.obj.2(query, pca.obj=ref@misc$pca_object,query.assay=query.assay)
    projected[["pca"]] <- CreateDimReducObject(embeddings = query.pca.proj,key = "PC_", assay = query.assay)
    print("DIRECTLY projecting query onto Reference UMAP space")
    query.umap.proj <- make.umap.predict(ref.umap=ref@misc$umap_object,pca.query.emb=query.pca.proj,fast.umap.predict=fast.umap.predict)
    projected[["umap"]] <- CreateDimReducObject(embeddings = query.umap.proj,key = "UMAP_", assay = query.assay)
    DefaultAssay(projected) <- query.assay
  }else{
    print(paste0("Aligning ", id, " to reference map for batch-correction..."))
    if(dim(ref@assays$integrated@scale.data)[2]==0){
      ref <- ScaleData(ref, do.center=FALSE, do.scale=FALSE, features = genes4integration)
    }
    query <- ScaleData(query, do.center=FALSE, do.scale=FALSE, features = genes4integration)
    ref <- RunPCA(ref, features = genes4integration,verbose = F)
    query <- RunPCA(query, features = genes4integration,verbose = F)
    proj.anchors <- FindIntegrationAnchors_local(
      object.list = c(ref, query),anchor.features = genes4integration,dims = 1:pca.dim,
      assay=c("integrated",query.assay),k.anchor=k.anchor,remove.thr=remove.thr,
      anchor.coverage=anchor.coverage,correction.scale=correction.scale,alpha=alpha
    )
    #Use all anchors for re-weighting - essentially disables local weighting
    if(k.weight == "max"){ k.weight <- length(unique(proj.anchors@anchors$cell2)) }
    #Do integration
    all.genes <- intersect(row.names(ref), row.names(query))
    proj.integrated <- IntegrateData(anchorset = proj.anchors, dims = 1:pca.dim,features.to.integrate = all.genes,
                                     k.weight = k.weight, preserve.order = T, verbose=F)
    #Subset query data from integrated space
    cells_query <- colnames(query)
    projected <- subset(proj.integrated, cells = cells_query)
    projected@meta.data <- query.metadata
    #Add anchor score to metadata
    aa <- proj.anchors@anchors
    aa.score <- aggregate(data=aa, x=aa$score, by=list(qcell = aa$cell2), FUN = mean)
    projected@meta.data[,"anchor.score"] <- 0
    projected@meta.data[aa.score$qcell, "anchor.score"] <- aa.score$x
    #Make PCA and UMAP projections
    cat("\nProjecting corrected query onto Reference PCA space\n")
    query.pca.proj <- apply.pca.obj.2(projected,pca.obj=ref@misc$pca_object,query.assay="integrated")
    projected[["pca"]] <- CreateDimReducObject(embeddings = query.pca.proj, key = "PC_", assay = "integrated")
    cat("\nProjecting corrected query onto Reference UMAP space\n")
    query.umap.proj <- make.umap.predict(ref.umap=ref@misc$umap_object,pca.query.emb=query.pca.proj,fast.umap.predict=fast.umap.predict)
    projected[["umap"]] <- CreateDimReducObject(embeddings = query.umap.proj, key = "UMAP_", assay = "integrated")
    DefaultAssay(projected) <- "integrated"
    if(retry.direct){
      projected <- query
      print("DIRECTLY projecting query onto Reference PCA space")
      query.pca.proj <- apply.pca.obj.2(query, pca.obj=ref@misc$pca_object,query.assay=query.assay)
      projected[["pca"]] <- CreateDimReducObject(embeddings = query.pca.proj,key = "PC_", assay = query.assay)
      print("DIRECTLY projecting query onto Reference UMAP space")
      query.umap.proj <- make.umap.predict(ref.umap=ref@misc$umap_object,pca.query.emb = query.pca.proj,fast.umap.predict=fast.umap.predict)
      projected[["umap"]] <- CreateDimReducObject(embeddings = query.umap.proj,key = "UMAP_", assay = query.assay)
      DefaultAssay(projected) <- query.assay
    }
  }
  if (!is.null(projected)) {
    projected@assays[[query.assay]]@var.features <- ref.var.features
    cellnames <- gsub("^Q_","",colnames(projected))  # remove prefix from cell names
    projected <- RenameCells(projected, new.names=cellnames)
  }
  return(projected)
}


apply.pca.obj.2 <- function(query, query.assay='RNA', pca.obj){
  newdata <- data.frame(t(as.matrix(query@assays[[query.assay]]@data)))
  newdata <- newdata[,order(names(newdata))]
  genes.use <- sort(intersect(colnames(newdata), names(pca.obj$center)))
  newdata.var <- newdata[,genes.use]
  center.use <- pca.obj$center[genes.use]
  scale.use <- pca.obj$scale[genes.use]
  rotation.use <- pca.obj$rotation[genes.use,]
  npca <- scale(newdata.var, center.use, scale.use) %*% rotation.use
  return(npca)
}

apply.ica.obj <- function(query, query.assay='RNA',ica.obj){
  newdata <- data.frame(t(as.matrix(query@assays[[query.assay]]@data)))
  newdata <- newdata[,order(names(newdata))]
  genes.use <- sort(intersect(colnames(newdata),names(ica.obj$center)))
  newdata.var <- newdata[,genes.use]
  center.use <- ica.obj$center[genes.use]
  scale.use <- ica.obj$scale[genes.use]
  npca <- scale(newdata.var,center.use,scale.use) %*% ica.obj$K[genes.use,] %*% ica.obj$W
  colnames(npca) <- colnames(ica.obj$S)
  return(npca)
}

# dispatch to UMAP prediction method (complete of fast)
make.umap.predict <- function(ref.umap, fast.umap.predict=F, ...){
  if(fast.umap.predict){
    nproj <- make.umap.predict.weighted.mean(ref.umap = ref.umap, ...)
  }else if(class(ref.umap)=='umap'){
    nproj <- make.umap.predict.2(ref.umap, method = 'umap', ...)
  }else if(!is.null(ref.umap$embedding)){
    nproj <- make.umap.predict.2(ref.umap = ref.umap, method = 'uwot', ...)
  }else{
    warning('No UMAP-predict model available. Using fast.umap.predict approximation.')
    nproj <- make.umap.predict.weighted.mean(ref.umap = ref.umap, ...)
  }
  return(nproj)
}

# UMAP predict using the umap package
make.umap.predict.2 <- function(ref.umap, query, query.assay='RNA', pca.obj,
                                pca.query.emb = NULL, method = 'uwot'){
  # if PCA query cell embeddings have been pre-calculated, read them from variable
  if(is.null(pca.query.emb)){
    pca.query.emb <- apply.pca.obj.2(query = query, query.assay = query.assay, pca.obj = pca.obj)
  }
  pca.dim <- dim(ref.umap$data)[2]
  if(method == 'umap'){
    nproj.umap <- umap:::predict.umap(ref.umap,pca.query.emb[,1:pca.dim])
  }else if(method == 'uwot'){
    nproj.umap <- uwot::umap_transform(pca.query.emb[,1:pca.dim], model = ref.umap)
  }else{
    stop('Unsupported UMAP method.')
  }
  return(nproj.umap)
}

# Fast projection mode: assign UMAP coordinates based on nearest neighbors in PCA space
make.umap.predict.weighted.mean <- function(ref.umap,query,query.assay='RNA',pca.obj,
                                            pca.query.emb=NULL, k = 8){
  if(is.null(pca.query.emb)){
    pca.query.emb <- apply.pca.obj.2(query = query, query.assay = query.assay, pca.obj = pca.obj)
  }
  ref.space <- ref.umap$data
  pca.dim <- ncol(ref.space)
  query.space <- pca.query.emb[,1:pca.dim]
  nn.ranked <- Seurat:::NNHelper(data = ref.space, query = query.space, k = k, method = 'rann')
  cellnames <- rownames(query.space)
  nproj.umap <- matrix(data = NA, nrow = length(cellnames), ncol = 2, dimnames = list(cellnames, c('UMAP_1','UMAP_2')))
  for(cell in 1:length(cellnames)){
    # calculate exp(-dist) as weights for nearest neighbors
    row <- exp(-nn.ranked@nn.dist[cell,])
    weights <- row/sum(row)
    # assign UMAP coordinates of (weighted) neighbors
    nproj.umap[cell,] <- weights %*% ref.umap$layout[nn.ranked@nn.idx[cell,],]
  }
  return(nproj.umap)
}

run.umap.2 <- function(pca.obj,ndim=NULL,n.neighbors=15,n.components=2,min.dist=0.3,metric='cosine',seed=1234){
  umap.config <- umap::umap.defaults
  umap.config$n_neighbors <- n.neighbors
  umap.config$min_dist <- min.dist
  umap.config$metric <- metric
  umap.config$n_components <- n.components
  umap.config$random_state <- seed
  umap.config$transform_state <- seed
  if(is.null(ndim)) ndim <- ncol(pca.obj$x)
  ref.umap <- umap::umap(pca.obj$x[,1:ndim], config = umap.config)
  colnames(ref.umap$layout) <- c('UMAP_1', 'UMAP_2')
  return(ref.umap)
}

run.ica <- function(object, assay='integrated', ndim=50){
  if(!require(fastICA)){install.packages('fastICA')}
  set.seed(1234)
  x <- scale(Matrix::t(object@assays[[assay]][object@assays[[assay]]@var.features,]))
  set.seed(1234)
  ref.ica <- fastICA::fastICA(x, n.comp=ndim, row.norm=T, maxit=1000, verbose=F, tol=1e-13, method='R')
  colnames(ref.ica$X) <- ref@assays$integrated@var.features
  rownames(ref.ica$X) <- colnames(ref)
  rownames(ref.ica$K) <- ref@assays$integrated@var.features
  colnames(ref.ica$A) <- colnames(ref.ica$X)
  colnames(ref.ica$S) <- paste0('ICA_', seq_len(ncol(ref.ica$S)))
  ref.ica$center <- attr(x, 'scaled:center')
  ref.ica$scale <- attr(x, 'scaled:scale')
  object[['ica']] <- CreateDimReducObject(embeddings=ref.ica$S,loadings=t(ref.ica$A),key='ICA_',assay=assay)
  object@misc$ica <- ref.ica
  return(object)
}

run.umap.uwot <- function(pca.obj,ndim=NULL,n.neighbors=15,n.components=2,min.dist=0.3,metric='cosine',seed=1234){
  if(is.null(ndim)) ndim <- ncol(pca.obj$x)
  set.seed(seed)
  ref.umap <- uwot::umap(pca.obj$x[,1:ndim],metric = metric, min_dist = min.dist,n_neighbors = n.neighbors, ret_model = T)
  colnames(ref.umap$embedding) <- c('UMAP_1', 'UMAP_2')
  return(ref.umap)
}

prcomp.seurat <- function(obj,assay=NULL,ndim=10,scale=T){
  if(is.null(assay)) assay <- DefaultAssay(obj)
  varfeat <- VariableFeatures(obj)
  refdata <- data.frame(t(as.matrix(obj@assays[[assay]]@data[varfeat,])))
  refdata <- refdata[,sort(colnames(refdata))]
  ref.pca <- prcomp(refdata, rank. = ndim, scale. = scale, center = T, retx = T)
  # save PCA rotation object
  obj@misc$pca_object <- ref.pca
  obj@reductions$pca@cell.embeddings <- ref.pca$x
  obj@reductions$pca@feature.loadings <- ref.pca$rotation
  colnames(obj@reductions$pca@cell.embeddings) <- gsub('PC(\\d+)','PC_\\1',colnames(ref.pca$x),perl = T)
  colnames(obj@reductions$pca@feature.loadings) <- gsub('PC(\\d+)','PC_\\1',colnames(ref.pca$rotation),perl = T)
  return(obj)
}

# calculate Silhouette coefficient between for cells in rows compared to set in columns with same labels
silhouette_2sets <- function(dist, labs.x, labs.y){
  labs.x <- as.character(labs.x); labs.y <- as.character(labs.y)
  ids <- sort(unique(c(labs.x, labs.y)))
  k <- length(ids)
  if(k <= 1) return(NA)
  if(nrow(dist) != length(labs.x)){
    stop(sprintf('Distance matrix has %i rows but %i row cluster labels are given',nrow(dist),length(labs.x)))
  }
  if(ncol(dist) != length(labs.y)){
    stop(sprintf('Distance matrix has %i columns but %i column cluster labels are given',ncol(dist),length(labs.y)))
  }
  res <- data.frame(matrix(NA, nrow(dist),2,dimnames = list(rownames(dist),c('cluster','sil_width'))))
  for(j in 1:k){
    lab <- ids[j]
    ix <- labs.x == lab; iy <- labs.y == lab
    Nx <- sum(ix); Ny <- sum(iy); Ny.n <- sum(!iy)
    if(Nx > 1){
      a.i <- rowSums(dist[ix, iy])/Ny
      b.i <- rowSums(dist[ix, !iy])/Ny.n
      s.i <- (b.i - a.i)/pmax(b.i, a.i)
      res[ix, 'cluster'] <- lab
      res[ix, 'sil_width'] <- s.i
    }
  }
  return(res)
}

circMarkCelltype <- function(seu, reduction, annot, sub, palette){
  library(Seurat)
  library(ggplot2)
  if(!require(ggunchull)) devtools::install_github('sajuukLyu/ggunchull',type = 'source')
  if(!require(tidydr)) install.packages('tidydr')
  dat <- as.data.frame(seu@reductions[[reduction]]@cell.embeddings)
  dat$celltype <- seu@meta.data[[annot]]
  ggplot(dat, aes(x = eval(str2expression(colnames(dat)[1])), y = eval(str2expression(colnames(dat)[2])), 
                  fill = celltype, color = celltype)) +
    stat_unchull(data = subset(dat,celltype==sub), alpha = 0.2, size = 1, 
                 show.legend = F, nbin = 300, nsm = 20,qval = 0.8,sfac = 1.5) +
    geom_point(size = 1) + theme_dr() + labs(x = colnames(dat)[1], y = colnames(dat)[2]) +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_manual(values = palette) + scale_color_manual(values = palette)
}