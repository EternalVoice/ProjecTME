#' @title Filter single-cell data by cell type
#' @description Apply scGate to filter specific cell types in a query dataset
#' @param data Seurat object containing a query data set - filtering will be applied to this object
#' @param model A single scGate model, or a list of scGate models. See Details for this format
#' @param pca.dim Number of dimensions for cluster analysis
#' @param assay Seurat assay to use
#' @param slot Data slot in Seurat object
#' @param pos.thr Minimum UCell score value for positive signatures
#' @param neg.thr Maximum UCell score value for negative signatures
#' @param maxRank Maximum number of genes that UCell will rank per cell
#' @param nfeatures Number of variable genes for dimensionality reduction
#' @param k.param Number of nearest neighbors for knn smoothing
#' @param reduction Dimensionality reduction to use for knn smoothing. By default, calculates a new reduction
#' based on the given \code{assay}; otherwise you may specify a precalculated dimensionality reduction (e.g.
#' in the case of an integrated dataset after batch-effect correction)
#' @param pca.dim Number of principal components for dimensionality reduction
#' @param resol Resolution for cluster analysis (if \code{by.knn=FALSE})
#' @param param_decay Controls decrease in parameter complexity at each iteration, between 0 and 1.
#' \code{param_decay == 0} gives no decay, increasingly higher \code{param_decay} gives increasingly stronger decay
#' @param ncores Number of processors for parallel processing
#' @param output.col.name Column name with 'pure/impure' annotation
#' @param min.cells Minimum number of cells to cluster or define cell types
#' @param additional.signatures A list of additional signatures, not included in the model, to be evaluated 
#' (e.g. a cycling signature). The scores for this list of signatures will be returned but not used for filtering.
#' @param save.levels Whether to save in metadata the filtering output for each gating model level
#' @param by.knn Perform k-nearest neighbor smoothing for UCell scores
#' @param keep.ranks Store UCell rankings in Seurat object. This will speed up calculations if the same object is 
#' applied again with new signatures.
#' @param genes.blacklist Genes blacklisted from variable features. The default loads the list of genes 
#' in \code{scGate::genes.blacklist.default}; you may deactivate blacklisting by setting \code{genes.blacklist=NULL}
#' @param multi.asNA How to label cells that are "Pure" for multiple annotations: "Multi" (FALSE) or NA (TRUE)
#' @param seed Integer seed for random number generator
#' @param verbose Verbose output
#' @return A new metadata column \code{is.pure} is added to the query Seurat object, indicating which cells correspond 
#' to the desidered cell type(s). The \code{active.ident} is also set to this variable.
#' @details Models for scGate are dataframes where each line is a signature for a given filtering level.
#' A database of models can be downloaded using the function \code{readRDS('data/models.rds')}.
#' You may directly use the models from the database, or edit one of these models to generate your own custom gating model.
#' Multiple models can also be evaluated at once, by running scGate with a list of models. Gating for each individual model is
#' returned as metadata, with a consensus annotation stored in \code{scGate_multi} metadata field. This allows using scGate as a
#' multi-class classifier, where only cells that are "Pure" for a single model are assigned a label, cells that are "Pure" for
#' more than one gating model are labeled as "Multi", all others cells are annotated as NA.
#' @examples
#' testing.datasets <- readRDS('data/testing.datasets.2k.rds')
#' seurat_object <- testing.datasets[["Satija"]]
#' # define basic gating model for B cells
#' my_scGate_model <- gating_model(name = "Bcell", signature = c("MS4A1"))
#' seurat_object <- scGate(data = seurat_object, model = my_scGate_model)
#' table(seurat_object$is.pure)
#' DimPlot(seurat_object)
#' # create a subsetted Seurat object with gated population
#' seurat_object_filtered <- subset(seurat_object, subset=is.pure=="Pure")
#' 
#' ############
#' # Using pre-defined models
#' scGate.model.db <- readRDS('data/models.rds')
#' seurat_object <- scGate(seurat_object, model=models$human$generic$PanBcell)
#' DimPlot(seurat_object)
#' seurat_object_filtered <- subset(seurat_object, subset=is.pure=="Pure")
#' 
#' ############
#' # Tree-like model visualization
#' models <- readRDS('data/models.rds')
#' plot_tree(models$human$generic$PanBcell) 
#' 
#' ############
#' # Using a manually edited model
#' my_scGate_model <- load_scGate_model("custom_model.tsv")
#' seurat_object <- scGate(seurat_object, model=my_scGate_model)
#' seurat_object_filtered <- subset(seurat_object, subset=is.pure=="Pure")
#' 
#' ############
#' # Run multiple models at once
#' models <- readRDS('data/models.rds')
#' model.list <- list("Bcell" = models$human$generic$Bcell,"Tcell" = models$human$generic$Tcell)
#' seurat_object <- scGate(seurat_object, model=model.list, verbose=T)
#' DimPlot(seurat_object, group.by = "scGate_multi")
#' 
#' ###########
#' # Run on pre-integrated space (e.g. using harmony)
#' models <- readRDS('data/models.rds')
#' seurat_object <- scGate(seurat_object, model=models$human$generic$Tcell, reduction="harmony")
#' DimPlot(seurat_object)
#' @import Seurat
#' @import ggplot2
#' @importFrom dplyr %>% distinct bind_rows
#' @export
# ---- scGate ----
scGate <- function(data,model,pos.thr=0.2,neg.thr=0.2,assay=NULL,slot="data",ncores=1,seed=123,keep.ranks=FALSE,
                   reduction=c("calculate","pca","umap","harmony","Liors_elephant"),min.cells=30,nfeatures=2000,
                   pca.dim=30,resol=3,param_decay=0.25,maxRank=1500,output.col.name='is.pure',
                   by.knn = TRUE,k.param=10,genes.blacklist="default",multi.asNA = FALSE,additional.signatures=NULL,
                   save.levels=FALSE,verbose=FALSE){
  set.seed(seed)
  if(!is.null(assay)){DefaultAssay(data) <- assay}
  assay <- DefaultAssay(data)
  if(assay == "integrated"){ #UCell should not run on integrated assay
    if('RNA' %in% Assays(data)){
      assay.ucell <- 'RNA'
    }else if('SCT' %in% Assays(data)){
      assay.ucell <- 'SCT'
    }else{
      stop("Cannot find assays with unintegrated data in this Seurat object")
    }
  }else{ assay.ucell <- assay }
  reduction <- reduction[1]
  if(is.null(reduction) || tolower(reduction)=="calculate"){
    reduction = "calculate"
  }else{
    if(!reduction %in% Reductions(data)){
      stop(sprintf("Could not find reduction %s in this object. Set reduction='calculate' to compute a new dimred", reduction))
    }
    pca.dim <- ncol(data@reductions[[reduction]])
  }
  #check gene blacklist
  if(!is.null(genes.blacklist)){
    if(length(genes.blacklist) == 1 && genes.blacklist == "default"){
      genes.blacklist <- readRDS('data/genes.blacklist.default.rds')
    }  
    if(is.list(genes.blacklist)){genes.blacklist <- unlist(genes.blacklist)}
    genes.blacklist <- unique(genes.blacklist)
  }
  #With single model, make a list of one element
  if(!inherits(model, "list")){ model <- list("Target" = model) }
  if(is.null(names(model))){ names(model) <- paste(output.col.name, seq_along(model), sep = ".") }
  # compute signature scores using UCell
  if(verbose){
    message(sprintf("Computing UCell scores for all signatures using %s assay...\n", assay.ucell))
  }
  data <- score.computing.for.scGate(
    data, model, ncores=ncores, assay = assay.ucell, slot = slot, maxRank=maxRank, 
    keep.ranks = keep.ranks, add.sign = additional.signatures
  )
  for(m in names(model)){
    col.id <- paste0(output.col.name, "_", m)
    data <- run_scGate_singlemodel(
      data, model = model[[m]], k.param = k.param, param_decay = param_decay, pca.dim = pca.dim, 
      resol = resol, nfeatures = nfeatures, by.knn = by.knn, min.cells = min.cells, assay = assay, 
      slot = slot, genes.blacklist = genes.blacklist, pos.thr = pos.thr, neg.thr = neg.thr, 
      verbose = verbose, reduction = reduction, colname = col.id, save.levels = save.levels
    )
    Idents(data) <- col.id
    n_rem <- sum(data[[col.id]]=="Impure")
    frac.to.rem <- n_rem/ncol(data)
    mess <- sprintf("\n### Detected a total of %i non-pure cells for %s (%.2f%% of total)", n_rem, m, 100*frac.to.rem)
    message(mess)
  }
  #Combine results from multiple model into single cell type annotation 
  data <- combine_scGate_multiclass(
    data, prefix = paste0(output.col.name,"_"), scGate_classes = names(m), multi.asNA = multi.asNA,
    min_cells = min.cells, out_column = "scGate_multi"
  )
  #Back-compatibility with previous versions
  if(names(model)[1] == 'Target'){
    cn <- paste0(output.col.name, "_Target")
    data@meta.data[,output.col.name] <- data@meta.data[,cn]
    data@meta.data[,cn] <- NULL
    if(save.levels){
      for(l in unique(model[[1]]$levels)){
        cn <- paste0(output.col.name, "_Target.",l)
        data@meta.data[,paste0(output.col.name,".",l)] <- data@meta.data[,cn]
        data@meta.data[,cn] <- NULL
      }
    }
  }
  return(data)
}


#' @title Plot model tree
#' @description View scGate model as a decision tree
#' @param model A scGate model to be visualized
#' @param box.size Box size
#' @param edge.text.size Edge text size
#' @return A plot of the model as a decision tree. At each level, green boxes indicate the 'positive' (accepted) 
#' cell types, red boxed indicate the 'negative' cell types (filtered out). The final Pure population is the 
#' bottom right subset in the tree.
#' @examples
#' models <- readRDS('data/models.rds')
#' plot_tree(model)
#' @export
plot_tree <- function(model, box.size = 8, edge.text.size = 4){
  if (!requireNamespace('ggparty', quietly = T)) {  #check whether ggparty is available
    stop("Please install and load package 'ggparty'")
  }else{ require(ggparty,quietly = T) }
  nlev <- length(unique(model$levels))
  if(nlev <= 1) stop("your model must contain at least two levels to be ploted as a tree")
  #restructure data for visualization
  level.list <- list()
  for(i in 1:nlev){
    level.list[[i]] <- list()
    sub <- subset(model, tolower(model$levels)==paste0("level",i))
    level.list[[i]][["positive"]] <- sub[sub$use_as=="positive","name"]
    level.list[[i]][["negative"]] <- sub[sub$use_as=="negative","name"]
  }
  #Initialize dataframe for tree
  df <- data.frame(matrix(ncol=nlev+1, nrow=nlev+1, data = 0))
  colnames(df) <- c(paste0("Level_", 1:nlev), "Pure")
  for(i in 2:nlev){ for(j in 1:(i-1)){ df[i,j] <- 1 } }
  df[nlev+1,] <- 1
  ##Construct tree structure
  pn <- list()
  #bottom level
  pn[[nlev]] <- partynode(
    nlev+1, split = partysplit(nlev, index=1:2, breaks = 0),
    kids = list(partynode(nlev+2),partynode(nlev+3))
  )
  for(i in (nlev-1):1){
    pn[[i]] <- partynode(
      i, split = partysplit(i, index=1:2, breaks=0),
      kids = list(partynode(i+1),pn[[i+1]])
    )
  }
  #first element in list has complete structure
  py <- party(pn[[1]], df)
  sign.annot <- vector(length=2*nlev+1)
  is.pos <- vector(length=2*nlev+1)
  sign.annot[1] <- "Root"
  is.pos[1] <- NA
  for(i in 1:nlev){
    sign.annot[2*i] <- paste0(level.list[[i]]$negative, collapse = "\n")
    sign.annot[2*i+1] <- paste0(level.list[[i]]$positive, collapse = "\n")
    is.pos[2*i] <- "Negative"
    is.pos[2*i+1] <- "Positive"
  }
  gg <- ggparty(py)
  gg$data$info <- sign.annot
  gg$data$p.value <- is.pos
  gg$data$breaks_label[grep("<=", gg$data$breaks_label)] <- "Negative"
  gg$data$breaks_label[grep(">", gg$data$breaks_label)] <- "Positive"
  gg <- gg + geom_edge() +
    geom_edge_label(size = edge.text.size) +
    geom_node_label(ids = "inner",mapping = aes(col = .data$p.value),
                    line_list = list(aes(label= .data$info)),line_gpar = list(list(size = box.size)))  +
    geom_node_label(ids = "terminal",mapping = aes(col = .data$p.value),nudge_y=0.01,
                    line_list = list(aes(label= .data$info)),line_gpar = list(list(size = box.size))) +
    scale_color_manual(values=c("#f60a0a", "#00ae60")) +
    theme(legend.position = "none", plot.margin = unit(c(1,1,1,1), "cm")) 
  return(gg)
}


#' @title Performance metrics
#' @description Evaluate model performance for binary tasks
#' @param actual Logical or numeric binary vector giving the actual cell labels. 
#' @param pred  Logical or numeric binary vector giving the predicted cell labels. 
#' @param return_contingency  Logical indicating if contingency table must be returned. Default is FALSE
#' @examples
#' results <- performance.metrics(actual= sample(c(1,0),20,replace =T),pred = sample(c(1,0),20,replace = T,prob = c(0.65,0.35)))
#' @export
performance.metrics <- function(actual,pred,return_contingency =F){
  actual <- as.numeric(actual + 0)
  pred <- as.numeric(pred +0)
  tp <- sum(actual&pred)
  tn <- sum((!actual)&(!pred))
  fn <- sum(actual&(!pred))
  fp <- sum((!actual)&pred)
  PREC <- tp/(tp +fp)
  REC <- tp/(tp + fn)
  #sqrt_ <- sqrt((tp + fp)*(tp+fn)*(tn+fp)*(tn+fn))
  sqrt_ <- exp(0.5* sum(log(c(tp+fp, tp+fn, tn+fp, tn+fn))))
  MCC <- (tp*tn - fp*fn) / sqrt_
  if(!return_contingency){  
    res.Summary <- c(PREC,REC,MCC); names(res.Summary) <- c("PREC","REC","MCC")
    return(res.Summary)
  }else{
    ct <- table(actual,pred)
    ## ordering contingency table, but avoiding errors when all predictions (or all actual cells) are equals
    nam.act <- unique(actual)%>%sort(decreasing = T)%>%as.character()  # 
    nam.pred <- unique(pred)%>%sort(decreasing = T)%>%as.character()
    ct <- ct[nam.act,nam.pred]  
    return(list('conting' = ct,'summary' = res.Summary ))
  }
}


#' @title Test your model
#' @description Wrapper for fast model testing on 3 sampled datasets 
#' @param model scGate model in data.frame format 
#' @param testing.version Character indicating the version of testing tatasets to be used. 
#' By default "hsa-latest" will be used. It will be ignored if custom.dataset variable is provied in Seurat format. 
#' Check available version in "https://figshare.com/XXXX/". 
#' @param custom.dataset Seurat object to be used as a testing dataset. For testing purposes, metadata seurat 
#' object must contain a column named 'cell_type' to be used as a gold standard. Also a set of positive targets 
#' must be provided in the target variable. 
#' @param target Positive target cell types. If default testing version is used this variable must be a character 
#' indicating one of the available target models 
#' ('immune','Lymphoid','Myeloid','Tcell','Bcell','CD8T','CD4T','NK','MoMacDC','Plasma_cell','PanBcell'). 
#' If a custom.dataset is provided in seurat format, this variable must be a vector of positive cell types in your data. 
#' The last case also require that such labels were named as in your cell_type meta.data column. 
#' @param plot Whether to return plots to device 
#' @examples
#' scGate.model.db <- readRDS('data/models.rds')
#' # Browse the list of models and select one:
#' model.panBcell <- scGate.model.db$human$generic$PanBcell
#' 
#' # Case 1: test the model with available testing datasets
#' panBcell.performance <- test_my_model(model.panBcell, testing.version = 'hsa.latest',target = "PanBcell")
#' model.Myeloid <-  scGate.model.db$human$generic$Myeloid
#' myeloid.performance <- test_my_model(model.Myeloid, testing.version = 'hsa.latest',target = "Myeloid")
#' 
#' # Case 2: use your own dataset for testing purposes. 
#' your.own.seurat.object <- readRDS(path.to.your.custom.dataset)
#' ## This make a copy of the cell_type reference column as required for scGate
#' your.own.seurat.object$cell_type <- your.own.seurat.object$your.gold.standard.column
#' # notice that 'Bcell' and 'PlasmaCell' must be celltypes present in your custom.dataset
#' performance.on.custom.dataset <- test_my_model(your.custom.PanBcell.model, 
#' custom.dataset = your.own.seurat.object, target = c("Bcell","PlasmaCell"))
#' @importFrom utils download.file
test_my_model <- function(model,testing.version='hsa.latest',custom.dataset=NULL,target=NULL, plot=T){
  performance.computation  <- ifelse(is.null(target), F, T)
  if(class(custom.dataset) == "Seurat"){
    testing.datasets <- list()
    testing.datasets$user.dataset <- custom.dataset
    custom = TRUE
  }else{ custom = FALSE }
  if(!custom){
    targets <- c('immune','Lymphoid','Myeloid','Tcell','Bcell','CD8T','CD4T','NK','MoMacDC','Plasma_cell','PanBcell')
    if(is.null(target)){
      message("warning: target cell_type not provided. Avoiding performance computation")  
      performace.computation = F
    }else if(!target %in% targets){
      stop(sprintf("target must be one of %s; or NULL for avoiding performance computation",paste(targets,collapse = "';'")))
    }
    ## check dataset version
    available.datasets = c("hsa.latest")
    if(!testing.version %in% available.datasets){
      stop("please, provide a valid testing.version paramter or provide a custom.dataset in seurat format")
    }
    # load testing datasets
    if(testing.version == "hsa.latest"){
      testing.datasets <- readRDS('data/testing.datasets.2k.rds')
    }
  }
  if(custom){
    if(!"cell_type" %in% colnames(custom.dataset@meta.data)){
      stop("please, provide a 'cell_type' column to be used as reference cell type")
    }
    if(is.null(target)){
      message("warning: target cell_type not provided. Avoiding performance computation")  
      performace.computation = F
    }else if(any(!target %in% custom.dataset$cell_type)){
      stop("all target celltypes must be included in cell_type metadata field. Otherwise, set target = NULL for avoiding performance computation")
    }
  }
  plt.out <- list()
  perf.out <- list()
  output <- list()
  # Drop is.pure cols if exists
  for(dset in names(testing.datasets)){
    obj <- testing.datasets[[dset]]
    plt <- list()
    dropcols = obj@meta.data %>% colnames() %>% grep("^is.pure",.,value =T) %>% unique()
    if(length(dropcols)>0){ for(col in dropcols){ obj[[col]] <- NULL } }
    ## scGate filtering
    obj <- scGate(obj, model = model, assay = DefaultAssay(obj))
    # add annotation plot
    nname <- sprintf("%s manual annot",dset)
    plt <- DimPlot(obj, group.by = "cell_type",label = T, repel =T,label.size = 3) + 
      ggtitle(nname) + NoLegend() +  theme(aspect.ratio = 1)
    # add one DimPlot by model level
    pure.plot <- DimPlot(obj, group.by = "is.pure", cols = list("Pure"="green","Impure"="gray")) +
      theme(aspect.ratio = 1)
    plt <- list("Annotation"=plt, "Gating"=pure.plot)
    #reserve plots of this dset
    plt.out[[dset]] <- patchwork::wrap_plots(plt,ncol = length(plt))
    if(performance.computation){
      if(!custom){    
        performance = performance.metrics(actual = obj@meta.data[,target], pred = obj$is.pure== "Pure")
      }else{
        performance = performance.metrics(actual = obj@cell_type %in% target, pred = obj$is.pure== "Pure")
      }
      perf.out[[dset]] <- performance 
    }
    output[[dset]] <- obj
  }
  if(performance.computation){
    perf <- Reduce(rbind,perf.out)
    rownames(perf) <- names(perf.out)
  }
  if(plot) {print(patchwork::wrap_plots(plt.out, ncol = 1))}
  if(performance.computation){
    return(list(performance = perf, plots = plt.out, objects = output))
  }else{
    return(list(plots = plt.out, objects = output))
  }
}


#' @title Plot scGate filtering results by level
#' @description Fast plotting of gating results over each model level.
#' @param obj Gated Seurat object output of scGate filtering function
#' @param pure.col Color code for pure category 
#' @param impure.col Color code for impure category 
#' @examples
#' scGate.model.db <- readRDS('data/models.rds')
#' # To see a specific model, browse the list of models:
#' scGate.model.db$human$generic$Myeloid
#' # Apply scGate with this model
#' query <- scGate(query, model=scGate.model.db$human$generic$Myeloid)
#' plot_levels(query)
#' @export
plot_levels <- function(obj, pure.col = "green", impure.col = "gray"){
  myCols <- grep("^is.pure", colnames(obj@meta.data),value = T)
  plots <- list()
  for(myCol in myCols){
    plots[[myCol]] <- DimPlot(obj, group.by = myCol, cols = list(Pure = pure.col,Impure = impure.col)) +
      theme(aspect.ratio = 1)
  }
  return(plots)
}


#' @title Plot UCell scores by level
#' @description Show distribution of UCell scores for each level of a given scGate model
#' @param obj Gated Seurat object (output of scGate)
#' @param model scGate model used to identify a target population in obj
#' @param pos.thr Threshold for positive signatures used in scGate model (set to NULL to disable)
#' @param neg.thr Threshold for negative signatures used in scGate model (set to NULL to disable) 
#' @param overlay Degree of overlay for ggridges
#' @param ncol Number of columns in output object (passed to wrap_plots)
#' @param combine Whether to combine plots into a single object, or to return a list of plots
#' @examples
#' scGate.model.db <- readRDS('data/models.rds')
#' model <- scGate.model.db$human$generic$Tcell
#' # Apply scGate with this model
#' query <- scGate(query, model=model, save.levels=T)
#' # View UCell score distribution
#' plot_UCell_scores(query, model)
#' @return Either a plot combined by patchwork (combine=T) or a list of plots (combine=F)
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom ggridges geom_density_ridges
#' @importFrom patchwork wrap_plots
#' @export
plot_UCell_scores <- function(obj,model,overlay=5,pos.thr=0.2,neg.thr=0.2,ncol=NULL,combine=T){
  u_cols <- grep('_UCell', colnames(obj@meta.data), value = T)
  levs <- unique(model$levels)
  pll <- list()
  palette <- c("#00fd0c","#f4340e")
  names(palette) <- c("Positive","Negative")
  if(sum(grepl("is.pure.level", colnames(obj@meta.data)))==0){
    obj$is.pure.level1 <- obj$is.pure
    if(length(levs)>1){
      warning("scGate levels were not stored in this object. Showing results only for top level.")
      levs <- "level1"
    }
  }
  for(l in seq_along(levs)){
    lev.name <- levs[l]
    sigs <- model[model$levels == lev.name, c("use_as","name")]
    col <- sprintf("%s_UCell", sigs$name)
    col <- col[col %in% u_cols]
    meta <- obj@meta.data
    if(l>1){ meta <- meta[meta[sprintf("is.pure.level%i",l-1)]=="Pure",] }
    ncells <- nrow(meta)
    stat <- table(meta[,sprintf("is.pure.level%i",l)])
    to.plot <- meta[,col, drop=FALSE]
    colnames(to.plot) <- gsub("_UCell","",colnames(to.plot))
    to.plot <- reshape2::melt(to.plot, id=NULL)
    colnames(to.plot) <- c("Signature","Score")
    to.plot$Class <- "Positive"
    to.plot$Class[to.plot$Signature %in% sigs[sigs$use_as =="negative","name"]] <- "Negative"
    #vertical lines (thresholds)
    to.plot$Thr <- NA
    if (!is.null(pos.thr)) { to.plot[to.plot$Class=="Positive","Thr"] <- pos.thr }
    if (!is.null(neg.thr)) { to.plot[to.plot$Class=="Negative","Thr"] <- neg.thr }
    #Make ggridges distribution plot
    pll[[l]] <- ggplot(to.plot, aes(x =.data$Score, y =.data$Signature, fill=.data$Class)) + 
      ggridges::geom_density_ridges(scale = overlay) +
      scale_fill_manual(values = palette) + theme_minimal() +
      theme(axis.title.y=element_blank()) + ggtitle(sprintf("%s - %i/%i pure ",lev.name, stat["Pure"], ncells)) 
    #Add threshold lines
    if(!is.null(pos.thr) | !is.null(neg.thr)){  
      pll[[l]] <- pll[[l]] + geom_segment(aes(x = .data$Thr, xend = .data$Thr,y = as.numeric(.data$Signature), 
                                              yend = as.numeric(.data$Signature)+0.9), linetype = "dashed")
    }
  }
  #Return combined plot or list of plots
  if(combine){
    return(patchwork::wrap_plots(pll, ncol=ncol))
  }else{ return(pll) }
}


#' @title Model creation and editing
#' @description Generate an scGate model from scratch or edit an existing one
#' @param model scGate model to be modified. When is NULL (default) a new model will be initialized.   
#' @param level integer. It refers to the hierarchical level of the model tree in which the signature will be added (level=1 by default)    
#' @param name Arbitrary signature name (i.e. Immune, Tcell, NK etc).   
#' @param signature character vector indicating gene symbols to be included in the signature (e.g. CD3D). If a minus sign is placed to the end of a gene name (e.g. "CD3D-"), this gene will be used as negative in UCell computing. See UCell documentation for details    
#' @param positive Logical indicating if the signature must be used as a positive signature in those model level. Default is TRUE. 
#' @param negative Same as `positive` but negated (negative=TRUE equals to positive=FALSE)
#' @param remove Whether to remove the given signature from the model
#' @examples
#' # create a simple gating model
#' my_model <- gating_model(model = my_model,level = 1, positive = T, name = "immune", signature = c("PTPRC"))
#' my_model <- gating_model(model = my_model,level = 1, positive = F, name = "Epithelial", signature = c("CDH1","FLT1") )
#' # Remove an existing signature
#' dropped_model <- gating_model(model = my_model, remove =TRUE, level = 1, name = "Epithelial")
#' @importFrom stats setNames
gating_model <- function(model=NULL,level= 1,name,signature,positive = T,negative = F,remove = F){
  template <- setNames(data.frame(matrix(ncol = 4,nrow = 0)),c("levels","use_as","name","signature"))
  if(negative) positive <-  F
  if(is.null(model)) model <- template 
  if(!remove){
    new.signature <- data.frame(
      levels = paste0("level",level), use_as = ifelse(positive, "positive","negative"),
      name = name,signature = ifelse(length(signature) > 1, paste(signature,collapse = ";"),signature)
    )
    model <- rbind(model,new.signature)
  }else{
    L <- paste0("level",level)
    model <- model[!((model$levels == L) & (model$name == name)),]
  }
  return(model)
}


#' @title Combine scGate annotations
#' @description If a single-cell dataset has precomputed results for multiple scGate models, combined them in multi-class annotation
#' @param obj Seurat object with scGate results for multiple models stored as metadata   
#' @param prefix Prefix in metadata column names for scGate result models
#' @param scGate_classes Vector of scGate model names. If NULL, use all columns that start with "prefix" above.
#' @param min_cells Minimum number of cells for a cell label to be considered
#' @param multi.asNA How to label cells that are "Pure" for multiple annotations: "Multi"(FALSE) or NA(TRUE)
#' @param out_column The name of the metadata column where to store the multi-class cell labels
#' @examples
#' obj <- combine_scGate_multiclass(obj)
#' @import Seurat
combine_scGate_multiclass <- function(obj,prefix="is.pure_",scGate_classes=NULL,min_cells=20,
                                      multi.asNA = FALSE,out_column="scGate_multi"){
  #Use all columns with given prefix
  if(is.null(scGate_classes)){  
    cols <- grep(prefix, colnames(obj@meta.data), value = T)
    cols <- grep("\\.level\\d+$", cols, invert=T, perl=T, value=T)
  }else{
    cols <- paste0(prefix, scGate_classes)
    cols <- cols[cols %in% colnames(obj@meta.data)]
  }
  if(is.null(cols)){
    stop("Could not find scGate annotations in this object metadata.")
  }
  meta <- obj@meta.data[,cols, drop=FALSE]
  meta[is.na(meta)] <- "Impure"  #Avoid NAs
  obj.logical <- meta=="Pure"
  label.sums <- apply(obj.logical,1,sum)
  obj.single <- obj.logical[label.sums==1, , drop=FALSE]
  obj.single.labels <- apply(obj.single,1,function(x) names(x)[x])
  #remove prefix
  if(!is.null(prefix)){ obj.single.labels <- gsub(prefix, "", obj.single.labels) }
  #Assign labels to uniquely identified cells
  labs <- rep(NA, ncol(obj))
  names(labs) <- colnames(obj)
  labs[names(obj.single.labels)] <- obj.single.labels
  #Set to NA classes with too few cells
  tt <- table(labs, useNA = "always")
  labs[labs %in% names(tt)[tt<min_cells]] <- NA
  if(multi.asNA){
    labs[names(label.sums[label.sums>1])] <- NA
  }else{ 
    labs[names(label.sums[label.sums>1])] <- "Multi"
  }
  obj@meta.data[[out_column]] <- labs
  obj
}


run_scGate_singlemodel <- function(data, model, pos.thr=0.2, neg.thr=0.2, assay=NULL, slot="data",
                                   reduction="calculate", nfeatures=2000, pca.dim=30, resol=3,
                                   param_decay=0.25, min.cells=30, by.knn = TRUE, k.param=10, 
                                   genes.blacklist="default", verbose=FALSE,
                                   colname="is.pure", save.levels = FALSE){
  if(!inherits(model, "data.frame")){
    stop("Invalid scGate model. Please check the format of your model")
  }
  list.model <- table.to.model(model)
  q <- data  #local copy to progressively remove cells
  tot.cells <- ncol(q)
  ## prepare output object (one pure/impure flag by level)
  output_by_level <- rep("Impure",length(list.model)*tot.cells)
  dim(output_by_level) <- c(tot.cells,length(list.model))
  colnames(output_by_level) <- names(list.model)
  output_by_level <- data.frame(output_by_level)
  rownames(output_by_level) <- data@meta.data %>% rownames()
  for(lev in 1:length(list.model)){
    if(verbose){message(sprintf("Running scGate on level %i...", lev))}
    pos.names <- sprintf("%s_UCell", names(list.model[[lev]]$positive))
    neg.names <- sprintf("%s_UCell", names(list.model[[lev]]$negative))
    ##Reduce parameter complexity at each iteration
    if(param_decay < 0 | param_decay > 1){
      stop("Parameter param_decay must be a number between 0 and 1")
    }
    if(reduction=="calculate"){
      pca.use <- round((1-param_decay)**(lev-1) * pca.dim)
      nfeat.use <- round((1-param_decay)**(lev-1) * nfeatures)
      res.use <- round((1-param_decay)**(lev-1) * resol)
    }else{pca.use <- pca.dim}
    q <- find.nn(q, by.knn=by.knn, assay=assay, slot=slot, min.cells=min.cells, nfeatures=nfeat.use, 
                 reduction=reduction, npca=pca.use, k.param=k.param, genes.blacklist=genes.blacklist)
    if(!by.knn){
      q  <- FindClusters(q, resolution = res.use, verbose = FALSE)
      q$clusterCT <- q@active.ident
    }
    q <- filter_bymean(q, positive=pos.names, negative=neg.names, pos.thr=pos.thr, assay=assay,
                       min.cells=min.cells, neg.thr=neg.thr, by.knn = by.knn)
    n_rem <- sum(q$is.pure=="Impure")
    frac.to.rem <- n_rem/tot.cells
    mess <- sprintf("scGate: Detected %i non-pure cells at level %i", n_rem, lev)
    if(verbose){ message(mess) }
    ## retain pure cells will give us info in case of all cell where filtered
    retain_pure_cells <- q$is.pure=="Pure"
    if(any(retain_pure_cells)){
      output_by_level[names(retain_pure_cells[retain_pure_cells==T]),lev] <- "Pure"  # save layer output
      q <- subset(q, subset=`is.pure`=="Pure")
    }else{
      break  # in case of all cells became filtered, we do not continue with the next layer
    }
  }
  #Add 'pure' labels to metadata
  data <- AddMetaData(data,col.name = colname, metadata = rep("Impure",tot.cells))
  if(any(retain_pure_cells)){  
    pure.cells <- colnames(q)
    data@meta.data[pure.cells, colname] <- "Pure"
  }else{
    message(sprintf("Warning, all cells were removed at level %i. Consider reviewing signatures or model layout...", lev))
  }
  data@meta.data[,colname] <- factor(data@meta.data[,colname], levels=c("Pure","Impure"))
  # Save output by level
  if(save.levels){
    for(name.lev in names(list.model)){
      combname <- paste0(colname,".",name.lev)
      data <- AddMetaData(data,col.name = combname, metadata = output_by_level[[name.lev]])
      data@meta.data[,combname] <- factor(data@meta.data[,combname], levels=c("Pure","Impure"))
    }
  }
  return(data)
}

find.nn <- function(q, assay = "RNA", slot="data", npca=30, nfeatures=2000, k.param=10,
                    min.cells=30,by.knn = F,reduction="calculate",genes.blacklist=NULL){
  DefaultAssay(q) <- assay
  ncells <- length(Cells(q))
  ngenes <- nrow(q)
  if(reduction=="calculate"){
    if(ncells < min.cells){
      q$clusterCT <- 0    #with very few cells, consider them as a single cluster
      return(q)
    }  
    if(ngenes < nfeatures){ nfeatures <- ngenes}
    if(ngenes/2 < npca){ npca <- ngenes/2 }
    if(slot=="counts"){ q <- NormalizeData(q, verbose = FALSE) }
    if(ngenes > 200){  #only perform this for high-dim data
      q <- FindVariableFeatures(q, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    }else{ q@assays[[assay]]@var.features <- rownames(q) }
    q@assays[[assay]]@var.features <- setdiff(q@assays[[assay]]@var.features, genes.blacklist)
    q <- ScaleData(q, verbose=FALSE)
    q <- RunPCA(q,features=q@assays[[assay]]@var.features,npcs=npca,verbose=FALSE,reduction.key="knnPCA_")
    red.use <- 'pca'
  }else{ red.use <- reduction }
  q <- suppressMessages(
    FindNeighbors(q, reduction = red.use, dims = 1:npca, k.param = k.param, verbose=FALSE,
                  return.neighbor = by.knn)
  )
  return(q)
}

## Filter by mean
filter_bymean <- function(q, positive, negative, pos.thr=0.1, neg.thr=0.2,min.cells=30,
                          assay="RNA", return_object = T, by.knn = F){
  DefaultAssay(q) <- assay
  ncells <- dim(q)[2]
  notfound <- c(positive[!positive %in% colnames(q@meta.data)], negative[!negative %in% colnames(q@meta.data)])
  if(length(notfound)>0){
    message(paste0("Warning: signatures not found: ", notfound))
  }
  positive <- positive[positive %in% colnames(q@meta.data)]
  negative <- negative[negative %in% colnames(q@meta.data)]
  cols <- c(positive, negative)
  means <- list()
  if(!by.knn){
    for(col in cols){
      means[[col]] <- sapply(levels(q$clusterCT), function(x) {
        mean(q@meta.data[q$clusterCT == x, col])
      })
    }
    meds <- Reduce(rbind, means)
    if(is.null(dim(meds))){dim(meds) <- c(1,length(meds))}
    rownames(meds) <- cols
    pos <- vector(length=dim(meds)[2])
    neg <- vector(length=dim(meds)[2])
    for(j in 1:dim(meds)[2]){
      pos[j] <- max(meds[positive,j])
      neg[j] <- max(meds[negative,j])
    }
    if(length(negative)>0){
      indices <- intersect(which(pos > pos.thr), which(neg < neg.thr))
    }else{
      indices <- which(pos > pos.thr) # case without negative signature
    }
    select.pures <- colnames(meds)[indices]
    ispure <- ifelse(q$clusterCT %in% select.pures,"Pure","Impure")
  }else{
    for(col in cols){
      meta.nn <- sprintf("%s.nn", assay)
      if(ncells < min.cells){   #very small dataset. Use all cells together
        neigs <- t(matrix(data = rep(1:ncells,ncells), nrow = ncells, ncol = ncells))
      }else{ neigs <- q@neighbors[[meta.nn]]@nn.idx }
      m <- q[[col]][[1]][neigs]
      dim(m) <- dim(neigs)
      means[[col]] <- apply(m,1,mean)
    }
    meds <- Reduce(rbind, means)
    if(is.null(dim(meds))){
      dim(meds) <- c(1,length(meds))
    }
    rownames(meds) <- cols
    if(length(positive)>1){
      pos <- meds[positive,]%>%apply(2,max)
    }else{
      pos <- meds[positive,]
    }
    if(length(negative)>1){
      neg <- meds[negative,]%>%apply(2,max)
    }else{
      neg<- meds[negative,]
    }
    ispure <- rep("Impure",dim(q)[2])
    if(length(negative)>0){
      ispure[(pos > pos.thr)&(neg < neg.thr)] <- "Pure"
    }else{
      ispure[pos > pos.thr] <- "Pure"  # case without negative signature
    }
  }
  q$is.pure <- ispure
  if(return_object) {return(q)}
  return(ispure)
}

score.computing.for.scGate <- function(data, model, ncores=1, assay="RNA", slot="data",
                                       add.sign=NULL, keep.ranks=FALSE, maxRank=1500){
  comb <- dplyr::bind_rows(model, .id = "Model_ID")
  # extract unique signatures
  model.uniq <- comb %>% dplyr::distinct(.data$name, .data$signature, .keep_all = T) 
  #Stop if there are signature with same name but different genes
  t <- table(model.uniq$name)
  dup <- names(t[t>1])
  if(length(dup)>0){
    s <- paste(dup,collapse=", ")
    stop(sprintf("Different gene sets have been assigned to signature with same name: %s", s))
  }
  ## generate list object to be used in computing stage
  signatures <- model.uniq$signature %>% strsplit("[,; ]+") %>% lapply(unlist)
  names(signatures) <- model.uniq$name
  if(!is.null(add.sign)){
    if (!inherits(add.sign, "list")){add.sign <- list("Additional_signature"=add.sign)}
    signatures <- append(signatures, add.sign)
  }
  data <- AddModuleScore_UCell(data, features = signatures, assay=assay, slot=slot,
                               ncores=ncores, storeRanks = keep.ranks, maxRank = maxRank)
  return(data)
}

model.to.table <- function(scGate.model){
  tab <- data.frame(levels=NULL,use_as = NULL,name = NULL,signature = NULL)
  ### Define first column: Levels
  levels <- scGate.model%>%names()
  lev <- rep(levels,rep(2,length(levels)))
  len_pos_neg <- lapply(scGate.model,function(x){ 
    res = lapply(x,function(y){length(y)})
    return(res)
  }) %>% unlist()
  extended.levels <- rep(lev,len_pos_neg)
  # second column: "use_as"
  useas <- rep(c("positive","negative"),length(levels))
  useas <- rep(useas , len_pos_neg)
  #third column: name
  signature.names <- lapply(scGate.model,function(x){ 
    res = lapply(x,function(y){names(y)})
    return(res)
  }) %>% unlist()
  ## Four column: signature
  signatures <- lapply(scGate.model,function(x){ 
    res = lapply(x,function(y){
      lapply(y, function(z){paste(z,collapse = ";")})})
    return(res)
  }) %>% unlist()
  tab <- data.frame("levels"=extended.levels,"use_as" = useas, "name" = signature.names,"signature" = signatures)
  return(tab)
}

table.to.model <- function(scGate.table){
  mod <- list()
  for(i in 1:nrow(scGate.table)){ 
    lev <- scGate.table$levels[i] 
    useas <- tolower(scGate.table$use_as[i])
    if(!useas %in% c("positive","negative")){
      message(sprintf("Error: row %i do not contain neither, 'positive' or 'negative' strings in 'use_as' column",i))
      return(NULL)
    }
    sign <- scGate.table$signature[i]
    name <- scGate.table$name[i]
    mod[[lev]][[useas]][[name]] <- strsplit(sign, "[,; ]+") %>% unlist()
  }
  return(mod)
}



#' @title Calculate module enrichment scores from single-cell data
#' Given a Seurat object, calculates module/signature enrichment scores at
#' single-cell level using the Mann-Whitney U statistic.
#' UCell scores are normalized U statistics (between 0 and 1), and they are
#' mathematically related to the Area under the ROC curve
#' (see [Mason and Graham](https://doi.org/10.1256/003590002320603584))
#' In contrast to Seurat's AddModuleScore, which is normalized by binning genes
#' of similar expression at the population level, UCell scores depend 
#' only on the gene expression ranks of individual cell, and therefore they are
#' robust across datasets regardless of dataset composition.
#' @param obj Seurat object
#' @param features A list of signatures, for example:
#' \code{list(Tcell_signature = c("CD2","CD3E","CD3D"),Myeloid_signature = c("SPI1","FCER1G","CSF1R"))}
#' You can also specify positive and negative gene sets by adding a + or -
#' sign to genes in the signature; see an example below
#' @param maxRank Maximum number of genes to rank per cell; above this rank, a given gene is considered as not expressed.
#' @param storeRanks Store ranks matrix in Seurat object ('UCellRanks' assay) for fast subsequent computations. 
#' This option may demand large amounts of RAM.
#' @param w_neg Weight on negative genes in signature. e.g. `w_neg=1` weighs equally up- and down-regulated genes,
#'`w_neg=0.5` gives 50% less importance to negative genes
#' @param assay Pull out data from this assay of the Seurat object (if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot of the Seurat object
#' @param chunk.size Number of cells to be processed simultaneously (lower size requires slightly more computation but reduces memory demands)
#' @param BPPARAM A [BiocParallel::bpparam()] object that tells UCell how to parallelize. If provided, it overrides the `ncores` parameter.     
#' @param ncores Number of processors to parallelize computation. If
#' \code{BPPARAM = NULL}, the function uses \code{BiocParallel::MulticoreParam(workers=ncores)}
#' @param ties.method How ranking ties should be resolved - passed on to [data.table::frank]
#' @param force.gc Explicitly call garbage collector to reduce memory footprint
#' @param name Name tag that will be appended at the end of each signature name, "_UCell" by 
#' default (e.g. signature score in meta data will be named: Myeloid_signature_UCell)
#' @return Returns a Seurat object with module/signature enrichment scores added to object meta data; 
#' each score is stored as the corresponding signature name provided in \code{features} followed by the tag given
#' in \code{name} (or "_UCell" by default )
#' @examples
#' gene.sets <- list(Tcell = c("CD2","CD3E","CD3D"), Myeloid = c("SPI1","FCER1G","CSF1R"))
#' load('data/sample.matrix.RData')
#' obj <- Seurat::CreateSeuratObject(sample.matrix)                
#' obj <- AddModuleScore_UCell(obj,features = gene.sets)
#' head(obj[[]])
#' 
#' ## Using positive and negative gene sets
#' gene.sets <- list()
#' gene.sets$Tcell_gd <- c("TRDC+","TRGC1+","TRGC2+","TRDV1+", "TRAC-","TRBC1-","TRBC2-")
#' gene.sets$NKcell <- c("FGFBP2+", "SPON2+", "KLRF1+", "FCGR3A+", "CD3E-","CD3G-")
#' obj <- AddModuleScore_UCell(obj, features = gene.sets, name=NULL)
#' head(obj$NKcell)
# ---- UCell ----
AddModuleScore_UCell <- function(obj,features,maxRank=1500,chunk.size=1000,BPPARAM=NULL,ncores=1,storeRanks=FALSE,
                                 w_neg=1,assay=NULL,slot="data",ties.method="average",force.gc=FALSE,name="_UCell"){
  features <- check_signature_names(features)
  if(is.null(assay)) assay <- Seurat::DefaultAssay(obj)
  # If rank matrix was pre-computed, evaluate the new signatures
  # from these ranks. Else, calculate new ranks to score signatures
  # (optionally storing ranks, takes up more memory but become very
  # fast to evaluate further signatures)
  if("UCellRanks" %in% Seurat::Assays(obj)){
    meta.list <- rankings2Uscore(
      Seurat::GetAssayData(obj, "counts", assay="UCellRanks"),
      features=features, chunk.size=chunk.size, w_neg=w_neg,
      ncores=ncores, BPPARAM=BPPARAM, force.gc=force.gc, name=name)
  }else{
    meta.list <- calculate_Uscore(
      Seurat::GetAssayData(obj, slot, assay=assay),
      features = features, maxRank = maxRank,chunk.size = chunk.size, w_neg = w_neg,
      ncores = ncores, BPPARAM = BPPARAM, ties.method = ties.method,
      force.gc = force.gc, storeRanks = storeRanks, name = name)
    #store ranks matrix?
    if (storeRanks==TRUE){
      cells_rankings.merge <- lapply(meta.list,function(x) rbind(x[["cells_rankings"]]))
      cells_rankings.merge <- Reduce(cbind, cells_rankings.merge)
      obj[["UCellRanks"]] <- Seurat::CreateAssayObject(cells_rankings.merge)
    }
  }
  meta.merge <- lapply(meta.list,function(x) rbind(x[["cells_AUC"]]))
  meta.merge <- Reduce(rbind, meta.merge)
  obj <- Seurat::AddMetaData(obj, as.data.frame(meta.merge))
  return(obj)
}


#' @title Calculate module enrichment scores from single-cell data
#' Given a gene vs. cell matrix, calculates module/signature enrichment scores
#' on single-cell level using Mann-Whitney U statistic.
#' UCell scores are normalized U statistics (between 0 and 1), and they are
#' mathematically related to the Area under the ROC curve (see
#' [Mason and Graham](https://doi.org/10.1256/003590002320603584))
#' These scores only depend on the gene expression ranks of individual cell, 
#' and therefore they are robust across datasets regardless of dataset composition.
#' @param matrix Input matrix, either stored in a [SingleCellExperiment] object or as a raw matrix. \code{dgCMatrix} format supported.
#' @param features A list of signatures, for example:
#' \code{list(Tcell_signature = c("CD2","CD3E","CD3D"), Myeloid_signature = c("SPI1","FCER1G","CSF1R"))}
#' You can also specify positive and negative gene sets by adding a + or - sign to genes in the signature; see an example below
#' @param precalc.ranks If you have pre-calculated ranks using
#' \code{\link{StoreRankings_UCell}}, you can specify the pre-calculated ranks instead of the gene vs. cell matrix.
#' @param maxRank Maximum number of genes to rank per cell; above this rank, a 
#' given gene is considered as not expressed. Note: this parameter is ignored if \code{precalc.ranks} are specified
#' @param assay The sce object assay where the data is to be found
#' @param chunk.size Number of cells to be processed simultaneously (lower size requires slightly more computation but reduces memory demands)
#' @param w_neg Weight on negative genes in signature. e.g. `w_neg=1` weighs equally up- and down-regulated genes, 
#' `w_neg=0.5` gives 50% less importance to negative genes
#' @param BPPARAM A [BiocParallel::bpparam()] object that tells UCell how to parallelize. If provided, it overrides the `ncores` parameter.     
#' @param ncores Number of processors to parallelize computation. If
#' \code{BPPARAM = NULL}, the function uses \code{BiocParallel::MulticoreParam(workers=ncores)}
#' @param ties.method How ranking ties should be resolved - passed on to [data.table::frank]
#' @param name Name suffix appended to signature names
#' @param force.gc Explicitly call garbage collector to reduce memory footprint
#' @return Returns input SingleCellExperiment object with UCell scores added to altExp
#' @examples
#' # Using sparse matrix
#' load('data/sample.matrix.RData')
#' gene.sets <- list(Tcell_signature = c("CD2","CD3E","CD3D"), Myeloid_signature = c("SPI1","FCER1G","CSF1R"))
#' scores <- ScoreSignatures_UCell(sample.matrix, features=gene.sets)
#' head(scores)
#' 
#' # Using sce object
#' library(SingleCellExperiment)
#' load('data/sample.matrix.RData')
#' my.sce <- SingleCellExperiment(list(counts=sample.matrix))
#' gene.sets <- list(Tcell_signature = c("CD2","CD3E","CD3D"), Myeloid_signature = c("SPI1","FCER1G","CSF1R"))
#' my.sce <- ScoreSignatures_UCell(my.sce, features=gene.sets)
#' altExp(my.sce, 'UCell')
#' @importFrom methods is 
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay assays SummarizedExperiment
#' @import Matrix
#' @export
ScoreSignatures_UCell <- function(matrix=NULL,features, precalc.ranks=NULL,maxRank=1500, w_neg=1, name="_UCell",
                                  assay="counts", chunk.size=1000,BPPARAM=NULL, ncores=1,ties.method="average", 
                                  force.gc=FALSE){
  features <- check_signature_names(features)
  #Check type of input
  if(is(matrix, "SingleCellExperiment")){ # sce object
    if(!assay %in% names(SummarizedExperiment::assays(matrix))){
      stop(sprintf("Assay %s not found in sce object.", assay))
    }
    m <- SummarizedExperiment::assay(matrix, assay) 
  }else if(is(matrix, "matrix") | is(matrix, "dgCMatrix") | is(matrix, "data.frame")) { 
    m <- matrix
  }else{
    m <- NULL
  }
  if(is.null(m) & is.null(precalc.ranks)){ stop("Unrecognized input format.") }
  #Run on pre-calculated ranks ('m' can be NULL)
  if(!is.null(precalc.ranks)){
    u.list <- rankings2Uscore(precalc.ranks, features=features,chunk.size=chunk.size,w_neg=w_neg,
                              ncores=ncores, BPPARAM=BPPARAM, force.gc=force.gc, name=name)
  }else{
    u.list <- calculate_Uscore(m, features=features, maxRank=maxRank,chunk.size=chunk.size, w_neg=w_neg,
                               ties.method=ties.method, ncores=ncores,BPPARAM=BPPARAM, force.gc=force.gc,name=name)
  }
  u.merge <- lapply(u.list,function(x) rbind(x[["cells_AUC"]]))
  u.merge <- Reduce(rbind, u.merge)
  if(is(matrix, "SingleCellExperiment")){
    SingleCellExperiment::altExp(matrix, "UCell") <- SummarizedExperiment(assays = list("UCell" = t(u.merge)))
    return(matrix)
  }else{
    return(u.merge)
  }
}


#' @title Calculate and store gene rankings for a single-cell dataset
#' Given a gene vs. cell matrix, calculates the rankings of expression for all
#' genes in each cell. 
#' While \code{\link{ScoreSignatures_UCell}} can be used 'on the fly' to
#' evaluate signatures in a query dataset, it requires recalculating gene
#' ranks at every execution. If you have a large dataset and plan to experiment
#' with multiple signatures, evaluating the same dataset multiple times,
#' this function allows you to store pre-calculated ranks so they do not have to
#' be recomputed every time. Pre-calculated ranks can then be applied to the
#' function \code{\link{ScoreSignatures_UCell}} to evaluate gene signatures in a
#' significantly faster way on successive iterations.
#' @param matrix Input matrix, either stored in a [SingleCellExperiment] object
#' or as a raw matrix. \code{dgCMatrix} format supported.
#' @param maxRank Maximum number of genes to rank per cell; above this rank, a given gene is considered as not expressed
#' @param assay Assay where the data is to be found (for input in 'sce' format)
#' @param chunk.size Number of cells to be processed simultaneously (lower size requires slightly more computation but reduces memory demands)
#' @param BPPARAM A [BiocParallel::bpparam()] object that tells UCell how to parallelize. If provided, 
#' it overrides the `ncores` parameter.     
#' @param ncores Number of processors to parallelize computation. If
#' \code{BPPARAM = NULL}, the function uses \code{BiocParallel::MulticoreParam(workers=ncores)}
#' @param ties.method How ranking ties should be resolved - passed on to [data.table::frank]
#' @param force.gc Explicitly call garbage collector to reduce memory footprint
#' @return Returns a sparse matrix of pre-calculated ranks that can be used multiple times to evaluate different signatures
#' @examples
#' load('data/sample.matrix.RData')
#' ranks <- StoreRankings_UCell(sample.matrix)
#' ranks[1:5,1:5]
#' gene.sets <- list(Tcell_signature = c("CD2","CD3E","CD3D"),Myeloid_signature = c("SPI1","FCER1G","CSF1R"))
#' scores <- ScoreSignatures_UCell(features=gene.sets, precalc.ranks=ranks)
#' head(scores)
#' @importFrom methods is 
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @import Matrix
#' @export
StoreRankings_UCell <- function(matrix,maxRank=1500,chunk.size=1000,BPPARAM=NULL,ncores=1,assay='counts',
                                ties.method="average", force.gc=FALSE){
  #Check type of input
  if(is(matrix, "SingleCellExperiment")){ # sce object
    if(!assay %in% names(SummarizedExperiment::assays(matrix))) {
      stop(sprintf("Assay %s not found in sce object.", assay))
    }
    m <- SummarizedExperiment::assay(matrix, assay)
  }else if(is(matrix, "matrix") | is(matrix, "dgCMatrix") | is(matrix, "data.frame")) { 
    m <- matrix
  }else{ stop("Unrecognized input format.") }
  features <- rownames(m)[1]  #placeholder signature
  meta.list <- calculate_Uscore(m, features=features, maxRank=maxRank,chunk.size=chunk.size, ncores=ncores, 
                                BPPARAM=BPPARAM, ties.method=ties.method, storeRanks=TRUE,force.gc=force.gc)
  ranks.all <- lapply(meta.list,function(x) rbind(x[["cells_rankings"]]))
  ranks.all <- Reduce(cbind, ranks.all)
  return(ranks.all)
}


#' @title Smooth signature scores by kNN
#' This function performs smoothing of single-cell scores by weighted
#' average of the k-nearest neighbors. It can be useful to 'impute' scores by
#' neighboring cells and partially correct data sparsity. While this function
#' has been designed to smooth UCell scores, it can be applied to any numerical
#' metadata contained in SingleCellExperiment or Seurat objects
#' @param obj Input object - either a [SingleCellExperiment] object or a Seurat object.
#' @param signature.names The names of the signatures (or any numeric metadata column) for which to calculate kNN-smoothed scores
#' @param reduction Which dimensionality reduction to use for kNN smoothing. It must be already present in the input object.
#' @param k Number of neighbors for kNN smoothing
#' @param BNPARAM A [BiocNeighborParam] object specifying the algorithm to use for kNN calculation.
#' @param BPPARAM A [BiocParallel::bpparam()] object for parallel computing, e.g. [MulticoreParam] or [SnowParam]
#' @param suffix Suffix to append to metadata columns for the new knn-smoothed scores  
#' @param assay For Seurat objects only - do smoothing on expression data from this assay. When NULL, only looks in metadata
#' @param slot For Seurat objects only - do smoothing on expression data from this slot
#' @param sce.expname For sce objects only - which experiment stores the signatures to be smoothed. Set to 'main' 
#' for smoothing gene expression stored in the main sce experiment.
#' @param sce.assay For sce objects only - pull data from this assay
#' @return An augmented \code{obj} with the smoothed signatures. If \code{obj}
#' is a Seurat object, smoothed signatures are added to metadata; if\code{obj} is a SingleCellExperiment object, 
#' smoothed signatures are returned in a new altExp. See the examples below.     
#' @examples
#' #### Using Seurat ####
#' library(Seurat)
#' gene.sets <- list(Tcell = c("CD2","CD3E","CD3D"), Myeloid = c("SPI1","FCER1G","CSF1R"))
#' load('data/sample.matrix.RData')
#' obj <- Seurat::CreateSeuratObject(sample.matrix)                
#' # Calculate UCell scores
#' obj <- AddModuleScore_UCell(obj,features = gene.sets, name=NULL)
#' # Run PCA
#' obj <- FindVariableFeatures(obj) |> ScaleData() |> RunPCA()
#' # Smooth signatures
#' obj <- SmoothKNN(obj, reduction="pca", signature.names=names(gene.sets))
#' head(obj[[]])
#' 
#' #### Using SingleCellExperiment ####
#' library(SingleCellExperiment)
#' library(scater)
#' load('data/sample.matrix.RData')
#' sce <- SingleCellExperiment(list(counts=sample.matrix))
#' gene.sets <- list( Tcell = c("CD2","CD3E","CD3D"), Myeloid = c("SPI1","FCER1G","CSF1R"))
#' # Calculate UCell scores
#' sce <- ScoreSignatures_UCell(sce, features=gene.sets, name=NULL)
#' # Run PCA
#' sce <- logNormCounts(sce)
#' sce <- runPCA(sce, scale=TRUE, ncomponents=20)
#' # Smooth signatures
#' sce <- SmoothKNN(sce, reduction="PCA", signature.names=names(gene.sets))
#' # See results
#' altExp(sce, 'UCell')
#' assays(altExp(sce, 'UCell'))
#' # Plot on UMAP
#' sce <- runUMAP(sce, dimred="PCA")
#' plotUMAP(sce, colour_by = "Tcell_kNN", by_exprs_values = "UCell_kNN")
#' 
#' @importFrom methods setMethod setGeneric is
#' @import BiocNeighbors
#' @importFrom Matrix t
#' @export SmoothKNN
SmoothKNN <- function(obj=NULL,signature.names=NULL,reduction="pca", k=10, BNPARAM=AnnoyParam(),BPPARAM=SerialParam(),
                      suffix="_kNN",assay=NULL,slot="data",sce.expname=c("UCell","main"),sce.assay=NULL){
  UseMethod("SmoothKNN")
}

#' @title Calculate Mann Whitney U from a vector of ranks
#' @param rank_value A vector of ranks
#' @param maxRank Max number of features to include in ranking
#' @param sparse Whether the vector of ranks is in sparse format
#' @return Normalized AUC (as U statistic) for the vector
u_stat <- function(rank_value, maxRank=1000, sparse=FALSE){
  if(sparse==TRUE){ rank_value[rank_value==0] <- maxRank+1 }
  insig <- rank_value > maxRank
  if(all(insig)){
    return(0L)
  }else{
    rank_value[insig] <- maxRank + 1
    rank_sum <- sum(rank_value)
    len_sig <- length(rank_value)
    u_value <- rank_sum - (len_sig * (len_sig + 1))/2
    auc <- 1 - u_value/(len_sig * maxRank)
    return(auc)
  }
}

#' @title Calculate U scores for a list of signatures, given a rank matrix
#' @param sig_list A list of signatures
#' @param ranks_matrix Matrix of pre-computed ranks
#' @param maxRank Max number of features to include in ranking, for u_stat function
#' @param sparse Whether the vector of ranks is in sparse format
#' @param w_neg Weight on negative signatures
#' @return A matrix of U scores
#' @import data.table
u_stat_signature_list <- function(sig_list, ranks_matrix, maxRank=1000,sparse=FALSE, w_neg=1){
  dim <- ncol(ranks_matrix)-1
  u_matrix <- vapply(sig_list, FUN.VALUE = numeric(dim), FUN=function(sig){
    sig_neg <- grep('-$', unlist(sig), perl=TRUE, value=TRUE)
    sig_pos <- setdiff(unlist(sig), sig_neg)
    if(length(sig_pos)>0){
      sig_pos <- gsub('\\+$','',sig_pos,perl=TRUE)
      u_p <- as.numeric(ranks_matrix[sig_pos,lapply(.SD, function(x) u_stat(x,maxRank = maxRank,sparse=sparse)), .SDcols=-1, on="rn"])
    }else{
      u_p <- rep(0, dim(ranks_matrix)[2]-1)
    }
    if(length(sig_neg)>0){
      sig_neg <- gsub('-$','',sig_neg,perl=TRUE)
      u_n <- as.numeric(ranks_matrix[sig_neg,lapply(.SD, function(x) u_stat(x,maxRank = maxRank,sparse=sparse)), .SDcols=-1, on="rn"])
    }else{
      u_n <- rep(0, dim(ranks_matrix)[2]-1)
    }
    diff <- u_p - w_neg*u_n # Subtract negative sets, if any
    diff[diff<0] <- 0
    return(diff)
  })
  if(is.vector(u_matrix)){  # Case of ncells=1
    u_matrix <- t(as.matrix(u_matrix))
  }
  rownames(u_matrix) <- colnames(ranks_matrix)[-1]
  return(u_matrix)
}

#' @title Calculate rankings and scores for query data and given signature set
#' @param matrix Input data matrix 
#' @param features List of signatures
#' @param maxRank Rank cutoff (1500) 
#' @param chunk.size Cells per sub-matrix (1000)
#' @param BPPARAM  A BioParallel object to instruct UCell how to parallelize  
#' @param ncores Number of cores to use for parallelization
#' @param w_neg  Weight on negative signatures
#' @param ties.method How to break ties, for data.table::frankv method ("average")
#' @param storeRanks Store ranks? (FALSE) 
#' @param force.gc Force garbage collection? (FALSE) 
#' @param name Suffix for metadata columns ("_UCell") 
#' @return  A list of signature scores
#' @importFrom methods is 
#' @import  Matrix
#' @import  BiocParallel
calculate_Uscore <- function(matrix,features,maxRank=1500,chunk.size=1000,BPPARAM = NULL,ncores=1,w_neg=1, 
                             ties.method="average",storeRanks=FALSE, force.gc=FALSE, name="_UCell"){
  #Make sure we have a sparse matrix
  if(!is(matrix, "dgCMatrix")){matrix <- Matrix::Matrix(as.matrix(matrix),sparse = TRUE)}
  #Check if all genes in signatures are present in the data matrix
  matrix <- check_genes(matrix, features)
  #Do not evaluate more genes than there are
  if(!is.numeric(maxRank)){ stop("Rank cutoff (maxRank) must be a number") }
  if(maxRank > nrow(matrix)){ maxRank <- nrow(matrix) }
  #Weight on neg signatures must be >=0
  if(is.null(w_neg)){w_neg <- 1}
  if(w_neg<0){stop("Weight on negative signatures (w_neg) must be >=0")}
  #Signatures cannot be larger than maxRank parameter
  sign.lgt <- lapply(features, length)
  if(any(sign.lgt > maxRank)){
    stop("One or more signatures contain more genes than maxRank parameter.
            Increase maxRank parameter or make shorter signatures")
  }
  #Split into manageable chunks
  split.data <- split_data.matrix(matrix=matrix, chunk.size=chunk.size)
  #Either take a BPPARAM object, or make one on the spot using 'ncores'
  if(is.null(BPPARAM)){
    BPPARAM <- BiocParallel::MulticoreParam(workers=ncores)
  }
  meta.list <- BiocParallel::bplapply(
    X = split.data, BPPARAM = BPPARAM, FUN = function(x){
      cells_rankings <- data_to_ranks_data_table(x,ties.method = ties.method)
      cells_AUC <- u_stat_signature_list(features,cells_rankings,maxRank = maxRank,sparse = F,w_neg = w_neg)
      colnames(cells_AUC) <- paste0(colnames(cells_AUC),name)
      if(storeRanks==T){
        gene.names <- as.character(as.matrix(cells_rankings[,1]))
        #make sparse
        cells_rankings[cells_rankings>maxRank] <- 0
        ranks.sparse <- Matrix::Matrix(as.matrix(cells_rankings[,-1]),sparse = TRUE)
        dimnames(ranks.sparse)[[1]] <- gene.names
        if(force.gc){
          cells_rankings <- NULL
          gc()
        }
        return(list(cells_rankings=ranks.sparse, cells_AUC=cells_AUC))
      }else{
        if(force.gc){
          cells_rankings <- NULL
          gc()
        }
        return(list(cells_AUC=cells_AUC))
      }
    }
  )
  return(meta.list)
}

#' @title Get signature scores from pre-computed rank matrix
#' @param ranks_matrix  A rank matrix
#' @param features List of signatures
#' @param chunk.size How many cells per matrix chunk
#' @param w_neg Weight on negative signatures
#' @param BPPARAM A BioParallel object to instruct UCell how to parallelize  
#' @param ncores How many cores to use for parallelization?
#' @param force.gc Force garbage collection to recover RAM? (FALSE)
#' @param name  Name suffix for metadata columns ("_UCell")
#' @return A list of signature scores
#' @import data.table
rankings2Uscore <- function(ranks_matrix, features, chunk.size=1000, w_neg=1,
                            BPPARAM = NULL,ncores=1, force.gc=FALSE, name="_UCell"){
  #Check if all genes in signatures are present in the stored signatures
  ranks_matrix <- check_genes(ranks_matrix, features)
  #Weight on neg signatures must be >=0
  if(is.null(w_neg)){w_neg <- 1}
  if(!is.numeric(w_neg) | w_neg < 0){
    stop("Weight on negative signatures (w_neg) must be >=0")}
  maxRank <- max(ranks_matrix)
  split.data <- split_data.matrix(matrix=ranks_matrix, chunk.size=chunk.size)
  rm(ranks_matrix)
  #Either take a BPPARAM object, or make one on the spot using 'ncores'
  if(is.null(BPPARAM)){
    BPPARAM <- BiocParallel::MulticoreParam(workers=ncores)
  }
  meta.list <- BiocParallel::bplapply(
    X = split.data, BPPARAM =  BPPARAM, FUN = function(x){
      dense <- as.matrix(x)
      dense <- data.table::as.data.table(dense, keep.rownames=TRUE)
      data.table::setkey(dense, "rn", physical=FALSE)
      cells_AUC <- u_stat_signature_list(features, dense, maxRank=maxRank, sparse=TRUE, w_neg=w_neg)
      colnames(cells_AUC) <- paste0(colnames(cells_AUC),name)
      if(force.gc){
        dense <- NULL
        gc()
      }
      return(list(cells_AUC=cells_AUC))
    }
  )
  return(meta.list)
}

#' @title Check genes
#' Check if all genes in signatures are found in data matrix - otherwise
#' add zero counts in data-matrix to complete it
#' @param matrix Input data matrix
#' @param features List of genes that must be present (otherwise they are added)
#' @return Same input matrix, extended to comprise any missing genes
check_genes <- function(matrix, features){
  features <- unlist(features)
  features <- gsub("[-+]$","",features,perl=TRUE)
  missing <- setdiff(features, rownames(matrix))
  ll <- length(missing)
  if(ll/length(features) > 0.5){
    n <- round(100*ll/length(features))
    mess <- sprintf("Over half of genes (%s%%)", n)
    mess <- paste(mess, "in specified signatures are missing from data.", "Check the integrity of your dataset.") 
    warning(mess, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
  }
  if(ll>0){
    dim1 <- length(missing)
    dim2 <- ncol(matrix)
    add.mat <-  Matrix::Matrix(data=min(matrix), nrow=dim1, ncol=dim2, sparse = TRUE)
    rownames(add.mat) <- missing
    matrix <- rbind(matrix, add.mat)
    missing.concatenate <- paste(missing, collapse=",")
    mess <- sprintf("The following genes were not found and will be imputed to exp=0:\n* %s",missing.concatenate)
    warning(mess, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
  }
  return(matrix)
}

#' @title Check signature names and add standard names is missing
#' @param features List of signatures for scoring
#' @return The input list of signatures, with standard names if provided un-named
check_signature_names <- function(features){
  defaultSigName <- paste0(rep("signature_",length(features)), seq_along(features))
  if(is.null(names(features))){
    names(features) <- defaultSigName
  }else{
    invalidNames <- names(features) == "" | duplicated(names(features))
    names(features)[invalidNames] <- defaultSigName[invalidNames]
  }
  return(features)
}

#' @title Calculate per-cell feature rankings
#' @param data  Expression data matrix 
#' @param ties.method How to break ties (passed on to data.table::frankv)
#' @return A data.table of ranks 
#' @import data.table
data_to_ranks_data_table <- function(data, ties.method="average"){
  dt <- data.table::as.data.table(as.matrix(data))
  rnaDT.ranks.dt <- dt[, lapply(.SD, function(x) data.table::frankv(x,ties.method=ties.method,order=c(-1L)))]
  rnaDT.ranks.rownames <- rownames(data)
  rnaDT.ranks.dt.rn <- cbind(rn=rnaDT.ranks.rownames, rnaDT.ranks.dt)
  data.table::setkey(rnaDT.ranks.dt.rn, "rn", physical = FALSE)
  return(rnaDT.ranks.dt.rn)
}

#' @title Split data matrix into smaller sub-matrices ('chunks')
#' @param matrix Input data matrix 
#' @param chunk.sizeHow many cells to include in each sub-matrix
#' @return  A list of sub-matrices, each with size {n_features x chunk_size}
split_data.matrix <- function(matrix, chunk.size=1000){
  ncols <- dim(matrix)[2]
  nchunks <- (ncols-1) %/% chunk.size + 1
  split.data <- list()
  min <- 1
  for(i in seq_len(nchunks)){
    if(i == nchunks-1){  #make last two chunks of equal size
      left <- ncols-(i-1)*chunk.size
      max <- min + round(left/2)-1
    }else{
      max <- min(i*chunk.size, ncols)
    }
    split.data[[i]] <- matrix[,min:max,drop=FALSE]
    min <- max+1    #for next chunk
  }
  return(split.data)
}

#' @title Smoothing scores by KNN
#' @param matrix Input data matrix 
#' @param nn A nearest neighbor object returned by [BiocNeighbors::findKNN]
#' @return A dataframe of knn-smoothed scores
knn_smooth_scores <- function(matrix=NULL,nn=NULL){
  sig.cols <- colnames(matrix)
  w.df <- vapply(sig.cols, FUN.VALUE = numeric(nrow(matrix)), FUN = function(s){
    ss.scores <- matrix[,s]
    weighted.scores <- vapply(seq_len(nrow(nn$index)),FUN.VALUE = numeric(1),FUN = function(x){
      r <- nn$index[x,]
      r <- c(x, r)
      d <- nn$distance[x,]
      d <- c(d[1], d)
      w <- 1/(0.001 + d)
      sum(w * ss.scores[r])/sum(w)
    })
  })
  rownames(w.df) <- rownames(matrix)
  as.data.frame(w.df)
}  

#' @rdname SmoothKNN
#' @method SmoothKNN Seurat
#' @export
SmoothKNN.Seurat <- function(obj=NULL,signature.names=NULL,reduction="pca",k=10,BNPARAM=BiocNeighbors::AnnoyParam(),
                             BPPARAM=BiocParallel::SerialParam(),suffix="_kNN",assay=NULL,slot="data", 
                             sce.expname=NULL,sce.assay=NULL){
  if(!reduction %in% Seurat::Reductions(obj)){
    stop(sprintf("Could not find reduction %s in this object", reduction))
  }
  if(is.null(signature.names)){
    stop("Please provide the metadata column names that you want to smooth")
  }
  if(is.null(assay)){  # Work on metadata
    found <- intersect(signature.names, colnames(obj[[]]))
    notfound <- setdiff(signature.names, found)
    if(length(found)==0){
      stop("Could not find any of the given signatures in this object")
    }
    if(length(notfound)>0){
      nf <- paste(notfound, collapse=",")
      mess <- sprintf("The following signature were found in metadata:\n* %s",nf)
      warning(mess, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
    }
    m <- obj[[found]]
  }else{  # Work directly on features
    exp <- Seurat::GetAssayData(obj, assay=assay, slot=slot)
    feats <- rownames(exp)
    found <- intersect(signature.names, feats)
    notfound <- setdiff(signature.names, found)
    if(length(found)==0){
      stop("Could not find any of the given features in this object")
    }
    if(length(notfound)>0){
      nf <- paste(notfound, collapse=",")
      mess <- sprintf("The following features were not found in assay %s:\n* %s",assay, nf)
      warning(mess, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
    }
    m <- t(exp[found, , drop=FALSE])
  }
  ncells <- ncol(obj)
  if(ncells <= k){
    k <- ncells-1
    warning("'k' capped at the number of observations minus 1")
  }
  if(ncells>1){
    # Find kNNs
    space <- Seurat::Embeddings(obj, reduction=reduction)
    nn <- findKNN(space, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
    # Do smoothing
    smooth.df <- knn_smooth_scores(matrix=m, nn=nn)  
  }else{
    smooth.df <- m
  }
  if(is.null(assay)){  #metadata
    colnames(smooth.df) <- paste0(colnames(smooth.df), suffix)
    obj <- Seurat::AddMetaData(obj, metadata = smooth.df)
  }else{  #new assay
    nas <- paste0(assay, suffix)
    obj[[nas]] <- Seurat::CreateAssayObject(data=t(smooth.df))
  }
  return(obj)
}

#' @rdname SmoothKNN
#' @method SmoothKNN SingleCellExperiment
#' @importFrom stats setNames
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay assays SummarizedExperiment assayNames
#' @export
SmoothKNN.SingleCellExperiment <- function(obj=NULL,signature.names=NULL,reduction="PCA",k=10,
                                           BNPARAM=BiocNeighbors::AnnoyParam(),
                                           BPPARAM=BiocParallel::SerialParam(),
                                           suffix="_kNN",assay=NULL,slot="data",
                                           sce.expname=c("UCell","main"), sce.assay=NULL){
  if(!reduction %in% reducedDimNames(obj)){
    stop(sprintf("Could not find reduction %s in this object", reduction))
  }
  if(is.null(signature.names)){
    stop("Please provide the metadata column names that you want to smooth")
  }
  sce.expname <- sce.expname[1]
  if(!sce.expname %in% c(altExpNames(obj), "main")){
    stop(sprintf("Cannot find summarized experiment name: %s", sce.expname))
  } 
  if(sce.expname == "main"){
    exp <- obj
  }else{
    exp <- altExp(obj, sce.expname)
  }
  found <- intersect(signature.names, rownames(exp))
  notfound <- setdiff(signature.names, found)
  if(length(found)==0){
    stop("Could not find any of the given signatures in this object")
  }
  if(length(notfound)>0){
    nf <- paste(notfound, collapse=",")
    mess <- sprintf("The following signatures were not found:\n* %s",nf)
    warning(mess, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
  }
  if(is.null(sce.assay)){
    sce.assay <- 1
  }else if(!sce.assay %in% assayNames(obj)) {
    stop(sprintf("Cannot find assay: %s", sce.assay))
  }
  m <- SummarizedExperiment::assay(exp, sce.assay)
  m <- t(m[found, ,drop=FALSE])
  ncells <- nrow(m)
  if(ncells <= k){
    k <- ncells-1
    warning("'k' capped at the number of observations minus 1")
  }
  if(ncells>1){
    # Find kNNs
    space <- reducedDim(obj, reduction)
    nn <- findKNN(space, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
    # Do smoothing
    m.smooth <- knn_smooth_scores(matrix=m, nn=nn) 
  }else{
    m.smooth <- m
  }
  #New experiment with smoothed scores
  colnames(m.smooth) <- paste0(colnames(m.smooth), suffix)
  sce.newexp <- paste0(sce.expname, suffix)
  l <- list("sce" = t(m.smooth))
  SingleCellExperiment::altExp(obj, sce.newexp) <- SummarizedExperiment(assays = setNames(l, sce.newexp))
  return(obj)
}



#' @title Calculate the expression entropy of each gene
#' @param expr The expression matrix. Rows should be genes and columns should be cells.
#' @param r A small fixed value to avoid log(0) of mean gene expression levels.
#' @return A tibble object with three columns 'Gene', 'mean.expr' and 'entropy'
# ---- ROGUE ----
Entropy <- function(expr, r = 1){
  tmp <- log(expr + 1)
  entropy <- Matrix::rowMeans(tmp)
  mean.expr <- log(Matrix::rowMeans(expr) + r)
  ent_res <- tibble(Gene = rownames(expr), mean.expr = mean.expr, entropy = entropy)
  return(ent_res)
}

#' @title Fit the relationship between expression entropy and mean gene expression
#' @description Fit the relationship between expression entropy and mean gene expression using loess regression
#' @param x A tibble object returned from the `Entropy` function.
#' @param span The parameter α which controls the degree of smoothing.
#' @param mt.method The multiple testing method used in p.adjust.
#' @return A tibble object with six columns
#' @examples 
#' ent.res <- Entropy(expr)
#' ent.res <- entropy_fit(ent.res, span = 0.3, mt.method = 'fdr')
entropy_fit <- function(x, span = 0.5, mt.method = 'fdr'){
  x <- x %>% dplyr::filter(is.finite(mean.expr)) %>% dplyr::filter(entropy > 0)
  fit <- loess(entropy~mean.expr, data = x, span = span)
  prd <- predict(fit, x$mean.expr)
  tmp <- x %>% dplyr::mutate(fit = prd) %>% dplyr::mutate(ds = fit - entropy) %>%
    dplyr::mutate(pv = 1 - pnorm(.$ds, mean = mean(.$ds), sd = sd(.$ds))) %>%
    dplyr::filter(pv > 0.1)
  fit <- loess(entropy~mean.expr, data = tmp, span = span)
  prd <- predict(fit, x$mean.expr)
  tmp <- x %>% dplyr::mutate(fit = prd) %>% dplyr::mutate(ds = fit - entropy) %>%
    dplyr::filter(is.finite(ds)) %>% dplyr::mutate(pv=1-pnorm(.$ds,mean=mean(.$ds),sd=sd(.$ds))) %>%
    dplyr::filter(pv > 0.1)
  fit <- loess(entropy~mean.expr, data = tmp, span = span)
  prd <- predict(fit, x$mean.expr)
  x <- x %>% dplyr::mutate(fit = prd) %>% dplyr::mutate(ds = fit - entropy) %>%
    dplyr::filter(is.finite(ds))
  x <- x %>% dplyr::mutate(p.value = 1-pnorm(x$ds,mean = mean(x$ds),sd = sd(x$ds)))
  p.adj <- p.adjust(x$p.value, method = mt.method)
  x <- x %>% dplyr::mutate(p.adj = p.adj) %>% dplyr::arrange(desc(ds))
  return(x)
}

#' @title ROGUE calculation
#' @description Using ROGUE to assess the purity of single cell population.
#' @param x A tibble object returned from the `Entropy` and `entropy_fit` function.
#' @param platform The platform ('UMI' or 'full-length') used for generating the tested dataset.
#' @param cutoff The threshold (adjusted p value) for identifying significant ds. Default: 0.05.
#' @param k The scaling factor for calculating ROGUE. The default value of k is set to 45 and
#' 500 for droplet-based ('UMI') and 'full-length' based datasets, respectively. When specifying
#' a custom k value, the 'platform' argument is redundant.
#' @param features Use these features to calcuate ROGUE.
#' @details By taking advantage of the wide applicability of S-E model to scRNA-seq data, we introduce the statistic ROGUE to measure the purity of a cell population as:
#' @details \deqn{ROGUE=1-∑sig.ds/(∑sig.ds+K)}
#' @details where K is an important parameter that constrains the ROGUE value between 0 and 1. A cell population with no significant ds for all genes will receive a ROGUE value of 1, while a population with maximum summarization of significant ds is supposed to yield a purity score of ~0.
#' @return A value of ROGUE.
#' calculateRogue(ent.res, platform = 'UMI')
#' calculateRogue(ent.res, k = 30)
calculateRogue <- function(x,platform=NULL,cutoff=0.05,k=NULL,features=NULL){
  if(is.null(k)){
    if(is.null(platform)){
      warning('Cannot find either "platform" or "k" value, calculating k value...')
      ent_res <- Entropy(expr = x, r = 1)
      ent_res <- entropy_fit(ent_res,span = 0.5, mt.method = 'fdr')
      ent_res <- ent_res %>% dplyr::mutate(p.adj = p.value)
      k <- ent_res %>% dplyr::filter(p.adj < 0.05) %>% dplyr::pull(ds) %>% sum()
      k <- k/2
    }else if(platform=='UMI'){
      k <- 45
    }else if(platform=='full-length'){
      k <- 500
    }else if(!is.null(platform) & !(platform %in% c('UMI','full-length'))){
      warning('Please provide valid platform argument: "UMI" or "full-length"...')
    }
  }else{k <- k}
  if(!is.null(features)){
    x <- x %>% dplyr::filter(Gene %in% features)
    sig_value <- sum(abs(x$ds))
    Rogue <- 1 - sig_value/(sig_value + k)
    return(Rogue)
  }else{
    sig_value <- abs(x$ds[x$p.adj < cutoff & x$p.value < cutoff])
    sig_value <- sum(sig_value)
    Rogue <- 1 - sig_value/(sig_value + k)
    return(Rogue)
  }
  
}


#' @title calculate Pathway Activity Score
#' @description Calculation seven tools
#' @param object Seurat object
#' @param gSets_path pathways/gene sets in classical GMT format
#' @param method select PAS tool
#' @param filter whether to filter for gene expressed in less than 5% cells
#' @param normalize normalization method
# ---- PAS ----
cal_all_tools <- function(object, gSets_path,
                          method = c('AUCell','VISION','GSVA','ssGSEA','plage','zscore'),
                          normalize = c('log','CLR','RC','scran','none'),
                          mat_type = 'counts', n_cores = 1, rand_seed = 123){
  # normalization
  if(normalize=='log'){
    object <- Seurat::NormalizeData(object, normalization.method = 'LogNormalize', assay = 'RNA')
  }else if(normalize=='CLR'){
    object <- Seurat::NormalizeData(object, normalization.method = 'CLR', assay = 'RNA')
  }else if(normalize=='RC'){
    object <- Seurat::NormalizeData(object, normalizaiton.method = 'RC', assay = 'RNA')
  }else if(normalize=='scran'){
    if(!require(scran)) BiocManager::install('scran')
    if(!require(scater)) BiocManager::install('scater')
    counts <- GetAssayData(object, assay = 'RNA', slot = 'counts')
    if(mat_type=='counts'){
      sc <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
      if(ncol(counts)>300){
        clusters <- scran::quickCluster(sc, assay.type = 'counts')
        sc <- scran::computeSumFactors(sc, clusters = clusters, assay.type = 'counts')
      }else{
        sc <- scran::computeSumFactors(sc, assay.type = 'counts')
      }
      sc <- scater::normalize(sc, exprs_values = 'counts')
      object@assays$RNA@data <- sc@assays$data$logcounts
    }else{
      sc <- SingleCellExperiment::SingleCellExperiment(assays = list(tpm = counts))
      if(ncol(counts)>300){
        clusters <- scran::quickCluster(sc, assay.type = 'tpm')
        sc <- scran::computeSumFactors(sc, clusters = clusters, assay.type = 'tpm')
      }else{
        sc <- scran::computeSumFactors(sc, assay.type = 'tpm')
      }
      sc <- scater::normalize(sc, exprs_values = 'tpm')
      object@assays$RNA@data <- sc@assays$data$logcounts
    }
  }else if(normalize=='sctransform'){
    object <- Seurat::SCTransform(object)
  }else if(normalize=='scnorm_5'){
    counts <- GetAssayData(object, assay = 'RNA', slot = 'counts')
    counts <- na.omit(counts)
    if(!require(SCnorm)) BiocManager::install('SCnorm')
    DataNorm <- SCnorm::SCnorm(counts[,],rep(c(1),ncol(counts)),K = 5,NCores = n_cores)
    counts <- DataNorm@assays$data$normcounts
    rm(DataNorm);gc()
  }else if(normalize=='scnorm_9'){
    counts <- GetAssayData(object, assay = 'RNA', slot = 'counts')
    counts <- na.omit(counts)
    if(!require(SCnorm)) BiocManager::install('SCnorm')
    DataNorm <- SCnorm::SCnorm(counts[,],rep(c(1),ncol(counts)),K = 9,NCores = n_cores)
    counts <- DataNorm@assays$data$normcounts
    rm(DataNorm);gc()
  }
  cat('normalize successfully!!\n')
  eval_tools <- vector(mode = 'list')
  gSets <- getGMT(gSets_path)
  for(i in 1:length(method)){
    method.tmp <- method[i]
    score <- switch(
      method.tmp,
      AUCell = cal_AUCell(GetAssayData(object,assay='RNA',slot='data'),gSets,n_cores),
      pagoda2 = cal_pagoda2(GetAssayData(object,assay='RNA',slot='data'),gSets,n_cores),
      fscLVM = cal_fscLVM(GetAssayData(object,assay='RNA',slot='data'),gSets_path,type = mat_type),
      VISION = cal_vision(GetAssayData(object,assay='RNA',slot='data'),gSets_path,n_cores),
      ROMA = cal_ROMA(GetAssayData(object,assay='RNA',slot='data'),gSets,n_cores),
      GSVA = GSVA::gsva(GetAssayData(object,assay='RNA',slot='data'),gSets,method='gsva',parallel.sz=n_cores),
      ssGSEA = GSVA::gsva(GetAssayData(object,assay='RNA',slot='data'),gSets,method='ssgsea',parallel.sz=n_cores),
      plage = GSVA::gsva(GetAssayData(object,assay='RNA',slot='data'),gSets,method='plage',parallel.sz=n_cores),
      zscore = GSVA::gsva(GetAssayData(object,assay='RNA',slot='data'),gSets,method='zscore',parallel.sz=n_cores),
      gene = counts
    )
    eval_tools[[i]] <- score
  }
  names(eval_tools) <- method
  return(eval_tools)
}

#' @title cal_AUCell
cal_AUCell <- function(counts,gSets,n_cores){
  if(!require(AUCell)) BiocManager::install('AUCell')
  ac_rankings <- AUCell::AUCell_buildRankings(counts,nCores = n_cores,plotStats = F,splitByBlocks = T)
  sc_AUC <- AUCell::AUCell_calcAUC(gSets,ac_rankings,normAUC = T,aucMaxRank = ceiling(0.05*nrow(ac_rankings)))
  score <- AUCell::getAUC(sc_AUC)
  gc()
  return(score)
}

#' @title cal_pagoda2
cal_pagoda2 <- function(counts,gSets,n_cores){
  if(!require(pagoda2)) install.packages('pagoda2')
  library(pagoda2)
  nPcs <- min(round(ncol(counts)/5),5)
  p2 <- Pagoda2$new(counts, n.cores = n_cores, log.scale = T)
  p2$adjustVariance(plot = F)
  p2$calculatePcaReduction(nPcs = nPcs, use.odgenes = F, fastpath = F)
  path_names <- c()
  env <- new.env(parent = globalenv())
  invisible(lapply(1:length(gSets), function(i){
    genes <- intersect(gSets[[i]], rownames(counts))
    name <- paste0(names(gSets[i]), i)
    if(length(genes)>3){
      assign(name, genes, envir = env)
      path_names <- c(path_names, name)
    }
  }))
  p2$testPathwayOverdispersion(setenv = env, verbose = T, recalculate.pca = T, min.pathway.size = 1)
  path_names <- names(p2$misc$pwpca)
  score <- matrix(NA, nrow = length(path_names), ncol = ncol(counts))
  rownames(score) <- path_names
  colnames(score) <- colnames(counts)
  for(i in 1:length(p2$misc$pwpca)){
    if(!is.null(p2$misc$pwpca[[i]]$xp$scores)){
      score[i,] <- as.numeric(p2$misc$pwpca[[i]]$xp$scores)
    }
  }
  gc()
  return(score)
}

#' @title cal_fscLVM
cal_fscLVM <- function(counts,gSets_path,n_cores,type='counts',
                       n_hidden=1,ifLog = T, rand_seed = 123){
  n_hidden <- ifelse(ncol(counts)>5000, ifelse(ncol(counts)>10000, 3, 2), 1)
  cat('n_hidden:', n_hidden)
  cat('type:',type)
  if(type=='counts'){
    if(ifLog){
      sc_mat <- SingleCellExperiment::SingleCellExperiment(assays=list(logcounts=log2(counts[,]+1)))
      slot <- 'logcounts'
    }else{
      sc_mat <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts[,]))
      slot <- 'counts'
    }
  }else if(type=='tpm'){
    if(ifLog){
      sc_mat <- SingleCellExperiment::SingleCellExperiment(assays=list(logtpm=log2(counts[,]+1)))
      slot <- 'logtpm'
    }else{
      sc_mat <- SingleCellExperiment::SingleCellExperiment(assays=list(tpm=counts[,]))
      slot <- 'tpm'
    }
  }
  print('slot:',slot)
  
  genesets <- GSEABase::getGmt(gSets_path)
  if(!require(slalom)) BiocManager::install('slalom')
  model <- slalom::newSlalomModel(
    sc_mat,genesets,n_hidden = n_hidden,min_genes = 1,assay_name = slot,
    prune_genes = T, design = NULL, anno_fpr = 0.01, anno_fnr = 0.01
  )
  cat('model construction ...\n')
  model <- slalom::initSlalom(model,noise_model = 'gauss',alpha_priors = NULL,epsilon_priors = NULL,
                              design = NULL, pi_prior = NULL, n_hidden = NULL, seed = rand_seed)
  cat('initial sucessfull ...\n')
  model <- slalom::trainSlalom(model,nIterations = 1000,minIterations = 300,shuffle = T,seed = rand_seed,
                               tolerance = 1e-05,forceIterations = F,pretrain = T,drop_factors = T)
  if(model@.xData$converged){
    score <- model$X_E1
    colnames(score) <- model$termNames
    rownames(score) <- colnames(counts)
    score <- score[, (n_hidden + 1):ncol(score)]
    score <- t(score)
    return(score)
  }else{ return('not converged') }
}

#' @title cal_vision
cal_vision <- function(counts,gSets_path,n_cores){
  if(!require(VISION)) devtools::install_github("YosefLab/VISION")
  vis <- VISION::Vision(counts, signatures = gSets_path, projection_methods = 'UMAP', sig_gene_threshold = 0)
  options(mc.cores = n_cores)
  vis <- VISION::analyze(vis)
  score <- t(vis@SigScores)
  return(score)
}

#' @title cal_ROMA
cal_ROMA <- function(counts,gSets,n_cores){
  if(!require(rROMA)) devtools::install_github("Albluca/rROMA")
  roma <- rRoma::rRoma.R(
    ExpressionMatrix = data.matrix(counts), ModuleList = gSets, ClusType = 'FORK',
    MinGenes = 1, UseParallel = 1, GeneOutThr = 1, OutGeneNumber = 1, nCores = n_cores,
    centerData = F, ExpFilter = F, UseWeigths = F, DefaultWeight = 1, OutGeneSpace = NULL,
    ApproxSamples = 5, nSamples = 100, Ncomp = 10, FixedCenter = T, GeneOutDetection = 'L1OutExpOut',
    GeneSelMode = 'All', SampleFilter = F, MoreInfo = F, PlotData = F, PCADims = 2,
    PCSignMode = 'none', PCSignThr = NULL, SamplingGeneWeights = NULL, Grouping = NULL,
    FullSampleInfo = F, GroupPCSign = F, CorMethod = 'pearson', PCAType = 'DimensionsAreGenes'
  )
  score <- roma$SampleMatrx
  return(score)
}

#' @title functions to convert data types
#' @description Reading a .gmt file containing pathways or gene sets.
#' @param gmt_file .gmt file
#' @param n_gene_thre minimum number of genes in gene sets
#' @examples
#' gsets <- getGMT('/path/kegg.gmt')
getGMT <- function(gmt_file, n_gene_thre = 0){
  paths <- readLines(gmt_file)
  gsets <- list()
  for(i in 1:length(paths)){
    t <- strsplit(paths[i], '\t')[[1]]
    genes <- t[3:length(t)]
    genes <- genes[which(genes != '')]
    genes <- unique(genes)
    if(length(genes) > n_gene_thre) gsets[[t[1]]] <- genes
  }
  return(gsets)
}


#' @title Compare pathways within a Seurat object.
#' @description This function takes a Seurat object as an input, and compare gene sets over specified
#' conditions/populations.
#' @param seu Seurat object with populations defined in the metadata.
#' @param group1 First comparison group as defined by metadata in Seurat object e.g. cell_type.
#' @param group1.population Populations within group1 to compare Seurat object e.g. hour.
#' @param group2 Second comparison group as defined by column names in Seurat object e.g. hour.
#' @param group2.population Population within group2 to compare e.g. 24.
#' @param pathways Pathway gene sets with each pathway in a separate list.
#' @param downsample Option to downsample cell numbers. Defaults to 500 cells per condition. If a
#' population has < 500 cells, all cells from that condition are used.
#' 
#' @examples 
#' pathways <- msigdbr::msigdbr(species = 'Homo sapiens', category = 'H') %>% format_pathways()
#' out <- compare_seurat(seu = seu, group1 = 'cell', group1_population = c("t_cell","b_cell"),
#' group2 = 'hour', group2_population = c('24'), pathways = pathways)
#' 
#' # Formatting pathway manually
#' gsets <- getGMT('reference/gmtFiles/mouse/reactome.gmt')
#' pathway <- list()
#' for(i in 1:length(gsets)){
#'   pathway[[i]] <- data.frame(Pathway = rep(names(gsets[i]),length(gsets[i])), Genes = gsets[[i]])
#' }
#' 
#' @return Statistical results. The qval should be the primary metric that is used to interpret pathway
#' differences i.e. a higher qval translates to larger pathway differences between conditions.
#' If only two samples are provided, a fold change (FC) enrichment score will also be calculated.
#' The FC output is generated from a running sum of mean changes in gene expression from all genes
#' of the pathway. It's calculated from average pathway expression in population1 - population2,
#' so a negative FC means the pathway is higher in population2.
# ---- Pseudotime ----
compare_seurat <- function(seu,group1=NULL,group1.population=NULL,group2=NULL,
                           group2.population=NULL,pathways,downsample=500){
  path.names <- sapply(pathways, function(x) unique(x$Pathway))
  if(is.null(group2)){
    samples <- list()
    for(i in group1.population){
      samples[[i]] <- seurat_extract(seu, meta1 = group1, value_meta1 = i)
    }
  }
  if(!is.null(group2)){
    samples <- list()
    for(i in group1.population){
      samples[[i]] <- seurat_extract(seu,meta1 = group1,value_meta1 = i,meta2 = group2,value_meta2 = group2.population)
    }
  }
  mcm_output <- compare_pathways(samples = samples, pathways = pathways)
}

#' @title Compare pathways within a SingleCellExperiment object.
#' @description This function takes a SingleCellExperiment object as an input, and compares gene sets
#' over specified conditions/populations.
#' @param sce SingleCellExperiment object with populations defined in the column data
#' @param assay.name Assay name to extract expression values from. Default to logcounts.
#' @param group1 First comparison group as defined by `colData()` columns of SingleCellExperiment object
#' e.g. cell_type.
#' @param group1.population Populations within group1 to compare e.g. c('t_cell', 'b_cell')
#' @param group2 Second comparison group as defined by `colData()` columns of SingleCellExperiment object
#' .g. hour
#' @param group2.population Popualtion within group2 to compare e.g. 24.
#' @param pathways Pathway gene sets with each pathway in a separate list.
#' @param downsample Option to downsample cell numbers. Defaults to 500 cells per condition. If a population
#' has < 500 cells, all cells from that condition are used.
#' 
#' @examples 
#' pathways <- msigdbr::msigdbr(species = 'Homo sapiens', category = 'H') %>% format_pathways()
#' out <- compare_sce(sce, group1='cell',group1_population = c("t_cell", "b_cell"),
#' group2 = "hour", group2_population = c("24"), pathways = pathways)
#' 
#' # Formatting pathway manually
#' gsets <- getGMT('reference/gmtFiles/mouse/reactome.gmt')
#' pathway <- list()
#' for(i in 1:length(gsets)){
#'   pathway[[i]] <- data.frame(Pathway = rep(names(gsets[i]),length(gsets[i])), Genes = gsets[[i]])
#' }
compare_sce <- function(sce,assay.name='logcounts',group1=NULL,group1.population=NULL,
                        group2=NULL,group2.population=NULL,pathways,downsample=500){
  path.names <- sapply(pathways, function(x) unique(x$Pathway))
  samples <- list()
  if(is.null(group2)){
    for(i in group1.population){
      samples[[i]] <- sce_extract(sce, assay.name = assay.name, meta1 = group1, value_meta1 = i)
    }
  }else{
    for(i in group1.population){
      samples[[i]] <- sce_extract(sce,assay.name = assay.name,meta1 = group1,value_meta1 = i,
                                  meta2 = group2,value_meta2 = group2.population)
    }
  }
  mcm_output <- compare_pathways(samples = samples, pathways = pathways)
  return(mcm_output)
}

#' @title Plot a heatmap of qvals
#' @description This function takes the output from `compare_pathway` and plots a heatmap of the qvals.
#' @param out Data frame that contains a "Pathway" column, and a column or multiple columns containing
#' qvals. This can be the direct output from `compare_pathways()`, or a custom data frame.
#' @param highlight.pathways Pathway or pathway to annotate on the heatmap, supplied as character vector.
#' If no argument is given, all pathways are shown.
#' @param row.fontsize Font size of pathway names.
#' @param column.fontsize Font size of the sample names.
#' @param column.names Option to supply names of the heatmap, column if more than one populations is present.
#' @param show.row.names Should row names be shown in the heatmap?
#' @param cluster.columns Should columns in the heatmap be clustered?
#' @param hm.colors Colors to be used in the heatmap?
#' @param scale.breaks Breaks to be used in the colors of the heatmap. Length must be equal to the number of colors.
#' @examples 
#' generate_heatmap(out = compare_pathway_results, highlight.pathways = 'mtorc')
generate_heatmap <- function(out,highlight.pathways=NULL,row.fontsize=6,column.fontsize=10,
                             column.names=colnames(out), show.row.names=T,cluster.columns=T,
                             hm.colors=NULL,scale.breaks=NULL){
  if(!is.null(hm.colors)){
    heatmap.colors <- hm.colors
  }else{
    heatmap.colors <- c('cornflowerblue', 'white', 'red')
  }
  scale <- out %>% dplyr::select(grep('qval',colnames(out),ignore.case = T,value = T)) %>% as.matrix()
  if(is.null(scale.breaks)){
    scale.breaks <- c(min(scale), mean(scale), max(scale))
  }else{
    scale.breaks <- scale.breaks
  }
  hm.col <- circlize::colorRamp2(colors = heatmap.colors, breaks = scale.breaks)
  if(is.null(highlight.pathways)==F){
    pathways <- out %>% dplyr::filter(grepl(paste(highlight.pathways,collapse = '|'),Pathway,ignore.case = T)) %>% dplyr::pull(Pathway)
    pathways <- gsub(pattern = '_', replacement = ' ', x = pathways)
    out <- out %>% dplyr::mutate(Pathway = gsub(pattern = '_', replacement = ' ',x = Pathway))
    position <- out %>% dplyr::mutate(position = 1:nrow(.)) %>%
      dplyr::filter(grepl(paste(pathways, collapse = '|'),Pathway, ignore.case = T)) %>%
      dplyr::pull(position)
    out <- out %>% tibble::remove_rownames() %>% tibble::column_to_rownames('Pathway') %>%
      dplyr::select(grep(pattern = 'qval',x = colnames(out),ignore.case = T,value = T))
    row.an <- ComplexHeatmap::rowAnnotation(
      Genes = ComplexHeatmap::anno_mark(at = which(rownames(out) %in% pathways),
                                        labels = rownames(out)[position],
                                        labels_gp = grid::gpar(fontsize = 7),
                                        link_width = grid::unit(2.5, 'mm'),
                                        padding = grid::unit(1, 'mm'),
                                        link_gp = grid::gpar(lwd = 0.5))
    )
    ComplexHeatmap::Heatmap(
      out, name = 'Qval', border = T, rect_gp = grid::gpar(col = 'white', lwd = 0.1),
      row_names_gp = grid::gpar(fontsize = row.fontsize),
      column_names_gp = grid::gpar(fontsize = column.fontsize),
      right_annotation = row.an, column_labels = column.names, col = hm.col,
      show_row_names = show.row.names, row_dend_width = grid::unit(6,'mm'),
      cluster_columns = cluster.columns
    )
  }else{
    out <- out %>% dplyr::mutate(Pathway = gsub(pattern = '_',replacement = ' ',x = Pathway)) %>%
      tibble::remove_rownames() %>% tibble::column_to_rownames('Pathway') %>%
      dplyr::select(grep(pattern = 'qval',x = colnames(out), ignore.case = T, value = T))
    ComplexHeatmap::Heatmap(
      out, name = 'Qval', border = T, rect_gp = grid::gpar(col = 'white', lwd = 0.1),
      row_names_gp = grid::gpar(fontsize = row.fontsize),
      column_names_gp = grid::gpar(fontsize = column.fontsize),
      show_row_names = show.row.names, column_labels = column.names, col = hm.col,
      row_dend_width = grid::unit(6,'mm'), cluster_columns = cluster.columns
    )
  }
}

#' @title Plot the rank of specific pathways
#' @description This function takes the output of `compare_pathway()` and plot the rank of a user
#' defined pathways.
#' @param out Data frame containing Pathways and qvals generated from `compare_pathways()`.
#' @param pathway Chosen pathway or pathways to plot the rank of. This can be specific e.g.
#' HALLMARK_GLYCOLYSIS to plot a specific result, or generic e.g. glycolysis to plot all
#' glycolysis pathways.
#' @param population.name Column name of the population to plot, if more than one population present
#' in the data frame.
#' @param base.point.size Size of the base points to plot on the graph.
#' @param base.point.color Color of the base points to plot on the graph
#' @param highlight.point.size Size of the highlighted points to plot on the graph.
#' @param highlight.point.color Color of the highlighted points to plot on the graph.
#' @param label.pathway Should the selected pathway be labelled?
#' @param label.size Text size of the pathway label.
#' @examples 
#' 
#' plot_rank(out = compare_pathway_result, pathway = 'interferon', population.name = 'qval')
plot_rank <- function(out,pathway,population.name='qval',base.point.size=2,base.point.color='grey70',
                      highlight.point.size=3,highlight.point.color='cornflowerblue',
                      label.pathway = T, label.size = 4){
  library(ggplot2)
  selected.pathways <- grep(paste(pathway,collapse = '|'),out$Pathway, ignore.case = T,value = T)
  pathway.ranking <- dplyr::arrange(out, dplyr::desc(population.name))
  pathway.ranking$path_rank <- dplyr::percent_rank(pathway.ranking[[population.name]])*100
  df.sub <- subset(pathway.ranking, pathway.ranking$Pathway %in% selected.pathways)
  if(label.pathway==T){
    pathway.label <- grep(paste(pathway,collapse = '|'),pathway.ranking$Pathway, value = T,ignore.case = T)
    pathway.label <- pathway.ranking[pathway.ranking$Pathway %in% pathway.label,]
    for(i in 1:length(pathway.label$Pathway)){
      pathway.label$Pathway[i] <- gsub(pathway.label$Pathway[i], pattern = '_', replacement = ' ')
    }
    pathway.label$Pathway <- stringr::str_to_title(pathway.label$Pathway)
    ggplot(pathway.ranking, aes(.data[[population.name]],path_rank)) +
      geom_hline(yintercept = c(0, 25, 50, 75, 100), linetype = 'dotted', lwd = 0.3, color = 'grey40') +
      geom_point(shape = 21,cex = base.point.size,color = 'black',fill = base.point.color,stroke = 0.05) +
      ggrepel::geom_label_repel(data = pathway.label, label = pathway.label$Pathway, size = label.size,
                                label.padding = unit(0.7,'mm'), label.r = unit(0.3,'mm'),nudge_x = -30) +
      geom_point(data = df.sub,shape = 21,cex = highlight.point.size,color = 'black',fill = highlight.point.color) +
      xlab('Qval') + ylab('Pathway rank') + scale_y_continuous(expand = c(0.03,0.03),breaks = c(0,25,50,75,100)) +
      scale_x_continuous(expand = c(0.2, 0.2)) +
      theme(panel.border = element_rect(fill = NA), panel.background = element_blank(),
            title = element_text(size = 9), axis.title = element_text(size = 11))
  }else{
    ggplot(pathway.ranking,aes(.data[[population.name]],path_rank)) +
      geom_hline(yintercept = c(0,25,50,75,100), linetype = 'dotted', lwd = 0.03, color='gray40') +
      geom_point(shape = 21, cex = base.point.size, color = 'black', fill = base.point.color, stroke = 0.05) +
      geom_point(data = df.sub,shape = 21,cex = highlight.point.size,color = 'black',fill = highlight.point.color) +
      xlab('Qval') + ylab('Pathway rank') + scale_y_continuous(expand = c(0.03,0.03),breaks = c(0,25,50,75,100)) +
      scale_x_continuous(expand = c(0.2,0.2)) +
      theme(panel.border = element_rect(fill = NA), panel.background = element_blank(),
            title = element_text(size = 9), axis.title = element_text(size = 11))
  }
}

#' @title Compare gene sets
#' @description This function takes an input of samples and pathways to compare gene set
#' perturbations over different conditions.
#' @param samples List of samples, each supplied as an expression matrix with cells in columns 
#' and genes in rows.
#' @param pathways Pathways and their genes with each pathway in a separate list.
#' @param downsample Option to downsample cell numbers. Defaults to 500 cells per condition. If a
#' population has < 500 cells, all cells from that condition are used.
#' @param min_genes Gene sets with fewer than this number of genes will be excluded.
#' @param max_genes Gene sets with more than this number of genes will excluded.
#' 
#' @examples 
#' # Case 1.
#' pathways <- msigdbr::msigdbr(species = 'Homo sapiens', category = 'H') %>% format_pathways()
#' result <- compare_pathways(list(sample1, sample2, sample3), pathways = pathways)
#' # Formatting pathway manually
#' gsets <- getGMT('reference/gmtFiles/mouse/reactome.gmt')
#' pathway <- list()
#' for(i in 1:length(gsets)){
#'   pathway[[i]] <- data.frame(Pathway = rep(names(gsets[i]),length(gsets[i])), Genes = gsets[[i]])
#' }
#' 
#' # Case 2. Compare pathway in two Seurat object
#' load('data/query_example_seurat.RData')
#' ref <- load.reference.map()
#' projected <- Run.ProjecTME(query_example_seurat, ref = ref, fast.umap.predict = T)
#' CD8_Tex <- seurat_extract(projected, meta1 = 'functional.cluster',value_meta1 = 'CD8_Tex')
#' CD8_Tpex <- seurat_extract(projected,meta1 = 'functional.cluster',value_meta1 = 'CD8_Tpex')
#' gsets <- getGMT('reference/gmtFiles/mouse/reactome.gmt')
#' pathways <- list()
#' for(i in 1:length(gsets)){
#'   pathways[[i]] <- data.frame(Pathway = rep(names(gsets[i]),length(gsets[i])), Genes = gsets[[i]])
#' }
#' Tex_Tpex <- compare_pathways(samples = list(CD8_Tex,CD8_Tpex),pathways = pathways)
#' # Visualizing results
#' Tex_Tpex <- Tex_Tpex %>%
#' mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
#'                          FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
#'                          FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
#'                          FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))
#' aa.path <- Tex_Tpex %>% filter(grepl('arachi',ignore.case = T,x = Pathway))
#' ggplot(Tex_Tpex, aes(-FC, qval)) +
#' geom_vline(xintercept = c(-2, 2), linetype = 'dashed', col = 'black', lwd = 0.3) +
#'   geom_point(cex = 4, shape = 21, fill = Tex_Tpex$color, stroke = 0.3) +
#'   geom_point(data = aa.path, shape = 21, cex = 5, fill = 'orangered2', color = 'black', stroke = 0.3) +
#'   xlim(-10, 10) + ylim(0, 4) + xlab('Enrichment') + ylab('Qval') +
#'   theme(panel.background = element_blank(),panel.border = element_rect(fill = NA),aspect.ratio = 1)
#' 
#' @return Statistical results from the analysis. The qval should be the primary metric that is used to
#' interpret pathway differences i.e. a higher qval translates to larger pathway differences between
#' conditions.
#' If only two samples are provided, a fold change (FC) enrichment score will also be calculated. The
#' FC statistic is generated from a running sum of mean changes in gene expression from all genes of
#' the pathway. It's calculated from average pathway expression in population1 - population2, so a 
#' negative FC means the pathway is higher in population2.
compare_pathways <- function(samples, pathways, downsample = 500, min_genes = 15, max_genes = 500){
  path.names <- sapply(pathways, function(x) unique(x$Pathway))
  # define the number of cells in each condition
  cell.number <- lapply(samples, function(x) ncol(x))
  cell.number <- sapply(cell.number, function(x) x[1])
  for(i in 1:length(cell.number)){
    message(paste('Cell numbers in population', i, '=', cell.number[i]))
  }
  message('- If greater than ', downsample, ' cells, these populations will be downsampled', '\n')
  # randomly sample cells
  for(i in 1:length(samples)){
    samples[[i]] <- samples[[i]][,sample(ncol(samples[[i]]),ifelse(cell.number[i]<500,cell.number[i],downsample))]
  }
  # Only take shared genes
  genes <- lapply(samples, function(x) rownames(x))
  genes <- table(unlist(genes))
  genes <- genes[genes == length(samples)]
  genes <- names(genes)
  samples <- lapply(samples, function(x) x[rownames(x) %in% genes,])
  # generate pathway matrices
  pop.paths <- vector(mode = 'list', length = length(samples))
  for(i in 1:length(pop.paths)){
    for(j in 1:length(pathways)){
      pop.paths[[i]][[j]] <- samples[[i]][rownames(samples[[i]]) %in% pathways[[j]]$Genes,]
    }
  }
  # set pathway names
  pop.paths <- lapply(pop.paths, function(x) purrr::set_names(x, path.names))
  # filter out pathways with < 15 genes or > 500
  filter.paths <- sapply(pop.paths[[1]], function(x) any(nrow(x) >= min_genes & nrow(x) <= max_genes))
  filtered.pathways <- names(filter.paths[filter.paths == 'FALSE'])
  if(length(filtered.pathways) > 0){
    message('Excluding ',length(filtered.pathways),' pathways based on min/max genes parameter: ',
            paste(head(filtered.pathways, 5),collapse = ', '), '...', '\n')
  }else{
    message('All ', length(pathways), ' pathways passed the min/max genes threshold','\n')
  }
  pop.paths <- lapply(pop.paths, function(x) x[unlist(filter.paths)])
  # transpose matrix
  pop.paths <- lapply(pop.paths, function(x) lapply(x, function(j) t(j)))
  # order columns
  pop.paths <- lapply(pop.paths, function(x) lapply(x, function(j) j[,sort(colnames(j))]))
  # FC calculation
  if(length(samples) == 2){
    message('Calculating pathway fold change...', '\n')
    avg.exp <- lapply(pop.paths, function(x) lapply(x, function(j) data.frame(colMeans(j))))
    samp.combined <- c()
    for(i in 1:length(pop.paths[[1]])){
      samp.combined[[i]] <- cbind(avg.exp[[1]][[i]], avg.exp[[2]][[i]])
    }
    samp.combined <- lapply(samp.combined, function(x) magrittr::set_colnames(x,c('Pop1','Pop2')))
    samp.combined <- lapply(samp.combined,function(x) cbind(x,logFC = x[,'Pop1'] - x[,'Pop2']))
    path.fc <- sapply(samp.combined, function(x) sum(x[,'logFC']))
  }
  if(length(samples) > 2){
    message('Performing a multisample analysis...')
    pb <- txtProgressBar(min = 0, max = length(pop.paths[[1]]), style = 3, width = 50)
    mcm_results <- list()
    for(i in 1:length(pop.paths[[1]])){
      Sys.sleep(0.1)
      setTxtProgressBar(pb, i)
      mcm_results[[i]] <- mcm(lapply(pop.paths,function(x) x[[i]]), level = 0.05)
    }
    close(pb)
    mcm_output <- data.frame(t(sapply(mcm_results, c)), stringsAsFactors = F)
    mcm_output$Pathway <- names(pop.paths[[1]])
    mcm_output$X2 <- NULL
    colnames(mcm_output)[1] <- c('Pval')
    mcm_output$Pval <- replace(mcm_output$Pval, mcm_output$Pval == 0, 10^-300)
    mcm_output$Pval <- as.numeric(mcm_output$Pval)
    mcm_output$adjPval <- p.adjust(mcm_output$Pval,method = 'bonferroni',n = nrow(mcm_output))
    mcm_output$qval <- sqrt(-log10(mcm_output$adjPval))
    mcm_output <- mcm_output[, c(2,1,3,4)]
    mcm_output <- mcm_output[order(-mcm_output$qval),]
  }else{
    message('Performing a two-sample analysis...')
    pb <- txtProgressBar(min = 0, max = length(pop.paths[[1]]), style = 3, width = 50)
    mcm_result <- list()
    for(i in 1:length(pop.paths[[1]])){
      Sys.sleep(0.1)
      setTxtProgressBar(pb, i)
      mcm_result[[i]] <- mcm(lapply(pop.paths,function(x) x[[i]]),level = 0.05)
    }
    close(pb)
    mcm_output <- data.frame(t(sapply(mcm_result,c)),stringsAsFactors = F)
    mcm_output$FC <- path.fc
    mcm_output$Pathway <- names(pop.paths[[1]])
    mcm_output$X2 <- NULL
    colnames(mcm_output)[1] <- c('Pval')
    mcm_output[mcm_output$Pval == 0] <- 10^-300
    mcm_output$Pval <- as.numeric(mcm_output$Pval)
    mcm_output$adjPval <- p.adjust(mcm_output$Pval, method = 'bonferroni',n = nrow(mcm_output))
    mcm_output$qval <- sqrt(-log10(mcm_output$adjPval))
    mcm_output <- mcm_output[, c(3,1,4,5,2)]
    mcm_output <- mcm_output[order(-mcm_output$qval),]
  }
  return(mcm_output)
}

#' @title Extract Normalized Data From Seurat
#' @description This function takes a Seurat object as an input, and returns an expression matrix based on
#' subsetting parameters. Either none, one, or two metadata features can be selected for a given input.
#' @param seu Seurat object containing normalized counts stored in `seu@assays$RNA@data`.
#' @param meta1 Metadata column to subset.
#' @param meta2 Metadata column to subset.
#' @param value_meta1 Value to select within `meta1` column.
#' @param value_meta2 Value to select within `meta2` column.
#' @param pseudocount Pseudocount to add to data. Defaults to 0.001.
#' @examples 
#' cd4 <- seurat_extract(seu, meta1 = 'Hour', value_meta1 = 12, meta2 = 'cell_type', value_meta2 = 'CD4')
seurat_extract <- function(seu, meta1 = NULL, value_meta1 = NULL, meta2 = NULL, 
                           value_meta2 = NULL, pseudocount = 0.001){
  if(is.null(meta1) && is.null(meta2)){
    message('No metadata selected. Converting whole Seurat object to matrix')
    seu <- as.matrix(seu@assays$RNA@data) + pseudocount
  }
  if(!is.null(meta1) && is.null(meta2)){
    message(paste0('Extracting cells where ', meta1, ' == ', value_meta1))
    met_1 <- Seurat::FetchData(object = seu, vars = meta1)
    seu <- seu[, which(met_1 == value_meta1)]
    seu <- as.matrix(seu@assays$RNA@data) + pseudocount
  }
  if(!is.null(meta1) && !is.null(meta2)){
    message(paste0('Extracting cells where ', meta1, ' == ', value_meta1, ' AND ', meta2, ' == ',value_meta2))
    met_1 <- Seurat::FetchData(object = seu, vars = meta1)
    meta_2 <- Seurat::FetchData(object = seu, vars = meta2)
    seu <- seu[, which(met_1 == value_meta1 & met_2 == value_meta2)]
    seu <- as.matrix(seu@assays$RNA@data) + pseudocount
  }
  return(seu)
}

#' @title Extract data from a SingleCellExperiment object
#' @description This function takes a SingleCellExperiment object as an input, and returns an expression matrix
#' based on subsetting parameters. Either none, one, or two metadata features can be selected for a given input.
#' @param sce SingleCellExperiment object containing expression data
#' @param assay.name Name of assay to pull from. Defaults to 'logcounts'.
#' @param meta1 Metadata column to subset.
#' @param value_meta1 Value to select within `meta1` column.
#' @param meta2 Metadata column to subset.
#' @param value_meta2 Value to select within `meta2` column.
#' @param pseudocount Pseudocount to add to data. Defaults to 0.001
#' @examples 
#' cd4 <- sce_extract(sce, meta1 = 'Hour', meta2 = 'cell_type', value_meta2 = 'CD4')
sce_extract <- function(sce, assay.name = 'logcounts', meta1 = NULL, value_meta1 = NULL,
                        meta2 = NULL, value_meta2 = NULL, pseudocount = 0.001){
  if(is.null(meta1) && is.null(meta2)){
    message('No metadata selected. Converting whole SCE object to matrix')
    sub.sce <- SummarizedExperiment::assay(sce, assay.name)
    sub.sce <- as.matrix(sub.sce) + pseudocount
  }
  if(!is.null(meta1) && is.null(meta2)){
    message(paste0('Extracting cells where ', meta1, ' == ', value_meta1))
    sub.sce <- sce[, sce[[meta1]] == value_meta1]
    sub.sce <- SummarizedExperiment::assay(sub.sce, assay.name)
    sub.sce <- as.matrix(sub.sce) + pseudocount
  }
  if(!is.null(meta1) && !is.null(meta2)){
    message(paste0('Extracting cells where ', meta1, ' == ',value_meta1, ' AND ', meta2, ' == ',value_meta2))
    sub.sce <- sce[,sce[[meta1]] == value_meta1 & sce[[meta2]] == value_meta2]
    sub.sce <- SummarizedExperiment::assay(sub.sce, assay.name)
    sub.sce <- as.matrix(sub.sce) + pseudocount
  }
  return(sub.sce)
}

#' @title Format the output of msigdbr.
#' @description This function takes the output of msigdbr and formats the pathways into something that
#' can be used.
#' @param msigdbr.output output from msigdbr
#' @examples 
#' library(magrittr)
#' library(dplyr)
#' pathways <- msigdbr::msigdbr(species = 'Homo sapiens', category = 'H') %>% format_pathways()
#' 
#' # Versatility in using msigdbr
#' ifn_pathways <- msigdbr::msigdbr('Homo sapiens') %>% filter(grepl('interferon',gs_name,ignore.case = T)) %>% format_pathways()
#' 
#' # or we could specify a combination of collections so we get all of the Hallmark,KEGG,and Reactome genesets.
#' pathways <- c('hallmark', 'kegg', 'reactome')
#' hkr_sets <- msigdbr::msigdbr('Homo sapiens') %>% 
#'   filter(grepl(paste(pathways,collapse = '|'),gs_name,ignore.case = T)) %>%
#'   format_pathways()
format_pathways <- function(msigdbr.output){
  msigdbr.output <- msigdbr.output[, c('gs_name', 'gene_symbol')]
  colnames(msigdbr.output) <- c('Pathway', 'Genes')
  msigdbr.output <- dplyr::group_split(msigdbr.output, Pathway)
  return(msigdbr.output)
}

#' @title Multisample generalization of Rosenbaum's crossmatch test
#' @description In this script, we present a framework inspired by Rosenbaum's crossmatch idea to tackle the
#' nonparametric, multisample problem wherein one is concerned with testing the equality of K unknown
#' multivariate probability distributions.
#' We implement two tests: the first is a multisample generalization of Rosenbaum's crossmatch (MCM), and the
#' other further introduces a Malahobis-type modification to the test (MMCM).
#' @param deta_list is list of multifeature matrices corresponding to the K different classes, so each element
#' of the list is a matrix, for a total of K matrices. Each matrix contains observations as the rows and 
#' features as the columns.
#' @param level is the level alpha for hypothesis testing.
#' @return The p-value corresponding to rejection of the alternative, along with the decision of the hypothesis
#' testing (Null being accepted versus rejected).
#' @examples 
#' # Simulation Example when the user wants to test whether K = 3 multivariate distribution are equal:
#' X1 = MASS::mvrnorm(10, rep(0,4), diag(2, 4), tol = 1e-6, empirical = F, EISPACK = F)
#' X2 = MASS::mvrnorm(10, rep(0,4), diag(1, 4), tol = 1e-6, empirical = F, EISPACK = F)
#' X3 = MASS::mvrnorm(10, rep(0,4), diag(3, 4), tol = 1e-6, empirical = F, EISPACK = F)
#' mcm(list(X1, X2, X3), 0.05)
mcm <- function(data_list, level){
  nvec <- rep(0, length(data_list))
  apmat <- c()
  for(i in (1:length(data_list))){
    nvec[i] <- nrow(data_list[[i]])
    apmat <- rbind(apmat, data_list[[i]])
  }
  K <- length(nvec)
  N <- sum(nvec)
  nvecsq <- nvec %*% t(nvec)
  sumninj <- sum(nvecsq[upper.tri(nvecsq, diag = F)])
  nnminus1 <- nvec * (nvec - 1)
  nnminus1vecsq <- nnminus1 %*% t(nnminus1)
  sumnnminus1 <- sum(nnminus1vecsq[upper.tri(nnminus1vecsq, diag = F)])
  s1 <- 0
  if(K >= 3){
    for(i in 1:(K - 2)){
      for(j in (i+1):(K-1)){
        for(k in (j+1):(K)){
          s1 <- s1+((nvec[i])*(nvec[j])*(nvec[k])*(nvec[i]-1))+((nvec[i])*(nvec[j])*(nvec[k])*(nvec[j]-1))+((nvec[i])*(nvec[j])*(nvec[k])*(nvec[k]-1))
        }
      }
    }
  }
  s2 <- 0
  if(K >= 4){
    for(i in 1:(K-3)){
      for(j in (i+1):(K-2)){
        for(k in (j+1):(K-1)){
          for(l in (k+1):K){
            s2 <- s2+((nvec[i])*(nvec[j])*(nvec[k])*(nvec[l]))
          }
        }
      }
    }
  }
  nullmean <- sumninj/(N - 1)
  nullvar <- (1/((N-1)*(N-3)))*(sumnnminus1+(2*s1)+(6*s2))-(((N-2)/(N*(N-1)^2))*sumninj^2)+((sumninj/(N-1))*(1-(2*sumninj/(N^2-N))))
  if(!require(nbpMatching)) install.packages('nbpMatching')
  smatch <- as.matrix(nbpMatching::nonbimatch(nbpMatching::distancematrix(as.matrix(dist(apmat,method = 'euclidean',diag = T,upper = T,p = 2))))$matches)
  multcm <- 0
  cs <- c(0, cumsum(nvec))
  for(k in 1:K){
    for(j in (cs[k]+1):(cs[k+1])){
      multcm <- multcm + ((as.numeric(smatch[j,4]) <= cs[k]) || (as.numeric(smatch[j,4]) > cs[k+1]))
    }
  }
  multcm <- multcm/2
  multstat <- (multcm - nullmean)/sqrt(nullvar)
  lowerpval <- pnorm(multstat)
  dec <- noquote('Accept')
  if(lowerpval < level) dec <- noquote('Reject')
  return(noquote(c(lowerpval, dec)))
}


#' @title Computes the BCS
#' @description This function compute the beyondcell score (BCS) and returns an object of class \code{\link[beyondcell]{beyondcell}}.
#' @param sc \code{\link[Seurat]{Seurat}} object or expression matrix.
#' @param gs \code{\link[beyondcell]{geneset}} object.
#' @param expr.thres Minimum fraction of signature genes that must be expressed in a cell to compute its BCS.
#' Cells with a number of expressed genes below this fraction will have a \code{NaN} BCS.
#' @return A \code{beyondcell} object.
# ---- beyondcell ----
bcScore <- function(sc, gs, expr.thres = 0.1){
  if('Seurat' %in% class(sc)){
    input <- 'Seurat object'
    default <- Seurat::DefaultAssay(sc)
    if(default %in% c('RNA', 'SCT')){
      if('data' %in% slotNames(sc@assays[[default]])){
        message(paste0('Using ', default, ' assay as input.'))
        expr.matrix <- as.matrix(Seurat::GetAssayData(sc, slot = 'data', assay = default))
      }else{
        stop('Default assay must include a normalized data (@data) slot.')
      }
    }else{
      stop('Seurat default assay must be either RNA or SCT.')
    }
  }else if('matrix' %in% class(sc) & is.numeric(sc)){
    input <- 'expression matrix'
    warning('Using count matrix as input. Please, check that this matrix is normalized and unscaled.')
    expr.matrix <- sc
    sc <- Seurat::CreateSeuratObject(expr.matrix)
  }else{ stop('sc must be either a Seurat object or a single-cell expression matrix.')}
  if(class(gs) != 'geneset') stop('gs must be a geneset object.')
  if(length(expr.thres) != 1 | expr.thres[1] < 0 | expr.thres[1] > 1){
    stop('expr.thres must be a positive number between 0 and 1.')
  }
  # check that gene names are in the same format
  sc.gene.case <- names(which.max(CaseFraction(rownames(expr.matrix))))
  gs.gene.case <- names(which.max(CaseFraction(unique(unlist(gs@genelist)))))
  if(sc.gene.case != gs.gene.case){
    warning(paste0('gs genes are ',sc.gene.case,' and sc genes are ',gs.gene.case,'. Please check your ', input,
                   ' and translate the genes if necessary.'))
  }
  
  bc <- beyondcell(expr.matrix = expr.matrix, meta.data = sc@meta.data,
                   SeuratInfo = list(reductions = sc@reductions),
                   regression = list(order = rep('',2),vars=NULL,order.background=rep('',2)),
                   n.genes = gs@n.genes, mode = gs@mode, thres = expr.thres)
  rownames(expr.matrix) <- tolower(rownames(expr.matrix))
  gs@genelist <- lapply(gs@genelist, function(x){
    if(is.null(x[['up']])) up <- NULL else up <- tolower(x[['up']])
    if(is.null(x[['down']])) down <- NULL else down <- tolower(x[['down']])
    list(up = up, down = down)
  })
  # Genes in expr.matrix
  genes <- rownames(expr.matrix)
  len.gs <- length(gs@genelist)
  total <- len.gs + (length(gs@mode) * len.gs)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  bins <- ceiling(total/100)
  below.thres <- t(
    sapply(1:len.gs, function(i){
      all.genes <- unique(unlist(gs@genelist[[i]]))
      sub.expr.matrix <- expr.matrix[genes %in% all.genes,,drop=F]
      n.expr.genes <- apply(sub.expr.matrix, 2, function(x) sum(x > 0))
      if(i %% bins == 0){
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = i)
      }
      return(n.expr.genes < (length(all.genes) * expr.thres))
    })
  )
  rownames(below.thres) <- names(gs@genelist)
  below.thres <- below.thres[,colnames(expr.matrix)]
  # If all cells are below the threshold, remove that signature and raise a warning
  nan.rows.idx <- which(rowSums(below.thres) == ncol(below.thres))
  if(length(nan.rows.idx) > 0){
    warning(paste0('The following signatures have no cells that pass the expr.thres and will be removed: ',
                   paste0(rownames(below.thres)[nan.rows.idx],collapse = ', '), '.'))
    below.thres <- below.thres[-nan.rows.idx,]
    gs@genelist <- gs@genelist[-nan.rows.idx]
  }
  # BCS
  bcs <- lapply(seq_along(gs@mode), function(j){
    score <- t(sapply(1:length(gs@genelist),function(k){
      # common genes between the expr.matrix and each signature
      common.genes <- unique(intersect(genes, gs@genelist[[k]][[gs@mode[j]]]))
      # subset expr.matrix with common.genes
      sub.expr.matrix <- expr.matrix[common.genes,,drop = F]
      # sum expression of those genes for each cell
      sum.expr <- colSums(sub.expr.matrix)
      # raw score (mean)
      raw <- colMeans(sub.expr.matrix)
      # stdev (for BCS normalization)
      sig.stdev <- apply(sub.expr.matrix, 2, sd)
      # normalized BCS
      norm.score <- raw * ((sum.expr - sig.stdev)/(raw + sig.stdev))
      # update the progress bar
      step <- len.gs + (j - 1) * len.gs + k + length(nan.rows.idx)
      if(step %% bins == 0 | step == total){
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = step)
      }
      return(norm.score)
    }))
    rownames(score) <- names(gs@genelist)
    return(score[, colnames(expr.matrix), drop = F])
  })
  names(bcs) <- gs@mode
  if(length(gs@mode)==2){
    nan.cells <- is.na(bcs[['up']]) & is.na(bcs[['down']])
    bcs[['up']][is.na(bcs[['up']])] <- 0
    bcs[['down']][is.na(bcs[['down']])] <- 0
    scoring.matrix <- bcs[['up']] - bcs[['down']]
    scoring.matrix[nan.cells | below.thres] <- NaN
  }else{
    scoring.matrix <- bcs[[gs@mode]]
    scoring.matrix[below.thres] <- NaN
    if(gs@mode == 'down') scoring.matrix <- 1 * scoring.matrix[,,drop = F]
  }
  # If genesets were obtained in a TREATED vs CONTROL comparison, invert BCS sign (exclude pathways)
  if(gs@comparison == 'treated_vs_control'){
    load('data/beyondcell/pathways.RData')
    not.paths <- which(!(rownames(scoring.matrix) %in% names(pathways)))
    paths <- which(rownames(scoring.matrix) %in% names(pathways))
    if(length(paths) > 0){
      scoring.matrix <- rbind((scoring.matrix[not.paths,,drop=F]* -1), scoring.matrix[paths,,drop=F])
    }else{
      scoring.matrix <- (-1) * scoring.matrix[,,drop = F]
    }
  }
  slot(bc, 'normalized') <- slot(bc, 'data') <- round(scoring.matrix, digits = 2)
  scaled.matrix <- t(apply(bc@normalized, 1, scales::rescale, to = c(0,1)))
  slot(bc, 'scaled') <- round(scaled.matrix, digits = 2)
  switch.point <- SwitchPoint(bc)
  slot(bc, 'switch.point') <- switch.point
  close(pb)
  return(bc)
}

#' @title Returns the fraction of each case type in the input
#' @description This function computes the fraction of each case type (uppercase, lowercase or capitalized)
#' in a character vector.
#' @param x Character vector
#' @return A named numeric vector with the fraction of each case type.
CaseFraction <- function(x){
  if(!is.character(x)) stop('x must be a character vector.')
  if(!require(useful)) install.packages('useful')
  perc <- sapply(c('upper', 'lower', 'capitalized'), function(case){
    if(case == 'upper' | case == 'lower'){
      p <- sum(useful::find.case(x, case = case))/length(x)
    }else{
      splitted <- strsplit(x, split = '')
      first.letter <- sapply(splitted, `[[`, 1)
      rest.letters <- sapply(splitted, function(y) paste0(y[2:length(y)], collapse = ''))
      p <- sum(useful::upper.case(first.letter) & useful::lower.case(rest.letters))/length(x)
    }
    return(round(p, digits = 2))
  })
  names(perc)[1:2] <- paste0('in ', names(perc)[1:2], 'case')
  return(perc)
}

#' @title Computes the switch point
#' @description This function computes the switch point of the signatures of a give \code{\link[beyondcell](beyondcell)}
#' object. The switch point is the (subsetted and/or regressed) scaled beyondcell score (BCS) that corresponds to the
#' point in which the normalized BCS in \code{beyondcell@data} switches from negative (insensitive) to positive (sensitive)
#' values. The closer to 0, the more sensitive are the cells to a given drug.
#' @param bc \code{beyondcell} object.
#' @return A named vector with the switch points of the signatures in \code{bc}.
SwitchPoint <- function(bc){
  if(class(bc) != 'beyondcell') stop('bc must be a beyondcell object.')
  if(any(!(bc@mode %in% c('up','down'))) | all(is.null(bc@mode))) stop('Incorrect mode.')
  bc@mode <- unique(bc@mode)
  
  sigs <- rownames(bc@normalized)
  cells <- colnames(bc@normalized)
  if(all(bc@mode == 'up')){
    switch.point <- rep(0, times = length(sigs))
  }else if(all(bc@mode == 'down')){
    switch.point <- rep(1, times = length(sigs))
  }else if(length(bc@mode) == 2){
    indexes <- lapply(sigs, function(x){
      m <- bc@data[x, cells]
      if(all(is.na(m))) return(NaN)
      m.nona <- na.omit(m)
      if(all(m.nona >= 0)){
        return(0)
      }else if(all(m.nona <= 0)){
        return(1)
      }else{
        exact.0 <- m.nona == 0
        if(any(exact.0)){
          return(rep(which(m == 0)[1], times = 2))
        }else{
          lower.bound <- which(m == max(m.nona[m.nona <= 0]))[1]
          upper.bound <- which(m == min(m.nona[m.nona >= 0]))[1]
          return(c(lower.bound, upper.bound))
        }
      }
    })
    names(indexes) <- sigs
    switch.point <- sapply(names(indexes),function(y){
      if(length(indexes[[y]]) == 2){
        return(round(sum(bc@scaled[y, indexes[[y]]])/2, digits = 2))
      }else{ return(indexes[[y]]) }
    })
  }
  names(switch.point) <- sigs
  return(switch.point)
}

#' @title Compute a UMAP projection and therapeutic clusters using the BCS
#' @description This function uses the beyondcell scores (BCS) to compute a UMAP projection of the data and to
#' cluster cells according to their sensitivity to the tested drugs (therapeutic clusters).
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param pc Number of principal components to use (which can be estimated from an elbow plot). If \code{pc = NULL}
#' (default), the function will stop prior to compute the UMAP projection and the therapeutic clusters.
#' @param k.neighbors (\code{\link[Seurat]{FindNeighbors}}' \code{k.param}). Defines \emph{k} for the k-Nearest
#' Neighbour algorithm.
#' @param res (\code{\link[Seurat]{FindClusters}}' \code{resolution}) Value of the resolution parameter, use a value
#' above/below 1.0 if you want to obtain a larger/smaller number of communities. Can be a single number of a numeric vector.
#' @param add.DSS Use background BCS computed with \code{DSS} signatures (\code{add.DSS = TRUE}) or just use the 
#' signatures included in the \code{bc} object (\code{add.DSS = FALSE}) to compute the UMAP projection and the therapeutic
#' clusters. If the number of drugs in \code{bc}.
#' (excluding pathways) is <= 20, it is recomended to set \code{add.DSS = TRUE}. Note that if \code{add.DSS = TRUE},
#' the regression and susbet steps that have been applied on \code{bc} will also be applied on the background BCS.
#' @param method (\code{\link[Seurat]{RunUMAP}}'s \code{umap.method}) UMAP implementation to run. Can be:
#' \itemize{
#'   \item{\code{uwot}}: Runs the UMAP via the \code{\link[uwot]{umap}} function of the \code{uwot} package.
#'   \item{\code{umap-learn}}: Runs the UMAP via the \code{Seurat} wrapper of the python \code{umap-learn} package.
#' }
#' @param return.model (\code{\link[Seurat]{RunUMAP}}'s \code{return.model}) Whether \code{RunUMAP} will return the
#' \code{uwot} model.
#' @details Thi function performs all the steps required to obtain a UMAP reduction of the data and cluster the cells
#' according to the BCS.
#' You will normally require to run the function twice:
#' \enumerate{
#'   \item Using \code{pc = NULL} to obtain the elbow plot.
#'   \item Specifying a value for the \code{pc} parameter according to this plot.
#'   This second time, the UMAP reduction and the therapeutic clusters will be computed.
#' }
#' Note that \code{add.DSS} must be the same in both runs, so the elbow plot obtained in 1 is still valid in 2. If
#' \code{add.DSS = TRUE}, the background BCS will be stored in the \code{bc} object and the function will skip this
#' step the second time.
#' @return A \code{beyondcell} object with the UMAP reduction in \code{@reductions} slot and the therapeutic clusters
#' for each \code{res} in \code{bc@meta.data}. Also, an elbow plot (\code{\link[ggplot2]{ggplot2}} object) is printed
#' (if \code{elbow.path = NULL}) or saved (if \code{elbow.path != NULL}).
bcUMAP <- function(bc, pc = NULL, k.neighbors = 20, res = 0.2, add.DSS = FALSE, 
                   method = 'uwot', return.model = FALSE){
  if(class(bc) != 'beyondcell') stop('bc must be a beyondcell object.')
  if(!is.null(pc)){if(length(pc) != 1 | pc[1] < 2) stop('pc must be an integer >= 2.')}
  if(length(k.neighbors) != 1 | k.neighbors < 1) stop('k.neighbors must be a positive integer.')
  if(any(sapply(res,function(x) !is.numeric(x))) | any(sapply(res,function(y) y < 0)))
    stop('res must be a vector of numbers >= 0.')
  # Check add.DSS
  load('data/beyondcell/pathways.RData')
  not.paths <- !(rownames(bc@normalized) %in% names(pathways))
  n.drugs <- sum(not.paths)
  if(length(add.DSS) != 1 | !is.logical(add.DSS)){
    stop('add.DSS must be TRUE or FALSE')
  }else if(!add.DSS){
    if(n.drugs <= 10){
      stop(paste('Only',n.drugs,'dug signatures (excluding pathways) are present in the bc object, please set add.DSS = TRUE.'))
    }else if(n.drugs <= 20){
      warning(paste('Computing an UMAP reduction for',n.drugs, 'drugs. We recommend to set add.DSS = TRUE',
                    'when the number of sigantures (excluding pathways) is below or equal to 20.'))
    }
  }
  if(length(method) != 1 | !(method[1] %in% c('uwot','umap-learn'))) stop('Incorrect method.')
  if(length(return.model) != 1 | !is.logical(return.model)) stop('return.model must be TRUE or FALSE')
  if(method=='umap-learn' & return.model == T){
    warning('return.model = TRUE is only valid when method = "umap-learn". Changing return.model to FALSE')
    return.model <- FALSE
  }
  
  cells <- colnames(bc@normalized)
  if(add.DSS){
    if(!identical(sort(rownames(bc@background),decreasing = F),sort(DSS[[1]]$sig_id,decreasing = F)) |
       !identical(sort(colnames(bc@background),decreasing = F),sort(cells, decreasing = F)) |
       !identical(bc@regression$order, bc@regression$order.background)){
      message('Computing background BCS using DSS signatures...')
      # Genesets
      gs.background <- GenerateGenesets(DSS, n.genes = bc@n.genes, mode = bc@mode, include.pathways = F)
      # BCS
      background <- bcScore(bc@expr.matrix, gs = gs.background, expr.thres = bc@thres)
      # Add metadata
      background@meta.data <- background@meta.data[, -c(1:ncol(background@meta.data))]
      background <- bcAddMetadata(background, metadata = bc@meta.data)
      # Subset and regress (if needed)
      if(bc@regression$order[1] == 'subset'){
        background <- bcSubset(background, cells = cells)
      }else if(bc@regression$order[1] == 'regression'){
        message('Regressing background BCS...')
        background <- bcRegressOut(background, vars.to.regress = bc@regression[['vars']])
      }
      if(bc@regression$order[2] == 'subset'){
        background <- bcSubset(background, cells = cells)
      }else if(bc@regression$order[2] == 'regression'){
        message('Regressing background BCS...')
        background <- bcRegressOut(background, vars.to.regress = bc@regression[['vars']])
      }
      # Add background@normalized to bc@background
      bc@background <- background@normalized
      # Add order.background to bc@regression
      bc@regression[['order.background']] <- bc@regression[['order']]
    }else{
      message('Background BCS already computed. Skipping this step.')
    }
    # Add background to bc.
    all.rows <- unique(c(rownames(bc@normalized), rownames(bc@background)))
    merged.score <- rbind(bc@normalized, bc@background[, cells])[all.rows,]
    # scale
    merged.score <- t(apply(merged.score, 1, scales::rescale, to = c(0,1)))
    bc.merged <- beyondcell(scaled = merged.score)
  }else{
    message('DSS background not computed. UMAP will be created just with the drugs (not pathways) in bc object.')
    bc.merged <- bc
  }
  sc <- Seurat::CreateSeuratObject(bc.merged@scaled[not.paths,,drop = F])
  sc <- Seurat::ScaleData(sc, features = rownames(sc), do.scale = F, do.center = F)
  sc <- Seurat::RunPCA(sc, features = rownames(sc), npcs = 100, maxit = 100000)
  if(!is.null(pc)){
    message('Obtaining therapeutic clusters...')
    sc <- Seurat::FindNeighbors(sc, dims = 1:pc, k.param = k.neighbors)
    sc <- Seurat::FindClusters(sc, resolution = res)
    message('Computing beyondcell\'s UMAP reduction...')
    sc <- Seurat::RunUMAP(sc, dims = 1:pc, umap.method = method, return.model = return.model,n.components = 2, verbose = F)
    bc@reductions <- sc@reductions
    # Therapeutic clusters
    message('Adding therapeutic clsuters to metadata...')
    meta <- sc@meta.data[,paste0('RNA_snn_res.',res),drop = F]
    colnames(meta) <- paste0('bc_clusters_res.', res)
    # bc@metadata without repeated beyondcell clusters
    repeated.res <- colnames(bc@meta.data) %in% colnames(meta)
    bc.metadata <- bc@meta.data[, !repeated.res, drop = F]
    # bc@meta.data rownames
    meta.cells <- rownames(bc@meta.data)
    # merge bc@meta.data (without repeated beyondcell clusters) and meta. If the bc object was subsetted and
    # some cells removed, their bc_clusters will be NA.
    new.meta <- transform(merge(bc.metadata, meta, by = 0, all.x = T),row.names = Row.names, Row.names = NULL)[meta.cells,]
    bc@meta.data <- new.meta
  }
  return(bc)
}

#' @title Ranks the signatures from most sensitive to least sensitive
#' @description This function computes the beyondcell score's (BCS) statistics of each signature and ranks them
#' according to the switch point and mean.
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param idents Name of the metadata column of interest. I \code{idents = NULL}, the function computes the ranks
#' using all cells. If \code{idents != NULL}, the signatures' ranks are computed for each level in \code{idents}.
#' @param extended If \code{extended = TRUE}, this function returns the switch point, mean, median, standard
#' deviaiton, variance, min, max, proportion of \code{NaN} and residuals' mean per signature. If \code{extended = F},
#' this function returns only the switch point, mean and residuals' mean.
#' @return A \code{beyondcell} object with the results in a new entry of \code{bc@ranks}: 
#' \code{bc@ranks[["general"]]} (if \code{idents = NULL}) or \code{bc@ranks[[idents]]} (if \code{idents != NULL}).
bcRanks <- function(bc, idents = NULL, extended = T){
  if(class(bc) != 'beyondcell') stop('bc must be a beyondcell object.')
  if(!is.null(idents)){
    if(length(idents) != 1){stop('Idents must be a single metadata column.')}
    if(idents %in% colnames(bc@meta.data)){
      if(idents %in% names(bc@ranks)){
        warning(paste0('$', idents, ' already exists in bc@ranks. Entry has been overwritten.'))
      }
      meta <- bc@meta.data[colnames(bc@normalized), idents, drop = F]
    }else{ stop('Idents not found.') }
  }else{
    if('general' %in% names(bc@ranks)){
      warning('$general already exists in bc@ranks. Entry has been overwritten.')
    }
  }
  if(length(extended) != 1 | !is.logical(extended[1])) stop('extended must be TRUE or FALSE')
  # Signatures in bc
  sigs <- rownames(bc@normalized)
  # n.rows <- length(sigs)
  cells <- colnames(bc@normalized)
  n <- ifelse(extended, 6, 1)
  if(is.null(idents)){
    total <- length(sigs) * n
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    idents <- 'general'
    order.col <- 'rank'
    final.stats <- GetStatistics(bc = bc, signatures = sigs, cells = cells, pb = pb, total = total, 
                                 i = 1, n.rows = length(sigs), extended = extended)
    sig.order <- order(-1 * final.stats$switch.point, final.stats$mean, decreasing = T)
    final.stats$rank[sig.order] <- 1:nrow(final.stats)
    final.stats <- final.stats[c('rank', colnames(final.stats)[-ncol(final.stats)])]
  }else{
    lvls <- sort(unique(as.factor(meta[[idents]])), decreasing = F)
    total <- length(sigs) * n * length(lvls)
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    order.col <- paste0('rank.',lvls[1])
    stats <- lapply(seq_along(lvls),function(i){
      group.cells <- rownames(bc@meta.data)[rownames(bc@meta.data) %in% cells & bc@meta.data[,idents]==lvls[i]]
      group.cells <- group.cells[!is.na(group.cells)]
      sub.bc <- bc
      sub.bc@regression <- list(order = rep('',2), vars = NULL, order.background = rep('',2))
      sub.bc@background <- matrix(ncol = 0, nrow = 0)
      sub.bc <- bcSubset(sub.bc, cells = group.cells)
      out <- GetStatistics(bc = sub.bc, signatures = sigs, cells = group.cells,pb = pb, total = total,
                           i = i, n.rows = length(sigs), extended = extended)
      # Add 4 squares group annotation
      # Get info about drugs (their corresponding name in bc, the preferred name used by beyondcell and the MoA).
      info <- FindDrugs(sub.bc, x = rownames(sub.bc@scaled), na.rm = F)
      sp <- data.frame(switch.point = bc@switch.point[info$bc.names], row.names = info$bc.names)
      res <- out[, 'residuals.mean', drop = F]
      df <- transform(merge(res, sp, by = 0), row.names = Row.names, Row.names = NULL)
      res.decil <- quantile(as.numeric(out$residuals.mean), prob = seq(from = 0, to = 1, length =11))
      # Group annotation
      sp_lower_01 <- as.numeric(df$switch.point) < 0.1
      sp_lower_06 <- as.numeric(df$switch.point) < 0.6
      sp_higher_04 <- as.numeric(df$switch.point) > 0.4
      sp_higher_09 <- as.numeric(df$switch.point) > 0.9
      res_lower_10 <- as.numeric(df$residuals.mean) < res.decil[['10%']]
      res_higher_90 <- as.numeric(df$residuals.mean) > res.decil[['90%']]
      df$group <- rep('', times = nrow(df))
      df$group[sp_lower_01 & res_higher_90] <- 'TOP-HighSensitivity'
      df$group[sp_higher_09 & res_lower_10] <- 'TOP-LowSensitivity'
      df$group[sp_higher_04 & sp_lower_06 & res_lower_10] <- 'TOP-Differential-LowSensitivity'
      df$group[sp_higher_04 & sp_lower_06 & res_higher_90] <- 'TOP-Differential-HighSensitivity'
      out$group <- df$group
      sig.order <- order(-1 * out$switch.point, out$mean, decreasing = T)
      out$rank[sig.order] <- 1:nrow(out)
      out <- out[c('rank', colnames(out)[-ncol(out)])]
      colnames(out) <- paste0(colnames(out), '.', lvls[i])
      return(out)
    })
    final.stats <- do.call(cbind.data.frame, args = stats)
  }
  # Add Drug name and MoA to final.stats
  cols <- colnames(final.stats)
  load('data/beyondcell/drugInfo.RData')
  info <- subset(drugInfo, subset = drugInfo$IDs %in% rownames(final.stats))
  if(dim(info)[1] > 0){
    info <- aggregate(.~IDs, data = info, na.action = NULL, function(x){ paste(na.omit(unique(x)),collapse = '; ')})
  }
  rownames(info) <- info$IDs
  info <- info[, c('drugs', 'preferred.drug.names', 'MoAs', 'targets', 'sources')]
  final.stats <- transform(merge(info,final.stats,by = 0, all.x = T), row.names = Row.names, Row.names = NULL)
  # final.stats <- final.stats[order(final.stats[,order.col],decreasing = F),
  #                            c('drugs','preferred.drug.names','MoAs','targets','sources',cols)]
  # Add to beyondcell object
  bc@ranks[[idents]] <- final.stats
  Sys.sleep(0.1)
  close(pb)
  return(bc)
}

#' @title Returns the first/last n ranked signatures
#' @description This function returns the top/bottom \code{n} signatures ranked by \code{\link[beyondcell]{bcRanks}}.
#' If the rank has not been previously computed, \code{rankSigs} performs the ranking itself.
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param idents Name of the metadata column of interest. If \code{idents = NULL}, the function uses the geenral rank
#' computed with all cells.
#' @param cond Level of \code{idents} to rank by the output vector. If \code{idents = NULL}, this parameter is deprecated.
#' @param n Number of signatures to return in the output vector.
#' @param decreasing Logical, return the top \code{n} signatures (default) or the bottom \code{n} signatures (\code{decreasing = F})
#' @return An ordered vector with the signature's names
rankSigs <- function(bc,idents = NULL, cond = NULL, n = 10, decreasing = T){
  if(class(bc) != 'beyondcell') stop('bc must be a beyondcell object.')
  if(!is.null(idents)){
    if(length(idents) != 1) stop('Idents must be a single metadata column.')
    if(!idents %in% colnames(bc@meta.data)) stop('Idents not found.')
    if(is.null(cond[1])) stop('Invalid cond.')
    meta <- idents
  }else{
    meta <- 'general'
    if(!is.null(cond[1])){
      warning('idents not specified, cond is deprecated.')
      cond <- NULL
    }
  }
  if(!is.null(cond)){
    if(length(cond) != 1) stop('cond must be a single idents level')
    if(!cond %in% levels(as.factor(bc@meta.data[,idents]))) stop(paste0(cond,' is not a level of ',idents,'.'))
  }
  if(length(n) != 1 (!is.null(n) & !is.character(n))) stop('n must be a single number of "all".')
  if(is.numeric(n) & (n < 1 | n %% 1 != 0)){
    stop('n must be an integer >0.')
  }else if(is.character(n)){
    if(n == 'all'){
      n <- nrow(bc@normalized)
    }else{
      stop('To select all signatures, please set n = "all".')
    }
  }
  if(length(decreasing) != 1 | !is.logical(decreasing)) stop('decreasing must be TRUE or FALSE.')
  
  # If ranks have not been computed, compute them now
  if(!meta %in% names(bc@ranks)){
    message('Computing ranks...')
    bc <- bcRanks(bc, idents = idents, extended = T)
  }
  # Get ranks for the specified idents
  df <- bc@ranks[[meta]]
  idx <- ifelse(decreasing, 1:n, nrow(bc@normalized):(nrow(bc@normalized)-n+1))
  # Return signatures whose rank == idx
  order.col <- ifelse(meta == 'general', 'rank', paste0('rank.',cond))
  ordered.df <- df[order(df[, order.col], decreasing = F),]
  sigs <- ordered.df[idx, 'Name']
  return(sigs)
}

#' @title Plot a UMAP reduction coloured by metadata information.
#' @description This fucntion returns a {\link[ggplot2]{ggplot2}} object with a UMAP reduction colored by
#' the specified metadata column.
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param idents Name of the metadata column to color by.
#' @param UMAP UMAP reduction to plot. Either \code{beyondcell}, computed using \code{\link[beyondcell]{bcUMAP}},
#' or \code{Seurat}, obtained using \code{Seurat}'s functions.
#' @param factor.col Logical indicating if \code{idents} column is a factor or not. Set \code{factor.col = FALSE}
#' if \code{idents} is a numeric column (such as \code{percent.mt} or \code{nFeature_RNA}).
#' @param ... Other arguments passed to \code{\link[Seurat]{DimPlot}}.
#' @return A \code{ggplot2} object with the UMAP reduction colored by \code{idents}.
bcClusters <- function(bc, idents, UMAP = 'beyondcell', factor.col = T, ...){
  if(class(bc) != 'beyondcell') stop('bc must be a beyondcell object.')
  if(length(idents) != 1) stop('Idents must be a single metadata column.')
  if(idents %in% colnames(bc@meta.data)){
    meta <- bc@meta.data[colnames(bc@scaled), idents, drop = F]
  }else{ stop('Idents not found.') }
  if(UMAP == 'beyondcell'){
    if(length(bc@reductions) == 0){
      stop('You must precompute beyondcell\'s UMAP projection using bcUMAP().')
    }
    sc <- Seurat::CreateSeuratObject(bc@scaled)
    reduction <- bc@reductions
  }else if(UMAP == 'Seurat'){
    if(length(bc@SeuratInfo$reductions)==0){
      stop('No UMAP projection available for your Seurat\'s object.')
    }
    sc <- Seurat::CreateSeuratObject(bc@expr.matrix[, colnames(bc@scaled)])
    reduction <- bc@SeuratInfo$reductions
  }else{
    stop('Incorrect UMAP argument. Please user either "Seurat" or "beyondcell".')
  }
  if(length(factor.col) != 1 | !is.logical(factor.col)) stop('factor.col must be TRUE or FALSE.')
  
  sc <- Seurat::AddMetaData(sc, metadata = meta)
  sc@reductions <- reduction
  if(factor.col){
    Seurat::Idents(sc) <- idents
    p <- Seurat::DimPlot(sc, reduction = 'umap', ...) + ggplot2::theme_minimal()
  }else{
    p <- Seurat::FeaturePlot(sc, reduction = 'umap', features = idents, ...) +
      ggplot2::theme_minimal() + ggplot2::labs(title = NULL)
  }
  return(p)
}

#' @title Plot a histogram with the BCS of the signature of interest
#' @description This function draw a histogram of beyondcell scores (BCS) for each signature of interest.
#' The plot can be a single histogram or a histogram for each level found in \code{idents}.
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param signatures Vector with the names of the signatures of interest. If \code{sigantures = "all"},
#' all signatures are selected.
#' @param idents Name of the metadata column of interest. If \code{idents = NULL}, a single histogram with all
#' BCS is drawn. On the other hand, if \code{idents != NULL} a histogram for each level found in \code{idents}
#' will be drawn.
#' @return A list of \code{\link[ggplot2]{ggplot2}} histograms, one for each signature of interest. In each
#' histogram, the median, mean and sd are reported. Also, the mean is indicated with a black dashed line and
#' the median with a red dashed line.
bcHistogram <- function(bc, signatures, idents = NULL){
  if(class(bc) != 'beyondcell') stop('bc must be a beyondcell object.')
  if(!is.character(signatures)) stop('Signatures must be a character vector.')
  if(length(signatures)==1 & signatures[1]=='all'){
    signatures <- rownames(bc@normalized)
    in.signatures <- rep(TRUE, times = nrow(bc@normalized))
  }else{
    in.signatures <- !is.null(signatures) & signatures %in% rownames(bc@normalized)
    if(all(!in.signatures)){
      stop('None of the specified signatures were found.')
    }else if(any(!in.signatures)){
      warning(paste0('There signatures were not found in bc: ',paste0(signatures[!in.signatures],collapse = ', '),'.'))
    }
  }
  if(!is.null(idents)){
    if(length(idents) != 1) stop('Idents must be a single metadata column.')
    if(idents %in% colnames(bc@meta.data)){
      meta <- paste(idents, '=', bc@meta.data[colnames(bc@normalized), idents,drop = T])
    }else{ stop('Idents not found.') }
  }else{
    meta <- rep('', times = ncol(bc@normalized))
  }
  
  lvls <- levels(as.factor(meta))
  sub.bc <- bc@data[signatures[in.signatures],,drop = F]
  # get maximum and minimum normalized BCS (for common x axis in all plots)
  limits <- c(min(as.vector(sub.bc),na.rm = T),max(as.vector(sub.bc),na.rm = T))
  # get info about drugs (their corresponding name in bc, the preferred name used by beyondcell and the MoA).
  info <- FindDrugs(bc, x = signatures[in.signatures])
  p <- lapply(signatures[in.signatures],function(x){
    sub.df <- na.omit(data.frame(bcscore = sub.bc[x,],condition = as.factor(meta),row.names = colnames(sub.bc)))
    stats.cond <- sapply(lvls, function(y){
      stats <- round(Mean.Med.SD(subset(sub.df$bcscore, subset = sub.df$condition == y)), digits = 2)
      return(stats)
    })
    stats.labels <- data.frame(label = apply(stats.cond,2,function(z){
      paste0(rownames(stats.cond), ' = ', z, collapse = '\n')
    }), mean = stats.cond['mean',], median = stats.cond['median',], condition = colnames(stats.cond))
    # Drug name and MoA
    if(x %in% info$IDs){
      drug.and.MoA <- info[which(info$IDs == x), c('preferred.and.sigs', 'MoAs')]
      drug.and.MoA[2] <- ifelse(drug.and.MoA[2] == 'NA', '', drug.and.MoA[2])
    }else{ drug.and.MoA <- c(x, '') }
    hist <- ggplot(sub.df, aes(x = bcscore, fill = condition)) +
      geom_histogram(data = transform(sub.df, condition = NULL),fill = 'grey85',alpha = 0.3,binwidth = 1) +
      facet_wrap(vars(condition), ncol = 1) +
      geom_histogram(color = '#F0F0F0',alpha = 0.6,position = 'identity',binwidth = 1) + 
      theme_minimal() + theme(legend.position = 'none') +
      geom_label(data = stats.labels, hjust = 1, vjust = 1, size = 3, x = Inf, y = Inf, 
                 label = stats.labels$label, fill = scales::hue_pal()(length(lvls))) +
      labs(title = drug.and.MoA[1], subtitle = drug.and.MoA[2]) +
      geom_vline(data = stats.labels, aes(xintercept = mean), linetype = 'dashed') +
      geom_vline(data = stats.labels, aes(xintercept = median), color = 'red', linetype = 'dashed')
    return(hist)
  })
  return(p)
}

#' @title Plot a UMAP reduction colored by BCS or gene expression values
#' @description This function returns a list of \code{\link[patchwork]{patchwork}} or \code{\link[ggplot2]{ggplot2}}
#' objects with the desired UMAP reduction colored by beyondcell scores (BCS) or gene expression values.
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param UMAP UMAP reduction to plot. Either \code{beyondcell}, computed using \code{\link[beyondcell]{bcUMAP}},
#' or \code{Seurat}, obtained using \code{Seurat}'s functions.
#' @param signatures List with plot parameters to color the UMAP by BCS:
#' \itmize{
#'   \item{\code{values}:} {Vector with the names of the signatures of interest. If \code{signatures[['values]]=='all'},
#'   all signatures are selected.}
#'   \item{\code{colorscale}:} {Either a \code{viridis}, \code{RColorBrewer} or a custom palette of 3 colors
#'   (low, medium and high) to color all signatures' plots. If \code{colorscale = NULL} (default), the plots
#'   are colored using \code{beyondcell}'s own palette.}
#'   \item{\code{alpha}:} {Transparency level between 1 (not transparent) and 0 (fully transparent).}
#'   \item{\code{na.value}:} {Color to use for missing values (\code{NA}s).}
#'   \item{\code{limits}:} {Vector with the desired limits for all signatures' plots.}
#'   \item{\code{center}:} {A single number indicating the centre of the \code{colorscale} for all signatures'
#'   plots. If \code{center = NULL} (default), the \code{center} for each signature is its switch point.}
#'   \item{\code{breaks}:} {A single number indicating the break size of the \code{colorscale}. Alternatively,
#'   it can be a vector with the desired breaks (which don't have to be symmetric or equally distributed).}
#' }
#' @param genes List with plot parameters to color the UMAP by gene expression values:
#' \itemize{
#'   \item{\code{values}:} {Vector with the names of the genes of interest. If \code{genes[["values"]]=="all"},
#'   all genes are selected.}
#'   \item{\code{limits}:} {Vector with the desired limits for all genes' plots. If \code{limits = c(NA,NA)}}
#'   (default), the \code{limits} are computed automatically.
#'   \item{\code{share.limits}:} {Logical argument. If \code{share.limits = TRUE}, all genes' plots will have
#'   the same \code{limits}. If \code{share.limits = FALE} (default), each gene plot will have its own 
#'   \code{limits}.}
#' }
#' @param merged If \code{merged != NULL}, two signatures will be superposed in the same plot. If \code{merged="direct"},
#' the signatures are assumed to have a direct relationship and the BCS will be added. On the other hand, if 
#' \code{merged = "indirect"}, the signatures are assumed to have an indirect relationship and their BCS will be
#' substracted.
#' @param blend Scale and blend expression values to visualize co-expression of two genes.
#' @param mfrow Numeric vector of the form \code{c(nr,nc)}. \code{nr} corresponds to the number of rows and \code{nc}
#' to the number of columns of the grid in which the plots will be drawn. If you want to draw the plots individually
#' set \code{mfrow = c(1,1)}.
#' @param ... Other arguments passed to \code{FeaturePlot}.
#' @details When \code{genes[['limits']] = c(NA, NA)}, \code{bcSignatures} computes the limits automatically.
#' You can make all plots share the same limits by specifying \code{genes[["share.limits"]] == TRUE}, or make
#' the function to compute the limits individually for each gene with \code{genes[["share.limits"]]=FALSE}.
#' Moreover, if you specify a value for \code{genes[["limits"]]}, \code{genes[["share.limits"]]} will
#' automatically set to \code{TRUE} and all plots will share those limits.
#' @return A list of \code{patchwork} (if \code{mfrow != c(1,1)}) or \code{ggplot2} objects (if \code{mfrow=c(1,1)})
#' of the desired UMAP reduction colored by the BCS (for signatures) or gene expression values (for genes).
bcSignatures <- function(bc, UMAP = "beyondcell",
                         signatures = list(values = NULL, colorscale = NULL,
                                           alpha = 0.7, na.value = "grey50",
                                           limits = c(0, 1), center = NULL,
                                           breaks = 0.1),
                         genes = list(values = NULL, limits = c(NA, NA),
                                      share.limits = FALSE),
                         merged = NULL, blend = FALSE, mfrow = c(1, 1), ...){
  if(UMAP == "beyondcell"){
    if(length(bc@reductions) == 0){
      stop('You must precompute beyondcell\'s UMAP projection using bcUMAP().')
    }
    reduction <- bc@reductions
    cells <- subset(rownames(bc@meta.data),subset = rownames(bc@meta.data) %in% colnames(bc@normalized))
  }else if(UMAP == "Seurat"){
    if (length(bc@SeuratInfo$reductions) == 0) {
      stop('No UMAP projection available for your Seurat\'s object.')
    }
    reduction <- bc@SeuratInfo$reductions
    cells <- rownames(bc@meta.data)
  }else{
    stop('Incorrect UMAP argument. Please use either "Seurat" or "beyondcell".')
  }
  # Check signatures' list values.
  default.sigs <- list(values = NULL, colorscale = NULL, alpha = 0.7,na.value = "grey", 
                       limits = c(0, 1), center = NULL,breaks = 0.1)
  selected.sigs <- names(signatures) %in% names(default.sigs)
  if(any(!selected.sigs)){
    warning(paste0('Incorrect entries in signatures: ',paste0(names(signatures)[!selected.sigs], collapse = ", ")))
  }
  signatures <- c(signatures,default.sigs[!(names(default.sigs) %in% names(signatures))])
  # Check genes' list values.
  default.genes <- list(values = NULL, limits = c(NA, NA), share.limits = FALSE)
  selected.genes <- names(genes) %in% names(default.genes)
  if(any(!selected.genes)){
    warning(paste0('Incorrect entries in genes: ',paste0(names(genes)[!selected.genes], collapse = ", ")))
  }
  genes <- c(genes, default.genes[!(names(default.genes) %in% names(genes))])
  # Check signatures and genes' values (features).
  if(is.null(signatures[["values"]]) & is.null(genes[["values"]])){
    stop('You must specify the signatures and/or genes of interest.')
  }
  # Check signatures' values.
  if(!is.null(signatures[["values"]])){
    if(length(signatures[["values"]]) == 1 & signatures[["values"]][1] == "all"){
      signatures[["values"]] <- rownames(bc@normalized)
      in.signatures <- rep(TRUE, times = nrow(bc@normalized))
    }else{
      in.signatures <- signatures[["values"]] %in% rownames(bc@normalized)
      if(all(!in.signatures)){
        stop('None of the specified signatures were found.')
      }else if (any(!in.signatures)){
        warning(paste0('These signatures were not found in bc: ',
                       paste0(signatures[["values"]][!in.signatures],collapse = ", "), '.'))
      }
    }
  }else{
    in.signatures <- NULL
  }
  # Check genes' values.
  if(!is.null(genes[["values"]])){
    if(length(genes[["values"]]) == 1 & genes[["values"]][1] == "all") {
      genes[["values"]] <- rownames(bc@expr.matrix)
      in.genes <- rep(TRUE, times = nrow(bc@expr.matrix))
    }else{
      in.genes <- toupper(genes[["values"]]) %in% toupper(rownames(bc@expr.matrix))
      if(all(!in.genes)){
        stop('None of the specified genes were found.')
      }else if (any(!in.genes)){
        warning(paste0('These genes were not found in bc@expr.matrix: ',
                       paste0(genes[["values"]][!in.genes], collapse = ", "), '.'))
      }
    }
  }else{
    in.genes <- NULL
  }
  # Sigs, gene and features.
  sigs <- unique(signatures[["values"]][in.signatures])
  gene <- unique(genes[["values"]][in.genes])
  features <- c(sigs, gene)
  # Check signature's colorscale.
  signatures[["colorscale"]] <- get_color_steps(signatures[["colorscale"]])
  # Check signatures' alpha, na.value and breaks -> inside center_scale_colour_stepsn().
  # Check signatures' limits.
  if(length(signatures[["limits"]]) != 2){
    stop('Signatures\' limits must be a vector of length 2.')
  }
  if(!is.numeric(signatures[["limits"]]) | any(signatures[["limits"]] < 0)){
    stop('Signatures\' limits must be numeric (>= 0).')
  }
  if(signatures[["limits"]][2] < signatures[["limits"]][1]){
    warning(paste('Signatures\' upper limit is smaller than lower limit. Sorting limits in increasing order.'))
    signatures[["limits"]] <- sort(signatures[["limits"]], decreasing = FALSE)
  }
  # Check signature's center.
  if(!is.null(signatures[["center"]])){
    if(length(signatures[["center"]]) != 1 | !is.numeric(signatures[["center"]])) {
      stop('Signatures\' center must be a single number or NULL.')
    }
  }
  # Check genes' limits.
  if(length(genes[["limits"]]) != 2){
    stop('Genes\' limits must be a vector of length 2.')
  }
  na.limits.genes <- is.na(genes[["limits"]])
  if(length(genes[["limits"]][!na.limits.genes]) > 0 &
     (!is.numeric(genes[["limits"]][!na.limits.genes]) |
      any(genes[["limits"]][!na.limits.genes] < 0))){
    stop('Genes\' limits must be numeric (>= 0) or NAs.')
  }
  if(all(!na.limits.genes) & genes[["limits"]][2] < genes[["limits"]][1]){
    warning(paste('Genes\' upper limit is smaller than lower limit. Sorting limits in increasing order.'))
    genes[["limits"]] <- sort(genes[["limits"]], decreasing = FALSE)
  }
  # Check genes' share.limits.
  if(length(genes[["share.limits"]]) != 1 | !is.logical(genes[["share.limits"]])){
    stop('genes$share.limits must be TRUE or FALSE.')
  }
  if(!genes[["share.limits"]] & !identical(genes[["limits"]], c(NA, NA))) {
    warning(paste('Genes\' limits were specified, setting genes[["share.limits"]] = TRUE.'))
  }
  # Check merged.
  if(!is.null(merged)){
    if(length(merged) != 1 | !(merged %in% c("direct", "indirect"))){
      stop(paste('Incorrect merged value: It must be either NULL, "direct" or "indirect".'))
    }
    if(length(features) != 2 | length(signatures[["values"]]) != 2){
      stop('When merged != NULL, the number of signatures must be exactly 2.')
    }else if (any(!(features %in% sigs))){
      stop(paste('The merged features must be signatures. For blending genes, please use blend = TRUE.'))
    }
    merged.symbol <- ifelse(test = merged == "direct", yes = " + ", no = " - ")
    merged.sigs <- paste0(sigs, collapse = merged.symbol)
    rm.NA <- FALSE
  }else{
    merged.symbol <- ""
    merged.sigs <- sigs
    rm.NA <- TRUE
  }
  # Check blend.
  if(length(blend) != 1 | !is.logical(blend)){
    stop('blend must be TRUE or FALSE.')
  }
  if(blend){
    if (length(features) != 2 | length(genes[["values"]]) != 2) {
      stop('When blend = TRUE, the number of genes must be exactly 2.')
    } else if (any(!(features %in% gene))) {
      stop(paste('The blended features must be genes. For merging signatures, please use merged argument.'))
    }
  }
  # Check mfrow.
  if(length(mfrow) != 2 | any(mfrow < 0) | any(mfrow%%1 != 0)){
    stop('mfrow must be a vector of two integers > 0.')
  }
  # --- Code ---
  if(blend){
    sc <- Seurat::CreateSeuratObject(bc@expr.matrix)
    sc@reductions <- reduction
    p <- Seurat::FeaturePlot(sc, features = gene, blend = TRUE, combine = FALSE)
  } else {
    # Get info about drugs (their corresponding name in bc, the preferred name
    # used by beyondcell and the MoA).
    info <- tryCatch(suppressWarnings(FindDrugs(bc, x = sigs, na.rm = rm.NA)),error = function(cond) data.frame())
    ### If we want to merge signatures, we must recompute the bc object using
    ### the added or substracted bc@normalized BCS.
    if (!is.null(merged)) {
      if (merged == "direct") {
        merged.bcscores <- colSums(bc@normalized[sigs, cells, drop = FALSE],na.rm = TRUE)
      } else if (merged == "indirect") {
        merged.bcscores <- colMinus(bc@normalized[sigs, cells, drop = FALSE],na.rm = TRUE)
      }
      bc@data <- matrix(merged.bcscores, nrow = 1, dimnames = list(merged.sigs, cells))
      bc <- suppressMessages(bcRecompute(bc, slot = "data"))
      features <- merged.sigs
    }
    ### Join scaled BCS and gene expression values for selected features.
    full.matrix <- rbind(bc@scaled[merged.sigs, cells, drop = FALSE],
                         bc@expr.matrix[gene, cells,drop = FALSE])[features, , drop = FALSE]
    ### Signature's center. If center == NULL, set center to switch.points.
    if (is.null(signatures[["center"]])) {
      center.sigs <- bc@switch.point[merged.sigs]
    } else {
      center.sigs <- setNames(rep(signatures[["center"]], times = length(merged.sigs)), merged.sigs)
    }
    ### Gene's colours.
    if((genes[["share.limits"]] | !identical(genes[["limits"]], c(NA, NA))) & !is.null(genes[["values"]])){
      if(any(na.limits.genes)){
        if (!is.null(merged)) range.genes <- pretty(full.matrix)
        else range.genes <- pretty(full.matrix[gene, ])
        genes[["limits"]][na.limits.genes] <- c(min(range.genes), max(range.genes))[na.limits.genes]
      }
      colors.genes <- ggplot2::scale_colour_gradient(low = "lightgrey", high = "blue",limits = genes[["limits"]])
    }else{
      colors.genes <- NULL
    }
    ### Create fake seurat object.
    sc <- Seurat::CreateSeuratObject(full.matrix)
    ### Add reductions.
    sc@reductions <- reduction
    ### Plot for each signtature/gene...
    p <- lapply(features, function(y){
      if(!is.null(merged)){
        ids <- unlist(strsplit(y, split = merged.symbol, fixed = TRUE))
      }else ids <- y
      ### Drug name and MoA.
      if(any(ids %in% info$IDs)){
        drug.and.MoA <- info[which(info$IDs %in% ids), c("preferred.and.sigs", "MoAs")]
        if (nrow(drug.and.MoA) > 1) {
          ### Paste both drug names and MoAs. If MoAs are the same, just print them one time
          drug.and.MoA <- as.data.frame(t(apply(drug.and.MoA, 2, function(z){
            paste0(unique(z), collapse = merged.symbol)
          })))
        }
        drug.and.MoA[, 2] <- BreakString(drug.and.MoA[, 2]) ### Format subtitle.
      }else{
        drug.and.MoA <- c(paste0(ids, collapse = merged.symbol), "")
      }
      ### Colours (depending on wether y is a signature or a gene).
      if(all(ids %in% sigs)){
        ### Binned colorscale centred around the switch point or the specified center.
        colors <- center_scale_colour_stepsn(
          full.matrix[y, ],colorscale = signatures[["colorscale"]],
          alpha = signatures[["alpha"]],na.value = signatures[["na.value"]],
          limits = signatures[["limits"]], center = center.sigs[y],
          breaks = signatures[["breaks"]]
        )
      }else{
        ### Continuous colorscale with default Seurat colours.
        colors <- colors.genes
      }
      fp <- Seurat::FeaturePlot(sc,features = gsub('_','-',y),combine = F,...)[[1]] +
        colors + ggplot2::labs(title = drug.and.MoA[1], subtitle = drug.and.MoA[2])
      return(fp)
    })
  }
  if(identical(mfrow, c(1, 1))){
    final.p <- p
  }else{
    ncol.nrow <- mfrow[1] * mfrow[2]
    final.p <- lapply(1:ceiling(length(p)/ncol.nrow), function(i) {
      start <- ((i - 1) * ncol.nrow) + 1
      end <- min(i * ncol.nrow, length(p))
      sub.p <- patchwork::wrap_plots(p[start:end], nrow = mfrow[1],ncol = mfrow[2])
      return(sub.p)
    })
  }
  return(final.p)
}

#' @title Plot a violinplot of the BCS grouped by cell cycle phase
#' @description This function draws for each signature of interest, a plot of the beyondcell scores (BCS) grouped
#' by cell cycle phase (G1, G2M or S). Note that this information must be present in \code{bc@meta.data} and
#' can be obtained using \code{\link[Seurat]{CellCycleScoring}}.
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param signatures Vector with the names of the signatures of interest. If \code{signatures = 'all'},
#' all signatures are selected.
#' @return A list of \code{\link[ggplot2]{ggplot2}} violindots, one for each signature of interst. In each violindot,
#' the BCS are grouped in G1, G2M or S phase groups.
bcCellCycle <- function(bc, signatures){
  if(class(bc) != 'beyondcell') stop('bc must be a beyondcell object.')
  if(!('Phase' %in% colnames(bc@meta.data))) stop('Cell cycle information not present in bc@meta.data')
  if(!is.character(signatures)) stop('Signatures must be a character vector.')
  if(length(signatures) == 1 & signatures[1] == 'all'){
    signatures <- rownames(bc@normalized)
    in.signatures <- rep(TRUE, times = nrow(bc@normalized))
  }else{
    in.signatures <- !is.null(signatures) & signatures %in% rownames(bc@normalized)
    if(all(!in.signatures)){
      stop('None of the specified signatures were found.')
    }else if(any(!in.signatures)){
      warning(paste0('These signatures were not found in bc: ',paste0(signatures[!in.signatures],collapse = ', '),'.'))
    }
  }
  
  cells <- subset(rownames(bc@meta.data), subset = rownames(bc@meta.data) %in% colnames(bc@normalized))
  info <- FindDrugs(bc, x = signatures[in.signatures])
  p <- lapply(signatures[in.signatures],function(x){
    sub.df <- na.omit(data.frame(bcscores = bc@data[x,cells],phase = bc@meta.data[cells,'Phase'],row.names = cells))
    if(x %in% info$IDs){
      drug.and.MoA <- info[which(info$IDs == x), c('preferred.and.sigs', 'MoAs')]
      drug.and.MoA[,2] <- BreakString(drug.and.MoA[,2])
    }else drug.and.MoA <- c(x, '')
    violindot <- ggplot(sub.df, aes(x = phase, y = bcscore, fill = phase)) +
      see::geom_violindot() + theme_minimal() + theme(legend.position = 'none') +
      labs(title = drug.and.MoA[1], subtitle = drug.and.MoA[2])
    return(violindot)
  })
  return(p)
}

#' @title Draw a 4 squares plot
#' @description This function draws a 4 square plot of the drug signatures present in a \code{\link[beyondcell]{beyondcell}}
#' object. A 4 squares plot consists in a scatter plot of the residuals' means (x axis) vs the switch points (y axis) of
#' a specified cluster (either a therapeutic cluster or a group defined by experimental condition or phenotype). 4 quadrants
#' are highlighted: the top-left nd bottom-right corners contain the drugs to which all selected cells are least/most
#' sensitive, respectively. The center quadrants show the drugs to which these cells are differentially insentitive or
#' sensitive when compared to the other clusters. 
#' x cut-offs: dirst and last deciles; y cut-offs: 0.1, 0.4, 0.6 and 0.9.
#' @param bc \code{beyondcell} object.
#' @param idents Name of the metadata column of interest.
#' @param lvl Character vector with the \code{idents}' level(s) of interest. If \code{lvl = NULL}, all levels will be plotted.
#' @param top Number of top drugs per quadrant to be labelled.
#' @param topnames Character vector with additional interesting drugs or pathways to be labeled (either their names or sig IDs).
#' @param force Force of repulsion between overlapping text labels. Default to 1.
#' @param alpha Transparency level between 1 (not transparent) and 0 (fully transparent).
#' @param pt.size Point size.
#' @param ... Other arguments passed to \code{\link[ggrepel]{geom_text_repel}}.
#' @details This function returns a list of \code{\link[ggplot2]{ggplot2}} objects, one per each \code{lvl}. Note that
#' residuals' means are different for each level while switch points are signature-specific. So, x axis will vary and
#' y axis will remain constant across all plots.
bc4Squares <- function(bc, idents, lvl = NULL, top = 3, topnames = NULL, force = 1, alpha = 0.7, pt.size = 3, ...){
  if(class(bc) != 'beyondcell') stop('bc must be a beyondcell object.')
  if(length(idents) != 1) stop('Idents must be a single metadata column.')
  if(!(idents %in% names(bc@ranks))){
    stop(paste0('$', idents, ' not found in bc@ranks.'))
  }else if(idents=='general'){
    stop('General rank can\'t be used in bc4Squares(). All residuals are 0.')
  }
  if(is.null(lvl)){
    lvl <- unique(bc@meta.data[, idents])
    in.lvl <- rep(TRUE, times = length(lvl))
  }else{
    in.lvl <- lvl %in% unique(bc@meta.data[, idents])
    if(all(!in.lvl)){
      stop(paste0('None of the specified levels were found in ', idents, '.'))
    }else if(any(!in.lvl)){
      warning(paste0('The following levels were not found in ', idents, ': ', paste0(lvl[!in.lvl],collapse = ', '),'.'))
    }
  }
  if(length(top) != 1 | top[1] %% 1 != 0 | top[1] < 0) stop('top must be a single integer >= 0.')
  if(!is.null(topnames)){
    in.topnames <- toupper(topnames) %in% drugInfo$drugs | tolower(topnames) %in% drugInfo$IDs | toupper(topnames) %in% toupper(rownames(bc@normalized))
    if(all(!in.topnames)){
      warning('None of the specified topname drugs were ound in bc.')
    }else if(any(!in.topnames)){
      warning(paste0('The following topname drugs were not found in bc: ',paste0(topnames[!in.topnames],collapse = ', '),'.'))
    }
  }else in.topnames <- NULL
  if(length(force) != 1 | force[1] < 0) stop('force must be a single number >= 0.')
  if(length(alpha) != 1 | alpha[1] < 0 | alpha[1] > 1) stop('alpha must be a single number between 0 and 1.')
  if(length(pt.size) != 1 | !is.numeric(pt.size)) stop('pt.size must be a single number.')
  
  info <- FindDrugs(bc, x = rownames(bc@scaled), na.rm = F)
  sp <- data.frame(switch.point = bc@switch.point[info$bc.names], row.names = info$bc.names)
  p4s <- lapply(lvl[in.lvl], function(l){
    res <- bc@ranks[[idents]][info$bc.names, paste0('residuals.mean.', l), drop = F]
    colnames(res) <- 'residuals.mean'
    df <- transform(merge(res, sp, by = 0), row.names = Row.names, Row.names = NULL)
    res.decil <- quantile(as.numeric(res$residuals.mean), prob = seq(from = 0, to = 1, length = 11), na.rm = T)
    # Drug annotation
    sp_lower_01 <- as.numeric(df$switch.point) < 0.1
    sp_lower_06 <- as.numeric(df$switch.point) < 0.6
    sp_higher_04 <- as.numeric(df$switch.point) > 0.4
    sp_higher_09 <- as.numeric(df$switch.point) > 0.9
    res_lower_10 <- as.numeric(df$residuals.mean) < res.decil[["10%"]]
    res_higher_90 <- as.numeric(df$residuals.mean) > res.decil[["90%"]]
    df$annotation <- rep("no", times = nrow(df))
    df$annotation[sp_lower_01 & res_higher_90] <- "TOP-HighSensitivityDrugs"
    df$annotation[sp_higher_09 & res_lower_10] <- "TOP-LowSensitivityDrugs"
    df$annotation[sp_higher_04 & sp_lower_06 & res_lower_10] <- "TOP-Differential-LowSensitivityDrugs"
    df$annotation[sp_higher_04 & sp_lower_06 & res_higher_90] <- "TOP-Differential-HighSensitivityDrugs"
    # Drug labels
    df$labels <- rep(NA, times = nrow(df))
    decreasing_order <- c('TOP-Differential-HighSensitivityDrugs', 'TOP-HighSensitivityDrugs')
    unique.annotations <- unique(df$annotation[df$annotation != 'no'])
    sel.labels <- unlist(sapply(unique.annotations, function(x){
      sub.df <- subset(df, subset = df$annotation == x)
      if(x %in% decreasing_order){
        sub.df <- sub.df[order(sub.df$residuals.mean, sub.df$switch.point, decreasing = T),]
      }else{
        sub.df <- sub.df[order(sub.df$residuals.mean, sub.df$switch.point, decreasing = F),]
      }
      return(rownames(sub.df)[1:min(top, nrow(sub.df))])
    }))
    df[sel.labels, 'labels'] <- info$preferred.and.sigs[match(sel.labels, info$bc.names)]
    # Topnames
    if(length(topnames[in.topnames]) > 0){
      topnames <- FindDrugs(bc, x = topnames[in.topnames], na.rm = F)
      df[match(topnames$bc.names, table = rownames(df)),'labels'] <- topnames$preferred.and.sigs
    }
    colors <- c("#1D61F2", "#DA0078", "orange", "#C7A2F5", "grey80", "black")
    names <- c("TOP-LowSensitivityDrugs", "TOP-HighSensitivityDrugs",
               "TOP-Differential-HighSensitivityDrugs",
               "TOP-Differential-LowSensitivityDrugs", "no", "black")
    # Circle's border colors
    df$borders <- df$annotation
    df$borderes[df$labels != ''] <- 'black'
      df <- rbind(subset(df, subset = df$borders != 'black'), subset(df, subset = df$borders == 'black'))
      p <- ggplot(df, aes(x = as.numeric(residuals.mean), y = as.numeric(switch.point), 
                          color = borders, fill = annotation)) +
        geom_point(shape = 21, alpha = alpha, size = pt.size) +
        scale_color_manual(values = setNames(colors, names)) +
        scale_fill_manual(values = setNames(colors, names), breaks = names[1:4], drop = F) + theme_classic() +
        geom_vline(xintercept = res.decil['10%'], linetype = 'dotted') +
        geom_vline(xintercept = res.decil['90%'], linetype = 'dotted') +
        geom_hline(yintercept = 0.9, linetype = 'dotted') +
        geom_hline(yintercept = 0.1, linetype = 'dotted') +
        geom_hline(yintercept = 0.4, linetype = 'dotted') +
        geom_hline(yintercept = 0.6, linetype = 'dotted') + ylim(0,1) +
        labs(title = paste(idents,'=',lvl), caption = paste0('x cut-offs: first and last deciles; 
                                                           y cut-offs: 0.1, 0.4, 0.6 and 0.9')) +
        xlab("Residual' Mean") + ylab("Switch Point") +
        ggrepel::geom_text_repel(label = df$labels, force = force, na.rm = T, ...) +
        guides(fill = guide_legend(title = 'Drug Annotation'), color = F) +
        cowplot::theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
      return(p)
  })
  return(p4s)
}

#' @title Creates a geneset object
#' @description This function creates a \code{\link[beyondcell]{geneset}} object.
#' @name GenerateGenesets
#' @param x A pre-loaded matrix, a ranked matrix or a path to a GMT file with custom gene sets.
#' @param n.genes Number of up and/or down-regulated genes used to compute each signature.
#' @param mode Whether the output \code{geneset} must contain up and/or down-regulated genes.
#' @param filters If \code{x} is a pre-loaded matrix, you can provide a list of filters to subset it.
#' You can specify which \code{drugs}, sig\code{IDs}, mechanisms of action (\code{MoAs}), \code{targets}
#' and/or \code{sources} you are interested in (cap insensitive). You can \code{\link[beyondcell]{ListFilters}}
#' to check all the available values for these filters. The signatures that pass \strong{ANY} of them are
#' included in the output.
#' @param comparison \code{"treated_vs_control"} or \code{"sensitive_vs_resistant"}.
#' @param include.pathways Logical, return \code{beyondcell}'s pre-computed signatures for functional pathways?
#' @details \code{x} can be:
#' \itemize{
#'   \item{A pre-loaded matrix:}{Either \code{PSc}, \code{SSc} pr \code{DSS}.}
#'   \item{A ranked matrix:}{A matrix with genes as rows and signatures columns that contains some type of
#'   numeric value such as t-stat or a LFC to rank the genes accordings.}
#'   \item{A path to a GMT file:}{A file that contains custom gene sets. Each gene set must have an "_UP"
#'   or "_DOWN" suffix.}
#' }
#' In addition, \code{mode} can be:
#' \itemize{
#'   \item{\code{"up"}:}{To compute the signatures using only up-regulated genes.}
#'   \item{\code{"down"}:}{To compute the signatures using only down-regulated genes.}
#'   \item{\code{c("up","down")}:}{To compute the signatures using both up and down-regulated genes.}
#' }
#' If \code{x} is path to a GMT file, \code{mode} is deprecated and the names of all gene sets must end
#' in "_UP" or "_DOWN" to indicate the \code{mode} of each one.
#' Finally, \code{comparison} can be:
#' \itemize{
#'   \item{\code{"treated_vs_control"}:}{(\code{PSc} and \code{DSS} like) When the numeric values or the
#'   gene sets in the GMT file were obtained from a comparison between drug treated and untreated cells.}
#'   \item{\code{"sensitive_vs_resistant"}:}{(\code{SSc} like) When the numeric values or the gene sets
#'   in the GMT file were obtained from a comparison between drug sensitive and resistant cells.}
#' }
#' When \code{x} is a pre-loaded matrix, \code{comparison} is set automatically.
#' @return A \code{geneset} object.
GenerateGenesets <- function(x, n.genes = 250, mode = c("up", "down"),
                             filters = list(drugs = NULL,IDs = NULL,MoAs = NULL,
                                            targets = NULL,sources = NULL),
                             comparison = NULL, include.pathways = TRUE){
  # Check if x is a pre-loaded matrix, a ranked matrix or a path to a GMT file.
  load('data/beyondcell/DSS.RData')
  load('data/beyondcell/PSc.RData')
  load('data/beyondcell/SSc.RData')
  load('data/beyondcell/drugInfo.RData')
  load('data/beyondcell/pathways.RData')
  is.D <- c(identical(x, PSc), identical(x, SSc), identical(x, DSS))
  if(any(is.D)){
    type <- "pre-loaded matrix"
    message(paste('Reading', c("PSc", "SSc", "DSS")[is.D], 'signatures...'))
    n.max <- 500
  }else if(is.matrix(x) | is.data.frame(x)){
    type <- "matrix"
    message('Reading input matrix...')
    ### Check if x is numeric.
    if(!is.numeric(x)) stop('x must be a numeric matrix.')
    x <- as.matrix(x)
    ### Check if there are missing values.
    if(sum(is.na(x)) > 0){
      warning('x contains NAs that will not be used to compute the geneset object.')
    }
    ### Check if there are duplicated values in the same column.
    dup.values <- apply(x, 2, function(y) any(duplicated(na.omit(y))))
    if(is.null(colnames(x))) colnames(x) <- 1:ncol(x)
    if(is.null(rownames(x))) rownames(x) <- 1:nrow(x)
    if(any(dup.values)){
      warning(paste0('The following columns contain duplicated values:',
                     paste0(colnames(x)[dup.values], collapse = ", "), '.'))
    }
    n.max <- max(apply(x, 2, function(z) length(na.omit(z))))
  }else{
    type <- "gmt"
    message('Reading gmt file...')
    gmt.file <- qusage::read.gmt(x)
    ### Check for duplicated gene sets.
    upper.gmt.names <- toupper(names(gmt.file))
    if(anyDuplicated(upper.gmt.names) != 0){
      duplicated.gene.sets <- unique(names(gmt.file)[duplicated(upper.gmt.names)])
      stop(paste0('The GMT file contains duplicated gene set\'s: ',
                  paste0(duplicated.gene.sets, collapse = ", "), '.'))
    }
  }
  # Check n.genes and mode.
  if(any(!(mode %in% c("up", "down")))) stop('Incorrect mode.')
  mode <- sort(unique(mode), decreasing = TRUE)
  if(type != "gmt"){
    ### Number of genes.
    if(length(n.genes) != 1 | n.genes[1]%%1 != 0 | n.genes[1] < 1){
      stop('n.genes must be a positive integer.')
    }else if(n.genes > n.max){
      stop(paste0('n.genes exceeds the maximum number of genes in signature (', n.max, ').'))
    }
  }else{
    ### Number of genes.
    if (!identical(n.genes, 250)) warning('x is a GMT file, n.genes is deprecated.')
    ### Mode in GMT files.
    n.up <- length(unique(grep(pattern = "_UP$", x = upper.gmt.names)))
    n.down <- length(unique(grep(pattern = "_DOWN$", x = upper.gmt.names)))
    if(n.up + n.down != length(names(gmt.file))){
      stop('All gene sets\' names in the GMT file must end in "_UP" or "_DOWN".')
    }else{
      if (n.up > 0 & n.down > 0) {
        if (!identical(mode, c("up", "down"))) {
          mode <- c("up", "down")
          warning(paste('The GMT file includes UP and DOWN gene sets. mode changed to c("up", "down").'))
        }
      }else if (n.up > 0){
        if(mode != "up"){
          mode <- "up"
          warning('The GMT file only includes UP gene sets. mode changed to "up".')
        }
      }else if (n.down > 0){
        if(mode != "down"){
          mode <- "down"
          warning('The GMT file only includes DOWN gene sets. mode changed to "down".')
        }
      }
    }
  }
  # Check filters and comparison.
  filters.names <- c("drugs", "IDs", "MoAs", "targets", "sources")
  selected.filters <- names(filters)
  if(any(!(selected.filters %in% filters.names))) stop('Invalid names in filters.')
  if(type != "pre-loaded matrix"){
    ### Filters.
    filters_class <- sapply(filters, is.null)
    if(any(!filters_class)){
      warning('x is not a pre-loaded matrix, filters is deprecated.')
    }
    ### Comparison.
    if(is.null(comparison)){
      stop(paste('Comparison must be either "treated_vs_control" or "control_vs_treated".'))
    }else if(length(comparison) != 1 | !(comparison[1] %in% c("treated_vs_control","control_vs_treated"))){
      stop(paste('Comparison must be either "treated_vs_control" or "control_vs_treated".'))
    }
  }else{
    ### Filters.
    filters_class <- sapply(filters, is.null) | sapply(filters, is.character)
    if(any(!filters_class)){
      stop(paste0('Incorrect value for filter\'s entry: "',
                  paste0(selected.filters[!filters_class], collapse = ", "),
                  '". You must provide a character vector.'))
    }
    selected.filters <- selected.filters[!sapply(filters, is.null)]
    ### Comparison.
    if(is.null(comparison)){
      comparison <- ifelse(is.D[2], 'control_vs_treated', 'treated_vs_control')
    }else{
      if(length(comparison) != 1 | !(comparison[1] %in% c("treated_vs_control","control_vs_treated"))){
        stop('Incorrect comparison.')
      }
      if(is.D[2] & comparison != "control_vs_treated"){
        comparison <- "control_vs_treated"
        warning('x = SSc, comparison changed to "control_vs_treated".')
      }else if(!is.D[2] & comparison != "treated_vs_control"){
        comparison <- "treated_vs_control"
        warning(paste0('x = ', c("PSc","SSc","DSS")[is.D],', comparison changed to "treated_vs_control".'))
      }
    }
  }
  # Check include.pathways.
  if(length(include.pathways) != 1 | !is.logical(include.pathways)){
    stop('include.pathways must be TRUE or FALSE.')
  }
  
  # If x is a pre-loaded matrix...
  if(type == "pre-loaded matrix"){
    ### sig IDs.
    if(is.D[1]){
      info <- subset(drugInfo, subset = drugInfo$sources == "LINCS")
    }else if(is.D[2]){
      info <- subset(drugInfo, subset = drugInfo$sources != "LINCS")
    }else if(is.D[3]){
      info <- subset(drugInfo, subset = drugInfo$IDs %in% DSS[[1]]$sig_id)
      x <- PSc # DSS is a subset of PSc
    }
    if(length(selected.filters) == 0){
      ids <- unique(info$IDs)
    }else{
      ids <- unique(unlist(lapply(selected.filters, function(y){
        tryCatch(suppressWarnings(GetIDs(values = filters[[y]], filter = y,df = info)),error = function(cond) character())
      })))
      warnings <- unlist(lapply(selected.filters, function(z){
        tryCatch(GetIDs(values = filters[[z]], filter = z, df = info),
                 error = function(cond){
                   err <- paste0(z, ": ", paste0(filters[[z]],collapse = ", "), ".\n")
                   return(err)
                 }, warning = function(cond){
                   warn <- as.character(cond)
                   warn.values <- strsplit(sapply(strsplit(warn, split = ": "),`[[`, 3),split = ", ")
                   return(paste0(z, ": ", warn.values))
                 })
      }))
      warnings <- warnings[!startsWith(warnings, prefix = "sig_")]
      if(length(ids) == 0){
        stop('Couldn\'t find signatures that matched any of the filters.')
      }else if(length(warnings) > 0){
        warning(paste('The following filters\' values yielded no results:\n',paste0("   - ", warnings, " ", collapse = "")))
      }
    }
    genes <- lapply(ids,function(sig){
      l <- list(up = x[["up"]][1:n.genes,sig], down = x[["down"]][1:n.genes, sig])
      return(l)
    })
    names(genes) <- ids
  }else if(type == "matrix"){
    genes <- apply(x, 2, function(sig){
      l <- list()
      if("up" %in% mode){
        up <- na.omit(rownames(x)[order(sig, decreasing = TRUE, na.last = NA)[1:n.genes]])
        l <- c(l, list(up = up))
      }
      if("down" %in% mode){
        down <- na.omit(rownames(x)[order(sig, decreasing = FALSE, na.last = NA)[1:n.genes]])
        l <- c(l, list(down = down))
      }
      return(l)
    })
  }else if(type == "gmt"){
    unique.gene.sets <- unique(gsub("_UP$|_DOWN$", "", names(gmt.file), ignore.case = TRUE))
    genes <- setNames(lapply(unique.gene.sets, function(sig){
      l <- list()
      if(toupper(paste0(sig, "_UP")) %in% upper.gmt.names){
        l <- c(l, list(up = gmt.file[[match(toupper(paste0(sig, "_UP")), table = upper.gmt.names)]]))
      }
      if(toupper(paste0(sig, "_DOWN")) %in% upper.gmt.names){
        l <- c(l, list(down = gmt.file[[match(toupper(paste0(sig, "_DOWN")), table = upper.gmt.names)]]))
      }
      return(l)
    }), unique.gene.sets)
  }
  if(type == "pre-loaded matrix"){
    info <- subset(info, subset = info$IDs %in% ids)
    info <- aggregate(.~ IDs, data = info, na.action = NULL, FUN = function(rw){
      paste(na.omit(unique(rw)), collapse = ", ")
    })
    info <- info[order(info$IDs, decreasing = FALSE), ]
  }else{
    info <- data.frame()
  }
  # Pathways.
  if(include.pathways){
    paths <- lapply(pathways, function(p) p[names(p)[mode %in% names(p)]])
  }else{
    paths <- list()
  }
  return(geneset(genelist = c(genes, paths), n.genes = n.genes,mode = mode, info = info, comparison = comparison))
}

#' @title Returns all the possible values for the specified filter
#' @description This function returns all the available values for \code{drugs},
#' \code{IDs}, \code{MoAs}, \code{targets} or \code{sources} filters in
#' \code{\link[beyondcell]{GenerateGenesets}} function.
#' @name ListFilters
#' @param entry Either \code{"drugs"}, \code{"IDs"}, \code{"MoAs"},
#' \code{"targets"} or \code{"sources"}.
#' @return All the possible values for the specified \code{entry}.
#' @export
ListFilters <- function(entry){
  if(entry == "drugs"){
    out <- sort(unique(drugInfo$drugs), decreasing = FALSE)
  }else if(entry == "IDs"){
    out <- sort(unique(drugInfo$IDs), decreasing = FALSE)
  }else if (entry == "MoAs"){
    out <- sort(unique(drugInfo$MoAs), decreasing = FALSE)
  }else if (entry == "targets"){
    out <- sort(unique(drugInfo$targets), decreasing = FALSE)
  }else if (entry == "sources"){
    out <- sort(unique(drugInfo$sources), decreasing = FALSE)
  }else{
    stop("Incorrect entry.")
  }
  return(out)
}

#' @title Subsets a beyondcell object
#' @description This function subsets a \code{\link[beyondcell]{beyondcell}} object based on signature
#' names, cells and/or the maximum proportion of \code{NaN} values desired in each signature and/or cell.
#' @param bc \code{beyondcell} object.
#' @param signatures Vector with the names of the signatures to subset by. If \code{signatures = NULL},
#' all signatures will be kept.
#' @param bg.signatures Vector with the names of the background signatures to subset by. If \code{bg.signatures = NULL},
#' all background signatures will be kept.
#' @param cells Vector with the names of the cells to subset by. If \code{cells = NULL}, all cells will be kept.
#' @param nan.sigs Maximum proportion of \code{NaN} values per cell in the output \code{beyondcell} object. All
#' cells with a proportion of \code{NaN > nan.cells} will be removed.
#' @details This function can subset a \code{beyondcell} object using its 5 parameters alone or in combination.
#' So, for example, if you specify \code{signatures} and \code{cells}, the resulting \code{beyondcell} object
#' (except the \code{@backround} slot) will be subsetted according to those vectors. The slot \code{@background}
#' will be only subsetted according to \code{cells}. If you want to subset it by signatures as well, you must
#' specify a value for \code{bg.signatures}.
#' 
#' On the other hand, if you specify \code{cells} and \code{nan.sigs}, the output \code{beyondcell} object will
#' keep the selected cells and those signatures with a proportion of \code{NaN} values below or equal to \code{nan.sigs}.
#' Note that \code{nan.sigs} and \code{nan.cells} arguments also subset the signatures and cells that meet the 
#' criteria in \code{@background} slot.
#' 
#' Finally, if you specify all parameters, the result will keep those signatures and cells of interest with a
#' proportion of \code{NaN} below or equal to \code{nan.sigs} and \code{nan.cells}, respectively.
bcSubset <- function(bc,signatures=NULL,bg.signatures=NULL,cells=NULL,nan.sigs=1,nan.cells=1){
  if(class(bc) != 'beyondcell') stop('bc must be a beyondcell object.')
  if(is.null(signatures)){
    pass.sigs <- rownames(bc@scaled)
  }else{
    in.signatures <- signatures %in% rownames(bc@scaled)
    if(all(!in.signatures)){
      stop('None of the specified signatures were found.')
    }else if(any(!in.signatures)){
      warning(paste0('These signatures were not found in bc: ',paste0(unique(signatures[!in.signatures]),collapse = ', '),'.'))
    }
    pass.sigs <- unique(signatures[in.signatures])
  }
  if(is.null(bg.signatures)){
    pass.sigs.bg <- rownames(bc@background)
  }else{
    in.signatures.bg <- bg.signatures %in% rownames(bc@background)
    if(all(!in.signatures.bg)){
      stop('None of the specified bg.signatures were found.')
    }else if(any(!in.signatures.bg)){
      warning(paste0('These bg.signatures were not found in bc: ',paste0(unique(bg.signatures[!in.signatures.bg]),collapse = ', '),'.'))
    }
    pass.sigs.bg <- unique(bg.signatures[in.signatures.bg])
  }
  if(is.null(cells)){
    pass.cells <- colnames(bc@scaled)
  }else{
    in.cells <- cells %in% colnames(bc@scaled)
    if(all(!in.cells)){
      stop('None of the specified cells were found.')
    }else if(any(!in.cells)){
      warning(paste0('These cells were not found in bc: ',paste0(unique(cells[!in.cells]),collapse = ', '),'.'))
    }
    pass.cells <- unique(cells[in.cells])
  }
  if(length(nan.sigs) != 1 | nan.sigs[1] < 0 | nan.sigs[1] > 1){
    stop('nan.sigs must be a positive number between 0 and 1.')
  }
  pass.nan.sigs <- rownames(bc@scaled)[apply(bc@scaled,1,function(x){
    sum(is.na(x)) <= ncol(bc@scaled) * nan.sigs
  })]
  pass.nan.sigs.bg <- rownames(bc@background)[apply(bc@background,1,function(y){sum(is.na(y)) <= ncol(bc@background)*nan.sigs})]
  if(length(nan.cells) != 1 | nan.cells[1] < 0 | nan.cells[1] > 1){
    stop('nan.cells must be a positive number between 0 and 1.')
  }
  pass.nan.cells <- colnames(bc@scaled)[apply(bc@scaled,2,function(y){sum(is.na(y)) <= nrow(bc@scaled)*nan.cells})]
  # Check regression and subset order
  reg.order <- bc@regression$order
  reg.order.bg <- bc@regression$order.background
  if(any(!(reg.order %in% c('regression','subset',''))) | reg.order[1] == '' & reg.order[2] != '' |
     length(reg.order) != 2 | (identical(reg.order[1],reg.order[2]) & reg.order[1] != '') |
     (!identical(reg.order.bg,c('','')) & !identical(reg.order,reg.order.bg))){
    warning(paste('Corrupt beyondcell object. Restoring original object before subsetting...'))
    bc <- bcRecompute(bc, slot = 'data')
    bc@background <- matrix(ncol = 0, nrow = 0)
    bc@regression <- list(order = c('subset',''), vars = NULL, order.background = rep('',times = 2))
    reg.order <- rep('', 2)
  }else{
    # If bc was previously subsetted and regressed, raise an error
    if(identical(reg.order, c('subset','regression'))){
      stop(paste('bc was previously subsetted and then regressed. Please run',
                 'bcSubset() on a new beyondcell object created with CreatebcObject(bc).'))
    }else if(identical(reg.order,rep('',2)) | identical(reg.order,c('regression',''))){
      bc@regression$order[match('',table = reg.order)] <- 'subset'
    }else if(tail(reg.order[which(reg.order != '')], n = 1) == 'subset'){
      warning('bc is an already subsetted object.')
    }
  }
  
  # Intersect pass.sig and pass.nan.sig (signatures that pass both filters.)
  final.sigs <- intersect(pass.sigs, pass.nan.sigs)
  # Intersect pass.sig.bg and pass.nan.sig.bg (bg.signatures that pass both filters).
  final.sigs.bg <- intersect(pass.sigs.bg, pass.nan.sigs.bg)
  # Intersect pass.cells and pass.nan.cells (cells that pass both filters)
  final.cells <- intersect(pass.cells, pass.nan.cells)
  if(!identical(sort(final.sigs, decreasing = F),sort(rownames(bc@normalized),decreasing = F)) |
     !identical(sort(final.sigs.bg,decreasing = F),sort(rownames(bc@background),decreasing = F)) |
     !identical(sort(final.cells,decreasing = F),sort(colnames(bc@normalized),decreasing = F)) |
     !identical(sort(final.cells,decreasing = F),sort(colnames(bc@background),decreasing = F))){
    if(length(final.sigs)==0){
      stop('No signature meet all filter criteria.')
    }else if(length(final.cells)==0){
      stop('No cell meet all filter criteria.')
    }else{
      bc@normalized <- round(bc@normalized[final.sigs, final.cells, drop = F], digits = 2)
      bc <- bcRecompute(bc, slot = 'normalized')
    }
    if(any(dim(bc@background) != 0)){
      if(length(final.sigs.bg)==0){
        stop('No background signature meet all filter criteria.')
      }else{
        bc@background <- bc@background[final.sigs.bg, final.cells, drop = F]
        bc@regression$order.background <- bc@regression$order
      }
    }
  }else{
    warning('bc was not subsetted.')
    bc@regression$order <- reg.order
  }
  return(bc)
}

#' @title Regresses out unwanted effects from BCS
#' @description This function regresses out unwanted effects from normalized beyondcell scores (BCS).
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param vars.to.regress Vector of metadata columns to regress out the BCS.
#' @return Returns a \code{beyondcell} object with regressed normalized BCS, regressed scaled BCS and
#' regressed switch points.
bcRegressOut <- function(bc, vars.to.regress){
  if(class(bc) != 'beyondcell') stop('bc must be a beyondcell object.')
  in.vars.to.regress <- !is.null(vars.to.regress) & vars.to.regress %in% colnames(bc@meta.data)
  if(all(!in.vars.to.regress)){
    stop('vars.to.regress not found.')
  }else if(any(!in.vars.to.regress)){
    stop(paste0('Some vars.to.regress not found: ',paste0(vars.to.regress[!in.vars.to.regress],collapse = ', '),'.'))
  }
  vars <- unique(vars.to.regress[in.vars.to.regress])
  reg.order <- bc@regression$order
  reg.order.bg <- bc@regression$order.background
  if(any(!(reg.order %in% c('regression','subset',''))) | reg.order[1] == '' & reg.order[2] != '' |
     length(reg.order) != 2 | (identical(reg.order[1],reg.order[2]) & reg.order[1] != '') |
     (!identical(reg.order.bg,c('','')) & !identical(reg.order,reg.order.bg))){
    warning(paste('Corrupt beyondcell object. Restoring original object before subsetting...'))
    bc <- bcRecompute(bc, slot = 'data')
    bc@background <- matrix(ncol = 0, nrow = 0)
    bc@regression <- list(order = c('regression',''),vars = vars, order.background = rep('',2))
  }else{
    if(identical(reg.order,c('regression','subset'))){
      stop(paste('bc was previously regressed and then subsetted. Please run',
                 'bcSubset() on a new beyondcell object created with CreatebcObject().'))
    }else if(identical(reg.order,rep('',2)) | identical(reg.order,c('subset',''))){
      bc@regression$order[match('',reg.order)] <- 'regression'
    }else if(tail(reg.order[which(reg.order != '')], n = 1) == 'regression'){
      warning('bc is an already regressed object.')
      vars <- unique(c(vars, bc@regression$vars))
      sigs <- rownames(bc@scaled)
      bg.sigs <- rownames(bc@background)
      cells <- colnames(bc@scaled)
      bc <- bcRecompute(bc, slot = 'data')
      bc@regression <- list(order = rep('',2), vars = NULL, order.background = rep('',2))
      if(any(dim(bc@background) != 0)){
        message('Restoring pre-regressed background matrix...')
        load('data/DSS.RData')
        gs.background <- GenerateGenesets(DSS, n.genes = bc@n.genes, mode = bc@mode, include.pathways = F)
        background <- bcScore(bc@expr.matrix,gs = gs.background, expr.thres = bc@thres)
        bc@background <- background@normalized
      }
      if('subset' %in% reg.order){
        bc <- bcSubset(bc, signatures = sigs, bg.signatures = bg.sigs, cells = cells)
      }
      bc@regression <- list(order = reg.order, order.background = reg.order)
    }
  }
  
  # Latent data
  latent.data <- bc@meta.data[colnames(bc@normalized), vars, drop = F]
  # Impute normalized BCS matrix
  message('Imputing normalized BCS...')
  if(!require(bnstruct)) install.packages('bnstruct')
  bc@normalized <- bnstruct::knn.impute(bc@normalized)
  # Limma formula
  fmla <- as.formula(object = paste('bcscore ~', paste(vars, collapse = '+')))
  # Compute regression and save it in bc@normalized
  message('Regressing scores...')
  total <- nrow(bc@normalized)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  bins <- ceiling(total/100)
  normalized.regressed <- t(
    apply(cbind(seq_len(nrow(bc@normalized)), bc@normalized), 1, function(x){
      regression.mat <- cbind(latent.data, x[-1])
      colnames(regression.mat) <- c(colnames(latent.data), 'bcscore')
      qr <- lm(fmla,data = regression.mat,qr = T,na.action = na.exclude)$qr
      resid <- qr.resid(qr = qr, y = x[-1])
      # Update the progress bar
      if(x[1] %% bins == 0 | x[1] == total){
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = x[1])
      }
      return(resid)
    })
  )
  bc@normalized <- round(normalized.regressed, digits = 2)
  # Close the progress bar
  Sys.sleep(0.1)
  close(pb)
  # Recompute the beyondcell object
  bc <- bcRecompute(bc, slot = 'normalized')
  # Add vars.to.regress to bc@regression$vars
  bc@regression$vars <- vars
  # Regress the background, if needed
  if(any(dim(bc@background) != 0)){
    message('Imputing background BCS...')
    bc@background <- bnstruct::knn.impute(bc@background)
    message('Regressing background BCS...')
    total.bg <- nrow(bc@background)
    pb.bg <- txtProgressBar(min = 0, max = total.bg, style = 3)
    bins.bg <- ceiling(total.bg/100)
    background.regressed <- t(
      apply(cbind(seq_len(nrow(bc@background)),bc@background),1,function(y){
        regression.mat <- cbind(latent.data, y[-1])
        colnames(regression.mat) <- c(colnames(latent.data), 'bcscore')
        qr <- lm(fmla,data = regression.mat,qr = T,na.action = na.exclude)$qr
        resid <- qr.resid(qr = qr, y = y[-1])
        # Update the progress bar
        if(y[1] %% bins == 0 | y[1] == total.bg){
          Sys.sleep(0.1)
          setTxtProgressBar(pb.bg, value = y[1])
        }
        return(resid)
      })
    )
    bc@background <- background.regressed
    bc@regression$order.background <- bc@regression$order
    # Close the background progress bar
    Sys.sleep(0.1)
    close(pb.bg)
  }
  return(bc)
}

#' @title Recomputes a beyondcell object
#' @description This function recompute a \code{\link[beyondcell]{beyondcell}} object using the matrix stored in the
#' slot \code{bc@data} (original scores) or \code{bc@normalized} (which can contain subsetted and/or regressed scores).
#' Columns added with \code{\link[beyondcell]{bcAddMetadata}} are preserved, except if they define therapeutic clusters.
#' Important: \code{bc@background} remains the same, while \code{bc@ranks} and \code{bc@reductions} are removed.
#' @param bc \code{beyondcell} object.
#' @param slot Score matrix to recompute the \code{beyondcell} object. Either \code{"data"} or \code{"normalized"}.
bcRecompute <- function(bc, slot = 'data'){
  if(slot == 'data'){
    bc@normalized <- bc@data
  }else if(slot == 'normalized'){
    bc@data <- bc@normalized <- round(bc@normalized, digits = 2)
  }
  # Recompute scaled BCS
  scaled <- round(t(apply(bc@normalized, 1, scales::rescale, to = c(0,1))),digits = 2)
  rownames(scaled) <- rownames(bc@normalized)
  colnames(scaled) <- colnames(bc@normalized)
  bc@scaled <- scaled
  # Recompute switch points
  bc@switch.point <- SwitchPoint(bc)
  # Remove ranks and reductions
  if(!is.null(names(bc@ranks))) message('Removing @ranks slot...')
  if(!is.null(names(bc@reductions))) message('Removing @reductions slot...')
  bc@ranks <- bc@reductions <- vector(mode = 'list', length = 0)
  # Remove therapeutic clusters from bc@meta.data
  therapeutic.clusters <- grep(pattern = 'bc_clusters_res.',x = colnames(bc@meta.data))
  if(length(therapeutic.clusters)>0){
    message('Removing therapeutic clsuters...')
    bc@meta.data <- bc@meta.data[, -c(therapeutic.clusters), drop = F]
  }
  return(bc)
}

#' @title Add new metadata to an existing beyondcell object
#' @description This function adds new metadata to an existing \code{\link[beyondcell]{beyondcell}} object.
#' @param bc \code{beyondcell} object.
#' @param metadata Matrix or dataframe with metadata to add. Rownames should be cell names and colnames should
#' not be already present in \code{bc@meta.data}.
#' @return Returns a \code{beyondcell} object with updated metadata.
bcAddMetadata <- function(bc, metadata){
  if(class(bc) != 'beyondcell') stop('bc must be a beyondcell object.')
  if(!is.matrix(metadata) & !is.data.frame(metadata)){
    stop('metadata must be a matrix or a data.frame')
  }
  # Check that the rownames of bc@meta.data and the new metadata are the same.
  if(!identical(sort(rownames(bc@meta.data),decreasing = F),sort(rownames(metadata),decreasing = F))){
    stop('metadata and bc@meta.data rownames are not the same.')
  }
  # Check that columns in metadata are different from the existing columns in bc@meta.data
  if(any(colnames(metadata) %in% colnames(bc@meta.data))){
    stop('Some metadata columns are already present in bc@meta.data slot.')
  }
  metadata <- metadata[rownames(bc@meta.data),,drop=F]
  bc@meta.data <- cbind(bc@meta.data, metadata)
  return(bc)
}

#' @title Merge two beyondcell objects
#' @description This function merges two \code{\link[beyondcell]{beyondcell}} objects obtained from the same 
#' single-cell matrix using the same \code{expr.thres}. It binds signatures, not cells.
#' @param bc1 First \code{beyondcell} object to merge
#' @param bc2 Second \code{beyondcell} object to merge
#' @return A merged \code{beyondcell} object
bcMerge <- function(bc1, bc2){
  if(class(bc1) != 'beyondcell') stop('bc1 must be a beyondcell object.')
  if(class(bc2) != 'beyondcell') stop('bc2 must be a beyondcell object.')
  if(!identical(bc1@thres, bc2@thres)){
    stop('bc objects weren\'t obtained using the same expression threshold.')
  }
  if(!identical(bc1@expr.matrix, bc2@expr.matrix) | !identical(bc1@SeuratInfo, bc2@SeuratInfo)){
    stop('bc objects weren\'t obtained from the same single-cell experiment.')
  }
  # Check for duplicated signatures
  duplicated.sigs <- intersect(rownames(bc1@data), rownames(bc2@data))
  if(length(duplicated.sigs)>0){
    identical.sigs <- sapply(duplicated.sigs,function(x){identical(bc1@data[x,], bc2@data[x,])})
    if(any(!identical.sigs)){
      stop(paste0('Duplicated signatures: ',paste0(duplicated.sigs[!identical.sigs],collapse = ', '),' without matching BCS in slot bc@data.'))
    }
  }
  if(!identical(bc1@regression, bc2@regression)){
    stop(paste('The two objects were not subsetted and/or regressed in the same order and with the same variables.'))
  }
  if(!identical(colnames(bc1@normalized),colnames(bc2@normalized))){
    stop('bc1 and bc2 do not contain the same cells.')
  }
  
  cells <- colnames(bc1@normalized)
  # Create a new beyondcell object
  bc <- beyondcell(expr.matrix = bc1@expr.matrix, SeuratInfo = bc1@SeuratInfo, regression = bc1@regression,
                   n.genes = unique(c(bc1@n.genes, bc2@n.genes)),
                   mode = unique(c(bc1@mode, bc2@mode)), thres = bc1@thres)
  # rbind scaled BCS
  bc@scaled <- unique(rbind(bc1@scaled, bc2@scaled[, cells]))[, cells]
  # rbind normalized BCS
  bc@normalized <- unique(rbind(bc1@normalized, bc2@normalized[,cells]))[,cells]
  # rbind data
  bc@data <- unique(rbind(bc1@data, bc2@data[,colnames(bc1@data)]))[,colnames(bc1@data)]
  # merge switch.points
  bc@switch.point <- c(bc1@switch.point, bc2@switch.point)[rownames(bc@scaled)]
  # merge meta.data
  bc@meta.data <- plyr::join(bc1@meta.data, bc2@meta.data)
  rownames(bc@meta.data) <- rownames(bc1@meta.data)
  # remove therapeutic clusters from bc@meta.data
  therapeutic.clusters <- grep('bc_clusters_res.', colnames(bc@meta.data))
  if(length(therapeutic.clusters) > 0){
    bc@meta.data <- bc@meta.data[, -c(therapeutic.clusters), drop = F]
  }
  # merge background
  bg <- list(bc1 = as.data.frame(bc1@background), bc2 = as.data.frame(bc2@background))
  is.empty.bg <- sapply(bg, function(x) dim(x)[2] == 0)
  if(all(is.empty.bg)){
    bc@background <- matrix(ncol = 0, nrow = 0)
  }else{
    background <- as.matrix(do.call('rbind', bg[!is.empty.bg]))
    rownames(background) <- gsub('bc[1|2]\\.','',rownames(background))
    bc@background <- background[unique(rownames(background)),]
  }
  return(bc)
}

#' @title Create a new beyondcell object
#' @description This function creates a new \code{\link[beyondcell]{beyondcell}} object.
#' @param bc \code{beyondcell} object.
#' @details This function creates a new \code{beyondcell} object by using the normalized BCS as the original
#' \code{bc@data}. Switch points are recomputed and \code{bc@regression} is restarted. The \code{bc@expr.matrix}
#' and \code{bc@meta.data} slots are subsetted to keep only those cells present in the new \code{bc@data} slot.
#' @return A new \code{beyondcell} object.
CreatebcObject <- function(bc){
  if(class(bc) != 'beyondcell') stop('bc must be a beyondcell object.')
  # make @data equal to @normalized
  bc@data <- bc@normalized
  # recompute switch points
  bc@switch.point <- SwitchPoint(bc)
  # subset those cells in @expr.matrix that have been removed from @data
  bc@expr.matrix <- bc@expr.matrix[, colnames(bc@data), drop = F]
  # subset those cells in @meta.data that have been removed from @data
  bc@meta.data <- bc@meta.data[colnames(bc@data),,drop=F]
  # restore @regression
  bc@regression <- list(order = rep('',2), vars = NULL, order.background = rep('',2))
  # restore @background
  bc@background <- matrix(ncol = 0, nrow = 0)
  return(bc)
}

#' @title Returns a vector of 5 colors
#' @description This function gets a color palette and returns a vector with 5 color values to be used in
#' \code{\link[beyondcell]{center_scale_colour_stepsn}}
#' @name get_colour_steps
#' @param colorscale Either a \code{viridis}, \code{RColorBrewer} or a custom palette of 3 colors (low, medium & and high).
#' If \code{colorsccale = NULL} (default), the function returns \code{beyondcell}'s own palette.
#' @return A vector with 5 color values.
get_color_steps <- function(colorscale = NULL){
  default <- c("#1D61F2", "#98B9FF", "#F7F7F7", "#FF9CBB", "#DA0078")
  if(is.null(colorscale)){
    colors <- default
  }else{
    guess.colors <- tryCatch(viridis::viridis(11,option = colorscale),error = function(cond) cond, warning = function(cond) cond)
    if(inherits(guess.colors,'error') | inherits(guess.colors,'warning')){
      guess.colors <- tryCatch(suppressWarnings(RColorBrewer::brewer.pal(12, colorscale)),error = function(cond) cond)
      if(inherits(guess.colors,'error')){
        guess.colors <- tryCatch(scale_colour_stepsn(colours = colorscale),error = function(cond) cond)
        if(inherits(guess.colors,'error')){
          warning('Colorscale not found. Setting default colorscale.')
          colors <- default
        }else{
          # If colorscale contains less than 3 values, set default colorscale.
          len.colors <- length(colorscale)
          if(len.colors < 3){
            warning(paste('Colorscale too short. It must contain 3 colors: high, medium and low. Setting default colorscale...'))
            colors <- default
          }else{
            # Else, construct a scale with 5 colors: the first, middle and last values in colorscale
            # and 2 intermediate colors between them (color.low and color.high, computed with colorRampPalette).
            color.middle <- colorscale[ceiling(len.colors/2)]
            color.low <- colorRampPalette(colors = c(colorscale[1],color.middle),space = 'Lab')(3)[2]
            color.high <- colorRampPalette(colors = c(colorscale[3],color.middle),space = 'Lab')(3)[2]
            colors <- c(colorscale[1],color.low,color.middle,color.high,colorscale[len.colors])
            if(len.colors > 3){
              warning(paste('Colorscale too long. It must contain 3 colors:','high medium and low. Colors chosen: ',
                            paste0(colors[c(1,3,5)],collapse = ', ')))
            }
          }
        }
      }else{
        len.guess <- length(guess.colors)
        idx.middle <- ceiling(len.guess/2)
        colors <- guess.colors[c(1,idx.middle - 1, idx.middle, idx.middle + 1, len.guess)]
      }
    }else{
      colors <- guess.colors[c(1,5,6,7,11)]
    }
  }
  return(colors)
}

#' @title Creates a centered sequential binned color gradient
#' @description This function creates a sequential binned color gradient (low-mid-high) centered around.
#' @name center_scale_colour_stepn
#' @param x A numeric vector. It can contain \code{NA}s.
#' @param colorscale A vector with 5 colors that can be obtained using \code{\link[beyondcell]{get_colour_steps}}.
#' @param alpha Transparency level between 1 (not transparent) and 0 (fully transparent).
#' @param na.value Colour to use for missing values.
#' @param limits Vector with the desired limits
#' @param center A single number indicating the centre of the \code{colorscale}. If \code{center = NULL} (default),
#' the centre is set to the middle point of \code{x}.
#' @param breaks A single number indicating the break size of the \code{colorscale}. Alternatively, it can be a vector
#' with the desired breaks (which don't have to be symmetric or equally distributed).
#' @return A centred sequential binned colour gradient that can be used to color \code{\link[ggplot2]{ggplot2}} objects.
center_scale_colour_stepsn <- function(x,colorscale,alpha=0.7,na.value='grey50',
                                       limits=c(NA,NA),center = NULL, breaks = 0.1){
  if(!is.numeric(x)) stop('x must be a numeric vector.')
  range.values <- pretty(x)
  # Check colorscale
  if(length(colorscale) != 5 | !tryCatch(is.matrix(col2rgb(colorscale)),error = function(cond) FALSE)){
    stop('colorscale must contain exactly 5 colours')
  }
  # Check alpha
  if(length(alpha) != 1 | alpha[1] < 0 | alpha[1] > 1){
    stop('alpha must be a positive number between 0 and 1.')
  }
  # Check na.value
  if(!tryCatch(is.matrix(col2rgb(na.value)),error = function(cond) FALSE)) stop('na.value is not a colour.')
  # Check limits
  if(length(limits) != 2) stop('limits must be a vector of length 2')
  na.limits <- is.na(limits)
  if(length(limits[!na.limits]) > 0 & !is.numeric(limits[!na.limits])){
    stop('limits must be numeric or NAs')
  }
  # If some limits are NAs, compute them
  if(any(na.limits)) limits[na.limits] <- c(min(range.values), max(range.values))[na.limits]
  # If limits are not sorted, sort them
  if(limits[2] < limits[1]){
    warning(paste('Upper limit is smaller than lower limit.','Sorting limits in increasing order.'))
    limits <- sort(limits, decreasing = F)
  }
  # Check center
  if(!is.null(center)){
    if(length(center) != 1 | !is.numeric(center)) stop('center must be a single number.')
    if(center < limits[1] | center > limits[2]){
      stop(paste('center =', center, 'outside of limits =',paste0(limits, collapse = ', ')))
    }
  }else{
    len.range <- length(range.values)
    # If len.range is odd, get the middle point
    if(len.range %% 2 == 1){
      center <- range.values[ceiling(len.range/2)]
    }else if(len.range %% 2 == 0){
      center <- round(sum(range.values[(len.range/2):((len.range/2)+1)])/2,digits = 2)
    }
  }
  # Check breaks
  if(!is.numeric(breaks)) stop('breaks must be numeric')
  # If breaks is a single number...
  if(length(breaks) == 1){
    if(breaks > abs(limits[1]-limits[2])){stop('breaks is bigger than the difference between limits.')}
  }else{
    if(any(breaks < limits[1]) | any(breaks > limits[2])){
      warning('Removing breaks outside the specified limits.')
      breaks <- breaks[which(breaks >= limits[1] & breaks <= limits[2])]
    }
  }
  # If breaks is not a vector
  if(length(breaks)==1){
    # Set a new center = center - breaks/2. The new center has to be at a minimum distance of breaks/2 from
    # the limits or be the upper limit itself.
    if(center < (limits[1] + (breaks/2))){
      original.center <- center
      center <- limits[1] + (breaks/2)
    }else if(center > (limits[2] - (breaks/2))){
      original.center <- center
      center <- limits[2]
    }
    center <- center - (breaks/2)
    # Compute brk.low (from the lower limit to be new center, by breaks)
    if(limits[1] < center){
      brk.low <- c(limits[1], seq(from = center, to = limits[1], by = -breaks))
      brk.low <- sort(unique(brk.low[which(brk.low >= limits[1])]),decreasing = F)
    }else{brk.low <- center}
    # compute brk.high (from the new center to the upper limit, by breaks)
    if(limits[2] > center){
      brk.high <- c(limits[2], seq(from = center, to = limits[2], by = breaks))
      brk.high <- sort(unique(brk.high[which(brk.high <= limits[2])]),decreasing = F)
    }else{brk.high <- center}
    # pseudo.center: the new center + breaks/2
    pseudo.center <- tail(brk.low, n = 1) + breaks/2
    # final breaks
    final.breaks <- brk.labels <- sort(unique(c(brk.low,pseudo.center,brk.high)),decreasing = F)
    # Remove all labels but the limits and the pseudo.center
    brk.labels[which(!(brk.labels %in% c(pseudo.center,limits)))] <- ''
    idx.pseudo.center <- which(brk.labels == pseudo.center)
    # If the original.center was at a distance < breaks/2 from the limits, modify the labels and the
    # pseudo.center.
    if(exists('original.center')){
      brk.labels[idx.pseudo.center] <- pseudo.center <- original.center
    }
    # If the pseudo.center == limit, remove the limits from the labels.
    if(pseudo.center %in% limits) brk.labels[-c(1,length(brk.labels))] <- ''
  }else{
    # Add limits to breaks
    breaks <- sort(unique(c(limits, breaks)), decreasing = F)
    pseudo.center <- center
    # Check which breaks element is the minimum value that is >= center
    idx.bigger.than.center <- which(cumsum(breaks >= center) == 1)
    # Set the new center to the previous element (or to the 1st element if idx.bigger.than.center == 1)
    center <- ifelse(idx.bigger.than.center > 1,breaks[idx.bigger.than.center - 1],breaks[idx.bigger.than.center])
    # brk.low (from the lower limit to the new center)
    brk.low <- breaks[1:which(breaks == center)]
    # brk.high (from the new center + 1 to the upper limits).
    brk.high <- breaks[(which(breaks == center) + 1):length(breaks)]
    # Final breaks and labels
    final.breaks <- brk.labels <- sort(unique(c(brk.low,pseudo.center,brk.high)),decreasing = F)
  }
  # COLOR
  # The color of the center (last break of brk.low and first break of brk.high)
  # and the pseudo.center is the same (so these three values form a single color break)
  rampcol.mid <- rep(colorscale[3], times = 3)
  # If brk.low is more than just the center, get a different colour for each break.
  if(length(brk.low) > 1){
    rampcol.low <- colorRampPalette(colors = colorscale[1:2],space='Lab')(length(brk.low)-1)
  }else rampcol.low <- character(0)
  # If brk.high is more than just the center and the limits[2], get a different colour for each break.
  if(length(brk.high) > 2){
    rampcol.high <- colorRampPalette(colors = colorscale[4:5],space = 'Lab')(length(brk.high)-1)
  }else rampcol.high <- character(0)
  # Rampcolors is the vector with the final colours
  rampcolors <- c(rampcol.low, rampcol.mid, rampcol.high)
  # Guide argument
  guide <- ggplot2::guide_colorsteps(even.steps = F,show.limits = F,title = 'Beyondcell')
  out <- ggplot2::scale_color_stepsn(
    colors = scales::alpha(rampcolors, alpha = alpha), breaks = final.breaks, labels = brk.labels,
    values = scales::rescale(final.breaks, to = c(0,1)), na.value = scales::alpha(na.value,alpha = alpha),
    limits = limits, guide = guide
  )
  return(out)
}

#' @title Breaks a string into several lines
#' @description This function breaks a string \code{x} formed by elements separated by \code{split} into lines of 
#' length \code{line.length}
#' @name BreakString
#' @param x String to be broken, formed by elements separated by \code{split}.
#' @param split Character that separates the elements of \code{x}.
#' @param line.length Length of the lines into which \code{x} will be broken.
#' @return A string with the same content as \code{x} broken in lines of length \code{line.length}
BreakString <- function(x, split = ', ',line.length = 50){
  if(length(x) != 1 | !is.character(x)) stop('x must be a single string.')
  if(length(split) != 1 | !is.character(split)) stop('split must be a single string.')
  if(length(line.length) != 1 | line.length[1] < 1 | line.length[1] %% 1 != 0){
    stop('line.length must be a single integer > 0')
  }
  
  if(nchar(x) <= line.length){
    final.x <- x
  }else{
    split.x <- unlist(strsplit(x, split = split))
    # Length of each element in x + length(split)
    n.char <- sapply(split.x, function(y) nchar(y) + length(split))
    # Line in which each element of x will be printed.
    n.line <- cumsum(n.char) %/% line.length + 1
    # Separate each element within a line with split; and eachline with split + "\n"
    final.x <- paste0(sapply(1:max(n.line),function(i){
      sub.x <- paste0(names(which(n.line == i)), collapse = split)
      return(sub.x)
    }),collapse = paste0(gsub(pattern = '^\\s+/\\s+$',replacement = '',x = split),'\n'))
  }
  return(final.x)
}

#' @title Substracts to the first element of a vector the rest of elements
#' @description This function substracts to the first element of a numeric vector \code{x[1]}
#' the result of elements of the same vector. \code{x[2:length(x)]}.
#' @name minus
#' @param x Numeric vector
#' @param na.rm Logical, should missing values (including \code{NaN} be removed?)
#' @return The result of the substraction
minus <- function(x, na.rm = F){
  if(length(x) == 1) x <- c(x, 0)
  # substract to the first element of x the result of elements
  out <- sum(x[1], na.rm = na.rm) - sum(x[2:length(x)], na.rm = na.rm)
  return(out)
}

#' @title Computes the column substraction
#' @description This function substracts to the first element of each column of a rectangular object
#' (\code{x[1, n]}) the rest of elements of the same column (\code{x[2:length(x), n]}).
#' @name colMinus
#' @param x A matrix or a dataframe
#' @param na.rm Logical, should missing values (including \code{NaN}) from rows \code{2:length(x)} be
#' omitted from the calculation?
#' @return A numeric rectangular object with the result of the substraction.
colMinus <- function(x, na.rm = F){
  # If x has a single row, append a row of zeros so we can run the next step
  if(dim(x)[1] == 1) x <- rbind(x, rep(0, times = ncol(x)))
  # subtract to the first row of x the rest of rows.
  first.row <- x[1, , drop = F]
  out <- first.row - colSums(x[2:nrow(x), , drop=F], na.rm = na.rm)
  return(out)
}

#' @title Computes the mean, the median and the sd of a vector
#' @description This function computes the mean, the median and the sd of a vector.
#' @name Mean.Med.SD
#' @param x Numeric vector.
#' @param na.rm Logical, should missing values (including NaN) be removed?
#' @return A named numeric vector with the mean, median and sd of \code{x}.
Mean.Med.SD <- function(x){
  stats.mean <- mean(x, na.rm = T)
  stats.median <- median(x, na.rm = T)
  stats.sd <- sd(x, na.rm = T)
  return(setNames(c(stats.mean, stats.median, stats.sd), c('mean', 'median', 'sd')))
}

#' @title Computes the BCS' statistics and ranks
#' @description This function computes the beyondcell scores' (BCS) statistics and ranks returned
#' by \code{\link[beyondcell]{bcRanks}}.
#' @name GetStatistics
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param signatures Vector with the names of the signatures of interest.
#' @param cells Vector with the name of the cells of interest.
#' @param pb \code{\link[utils]{txtProgressBar}}
#' @param total Number of iterations to complete the \code{pb}.
#' @param i Iteration number, used to increase the \code{pb}.
#' @param n.rows Number of signatures, used to increase the \code{pb}.
#' @param extended If \code{extended = TRUE}, this function returns the switch point, mean, median, sd,
#' variance, min, max, proportion of \code{NaN} and residuals' mean or per siganture. If \code{extended = F},
#' this function returns only the switch point, mean and residuals' mean.
#' @return A \code{data.frame} with the BCS' statistics and ranks
GetStatistics <- function(bc, signatures, cells, pb, total, i, n.rows, extended){
  # Check that bc is a beyondcell object
  if(class(bc) != 'beyondcell') stop('bc must be a beyondcell object.')
  # Check signatures
  in.signatures <- !is.null(signatures) & signatures %in% rownames(bc@normalized)
  if(all(!in.signatures)){
    stop('None of the specified signatures were found.')
  }else if(any(!in.signatures)){
    warning(paste0('These signatures were not found in bc: ',paste0(signatures[!in.signatures],collapse = ', '),'.'))
  }
  # Check cells
  in.cells <- !is.null(cells) & cells %in% colnames(bc@normalized)
  if(all(!in.cells)){
    stop('None of the specified cells were found.')
  }else if(any(!in.cells)){
    warning(paste0('These cells were not found in bc: ',paste0(cells[!in.cells],collapse = ', '),'.'))
  }
  # Check pb
  if(class(pb) != 'txtProgressBar') stop('pb must be a txtProgressBar object.')
  # Check total
  if(length(total) != 1 | total[1] <= 0 | total[1] %% 1 != 0){
    stop('total must be a positive integer.')
  }
  # Check i
  if(length(i) != 1 | i[1] <= 0 | i[1] %% 1 != 0){ stop('i must be a positive integer.') }
  # Check n.rows
  if(length(n.rows) != 1 | n.rows[1] <= 0 | n.rows[1] %% 1 != 0){
    stop('n.rows must be a positive integer.')
  }
  # Check extended
  if(length(extended) != 1 | !is.logical(extended[1])) stop('extended must be TRUE or FALSE.')
  
  # signatures & cells
  signatures <- signatures[in.signatures];cells <- cells[in.cells]
  # Bins: Number of iterations to complete each 1% of the progress bar.
  bins <- ceiling(total/100)
  # switch points per signature
  switch.p <- bc@switch.point[signatures]
  if(extended){
    # Dataframe with mean, median and sd per signature
    data <- cbind(seq_len(n.rows) + (n.rows * 6) * (i - 1), bc@data[signatures, cells])
    mean.med.sd <- as.data.frame(t(apply(data, 1, function(u){
      mms <- round(Mean.Med.SD(u[-1]), digits = 2)
      # Update the progress bar
      if(u[1] %% bins == 0){
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = u[1])
      }
      return(mms)
    })))
    # Variance per signature
    data[,1] <- data[,1] + n.rows
    variance.bcscore <- apply(data, 1, function(v){
      variance <- var(v[-1], na.rm = T)
      # Update the progress bar
      if(v[1] %% bins == 0){
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = v[1])
      }
      return(variance)
    })
    # Min normalized BCS per signature
    data[,1] <- data[,1] + n.rows
    min.bcscore <- apply(data, 1, function(w){
      min.bcs <- min(w[-1], na.rm = T)
      # Update the progress bar
      if(w[1] %% bins == 0){
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = w[1])
      }
      return(min.bcs)
    })
    # Max normalized BCS per signature
    data[,1] <- data[,1] + n.rows
    max.bcscore <- apply(data, 1, function(x){
      max.bcs <- max(x[-1], na.rm = T)
      # Update the progress bar
      if(x[1] %% bins == 0){
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = x[1])
      }
      return(max.bcs)
    })
    # NA proportion per signature
    data[,1] <- data[,1] + n.rows
    prop.na <- apply(data, 1, function(y){
      nas <- round(sum(is.na(y[-1]))/length(y[-1]), digits = 2)
      # Update the progress bar
      if(y[1] %% bins == 0){
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = y[1])
      }
      return(nas)
    })
    # Residuals
    normalized <- cbind(seq_len(n.rows) + (n.rows*i*5) + (n.rows*(i-1)),bc@normalized[signatures,cells])
    # Create dataframe
    stats <- data.frame(
      switch.point = switch.p, mean = mean.med.sd$mean, median = mean.med.sd$median,
      sd = mean.med.sd$sd, variance = variance.bcscore, min = min.bcscore,
      max = max.bcscore, prop.na = prop.na, row.names = signatures
    )
  }else{
    # Mean BCS per signature
    mean.bc <- round(rowMeans(bc@data[signatures, cells], na.rm = T),digits = 2)
    # Residuals
    normalized <- cbind(seq_len(n.rows) + (n.rows * (i-1)),bc@normalized[signatures,cells])
    # Create dataframe
    stats <- data.frame(switch.point = switch.p, mean = mean.bc, row.names = signatures)
  }
  # Residuals' mean
  resid <- apply(normalized, 1, function(z){
    res <- round(mean(z[-1],na.rm = T), digits = 2)
    # Update the progress bar
    if(z[1] %% bins == 0 | z[1] == total){
      Sys.sleep(0.1)
      setTxtProgressBar(pb, value = z[1])
    }
    return(res)
  })
  # Update dataframe
  stats <- cbind(stats,data.frame(residuals.mean = resid, row.names = signatures))
  return(stats)
}

#' @title Returns the IDs that match the specified filter's values
#' @description This function subsets \code{df} to select only the entries that match the specified 
#' \code{filter}'s \code{values} and returns the corresponding IDs.
#' @name GetIDs
#' @param values User-supplied filtering vector.
#' @param filter Column name or number to subset by.
#' @param df \code{data.frame} with drug information. It must contain, at least, two columns:
#' \code{"IDs"} and \code{filter}.
#' @return A vector with the IDs that match the \code{filter}'s values
GetIDs <- function(values, filter, df = drugInfo){
  if(length(values) < 1 | !is.character(values)) stop('values must be a character vector.')
  if(length(filter) != 1) stop('You must specify a single filter.')
  if(!is.character(filter) & !(filter %in% colnames(df))){
    stop(paste('filter = ', filter, ' is not a column of df.'))
  }
  if(is.numeric(filter) & (filter < 1 | filter > ncol(df))){
    stop(paste('filter = ', filter, ' is out of range.'))
  }
  if(class(df) != 'data.frame') stop('df must be a data.frame')
  if(!('IDs' %in% colnames(df))) stop('df must contain an "IDs" column.')
  
  upper.values <- toupper(values)
  selected <- subset(df, subset = toupper(df[[filter]]) %in% upper.values)
  if(filter == 'drugs' & 'preferred.drug.names' %in% colnames(df)){
    synonyms <- subset(df, subset = toupper(df[['preferred.drug.names']]) %in% unique(toupper(selected[['preferred.drug.names']])))
    selected <- unique(rbind(selected, synonyms))
  }
  ids <- unique(selected$IDs)
  not.found <- values[!(upper.values %in% toupper(df[[filter]]))]
  if(all(values %in% not.found)){
    stop('No sig ID was found for any of the elements in values.')
  }else if(length(not.found) > 0){
    filtername <- gsub(pattern = '"', replacement = '',x = deparse(substitute(filter)))
    warning(paste0('sig IDs were not found for ',length(not.found),' out of ',length(values), ' ',
                   filtername,': ', paste0(not.found, collapse = ', '),'.'))
  }
  return(ids)
}

#' @title Returns a dataframe with information about the input drugs
#' @description This function searches the input drugs in the pre-loaded \code{beyondcell} matrices and
#' returns a dataframe with drug information, including drug synonyms and MoAs.
#' @name FindDrugs
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param x A character vector with drug names and/or sig IDs.
#' @param na.rm Logical, should \code{x} entries with no available drug information be removed from the
#' final output?
#' @details The output \code{data.frame} has the following columns:
#' \itemize{
#'   \item{\code{original.names}}: Input drug names.
#'   \item{\code{bc.names}}: Drug names used in \code{bc}.
#'   \item{\code{preferred.drug.names}}: Standard drug names.
#'   \item{\code{drugs}}: Other drug names.
#'   \item{\code{IDs}}: Signature IDs.
#'   \item{\code{preferred.and.sigs}}: \code{preferred,drug.names} (or alternatively \code{bc.names})
#'   and \code{IDs}. Used as title in \code{beyondcell} plots.
#'   \item{\code{MoAs}}: Mechanism(s) of action.
#' }
#' @return A \code{data.frame}
FindDrugs <- function(bc, x, na.rm = T){
  if(class(bc) != 'beyondcell') stop('bc must be a beyondcell object.')
  if(!is.character(x)) stop('x must be a character vector.')
  if(length(na.rm) != 1 | !is.logical(na.rm)) stop('na.rm must be TRUE or FALSE.')
  sigs <- rownames(bc@normalized)
  load('data/beyondcell/drugInfo.RData')
  # Match x with bc signatures and get the indexes of matching elements.
  indexes <- lapply(x, function(y){
    idx <- match(toupper(y), table = toupper(sigs), nomatch = 0)
    if(idx == 0){
      idx <- unique(match(drugInfo$IDs[drugInfo$drugs == toupper(y)], table = sigs))
    }
    return(idx[!is.na(idx)])
  })
  # Original names (x) and bc names (sigs)
  df <- data.frame(original.names = unlist(sapply(seq_along(x),function(i){
    rep(x[i], times = length(indexes[[i]]))
  })), IDs = unlist(sapply(indexes, function(z) sigs[z])))
  df.not.found <- !(x %in% df$original.names)
  if(any(df.not.found)){
    empty.df <- data.frame(original.names = x[df.not.found],IDs = rep(NA,sum(df.not.found)))
    df <- rbind(df, empty.df)
  }
  # Get the names and pathways of the selected signatures
  info <- subset(drugInfo, subset = IDs %in% df$IDs)
  if(all(dim(info) != 0)){
    info <- aggregate(.~IDs, data = info, na.action = NULL,FUN = function(w){
      paste(na.omit(unique(w)), collapse = ', ')
    })
  }
  info.not.found <- !(df$IDs %in% drugInfo$IDs)
  if(any(info.not.found)){
    empty.info <- matrix(rep(NA, times = sum(info.not.found)*6), ncol = 6,
                         dimnames = list(1:sum(info.not.found),colnames(drugInfo)))
    info <- rbind(info, as.data.frame(empty.info))
  }
  # Merge df and info
  df <- unique(merge(df, info[,c('IDs','drugs','preferred.drug.names','MoAs')],by = 'IDs', all.x = T))
  # Add bc.names column and remove names that are not sig IDs from sig_id column.
  df$bc.names <- df$IDs
  df$IDs[!startsWith(df$IDs, prefix = 'sig_')] <- NA
  # Create preferred.and.sigs column: Preferred name and sig_id
  df$preferred.and.sigs <- sapply(1:nrow(df),function(j){
    return(
      ifelse(test = !is.na(df$preferred.drug.names[j]),
             yes = paste0(df$preferred.drug.names[j],paste0(' (',df$IDs[j],')')),
             no = df$bc.names[j])
    )
  })
  # Reorder df
  rows <- unlist(lapply(x,function(entry) which(df$original.names == entry)))
  cols <- c('original.names','bc.names','preferred.drug.names','drugs','IDs','preferred.and.sigs','MoAs')
  df <- df[rows, cols]
  # If na.rm = TRUE, remove rows with NAs in "preferred.drug.names" and "drugs" fields
  if(na.rm) df <- df[rowSums(is.na(df[, 3:4])) < 2,]
  return(df)
}

#' @title Geneset class
#' @description An object to represent signatures
#' @slot genelist A list of drug signatures and functional \code{pathways} (if specified) with up and/or
#' down-regulated genes.
#' @slot n.genes Argument passed to \code{\link[beyondcell]{GenerateGenesets}}. Number of up and/or 
#' down-regulated genes per signature.
#' @slot mode Argument passed to \code{GenerateGenesets}. Whether the \code{geneset} contains up and/or 
#' down-regulated genes.
#' @slot info Dataframe with drug signatures information, including sig IDs, drug names, NoAs, target genes and
#' data sources (LINCS, CTRP, GDSC or CCLE). This slot is only filled if \code{GenerateGenesets}' input is a
#' pre-loaded matrix.
#' @slot comparison Argument passed to \code{GenerateGenesets}. Either \code{"treated_vs_control"} or
#' \code{"control_vs_treated"}.
geneset <- setClass(
  'geneset',
  slots = list(genelist='list', n.genes='numeric', mode='character', info='data.frame', comparison='character')
)

#' @title Beyondcell class
#' @description An object to represent the beyondcell scores (BCS) for each cell and signature.
#' @slot scaled (Subsetted and/or regressed) scaled BCS.
#' @slot normalized (Subsetted and/or regressed) normalized BCS.
#' @slot data Original normalized BCS, without subsetting or regression.
#' @slot switch.point (Subsetted and/or regressed) scaled BCS for which the normalized score in \code{bc@data} is
#' 0 (one switch point per signature).
#' @slot ranks List of dataframes with the BCS' statistics and ranks returned by \code{\link[beyondcell]{bcRanks}}.
#' @slot expr.matrix Single-cell expression matrix used to compute the BCS.
#' @slot meta.data Dataframe that contains information about each cell (including the therapeutic clusters and
#' \code{\link[Seurat]{Seurat}}'s \code{@meta.data}).
#' @slot SeuratInfo List with information about the input \code{Seurat} object, including the \code{@reductions}.
#' @slot background (Subsetted and/or regressed) normalized BCS obtained using DSS signatures. Useful to compute
#' \code{beyondcell}'s UMAP reduction and the therapeutic clusters when the number of drug signatures is low.
#' @slot reductions A list of dimensional reductions for this object.
#' @slot regression A list with the order of subset and regression steps performed on the \code{beyondcell} object
#' and the variables used for regression.
#' @slot n.genes Argument passed to \code{\link[beyondcell]{GenerateGenesets}}. Number of up and/or down-regulated
#' genes for signature.
#' @slot mode Argument passed to \code{GenerateGenesets}. Whether the \code{geneset} contains up and/or down-regulated genes.
#' @slot thres Argument \code{expr.thres} passed to \code{\link[beyondcell]{bcScore}}. Minimum fraction of signature
#' genes that must be expressed in a cell to compute its BCS.
beyondcell <- setClass(
  'beyondcell',
  slots = list(
    scaled = 'matrix', normalized = 'matrix', data = 'matrix', switch.point = 'numeric', ranks = 'list',
    expr.matrix = 'matrix', meta.data = 'data.frame', SeuratInfo = 'list', background = 'matrix',
    reductions = 'list', regression = 'list', n.genes = 'numeric', mode = 'character', thres = 'numeric'
  )
)