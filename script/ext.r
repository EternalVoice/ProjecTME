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

