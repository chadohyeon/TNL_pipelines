#!/usr/bin/Rscript
runSeuratPCA=function(obj, ndim=10, regress=c(), normalization.method="LogNormalize", sct_method="glmGamPoi"){
  require(dplyr); require(Seurat)
  if (normalization.method=="LogNormalize"){
    DefaultAssay(obj)="RNA"
    obj = obj %>% NormalizeData(.) %>% FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>% ScaleData(., vars.to.regress=regress)
  }else if(normalization.method=="SCT"){
    require(glmGamPoi)
    obj = obj %>% SCTransform(., method=sct_method, vars.to.regress=regress)
  }else if(normalization.method=="None"){
    obj = obj %>% ScaleData(., vars.to.regress=regress)
  }else{
    DefaultAssay(obj)="RNA"
    obj = obj %>% NormalizeData(.) %>% FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>% ScaleData(., vars.to.regress=regress)
  }
  obj = obj %>% RunPCA(., npcs=ndim)
}


runSeurat.umap.louvain=function(obj, ndim=10, louvainResolution=1.0 , reductionName="pca", tsne=FALSE){
  require(dplyr); require(Seurat)
  obj = obj %>% RunUMAP(reduction = reductionName, dims = 1:ndim) %>% 
    FindNeighbors(reduction = reductionName, dims = 1:ndim) %>% 
    FindClusters(resolution = louvainResolution)
  
  if (tsne){
    obj = obj %>% RunTSNE(reduction = reductionName, dims = 1:ndim)
  }
  return(obj)
}


wrapperRunHarmony=function(inputData,project_name="harmonyMerged", anchorFeature="orig.ident", ndim=10, louvainResolution=1.0, regress=c(), SCT_normalization=FALSE){
  require(harmony); require(dplyr); require(Seurat)
  if (mode(inputData)=="S4"){
    print("[Seurat S4] detected")
    data_vector=SplitObject(inputData, anchorFeature)
  }else if(mode(inputData)=="list" & unique(as.character(lapply(inputData, function(x) mode(x))))=="S4"){
    print("[list of Seurat S4] detected")
  }else{
    print("Please specify input data as [Seurat S4] or [list of Seurat S4s]")
    return(NULL)
  }
  cellIDs = data_vector %>% lapply(., function(x) x@project.name) %>% as.character()
  obj1=data_vector[[1]]
  obj2=data_vector[2:length(data_vector)]
  merged.obj = merge(obj1, y = obj2, add.cell.ids = cellIDs, project = project_name)
  merged.obj[["toAnchor"]] = merged.obj[[anchorFeature]]
  
  if(SCT_normalization){
    merged.obj = runSeuratPCA(merged.obj, ndim, regress, "SCT")
    merged.obj = merged.obj %>% RunHarmony("toAnchor", plot_convergence = TRUE, assay.use="SCT")
  }else{
    merged.obj = runSeuratPCA(merged.obj, ndim, regress, "LogNormalize")
    merged.obj = merged.obj %>% RunHarmony("toAnchor", plot_convergence = TRUE, assay.use="RNA")
  }
  
  merged.obj %>% DimPlot(object = ., reduction = "pca", group.by = "toAnchor", pt.size = .1) %>% print(.)
  merged.obj %>% VlnPlot(object = ., features = "PC_1", group.by = "toAnchor", pt.size = 0) %>% print(.)
  
  merged.obj %>% DimPlot(object = ., reduction = "harmony", group.by = "toAnchor", pt.size = .1) %>% print(.)
  merged.obj %>% VlnPlot(object = ., features = "harmony_1", group.by = "toAnchor", pt.size = 0) %>% print(.)
  
  # Downstream analysis
  merged.obj = runSeurat.umap.louvain(merged.obj, ndim, louvainResolution, reductionName="harmony")
  
  merged.obj %>% DimPlot(., reduction = "umap",label = TRUE, pt.size = .1) %>% print(.)
  merged.obj %>% DimPlot(., reduction = "umap",group.by = "toAnchor",label = TRUE, pt.size = .1) %>% print(.)
  return(merged.obj)
}

wrapperRunSeuratCCA=function(inputData,project_name="CCA_Merged", anchorFeature="orig.ident", ndim=30, louvainResolution=1.0, regress=c(), SCT_normalization=FALSE, k.weights=100, memMB=20000){
  require(dplyr); require(Seurat); require(glmGamPoi)
  options(future.globals.maxSize=memMB * 1024^2)
  if (mode(inputData)=="S4"){
    print("[Seurat S4] detected")
    data_vector=SplitObject(inputData, anchorFeature)
  }else if(mode(inputData)=="list" & unique(as.character(lapply(inputData, function(x) mode(x))))=="S4"){
    print("[list of Seurat S4] detected")
  }else{
    print("Please specify input data as [Seurat S4] or [list of Seurat S4s]")
    return(NULL)
  }
  cellIDs = data_vector %>% lapply(., function(x) x@project.name) %>% as.character()
  obj1=data_vector[[1]]
  obj2=data_vector[2:length(data_vector)]
  merged.obj = merge(obj1, y = obj2, add.cell.ids = cellIDs, project = project_name)
  
  individual_list = SplitObject(merged.obj, split.by = anchorFeature)
  minCellNum=individual_list %>% lapply(., function(x) length(Cells(x))) %>% unlist(.) %>% min(.)
  k.weights=min(k.weights, minCellNum)
  
  if(SCT_normalization){
    individual_list = lapply(X = individual_list, function(x) SCTransform(x, method="glmGamPoi", vars.to.regress = regress))
    SCT_features=SelectIntegrationFeatures(object.list = individual_list, nfeatures=3000) ### Seurat suggestion
    individual_list=PrepSCTIntegration(object.list = individual_list, anchor.features = SCT_features)
    CCA_anchors = FindIntegrationAnchors(object.list = individual_list, dims = 1:ndim, normalization.method = "SCT", anchor.features = SCT_features)
    obj_integrated = IntegrateData(anchorset = CCA_anchors, dims = 1:ndim, normalization.method = "SCT", k.weight = k.weights)
    DefaultAssay(obj_integrated) = "integrated"
    obj_integrated = obj_integrated %>% RunPCA(., npcs=ndim)
  }else{
    individual_list = lapply(X = individual_list, FUN = NormalizeData)
    individual_list = lapply(X = individual_list, FUN = FindVariableFeatures)
    features = SelectIntegrationFeatures(object.list = individual_list)
    CCA_anchors = FindIntegrationAnchors(object.list = individual_list, dims = 1:ndim, anchor.features = features)
    obj_integrated = IntegrateData(anchorset = CCA_anchors, dims = 1:ndim, k.weight = k.weights)
    DefaultAssay(obj_integrated) = "integrated"
    obj_integrated = obj_integrated %>% ScaleData(.) %>% RunPCA(., npcs=ndim)
  }
  
  obj_integrated = obj_integrated %>% runSeurat.umap.louvain(., ndim=ndim, louvainResolution=louvainResolution)
  return(obj_integrated)
}

wrapperRunSeuratRPCA=function(inputData,project_name="RPCA_Merged", anchorFeature="orig.ident", ndim=30, louvainResolution=1.0, regress=c(), SCT_normalization=FALSE, k.anchors=5, memMB=20000){
  require(dplyr); require(Seurat); require(glmGamPoi)
  options(future.globals.maxSize=memMB * 1024^2)
  if (mode(inputData)=="S4"){
    print("[Seurat S4] detected")
    data_vector=SplitObject(inputData, anchorFeature)
  }else if(mode(inputData)=="list" & unique(as.character(lapply(inputData, function(x) mode(x))))=="S4"){
    print("[list of Seurat S4] detected")
  }else{
    print("Please specify input data as [Seurat S4] or [list of Seurat S4s]")
    return(NULL)
  }
  cellIDs = data_vector %>% lapply(., function(x) x@project.name) %>% as.character()
  obj1=data_vector[[1]]
  obj2=data_vector[2:length(data_vector)]
  merged.obj = merge(obj1, y = obj2, add.cell.ids = cellIDs, project = project_name)
  
  individual_list = SplitObject(merged.obj, split.by = anchorFeature)
  #minCellNum=individual_list %>% lapply(., function(x) length(Cells(x))) %>% unlist(.) %>% min(.)
  #k.weights=min(k.weights, minCellNum)
  
  if(SCT_normalization){
    individual_list = lapply(X = individual_list, function(x) runSeuratPCA(x, ndim=ndim, normalization.method="SCT", regress = regress))
    SCT_features=SelectIntegrationFeatures(object.list = individual_list, nfeatures=3000) ### Seurat suggestion
    individual_list=PrepSCTIntegration(object.list = individual_list, anchor.features = SCT_features)
    RPCA_anchors = FindIntegrationAnchors(object.list = individual_list, dims = 1:ndim, normalization.method = "SCT", anchor.features = SCT_features, k.anchor = k.anchors, reduction="rpca")
    obj_integrated = IntegrateData(anchorset = RPCA_anchors, dims = 1:ndim, normalization.method = "SCT")
    DefaultAssay(obj_integrated) = "integrated"
    obj_integrated = obj_integrated %>% RunPCA(., npcs=ndim)
  }else{
    individual_list = lapply(X = individual_list, function(x) runSeuratPCA(x, ndim=ndim, regress = regress))
    features = SelectIntegrationFeatures(object.list = individual_list)
    RPCA_anchors = FindIntegrationAnchors(object.list = individual_list, dims = 1:ndim, anchor.features = features, k.anchor = k.anchors, reduction="rpca")
    obj_integrated = IntegrateData(anchorset = RPCA_anchors, dims = 1:ndim)
    DefaultAssay(obj_integrated) = "integrated"
    obj_integrated = obj_integrated %>% ScaleData(.) %>% RunPCA(., npcs=ndim)
  }
  
  obj_integrated = obj_integrated %>% runSeurat.umap.louvain(., ndim=ndim, louvainResolution=louvainResolution)
  return(obj_integrated)
}


cca_anchor_transfer=function(query_obj, ref_obj, nDim=30, mnn=FALSE, mnn_feature="batch"){
  require(Seurat); require(dplyr)
  ref_obj=ref_obj %>% runSeuratPCA(., ndim=nDim)
  cca.anchors = FindTransferAnchors(reference = ref_obj, query = query_obj, features = VariableFeatures(object = ref_obj),
                                    reference.assay = "RNA", query.assay = "RNA", reduction = "cca")
  predictions = TransferData(anchorset = cca.anchors, refdata = ref_obj$celltype,
                             dims = 1:nDim, weight.reduction="cca")
  query_obj$celltype_prev=query_obj$celltype
  query_obj$celltype=predictions[["predicted.id"]]
  
  if(mnn){
    ref_obj = ref_obj %>% reducedMNN_seurat(., nDim=nDim, anchorFeature=mnn_feature)
    query_obj = MapQuery(anchorset=cca.anchors, reference=ref_obj, query=query_obj,
                         refdata = list(celltype = "celltype"), reference.reduction = "reducedmnn", reduction.model = "umap") 
  }else{
    ref_obj = ref_obj %>% RunUMAP(., dims=1:nDim, reduction="pca", return.model = TRUE)
    query_obj = MapQuery(anchorset=cca.anchors, reference=ref_obj, query=query_obj,
                         refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap") 
  }
  resList=list(query_obj, ref_obj)
  names(resList)=c("query", "ref")
  return(resList)
}

transferLatentSpace=function(query_obj, ref_obj, nDim=30, latentSpace="pca"){
  require(Seurat); require(dplyr)
  if(unique(dim(Loadings(ref_obj[[latentSpace]]))==0)){
    ref_obj[[latentSpace]]@feature.loadings=(ref_obj[[ref_obj[[latentSpace]]@assay.used]]@scale.data) %*% MASS::ginv(t(ref_obj[[latentSpace]]@cell.embeddings))
    colnames(ref_obj[[latentSpace]]@feature.loadings)=colnames(ref_obj[[latentSpace]]@cell.embeddings)
  }
  transfer.anchors = FindTransferAnchors(reference = ref_obj, query = query_obj, dims = 1:nDim, reference.reduction = latentSpace)
  predictions = TransferData(anchorset = transfer.anchors, refdata = ref_obj$celltype, dims = 1:nDim)
  query_obj$celltype_prev=query_obj$celltype
  query_obj$celltype=predictions[["predicted.id"]]
  resList=list(query_obj, ref_obj)
  names(resList)=c("query", "ref")
  return(resList)
}


reducedMNN_seurat=function(obj, nDim=10,  anchorFeature="batch",  clusterResolution=1.0, SEED=100, reduction="pca"){
  require(batchelor); require(Seurat); require(dplyr)
  set.seed(SEED)
  seuratPCs=obj@reductions[[reduction]]@cell.embeddings[,1:nDim]
  reducedMNN.obj = reducedMNN(seuratPCs, batch=obj@meta.data[[anchorFeature]])
  
  batchelor_reducedMNN_correctedPCs=reducedMNN.obj$corrected
  colnames(batchelor_reducedMNN_correctedPCs) = 1:ncol(batchelor_reducedMNN_correctedPCs) %>% paste0("reducedMNN_", .)
  obj[["reducedmnn"]] = CreateDimReducObject(embeddings = batchelor_reducedMNN_correctedPCs, key = "reducedMNN_", assay = DefaultAssay(obj))
  obj = obj %>% runSeurat.umap.louvain(., ndim=nDim, reductionName="reducedmnn", louvainResolution = clusterResolution)
  
  return(obj)
}

fastMNN_seurat=function(obj, nDim=10, anchorFeature="batch", clusterResolution=1.0, SEED=100, SCT_normalization=FALSE){
  require(batchelor); require(Seurat); require(dplyr); require(glmGamPoi)
  set.seed(SEED)
  if(SCT_normalization){
    obj = obj %>% SCTransform(., method = "glmGamPoi")
  }else{
    DefaultAssay(obj)="RNA"
    obj = obj %>% NormalizeData(.) %>% FindVariableFeatures(.)
  }
  chosen.hvgs=obj[[DefaultAssay(obj)]]@var.features
  seuratCount=obj[[DefaultAssay(obj)]]@counts
  seuratLogCount=obj[[DefaultAssay(obj)]]@data
  sce= SingleCellExperiment(list(counts=seuratCount, logcounts=seuratLogCount),
                            colData=obj@meta.data)
  
  sce = fastMNN(sce, batch=sce[[anchorFeature]], subset.row=chosen.hvgs, d=nDim)
  
  batchelor_fastMNN_correctedPCs=sce@int_colData$reducedDims$corrected
  colnames(batchelor_fastMNN_correctedPCs) = 1:ncol(batchelor_fastMNN_correctedPCs) %>% paste0("fastMNN_", .)
  obj[["fastmnn"]] = CreateDimReducObject(embeddings = batchelor_fastMNN_correctedPCs, key = "fastMNN_", assay = DefaultAssay(obj))
  obj = obj %>% runSeurat.umap.louvain(., ndim=nDim, reductionName="fastmnn", louvainResolution = clusterResolution)
  
  return(obj)
}



args=commandArgs(trailingOnly = TRUE)
sctNorm=args[1]
k_weights=100

if(sctNorm=="y"){
  sctNorm=TRUE
  prep="SCT"
}else{
  sctNorm=FALSE
  prep="logNorm"}

obj=readRDS(paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/obj.fastmnn_newAnnot.rds"))

rpcaRes = wrapperRunSeuratRPCA(obj, anchorFeature = "batch", SCT_normalization = sctNorm, louvainResolution = 3.0)
saveRDS(rpcaRes, paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/obj_",prep, "_rpca.rds"))
