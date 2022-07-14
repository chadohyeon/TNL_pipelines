#!/usr/bin/Rscript
require(dplyr)

runMAST_time=function(file_name, obj, filterThreshold=0.1, downsampling=TRUE){
  require(MAST); require(Seurat); require(dplyr); require(lme4)
  set.seed(100)
  if(downsampling){
    cells_downsampled=c()
    minCellCount=obj@meta.data[["timePoint"]] %>% table(.) %>% as.numeric(.) %>% min(.)
    for (i in obj@meta.data[["timePoint"]] %>% unique(.)){
      cells_downsampled=c(cells_downsampled, obj %>% subset(., subset=timePoint==i) %>% .@meta.data %>% rownames(.) %>% sample(.,size=minCellCount, replace=FALSE))}
    obj=obj %>% subset(., cells=cells_downsampled)
  }
  obj.data=obj@assays$RNA@data
  
  obj.sce= SingleCellExperiment(list(logNorm=obj.data),
                                colData=obj@meta.data)
  obj.sca=obj.sce %>% SceToSingleCellAssay(.)
  colData(obj.sca)$cngeneson = scale(colSums(assay(obj.sca)>0))
  
  obj.sca = filterLowExpressedGenes(obj.sca, threshold=filterThreshold)
  zlmCond = zlm(formula = as.formula("~timePoint+cngeneson+(1|batch)+celltype+percent.mt+percent.ribo+S.Score+G2M.Score"), sca=obj.sca,
                method="glmer", ebayes=FALSE, strictConvergence=FALSE, fitArgsD = list(nAGQ = 0))
  
  summaryCond = summary(zlmCond, doLRT="timePoint")
  summaryDt = summaryCond$datatable
  dt1 = summaryDt[contrast=="timePoint" & component=="H", .(primerid, `Pr(>Chisq)`)]
  dt2 = summaryDt[contrast=="timePoint" & component=="logFC", .(primerid, coef, ci.hi, ci.lo)]
  de_res = merge(dt1, dt2, by="primerid")
  colnames(de_res) = c("gene", "timePoint.H_p", "timePoint.logFC", 'timePoint.logFC_ci.hi', 'timePoint.logFC_ci.lo')
  de_res$timePoint.H_fdr = p.adjust(de_res$timePoint.H_p, "fdr")
  de_res=de_res %>% na.omit(.)
  if(downsampling){
    saveRDS(obj, file=gsub(".rds", ".downsampled.MAST.obj.rds", fullFn))}
  return(de_res)
}


args=commandArgs(trailingOnly = TRUE)
file_name=args[1]

path="/home/dcha/02.glioblastoma_scRNAseq/rdata/mast_data/"
fullFn=paste0(path, file_name)
obj=readRDS(fullFn)

mastRes=runMAST_time(fullFn, obj, downsampling=FALSE) %>% arrange(timePoint.H_fdr)
write.table(mastRes, gsub(".rds", ".mastRes.txt", fullFn), row.names = F, quote=F)
