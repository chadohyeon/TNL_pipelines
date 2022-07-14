#!/usr/bin/Rscript
require(magrittr); require(slingshot); require(PseudotimeDE); require(scales); require(dplyr); require(tradeSeq); require(monocle3); require(irlba)
path="/home/dcha/02.glioblastoma_scRNAseq/rdata/slingshot_data_pca/"
args=commandArgs(trailingOnly = TRUE)
fn=args[1]
start_clus=args[2]
npcs=as.numeric(args[3])
cores=as.numeric(args[4])
subsampleNum=as.numeric(args[5])

sling_sce=readRDS(paste0(path,fn))
slingPAG=readRDS(paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/slingshot_data_pca/mc3.slingTradePT.pc20.nk6.nb.sample100.",fn))

set.seed(100)

sling_sce@colData$slingshot=NULL; sling_sce@colData$slingPseudotime_1=NULL
FQnorm = function(counts){
  rk = apply(counts,2,rank,ties.method='min')
  counts.sort = apply(counts,2,sort)
  refdist = apply(counts.sort,1,median)
  norm = apply(rk,2,function(r){ refdist[r] })
  rownames(norm) = rownames(counts)
  return(norm)
}
assays(sling_sce)$norm = log1p(FQnorm(assays(sling_sce)$counts))

#rd = prcomp(t(assays(sling_sce)$norm), scale. = FALSE)$x[,1:npcs]
rd = prcomp_irlba(t(assays(sling_sce)$norm), scale. = FALSE, n=npcs)$x; rownames(rd)= colnames(sling_sce)

reducedDims(sling_sce) = SimpleList(PCA = rd)
fit_ori = slingshot(sling_sce, reducedDim = 'PCA', clusterLabels = "celltype", start.clus=start_clus)

sling_ori_tbl=tibble(cell = colnames(sling_sce), pseudotime = rescale(colData(fit_ori)$slingPseudotime_1))

options(mc.cores = cores)
sling_index = mclapply(seq_len(subsampleNum), function(x) {
  sample(x = c(1:ncol(sling_sce)), size = 0.8*ncol(sling_sce), replace = FALSE)
})

sling_sub_tbl = mclapply(sling_index, function(x, sce) {
  sce = sce[, x]
  #rd = prcomp(t(assays(sce)$norm), scale. = FALSE)$x[,1:npcs]
  rd = prcomp_irlba(t(assays(sce)$norm), scale. = FALSE, n=npcs)$x; rownames(rd)= colnames(sce)
  reducedDims(sce) = SimpleList(PCA = rd)
  
  fit = slingshot(sce, reducedDim = 'PCA', clusterLabels = "celltype", start.clus=start_clus)
  tbl = tibble(cell = colnames(sce), pseudotime = rescale(colData(fit)$slingPseudotime_1))
  
  ## Make sure the direction is the same as the original pseudotime
  merge.tbl = left_join(tbl, sling_ori_tbl, by = "cell")
  
  if(cor(merge.tbl$pseudotime.x, merge.tbl$pseudotime.y) < 0) {
    tbl = dplyr::mutate(tbl, pseudotime = 1-pseudotime)
  }
  
  tbl
}, sce = sling_sce)

pcIterativePseudotime=colnames(sling_sce) %>% as.data.frame(.) %>% setNames("cell")
for(i in 1:length(sling_sub_tbl)){
  pcIterativePseudotime=merge(x=pcIterativePseudotime, y=as.data.frame(sling_sub_tbl[[i]]), by="cell", all.x=T)
}
pcIterativePseudotime = pcIterativePseudotime %>% tibble::column_to_rownames("cell") %>% setNames(paste0("pseudotime",1:length(sling_sub_tbl)))

avgSlingPC_pseudotime=(pcIterativePseudotime %>% rowMeans(.,na.rm = T))[colnames(sling_sce)] ### Averaged pseudotime

sling_avg_tbl=tibble(cell = colnames(sling_sce), pseudotime = avgSlingPC_pseudotime)

firstSlingPC_pseudotime=fit_ori@colData$slingPseudotime_1
  
sling.avgPT.sce = fit_ori
sling.avgPT.sce@colData$slingPseudotime_1=avgSlingPC_pseudotime
avgPT=avgSlingPC_pseudotime %>% as.matrix(.); colnames(avgPT)="Lineage1"
sling.avgPT.sce@colData$slingshot@assays@data$pseudotime=avgPT
cellWeights=sling.avgPT.sce@colData$slingshot@assays@data$weights

#fitGAM_sce = fitGAM(fit_ori,nknots=nknots)
#assoRes = associationTest(fitGAM_sce)
#assoRes$fdr_q=p.adjust(assoRes$pvalue)


slingPAG$sce=sling.avgPT.sce
slingPAG$firstPT=firstSlingPC_pseudotime
slingPAG$avgPT=avgSlingPC_pseudotime
saveRDS(slingPAG, paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/slingshot_data_pca/pt_add.mc3.slingTradePT.pc20.nk6.nb.sample100.",fn))
