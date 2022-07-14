#!/usr/bin/Rscript
require(cFIT); require(dplyr); require(Seurat)

args=commandArgs(trailingOnly = TRUE)
toCFIT=args[1]
nCores=args[2]

obj=readRDS("/home/dcha/02.glioblastoma_scRNAseq/rdata/merged_hart.rds")

mouseTerms=obj@meta.data %>% filter(species=="mouse") %>% .[[toCFIT]] %>% unique(.)
humanTerms=obj@meta.data %>% filter(species=="human") %>% .[[toCFIT]] %>% unique(.)

data.list = split_dataset_by_batch(X=t(as.matrix(obj@assays$RNA@counts)), 
                                   batch = obj@meta.data[[toCFIT]], 
                                   labels = obj@meta.data[["celltype"]],
                                   metadata = obj@meta.data)

cfitHVG = select_genes(data.list$X.list, ngenes=2000, verbose=F)
exprs.list = preprocess_for_integration(data.list$X.list, cfitHVG, scale.factor=10^4, scale=T, center=F)

ref.data=which(names(exprs.list) %in% mouseTerms)
query.data=which(names(exprs.list) %in% humanTerms)

int.out = CFITIntegrate_sketched(X.list=exprs.list[ref.data], r=15, subsample.prop = 0.1, max.niter=100, tol=1e-8, n.cores=nCores, seed=42, verbose=F)
cFitOut=list("int.out"=int.out, "exprs.list"=exprs.list, "data.list"=data.list)

saveRDS(cFitOut, paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/cfit_",toCFIT, ".rds"))

queryTarget=exprs.list[[query.data[1]]]
if (length(query.data)>1){
  for (i in 2:length(query.data)){
    queryTarget=rbind(queryTarget, exprs.list[[query.data[i]]])
  }
}

tf.out = CFITTransfer(Xtarget=queryTarget, Wref=int.out$W, max.niter = 100, seed=0, verbose=F, n.cores = nCores)
query_norm = rbind(do.call(rbind, int.out$H.list), tf.out$H) %*% diag(colSums(int.out$W))

est.labels = asign_labels(exprs.source=do.call(rbind, int.out$H.list), 
                          exprs.target=tf.out$H, 
                          labels.source=do.call(c, data.list$labels.list[ref.data]))

cFitOut=list("int.out"=int.out, "exprs.list"=exprs.list, "data.list"=data.list, "tf.out"=tf.out, "est.labels"=est.labels)

saveRDS(cFitOut, paste0("/home/dcha/02.glioblastoma_scRNAseq/rdata/cfit_",toCFIT, ".rds"))

#cFitOut=readRDS("/home/dcha/02.glioblastoma_scRNAseq/rdata/cfit_species.rds")
#toCFIT="batch"

#int.out=cFitOut$int.out
#exprs.list=cFitOut$exprs.list
#data.list=cFitOut$data.list
#remove(cFitOut)