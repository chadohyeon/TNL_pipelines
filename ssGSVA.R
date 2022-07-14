#!/usr/bin/Rscript

require(GSVA); require(msigdbr); require(dplyr)
setwd("/home/dcha/02.glioblastoma_scRNAseq/rdata/")

ncores=20

h_gene_sets = msigdbr(species = "Mus musculus", category = "H")
kegg_gene_sets = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
reac_gene_sets = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
#wp_gene_sets = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:WIKIPATHWAYS")

h_list=h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)
kegg_list=kegg_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)
reac_list=reac_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)
#wp_list=wp_gene_sets %>% split(x=.$gene_symbol, f=.$gs_name)

gs_list=c(h_list, kegg_list, reac_list)#, wp_list)

obj.fastmnn_finalAnnot=readRDS("/home/dcha/02.glioblastoma_scRNAseq/rdata/obj.fastmnn_finalAnnot.rds")
opcRoot.obj = obj.fastmnn_finalAnnot %>% subset(., subset=celltype%in%c("OPC","pre-CC1","CC1-A", "CC1-B1", "CC1-B2"))

opc_pag=readRDS("/home/dcha/02.glioblastoma_scRNAseq/rdata/opc_pag.rds")

gsva_obj=opcRoot.obj[["RNA"]]@data
gsvaOut = GSVA::gsva(expr = gsva_obj, gset.idx.list = gs_list,
                     min.sz=5, kcdf="Gaussian", method="gsva")

saveRDS(gsvaOut, "opc.all.gsva.3set.rds")
