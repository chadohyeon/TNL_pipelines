#!/usr/bin/Rscript
require(SCENIC); require(Seurat); require(ggplot2); require(dplyr)
set.seed(100)
cancerCols=c("pre-CC1"= '#F8766D', "CC1"= '#D89000', "CC2"='#A3A500', "CC3"= '#39B600', "CC4"='#00BF7D', 
             "CC5"='#00BFC4',  "CC6"= '#00B0F6', "pre-CC2"='#9590FF', "CC7"= '#E76BF3', "CC8"='#FF62BC')
newCancerCols=c("pre-CC1"= '#F8766D', "pre-CC2"= '#ad180e', "CC1-A"='#cf9c63', "CC2"= '#de7b37', "CC1-B1"='#d1c630', 
                "CC1-B2"='#8FAA00')
on_cols = c("AC"=  '#edeb53', "TAC"= '#66d99c', "NB"=  '#00bc59', "Neuron"= '#07753b',"OPC"= '#85c9f2',  "COP"=  '#1095e6', "OL"='#2f3bbd', "EC"='#9b4db8', "PC"= '#e06360', "VSMC"='#f2a274', "BAM"='#f5bfd0', "MG"='#d1889f', "TAM"= '#ad4565', "Meninges"='#63110b', "Lymphocyte"="#ff3874")
newAnnotCols=c(newCancerCols, on_cols)


setwd("/home/dcha/02.glioblastoma_scRNAseq/rdata/scenic/")
obj.fastmnn_finalAnnot=readRDS("/home/dcha/02.glioblastoma_scRNAseq/rdata/obj.fastmnn_finalAnnot.rds")
opcRoot.obj = obj.fastmnn_finalAnnot %>% subset(., subset=celltype%in%c("OPC","pre-CC1","CC1-A", "CC1-B1", "CC1-B2"))

geneToUse=readRDS("/home/dcha/02.glioblastoma_scRNAseq/rdata/opc_pag.rds")
obj=opcRoot.obj

exprMat=obj@assays$RNA@counts[geneToUse,] %>% as.matrix(.)
exprMat_log=obj@assays$RNA@data[geneToUse,] %>% as.matrix(.)

cellInfo=obj@meta.data %>% dplyr::select(c(celltype,nFeature_RNA,nCount_RNA)) %>% setNames(c("CellType","nGene","nUMI"))


saveRDS(cellInfo, file="int/cellInfo.Rds")
colVars=list(); colVars[["CellType"]]=newAnnotCols
colVars$CellType = colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")

org = "mgi" # or hgnc, or dmel
dbDir = "cisTarget_databases" 
myDatasetTitle = "SCENIC on mouse GBM"
data(defaultDbNames)
dbs = defaultDbNames[[org]]
scenicOptions = initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=25) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

genesKept = geneFiltering(exprMat, scenicOptions=scenicOptions,
                          minCountsPerGene=3*.01*ncol(exprMat),
                          minSamples=ncol(exprMat)*.01)

exprMat_filtered = exprMat[genesKept, ]
exprMat_filtered_log=exprMat_log[genesKept,]

runCorrelation(exprMat_filtered, scenicOptions)
runGenie3(exprMat_filtered_log, scenicOptions, resumePreviousRun = TRUE)

scenicOptions = runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions = runSCENIC_2_createRegulons(scenicOptions)
scenicOptions = runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

aucellApp = plotTsne_AUCellApp(scenicOptions, exprMat_log)





library(GENIE3)
library(igraph)
library(RCy3)
library(Rgraphviz)

expr_data <- as.matrix(read.table("upanddowngene.csv", header = T, sep = ",",row.names = 1))
set.seed(expr_data)

weightMat <- GENIE3(expr_data)
linkList <- getLinkList(weightMat)

edge_listsi <- linkList
Gsi <- graph.data.frame(edge_listsi,directed = F)
Asi <- get.adjacency(Gsi,sparse = F,attr = "weight",type = "both")

g_arasi <- graph.adjacency(Asi,mode = "undirected",weighted = T)
g.cyto <- igraph.to.graphNEL(g_arasi)

cw = createNetworkFromGraph("net", graph=g.cyto)
displayGraph (cw)