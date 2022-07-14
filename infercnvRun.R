#!/usr/bin/Rscript
setwd("/home/dcha/02.glioblastoma_scRNAseq/rdata/infercnv")
args=commandArgs(trailingOnly = TRUE)
sample=args[1]
require(infercnv)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(sample, ".count.txt"),
                                    annotations_file=paste0(sample, ".celltype.txt"),
                                    delim="\t",
                                    gene_order_file="human_gene_pos.canonical.txt",
                                    ref_group_names=c("Immune"), chr_exclude=c("chrX", "chrY", "chrM"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=paste0("output_dir.",sample), 
                             cluster_by_groups=TRUE, 
                              denoise=TRUE, num_threads = 10,
                             HMM=TRUE)
saveRDS(infercnv_obj, paste0(sample,"_inferCnvObj.rds"))
