#setwd("/mnt/ceph/mbp/servers/bioinformatics-platform/home/sarah.williams/projects/cell_groupings/analysis/18_bigpbmc")



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reran main_process_pure_pbmc.R locally, and saved the environment.
# 
# Get important parts from that .rdata file, and save them in someting smaller!
# setwd("~/data/projects/cell_groupings/ref_data/single-cell-3prime-paper/pbmc68k_data")
# load("checkpoint5_main_process_pure_pbmc_run.rdata")
# saveRDS(sub_idx,        "sub_idx.rds")
# saveRDS(pure_select_11, 'pure_select_11.rds')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Also, saved specific subsets of the reference datasets (all of them?)
all_pure_pbmc_data      <- readRDS("all_pure_pbmc_data.rds") # Big, but loadable.
all_pure_select_11types <- readRDS("all_pure_select_11types.rds")
pbmc68k_data            <- readRDS("pbmc68k_data.rds")


genes <-all_pure_pbmc_data$all_data[[1]]$hg19$gene_symbols
saveRDS(genes, "genes.rds")


## all_json <- all_pure_pbmc_data$all_json
## all_metrics <- all_pure_pbmc_data$all_metrics
## all_mol_info <- all_pure_pbmc_data$all_mol_info

# 
## saveRDS(all_pure_pbmc_data$all_data[[10]]$hg19, "data_sample10_hg19.rds")
## saveRDS(all_json,   "all_json.rds")
## saveRDS(all_metrics, "all_metrics.rds")
## saveRDS(all_mol_info, "all_mol_info.rds") #  This is the large part.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# In a new session, using the objects as saved above - 
# process data with celaref.
library(tidyverse)
library(Matrix)
library(celaref)
library(SingleCellExperiment)
library(HDF5Array)

genes          <- readRDS("genes.rds") 
pure_select_11 <- readRDS('pure_select_11.rds')



# Order from supplied scripts. Order matters!
pure_id<-c("CD34+","CD56+ NK","CD4+/CD45RA+/CD25- Naive T", 
           "CD4+/CD25 T Reg","CD8+/CD45RA+ Naive Cytotoxic",
           "CD4+/CD45RO+ Memory","CD8+ Cytotoxic T","CD19+ B",
           "CD4+ T Helper2","CD14+ Monocyte","Dendritic")


## Get Gene info table
gene_info_table <- tibble(GeneIndex  = seq_len(length(genes)),
                          GeneSymbol = genes)  #NB: genes isn't unique.

## Cell info table

# Make up arbitrary cell IDs
one_cluster_cell_info_table <- function(num_cells, cell_type, 
                                        cell_type_id, the_library){
   tibble (CellId   = paste0(cell_type_id,"_", seq_len(num_cells)), 
           CellType = as.character(cell_type),
           Library  = as.character(the_library))
}
cells_per_type         <- sapply(FUN=nrow, pure_select_11)
cell_info_table        <- bind_rows(mapply(
   FUN=one_cluster_cell_info_table, 
   cells_per_type, pure_id, seq_len(11), c(1:10,10), # 10 libs, 11 cell types
   SIMPLIFY = FALSE )) 
cell_info_table$CellType <- factor(cell_info_table$CellType, levels = pure_id)
#CellId Cluster Library
#<chr>  <fct>   <chr>  
#1 1_1    CD34+   1      
#2 1_2    CD34+   1     


## Big counts matrix -  SHOULD SAVE Hdt5summarisedExperiments here!
# Unlabelled, and also want transposed.
counts_matrix <-  t(pure_select_11[[1]])
for (i in 2:11) {
   counts_matrix <- cbind(counts_matrix, t(pure_select_11[[i]]))
}
colnames(counts_matrix) <- cell_info_table$CellId
rownames(counts_matrix) <- gene_info_table$GeneIndex


## ALT.

save_cell_type_hdf5 <- function(cell_type_id, the_library) {
   
   cell_type          <- pure_id[cell_type_id]
   num_cells_for_type <- nrow(pure_select_11[cell_type])
   
   cluster_cell_info_table <- tibble (
      CellId   = paste0(cell_type_id,"_", seq_len(num_cells_for_type)), 
      CellType = as.character(cell_type),
      Library  = as.character(the_library)
   )
   
   sce <- SingleCellExperiment(
      assays = list(counts = t(pure_select_11[[cell_type]]) ) ,
      colData = cluster_cell_info_table,
      rowData = gene_info_table
   )
   return(sce)
}

# NB: 10 libs, 11 cell types - last one has two!
sce.purePBMC <- cbind(mapply(
   FUN=save_cell_type_hdf5 , c(seq_len(10),10), seq_len(11)))

sce.purePBM$CellType <- factor(sce.purePBMC$CellType, levels = pure_id)

saveHDF5SummarizedExperiment(sce.purePBMC, "sce_Zheng2017_purePBMC")







#----------------
# Create a summarised experiment objct.
dataset_se.pbmc68kpure  <- load_se_from_tables(counts_matrix, 
                                               cell_info_table = cell_info_table,
                                               gene_info_table = gene_info_table,
                                               group_col_name  = 'CellType')

rowData(dataset_se.pbmc68kpure.subset)$total_count <- Matrix::rowSums(assay(dataset_se.pbmc68kpure.subset))
convert_se_gene_ids( dataset_se.pbmc68kpure.subset,  new_id='GeneSymbol', eval_col='total_count')

#GeneIndex






dataset_se.pbmc68kpure.subset <- subset_se_cells_by_group(dataset_se.pbmc68kpure)

# Use gene symbols for ids (this does remove some genes)
rowData(dataset_se.pbmc68kpure.subset)$total_count <- Matrix::rowSums(assay(dataset_se.pbmc68kpure.subset))
dataset_se.pbmc68kpure.subset <-  convert_se_gene_ids( dataset_se.pbmc68kpure.subset,  new_id='GeneSymbol', eval_col='total_count')

# Useful filtering  - remove too-small groups, unwanted cells
dataset_se.pbmc68kpure.subset <- trim_small_groups_and_low_expression_genes(dataset_se.pbmc68kpure.subset)

