library(tidyverse)
library(Matrix)
library(celaref)
library(SingleCellExperiment)
library(HDF5Array)


# Downloaded all_pure_pbmc_data.rds includes per 10x library data, 
# with gene info. And lots of other stuff. 
# Only useing gene info - data not split into cell types in this yet.
if (! file.exists("all_data.rds")) {
   all_pure_pbmc_data <- readRDS("all_pure_pbmc_data.rds") 
   all_data <- all_pure_pbmc_data$all_data
   saveRDS(all_data, "all_data.rds")
   rm(all_pure_pbmc_data)
} else {
   all_data  <- readRDS("all_data.rds")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reran main_process_pure_pbmc.R locally, 
# and saved the final environment as "main_process_pure_pbmc.rdata"
# Thne extracted counts matricies that have per 11x cell-type data 
if (! file.exists("pure_select_11.rds")) {
   all_pure_pbmc_data <- load("main_process_pure_pbmc.rdata")
   saveRDS(pure_select_11, 'pure_select_11.rds')
} else {
   pure_select_11 <- readRDS("pure_select_11.rds")
}

# Order from supplied scripts. Order matters!
pure_id<-c("CD34+","CD56+ NK","CD4+/CD45RA+/CD25- Naive T", 
           "CD4+/CD25 T Reg","CD8+/CD45RA+ Naive Cytotoxic",
           "CD4+/CD45RO+ Memory","CD8+ Cytotoxic T","CD19+ B",
           "CD4+ T Helper2","CD14+ Monocyte","Dendritic")

# Genes can be lifted from all_data (see main_process_pure_pbmc.R)
gene_info_table <- tibble(
   ID         = all_data[[1]][1]$hg19$genes,
   EnsemblID  = all_data[[1]][1]$hg19$genes,
   GeneSymbol = all_data[[1]][1]$hg19$gene_symbols
)

save_cell_type_hdf5 <- function(cell_type_id, the_library) {
   
   cell_type          <- pure_id[cell_type_id]
   num_cells_for_type <- nrow(pure_select_11[[cell_type_id]])
   
   cluster_cell_info_table <- tibble (
      cell_sample = paste0(cell_type_id,"_", seq_len(num_cells_for_type)), 
      CellType = as.character(cell_type),
      Library  = as.character(the_library)
   )
   
   counts_matrix <-  t(pure_select_11[[cell_type_id]])
   colnames(counts_matrix) <- cluster_cell_info_table$cell_sample
   rownames(counts_matrix) <- gene_info_table$ID
   
   sce <- SingleCellExperiment(
      assays = list(counts = counts_matrix) ,
      colData = cluster_cell_info_table,
      rowData = gene_info_table
   )
   return(sce)
}






# NB: 10 libs, 11 cell types - last one has two!
sce.purePBMC <- do.call(cbind, 
        mapply(FUN=save_cell_type_hdf5 ,  seq_len(11), c(seq_len(10),10))) 

# Tidy, use GeneSymbol, filter and save.
sce.purePBMC$CellType <- factor(sce.purePBMC$CellType, levels = pure_id)
sce.purePBMC$group    <- sce.purePBMC$CellType
rowData(sce.purePBMC)$total_count <- Matrix::rowSums(counts(sce.purePBMC))

sce.purePBMC <- convert_se_gene_ids( sce.purePBMC,  
                                     new_id   = 'GeneSymbol', 
                                     eval_col = 'total_count')
sce.purePBMC <- trim_small_groups_and_low_expression_genes(sce.purePBMC)
sce.purePBMC <- saveHDF5SummarizedExperiment(sce.purePBMC, 
                                             "sce_Zheng2017_purePBMC")



#--------------------------------
# Summarised experiment object created - now process for celaref reference.
# pure.purePBMC <- loadHDF5SummarizedExperiment(sce.purePBMC, 
#                                               "sce_Zheng2017_purePBMC")

# Don't need all cells, 1000 from each group should be plenty for DE!
# Because this is cell sorted data, cell proportions are meaningless.
set.seed(12)
sce.purePBMC.1000 <- subset_cells_by_group(sce.purePBMC, n.group = 1000)
de_table.Zheng2017purePBMC <- contrast_each_group_to_the_rest(
   sce.purePBMC.1000,  
   dataset_name="Zheng2017purePBMC", 
   num_cores = 4)

saveRDS(de_table.Zheng2017purePBMC,  
        file.path(output_dir, "de_table_Zheng2017purePBMC.rds"))

print('done')


