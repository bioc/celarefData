# make-data.R
# Refer to celaref package vignette for more info.
library(tidyverse)
library(celaref)

###############################################################################
# Paths and config
datasets_dir    <- file.path('core_datasets', 'ref_data','datasets')
output_dir      <- 'celarefData'
dir.create(file.path(output_dir), showWarnings = FALSE)


################################################################################
# Process and Save

#-------------------------------------------------------------------------------
# 10X_pbmc4k_k7
# 10X genomics has several datasets available to download from their website,
# including the pbmc4k dataset, which contains PBMCs derived from a healthy
# individual.
# This was processed directly from the output files as per code in the vignette.
#
# pbmc4k Dataset available here:
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k

dataset_se.10X_pbmc4k_k7 <- load_dataset_10Xdata(
  dataset_path   = file.path(datasets_dir,'10X_pbmc4k'),
  dataset_genome = "GRCh38",
  clustering_set = "kmeans_7_clusters",
  id_to_use      = "GeneSymbol")
dataset_se.10X_pbmc4k_k7 <- trim_small_groups_and_low_expression_genes(dataset_se.10X_pbmc4k_k7)

#Then prepare the datasets with the within-experiment comparisons. Setting the num-cores to 7 to let each group run in parallel.
de_table.10X_pbmc4k_k7   <- contrast_each_group_to_the_rest(dataset_se.10X_pbmc4k_k7, dataset_name="10X_pbmc4k_k7", num_cores=7)

saveRDS(de_table.10X_pbmc4k_k7,  file.path(output_dir, "de_table_10X_pbmc4k_k7.rds"))

#-------------------------------------------------------------------------------
# Watkins2009PBMCs
# NB: the microarray processing step isn't long
#The ('Watkins2009') ‘HaemAtlas’ microarray dataset of purified PBMC cell types
# was downloaded as a normalised table from the 'haemosphere' website:
# http://haemosphere.org/datasets/show
# Processing for those data files is described in the vignette.
#
# CITATIONS:
# Watkins, Nicholas a, Arief Gusnanto, Bernard de Bono, Subhajyoti De, Diego Miranda-Saavedra, Debbie L Hardie, Will G J Angenent, et al. 2009. “A HaemAtlas: characterizing gene expression in differentiated human blood cells.” Blood 113 (19): e1–9. doi:10.1182/blood-2008-06-162958.
# Graaf, Carolyn A. de, Jarny Choi, Tracey M. Baldwin, Jessica E. Bolden, Kirsten A. Fairfax, Aaron J. Robinson, Christine Biben, et al. 2016. “Haemopedia: An Expression Atlas of Murine Hematopoietic Cells.” Stem Cell Reports 7 (3): 571–82. doi:10.1016/j.stemcr.2016.07.007.

library("illuminaHumanv2.db")


this_dataset_dir     <- file.path(datasets_dir,     'haemosphere_datasets','watkins')
norm_expression_file <- file.path(this_dataset_dir, "watkins_expression.txt")
samples_file         <- file.path(this_dataset_dir, "watkins_samples.txt")

norm_expression_table.full <- read.table(norm_expression_file, sep="\t", header=TRUE, quote="", comment.char="", row.names=1, check.names=FALSE)

samples_table              <- read_tsv(samples_file, col_types = cols())
samples_table$description  <- make.names( samples_table$description) # Avoid group or extra_factor names starting with numbers, for microarrays

#From the sample table, can see that this dataset includes other tissues, but as a PBMC reference, we only want to consider the peripheral blood samples. Like the other data loading functions, to remove a sample (or cell) from the analysis, it is enough to remove it from the sample table.
samples_table        <- samples_table[samples_table$tissue == "Peripheral Blood",]


# Setup id mappings
probes_with_gene_symbol_and_with_data <- intersect(keys(illuminaHumanv2SYMBOL),rownames(norm_expression_table.full))

# Get mappings - non NA
probe_to_symbol <- select(illuminaHumanv2.db, keys=rownames(norm_expression_table.full), columns=c("SYMBOL"), keytype="PROBEID")
probe_to_symbol <- unique(probe_to_symbol[! is.na(probe_to_symbol$SYMBOL),])
# no multimapping probes
genes_per_probe <- table(probe_to_symbol$PROBEID) # How many genes a probe is annotated against?
multimap_probes <- names(genes_per_probe)[genes_per_probe  > 1]
probe_to_symbol <- probe_to_symbol[!probe_to_symbol$PROBEID %in% multimap_probes, ]


convert_expression_table_ids<- function(expression_table, the_probes_table, old_id_name, new_id_name){

  the_probes_table <- the_probes_table[,c(old_id_name, new_id_name)]
  colnames(the_probes_table) <- c("old_id", "new_id")

  # Before DE, just pick the top expresed probe to represent the gene
  # Not perfect, but this is a ranking-based analysis.
  # hybridisation issues aside, would expect higher epressed probes to be more relevant to Single cell data anyway.
  probe_expression_levels <- rowSums(expression_table)
  the_probes_table$avgexpr <- probe_expression_levels[as.character(the_probes_table$old_id)]

  the_genes_table <-  the_probes_table %>%
    group_by(new_id) %>%
    top_n(1, avgexpr)

  expression_table <- expression_table[the_genes_table$old_id,]
  rownames(expression_table) <- the_genes_table$new_id

  return(expression_table)
}

# Just the most highly expressed probe foreach gene.
norm_expression_table.genes <- convert_expression_table_ids(norm_expression_table.full,
                                                            probe_to_symbol, old_id_name="PROBEID", new_id_name="SYMBOL")

de_table.Watkins2009PBMCs <- contrast_each_group_to_the_rest_for_norm_ma_with_limma(
  norm_expression_table = norm_expression_table.genes,
  sample_sheet_table    = samples_table,
  dataset_name          = "Watkins2009PBMCs",
  extra_factor_name     = 'description',
  sample_name           = "sampleId",
  group_name            = 'celltype')


saveRDS(de_table.Watkins2009PBMCs, file.path(output_dir, "de_table_Watkins2009_pbmcs.rds"))





#-------------------------------------------------------------------------------
# Zeisel2015
#
# In their paper 'Cell types in the mouse cortex and hippocampus revealed by
# single-cell RNA-seq' Zeisel et al. (2015) performed single cell RNA sequencing
# in mouse, in two tissues (sscortex and ca1hippocampus).
#
# This data was download from the link provided in the paper: http://linnarssonlab.org/cortex
#
# Specfically, both counts and cell annotations were parsed from this file:
# mRNA : https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt
# NB: The cell info/experimental details were parsed out of the first 10 lines
# of this file manually. The rest is counts data.
#
# CITATIONS:
# Zeisel, A., A. B. M. Manchado, S. Codeluppi, P. Lonnerberg, G. La Manno, A. Jureus, S. Marques, et al. 2015. “Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq.” Science 347 (6226): 1138–42. doi:10.1126/science.aaa1934.

zeisel_cell_info_file <- file.path(datasets_dir, "zeisel2015", "zeisel2015_mouse_scs_detail.tab")
zeisel_counts_file    <- file.path(datasets_dir, "zeisel2015", "zeisel2015_mouse_scs_counts.tab")

dataset_se.zeisel <- load_se_from_files(zeisel_counts_file, zeisel_cell_info_file,
                                        group_col_name = "level1class",
                                        cell_col_name  = "cell_id" )


# Subset the summarizedExperiment object into two tissue-specific objects
dataset_se.cortex <- dataset_se.zeisel[,dataset_se.zeisel$tissue == "sscortex"]
dataset_se.hippo  <- dataset_se.zeisel[,dataset_se.zeisel$tissue == "ca1hippocampus"]

# And filter them
dataset_se.cortex  <- trim_small_groups_and_low_expression_genes(dataset_se.cortex )
dataset_se.hippo   <- trim_small_groups_and_low_expression_genes(dataset_se.hippo )

de_table.zeisel.cortex <- contrast_each_group_to_the_rest(dataset_se.cortex, dataset_name="zeisel_sscortex",       num_cores=6)
de_table.zeisel.hippo  <- contrast_each_group_to_the_rest(dataset_se.hippo,  dataset_name="zeisel_ca1hippocampus", num_cores=6)


saveRDS(de_table.zeisel.cortex,  file.path(output_dir, "de_table_Zeisel2015_cortex.rds"))
saveRDS(de_table.zeisel.hippo,   file.path(output_dir, "de_table_Zeisel2015_hc.rds"))


#-------------------------------------------------------------------------------
# Farmer 2017

# Farmer et al. (2017) have published a survey of cell types in the mouse
# lacrimal gland at two developmental stages.
# Loading only P4 timepoint here.
# For the two datasets 'barcodes', 'genes' and 'matrix' files were downloaded
# from the GEO repository via accession codes GSM2671415 (E16) and GSM2671416
# (P4). These were processed as described in the vignette.
# Cell and cluster information were downloaded from the supplementary
# information (Tables S9, S10, S11 and S12) and formatted manually.
# CITATIONS:
# Farmer, D’Juan T., Sara Nathan, Jennifer K. Finley, Kevin Shengyang Yu, Elaine Emmerson, Lauren E. Byrnes, Julie B. Sneddon, Michael T. McManus, Aaron D. Tward, and Sarah M. Knox. 2017. “Defining epithelial cell dynamics and lineage relationships in the developing lacrimal gland.” Development 144 (13): 2517–28. doi:10.1242/dev.150789.

library(Matrix)
Farmer2017lacrimal_dir  <- file.path(datasets_dir, "Farmer2017_lacrimal", "GSM2671416_P4")

# Counts matrix
Farmer2017lacrimal_matrix_file   <- file.path(Farmer2017lacrimal_dir, "GSM2671416_P4_matrix.mtx")
Farmer2017lacrimal_barcodes_file <- file.path(Farmer2017lacrimal_dir, "GSM2671416_P4_barcodes.tsv")
Farmer2017lacrimal_genes_file    <- file.path(Farmer2017lacrimal_dir, "GSM2671416_P4_genes.tsv")

counts_matrix <- readMM(Farmer2017lacrimal_matrix_file)
counts_matrix <- as.matrix(counts_matrix)
storage.mode(counts_matrix) <- "integer"

genes <- read.table(Farmer2017lacrimal_genes_file,    sep="", stringsAsFactors = FALSE)[,1]
cells <- read.table(Farmer2017lacrimal_barcodes_file, sep="", stringsAsFactors = FALSE)[,1]
rownames(counts_matrix) <- genes
colnames(counts_matrix) <- cells


# Gene info table
gene_info_table.Farmer2017lacrimal <- as.data.frame(read.table(Farmer2017lacrimal_genes_file, sep="", stringsAsFactors = FALSE), stringsAsFactors = FALSE)
colnames(gene_info_table.Farmer2017lacrimal) <- c("ensemblID","GeneSymbol") # ensemblID is first, will become ID

## Cell/sample info
Farmer2017lacrimal_cells2groups_file  <- file.path(datasets_dir, "Farmer2017_lacrimal", "Farmer2017_supps", paste0("P4_cellinfo.tab"))
Farmer2017lacrimal_clusterinfo_file   <- file.path(datasets_dir, "Farmer2017_lacrimal", "Farmer2017_supps", paste0("Farmer2017_clusterinfo_P4.tab"))

# Cells to cluster number (just a number)
Farmer2017lacrimal_cells2groups_table <- read_tsv(Farmer2017lacrimal_cells2groups_file, col_types=cols())
# Cluster info - number to classification
Farmer2017lacrimal_clusterinfo_table <- read_tsv(Farmer2017lacrimal_clusterinfo_file, col_types=cols())
# Add in cluster info
Farmer2017lacrimal_cells2groups_table <- merge(x=Farmer2017lacrimal_cells2groups_table, y=Farmer2017lacrimal_clusterinfo_table, by.x="cluster", by.y="ClusterNum")

# Cell sample2group
cell_sample_2_group.Farmer2017lacrimal <- Farmer2017lacrimal_cells2groups_table[,c("Cell identity","ClusterID", "nGene", "nUMI")]
colnames(cell_sample_2_group.Farmer2017lacrimal) <- c("cell_sample", "group", "nGene", "nUMI")
# Add -1 onto each of the names, that seems to be in the counts
cell_sample_2_group.Farmer2017lacrimal$cell_sample <- paste0(cell_sample_2_group.Farmer2017lacrimal$cell_sample, "-1")

# Create a summarised experiment object.
dataset_se.P4  <- load_se_from_tables(counts_matrix,
                                      cell_info_table = cell_sample_2_group.Farmer2017lacrimal,
                                      gene_info_table = gene_info_table.Farmer2017lacrimal )

# For doco only.
# dataset_se.P4 <- dataset_se.P4[1:10,1:10]

# Change id
rowData(dataset_se.P4)$total_count <- rowSums(assay(dataset_se.P4))
dataset_se.P4  <-  convert_se_gene_ids( dataset_se.P4,  new_id='GeneSymbol', eval_col='total_count')

# filter
dataset_se.P4 <- trim_small_groups_and_low_expression_genes(dataset_se.P4)
de_table.Farmer2017lacrimalP4  <- contrast_each_group_to_the_rest(dataset_se.P4,  dataset_name="Farmer2017lacrimalP4", num_cores = 4)

saveRDS(de_table.Farmer2017lacrimalP4,  file.path(output_dir, "de_table_Farmer2017_lacrimalP4.rds"))

print('done')



#-------------------------------------------------------------------------------
# Purified PBMCs
#
#Zheng2017
#Massively parallel digital transcriptional profiling of single cells

# Data and scripts obtained from 
# https://github.com/10XGenomics/single-cell-3prime-paper/tree/master/pbmc68k_analysis
# 
# Note that the cellcluster labels were obtained by re-running the 
# scripts (specifically 'main_process_pure_pbmc.R' ) also provided by the authors at:
# https://github.com/10XGenomics/single-cell-3prime-paper/tree/master/pbmc68k_analysis
# 
# The environment from running main_process_pure_pbmc.R was saved....
#



