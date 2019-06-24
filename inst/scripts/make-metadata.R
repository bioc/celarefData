# make-metadata.R
# based on https://github.com/Bioconductor/GSE62944/blob/master/inst/scripts/make-metadata.R

data_title <- c('de_table_10X_pbmc4k_k7',
                'de_table_Watkins2009_pbmcs',
                'de_table_Zeisel2015_cortex',
                'de_table_Zeisel2015_hc',
                'de_table_Farmer2017_lacrimalP4',
                'de_table_Zheng2017purePBMC',
                'de_table_Zheng2017purePBMC_ensembl')

the_species <- c('human', 'human', 'mouse', 'mouse', 'mouse', 'human', 'human')
the_taxid   <- c(9606,     9606,    10090,   10090,   10090,   9606  ,  9606    )

the_data_provider <- c("10X",
                       "Haemosphere",
                       "linnarssonlab.org",
                       "linnarssonlab.org",
                       "GEO",
                       "10X",
                       "10X"
                       )

the_urls  <- c('https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k',
               'http://haemosphere.org/datasets/show',
               'https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt',
               'https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt',
               'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2671416',
               'https://github.com/10XGenomics/single-cell-3prime-paper/tree/master/pbmc68k_analysis',
               'https://github.com/10XGenomics/single-cell-3prime-paper/tree/master/pbmc68k_analysis'
               )

the_source_version <- c("March 2018", # from dl dates
                        "Jan 2018",
                        "17-Aug-2014",
                        "17-Aug-2014",
                        "Oct 2018",
                        "Jul 26, 2016", # script edit date.
                        "Jul 26, 2016" 
                        )

meta <- data.frame(
  Title              = data_title,
  Description        = paste(data_title, "object pre-processed by celaref. See celaref vignette for citation and details."),
  BiocVersion        = "3.8",
  Genome             = "", #irrelevant, really.
  SourceType         = "txt",
  SourceUrl          = the_urls,
  SourceVersion      = the_source_version,
  Species            = the_species,
  TaxonomyId         = the_taxid,
  Coordinate_1_based = TRUE, #irrelevant.
  DataProvider       = the_data_provider,
  Maintainer         = "Sarah Williams <sarah.williams1@monash.edu>",
  RDataClass         = "data.frame" ,
  DispatchClass      = "Rds",
  RDataPath          = paste0('celarefData/',data_title,".rds"),
  Tags               = "",
  Notes              = ""
)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
