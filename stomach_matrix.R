library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(UCSCXenaTools)


# 1.- Loading data --------------------------------------------------------

# 1.1 Query con base a RNA-Seq y a 3 casos y 3 controles

tcga_rna <- GDCquery("TCGA-STAD",
                     data.category = "Transcriptome Profiling",
                     access = "open",
                     experimental.strategy = "RNA-Seq",
                     workflow.type = "STAR - Counts"
)

#GDCdownload(tcga_rna, method = "api", files.per.chunk = 5, directory = "D:/Matrix/TCGA/Stomach")

# 1.2 Prepare data for usage

tcga_os_data <- GDCprepare(
  tcga_rna,
  summarizedExperiment = TRUE,
  directory = "D:/Matrix/TCGA/Stomach"
  )

# 1.3 Count matrix

stom_matrix <- assay(tcga_os_data, "fpkm_unstrand")

# 1.4 Convert to data frmae

stom_data <- stom_matrix %>% 
  as.data.frame()

# 1.4.2 Convert the dots to slashes

colnames(stom_data) <- gsub("\\.", "-", colnames(stom_data))


# 1.5.- Eliminate duplicates ----------------------------------------------


# 1.5 Extract sample type from TCGA barcode

# 1.5.2 Select the 14th - 16th value which correspond to sample type codes https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
# this to then select only the samples that correspond to primary tumors

sample_type_full <- substr(colnames(stom_data), 14, 16)

# 1.5.3 Maintain only primary tumor samples so as to avoid duplicates

# 1.5.3.1 Keep only counts that correspond to primary tumor

stom_data <- stom_data[, sample_type_full == "01A"]

# 1.5.4 Assign to new object that will be modified to eliminate duplicates

stom_data2 <- stom_data

# 1.5.4.2 Keep the names up until the -01 so as to have it in the same nomenclature as the metadata

colnames(stom_data2) <- substr(colnames(stom_data2), 1, 15)

# 1.5.4.3 We can see that there are no duplicates

names(stom_data)[substr(colnames(stom_data), 1, 15) %in% names(stom_data2)[duplicated(names(stom_data2))]]


# 2.- Metadata ------------------------------------------------------------

# 2.1 Generate and Query

data_query <- XenaGenerate(subset = XenaDatasets == "TCGA.STAD.sampleMap/STAD_clinicalMatrix") %>% 
  XenaQuery()

# 2.2 Download 

xe_download <- XenaDownload(data_query, destdir = "D:/Matrix/TCGA/Stomach/Metadata")

# 2.3 Prepare (Load) the data

stom_clinical <- XenaPrepare(xe_download)

# 2.4 Observe that there are no control patients

table(stom_clinical$`_primary_disease`)


# 3.- Deleting duplicates and asigning ensembl as rownames ----------------

# 3.1 Deleting the version of ensembl and keeping only the full name

counts_data_duplicates <- 
  stom_data2 %>% 
  rownames_to_column("ensembl_version") %>% 
  mutate(ensembl = gsub("\\..*", "", ensembl_version)) %>% 
  column_to_rownames("ensembl_version")

# 3.2 Asign symbol as a column

counts_data_duplicates$symbol <- mapIds(
  org.Hs.eg.db,
  keys = counts_data_duplicates$ensembl, # A donde va a buscar
  column = "SYMBOL",  # Nueva columna con ese formato
  keytype = "ENSEMBL",  # Que formato va a buscar en keys
  multiVals = "first") # Que hacer si hay varios del mismo Key


# 3.2 Variance

numeric_data <- counts_data_duplicates %>%
  dplyr::select(where(is.numeric))

counts_data_duplicates$variance <- apply(numeric_data, 1, var)


# 3.2.2 Only mantain the version of the gene duplicate with higher variance

counts_data.tcga <- counts_data_duplicates %>% # Initial data
  group_by(symbol) %>% # Group by ensembl
  slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% # Order by variance and keep the highest
  ungroup() %>% 
  rownames_to_column("genes") %>% 
  dplyr::select( - variance, 
                 - ensembl,
                 - genes) %>% # Delete variance, and both of the ensembl ids columns
  filter(!is.na(symbol)) %>%  # Delete those that had a NA in symbol
  column_to_rownames("symbol")



# Connect to Ensembl human dataset

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query the gene types

gene_types <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype"),
  filters = "hgnc_symbol",
  values = rownames(counts_data.tcga),
  mart = mart
)

# Maintain only the protein coding genes in a list

protein_cod <- gene_types[gene_types$gene_biotype == "protein_coding",]

# Keep only protein coding genes in the counts

counts_data.tcga <- counts_data.tcga[rownames(counts_data.tcga) %in% protein_cod$hgnc_symbol ,]

# Transpose

tcga_transposed <- 
  t(counts_data.tcga)

# 5.2 Log2 Transform 

tcga_log <- 
  log2(tcga_transposed + 1)
