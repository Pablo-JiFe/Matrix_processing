# Script to perform normalization, annotation and batch correction for PCA and cluster analysis

library(tidyverse)
library(limma)
library(GEOquery)
library(AnnotationDbi)
library(illuminaHumanv4.db)

# 1.- Load data and metadata ----------------------------------------------

# 1.1 Get supplementary files

getGEOSuppFiles("GSE138458", baseDir = "D:/Matrix")


# 1.1.2 Read tsv

raw_data_gse138458 <- read_tsv("D:/Matrix/GSE138458/GSE138458_non-normalized.txt/GSE138458_non-normalized.txt")

# 1.1.3 Id probes to rownbames

raw_data_gse138458 <- raw_data_gse138458 %>% 
  column_to_rownames("ID_REF")

# 1.1.4 Dont select detect p vals

data_nopval <- raw_data_gse138458 %>% 
  select(-contains("Detection")) %>%
  as.matrix()


# 1.2 Metadata query

gse_gse138458 <- getGEO("GSE138458", destdir = "D:/Matrix/GSE138458/Metadata/", getGPL = FALSE)

# 1.2.2 Extract metadata (pData)

metadata_gse138458 <- pData(phenoData(gse_gse138458[[1]]))

# 3.- Preprocessing metadata --------------------------------------------------

#3.1 Clean metadata

metadata_gse138458 <- metadata_gse138458 %>%
  mutate(
    case = ifelse(`case/control:ch1` == "SLE Case",
    "Case",
    "Control"))


# 4.- Preprocess data -----------------------------------------------------

# 4.1 log2 transformation

log2_matrix <- log2(data_nopval + 1)

# 4.2 Quantile Normalize across arrays

normalized_matrix <- normalizeBetweenArrays(log2_matrix, method = "quantile")


# 5.- Probe ID to symbol --------------------------------------------------

# 5.1 Change probe IDs to gene names

probe_ids <- rownames(normalized_matrix)

# 5.2 Map IDs to Gene Symbols

gene_symbols <- mapIds(
  illuminaHumanv4.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first" 
)

# 5.3 Convert to a data frame 

anno_df <- data.frame(
  ID_REF = names(gene_symbols),
  Gene_Symbol = gene_symbols,
  stringsAsFactors = FALSE
)

# 5.3.2 Add the column with symbols

annotated_data <- normalized_matrix %>%
  as.data.frame() %>% 
  rownames_to_column("ID_REF") %>%
  inner_join(anno_df, by = "ID_REF") %>% 
  dplyr::select(-ID_REF)

# 5.4 Object to calculate variance

numeric_data <- annotated_data %>%
  dplyr::select(where(is.numeric))

# 5.4.2 Calculate variane

annotated_data$variance <- apply(numeric_data, 1, var)

# 5.5 Maintain gene symbol with highest variance

counts_data <- annotated_data %>% 
  group_by(Gene_Symbol) %>% 
  slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% # Order by variance and keep the highest
  ungroup() %>% 
  filter(!is.na(Gene_Symbol) & Gene_Symbol != "") %>%
  column_to_rownames("Gene_Symbol") %>% 
  dplyr::select(-variance)

# 6. Metadata with only cases

metadata_gse138458_case <- metadata_gse138458 %>% 
  filter(case == "Case")

# Maintain only high flares

metadata_uniquegse138458_case <- metadata_gse138458_case %>% 
  filter(`sledai group:ch1` == "High") %>% 
  group_by(`subject id:ch1`) %>% 
  distinct(`subject id:ch1`, .keep_all = TRUE)

# Count data of only cases

counts_data_lupus <- counts_data[, colnames(counts_data) %in% metadata_uniquegse138458_case$description]

# CSV

write.csv(counts_data_lupus, "Matrix_res/lupus_matrix.csv")
