#Hannah Crohns Code
# last updated 6/16/22

library(Matrix)
library(readr)

C2matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976992_2C_matrix.mtx.gz")
C2genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976992_2C_genes.tsv.gz", col_names = FALSE)
C2cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976992_2C_barcodes.tsv.gz", col_names = FALSE)$X1
C2gene_ids <- C2genes$X1
rownames(C2matrix) <- C2gene_ids
colnames(C2matrix) <- C2cell_ids

P2matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976993_2P_matrix.mtx.gz")
P2genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976993_2P_genes.tsv.gz", col_names = FALSE)
P2cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976993_2P_barcodes.tsv.gz", col_names = FALSE)$X1
P2gene_ids <- P2genes$X1
rownames(P2matrix) <- P2gene_ids
colnames(P2matrix) <- P2cell_ids

C3matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976994_3C_matrix.mtx.gz")
C3genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976994_3C_genes.tsv.gz", col_names = FALSE)
C3cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976994_3C_barcodes.tsv.gz", col_names = FALSE)$X1
C3gene_ids <- P2genes$X1
rownames(C3matrix) <- C3gene_ids
colnames(C3matrix) <- C3cell_ids

P3matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976995_3P_matrix.mtx.gz")
P3genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976995_3P_genes.tsv.gz", col_names = FALSE)
P3cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976995_3P_barcodes.tsv.gz", col_names = FALSE)$X1
P3gene_ids <- P2genes$X1
rownames(P3matrix) <- P3gene_ids
colnames(P3matrix) <- P3cell_ids

C5matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976996_5C_matrix.mtx.gz")
C5genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976996_5C_genes.tsv.gz", col_names = FALSE)
C5cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976996_5C_barcodes.tsv.gz", col_names = FALSE)$X1
C5gene_ids <- C5_genes$X1
rownames(C5matrix) <- C5gene_ids
colnames(C5matrix) <- C5cell_ids

P5matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976997_5P_matrix.mtx.gz")
P5genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976997_5P_genes.tsv.gz", col_names = FALSE)
P5cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976997_5P_barcodes.tsv.gz", col_names = FALSE)$X1
P5gene_ids <- P5genes$X1
rownames(P5matrix) <- P5gene_ids
colnames(P5matrix) <- P5cell_ids

C7matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976998_7C_matrix.mtx.gz")
C7genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976998_7C_genes.tsv.gz", col_names = FALSE)
C7cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976998_7C_barcodes.tsv.gz", col_names = FALSE)$X1
C7gene_ids <- C7genes$X1
rownames(C7matrix) <- C7gene_ids
colnames(C7matrix) <- C7cell_ids

P7matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976999_7P_matrix.mtx.gz")
P7genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976999_7P_genes.tsv.gz", col_names = FALSE)
P7cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4976999_7P_barcodes.tsv.gz", col_names = FALSE)$X1
P7gene_ids <- P2genes$X1
rownames(P7matrix) <- P7gene_ids
colnames(P7matrix) <- P7cell_ids

C21matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977000_4_21C_matrix.mtx.gz")
C21genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977000_4_21C_genes.tsv.gz", col_names = FALSE)
C21cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977000_4_21C_barcodes.tsv.gz", col_names = FALSE)$X1
C21gene_ids <- C21genes$X1
rownames(C21matrix) <- C21gene_ids
colnames(C21matrix) <- C21cell_ids

B21matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977001_3_21B_matrix.mtx.gz")
B21genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977001_3_21B_genes.tsv.gz", col_names = FALSE)
B21cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977001_3_21B_barcodes.tsv.gz", col_names = FALSE)$X1
B21gene_ids <- B21genes$X1
rownames(B21matrix) <- B21gene_ids
colnames(B21matrix) <- B21cell_ids

C23matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977002_6_23C_matrix.mtx.gz")
C23genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977002_6_23C_genes.tsv.gz", col_names = FALSE)
C23cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977002_6_23C_barcodes.tsv.gz", col_names = FALSE)$X1
C23gene_ids <- C23genes$X1
rownames(C23matrix) <- C23gene_ids
colnames(C23matrix) <- C23cell_ids

B23matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977003_5_23B_matrix.mtx.gz")
B23genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977003_5_23B_genes.tsv.gz", col_names = FALSE)
B23cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977003_5_23B_barcodes.tsv.gz", col_names = FALSE)$X1
B23gene_ids <- B23genes$X1
rownames(B23matrix) <- B23gene_ids
colnames(B23matrix) <- B23cell_ids

C27matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977004_8_27C_matrix.mtx.gz")
C27genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977004_8_27C_genes.tsv.gz", col_names = FALSE)
C27cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977004_8_27C_barcodes.tsv.gz", col_names = FALSE)$X1
C27gene_ids <- C27genes$X1
rownames(C27matrix) <- C27gene_ids
colnames(C27matrix) <- C27cell_ids

B27matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977005_7_27B_matrix.mtx.gz")
B27genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977005_7_27B_genes.tsv.gz", col_names = FALSE)
B27cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977005_7_27B_barcodes.tsv.gz", col_names = FALSE)$X1
B27gene_ids <- B27genes$X1
rownames(B27matrix) <- B27gene_ids
colnames(B27matrix) <- B27cell_ids

C33matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977006_2_33C_matrix.mtx.gz")
C33genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977006_2_33C_genes.tsv.gz", col_names = FALSE)
C33cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977006_2_33C_barcodes.tsv.gz", col_names = FALSE)$X1
C33gene_ids <- C33genes$X1
rownames(C33matrix) <- C33gene_ids
colnames(C33matrix) <- C33cell_ids

B33matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977007_1_33B_matrix.mtx.gz")
B33genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977007_1_33B_genes.tsv.gz", col_names = FALSE)
B33cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_RAW/GSM4977007_1_33B_barcodes.tsv.gz", col_names = FALSE)$X1
B33gene_ids <- B33genes$X1
rownames(B33matrix) <- B33gene_ids
colnames(B33matrix) <- B33cell_ids

colnames(C2matrix) <- paste("C2_", colnames(C2matrix))
colnames(P2matrix) <- paste("P2_", colnames(P2matrix))
colnames(C3matrix) <- paste("C3_", colnames(C3matrix))
colnames(P3matrix) <- paste("P3_", colnames(P3matrix))
colnames(C5matrix) <- paste("C5_", colnames(C5matrix))
colnames(P5matrix) <- paste("P5_", colnames(P5matrix))
colnames(C7matrix) <- paste("C7_", colnames(C7matrix))
colnames(P7matrix) <- paste("P7_", colnames(P7matrix))
colnames(C21matrix) <- paste("C21_", colnames(C21matrix))
colnames(B21matrix) <- paste("B21_", colnames(B21matrix))
colnames(C23matrix) <- paste("C23_", colnames(C23matrix))
colnames(B23matrix) <- paste("B23_", colnames(B23matrix))
colnames(C27matrix) <- paste("C27_", colnames(C27matrix))
colnames(B27matrix) <- paste("B27_", colnames(B27matrix))
colnames(C33matrix) <- paste("C33_", colnames(C33matrix))
colnames(B33matrix) <- paste("B33_", colnames(B33matrix))


crohns_combined <- cbind(C2matrix, P2matrix, C3matrix, P3matrix, C5matrix, 
                         P5matrix, C7matrix, P7matrix, C21matrix, B21matrix, 
                         C23matrix, B23matrix, C27matrix, B27matrix, C33matrix, 
                         B33matrix)

#creating metadata for crohns, (sample id, sample type, disease)

crohns_meta <- matrix(nrow = 62166, ncol = 3)
rownames(crohns_meta) <- colnames(crohns_combined)

colnames(crohns_meta) <- c("Sample", "Type","Disease")

#rownames(crohns_meta) <- sub(".*_ ", "", rownames(crohns_meta))
#rownames(crohns_meta) <- sub("*..", "", rownames(crohns_meta))

crohns_meta[,1] <- sub("_.*", "", colnames(crohns_combined))

Co<-0
Bl<-0
for (row in 1:nrow(crohns_meta)) {
  print(grepl("C", crohns_meta[row,1], fixed = TRUE))
  if(grepl("C", crohns_meta[row,1], fixed = TRUE)){
    crohns_meta[row,2] <- ("Colon")
    Co<- Co+1
  }  else {crohns_meta[row,2] <- ("Blood")
    Bl <- Bl + 1
  }
}
print(paste0("blood:", Bl))
print(paste0("Colon:", Co))

CDAS<- 0
CD<-0
AS<-0
control<-0
for (row in 1:nrow(crohns_meta)) {
  if("B2" %in% crohns_meta[row,1] || "C2" %in% crohns_meta[row,1] || (grepl("33", crohns_meta[row,1], fixed = TRUE))){
    crohns_meta[row,3] <- ("CDAS")
    CDAS<-CDAS+1
  } else {
    if("B3" %in% crohns_meta[row,1] || "C3" %in% crohns_meta[row,1] || (grepl("27", crohns_meta[row,1], fixed = TRUE))){
      crohns_meta[row,3] <- ("AS")
      AS<-AS+1
    } else {
      if("B7" %in% crohns_meta[row,1] || "C7" %in% crohns_meta[row,1] || (grepl("23", crohns_meta[row,1], fixed = TRUE))){
        crohns_meta[row,3] <- ("CD")
        CD<-CD+1
      } else {crohns_meta[row,3] <- ("Control")
        control<-control+1
        }
      }
    }
  }
print(paste0("AS:", AS))
print(paste0("CD:", CD))
print(paste0("CDAS:", CDAS))
print(paste0("control:", control))

saveRDS(crohns_meta, "crohns_kuhn.rds")
saveRDS(crohns_combined, "crohns_exp_kuhn.rds")

# import metadata

metadata <- read.csv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_All.combined.metadata.csv.gz")

# add column for patient id and patient type

patient_id <- matrix(0:0, nrow = 60194, ncol = 1)
metadata <- cbind(metadata, patient_id)

#for(row in 1:nrow(metadata)) {
#create a data column for patient#
metadata[,11] <- metadata[1:60194,2] 
#}

#for(row in 1:nrow(metadata)) {
#metadata[row,11] <- metadata[row,2] }

metadata[1:60194,11] <- sub("_Control", "", metadata[1:60194,11])
metadata[1:60194,11] <- sub("_CDAS", "", metadata[1:60194,11])
metadata[1:60194,11] <- sub("_CD", "", metadata[1:60194,11])
metadata[1:60194,11] <- sub("_AS", "", metadata[1:60194,11])
metadata[1:60194,11] <- sub("AS", "", metadata[1:60194,11])

patient_type <- matrix(0:0, nrow=60194, ncol=1)
metadata <- cbind(metadata,patient_type)

#create a column for data type
metadata[,12] <- metadata[1:60194,2] 

for(row in 1:nrow(metadata)) {
  metadata[row,12] <- sub(paste0(metadata[row,11],"_"), "", metadata[row,12])
}

metadata[1:60194,12] <- sub("_", "", metadata[1:60194,12])
