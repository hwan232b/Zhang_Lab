# Hannah HIV code
# last updated 6/15/22

library(Matrix)
library(readr)

C1matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775588_C1matrix.mtx.gz")
C1genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775588_C1genes.tsv.gz", col_names = FALSE)
C1cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775588_C1barcodes.tsv.gz", col_names = FALSE)$X1
C1gene_ids <- C1genes$X1
rownames(C1matrix) <- C1gene_ids
colnames(C1matrix) <- C1cell_ids


Q1matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775591_Q1matrix.mtx.gz")
Q1genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775589_Q2genes.tsv.gz", col_names = FALSE)
Q1cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775591_Q1barcodes.tsv.gz", col_names = FALSE)$X1
Q1gene_ids <- Q1genes$X1
rownames(Q1matrix) <- Q1gene_ids
colnames(Q1matrix) <- Q1cell_ids


Q2matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775589_Q2matrix.mtx.gz")
Q2genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775589_Q2genes.tsv.gz", col_names = FALSE)
Q2gene_ids <- Q2genes$X1
Q2cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775589_Q2barcodes.tsv.gz", col_names = FALSE)$X1
rownames(Q2matrix) <- Q2gene_ids
colnames(Q2matrix) <- Q2cell_ids


Q3matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775590_Q3matrix.mtx.gz")
Q3genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775590_Q3genes.tsv.gz", col_names = FALSE)
Q3gene_ids <- Q3genes$X1
Q3cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775590_Q3barcodes.tsv.gz", col_names = FALSE)$X1
rownames(Q3matrix) <- Q3gene_ids
colnames(Q3matrix) <- Q3cell_ids


Q4matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775592_Q4matrix.mtx.gz")
Q4genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775592_Q4genes.tsv.gz", col_names = FALSE)
Q4gene_ids <- Q4genes$X1
Q4cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775592_Q4barcodes.tsv.gz", col_names = FALSE)$X1
rownames(Q4matrix) <- Q4gene_ids
colnames(Q4matrix) <- Q4cell_ids


Q5matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775593_Q5matrix.mtx.gz")
Q5genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775593_Q5genes.tsv.gz", col_names = FALSE)
Q5gene_ids <- Q5genes$X1
Q5cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775593_Q5barcodes.tsv.gz", col_names = FALSE)$X1
rownames(Q5matrix) <- Q5gene_ids
colnames(Q5matrix) <- Q5cell_ids

Q7matrix <- readMM("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775594_Q7matrix.mtx.gz")
Q7genes <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775594_Q7genes.tsv.gz", col_names = FALSE)
Q7gene_ids <- Q7genes$X1
Q7cell_ids <- read_tsv("/Users/hannahwang/github/Projects/Zhang_lab/GSE157829_RAW/GSM4775594_Q7barcodes.tsv.gz", col_names = FALSE)$X1
rownames(Q7matrix) <- Q7gene_ids
colnames(Q7matrix) <- Q7cell_ids

rownames(Q7matrix) <- sub("^hg19_", "", rownames(Q7matrix))

colnames(C1matrix) <- paste("C1_", colnames(C1matrix))
colnames(Q1matrix) <- paste("Q1_", colnames(Q1matrix))
colnames(Q2matrix) <- paste("Q2_", colnames(Q2matrix))
colnames(Q3matrix) <- paste("Q3_", colnames(Q3matrix))
colnames(Q4matrix) <- paste("Q4_", colnames(Q4matrix))
colnames(Q5matrix) <- paste("Q5_", colnames(Q5matrix))
colnames(Q7matrix) <- paste("Q7_", colnames(Q7matrix))


rownames1 <- (intersect(rownames(Q7matrix), rownames(Q5matrix)))

Q7matrix_part <- Q7matrix[rownames(Q7matrix) %in% rownames1, ]

length(intersect(rownames(Q7matrix), rownames(Q5matrix)))
length(intersect(rownames(Q7matrix_part), rownames(Q5matrix)))


HIV_combined <- cbind(C1matrix, Q1matrix, Q2matrix, Q3matrix, Q4matrix, Q5matrix, Q7matrix_part)

#max value
which(HIV_combined==max(HIV_combined), arr.ind=TRUE)

# meta: cell, sample, disease (HIV or control)

HIV_meta <- matrix(nrow = 35750, ncol = 2)
rownames(HIV_meta) <- colnames(HIV_combined)

colnames(HIV_meta) <- c("Sample", "Disease")

#rownames(HIV_meta) <- sub(".*_ ", "", rownames(HIV_meta))
#rownames(HIV_meta) <- sub("*..", "", rownames(HIV_meta))

HIV_meta[,1] <- sub("_.*", "", colnames(HIV_combined))

for (row in 1:nrow(HIV_meta)) {
if("C1" %in% HIV_meta[row,1]){
  HIV_meta[row,2] <- ("Control")
} else {HIV_meta[row,2] <- ("HIV")
}
}

saveRDS(HIV_meta, "HIV_ShaoboWang.rds")
saveRDS(HIV_combined, "HIV_exp_ShaoboWang.rds")