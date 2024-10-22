# Kfold_CrossValidation.R
# This script processes AffiXcan K-fold Cross Validation output to extract gene IDs and R2 values.

library(data.table)
library(reshape2)
library(biomaRt)

# Load AffiXcan 5-fold CV data
affixcan.5k <- readRDS(snakemake@input[[1]])

# Function to process each fold's data
process_fold <- function(data) {
    df <- as.data.table(data, keep.rownames = TRUE)
    melted_df <- melt(t(df)[, 2])  # Transpose and melt
    return(melted_df)
}

# Process all five folds
affix_r2_melts <- lapply(affixcan.5k, process_fold)

# Merge all folds
merged_data <- Reduce(function(x, y) merge(x, y, by = "ID", all.x = TRUE), affix_r2_melts)

# Rename columns
colnames(merged_data) <- c("ID", paste0("r.sq", 1:5))

# Calculate mean R2 for genes with records in all five folds
merged_data$R2 <- rowMeans(merged_data[, 2:6], na.rm = TRUE)

# Filter out genes without R2 values
final_data <- merged_data[!is.na(R2), c("ID", "R2")]

# Write output to file
write.table(final_data, "grex/Cells_EBV_transformed_lymphocytes/gtex/eur/r2.affixcan.5KCV.no.outer.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

###################
### abc_model #####
###################

# Load GTEx bed file for chromosome 1
bed.chr1.gtex <- read.table("/bioinfo/prj/compute_tba/dataset/bed/chr1.ir.hg19.gencode_v19.bed.gz",
                            header = FALSE, sep = "\t")
bed.chr1.gtex$V4 <- gsub("\\.\\d+$", "", bed.chr1.gtex$V4)
colnames(bed.chr1.gtex) <- c("chr", "start", "end", "ID")

# Merge bed file with final data
merged_bed <- merge(bed.chr1.gtex, final_data, by = "ID")

# Select unique genes randomly
uniq_genes <- unique(merged_bed$ID)
gen_150_random <- sample(uniq_genes, 300)

# Convert Ensembl IDs to gene symbols
genes <- data.table(uniq_genes = gen_150_random)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol"),
                values = genes$uniq_genes, mart = mart)

# Remove genes without registered HGNC symbols
G_list <- G_list[G_list$hgnc_symbol != "", ]
colnames(G_list) <- c("ID", "TargetGene")

# Load predictions and filter for specific cell type
allPred_abc_paper <- read.table("abc_model/abc_scores/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz",
                                header = TRUE, sep = "\t")
allPred_lcl <- allPred_abc_paper[allPred_abc_paper$CellType == "GM12878-Roadmap", c(1, 2, 3, 7, 24)]

# Merge predictions with gene list
merge_gtex_paper_genes <- merge(allPred_lcl, G_list, by = "TargetGene")
bed_abc <- merge_gtex_paper_genes[, c(2, 3, 4, 6)]

# Write merged bed file
write.table(bed_abc, "abc_model/compute_tba/gtex/bed/bed.merged.chr1.txt", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

################
## R2 Difference
################

# Merge R2 data for comparison
aff_cv_lcl_gtex_abc_genes_x <- merge(final_data, merge.5fold_final_abc_random_x, by = "ID")
colnames(aff_cv_lcl_gtex_abc_genes_x) <- c("ID", "R2.all", "R2.abc")
aff_cv_lcl_gtex_abc_genes_x$diff_all_abc <- aff_cv_lcl_gtex_abc_genes_x$R2.all - aff_cv_lcl_gtex_abc_genes_x$R2.abc

# Perform t-test
t_test_result <- t.test(aff_cv_lcl_gtex_abc_genes_x$R2.all, aff_cv_lcl_gtex_abc_genes_x$R2.abc, 
                            alternative = "two.sided", paired = TRUE)

# Count rows with positive differences
positive_diff_count <- sum(aff_cv_lcl_gtex_abc_genes_x$diff_all_abc < 0)

# Write results to file
write.table(aff_cv_lcl_gtex_abc_genes_x, "r2.abc.random", sep = "\t", row.names = FALSE, quote = FALSE)