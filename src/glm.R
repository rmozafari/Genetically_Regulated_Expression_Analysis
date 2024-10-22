# glm.R
# This script performs logistic regression analysis on ADNI data.

library(dplyr)
library(data.table)
library(qqman)

# Load phenotypic and covariate data
adni.pheno.data <- read.table("/bioinfo/prj/expr_reg_pred/data/adni/adnimerge.tsv", 
                                header = TRUE, sep = "\t")

# Filter for White participants and select relevant columns
adni.white <- adni.pheno.data[adni.pheno.data$PTRACCAT == 'White',]
adni.sub <- adni.white[, c("PTID", "DX.bl", "AGE", "PTGENDER", "PTRACCAT", "PTEDUCAT", "PTMARRY", "APOE4")]
colnames(adni.sub) <- c("IID", "DX", "AGE", "GENDER", "RACE", "EDUCAT", "MARRY", "APOE4")

# Remove duplicate rows, keeping the first occurrence
adni.uniq <- adni.sub %>% distinct(IID, .keep_all = TRUE)

# Recode DX and GENDER
adni.uniq$DX <- as.factor(ifelse(adni.uniq$DX == "CN", "CN", "AD"))
adni.uniq$GENDER <- as.factor(ifelse(adni.uniq$GENDER == "Male", "M", "F"))
adni.uniq$EDUCAT <- as.numeric(adni.uniq$EDUCAT)
adni.uniq$APOE4 <- as.factor(adni.uniq$APOE4)

# Load GReX data and merge with phenotypic data
grex <- read.table("grex/Brain_Cortex/adni/affixcan.imputing_GReX.Brain_Cortex.gz", header = TRUE, sep = "\t")
merge.id <- merge(adni.uniq, grex, by = "IID")
merge.id$RACE <- NULL

# Logistic regression for each gene
glm.output <- lapply(merge.id[, 8:ncol(merge.id)], function(x) {
    summary(glm(merge.id$DX ~ x + merge.id$AGE + merge.id$GENDER + merge.id$EDUCAT + merge.id$MARRY + merge.id$APOE4, family = "binomial"))
})

# Extract coefficients
gene.list <- lapply(glm.output, function(gene) gene$coefficients[rownames(gene$coefficients) == "x"])
gene.df <- do.call(rbind, gene.list)
gene.list.df <- as.data.frame(gene.df)

# Convert row names into a column
setDT(gene.list.df, keep.rownames = TRUE)[]
colnames(gene.list.df) <- c("IID", "Estimate", "Std.Error", "z.value", "p.val")

# Load GTEx annotation and merge with gene list
gtex.annot <- read.table("/bioinfo/data/gtex_analysis_v7/reference/gencode.v19.genes.v7.patched_contigs.gtf", 
                            header = FALSE, sep = "\t")
gtex.chr.tss.ids <- gtex.annot[gtex.annot$V3 == "gene", c("V1", "V4", "V9")]
colnames(gtex.chr.tss.ids) <- c("CHR", "TSS", "IID")

# Clean up gene IDs
gtex.chr.tss.ids$IID <- gsub("\\;.*|gene_id |\\.\\d+$", "", gtex.chr.tss.ids$IID)

# Merge gene data with chromosome info
gene.merge <- merge(gtex.chr.tss.ids, gene.list.df, by = "IID")
gene.merge$CHR <- as.numeric(gene.merge$CHR)  # Ensure CHR is numeric

# Create Manhattan plot
manhattan(gene.merge, chr = "CHR", bp = "TSS", p = "p.val", snp = "IID", col = c("gray10", "skyblue"))