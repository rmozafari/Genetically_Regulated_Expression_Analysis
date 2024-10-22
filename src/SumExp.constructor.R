#!/bioinfo/miniconda3/envs/tba_grex/bin/Rscript

# The script reads a table with the expression data, transposes it, and creates a SummarizedExperiment object.
# The script is called by the following bash script:
# se_constructor.smk    

library("SummarizedExperiment")
args <- commandArgs(TRUE)
expr <- read.table(args[1],header=TRUE,row.names=1)
t.expr <- t(as.data.frame(expr))
exprMatrix <- SummarizedExperiment(assays=list(values=t.expr))
saveRDS(exprMatrix,file=args[2])