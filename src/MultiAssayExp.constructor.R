#!/usr/bin/env Rscript

#  The script reads a table with the binding affinity data, splits it by the REGION column, and creates a MultiAssayExperiment object.
#  The script is called by the following bash script:
#  mae_constructor.smk

library("MultiAssayExperiment")
library("plyr")
args <- commandArgs(TRUE)
tba <- read.table(file=args[1],header=TRUE)
tba_by_region <- dlply(.data=tba, .(REGION), .fun=function(data){
	data$REGION<-NULL;rownames(data)<-data$IID;data$IID<-NULL;return(data)
	}, 
	.parallel=FALSE)
totalBindingAffinity <- MultiAssayExperiment(experiments=tba_by_region)
saveRDS(totalBindingAffinity,file=args[2])


