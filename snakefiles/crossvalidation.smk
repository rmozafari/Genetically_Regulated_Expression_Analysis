# Author: Reza Mozafari
# Date: 2020-09-01
# Email: rezza.mozafari@gmail.com
# Description: Snakemake pipeline for processing genomic data, including 
#              steps for data preparation, variant filtering, 
#              expression analysis, pruning and imputation of gene expression 
#              (GReX) using the modified version of AffiXcan package.


# Define the final targets: all the others depend on them
ALL = ["path/to/output/foo.bed", "path/to/output/bar.png"]

# Run the entire pipeline and copy the files to a mountable directory
rule docker:
    input:
        ALL
    shell:
        "cp -R path/to/dataset/.?* path/to/results || exit 1"  # Generalized paths

################################################
################################################

# Set the config file
configfile: "/path/to/snakemake_config.yml"

# MultiAssayExperiment from tba matrix
rule mae_tba:
    input: "/bioinfo/prj/tba_grex/dataset/1/tba/Cells_-_EBV-transformed_lymphocytes/GEUVADIS_EUR/ir/hocomoco_mono_v10/all/chr{chr}.tba.matrix.gz"
    output: "tba/Cells_EBV-transformed_lymphocytes/GEUVADIS/EUR/chr{chr}.tba.rds"
    shell:
        "mae_constructor {input} {output} || exit 1"  

rule all_mae_tba:
    input: expand("tba/{tissue}/{dataset}/{race}/chr{chr}.tba.rds", tissue="Cells_-_EBV-transformed_lymphocytes", dataset="GEUVADIS_EUR", chr=config["chr"])

# SummarizedExperiment from expression matrix
rule trim_expr_matrix:
    input:
        xpr="/bioinfo/data/geuvadis/E-GEUV-1/analysis_results/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz",
        iid_wanted_sam="selected_samples/all_samples"
    output:
        "expr/{tissue}/{dataset}/{race}/transposed.expr.matrix.gz"
    shell:
        "zcat {input.xpr} | cut -f2,5- | sed 's/Gene_Symbol/IID/' | transpose"
        " | filter_1col 1 <(cat <(echo \"IID\") {input.iid_wanted_sam}) | gzip > {output} || exit 1"  

rule se_expr:
    input:
        "/bioinfo/prj/test_reza/database/MergedVcf/chr{chr}.indexed_GEUVADIS.vcf.gz",
        "expr/{tissue}/{dataset}/{race}/transposed.expr.matrix.gz"
    output:
        "expr/{tissue}/{dataset}/{race}/expr.se.rds"
    shell:
        "se_constructor {input} {output} || exit 1"  

# Expressed genes and regulatory regions associations 
rule regulation_association:
    input:
        xpr_wanted_sam="expr/{tissue}/{dataset}/{race}/transposed.expr.matrix.gz",  # race=ALL
        tba=expand("/bioinfo/prj/tba_grex/dataset/1/tba/Cells_-_EBV-transformed_lymphocytes/GEUVADIS_EUR/ir/hocomoco_mono_v10/all/chr{chr}.tba.matrix.gz", chr=config["chr"])
    output:
        "regions_association/{tissue}/{dataset}/{race}/reg_assoc.gz"
    shell:
        "cat <(echo \"EXPRESSED_REGION\tREGULATORY_REGION\")"
        " <(zcat {input.tba} | cut -f1 | sed '1d' | sed 's/@/\t/g' | bsort"
        " | uniq | filter_1col 1"
        " <(zcat {input.xpr_wanted_sam} | grep IID | cut -f2- | transpose)"
        " | bawk '{{if($2==\"\"){{print $1,$1}}else{{print $0}}}}')"
        " | gzip > {output} || exit 1"  

rule associations_save_rds:
    input:
        "regions_association/{tissue}/{dataset}/{race}/reg_assoc.gz"
    output:
        "regions_association/{tissue}/{dataset}/{race}/reg_assoc.rds"
    shell:
        "R -e 'regionAssoc <- read.table(\"{input}\",header=TRUE);"
        "saveRDS(regionAssoc,file=\"{output}\")' || exit 1"  

##################################################
##### Covariates of the population structure #####

rule list_samples:
    input:
        wanted_sam=lambda wildcards: config["samples"][wildcards.dataset],
        all_sam=lambda wildcards: config["all_samples"][wildcards.dataset]
    output:
        "covariates/{tissue}/{dataset}/{race}/samples"
    shell:
        "filter_1col 1 {input.wanted_sam} < {input.all_sam} > {output} || exit 1"  

rule add_unknown_sex:
    input:
        "covariates/{tissue}/{dataset}/{race}/samples"
    output:
        "covariates/{tissue}/{dataset}/{race}/eigenstrat.ind"
    shell:
        "bawk '{{print $1,\"U\",$2}}' {input} > {output} || exit 1"  

# Selection of individuals and some filters (MAF and LRLD regions)
rule filter_vcf:
    input:
        vcf="/bioinfo/prj/compute_tba/dataset/geuvadis_eur/phased_vcf/chr{chr}.vcf.gz",
        sam="covariates/{tissue}/{dataset}/{race}/samples"
    output:
        "covariates/{tissue}/{dataset}/{race}/chr{chr}.filtered.vcf.gz"
    params:
        MAF_CUTOFF=0.05,
        LRLD=lambda wildcards: config["LRLD_regions"][wildcards.dataset]
    threads: 4  # Specify the number of threads to use
    shell:
        "vcftools --gzvcf {input.vcf} --keep <(cut -f1 {input.sam})"
        " --maf {params.MAF_CUTOFF} --exclude-bed {params.LRLD} --recode"
        " --stdout | gzip > {output}; rm out.log || exit 1"  

# Automatic execution of rule filter_vcf for {chr} and fixed wildcards
rule all_filter_vcf:
    input: expand("covariates/{tissue}/{dataset}/{race}/chr{chr}.filtered.vcf.gz", tissue="Cells_EBV-transformed_lymphocytes", dataset="GEUVADIS", race="EUR", chr=config["chr"])

# Pruning of genetic variants according to the LD pattern
rule ld_pruning:
    input:
        "covariates/{tissue}/{dataset}/{race}/chr{chr}.filtered.vcf.gz"
    output:
        "covariates/{tissue}/{dataset}/{race}/chr{chr}.filtered.prune.in"
    params:
        WINDOW_SIZE=50,
        STEP_SIZE=5,
        RSQUARED_CUTOFF=0.80,
        prefix=lambda wildcards: expand("covariates/{tissue}/{dataset}/{race}/chr{chr}",
            tissue=wildcards.tissue, dataset=wildcards.dataset, race=wildcards.race, chr=wildcards.chr)
    shell:
        "plink --vcf {input} --indep-pairwise {params.WINDOW_SIZE}"
        " {params.STEP_SIZE} {params.RSQUARED_CUTOFF}"
        " --out {params.prefix}.filtered;"
        " rm {params.prefix}.filtered.prune.out {params.prefix}.filtered.nosex"
        " {params.prefix}.filtered.log || exit 1"  

# Filter once again the VCF files, including only the genetic variants
# selected after the LD pruning
rule refilter_vcf:
    input:
        vcf="covariates/{tissue}/{dataset}/{race}/chr{chr}.filtered.vcf.gz",
        snp="covariates/{tissue}/{dataset}/{race}/chr{chr}.filtered.prune.in"
    output:
        "covariates/{tissue}/{dataset}/{race}/chr{chr}.filtered.pruned.vcf.gz"
    shell:
        "vcftools --gzvcf {input.vcf} --snps {input.snp} --recode --stdout"
        " | gzip > {output}; rm out.log || exit 1"  

# The information present in VCF files must be split into multiple input files
# in order to run EIGENSTRAT
# The name of the output files is imposed by fromVcftoEigenstrat
rule execute_fromVcftoEigenstrat:
    input:
        "covariates/{tissue}/{dataset}/{race}/chr{chr}.filtered.pruned.vcf.gz"
    output:
        "covariates/{tissue}/{dataset}/{race}/chr{chr}.filtered.pruned.snps.gz"
    params:
        prefix="covariates/{tissue}/{dataset}/{race}/chr{chr}"
    shell:
        "zcat {input} | fromVcftoEigenstrat -p {params.prefix}.filtered.pruned || exit 1"  

# Automatic execution of rule execute_fromVcftoEigenstrat for {chr} and fixed wildcards
rule all_execute_fromVcftoEigenstrat:
    input: expand("covariates/{tissue}/{dataset}/{race}/chr{chr}.filtered.pruned.snps.gz", tissue="Cells_EBV-transformed_lymphocytes", dataset="GEUVADIS", race="EUR", chr=config["chr"])

rule cat_snps:
    input:
        expand("covariates/Cells_EBV-transformed_lymphocytes/GEUVADIS/EUR/chr{chr}.filtered.pruned.snps.gz", chr=config["chr"])  # should be used the [wildcards.dataset] and so on to expand the wildcards for all chromosomes!!!!
    output:
        "covariates/{tissue}/{dataset}/{race}/all.filtered.pruned.snps.gz"
    shell:
        "zcat {input} | gzip -c > {output} || exit 1" 

rule cat_samples:
    input:
        "covariates/{tissue}/{dataset}/{race}/all.filtered.pruned.snps.gz"
    output:
        "covariates/{tissue}/{dataset}/{race}/all.filtered.pruned.samples.gz"
    params:
        "covariates/{tissue}/{dataset}/{race}/"
    shell:
        # chr*.filtered.pruned.samples.gz are generated by fromVcftoEigenstrat
        "zcat {params}chr*.filtered.pruned.samples.gz | gzip -c > {output} || exit 1"  

# The main output is the one with autovectors/PCs
rule execute_smartpca:
    input:
        pruned_snps="covariates/{tissue}/{dataset}/{race}/all.filtered.pruned.snps.gz",
        pruned_samples="covariates/{tissue}/{dataset}/{race}/all.filtered.pruned.samples.gz",
        eigenstrat="covariates/{tissue}/{dataset}/{race}/eigenstrat.ind"
    output:
        "covariates/{tissue}/{dataset}/{race}/all.filtered.pruned.evec"
    threads: 4  # Specify the number of threads to use
    shell:
        "LD_LIBRARY_PATH=/bioinfo/miniconda3/envs/test_reza/lib; export LD_LIBRARY_PATH;"
        " zcat {input.pruned_snps} > {output}.snp;"
        " zcat {input.pruned_samples} > {output}.geno;"
        " PREFIX=`echo {output} | sed 's/.evec//1`; "
        " smartpca.perl -i {output}.geno -a {output}.snp -b {input.eigenstrat} || exit 1"  

# Refined output; only PC1, PC2, and PC3 are selected
rule cov_refinement:
    input:
        "covariates/{tissue}/{dataset}/{race}/all.filtered.pruned.evec"
    output:
        "covariates/{tissue}/{dataset}/{race}/covariates"   
    shell:
        "cat"
        " <(echo -e 'id\tPC1\tPC2\tPC3')"
        " <(tr -s \" \" < {input} | sed 's/ /\t/g' | cut -f2,3,4,5"
        " | grep -v \"#\") > {output} || exit 1"  

# Final file used by AffiXcan:
rule cov_rds:
    input:
        cov="covariates/{tissue}/{dataset}/{race}/covariates",
        iid="selected_samples/all_samples"
    output:
        "covariates/{tissue}/{dataset}/{race}/{samples}.cov.rds"
    params:
        "covariates/{tissue}/{dataset}/{race}/"
    shell:
        "cat <(grep id < {input.cov} | sed 's/id/IID/1')"
        " <(filter_1col 1 {input.iid} < {input.cov})"
        " > {params}training.cov.tmp; R -e 'trainingCovariates <- read.table("
        " \"{params}training.cov.tmp\",header=TRUE,row.names=1);"
        " saveRDS(trainingCovariates,file=\"{output}\")' || exit 1"  

#################################
### Impute GReX with AffiXcan ###

rule training_the_training_dataset:
    input:
        training_tba=expand("tba/{{tissue}}/{{training_dataset}}/{{race}}/chr{chr}.tba.rds", chr=config["chr"]),
        training_expr="expr/Cells_EBV-transformed_lymphocytes/GEUVADIS/ALL/expr.se.rds",
        region_assoc="regions_association/Cells_EBV-transformed_lymphocytes/GEUVADIS/EUR/reg_assoc.rds",
        training_cov="covariates/Cells_EBV-transformed_lymphocytes/GEUVADIS/EUR/all_samples.cov.rds"
    output:
        "grex/{tissue}/{training_dataset}/{race}/trained_datase_on_all_samples.rds"
    params:
        training_tba_path="tba/{tissue}/{training_dataset}/"
    shell:
        "R -e"
        # load needed data inside R session
        " 'trainingTbaPaths <- list.files(\"{params.training_tba_path}\", full.names=TRUE);"
        " exprMatrix <- readRDS(\"{input.training_expr}\");"
        " regionAssoc <- readRDS(\"{input.region_assoc}\");"
        " cov <- readRDS(\"{input.training_cov}\");"
        # Load AffiXcan functions
        " library(\"AffiXcan\");"
        # Set up BiocParallelParam
        " library(BiocParallel);"
        " BPPARAM <- SnowParam(workers=4);"
        # Execute AffiXcan training phase
        " training <- affiXcanTrain(exprMatrix=exprMatrix,"
            " assay=\"values\",tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc,"
            " cov=cov, varExplained=80, scale=TRUE, BPPARAM=BPPARAM);"
        " saveRDS(training,file=\"{output}\")' || exit 1"  

rule impute_grex:
    # The training dataset used as a temporary solution! 
    input:
        trained_dataset="grex/{tissue}/{training_dataset}/{race}/trained_datase_on_all_samples.rds",
        testing_tba=expand("tba/Cells_EBV-transformed_lymphocytes/GEUVADIS/EUR/chr{chr}.tba.rds", chr=config["chr"])
    output:
        "grex/{tissue}/{training_dataset}/{race}/imputed_GReX.trained_on_all_samples.gz"
    params:
        testing_tba_path="tba/{tissue}/{training_dataset}/",
        grex_prefix="grex/{tissue}/{training_dataset}/{race}/"
    shell:
        "R -e"
        # Load needed data inside R session
        " 'testingTbaPaths <- list.files(\"{params.testing_tba_path}\", full.names=TRUE);"
        " training <- readRDS(\"{input.trained_dataset}\");"
        # Load AffiXcan functions
        " library(\"AffiXcan\");"
        # Set up BiocParallelParam
        " library(BiocParallel);"
        " BPPARAM <- SnowParam(workers=4);"
        # Execute AffiXcan testing phase
        " library(SummarizedExperiment);"
        " exprmatrix <- affiXcanImpute(tbaPaths=testingTbaPaths, affiXcanTraining=training, scale=TRUE, BPPARAM=BPPARAM);"
        # Write the matrix with imputed GReX values
        " GReX <- assays(exprmatrix)$GReX;"
        " t_GReX <- t(GReX);"
        " IID <- c(rownames(t_GReX));"
        " t_GReX_annot <- cbind(IID,t_GReX);"
        " write.table(x=as.data.frame(t_GReX_annot),"
        " file=\"{params.grex_prefix}imputed_GReX.tmp\","
        " quote=FALSE, row.names=FALSE)';"
        # save the imputed GReX in .gz format
        " sed 's/ /\\t/g' < {params.grex_prefix}imputed_GReX.tmp"
        " | gzip > {output}; rm {params.grex_prefix}imputed_GReX.tmp || exit 1" 