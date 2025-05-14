#### Package Installation (if neccessary)####
if(!require(optparse)) { install.packages("optparse") }
if(!require(ggplot2)) { install.packages('ggplot2')}
if(!require(dplyr)) { install.packages('dplyr')}
if(!require(vegan)) { install.packages('vegan')}
if(!require(cgwtools)) { install.packages('cgwtools')}
if(!require(patchwork)) { install.packages('patchwork')}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require(dada2)) {BiocManager::install("dada2", version = "3.20")}
if(!require(phyloseq)) {BiocManager::install("phyloseq", version = "3.20")}
if(!require(microbiome)) {BiocManager::install("microbiome")}
if(!require(DECIPHER)) {BiocManager::install("DECIPHER")}
if(!require(Maaslin2)) {BiocManager::install("Maaslin2")}
if(!require(microbiomeutilities)) {remotes::install_github("microsud/microbiomeutilities")}
if(!require(msa)) {BiocManager::install("msa")}
if(!require(phangorn)) {BiocManager::install("phangorn")}

option_list("--raw_reads", type = "character", help = "filepath containing the raw, untrimmed reads",
            "--metadata", type = "character", help = "filepath that contains the Comma Separated Values (csv) file of the metadata",
            "--pheno", type = "character", help = "filepath that contains the Comma Separated Values (csv) file of the phenotype data (biomass and nodule counts)",
            "--cutadapt", type = "character", help = "filepath to the usable version of cutadapt")


