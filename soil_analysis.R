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
if(!require(Maaslin2)) {BiocManager::install("Maaslin2")}
if(!require(microbiomeutilities)) {remotes::install_github("microsud/microbiomeutilities")}
if(!require(msa)) {BiocManager::install("msa")}
if(!require(phangorn)) {BiocManager::install("phangorn")}
if(!require(Biostrings)) {BiocManager::install("Biostrings")}
if(!require(ShortRead)) {BiocManager::install("ShortRead")}


# Rscript ~/PSF_MBIOME/soil_analysis.R --raw_soil ~/test/reads/raw/soil_reads --soil_metadata ~/test/metadata/soil_metadata.csv --pheno ~/test/nodnbio.csv --reference ~/test/reference/rdp_19_toGenus_trainset.fa.gz | cat > PSF_log.txt #

#### Argument Parsing ####
library(optparse); packageVersion("optparse")
option_list <- list(
  make_option("--raw_soil", type = "character", help = "filepath containing the raw, untrimmed reads for the reads generated from the v4 primers (bulk soil, rhizosphere, and some nodule samples)"),
  make_option("--soil_metadata", type = "character", help = "filepath that contains the Comma Separated Values (csv) file of the soil metadata"),
  make_option("--pheno", type = "character", help = "filepath that contains the Comma Separated Values (csv) file of the phenotype data (biomass and nodule counts)"),
  make_option("--reference", type = "character", help = "filepath that contains the reference database to assign taxonomy to the reads"))

opt <- parse_args(OptionParser(option_list=option_list))
soil.dir <- opt$raw_soil
soil.met <- opt$soil_metadata
nodnbio <- opt$pheno
reference <- opt$reference

load('~/test.RData')
#### dada2 Implementation for the Soil Samples ####
# Learn the error rates that are specific to your data #
soil_for.er <- learnErrors(post_soil.ffp, multithread=FALSE, verbose = TRUE)
soil_rev.er <- learnErrors(post_soil.rfp, multithread=FALSE, verbose = TRUE)

# Plot the learned errors #
plotErrors(soil_for.er, nominalQ=TRUE)
plotErrors(soil_rev.er, nominalQ=TRUE)

# Dereplicate the reads #
soil.fderep <- derepFastq(post_soil.ffp, verbose=TRUE)
soil.rderep <- derepFastq(post_soil.rfp, verbose=TRUE)

# Construct the dada-class object #
soil.fdada <- dada(soil.fderep, err=soil_for.er, multithread=FALSE, verbose = TRUE)
soil.rdada <- dada(soil.rderep, err=soil_rev.er, multithread=FALSE, verbose = TRUE)

# Merge the denoised forward and reversed reads #
soil.remerged <- mergePairs(soil.fdada, post_soil.ffp, soil.rdada, post_soil.rfp, verbose=TRUE)

# Construct Sequence (ASV) Table #
soil.st <- makeSequenceTable(soil.remerged)
dim(soil.st)
table(nchar(getSequences(soil.st)))

# Remove chimeras #
soil_nochim.st <- removeBimeraDenovo(soil.st, method="consensus", multithread=FALSE, verbose=TRUE)
dim(soil_nochim.st)

# Determine the ratio of non-chimeras to all reads #
sum(soil_nochim.st)/sum(soil.st)
soil_nochim.st <- t(soil_nochim.st)

# track reads through the pipeline #
getN <- function(x) sum(getUniques(x))
soil_final.track <- cbind(soil_prefilt.track[,1], soil_prefilt.track[,2], soil_postfilt.track[,2], sapply(soil.fdada, getN), sapply(soil.rdada, getN), sapply(soil.remerged, getN), rowSums(soil_nochim.st))
colnames(soil_final.track) <- c("pre-cutadapt", "post-cutadapt", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(soil_final.track) <- soil.names
soil_final.track <- as.data.frame(soil_final.track)

# Assign Taxonomy #
soil_rdp.taxa <- assignTaxonomy(rownames(soil_nochim.st), refFasta = reference, verbose = TRUE)
soil_rdp.taxa <- as.matrix(soil_rdp.taxa)

# Load the metadata #
soil_raw.met <- read.csv2(soil_metadata, sep = ',')

save.image("~/test.RData")