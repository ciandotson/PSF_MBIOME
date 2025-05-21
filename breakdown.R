#### Argument Parsing ####
library(optparse); packageVersion("optparse")
option_list <- list(
  make_option("--raw_soil", type = "character", help = "filepath containing the raw, untrimmed reads for the reads generated from the v4 primers (bulk soil, rhizosphere, and some nodule samples)"),
  make_option("--soil_metadata", type = "character", help = "filepath that contains the Comma Separated Values (csv) file of the soil metadata"),
  make_option("--pheno", type = "character", help = "filepath that contains the Comma Separated Values (csv) file of the phenotype data (biomass and nodule counts)"),
  make_option("--reference", type = "character", help = "filepath that contains the reference database to assign taxonomy to the reads"),
  make_option("--raw_root", type = "character", help = "filepath containing the raw, untrimmed reads for the reads generated from the v5-v7 primers( root endosphere and some nodule samples)"),
  make_option("--root_metadata", type = "character", help = "filepath that contains the Comma Separated Values (csv) file of the root metadata"),
  make_option("--zymo", type = "character", help = "filepath to where the reference fasta files are for the ZYMO mock culture"))

opt <- parse_args(OptionParser(option_list=option_list))

soil.dir <- opt$raw_soil
soil.met <- opt$soil_metadata
nodnbio <- opt$pheno
reference <- opt$reference
root.dir <- opt$raw_root
root.met <- opt$root_metadata
reference <- opt$reference
zymo <- opt$zymo


#### root Primer Removal ####
# Ensure you have the right files #
list.files(root.dir)

# Create a list of the files that are corresponding to the forward and reverse reads #
raw_root.ffp <- sort(list.files(root.dir, pattern = "_R1_001.fastq.gz", full.names = TRUE))
raw_root.rfp <- sort(list.files(root.dir, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Save the names based on the file names #
root.names <- sapply(strsplit(basename(raw_root.ffp), "_"), `[`, 1)
root.names <- as.character(root.names)

# Find all orientations of each primer for primer trimming #
library(Biostrings); packageVersion('Biostrings')
allOrients <- function(primer) {
  # Create all orientations of the input sequence #
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

root.fprimer <- "AACMGGATTAGATACCCKG"
root.rprimer <- "ACGTCATCCCCACCTTCC"

root.fori <- allOrients(root.fprimer)
root.rori <- allOrients(root.fprimer)

# Make filepaths for pretrimmed fastqs #
system('mkdir ./reads/pretrim/root_pretrim')
pre_root.ffp <- file.path('./reads/pretrim/root_pretrim', paste0(root.names, '_pretrim_R1.fastq.gz'))
pre_root.rfp <- file.path('./reads/pretrim/root_pretrim', paste0(root.names, '_pretrim_R2.fastq.gz'))

# Filter reads less than 75 bp and save the filtered fastqs to the pretrim filepaths #
library(dada2); packageVersion('dada2')
root_prefilt.track <- filterAndTrim(raw_root.ffp, pre_root.ffp, raw_root.rfp, pre_root.rfp, minLen = 75,
                                    compress = TRUE, multithread = TRUE)

# Take reverse complements of the forward and reverse primers #
root_fprim.rc <- dada2::rc(root.fprimer)
root_rprim.rc <- dada2::rc(root.fprimer)

# Create filepaths for the forward and reverse fastqs with primers trimmed #
system('mkdir ./reads/ptrimmed/root_ptrimmed')
pt_root.ffp <- file.path("./reads/ptrimmed/root_ptrimmed", paste0(root.names, "_ptrimmed_R1.fastq.gz"))
pt_root.rfp <- file.path("./reads/ptrimmed/root_ptrimmed", paste0(root.names, "_ptrimmed_R2.fastq.gz"))

# Perform primmer trimming using cutadapt #
for(i in seq_along(pre_root.ffp)){
  system(paste0('cutadapt -g ', root.fprimer, ' -a ', root_rprim.rc, ' -G ', root.rprimer, ' -A ', root_fprim.rc, ' --pair-filter=any -m 50:50 -o ', pt_root.ffp[i], ' -p ', pt_root.rfp[i], ' ', pre_root.ffp[i], ' ', pre_root.rfp[i]))
}

# Perform additional filtering prior to dada2 algorithm #
system('mkdir ./reads/postfilt/root_postfilt')
post_root.ffp <- file.path('./reads/postfilt/root_postfilt', paste0(root.names, '_filt_R1.fastq.gz'))
post_root.rfp <- file.path('./reads/postfilt/root_postfilt', paste0(root.names, '_filt_R2.fastq.gz'))

root_postfilt.track <- filterAndTrim(pt_root.ffp, post_root.ffp, pt_root.rfp, post_root.rfp, truncLen=c(230,200),
                                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                                     compress=TRUE, multithread=TRUE)

save.image("./root.RData")
#### dada2 Implementation for the root Samples ####
# Learn the error rates that are specific to your data #
root_for.er <- learnErrors(post_root.ffp, multithread=TRUE, verbose = TRUE)
root_rev.er <- learnErrors(post_root.rfp, multithread=TRUE, verbose = TRUE)

# Plot the learned errors #
plotErrors(root_for.er, nominalQ=TRUE)
plotErrors(root_rev.er, nominalQ=TRUE)

# Dereplicate the reads #
root.fderep <- derepFastq(post_root.ffp, verbose=TRUE)
root.rderep <- derepFastq(post_root.rfp, verbose=TRUE)

# Construct the dada-class object #
root.fdada <- dada(root.fderep, err=root_for.er, multithread=TRUE, verbose = TRUE)
root.rdada <- dada(root.rderep, err=root_rev.er, multithread=TRUE, verbose = TRUE)

save.image("./root.RData")
# Merge the denoised forward and reversed reads #
root.remerged <- mergePairs(root.fdada, post_root.ffp, root.rdada, post_root.rfp, verbose=TRUE)

# Construct Sequence (ASV) Table #
root.st <- makeSequenceTable(root.remerged)
dim(root.st)
table(nchar(getSequences(root.st)))

# Remove chimeras #
root_nochim.st <- removeBimeraDenovo(root.st, method="consensus", multithread=TRUE, verbose=TRUE)
dim(root_nochim.st)

# Determine the ratio of non-chimeras to all reads #
sum(root_nochim.st)/sum(root.st)
root_nochim.st <- t(root_nochim.st)
colnames(root_nochim.st) <- root.names

save.image("./root.RData")