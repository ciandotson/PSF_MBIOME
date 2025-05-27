# Rscript ~/PSF_MBIOME/analysis.R --raw_soil ~/test_PSF/reads/soil_reads --raw_root ~/test_PSF/reads/endo_reads --soil_metadata ~/PSF_MBIOME/metadata/soil_metadata.csv --pheno ~/PSF_MBIOME/nodnbio.csv --reference ~/PSF_MBIOME/reference/rdp_19_toGenus_trainset.fa.gz | cat > PSF_log.txt #

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
zymo <- opt$zymo

#### Nodule Count and Biomass Data Visualization ####
# Read in the phenotypic data and clean it for analysis #
nodnbio.data <- read.csv2(nodnbio, sep = ',')
colnames(nodnbio.data) <- c('Plant_Sample', 'Soil_Treatment', 'Nodule_Count', 'Aboveground_Biomass')
nodnbio.data$Aboveground_Biomass <- as.numeric(nodnbio.data$Aboveground_Biomass)
for(i in 1:nrow(nodnbio.data)){
  nodnbio.data$Grouped[i] <- paste0(nodnbio.data$Plant_Sample[i], '; ', nodnbio.data$Soil_Treatment[i])
}

# Group all observations by Plant Species and Soil Treatment and find the group mean and standard error # 
library(dplyr); packageVersion('dplyr')
nodnbio.mnsd <- nodnbio.data %>%
  group_by(Plant_Sample, Soil_Treatment) %>%
  summarize(
    nod.mean = mean(Nodule_Count),
    bio.mean = mean(Aboveground_Biomass),
    nod.sd = sd(Nodule_Count),
    bio.sd = sd(Aboveground_Biomass),
    .groups = "drop" # Prevent grouping in the result
  )

# Clean up the names by removing dashes #
nodnbio.data$Soil_Treatment <- gsub('Non-PSF Soil', 'Non PSF Soil', nodnbio.data$Soil_Treatment)

# subset each set of values by plant species #
fb_nod.data <- filter(nodnbio.data, Plant_Sample == 'S. helvola')
cc_nod.data <- filter(nodnbio.data, Plant_Sample == 'C. fasciculata')
ds_nod.data <- filter(nodnbio.data, Plant_Sample == 'D. illinoense')
hp_nod.data <- filter(nodnbio.data, Plant_Sample == 'A. braceteata')
cl_nod.data <- filter(nodnbio.data, Plant_Sample == 'T. repens')
md_nod.data <- filter(nodnbio.data, Plant_Sample == 'M. truncatula')

# Perform ANOVA on nodule counts for each plant group #
fb_nod.aov <- aov(Nodule_Count~Soil_Treatment, fb_nod.data)
cc_nod.aov <- aov(Nodule_Count~Soil_Treatment, cc_nod.data)
ds_nod.aov <- aov(Nodule_Count~Soil_Treatment, ds_nod.data)
hp_nod.aov <- aov(Nodule_Count~Soil_Treatment, hp_nod.data)
cl_nod.aov <- aov(Nodule_Count~Soil_Treatment, cl_nod.data)
md_nod.aov <- aov(Nodule_Count~Soil_Treatment, md_nod.data)

# Look at summaries of each ANOVA of noduel values #
summary(fb_nod.aov)
summary(cc_nod.aov)
summary(ds_nod.aov)
summary(hp_nod.aov)
summary(cl_nod.aov)
summary(md_nod.aov)

# Perform Tukey Honest Significant Difference for each test #
fb_nod.hsd <- TukeyHSD(fb_nod.aov)
cc_nod.hsd <- TukeyHSD(cc_nod.aov)
ds_nod.hsd <- TukeyHSD(ds_nod.aov)
hp_nod.hsd <- TukeyHSD(hp_nod.aov)
cl_nod.hsd <- TukeyHSD(cl_nod.aov)
md_nod.hsd <- TukeyHSD(md_nod.aov)

# Assign letters based on the Tukey HSD results #
library(multcompView); packageVersion('multcompView')
fb_nod.let <- multcompLetters4(fb_nod.aov, fb_nod.hsd)
cc_nod.let <- multcompLetters4(cc_nod.aov, cc_nod.hsd)
ds_nod.let <- multcompLetters4(ds_nod.aov, ds_nod.hsd)
hp_nod.let <- multcompLetters4(hp_nod.aov, hp_nod.hsd)
cl_nod.let <- multcompLetters4(cl_nod.aov, cl_nod.hsd)
md_nod.let <- multcompLetters4(md_nod.aov, md_nod.hsd)

# Save a master list of the letters assigned to each comparison #
nod.let <- c(hp_nod.let$Soil_Treatment$Letters,
             cc_nod.let$Soil_Treatment$Letters,
             ds_nod.let$Soil_Treatment$Letters,
             md_nod.let$Soil_Treatment$Letters,
             fb_nod.let$Soil_Treatment$Letters,
             cl_nod.let$Soil_Treatment$Letters)

# Create a data.frame that has both mean/ses and nodule comparisons #
nod.data <- cbind(nodnbio.mnsd, nod.let)

# Clean up the letters to match their tests and clean for plotting #
nod.data$nod.let[8] <- 'a'
nod.data$nod.let[14] <- 'a'
nod.data$nod.let[16] <- 'ab'
nod.data$nod.let[17] <- 'a'
nod.data$nod.let[3] <- 'b'
nod.data[1,3:6] <- 0
nod.data$Group <- factor(nod.data$Plant_Sample, levels = c('T. repens', 'M. truncatula', 'S. helvola', 'C. fasciculata', 'D. illinoense', 'A. braceteata'))

# Plot the results #
library(ggplot2); packageVersion('ggplot2')
library(ggprism); packageVersion('ggprism')
nod.plot <- ggplot(nod.data, aes(x = Group, y = nod.mean, fill = Soil_Treatment, color = Soil_Treatment)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = nod.mean, ymax = nod.mean + nod.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = nod.let, y = nod.mean + nod.sd + 1), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 10) +
  ylab('Nodule Count') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,65),expand = expansion(mult = c(0, 0.05)), sec.axis = dup_axis()) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24, family = 'Liberation Sans', face = 'bold'),
        axis.text.x = element_text(size = 22, face = c('bold.italic')),
        axis.title.y = element_text(size = 24),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm')) +
  labs(tag = "A.")

# Do the same for biomass as seen in lines 74-159 #
fb_bio.data <- filter(nodnbio.data, Plant_Sample == 'S. helvola')
cc_bio.data <- filter(nodnbio.data, Plant_Sample == 'C. fasciculata')
ds_bio.data <- filter(nodnbio.data, Plant_Sample == 'D. illinoense')
hp_bio.data <- filter(nodnbio.data, Plant_Sample == 'A. braceteata')
cl_bio.data <- filter(nodnbio.data, Plant_Sample == 'T. repens')
md_bio.data <- filter(nodnbio.data, Plant_Sample == 'M. truncatula')

fb_bio.aov <- aov(Nodule_Count~Soil_Treatment, fb_bio.data)
cc_bio.aov <- aov(Nodule_Count~Soil_Treatment, cc_bio.data)
ds_bio.aov <- aov(Nodule_Count~Soil_Treatment, ds_bio.data)
hp_bio.aov <- aov(Nodule_Count~Soil_Treatment, hp_bio.data)
cl_bio.aov <- aov(Nodule_Count~Soil_Treatment, cl_bio.data)
md_bio.aov <- aov(Nodule_Count~Soil_Treatment, md_bio.data)

summary(fb_bio.aov)
summary(cc_bio.aov)
summary(ds_bio.aov)
summary(hp_bio.aov)
summary(cl_bio.aov)
summary(md_bio.aov)

fb_bio.hsd <- TukeyHSD(fb_bio.aov)
cc_bio.hsd <- TukeyHSD(cc_bio.aov)
ds_bio.hsd <- TukeyHSD(ds_bio.aov)
hp_bio.hsd <- TukeyHSD(hp_bio.aov)
cl_bio.hsd <- TukeyHSD(cl_bio.aov)
md_bio.hsd <- TukeyHSD(md_bio.aov)

fb_bio.let <- multcompLetters4(fb_bio.aov, fb_bio.hsd)
cc_bio.let <- multcompLetters4(cc_bio.aov, cc_bio.hsd)
ds_bio.let <- multcompLetters4(ds_bio.aov, ds_bio.hsd)
hp_bio.let <- multcompLetters4(hp_bio.aov, hp_bio.hsd)
cl_bio.let <- multcompLetters4(cl_bio.aov, cl_bio.hsd)
md_bio.let <- multcompLetters4(md_bio.aov, md_bio.hsd)

bio.let <- c(hp_bio.let$Soil_Treatment$Letters,
             cc_bio.let$Soil_Treatment$Letters,
             ds_bio.let$Soil_Treatment$Letters,
             md_bio.let$Soil_Treatment$Letters,
             fb_bio.let$Soil_Treatment$Letters,
             cl_bio.let$Soil_Treatment$Letters)

bio.data <- cbind(nodnbio.mnsd, bio.let)

bio.data$bio.let[8] <- 'a'
bio.data$bio.let[14] <- 'a'
bio.data$bio.let[16] <- 'ab'
bio.data$bio.let[17] <- 'a'
bio.data$bio.let[3] <- 'b'

bio.data[1,3:6] <- 0

bio.data$Group <- factor(bio.data$Plant_Sample, levels = c('T. repens', 'M. truncatula', 'S. helvola', 'C. fasciculata', 'D. illinoense', 'A. braceteata'))
bio.plot <- ggplot(bio.data, aes(x = Group, y = bio.mean, fill = Soil_Treatment, color = Soil_Treatment)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = bio.mean, ymax = bio.mean + bio.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = bio.let, y = bio.mean + bio.sd + 0.01), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 10) +
  ylab('Aboveground Biomass (grams)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,1.6),expand = expansion(mult = c(0, 0.05)), sec.axis = dup_axis()) +
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 24, family = 'Liberation Sans', face = 'bold'),
        axis.text.x = element_text(size = 22, face = c('bold.italic')),
        axis.title.y.left = element_text(size = 20, vjust = 0.5),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm')) +
  labs(tag = "B.")

# Join the nodule count and biomass plots into one plot #
library(patchwork); packageVersion('patchwork')
nodnbio.plot <- (nod.plot) /
  (bio.plot) +
  plot_layout(guides = 'keep') &
  theme(plot.tag = element_text(size = 20))

#### FastQC on the reads ####
# Call fastqc to do quality control from the command line #
system('mkdir QC')
system('mkdir ./QC/raw_qc')
system('mkdir ./QC/raw_qc/raw_soil_qc')
system(paste0("fastqc --noextract ",soil.dir, "*fastq.gz -o ./QC/raw_qc/raw_soil_qc"))
system('mkdir ./QC/raw_qc/raw_root_qc')
system(paste0("fastqc --noextract ",root.dir, "*fastq.gz -o ./QC/raw_qc/raw_root_qc"))
system("rm ./QC/raw_qc/raw_soil_qc/*.zip")
system("rm ./QC/raw_qc/raw_root_qc/*.zip")

#### Soil Primer Removal ####
# Ensure you have the right files #
list.files(soil.dir)

# Create a list of the files that are corresponding to the forward and reverse reads #
raw_soil.ffp <- sort(list.files(soil.dir, pattern = "_R1_001.fastq.gz", full.names = TRUE))
raw_soil.rfp <- sort(list.files(soil.dir, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Save the names based on the file names #
soil.names <- sapply(strsplit(basename(raw_soil.ffp), "_"), `[`, 1)
soil.names <- as.character(soil.names)

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

soil.fprimer <- "GTGCCAGCMGCCGCGGTAA"
soil.rprimer <- "GGACTACHVGGGTWTCTAAT"

soil.fori <- allOrients(soil.fprimer)
soil.rori <- allOrients(soil.rprimer)

# Make filepaths for pretrimmed fastqs #
system('mkdir ./reads/pretrim')
system('mkdir ./reads/pretrim/soil_pretrim')
pre_soil.ffp <- file.path('./reads/pretrim/soil_pretrim', paste0(soil.names, '_pretrim_R1.fastq.gz'))
pre_soil.rfp <- file.path('./reads/pretrim/soil_pretrim', paste0(soil.names, '_pretrim_R2.fastq.gz'))

# Filter reads less than 75 bp and save the filtered fastqs to the pretrim filepaths #
library(dada2); packageVersion('dada2')
soil_prefilt.track <- filterAndTrim(raw_soil.ffp, pre_soil.ffp, raw_soil.rfp, pre_soil.rfp, minLen = 75,
                                    compress = TRUE, multithread = TRUE)

# Take reverse complements of the forward and reverse primers #
soil_fprim.rc <- dada2::rc(soil.fprimer)
soil_rprim.rc <- dada2::rc(soil.fprimer)

# Create filepaths for the forward and reverse fastqs with primers trimmed #
system('mkdir ./reads/ptrimmed')
system('mkdir ./reads/ptrimmed/soil_ptrimmed')
pt_soil.ffp <- file.path("./reads/ptrimmed/soil_ptrimmed", paste0(soil.names, "_ptrimmed_R1.fastq.gz"))
pt_soil.rfp <- file.path("./reads/ptrimmed/soil_ptrimmed", paste0(soil.names, "_ptrimmed_R2.fastq.gz"))

# Perform primmer trimming using cutadapt #
for(i in seq_along(pre_soil.ffp)){
  system(paste0('cutadapt -g ', soil.fprimer, ' -a ', soil_rprim.rc, ' -G ', soil.rprimer, ' -A ', soil_fprim.rc, ' --pair-filter=any -m 50:50 -o ', pt_soil.ffp[i], ' -p ', pt_soil.rfp[i], ' ', pre_soil.ffp[i], ' ', pre_soil.rfp[i]))
}

# Perform additional filtering prior to dada2 algorithm #
system('mkdir ./reads/postfilt')
system('mkdir ./reads/postfilt/soil_postfilt')
post_soil.ffp <- file.path('./reads/postfilt/soil_postfilt', paste0(soil.names, '_filt_R1.fastq.gz'))
post_soil.rfp <- file.path('./reads/postfilt/soil_postfilt', paste0(soil.names, '_filt_R2.fastq.gz'))

soil_postfilt.track <- filterAndTrim(pt_soil.ffp, post_soil.ffp, pt_soil.rfp, post_soil.rfp, truncLen=c(230,200),
                                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                                     compress=TRUE, multithread=TRUE)

save.image("./test.RData")
#### dada2 Implementation for the Soil Samples ####
# Learn the error rates that are specific to your data #
soil_for.er <- learnErrors(post_soil.ffp, multithread=TRUE, verbose = TRUE)
soil_rev.er <- learnErrors(post_soil.rfp, multithread=TRUE, verbose = TRUE)

# Plot the learned errors #
plotErrors(soil_for.er, nominalQ=TRUE)
plotErrors(soil_rev.er, nominalQ=TRUE)

# Dereplicate the reads #
soil.fderep <- derepFastq(post_soil.ffp, verbose=TRUE)
soil.rderep <- derepFastq(post_soil.rfp, verbose=TRUE)

# Construct the dada-class object #
soil.fdada <- dada(soil.fderep, err=soil_for.er, multithread=TRUE, verbose = TRUE)
soil.rdada <- dada(soil.rderep, err=soil_rev.er, multithread=TRUE, verbose = TRUE)

save.image("./test.RData")
# Merge the denoised forward and reversed reads #
soil.remerged <- mergePairs(soil.fdada, post_soil.ffp, soil.rdada, post_soil.rfp, verbose=TRUE)

# Construct Sequence (ASV) Table #
soil.st <- makeSequenceTable(soil.remerged)
dim(soil.st)
table(nchar(getSequences(soil.st)))

# Remove chimeras #
soil_nochim.st <- removeBimeraDenovo(soil.st, method="consensus", multithread=TRUE, verbose=TRUE)
dim(soil_nochim.st)

# Determine the ratio of non-chimeras to all reads #
sum(soil_nochim.st)/sum(soil.st)
soil_nochim.st <- t(soil_nochim.st)
colnames(soil_nochim.st) <- soil.names

# track reads through the pipeline #
getN <- function(x) sum(getUniques(x))
soil_final.track <- cbind(soil_prefilt.track[,1], soil_prefilt.track[,2], soil_postfilt.track[,2], sapply(soil.fdada, getN), sapply(soil.rdada, getN), sapply(soil.remerged, getN), colSums(soil_nochim.st))
colnames(soil_final.track) <- c("pre-cutadapt", "post-cutadapt", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(soil_final.track) <- soil.names
soil_final.track <- as.data.frame(soil_final.track)

# Assign Taxonomy #
soil_rdp.taxa <- assignTaxonomy(rownames(soil_nochim.st), refFasta = reference, multithread = TRUE, verbose = TRUE)
soil_rdp.taxa <- as.matrix(soil_rdp.taxa)

# Load the metadata #
soil_raw.met <- read.csv2(soil.met, sep = ',')
rownames(soil_raw.met) <- soil_raw.met$Sample
soil_raw.met <- soil_raw.met[,c('Sample', 'Plant', 'Soil_Treatment', 'Compartment')]

# Check to see if ASVs were denoise properly #
unqs.mock <- soil_nochim.st[,"ZymoMockDNA"]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
mock.ref <- getSequences(zymo)
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

#### Phyloseq Object Construction and Filtering for Soils ####
library(phyloseq)
raw_soil.ps <- phyloseq(otu_table(soil_nochim.st, taxa_are_rows = TRUE),
                        sample_data(soil_raw.met),
                        tax_table(soil_rdp.taxa))

raw_soil.dna <- Biostrings::DNAStringSet(taxa_names(raw_soil.ps))
names(raw_soil.dna) <- taxa_names(raw_soil.ps)
raw_soil.ps <- merge_phyloseq(raw_soil.ps, raw_soil.dna)

raw_soil.ps <- subset_taxa(raw_soil.ps, Phylum != "Plantae")
raw_soil.ps <- subset_taxa(raw_soil.ps, Phylum != "Cyanobacteriota")

raw_soil.ps <- subset_samples(raw_soil.ps, Plant != 'N/A')
raw_soil.ps <- subset_samples(raw_soil.ps, Compartment != 'Leaf')
raw_soil.ps <- subset_samples(raw_soil.ps, Plant != 'Lupine')

raw_soil.ps <- subset_taxa(raw_soil.ps, taxa_sums(raw_soil.ps) > 1000)
colnames(tax_table(raw_soil.ps)) <- c('rdp_Kingdom', 'rdp_Phylum', 'rdp_Class', 'rdp_Order', 'rdp_Family', 'rdp_Genus')

taxa_names(raw_soil.ps) <- paste0("ASV", seq(ntaxa(raw_soil.ps)))

decompose_ps <- function(ps, label){
  # function that decomposes a phyloseq object into separate data.frame and refseq objects (does not include tree) #
  tax.tab <- as.data.frame(tax_table(ps))
  otu.tab <- as.data.frame(otu_table(ps))
  met.tab <- as(sample_data(ps), 'data.frame')
  dna.tab <- refseq(ps)
  fra.tab <- cbind(tax.tab, otu.tab)
    decomposed = list(
    tax = tax.tab,
    otu = otu.tab,
    met = met.tab,
    dna = dna.tab,
    fra = fra.tab
    )
  assign(label, decomposed, envir = .GlobalEnv)
  invisible(decomposed)
}

decompose_ps(raw_soil.ps, 'raw_soil')
#### Cross-Validation of Soil Reads Using BLAST ####
library(rBLAST)

# Create local blast database from the 16S rRNA database using rBLAST #
blast.tar <- blast_db_get("16S_ribosomal_RNA.tar.gz", baseURL = 'https://ftp.ncbi.nlm.nih.gov/blast/db/', check_update = TRUE)
untar(blast.tar, exdir = './reference/16S_database')
list.files('./reference/16S_database')
blast.db <- blast(db = './reference/16S_database/16S_ribosomal_RNA')

# Performs the blast for each read and returns the best hit # 
soil.hits <- matrix(nrow = nrow(raw_soil$tax), ncol = 12)
soil.hits <- as.data.frame(soil.hits) 
hold <- c()
for(i in 1:length(raw_soil$dna)){
  hold <- predict(blast.db, raw_soil$dna[i])
  soil.hits[i,] <- hold[1,]
  raw_soil$tax$Best_Hit[i] <- hold[1, 2]
}

# Filter out reads that do not correspond to a NCBI entry #
library(dplyr); packageVersion('dplyr')
filt_soil.tax <- filter(raw_soil$tax, !is.na(raw_soil$tax$Best_Hit))

# Output the resulting NCBI entry names to a list #
if(!dir.exists("./blast_hits")){
  dir.create('./blast_hits')
}
write.table(filt_soil.tax$Best_Hit, './blast_hits/soil_blast_hits.txt')

# Call the python script to retrieve the taxonomies of the matched entries #
system('python3 ~/PSF_MBIOME/rRNA_BLAST.py -i ./blast_hits/soil_blast_hits.txt -o ./blast_hits/soil_ncbi_hits.csv')

# Read in the output from the python script and make new taxonomy table "soil_ncbi_fin.tax" #
soil_ncbi.taxa <- read.csv2('./blast_hits/soil_ncbi_hits.csv', header = FALSE, fill = TRUE)
soil_ncbi.int <- strsplit(as.character(soil_ncbi.taxa$V1), ",")
soil_ncbi_fin.tax <- do.call(rbind, lapply(soil_ncbi.int, function(x) { length(x) <- max(sapply(soil_ncbi.int, length)); x }))
soil_ncbi_fin.tax <- as.data.frame(soil_ncbi_fin.tax, stringsAsFactors = FALSE)
rownames(soil_ncbi.taxa) <- rownames(filt_soil.tax)
colnames(soil_ncbi_fin.tax) <- c('Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'hold')
for(i in 1:nrow(soil_ncbi_fin.tax)){
  if(!is.na(soil_ncbi_fin.tax$hold[i])){
    soil_ncbi_fin.tax$Genus[i] <- soil_ncbi_fin.tax$hold[i]
  }
}

# Filter otu table and refseq object such that all reads without a BLAST assignment are removed #
soil_ncbi_fin.tax <- soil_ncbi_fin.tax[,1:7]
filter_soil.tax <- cbind(filt_soil.tax, soil_ncbi_fin.tax)

decompose_ps(raw_soil.ps, 'filt_soil')

filt_soil$tax <- filter_soil.tax
filt_soil$otu <- filter(filt_soil$otu, rownames(filt_soil$otu) %in% rownames(filt_soil$tax))
soil_dna.df <- as.data.frame(filt_soil$dna)
soil_dna.df <- filter(soil_dna.df, rownames(soil_dna.df) %in% rownames(filt_soil$tax))
filt_soil$dna <- DNAStringSet(soil_dna.df$x)
names(filt_soil$dna) <- rownames(filt_soil$tax)
filt_soil$tax <- as.matrix(filt_soil$tax)

# Make phyloseq object with filtered tables #
soil.ps <- phyloseq(otu_table(filt_soil$otu, taxa_are_rows = TRUE),
                        sample_data(filt_soil$met),
                        tax_table(filt_soil$tax),
                        refseq(filt_soil$dna))

soil.ps <- subset_taxa(soil.ps, taxa_sums(soil.ps) > 1000)

# Change the taxa names to represent comparative abundance and lowest identification level #
taxa_names(soil.ps) <- paste0('ASV', seq(ntaxa(soil.ps)))
decompose_ps(soil.ps, 'soil')
for(i in 1:nrow(soil$tax)){
  if(!is.na(soil$tax$Genus[i])){
    taxa_names(soil.ps)[i] = paste0(taxa_names(soil.ps)[i], '(', soil$tax$Genus[i], ')')
  }else if(!is.na(soil$tax$Family[i])){
    taxa_names(soil.ps)[i] = paste0(taxa_names(soil.ps)[i], '(', soil$tax$Family[i], ')')
  }else if(!is.na(soil$tax$Order[i])){
    taxa_names(soil.ps)[i] = paste0(taxa_names(soil.ps)[i], '(', soil$tax$Order[i], ')')
  }else{
    taxa_names(soil.ps)[i] = paste0(taxa_names(soil.ps)[i], '(NA)')
  }
}

# produce final decomposed phyloseq object
decompose_ps(soil.ps, 'soil')

#### Phylogenetic Tree Construction for Soils ####
# Output the reads into a fasta file #
writeXStringSet(soil$dna, "./reads/soil_input.fasta", use.names = TRUE)

# Perform a multiple sequence alignment using MAFFT #
system('mafft --auto --thread -1 ./reads/soil_input.fasta > ./reads/soils_aligned.fasta')

# Construct a tree using IQTree with a general time reversible model with a gamma distribution and invariant site copies #
system('iqtree -s ./reads/soils_aligned.fasta -m GTR+G+I -nt AUTO')

# Read the tree using ape and check that the tip labels match #
library(ape); packageVersion('ape')
soil.tre <- read.tree("./reads/soils_aligned.fasta.treefile")
soil.tre$tip.label <- sub("^(ASV[0-9]+)_([^_]+)_$", "\\1(\\2)", soil.tre$tip.label)
soil.tre$tip.label <- gsub("ASV47_Enterobacter", "ASV47(Enterobacter cloacae complex)", soil.tre$tip.label)
soil.tre$tip.label <- gsub("ASV32_Streptomyces", "ASV32(Streptomyces aurantiacus group)", soil.tre$tip.label)
soil.tre$tip.label <- gsub("ASV53_Streptomyces", "ASV53(Streptomyces aurantiacus group)", soil.tre$tip.label)

# Combine final phyloseq object for soil samples #
decompose_ps(soil.ps, 'soil')
soil$tax$ASV <- rownames(soil$tax)
soil$tax <- as.matrix(soil$tax)
soil.ps <- phyloseq(otu_table(soil$otu, taxa_are_rows = TRUE),
                    sample_data(soil$met),
                    tax_table(soil$tax),
                    refseq(soil$dna),
                    phy_tree(soil.tre))

#### Soil Nodule Stacked Histograms ####
# Create a nodule specifc phyloseq object #
soil_nod.ps <- subset_samples(soil.ps, Compartment == "Nodule")
soil_nod.ps <- subset_taxa(soil_nod.ps, taxa_sums(soil_nod.ps) > 0)
decompose_ps(soil_nod.ps, "soil_nod")

# Create a color palette for each ASV #
library(Polychrome)
soil_nod.colr <- createPalette(ntaxa(soil_nod.ps),  c("#ff0000", "#00ff00", "#0000ff"))
soil_nod.colr <- as.data.frame(soil_nod.colr)
rownames(soil_nod.colr) <- rownames(soil_nod$tax)
soil_nod.colr[118,] <- "#D4D4D4" 
rownames(soil_nod.colr)[118] <- "Other" 

# Create a tree for ASVs that are found in the nodules and are in the Order Hyphomicrobiales #
soil_hyph.ps <- subset_taxa(soil_nod.ps, Order == "Hyphomicrobiales")
plot_tree(soil_hyph.ps, label.tips = "taxa_names", ladderize = TRUE, color = "Family")

# save a new phyloseq object without the nodules #
soil.ps <- subset_samples(soil.ps, Compartment != "Nodule")
soil.ps <- subset_taxa(soil.ps, taxa_sums(soil.ps) > 0)
decompose_ps(soil.ps, 'soil')

# Construct a phyloseq object for each individual plant taxon's nodule community # 
## Fuzzy bean ##
fb_soil_nod.ps <- subset_samples(soil_nod.ps, Plant == "S. helvola")
fb_soil_nod.ps <- subset_taxa(fb_soil_nod.ps, taxa_sums(fb_soil_nod.ps) > 0)
fb_soil_nod.ps <- aggregate_top_taxa2(fb_soil_nod.ps, 8, "ASV")
fb_nod_soil.name <- names(sort(taxa_sums(fb_soil_nod.ps), decreasing = TRUE))
fb_soil_nod.colr <- soil_nod.colr[fb_nod_soil.name,]
fb_soil_nod.df <- psmelt(fb_soil_nod.ps)
fb_soil_nod.df$ASVs <- factor(fb_soil_nod.df$ASV, levels = fb_nod_soil.name)
fb_soil_nod.df$Soil <- factor(fb_soil_nod.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

fb_soil_nod.plot <- ggplot(fb_soil_nod.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = fb_soil_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "S. helvola")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "A.")
fb_soil_nod.plot

## Chamaecrista ##
cc_soil_nod.ps <- subset_samples(soil_nod.ps, Plant == "C. fasciculata")
cc_soil_nod.ps <- subset_taxa(cc_soil_nod.ps, taxa_sums(cc_soil_nod.ps) > 0)
cc_soil_nod.ps <- aggregate_top_taxa2(cc_soil_nod.ps, 8, "ASV")
cc_nod_soil.name <- names(sort(taxa_sums(cc_soil_nod.ps), decreasing = TRUE))
cc_soil_nod.colr <- soil_nod.colr[cc_nod_soil.name,]
cc_soil_nod.df <- psmelt(cc_soil_nod.ps)
cc_soil_nod.df$ASVs <- factor(cc_soil_nod.df$ASV, levels = cc_nod_soil.name)
cc_soil_nod.df$Soil <- factor(cc_soil_nod.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

cc_soil_nod.plot <- ggplot(cc_soil_nod.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = cc_soil_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "C. fasciculata")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "B.")
cc_soil_nod.plot

## Desmodium ##
ds_soil_nod.ps <- subset_samples(soil_nod.ps, Plant == "D. illinoense")
ds_soil_nod.ps <- subset_taxa(ds_soil_nod.ps, taxa_sums(ds_soil_nod.ps) > 0)
ds_nod_soil.name <- names(sort(taxa_sums(ds_soil_nod.ps), decreasing = TRUE))
ds_soil_nod.colr <- soil_nod.colr[ds_nod_soil.name,]
ds_soil_nod.df <- psmelt(ds_soil_nod.ps)
ds_soil_nod.df$ASVs <- factor(ds_soil_nod.df$ASV, levels = ds_nod_soil.name)
ds_soil_nod.df$Soil <- factor(ds_soil_nod.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

ds_soil_nod.plot <- ggplot(ds_soil_nod.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = ds_soil_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "D. illinoense")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "C.")
ds_soil_nod.plot

## Hog Peanut ##
hp_soil_nod.ps <- subset_samples(soil_nod.ps, Plant == "A. bracteata")
hp_soil_nod.ps <- subset_taxa(hp_soil_nod.ps, taxa_sums(hp_soil_nod.ps) > 0)
hp_nod_soil.name <- names(sort(taxa_sums(hp_soil_nod.ps), decreasing = TRUE))
hp_soil_nod.colr <- soil_nod.colr[hp_nod_soil.name,]
hp_soil_nod.df <- psmelt(hp_soil_nod.ps)
hp_soil_nod.df$ASVs <- factor(hp_soil_nod.df$ASV, levels = hp_nod_soil.name)
hp_soil_nod.df$Soil <- factor(hp_soil_nod.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

hp_soil_nod.plot <- ggplot(hp_soil_nod.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = hp_soil_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "A. bracteata")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "D.")
hp_soil_nod.plot

## Clover ##
cl_soil_nod.ps <- subset_samples(soil_nod.ps, Plant == "T. repens")
cl_soil_nod.ps <- subset_taxa(cl_soil_nod.ps, taxa_sums(cl_soil_nod.ps) > 0)
cl_soil_nod.ps <- aggregate_top_taxa2(cl_soil_nod.ps, 8, "ASV")
cl_nod_soil.name <- names(sort(taxa_sums(cl_soil_nod.ps), decreasing = TRUE))
cl_soil_nod.colr <- soil_nod.colr[cl_nod_soil.name,]
cl_soil_nod.df <- psmelt(cl_soil_nod.ps)
cl_soil_nod.df$ASVs <- factor(cl_soil_nod.df$ASV, levels = cl_nod_soil.name)
cl_soil_nod.df$Soil <- factor(cl_soil_nod.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

cl_soil_nod.plot <- ggplot(cl_soil_nod.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = cl_soil_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "T. repens")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "F.")
cl_soil_nod.plot

## Medicago ##
md_soil_nod.ps <- subset_samples(soil_nod.ps, Plant == "M. truncatula")
md_soil_nod.ps <- subset_taxa(md_soil_nod.ps, taxa_sums(md_soil_nod.ps) > 0)
md_soil_nod.ps <- aggregate_top_taxa2(md_soil_nod.ps, 8, "ASV")
md_nod_soil.name <- names(sort(taxa_sums(md_soil_nod.ps), decreasing = TRUE))
md_soil_nod.colr <- soil_nod.colr[md_nod_soil.name,]
md_soil_nod.df <- psmelt(md_soil_nod.ps)
md_soil_nod.df$ASVs <- factor(md_soil_nod.df$ASV, levels = md_nod_soil.name)
md_soil_nod.df$Soil <- factor(md_soil_nod.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

md_soil_nod.plot <- ggplot(md_soil_nod.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = md_soil_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "M. truncatula")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "G.")
md_soil_nod.plot

#### Root Primer Removal ####
# Ensure you have the right files #
list.files(root.dir)

# Create a list of the files that are corresponding to the forward and reverse reads #
raw_root.ffp <- sort(list.files(root.dir, pattern = "_R1_001.fastq.gz", full.names = TRUE))
raw_root.rfp <- sort(list.files(root.dir, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Save the names based on the file names #
root.names <- sub("^([^_]+_[^_]+)_.*$", "\\1", basename(raw_root.ffp))
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
root.rori <- allOrients(root.rprimer)

# Make filepaths for pretrimmed fastqs #
system('mkdir ./reads/pretrim/root_pretrim')
pre_root.ffp <- file.path('./reads/pretrim/root_pretrim', paste0(root.names, '_pretrim_R1.fastq.gz'))
pre_root.rfp <- file.path('./reads/pretrim/root_pretrim', paste0(root.names, '_pretrim_R2.fastq.gz'))

# Filter reads less than 75 bp and save the filtered fastqs to the pretrim filepaths #
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

save.image("./test.RData")
#### dada2 Implementation for the Root Samples ####
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

save.image("./test.RData")
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

# track reads through the pipeline #
getN <- function(x) sum(getUniques(x))
root_final.track <- cbind(root_prefilt.track[,1], root_prefilt.track[,2], root_postfilt.track[,2], sapply(root.fdada, getN), sapply(root.rdada, getN), sapply(root.remerged, getN), colSums(root_nochim.st))
colnames(root_final.track) <- c("pre-cutadapt", "post-cutadapt", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(root_final.track) <- root.names
root_final.track <- as.data.frame(root_final.track)

# Assign Taxonomy #
root_rdp.taxa <- assignTaxonomy(rownames(root_nochim.st), refFasta = reference, multithread = TRUE, verbose = TRUE)
root_rdp.taxa <- as.matrix(root_rdp.taxa)

# Load the metadata #
root_raw.met <- read.csv2(root.met, sep = ',')
rownames(root_raw.met) <- root_raw.met$Sample
root_raw.met <- root_raw.met[,c('Sample.Name', 'Plant.Species', 'Soil.Origin', 'Compartment')]
rownames(root_raw.met) <- sub("^([^_]+_[^_]+)_.*$", "\\1", rownames(root_raw.met))
colnames(root_nochim.st) <- rownames(root_raw.met)

#### Phyloseq Object Construction and Filtering for roots ####
root_rdp.taxa <- as.matrix(root_rdp.taxa)
raw_root.ps <- phyloseq(otu_table(root_nochim.st, taxa_are_rows = TRUE),
                        sample_data(root_raw.met),
                        tax_table(root_rdp.taxa))

raw_root.dna <- Biostrings::DNAStringSet(taxa_names(raw_root.ps))
names(raw_root.dna) <- taxa_names(raw_root.ps)
raw_root.ps <- merge_phyloseq(raw_root.ps, raw_root.dna)

raw_root.ps <- subset_taxa(raw_root.ps, Phylum != "Cyanobacteriota")
raw_root.ps <- subset_taxa(raw_root.ps, Phylum != "Plantae")

raw_root.ps <- subset_samples(raw_root.ps, Compartment != 'Leaf Endosphere')

raw_root.ps <- subset_taxa(raw_root.ps, taxa_sums(raw_root.ps) > 100)
colnames(tax_table(raw_root.ps)) <- c('rdp_Kingdom', 'rdp_Phylum', 'rdp_Class', 'rdp_Order', 'rdp_Family', 'rdp_Genus')

taxa_names(raw_root.ps) <- paste0("ASV", seq(ntaxa(raw_root.ps)))
decompose_ps(raw_root.ps, 'raw_root')
#### Cross-Validation of root Reads Using BLAST ####

# Performs the blast for each read and returns the best hit # 
root.hits <- matrix(nrow = nrow(raw_root$tax), ncol = 12)
root.hits <- as.data.frame(root.hits) 
hold <- c()
for(i in 1:length(raw_root$dna)){
  hold <- predict(blast.db, raw_root$dna[i])
  root.hits[i,] <- hold[1,]
  raw_root$tax$Best_Hit[i] <- hold[1, 2]
}

# Filter out reads that do not correspond to a NCBI entry #
filt_root.tax <- filter(raw_root$tax, !is.na(raw_root$tax$Best_Hit))

# Output the resulting NCBI entry names to a list #
if(!dir.exists("./blast_hits")){
  dir.create('./blast_hits')
}
write.table(filt_root.tax$Best_Hit, './blast_hits/root_blast_hits.txt')

# Call the python script to retrieve the taxonomies of the matched entries #
system('python3 ~/PSF_MBIOME/rRNA_BLAST.py -i ./blast_hits/root_blast_hits.txt -o ./blast_hits/root_ncbi_hits.csv')

# Read in the output from the python script and make new taxonomy table "root_ncbi_fin.tax" #
root_ncbi.taxa <- read.csv2('./blast_hits/root_ncbi_hits.csv', header = FALSE, fill = TRUE)
root_ncbi.int <- strsplit(as.character(root_ncbi.taxa$V1), ",")
root_ncbi_fin.tax <- do.call(rbind, lapply(root_ncbi.int, function(x) { length(x) <- max(sapply(root_ncbi.int, length)); x }))
root_ncbi_fin.tax <- as.data.frame(root_ncbi_fin.tax, stringsAsFactors = FALSE)
rownames(root_ncbi.taxa) <- rownames(filt_root.tax)
colnames(root_ncbi_fin.tax) <- c('Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'hold')
for(i in 1:nrow(root_ncbi_fin.tax)){
  if(!is.na(root_ncbi_fin.tax$hold[i])){
    root_ncbi_fin.tax$Genus[i] <- root_ncbi_fin.tax$hold[i]
  }
}

# Filter otu table and refseq object such that all reads without a BLAST assignment are removed #
root_ncbi_fin.tax <- root_ncbi_fin.tax[,1:7]
filter_root.tax <- cbind(filt_root.tax, root_ncbi_fin.tax)

decompose_ps(raw_root.ps, 'filt_root')

filt_root$tax <- filter_root.tax
filt_root$otu <- filter(filt_root$otu, rownames(filt_root$otu) %in% rownames(filt_root$tax))
root_dna.df <- as.data.frame(filt_root$dna)
root_dna.df <- filter(root_dna.df, rownames(root_dna.df) %in% rownames(filt_root$tax))
filt_root$dna <- DNAStringSet(root_dna.df$x)
names(filt_root$dna) <- rownames(filt_root$tax)
filt_root$tax <- as.matrix(filt_root$tax)

# Make phyloseq object with filtered tables #
root.ps <- phyloseq(otu_table(filt_root$otu, taxa_are_rows = TRUE),
                    sample_data(filt_root$met),
                    tax_table(filt_root$tax),
                    refseq(filt_root$dna))

root.ps <- subset_taxa(root.ps, taxa_sums(root.ps) > 100)

# Change the taxa names to represent comparatiove abundance and lowest identification level #
decompose_ps(root.ps, 'root')
taxa_names(root.ps) <- paste0('ASV', seq(ntaxa(root.ps)))
for(i in 1:nrow(root$tax)){
  if(!is.na(root$tax$Genus[i])){
    taxa_names(root.ps)[i] = paste0(taxa_names(root.ps)[i], '(', root$tax$Genus[i], ')')
  }else if(!is.na(root$tax$Family[i])){
    taxa_names(root.ps)[i] = paste0(taxa_names(root.ps)[i], '(', root$tax$Family[i], ')')
  }else if(!is.na(root$tax$Order[i])){
    taxa_names(root.ps)[i] = paste0(taxa_names(root.ps)[i], '(', root$tax$Order[i], ')')
  }else{
    taxa_names(root.ps)[i] = paste0(taxa_names(root.ps)[i], '(NA)')
  }
}

# produce final decomposed phyloseq object
decompose_ps(root.ps, 'root')

#### Phylogenetic Tree Construction for roots ####
# Output the reads into a fasta file #
writeXStringSet(as.character(root$dna, "./reads/root_input.fasta", use.names = TRUE))

# Perform a multiple sequence alignment using MAFFT #
system('mafft --auto --thread -1 ./reads/root_input.fasta > ./reads/roots_aligned.fasta')

# Construct a tree using IQTree with a general time reversible model with a gamma distribution and invariant site copies #
system('iqtree -s ./reads/roots_aligned.fasta -m GTR+G+I -nt AUTO')

# Read the tree using ape and check that the tip labels match #
root.tre <- read.tree("./reads/roots_aligned.fasta.treefile")
root.tre$tip.label <- sub("^(ASV[0-9]+)_([^_]+)_$", "\\1(\\2)", root.tre$tip.label)
root.tre$tip.label <- gsub("_Enterobacter", "(Enterobacter cloacae complex)", root.tre$tip.label)
root.tre$tip.label <- gsub("ASV110_Streptomyces", "ASV110(Streptomyces phaeochromogenes group)", root.tre$tip.label)
root.tre$tip.label <- gsub("_Streptomyces", "(Streptomyces aurantiacus group)", root.tre$tip.label)
root.tre$tip.label <- gsub("_Bacillus", "(Bacillus cereus group)", root.tre$tip.label)
root.tre$tip.label <- gsub("_Citrobacter", "(Citrobacter freundii complex)", root.tre$tip.label)
root.tre$tip.label <- gsub("_Pantoea", "(Pantoea agglomerans group)", root.tre$tip.label)
root.tre$tip.label <- gsub("_Candidatus", "(Candidatus Protochlamydia)", root.tre$tip.label)



# Combine final phyloseq object for root samples #
decompose_ps(root.ps, 'root')
root$tax$ASV <- rownames(root$tax)
root$tax <- as.matrix(root$tax)
root.ps <- phyloseq(otu_table(root$otu, taxa_are_rows = TRUE),
                    sample_data(root$met),
                    tax_table(root$tax),
                    refseq(root$dna),
                    phy_tree(root.tre))

#### Root Nodule Stacked Histograms ####
# Create a nodule specifc phyloseq object #
root_nod.ps <- subset_samples(root.ps, Compartment == "Nodule")
root_nod.ps <- subset_taxa(root_nod.ps, taxa_sums(root_nod.ps) > 0)
decompose_ps(root_nod.ps, "root_nod")

# Create a color palette for each ASV #
set.seed(248)
root_nod.colr <- createPalette(ntaxa(root_nod.ps),  c("#ff0000", "#00ff00", "#0000ff"))
root_nod.colr <- as.data.frame(root_nod.colr)
rownames(root_nod.colr) <- rownames(root_nod$tax)
root_nod.colr[252,] <- "#D4D4D4" 
rownames(root_nod.colr)[252] <- "Other" 

# Create a tree for ASVs that are found in the nodules and are in the Order Hyphomicrobiales #
root_hyph.ps <- subset_taxa(root_nod.ps, Order == "Hyphomicrobiales")
plot_tree(root_hyph.ps, label.tips = "taxa_names", ladderize = TRUE, color = "Family")

# save a new phyloseq object without the nodules #
root.ps <- subset_samples(root.ps, Compartment != "Nodule")
root.ps <- subset_taxa(root.ps, taxa_sums(root.ps) > 0)
decompose_ps(root.ps, 'root')

# Construct a phyloseq object for each individual plant taxon's nodule community # 
## Fuzzy bean ##
fb_root_nod.ps <- subset_samples(root_nod.ps, Plant.Species == "S. helvola")
fb_root_nod.ps <- subset_taxa(fb_root_nod.ps, taxa_sums(fb_root_nod.ps) > 0)
fb_root_nod.ps <- aggregate_top_taxa2(fb_root_nod.ps, 8, "ASV")
fb_nod_root.name <- names(sort(taxa_sums(fb_root_nod.ps), decreasing = TRUE))
fb_root_nod.colr <- root_nod.colr[fb_nod_root.name,]
fb_root_nod.df <- psmelt(fb_root_nod.ps)
fb_root_nod.df$ASVs <- factor(fb_root_nod.df$ASV, levels = fb_nod_root.name)
fb_root_nod.df$Soil <- factor(fb_root_nod.df$Soil.Origin, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

fb_root_nod.plot <- ggplot(fb_root_nod.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = fb_root_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "S. helvola")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "A.")
fb_root_nod.plot

(fb_soil_nod.plot | fb_root_nod.plot)

## Chamaecrista ##
cc_root_nod.ps <- subset_samples(root_nod.ps, Plant.Species == "C. fasciculata")
cc_root_nod.ps <- subset_taxa(cc_root_nod.ps, taxa_sums(cc_root_nod.ps) > 0)
cc_root_nod.ps <- aggregate_top_taxa2(cc_root_nod.ps, 8, "ASV")
cc_nod_root.name <- names(sort(taxa_sums(cc_root_nod.ps), decreasing = TRUE))
cc_root_nod.colr <- root_nod.colr[cc_nod_root.name,]
cc_root_nod.df <- psmelt(cc_root_nod.ps)
cc_root_nod.df$ASVs <- factor(cc_root_nod.df$ASV, levels = cc_nod_root.name)
cc_root_nod.df$Soil <- factor(cc_root_nod.df$Soil.Origin, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

cc_root_nod.plot <- ggplot(cc_root_nod.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = cc_root_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "C. fasciculata")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "B.")
cc_root_nod.plot

(cc_soil_nod.plot | cc_root_nod.plot)

## Desmodium ##
ds_root_nod.ps <- subset_samples(root_nod.ps, Plant.Species == "D. illinoense")
ds_root_nod.ps <- subset_taxa(ds_root_nod.ps, taxa_sums(ds_root_nod.ps) > 0)
ds_root_nod.ps <- aggregate_top_taxa2(ds_root_nod.ps, 8, "ASV")
ds_nod_root.name <- names(sort(taxa_sums(ds_root_nod.ps), decreasing = TRUE))
ds_root_nod.colr <- root_nod.colr[ds_nod_root.name,]
ds_root_nod.df <- psmelt(ds_root_nod.ps)
ds_root_nod.df$ASVs <- factor(ds_root_nod.df$ASV, levels = ds_nod_root.name)
ds_root_nod.df$Soil <- factor(ds_root_nod.df$Soil.Origin, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

ds_root_nod.plot <- ggplot(ds_root_nod.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = ds_root_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "D. illinoense")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "C.")
ds_root_nod.plot

(ds_soil_nod.plot | ds_root_nod.plot)

## Hog Peanut ##
hp_root_nod.ps <- subset_samples(root_nod.ps, Plant.Species == "A. bracteata")
hp_root_nod.ps <- subset_taxa(hp_root_nod.ps, taxa_sums(hp_root_nod.ps) > 0)
hp_nod_root.name <- names(sort(taxa_sums(hp_root_nod.ps), decreasing = TRUE))
hp_root_nod.colr <- root_nod.colr[hp_nod_root.name,]
hp_root_nod.df <- psmelt(hp_root_nod.ps)
hp_root_nod.df$ASVs <- factor(hp_root_nod.df$ASV, levels = hp_nod_root.name)
hp_root_nod.df$Soil <- factor(hp_root_nod.df$Soil.Origin, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

hp_root_nod.plot <- ggplot(hp_root_nod.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = hp_root_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "A. bracteata")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "D.")
hp_root_nod.plot

(hp_soil_nod.plot | hp_root_nod.plot)

## Clover ##
cl_root_nod.ps <- subset_samples(root_nod.ps, Plant.Species == "T. repens")
cl_root_nod.ps <- subset_taxa(cl_root_nod.ps, taxa_sums(cl_root_nod.ps) > 0)
cl_root_nod.ps <- aggregate_top_taxa2(cl_root_nod.ps, 8, "ASV")
cl_nod_root.name <- names(sort(taxa_sums(cl_root_nod.ps), decreasing = TRUE))
cl_root_nod.colr <- root_nod.colr[cl_nod_root.name,]
cl_root_nod.df <- psmelt(cl_root_nod.ps)
cl_root_nod.df$ASVs <- factor(cl_root_nod.df$ASV, levels = cl_nod_root.name)
cl_root_nod.df$Soil<- factor(cl_root_nod.df$Soil.Origin, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

cl_root_nod.plot <- ggplot(cl_root_nod.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = cl_root_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "T. repens")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "F.")
cl_root_nod.plot

(cl_soil_nod.plot | cl_root_nod.plot)

## Medicago ##
md_root_nod.ps <- subset_samples(root_nod.ps, Plant.Species == "M. truncatula")
md_root_nod.ps <- subset_taxa(md_root_nod.ps, taxa_sums(md_root_nod.ps) > 0)
md_root_nod.ps <- aggregate_top_taxa2(md_root_nod.ps, 8, "ASV")
md_nod_root.name <- names(sort(taxa_sums(md_root_nod.ps), decreasing = TRUE))
md_root_nod.colr <- root_nod.colr[md_nod_root.name,]
md_root_nod.df <- psmelt(md_root_nod.ps)
md_root_nod.df$ASVs <- factor(md_root_nod.df$ASV, levels = md_nod_root.name)
md_root_nod.df$Soil <- factor(md_root_nod.df$Soil.Origin, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

md_root_nod.plot <- ggplot(md_root_nod.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = md_root_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "M. truncatula")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "G.")
md_root_nod.plot

(md_soil_nod.plot | md_root_nod.plot)

#### Alpha Diversity Analysis and Visualization ####
# Add additional variables as factors for Soils #
soil$met$Soils <- factor(soil$met$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))
soil$met$Soils <- gsub("Non-PSF Soil", "Non PSF Soil", soil$met$Soils)
soil$met$Comps <- factor(soil$met$Compartment, levels = c("Bulk Soil", "Rhizosphere"))
soil$met$Plants <- factor(soil$met$Plant, levels = c("T. repens", "M. truncatula", "S. helvola", "C. fasciculata", "D. illinoense", "A. bracteata"))
for(i in 1:nrow(soil$met)){
  soil$met$Tri[i] <- paste0(substr(soil$met$Plant[i],1,1), substr(soil$met$Soil_Treatment[i],1,1), substr(soil$met$Compartment[i],1,2)) 
}

# Add additional variables as factors for Soils #
root$met$Soils <- factor(root$met$Soil.Origin, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))
root$met$Soils <- gsub("Non-PSF Soil", "Non PSF Soil", root$met$Soils)
root$met$Comps <- factor(root$met$Compartment, levels = "Root Endosphere")
root$met$Plants <- factor(root$met$Plant.Species, levels = c("T. repens", "M. truncatula", "S. helvola", "C. fasciculata", "D. illinoense", "A. bracteata"))
for(i in 1:nrow(root$met)){
  root$met$Tri[i] <- paste0(substr(root$met$Plant.Species[i],1,1), substr(root$met$Soil.Origin[i],1,1), substr(root$met$Compartment[i],1,2)) 
}
soil$met <- soil$met[,c("Plant", "Soil_Treatment", "Compartment", "Soils", "Comps", "Plants", "Tri")]
root$met <- root$met[,c("Plant.Species", "Soil.Origin", "Compartment", "Soils", "Comps", "Plants", "Tri")]
colnames(root$met) <- colnames(soil$met)

# Calculate the diversitry metrics #
soil.rich <- estimate_richness(soil.ps)
soil.rich <- as.data.frame(soil.rich)
soil.rich <- cbind(soil$met, soil.rich)
root.rich <- estimate_richness(root.ps)
root.rich <- as.data.frame(root.rich)
root.rich <- cbind(root$met, root.rich)
all.rich <- rbind(soil.rich, root.rich)

# Calculate Shannon Evenness by dividing the Shannon Diversity value by the  #
for(i in 1:nrow(all.rich)){
  all.rich$ShaEvn[i] <- all.rich$Shannon[i]/log(all.rich$Chao1[i], base = 2.718) 
}


# Perform three way anova tests for plant species, compartment, and soil origin #
all_cha.aov <- aov(Chao1~Soils*Plants*Comps, all.rich)
summary(all_cha.aov)
all_evn.aov <- aov(ShaEvn~Soils*Plants*Comps, all.rich)
summary(all_evn.aov)
all_sha.aov <- aov(Shannon~Soils*Plants*Comps, all.rich)
summary(all_sha.aov)

# Find the mean and standard error for each unique tripartite grouping
all_rich.mnsd <- all.rich %>%
  group_by(Plant, Soil_Treatment, Compartment) %>%
  summarize(
    sha.mean = mean(Shannon),
    evn.mean = mean(ShaEvn),
    cha.mean = mean(Chao1),
    sha.sd = sd(Shannon),
    evn.sd = sd(ShaEvn),
    cha.sd = sd(Chao1),
    .groups = "drop" # Prevent grouping in the result
  )

all_rich.mnsd <- as.data.frame(all_rich.mnsd)
all_rich.mnsd$`Plant Species` <- factor(all_rich.mnsd$Plant, levels = c('A. bracteata', 'D. illinoense', 'C. fasciculata', 'S. helvola', 'T. repens', 'M. truncatula'))
all_rich.mnsd$`Soil Type` <- factor(all_rich.mnsd$Soil_Treatment, levels = c('Common Soil', 'Non-PSF Soil', 'PSF Soil'))
all_rich.mnsd$`Soil Type` <- gsub("Non-PSF Soil", "Non PSF Soil", all_rich.mnsd$`Soil Type`)

# All Rhizosphere Samples #
rhiz.rich <- filter(all_rich.mnsd, Compartment == "Rhizosphere")
rhiz.richraw <- filter(all.rich, Compartment == 'Rhizosphere')

## Rhizosphere Shannon Diversity ##
### Strophostyles ###
fb_rhiz.richraw <- filter(rhiz.richraw, Plant == 'S. helvola')
fb_rhiz_sha.aov <- aov(Shannon ~ Soils, data = fb_rhiz.richraw)
summary(fb_rhiz_sha.aov)

fb_rhiz_sha.hsd <- TukeyHSD(fb_rhiz_sha.aov)
fb_rhiz_sha.hsd

fb_rhiz_sha.let <- multcompLetters4(fb_rhiz_sha.aov, fb_rhiz_sha.hsd)
fb_rhiz_sha.let <- fb_rhiz_sha.let$Soils$Letters[sort(names(fb_rhiz_sha.let$Soils$Letters))]

## Chamecrista ##
cc_rhiz.richraw <- filter(rhiz.richraw, Plant == 'C. fasciculata')
cc_rhiz_sha.aov <- aov(Shannon ~ Soils, data = cc_rhiz.richraw)
summary(cc_rhiz_sha.aov)

cc_rhiz_sha.hsd <- TukeyHSD(cc_rhiz_sha.aov)
cc_rhiz_sha.hsd

cc_rhiz_sha.let <- multcompLetters4(cc_rhiz_sha.aov, cc_rhiz_sha.hsd)
cc_rhiz_sha.let <- cc_rhiz_sha.let$Soils$Letters[sort(names(cc_rhiz_sha.let$Soils$Letters))]


## Desmodium ##
ds_rhiz.richraw <- filter(rhiz.richraw, Plant == 'D. illinoense')
ds_rhiz_sha.aov <- aov(Shannon ~ Soils, data = ds_rhiz.richraw)
summary(ds_rhiz_sha.aov)

ds_rhiz_sha.hsd <- TukeyHSD(ds_rhiz_sha.aov)
ds_rhiz_sha.hsd

ds_rhiz_sha.let <- multcompLetters4(ds_rhiz_sha.aov, ds_rhiz_sha.hsd)
ds_rhiz_sha.let <- ds_rhiz_sha.let$Soils$Letters[sort(names(ds_rhiz_sha.let$Soils$Letters))]

## Amphicarpaea ##
hp_rhiz.richraw <- filter(rhiz.richraw, Plant == 'A. bracteata')
hp_rhiz_sha.aov <- aov(Shannon ~ Soils, data = hp_rhiz.richraw)
summary(hp_rhiz_sha.aov)

hp_rhiz_sha.hsd <- TukeyHSD(hp_rhiz_sha.aov)
hp_rhiz_sha.hsd

hp_rhiz_sha.let <- multcompLetters4(hp_rhiz_sha.aov, hp_rhiz_sha.hsd)
hp_rhiz_sha.let <- hp_rhiz_sha.let$Soils$Letters[sort(names(hp_rhiz_sha.let$Soils$Letters))]

## Trifolium ##
cl_rhiz.richraw <- filter(rhiz.richraw, Plant == 'T. repens')
cl_rhiz_sha.aov <- aov(Shannon ~ Soils, data = cl_rhiz.richraw)
summary(cl_rhiz_sha.aov)

cl_rhiz_sha.hsd <- TukeyHSD(cl_rhiz_sha.aov)
cl_rhiz_sha.hsd

cl_rhiz_sha.let <- multcompLetters4(cl_rhiz_sha.aov, cl_rhiz_sha.hsd)
cl_rhiz_sha.let <- cl_rhiz_sha.let$Soils$Letters[sort(names(cl_rhiz_sha.let$Soils$Letters))]

## Medicago ##
md_rhiz.richraw <- filter(rhiz.richraw, Plant == 'M. truncatula')
md_rhiz_sha.aov <- aov(Shannon ~ Soils, data = md_rhiz.richraw)
summary(md_rhiz_sha.aov)

md_rhiz_sha.hsd <- TukeyHSD(md_rhiz_sha.aov)
md_rhiz_sha.hsd

md_rhiz_sha.let <- multcompLetters4(md_rhiz_sha.aov, md_rhiz_sha.hsd)
md_rhiz_sha.let <- md_rhiz_sha.let$Soils$Letters[sort(names(md_rhiz_sha.let$Soils$Letters))]

## Adding Letters ##
sha_rhiz.let <- c(hp_rhiz_sha.let,
                  cc_rhiz_sha.let,
                  ds_rhiz_sha.let,
                  md_rhiz_sha.let,
                  fb_rhiz_sha.let,
                  cl_rhiz_sha.let)
rhiz.rich <- cbind(rhiz.rich, sha_rhiz.let)

### Final Plot ###
sha_rhiz.plot <- ggplot(rhiz.rich, aes(x = `Plant Species`, y = sha.mean, fill = `Soil Type`, color = `Soil Type`)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = sha.mean - sha.sd, ymax = sha.mean + sha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = sha_rhiz.let, y = sha.mean + sha.sd + 0.1), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 4) +
  ylab('Shannon Diversity (H)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,6),expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ .,name = "Rhizosphere")) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm'))

sha_rhiz.plot

# Rhizosphere Shannon Evenness #
## Strophostyles ##
fb_rhiz.richraw <- filter(rhiz.richraw, Plant == 'S. helvola')
fb_rhiz_evn.aov <- aov(ShaEvn ~ Soils, data = fb_rhiz.richraw)
summary(fb_rhiz_evn.aov)

fb_rhiz_evn.hsd <- TukeyHSD(fb_rhiz_evn.aov)
fb_rhiz_evn.hsd

fb_rhiz_evn.let <- multcompLetters4(fb_rhiz_evn.aov, fb_rhiz_evn.hsd)
fb_rhiz_evn.let <- fb_rhiz_evn.let$Soils$Letters[sort(names(fb_rhiz_evn.let$Soils$Letters))]

## Chamecrista ##
cc_rhiz.richraw <- filter(rhiz.richraw, Plant == 'C. fasciculata')
cc_rhiz_evn.aov <- aov(ShaEvn ~ Soils, data = cc_rhiz.richraw)
summary(cc_rhiz_evn.aov)

cc_rhiz_evn.hsd <- TukeyHSD(cc_rhiz_evn.aov)
cc_rhiz_evn.hsd

cc_rhiz_evn.let <- multcompLetters4(cc_rhiz_evn.aov, cc_rhiz_evn.hsd)
cc_rhiz_evn.let <- cc_rhiz_evn.let$Soils$Letters[sort(names(cc_rhiz_evn.let$Soils$Letters))]

## Desmodium ##
ds_rhiz.richraw <- filter(rhiz.richraw, Plant == 'D. illinoense')
ds_rhiz_evn.aov <- aov(ShaEvn ~ Soils, data = ds_rhiz.richraw)
summary(ds_rhiz_evn.aov)

ds_rhiz_evn.hsd <- TukeyHSD(ds_rhiz_evn.aov)
ds_rhiz_evn.hsd

ds_rhiz_evn.let <- multcompLetters4(ds_rhiz_evn.aov, ds_rhiz_evn.hsd)
ds_rhiz_evn.let <- ds_rhiz_evn.let$Soils$Letters[sort(names(ds_rhiz_evn.let$Soils$Letters))]

## Amphicarpaea ##
hp_rhiz.richraw <- filter(rhiz.richraw, Plant == 'A. bracteata')
hp_rhiz_evn.aov <- aov(ShaEvn ~ Soils, data = hp_rhiz.richraw)
summary(hp_rhiz_evn.aov)

hp_rhiz_evn.hsd <- TukeyHSD(hp_rhiz_evn.aov)
hp_rhiz_evn.hsd

hp_rhiz_evn.let <- multcompLetters4(hp_rhiz_evn.aov, hp_rhiz_evn.hsd)
hp_rhiz_evn.let <- hp_rhiz_evn.let$Soils$Letters[sort(names(hp_rhiz_evn.let$Soils$Letters))]

## Trifolium ##
cl_rhiz.richraw <- filter(rhiz.richraw, Plant == 'T. repens')
cl_rhiz_evn.aov <- aov(ShaEvn ~ Soils, data = cl_rhiz.richraw)
summary(cl_rhiz_evn.aov)

cl_rhiz_evn.hsd <- TukeyHSD(cl_rhiz_evn.aov)
cl_rhiz_evn.hsd

cl_rhiz_evn.let <- multcompLetters4(cl_rhiz_evn.aov, cl_rhiz_evn.hsd)
cl_rhiz_evn.let <- cl_rhiz_evn.let$Soils$Letters[sort(names(cl_rhiz_evn.let$Soils$Letters))]

## Medicago ##
md_rhiz.richraw <- filter(rhiz.richraw, Plant == 'M. truncatula')
md_rhiz_evn.aov <- aov(ShaEvn ~ Soils, data = md_rhiz.richraw)
summary(md_rhiz_evn.aov)

md_rhiz_evn.hsd <- TukeyHSD(md_rhiz_evn.aov)
md_rhiz_evn.hsd

md_rhiz_evn.let <- multcompLetters4(md_rhiz_evn.aov, md_rhiz_evn.hsd)
md_rhiz_evn.let <- md_rhiz_evn.let$Soils$Letters[sort(names(md_rhiz_evn.let$Soils$Letters))]

## Adding Letters ##
evn_rhiz.let <- c(hp_rhiz_evn.let,
                  cc_rhiz_evn.let,
                  ds_rhiz_evn.let,
                  md_rhiz_evn.let,
                  fb_rhiz_evn.let,
                  cl_rhiz_evn.let)
rhiz.rich <- cbind(rhiz.rich, evn_rhiz.let)

### Final Plot ###
evn_rhiz.plot <- ggplot(rhiz.rich, aes(x = `Plant Species`, y = evn.mean, fill = `Soil Type`, color = `Soil Type`)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = evn.mean - evn.sd, ymax = evn.mean + evn.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = evn_rhiz.let, y = evn.mean + evn.sd + 0.02), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 4) +
  ylab('Shannon Evenness (H/ln(S))') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,1),expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = "")) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm'))

evn_rhiz.plot
# Chao1 Observed ASVs #
## Strophostyles ##
fb_rhiz.richraw <- filter(rhiz.richraw, Plant == 'S. helvola')
fb_rhiz_cha.aov <- aov(Chao1 ~ Soils, data = fb_rhiz.richraw)
summary(fb_rhiz_cha.aov)

fb_rhiz_cha.hsd <- TukeyHSD(fb_rhiz_cha.aov)
fb_rhiz_cha.hsd

fb_rhiz_cha.let <- multcompLetters4(fb_rhiz_cha.aov, fb_rhiz_cha.hsd)
fb_rhiz_cha.let <- fb_rhiz_cha.let$Soils$Letters[sort(names(fb_rhiz_cha.let$Soils$Letters))]

## Chamecrista ##
cc_rhiz.richraw <- filter(rhiz.richraw, Plant == 'C. fasciculata')
cc_rhiz_cha.aov <- aov(Chao1~ Soils, data = cc_rhiz.richraw)
summary(cc_rhiz_cha.aov)

cc_rhiz_cha.hsd <- TukeyHSD(cc_rhiz_cha.aov)
cc_rhiz_cha.hsd

cc_rhiz_cha.let <- multcompLetters4(cc_rhiz_cha.aov, cc_rhiz_cha.hsd)
cc_rhiz_cha.let <- cc_rhiz_cha.let$Soils$Letters[sort(names(cc_rhiz_cha.let$Soils$Letters))]

## Desmodium ##
ds_rhiz.richraw <- filter(rhiz.richraw, Plant == 'D. illinoense')
ds_rhiz_cha.aov <- aov(Chao1 ~ Soils, data = ds_rhiz.richraw)
summary(ds_rhiz_cha.aov)

ds_rhiz_cha.hsd <- TukeyHSD(ds_rhiz_cha.aov)
ds_rhiz_cha.hsd

ds_rhiz_cha.let <- multcompLetters4(ds_rhiz_cha.aov, ds_rhiz_cha.hsd)
ds_rhiz_cha.let <- ds_rhiz_cha.let$Soils$Letters[sort(names(ds_rhiz_cha.let$Soils$Letters))]

## Amphicarpaea ##
hp_rhiz.richraw <- filter(rhiz.richraw, Plant == 'A. bracteata')
hp_rhiz_cha.aov <- aov(Chao1 ~ Soils, data = hp_rhiz.richraw)
summary(hp_rhiz_cha.aov)

hp_rhiz_cha.hsd <- TukeyHSD(hp_rhiz_cha.aov)
hp_rhiz_cha.hsd

hp_rhiz_cha.let <- multcompLetters4(hp_rhiz_cha.aov, hp_rhiz_cha.hsd)
hp_rhiz_cha.let <- hp_rhiz_cha.let$Soils$Letters[sort(names(hp_rhiz_cha.let$Soils$Letters))]

## Trifolium ##
cl_rhiz.richraw <- filter(rhiz.richraw, Plant == 'T. repens')
cl_rhiz_cha.aov <- aov(Chao1 ~ Soils, data = cl_rhiz.richraw)
summary(cl_rhiz_cha.aov)

cl_rhiz_cha.hsd <- TukeyHSD(cl_rhiz_cha.aov)
cl_rhiz_cha.hsd

cl_rhiz_cha.let <- multcompLetters4(cl_rhiz_cha.aov, cl_rhiz_cha.hsd)
cl_rhiz_cha.let <- cl_rhiz_cha.let$Soils$Letters[sort(names(cl_rhiz_cha.let$Soils$Letters))]

## Medicago ##
md_rhiz.richraw <- filter(rhiz.richraw, Plant == 'M. truncatula')
md_rhiz_cha.aov <- aov(Chao1 ~ Soils, data = md_rhiz.richraw)
summary(md_rhiz_cha.aov)

md_rhiz_cha.hsd <- TukeyHSD(md_rhiz_cha.aov)
md_rhiz_cha.hsd

md_rhiz_cha.let <- multcompLetters4(md_rhiz_cha.aov, md_rhiz_cha.hsd)
md_rhiz_cha.let <- md_rhiz_cha.let$Soils$Letters[sort(names(md_rhiz_cha.let$Soils$Letters))]

## Adding Letters ##
cha_rhiz.let <- c(hp_rhiz_cha.let,
                  cc_rhiz_cha.let,
                  ds_rhiz_cha.let,
                  md_rhiz_cha.let,
                  fb_rhiz_cha.let,
                  cl_rhiz_cha.let)
rhiz.rich <- cbind(rhiz.rich, cha_rhiz.let)

### Final Plot ###
cha_rhiz.plot <- ggplot(rhiz.rich, aes(x = `Plant Species`, y = cha.mean, fill = `Soil Type`, color = `Soil Type`)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = cha.mean - cha.sd, ymax = cha.mean + cha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = cha_rhiz.let, y = cha.mean + cha.sd + 10), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 4) +
  ylab('Observed ASV Richness (S)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,500),expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = "")) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm'))

cha_rhiz.plot

# All Bulk Soil Samples #
bulk.rich <- filter(all_rich.mnsd, Compartment == "Bulk Soil")
bulk.richraw <- filter(all.rich, Compartment == 'Bulk Soil')
bulk.richraw$Soils <- gsub('Non-PSF Soil', 'Non PSF Soil', bulk.richraw$Soils)
shapiro.test(bulk.richraw$Shannon)

comm_bulk.rich <- bulk.rich[7,]
comm_bulk.rich$`Plant Species` <- gsub("S. helvola", "C. fasciculata", comm_bulk.rich$`Plant Species`)
bulk.rich <- rbind(bulk.rich, comm_bulk.rich)
comm_bulk.rich$`Plant Species` <- gsub("C. fasciculata", "D. illinoense", comm_bulk.rich$`Plant Species`)
bulk.rich <- rbind(bulk.rich, comm_bulk.rich)
comm_bulk.rich$`Plant Species` <- gsub("D. illinoense", "A. bracteata", comm_bulk.rich$`Plant Species`)
bulk.rich <- rbind(bulk.rich, comm_bulk.rich)
comm_bulk.rich$`Plant Species` <- gsub("A. bracteata", "T. repens", comm_bulk.rich$`Plant Species`)
bulk.rich <- rbind(bulk.rich, comm_bulk.rich)
comm_bulk.rich$`Plant Species` <- gsub("T. repens", "M. truncatula", comm_bulk.rich$`Plant Species`)
bulk.rich <- rbind(bulk.rich, comm_bulk.rich)
npsf.rich <- bulk.rich[2,]
npsf.rich$`Plant Species` <- gsub("C. fasciculata", "S. helvola", npsf.rich$`Plant Species`)
bulk.rich <- rbind(bulk.rich, npsf.rich)
npsf.rich <- bulk.rich[4,]
npsf.rich$`Plant Species` <- gsub("D. illinoense", "A. bracteata", npsf.rich$`Plant Species`)
bulk.rich <- rbind(bulk.rich, npsf.rich)
npsf.rich <- bulk.rich[9,]
npsf.rich$`Plant Species` <- gsub("T. repens", "M. truncatula", npsf.rich$`Plant Species`)
bulk.rich <- rbind(bulk.rich, npsf.rich)

## Strophostyles ##
fb_bulk.richraw <- filter(bulk.richraw, Plant == 'S. helvola')
fb_bulk.richraw <- rbind(fb_bulk.richraw, bulk.richraw[c('5076', '5077', '5078'),])

fb_bulk.rich <- filter(bulk.rich, `Plant Species` == 'S. helvola')
fb_bulk_sha.aov <- aov(Shannon ~ Soils, data = fb_bulk.richraw)
summary(fb_bulk_sha.aov)

fb_bulk_sha.hsd <- TukeyHSD(fb_bulk_sha.aov)
fb_bulk_sha.hsd

fb_bulk_sha.let <- multcompLetters4(fb_bulk_sha.aov, fb_bulk_sha.hsd)
fb_bulk_sha.let <- fb_bulk_sha.let$Soils$Letters[sort(names(fb_bulk_sha.let$Soils$Letters))]

## Chamaecrista ##
cc_bulk.richraw <- filter(bulk.richraw, Plants == 'C. fasciculata')
cc_bulk.richraw <- rbind(cc_bulk.richraw, bulk.richraw[c('5044', '5045', '5046'),])

cc_bulk_sha.aov <- aov(Shannon ~ Soils, data = cc_bulk.richraw)
summary(cc_bulk_sha.aov)

cc_bulk_sha.hsd <- TukeyHSD(cc_bulk_sha.aov)
cc_bulk_sha.hsd

cc_bulk_sha.let <- multcompLetters4(cc_bulk_sha.aov, cc_bulk_sha.hsd)
cc_bulk_sha.let <- cc_bulk_sha.let$Soils$Letters[sort(names(cc_bulk_sha.let$Soils$Letters))]

## Desmodium ##
ds_bulk.richraw <- filter(bulk.richraw, Plant == 'D. illinoense')
ds_bulk.richraw <- rbind(ds_bulk.richraw, bulk.richraw[c('5044', '5045', '5046'),])

ds_bulk_sha.aov <- aov(Shannon ~ Soils, data = ds_bulk.richraw)
summary(ds_bulk_sha.aov)

ds_bulk_sha.hsd <- TukeyHSD(ds_bulk_sha.aov)
ds_bulk_sha.hsd

ds_bulk_sha.let <- multcompLetters4(ds_bulk_sha.aov, ds_bulk_sha.hsd)
ds_bulk_sha.let <- ds_bulk_sha.let$Soils$Letters[sort(names(ds_bulk_sha.let$Soils$Letters))]

## Amphicarpaea ##
hp_bulk.richraw <- filter(bulk.richraw, Plant == 'A. bracteata')
hp_bulk.richraw <- rbind(hp_bulk.richraw, bulk.richraw[c('5044', '5045','5046', '5129-1', '5129-2', '5129-3'),])

hp_bulk_sha.aov <- aov(Shannon ~ Soils, data = hp_bulk.richraw)
summary(hp_bulk_sha.aov)

hp_bulk_sha.hsd <- TukeyHSD(hp_bulk_sha.aov)
hp_bulk_sha.hsd

hp_bulk_sha.let <- multcompLetters4(hp_bulk_sha.aov, hp_bulk_sha.hsd)
hp_bulk_sha.let <- hp_bulk_sha.let$Soils$Letters[sort(names(hp_bulk_sha.let$Soils$Letters))]

## Trifolium ##
cl_bulk.richraw <- filter(bulk.richraw, Plant == 'T. repens')
cl_bulk.richraw <- rbind(cl_bulk.richraw, bulk.richraw[c('5044', '5045','5046'),])

cl_bulk_sha.aov <- aov(Shannon ~ Soils, data = cl_bulk.richraw)
summary(cl_bulk_sha.aov)

cl_bulk_sha.hsd <- TukeyHSD(cl_bulk_sha.aov)
cl_bulk_sha.hsd

cl_bulk_sha.let <- multcompLetters4(cl_bulk_sha.aov, cl_bulk_sha.hsd)
cl_bulk_sha.let <- cl_bulk_sha.let$Soils$Letters[sort(names(cl_bulk_sha.let$Soils$Letters))]

## Medicago ##
md_bulk.richraw <- filter(bulk.richraw, Plant == 'M. truncatula')
md_bulk.richraw <- rbind(md_bulk.richraw, bulk.richraw[c('5044', '5045','5046', '5210', '5211', '5212'),])

md_bulk_sha.aov <- aov(Shannon ~ Soils, data = md_bulk.richraw)
summary(md_bulk_sha.aov)

md_bulk_sha.hsd <- TukeyHSD(md_bulk_sha.aov)
md_bulk_sha.hsd

md_bulk_sha.let <- multcompLetters4(md_bulk_sha.aov, md_bulk_sha.hsd)
md_bulk_sha.let <- md_bulk_sha.let$Soils$Letters[sort(names(md_bulk_sha.let$Soils$Letters))]

## Adding Letters ##
sha_bulk.let <- c(hp_bulk_sha.let,
                  ds_bulk_sha.let,
                  cc_bulk_sha.let,
                  fb_bulk_sha.let,
                  cl_bulk_sha.let,
                  md_bulk_sha.let)
bulk.rich <- arrange(bulk.rich, `Plant Species`, `Soil Type`)
bulk.rich <- cbind(bulk.rich, sha_bulk.let)

### Final Plot ###
sha_bulk.plot <- ggplot(bulk.rich, aes(x = `Plant Species`, y = sha.mean, fill = `Soil Type`, color = `Soil Type`)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = sha.mean - sha.sd, ymax = sha.mean + sha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = sha_bulk.let, y = sha.mean + sha.sd + 0.1), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 4) +
  ylab('Shannon Diversity (H)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,7.5),expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = "Bulk Soil")) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm'))

sha_bulk.plot

# Bulk Soil Shannon Evenness #
## Strophostyles ##
fb_bulk_evn.aov <- aov(ShaEvn ~ Soils, data = fb_bulk.richraw)
summary(fb_bulk_evn.aov)

fb_bulk_evn.hsd <- TukeyHSD(fb_bulk_evn.aov)
fb_bulk_evn.hsd

fb_bulk_evn.let <- multcompLetters4(fb_bulk_evn.aov, fb_bulk_evn.hsd)
fb_bulk_evn.let <- fb_bulk_evn.let$Soils$Letters[sort(names(fb_bulk_evn.let$Soils$Letters))]

## Chamecrista ##
cc_bulk_evn.aov <- aov(ShaEvn ~ Soils, data = cc_bulk.richraw)
summary(cc_bulk_evn.aov)

cc_bulk_evn.hsd <- TukeyHSD(cc_bulk_evn.aov)
cc_bulk_evn.hsd

cc_bulk_evn.let <- multcompLetters4(cc_bulk_evn.aov, cc_bulk_evn.hsd)
cc_bulk_evn.let <- cc_bulk_evn.let$Soils$Letters[sort(names(cc_bulk_evn.let$Soils$Letters))]

## Desmodium ##
ds_bulk_evn.aov <- aov(ShaEvn ~ Soils, data = ds_bulk.richraw)
summary(ds_bulk_evn.aov)

ds_bulk_evn.hsd <- TukeyHSD(ds_bulk_evn.aov)
ds_bulk_evn.hsd

ds_bulk_evn.let <- multcompLetters4(ds_bulk_evn.aov, ds_bulk_evn.hsd)
ds_bulk_evn.let <- ds_bulk_evn.let$Soils$Letters[sort(names(ds_bulk_evn.let$Soils$Letters))]

## Amphicarpaea ##
hp_bulk_evn.aov <- aov(ShaEvn ~ Soils, data = hp_bulk.richraw)
summary(hp_bulk_evn.aov)

hp_bulk_evn.hsd <- TukeyHSD(hp_bulk_evn.aov)
hp_bulk_evn.hsd

hp_bulk_evn.let <- multcompLetters4(hp_bulk_evn.aov, hp_bulk_evn.hsd)
hp_bulk_evn.let <- hp_bulk_evn.let$Soils$Letters[sort(names(hp_bulk_evn.let$Soils$Letters))]

## Trifolium ##
cl_bulk_evn.aov <- aov(ShaEvn ~ Soils, data = cl_bulk.richraw)
summary(cl_bulk_evn.aov)

cl_bulk_evn.hsd <- TukeyHSD(cl_bulk_evn.aov)
cl_bulk_evn.hsd

cl_bulk_evn.let <- multcompLetters4(cl_bulk_evn.aov, cl_bulk_evn.hsd)
cl_bulk_evn.let <- cl_bulk_evn.let$Soils$Letters[sort(names(cl_bulk_evn.let$Soils$Letters))]

## Medicago ##
md_bulk_evn.aov <- aov(ShaEvn ~ Soils, data = md_bulk.richraw)
summary(md_bulk_evn.aov)

md_bulk_evn.hsd <- TukeyHSD(md_bulk_evn.aov)
md_bulk_evn.hsd

md_bulk_evn.let <- multcompLetters4(md_bulk_evn.aov, md_bulk_evn.hsd)
md_bulk_evn.let <- md_bulk_evn.let$Soils$Letters[sort(names(md_bulk_evn.let$Soils$Letters))]

## Adding Letters ##
evn_bulk.let <- c(hp_bulk_evn.let,
                  ds_bulk_evn.let,
                  cc_bulk_evn.let,
                  fb_bulk_evn.let,
                  cl_bulk_evn.let,
                  md_bulk_evn.let)
bulk.rich <- cbind(bulk.rich, evn_bulk.let)


### Final Plot ###
evn_bulk.plot <- ggplot(bulk.rich, aes(x = `Plant Species`, y = evn.mean, fill = `Soil Type`, color = `Soil Type`)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = evn.mean - evn.sd, ymax = evn.mean + evn.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = evn_bulk.let, y = evn.mean + evn.sd + 0.02), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 4) +
  ylab('Shannon Evenness (H/ln(S))') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,1),expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ .,name = "")) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm'))

evn_bulk.plot
# Chao1 Observed ASVs #
## Strophostyles ##
fb_bulk_cha.aov <- aov(Chao1 ~ Soils, data = fb_bulk.richraw)
summary(fb_bulk_cha.aov)

fb_bulk_cha.hsd <- TukeyHSD(fb_bulk_cha.aov)
fb_bulk_cha.hsd

fb_bulk_cha.let <- multcompLetters4(fb_bulk_cha.aov, fb_bulk_cha.hsd)
fb_bulk_cha.let <- fb_bulk_cha.let$Soils$Letters[sort(names(fb_bulk_cha.let$Soils$Letters))]

## Chamecrista ##
cc_bulk_cha.aov <- aov(Chao1~ Soils, data = cc_bulk.richraw)
summary(cc_bulk_cha.aov)

cc_bulk_cha.hsd <- TukeyHSD(cc_bulk_cha.aov)
cc_bulk_cha.hsd

cc_bulk_cha.let <- multcompLetters4(cc_bulk_cha.aov, cc_bulk_cha.hsd)
cc_bulk_cha.let <- cc_bulk_cha.let$Soils$Letters[sort(names(cc_bulk_cha.let$Soils$Letters))]

## Desmodium ##
ds_bulk_cha.aov <- aov(Chao1 ~ Soils, data = ds_bulk.richraw)
summary(ds_bulk_cha.aov)

ds_bulk_cha.hsd <- TukeyHSD(ds_bulk_cha.aov)
ds_bulk_cha.hsd

ds_bulk_cha.let <- multcompLetters4(ds_bulk_cha.aov, ds_bulk_cha.hsd)
ds_bulk_cha.let <- ds_bulk_cha.let$Soils$Letters[sort(names(ds_bulk_cha.let$Soils$Letters))]

## Amphicarpaea ##
hp_bulk_cha.aov <- aov(Chao1 ~ Soils, data = hp_bulk.richraw)
summary(hp_bulk_cha.aov)

hp_bulk_cha.hsd <- TukeyHSD(hp_bulk_cha.aov)
hp_bulk_cha.hsd

hp_bulk_cha.let <- multcompLetters4(hp_bulk_cha.aov, hp_bulk_cha.hsd)
hp_bulk_cha.let <- hp_bulk_cha.let$Soils$Letters[sort(names(hp_bulk_cha.let$Soils$Letters))]

## Trifolium ##
cl_bulk_cha.aov <- aov(Chao1 ~ Soils, data = cl_bulk.richraw)
summary(cl_bulk_cha.aov)

cl_bulk_cha.hsd <- TukeyHSD(cl_bulk_cha.aov)
cl_bulk_cha.hsd

cl_bulk_cha.let <- multcompLetters4(cl_bulk_cha.aov, cl_bulk_cha.hsd)
cl_bulk_cha.let <- cl_bulk_cha.let$Soils$Letters[sort(names(cl_bulk_cha.let$Soils$Letters))]

## Medicago ##
md_bulk_cha.aov <- aov(Chao1 ~ Soils, data = md_bulk.richraw)
summary(md_bulk_cha.aov)

md_bulk_cha.hsd <- TukeyHSD(md_bulk_cha.aov)
md_bulk_cha.hsd

md_bulk_cha.let <- multcompLetters4(md_bulk_cha.aov, md_bulk_cha.hsd)
md_bulk_cha.let <- md_bulk_cha.let$Soils$Letters[sort(names(md_bulk_cha.let$Soils$Letters))]

## Adding Letters ##
cha_bulk.let <- c(hp_bulk_cha.let,
                  ds_bulk_cha.let,
                  cc_bulk_cha.let,
                  fb_bulk_cha.let,
                  cl_bulk_cha.let,
                  md_bulk_cha.let)
bulk.rich <- cbind(bulk.rich, cha_bulk.let)

### Final Plot ###
cha_bulk.plot <- ggplot(bulk.rich, aes(x = `Plant Species`, y = cha.mean, fill = `Soil Type`, color = `Soil Type`)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = cha.mean - cha.sd, ymax = cha.mean + cha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = cha_bulk.let, y = cha.mean + cha.sd + 10), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 4) +
  ylab('Observed ASV Richness (S)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,400),expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = "")) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm'))

cha_bulk.plot

# All Root Samples #
root.rich <- filter(all_rich.mnsd, Compartment == "Root Endosphere")
root.richraw <- filter(all.rich, Compartment == 'Root Endosphere')
root.richraw$Soils <- gsub('Non-PSF Soil', 'Non PSF Soil', root.richraw$Soils)
shapiro.test(root.richraw$Shannon)

# Root Shannon Diversity #

## Strophostyles ##
fb_root.richraw <- filter(root.richraw, Plant == 'S. helvola')
fb_root_sha.aov <- aov(Shannon ~ Soils, data = fb_root.richraw)
summary(fb_root_sha.aov)

fb_root_sha.hsd <- TukeyHSD(fb_root_sha.aov)
fb_root_sha.hsd

fb_root_sha.let <- multcompLetters4(fb_root_sha.aov, fb_root_sha.hsd)
fb_root_sha.let <- fb_root_sha.let$Soils$Letters[sort(names(fb_root_sha.let$Soils$Letters))]

## Chamecrista ##
cc_root.richraw <- filter(root.richraw, Plant == 'C. fasciculata')
cc_root_sha.aov <- aov(Shannon ~ Soils, data = cc_root.richraw)
summary(cc_root_sha.aov)

cc_root_sha.hsd <- TukeyHSD(cc_root_sha.aov)
cc_root_sha.hsd

cc_root_sha.let <- multcompLetters4(cc_root_sha.aov, cc_root_sha.hsd)
cc_root_sha.let <- cc_root_sha.let$Soils$Letters[sort(names(cc_root_sha.let$Soils$Letters))]

## Desmodium ##
ds_root.richraw <- filter(root.richraw, Plant == 'D. illinoense')
ds_root_sha.aov <- aov(Shannon ~ Soils, data = ds_root.richraw)
summary(ds_root_sha.aov)

ds_root_sha.hsd <- TukeyHSD(ds_root_sha.aov)
ds_root_sha.hsd

ds_root_sha.let <- multcompLetters4(ds_root_sha.aov, ds_root_sha.hsd)
ds_root_sha.let <- ds_root_sha.let$Soils$Letters[sort(names(ds_root_sha.let$Soils$Letters))]

## Amphicarpaea ##
hp_root.richraw <- filter(root.richraw, Plant == 'A. bracteata')
hp_root_sha.aov <- aov(Shannon ~ Soils, data = hp_root.richraw)
summary(hp_root_sha.aov)

hp_root_sha.hsd <- TukeyHSD(hp_root_sha.aov)
hp_root_sha.hsd

hp_root_sha.let <- multcompLetters4(hp_root_sha.aov, hp_root_sha.hsd)
hp_root_sha.let <- hp_root_sha.let$Soils$Letters[sort(names(hp_root_sha.let$Soils$Letters))]

## Trifolium ##
cl_root.richraw <- filter(root.richraw, Plant == 'T. repens')
cl_root_sha.aov <- aov(Shannon ~ Soils, data = cl_root.richraw)
summary(cl_root_sha.aov)

cl_root_sha.hsd <- TukeyHSD(cl_root_sha.aov)
cl_root_sha.hsd

cl_root_sha.let <- multcompLetters4(cl_root_sha.aov, cl_root_sha.hsd)
cl_root_sha.let <- cl_root_sha.let$Soils$Letters[sort(names(cl_root_sha.let$Soils$Letters))]

## Medicago ##
md_root.richraw <- filter(root.richraw, Plant == 'M. truncatula')
md_root_sha.aov <- aov(Shannon ~ Soils, data = md_root.richraw)
summary(md_root_sha.aov)

md_root_sha.hsd <- TukeyHSD(md_root_sha.aov)
md_root_sha.hsd

md_root_sha.let <- multcompLetters4(md_root_sha.aov, md_root_sha.hsd)
md_root_sha.let <- md_root_sha.let$Soils$Letters[sort(names(md_root_sha.let$Soils$Letters))]

## Adding Letters ##
sha_root.let <- c(hp_root_sha.let,
                  cc_root_sha.let,
                  ds_root_sha.let,
                  md_root_sha.let,
                  fb_root_sha.let,
                  cl_root_sha.let)
root.rich <- cbind(root.rich, sha_root.let)

### Final Plot ###
sha_root.plot <- ggplot(root.rich, aes(x = `Plant Species`, y = sha.mean, fill = `Soil Type`, color = `Soil Type`)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = sha.mean - sha.sd, ymax = sha.mean + sha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = sha_root.let, y = sha.mean + sha.sd + 0.1), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 4) +
  ylab('Shannon Diversity (H)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,4.5), expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = 'Root Endosphere')) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24),
        axis.text.x = element_text(size = 12, face = c('bold.italic'), angle = 45, hjust = 1, vjust = 1),
        axis.title.y.right = element_text(size = 12, angle = -90, hjust = 0.5),
        axis.title.y.left = element_text(size =12, angle = 90, hjust = 0.5),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm'))

sha_root.plot
# Root Shannon Evenness #
## Strophostyles ##
fb_root.richraw <- filter(root.richraw, Plant == 'S. helvola')
fb_root_evn.aov <- aov(ShaEvn ~ Soils, data = fb_root.richraw)
summary(fb_root_evn.aov)

fb_root_evn.hsd <- TukeyHSD(fb_root_evn.aov)
fb_root_evn.hsd

fb_root_evn.let <- multcompLetters4(fb_root_evn.aov, fb_root_evn.hsd)
fb_root_evn.let <- fb_root_evn.let$Soils$Letters[sort(names(fb_root_evn.let$Soils$Letters))]

## Chamecrista ##
cc_root.richraw <- filter(root.richraw, Plant == 'C. fasciculata')
cc_root_evn.aov <- aov(ShaEvn ~ Soils, data = cc_root.richraw)
summary(cc_root_evn.aov)

cc_root_evn.hsd <- TukeyHSD(cc_root_evn.aov)
cc_root_evn.hsd

cc_root_evn.let <- multcompLetters4(cc_root_evn.aov, cc_root_evn.hsd)
cc_root_evn.let <- cc_root_evn.let$Soils$Letters[sort(names(cc_root_evn.let$Soils$Letters))]

## Desmodium ##
ds_root.richraw <- filter(root.richraw, Plant == 'D. illinoense')
ds_root_evn.aov <- aov(ShaEvn ~ Soils, data = ds_root.richraw)
summary(ds_root_evn.aov)

ds_root_evn.hsd <- TukeyHSD(ds_root_evn.aov)
ds_root_evn.hsd

ds_root_evn.let <- multcompLetters4(ds_root_evn.aov, ds_root_evn.hsd)
ds_root_evn.let <- ds_root_evn.let$Soils$Letters[sort(names(ds_root_evn.let$Soils$Letters))]

## Amphicarpaea ##
hp_root.richraw <- filter(root.richraw, Plant == 'A. bracteata')
hp_root_evn.aov <- aov(ShaEvn ~ Soils, data = hp_root.richraw)
summary(hp_root_evn.aov)

hp_root_evn.hsd <- TukeyHSD(hp_root_evn.aov)
hp_root_evn.hsd

hp_root_evn.let <- multcompLetters4(hp_root_evn.aov, hp_root_evn.hsd)
hp_root_evn.let <- hp_root_evn.let$Soils$Letters[sort(names(hp_root_evn.let$Soils$Letters))]

## Trifolium ##
cl_root.richraw <- filter(root.richraw, Plant == 'T. repens')
cl_root_evn.aov <- aov(ShaEvn ~ Soils, data = cl_root.richraw)
summary(cl_root_evn.aov)

cl_root_evn.hsd <- TukeyHSD(cl_root_evn.aov)
cl_root_evn.hsd

cl_root_evn.let <- multcompLetters4(cl_root_evn.aov, cl_root_evn.hsd)
cl_root_evn.let <- cl_root_evn.let$Soils$Letters[sort(names(cl_root_evn.let$Soils$Letters))]

## Medicago ##
md_root.richraw <- filter(root.richraw, Plant == 'M. truncatula')
md_root_evn.aov <- aov(ShaEvn ~ Soils, data = md_root.richraw)
summary(md_root_evn.aov)

md_root_evn.hsd <- TukeyHSD(md_root_evn.aov)
md_root_evn.hsd

md_root_evn.let <- multcompLetters4(md_root_evn.aov, md_root_evn.hsd)
md_root_evn.let <- md_root_evn.let$Soils$Letters[sort(names(md_root_evn.let$Soils$Letters))]

## Adding Letters ##
evn_root.let <- c(hp_root_evn.let,
                  cc_root_evn.let,
                  ds_root_evn.let,
                  md_root_evn.let,
                  fb_root_evn.let,
                  cl_root_evn.let)
root.rich <- cbind(root.rich, evn_root.let)

### Final Plot ###
evn_root.plot <- ggplot(root.rich, aes(x = `Plant Species`, y = evn.mean, fill = `Soil Type`, color = `Soil Type`)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = evn.mean - evn.sd, ymax = evn.mean + evn.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = evn_root.let, y = evn.mean + evn.sd + 0.02), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 4) +
  ylab('Shannon Evenness (H/ln(S))') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,1),expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = "")) +
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 24, face = 'bold', family = 'Liberation Sans'),
        axis.text.x = element_text(size = 12, face = c('bold.italic'), angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 12),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm'))

evn_root.plot
# Chao1 Observed ASVs #
## Strophostyles ##
fb_root.richraw <- filter(root.richraw, Plant == 'S. helvola')
fb_root_cha.aov <- aov(Chao1 ~ Soils, data = fb_root.richraw)
summary(fb_root_cha.aov)

fb_root_cha.hsd <- TukeyHSD(fb_root_cha.aov)
fb_root_cha.hsd

fb_root_cha.let <- multcompLetters4(fb_root_cha.aov, fb_root_cha.hsd)
fb_root_cha.let <- fb_root_cha.let$Soils$Letters[sort(names(fb_root_cha.let$Soils$Letters))]

## Chamecrista ##
cc_root.richraw <- filter(root.richraw, Plant == 'C. fasciculata')
cc_root_cha.aov <- aov(Chao1~ Soils, data = cc_root.richraw)
summary(cc_root_cha.aov)

cc_root_cha.hsd <- TukeyHSD(cc_root_cha.aov)
cc_root_cha.hsd

cc_root_cha.let <- multcompLetters4(cc_root_cha.aov, cc_root_cha.hsd)
cc_root_cha.let <- cc_root_cha.let$Soils$Letters[sort(names(cc_root_cha.let$Soils$Letters))]

## Desmodium ##
ds_root.richraw <- filter(root.richraw, Plant == 'D. illinoense')
ds_root_cha.aov <- aov(Chao1 ~ Soils, data = ds_root.richraw)
summary(ds_root_cha.aov)

ds_root_cha.hsd <- TukeyHSD(ds_root_cha.aov)
ds_root_cha.hsd

ds_root_cha.let <- multcompLetters4(ds_root_cha.aov, ds_root_cha.hsd)
ds_root_cha.let <- ds_root_cha.let$Soils$Letters[sort(names(ds_root_cha.let$Soils$Letters))]

## Amphicarpaea ##
hp_root.richraw <- filter(root.richraw, Plant == 'A. bracteata')
hp_root_cha.aov <- aov(Chao1 ~ Soils, data = hp_root.richraw)
summary(hp_root_cha.aov)

hp_root_cha.hsd <- TukeyHSD(hp_root_cha.aov)
hp_root_cha.hsd

hp_root_cha.let <- multcompLetters4(hp_root_cha.aov, hp_root_cha.hsd)
hp_root_cha.let <- hp_root_cha.let$Soils$Letters[sort(names(hp_root_cha.let$Soils$Letters))]

## Trifolium ##
cl_root.richraw <- filter(root.richraw, Plant == 'T. repens')
cl_root_cha.aov <- aov(Chao1 ~ Soils, data = cl_root.richraw)
summary(cl_root_cha.aov)

cl_root_cha.hsd <- TukeyHSD(cl_root_cha.aov)
cl_root_cha.hsd

cl_root_cha.let <- multcompLetters4(cl_root_cha.aov, cl_root_cha.hsd)
cl_root_cha.let <- cl_root_cha.let$Soils$Letters[sort(names(cl_root_cha.let$Soils$Letters))]

## Medicago ##
md_root.richraw <- filter(root.richraw, Plant == 'M. truncatula')
md_root_cha.aov <- aov(Chao1 ~ Soils, data = md_root.richraw)
summary(md_root_cha.aov)

md_root_cha.hsd <- TukeyHSD(md_root_cha.aov)
md_root_cha.hsd

md_root_cha.let <- multcompLetters4(md_root_cha.aov, md_root_cha.hsd)
md_root_cha.let <- md_root_cha.let$Soils$Letters[sort(names(md_root_cha.let$Soils$Letters))]

## Adding Letters ##
cha_root.let <- c(hp_root_cha.let,
                  cc_root_cha.let,
                  ds_root_cha.let,
                  md_root_cha.let,
                  fb_root_cha.let,
                  cl_root_cha.let)
root.rich <- cbind(root.rich, cha_root.let)

### Final Plot ###
cha_root.plot <- ggplot(root.rich, aes(x = `Plant Species`, y = cha.mean, fill = `Soil Type`, color = `Soil Type`)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = cha.mean - cha.sd, ymax = cha.mean + cha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = cha_root.let, y = cha.mean + cha.sd + 2), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 4) +
  ylab('Observed ASV Richness (S)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,120), expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = "")) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24),
        axis.text.x = element_text(size = 12, face = c('bold.italic'), angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 12),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm'))


cha_root.plot

## Putting all figures together ##
(cha_bulk.plot | evn_bulk.plot | sha_bulk.plot) /
  (cha_rhiz.plot | evn_rhiz.plot | sha_rhiz.plot) /
  (cha_root.plot | evn_root.plot | sha_root.plot) +
  plot_layout(guides = 'auto')

