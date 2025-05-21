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
reference <- opt$reference
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
soil.rori <- allOrients(soil.fprimer)

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

save.image("./test.RData")
# track reads through the pipeline #
getN <- function(x) sum(getUniques(x))
soil_final.track <- cbind(soil_prefilt.track[,1], soil_prefilt.track[,2], soil_postfilt.track[,2], sapply(soil.fdada, getN), sapply(soil.rdada, getN), sapply(soil.remerged, getN), rowSums(soil_nochim.st))
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
unqs.mock <- soil_nochim.st["ZymoMockDNA_1870_S132",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
mock.ref <- getSequences(zymo)
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

#### Phyloseq Object Construction and Filtering for Soils ####
raw_soil.ps <- phyloseq(otu_table(soil_nochim.st, taxa_are_rows = TRUE),
                        sample_data(soil_raw.met),
                        tax_table(soil_rdp.taxa))

raw_soil.dna <- Biostrings::DNAStringSet(taxa_names(raw_soil.ps))
names(raw_soil.dna) <- taxa_names(raw_soil.ps)
raw_soil.ps <- merge_phyloseq(raw_soil.ps, raw_soil.dna)

raw_soil.ps <- subset_taxa(raw_soil.ps, Phylum != "Plantae")
raw_soil.ps <- subset_taxa(raw_soil.ps, Phylum != "Cyanobacteriota")

raw_soil.ps <- subset_samples(raw_soil.ps, Plant != 'Imma')
raw_soil.ps <- subset_samples(raw_soil.ps, Compartment != 'Leaf')
raw_soil.ps <- subset_samples(raw_soil.ps, Plant != 'Lupine')

raw_soil.ps <- subset_taxa(raw_soil.ps, taxa_sums(raw_soil.ps) > 1000)
colnames(tax_table(raw_soil.ps)) <- c('rdp_Kingdom', 'rdp_Phylum', 'rdp_Class', 'rdp_Order', 'rdp_Family', 'rdp_Genus')

decompose_ps <- function(ps, label){
  # function that decomposes a phyloseq object into separate data.frame and refseq objects (does not include tree) #
  assign(paste0(label, '.tax'), as.data.frame(tax_table(ps)))
  assign(paste0(label, '.otu'), as.data.frame(otu_table(ps)))
  assign(paste0(label, '.met'), as(sample_data(ps), 'data.frame'))
  assign(paste0(label, '.dna'), refseq(ps))
  assign(paste0(label, '.fra'), cbind(get(paste0(label, '.tax')), get(paste0(label, '.otu'))))
}

decompose_ps(raw_soil.ps, raw_soil)
#### Cross-Validation of Soil Reads Using BLAST ####
library(rBLAST)

# Create local blast database from the 16S rRNA database using rBLAST #
blast.tar <- blast_db_get("16S_ribosomal_RNA.tar.gz", baseURL = 'https://ftp.ncbi.nlm.nih.gov/blast/db/', check_update = TRUE)
untar(blast.tar, exdir = './reference/16S_database')
list.files('./reference/16S_database')
blast.db <- blast(db = './reference/16S_database/16S_ribosomal_RNA')

# Performs the blast for each read and returns the best hit # 
soil.hits <- matrix(nrow = nrow(raw_soil.tax), ncol = 12)
soil.hits <- as.data.frame(soil.hits) 
hold <- c()
for(i in 1:length(soil.dna)){
  hold <- predict(blast.db, raw_soil.dna[i])
  soil.hits[i,] <- hold[1,]
  raw_soil.tax$Best_Hit[i] <- hold[1, 2]
}

# Filter out reads that do not correspond to a NCBI entry #
library(dplyr); packageVersion('dplyr')
filt_soil.tax <- filter(raw_soil.tax, !is.na(raw_soil.tax$Best_Hit))

# Output the resulting NCBI entry names to a list #
if(!dir.exists("./blast_hits")){
  dir.create('./blast_hits')
}
write.table(filt_soil.tax$Best_Hit, './blast_hits/soil_blast_hits.txt')

# Call the python script to retrieve the taxonomies of the matched entries #
system('python3 rRNA_BLAST.py -i ./blast_hits/soil_blast_hits.txt -o ./blast_hits/soil_ncbi_hits.csv')

# Read in the output from the python script and make new taxonomy table "soil_ncbi_fin.tax" #
soil_ncbi.taxa <- read.csv2('./blast_hits/soil_ncbi_hits.csv', header = FALSE, fill = TRUE)
soil_ncbi.int <- strsplit(as.character(soil_ncbi.taxa$V1), ",")
soil_ncbi_fin.tax <- do.call(rbind, lapply(soil_ncbi.int, function(x) { length(x) <- max(sapply(soil_ncbi.int, length)); x }))
soil_ncbi_fin.tax <- as.data.frame(soil_ncbi_fin.tax, stringsAsFactors = FALSE)
rownames(ncbi.taxa) <- rownames(filt_soil.tax)
colnames(soil_ncbi_fin.tax) <- c('Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'hold')
for(i in 1:nrow(soil_ncbi_fin.tax)){
  if(!is.na(soil_ncbi_fin.tax$hold[i])){
    soil_ncbi_fin.tax$Genus[i] <- soil_ncbi_fin.tax$hold[i]
  }
}

# Filter otu table and refseq object such that all reads without a BLAST assignment are removed #
soil_ncbi_fin.tax <- soil_ncbi_fin.tax[,1:7]
filter_soil.tax <- cbind(filt_soil.tax, soil_ncbi_fin.tax)

decompose_ps(raw_soil.ps, filt_soil)

filt_soil.tax <- filter_soil.tax
filt_soil.otu <- filter(filt_soil.otu, rownames(filt_soil.otu) %in% rownames(filt_soil.tax))
soil_dna.df <- as.data.frame(filt_soil.dna)
soil_dna.df <- filter(soil_dna.df, rownames(soil_dna.df) %in% rownames(filt_soil.tax))
filt_soil.dna <- DNAStringSet(soil_dna.df$x)
names(filt_soil.dna) <- rownames(filt_soil.tax)
filt_soil.tax <- as.matrix(filt_soil.tax)

# Make phyloseq object with filtered tables #
soil.ps <- phyloseq(otu_table(filt_soil.otu, taxa_are_rows = TRUE),
                        sample_data(filt_soil.met),
                        tax_table(filt_soil.tax),
                        refseq(filt_soil.dna))

soil.ps <- subset_taxa(soil.ps, taxa_sums(soil.ps) > 1000)

# Change the taxa names to represent comparatiove abundance and lowest identification level #
taxa_names(soil.ps) <- paste0('ASV', seq(ntaxa(soil.ps)))
soil.tax <- as.data.frame(tax_table(soil.ps))
for(i in 1:nrow(soil.tax)){
  if(!is.na(soil.tax$Genus[i])){
    taxa_names(soil.ps)[i] = paste0(taxa_names(soil.ps)[i], '(', soil.tax$Genus[i], ')')
  }else if(!is.na(soil.tax$Family[i])){
    taxa_names(soil.ps)[i] = paste0(taxa_names(soil.ps)[i], '(', soil.tax$Family[i], ')')
  }else if(!is.na(soil.tax$Order[i])){
    taxa_names(soil.ps)[i] = paste0(taxa_names(soil.ps)[i], '(', soil.tax$Order[i], ')')
  }else{
    taxa_names(soil.ps)[i] = paste0(taxa_names(soil.ps)[i], '(NA)')
  }
}

# produce final decomposed phyloseq object
decompose_ps(soil.ps, soil)
save.image("./test.RData")

#### Phylogenetic Tree Construction for Soils ####
# Output the reads into a fasta file #
writeXStringSet(as.character(soil.dna, "./reads/soil_input.fasta"))

# Perform a multiple sequence alignment using MAFFT #
system('mafft --auto --thread -1 ./reads/soil_input.fasta > ./reads/soils_aligned.fasta')

# Construct a tree using IQTree with a general time reversible model with a gamma distribution and invariant site copies #
system('iqtree -s ./reads/soils_aligned.fasta -m GTR+G+I -nt AUTO')

# Read the tree using ape and check that the tip labels match #
library(ape); packageVersion('ape')
soil.tre <- read.tree("./reads/soils_aligned.fasta.treefile")

# Combine final phyloseq object for soil samples #
soil.ps <- phyloseq(otu_table(soil.otu, taxa_are_rows = TRUE),
                    sample_data(soil.met),
                    tax_table(soil.tax),
                    refseq(soil.dna),
                    phy_tree(soil.tre))

#### Root Primer Removal ####
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

root.fprimer <- "GTGCCAGCMGCCGCGGTAA"
root.rprimer <- "GGACTACHVGGGTWTCTAAT"

root.fori <- allOrients(root.fprimer)
root.rori <- allOrients(root.fprimer)

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
save.image("./test.RData")

# track reads through the pipeline #
getN <- function(x) sum(getUniques(x))
root_final.track <- cbind(root_prefilt.track[,1], root_prefilt.track[,2], root_postfilt.track[,2], sapply(root.fdada, getN), sapply(root.rdada, getN), sapply(root.remerged, getN), rowSums(root_nochim.st))
colnames(root_final.track) <- c("pre-cutadapt", "post-cutadapt", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(root_final.track) <- root.names
root_final.track <- as.data.frame(root_final.track)

# Assign Taxonomy #
root_rdp.taxa <- assignTaxonomy(rownames(root_nochim.st), refFasta = reference, multithread = TRUE, verbose = TRUE)
root_rdp.taxa <- as.matrix(root_rdp.taxa)

# Load the metadata #
root_raw.met <- read.csv2(root_metadata, sep = ',', row.names = TRUE)
rownames(root_raw.met) <- root_raw.met$Sample
root_raw.met <- root_raw.met[,c('Sample', 'Plant', 'Soil_Treatment', 'Compartment')]

#### Phyloseq Object Construction and Filtering for roots ####
raw_root.ps <- phyloseq(otu_table(root_nochim.st, taxa_are_rows = TRUE),
                        sample_data(root_raw.met),
                        tax_table(root_rdp.taxa))

raw_root.dna <- Biostrings::DNAStringSet(taxa_names(raw_root.ps))
names(raw_root.dna) <- taxa_names(raw_root.ps)
raw_root.ps <- merge_phyloseq(raw_root.ps, raw_root.dna)

raw_root.ps <- subset_taxa(raw_root.ps, Phylum != "Plantae")
raw_root.ps <- subset_taxa(raw_root.ps, Phylum != "Cyanobacteriota")

raw_root.ps <- subset_samples(raw_root.ps, Plant != 'Imma')
raw_root.ps <- subset_samples(raw_root.ps, Compartment != 'Leaf')
raw_root.ps <- subset_samples(raw_root.ps, Plant != 'Lupine')

raw_root.ps <- subset_taxa(raw_root.ps, taxa_sums(raw_root.ps) > 1000)
colnames(tax_table(raw_root.ps)) <- c('rdp_Kingdom', 'rdp_Phylum', 'rdp_Class', 'rdp_Order', 'rdp_Family', 'rdp_Genus')

decompose_ps(raw_root.ps, raw_root)
#### Cross-Validation of root Reads Using BLAST ####

# Performs the blast for each read and returns the best hit # 
root.hits <- matrix(nrow = nrow(raw_root.tax), ncol = 12)
root.hits <- as.data.frame(root.hits) 
hold <- c()
for(i in 1:length(root.dna)){
  hold <- predict(blast.db, raw_root.dna[i])
  root.hits[i,] <- hold[1,]
  raw_root.tax$Best_Hit[i] <- hold[1, 2]
}

# Filter out reads that do not correspond to a NCBI entry #
filt_root.tax <- filter(raw_root.tax, !is.na(raw_root.tax$Best_Hit))

# Output the resulting NCBI entry names to a list #
if(!dir.exists("./blast_hits")){
  dir.create('./blast_hits')
}
write.table(filt_root.tax$Best_Hit, './blast_hits/root_blast_hits.txt')

# Call the python script to retrieve the taxonomies of the matched entries #
system('python3 rRNA_BLAST.py -i ./blast_hits/root_blast_hits.txt -o ./blast_hits/root_ncbi_hits.csv')

# Read in the output from the python script and make new taxonomy table "root_ncbi_fin.tax" #
root_ncbi.taxa <- read.csv2('./blast_hits/root_ncbi_hits.csv', header = FALSE, fill = TRUE)
root_ncbi.int <- strsplit(as.character(root_ncbi.taxa$V1), ",")
root_ncbi_fin.tax <- do.call(rbind, lapply(root_ncbi.int, function(x) { length(x) <- max(sapply(root_ncbi.int, length)); x }))
root_ncbi_fin.tax <- as.data.frame(root_ncbi_fin.tax, stringsAsFactors = FALSE)
rownames(ncbi.taxa) <- rownames(filt_root.tax)
colnames(root_ncbi_fin.tax) <- c('Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'hold')
for(i in 1:nrow(root_ncbi_fin.tax)){
  if(!is.na(root_ncbi_fin.tax$hold[i])){
    root_ncbi_fin.tax$Genus[i] <- root_ncbi_fin.tax$hold[i]
  }
}

# Filter otu table and refseq object such that all reads without a BLAST assignment are removed #
root_ncbi_fin.tax <- root_ncbi_fin.tax[,1:7]
filter_root.tax <- cbind(filt_root.tax, root_ncbi_fin.tax)

decompose_ps(raw_root.ps, filt_root)

filt_root.tax <- filter_root.tax
filt_root.otu <- filter(filt_root.otu, rownames(filt_root.otu) %in% rownames(filt_root.tax))
root_dna.df <- as.data.frame(filt_root.dna)
root_dna.df <- filter(root_dna.df, rownames(root_dna.df) %in% rownames(filt_root.tax))
filt_root.dna <- DNAStringSet(root_dna.df$x)
names(filt_root.dna) <- rownames(filt_root.tax)
filt_root.tax <- as.matrix(filt_root.tax)

# Make phyloseq object with filtered tables #
root.ps <- phyloseq(otu_table(filt_root.otu, taxa_are_rows = TRUE),
                    sample_data(filt_root.met),
                    tax_table(filt_root.tax),
                    refseq(filt_root.dna))

root.ps <- subset_taxa(root.ps, taxa_sums(root.ps) > 1000)

# Change the taxa names to represent comparatiove abundance and lowest identification level #
taxa_names(root.ps) <- paste0('ASV', seq(ntaxa(root.ps)))
root.tax <- as.data.frame(tax_table(root.ps))
for(i in 1:nrow(root.tax)){
  if(!is.na(root.tax$Genus[i])){
    taxa_names(root.ps)[i] = paste0(taxa_names(root.ps)[i], '(', root.tax$Genus[i], ')')
  }else if(!is.na(root.tax$Family[i])){
    taxa_names(root.ps)[i] = paste0(taxa_names(root.ps)[i], '(', root.tax$Family[i], ')')
  }else if(!is.na(root.tax$Order[i])){
    taxa_names(root.ps)[i] = paste0(taxa_names(root.ps)[i], '(', root.tax$Order[i], ')')
  }else{
    taxa_names(root.ps)[i] = paste0(taxa_names(root.ps)[i], '(NA)')
  }
}

# produce final decomposed phyloseq object
decompose_ps(root.ps, root)
save.image("./test.RData")

#### Phylogenetic Tree Construction for roots ####
# Output the reads into a fasta file #
writeXStringSet(as.character(root.dna, "./reads/root_input.fasta"))

# Perform a multiple sequence alignment using MAFFT #
system('mafft --auto --thread -1 ./reads/root_input.fasta > ./reads/roots_aligned.fasta')

# Construct a tree using IQTree with a general time reversible model with a gamma distribution and invariant site copies #
system('iqtree -s ./reads/roots_aligned.fasta -m GTR+G+I -nt AUTO')

# Read the tree using ape and check that the tip labels match #
library(ape); packageVersion('ape')
root.tre <- read.tree("./reads/roots_aligned.fasta.treefile")

# Combine final phyloseq object for root samples #
root.ps <- phyloseq(otu_table(root.otu, taxa_are_rows = TRUE),
                    sample_data(root.met),
                    tax_table(root.tax),
                    refseq(root.dna),
                    phy_tree(root.tre))



save.image("./test.RData")
