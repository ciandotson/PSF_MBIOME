# Rscript ~/PSF_MBIOME/analysis.R --raw_soil ~/test_PSF/reads/soil_reads --raw_root ~/test_PSF/reads/endo_reads | cat > PSF_log.txt #

#### Argument Parsing ####
if(!requireNamespace('optparse', quietly = TRUE)) install.packages("optparse")
library(optparse); packageVersion("optparse")
option_list <- list(
  make_option("--raw_soil", type = "character", help = "filepath containing the raw, untrimmed reads for the reads generated from the v4 primers (bulk soil, rhizosphere, and some nodule samples)"),
  make_option("--raw_root", type = "character", help = "filepath containing the raw, untrimmed reads for the reads generated from the v5-v7 primers( root endosphere and some nodule samples)"))

opt <- parse_args(OptionParser(option_list=option_list))

soil.dir <- opt$raw_soil
root.dir <- opt$raw_root

#### Nodule Count and Biomass Data Visualization ####
# Read in the phenotypic data and clean it for analysis #
nodnbio.data <- read.csv2("./nodnbio.csv", sep = ',')
colnames(nodnbio.data) <- c('Plant_Sample', 'Soil_Treatment', 'Nodule_Count', 'Aboveground_Biomass')
nodnbio.data$Aboveground_Biomass <- as.numeric(nodnbio.data$Aboveground_Biomass)
for(i in 1:nrow(nodnbio.data)){
  nodnbio.data$Grouped[i] <- paste0(nodnbio.data$Plant_Sample[i], '; ', nodnbio.data$Soil_Treatment[i])
}

# Group all observations by Plant Species and Soil Treatment and find the group mean and standard error # 
if(!requireNamespace('dplyr')) installed.packages('dplyr')
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
hp_nod.data <- filter(nodnbio.data, Plant_Sample == 'A. bracteata')
cl_nod.data <- filter(nodnbio.data, Plant_Sample == 'T. repens')
md_nod.data <- filter(nodnbio.data, Plant_Sample == 'M. truncatula')

# Perform ANOVA on nodule counts for each plant group #
fb_nod.aov <- aov(Nodule_Count~Soil_Treatment, fb_nod.data)
cc_nod.aov <- aov(Nodule_Count~Soil_Treatment, cc_nod.data)
ds_nod.aov <- aov(Nodule_Count~Soil_Treatment, ds_nod.data)
hp_nod.aov <- aov(Nodule_Count~Soil_Treatment, hp_nod.data)
cl_nod.aov <- aov(Nodule_Count~Soil_Treatment, cl_nod.data)
md_nod.aov <- aov(Nodule_Count~Soil_Treatment, md_nod.data)

# Look at summaries of each ANOVA of nodule values #
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
fb_nod.let <- fb_nod.let$Soil_Treatment$Letters[sort(names(fb_nod.let$Soil_Treatment$Letters))]
cc_nod.let <- multcompLetters4(cc_nod.aov, cc_nod.hsd)
cc_nod.let <- cc_nod.let$Soil_Treatment$Letters[sort(names(cc_nod.let$Soil_Treatment$Letters))]
ds_nod.let <- multcompLetters4(ds_nod.aov, ds_nod.hsd)
ds_nod.let <- ds_nod.let$Soil_Treatment$Letters[sort(names(ds_nod.let$Soil_Treatment$Letters))]
hp_nod.let <- multcompLetters4(hp_nod.aov, hp_nod.hsd)
hp_nod.let <- hp_nod.let$Soil_Treatment$Letters[sort(names(hp_nod.let$Soil_Treatment$Letters))]
cl_nod.let <- multcompLetters4(cl_nod.aov, cl_nod.hsd)
cl_nod.let <- cl_nod.let$Soil_Treatment$Letters[sort(names(cl_nod.let$Soil_Treatment$Letters))]
md_nod.let <- multcompLetters4(md_nod.aov, md_nod.hsd)
md_nod.let <- md_nod.let$Soil_Treatment$Letters[sort(names(md_nod.let$Soil_Treatment$Letters))]

# Save a master list of the letters assigned to each comparison #
nod.let <- c(hp_nod.let,
             cc_nod.let,
             ds_nod.let,
             md_nod.let,
             fb_nod.let,
             cl_nod.let)

# Create a data.frame that has both mean/ses and nodule comparisons #
nod.data <- cbind(nodnbio.mnsd, nod.let)

# Clean up the letters to match their tests and clean for plotting #
nod.data$nod.let[3] <- 'a'
nod.data$nod.let[1:2] <- 'b'
nod.data[1,3:6] <- 0
nod.data$Group <- factor(nod.data$Plant_Sample, levels = c('T. repens', 'M. truncatula', 'S. helvola', 'C. fasciculata', 'D. illinoense', 'A. bracteata'))

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
hp_bio.data <- filter(nodnbio.data, Plant_Sample == 'A. bracteata')
cl_bio.data <- filter(nodnbio.data, Plant_Sample == 'T. repens')
md_bio.data <- filter(nodnbio.data, Plant_Sample == 'M. truncatula')

fb_bio.aov <- aov(Aboveground_Biomass~Soil_Treatment, fb_bio.data)
cc_bio.aov <- aov(Aboveground_Biomass~Soil_Treatment, cc_bio.data)
ds_bio.aov <- aov(Aboveground_Biomass~Soil_Treatment, ds_bio.data)
hp_bio.aov <- aov(Aboveground_Biomass~Soil_Treatment, hp_bio.data)
cl_bio.aov <- aov(Aboveground_Biomass~Soil_Treatment, cl_bio.data)
md_bio.aov <- aov(Aboveground_Biomass~Soil_Treatment, md_bio.data)

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
fb_bio.let <- fb_bio.let$Soil_Treatment$Letters[sort(names(fb_bio.let$Soil_Treatment$Letters))]
cc_bio.let <- multcompLetters4(cc_bio.aov, cc_bio.hsd)
cc_bio.let <- cc_bio.let$Soil_Treatment$Letters[sort(names(cc_bio.let$Soil_Treatment$Letters))]
ds_bio.let <- multcompLetters4(ds_bio.aov, ds_bio.hsd)
ds_bio.let <- ds_bio.let$Soil_Treatment$Letters[sort(names(ds_bio.let$Soil_Treatment$Letters))]
hp_bio.let <- multcompLetters4(hp_bio.aov, hp_bio.hsd)
hp_bio.let <- hp_bio.let$Soil_Treatment$Letters[sort(names(hp_bio.let$Soil_Treatment$Letters))]
cl_bio.let <- multcompLetters4(cl_bio.aov, cl_bio.hsd)
cl_bio.let <- cl_bio.let$Soil_Treatment$Letters[sort(names(cl_bio.let$Soil_Treatment$Letters))]
md_bio.let <- multcompLetters4(md_bio.aov, md_bio.hsd)
md_bio.let <- md_bio.let$Soil_Treatment$Letters[sort(names(md_bio.let$Soil_Treatment$Letters))]

# Save a master list of the letters assigned to each comparison #
bio.let <- c(hp_bio.let,
             cc_bio.let,
             ds_bio.let,
             md_bio.let,
             fb_bio.let,
             cl_bio.let)


bio.data <- cbind(nodnbio.mnsd, bio.let)

bio.data[1,3:6] <- 0
bio.data[3,7] <- 'a'
bio.data[1:2,7] <- 'b'

bio.data$Group <- factor(bio.data$Plant_Sample, levels = c('T. repens', 'M. truncatula', 'S. helvola', 'C. fasciculata', 'D. illinoense', 'A. bracteata'))
bio.plot <- ggplot(bio.data, aes(x = Group, y = bio.mean, fill = Soil_Treatment, color = Soil_Treatment)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = bio.mean, ymax = bio.mean + bio.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = bio.let, y = bio.mean + bio.sd + 0.01), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 10) +
  ylab('Aboveground Biomass (grams)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(labels = c("Common Soil", "Non-PSF Soil", "PSF Soil"), values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(labels = c("Common Soil", "Non-PSF Soil", "PSF Soil"), values = c('black', 'black', 'black')) +
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
if(!requireNamespace("Biostrings")) install.packages("Biostrings")
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
if(!requireNamespace("dada2")) BiocManager::install("dada2")
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
soil_rdp.taxa <- assignTaxonomy(rownames(soil_nochim.st), refFasta = './reference/rdp_19_toGenus_trainset.fa.gz', multithread = TRUE, verbose = TRUE)
soil_rdp.taxa <- as.matrix(soil_rdp.taxa)

# Load the metadata #
soil_raw.met <- read.csv2('./metadata/soil_metadata.csv', sep = ',')
rownames(soil_raw.met) <- soil_raw.met$Sample
soil_raw.met <- soil_raw.met[,c('Sample', 'Plant', 'Soil_Treatment', 'Compartment')]

# Check to see if ASVs were denoise properly #
unqs.mock <- soil_nochim.st[,"ZymoMockDNA"]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
mock.ref <- getSequences('./reference/ssrRNAs/')
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
save(raw_soil.ps, file = './psf_abirdged.RData'
#### Cross-Validation of Soil Reads Using BLAST ####
if(!requireNamespace('rBLAST')) BiocManager::install('rBLAST')
library(rBLAST);packageVersion('rBLAST')

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
if(!requireNamespace('ape')) install.packages('ape')
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
fb_soil_nod_raw.ps <- subset_samples(soil_nod.ps, Plant == "S. helvola")
fb_soil_nod_raw.ps <- subset_taxa(fb_soil_nod_raw.ps, taxa_sums(fb_soil_nod_raw.ps) > 0)
fb_soil_nod.ps <- aggregate_top_taxa2(fb_soil_nod_raw.ps, 8, "ASV")
fb_nod_soil.name <- names(sort(taxa_sums(fb_soil_nod.ps), decreasing = TRUE))
fb_soil_nod.colr <- soil_nod.colr[fb_nod_soil.name,]
fb_soil_nod.df <- psmelt(fb_soil_nod.ps)
fb_soil_nod.df$ASVs <- factor(fb_soil_nod.df$ASV, levels = fb_nod_soil.name)
fb_soil_nod.df$Soil <- factor(fb_soil_nod.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

fb_soil_nod.plot <- ggplot(fb_soil_nod.df, aes(x = Soil, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(name = "Soil ASV", values = fb_soil_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "S. helvola")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(color = "black", size = 18, family = "Liberation Sans", angle = -45, vjust = 0.6, hjust = 0.1),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "B.")
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
  scale_fill_manual(name = "Soil ASV", values = cc_soil_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "C. fasciculata")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(size = 18, family = "Liberation Sans", angle = -45, hjust= 0.1, vjust = 0.6),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
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
  scale_fill_manual(name = "Soil ASV", values = ds_soil_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "D. illinoense")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(size = 18, family = "Liberation Sans", angle = -45, hjust= 0.1, vjust = 0.6),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "B.")
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
  scale_fill_manual(name = "Soil ASV", values = hp_soil_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "A. bracteata")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(size = 18, family = "Liberation Sans", angle = -45, hjust= 0.1, vjust = 0.6),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "B.")
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
  scale_fill_manual(name = "Soil ASV", values = cl_soil_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "T. repens")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(size = 18, family = "Liberation Sans", angle = -45, hjust= 0.1, vjust = 0.6),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "B.")
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
  scale_fill_manual(name = "Soil ASV", values = md_soil_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "M. truncatula")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(size = 18, family = "Liberation Sans", angle = -45, hjust= 0.1, vjust = 0.6),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "B.")
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
root_rdp.taxa <- assignTaxonomy(rownames(root_nochim.st), refFasta = './reference/rdp_19_toGenus_trainset.fa.gz', multithread = TRUE, verbose = TRUE)
root_rdp.taxa <- as.matrix(root_rdp.taxa)

# Load the metadata #
root_raw.met <- read.csv2('./metadata/endo_metadata.csv', sep = ',')
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
fb_root_nod_raw.ps <- subset_samples(root_nod.ps, Plant.Species == "S. helvola")
fb_root_nod_raw.ps <- subset_taxa(fb_root_nod_raw.ps, taxa_sums(fb_root_nod_raw.ps) > 0)
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
  scale_fill_manual(name = "Root ASV", values = fb_root_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "S. helvola")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(size = 18, family = "Liberation Sans", angle = -45, hjust= 0.1, vjust = 0.6),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "C.")
fb_root_nod.plot

fb_nod.plot <- (fb_soil_nod.plot | fb_root_nod.plot)

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
  scale_fill_manual(name = "Root ASV", values = cc_root_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "C. fasciculata")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(size = 18, family = "Liberation Sans", angle = -45, hjust= 0.1, vjust = 0.6),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right')+
  labs(tag = "C.")
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
  scale_fill_manual(name = "Root ASV", values = ds_root_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "D. illinoense")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(size = 18, family = "Liberation Sans", angle = -45, hjust= 0.1, vjust = 0.6),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
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
  scale_fill_manual(name = "Root ASV",values = hp_root_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "A. bracteata")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(size = 18, family = "Liberation Sans", angle = -45, hjust= 0.1, vjust = 0.6),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "C.")
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
  scale_fill_manual(name = "Root ASV", values = cl_root_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "T. repens")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(size = 18, family = "Liberation Sans", angle = -45, hjust= 0.1, vjust = 0.6),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "C.")
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
  scale_fill_manual(name = "Root ASV", values = md_root_nod.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = "M. truncatula")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18, family = "Liberation Sans"),
        axis.text.x.bottom = element_text(size = 18, family = "Liberation Sans", angle = -45, hjust= 0.1, vjust = 0.6),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size =18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.title = element_text(size = 18, face = "bold", family = "Liberation Sans"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'right') +
  labs(tag = "C.")
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

all.rich$Soils <- factor(all.rich$Soil_Treatment, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))
# Perform three way anova tests for plant species, compartment, and soil origin #
library(car); packageVersion("car")
all_cha.lm <- lm(Chao1~Soils*Plants*Comps, all.rich)
Anova(all_cha.lm)
all_evn.lm <- lm(ShaEvn~Soils*Plants*Comps, all.rich)
Anova(all_evn.lm)
all_sha.lm <- lm(Shannon~Soils*Plants*Comps, all.rich)
Anova(all_sha.lm)

# Do the same except for just the bulk soil and rhizosphere data #
soil.rich <- filter(all.rich, Compartment != "Root Endosphere") 
soil_cha.lm <- lm(Chao1~Soils*Plants*Comps, soil.rich)
Anova(soil_cha.lm)
soil_evn.lm <- lm(ShaEvn~Soils*Plants*Comps, soil.rich)
Anova(soil_evn.lm)
soil_sha.lm <- lm(Shannon~Soils*Plants*Comps, soil.rich)
Anova(soil_sha.lm)

# Do the same with only the native plant data #
naty.rich <- filter(soil.rich, Plant != "M. truncatula" & Plant != "T. repens")
naty.bin <- filter(naty.rich, Compartment == "Bulk Soil" & Soil_Treatment == "Non-PSF Soil")
naty.bin$Plant <- gsub("C. fasciculata", "S. helvola", naty.bin$Plant)
naty.bin$Plant <- gsub("D. illinoense", "A. bracteata", naty.bin$Plant)
naty_comm.bin <- filter(soil.rich, Tri == "SCBu")
naty_comm.bin$Plant <- gsub("S. helvola", "C. fasciculata", naty_comm.bin$Plant)
naty.rich <- rbind(naty.rich, naty.bin, naty_comm.bin)
naty_comm.bin$Plant <- gsub("C. fasciculata", "D. illinoense", naty_comm.bin$Plant)
naty.rich <- rbind(naty.rich, naty_comm.bin)
naty_comm.bin$Plant <- gsub("D. illinoense", "A. bracteata", naty_comm.bin$Plant)
naty.rich <- rbind(naty.rich, naty_comm.bin)

naty.rich$Soils <- gsub("Non-PSF Soil", "Non PSF Soil", naty.rich$Soils)
naty.rich$Soils <- factor(naty.rich$Soils, levels = c("Common Soil", "Non PSF Soil", "PSF Soil"))
naty.rich$Comps <- factor(naty.rich$Compartment, levels = c("Bulk Soil", "Rhizosphere"))
naty.rich$Plants <- factor(naty.rich$Plant, levels = c("S. helvola", "C. fasciculata", "D. illinoense", "A. bracteata"))
naty_cha.lm <- lm(Chao1~Soils*Plants*Comps, naty.rich)
Anova(naty_cha.lm)
naty_evn.lm <- lm(ShaEvn~Soils*Plants*Comps, naty.rich)
Anova(naty_evn.lm)
naty_sha.lm <- lm(Shannon~Soils*Plants*Comps, naty.rich)
Anova(naty_sha.lm)

# Do the same with only the non-native plant data #
nnat.rich <- filter(soil.rich, Plant == "M. truncatula" | Plant == "T. repens" | Tri == "SCBu")
nnat.rich$Soils <- gsub("Non-PSF Soil", "Non PSF Soil", nnat.rich$Soil_Treatment)
nnat.rich$Plants <- gsub("S. helvola", "T. repens", nnat.rich$Plant)
nnat.bin <- filter(nnat.rich, Tri == "TNBu" | Tri == "SCBu")
nnat.bin$Plants <- gsub("T. repens", "M. truncatula", nnat.bin$Plant)
nnat.rich <- rbind(nnat.rich, nnat.bin)
nnat.rich$Plants <- gsub("S. helvola", "M. truncatula", nnat.rich$Plants)

nnat.rich$Soils <- factor(nnat.rich$Soils, levels = c("Common Soil", "Non PSF Soil", "PSF Soil"))
nnat.rich$Comps <- factor(nnat.rich$Compartment, levels = c("Bulk Soil", "Rhizosphere"))
nnat.rich$Plants <- factor(nnat.rich$Plants, levels = c("T. repens", "M. truncatula"))
nnat_cha.lm <- lm(Chao1~Soils*Plants*Comps, nnat.rich)
Anova(nnat_cha.lm, type = "II")
nnat_evn.lm <- lm(ShaEvn~Soils*Plants*Comps, nnat.rich)
Anova(nnat_evn.lm)
nnat_sha.lm <- lm(Shannon~Soils*Plants*Comps, nnat.rich)
Anova(nnat_sha.lm)

# Now do the same for the roots #
root.rich <- filter(all.rich, Compartment == "Root Endosphere")
root_cha.lm <- lm(Chao1~Soils*Plants, root.rich)
Anova(root_cha.lm, type = "II")
root_evn.lm <- lm(ShaEvn~Soils*Plants, root.rich)
Anova(root_evn.lm)
root_sha.lm <- lm(Shannon~Soils*Plants, root.rich)
Anova(root_sha.lm)

# Now do the same for the native roots #
raty.rich <- filter(root.rich, Plant != "T. repens" & Plant != "M. truncatula")
raty_cha.lm <- lm(Chao1~Soils*Plants, raty.rich)
Anova(raty_cha.lm, type = "II")
raty_evn.lm <- lm(ShaEvn~Soils*Plants, raty.rich)
Anova(raty_evn.lm)
raty_sha.lm <- lm(Shannon~Soils*Plants, raty.rich)
Anova(raty_sha.lm)

# Now do the same for the non-native roots #
rnat.rich <- filter(root.rich, Plant == "T. repens" | Plant == "M. truncatula")
rnat_cha.lm <- lm(Chao1~Soils*Plants, rnat.rich)
Anova(rnat_cha.lm, type = "II")
rnat_evn.lm <- lm(ShaEvn~Soils*Plants, rnat.rich)
Anova(rnat_evn.lm)
rnat_sha.lm <- lm(Shannon~Soils*Plants, rnat.rich)
Anova(rnat_sha.lm)

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
  geom_text(aes(label = sha_rhiz.let, y = sha.mean + sha.sd + 0.1), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 8) +
  ylab('Shannon Diversity (H)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(labels = c("Common Soil", "Non-PSF Soil", "PSF Soil"), values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(labels = c("Common Soil", "Non-PSF Soil", "PSF Soil"), values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,6.3), breaks = seq(0,6.25, by = 1.25), expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ .,name = "Rhizosphere")) +
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 24, face = 'bold', family = "Liberation Sans"),
        axis.text.x = element_text(size = 24, face = 'bold.italic', family = "Liberation Sans"),
        axis.title.y = element_text(size = 12),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm')) +
    labs(tag = "C.")

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
  geom_text(aes(label = evn_rhiz.let, y = evn.mean + evn.sd + 0.02), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 8) +
  ylab('Shannon Evenness (H/ln(S))') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(labels = c("Common Soil", "Non-PSF Soil", "PSF Soil"), values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(labels = c("Common Soil", "Non-PSF Soil", "PSF Soil"), values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,1.1),breaks = seq(0,1,by = 0.25), expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = "Rhizosphere")) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm')) +
  labs(tag = "B.")
  

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
  geom_text(aes(label = cha_rhiz.let, y = cha.mean + cha.sd + 10), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 8) +
  ylab('Observed ASV Richness (S)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(labels = c("Common Soil", "Non-PSF Soil", "PSF Soil"), values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(labels = c("Common Soil", "Non-PSF Soil", "PSF Soil"), values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,500),expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = "Rhizosphere")) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm')) +
  labs(tag = "A.")

cha_rhiz.plot

# All Rhizosphere Sample Plots #
(cha_rhiz.plot) / 
  (evn_rhiz.plot) / 
  (sha_rhiz.plot) &
  theme(plot.tag = element_text(size = 22))

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
  geom_text(aes(label = sha_bulk.let, y = sha.mean + sha.sd + 0.1), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 8) +
  ylab('Shannon Diversity (H)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(labels = c("Common Soil", "Non-PSF Soil","PSF Soil"), values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(labels = c("Common Soil", "Non-PSF Soil","PSF Soil"), values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, by = 2), expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = "Bulk Soil")) +
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 24, face = 'bold', family = "Liberation Sans"),
        axis.text.x = element_text(size = 24, face = 'bold.italic', family = "Liberation Sans"),
        axis.title.y.left = element_text(size = 14, hjust = 0.5),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm')) +
  labs(tag = "C.")

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
  geom_text(aes(label = evn_bulk.let, y = evn.mean + evn.sd + 0.02), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 8) +
  ylab('Shannon Evenness (H/ln(S))') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by = 0.25), expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ .,name = "Bulk Soil")) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size = 14, hjust = 1),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm')) +
  labs(tag = "B.")

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
  geom_text(aes(label = cha_bulk.let, y = cha.mean + cha.sd + 10), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 8) +
  ylab('Observed ASV Richness (S)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,400), breaks = seq(0,400, by = 100), expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = "Bulk Soil")) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size = 14, hjust = 1),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm')) +
  labs(tag = "A.")

cha_bulk.plot

# All Bulk Soil Sample Plots #
(cha_bulk.plot) / 
(evn_bulk.plot) / 
(sha_bulk.plot) &
  theme(plot.tag = element_text(size = 22))
  

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

#### Beta Diversity Measurements and Visualizations ####
# All Soil Samples #
# Construct a Weighted Unifrac distance matrix and PCoA ordination #
soil_prop.ps <- transform_sample_counts(soil.ps, function(x) x/sum(x))

set.seed(248)
soil.wuni <- phyloseq::distance(soil_prop.ps, method = 'wunifrac')
soil.pcoa <- phyloseq::ordinate(soil_prop.ps, 'PCoA', distance = soil.wuni)

# Perform an NMDS analysis using the weighted Unifrac distance matrix, with the PCoA ordination as the starting ordination # 
soil.nmds <- metaMDS(soil.wuni, 
                     k = 5, try = 100, trymax = 1000, maxit = 999,
                     model = 'global', 
                     autotransform = FALSE, previous.best = soil.pcoa$vectors[,1:5])

# Save the loading scores for all axes and make a distance matrix from these scores #
soil_nmds.scores <- scores(soil.nmds, display = 'sites')
soil_nmds.dist <- dist(soil_nmds.scores)

# Fit a linear model using the vectorized weighted Unifrac ditsance matrix as a predictor of the vectorized loading score distance matrix to fine total R^2 of the model # 
soil_nmds.ffit <- lm(as.vector(soil_nmds.dist) ~ as.vector(soil.wuni))
summary(soil_nmds.ffit)
soil_nmds.totr2 <- summary(soil_nmds.ffit)$r.squared

# Axes Variance Calculation #
# Fit linear models as before expect to preidtc the distance matrix of each individual axis #
soil_nmds.dist1 <- dist(soil_nmds.scores[,1])
soil_nmds.fit1 <- lm(as.vector(soil_nmds.dist1)~as.vector(soil.wuni))
soil_nmds.r1 <- summary(soil_nmds.fit1)$r.squared

soil_nmds.dist2 <- dist(soil_nmds.scores[,2])
soil_nmds.fit2 <- lm(as.vector(soil_nmds.dist2)~as.vector(soil.wuni))
soil_nmds.r2 <- summary(soil_nmds.fit2)$r.squared

soil_nmds.dist3 <- dist(soil_nmds.scores[,3])
soil_nmds.fit3 <- lm(as.vector(soil_nmds.dist3)~as.vector(soil.wuni))
soil_nmds.r3 <- summary(soil_nmds.fit3)$r.squared

soil_nmds.dist4 <- dist(soil_nmds.scores[,4])
soil_nmds.fit4 <- lm(as.vector(soil_nmds.dist4)~as.vector(soil.wuni))
soil_nmds.r4 <- summary(soil_nmds.fit4)$r.squared

soil_nmds.dist5 <- dist(soil_nmds.scores[,5])
soil_nmds.fit5 <- lm(as.vector(soil_nmds.dist5)~as.vector(soil.wuni))
soil_nmds.r5 <- summary(soil_nmds.fit5)$r.squared

# Take the sum the R^2 value from each axis #
soil_nmds.comb <- soil_nmds.r1 + soil_nmds.r2 + soil_nmds.r3 + soil_nmds.r4 + soil_nmds.r5

# Divide each axis R^2 by the total of all axes and then multiply by the variation explained by the whole model
soil_nmds.axisr <- c()
for(i in 1:ncol(soil_nmds.scores)){
  soil_nmds.axisr[i] <- (get(paste0('soil_nmds.r', i)) / soil_nmds.comb) * soil_nmds.totr2 
}

### Calculating Variance Components ###
library(lme4); packageVersion('lme4')

# Construct a data.frame that has sample info and their loading scores #
decompose_ps(soil.ps, 'soil')
for(i in 1:nrow(soil$met)){
  soil$met$PSC[i] <- paste0(substr(soil$met$Plant[i],1,1), substr(soil$met$Soil_Treatment[i],1,1), substr(soil$met$Compartment[i],1,2))
}
soil$met$Plants <- factor(soil$met$Plant, levels = c('S. helvola', 'C. fasciculata', 'D. illinoense', 'A. bracteata', 'T. repens', 'M. truncatula', 'Common Soil'))
soil$met$Comps <- factor(soil$met$Compartment, levels = c('Bulk Soil', 'Rhizosphere'))
soil$met$Soils <- factor(soil$met$Soil_Treatment, levels = c('Common Soil', "Non-PSF Soil", "PSF Soil"))

sample_data(soil.ps) <- soil$met
soil_nmds.load <- cbind(soil$met, soil_nmds.scores)

# Test the mixed linear model on the first NMDS axis #
soil_nmds.vfit1 <- lmer(NMDS1 ~ (1|Soils) + (1|Plants) + (1|Comps) + (1|Soils:Plants) + (1|Soils:Comps) + (1|Plants:Comps) + (1|Soils:Plants:Comps), data = soil_nmds.load, REML = TRUE)
summary(soil_nmds.vfit1)
soil_nmds.vca1 <- as.data.frame(VarCorr(soil_nmds.vfit1))

# Using Loop to do each NMDS axis #
soil_nmds.vca <- matrix(nrow = 8, ncol = ncol(soil_nmds.scores))
hold <- c()
for(i in 1:ncol(soil_nmds.scores)){
  hold <- lmer(soil_nmds.scores[,i] ~ (1|Soils) + (1|Plants) + (1|Comps) + (1|Soils:Plants) + (1|Soils:Comps) + (1|Plants:Comps) + (1|Soils:Plants:Comps), data = soil_nmds.load, REML = TRUE)
  hold <- as.data.frame(VarCorr(hold))
  soil_nmds.vca[1,i] <- hold[1,4]
  soil_nmds.vca[2,i] <- hold[2,4]
  soil_nmds.vca[3,i] <- hold[3,4]
  soil_nmds.vca[4,i] <- hold[4,4]
  soil_nmds.vca[5,i] <- hold[5,4]
  soil_nmds.vca[6,i] <- hold[6,4]
  soil_nmds.vca[7,i] <- hold[7,4]
  soil_nmds.vca[8,i] <- hold[8,4]
}

# Save the variance components to their assigned variable/variable interaction and their NMDS loading axis #
rownames(soil_nmds.vca) <- c('Soil x Plant X Comp', 'Soil x Plant', 'Plant x Comp', 'Soil x Comp', 'Plant', 'Soil', 'Comp', 'Residual')
colnames(soil_nmds.vca) <- colnames(soil_nmds.scores)

# Calculate the total variance of each variance component#
soil_nmds.vtot <- colSums(soil_nmds.vca)

# Weight each variance component by the amount of variation each axis explains #
soil_nmds.wvca <- matrix(nrow = nrow(soil_nmds.vca), ncol = length(soil_nmds.axisr))
for(i in 1:length(soil_nmds.axisr)){
  for(j in 1:nrow(soil_nmds.vca)){
    soil_nmds.wvca[j,i] <- soil_nmds.vca[j,i]*soil_nmds.axisr[i] 
  }
}
# Take the total variance explained by each predictor and take the sum of those values #
rownames(soil_nmds.wvca) <- rownames(soil_nmds.vca); colnames(soil_nmds.wvca) <- colnames(soil_nmds.vca)
soil_nmds.tvca <- rowSums(soil_nmds.wvca)
soil_nmds.tvca <- as.data.frame(soil_nmds.tvca)
soil_nmds.ptot <- colSums(soil_nmds.tvca)

# Take the variance explained by each predictor and divide by the total variance explained and multiply by 100% #
soil_nmds.pvca <- matrix(nrow = nrow(soil_nmds.tvca), ncol = 1)
for(i in 1:nrow(soil_nmds.tvca)){
  soil_nmds.pvca[i,1] <- soil_nmds.tvca[i,1] / soil_nmds.ptot * 100
}

# Save the variation explained percentages of each predictor/predictor interaction #
rownames(soil_nmds.pvca) <- rownames(soil_nmds.vca); colnames(soil_nmds.pvca) <- 'Variance Explained'
soil_nmds.pvca

# Perform PermANOVA using all samples #
soil.adon <- adonis2(soil.wuni~Soils*Plants*Comps, soil$met, permutations = 999)
soil.adon_by <- adonis2(soil.wuni~Soils*Plants*Comps, soil$met, permutations = 999, by = 'terms')
soil.adon
soil.adon_by

soil.man <- manova(cbind(NMDS1, NMDS2) ~ Soils*Plants*Comps, soil_nmds.load)

# Perform PermDISP using all samples # 
soil.bdis <- betadisper(soil.wuni, group = soil$met$PS)
anova(soil.bdis)
TukeyHSD(soil.bdis)

## Native vs Non-Native ##
# Native #
# Construct a Weighted Unifrac distance matrix and PCoA ordination #
naty.ps <- subset_samples(soil.ps, Plant != "M. truncatula")
naty.ps <- subset_samples(naty.ps, Plant != "T. repens")
naty_prop.ps <- transform_sample_counts(naty.ps, function(x) x/sum(x))

set.seed(248)
naty.wuni <- phyloseq::distance(naty_prop.ps, method = 'wunifrac')
naty.pcoa <- phyloseq::ordinate(naty_prop.ps, 'PCoA', distance = naty.wuni)

# Perform an NMDS analysis using the weighted Unifrac distance matrix, with the PCoA ordination as the starting ordination # 
naty.nmds <- metaMDS(naty.wuni, 
                     k = 5, try = 100, trymax = 1000, maxit = 999,
                     model = 'global', 
                     autotransform = FALSE, previous.best = naty.pcoa$vectors[,1:5])

# Save the loading scores for all axes and make a distance matrix from these scores #
naty_nmds.scores <- scores(naty.nmds, display = 'sites')
naty_nmds.dist <- dist(naty_nmds.scores)

# Fit a linear model using the vectorized weighted Unifrac ditsance matrix as a predictor of the vectorized loading score distance matrix to fine total R^2 of the model # 
naty_nmds.ffit <- lm(as.vector(naty_nmds.dist) ~ as.vector(naty.wuni))
summary(naty_nmds.ffit)
naty_nmds.totr2 <- summary(naty_nmds.ffit)$r.squared

# Axes Variance Calculation #
# Fit linear models as before expect to preidtc the distance matrix of each individual axis #
naty_nmds.dist1 <- dist(naty_nmds.scores[,1])
naty_nmds.fit1 <- lm(as.vector(naty_nmds.dist1)~as.vector(naty.wuni))
naty_nmds.r1 <- summary(naty_nmds.fit1)$r.squared

naty_nmds.dist2 <- dist(naty_nmds.scores[,2])
naty_nmds.fit2 <- lm(as.vector(naty_nmds.dist2)~as.vector(naty.wuni))
naty_nmds.r2 <- summary(naty_nmds.fit2)$r.squared

naty_nmds.dist3 <- dist(naty_nmds.scores[,3])
naty_nmds.fit3 <- lm(as.vector(naty_nmds.dist3)~as.vector(naty.wuni))
naty_nmds.r3 <- summary(naty_nmds.fit3)$r.squared

naty_nmds.dist4 <- dist(naty_nmds.scores[,4])
naty_nmds.fit4 <- lm(as.vector(naty_nmds.dist4)~as.vector(naty.wuni))
naty_nmds.r4 <- summary(naty_nmds.fit4)$r.squared

naty_nmds.dist5 <- dist(naty_nmds.scores[,5])
naty_nmds.fit5 <- lm(as.vector(naty_nmds.dist5)~as.vector(naty.wuni))
naty_nmds.r5 <- summary(naty_nmds.fit5)$r.squared

# Take the sum the R^2 value from each axis #
naty_nmds.comb <- naty_nmds.r1 + naty_nmds.r2 + naty_nmds.r3 + naty_nmds.r4 + naty_nmds.r5

# Divide each axis R^2 by the total of all axes and then multiply by the variation explained by the whole model
naty_nmds.axisr <- c()
for(i in 1:ncol(naty_nmds.scores)){
  naty_nmds.axisr[i] <- (get(paste0('naty_nmds.r', i)) / naty_nmds.comb) * naty_nmds.totr2 
}

### Calculating Variance Components ###
library(lme4); packageVersion('lme4')

# Construct a data.frame that has sample info and their loading scores #
decompose_ps(naty.ps, 'naty')
for(i in 1:nrow(naty$met)){
  naty$met$PSC[i] <- paste0(substr(naty$met$Plant[i],1,1), substr(naty$met$Soil_Treatment[i],1,1), substr(naty$met$Compartment[i],1,2))
}
naty$met$Plants <- factor(naty$met$Plant, levels = c('S. helvola', 'C. fasciculata', 'D. illinoense', 'A. bracteata', 'T. repens', 'M. truncatula'))
naty$met$Comps <- factor(naty$met$Compartment, levels = c('Bulk Soil', 'Rhizosphere'))
naty$met$Soils <- factor(naty$met$Soil_Treatment, levels = c('Common Soil', "Non-PSF Soil", "PSF Soil"))

sample_data(naty.ps) <- naty$met
naty_nmds.load <- cbind(naty$met, naty_nmds.scores)

# Test the mixed linear model on the first NMDS axis #
naty_nmds.vfit1 <- lmer(NMDS1 ~ (1|Soils) + (1|Plants) + (1|Comps) + (1|Soils:Plants) + (1|Soils:Comps) + (1|Plants:Comps) + (1|Soils:Plants:Comps), data = naty_nmds.load, REML = TRUE)
summary(naty_nmds.vfit1)
naty_nmds.vca1 <- as.data.frame(VarCorr(naty_nmds.vfit1))

# Using Loop to do each NMDS axis #
naty_nmds.vca <- matrix(nrow = 8, ncol = ncol(naty_nmds.scores))
hold <- c()
for(i in 1:ncol(naty_nmds.scores)){
  hold <- lmer(naty_nmds.scores[,i] ~ (1|Soils) + (1|Plants) + (1|Comps) + (1|Soils:Plants) + (1|Soils:Comps) + (1|Plants:Comps) + (1|Soils:Plants:Comps), data = naty_nmds.load, REML = TRUE)
  hold <- as.data.frame(VarCorr(hold))
  naty_nmds.vca[1,i] <- hold[1,4]
  naty_nmds.vca[2,i] <- hold[2,4]
  naty_nmds.vca[3,i] <- hold[3,4]
  naty_nmds.vca[4,i] <- hold[4,4]
  naty_nmds.vca[5,i] <- hold[5,4]
  naty_nmds.vca[6,i] <- hold[6,4]
  naty_nmds.vca[7,i] <- hold[7,4]
  naty_nmds.vca[8,i] <- hold[8,4]
}

# Save the variance components to their assigned variable/variable interaction and their NMDS loading axis #
rownames(naty_nmds.vca) <- c('Soil x Plant X Comp', 'Soil x Plant', 'Plant x Comp', 'Soil x Comp', 'Plant', 'Soil', 'Comp', 'Residual')
colnames(naty_nmds.vca) <- colnames(naty_nmds.scores)

# Calculate the total variance of each variance component#
naty_nmds.vtot <- colSums(naty_nmds.vca)

# Weight each variance component by the amount of variation each axis explains #
naty_nmds.wvca <- matrix(nrow = nrow(naty_nmds.vca), ncol = length(naty_nmds.axisr))
for(i in 1:length(naty_nmds.axisr)){
  for(j in 1:nrow(naty_nmds.vca)){
    naty_nmds.wvca[j,i] <- naty_nmds.vca[j,i]*naty_nmds.axisr[i] 
  }
}
# Take the total variance explained by each predictor and take the sum of those values #
rownames(naty_nmds.wvca) <- rownames(naty_nmds.vca); colnames(naty_nmds.wvca) <- colnames(naty_nmds.vca)
naty_nmds.tvca <- rowSums(naty_nmds.wvca)
naty_nmds.tvca <- as.data.frame(naty_nmds.tvca)
naty_nmds.ptot <- colSums(naty_nmds.tvca)

# Take the variance explained by each predictor and divide by the total variance explained and multiply by 100% #
naty_nmds.pvca <- matrix(nrow = nrow(naty_nmds.tvca), ncol = 1)
for(i in 1:nrow(naty_nmds.tvca)){
  naty_nmds.pvca[i,1] <- naty_nmds.tvca[i,1] / naty_nmds.ptot * 100
}

# Save the variation explained percentages of each predictor/predictor interaction #
rownames(naty_nmds.pvca) <- rownames(naty_nmds.vca); colnames(naty_nmds.pvca) <- 'Variance Explained'
naty_nmds.pvca

# Perform PermANOVA using all samples #
naty.adon <- adonis2(naty.wuni~Soils*Plants*Comps, naty$met, permutations = 999)
naty.adon_by <- adonis2(naty.wuni~Soils*Plants*Comps, naty$met, permutations = 999, by = 'terms')
naty.adon
naty.adon_by

# Perform PermDISP using all samples # 
naty.bdis <- betadisper(naty.wuni, group = naty$met$PS)
anova(naty.bdis)
TukeyHSD(naty.bdis)

# Non-Native #
# Construct a Weighted Unifrac distance matrix and PCoA ordination #
nnat.ps <- subset_samples(soil.ps, Plant == "M. truncatula" | Plant == "T. repens" | PSC == "SCBu")

nnat_prop.ps <- transform_sample_counts(nnat.ps, function(x) x/sum(x))

set.seed(248)
nnat.wuni <- phyloseq::distance(nnat_prop.ps, method = 'wunifrac')
nnat.pcoa <- phyloseq::ordinate(nnat_prop.ps, 'PCoA', distance = nnat.wuni)

# Perform an NMDS analysis using the weighted Unifrac distance matrix, with the PCoA ordination as the starting ordination # 
nnat.nmds <- metaMDS(nnat.wuni, 
                     k = 7, try = 100, trymax = 1000, maxit = 999,
                     model = 'global', 
                     autotransform = FALSE, previous.best = nnat.pcoa$vectors[,1:7])

# Save the loading scores for all axes and make a distance matrix from these scores #
nnat_nmds.scores <- scores(nnat.nmds, display = 'sites')
nnat_nmds.dist <- dist(nnat_nmds.scores)

# Fit a linear model using the vectorized weighted Unifrac ditsance matrix as a predictor of the vectorized loading score distance matrix to fine total R^2 of the model # 
nnat_nmds.ffit <- lm(as.vector(nnat_nmds.dist) ~ as.vector(nnat.wuni))
summary(nnat_nmds.ffit)
nnat_nmds.totr2 <- summary(nnat_nmds.ffit)$r.squared

# Axes Variance Calculation #
# Fit linear models as before expect to preidtc the distance matrix of each individual axis #
nnat_nmds.dist1 <- dist(nnat_nmds.scores[,1])
nnat_nmds.fit1 <- lm(as.vector(nnat_nmds.dist1)~as.vector(nnat.wuni))
nnat_nmds.r1 <- summary(nnat_nmds.fit1)$r.squared

nnat_nmds.dist2 <- dist(nnat_nmds.scores[,2])
nnat_nmds.fit2 <- lm(as.vector(nnat_nmds.dist2)~as.vector(nnat.wuni))
nnat_nmds.r2 <- summary(nnat_nmds.fit2)$r.squared

nnat_nmds.dist3 <- dist(nnat_nmds.scores[,3])
nnat_nmds.fit3 <- lm(as.vector(nnat_nmds.dist3)~as.vector(nnat.wuni))
nnat_nmds.r3 <- summary(nnat_nmds.fit3)$r.squared

nnat_nmds.dist4 <- dist(nnat_nmds.scores[,4])
nnat_nmds.fit4 <- lm(as.vector(nnat_nmds.dist4)~as.vector(nnat.wuni))
nnat_nmds.r4 <- summary(nnat_nmds.fit4)$r.squared

nnat_nmds.dist5 <- dist(nnat_nmds.scores[,5])
nnat_nmds.fit5 <- lm(as.vector(nnat_nmds.dist5)~as.vector(nnat.wuni))
nnat_nmds.r5 <- summary(nnat_nmds.fit5)$r.squared

nnat_nmds.dist6 <- dist(nnat_nmds.scores[,6])
nnat_nmds.fit6 <- lm(as.vector(nnat_nmds.dist6)~as.vector(nnat.wuni))
nnat_nmds.r6 <- summary(nnat_nmds.fit6)$r.squared

nnat_nmds.dist7 <- dist(nnat_nmds.scores[,7])
nnat_nmds.fit7 <- lm(as.vector(nnat_nmds.dist7)~as.vector(nnat.wuni))
nnat_nmds.r7 <- summary(nnat_nmds.fit7)$r.squared

# Take the sum the R^2 value from each axis #
nnat_nmds.comb <- nnat_nmds.r1 + nnat_nmds.r2 + nnat_nmds.r3 + nnat_nmds.r4 + nnat_nmds.r5 + nnat_nmds.r6 + nnat_nmds.r7

# Divide each axis R^2 by the total of all axes and then multiply by the variation explained by the whole model
nnat_nmds.axisr <- c()
for(i in 1:ncol(nnat_nmds.scores)){
  nnat_nmds.axisr[i] <- (get(paste0('nnat_nmds.r', i)) / nnat_nmds.comb) * nnat_nmds.totr2 
}

### Calculating Variance Components ###
library(lme4); packageVersion('lme4')

# Construct a data.frame that has sample info and their loading scores #
decompose_ps(nnat.ps, 'nnat')
for(i in 1:nrow(nnat$met)){
  nnat$met$PSC[i] <- paste0(substr(nnat$met$Plant[i],1,1), substr(nnat$met$Soil_Treatment[i],1,1), substr(nnat$met$Compartment[i],1,2))
}
nnat$met$Plants <- factor(nnat$met$Plant, levels = c('S. helvola', 'C. fasciculata', 'D. illinoense', 'A. bracteata', 'T. repens', 'M. truncatula'))
nnat$met$Comps <- factor(nnat$met$Compartment, levels = c('Bulk Soil', 'Rhizosphere'))
nnat$met$Soils <- factor(nnat$met$Soil_Treatment, levels = c('Common Soil', "Non-PSF Soil", "PSF Soil"))

sample_data(nnat.ps) <- nnat$met
nnat_nmds.load <- cbind(nnat$met, nnat_nmds.scores)

# Test the mixed linear model on the first NMDS axis #
nnat_nmds.vfit1 <- lmer(NMDS1 ~ (1|Soils) + (1|Plants) + (1|Comps) + (1|Soils:Plants) + (1|Soils:Comps) + (1|Plants:Comps) + (1|Soils:Plants:Comps), data = nnat_nmds.load, REML = TRUE)
summary(nnat_nmds.vfit1)
nnat_nmds.vca1 <- as.data.frame(VarCorr(nnat_nmds.vfit1))

# Using Loop to do each NMDS axis #
nnat_nmds.vca <- matrix(nrow = 8, ncol = ncol(nnat_nmds.scores))
hold <- c()
for(i in 1:ncol(nnat_nmds.scores)){
  hold <- lmer(nnat_nmds.scores[,i] ~ (1|Soils) + (1|Plants) + (1|Comps) + (1|Soils:Plants) + (1|Soils:Comps) + (1|Plants:Comps) + (1|Soils:Plants:Comps), data = nnat_nmds.load, REML = TRUE)
  hold <- as.data.frame(VarCorr(hold))
  nnat_nmds.vca[1,i] <- hold[1,4]
  nnat_nmds.vca[2,i] <- hold[2,4]
  nnat_nmds.vca[3,i] <- hold[3,4]
  nnat_nmds.vca[4,i] <- hold[4,4]
  nnat_nmds.vca[5,i] <- hold[5,4]
  nnat_nmds.vca[6,i] <- hold[6,4]
  nnat_nmds.vca[7,i] <- hold[7,4]
  nnat_nmds.vca[8,i] <- hold[8,4]
}

# Save the variance components to their assigned variable/variable interaction and their NMDS loading axis #
rownames(nnat_nmds.vca) <- c('Soil x Plant X Comp', 'Soil x Plant', 'Plant x Comp', 'Soil x Comp', 'Plant', 'Soil', 'Comp', 'Residual')
colnames(nnat_nmds.vca) <- colnames(nnat_nmds.scores)

# Calculate the total variance of each variance component#
nnat_nmds.vtot <- colSums(nnat_nmds.vca)

# Weight each variance component by the amount of variation each axis explains #
nnat_nmds.wvca <- matrix(nrow = nrow(nnat_nmds.vca), ncol = length(nnat_nmds.axisr))
for(i in 1:length(nnat_nmds.axisr)){
  for(j in 1:nrow(nnat_nmds.vca)){
    nnat_nmds.wvca[j,i] <- nnat_nmds.vca[j,i]*nnat_nmds.axisr[i] 
  }
}
# Take the total variance explained by each predictor and take the sum of those values #
rownames(nnat_nmds.wvca) <- rownames(nnat_nmds.vca); colnames(nnat_nmds.wvca) <- colnames(nnat_nmds.vca)
nnat_nmds.tvca <- rowSums(nnat_nmds.wvca)
nnat_nmds.tvca <- as.data.frame(nnat_nmds.tvca)
nnat_nmds.ptot <- colSums(nnat_nmds.tvca)

# Take the variance explained by each predictor and divide by the total variance explained and multiply by 100% #
nnat_nmds.pvca <- matrix(nrow = nrow(nnat_nmds.tvca), ncol = 1)
for(i in 1:nrow(nnat_nmds.tvca)){
  nnat_nmds.pvca[i,1] <- nnat_nmds.tvca[i,1] / nnat_nmds.ptot * 100
}

# Save the variation explained percentages of each predictor/predictor interaction #
rownames(nnat_nmds.pvca) <- rownames(nnat_nmds.vca); colnames(nnat_nmds.pvca) <- 'Variance Explained'
nnat_nmds.pvca

# Perform PermANOVA using all samples #
nnat.adon <- adonis2(nnat.wuni~Soils*Plants*Comps, nnat$met, permutations = 999)
nnat.adon_by <- adonis2(nnat.wuni~Soils*Plants*Comps, nnat$met, permutations = 999, by = 'terms')
nnat.adon
nnat.adon_by

# Perform PermDISP using all samples # 
nnat.bdis <- betadisper(nnat.wuni, group = nnat$met$PS)
anova(nnat.bdis)
TukeyHSD(nnat.bdis)


# Bulk Soils #
# Construct a Weighted Unifrac distance matrix and PCoA ordination #
bulk.ps <- subset_samples(soil.ps, Compartment == "Bulk Soil")
bulk.ps <- subset_taxa(bulk.ps, taxa_sums(bulk.ps) > 0)

bulk_prop.ps <- transform_sample_counts(bulk.ps, function(x) x/sum(x))

set.seed(248)
bulk.wuni <- phyloseq::distance(bulk_prop.ps, method = 'wunifrac')
bulk.pcoa <- phyloseq::ordinate(bulk_prop.ps, 'PCoA', distance = bulk.wuni)

# Perform an NMDS analysis using the weighted Unifrac distance matrix, with the PCoA ordination as the starting ordination # 
bulk.nmds <- metaMDS(bulk.wuni, 
                    k = 4, try = 20, trymax = 1000, maxit = 999,
                    model = 'global', 
                    autotransform = FALSE, previous.best = bulk.pcoa$vectors[,1:5])

# Save the loading scores for all axes and make a distance matrix from these scores #
bulk_nmds.scores <- scores(bulk.nmds, display = 'sites')
bulk_nmds.dist <- dist(bulk_nmds.scores)

# Fit a linear model using the vectorized weighted Unifrac ditsance matrix as a predictor of the vectorized loading score distance matrix to fine total R^2 of the model # 
bulk_nmds.ffit <- lm(as.vector(bulk_nmds.dist) ~ as.vector(bulk.wuni))
summary(bulk_nmds.ffit)
bulk_nmds.totr2 <- summary(bulk_nmds.ffit)$r.squared

# Axes Variance Calculation #
# Fit linear models as before expect to preidtc the distance matrix of each individual axis #
bulk_nmds.dist1 <- dist(bulk_nmds.scores[,1])
bulk_nmds.fit1 <- lm(as.vector(bulk_nmds.dist1)~as.vector(bulk.wuni))
bulk_nmds.r1 <- summary(bulk_nmds.fit1)$r.squared

bulk_nmds.dist2 <- dist(bulk_nmds.scores[,2])
bulk_nmds.fit2 <- lm(as.vector(bulk_nmds.dist2)~as.vector(bulk.wuni))
bulk_nmds.r2 <- summary(bulk_nmds.fit2)$r.squared

bulk_nmds.dist3 <- dist(bulk_nmds.scores[,3])
bulk_nmds.fit3 <- lm(as.vector(bulk_nmds.dist3)~as.vector(bulk.wuni))
bulk_nmds.r3 <- summary(bulk_nmds.fit3)$r.squared

bulk_nmds.dist4 <- dist(bulk_nmds.scores[,4])
bulk_nmds.fit4 <- lm(as.vector(bulk_nmds.dist4)~as.vector(bulk.wuni))
bulk_nmds.r4 <- summary(bulk_nmds.fit4)$r.squared

# Take the sum the R^2 value from each axis #
bulk_nmds.comb <- bulk_nmds.r1 + bulk_nmds.r2 + bulk_nmds.r3 + bulk_nmds.r4

# Divide each axis R^2 by the total of all axes and then multiply by the variation explained by the whole model
bulk_nmds.axisr <- c()
for(i in 1:ncol(bulk_nmds.scores)){
  bulk_nmds.axisr[i] <- (get(paste0('bulk_nmds.r', i)) / bulk_nmds.comb) * bulk_nmds.totr2 
}

### Calculating Variance Components ###
library(lme4); packageVersion('lme4')

# Construct a data.frame that has sample info and their loading scores #
decompose_ps(bulk.ps, 'bulk')
for(i in 1:nrow(bulk$met)){
  bulk$met$PS[i] <- paste0(substr(bulk$met$Plant[i],1,1), substr(bulk$met$Soil_Treatment[i],1,1))
}
for(i in 1:nrow(bulk$met)){
  if(bulk$met$PS[i] == "SC"){
    bulk$met$Plant[i] = "Common Soil"
  }
}
bulk$met$Plants <- factor(bulk$met$Plant, levels = c('S. helvola', 'C. fasciculata', 'D. illinoense', 'A. bracteata', 'T. repens', 'M. truncatula', 'Common Soil'))
bulk$met$Comps <- factor(bulk$met$Compartment, c('Bulk Soil'))
bulk$met$Soils <- factor(bulk$met$Soil_Treatment, levels = c('Common Soil', "Non-PSF Soil", "PSF Soil"))

sample_data(bulk.ps) <- bulk$met
bulk_nmds.load <- cbind(bulk$met, bulk_nmds.scores)

# Test the mixed linear model on the first NMDS axis #
bulk_nmds.vfit1 <- lmer(NMDS1 ~ (1|Soils) + (1|Plants) +(1|Soils:Plants), data = bulk_nmds.load, REML = TRUE)
summary(bulk_nmds.vfit1)
bulk_nmds.vca1 <- as.data.frame(VarCorr(bulk_nmds.vfit1))

# Using Loop to do each NMDS axis #
bulk_nmds.vca <- matrix(nrow = 4, ncol = ncol(bulk_nmds.scores))
hold <- c()
for(i in 1:ncol(bulk_nmds.scores)){
  hold <- lmer(bulk_nmds.scores[,i] ~  (1|Soils) + (1|Plants) +(1|Soils:Plants), data = bulk_nmds.load, REML = TRUE)
  hold <- as.data.frame(VarCorr(hold))
  bulk_nmds.vca[1,i] <- hold[1,4]
  bulk_nmds.vca[2,i] <- hold[2,4]
  bulk_nmds.vca[3,i] <- hold[3,4]
  bulk_nmds.vca[4,i] <- hold[4,4]
}

# Save the variance components to their assigned variable/variable interaction and their NMDS loading axis #
rownames(bulk_nmds.vca) <- c('PSF History x Plant', 'Plant', 'PSF History', 'Residual')
colnames(bulk_nmds.vca) <- colnames(bulk_nmds.scores)

# Calculate the total variance of each variance component#
bulk_nmds.vtot <- colSums(bulk_nmds.vca)

# Weight each variance component by the amount of variation each axis explains #
bulk_nmds.wvca <- matrix(nrow = nrow(bulk_nmds.vca), ncol = length(bulk_nmds.axisr))
for(i in 1:length(bulk_nmds.axisr)){
  for(j in 1:nrow(bulk_nmds.vca)){
    bulk_nmds.wvca[j,i] <- bulk_nmds.vca[j,i]*bulk_nmds.axisr[i] 
  }
}
# Take the total variance explained by each predictor and take the sum of those values #
rownames(bulk_nmds.wvca) <- rownames(bulk_nmds.vca); colnames(bulk_nmds.wvca) <- colnames(bulk_nmds.vca)
bulk_nmds.tvca <- rowSums(bulk_nmds.wvca)
bulk_nmds.tvca <- as.data.frame(bulk_nmds.tvca)
bulk_nmds.ptot <- colSums(bulk_nmds.tvca)

# Take the variance explained by each predictor and divide by the total variance explained and multiply by 100% #
bulk_nmds.pvca <- matrix(nrow = nrow(bulk_nmds.tvca), ncol = 1)
for(i in 1:nrow(bulk_nmds.tvca)){
  bulk_nmds.pvca[i,1] <- bulk_nmds.tvca[i,1] / bulk_nmds.ptot * 100
}

# Save the variation explained percentages of each predictor/predictor interaction #
rownames(bulk_nmds.pvca) <- rownames(bulk_nmds.vca); colnames(bulk_nmds.pvca) <- 'Variance Explained'
bulk_nmds.pvca

# Perform PermANOVA using all samples #
bulk.adon <- adonis2(bulk.wuni~Soils*Plants, bulk$met, permutations = 999)
bulk.adon_by <- adonis2(bulk.wuni~Soils*Plants, bulk$met, permutations = 999, by = 'terms')
bulk.adon
bulk.adon_by

# Perform PermDISP using all samples # 
bulk.bdis <- betadisper(bulk.wuni, group = bulk$met$PS)
anova(bulk.bdis)
TukeyHSD(bulk.bdis)

# Make an object with the data to be plotted #
bulk_nmds.load <- cbind(bulk$met, bulk_nmds.scores)
library(ggplot2);packageVersion('ggplot2')

# plot the NMDS ordination of all samples that will be used to make a patchwork plot #
bulk_nmds.plot <- ggplot(bulk_nmds.load, aes(NMDS1, NMDS2, color = Plants, shape = Soils)) +
  geom_point(size = 12) +
  scale_color_manual(name = "Bulk Soil Community Source", labels = c(expression(italic('S. helvola')), expression(italic('C. fasciculata')), expression(italic('D. illinoense')), expression(italic('A. bracteata')), expression(italic('T. repens')), expression(italic('M. truncatula')), "Common Soil"), values = c("#A6CEE3","#1F78B4","#FDBF6F", "#FF7F00","#CAB2D6", "#6A3D9A", "#B15928")) +
  scale_shape_manual(name = "PSF History", labels = c("Common Soil", "Non-PSF Soil", "PSF Soil"), values = c(15,17,16)) +
  xlab(expression("NMDS1 (R"^2 ~ "= 0.8167)")) +
  ylab(expression("NMDS2 (R"^2 ~ "= 0.1551)")) +
  theme_bw() +
  theme(legend.position = 'right',
        panel.grid = element_blank(),
        legend.text = element_text(size = 20, family = "Liberation Sans"),
        legend.title = element_text(size = 24, face = 'bold', family = "Liberation Sans"),
        axis.title = element_text(size= 24, face = 'bold', family = "Liberation Sans"),
        axis.text = element_text(color = "black", size = 8, family = "Liberation Sans")) +
  annotate('text', x = 0.06, y = -0.15,
           label = expression("(PERMANOVA) F"["9,29"] ~ "= 39.918, P < 0.001;  (PERMDISP) F"["9,29"] ~  "= 0.3344, P = 0.9526; 3D Stress = 0.0127"),
           size = 6, family = 'Liberation Sans') +
  coord_cartesian(ylim = c(-0.145,0.2))

bulk_nmds.plot

## Fuzzy Bean ##
# Create a phyloseq object that has all samples of the specific plant species and compartment #
fb_bulk.ps <- subset_samples(bulk.ps, Plant == "S. helvola" | Plant == "Common Soil" | Plant == "C. fasciculata" & Soil_Treatment == "Non-PSF Soil")
fb_bulk.ps <- subset_taxa(fb_bulk.ps, taxa_sums(fb_bulk.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(fb_bulk.ps, 'fb_bulk')
fb_bulk_prop.ps <- transform_sample_counts(fb_bulk.ps, function(x) x/sum(x))
set.seed(248)
fb_bulk.wuni <- phyloseq::distance(fb_bulk_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
fb_bulk.bdis <- betadisper(fb_bulk.wuni, group = fb_bulk$met$Soils)
anova(fb_bulk.bdis)
TukeyHSD(fb_bulk.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
fb_bulk.met <- filter(bulk_nmds.load, Plant == "S. helvola" | Plant == "Common Soil" | Plant == "C. fasciculata" & Soil_Treatment == "Non-PSF Soil")

# Perform MANOVA on all samples #
fb_bulk.man <- manova(cbind(NMDS1,NMDS2)~Soils, fb_bulk.met)
summary(fb_bulk.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
fb_bulk_pvsn.met <- filter(fb_bulk.met, Soil_Treatment != "Common Soil")
fb_bulk_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, fb_bulk_pvsn.met)
summary(fb_bulk_pvsn.man)

## PSF vs. Common ##
fb_bulk_pvsc.met <- filter(fb_bulk.met, Soil_Treatment != "Non-PSF Soil")
fb_bulk_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, fb_bulk_pvsc.met)
summary(fb_bulk_pvsc.man)

## Non-PSF vs. Common
fb_bulk_nvsc.met <- filter(fb_bulk.met, Soil_Treatment != "PSF Soil")
fb_bulk_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, fb_bulk_nvsc.met)
summary(fb_bulk_nvsc.man)

## Chamaecrista ##
# Create a phyloseq object that has all samples of the specific plant species and compartment #
cc_bulk.ps <- subset_samples(bulk.ps, Plant == "C. fasciculata" | Plant == "Common Soil")
cc_bulk.ps <- subset_taxa(cc_bulk.ps, taxa_sums(cc_bulk.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(cc_bulk.ps, 'cc_bulk')
cc_bulk_prop.ps <- transform_sample_counts(cc_bulk.ps, function(x) x/sum(x))
set.seed(248)
cc_bulk.wuni <- phyloseq::distance(cc_bulk_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
cc_bulk.bdis <- betadisper(cc_bulk.wuni, group = cc_bulk$met$Soils)
anova(cc_bulk.bdis)
TukeyHSD(cc_bulk.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
cc_bulk.met <- filter(bulk_nmds.load,  Plant == "C. fasciculata" | Plant == "Common Soil")

# Perform MANOVA on all samples #
cc_bulk.man <- manova(cbind(NMDS1,NMDS2)~Soils, cc_bulk.met)
summary(cc_bulk.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
cc_bulk_pvsn.met <- filter(cc_bulk.met, Soil_Treatment != "Common Soil")
cc_bulk_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, cc_bulk_pvsn.met)
summary(cc_bulk_pvsn.man)

## PSF vs. Common ##
cc_bulk_pvsc.met <- filter(cc_bulk.met, Soil_Treatment != "Non-PSF Soil")
cc_bulk_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, cc_bulk_pvsc.met)
summary(cc_bulk_pvsc.man)

## Non-PSF vs. Common
cc_bulk_nvsc.met <- filter(cc_bulk.met, Soil_Treatment != "PSF Soil")
cc_bulk_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, cc_bulk_nvsc.met)
summary(cc_bulk_nvsc.man)

## Desmodium ##
ds_bulk.ps <- subset_samples(bulk.ps, Plant == "D. illinoense" | Plant == "Common Soil")
ds_bulk.ps <- subset_taxa(ds_bulk.ps, taxa_sums(ds_bulk.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(ds_bulk.ps, 'ds_bulk')
ds_bulk_prop.ps <- transform_sample_counts(ds_bulk.ps, function(x) x/sum(x))
set.seed(248)
ds_bulk.wuni <- phyloseq::distance(ds_bulk_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
ds_bulk.bdis <- betadisper(ds_bulk.wuni, group = ds_bulk$met$Soils)
anova(ds_bulk.bdis)
TukeyHSD(ds_bulk.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
ds_bulk.met <- filter(bulk_nmds.load,  Plant == "D. illinoense" | Plant == "Common Soil")

# Perform MANOVA on all samples #
ds_bulk.man <- manova(cbind(NMDS1,NMDS2)~Soils, ds_bulk.met)
summary(ds_bulk.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
ds_bulk_pvsn.met <- filter(ds_bulk.met, Soil_Treatment != "Common Soil")
ds_bulk_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, ds_bulk_pvsn.met)
summary(ds_bulk_pvsn.man)

## PSF vs. Common ##
ds_bulk_pvsc.met <- filter(ds_bulk.met, Soil_Treatment != "Non-PSF Soil")
ds_bulk_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, ds_bulk_pvsc.met)
summary(ds_bulk_pvsc.man)

## Non-PSF vs. Common
ds_bulk_nvsc.met <- filter(ds_bulk.met, Soil_Treatment != "PSF Soil")
ds_bulk_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, ds_bulk_nvsc.met)
summary(ds_bulk_nvsc.man)

## Hog Peanut ##
hp_bulk.ps <- subset_samples(bulk.ps, Plant == "A. bracteata" | Plant == "Common Soil" | Plant == "D. illinoense" & Soil_Treatment == "Non-PSF Soil")
hp_bulk.ps <- subset_taxa(hp_bulk.ps, taxa_sums(hp_bulk.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(hp_bulk.ps, 'hp_bulk')
hp_bulk_prop.ps <- transform_sample_counts(hp_bulk.ps, function(x) x/sum(x))
set.seed(248)
hp_bulk.wuni <- phyloseq::distance(hp_bulk_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
hp_bulk.bdis <- betadisper(hp_bulk.wuni, group = hp_bulk$met$Soils)
anova(hp_bulk.bdis)
TukeyHSD(hp_bulk.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
hp_bulk.met <- filter(bulk_nmds.load, Plant == "A. bracteata" | Plant == "Common Soil" | Plant == "D. illinoense" & Soil_Treatment == "Non-PSF Soil")

# Perform MANOVA on all samples #
hp_bulk.man <- manova(cbind(NMDS1,NMDS2)~Soils, hp_bulk.met)
summary(hp_bulk.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
hp_bulk_pvsn.met <- filter(hp_bulk.met, Soil_Treatment != "Common Soil")
hp_bulk_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, hp_bulk_pvsn.met)
summary(hp_bulk_pvsn.man)

## PSF vs. Common ##
hp_bulk_pvsc.met <- filter(hp_bulk.met, Soil_Treatment != "Non-PSF Soil")
hp_bulk_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, hp_bulk_pvsc.met)
summary(hp_bulk_pvsc.man)

## Non-PSF vs. Common
hp_bulk_nvsc.met <- filter(hp_bulk.met, Soil_Treatment != "PSF Soil")
hp_bulk_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, hp_bulk_nvsc.met)
summary(hp_bulk_nvsc.man)

## Clover ##
cl_bulk.ps <- subset_samples(bulk.ps, Plant == "T. repens" | Plant == "Common Soil")
cl_bulk.ps <- subset_taxa(cl_bulk.ps, taxa_sums(cl_bulk.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(cl_bulk.ps, 'cl_bulk')
cl_bulk_prop.ps <- transform_sample_counts(cl_bulk.ps, function(x) x/sum(x))
set.seed(248)
cl_bulk.wuni <- phyloseq::distance(cl_bulk_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
cl_bulk.bdis <- betadisper(cl_bulk.wuni, group = cl_bulk$met$Soils)
anova(cl_bulk.bdis)
TukeyHSD(cl_bulk.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
cl_bulk.met <- filter(bulk_nmds.load, Plant == "T. repens" | Plant == "Common Soil")

# Perform MANOVA on all samples #
cl_bulk.man <- manova(cbind(NMDS1,NMDS2)~Soils, cl_bulk.met)
summary(cl_bulk.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
cl_bulk_pvsn.met <- filter(cl_bulk.met, Soil_Treatment != "Common Soil")
cl_bulk_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, cl_bulk_pvsn.met)
summary(cl_bulk_pvsn.man)

## PSF vs. Common ##
cl_bulk_pvsc.met <- filter(cl_bulk.met, Soil_Treatment != "Non-PSF Soil")
cl_bulk_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, cl_bulk_pvsc.met)
summary(cl_bulk_pvsc.man)

## Non-PSF vs. Common
cl_bulk_nvsc.met <- filter(cl_bulk.met, Soil_Treatment != "PSF Soil")
cl_bulk_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, cl_bulk_nvsc.met)
summary(cl_bulk_nvsc.man)

## Medicago ##
md_bulk.ps <- subset_samples(bulk.ps, Plant == "M. truncatula" | Plant == "Common Soil" | Plant == "T. repens" & Soil_Treatment == "Non-PSF Soil")
md_bulk.ps <- subset_taxa(md_bulk.ps, taxa_sums(md_bulk.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(md_bulk.ps, 'md_bulk')
md_bulk_prop.ps <- transform_sample_counts(md_bulk.ps, function(x) x/sum(x))
set.seed(248)
md_bulk.wuni <- phyloseq::distance(md_bulk_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
md_bulk.bdis <- betadisper(md_bulk.wuni, group = md_bulk$met$Soils)
anova(md_bulk.bdis)
TukeyHSD(md_bulk.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
md_bulk.met <- filter(bulk_nmds.load, Plant == "M. truncatula" | Plant == "Common Soil" | Plant == "T. repens" & Soil_Treatment == "Non-PSF Soil")

# Perform MANOVA on all samples #
md_bulk.man <- manova(cbind(NMDS1,NMDS2)~Soils, md_bulk.met)
summary(md_bulk.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
md_bulk_pvsn.met <- filter(md_bulk.met, Soil_Treatment != "Common Soil")
md_bulk_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, md_bulk_pvsn.met)
summary(md_bulk_pvsn.man)

## PSF vs. Common ##
md_bulk_pvsc.met <- filter(md_bulk.met, Soil_Treatment != "Non-PSF Soil")
md_bulk_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, md_bulk_pvsc.met)
summary(md_bulk_pvsc.man)

## Non-PSF vs. Common
md_bulk_nvsc.met <- filter(md_bulk.met, Soil_Treatment != "PSF Soil")
md_bulk_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, md_bulk_nvsc.met)
summary(md_bulk_nvsc.man)

# Rhizosphere #
# Construct a Weighted Unifrac distance matrix and PCoA ordination #
rhiz.ps <- subset_samples(soil.ps, Compartment == "Rhizosphere")
rhiz.ps <- subset_taxa(rhiz.ps, taxa_sums(rhiz.ps) > 0)

rhiz_prop.ps <- transform_sample_counts(rhiz.ps, function(x) x/sum(x))

set.seed(248)
rhiz.wuni <- phyloseq::distance(rhiz_prop.ps, method = 'wunifrac')
rhiz.pcoa <- phyloseq::ordinate(rhiz_prop.ps, 'PCoA', distance = rhiz.wuni)

# Perform an NMDS analysis using the weighted Unifrac distance matrix, with the PCoA ordination as the starting ordination # 
rhiz.nmds <- metaMDS(rhiz.wuni, 
                     k = 4, try = 20, trymax = 1000, maxit = 999,
                     model = 'global', 
                     autotransform = FALSE, previous.best = rhiz.pcoa$vectors[,1:5])

# Save the loading scores for all axes and make a distance matrix from these scores #
rhiz_nmds.scores <- scores(rhiz.nmds, display = 'sites')
rhiz_nmds.dist <- dist(rhiz_nmds.scores)

# Fit a linear model using the vectorized weighted Unifrac ditsance matrix as a predictor of the vectorized loading score distance matrix to fine total R^2 of the model # 
rhiz_nmds.ffit <- lm(as.vector(rhiz_nmds.dist) ~ as.vector(rhiz.wuni))
summary(rhiz_nmds.ffit)
rhiz_nmds.totr2 <- summary(rhiz_nmds.ffit)$r.squared

# Axes Variance Calculation #
# Fit linear models as before expect to preidtc the distance matrix of each individual axis #
rhiz_nmds.dist1 <- dist(rhiz_nmds.scores[,1])
rhiz_nmds.fit1 <- lm(as.vector(rhiz_nmds.dist1)~as.vector(rhiz.wuni))
rhiz_nmds.r1 <- summary(rhiz_nmds.fit1)$r.squared

rhiz_nmds.dist2 <- dist(rhiz_nmds.scores[,2])
rhiz_nmds.fit2 <- lm(as.vector(rhiz_nmds.dist2)~as.vector(rhiz.wuni))
rhiz_nmds.r2 <- summary(rhiz_nmds.fit2)$r.squared

rhiz_nmds.dist3 <- dist(rhiz_nmds.scores[,3])
rhiz_nmds.fit3 <- lm(as.vector(rhiz_nmds.dist3)~as.vector(rhiz.wuni))
rhiz_nmds.r3 <- summary(rhiz_nmds.fit3)$r.squared

rhiz_nmds.dist4 <- dist(rhiz_nmds.scores[,4])
rhiz_nmds.fit4 <- lm(as.vector(rhiz_nmds.dist4)~as.vector(rhiz.wuni))
rhiz_nmds.r4 <- summary(rhiz_nmds.fit4)$r.squared

# Take the sum the R^2 value from each axis #
rhiz_nmds.comb <- rhiz_nmds.r1 + rhiz_nmds.r2 + rhiz_nmds.r3 + rhiz_nmds.r4

# Divide each axis R^2 by the total of all axes and then multiply by the variation explained by the whole model
rhiz_nmds.axisr <- c()
for(i in 1:ncol(rhiz_nmds.scores)){
  rhiz_nmds.axisr[i] <- (get(paste0('rhiz_nmds.r', i)) / rhiz_nmds.comb) * rhiz_nmds.totr2 
}

### Calculating Variance Components ###

# Construct a data.frame that has sample info and their loading scores #
decompose_ps(rhiz.ps, 'rhiz')
for(i in 1:nrow(rhiz$met)){
  rhiz$met$PS[i] <- paste0(substr(rhiz$met$Plant[i],1,1), substr(rhiz$met$Soil_Treatment[i],1,1))
}
rhiz$met$Plants <- factor(rhiz$met$Plant, levels = c('S. helvola', 'C. fasciculata', 'D. illinoense', 'A. bracteata', 'T. repens', 'M. truncatula'))
rhiz$met$Comps <- factor(rhiz$met$Compartment, c('Rhizosphere'))
rhiz$met$Soils <- factor(rhiz$met$Soil_Treatment, levels = c('Common Soil', "Non-PSF Soil", "PSF Soil"))

sample_data(rhiz.ps) <- rhiz$met
rhiz_nmds.load <- cbind(rhiz$met, rhiz_nmds.scores)

# Test the mixed linear model on the first NMDS axis #
rhiz_nmds.vfit1 <- lmer(NMDS1 ~ (1|Soils) + (1|Plants) +(1|Soils:Plants), data = rhiz_nmds.load, REML = TRUE)
summary(rhiz_nmds.vfit1)
rhiz_nmds.vca1 <- as.data.frame(VarCorr(rhiz_nmds.vfit1))

# Using Loop to do each NMDS axis #
rhiz_nmds.vca <- matrix(nrow = 4, ncol = ncol(rhiz_nmds.scores))
hold <- c()
for(i in 1:ncol(rhiz_nmds.scores)){
  hold <- lmer(rhiz_nmds.scores[,i] ~  (1|Soils) + (1|Plants) +(1|Soils:Plants), data = rhiz_nmds.load, REML = TRUE)
  hold <- as.data.frame(VarCorr(hold))
  rhiz_nmds.vca[1,i] <- hold[1,4]
  rhiz_nmds.vca[2,i] <- hold[2,4]
  rhiz_nmds.vca[3,i] <- hold[3,4]
  rhiz_nmds.vca[4,i] <- hold[4,4]
}

# Save the variance components to their assigned variable/variable interaction and their NMDS loading axis #
rownames(rhiz_nmds.vca) <- c('PSF History x Plant', 'Plant', 'PSF History', 'Residual')
colnames(rhiz_nmds.vca) <- colnames(rhiz_nmds.scores)

# Calculate the total variance of each variance component#
rhiz_nmds.vtot <- colSums(rhiz_nmds.vca)

# Weight each variance component by the amount of variation each axis explains #
rhiz_nmds.wvca <- matrix(nrow = nrow(rhiz_nmds.vca), ncol = length(rhiz_nmds.axisr))
for(i in 1:length(rhiz_nmds.axisr)){
  for(j in 1:nrow(rhiz_nmds.vca)){
    rhiz_nmds.wvca[j,i] <- rhiz_nmds.vca[j,i]*rhiz_nmds.axisr[i] 
  }
}
# Take the total variance explained by each predictor and take the sum of those values #
rownames(rhiz_nmds.wvca) <- rownames(rhiz_nmds.vca); colnames(rhiz_nmds.wvca) <- colnames(rhiz_nmds.vca)
rhiz_nmds.tvca <- rowSums(rhiz_nmds.wvca)
rhiz_nmds.tvca <- as.data.frame(rhiz_nmds.tvca)
rhiz_nmds.ptot <- colSums(rhiz_nmds.tvca)

# Take the variance explained by each predictor and divide by the total variance explained and multiply by 100% #
rhiz_nmds.pvca <- matrix(nrow = nrow(rhiz_nmds.tvca), ncol = 1)
for(i in 1:nrow(rhiz_nmds.tvca)){
  rhiz_nmds.pvca[i,1] <- rhiz_nmds.tvca[i,1] / rhiz_nmds.ptot * 100
}

# Save the variation explained percentages of each predictor/predictor interaction #
rownames(rhiz_nmds.pvca) <- rownames(rhiz_nmds.vca); colnames(rhiz_nmds.pvca) <- 'Variance Explained'
rhiz_nmds.pvca

# Perform PermANOVA using all samples #
rhiz.adon <- adonis2(rhiz.wuni~Plants*Soils, rhiz$met, permutations = 999)
rhiz.adon_by <- adonis2(rhiz.wuni~Plants*Soils, rhiz$met, permutations = 999, by = 'terms')
rhiz.adon
rhiz.adon_by

# Perform PermDISP using all samples # 
rhiz.bdis <- betadisper(rhiz.wuni, group = rhiz$met$PS)
anova(rhiz.bdis)
TukeyHSD(rhiz.bdis)

# Make an object with the data to be plotted #
rhiz_nmds.load <- cbind(rhiz$met, rhiz_nmds.scores)

# plot the NMDS ordination of all samples that will be used to make a patchwork plot #
rhiz_nmds.plot <- ggplot(rhiz_nmds.load, aes(NMDS1, NMDS2, color = Plants, shape = Soils)) +
  geom_point(size = 8) +
  scale_color_manual(name = "Rhizosphere Community Source", labels = c(expression(italic('S. helvola')), expression(italic('C. fasciculata')), expression(italic('D. illinoense')), expression(italic('A. bracteata')), expression(italic('T. repens')), expression(italic('M. truncatula')), "Common Soil"), values = c("#A6CEE3","#1F78B4","#FDBF6F", "#FF7F00","#CAB2D6", "#6A3D9A", "#B15928")) +
  scale_shape_manual(name = "PSF History", labels = c("Common Soil", "Non-PSF Soil", "PSF Soil"), values = c(15,17,16)) +
  xlab(expression("NMDS1 (R"^2 ~ "= 0.5386)")) +
  ylab(expression("NMDS2 (R"^2 ~ "= 0.2525)")) +
  theme_bw() +
  theme(legend.position = 'right',
        panel.grid = element_blank(),
        legend.text = element_text(size = 20, family = "Liberation Sans"),
        legend.title = element_text(size = 24, face = 'bold', family = "Liberation Sans"),
        axis.title = element_text(size= 24, face = 'bold', family = "Liberation Sans"),
        axis.text = element_text(color = "black", size = 8, family = "Liberation Sans")) +
  annotate('text', x = -0.1, y = -0.29,
           label = expression("(PERMANOVA) F"["16,70"] ~ "= 25.443, P = 0.001;  (PERMDISP) F"["16,70"] ~  "= 4.9933, P < 0.001; 3D Stress = 0.0224"),
           size = 6, family = 'Liberation Sans') +
  coord_cartesian(ylim = c(-0.28, 0.15))

rhiz_nmds.plot

# Separate samples by plant species #
## Fuzzy Bean ##
fb_rhiz.ps <- subset_samples(rhiz.ps, Plant == "S. helvola")
fb_rhiz.ps <- subset_taxa(fb_rhiz.ps, taxa_sums(fb_rhiz.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(fb_rhiz.ps, 'fb_rhiz')
fb_rhiz_prop.ps <- transform_sample_counts(fb_rhiz.ps, function(x) x/sum(x))
set.seed(248)
fb_rhiz.wuni <- phyloseq::distance(fb_rhiz_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
fb_rhiz.bdis <- betadisper(fb_rhiz.wuni, group = fb_rhiz$met$Soils)
anova(fb_rhiz.bdis)
TukeyHSD(fb_rhiz.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
fb_rhiz.met <- filter(rhiz_nmds.load, Plant == "S. helvola")

# Perform MANOVA on all samples #
fb_rhiz.man <- manova(cbind(NMDS1,NMDS2)~Soils, fb_rhiz.met)
summary(fb_rhiz.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
fb_rhiz_pvsn.met <- filter(fb_rhiz.met, Soil_Treatment != "Common Soil")
fb_rhiz_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, fb_rhiz_pvsn.met)
summary(fb_rhiz_pvsn.man)

## PSF vs. Common ##
fb_rhiz_pvsc.met <- filter(fb_rhiz.met, Soil_Treatment != "Non-PSF Soil")
fb_rhiz_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, fb_rhiz_pvsc.met)
summary(fb_rhiz_pvsc.man)

## Non-PSF vs. Common ##
fb_rhiz_nvsc.met <- filter(fb_rhiz.met, Soil_Treatment != "PSF Soil")
fb_rhiz_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, fb_rhiz_nvsc.met)
summary(fb_rhiz_nvsc.man)

## Chamaecrista ##
cc_rhiz.ps <- subset_samples(rhiz.ps, Plant == "C. fasciculata")
cc_rhiz.ps <- subset_taxa(cc_rhiz.ps, taxa_sums(cc_rhiz.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(cc_rhiz.ps, 'cc_rhiz')
cc_rhiz_prop.ps <- transform_sample_counts(cc_rhiz.ps, function(x) x/sum(x))
set.seed(248)
cc_rhiz.wuni <- phyloseq::distance(cc_rhiz_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
cc_rhiz.bdis <- betadisper(cc_rhiz.wuni, group = cc_rhiz$met$Soils)
anova(cc_rhiz.bdis)
TukeyHSD(cc_rhiz.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
cc_rhiz.met <- filter(rhiz_nmds.load, Plant == "C. fasciculata")

# Perform MANOVA on all samples #
cc_rhiz.man <- manova(cbind(NMDS1,NMDS2)~Soils, cc_rhiz.met)
summary(cc_rhiz.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
cc_rhiz_pvsn.met <- filter(cc_rhiz.met, Soil_Treatment != "Common Soil")
cc_rhiz_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, cc_rhiz_pvsn.met)
summary(cc_rhiz_pvsn.man)

## PSF vs. Common ##
cc_rhiz_pvsc.met <- filter(cc_rhiz.met, Soil_Treatment != "Non-PSF Soil")
cc_rhiz_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, cc_rhiz_pvsc.met)
summary(cc_rhiz_pvsc.man)

## Non-PSF vs. Common ##
cc_rhiz_nvsc.met <- filter(cc_rhiz.met, Soil_Treatment != "PSF Soil")
cc_rhiz_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, cc_rhiz_nvsc.met)
summary(cc_rhiz_nvsc.man)

## Desmodium ##
ds_rhiz.ps <- subset_samples(rhiz.ps, Plant == "D. illinoense")
ds_rhiz.ps <- subset_taxa(ds_rhiz.ps, taxa_sums(ds_rhiz.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(ds_rhiz.ps, 'ds_rhiz')
ds_rhiz_prop.ps <- transform_sample_counts(ds_rhiz.ps, function(x) x/sum(x))
set.seed(248)
ds_rhiz.wuni <- phyloseq::distance(ds_rhiz_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
ds_rhiz.bdis <- betadisper(ds_rhiz.wuni, group = ds_rhiz$met$Soils)
anova(ds_rhiz.bdis)
TukeyHSD(ds_rhiz.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
ds_rhiz.met <- filter(rhiz_nmds.load, Plant == "D. illinoense")

# Perform MANOVA on all samples #
ds_rhiz.man <- manova(cbind(NMDS1,NMDS2)~Soils, ds_rhiz.met)
summary(ds_rhiz.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
ds_rhiz_pvsn.met <- filter(ds_rhiz.met, Soil_Treatment != "Common Soil")
ds_rhiz_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, ds_rhiz_pvsn.met)
summary(ds_rhiz_pvsn.man)

## PSF vs. Common ##
ds_rhiz_pvsc.met <- filter(ds_rhiz.met, Soil_Treatment != "Non-PSF Soil")
ds_rhiz_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, ds_rhiz_pvsc.met)
summary(ds_rhiz_pvsc.man)

## Non-PSF vs. Common ##
ds_rhiz_nvsc.met <- filter(ds_rhiz.met, Soil_Treatment != "PSF Soil")
ds_rhiz_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, ds_rhiz_nvsc.met)
summary(ds_rhiz_nvsc.man)

## Hog Peanut ##
hp_rhiz.ps <- subset_samples(rhiz.ps, Plant == "A. bracteata")
hp_rhiz.ps <- subset_taxa(hp_rhiz.ps, taxa_sums(hp_rhiz.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(hp_rhiz.ps, 'hp_rhiz')
hp_rhiz_prop.ps <- transform_sample_counts(hp_rhiz.ps, function(x) x/sum(x))
set.seed(248)
hp_rhiz.wuni <- phyloseq::distance(hp_rhiz_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
hp_rhiz.bdis <- betadisper(hp_rhiz.wuni, group = hp_rhiz$met$Soils)
anova(hp_rhiz.bdis)
TukeyHSD(hp_rhiz.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
hp_rhiz.met <- filter(rhiz_nmds.load, Plant == "A. bracteata")

# Perform MANOVA on all samples #
hp_rhiz.man <- manova(cbind(NMDS1,NMDS2)~Soils, hp_rhiz.met)
summary(hp_rhiz.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
hp_rhiz_pvsn.met <- filter(hp_rhiz.met, Soil_Treatment != "Common Soil")
hp_rhiz_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, hp_rhiz_pvsn.met)
summary(hp_rhiz_pvsn.man)

## PSF vs. Common ##
hp_rhiz_pvsc.met <- filter(hp_rhiz.met, Soil_Treatment != "Non-PSF Soil")
hp_rhiz_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, hp_rhiz_pvsc.met)
summary(hp_rhiz_pvsc.man)

## Non-PSF vs. Common ##
hp_rhiz_nvsc.met <- filter(hp_rhiz.met, Soil_Treatment != "PSF Soil")
hp_rhiz_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, hp_rhiz_nvsc.met)
summary(hp_rhiz_nvsc.man)

## Clover ##
cl_rhiz.ps <- subset_samples(rhiz.ps, Plant == "T. repens")
cl_rhiz.ps <- subset_taxa(cl_rhiz.ps, taxa_sums(cl_rhiz.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(cl_rhiz.ps, 'cl_rhiz')
cl_rhiz_prop.ps <- transform_sample_counts(cl_rhiz.ps, function(x) x/sum(x))
set.seed(248)
cl_rhiz.wuni <- phyloseq::distance(cl_rhiz_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
cl_rhiz.bdis <- betadisper(cl_rhiz.wuni, group = cl_rhiz$met$Soils)
anova(cl_rhiz.bdis)
TukeyHSD(cl_rhiz.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
cl_rhiz.met <- filter(rhiz_nmds.load, Plant == "T. repens")

# Perform MANOVA on all samples #
cl_rhiz.man <- manova(cbind(NMDS1,NMDS2)~Soils, cl_rhiz.met)
summary(cl_rhiz.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Common ##
cl_rhiz_pvsc.met <- filter(cl_rhiz.met, Soil_Treatment != "Non-PSF Soil")
cl_rhiz_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, cl_rhiz_pvsc.met)
summary(cl_rhiz_pvsc.man)

## Medicago ##
md_rhiz.ps <- subset_samples(rhiz.ps, Plant == "M. truncatula")
md_rhiz.ps <- subset_taxa(md_rhiz.ps, taxa_sums(md_rhiz.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(md_rhiz.ps, 'md_rhiz')
md_rhiz_prop.ps <- transform_sample_counts(md_rhiz.ps, function(x) x/sum(x))
set.seed(248)
md_rhiz.wuni <- phyloseq::distance(md_rhiz_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
md_rhiz.bdis <- betadisper(md_rhiz.wuni, group = md_rhiz$met$Soils)
anova(md_rhiz.bdis)
TukeyHSD(md_rhiz.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
md_rhiz.met <- filter(rhiz_nmds.load, Plant == "M. truncatula")

# Perform MANOVA on all samples #
md_rhiz.man <- manova(cbind(NMDS1,NMDS2)~Soils, md_rhiz.met)
summary(md_rhiz.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
md_rhiz_pvsn.met <- filter(md_rhiz.met, Soil_Treatment != "Common Soil")
md_rhiz_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, md_rhiz_pvsn.met)
summary(md_rhiz_pvsn.man)

## PSF vs. Common ##
md_rhiz_pvsc.met <- filter(md_rhiz.met, Soil_Treatment != "Non-PSF Soil")
md_rhiz_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, md_rhiz_pvsc.met)
summary(md_rhiz_pvsc.man)

## Non-PSF vs. Common ##
md_rhiz_nvsc.met <- filter(md_rhiz.met, Soil_Treatment != "PSF Soil")
md_rhiz_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, md_rhiz_nvsc.met)
summary(md_rhiz_nvsc.man)


# Root Endosphere #
# Construct a Weighted Unifrac distance matrix and PCoA ordination #
root_prop.ps <- transform_sample_counts(root.ps, function(x) x/sum(x))

set.seed(248)
root.wuni <- phyloseq::distance(root_prop.ps, method = 'wunifrac')
root.pcoa <- phyloseq::ordinate(root_prop.ps, 'PCoA', distance = root.wuni)

# Perform an NMDS analysis using the weighted Unifrac distance matrix, with the PCoA ordination as the starting ordination # 
root.nmds <- metaMDS(root.wuni, 
                     k = 11, try = 20, trymax = 1000, maxit = 999,
                     model = 'global', 
                     autotransform = FALSE, previous.best = root.pcoa$vectors[,1:5])

# Save the loading scores for all axes and make a distance matrix from these scores #
root_nmds.scores <- scores(root.nmds, display = 'sites')
root_nmds.dist <- dist(root_nmds.scores)

# Fit a linear model using the vectorized weighted Unifrac ditsance matrix as a predictor of the vectorized loading score distance matrix to fine total R^2 of the model # 
root_nmds.ffit <- lm(as.vector(root_nmds.dist) ~ as.vector(root.wuni))
summary(root_nmds.ffit)
root_nmds.totr2 <- summary(root_nmds.ffit)$r.squared

# Axes Variance Calculation #
# Fit linear models as before expect to predict the distance matrix of each individual axis #
root_nmds.dist1 <- dist(root_nmds.scores[,1])
root_nmds.fit1 <- lm(as.vector(root_nmds.dist1)~as.vector(root.wuni))
root_nmds.r1 <- summary(root_nmds.fit1)$r.squared

root_nmds.dist2 <- dist(root_nmds.scores[,2])
root_nmds.fit2 <- lm(as.vector(root_nmds.dist2)~as.vector(root.wuni))
root_nmds.r2 <- summary(root_nmds.fit2)$r.squared

root_nmds.dist3 <- dist(root_nmds.scores[,3])
root_nmds.fit3 <- lm(as.vector(root_nmds.dist3)~as.vector(root.wuni))
root_nmds.r3 <- summary(root_nmds.fit3)$r.squared

root_nmds.dist4 <- dist(root_nmds.scores[,4])
root_nmds.fit4 <- lm(as.vector(root_nmds.dist4)~as.vector(root.wuni))
root_nmds.r4 <- summary(root_nmds.fit4)$r.squared

root_nmds.dist5 <- dist(root_nmds.scores[,5])
root_nmds.fit5 <- lm(as.vector(root_nmds.dist5)~as.vector(root.wuni))
root_nmds.r5 <- summary(root_nmds.fit5)$r.squared

root_nmds.dist6 <- dist(root_nmds.scores[,6])
root_nmds.fit6 <- lm(as.vector(root_nmds.dist6)~as.vector(root.wuni))
root_nmds.r6 <- summary(root_nmds.fit6)$r.squared

root_nmds.dist7 <- dist(root_nmds.scores[,7])
root_nmds.fit7 <- lm(as.vector(root_nmds.dist7)~as.vector(root.wuni))
root_nmds.r7 <- summary(root_nmds.fit7)$r.squared

root_nmds.dist8 <- dist(root_nmds.scores[,8])
root_nmds.fit8 <- lm(as.vector(root_nmds.dist8)~as.vector(root.wuni))
root_nmds.r8 <- summary(root_nmds.fit8)$r.squared

root_nmds.dist9 <- dist(root_nmds.scores[,9])
root_nmds.fit9 <- lm(as.vector(root_nmds.dist9)~as.vector(root.wuni))
root_nmds.r9 <- summary(root_nmds.fit9)$r.squared

root_nmds.dist10 <- dist(root_nmds.scores[,10])
root_nmds.fit10 <- lm(as.vector(root_nmds.dist10)~as.vector(root.wuni))
root_nmds.r10 <- summary(root_nmds.fit10)$r.squared

root_nmds.dist11 <- dist(root_nmds.scores[,11])
root_nmds.fit11 <- lm(as.vector(root_nmds.dist11)~as.vector(root.wuni))
root_nmds.r11 <- summary(root_nmds.fit11)$r.squared

# Take the sum the R^2 value from each axis #
root_nmds.comb <- root_nmds.r1 + root_nmds.r2 + root_nmds.r3 + root_nmds.r4 + root_nmds.r5 + root_nmds.r6 + root_nmds.r7 + root_nmds.r8 + root_nmds.r9 + root_nmds.r10 + root_nmds.r11

# Divide each axis R^2 by the total of all axes and then multiply by the variation explained by the whole model
root_nmds.axisr <- c()
for(i in 1:ncol(root_nmds.scores)){
  root_nmds.axisr[i] <- (get(paste0('root_nmds.r', i)) / root_nmds.comb) * root_nmds.totr2 
}

### Calculating Variance Components ###
# Construct a data.frame that has sample info and their loading scores #
decompose_ps(root.ps, 'root')
for(i in 1:nrow(root$met)){
  root$met$PS[i] <- paste0(substr(root$met$Plant.Species[i],1,1), substr(root$met$Soil.Origin[i],1,1))
}
root$met$Plants <- factor(root$met$Plant.Species, levels = c('S. helvola', 'C. fasciculata', 'D. illinoense', 'A. bracteata', 'T. repens', 'M. truncatula'))
root$met$Comps <- factor(root$met$Compartment, c('Root Endosphere'))
root$met$Soils <- factor(root$met$Soil.Origin, levels = c('Common Soil', "Non-PSF Soil", "PSF Soil"))

sample_data(root.ps) <- root$met
root_nmds.load <- cbind(root$met, root_nmds.scores)

# Test the mixed linear model on the first NMDS axis #
root_nmds.vfit1 <- lmer(NMDS1 ~ (1|Soils) + (1|Plants) +(1|Soils:Plants), data = root_nmds.load, REML = TRUE)
summary(root_nmds.vfit1)
root_nmds.vca1 <- as.data.frame(VarCorr(root_nmds.vfit1))

# Using Loop to do each NMDS axis #
root_nmds.vca <- matrix(nrow = 4, ncol = ncol(root_nmds.scores))
hold <- c()
for(i in 1:ncol(root_nmds.scores)){
  hold <- lmer(root_nmds.scores[,i] ~  (1|Soils) + (1|Plants) +(1|Soils:Plants), data = root_nmds.load, REML = TRUE)
  hold <- as.data.frame(VarCorr(hold))
  root_nmds.vca[1,i] <- hold[1,4]
  root_nmds.vca[2,i] <- hold[2,4]
  root_nmds.vca[3,i] <- hold[3,4]
  root_nmds.vca[4,i] <- hold[4,4]
}

# Save the variance components to their assigned variable/variable interaction and their NMDS loading axis #
rownames(root_nmds.vca) <- c('PSF History x Plant', 'Plant', 'PSF History', 'Residual')
colnames(root_nmds.vca) <- colnames(root_nmds.scores)

# Calculate the total variance of each variance component#
root_nmds.vtot <- colSums(root_nmds.vca)

# Weight each variance component by the amount of variation each axis explains #
root_nmds.wvca <- matrix(nrow = nrow(root_nmds.vca), ncol = length(root_nmds.axisr))
for(i in 1:length(root_nmds.axisr)){
  for(j in 1:nrow(root_nmds.vca)){
    root_nmds.wvca[j,i] <- root_nmds.vca[j,i]*root_nmds.axisr[i] 
  }
}
# Take the total variance explained by each predictor and take the sum of those values #
rownames(root_nmds.wvca) <- rownames(root_nmds.vca); colnames(root_nmds.wvca) <- colnames(root_nmds.vca)
root_nmds.tvca <- rowSums(root_nmds.wvca)
root_nmds.tvca <- as.data.frame(root_nmds.tvca)
root_nmds.ptot <- colSums(root_nmds.tvca)

# Take the variance explained by each predictor and divide by the total variance explained and multiply by 100% #
root_nmds.pvca <- matrix(nrow = nrow(root_nmds.tvca), ncol = 1)
for(i in 1:nrow(root_nmds.tvca)){
  root_nmds.pvca[i,1] <- root_nmds.tvca[i,1] / root_nmds.ptot * 100
}

# Save the variation explained percentages of each predictor/predictor interaction #
rownames(root_nmds.pvca) <- rownames(root_nmds.vca); colnames(root_nmds.pvca) <- 'Variance Explained'
root_nmds.pvca

# Perform PermANOVA using all samples #
root.adon <- adonis2(root.wuni~Plants*Soils, root$met, permutations = 999)
root.adon_by <- adonis2(root.wuni~Plants*Soils, root$met, permutations = 999, by = 'terms')
root.adon
root.adon_by

# Perform PermDISP using all samples # 
root.bdis <- betadisper(root.wuni, group = root$met$PS)
anova(root.bdis)
TukeyHSD(root.bdis)

# Make an object with the data to be plotted #
root_nmds.load <- cbind(root$met, root_nmds.scores)

# plot the NMDS ordination of all samples that will be used to make a patchwork plot #
root_nmds.plot <- ggplot(root_nmds.load, aes(NMDS1, NMDS2, color = Plants, shape = Soils)) +
  geom_point(size = 8) +
  scale_color_manual(name = "Root Endosphere Community Source", labels = c(expression(italic('S. helvola')), expression(italic('C. fasciculata')), expression(italic('D. illinoense')), expression(italic('A. bracteata')), expression(italic('T. repens')), expression(italic('M. truncatula')), "Common Soil"), values = c("#A6CEE3","#1F78B4","#FDBF6F", "#FF7F00","#CAB2D6", "#6A3D9A", "#B15928")) +
  scale_shape_manual(name = "PSF History", labels = c("Common Soil", "Non-PSF Soil", "PSF Soil"), values = c(15,17,16)) +
  xlab(expression("NMDS1 (R"^2 ~ "= 0.5695)")) +
  ylab(expression("NMDS2 (R"^2 ~ "= 0.1246)")) +
  theme_bw() +
  theme(legend.position = 'right',
        panel.grid = element_blank(),
        legend.text = element_text(size = 12, family = "Liberation Sans"),
        legend.title = element_text(size = 16, family = "Liberation Sans"),
        axis.title = element_text(size=20, family = "Liberation Sans"),
        axis.text = element_text(color = "black", size = 8, family = "Liberation Sans")) +
  annotate('text', x = 0.08, y = -0.3,
           label = expression("(PERMANOVA) F"["17,62"] ~ "= 6.4871, P = 0.001;  (PERMDISP) F"["17,62"] ~  "= 0.9935, P = 0.4821; 3D Stress = 0.0188"),
           size = 8, family = 'Liberation Sans') +
  coord_cartesian(ylim = c(-0.3, 0.25))

root_nmds.plot

# Separate samples by plant species #
## Fuzzy Bean ##
fb_root.ps <- subset_samples(root.ps, Plant.Species == "S. helvola")
fb_root.ps <- subset_taxa(fb_root.ps, taxa_sums(fb_root.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(fb_root.ps, 'fb_root')
fb_root_prop.ps <- transform_sample_counts(fb_root.ps, function(x) x/sum(x))
set.seed(248)
fb_root.wuni <- phyloseq::distance(fb_root_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
fb_root.bdis <- betadisper(fb_root.wuni, group = fb_root$met$Soils)
anova(fb_root.bdis)
TukeyHSD(fb_root.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
fb_root.met <- filter(root_nmds.load, Plant.Species == "S. helvola")

# Perform MANOVA on all samples #
fb_root.man <- manova(cbind(NMDS1,NMDS2)~Soils, fb_root.met)
summary(fb_root.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
fb_root_pvsn.met <- filter(fb_root.met, Soil.Origin != "Common Soil")
fb_root_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, fb_root_pvsn.met)
summary(fb_root_pvsn.man)

## PSF vs. Common ##
fb_root_pvsc.met <- filter(fb_root.met, Soil.Origin != "Non-PSF Soil")
fb_root_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, fb_root_pvsc.met)
summary(fb_root_pvsc.man)

## Non-PSF vs. Common ##
fb_root_nvsc.met <- filter(fb_root.met, Soil.Origin != "PSF Soil")
fb_root_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, fb_root_nvsc.met)
summary(fb_root_nvsc.man)

## Chamaecrista ##
cc_root.ps <- subset_samples(root.ps, Plant.Species == "C. fasciculata")
cc_root.ps <- subset_taxa(cc_root.ps, taxa_sums(cc_root.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(cc_root.ps, 'cc_root')
cc_root_prop.ps <- transform_sample_counts(cc_root.ps, function(x) x/sum(x))
set.seed(248)
cc_root.wuni <- phyloseq::distance(cc_root_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
cc_root.bdis <- betadisper(cc_root.wuni, group = cc_root$met$Soils)
anova(cc_root.bdis)
TukeyHSD(cc_root.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
cc_root.met <- filter(root_nmds.load, Plant.Species == "C. fasciculata")

# Perform MANOVA on all samples #
cc_root.man <- manova(cbind(NMDS1,NMDS2)~Soils, cc_root.met)
summary(cc_root.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
cc_root_pvsn.met <- filter(cc_root.met, Soil.Origin != "Common Soil")
cc_root_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, cc_root_pvsn.met)
summary(cc_root_pvsn.man)

## PSF vs. Common ##
cc_root_pvsc.met <- filter(cc_root.met, Soil.Origin != "Non-PSF Soil")
cc_root_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, cc_root_pvsc.met)
summary(cc_root_pvsc.man)

## Non-PSF vs. Common ##
cc_root_nvsc.met <- filter(cc_root.met, Soil.Origin != "PSF Soil")
cc_root_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, cc_root_nvsc.met)
summary(cc_root_nvsc.man)

## Desmodium ##
ds_root.ps <- subset_samples(root.ps, Plant.Species == "D. illinoense")
ds_root.ps <- subset_taxa(ds_root.ps, taxa_sums(ds_root.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(ds_root.ps, 'ds_root')
ds_root_prop.ps <- transform_sample_counts(ds_root.ps, function(x) x/sum(x))
set.seed(248)
ds_root.wuni <- phyloseq::distance(ds_root_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
ds_root.bdis <- betadisper(ds_root.wuni, group = ds_root$met$Soils)
anova(ds_root.bdis)
TukeyHSD(ds_root.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
ds_root.met <- filter(root_nmds.load, Plant.Species == "D. illinoense")

# Perform MANOVA on all samples #
ds_root.man <- manova(cbind(NMDS1,NMDS2)~Soils, ds_root.met)
summary(ds_root.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
ds_root_pvsn.met <- filter(ds_root.met, Soil.Origin != "Common Soil")
ds_root_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, ds_root_pvsn.met)
summary(ds_root_pvsn.man)

## PSF vs. Common ##
ds_root_pvsc.met <- filter(ds_root.met, Soil.Origin != "Non-PSF Soil")
ds_root_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, ds_root_pvsc.met)
summary(ds_root_pvsc.man)

## Non-PSF vs. Common ##
ds_root_nvsc.met <- filter(ds_root.met, Soil.Origin != "PSF Soil")
ds_root_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, ds_root_nvsc.met)
summary(ds_root_nvsc.man)

## Hog Peanut ##
hp_root.ps <- subset_samples(root.ps, Plant.Species == "A. bracteata")
hp_root.ps <- subset_taxa(hp_root.ps, taxa_sums(hp_root.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(hp_root.ps, 'hp_root')
hp_root_prop.ps <- transform_sample_counts(hp_root.ps, function(x) x/sum(x))
set.seed(248)
hp_root.wuni <- phyloseq::distance(hp_root_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
hp_root.bdis <- betadisper(hp_root.wuni, group = hp_root$met$Soils)
anova(hp_root.bdis)
TukeyHSD(hp_root.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
hp_root.met <- filter(root_nmds.load, Plant.Species == "A. bracteata")

# Perform MANOVA on all samples #
hp_root.man <- manova(cbind(NMDS1,NMDS2)~Soils, hp_root.met)
summary(hp_root.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
hp_root_pvsn.met <- filter(hp_root.met, Soil.Origin != "Common Soil")
hp_root_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, hp_root_pvsn.met)
summary(hp_root_pvsn.man)

## PSF vs. Common ##
hp_root_pvsc.met <- filter(hp_root.met, Soil.Origin != "Non-PSF Soil")
hp_root_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, hp_root_pvsc.met)
summary(hp_root_pvsc.man)

## Non-PSF vs. Common ##
hp_root_nvsc.met <- filter(hp_root.met, Soil.Origin != "PSF Soil")
hp_root_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, hp_root_nvsc.met)
summary(hp_root_nvsc.man)

## Clover ##
cl_root.ps <- subset_samples(root.ps, Plant.Species == "T. repens")
cl_root.ps <- subset_taxa(cl_root.ps, taxa_sums(cl_root.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(cl_root.ps, 'cl_root')
cl_root_prop.ps <- transform_sample_counts(cl_root.ps, function(x) x/sum(x))
set.seed(248)
cl_root.wuni <- phyloseq::distance(cl_root_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
cl_root.bdis <- betadisper(cl_root.wuni, group = cl_root$met$Soils)
anova(cl_root.bdis)
TukeyHSD(cl_root.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
cl_root.met <- filter(root_nmds.load, Plant.Species == "T. repens")

# Perform MANOVA on all samples #
cl_root.man <- manova(cbind(NMDS1,NMDS2)~Soils, cl_root.met)
summary(cl_root.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
cl_root_pvsn.met <- filter(cl_root.met, Soil.Origin != "Common Soil")
cl_root_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, cl_root_pvsn.met)
summary(cl_root_pvsn.man)

## PSF vs. Common ##
cl_root_pvsc.met <- filter(cl_root.met, Soil.Origin != "Non-PSF Soil")
cl_root_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, cl_root_pvsc.met)
summary(cl_root_pvsc.man)

## Non-PSF vs. Common ##
cl_root_nvsc.met <- filter(cl_root.met, Soil.Origin != "PSF Soil")
cl_root_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, cl_root_nvsc.met)
summary(cl_root_nvsc.man)

## Medicago ##
md_root.ps <- subset_samples(root.ps, Plant.Species == "M. truncatula")
md_root.ps <- subset_taxa(md_root.ps, taxa_sums(md_root.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(md_root.ps, 'md_root')
md_root_prop.ps <- transform_sample_counts(md_root.ps, function(x) x/sum(x))
set.seed(248)
md_root.wuni <- phyloseq::distance(md_root_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
md_root.bdis <- betadisper(md_root.wuni, group = md_root$met$Soils)
anova(md_root.bdis)
TukeyHSD(md_root.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
md_root.met <- filter(root_nmds.load, Plant.Species == "M. truncatula")

# Perform MANOVA on all samples #
md_root.man <- manova(cbind(NMDS1,NMDS2)~Soils, md_root.met)
summary(md_root.man)

# Subset the data for pairwise MANOVA tests #
## PSF vs. Non-PSF ##
md_root_pvsn.met <- filter(md_root.met, Soil.Origin != "Common Soil")
md_root_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Soils, md_root_pvsn.met)
summary(md_root_pvsn.man)

## PSF vs. Common ##
md_root_pvsc.met <- filter(md_root.met, Soil.Origin != "Non-PSF Soil")
md_root_pvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, md_root_pvsc.met)
summary(md_root_pvsc.man)

## Non-PSF vs. Common ##
md_root_nvsc.met <- filter(md_root.met, Soil.Origin != "PSF Soil")
md_root_nvsc.man <- manova(cbind(NMDS1,NMDS2)~Soils, md_root_nvsc.met)
summary(md_root_nvsc.man)

# Root Endosphere for Natives #
# Construct a Weighted Unifrac distance matrix and PCoA ordination #
raty.ps <- subset_samples(root.ps, Plant.Species != "T. repens")
raty.ps <- subset_samples(raty.ps, Plant.Species != "M. truncatula") 
raty_prop.ps <- transform_sample_counts(raty.ps, function(x) x/sum(x))

set.seed(248)
raty.wuni <- phyloseq::distance(raty_prop.ps, method = 'wunifrac')
raty.pcoa <- phyloseq::ordinate(raty_prop.ps, 'PCoA', distance = raty.wuni)

# Perform an NMDS analysis using the weighted Unifrac distance matrix, with the PCoA ordination as the starting ordination # 
raty.nmds <- metaMDS(raty.wuni, 
                     k = 8, try = 100, trymax = 1000, maxit = 999,
                     model = 'global', 
                     autotransform = FALSE, previous.best = raty.pcoa$vectors[,1:8])

# Save the loading scores for all axes and make a distance matrix from these scores #
raty_nmds.scores <- scores(raty.nmds, display = 'sites')
raty_nmds.dist <- dist(raty_nmds.scores)

# Fit a linear model using the vectorized weighted Unifrac ditsance matrix as a predictor of the vectorized loading score distance matrix to fine total R^2 of the model # 
raty_nmds.ffit <- lm(as.vector(raty_nmds.dist) ~ as.vector(raty.wuni))
summary(raty_nmds.ffit)
raty_nmds.totr2 <- summary(raty_nmds.ffit)$r.squared

# Axes Variance Calculation #
# Fit linear models as before expect to predict the distance matrix of each individual axis #
raty_nmds.dist1 <- dist(raty_nmds.scores[,1])
raty_nmds.fit1 <- lm(as.vector(raty_nmds.dist1)~as.vector(raty.wuni))
raty_nmds.r1 <- summary(raty_nmds.fit1)$r.squared

raty_nmds.dist2 <- dist(raty_nmds.scores[,2])
raty_nmds.fit2 <- lm(as.vector(raty_nmds.dist2)~as.vector(raty.wuni))
raty_nmds.r2 <- summary(raty_nmds.fit2)$r.squared

raty_nmds.dist3 <- dist(raty_nmds.scores[,3])
raty_nmds.fit3 <- lm(as.vector(raty_nmds.dist3)~as.vector(raty.wuni))
raty_nmds.r3 <- summary(raty_nmds.fit3)$r.squared

raty_nmds.dist4 <- dist(raty_nmds.scores[,4])
raty_nmds.fit4 <- lm(as.vector(raty_nmds.dist4)~as.vector(raty.wuni))
raty_nmds.r4 <- summary(raty_nmds.fit4)$r.squared

raty_nmds.dist5 <- dist(raty_nmds.scores[,5])
raty_nmds.fit5 <- lm(as.vector(raty_nmds.dist5)~as.vector(raty.wuni))
raty_nmds.r5 <- summary(raty_nmds.fit5)$r.squared

raty_nmds.dist6 <- dist(raty_nmds.scores[,6])
raty_nmds.fit6 <- lm(as.vector(raty_nmds.dist6)~as.vector(raty.wuni))
raty_nmds.r6 <- summary(raty_nmds.fit6)$r.squared

raty_nmds.dist7 <- dist(raty_nmds.scores[,7])
raty_nmds.fit7 <- lm(as.vector(raty_nmds.dist7)~as.vector(raty.wuni))
raty_nmds.r7 <- summary(raty_nmds.fit7)$r.squared

raty_nmds.dist8 <- dist(raty_nmds.scores[,8])
raty_nmds.fit8 <- lm(as.vector(raty_nmds.dist8)~as.vector(raty.wuni))
raty_nmds.r8 <- summary(raty_nmds.fit8)$r.squared

# Take the sum the R^2 value from each axis #
raty_nmds.comb <- raty_nmds.r1 + raty_nmds.r2 + raty_nmds.r3 + raty_nmds.r4 + raty_nmds.r5 + raty_nmds.r6 + raty_nmds.r7 + raty_nmds.r8

# Divide each axis R^2 by the total of all axes and then multiply by the variation explained by the whole model
raty_nmds.axisr <- c()
for(i in 1:ncol(raty_nmds.scores)){
  raty_nmds.axisr[i] <- (get(paste0('raty_nmds.r', i)) / raty_nmds.comb) * raty_nmds.totr2 
}

### Calculating Variance Components ###
# Construct a data.frame that has sample info and their loading scores #
decompose_ps(raty.ps, 'raty')
for(i in 1:nrow(raty$met)){
  raty$met$PS[i] <- paste0(substr(raty$met$Plant.Species[i],1,1), substr(raty$met$Soil.Origin[i],1,1))
}
raty$met$Plants <- factor(raty$met$Plant.Species, levels = c('S. helvola', 'C. fasciculata', 'D. illinoense', 'A. bracteata', 'T. repens', 'M. truncatula'))
raty$met$Comps <- factor(raty$met$Compartment, c('Root Endosphere'))
raty$met$Soils <- factor(raty$met$Soil.Origin, levels = c('Common Soil', "Non-PSF Soil", "PSF Soil"))

sample_data(raty.ps) <- raty$met
raty_nmds.load <- cbind(raty$met, raty_nmds.scores)

# Test the mixed linear model on the first NMDS axis #
raty_nmds.vfit1 <- lmer(NMDS1 ~ (1|Soils) + (1|Plants) +(1|Soils:Plants), data = raty_nmds.load, REML = TRUE)
summary(raty_nmds.vfit1)
raty_nmds.vca1 <- as.data.frame(VarCorr(raty_nmds.vfit1))

# Using Loop to do each NMDS axis #
raty_nmds.vca <- matrix(nrow = 4, ncol = ncol(raty_nmds.scores))
hold <- c()
for(i in 1:ncol(raty_nmds.scores)){
  hold <- lmer(raty_nmds.scores[,i] ~  (1|Soils) + (1|Plants) +(1|Soils:Plants), data = raty_nmds.load, REML = TRUE)
  hold <- as.data.frame(VarCorr(hold))
  raty_nmds.vca[1,i] <- hold[1,4]
  raty_nmds.vca[2,i] <- hold[2,4]
  raty_nmds.vca[3,i] <- hold[3,4]
  raty_nmds.vca[4,i] <- hold[4,4]
}

# Save the variance components to their assigned variable/variable interaction and their NMDS loading axis #
rownames(raty_nmds.vca) <- c('PSF History x Plant', 'Plant', 'PSF History', 'Residual')
colnames(raty_nmds.vca) <- colnames(raty_nmds.scores)

# Calculate the total variance of each variance component#
raty_nmds.vtot <- colSums(raty_nmds.vca)

# Weight each variance component by the amount of variation each axis explains #
raty_nmds.wvca <- matrix(nrow = nrow(raty_nmds.vca), ncol = length(raty_nmds.axisr))
for(i in 1:length(raty_nmds.axisr)){
  for(j in 1:nrow(raty_nmds.vca)){
    raty_nmds.wvca[j,i] <- raty_nmds.vca[j,i]*raty_nmds.axisr[i] 
  }
}
# Take the total variance explained by each predictor and take the sum of those values #
rownames(raty_nmds.wvca) <- rownames(raty_nmds.vca); colnames(raty_nmds.wvca) <- colnames(raty_nmds.vca)
raty_nmds.tvca <- rowSums(raty_nmds.wvca)
raty_nmds.tvca <- as.data.frame(raty_nmds.tvca)
raty_nmds.ptot <- colSums(raty_nmds.tvca)

# Take the variance explained by each predictor and divide by the total variance explained and multiply by 100% #
raty_nmds.pvca <- matrix(nrow = nrow(raty_nmds.tvca), ncol = 1)
for(i in 1:nrow(raty_nmds.tvca)){
  raty_nmds.pvca[i,1] <- raty_nmds.tvca[i,1] / raty_nmds.ptot * 100
}

# Save the variation explained percentages of each predictor/predictor interaction #
rownames(raty_nmds.pvca) <- rownames(raty_nmds.vca); colnames(raty_nmds.pvca) <- 'Variance Explained'
raty_nmds.pvca

# Perform PermANOVA using all samples #
raty.adon <- adonis2(raty.wuni~Plants*Soils, raty$met, permutations = 999)
raty.adon_by <- adonis2(raty.wuni~Plants*Soils, raty$met, permutations = 999, by = 'terms')
raty.adon
raty.adon_by

# Perform PermDISP using all samples # 
raty.bdis <- betadisper(raty.wuni, group = raty$met$PS)
anova(raty.bdis)
TukeyHSD(raty.bdis)

# Root Endosphere for Non-Natives #
# Construct a Weighted Unifrac distance matrix and PCoA ordination #
rnat.ps <- subset_samples(root.ps, Plant.Species == "T. repens" | Plant.Species == "M. truncatula")
rnat_prop.ps <- transform_sample_counts(rnat.ps, function(x) x/sum(x))

set.seed(248)
rnat.wuni <- phyloseq::distance(rnat_prop.ps, method = 'wunifrac')
rnat.pcoa <- phyloseq::ordinate(rnat_prop.ps, 'PCoA', distance = rnat.wuni)

# Perform an NMDS analysis using the weighted Unifrac distance matrix, with the PCoA ordination as the starting ordination # 
rnat.nmds <- metaMDS(rnat.wuni, 
                     k = 5, try = 100, trymax = 1000, maxit = 999,
                     model = 'global', 
                     autotransform = FALSE, previous.best = rnat.pcoa$vectors[,1:5])

# Save the loading scores for all axes and make a distance matrix from these scores #
rnat_nmds.scores <- scores(rnat.nmds, display = 'sites')
rnat_nmds.dist <- dist(rnat_nmds.scores)

# Fit a linear model using the vectorized weighted Unifrac ditsance matrix as a predictor of the vectorized loading score distance matrix to fine total R^2 of the model # 
rnat_nmds.ffit <- lm(as.vector(rnat_nmds.dist) ~ as.vector(rnat.wuni))
summary(rnat_nmds.ffit)
rnat_nmds.totr2 <- summary(rnat_nmds.ffit)$r.squared

# Axes Variance Calculation #
# Fit linear models as before expect to predict the distance matrix of each individual axis #
rnat_nmds.dist1 <- dist(rnat_nmds.scores[,1])
rnat_nmds.fit1 <- lm(as.vector(rnat_nmds.dist1)~as.vector(rnat.wuni))
rnat_nmds.r1 <- summary(rnat_nmds.fit1)$r.squared

rnat_nmds.dist2 <- dist(rnat_nmds.scores[,2])
rnat_nmds.fit2 <- lm(as.vector(rnat_nmds.dist2)~as.vector(rnat.wuni))
rnat_nmds.r2 <- summary(rnat_nmds.fit2)$r.squared

rnat_nmds.dist3 <- dist(rnat_nmds.scores[,3])
rnat_nmds.fit3 <- lm(as.vector(rnat_nmds.dist3)~as.vector(rnat.wuni))
rnat_nmds.r3 <- summary(rnat_nmds.fit3)$r.squared

rnat_nmds.dist4 <- dist(rnat_nmds.scores[,4])
rnat_nmds.fit4 <- lm(as.vector(rnat_nmds.dist4)~as.vector(rnat.wuni))
rnat_nmds.r4 <- summary(rnat_nmds.fit4)$r.squared

rnat_nmds.dist5 <- dist(rnat_nmds.scores[,5])
rnat_nmds.fit5 <- lm(as.vector(rnat_nmds.dist5)~as.vector(rnat.wuni))
rnat_nmds.r5 <- summary(rnat_nmds.fit5)$r.squared

# Take the sum the R^2 value from each axis #
rnat_nmds.comb <- rnat_nmds.r1 + rnat_nmds.r2 + rnat_nmds.r3 + rnat_nmds.r4 + rnat_nmds.r5

# Divide each axis R^2 by the total of all axes and then multiply by the variation explained by the whole model
rnat_nmds.axisr <- c()
for(i in 1:ncol(rnat_nmds.scores)){
  rnat_nmds.axisr[i] <- (get(paste0('rnat_nmds.r', i)) / rnat_nmds.comb) * rnat_nmds.totr2 
}

### Calculating Variance Components ###
# Construct a data.frame that has sample info and their loading scores #
decompose_ps(rnat.ps, 'rnat')
for(i in 1:nrow(rnat$met)){
  rnat$met$PS[i] <- paste0(substr(rnat$met$Plant.Species[i],1,1), substr(rnat$met$Soil.Origin[i],1,1))
}
rnat$met$Plants <- factor(rnat$met$Plant.Species, levels = c('S. helvola', 'C. fasciculata', 'D. illinoense', 'A. bracteata', 'T. repens', 'M. truncatula'))
rnat$met$Comps <- factor(rnat$met$Compartment, c('Root Endosphere'))
rnat$met$Soils <- factor(rnat$met$Soil.Origin, levels = c('Common Soil', "Non-PSF Soil", "PSF Soil"))

sample_data(rnat.ps) <- rnat$met
rnat_nmds.load <- cbind(rnat$met, rnat_nmds.scores)

# Test the mixed linear model on the first NMDS axis #
rnat_nmds.vfit1 <- lmer(NMDS1 ~ (1|Soils) + (1|Plants) +(1|Soils:Plants), data = rnat_nmds.load, REML = TRUE)
summary(rnat_nmds.vfit1)
rnat_nmds.vca1 <- as.data.frame(VarCorr(rnat_nmds.vfit1))

# Using Loop to do each NMDS axis #
rnat_nmds.vca <- matrix(nrow = 4, ncol = ncol(rnat_nmds.scores))
hold <- c()
for(i in 1:ncol(rnat_nmds.scores)){
  hold <- lmer(rnat_nmds.scores[,i] ~  (1|Soils) + (1|Plants) +(1|Soils:Plants), data = rnat_nmds.load, REML = TRUE)
  hold <- as.data.frame(VarCorr(hold))
  rnat_nmds.vca[1,i] <- hold[1,4]
  rnat_nmds.vca[2,i] <- hold[2,4]
  rnat_nmds.vca[3,i] <- hold[3,4]
  rnat_nmds.vca[4,i] <- hold[4,4]
}

# Save the variance components to their assigned variable/variable interaction and their NMDS loading axis #
rownames(rnat_nmds.vca) <- c('PSF History x Plant', 'PSF History', 'Plant', 'Residual')
colnames(rnat_nmds.vca) <- colnames(rnat_nmds.scores)

# Calculate the total variance of each variance component#
rnat_nmds.vtot <- colSums(rnat_nmds.vca)

# Weight each variance component by the amount of variation each axis explains #
rnat_nmds.wvca <- matrix(nrow = nrow(rnat_nmds.vca), ncol = length(rnat_nmds.axisr))
for(i in 1:length(rnat_nmds.axisr)){
  for(j in 1:nrow(rnat_nmds.vca)){
    rnat_nmds.wvca[j,i] <- rnat_nmds.vca[j,i]*rnat_nmds.axisr[i] 
  }
}
# Take the total variance explained by each predictor and take the sum of those values #
rownames(rnat_nmds.wvca) <- rownames(rnat_nmds.vca); colnames(rnat_nmds.wvca) <- colnames(rnat_nmds.vca)
rnat_nmds.tvca <- rowSums(rnat_nmds.wvca)
rnat_nmds.tvca <- as.data.frame(rnat_nmds.tvca)
rnat_nmds.ptot <- colSums(rnat_nmds.tvca)

# Take the variance explained by each predictor and divide by the total variance explained and multiply by 100% #
rnat_nmds.pvca <- matrix(nrow = nrow(rnat_nmds.tvca), ncol = 1)
for(i in 1:nrow(rnat_nmds.tvca)){
  rnat_nmds.pvca[i,1] <- rnat_nmds.tvca[i,1] / rnat_nmds.ptot * 100
}

# Save the variation explained percentages of each predictor/predictor interaction #
rownames(rnat_nmds.pvca) <- rownames(rnat_nmds.vca); colnames(rnat_nmds.pvca) <- 'Variance Explained'
rnat_nmds.pvca

# Perform PermANOVA using all samples #
rnat.adon <- adonis2(rnat.wuni~Plants*Soils, rnat$met, permutations = 999)
rnat.adon_by <- adonis2(rnat.wuni~Plants*Soils, rnat$met, permutations = 999, by = 'terms')
rnat.adon
rnat.adon_by

# Perform PermDISP using all samples # 
rnat.bdis <- betadisper(rnat.wuni, group = rnat$met$PS)
anova(rnat.bdis)
TukeyHSD(rnat.bdis)

# Bulk Soil vs Rhizosphere No Plant #
# Construct a phyloseq object containing only samples in the Non-PSF Soil #
npsf.ps <- subset_samples(soil.ps, Soil_Treatment == "Non-PSF Soil")

# Construct a Weighted Unifrac distance matrix and PCoA ordination #
npsf_prop.ps <- transform_sample_counts(npsf.ps, function(x) x/sum(x))

set.seed(248)
npsf.wuni <- phyloseq::distance(npsf_prop.ps, method = 'wunifrac')
npsf.pcoa <- phyloseq::ordinate(npsf_prop.ps, 'PCoA', distance = npsf.wuni)

# Perform an NMDS analysis using the weighted Unifrac distance matrix, with the PCoA ordination as the starting ordination # 
npsf.nmds <- metaMDS(npsf.wuni, 
                     k = 4, try = 20, trymax = 1000, maxit = 999,
                     model = 'global', 
                     autotransform = FALSE, previous.best = npsf.pcoa$vectors[,1:5])

# Save the loading scores for all axes and make a distance matrix from these scores #
npsf_nmds.scores <- scores(npsf.nmds, display = 'sites')
npsf_nmds.dist <- dist(npsf_nmds.scores)

# Fit a linear model using the vectorized weighted Unifrac ditsance matrix as a predictor of the vectorized loading score distance matrix to fine total R^2 of the model # 
npsf_nmds.ffit <- lm(as.vector(npsf_nmds.dist) ~ as.vector(npsf.wuni))
summary(npsf_nmds.ffit)
npsf_nmds.totr2 <- summary(npsf_nmds.ffit)$r.squared

# Axes Variance Calculation #
# Fit linear models as before expect to predict the distance matrix of each individual axis #
npsf_nmds.dist1 <- dist(npsf_nmds.scores[,1])
npsf_nmds.fit1 <- lm(as.vector(npsf_nmds.dist1)~as.vector(npsf.wuni))
npsf_nmds.r1 <- summary(npsf_nmds.fit1)$r.squared

npsf_nmds.dist2 <- dist(npsf_nmds.scores[,2])
npsf_nmds.fit2 <- lm(as.vector(npsf_nmds.dist2)~as.vector(npsf.wuni))
npsf_nmds.r2 <- summary(npsf_nmds.fit2)$r.squared

npsf_nmds.dist3 <- dist(npsf_nmds.scores[,3])
npsf_nmds.fit3 <- lm(as.vector(npsf_nmds.dist3)~as.vector(npsf.wuni))
npsf_nmds.r3 <- summary(npsf_nmds.fit3)$r.squared

npsf_nmds.dist4 <- dist(npsf_nmds.scores[,4])
npsf_nmds.fit4 <- lm(as.vector(npsf_nmds.dist4)~as.vector(npsf.wuni))
npsf_nmds.r4 <- summary(npsf_nmds.fit4)$r.squared

# Take the sum the R^2 value from each axis #
npsf_nmds.comb <- npsf_nmds.r1 + npsf_nmds.r2 + npsf_nmds.r3 + npsf_nmds.r4

# Divide each axis R^2 by the total of all axes and then multiply by the variation explained by the whole model
npsf_nmds.axisr <- c()
for(i in 1:ncol(npsf_nmds.scores)){
  npsf_nmds.axisr[i] <- (get(paste0('npsf_nmds.r', i)) / npsf_nmds.comb) * npsf_nmds.totr2 
}

### Calculating Variance Components ###
# Construct a data.frame that has sample info and their loading scores #
decompose_ps(npsf.ps, 'npsf')
for(i in 1:nrow(npsf$met)){
  npsf$met$PC[i] <- paste0(substr(npsf$met$Plant[i],1,1), substr(npsf$met$Compartment[i],1,2))
}
npsf$met$Plants <- factor(npsf$met$Plant, levels = c('S. helvola', 'C. fasciculata', 'D. illinoense', 'A. bracteata', 'T. repens', 'M. truncatula'))
npsf$met$Comps <- factor(npsf$met$Compartment, c('Bulk Soil', 'Rhizosphere'))
npsf$met$Soils <- factor(npsf$met$Soil_Treatment, levels = c("Non-PSF Soil"))

sample_data(npsf.ps) <- npsf$met
npsf_nmds.load <- cbind(npsf$met, npsf_nmds.scores)

# Test the mixed linear model on the first NMDS axis #
npsf_nmds.vfit1 <- lmer(NMDS1 ~ (1|Comps) + (1|Plants) +(1|Comps:Plants), data = npsf_nmds.load, REML = TRUE)
summary(npsf_nmds.vfit1)
npsf_nmds.vca1 <- as.data.frame(VarCorr(npsf_nmds.vfit1))

# Using Loop to do each NMDS axis #
npsf_nmds.vca <- matrix(nrow = 4, ncol = ncol(npsf_nmds.scores))
hold <- c()
for(i in 1:ncol(npsf_nmds.scores)){
  hold <- lmer(npsf_nmds.scores[,i] ~  (1|Comps) + (1|Plants) +(1|Comps:Plants), data = npsf_nmds.load, REML = TRUE)
  hold <- as.data.frame(VarCorr(hold))
  npsf_nmds.vca[1,i] <- hold[1,4]
  npsf_nmds.vca[2,i] <- hold[2,4]
  npsf_nmds.vca[3,i] <- hold[3,4]
  npsf_nmds.vca[4,i] <- hold[4,4]
}

# Save the variance components to their assigned variable/variable interaction and their NMDS loading axis #
rownames(npsf_nmds.vca) <- c('PSF History x Plant', 'Plant', 'PSF History', 'Residual')
colnames(npsf_nmds.vca) <- colnames(npsf_nmds.scores)

# Calculate the total variance of each variance component#
npsf_nmds.vtot <- colSums(npsf_nmds.vca)

# Weight each variance component by the amount of variation each axis explains #
npsf_nmds.wvca <- matrix(nrow = nrow(npsf_nmds.vca), ncol = length(npsf_nmds.axisr))
for(i in 1:length(npsf_nmds.axisr)){
  for(j in 1:nrow(npsf_nmds.vca)){
    npsf_nmds.wvca[j,i] <- npsf_nmds.vca[j,i]*npsf_nmds.axisr[i] 
  }
}
# Take the total variance explained by each predictor and take the sum of those values #
rownames(npsf_nmds.wvca) <- rownames(npsf_nmds.vca); colnames(npsf_nmds.wvca) <- colnames(npsf_nmds.vca)
npsf_nmds.tvca <- rowSums(npsf_nmds.wvca)
npsf_nmds.tvca <- as.data.frame(npsf_nmds.tvca)
npsf_nmds.ptot <- colSums(npsf_nmds.tvca)

# Take the variance explained by each predictor and divide by the total variance explained and multiply by 100% #
npsf_nmds.pvca <- matrix(nrow = nrow(npsf_nmds.tvca), ncol = 1)
for(i in 1:nrow(npsf_nmds.tvca)){
  npsf_nmds.pvca[i,1] <- npsf_nmds.tvca[i,1] / npsf_nmds.ptot * 100
}

# Save the variation explained percentages of each predictor/predictor interaction #
rownames(npsf_nmds.pvca) <- rownames(npsf_nmds.vca); colnames(npsf_nmds.pvca) <- 'Variance Explained'
npsf_nmds.pvca

# Perform PermANOVA using all samples #
npsf.adon <- adonis2(npsf.wuni~Comps*Plants, npsf$met, permutations = 999)
npsf.adon_by <- adonis2(npsf.wuni~Comps*Plants, npsf$met, permutations = 999, by = 'terms')
npsf.adon
npsf.adon_by

# Perform PermDISP using all samples # 
npsf.bdis <- betadisper(npsf.wuni, group = npsf$met$PC)
anova(npsf.bdis)
TukeyHSD(npsf.bdis)

# Make an object with the data to be plotted #
npsf_nmds.load <- cbind(npsf$met, npsf_nmds.scores)
npsf_nmds.load$Plants <- gsub("T. repens", "M. truncatula", npsf_nmds.load$Plants)
npsf_nmds.load$Plants <- factor(npsf_nmds.load$Plants, levels = c("S. helvola", "C. fasciculata", "D. illinoense", "A. bracteata", "M. truncatula"))

# plot the NMDS ordination of all samples that will be used to make a patchwork plot #
npsf_nmds.plot <- ggplot(npsf_nmds.load, aes(NMDS1, NMDS2, color = Plants, shape = Comps)) +
  geom_point(size = 8) +
  scale_color_manual(name = "Community Source", labels = c(expression(italic('S. helvola')), expression(italic('C. fasciculata')), expression(italic('D. illinoense')), expression(italic('A. bracteata')), expression(italic('M. truncatula'))), values = c("#A6CEE3","#1F78B4","#FDBF6F", "#FF7F00", "#6A3D9A")) +
  scale_shape_manual(name = "Compartment", labels = c("Bulk Soil", "Rhizosphere", "PSF Soil"), values = c(15,16)) +
  xlab(expression("NMDS1 (R"^2 ~ "= 0.6251)")) +
  ylab(expression("NMDS2 (R"^2 ~ "= 0.2766)")) +
  theme_bw() +
  theme(legend.position = 'right',
        panel.grid = element_blank(),
        legend.text = element_text(size = 20, family = "Liberation Sans"),
        legend.title = element_text(size = 24, face = 'bold', family = "Liberation Sans"),
        axis.title = element_text(size= 24, face = 'bold', family = "Liberation Sans"),
        axis.text = element_text(color = "black", size = 8, family = "Liberation Sans")) +  
  annotate('text', x = -0.05, y = -0.3,
           label = expression("(PERMANOVA) F"["7,27"] ~ "= 26.094, P = 0.001;  (PERMDISP) F"["7,27"] ~  "= 6.0835, P < 0.001; 3D Stress = 0.0083"),
           size = 8, family = 'Liberation Sans') +
  coord_cartesian(ylim = c(-0.31,0.23))+
  labs(tag = "B.")

npsf_nmds.plot

# Separate samples by plant species #
## Fuzzy Bean ##
fb_npsf.ps <- subset_samples(npsf.ps, Plant == "S. helvola" | Plant == "C. fasciculata" & Compartment == "Bulk Soil")
fb_npsf.ps <- subset_taxa(fb_npsf.ps, taxa_sums(fb_npsf.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(fb_npsf.ps, 'fb_npsf')
fb_npsf_prop.ps <- transform_sample_counts(fb_npsf.ps, function(x) x/sum(x))
set.seed(248)
fb_npsf.wuni <- phyloseq::distance(fb_npsf_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
fb_npsf.bdis <- betadisper(fb_npsf.wuni, group = fb_npsf$met$Comps)
anova(fb_npsf.bdis)
TukeyHSD(fb_npsf.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
fb_npsf.met <- filter(npsf_nmds.load, Plant == "S. helvola" | Plant == "C. fasciculata" & Compartment == "Bulk Soil")

# Perform MANOVA on all samples #
fb_npsf.man <- manova(cbind(NMDS1,NMDS2)~Comps, fb_npsf.met)
summary(fb_npsf.man)

# Subset the data for pairwise MANOVA tests #
## Bulk Soil vs. Rhizosphere ##
fb_npsf_pvsn.met <- filter(fb_npsf.met, Soil_Treatment != "Common Soil")
fb_npsf_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Comps, fb_npsf_pvsn.met)
summary(fb_npsf_pvsn.man)

## Chamaecrista ##
cc_npsf.ps <- subset_samples(npsf.ps, Plant == "C. fasciculata")
cc_npsf.ps <- subset_taxa(cc_npsf.ps, taxa_sums(cc_npsf.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(cc_npsf.ps, 'cc_npsf')
cc_npsf_prop.ps <- transform_sample_counts(cc_npsf.ps, function(x) x/sum(x))
set.seed(248)
cc_npsf.wuni <- phyloseq::distance(cc_npsf_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
cc_npsf.bdis <- betadisper(cc_npsf.wuni, group = cc_npsf$met$Comps)
anova(cc_npsf.bdis)
TukeyHSD(cc_npsf.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
cc_npsf.met <- filter(npsf_nmds.load, Plant == "C. fasciculata")

# Perform MANOVA on all samples #
cc_npsf.man <- manova(cbind(NMDS1,NMDS2)~Comps, cc_npsf.met)
summary(cc_npsf.man)

# Subset the data for pairwise MANOVA tests #
## Bulk Soil vs. Rhizosphere ##
cc_npsf_pvsn.met <- filter(cc_npsf.met, Soil_Treatment != "Common Soil")
cc_npsf_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Comps, cc_npsf_pvsn.met)
summary(cc_npsf_pvsn.man)

## Desmodium ##
ds_npsf.ps <- subset_samples(npsf.ps, Plant == "D. illinoense")
ds_npsf.ps <- subset_taxa(ds_npsf.ps, taxa_sums(ds_npsf.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(ds_npsf.ps, 'ds_npsf')
ds_npsf_prop.ps <- transform_sample_counts(ds_npsf.ps, function(x) x/sum(x))
set.seed(248)
ds_npsf.wuni <- phyloseq::distance(ds_npsf_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
ds_npsf.bdis <- betadisper(ds_npsf.wuni, group = ds_npsf$met$Comps)
anova(ds_npsf.bdis)
TukeyHSD(ds_npsf.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
ds_npsf.met <- filter(npsf_nmds.load, Plant == "D. illinoense")

# Perform MANOVA on all samples #
ds_npsf.man <- manova(cbind(NMDS1,NMDS2)~Comps, ds_npsf.met)
summary(ds_npsf.man)

# Subset the data for pairwise MANOVA tests #
## Bulk Soil vs. Rhizosphere ##
ds_npsf_pvsn.met <- filter(ds_npsf.met, Soil_Treatment != "Common Soil")
ds_npsf_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Comps, ds_npsf_pvsn.met)
summary(ds_npsf_pvsn.man)

## Hog Peanut ##
hp_npsf.ps <- subset_samples(npsf.ps, Plant == "A. bracteata" | Plant == "D. illinoense" & Compartment == "Bulk Soil")
hp_npsf.ps <- subset_taxa(hp_npsf.ps, taxa_sums(hp_npsf.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(hp_npsf.ps, 'hp_npsf')
hp_npsf_prop.ps <- transform_sample_counts(hp_npsf.ps, function(x) x/sum(x))
set.seed(248)
hp_npsf.wuni <- phyloseq::distance(hp_npsf_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
hp_npsf.bdis <- betadisper(hp_npsf.wuni, group = hp_npsf$met$Comps)
anova(hp_npsf.bdis)
TukeyHSD(hp_npsf.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
hp_npsf.met <- filter(npsf_nmds.load, Plant == "A. bracteata" | Plant == "D. illinoense" & Compartment == "Bulk Soil")

# Perform MANOVA on all samples #
hp_npsf.man <- manova(cbind(NMDS1,NMDS2)~Comps, hp_npsf.met)
summary(hp_npsf.man)

# Subset the data for pairwise MANOVA tests #
## Bulk Soil vs. Rhizosphere ##
hp_npsf_pvsn.met <- filter(hp_npsf.met, Soil_Treatment != "Common Soil")
hp_npsf_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Comps, hp_npsf_pvsn.met)
summary(hp_npsf_pvsn.man)

## Medicago ##
md_npsf.ps <- subset_samples(npsf.ps, Plant == "M. truncatula" | Plant == "T. repens" & Compartment == "Bulk Soil")
md_npsf.ps <- subset_taxa(md_npsf.ps, taxa_sums(md_npsf.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(md_npsf.ps, 'md_npsf')
md_npsf_prop.ps <- transform_sample_counts(md_npsf.ps, function(x) x/sum(x))
set.seed(248)
md_npsf.wuni <- phyloseq::distance(md_npsf_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
md_npsf.bdis <- betadisper(md_npsf.wuni, group = md_npsf$met$Comps)
anova(md_npsf.bdis)
TukeyHSD(md_npsf.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
md_npsf.met <- filter(npsf_nmds.load, Plants == "M. truncatula")

# Perform MANOVA on all samples #
md_npsf.man <- manova(cbind(NMDS1,NMDS2)~Comps, md_npsf.met)
summary(md_npsf.man)

# Subset the data for pairwise MANOVA tests #
## Bulk Soil vs. Rhizosphere ##
md_npsf_pvsn.met <- filter(md_npsf.met, Soil_Treatment != "Common Soil")
md_npsf_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Comps, md_npsf_pvsn.met)
summary(md_npsf_pvsn.man)

# Bulk Soil vs Rhizosphere PSF Soil #
# Construct a phyloseq object containing only samples in the Non-PSF Soil #
wpsf.ps <- subset_samples(soil.ps, Soil_Treatment == "PSF Soil")

# Construct a Weighted Unifrac distance matrix and PCoA ordination #
wpsf_prop.ps <- transform_sample_counts(wpsf.ps, function(x) x/sum(x))

set.seed(248)
wpsf.wuni <- phyloseq::distance(wpsf_prop.ps, method = 'wunifrac')
wpsf.pcoa <- phyloseq::ordinate(wpsf_prop.ps, 'PCoA', distance = wpsf.wuni)

# Perform an NMDS analysis using the weighted Unifrac distance matrix, with the PCoA ordination as the starting ordination # 
wpsf.nmds <- metaMDS(wpsf.wuni, 
                     k = 5, try = 100, trymax = 1000, maxit = 999,
                     model = 'global', 
                     autotransform = FALSE, previous.best = wpsf.pcoa$vectors[,1:5])

# Save the loading scores for all axes and make a distance matrix from these scores #
wpsf_nmds.scores <- scores(wpsf.nmds, display = 'sites')
wpsf_nmds.dist <- dist(wpsf_nmds.scores)

# Fit a linear model using the vectorized weighted Unifrac ditsance matrix as a predictor of the vectorized loading score distance matrix to fine total R^2 of the model # 
wpsf_nmds.ffit <- lm(as.vector(wpsf_nmds.dist) ~ as.vector(wpsf.wuni))
summary(wpsf_nmds.ffit)
wpsf_nmds.totr2 <- summary(wpsf_nmds.ffit)$r.squared

# Axes Variance Calculation #
# Fit linear models as before expect to predict the distance matrix of each individual axis #
wpsf_nmds.dist1 <- dist(wpsf_nmds.scores[,1])
wpsf_nmds.fit1 <- lm(as.vector(wpsf_nmds.dist1)~as.vector(wpsf.wuni))
wpsf_nmds.r1 <- summary(wpsf_nmds.fit1)$r.squared

wpsf_nmds.dist2 <- dist(wpsf_nmds.scores[,2])
wpsf_nmds.fit2 <- lm(as.vector(wpsf_nmds.dist2)~as.vector(wpsf.wuni))
wpsf_nmds.r2 <- summary(wpsf_nmds.fit2)$r.squared

wpsf_nmds.dist3 <- dist(wpsf_nmds.scores[,3])
wpsf_nmds.fit3 <- lm(as.vector(wpsf_nmds.dist3)~as.vector(wpsf.wuni))
wpsf_nmds.r3 <- summary(wpsf_nmds.fit3)$r.squared

wpsf_nmds.dist4 <- dist(wpsf_nmds.scores[,4])
wpsf_nmds.fit4 <- lm(as.vector(wpsf_nmds.dist4)~as.vector(wpsf.wuni))
wpsf_nmds.r4 <- summary(wpsf_nmds.fit4)$r.squared

wpsf_nmds.dist5 <- dist(wpsf_nmds.scores[,5])
wpsf_nmds.fit5 <- lm(as.vector(wpsf_nmds.dist5)~as.vector(wpsf.wuni))
wpsf_nmds.r5 <- summary(wpsf_nmds.fit5)$r.squared

# Take the sum the R^2 value from each axis #
wpsf_nmds.comb <- wpsf_nmds.r1 + wpsf_nmds.r2 + wpsf_nmds.r3 + wpsf_nmds.r4 + wpsf_nmds.r5

# Divide each axis R^2 by the total of all axes and then multiply by the variation explained by the whole model
wpsf_nmds.axisr <- c()
for(i in 1:ncol(wpsf_nmds.scores)){
  wpsf_nmds.axisr[i] <- (get(paste0('wpsf_nmds.r', i)) / wpsf_nmds.comb) * wpsf_nmds.totr2 
}

### Calculating Variance Components ###
# Construct a data.frame that has sample info and their loading scores #
decompose_ps(wpsf.ps, 'wpsf')
for(i in 1:nrow(wpsf$met)){
  wpsf$met$PC[i] <- paste0(substr(wpsf$met$Plant[i],1,1), substr(wpsf$met$Compartment[i],1,2))
}
wpsf$met$Plants <- factor(wpsf$met$Plant, levels = c('S. helvola', 'C. fasciculata', 'D. illinoense', 'A. bracteata', 'T. repens', 'M. truncatula'))
wpsf$met$Comps <- factor(wpsf$met$Compartment, c('Bulk Soil', 'Rhizosphere'))
wpsf$met$Soils <- factor(wpsf$met$Soil_Treatment, levels = c("PSF Soil"))

sample_data(wpsf.ps) <- wpsf$met
wpsf_nmds.load <- cbind(wpsf$met, wpsf_nmds.scores)

# Test the mixed linear model on the first NMDS axis #
wpsf_nmds.vfit1 <- lmer(NMDS1 ~ (1|Comps) + (1|Plants) +(1|Comps:Plants), data = wpsf_nmds.load, REML = TRUE)
summary(wpsf_nmds.vfit1)
wpsf_nmds.vca1 <- as.data.frame(VarCorr(wpsf_nmds.vfit1))

# Using Loop to do each NMDS axis #
wpsf_nmds.vca <- matrix(nrow = 4, ncol = ncol(wpsf_nmds.scores))
hold <- c()
for(i in 1:ncol(wpsf_nmds.scores)){
  hold <- lmer(wpsf_nmds.scores[,i] ~  (1|Comps) + (1|Plants) +(1|Comps:Plants), data = wpsf_nmds.load, REML = TRUE)
  hold <- as.data.frame(VarCorr(hold))
  wpsf_nmds.vca[1,i] <- hold[1,4]
  wpsf_nmds.vca[2,i] <- hold[2,4]
  wpsf_nmds.vca[3,i] <- hold[3,4]
  wpsf_nmds.vca[4,i] <- hold[4,4]
}

# Save the variance components to their assigned variable/variable interaction and their NMDS loading axis #
rownames(wpsf_nmds.vca) <- c('Comp x Plant', 'Plant', 'Comp', 'Residual')
colnames(wpsf_nmds.vca) <- colnames(wpsf_nmds.scores)

# Calculate the total variance of each variance component#
wpsf_nmds.vtot <- colSums(wpsf_nmds.vca)

# Weight each variance component by the amount of variation each axis explains #
wpsf_nmds.wvca <- matrix(nrow = nrow(wpsf_nmds.vca), ncol = length(wpsf_nmds.axisr))
for(i in 1:length(wpsf_nmds.axisr)){
  for(j in 1:nrow(wpsf_nmds.vca)){
    wpsf_nmds.wvca[j,i] <- wpsf_nmds.vca[j,i]*wpsf_nmds.axisr[i] 
  }
}
# Take the total variance explained by each predictor and take the sum of those values #
rownames(wpsf_nmds.wvca) <- rownames(wpsf_nmds.vca); colnames(wpsf_nmds.wvca) <- colnames(wpsf_nmds.vca)
wpsf_nmds.tvca <- rowSums(wpsf_nmds.wvca)
wpsf_nmds.tvca <- as.data.frame(wpsf_nmds.tvca)
wpsf_nmds.ptot <- colSums(wpsf_nmds.tvca)

# Take the variance explained by each predictor and divide by the total variance explained and multiply by 100% #
wpsf_nmds.pvca <- matrix(nrow = nrow(wpsf_nmds.tvca), ncol = 1)
for(i in 1:nrow(wpsf_nmds.tvca)){
  wpsf_nmds.pvca[i,1] <- wpsf_nmds.tvca[i,1] / wpsf_nmds.ptot * 100
}

# Save the variation explained percentages of each predictor/predictor interaction #
rownames(wpsf_nmds.pvca) <- rownames(wpsf_nmds.vca); colnames(wpsf_nmds.pvca) <- 'Variance Explained'
wpsf_nmds.pvca

# Perform PermANOVA using all samples #
wpsf.adon <- adonis2(wpsf.wuni~Comps*Plants, wpsf$met, permutations = 999)
wpsf.adon_by <- adonis2(wpsf.wuni~Comps*Plants, wpsf$met, permutations = 999, by = 'terms')
wpsf.adon
wpsf.adon_by

# Perform PermDISP using all samples # 
wpsf.bdis <- betadisper(wpsf.wuni, group = wpsf$met$PC)
anova(wpsf.bdis)
TukeyHSD(wpsf.bdis)

# Make an object with the data to be plotted #
wpsf_nmds.load <- cbind(wpsf$met, wpsf_nmds.scores)

# plot the NMDS ordination of all samples that will be used to make a patchwork plot #
wpsf_nmds.plot <- ggplot(wpsf_nmds.load, aes(NMDS1, NMDS2, color = Plants, shape = Comps)) +
  geom_point(size = 8) +
  scale_color_manual(name = "Community Source", labels = c(expression(italic('S. helvola')), expression(italic('C. fasciculata')), expression(italic('D. illinoense')), expression(italic('A. bracteata')), expression(italic('T. repens')), expression(italic('M. truncatula')), "Common Soil"), values = c("#A6CEE3","#1F78B4","#FDBF6F", "#FF7F00","#CAB2D6", "#6A3D9A", "#B15928")) +
  scale_shape_manual(name = "Compartment", labels = c("Bulk Soil", "Rhizosphere", "PSF Soil"), values = c(15,16)) +
  xlab(expression("NMDS1 (R"^2 ~ "= 0.4873)")) +
  ylab(expression("NMDS2 (R"^2 ~ "= 0.3939)")) +
  theme_bw() +
  theme(legend.position = 'right',
        panel.grid = element_blank(),
        legend.text = element_text(size = 20, family = "Liberation Sans"),
        legend.title = element_text(size = 24, face = 'bold', family = "Liberation Sans"),
        axis.title = element_text(size= 24, face = 'bold', family = "Liberation Sans"),
        axis.text = element_text(color = "black", size = 8, family = "Liberation Sans")) +  
  annotate('text', x = -0.025, y = -0.3,
           label = expression("(PERMANOVA) F"["11,45"] ~ "= 31.072, P = 0.001;  (PERMDISP) F"["11,45"] ~  "= 9.1545, P < 0.001; 3D Stress = 0.0143"),
           size = 8, family = 'Liberation Sans') +
  coord_cartesian(ylim = c(-0.31,0.23)) +
  labs(tag = "A.")

wpsf_nmds.plot

# Separate samples by plant species #
## Fuzzy Bean ##
fb_wpsf.ps <- subset_samples(wpsf.ps, Plant == "S. helvola")
fb_wpsf.ps <- subset_taxa(fb_wpsf.ps, taxa_sums(fb_wpsf.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(fb_wpsf.ps, 'fb_wpsf')
fb_wpsf_prop.ps <- transform_sample_counts(fb_wpsf.ps, function(x) x/sum(x))
set.seed(248)
fb_wpsf.wuni <- phyloseq::distance(fb_wpsf_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
fb_wpsf.bdis <- betadisper(fb_wpsf.wuni, group = fb_wpsf$met$Comps)
anova(fb_wpsf.bdis)
TukeyHSD(fb_wpsf.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
fb_wpsf.met <- filter(wpsf_nmds.load, Plant == "S. helvola")

# Perform MANOVA on all samples #
fb_wpsf.man <- manova(cbind(NMDS1,NMDS2)~Comps, fb_wpsf.met)
summary(fb_wpsf.man)

# Subset the data for pairwise MANOVA tests #
## Bulk Soil vs. Rhizosphere ##
fb_wpsf_pvsn.met <- filter(fb_wpsf.met, Soil_Treatment != "Common Soil")
fb_wpsf_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Comps, fb_wpsf_pvsn.met)
summary(fb_wpsf_pvsn.man)

## Chamaecrista ##
cc_wpsf.ps <- subset_samples(wpsf.ps, Plant == "C. fasciculata")
cc_wpsf.ps <- subset_taxa(cc_wpsf.ps, taxa_sums(cc_wpsf.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(cc_wpsf.ps, 'cc_wpsf')
cc_wpsf_prop.ps <- transform_sample_counts(cc_wpsf.ps, function(x) x/sum(x))
set.seed(248)
cc_wpsf.wuni <- phyloseq::distance(cc_wpsf_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
cc_wpsf.bdis <- betadisper(cc_wpsf.wuni, group = cc_wpsf$met$Comps)
anova(cc_wpsf.bdis)
TukeyHSD(cc_wpsf.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
cc_wpsf.met <- filter(wpsf_nmds.load, Plant == "C. fasciculata")

# Perform MANOVA on all samples #
cc_wpsf.man <- manova(cbind(NMDS1,NMDS2)~Comps, cc_wpsf.met)
summary(cc_wpsf.man)

# Subset the data for pairwise MANOVA tests #
## Bulk Soil vs. Rhizosphere ##
cc_wpsf_pvsn.met <- filter(cc_wpsf.met, Soil_Treatment != "Common Soil")
cc_wpsf_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Comps, cc_wpsf_pvsn.met)
summary(cc_wpsf_pvsn.man)


## Desmodium ##
ds_wpsf.ps <- subset_samples(wpsf.ps, Plant == "D. illinoense")
ds_wpsf.ps <- subset_taxa(ds_wpsf.ps, taxa_sums(ds_wpsf.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(ds_wpsf.ps, 'ds_wpsf')
ds_wpsf_prop.ps <- transform_sample_counts(ds_wpsf.ps, function(x) x/sum(x))
set.seed(248)
ds_wpsf.wuni <- phyloseq::distance(ds_wpsf_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
ds_wpsf.bdis <- betadisper(ds_wpsf.wuni, group = ds_wpsf$met$Comps)
anova(ds_wpsf.bdis)
TukeyHSD(ds_wpsf.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
ds_wpsf.met <- filter(wpsf_nmds.load, Plant == "D. illinoense")

# Perform MANOVA on all samples #
ds_wpsf.man <- manova(cbind(NMDS1,NMDS2)~Comps, ds_wpsf.met)
summary(ds_wpsf.man)

# Subset the data for pairwise MANOVA tests #
## Bulk Soil vs. Rhizosphere ##
ds_wpsf_pvsn.met <- filter(ds_wpsf.met, Soil_Treatment != "Common Soil")
ds_wpsf_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Comps, ds_wpsf_pvsn.met)
summary(ds_wpsf_pvsn.man)


## Hog Peanut ##
hp_wpsf.ps <- subset_samples(wpsf.ps, Plant == "A. bracteata")
hp_wpsf.ps <- subset_taxa(hp_wpsf.ps, taxa_sums(hp_wpsf.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(hp_wpsf.ps, 'hp_wpsf')
hp_wpsf_prop.ps <- transform_sample_counts(hp_wpsf.ps, function(x) x/sum(x))
set.seed(248)
hp_wpsf.wuni <- phyloseq::distance(hp_wpsf_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
hp_wpsf.bdis <- betadisper(hp_wpsf.wuni, group = hp_wpsf$met$Comps)
anova(hp_wpsf.bdis)
TukeyHSD(hp_wpsf.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
hp_wpsf.met <- filter(wpsf_nmds.load, Plant == "A. bracteata")

# Perform MANOVA on all samples #
hp_wpsf.man <- manova(cbind(NMDS1,NMDS2)~Comps, hp_wpsf.met)
summary(hp_wpsf.man)

# Subset the data for pairwise MANOVA tests #
## Bulk Soil vs. Rhizosphere ##
hp_wpsf_pvsn.met <- filter(hp_wpsf.met, Soil_Treatment != "Common Soil")
hp_wpsf_pvsn.man <- manova(cbind(NMDS1,NMDS2)~Comps, hp_wpsf_pvsn.met)
summary(hp_wpsf_pvsn.man)

## Clover ##
cl_wpsf.ps <- subset_samples(wpsf.ps, Plant == "T. repens")
cl_wpsf.ps <- subset_taxa(cl_wpsf.ps, taxa_sums(cl_wpsf.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(cl_wpsf.ps, 'cl_wpsf')
cl_wpsf_prop.ps <- transform_sample_counts(cl_wpsf.ps, function(x) x/sum(x))
set.seed(248)
cl_wpsf.wuni <- phyloseq::distance(cl_wpsf_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
cl_wpsf.bdis <- betadisper(cl_wpsf.wuni, group = cl_wpsf$met$Comps)
anova(cl_wpsf.bdis)
TukeyHSD(cl_wpsf.bdis)

# Subset the NMDS plot metadata and scores to only include samples of the specified compartment and plant species #
cl_wpsf.met <- filter(wpsf_nmds.load, Plant == "T. repens")

# Perform MANOVA on all samples #
cl_wpsf.man <- manova(cbind(NMDS1,NMDS2)~Comps, cl_wpsf.met)
summary(cl_wpsf.man)

## Medicago ##
md_wpsf.ps <- subset_samples(wpsf.ps, Plant == "M. truncatula")
md_wpsf.ps <- subset_taxa(md_wpsf.ps, taxa_sums(md_wpsf.ps) > 0)

# Construct a distance matrix based on weighted unifrac distances for the samples of the given plant species and compartment #
decompose_ps(md_wpsf.ps, 'md_wpsf')
md_wpsf_prop.ps <- transform_sample_counts(md_wpsf.ps, function(x) x/sum(x))
set.seed(248)
md_wpsf.wuni <- phyloseq::distance(md_wpsf_prop.ps, method = 'wunifrac')

# Perform a PermDISP with TukeyHSD to test the homogeneity of the variance #
md_wpsf.bdis <- betadisper(md_wpsf.wuni, group = md_wpsf$met$Comps)
anova(md_wpsf.bdis)
TukeyHSD(md_wpsf.bdis)

# Subset the NMDS plot metadata and scores to only inmdude samples of the specified compartment and plant species #
md_wpsf.met <- filter(wpsf_nmds.load, Plant == "M. truncatula")

# Perform MANOVA on all samples #
md_wpsf.man <- manova(cbind(NMDS1,NMDS2)~Comps, md_wpsf.met)
summary(md_wpsf.man)

# PSF vs. Non-PSF Soil Plot #
(wpsf_nmds.plot) / 
(npsf_nmds.plot) +
  plot_layout(guides = "collect") &
  theme(plot.tag = element_text(size = 22, face = "bold", family = "Liberation Sans"))

#### Non-Nodule Stacked Histograms ####
library(microbiome); packageVersion('microbiome')
library(microbiomeutilities); packageVersion("microbiomeutilities")

# Create a phyloseq object that contains all non-nodule reads #
decompose_ps(soil.ps, 'soil_temp')
for(i in 1:nrow(soil_temp$met)){
  soil_temp$met$PS[i] <- paste0(substr(soil_temp$met$Plant[i],1,1), substr(soil_temp$met$Soil_Treatment[i], 1,1))
}
soil_temp$met$Plants <- factor(soil_temp$met$Plant, levels = c("S. helvola", "C. fasciculata", "D. illinoense", "A. bracteata", "T. repens", "M. truncatula"))
soil_temp$met$Comps <- factor(soil_temp$met$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Root Endosphere"))
soil_temp$met$Soils <- factor(soil_temp$met$Soil_Treatment, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

decompose_ps(root.ps, 'root_temp')
root_temp$met$Comps <- factor(root_temp$met$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Root Endosphere"))
colnames(root_temp$met) <- colnames(soil_temp$met)

soil_temp.ps <- phyloseq(otu_table(soil_temp$otu, taxa_are_rows = TRUE),
                         tax_table(as.matrix(soil_temp$tax[,c("Family", "Genus")])),
                         sample_data(soil_temp$met))
soil_temp.ps <- tax_glom(soil_temp.ps, taxrank = "Genus")
taxa_names(soil_temp.ps) <- tax_table(soil_temp.ps)[,"Genus"]
root_temp.ps <- phyloseq(otu_table(root_temp$otu, taxa_are_rows = TRUE),
                         tax_table(as.matrix(root_temp$tax[,c("Family", "Genus")])),
                         sample_data(root_temp$met))
root_temp.ps <- tax_glom(root_temp.ps, taxrank = "Genus")
taxa_names(root_temp.ps) <- tax_table(root_temp.ps)[,"Genus"]

all.ps <- merge_phyloseq(soil_temp.ps, root_temp.ps)

# Create a palette for each genus represented in the phyloseq object #
library(Polychrome); packageVersion("Polychrome")
all.colr <- createPalette(ntaxa(all.ps), c("#ff0000", "#00ff00", "#0000ff"))
all.colr <- as.data.frame(all.colr)
rownames(all.colr) <- sort(taxa_names(all.ps))
all.colr[350,] <- "#D4D4D4"
rownames(all.colr)[350] <- "Other"

# Separate by plant species # 
## Fuzzy Bean ##
fb_all.ps <- subset_samples(all.ps, Plant == "S. helvola")
fb_all.ps <- merge_phyloseq(fb_all.ps, subset_samples(all.ps, Plant == "C. fasciculata" & Soil_Treatment == "Non-PSF Soil" & Compartment == "Bulk Soil"))

fb_top.ps <- aggregate_top_taxa2(fb_all.ps, 19, "Genus")
fb_top.name <- names(sort(taxa_sums(fb_top.ps), decreasing = TRUE))
fb_top.colr <- all.colr[fb_top.name,]

fb_top.df <- psmelt(fb_top.ps)
fb_top.df$Genera <- factor(fb_top.df$Genus, levels = fb_top.name)
fb_top.df$Group <- factor(fb_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

fb_top.plot <- ggplot(fb_top.df, aes(x = Group, y = Abundance, fill = Genera)) +
  geom_bar(stat='identity', position = 'fill') +
  facet_wrap(~Comps, scales = "fixed") +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = fb_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('S. helvola')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18, angle = -45, hjust = 0),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90)) +
  labs(tags = "A.")
fb_top.plot

# Joining histograms #
(fb_top.plot) /
(fb_soil_nod.plot | fb_root_nod.plot) +
  plot_layout(guides = "keep") &
  theme(plot.tag = element_text(size = 22, face = "bold", family = "Liberation Sans"))

## Chamaecrista ##
cc_all.ps <- subset_samples(all.ps, Plant == "C. fasciculata")
cc_all.ps <- merge_phyloseq(cc_all.ps, subset_samples(all.ps, Plant == "S. helvola" & Soil_Treatment == "Common Soil" & Compartment == "Bulk Soil"))

cc_top.ps <- aggregate_top_taxa2(cc_all.ps, 19, "Genus")
cc_top.name <- names(sort(taxa_sums(cc_top.ps), decreasing = TRUE))
cc_top.colr <- all.colr[cc_top.name,]

cc_top.df <- psmelt(cc_top.ps)
cc_top.df$Genera <- factor(cc_top.df$Genus, levels = cc_top.name)
cc_top.df$Group <- factor(cc_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

cc_top.plot <- ggplot(cc_top.df, aes(x = Group, y = Abundance, fill = Genera)) +
  geom_bar(stat='identity', position = 'fill') +
  facet_wrap(~Comps, scales = "fixed") +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = cc_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('C. fasciculata')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18, angle = -45, hjust = 0),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90)) +
  labs(tags = "A.")
cc_top.plot

# Joining histograms #
(cc_top.plot) /
  (cc_soil_nod.plot | cc_root_nod.plot) +
  plot_layout(guides = "keep") &
  theme(plot.tag = element_text(size = 22, face = "bold", family = "Liberation Sans"))

## Desmodium ##
ds_all.ps <- subset_samples(all.ps, Plant == "D. illinoense")
ds_all.ps <- merge_phyloseq(ds_all.ps, subset_samples(all.ps, Plant == "S. helvola" & Soil_Treatment == "Common Soil" & Compartment == "Bulk Soil"))

ds_top.ps <- aggregate_top_taxa2(ds_all.ps, 19, "Genus")
ds_top.name <- names(sort(taxa_sums(ds_top.ps), decreasing = TRUE))
ds_top.colr <- all.colr[ds_top.name,]

ds_top.df <- psmelt(ds_top.ps)
ds_top.df$Genera <- factor(ds_top.df$Genus, levels = ds_top.name)
ds_top.df$Group <- factor(ds_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

ds_top.plot <- ggplot(ds_top.df, aes(x = Group, y = Abundance, fill = Genera)) +
  geom_bar(stat='identity', position = 'fill') +
  facet_wrap(~Comps, scales = "fixed") +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = ds_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('D. illinoense')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18, angle = -45, hjust = 0),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90)) +
  labs(tags = "A.")
ds_top.plot

# Joining histograms #
(ds_top.plot) /
  (ds_soil_nod.plot | ds_root_nod.plot) +
  plot_layout(guides = "keep") &
  theme(plot.tag = element_text(size = 22, face = "bold", family = "Liberation Sans"))

## Hog Peanut ##
hp_all.ps <- subset_samples(all.ps, Plant == "A. bracteata")
hp_all.ps <- merge_phyloseq(hp_all.ps, subset_samples(all.ps, Plant == "S. helvola" & Soil_Treatment == "Common Soil" & Compartment == "Bulk Soil"))
hp_all.ps <- merge_phyloseq(hp_all.ps, subset_samples(all.ps, Plant == "D. illinoense" & Soil_Treatment == "Non-PSF Soil" & Compartment == "Bulk Soil"))


hp_top.ps <- aggregate_top_taxa2(hp_all.ps, 19, "Genus")
hp_top.name <- names(sort(taxa_sums(hp_top.ps), decreasing = TRUE))
hp_top.colr <- all.colr[hp_top.name,]

hp_top.df <- psmelt(hp_top.ps)
hp_top.df$Genera <- factor(hp_top.df$Genus, levels = hp_top.name)
hp_top.df$Group <- factor(hp_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

hp_top.plot <- ggplot(hp_top.df, aes(x = Group, y = Abundance, fill = Genera)) +
  geom_bar(stat='identity', position = 'fill') +
  facet_wrap(~Comps, scales = "fixed") +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = hp_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('A. bracteata')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18, angle = -45, hjust = 0),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90)) +
  labs(tags = "A.")
hp_top.plot

# Joining histograms #
(hp_top.plot) /
  (hp_soil_nod.plot | hp_root_nod.plot) +
  plot_layout(guides = "keep") &
  theme(plot.tag = element_text(size = 22, face = "bold", family = "Liberation Sans"))

## Clover ##
cl_all.ps <- subset_samples(all.ps, Plant == "T. repens")
cl_all.ps <- merge_phyloseq(cl_all.ps, subset_samples(all.ps, Plant == "S. helvola" & Soil_Treatment == "Common Soil" & Compartment == "Bulk Soil"))


cl_top.ps <- aggregate_top_taxa2(cl_all.ps, 19, "Genus")
cl_top.name <- names(sort(taxa_sums(cl_top.ps), decreasing = TRUE))
cl_top.name <- c("Other", cl_top.name[-2])
cl_top.colr <- all.colr[cl_top.name,]

cl_top.df <- psmelt(cl_top.ps)
cl_top.df$Genera <- factor(cl_top.df$Genus, levels = cl_top.name)
cl_top.df$Group <- factor(cl_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

cl_top.plot <- ggplot(cl_top.df, aes(x = Group, y = Abundance, fill = Genera)) +
  geom_bar(stat='identity', position = 'fill') +
  facet_wrap(~Comps, scales = "fixed") +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = cl_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('T. repens')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18, angle = -45, hjust = 0),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90)) +
  labs(tags = "A.")
cl_top.plot

# Joining histograms #
(cl_top.plot) /
  (cl_soil_nod.plot | cl_root_nod.plot) +
  plot_layout(guides = "keep") &
  theme(plot.tag = element_text(size = 22, face = "bold", family = "Liberation Sans"))

## Medicago ##
md_all.ps <- subset_samples(all.ps, Plant == "M. truncatula")
md_all.ps <- merge_phyloseq(md_all.ps, subset_samples(all.ps, Plant == "S. helvola" & Soil_Treatment == "Common Soil" & Compartment == "Bulk Soil"))
md_all.ps <- merge_phyloseq(md_all.ps, subset_samples(all.ps, Plant == "T. repens" & Soil_Treatment == "Non-PSF Soil" & Compartment == "Bulk Soil"))


md_top.ps <- aggregate_top_taxa2(md_all.ps, 19, "Genus")
md_top.name <- names(sort(taxa_sums(md_top.ps), decreasing = TRUE))
md_top.colr <- all.colr[md_top.name,]

md_top.df <- psmelt(md_top.ps)
md_top.df$Genera <- factor(md_top.df$Genus, levels = md_top.name)
md_top.df$Group <- factor(md_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

md_top.plot <- ggplot(md_top.df, aes(x = Group, y = Abundance, fill = Genera)) +
  geom_bar(stat='identity', position = 'fill') +
  facet_wrap(~Comps, scales = "fixed") +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = md_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('M. truncatula')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18, angle = -45, hjust = 0),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90)) +
  labs(tags = "A.")
md_top.plot

# Joining histograms #
(md_top.plot) /
  (md_soil_nod.plot | md_root_nod.plot) +
  plot_layout(guides = "keep") &
  theme(plot.tag = element_text(size = 22, face = "bold", family = "Liberation Sans"))

#### Differential Abundance ####
library(Maaslin2); packageVersion("Maaslin2")

# Make a decomposed phyloseq object and at a variable for each unique tripartite group of the soil and root samples #
decompose_ps(soil.ps, "soil_maas")
for(i in 1:nrow(soil_maas$met)){
  soil_maas$met$PSC[i] <- paste0(substr(soil_maas$met$Plant[i], 1, 1), substr(soil_maas$met$Soil_Treatment[i], 1, 1), substr(soil_maas$met$Compartment[i],1,2)) 
}
decompose_ps(root.ps, "root_maas")
for(i in 1:nrow(root_maas$met)){
  root_maas$met$PSC[i] <- paste0(substr(root_maas$met$Plant.Species[i], 1, 1), substr(root_maas$met$Soil.Origin[i], 1, 1), substr(root_maas$met$Compartment[i],1,2)) 
}

# Make a directory for the differential abundance plots #
if(!dir.exists('./maas_results')){
  system('mkdir ./maas_results')
}

### Perform Maaslin2 analyses on a per plant basis ###
## Fuzzy Bean ##
# Organize ASVs by number #
fb_bulk.sort <- as.numeric(sub("ASV([0-9]+).*", "\\1", rownames(fb_bulk$otu)))
fb_bulk$otu <- fb_bulk$otu[order(fb_bulk.sort),]
# Maaslin2 analysis for the rhiz soil data (PSF Reference) # 
fb_bulk_wpsf.maas <- Maaslin2(input_data = fb_bulk$otu,
                               input_metadata = fb_bulk$met,
                               output = "./maas_results/fb_bulk_wpsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within bulk soils of Fuzzy Bean #
fb_bulk_wpsf.res <- fb_bulk_wpsf.maas$results
fb_bulk_wpsf.dares <- c()
for(i in 1:nrow(fb_bulk_wpsf.res)){
  if(fb_bulk_wpsf.res$qval[i] < 0.05){
    fb_bulk_wpsf.dares <- rbind(fb_bulk_wpsf.dares, fb_bulk_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the bulk soil data (PSF Reference) # 
fb_bulk_npsf.maas <- Maaslin2(input_data = fb_bulk$otu,
                               input_metadata = fb_bulk$met,
                               output = "./maas_results/fb_bulk_npsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,Non-PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within bulk soils of Fuzzy Bean #
fb_bulk_npsf.res <- fb_bulk_npsf.maas$results
fb_bulk_npsf.dares <- c()
for(i in 1:nrow(fb_bulk_npsf.res)){
  if(fb_bulk_npsf.res$qval[i] < 0.05){
    fb_bulk_npsf.dares <- rbind(fb_bulk_npsf.dares, fb_bulk_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Organize ASVs by number #
fb_rhiz.sort <- as.numeric(sub("ASV([0-9]+).*", "\\1", rownames(fb_rhiz$otu)))
fb_rhiz$otu <- fb_rhiz$otu[order(fb_rhiz.sort),]

# Maaslin2 analysis for the bulk soil data (PSF Reference) # 
fb_rhiz_wpsf.maas <- Maaslin2(input_data = fb_rhiz$otu,
                               input_metadata = fb_rhiz$met,
                               output = "./maas_results/fb_rhiz_wpsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within rhiz soils of Fuzzy Bean #
fb_rhiz_wpsf.res <- fb_rhiz_wpsf.maas$results
fb_rhiz_wpsf.dares <- c()
for(i in 1:nrow(fb_rhiz_wpsf.res)){
  if(fb_rhiz_wpsf.res$qval[i] < 0.05){
    fb_rhiz_wpsf.dares <- rbind(fb_rhiz_wpsf.dares, fb_rhiz_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the rhiz soil data (PSF Reference) # 
fb_rhiz_npsf.maas <- Maaslin2(input_data = fb_rhiz$otu,
                               input_metadata = fb_rhiz$met,
                               output = "./maas_results/fb_rhiz_npsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,Non-PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within rhiz soils of Fuzzy Bean #
fb_rhiz_npsf.res <- fb_rhiz_npsf.maas$results
fb_rhiz_npsf.dares <- c()
for(i in 1:nrow(fb_rhiz_npsf.res)){
  if(fb_rhiz_npsf.res$qval[i] < 0.05){
    fb_rhiz_npsf.dares <- rbind(fb_rhiz_npsf.dares, fb_rhiz_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}
# Maaslin2 analysis for the root endosphere data (PSF Soil Reference)# 
fb_root_wpsf.maas <- Maaslin2(input_data = root_maas$otu,
                              input_metadata = root_maas$met,
                              output = "./maas_results/fb_root_wpsf.maas",
                              fixed_effects = c("PSC"),
                              analysis_method = "LM",
                              normalization = "CSS",
                              transform = "NONE",
                              min_prevalence = 0.05,
                              correction = "BH",
                              max_significance = 0.05,
                              reference = c("PSC,SPRo"), 
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE,
                              save_scatter = FALSE)

# Save only the pairwise comparisons within root endosphere of Fuzzy Bean #
fb_root_wpsf.res <- fb_root_wpsf.maas$results
fb_root_wpsf.dares <- c()
for(i in 1:nrow(fb_root_wpsf.res)){
  if(fb_root_wpsf.res$value[i] == "SNRo" | fb_root_wpsf.res$value[i] == "SCRo"){
    if(fb_root_wpsf.res$qval[i] < 0.05){
      fb_root_wpsf.dares <- rbind(fb_root_wpsf.dares, fb_root_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the root endosphere data (Non-PSF Soil Reference)# 
fb_root_npsf.maas <- Maaslin2(input_data = root_maas$otu,
                              input_metadata = root_maas$met,
                              output = "./maas_results/fb_root_npsf.maas",
                              fixed_effects = c("PSC"),
                              analysis_method = "LM",
                              normalization = "CSS",
                              transform = "NONE",
                              min_prevalence = 0.05,
                              correction = "BH",
                              max_significance = 0.05,
                              reference = c("PSC,SNRo"), 
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE,
                              save_scatter = FALSE)

# Save only the pairwise comparisons within root endosphere of Fuzzy Bean #
fb_root_npsf.res <- fb_root_npsf.maas$results
fb_root_npsf.dares <- c()
for(i in 1:nrow(fb_root_npsf.res)){
  if(fb_root_npsf.res$value[i] == "SPRo" | fb_root_npsf.res$value[i] == "SCRo"){
    if(fb_root_npsf.res$qval[i] < 0.05){
      fb_root_npsf.dares <- rbind(fb_root_npsf.dares, fb_root_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data # 
fb_npsf.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/fb_npsf.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,SNRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)
# Save only the pairwise comparisons within Non-PSF treatments of Fuzzy Bean #
fb_npsf.res <- fb_npsf.maas$results
fb_npsf.dares <- c()
for(i in 1:nrow(fb_npsf.res)){
  if(fb_npsf.res$value[i] == "CNBu"){
    if(fb_npsf.res$qval[i] < 0.05){
      fb_npsf.dares <- rbind(fb_npsf.dares, fb_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data #
fb_comm.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/fb_comm.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,SCRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)
# Save only the pairwise comparisons within Common Soil treatments of Fuzzy Bean #
fb_comm.res <- fb_comm.maas$results
fb_comm.dares <- c()
for(i in 1:nrow(fb_comm.res)){
  if(fb_comm.res$value[i] == "SCBu"){
    if(fb_comm.res$qval[i] < 0.05){
      fb_comm.dares <- rbind(fb_comm.dares, fb_comm.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data #
fb_wpsf.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/fb_wpsf.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,SPRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)

# Save only the pairwise comparisons within PSF Soil treatments of Fuzzy Bean #
fb_wpsf.res <- fb_wpsf.maas$results
fb_wpsf.dares <- c()
for(i in 1:nrow(fb_wpsf.res)){
  if(fb_wpsf.res$value[i] == "SPBu"){
    if(fb_wpsf.res$qval[i] < 0.05){
      fb_wpsf.dares <- rbind(fb_wpsf.dares, fb_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

## Chamaecrista ##
# Organize ASVs by number #
cc_bulk.sort <- as.numeric(sub("ASV([0-9]+).*", "\\1", rownames(cc_bulk$otu)))
cc_bulk$otu <- cc_bulk$otu[order(cc_bulk.sort),]
# Maaslin2 analysis for the rhiz soil data (PSF Reference) # 
cc_bulk_wpsf.maas <- Maaslin2(input_data = cc_bulk$otu,
                               input_metadata = cc_bulk$met,
                               output = "./maas_results/cc_bulk_wpsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within bulk soils of Fuzzy Bean #
cc_bulk_wpsf.res <- cc_bulk_wpsf.maas$results
cc_bulk_wpsf.dares <- c()
for(i in 1:nrow(cc_bulk_wpsf.res)){
  if(cc_bulk_wpsf.res$qval[i] < 0.05){
    cc_bulk_wpsf.dares <- rbind(cc_bulk_wpsf.dares, cc_bulk_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the bulk soil data (PSF Reference) # 
cc_bulk_npsf.maas <- Maaslin2(input_data = cc_bulk$otu,
                               input_metadata = cc_bulk$met,
                               output = "./maas_results/cc_bulk_npsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,Non-PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within bulk soils of Fuzzy Bean #
cc_bulk_npsf.res <- cc_bulk_npsf.maas$results
cc_bulk_npsf.dares <- c()
for(i in 1:nrow(cc_bulk_npsf.res)){
  if(cc_bulk_npsf.res$qval[i] < 0.05){
    cc_bulk_npsf.dares <- rbind(cc_bulk_npsf.dares, cc_bulk_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Organize ASVs by number #
cc_rhiz.sort <- as.numeric(sub("ASV([0-9]+).*", "\\1", rownames(cc_rhiz$otu)))
cc_rhiz$otu <- cc_rhiz$otu[order(cc_rhiz.sort),]

# Maaslin2 analysis for the bulk soil data (PSF Reference) # 
cc_rhiz_wpsf.maas <- Maaslin2(input_data = cc_rhiz$otu,
                               input_metadata = cc_rhiz$met,
                               output = "./maas_results/cc_rhiz_wpsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within rhiz soils of Fuzzy Bean #
cc_rhiz_wpsf.res <- cc_rhiz_wpsf.maas$results
cc_rhiz_wpsf.dares <- c()
for(i in 1:nrow(cc_rhiz_wpsf.res)){
  if(cc_rhiz_wpsf.res$qval[i] < 0.05){
    cc_rhiz_wpsf.dares <- rbind(cc_rhiz_wpsf.dares, cc_rhiz_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the rhiz soil data (PSF Reference) # 
cc_rhiz_npsf.maas <- Maaslin2(input_data = cc_rhiz$otu,
                               input_metadata = cc_rhiz$met,
                               output = "./maas_results/cc_rhiz_npsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,Non-PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within rhiz soils of Fuzzy Bean #
cc_rhiz_npsf.res <- cc_rhiz_npsf.maas$results
cc_rhiz_npsf.dares <- c()
for(i in 1:nrow(cc_rhiz_npsf.res)){
  if(cc_rhiz_npsf.res$qval[i] < 0.05){
    cc_rhiz_npsf.dares <- rbind(cc_rhiz_npsf.dares, cc_rhiz_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the root endosphere data (PSF Soil Reference)# 
cc_root_wpsf.maas <- Maaslin2(input_data = root_maas$otu,
                              input_metadata = root_maas$met,
                              output = "./maas_results/cc_root_wpsf.maas",
                              fixed_effects = c("PSC"),
                              analysis_method = "LM",
                              normalization = "CSS",
                              transform = "NONE",
                              min_prevalence = 0.05,
                              correction = "BH",
                              max_significance = 0.05,
                              reference = c("PSC,CPRo"), 
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE,
                              save_scatter = FALSE)

# Save only the pairwise comparisons within root endosphere of Fuzzy Bean #
cc_root_wpsf.res <- cc_root_wpsf.maas$results
cc_root_wpsf.dares <- c()
for(i in 1:nrow(cc_root_wpsf.res)){
  if(cc_root_wpsf.res$value[i] == "CNRo" | cc_root_wpsf.res$value[i] == "CCRo"){
    if(cc_root_wpsf.res$qval[i] < 0.05){
      cc_root_wpsf.dares <- rbind(cc_root_wpsf.dares, cc_root_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the root endosphere data (Non-PSF Soil Reference)# 
cc_root_npsf.maas <- Maaslin2(input_data = root_maas$otu,
                              input_metadata = root_maas$met,
                              output = "./maas_results/cc_root_npsf.maas",
                              fixed_effects = c("PSC"),
                              analysis_method = "LM",
                              normalization = "CSS",
                              transform = "NONE",
                              min_prevalence = 0.05,
                              correction = "BH",
                              max_significance = 0.05,
                              reference = c("PSC,CNRo"), 
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE,
                              save_scatter = FALSE)

# Save only the pairwise comparisons within root endosphere of Fuzzy Bean #
cc_root_npsf.res <- cc_root_npsf.maas$results
cc_root_npsf.dares <- c()
for(i in 1:nrow(cc_root_npsf.res)){
  if(cc_root_npsf.res$value[i] == "CPRo" | cc_root_npsf.res$value[i] == "CCRo"){
    if(cc_root_npsf.res$qval[i] < 0.05){
      cc_root_npsf.dares <- rbind(cc_root_npsf.dares, cc_root_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data # 
cc_npsf.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/cc_npsf.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,CNRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)
# Save only the pairwise comparisons within Non-PSF treatments of Fuzzy Bean #
cc_npsf.res <- cc_npsf.maas$results
cc_npsf.dares <- c()
for(i in 1:nrow(cc_npsf.res)){
  if(cc_npsf.res$value[i] == "CNBu"){
    if(cc_npsf.res$qval[i] < 0.05){
      cc_npsf.dares <- rbind(cc_npsf.dares, cc_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data #
cc_comm.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/cc_comm.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,CCRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)
# Save only the pairwise comparisons within Common Soil treatments of Fuzzy Bean #
cc_comm.res <- cc_comm.maas$results
cc_comm.dares <- c()
for(i in 1:nrow(cc_comm.res)){
  if(cc_comm.res$value[i] == "SCBu"){
    if(cc_comm.res$qval[i] < 0.05){
      cc_comm.dares <- rbind(cc_comm.dares, cc_comm.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data #
cc_wpsf.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/cc_wpsf.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,CPRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)

# Save only the pairwise comparisons within PSF Soil treatments of Fuzzy Bean #
cc_wpsf.res <- cc_wpsf.maas$results
cc_wpsf.dares <- c()
for(i in 1:nrow(cc_wpsf.res)){
  if(cc_wpsf.res$value[i] == "CPBu"){
    if(cc_wpsf.res$qval[i] < 0.05){
      cc_wpsf.dares <- rbind(cc_wpsf.dares, cc_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

## Desmodium ##
# Organize ASVs by number #
ds_bulk.sort <- as.numeric(sub("ASV([0-9]+).*", "\\1", rownames(ds_bulk$otu)))
ds_bulk$otu <- ds_bulk$otu[order(ds_bulk.sort),]
# Maaslin2 analysis for the rhiz soil data (PSF Reference) # 
ds_bulk_wpsf.maas <- Maaslin2(input_data = ds_bulk$otu,
                               input_metadata = ds_bulk$met,
                               output = "./maas_results/ds_bulk_wpsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within bulk soils of Fuzzy Bean #
ds_bulk_wpsf.res <- ds_bulk_wpsf.maas$results
ds_bulk_wpsf.dares <- c()
for(i in 1:nrow(ds_bulk_wpsf.res)){
  if(ds_bulk_wpsf.res$qval[i] < 0.05){
    ds_bulk_wpsf.dares <- rbind(ds_bulk_wpsf.dares, ds_bulk_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the bulk soil data (PSF Reference) # 
ds_bulk_npsf.maas <- Maaslin2(input_data = ds_bulk$otu,
                               input_metadata = ds_bulk$met,
                               output = "./maas_results/ds_bulk_npsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,Non-PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within bulk soils of Fuzzy Bean #
ds_bulk_npsf.res <- ds_bulk_npsf.maas$results
ds_bulk_npsf.dares <- c()
for(i in 1:nrow(ds_bulk_npsf.res)){
  if(ds_bulk_npsf.res$qval[i] < 0.05){
    ds_bulk_npsf.dares <- rbind(ds_bulk_npsf.dares, ds_bulk_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Organize ASVs by number #
ds_rhiz.sort <- as.numeric(sub("ASV([0-9]+).*", "\\1", rownames(ds_rhiz$otu)))
ds_rhiz$otu <- ds_rhiz$otu[order(ds_rhiz.sort),]

# Maaslin2 analysis for the bulk soil data (PSF Reference) # 
ds_rhiz_wpsf.maas <- Maaslin2(input_data = ds_rhiz$otu,
                               input_metadata = ds_rhiz$met,
                               output = "./maas_results/ds_rhiz_wpsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within rhiz soils of Fuzzy Bean #
ds_rhiz_wpsf.res <- ds_rhiz_wpsf.maas$results
ds_rhiz_wpsf.dares <- c()
for(i in 1:nrow(ds_rhiz_wpsf.res)){
  if(ds_rhiz_wpsf.res$qval[i] < 0.05){
    ds_rhiz_wpsf.dares <- rbind(ds_rhiz_wpsf.dares, ds_rhiz_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the rhiz soil data (PSF Reference) # 
ds_rhiz_npsf.maas <- Maaslin2(input_data = ds_rhiz$otu,
                               input_metadata = ds_rhiz$met,
                               output = "./maas_results/ds_rhiz_npsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,Non-PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within rhiz soils of Fuzzy Bean #
ds_rhiz_npsf.res <- ds_rhiz_npsf.maas$results
ds_rhiz_npsf.dares <- c()
for(i in 1:nrow(ds_rhiz_npsf.res)){
  if(ds_rhiz_npsf.res$qval[i] < 0.05){
    ds_rhiz_npsf.dares <- rbind(ds_rhiz_npsf.dares, ds_rhiz_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the root endosphere data (PSF Soil Reference)# 
ds_root_wpsf.maas <- Maaslin2(input_data = root_maas$otu,
                              input_metadata = root_maas$met,
                              output = "./maas_results/ds_root_wpsf.maas",
                              fixed_effects = c("PSC"),
                              analysis_method = "LM",
                              normalization = "CSS",
                              transform = "NONE",
                              min_prevalence = 0.05,
                              correction = "BH",
                              max_significance = 0.05,
                              reference = c("PSC,DPRo"), 
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE,
                              save_scatter = FALSE)

# Save only the pairwise comparisons within root endosphere of Fuzzy Bean #
ds_root_wpsf.res <- ds_root_wpsf.maas$results
ds_root_wpsf.dares <- c()
for(i in 1:nrow(ds_root_wpsf.res)){
  if(ds_root_wpsf.res$value[i] == "DNRo" | ds_root_wpsf.res$value[i] == "DCRo"){
    if(ds_root_wpsf.res$qval[i] < 0.05){
      ds_root_wpsf.dares <- rbind(ds_root_wpsf.dares, ds_root_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the root endosphere data (Non-PSF Soil Reference)# 
ds_root_npsf.maas <- Maaslin2(input_data = root_maas$otu,
                              input_metadata = root_maas$met,
                              output = "./maas_results/ds_root_npsf.maas",
                              fixed_effects = c("PSC"),
                              analysis_method = "LM",
                              normalization = "CSS",
                              transform = "NONE",
                              min_prevalence = 0.05,
                              correction = "BH",
                              max_significance = 0.05,
                              reference = c("PSC,DNRo"), 
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE,
                              save_scatter = FALSE)

# Save only the pairwise comparisons within root endosphere of Fuzzy Bean #
ds_root_npsf.res <- ds_root_npsf.maas$results
ds_root_npsf.dares <- c()
for(i in 1:nrow(ds_root_npsf.res)){
  if(ds_root_npsf.res$value[i] == "DPRo" | ds_root_npsf.res$value[i] == "DCRo"){
    if(ds_root_npsf.res$qval[i] < 0.05){
      ds_root_npsf.dares <- rbind(ds_root_npsf.dares, ds_root_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data # 
ds_npsf.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/ds_npsf.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,DNRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)
# Save only the pairwise comparisons within Non-PSF treatments of Fuzzy Bean #
ds_npsf.res <- ds_npsf.maas$results
ds_npsf.dares <- c()
for(i in 1:nrow(ds_npsf.res)){
  if(ds_npsf.res$value[i] == "DNBu"){
    if(ds_npsf.res$qval[i] < 0.05){
      ds_npsf.dares <- rbind(ds_npsf.dares, ds_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data #
ds_comm.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/ds_comm.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,DCRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)
# Save only the pairwise comparisons within Common Soil treatments of Fuzzy Bean #
ds_comm.res <- ds_comm.maas$results
ds_comm.dares <- c()
for(i in 1:nrow(ds_comm.res)){
  if(ds_comm.res$value[i] == "SCBu"){
    if(ds_comm.res$qval[i] < 0.05){
      ds_comm.dares <- rbind(ds_comm.dares, ds_comm.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data #
ds_wpsf.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/ds_wpsf.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,DPRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)

# Save only the pairwise comparisons within PSF Soil treatments of Fuzzy Bean #
ds_wpsf.res <- ds_wpsf.maas$results
ds_wpsf.dares <- c()
for(i in 1:nrow(ds_wpsf.res)){
  if(ds_wpsf.res$value[i] == "DPBu"){
    if(ds_wpsf.res$qval[i] < 0.05){
      ds_wpsf.dares <- rbind(ds_wpsf.dares, ds_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

## Hog Peanut ##
# Organize ASVs by number #
hp_bulk.sort <- as.numeric(sub("ASV([0-9]+).*", "\\1", rownames(hp_bulk$otu)))
hp_bulk$otu <- hp_bulk$otu[order(hp_bulk.sort),]
# Maaslin2 analysis for the rhiz soil data (PSF Reference) # 
hp_bulk_wpsf.maas <- Maaslin2(input_data = hp_bulk$otu,
                               input_metadata = hp_bulk$met,
                               output = "./maas_results/hp_bulk_wpsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within bulk soils of Fuzzy Bean #
hp_bulk_wpsf.res <- hp_bulk_wpsf.maas$results
hp_bulk_wpsf.dares <- c()
for(i in 1:nrow(hp_bulk_wpsf.res)){
  if(hp_bulk_wpsf.res$qval[i] < 0.05){
    hp_bulk_wpsf.dares <- rbind(hp_bulk_wpsf.dares, hp_bulk_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the bulk soil data (PSF Reference) # 
hp_bulk_npsf.maas <- Maaslin2(input_data = hp_bulk$otu,
                               input_metadata = hp_bulk$met,
                               output = "./maas_results/hp_bulk_npsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,Non-PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within bulk soils of Fuzzy Bean #
hp_bulk_npsf.res <- hp_bulk_npsf.maas$results
hp_bulk_npsf.dares <- c()
for(i in 1:nrow(hp_bulk_npsf.res)){
  if(hp_bulk_npsf.res$qval[i] < 0.05){
    hp_bulk_npsf.dares <- rbind(hp_bulk_npsf.dares, hp_bulk_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Organize ASVs by number #
hp_rhiz.sort <- as.numeric(sub("ASV([0-9]+).*", "\\1", rownames(hp_rhiz$otu)))
hp_rhiz$otu <- hp_rhiz$otu[order(hp_rhiz.sort),]

# Maaslin2 analysis for the bulk soil data (PSF Reference) # 
hp_rhiz_wpsf.maas <- Maaslin2(input_data = hp_rhiz$otu,
                               input_metadata = hp_rhiz$met,
                               output = "./maas_results/hp_rhiz_wpsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within rhiz soils of Fuzzy Bean #
hp_rhiz_wpsf.res <- hp_rhiz_wpsf.maas$results
hp_rhiz_wpsf.dares <- c()
for(i in 1:nrow(hp_rhiz_wpsf.res)){
  if(hp_rhiz_wpsf.res$qval[i] < 0.05){
    hp_rhiz_wpsf.dares <- rbind(hp_rhiz_wpsf.dares, hp_rhiz_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the rhiz soil data (PSF Reference) # 
hp_rhiz_npsf.maas <- Maaslin2(input_data = hp_rhiz$otu,
                               input_metadata = hp_rhiz$met,
                               output = "./maas_results/hp_rhiz_npsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,Non-PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within rhiz soils of Fuzzy Bean #
hp_rhiz_npsf.res <- hp_rhiz_npsf.maas$results
hp_rhiz_npsf.dares <- c()
for(i in 1:nrow(hp_rhiz_npsf.res)){
  if(hp_rhiz_npsf.res$qval[i] < 0.05){
    hp_rhiz_npsf.dares <- rbind(hp_rhiz_npsf.dares, hp_rhiz_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the root endosphere data (PSF Soil Reference)# 
hp_root_wpsf.maas <- Maaslin2(input_data = root_maas$otu,
                              input_metadata = root_maas$met,
                              output = "./maas_results/hp_root_wpsf.maas",
                              fixed_effects = c("PSC"),
                              analysis_method = "LM",
                              normalization = "CSS",
                              transform = "NONE",
                              min_prevalence = 0.05,
                              correction = "BH",
                              max_significance = 0.05,
                              reference = c("PSC,APRo"), 
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE,
                              save_scatter = FALSE)

# Save only the pairwise comparisons within root endosphere of Fuzzy Bean #
hp_root_wpsf.res <- hp_root_wpsf.maas$results
hp_root_wpsf.dares <- c()
for(i in 1:nrow(hp_root_wpsf.res)){
  if(hp_root_wpsf.res$value[i] == "ANRo" | hp_root_wpsf.res$value[i] == "ACRo"){
    if(hp_root_wpsf.res$qval[i] < 0.05){
      hp_root_wpsf.dares <- rbind(hp_root_wpsf.dares, hp_root_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the root endosphere data (Non-PSF Soil Reference)# 
hp_root_npsf.maas <- Maaslin2(input_data = root_maas$otu,
                              input_metadata = root_maas$met,
                              output = "./maas_results/hp_root_npsf.maas",
                              fixed_effects = c("PSC"),
                              analysis_method = "LM",
                              normalization = "CSS",
                              transform = "NONE",
                              min_prevalence = 0.05,
                              correction = "BH",
                              max_significance = 0.05,
                              reference = c("PSC,ANRo"), 
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE,
                              save_scatter = FALSE)

# Save only the pairwise comparisons within root endosphere of Fuzzy Bean #
hp_root_npsf.res <- hp_root_npsf.maas$results
hp_root_npsf.dares <- c()
for(i in 1:nrow(hp_root_npsf.res)){
  if(hp_root_npsf.res$value[i] == "APRo" | hp_root_npsf.res$value[i] == "ACRo"){
    if(hp_root_npsf.res$qval[i] < 0.05){
      hp_root_npsf.dares <- rbind(hp_root_npsf.dares, hp_root_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data # 
hp_npsf.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/hp_npsf.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,ANRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)
# Save only the pairwise comparisons within Non-PSF treatments of Fuzzy Bean #
hp_npsf.res <- hp_npsf.maas$results
hp_npsf.dares <- c()
for(i in 1:nrow(hp_npsf.res)){
  if(hp_npsf.res$value[i] == "DNBu"){
    if(hp_npsf.res$qval[i] < 0.05){
      hp_npsf.dares <- rbind(hp_npsf.dares, hp_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data #
hp_comm.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/hp_comm.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,ACRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)
# Save only the pairwise comparisons within Common Soil treatments of Fuzzy Bean #
hp_comm.res <- hp_comm.maas$results
hp_comm.dares <- c()
for(i in 1:nrow(hp_comm.res)){
  if(hp_comm.res$value[i] == "SCBu"){
    if(hp_comm.res$qval[i] < 0.05){
      hp_comm.dares <- rbind(hp_comm.dares, hp_comm.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data #
hp_wpsf.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/hp_wpsf.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,APRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)

# Save only the pairwise comparisons within PSF Soil treatments of Fuzzy Bean #
hp_wpsf.res <- hp_wpsf.maas$results
hp_wpsf.dares <- c()
for(i in 1:nrow(hp_wpsf.res)){
  if(hp_wpsf.res$value[i] == "APBu"){
    if(hp_wpsf.res$qval[i] < 0.05){
      hp_wpsf.dares <- rbind(hp_wpsf.dares, hp_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

## Clover ##
# Organize ASVs by number #
cl_bulk.sort <- as.numeric(sub("ASV([0-9]+).*", "\\1", rownames(cl_bulk$otu)))
cl_bulk$otu <- cl_bulk$otu[order(cl_bulk.sort),]
# Maaslin2 analysis for the rhiz soil data (PSF Reference) # 
cl_bulk_wpsf.maas <- Maaslin2(input_data = cl_bulk$otu,
                               input_metadata = cl_bulk$met,
                               output = "./maas_results/cl_bulk_wpsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within bulk soils of Fuzzy Bean #
cl_bulk_wpsf.res <- cl_bulk_wpsf.maas$results
cl_bulk_wpsf.dares <- c()
for(i in 1:nrow(cl_bulk_wpsf.res)){
  if(cl_bulk_wpsf.res$qval[i] < 0.05){
    cl_bulk_wpsf.dares <- rbind(cl_bulk_wpsf.dares, cl_bulk_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the bulk soil data (PSF Reference) # 
cl_bulk_npsf.maas <- Maaslin2(input_data = cl_bulk$otu,
                               input_metadata = cl_bulk$met,
                               output = "./maas_results/cl_bulk_npsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,Non-PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within bulk soils of Fuzzy Bean #
cl_bulk_npsf.res <- cl_bulk_npsf.maas$results
cl_bulk_npsf.dares <- c()
for(i in 1:nrow(cl_bulk_npsf.res)){
  if(cl_bulk_npsf.res$qval[i] < 0.05){
    cl_bulk_npsf.dares <- rbind(cl_bulk_npsf.dares, cl_bulk_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Organize ASVs by number #
cl_rhiz.sort <- as.numeric(sub("ASV([0-9]+).*", "\\1", rownames(cl_rhiz$otu)))
cl_rhiz$otu <- cl_rhiz$otu[order(cl_rhiz.sort),]

# Maaslin2 analysis for the bulk soil data (PSF Reference) # 
cl_rhiz_wpsf.maas <- Maaslin2(input_data = cl_rhiz$otu,
                               input_metadata = cl_rhiz$met,
                               output = "./maas_results/cl_rhiz_wpsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within rhiz soils of Fuzzy Bean #
cl_rhiz_wpsf.res <- cl_rhiz_wpsf.maas$results
cl_rhiz_wpsf.dares <- c()
for(i in 1:nrow(cl_rhiz_wpsf.res)){
  if(cl_rhiz_wpsf.res$qval[i] < 0.05){
    cl_rhiz_wpsf.dares <- rbind(cl_rhiz_wpsf.dares, cl_rhiz_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the root endosphere data (PSF Soil Reference)# 
cl_root_wpsf.maas <- Maaslin2(input_data = root_maas$otu,
                              input_metadata = root_maas$met,
                              output = "./maas_results/cl_root_wpsf.maas",
                              fixed_effects = c("PSC"),
                              analysis_method = "LM",
                              normalization = "CSS",
                              transform = "NONE",
                              min_prevalence = 0.05,
                              correction = "BH",
                              max_significance = 0.05,
                              reference = c("PSC,TPRo"), 
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE,
                              save_scatter = FALSE)

# Save only the pairwise comparisons within root endosphere of Fuzzy Bean #
cl_root_wpsf.res <- cl_root_wpsf.maas$results
cl_root_wpsf.dares <- c()
for(i in 1:nrow(cl_root_wpsf.res)){
  if(cl_root_wpsf.res$value[i] == "TNRo" | cl_root_wpsf.res$value[i] == "TCRo"){
    if(cl_root_wpsf.res$qval[i] < 0.05){
      cl_root_wpsf.dares <- rbind(cl_root_wpsf.dares, cl_root_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the root endosphere data (Non-PSF Soil Reference)# 
cl_root_npsf.maas <- Maaslin2(input_data = root_maas$otu,
                              input_metadata = root_maas$met,
                              output = "./maas_results/cl_root_npsf.maas",
                              fixed_effects = c("PSC"),
                              analysis_method = "LM",
                              normalization = "CSS",
                              transform = "NONE",
                              min_prevalence = 0.05,
                              correction = "BH",
                              max_significance = 0.05,
                              reference = c("PSC,TNRo"), 
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE,
                              save_scatter = FALSE)

# Save only the pairwise comparisons within root endosphere of Fuzzy Bean #
cl_root_npsf.res <- cl_root_npsf.maas$results
cl_root_npsf.dares <- c()
for(i in 1:nrow(cl_root_npsf.res)){
  if(cl_root_npsf.res$value[i] == "TPRo" | cl_root_npsf.res$value[i] == "TCRo"){
    if(cl_root_npsf.res$qval[i] < 0.05){
      cl_root_npsf.dares <- rbind(cl_root_npsf.dares, cl_root_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data #
cl_comm.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/cl_comm.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,TCRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)
# Save only the pairwise comparisons within Common Soil treatments of Fuzzy Bean #
cl_comm.res <- cl_comm.maas$results
cl_comm.dares <- c()
for(i in 1:nrow(cl_comm.res)){
  if(cl_comm.res$value[i] == "SCBu"){
    if(cl_comm.res$qval[i] < 0.05){
      cl_comm.dares <- rbind(cl_comm.dares, cl_comm.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data #
cl_wpsf.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/cl_wpsf.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,TPRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)

# Save only the pairwise comparisons within PSF Soil treatments of Fuzzy Bean #
cl_wpsf.res <- cl_wpsf.maas$results
cl_wpsf.dares <- c()
for(i in 1:nrow(cl_wpsf.res)){
  if(cl_wpsf.res$value[i] == "TPBu"){
    if(cl_wpsf.res$qval[i] < 0.05){
      cl_wpsf.dares <- rbind(cl_wpsf.dares, cl_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

## Medicago ##
# Organize ASVs by number #
md_bulk.sort <- as.numeric(sub("ASV([0-9]+).*", "\\1", rownames(md_bulk$otu)))
md_bulk$otu <- md_bulk$otu[order(md_bulk.sort),]
# Maaslin2 analysis for the rhiz soil data (PSF Reference) # 
md_bulk_wpsf.maas <- Maaslin2(input_data = md_bulk$otu,
                               input_metadata = md_bulk$met,
                               output = "./maas_results/md_bulk_wpsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within bulk soils of Fuzzy Bean #
md_bulk_wpsf.res <- md_bulk_wpsf.maas$results
md_bulk_wpsf.dares <- c()
for(i in 1:nrow(md_bulk_wpsf.res)){
  if(md_bulk_wpsf.res$qval[i] < 0.05){
    md_bulk_wpsf.dares <- rbind(md_bulk_wpsf.dares, md_bulk_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the bulk soil data (PSF Reference) # 
md_bulk_npsf.maas <- Maaslin2(input_data = md_bulk$otu,
                               input_metadata = md_bulk$met,
                               output = "./maas_results/md_bulk_npsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,Non-PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within bulk soils of Fuzzy Bean #
md_bulk_npsf.res <- md_bulk_npsf.maas$results
md_bulk_npsf.dares <- c()
for(i in 1:nrow(md_bulk_npsf.res)){
  if(md_bulk_npsf.res$qval[i] < 0.05){
    md_bulk_npsf.dares <- rbind(md_bulk_npsf.dares, md_bulk_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Organize ASVs by number #
md_rhiz.sort <- as.numeric(sub("ASV([0-9]+).*", "\\1", rownames(md_rhiz$otu)))
md_rhiz$otu <- md_rhiz$otu[order(md_rhiz.sort),]

# Maaslin2 analysis for the bulk soil data (PSF Reference) # 
md_rhiz_wpsf.maas <- Maaslin2(input_data = md_rhiz$otu,
                               input_metadata = md_rhiz$met,
                               output = "./maas_results/md_rhiz_wpsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within rhiz soils of Fuzzy Bean #
md_rhiz_wpsf.res <- md_rhiz_wpsf.maas$results
md_rhiz_wpsf.dares <- c()
for(i in 1:nrow(md_rhiz_wpsf.res)){
  if(md_rhiz_wpsf.res$qval[i] < 0.05){
    md_rhiz_wpsf.dares <- rbind(md_rhiz_wpsf.dares, md_rhiz_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the rhiz soil data (PSF Reference) # 
md_rhiz_npsf.maas <- Maaslin2(input_data = md_rhiz$otu,
                               input_metadata = md_rhiz$met,
                               output = "./maas_results/md_rhiz_npsf.maas",
                               fixed_effects = c("Soils"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.32,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Soils,Non-PSF Soil"), 
                               plot_heatmap = FALSE,
                               plot_scatter = FALSE,
                               save_scatter = FALSE)

# Save only the pairwise comparisons within rhiz soils of Fuzzy Bean #
md_rhiz_npsf.res <- md_rhiz_npsf.maas$results
md_rhiz_npsf.dares <- c()
for(i in 1:nrow(md_rhiz_npsf.res)){
  if(md_rhiz_npsf.res$qval[i] < 0.05){
    md_rhiz_npsf.dares <- rbind(md_rhiz_npsf.dares, md_rhiz_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
  }
}

# Maaslin2 analysis for the root endosphere data (PSF Soil Reference)# 
md_root_wpsf.maas <- Maaslin2(input_data = root_maas$otu,
                              input_metadata = root_maas$met,
                              output = "./maas_results/md_root_wpsf.maas",
                              fixed_effects = c("PSC"),
                              analysis_method = "LM",
                              normalization = "CSS",
                              transform = "NONE",
                              min_prevalence = 0.05,
                              correction = "BH",
                              max_significance = 0.05,
                              reference = c("PSC,MPRo"), 
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE,
                              save_scatter = FALSE)

# Save only the pairwise comparisons within root endosphere of Fuzzy Bean #
md_root_wpsf.res <- md_root_wpsf.maas$results
md_root_wpsf.dares <- c()
for(i in 1:nrow(md_root_wpsf.res)){
  if(md_root_wpsf.res$value[i] == "MNRo" | md_root_wpsf.res$value[i] == "MCRo"){
    if(md_root_wpsf.res$qval[i] < 0.05){
      md_root_wpsf.dares <- rbind(md_root_wpsf.dares, md_root_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the root endosphere data (Non-PSF Soil Reference)# 
md_root_npsf.maas <- Maaslin2(input_data = root_maas$otu,
                              input_metadata = root_maas$met,
                              output = "./maas_results/md_root_npsf.maas",
                              fixed_effects = c("PSC"),
                              analysis_method = "LM",
                              normalization = "CSS",
                              transform = "NONE",
                              min_prevalence = 0.05,
                              correction = "BH",
                              max_significance = 0.05,
                              reference = c("PSC,MNRo"), 
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE,
                              save_scatter = FALSE)

# Save only the pairwise comparisons within root endosphere of Fuzzy Bean #
md_root_npsf.res <- md_root_npsf.maas$results
md_root_npsf.dares <- c()
for(i in 1:nrow(md_root_npsf.res)){
  if(md_root_npsf.res$value[i] == "MPRo" | md_root_npsf.res$value[i] == "MCRo"){
    if(md_root_npsf.res$qval[i] < 0.05){
      md_root_npsf.dares <- rbind(md_root_npsf.dares, md_root_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data # 
md_npsf.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/md_npsf.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,MNRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)
# Save only the pairwise comparisons within Non-PSF treatments of Fuzzy Bean #
md_npsf.res <- md_npsf.maas$results
md_npsf.dares <- c()
for(i in 1:nrow(md_npsf.res)){
  if(md_npsf.res$value[i] == "TNBu"){
    if(md_npsf.res$qval[i] < 0.05){
      md_npsf.dares <- rbind(md_npsf.dares, md_npsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data #
md_comm.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/md_comm.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,MCRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)
# Save only the pairwise comparisons within Common Soil treatments of Fuzzy Bean #
md_comm.res <- md_comm.maas$results
md_comm.dares <- c()
for(i in 1:nrow(md_comm.res)){
  if(md_comm.res$value[i] == "SCBu"){
    if(md_comm.res$qval[i] < 0.05){
      md_comm.dares <- rbind(md_comm.dares, md_comm.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

# Maaslin2 analysis for the Non-PSF Soil data #
md_wpsf.maas <- Maaslin2(input_data = soil_maas$otu,
                         input_metadata = soil_maas$met,
                         output = "./maas_results/md_wpsf.maas",
                         fixed_effects = c("PSC"),
                         analysis_method = "LM",
                         normalization = "CSS",
                         transform = "NONE",
                         min_prevalence = 0.05,
                         correction = "BH",
                         max_significance = 0.05,
                         reference = c("PSC,MPRh"), 
                         plot_heatmap = FALSE,
                         plot_scatter = FALSE,
                         save_scatter = FALSE)

# Save only the pairwise comparisons within PSF Soil treatments of Fuzzy Bean #
md_wpsf.res <- md_wpsf.maas$results
md_wpsf.dares <- c()
for(i in 1:nrow(md_wpsf.res)){
  if(md_wpsf.res$value[i] == "MPBu"){
    if(md_wpsf.res$qval[i] < 0.05){
      md_wpsf.dares <- rbind(md_wpsf.dares, md_wpsf.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

#### Visualizing Differential Abundance ####
### Fuzzy Bean Bulk Soil ###
## Soil Nodule ASV: ASV9(Mesorhizobium) ##
# Transform the abundance data using total sum scaling #
fb_bulk_prop.ps <- transform_sample_counts(fb_bulk.ps, function(x) x/sum(x))
decompose_ps(fb_bulk_prop.ps, "fb_bulk_prop")

# Make a data frame of just the bulk soil samples as well as the abundances of rhizobia of the PSF nodule #
fb_bulk.meso <- cbind(fb_bulk_prop$met, t(fb_bulk_prop$otu["ASV9(Mesorhizobium)",]))

# Make soils a factor with the proper levels #
fb_bulk.meso$Soils <- factor(fb_bulk.meso$Soil_Treatment, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
fb_bulk.pvsn <- "*"
fb_bulk.pvsc <- "**"
fb_bulk.nvsc <- "ns"

# Plot the differential abundance #
fb_bulk_meso.plot <- ggplot(fb_bulk.meso, aes(x = Comps, y = `ASV9(Mesorhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  coord_cartesian(ylim = c(0, 0.01)) +
  scale_y_continuous(limits = c(0,0.01), breaks = seq(0,0.01, by = 0.002), sec.axis = dup_axis(name = "")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Bulk Soil<br>ASV9(*Mesorhizobium*)") +
  ylab('Relative Abundance') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'none') +
  geom_signif(annotation = fb_bulk.pvsn,
              y_position = 0.0065, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = fb_bulk.pvsc,
              y_position = 0.007, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = fb_bulk.nvsc,
              y_position = 0.0025, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  labs(tag = "A.")
fb_bulk_meso.plot

### Fuzzy Bean Rhizosphere ###
## Soil Nodule ASV: ASV9(Mesorhizobium) ##
# Transform the abundance data using total sum scaling #
fb_rhiz_prop.ps <- transform_sample_counts(fb_rhiz.ps, function(x) x/sum(x))
decompose_ps(fb_rhiz_prop.ps, "fb_rhiz_prop")

# Make a data frame of just the rhizosphere samples as well as the abundances of rhizobia of the PSF nodule #
fb_rhiz.meso <- cbind(fb_rhiz_prop$met, t(fb_rhiz_prop$otu["ASV9(Mesorhizobium)",]))

# Make soils a factor with the proper levels #
fb_rhiz.meso$Soils <- factor(fb_rhiz.meso$Soil_Treatment, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
fb_rhiz.pvsn <- "ns"
fb_rhiz.pvsc <- "ns"
fb_rhiz.nvsc <- "*"

# Plot the differential abundance #
fb_rhiz_meso.plot <- ggplot(fb_rhiz.meso, aes(x = Comps, y = `ASV9(Mesorhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  coord_cartesian(ylim = c(0, 0.01)) +
  scale_y_continuous(limits = c(0,0.01), breaks = seq(0,0.01, by = 0.002), sec.axis = dup_axis(name = "")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Rhizosphere<br>ASV9(*Mesorhizobium*)") +
  ylab('') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'bottom') +
  geom_signif(annotation = fb_rhiz.pvsn,
              y_position = 0.0075, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = fb_rhiz.pvsc,
              y_position = 0.0085, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = fb_rhiz.nvsc,
              y_position = 0.0055, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  labs(tag = "B.")
fb_rhiz_meso.plot

### Fuzzy Bean Root Endosphere ###
## Root Nodule ASV: ASV12(Mesorhizobium) ##
# Transform the abundance data using total sum scaling #
fb_root_prop.ps <- transform_sample_counts(fb_root.ps, function(x) x/sum(x))
decompose_ps(fb_root_prop.ps, "fb_root_prop")

# Make a data frame of just the root endosphere samples as well as the abundances of rhizobia of the PSF nodule #
fb_root.meso <- cbind(fb_root_prop$met, t(fb_root_prop$otu["ASV12(Mesorhizobium)",]))

# Make soils a factor with the proper levels #
fb_root.meso$Soils <- factor(fb_root.meso$Soil.Origin, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
fb_root.pvsn <- "**"
fb_root.pvsc <- "**"
fb_root.nvsc <- "ns"

# Plot the differential abundance #
fb_root_meso.plot <- ggplot(fb_root.meso, aes(x = Comps, y = `ASV12(Mesorhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  scale_y_continuous(limits = c(0,0.35), breaks = seq(0,0.35, by = 0.05), sec.axis = dup_axis(name = "S. helvola")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Root Endosphere<br>ASV12(*Mesorhizobium*)") +
  ylab('') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'none') +
  geom_signif(annotation = fb_root.pvsn,
              y_position = 0.31, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = fb_root.pvsc,
              y_position = 0.33, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = fb_root.nvsc,
              y_position = 0.01, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  labs(tag = "C.")
fb_root_meso.plot

# Make a patchwork plot for all compartments #
(fb_bulk_meso.plot | fb_rhiz_meso.plot | fb_root_meso.plot) &
  theme(plot.tag = element_text(size = 22, face = 'bold', family = "Liberation Sans"))


### Chamaecrista Bulk Soil ###
## Soil Nodule ASV: ASV1(Bradyrhizobium) ##
# Transform the abundance data using total sum scaling #
cc_bulk_prop.ps <- transform_sample_counts(cc_bulk.ps, function(x) x/sum(x))
decompose_ps(cc_bulk_prop.ps, "cc_bulk_prop")

# Make a data frame of just the bulk soil samples as well as the abundances of rhizobia of the PSF nodule #
cc_bulk.meso <- cbind(cc_bulk_prop$met, t(cc_bulk_prop$otu["ASV1(Bradyrhizobium)",]))

# Make soils a factor with the proper levels #
cc_bulk.meso$Soils <- factor(cc_bulk.meso$Soil_Treatment, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
cc_bulk.pvsn <- "**"
cc_bulk.pvsc <- "ns"
cc_bulk.nvsc <- "**"

# Plot the differential abundance #
cc_bulk_meso.plot <- ggplot(cc_bulk.meso, aes(x = Comps, y = `ASV1(Bradyrhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  coord_cartesian(ylim = c(0, 0.1)) +
  scale_y_continuous(limits = c(0,0.11), breaks = seq(0,0.1, by = 0.02), sec.axis = dup_axis(name = "")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Bulk Soil<br>ASV1(*Bradyrhizobium*)") +
  ylab('Relative Abundance') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'none') +
  geom_signif(annotation = cc_bulk.pvsn,
              y_position = 0.09, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = cc_bulk.pvsc,
              y_position = 0.1, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = cc_bulk.nvsc,
              y_position = 0.095, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  labs(tag = "A.")
cc_bulk_meso.plot

### Chamaecrista Rhizosphere ###
## Soil Nodule ASV: ASV9(Mesorhizobium) ##
# Transform the abundance data using total sum scaling #
cc_rhiz_prop.ps <- transform_sample_counts(cc_rhiz.ps, function(x) x/sum(x))
decompose_ps(cc_rhiz_prop.ps, "cc_rhiz_prop")

# Make a data frame of just the rhizosphere samples as well as the abundances of rhizobia of the PSF nodule #
cc_rhiz.meso <- cbind(cc_rhiz_prop$met, t(cc_rhiz_prop$otu["ASV1(Bradyrhizobium)",]))

# Make soils a factor with the proper levels #
cc_rhiz.meso$Soils <- factor(cc_rhiz.meso$Soil_Treatment, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
cc_rhiz.pvsn <- "*"
cc_rhiz.pvsc <- "ns"
cc_rhiz.nvsc <- "**"

# Plot the differential abundance #
cc_rhiz_meso.plot <- ggplot(cc_rhiz.meso, aes(x = Comps, y = `ASV1(Bradyrhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  coord_cartesian(ylim = c(0, 0.8)) +
  scale_y_continuous(limits = c(0,0.8), breaks = seq(0,0.8, by = 0.1), sec.axis = dup_axis(name = "")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Rhizosphere<br>ASV1(*Bradyrhizobium*)") +
  ylab('') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'bottom') +
  geom_signif(annotation = cc_rhiz.pvsn,
              y_position = 0.65, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = cc_rhiz.pvsc,
              y_position = 0.75, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = cc_rhiz.nvsc,
              y_position = 0.7, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  labs(tag = "B.")
cc_rhiz_meso.plot

### Chamaecrista Root Endosphere ###
## Root Nodule ASV: ASV3(Bradyrhizobium) ##
# Transform the abundance data using total sum scaling #
cc_root_prop.ps <- transform_sample_counts(cc_root.ps, function(x) x/sum(x))
decompose_ps(cc_root_prop.ps, "cc_root_prop")

# Make a data frame of just the root endosphere samples as well as the abundances of rhizobia of the PSF nodule #
cc_root.meso <- cbind(cc_root_prop$met, t(cc_root_prop$otu[c("ASV3(Bradyrhizobium)", "ASV6(Bradyrhizobium)", "ASV14(Bradyrhizobium"),]))

# Make soils a factor with the proper levels #
cc_root.meso$Soils <- factor(cc_root.meso$Soil.Origin, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
cc_root.pvsn <- "ns"
cc_root.pvsc <- "**"
cc_root.nvsc <- "**"

# Plot the differential abundance #
cc_root_meso.plot <- ggplot(cc_root.meso, aes(x = Comps, y = `ASV3(Bradyrhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.2), sec.axis = dup_axis(name = "C. fasciculata")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Root Endosphere<br>ASV3(*Bradyrhizobium*)") +
  ylab('') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'none') +
  geom_signif(annotation = cc_root.pvsn,
              y_position = 0.85, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = cc_root.pvsc,
              y_position = 0.95, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = cc_root.nvsc,
              y_position = 0.9, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  labs(tag = "C.")
cc_root_meso.plot

# Make a patchwork plot for all compartments #
(cc_bulk_meso.plot | cc_rhiz_meso.plot | cc_root_meso.plot) &
  theme(plot.tag = element_text(size = 22, face = 'bold', family = "Liberation Sans"))

## Bonus Plots for other noduel rhizobia in roots ##
# Plot the differential abundance #
cc_root_meso6.plot <- ggplot(cc_root.meso, aes(x = Comps, y = `ASV6(Bradyrhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  scale_y_continuous(limits = c(0,0.4), breaks = seq(0,0.4, by = 0.05), sec.axis = dup_axis(name = "C. fasciculata")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Root Endosphere<br>ASV6(*Bradyrhizobium*)") +
  ylab('') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'none') +
  geom_signif(annotation = "**",
              y_position = 0.325, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = "ns",
              y_position = 0.375, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = "**",
              y_position = 0.35, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  labs(tag = "A.")
cc_root_meso6.plot

# Plot the differential abundance #
cc_root_meso14.plot <- ggplot(cc_root.meso, aes(x = Comps, y = `ASV14(Bradyrhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  scale_y_continuous(limits = c(0,0.4), breaks = seq(0,0.4, by = 0.05), sec.axis = dup_axis(name = "C. fasciculata")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Root Endosphere<br>ASV14(*Bradyrhizobium*)") +
  ylab('') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'none') +
  geom_signif(annotation = "***",
              y_position = 0.325, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = "ns",
              y_position = 0.375, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = "***",
              y_position = 0.35, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  lab(tags = "B.")
cc_root_meso14.plot

(cc_root_meso.plot | cc_root_meso6.plot | cc_root_meso14.plot)

### Desmodium Bulk Soil ###
## Soil Nodule ASV: ASV1(Bradyrhizobium) ##
# Transform the abundance data using total sum scaling #
ds_bulk_prop.ps <- transform_sample_counts(ds_bulk.ps, function(x) x/sum(x))
decompose_ps(ds_bulk_prop.ps, "ds_bulk_prop")

# Make a data frame of just the bulk soil samples as well as the abundances of rhizobia of the PSF nodule #
ds_bulk.meso <- cbind(ds_bulk_prop$met, t(ds_bulk_prop$otu["ASV1(Bradyrhizobium)",]))

# Make soils a factor with the proper levels #
ds_bulk.meso$Soils <- factor(ds_bulk.meso$Soil_Treatment, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
ds_bulk.pvsn <- "***"
ds_bulk.pvsc <- "***"
ds_bulk.nvsc <- "**"

# Plot the differential abundance #
ds_bulk_meso.plot <- ggplot(ds_bulk.meso, aes(x = Comps, y = `ASV1(Bradyrhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  coord_cartesian(ylim = c(0, 0.08)) +
  scale_y_continuous(limits = c(0,0.08), breaks = seq(0,0.08, by = 0.01), sec.axis = dup_axis(name = "")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Bulk Soil<br>ASV1(*Bradyrhizobium*)") +
  ylab('Relative Abundance') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'none') +
  geom_signif(annotation = ds_bulk.pvsn,
              y_position = 0.065, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = ds_bulk.pvsc,
              y_position = 0.075, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = ds_bulk.nvsc,
              y_position = 0.07, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  labs(tag = "A.")
ds_bulk_meso.plot

### Desmodium Rhizosphere ###
## Soil Nodule ASV: ASV1(Bradyrhizobium) ##
# Transform the abundance data using total sum scaling #
ds_rhiz_prop.ps <- transform_sample_counts(ds_rhiz.ps, function(x) x/sum(x))
decompose_ps(ds_rhiz_prop.ps, "ds_rhiz_prop")

# Make a data frame of just the rhizosphere samples as well as the abundances of rhizobia of the PSF nodule #
ds_rhiz.meso <- cbind(ds_rhiz_prop$met, t(ds_rhiz_prop$otu["ASV1(Bradyrhizobium)",]))

# Make soils a factor with the proper levels #
ds_rhiz.meso$Soils <- factor(ds_rhiz.meso$Soil_Treatment, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
ds_rhiz.pvsn <- "ns"
ds_rhiz.pvsc <- "***"
ds_rhiz.nvsc <- "***"

# Plot the differential abundance #
ds_rhiz_meso.plot <- ggplot(ds_rhiz.meso, aes(x = Comps, y = `ASV1(Bradyrhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  coord_cartesian(ylim = c(0, 0.1)) +
  scale_y_continuous(limits = c(0,0.1), breaks = seq(0,0.1, by = 0.025), sec.axis = dup_axis(name = "")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Rhizosphere<br>ASV1(*Bradyrhizobium*)") +
  ylab('') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'bottom') +
  geom_signif(annotation = ds_rhiz.pvsn,
              y_position = 0.085, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = ds_rhiz.pvsc,
              y_position = 0.095, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = ds_rhiz.nvsc,
              y_position = 0.09, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  labs(tag = "B.")
ds_rhiz_meso.plot

### Desmodium Root Endosphere ###
## Root Nodule ASV: ASV6(Bradyrhizobium) ##
# Transform the abundance data using total sum scaling #
ds_root_prop.ps <- transform_sample_counts(ds_root.ps, function(x) x/sum(x))
decompose_ps(ds_root_prop.ps, "ds_root_prop")

# Make a data frame of just the root endosphere samples as well as the abundances of rhizobia of the PSF nodule #
ds_root.meso <- cbind(ds_root_prop$met, t(ds_root_prop$otu["ASV7(Bradyrhizobium)",]))

# Make soils a factor with the proper levels #
ds_root.meso$Soils <- factor(ds_root.meso$Soil.Origin, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
ds_root.pvsn <- "ns"
ds_root.pvsc <- "ns"
ds_root.nvsc <- "ns"

# Plot the differential abundance #
ds_root_meso.plot <- ggplot(ds_root.meso, aes(x = Comps, y = `ASV7(Bradyrhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  scale_y_continuous(limits = c(0,0.4), breaks = seq(0,0.4, by = 0.1), sec.axis = dup_axis(name = "D. illinoense")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Root Endosphere<br>ASV7(*Bradyrhizobium*)") +
  ylab('') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'none') +
  geom_signif(annotation = ds_root.pvsn,
              y_position = 0.3, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = ds_root.pvsc,
              y_position = 0.36, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = ds_root.nvsc,
              y_position = 0.33, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  labs(tag = "C.")
ds_root_meso.plot

# Make a patchwork plot for all compartments #
(ds_bulk_meso.plot | ds_rhiz_meso.plot | ds_root_meso.plot) &
  theme(plot.tag = element_text(size = 22, face = 'bold', family = "Liberation Sans"))

### Hog Peanut Bulk Soil ###
## Soil Nodule ASV: ASV1(Bradyrhizobium) ##
# Transform the abundance data using total sum scaling #
hp_bulk_prop.ps <- transform_sample_counts(hp_bulk.ps, function(x) x/sum(x))
decompose_ps(hp_bulk_prop.ps, "hp_bulk_prop")

# Make a data frame of just the bulk soil samples as well as the abundances of rhizobia of the PSF nodule #
hp_bulk.meso <- cbind(hp_bulk_prop$met, t(hp_bulk_prop$otu["ASV1(Bradyrhizobium)",]))

# Make soils a factor with the proper levels #
hp_bulk.meso$Soils <- factor(hp_bulk.meso$Soil_Treatment, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
hp_bulk.pvsn <- "***"
hp_bulk.pvsc <- "ns"
hp_bulk.nvsc <- "***"

# Plot the differential abundance #
hp_bulk_meso.plot <- ggplot(hp_bulk.meso, aes(x = Comps, y = `ASV1(Bradyrhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  coord_cartesian(ylim = c(0, 0.05)) +
  scale_y_continuous(limits = c(0,0.05), breaks = seq(0,0.05, by = 0.01), sec.axis = dup_axis(name = "")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Bulk Soil<br>ASV1(*Bradyrhizobium*)") +
  ylab('Relative Abundance') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'none') +
  geom_signif(annotation = hp_bulk.pvsn,
              y_position = 0.038, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = hp_bulk.pvsc,
              y_position = 0.046, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = hp_bulk.nvsc,
              y_position = 0.042, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  labs(tag = "A.")
hp_bulk_meso.plot

### Hog Peanut Rhizosphere ###
## Soil Nodule ASV: ASV1(Bradyrhizobium) ##
# Transform the abundance data using total sum scaling #
hp_rhiz_prop.ps <- transform_sample_counts(hp_rhiz.ps, function(x) x/sum(x))
decompose_ps(hp_rhiz_prop.ps, "hp_rhiz_prop")

# Make a data frame of just the rhizosphere samples as well as the abundances of rhizobia of the PSF nodule #
hp_rhiz.meso <- cbind(hp_rhiz_prop$met, t(hp_rhiz_prop$otu["ASV1(Bradyrhizobium)",]))

# Make soils a factor with the proper levels #
hp_rhiz.meso$Soils <- factor(hp_rhiz.meso$Soil_Treatment, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
hp_rhiz.pvsn <- "ns"
hp_rhiz.pvsc <- "ns"
hp_rhiz.nvsc <- "ns"

# Plot the differential abundance #
hp_rhiz_meso.plot <- ggplot(hp_rhiz.meso, aes(x = Comps, y = `ASV1(Bradyrhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  coord_cartesian(ylim = c(0, 0.1)) +
  scale_y_continuous(limits = c(0,0.1), breaks = seq(0,0.1, by = 0.025), sec.axis = dup_axis(name = "")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Rhizosphere<br>ASV1(*Bradyrhizobium*)") +
  ylab('') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'bottom') +
  geom_signif(annotation = hp_rhiz.pvsn,
              y_position = 0.085, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = hp_rhiz.pvsc,
              y_position = 0.096, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = hp_rhiz.nvsc,
              y_position = 0.09, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  labs(tag = "B.")
hp_rhiz_meso.plot

### Hog Peanut Root Endosphere ###
## Root Nodule ASV: ASV6(Bradyrhizobium) ##
# Transform the abundance data using total sum scaling #
hp_root_prop.ps <- transform_sample_counts(hp_root.ps, function(x) x/sum(x))
decompose_ps(hp_root_prop.ps, "hp_root_prop")

# Make a data frame of just the root endosphere samples as well as the abundances of rhizobia of the PSF nodule #
hp_root.meso <- cbind(hp_root_prop$met, t(hp_root_prop$otu["ASV7(Bradyrhizobium)",]))

# Make soils a factor with the proper levels #
hp_root.meso$Soils <- factor(hp_root.meso$Soil.Origin, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
hp_root.pvsn <- "ns"
hp_root.pvsc <- "ns"
hp_root.nvsc <- "ns"

# Since ASV9 was not found in the root endosphere samples, just make a column of zeros #
hp_root.meso$`ASV9(Bradyrhizobium)` <- 0

# Plot the differential abundance #
hp_root_meso.plot <- ggplot(hp_root.meso, aes(x = Comps, y = `ASV9(Bradyrhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  scale_y_continuous(limits = c(0,0.8), breaks = seq(0,0.8, by = 0.2), sec.axis = dup_axis(name = "A. bracteata")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Root Endosphere<br>ASV9(*Bradyrhizobium*)") +
  ylab('') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'none') +
  labs(tag = "C.")
hp_root_meso.plot

# Make a patchwork plot for all compartments #
(hp_bulk_meso.plot | hp_rhiz_meso.plot | hp_root_meso.plot) &
  theme(plot.tag = element_text(size = 22, face = 'bold', family = "Liberation Sans"))

### Clover Bulk Soil ###
## Soil Nodule ASV: ASV1(Bradyrhizobium) ##
# Transform the abundance data using total sum scaling #
cl_bulk_prop.ps <- transform_sample_counts(cl_bulk.ps, function(x) x/sum(x))
decompose_ps(cl_bulk_prop.ps, "cl_bulk_prop")

# Make a data frame of just the bulk soil samples as well as the abundances of rhizobia of the PSF nodule #
cl_bulk.meso <- cbind(cl_bulk_prop$met, t(cl_bulk_prop$otu["ASV2(Agrobacterium)",]))

# Make soils a factor with the proper levels #
cl_bulk.meso$Soils <- factor(cl_bulk.meso$Soil_Treatment, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
cl_bulk.pvsn <- "ns"
cl_bulk.pvsc <- "ns"
cl_bulk.nvsc <- "ns"

# Plot the differential abundance #
cl_bulk_meso.plot <- ggplot(cl_bulk.meso, aes(x = Comps, y = `ASV2(Agrobacterium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  coord_cartesian(ylim = c(0, 0.01)) +
  scale_y_continuous(limits = c(0,0.01), breaks = seq(0,0.01, by = 0.002), sec.axis = dup_axis(name = "")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Bulk Soil<br>ASV2(*Agrobacterium*)") +
  ylab('Relative Abundance') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'bottom') +
  geom_signif(annotation = cl_bulk.pvsn,
              y_position = 0.038, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = cl_bulk.pvsc,
              y_position = 0.046, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = cl_bulk.nvsc,
              y_position = 0.042, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  labs(tag = "A.")
cl_bulk_meso.plot

### Clover Rhizosphere ###
## Soil Nodule ASV: ASV1(Bradyrhizobium) ##
# Transform the abundance data using total sum scaling #
cl_rhiz_prop.ps <- transform_sample_counts(cl_rhiz.ps, function(x) x/sum(x))
decompose_ps(cl_rhiz_prop.ps, "cl_rhiz_prop")

# Make a data frame of just the rhizosphere samples as well as the abundances of rhizobia of the PSF nodule #
cl_rhiz.meso <- cbind(cl_rhiz_prop$met, t(cl_rhiz_prop$otu["ASV2(Agrobacterium)",]))

# Make soils a factor with the proper levels #
cl_rhiz.meso$Soils <- factor(cl_rhiz.meso$Soil_Treatment, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
cl_rhiz.pvsn <- "ns"
cl_rhiz.pvsc <- "ns"
cl_rhiz.nvsc <- "ns"

# Plot the differential abundance #
cl_rhiz_meso.plot <- ggplot(cl_rhiz.meso, aes(x = Comps, y = `ASV2(Agrobacterium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  coord_cartesian(ylim = c(0, 0.01)) +
  scale_y_continuous(limits = c(0,0.01), breaks = seq(0,0.01, by = 0.0025), sec.axis = dup_axis(name = "")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Rhizosphere<br>ASV2(*Agrobacterium*)") +
  ylab('') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'none') +
  labs(tag = "B.")
cl_rhiz_meso.plot

### Hog Peanut Root Endosphere ###
## Root Nodule ASV: ASV6(Bradyrhizobium) ##
# Transform the abundance data using total sum scaling #
cl_root_prop.ps <- transform_sample_counts(cl_root.ps, function(x) x/sum(x))
decompose_ps(cl_root_prop.ps, "cl_root_prop")

# Make a data frame of just the root endosphere samples as well as the abundances of rhizobia of the PSF nodule #
cl_root.meso <- cbind(cl_root_prop$met, t(cl_root_prop$otu["ASV1(Rhizobium)",]))

# Make soils a factor with the proper levels #
cl_root.meso$Soils <- factor(cl_root.meso$Soil.Origin, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
cl_root.pvsn <- "ns"
cl_root.pvsc <- "ns"
cl_root.nvsc <- "ns"

# Plot the differential abundance #
cl_root_meso.plot <- ggplot(cl_root.meso, aes(x = Comps, y = `ASV1(Rhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = 0.2), sec.axis = dup_axis(name = "T. repens")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Root Endosphere<br>ASV1(*Rhizobium*)") +
  ylab('') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'none') +
  geom_signif(annotation = cl_root.pvsn,
              y_position = 0.85, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = cl_root.pvsc,
              y_position = 0.96, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = cl_root.nvsc,
              y_position = 0.9, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  labs(tag = "C.")
cl_root_meso.plot

# Make a patchwork plot for all compartments #
(cl_bulk_meso.plot + plot_layout(guides = "keep") & theme(legend.position = "bottom", legend.justification = c(-0.85,0)) | cl_rhiz_meso.plot   | cl_root_meso.plot) &
  theme(plot.tag = element_text(size = 22, face = 'bold', family = "Liberation Sans"))


### Medicago Bulk Soil ###
## Soil Nodule ASV: ASV1(Bradyrhizobium) ##
# Transform the abundance data using total sum scaling #
md_bulk_prop.ps <- transform_sample_counts(md_bulk.ps, function(x) x/sum(x))
decompose_ps(md_bulk_prop.ps, "md_bulk_prop")

# Make a data frame of just the bulk soil samples as well as the abundances of rhizobia of the PSF nodule #
md_bulk.meso <- cbind(md_bulk_prop$met, t(md_bulk_prop$otu["ASV2(Agrobacterium)",]))

# Make soils a factor with the proper levels #
md_bulk.meso$Soils <- factor(md_bulk.meso$Soil_Treatment, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
md_bulk.pvsn <- "ns"
md_bulk.pvsc <- "ns"
md_bulk.nvsc <- "ns"

# Plot the differential abundance #
md_bulk_meso.plot <- ggplot(md_bulk.meso, aes(x = Comps, y = `ASV2(Agrobacterium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  coord_cartesian(ylim = c(0, 0.01)) +
  scale_y_continuous(limits = c(0,0.01), breaks = seq(0,0.01, by = 0.0025), sec.axis = dup_axis(name = "")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Bulk Soil<br>ASV2(*Agrobacterium*)") +
  ylab('Relative Abundance') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'none') +
  geom_signif(annotation = md_bulk.pvsn,
              y_position = 0.038, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = md_bulk.pvsc,
              y_position = 0.046, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = md_bulk.nvsc,
              y_position = 0.042, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  labs(tag = "A.")
md_bulk_meso.plot

### Medicago Rhizosphere ###
## Soil Nodule ASV: ASV1(Bradyrhizobium) ##
# Transform the abundance data using total sum scaling #
md_rhiz_prop.ps <- transform_sample_counts(md_rhiz.ps, function(x) x/sum(x))
decompose_ps(md_rhiz_prop.ps, "md_rhiz_prop")

# Make a data frame of just the rhizosphere samples as well as the abundances of rhizobia of the PSF nodule #
md_rhiz.meso <- cbind(md_rhiz_prop$met, t(md_rhiz_prop$otu["ASV2(Agrobacterium)",]))

# Make soils a factor with the proper levels #
md_rhiz.meso$Soils <- factor(md_rhiz.meso$Soil_Treatment, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
md_rhiz.pvsn <- "ns"
md_rhiz.pvsc <- "ns"
md_rhiz.nvsc <- "ns"

# Plot the differential abundance #
md_rhiz_meso.plot <- ggplot(md_rhiz.meso, aes(x = Comps, y = `ASV2(Agrobacterium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  coord_cartesian(ylim = c(0, 0.01)) +
  scale_y_continuous(limits = c(0,0.01), breaks = seq(0,0.01, by = 0.0025), sec.axis = dup_axis(name = "")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Rhizosphere<br>ASV2(*Agrobacterium*)") +
  ylab('') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'bottom') +
  geom_signif(annotation = md_rhiz.pvsn,
              y_position = 0.085, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = md_rhiz.pvsc,
              y_position = 0.096, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = md_rhiz.nvsc,
              y_position = 0.09, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  labs(tag = "B.")
md_rhiz_meso.plot

### Medicago Root Endosphere ###
## Root Nodule ASV: ASV6(Bradyrhizobium) ##
# Transform the abundance data using total sum scaling #
md_root_prop.ps <- transform_sample_counts(md_root.ps, function(x) x/sum(x))
decompose_ps(md_root_prop.ps, "md_root_prop")

# Make a data frame of just the root endosphere samples as well as the abundances of rhizobia of the PSF nodule #
md_root.meso <- cbind(md_root_prop$met, t(md_root_prop$otu[c("ASV4(Sinorhizobium)", "ASV2(Sinorhizobium)"),]))

# Make soils a factor with the proper levels #
md_root.meso$Soils <- factor(md_root.meso$Soil.Origin, levels = c("Common Soil", "Non-PSF Soil", "PSF Soil"))

# Make annotations for significant differences #
md_root.pvsn <- "***"
md_root.pvsc <- "***"
md_root.nvsc <- "ns"

# Plot the differential abundance #
md_root_meso.plot <- ggplot(md_root.meso, aes(x = Comps, y = `ASV4(Sinorhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  coord_cartesian(ylim = c(0,0.65)) +
  scale_y_continuous(limits = c(0,0.6), breaks = seq(0,0.6, by = 0.15), sec.axis = dup_axis(name = "M. truncatula")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Root Endosphere<br>ASV4(*Sinorhizobium*)") +
  ylab('') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'none') +
  geom_signif(annotation = md_root.pvsn,
              y_position = 0.50, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = md_root.pvsc,
              y_position = 0.6, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = md_root.nvsc,
              y_position = 0.55, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  labs(tag = "C.")
md_root_meso.plot

# Make a patchwork plot for all compartments #
(md_bulk_meso.plot | md_rhiz_meso.plot | md_root_meso.plot) &
  theme(plot.tag = element_text(size = 22, face = 'bold', family = "Liberation Sans"))

# Plot the differential abundance #
md_root_meso2.plot <- ggplot(md_root.meso, aes(x = Comps, y = `ASV2(Sinorhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Soils`)) +
  coord_cartesian(ylim = c(0,0.95)) +
  scale_y_continuous(limits = c(0,0.9), breaks = seq(0,0.9, by = 0.15), sec.axis = dup_axis(name = "M. truncatula")) +
  scale_fill_manual(name = "PSF History", values = c("Common Soil" = "white", "Non-PSF Soil" = "gray", "PSF Soil" = "#4D4D4D")) +
  scale_x_discrete(labels = "Root Endosphere<br>ASV2(*Sinorhizobium*)") +
  ylab('') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = ggtext::element_markdown(color = 'black', size = 32),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 22),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90, face = "bold.italic"),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'bottom') +
  geom_signif(annotation = "**",
              y_position = 0.72, xmin = 1, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = "***",
              y_position = 0.86, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65) +
  geom_signif(annotation = "*",
              y_position = 0.81, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 20,
              vjust = 0.65)
md_root_meso2.plot

(md_root_meso.plot | md_root_meso2.plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

#### Alpha Rarefaction Curve ####
library(amplicon); packageVersion('amplicon')

# Perform the rarefaction using all soil samples # 
soild <- alpha_rare_all(otu = soil$otu, map = soil$met, method = "diversity_shannon", start = 0, step = 2500, group = "Plants", ps = soil.ps)
soil.radf <- soil.rare[[2]]

# Join the metadata to the rarefaction dataframe joined based on ID of the sample #
soil.rfdf <- c()
for(i in 1:nrow(soil.radf)){
  if(soil.radf$ID[i] %in% rownames(soil$met)){
    hold <- cbind(soil.radf[i,], soil$met[soil.radf$ID[i],])
    soil.rfdf <- rbind(soil.rfdf, hold)
    hold <- c()
  }
}
soil.radf <- soil.rfdf

# Plot the Alpha-rarefaction Curve #
ggplot(soil.radf, aes(x = i, y = index, group = ID, color = Group)) +
  geom_smooth(span = 0.7, se = FALSE, method = "loess") +
  theme_bw() +
  xlab('Sequencing Coverage') +
  ylab("Shannon Diversity") +
  scale_color_manual(name = "Plant Species", labels = c(expression(italic('S. helvola')), expression(italic('C. fasciculata')), expression(italic('D. illinoense')), expression(italic('A. bracteata')), expression(italic('T. repens')), expression(italic('M. truncatula')), "Common Soil"), values = c("#A6CEE3","#1F78B4","#FDBF6F", "#FF7F00","#CAB2D6", "#6A3D9A")) +
  theme(axis.title = element_text(size = 22, color = 'black', face = 'bold', family = 'Liberation Sans'),
        axis.text = element_text(size = 18, color = 'black', face = 'bold', family = 'Liberation Sans'),
        legend.title = element_text(size = 22, color = 'black', face = 'bold', family = 'Liberation Sans'),
        legend.text = element_text(size = 18, color = 'black', face = 'bold.italic'))
  
# Perform the rarefaction using all root endosphere samples # 
root.rare <- alpha_rare_all(otu = root$otu, map = root$met, method = "diversity_shannon", start = 0, step = 2500, group = "Plants", ps = root.ps)
root.radf <- root.rare[[2]]

# Join the metadata to the rarefaction dataframe joined based on ID of the sample #
root.rfdf <- c()
for(i in 1:nrow(root.radf)){
  if(root.radf$ID[i] %in% rownames(root$met)){
    hold <- cbind(root.radf[i,], root$met[root.radf$ID[i],])
    root.rfdf <- rbind(root.rfdf, hold)
    hold <- c()
  }
}
root.radf <- root.rfdf

# Plot the Alpha-rarefaction Curve #
ggplot(root.radf, aes(x = i, y = index, group = ID, color = Group)) +
  geom_smooth(span = 0.7, se = FALSE, method = "loess") +
  theme_bw() +
  xlab('Sequencing Coverage') +
  ylab("Shannon Diversity") +
  scale_color_manual(name = "Plant Species", labels = c(expression(italic('S. helvola')), expression(italic('C. fasciculata')), expression(italic('D. illinoense')), expression(italic('A. bracteata')), expression(italic('T. repens')), expression(italic('M. truncatula')), "Common Soil"), values = c("#A6CEE3","#1F78B4","#FDBF6F", "#FF7F00","#CAB2D6", "#6A3D9A")) +
  theme(axis.title = element_text(size = 22, color = 'black', face = 'bold', family = 'Liberation Sans'),
        axis.text = element_text(size = 18, color = 'black', face = 'bold', family = 'Liberation Sans'),
        legend.title = element_text(size = 22, color = 'black', face = 'bold', family = 'Liberation Sans'),
        legend.text = element_text(size = 18, color = 'black', face = 'bold.italic'))

save.image("PSF.RData")
