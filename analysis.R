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

# Rscript ~/PSF_MBIOME/analysis.R --raw_soil ~/test/reads/raw/soil_reads/ --raw_root ~/test/read/raw/endo_reads/ --soil_metadata ~/test/metadata/soil_metadata.csv --root_metadata ~/test/metadata/endo_metadata.csv --pheno ~/test/nodnbio.csv #

#### Argument Parsing ####
library(optparse)
option_list <- list(
            make_option("--raw_soil", type = "character", help = "filepath containing the raw, untrimmed reads for the reads generated from the v4 primers (bulk soil, rhizosphere, and some nodule samples)"),
            make_option("--raw_root", type = "character", help = "filepath containing the raw, untrimmed reads for the reads generated from the v5-v7 primers( root endosphere and some nodule samples)"),
            make_option("--soil_metadata", type = "character", help = "filepath that contains the Comma Separated Values (csv) file of the soil metadata"),
            make_option("--root_metadata", type = "character", help = "filepath that contains the Comma Separated Values (csv) file of the root metadata"),
            make_option("--pheno", type = "character", help = "filepath that contains the Comma Separated Values (csv) file of the phenotype data (biomass and nodule counts)"))

opt <- parse_args(OptionParser(option_list=option_list))
soil.dir <- opt$raw_soil
root.dir <- opt$raw_root
metadata <- opt$metadata
nodnbio <- opt$pheno
cutadapt <- opt$cutadapt

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
system(paste0("fastqc --noextract ",soil.dir, "*fastq.gz -o ./QC/raw_qc/raw_root_qc"))
system("rm ./QC/raw_qc/raw_soil_qc/*.zip")
system("rm ./QC/raw_qc/raw_root_qc/*.zip")

save.image("./test.RData")
