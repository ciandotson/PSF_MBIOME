setwd("~/PSF")

#### Package Installation ####
install.packages('remotes')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.20")

BiocManager::install("phyloseq", version = "3.20")

BiocManager::install('microbiome')

BiocManager::install('DECIPHER')

install.packages('vegan')

install.packages('dplyr')

install.packages('ggplot2')

BiocManager::install('Maaslin2')

remotes::install_github("microsud/microbiomeutilities")

BiocManager::install('NetCoMi')

remotes::install_github("zdk123/SpiecEasi")

remotes::install_github('GraceYoon/SPRING')

remotes::install_github("stefpeschel/NetCoMi", 
                        repos = c("https://cloud.r-project.org/",
                                  BiocManager::repositories()))
install.packages('cgwtools')

install.packages('multcompView')

install.packages('patchwork')

BiocManager::install('msa')

#### Package Depot ####
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(microbiome); packageVersion('microbiome')
library(vegan); packageVersion('vegan')
library(dplyr); packageVersion('dplyr')
library(microbiomeutilities); packageVersion('microbiomeutilities')
library(ggplot2); packageVersion('ggplot2')
library(Maaslin2); packageVersion('Maaslin2')
library(NetCoMi); packageVersion('NetCoMi')
library(Polychrome); packageVersion('Polychrome')
library(ggprism); packageVersion('ggprism')
library(cgwtools); packageVersion('cgwtools')
library(multcompView); packageVersion('multcompView')
library(patchwork); packageVersion('patchwork')
library(msa); packageVersion('msa')

#### Nodule Count and Biomass Data Visualization ####
nodnbio.data <- read.csv2('~/Cian_Legume16S/nodnbio.csv', sep = ',')
colnames(nodnbio.data) <- c('Plant_Sample', 'Soil_Treatment', 'Nodule_Count', 'Aboveground_Biomass')

nodnbio.data$Aboveground_Biomass <- as.numeric(nodnbio.data$Aboveground_Biomass)

for(i in 1:nrow(nodnbio.data)){
  nodnbio.data$Grouped[i] <- paste0(nodnbio.data$Plant_Sample[i], '; ', nodnbio.data$Soil_Treatment[i])
}

nodnbio.mnsd <- nodnbio.data %>%
  group_by(Plant_Sample, Soil_Treatment) %>%
  summarize(
    nod.mean = mean(Nodule_Count),
    bio.mean = mean(Aboveground_Biomass),
    nod.sd = sd(Nodule_Count),
    bio.sd = sd(Aboveground_Biomass),
    .groups = "drop" # Prevent grouping in the result
  )

nodnbio.data$Soil_Treatment <- gsub('Non-PSF Soil', 'Non PSF Soil', nodnbio.data$Soil_Treatment)

fb_nod.data <- filter(nodnbio.data, Plant_Sample == 'S. helvola')
cc_nod.data <- filter(nodnbio.data, Plant_Sample == 'C. fasciculata')
ds_nod.data <- filter(nodnbio.data, Plant_Sample == 'D. illinoense')
hp_nod.data <- filter(nodnbio.data, Plant_Sample == 'A. braceteata')
cl_nod.data <- filter(nodnbio.data, Plant_Sample == 'T. repens')
md_nod.data <- filter(nodnbio.data, Plant_Sample == 'M. truncatula')

fb_nod.aov <- aov(Nodule_Count~Soil_Treatment, fb_nod.data)
cc_nod.aov <- aov(Nodule_Count~Soil_Treatment, cc_nod.data)
ds_nod.aov <- aov(Nodule_Count~Soil_Treatment, ds_nod.data)
hp_nod.aov <- aov(Nodule_Count~Soil_Treatment, hp_nod.data)
cl_nod.aov <- aov(Nodule_Count~Soil_Treatment, cl_nod.data)
md_nod.aov <- aov(Nodule_Count~Soil_Treatment, md_nod.data)

summary(fb_nod.aov)
summary(cc_nod.aov)
summary(ds_nod.aov)
summary(hp_nod.aov)
summary(cl_nod.aov)
summary(md_nod.aov)

fb_nod.hsd <- TukeyHSD(fb_nod.aov)
cc_nod.hsd <- TukeyHSD(cc_nod.aov)
ds_nod.hsd <- TukeyHSD(ds_nod.aov)
hp_nod.hsd <- TukeyHSD(hp_nod.aov)
cl_nod.hsd <- TukeyHSD(cl_nod.aov)
md_nod.hsd <- TukeyHSD(md_nod.aov)

fb_nod.let <- multcompLetters4(fb_nod.aov, fb_nod.hsd)
cc_nod.let <- multcompLetters4(cc_nod.aov, cc_nod.hsd)
ds_nod.let <- multcompLetters4(ds_nod.aov, ds_nod.hsd)
hp_nod.let <- multcompLetters4(hp_nod.aov, hp_nod.hsd)
cl_nod.let <- multcompLetters4(cl_nod.aov, cl_nod.hsd)
md_nod.let <- multcompLetters4(md_nod.aov, md_nod.hsd)

nod.let <- c(hp_nod.let$Soil_Treatment$Letters,
                 cc_nod.let$Soil_Treatment$Letters,
                 ds_nod.let$Soil_Treatment$Letters,
                 md_nod.let$Soil_Treatment$Letters,
                 fb_nod.let$Soil_Treatment$Letters,
                 cl_nod.let$Soil_Treatment$Letters)

nod.data <- cbind(nodnbio.mnsd, nod.let)

nod.data$nod.let[8] <- 'a'
nod.data$nod.let[14] <- 'a'
nod.data$nod.let[16] <- 'ab'
nod.data$nod.let[17] <- 'a'
nod.data$nod.let[2:3] <- ''


nod.data$nod.sd[3] <- 0

nod.data$Group <- factor(nod.data$Plant_Sample, levels = c('T. repens', 'M. truncatula', 'S. helvola', 'C. fasciculata', 'D. illinoense', 'A. braceteata'))
nod.data <- filter(nod.data, Plant_Sample != 'A. braceteata')
ggplot(nod.data, aes(x = Group, y = nod.mean, fill = Soil_Treatment, color = Soil_Treatment)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = nod.mean, ymax = nod.mean + nod.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = nod.let, y = nod.mean + nod.sd + 1), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 10) +
  ylab('Nodule Count') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  theme(legend.position = 'top',
        legend.text = element_text(size = 24),
        axis.text.x = element_text(size = 22, face = c('bold.italic')),
        axis.title.y = element_text(size = 24),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm')) +
  scale_y_continuous(limits = c(0,65),expand = expansion(mult = c(0, 0.05)))

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
bio.data$bio.let[2:3] <- ''

bio.data$bio.sd[3] <- 0

bio.data$Group <- factor(bio.data$Plant_Sample, levels = c('T. repens', 'M. truncatula', 'S. helvola', 'C. fasciculata', 'D. illinoense', 'A. braceteata'))
bio.data <- filter(bio.data, Plant_Sample != 'A. braceteata')

ggplot(bio.data, aes(x = Group, y = bio.mean, fill = Soil_Treatment, color = Soil_Treatment)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = bio.mean, ymax = bio.mean + bio.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = bio.let, y = bio.mean + bio.sd + 0.01), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 10) +
  ylab('Aboveground Biomass (grams)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  theme(legend.position = 'top',
        legend.text = element_text(size = 24),
        axis.text.x = element_text(size = 22, face = c('bold.italic')),
        axis.title.y = element_text(size = 24),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.spacing.x = unit(6, 'cm')) +
  scale_y_continuous(limits = c(0,1.5),expand = expansion(mult = c(0, 0.05)))

#### DADA2 Pipeline for Soil Communities ####
setwd("~/Cian_Legume16S")

path <- "~/Cian_Legume16S/ptrimmed_reads"
list.files(path)

# separate forward and reverse reads #
fnFs <- sort(list.files(path, pattern="_R1_001.trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.trimmed.fastq", full.names = TRUE))

# create a pdf with all of the quality profiles prior to filtering #
pre.qual <- matrix(nrow = length(fnFs), ncol = 2)
for(i in 1:nrow(pre.qual)){
  pre.qual[i,1] <- fnFs[i]
  pre.qual[i,2] <- fnRs[i]
}
prequal.list <- c()

# create a pdf with all of the quality profiles prior to filtering #
for(i in 1:nrow(pre.qual)){
  hold <- plotQualityProfile(pre.qual[i,])
  prequal.list[[i]] <- hold
  dev.off()
}
pdf('prequal.pdf')

# Assigning each individual submission by its sample name #
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# placing filtered files in subdirectory #
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filtering #
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread = TRUE)
head(out)
out_diff <- matrix(nrow = nrow(out), ncol = 1)
for(i in 1:nrow(out)){
  out_diff[i,1] <- out[i,1] - out[i,2] 
}
# Learning error rates # 
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Plotting Errors #
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Constructing dada-class objects #
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

# merging denoised forward and reversed reads #
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Constructing Sequence Table #
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# removing chimeras #
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# determine ratio of non-chimeras to all reads #
sum(seqtab.nochim)/sum(seqtab)

# track reads through the pipeline #
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Assign Taxonomy #
taxa <- assignTaxonomy(seqtab.nochim, "~/Cian_Legume16S/rdp_19_toGenus_trainset.fa.gz", multithread = FALSE)
soil.taxa <- taxa

unqs.mock <- seqtab.nochim["ZymoMockDNA",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE)
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
mock.ref <- getSequences("~/Cian_Legume16S/ZymoBIOMICS.STD.refseq.v3/ssrRNAs")
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

samples.out <- rownames(seqtab.nochim)
metadata <- read.csv2("~/Cian_Legume16S/soil_metadata.csv", sep = ',')
rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]

soil.ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(soil.taxa))

soil.dna <- Biostrings::DNAStringSet(taxa_names(soil.ps))
names(dna) <- taxa_names(soil.ps)
soil.ps <- merge_phyloseq(soil.ps, dna)
taxa_names(soil.ps) <- paste0("ASV", seq(ntaxa(soil.ps)))
soil.ps

soil.ps <- subset_taxa(soil.ps, Phylum != "Plantae")
soil.ps <- subset_taxa(soil.ps, Phylum != "Cyanobacteriota")

soil.ps <- subset_samples(soil.ps, Plant != 'Imma')
soil.ps <- subset_samples(soil.ps, Compartment != 'Leaf')
soil.ps <- subset_samples(soil.ps, Plant != 'Lupine')

soil.ps <- subset_taxa(soil.ps, taxa_sums(soil.ps) > 100)

colnames(tax_table(soil.ps)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')


soil.tax <- as.data.frame(tax_table(soil.ps))

for(i in 1:nrow(soil.tax)){
  if(!is.na(soil.tax$Genus[i])){
  rownames(soil.tax)[i] = paste0(rownames(soil.tax)[i], '(', soil.tax$Genus[i], ')')
  } else if(!is.na(soil.tax$Family[i])){
    rownames(soil.tax)[i] = paste0(rownames(soil.tax)[i], '(', soil.tax$Family[i], ')')
  } else if(!is.na(soil.tax$Order[i])){
    rownames(soil.tax)[i] = paste0(rownames(soil.tax)[i], '(', soil.tax$Order[i], ')')
  } else if (!is.na(soil.tax$Class[i])){
    rownames(soil.tax)[i] = paste0(rownames(soil.tax)[i], '(', soil.tax$Class[i], ')')
  } else if (is.na(soil.tax$Phylum[i])){
    rownames(soil.tax)[i] = paste0(rownames(soil.tax)[i], '(', soil.tax$Phylum[i], ')')
  } else{
    rownames(soil.tax)[i] = paste0(rownames(soil.tax)[i], '(NA)')
  }
}

soil.tax <- as.matrix(soil.tax)  
soil.otu <- as.data.frame(otu_table(soil.ps))
rownames(soil.otu) <- rownames(soil.tax)
soil.met <- as(sample_data(soil.ps), 'data.frame')
soil.dna <- refseq(soil.ps)
names(soil.dna) <- rownames(soil.tax)

soil.ps <- phyloseq(otu_table(soil.otu, taxa_are_rows = TRUE),
                    tax_table(soil.tax),
                    sample_data(soil.met),
                    refseq(soil.dna))

nodsoil.ps <- subset_samples(soil.ps, Compartment == 'Nodule')
nodsoil.ps <- subset_taxa(nodsoil.ps, taxa_sums(nodsoil.ps) > 0)

# Constructing MSA and a Phylogenetic Tree #
## See msa.R file ##
soil.ps <- phyloseq(otu_table(soil.otu, taxa_are_rows = TRUE),
                    tax_table(soil.tax),
                    sample_data(soil.met),
                    refseq(soil.dna),
                    phy_tree(soil.gtrfit$tree))

#### DADA2 Pipeline for Root Endosphere Communities ####
end_taxa <- assignTaxonomy(true_endo_seqs, "~/Cian_Legume16S/rdp_19_toGenus_trainset.fa.gz", multithread = FALSE)

endo.met <- read.csv2('~/Cian_Legume16S/endo_metadata.csv', sep = ',')
rownames(endo.met) <- endo.met[,1]
endo.met <- endo.met[,-1]
true_endo_seqs <- t(true_endo_seqs)
endo.otu <- matrix(nrow=nrow(true_endo_seqs), ncol = nrow(endo.met))
endo.otu <- as.data.frame(endo.otu)

count = 1
for(i in 1:ncol(true_endo_seqs)){
  if(colnames(true_endo_seqs)[i] %in% rownames(endo.met)){
    endo.otu[,count] <- true_endo_seqs[,i]
    colnames(endo.otu)[count] <- colnames(true_endo_seqs)[i]
    count = count + 1
  }
}
rownames(endo.otu) <-rownames(end_taxa)
end.ps <- phyloseq(otu_table(endo.otu, taxa_are_rows = TRUE),
                   sample_data(endo.met),
                   tax_table(end_taxa))

end.dna <- Biostrings::DNAStringSet(rownames(endo.otu))
names(end.dna) <- taxa_names(end.ps)
end.ps <- merge_phyloseq(end.ps, end.dna)
taxa_names(end.ps) <- paste0("ASV", seq(ntaxa(end.ps)))
end.ps

end.ps <- subset_taxa(end.ps, Order != 'Chloroplast')

end.ps <- subset_taxa(end.ps, taxa_sums(end.ps) > 100)

taxa_names(end.ps) <- paste0("ASV", seq(ntaxa(end.ps)))
end.tax <- as.data.frame(tax_table(end.ps))

for(i in 1:nrow(end.tax)){
  if(!is.na(end.tax$Genus[i])){
    rownames(end.tax)[i] = paste0(rownames(end.tax)[i], '(', end.tax$Genus[i], ')')
  } else if(!is.na(end.tax$Family[i])){
    rownames(end.tax)[i] = paste0(rownames(end.tax)[i], '(', end.tax$Family[i], ')')
  } else if(!is.na(end.tax$Order[i])){
    rownames(end.tax)[i] = paste0(rownames(end.tax)[i], '(', end.tax$Order[i], ')')
  } else if (!is.na(end.tax$Class[i])){
    rownames(end.tax)[i] = paste0(rownames(end.tax)[i], '(', end.tax$Class[i], ')')
  } else if (is.na(end.tax$Phylum[i])){
    rownames(end.tax)[i] = paste0(rownames(end.tax)[i], '(', end.tax$Phylum[i], ')')
  } else{
    rownames(end.tax)[i] = paste0(rownames(end.tax)[i], '(NA)')
  }
}

end.tax <- as.matrix(end.tax)  
end.otu <- as.data.frame(otu_table(end.ps))
rownames(end.otu) <- rownames(end.tax)
end.met <- as(sample_data(end.ps), 'data.frame')
end.dna <- refseq(end.ps)
names(end.dna) <- rownames(end.tax)

end.ps <- phyloseq(otu_table(end.otu, taxa_are_rows = TRUE),
                    tax_table(end.tax),
                    sample_data(end.met),
                    refseq(end.dna),
                   phy_tree(end.gtrfit$tree))

save(end.dna, file = '~/Cian_Legume16S/endserver.RData')

# Constructing MSA and Phylogenetic Tree #
## See msa.R file ##

#### Across Species Analysis ####
# phyloseq and dataframe of all ASVs grouped by genera #
all.glom <- tax_glom(all.ps, taxrank = "Genus")
all.df <- psmelt(all.glom)
all.df

# making a color pallette that remains the same across figures #
set.seed(4)
all.color <- createPalette(ntaxa(all.glom),  c("#ff0000", "#00ff00", "#0000ff"))
all.color <- as.data.frame(all.color)
rownames(all.color) <- unique(tax_table())

#### Alpha Diversity Measures and Figures ####
### Soils First ###
soil.rich <- estimate_richness(soil.ps)
soil.rich <- as.data.frame(soil.rich)
soil.rich <- cbind(sample_data(soil.ps), soil.rich)

end.rich <- estimate_richness(end.ps)
end.rich <- as.data.frame(end.rich)
end.rich <- cbind(sample_data(end.ps), end.rich)

all.rich <- rbind(soil.rich, end.rich)

for(i in 1:nrow(all.rich)){
  all.rich$ShaEvn[i] <- all.rich$Shannon[i]/log(all.rich$Chao1[i], base = 2.718) 
}

all.rich$`Plant Species` <- factor(all.rich$Plant, levels = c('A. bracteata', 'D. illinoense', 'C. fasciculata', 'S. helvola', 'T. repens', 'M. truncatula'))
all.rich$`Soil Type` <- factor(all.rich$Soil_Treatment, levels = c('Common Soil', 'Non-PSF Soil', 'PSF Soil'))

## ANOVA of Alpha Diversity in Soils ##
soils.rich <- filter(all.rich, Compartment != 'Nodule')
soils.rich <- filter(soils.rich, Compartment != 'Root')
soils.rich$Plants <- factor(soils.rich$Plant)
soils.rich$Soil <- factor(soils.rich$Soil_Treatment)
soils.rich$Comps <- factor(soils.rich$Compartment)

# Tripartite Interaction #
for(i in 1:nrow(soils.rich)){
  soils.rich$PSC[i] <- paste0(soils.rich$Plant[i], ';', soils.rich$Soil_Treatment[i], ';', soils.rich$Compartment[i])
}
soil_psc.sha <- aov(Shannon~PSC, soils.rich)
summary(soil_psc.sha)
soil_psc.evn <- aov(ShaEvn~PSC, soils.rich)
summary(soil_psc.evn)
soil_psc.cha <- aov(Chao1~PSC, soils.rich)
summary(soil_psc.cha)

# Soil x Compartment #
for(i in 1:nrow(soils.rich)){
  soils.rich$SC[i] <- paste0(soils.rich$Soil_Treatment[i], ';', soils.rich$Compartment[i])
}
soil_sc.sha <- aov(Shannon~SC, soils.rich)
summary(soil_sc.sha)
soil_sc.evn <- aov(ShaEvn~SC, soils.rich)
summary(soil_sc.evn)
soil_sc.cha <- aov(Chao1~SC, soils.rich)
summary(soil_sc.cha)

# Plant x Compartment #
for(i in 1:nrow(soils.rich)){
  soils.rich$PC[i] <- paste0(soils.rich$Plant[i], ';', soils.rich$Compartment[i])
}
soil_pc.sha <- aov(Shannon~PC, soils.rich)
summary(soil_pc.sha)
soil_pc.evn <- aov(ShaEvn~PC, soils.rich)
summary(soil_pc.evn)
soil_pc.cha <- aov(Chao1~PC, soils.rich)
summary(soil_pc.cha)

# Plant x Soil #
for(i in 1:nrow(soils.rich)){
  soils.rich$PS[i] <- paste0(soils.rich$Plant[i], ';', soils.rich$Soil_Treatment[i])
}
soil_ps.sha <- aov(Shannon~PS, soils.rich)
summary(soil_ps.sha)
soil_ps.evn <- aov(ShaEvn~PS, soils.rich)
summary(soil_ps.evn)
soil_ps.cha <- aov(Chao1~PS, soils.rich)
summary(soil_ps.cha)

# Compartment #
soil_c.sha <- aov(Shannon~Compartment, soils.rich)
summary(soil_c.sha)
soil_c.evn <- aov(ShaEvn~Compartment, soils.rich)
summary(soil_c.evn)
soil_c.cha <- aov(Chao1~Compartment, soils.rich)
summary(soil_c.cha)

# Soil #
soil_s.sha <- aov(Shannon~Soil_Treatment, soils.rich)
summary(soil_s.sha)
soil_s.evn <- aov(ShaEvn~Soil_Treatment, soils.rich)
summary(soil_s.evn)
soil_s.cha <- aov(Chao1~Soil_Treatment, soils.rich)
summary(soil_s.cha)

# Soil #
soil_p.sha <- aov(Shannon~Plant, soils.rich)
summary(soil_p.sha)
soil_p.evn <- aov(ShaEvn~Plant, soils.rich)
summary(soil_p.evn)
soil_p.cha <- aov(Chao1~Plant, soils.rich)
summary(soil_p.cha)

all.sha <- aov(Shannon~Plant*Soil_Treatment*Compartment, data = all.rich) 
summary(all.sha)
all.evn <- aov(ShaEvn~Plant*Soil_Treatment*Compartment, data = all.rich) 
summary(all.evn)
all.cha <- aov(Chao1~Plant*Soil_Treatment*Compartment, data = all.rich) 
summary(all.cha)

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

# All Rhizosphere Samples #
rhiz.rich <- filter(all_rich.mnsd, Compartment == "Rhizosphere")
rhiz.richraw <- filter(all.rich, Compartment == 'Rhizosphere')

shapiro.test(rhiz.richraw$Shannon)

# Rhizosphere Shannon Diversity #
rhiz.rich <- filter(all_rich.mnsd, Compartment == "Rhizosphere")
rhiz.richraw <- filter(all.rich, Compartment == 'Rhizosphere')
rhiz.richraw$Soil_Treatment <- gsub('Non-PSF Soil', 'Non PSF Soil', rhiz.richraw$Soil_Treatment)

## Strophostyles ##
fb_rhiz.richraw <- filter(rhiz.richraw, Plant == 'S. helvola')
fb_rhiz_sha.aov <- aov(Shannon ~ Soil_Treatment, data = fb_rhiz.richraw)
summary(fb_rhiz_sha.aov)

fb_rhiz_sha.hsd <- TukeyHSD(fb_rhiz_sha.aov)
fb_rhiz_sha.hsd

fb_rhiz_sha.let <- multcompLetters4(fb_rhiz_sha.aov, fb_rhiz_sha.hsd)
fb_rhiz_sha.let

## Chamecrista ##
cc_rhiz.richraw <- filter(rhiz.richraw, Plant == 'C. fasciculata')
cc_rhiz_sha.aov <- aov(Shannon ~ Soil_Treatment, data = cc_rhiz.richraw)
summary(cc_rhiz_sha.aov)

cc_rhiz_sha.hsd <- TukeyHSD(cc_rhiz_sha.aov)
cc_rhiz_sha.hsd

cc_rhiz_sha.let <- multcompLetters4(cc_rhiz_sha.aov, cc_rhiz_sha.hsd)
cc_rhiz_sha.let

## Desmodium ##
ds_rhiz.richraw <- filter(rhiz.richraw, Plant == 'D. illinoense')
ds_rhiz_sha.aov <- aov(Shannon ~ Soil_Treatment, data = ds_rhiz.richraw)
summary(ds_rhiz_sha.aov)

ds_rhiz_sha.hsd <- TukeyHSD(ds_rhiz_sha.aov)
ds_rhiz_sha.hsd

ds_rhiz_sha.let <- multcompLetters4(ds_rhiz_sha.aov, ds_rhiz_sha.hsd)
ds_rhiz_sha.let

## Amphicarpaea ##
hp_rhiz.richraw <- filter(rhiz.richraw, Plant == 'A. bracteata')
hp_rhiz_sha.aov <- aov(Shannon ~ Soil_Treatment, data = hp_rhiz.richraw)
summary(hp_rhiz_sha.aov)

hp_rhiz_sha.hsd <- TukeyHSD(hp_rhiz_sha.aov)
hp_rhiz_sha.hsd

hp_rhiz_sha.let <- multcompLetters4(hp_rhiz_sha.aov, hp_rhiz_sha.hsd)
hp_rhiz_sha.let

## Trifolium ##
cl_rhiz.richraw <- filter(rhiz.richraw, Plant == 'T. repens')
cl_rhiz_sha.aov <- aov(Shannon ~ Soil_Treatment, data = cl_rhiz.richraw)
summary(cl_rhiz_sha.aov)

cl_rhiz_sha.hsd <- TukeyHSD(cl_rhiz_sha.aov)
cl_rhiz_sha.hsd

cl_rhiz_sha.let <- multcompLetters4(cl_rhiz_sha.aov, cl_rhiz_sha.hsd)
cl_rhiz_sha.let

## Medicago ##
md_rhiz.richraw <- filter(rhiz.richraw, Plant == 'M. truncatula')
md_rhiz_sha.aov <- aov(Shannon ~ Soil_Treatment, data = md_rhiz.richraw)
summary(md_rhiz_sha.aov)

md_rhiz_sha.hsd <- TukeyHSD(md_rhiz_sha.aov)
md_rhiz_sha.hsd

md_rhiz_sha.let <- multcompLetters4(md_rhiz_sha.aov, md_rhiz_sha.hsd)
md_rhiz_sha.let

## Adding Letters ##
sha_rhiz.let <- c(hp_rhiz_sha.let$Soil_Treatment$Letters,
                      cc_rhiz_sha.let$Soil_Treatment$Letters,
                      ds_rhiz_sha.let$Soil_Treatment$Letters,
                      md_rhiz_sha.let$Soil_Treatment$Letters,
                      fb_rhiz_sha.let$Soil_Treatment$Letters,
                      cl_rhiz_sha.let$Soil_Treatment$Letters)
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
  scale_y_continuous(limits = c(0,7.5),expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ .,name = "Rhizosphere")) +
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
fb_rhiz_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = fb_rhiz.richraw)
summary(fb_rhiz_evn.aov)

fb_rhiz_evn.hsd <- TukeyHSD(fb_rhiz_evn.aov)
fb_rhiz_evn.hsd

fb_rhiz_evn.let <- multcompLetters4(fb_rhiz_evn.aov, fb_rhiz_evn.hsd)
fb_rhiz_evn.let

## Chamecrista ##
cc_rhiz.richraw <- filter(rhiz.richraw, Plant == 'C. fasciculata')
cc_rhiz_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = cc_rhiz.richraw)
summary(cc_rhiz_evn.aov)

cc_rhiz_evn.hsd <- TukeyHSD(cc_rhiz_evn.aov)
cc_rhiz_evn.hsd

cc_rhiz_evn.let <- multcompLetters4(cc_rhiz_evn.aov, cc_rhiz_evn.hsd)
cc_rhiz_evn.let

## Desmodium ##
ds_rhiz.richraw <- filter(rhiz.richraw, Plant == 'D. illinoense')
ds_rhiz_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = ds_rhiz.richraw)
summary(ds_rhiz_evn.aov)

ds_rhiz_evn.hsd <- TukeyHSD(ds_rhiz_evn.aov)
ds_rhiz_evn.hsd

ds_rhiz_evn.let <- multcompLetters4(ds_rhiz_evn.aov, ds_rhiz_evn.hsd)
ds_rhiz_evn.let

## Amphicarpaea ##
hp_rhiz.richraw <- filter(rhiz.richraw, Plant == 'A. bracteata')
hp_rhiz_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = hp_rhiz.richraw)
summary(hp_rhiz_evn.aov)

hp_rhiz_evn.hsd <- TukeyHSD(hp_rhiz_evn.aov)
hp_rhiz_evn.hsd

hp_rhiz_evn.let <- multcompLetters4(hp_rhiz_evn.aov, hp_rhiz_evn.hsd)
hp_rhiz_evn.let

## Trifolium ##
cl_rhiz.richraw <- filter(rhiz.richraw, Plant == 'T. repens')
cl_rhiz_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = cl_rhiz.richraw)
summary(cl_rhiz_evn.aov)

cl_rhiz_evn.hsd <- TukeyHSD(cl_rhiz_evn.aov)
cl_rhiz_evn.hsd

cl_rhiz_evn.let <- multcompLetters4(cl_rhiz_evn.aov, cl_rhiz_evn.hsd)
cl_rhiz_evn.let

## Medicago ##
md_rhiz.richraw <- filter(rhiz.richraw, Plant == 'M. truncatula')
md_rhiz_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = md_rhiz.richraw)
summary(md_rhiz_evn.aov)

md_rhiz_evn.hsd <- TukeyHSD(md_rhiz_evn.aov)
md_rhiz_evn.hsd

md_rhiz_evn.let <- multcompLetters4(md_rhiz_evn.aov, md_rhiz_evn.hsd)
md_rhiz_evn.let

## Adding Letters ##
evn_rhiz.let <- c(hp_rhiz_evn.let$Soil_Treatment$Letters,
                  cc_rhiz_evn.let$Soil_Treatment$Letters,
                  ds_rhiz_evn.let$Soil_Treatment$Letters,
                  md_rhiz_evn.let$Soil_Treatment$Letters,
                  fb_rhiz_evn.let$Soil_Treatment$Letters,
                  cl_rhiz_evn.let$Soil_Treatment$Letters)
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
fb_rhiz_cha.aov <- aov(Chao1 ~ Soil_Treatment, data = fb_rhiz.richraw)
summary(fb_rhiz_cha.aov)

fb_rhiz_cha.hsd <- TukeyHSD(fb_rhiz_cha.aov)
fb_rhiz_cha.hsd

fb_rhiz_cha.let <- multcompLetters4(fb_rhiz_cha.aov, fb_rhiz_cha.hsd)
fb_rhiz_cha.let

## Chamecrista ##
cc_rhiz.richraw <- filter(rhiz.richraw, Plant == 'C. fasciculata')
cc_rhiz_cha.aov <- aov(Chao1~ Soil_Treatment, data = cc_rhiz.richraw)
summary(cc_rhiz_cha.aov)

cc_rhiz_cha.hsd <- TukeyHSD(cc_rhiz_cha.aov)
cc_rhiz_cha.hsd

cc_rhiz_cha.let <- multcompLetters4(cc_rhiz_cha.aov, cc_rhiz_cha.hsd)
cc_rhiz_cha.let

## Desmodium ##
ds_rhiz.richraw <- filter(rhiz.richraw, Plant == 'D. illinoense')
ds_rhiz_cha.aov <- aov(Chao1 ~ Soil_Treatment, data = ds_rhiz.richraw)
summary(ds_rhiz_cha.aov)

ds_rhiz_cha.hsd <- TukeyHSD(ds_rhiz_cha.aov)
ds_rhiz_cha.hsd

ds_rhiz_cha.let <- multcompLetters4(ds_rhiz_cha.aov, ds_rhiz_cha.hsd)
ds_rhiz_cha.let

## Amphicarpaea ##
hp_rhiz.richraw <- filter(rhiz.richraw, Plant == 'A. bracteata')
hp_rhiz_cha.aov <- aov(Chao1 ~ Soil_Treatment, data = hp_rhiz.richraw)
summary(hp_rhiz_cha.aov)

hp_rhiz_cha.hsd <- TukeyHSD(hp_rhiz_cha.aov)
hp_rhiz_cha.hsd

hp_rhiz_cha.let <- multcompLetters4(hp_rhiz_cha.aov, hp_rhiz_cha.hsd)
hp_rhiz_cha.let

## Trifolium ##
cl_rhiz.richraw <- filter(rhiz.richraw, Plant == 'T. repens')
cl_rhiz_cha.aov <- aov(Chao1 ~ Soil_Treatment, data = cl_rhiz.richraw)
summary(cl_rhiz_cha.aov)

cl_rhiz_cha.hsd <- TukeyHSD(cl_rhiz_cha.aov)
cl_rhiz_cha.hsd

cl_rhiz_cha.let <- multcompLetters4(cl_rhiz_cha.aov, cl_rhiz_cha.hsd)
cl_rhiz_cha.let

## Medicago ##
md_rhiz.richraw <- filter(rhiz.richraw, Plant == 'M. truncatula')
md_rhiz_cha.aov <- aov(Chao1 ~ Soil_Treatment, data = md_rhiz.richraw)
summary(md_rhiz_cha.aov)

md_rhiz_cha.hsd <- TukeyHSD(md_rhiz_cha.aov)
md_rhiz_cha.hsd

md_rhiz_cha.let <- multcompLetters4(md_rhiz_cha.aov, md_rhiz_cha.hsd)
md_rhiz_cha.let

## Adding Letters ##
cha_rhiz.let <- c(hp_rhiz_cha.let$Soil_Treatment$Letters,
                  cc_rhiz_cha.let$Soil_Treatment$Letters,
                  ds_rhiz_cha.let$Soil_Treatment$Letters,
                  md_rhiz_cha.let$Soil_Treatment$Letters,
                  fb_rhiz_cha.let$Soil_Treatment$Letters,
                  cl_rhiz_cha.let$Soil_Treatment$Letters)
rhiz.rich <- cbind(rhiz.rich, cha_rhiz.let)
rhiz.rich$cha_rhiz.let[13] <- 'ab'
rhiz.rich$cha_rhiz.let[14] <- 'b'
rhiz.rich$cha_rhiz.let[15] <- 'a'
rhiz.rich$cha_rhiz.let[5] <- 'b'
rhiz.rich$cha_rhiz.let[6] <- 'ab'



### Final Plot ###
cha_rhiz.plot <- ggplot(rhiz.rich, aes(x = `Plant Species`, y = cha.mean, fill = `Soil Type`, color = `Soil Type`)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = cha.mean - cha.sd, ymax = cha.mean + cha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = cha_rhiz.let, y = cha.mean + cha.sd + 25), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 4) +
  ylab('Observed ASV Richness (S)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,1700),expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = "")) +
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
bulk.richraw$Soil_Treatment <- gsub('Non-PSF Soil', 'Non PSF Soil', bulk.richraw$Soil_Treatment)
shapiro.test(bulk.richraw$Shannon)

bulk.rich[11,] <- c('A. bracteata', 'Non-PSF Soil', 'Bulk Soil', 6.38695563671935, 0.916637445922254, 1064.33333333333, 0.090210043413787, 0.0014385357114979, 99.5707453689754, 'A. bracteata', 'Non-PSF Soil')
bulk.rich[12,] <- c('S. helvola', 'Non-PSF Soil', 'Bulk Soil', 5.81122621101946, 0.860055149110637, 861.333333333333, 0.0885416648141404,  0.0134357294023773, 68.2959247198055, 'S. helvola', 'Non-PSF Soil')
bulk.rich[13,] <- c('C. fasciculata', 'Common Soil', 'Bulk Soil', 6.20898577865338, 0.909855472564764, 920.333333333333, 0.0692706178522534, 0.00707817919502481, 56.0832714214616, 'C. fasciculata', 'Common Soil')
bulk.rich[14,] <- c('D. illinoense', 'Common Soil', 'Bulk Soil', 6.20898577865338, 0.909855472564764, 920.333333333333, 0.0692706178522534, 0.00707817919502481, 56.0832714214616, 'D. illinoense', 'Common Soil')
bulk.rich[15,] <- c('A. bracteata', 'Common Soil', 'Bulk Soil', 6.20898577865338, 0.909855472564764, 920.333333333333, 0.0692706178522534, 0.00707817919502481, 56.0832714214616, 'A. bracteata', 'Common Soil')
bulk.rich[16,] <- c('T. repens', 'Common Soil', 'Bulk Soil', 6.20898577865338, 0.909855472564764, 920.333333333333, 0.0692706178522534, 0.00707817919502481, 56.0832714214616, 'T. repens', 'Common Soil')
bulk.rich[17,] <- c('M. truncatula', 'Common Soil', 'Bulk Soil', 6.20898577865338, 0.909855472564764, 920.333333333333, 0.0692706178522534, 0.00707817919502481, 56.0832714214616, 'M. truncatula', 'Common Soil')
bulk.rich[18,] <- c('M. truncatula', 'Common Soil', 'Bulk Soil', 5.95648139950243, 0.898640165840898, 814.333333333333,  0.388295980568504, 0.00896487023256766, 333.008007911722,'M. truncatula', 'Non-PSF Soil')

bulk.rich$sha.mean <- as.numeric(bulk.rich$sha.mean)
bulk.rich$sha.sd <- as.numeric(bulk.rich$sha.sd)
bulk.rich$evn.mean <- as.numeric(bulk.rich$evn.mean)
bulk.rich$evn.sd <- as.numeric(bulk.rich$evn.sd)
bulk.rich$cha.mean <- as.numeric(bulk.rich$cha.mean)
bulk.rich$cha.sd <- as.numeric(bulk.rich$cha.sd)

bulk.richraw$Chao1 <- as.numeric(bulk.richraw$Chao1) 
bulk.richraw$Shannon <- as.numeric(bulk.richraw$Shannon)
bulk.richraw$ShaEvn <- as.numeric(bulk.richraw$ShaEvn)

## Strophostyles ##
fb_bulk.richraw <- filter(bulk.richraw, Plant == 'S. helvola')
fb_bulk.richraw[7,] <- c(5076, 'S. helvola', 'Non PSF Soil', 'Bulk Soil', 938, 938, 0, 938, 14.39472, 5.870144, 0.9917414, 121.08588, 172.2512, 0.8576490, 'S. helvola', 'Non-PSF Soil')
fb_bulk.richraw[8,] <- c(5077, 'S. helvola', 'Non PSF Soil', 'Bulk Soil', 839, 839, 0, 839, 13.98995, 5.709406, 0.9891960, 92.55843, 160.3916, 0.8479851, 'S. helvola', 'Non-PSF Soil')
fb_bulk.richraw[9,] <- c(5078, 'S. helvola', 'Non PSF Soil', 'Bulk Soil', 807, 807, 0, 807, 13.24827, 5.854129, 0.9931115, 145.16855, 153.1863, 0.8745314, 'S. helvola', 'Non-PSF Soil')

fb_bulk.richraw$Chao1 <- as.numeric(fb_bulk.richraw$Chao1) 
fb_bulk.richraw$Shannon <- as.numeric(fb_bulk.richraw$Shannon)
fb_bulk.richraw$ShaEvn <- as.numeric(fb_bulk.richraw$ShaEvn)

fb_bulk.rich <- filter(bulk.rich, Plant == 'S. helvola')
fb_bulk_sha.aov <- aov(Shannon ~ Soil_Treatment, data = fb_bulk.richraw)
summary(fb_bulk_sha.aov)

fb_bulk_sha.hsd <- TukeyHSD(fb_bulk_sha.aov)
fb_bulk_sha.hsd

fb_bulk_sha.let <- multcompLetters4(fb_bulk_sha.aov, fb_bulk_sha.hsd)
fb_bulk_sha.let

## Chamecrista ##
cc_bulk.richraw <- filter(bulk.richraw, Plant == 'C. fasciculata')
cc_bulk.richraw[7,] <- fb_bulk.richraw[4,]
cc_bulk.richraw[8,] <- fb_bulk.richraw[5,]
cc_bulk.richraw[9,] <- fb_bulk.richraw[6,]

cc_bulk.richraw$Chao1 <- as.numeric(cc_bulk.richraw$Chao1) 
cc_bulk.richraw$Shannon <- as.numeric(cc_bulk.richraw$Shannon)
cc_bulk.richraw$ShaEvn <- as.numeric(cc_bulk.richraw$ShaEvn)

cc_bulk_sha.aov <- aov(Shannon ~ Soil_Treatment, data = cc_bulk.richraw)
summary(cc_bulk_sha.aov)

cc_bulk_sha.hsd <- TukeyHSD(cc_bulk_sha.aov)
cc_bulk_sha.hsd

cc_bulk_sha.let <- multcompLetters4(cc_bulk_sha.aov, cc_bulk_sha.hsd)
cc_bulk_sha.let

## Desmodium ##
ds_bulk.richraw <- filter(bulk.richraw, Plant == 'D. illinoense')
ds_bulk.richraw[7,] <- fb_bulk.richraw[4,]
ds_bulk.richraw[8,] <- fb_bulk.richraw[5,]
ds_bulk.richraw[9,] <- fb_bulk.richraw[6,]

ds_bulk.richraw$Chao1 <- as.numeric(ds_bulk.richraw$Chao1) 
ds_bulk.richraw$Shannon <- as.numeric(ds_bulk.richraw$Shannon)
ds_bulk.richraw$ShaEvn <- as.numeric(ds_bulk.richraw$ShaEvn)

ds_bulk_sha.aov <- aov(Shannon ~ Soil_Treatment, data = ds_bulk.richraw)
summary(ds_bulk_sha.aov)

ds_bulk_sha.hsd <- TukeyHSD(ds_bulk_sha.aov)
ds_bulk_sha.hsd

ds_bulk_sha.let <- multcompLetters4(ds_bulk_sha.aov, ds_bulk_sha.hsd)
ds_bulk_sha.let

## Amphicarpaea ##
hp_bulk.richraw <- filter(bulk.richraw, Plant == 'A. bracteata')
hp_bulk.richraw[4,] <- ds_bulk.richraw[4,]
hp_bulk.richraw[5,] <- ds_bulk.richraw[5,]
hp_bulk.richraw[6,] <- ds_bulk.richraw[6,]
hp_bulk.richraw[7,] <- fb_bulk.richraw[4,]
hp_bulk.richraw[8,] <- fb_bulk.richraw[5,]
hp_bulk.richraw[9,] <- fb_bulk.richraw[6,]

hp_bulk.richraw$Chao1 <- as.numeric(hp_bulk.richraw$Chao1) 
hp_bulk.richraw$Shannon <- as.numeric(hp_bulk.richraw$Shannon)
hp_bulk.richraw$ShaEvn <- as.numeric(hp_bulk.richraw$ShaEvn)

hp_bulk_sha.aov <- aov(Shannon ~ Soil_Treatment, data = hp_bulk.richraw)
summary(hp_bulk_sha.aov)

hp_bulk_sha.hsd <- TukeyHSD(hp_bulk_sha.aov)
hp_bulk_sha.hsd

hp_bulk_sha.let <- multcompLetters4(hp_bulk_sha.aov, hp_bulk_sha.hsd)
hp_bulk_sha.let

## Trifolium ##
cl_bulk.richraw <- filter(bulk.richraw, Plant == 'T. repens')
cl_bulk.richraw[7,] <- fb_bulk.richraw[4,]
cl_bulk.richraw[8,] <- fb_bulk.richraw[5,]
cl_bulk.richraw[9,] <- fb_bulk.richraw[6,]

cl_bulk.richraw$Chao1 <- as.numeric(cl_bulk.richraw$Chao1) 
cl_bulk.richraw$Shannon <- as.numeric(cl_bulk.richraw$Shannon)
cl_bulk.richraw$ShaEvn <- as.numeric(cl_bulk.richraw$ShaEvn)

cl_bulk_sha.aov <- aov(Shannon ~ Soil_Treatment, data = cl_bulk.richraw)
summary(cl_bulk_sha.aov)

cl_bulk_sha.hsd <- TukeyHSD(cl_bulk_sha.aov)
cl_bulk_sha.hsd

cl_bulk_sha.let <- multcompLetters4(cl_bulk_sha.aov, cl_bulk_sha.hsd)
cl_bulk_sha.let

## Medicago ##
md_bulk.richraw <- filter(bulk.richraw, Plant == 'M. truncatula')
md_bulk.richraw[4,] <- cl_bulk.richraw[4,]
md_bulk.richraw[5,] <- cl_bulk.richraw[5,]
md_bulk.richraw[6,] <- cl_bulk.richraw[6,]
md_bulk.richraw[7,] <- fb_bulk.richraw[4,]
md_bulk.richraw[8,] <- fb_bulk.richraw[5,]
md_bulk.richraw[9,] <- fb_bulk.richraw[6,]

md_bulk.richraw$Chao1 <- as.numeric(md_bulk.richraw$Chao1) 
md_bulk.richraw$Shannon <- as.numeric(md_bulk.richraw$Shannon)
md_bulk.richraw$ShaEvn <- as.numeric(md_bulk.richraw$ShaEvn)

md_bulk_sha.aov <- aov(Shannon ~ Soil_Treatment, data = md_bulk.richraw)
summary(md_bulk_sha.aov)

md_bulk_sha.hsd <- TukeyHSD(md_bulk_sha.aov)
md_bulk_sha.hsd

md_bulk_sha.let <- multcompLetters4(md_bulk_sha.aov, md_bulk_sha.hsd)
md_bulk_sha.let

## Adding Letters ##
sha_bulk.let <- c(hp_bulk_sha.let$Soil_Treatment$Letters,
                  cc_bulk_sha.let$Soil_Treatment$Letters,
                  ds_bulk_sha.let$Soil_Treatment$Letters,
                  md_bulk_sha.let$Soil_Treatment$Letters,
                  fb_bulk_sha.let$Soil_Treatment$Letters,
                  cl_bulk_sha.let$Soil_Treatment$Letters)
bulk.rich <- cbind(bulk.rich, sha_bulk.let)
for(i in 1:nrow(bulk.rich)){
  bulk.rich$sha_bulk.let[i] <- 'a'
}

bulk.rich$sha_bulk.let[12] <- 'b'
bulk.rich$sha_bulk.let[2] <- 'b'

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
fb_bulk_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = fb_bulk.richraw)
summary(fb_bulk_evn.aov)

fb_bulk_evn.hsd <- TukeyHSD(fb_bulk_evn.aov)
fb_bulk_evn.hsd

fb_bulk_evn.let <- multcompLetters4(fb_bulk_evn.aov, fb_bulk_evn.hsd)
fb_bulk_evn.let

## Chamecrista ##
cc_bulk_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = cc_bulk.richraw)
summary(cc_bulk_evn.aov)

cc_bulk_evn.hsd <- TukeyHSD(cc_bulk_evn.aov)
cc_bulk_evn.hsd

cc_bulk_evn.let <- multcompLetters4(cc_bulk_evn.aov, cc_bulk_evn.hsd)
cc_bulk_evn.let

## Desmodium ##
ds_bulk_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = ds_bulk.richraw)
summary(ds_bulk_evn.aov)

ds_bulk_evn.hsd <- TukeyHSD(ds_bulk_evn.aov)
ds_bulk_evn.hsd

ds_bulk_evn.let <- multcompLetters4(ds_bulk_evn.aov, ds_bulk_evn.hsd)
ds_bulk_evn.let

## Amphicarpaea ##
hp_bulk_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = hp_bulk.richraw)
summary(hp_bulk_evn.aov)

hp_bulk_evn.hsd <- TukeyHSD(hp_bulk_evn.aov)
hp_bulk_evn.hsd

hp_bulk_evn.let <- multcompLetters4(hp_bulk_evn.aov, hp_bulk_evn.hsd)
hp_bulk_evn.let

## Trifolium ##
cl_bulk_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = cl_bulk.richraw)
summary(cl_bulk_evn.aov)

cl_bulk_evn.hsd <- TukeyHSD(cl_bulk_evn.aov)
cl_bulk_evn.hsd

cl_bulk_evn.let <- multcompLetters4(cl_bulk_evn.aov, cl_bulk_evn.hsd)
cl_bulk_evn.let

## Medicago ##
md_bulk_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = md_bulk.richraw)
summary(md_bulk_evn.aov)

md_bulk_evn.hsd <- TukeyHSD(md_bulk_evn.aov)
md_bulk_evn.hsd

md_bulk_evn.let <- multcompLetters4(md_bulk_evn.aov, md_bulk_evn.hsd)
md_bulk_evn.let

## Adding Letters ##
evn_bulk.let <- c(hp_bulk_evn.let$Soil_Treatment$Letters,
                  cc_bulk_evn.let$Soil_Treatment$Letters,
                  ds_bulk_evn.let$Soil_Treatment$Letters,
                  md_bulk_evn.let$Soil_Treatment$Letters,
                  fb_bulk_evn.let$Soil_Treatment$Letters,
                  cl_bulk_evn.let$Soil_Treatment$Letters)
bulk.rich <- cbind(bulk.rich, evn_bulk.let)
for(i in 1:nrow(bulk.rich)){
  bulk.rich$evn_bulk.let[i] <- 'a'
}

bulk.rich$evn_bulk.let[12] <- 'b'
bulk.rich$evn_bulk.let[2] <- 'b'
bulk.rich$evn_bulk.let[5] <- 'b'


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
fb_bulk_cha.aov <- aov(Chao1 ~ Soil_Treatment, data = fb_bulk.richraw)
summary(fb_bulk_cha.aov)

fb_bulk_cha.hsd <- TukeyHSD(fb_bulk_cha.aov)
fb_bulk_cha.hsd

fb_bulk_cha.let <- multcompLetters4(fb_bulk_cha.aov, fb_bulk_cha.hsd)
fb_bulk_cha.let

## Chamecrista ##
cc_bulk_cha.aov <- aov(Chao1~ Soil_Treatment, data = cc_bulk.richraw)
summary(cc_bulk_cha.aov)

cc_bulk_cha.hsd <- TukeyHSD(cc_bulk_cha.aov)
cc_bulk_cha.hsd

cc_bulk_cha.let <- multcompLetters4(cc_bulk_cha.aov, cc_bulk_cha.hsd)
cc_bulk_cha.let

## Desmodium ##
ds_bulk_cha.aov <- aov(Chao1 ~ Soil_Treatment, data = ds_bulk.richraw)
summary(ds_bulk_cha.aov)

ds_bulk_cha.hsd <- TukeyHSD(ds_bulk_cha.aov)
ds_bulk_cha.hsd

ds_bulk_cha.let <- multcompLetters4(ds_bulk_cha.aov, ds_bulk_cha.hsd)
ds_bulk_cha.let

## Amphicarpaea ##
hp_bulk_cha.aov <- aov(Chao1 ~ Soil_Treatment, data = hp_bulk.richraw)
summary(hp_bulk_cha.aov)

hp_bulk_cha.hsd <- TukeyHSD(hp_bulk_cha.aov)
hp_bulk_cha.hsd

hp_bulk_cha.let <- multcompLetters4(hp_bulk_cha.aov, hp_bulk_cha.hsd)
hp_bulk_cha.let

## Trifolium ##
cl_bulk_cha.aov <- aov(Chao1 ~ Soil_Treatment, data = cl_bulk.richraw)
summary(cl_bulk_cha.aov)

cl_bulk_cha.hsd <- TukeyHSD(cl_bulk_cha.aov)
cl_bulk_cha.hsd

cl_bulk_cha.let <- multcompLetters4(cl_bulk_cha.aov, cl_bulk_cha.hsd)
cl_bulk_cha.let

## Medicago ##
md_bulk_cha.aov <- aov(Chao1 ~ Soil_Treatment, data = md_bulk.richraw)
summary(md_bulk_cha.aov)

md_bulk_cha.hsd <- TukeyHSD(md_bulk_cha.aov)
md_bulk_cha.hsd

md_bulk_cha.let <- multcompLetters4(md_bulk_cha.aov, md_bulk_cha.hsd)
md_bulk_cha.let

## Adding Letters ##
cha_bulk.let <- c(hp_bulk_cha.let$Soil_Treatment$Letters,
                  cc_bulk_cha.let$Soil_Treatment$Letters,
                  ds_bulk_cha.let$Soil_Treatment$Letters,
                  md_bulk_cha.let$Soil_Treatment$Letters,
                  fb_bulk_cha.let$Soil_Treatment$Letters,
                  cl_bulk_cha.let$Soil_Treatment$Letters)
bulk.rich <- cbind(bulk.rich, cha_bulk.let)

for(i in 1:nrow(bulk.rich)){
  bulk.rich$cha_bulk.let[i] <- 'a'
}

bulk.rich$cha_bulk.let[] <- 'b'

### Final Plot ###
cha_bulk.plot <- ggplot(bulk.rich, aes(x = `Plant Species`, y = cha.mean, fill = `Soil Type`, color = `Soil Type`)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = cha.mean - cha.sd, ymax = cha.mean + cha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = cha_bulk.let, y = cha.mean + cha.sd + 25), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 4) +
  ylab('Observed ASV Richness (S)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,1300),expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = "")) +
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
root.rich <- filter(all_rich.mnsd, Compartment == "Root")
root.richraw <- filter(all.rich, Compartment == 'Root')
root.richraw$Soil_Treatment <- gsub('Non-PSF Soil', 'Non PSF Soil', root.richraw$Soil_Treatment)
shapiro.test(root.richraw$Shannon)

# Root Shannon Diversity #

## Strophostyles ##
fb_root.richraw <- filter(root.richraw, Plant == 'S. helvola')
fb_root_sha.aov <- aov(Shannon ~ Soil_Treatment, data = fb_root.richraw)
summary(fb_root_sha.aov)

fb_root_sha.hsd <- TukeyHSD(fb_root_sha.aov)
fb_root_sha.hsd

fb_root_sha.let <- multcompLetters4(fb_root_sha.aov, fb_root_sha.hsd)
fb_root_sha.let

## Chamecrista ##
cc_root.richraw <- filter(root.richraw, Plant == 'C. fasciculata')
cc_root_sha.aov <- aov(Shannon ~ Soil_Treatment, data = cc_root.richraw)
summary(cc_root_sha.aov)

cc_root_sha.hsd <- TukeyHSD(cc_root_sha.aov)
cc_root_sha.hsd

cc_root_sha.let <- multcompLetters4(cc_root_sha.aov, cc_root_sha.hsd)
cc_root_sha.let

## Desmodium ##
ds_root.richraw <- filter(root.richraw, Plant == 'D. illinoense')
ds_root_sha.aov <- aov(Shannon ~ Soil_Treatment, data = ds_root.richraw)
summary(ds_root_sha.aov)

ds_root_sha.hsd <- TukeyHSD(ds_root_sha.aov)
ds_root_sha.hsd

ds_root_sha.let <- multcompLetters4(ds_root_sha.aov, ds_root_sha.hsd)
ds_root_sha.let

## Amphicarpaea ##
hp_root.richraw <- filter(root.richraw, Plant == 'A. bracteata')
hp_root_sha.aov <- aov(Shannon ~ Soil_Treatment, data = hp_root.richraw)
summary(hp_root_sha.aov)

hp_root_sha.hsd <- TukeyHSD(hp_root_sha.aov)
hp_root_sha.hsd

hp_root_sha.let <- multcompLetters4(hp_root_sha.aov, hp_root_sha.hsd)
hp_root_sha.let

## Trifolium ##
cl_root.richraw <- filter(root.richraw, Plant == 'T. repens')
cl_root_sha.aov <- aov(Shannon ~ Soil_Treatment, data = cl_root.richraw)
summary(cl_root_sha.aov)

cl_root_sha.hsd <- TukeyHSD(cl_root_sha.aov)
cl_root_sha.hsd

cl_root_sha.let <- multcompLetters4(cl_root_sha.aov, cl_root_sha.hsd)
cl_root_sha.let

## Medicago ##
md_root.richraw <- filter(root.richraw, Plant == 'M. truncatula')
md_root_sha.aov <- aov(Shannon ~ Soil_Treatment, data = md_root.richraw)
summary(md_root_sha.aov)

md_root_sha.hsd <- TukeyHSD(md_root_sha.aov)
md_root_sha.hsd

md_root_sha.let <- multcompLetters4(md_root_sha.aov, md_root_sha.hsd)
md_root_sha.let

## Adding Letters ##
sha_root.let <- c(hp_root_sha.let$Soil_Treatment$Letters,
                  cc_root_sha.let$Soil_Treatment$Letters,
                  ds_root_sha.let$Soil_Treatment$Letters,
                  md_root_sha.let$Soil_Treatment$Letters,
                  fb_root_sha.let$Soil_Treatment$Letters,
                  cl_root_sha.let$Soil_Treatment$Letters)
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
  scale_y_continuous(limits = c(0,5.5), expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = 'Root Endosphere')) +
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
# Rhizosphere Shannon Evenness #
## Strophostyles ##
fb_root.richraw <- filter(root.richraw, Plant == 'S. helvola')
fb_root_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = fb_root.richraw)
summary(fb_root_evn.aov)

fb_root_evn.hsd <- TukeyHSD(fb_root_evn.aov)
fb_root_evn.hsd

fb_root_evn.let <- multcompLetters4(fb_root_evn.aov, fb_root_evn.hsd)
fb_root_evn.let

## Chamecrista ##
cc_root.richraw <- filter(root.richraw, Plant == 'C. fasciculata')
cc_root_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = cc_root.richraw)
summary(cc_root_evn.aov)

cc_root_evn.hsd <- TukeyHSD(cc_root_evn.aov)
cc_root_evn.hsd

cc_root_evn.let <- multcompLetters4(cc_root_evn.aov, cc_root_evn.hsd)
cc_root_evn.let

## Desmodium ##
ds_root.richraw <- filter(root.richraw, Plant == 'D. illinoense')
ds_root_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = ds_root.richraw)
summary(ds_root_evn.aov)

ds_root_evn.hsd <- TukeyHSD(ds_root_evn.aov)
ds_root_evn.hsd

ds_root_evn.let <- multcompLetters4(ds_root_evn.aov, ds_root_evn.hsd)
ds_root_evn.let

## Amphicarpaea ##
hp_root.richraw <- filter(root.richraw, Plant == 'A. bracteata')
hp_root_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = hp_root.richraw)
summary(hp_root_evn.aov)

hp_root_evn.hsd <- TukeyHSD(hp_root_evn.aov)
hp_root_evn.hsd

hp_root_evn.let <- multcompLetters4(hp_root_evn.aov, hp_root_evn.hsd)
hp_root_evn.let

## Trifolium ##
cl_root.richraw <- filter(root.richraw, Plant == 'T. repens')
cl_root_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = cl_root.richraw)
summary(cl_root_evn.aov)

cl_root_evn.hsd <- TukeyHSD(cl_root_evn.aov)
cl_root_evn.hsd

cl_root_evn.let <- multcompLetters4(cl_root_evn.aov, cl_root_evn.hsd)
cl_root_evn.let

## Medicago ##
md_root.richraw <- filter(root.richraw, Plant == 'M. truncatula')
md_root_evn.aov <- aov(ShaEvn ~ Soil_Treatment, data = md_root.richraw)
summary(md_root_evn.aov)

md_root_evn.hsd <- TukeyHSD(md_root_evn.aov)
md_root_evn.hsd

md_root_evn.let <- multcompLetters4(md_root_evn.aov, md_root_evn.hsd)
md_root_evn.let

## Adding Letters ##
evn_root.let <- c(hp_root_evn.let$Soil_Treatment$Letters,
                  cc_root_evn.let$Soil_Treatment$Letters,
                  ds_root_evn.let$Soil_Treatment$Letters,
                  md_root_evn.let$Soil_Treatment$Letters,
                  fb_root_evn.let$Soil_Treatment$Letters,
                  cl_root_evn.let$Soil_Treatment$Letters)
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
  scale_y_continuous(limits = c(0,0.85),expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = "")) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24),
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
fb_root_cha.aov <- aov(Chao1 ~ Soil_Treatment, data = fb_root.richraw)
summary(fb_root_cha.aov)

fb_root_cha.hsd <- TukeyHSD(fb_root_cha.aov)
fb_root_cha.hsd

fb_root_cha.let <- multcompLetters4(fb_root_cha.aov, fb_root_cha.hsd)
fb_root_cha.let

## Chamecrista ##
cc_root.richraw <- filter(root.richraw, Plant == 'C. fasciculata')
cc_root_cha.aov <- aov(Chao1~ Soil_Treatment, data = cc_root.richraw)
summary(cc_root_cha.aov)

cc_root_cha.hsd <- TukeyHSD(cc_root_cha.aov)
cc_root_cha.hsd

cc_root_cha.let <- multcompLetters4(cc_root_cha.aov, cc_root_cha.hsd)
cc_root_cha.let

## Desmodium ##
ds_root.richraw <- filter(root.richraw, Plant == 'D. illinoense')
ds_root_cha.aov <- aov(Chao1 ~ Soil_Treatment, data = ds_root.richraw)
summary(ds_root_cha.aov)

ds_root_cha.hsd <- TukeyHSD(ds_root_cha.aov)
ds_root_cha.hsd

ds_root_cha.let <- multcompLetters4(ds_root_cha.aov, ds_root_cha.hsd)
ds_root_cha.let

## Amphicarpaea ##
hp_root.richraw <- filter(root.richraw, Plant == 'A. bracteata')
hp_root_cha.aov <- aov(Chao1 ~ Soil_Treatment, data = hp_root.richraw)
summary(hp_root_cha.aov)

hp_root_cha.hsd <- TukeyHSD(hp_root_cha.aov)
hp_root_cha.hsd

hp_root_cha.let <- multcompLetters4(hp_root_cha.aov, hp_root_cha.hsd)
hp_root_cha.let

## Trifolium ##
cl_root.richraw <- filter(root.richraw, Plant == 'T. repens')
cl_root_cha.aov <- aov(Chao1 ~ Soil_Treatment, data = cl_root.richraw)
summary(cl_root_cha.aov)

cl_root_cha.hsd <- TukeyHSD(cl_root_cha.aov)
cl_root_cha.hsd

cl_root_cha.let <- multcompLetters4(cl_root_cha.aov, cl_root_cha.hsd)
cl_root_cha.let

## Medicago ##
md_root.richraw <- filter(root.richraw, Plant == 'M. truncatula')
md_root_cha.aov <- aov(Chao1 ~ Soil_Treatment, data = md_root.richraw)
summary(md_root_cha.aov)

md_root_cha.hsd <- TukeyHSD(md_root_cha.aov)
md_root_cha.hsd

md_root_cha.let <- multcompLetters4(md_root_cha.aov, md_root_cha.hsd)
md_root_cha.let

## Adding Letters ##
cha_root.let <- c(hp_root_cha.let$Soil_Treatment$Letters,
                  cc_root_cha.let$Soil_Treatment$Letters,
                  ds_root_cha.let$Soil_Treatment$Letters,
                  md_root_cha.let$Soil_Treatment$Letters,
                  fb_root_cha.let$Soil_Treatment$Letters,
                  cl_root_cha.let$Soil_Treatment$Letters)
root.rich <- cbind(root.rich, cha_root.let)
root.rich$cha_root.let[1] <- 'ab'
root.rich$cha_root.let[2] <- 'b'
root.rich$cha_root.let[3] <- 'a'

### Final Plot ###
cha_root.plot <- ggplot(root.rich, aes(x = `Plant Species`, y = cha.mean, fill = `Soil Type`, color = `Soil Type`)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = cha.mean - cha.sd, ymax = cha.mean + cha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = cha_root.let, y = cha.mean + cha.sd + 25), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 4) +
  ylab('Observed ASV Richness (S)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(values = c('white', "gray", '#4D4D4D')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,1700), expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ ., name = "")) +
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
(cha_root.plot | evn_root.plot | sha_root.plot)

#### Beta Diversity Measure and Figures ####
### Generating Weighted Unifrac Matrix and PCoA ordinations for NMDS ###
## Soils First ##
soil_nono.ps <- subset_samples(soil.ps, Compartment != 'Nodule') 
soil_nono.ps <- subset_taxa(soil_nono.ps, taxa_sums(soil_nono.ps) > 0)
prop_soil.ps <- transform_sample_counts(soil_nono.ps, function(x) x/sum(x))

set.seed(248)
soil.wuni <- phyloseq::distance(prop_soil.ps, method = 'wunifrac')
soil.pcoa <- phyloseq::ordinate(prop_soil.ps, 'PCoA', distance = 'bray')

### Performing the NMDS ###
## Play around with the k values to get the lowest K that explains => 99% of variation in the linear model ##
library(vegan); packageVersion('vegan')
soil.nmds <- metaMDS(soil.wuni, 
                     k = 7, try = 20, trymax = 1000, maxit = 999,
                     model = 'global', 
                     autotransform = FALSE, previous.best = soil.pcoa$vectors[,1:7])
soil_nmds.scores <- scores(soil.nmds, display = 'sites')
soil_nmds.dist <- dist(soil_nmds.scores)
### Calculating the Variance of the Principal Components (or NMDS axes in our case) ###
## Total Model Variance Calculation ##
soil_nmds.ffit <- lm(as.vector(soil_nmds.dist) ~ as.vector(soil.wuni))
soil_nmds.totr2 <- summary(soil_nmds.ffit)$r.squared

# Axes Variance Calculation #
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

soil_nmds.dist6 <- dist(soil_nmds.scores[,6])
soil_nmds.fit6 <- lm(as.vector(soil_nmds.dist6)~as.vector(soil.wuni))
soil_nmds.r6 <- summary(soil_nmds.fit6)$r.squared

soil_nmds.dist7 <- dist(soil_nmds.scores[,7])
soil_nmds.fit7 <- lm(as.vector(soil_nmds.dist7)~as.vector(soil.wuni))
soil_nmds.r7 <- summary(soil_nmds.fit7)$r.squared

soil_nmds.comb <- soil_nmds.r1 + soil_nmds.r2 + soil_nmds.r3 + soil_nmds.r4 + soil_nmds.r5 + soil_nmds.r6 + soil_nmds.r7

soil_nmds.axisr <- c()
for(i in 1:ncol(soil_nmds.scores)){
  soil_nmds.axisr[i] <- (get(paste0('soil_nmds.r', i)) / soil_nmds.comb) * soil_nmds.totr2 
}

### Calculating Variance Components ###
library(lme4); packageVersion('lme4')

## Constructing a data.frame that has sample info and their loading scores ##
soil_nono.met <- as(sample_data(soil_nono.ps), 'data.frame')
soil_nono.met$Plants <- factor(soil_nono.met$Plant)
soil_nono.met$Comps <- factor(soil_nono.met$Compartment)
soil_nono.met$Soils <- factor(soil_nono.met$Soil_Treatment)
soil_nmds.load <- cbind(soil_nono.met, soil_nmds.scores)

## Testing the mixed linear model on the first NMDS axis
soil_nmds.vfit1 <- lmer(NMDS1 ~  (1|Plants) + (1|Comps) + (1|Soils) +(1|Plants:Comps) + (1|Plants:Soils) + (1|Comps:Soils) + (1|Comps:Plants:Soils), data = soil_nmds.load)
summary(soil_nmds.vfit1)
soil_nmds.vca1 <- as.data.frame(VarCorr(soil_nmds.vfit1))
## Using Loop to do each NMDS axis ##
soil_nmds.vca <- matrix(nrow = 8, ncol = ncol(soil_nmds.scores))
hold <- c()
for(i in 1:ncol(soil_nmds.scores)){
  hold <- lmer(soil_nmds.scores[,i] ~ (1|Plants) + (1|Comps) + (1|Soils) +(1|Plants:Comps) + (1|Plants:Soils) + (1|Comps:Soils) + (1|Comps:Plants:Soils), data = soil_nmds.load)
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


rownames(soil_nmds.vca) <- c('Compartment x Plant x Soil', 'Plant x Soil', 'Plant x Compartment', 'Soil x Compartment', 'Plant', 'Soil', 'Compartment', 'Residuals')
colnames(soil_nmds.vca) <- colnames(soil_nmds.scores)
soil_nmds.vtot <- colSums(soil_nmds.vca)

### Putting original variance and variance components together for PVCA ###
## weighting each variance component by the variance explained by each NMDS axis ##
soil_nmds.wvca <- matrix(nrow = nrow(soil_nmds.vca), ncol = length(soil_nmds.axisr))
for(i in 1:length(soil_nmds.axisr)){
  for(j in 1:nrow(soil_nmds.vca)){
    soil_nmds.wvca[j,i] <- soil_nmds.vca[j,i]*soil_nmds.axisr[i] 
  }
}

rownames(soil_nmds.wvca) <- rownames(soil_nmds.vca); colnames(soil_nmds.wvca) <- colnames(soil_nmds.vca)
soil_nmds.tvca <- rowSums(soil_nmds.wvca)
soil_nmds.tvca <- as.data.frame(soil_nmds.tvca)
soil_nmds.ptot <- colSums(soil_nmds.tvca)

soil_nmds.pvca <- matrix(nrow = nrow(soil_nmds.tvca), ncol = 1)
for(i in 1:nrow(soil_nmds.tvca)){
  soil_nmds.pvca[i,1] <- soil_nmds.tvca[i,1] / soil_nmds.ptot * 100
}

rownames(soil_nmds.pvca) <- rownames(soil_nmds.vca); colnames(soil_nmds.pvca) <- 'Variance Explained'

## Roots Next ##
end_nono.ps <- subset_samples(end.ps, Compartment != 'Nodule')
end_nono.ps <- subset_taxa(end_nono.ps, taxa_sums(end_nono.ps) > 0)
prop_end.ps <- transform_sample_counts(end_nono.ps, function(x) x/sum(x))
end.wuni <- phyloseq::distance(prop_end.ps, method = 'wunifrac')
soil.pcoa <- phyloseq::ordinate(prop_soil.ps, 'PCoA', distance = 'bray')
soil.tpcoa <- soil.pcoa$vectors[,1:2]

### Performing the NMDS ###
set.seed(248)


## Play around with the k values to get the lowest K that explains => 99% of variation in the linear model ##
soil.nmds <- metaMDS(soil.wuni, 
                     k = 7, try = 20, trymax = 1000, maxit = 999,
                     model = 'global', 
                     autotransform = FALSE, previous.best = soil.pcoa$vectors[,1:7])

soil_nmds.dist <- dist(soil_nmds.scores)
### Calculating the Variance of the Principal Components (or NMDS axes in our case) ###
## Total Model Variance Calculation ##
soil_nmds.scores <- scores(soil.nmds, display = 'sites')
soil_nmds.ffit <- lm(as.vector(soil_nmds.dist) ~ as.vector(soil.wuni))
soil_nmds.totr2 <- summary(soil_nmds.ffit)$r.squared

# Axes Variance Calculation #
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

soil_nmds.dist6 <- dist(soil_nmds.scores[,6])
soil_nmds.fit6 <- lm(as.vector(soil_nmds.dist6)~as.vector(soil.wuni))
soil_nmds.r6 <- summary(soil_nmds.fit6)$r.squared

soil_nmds.dist7 <- dist(soil_nmds.scores[,7])
soil_nmds.fit7 <- lm(as.vector(soil_nmds.dist7)~as.vector(soil.wuni))
soil_nmds.r7 <- summary(soil_nmds.fit7)$r.squared

soil_nmds.comb <- soil_nmds.r1 + soil_nmds.r2 + soil_nmds.r3 + soil_nmds.r4 + soil_nmds.r5 + soil_nmds.r6 + soil_nmds.r7

soil_nmds.axisr <- c()
for(i in 1:ncol(soil_nmds.scores)){
  soil_nmds.axisr[i] <- (get(paste0('soil_nmds.r', i)) / soil_nmds.comb) * soil_nmds.totr2 
}

### Calculating Variance Components ###
library(lme4); packageVersion('lme4')

## Constructing a data.frame that has sample info and their loading scores ##
soil_nmds.load <- cbind(soil.met, soil_nmds.scores)

## Testing the mixed linear model on the first NMDS axis
soil_nmds.vfit1 <- lmer(NMDS1 ~  (1|Plants) + (1|Comps) + (1|Soils) +(1|Plants:Comps) + (1|Plants:Soils) + (1|Comps:Soils) + (1|Comps:Plants:Soils), data = soil_nmds.load)
summary(soil_nmds.vfit1)
soil_nmds.vca1 <- as.data.frame(VarCorr(soil_nmds.vfit1))
## Using Loop to do each NMDS axis ##
soil_nmds.vca <- matrix(nrow = 8, ncol = ncol(soil_nmds.scores))
hold <- c()
for(i in 1:ncol(soil_nmds.scores)){
  hold <- lmer(soil_nmds.scores[,i] ~ (1|Plants) + (1|Comps) + (1|Soils) +(1|Plants:Comps) + (1|Plants:Soils) + (1|Comps:Soils) + (1|Comps:Plants:Soils), data = soil_nmds.load)
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


rownames(soil_nmds.vca) <- c('Compartment x Plant x Soil', 'Plant x Soil', 'Plant x Compartment', 'Soil x Compartment', 'Plant', 'Soil', 'Compartment', 'Residuals')
colnames(soil_nmds.vca) <- colnames(soil_nmds.scores)
soil_nmds.vtot <- colSums(soil_nmds.vca)

### Putting original variance and variance components together for PVCA ###
## weighting each variance component by the variance explained by each NMDS axis ##
soil_nmds.wvca <- matrix(nrow = nrow(soil_nmds.vca), ncol = length(soil_nmds.axisr))
for(i in 1:length(soil_nmds.axisr)){
  for(j in 1:nrow(soil_nmds.vca)){
    soil_nmds.wvca[j,i] <- soil_nmds.vca[j,i]*soil_nmds.axisr[i] 
  }
}

rownames(soil_nmds.wvca) <- rownames(soil_nmds.vca); colnames(soil_nmds.wvca) <- colnames(soil_nmds.vca)
soil_nmds.tvca <- rowSums(soil_nmds.wvca)
soil_nmds.tvca <- as.data.frame(soil_nmds.tvca)
soil_nmds.ptot <- colSums(soil_nmds.tvca)

soil_nmds.pvca <- matrix(nrow = nrow(soil_nmds.tvca), ncol = 1)
for(i in 1:nrow(soil_nmds.tvca)){
  soil_nmds.pvca[i,1] <- soil_nmds.tvca[i,1] / soil_nmds.ptot * 100
}

rownames(soil_nmds.pvca) <- rownames(soil_nmds.vca); colnames(soil_nmds.pvca) <- 'Variance Explained'
#### Maaslin2 Differential Abundance ####
### Bulk Soil Samples ###

soil.otu <- as.data.frame(otu_table(soil.ps))
soil.met <- as(sample_data(soil.ps), 'data.frame')
soil.mett <- t(soil.met)

## Fuzzy Bean ##
bulkfb.met <- filter(soil.met, Plant == 'S. helvola' & Compartment == 'Bulk Soil')
bulkfb.met <- rbind(bulkfb.met, soil.met[98:100,])
bulkfb.otu <- 0
for(i in 1:ncol(soil.otu)){
  if(colnames(soil.otu)[i] %in% rownames(bulkfb.met)){
    bulkfb.otu <- cbind(bulkfb.otu, soil.otu[,i])
  }
}
bulkfb.otu <- bulkfb.otu[,-1]
rownames(bulkfb.otu) <- rownames(soil.otu)
colnames(bulkfb.otu) <- rownames(bulkfb.met)
bulkfb.otu <- as.data.frame(bulkfb.otu)

bulkfb.met$st <- factor(bulkfb.met$Soil_Treatment, levels = c('PSF Soil', 'Non-PSF Soil', 'Common Soil'))
bulkpfb.maas <- Maaslin2(input_data = bulkfb.otu,
                        input_metadata = bulkfb.met,
                        output = '~/Cian_Legume16S/maas/bulkpfb.maas',
                        fixed_effects = c('st'),
                        analysis_method = 'LM',
                        normalization = 'TSS',
                        transform = 'NONE',
                        min_abundance = 0.01,
                        min_prevalence = 0.32,
                        max_significance = 0.05)
