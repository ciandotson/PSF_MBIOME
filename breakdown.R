#### Argument Parsing ####
if(!requireNamespace('optparse', quietly = TRUE)) install.packages("optparse")
library(optparse); packageVersion("optparse")
option_list <- list(
  make_option("--raw_soil", type = "character", help = "filepath containing the raw, untrimmed reads for the reads generated from the v4 primers (bulk soil, rhizosphere, and some nodule samples)"),
  make_option("--raw_root", type = "character", help = "filepath containing the raw, untrimmed reads for the reads generated from the v5-v7 primers( root endosphere and some nodule samples)"))

opt <- parse_args(OptionParser(option_list=option_list))

soil.dir <- opt$raw_soil
root.dir <- opt$raw_root

# Save the outputs and messages to a log file # 
sink(file = './psf.log', append = TRUE, type = c("output", "message")) # Redirects stdout (e.g., print, cat) and stderr (e.g., warnings, message) #
cat("## Script started at", Sys.time(), "\n\n")

# Set the working directory to the cloned github repository #
setwd("~/PSF_MBIOME")

load('./test.RData')

library(dada2)
library(phyloseq)
library(ggplot2)
library(microbiome)
library(patchwork)
library(cgwtools)
library(ape)
library(dplyr)
library(Biostrings)
library(ggprism)
library(rBLAST)
library(Polychrome)

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
save.image("./test.RData")
raw_root.dna <- Biostrings::DNAStringSet(taxa_names(raw_root.ps))
names(raw_root.dna) <- taxa_names(raw_root.ps)
raw_root.ps <- merge_phyloseq(raw_root.ps, raw_root.dna)

raw_root.ps <- subset_taxa(raw_root.ps, Phylum != "Cyanobacteriota")
raw_root.ps <- subset_taxa(raw_root.ps, Phylum != "Plantae")

raw_root.ps <- subset_taxa(raw_root.ps, taxa_sums(raw_root.ps) > 100)
colnames(tax_table(raw_root.ps)) <- c('rdp_Kingdom', 'rdp_Phylum', 'rdp_Class', 'rdp_Order', 'rdp_Family', 'rdp_Genus')

taxa_names(raw_root.ps) <- paste0("ASV", seq(ntaxa(raw_root.ps)))
decompose_ps(raw_root.ps, 'raw_root')
resave(raw_root.ps, file = './psf_abridged.RData')

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
save.image("./test.RData")

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

# Change the taxa names to represent comparative abundance and lowest identification level #
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

save.image("./test.RData")
# produce final decomposed phyloseq object
decompose_ps(root.ps, 'root')

#### Phylogenetic Tree Construction for roots ####
# Output the reads into a fasta file #
writeXStringSet(root$dna, "./reads/root_input.fasta")

# Perform a multiple sequence alignment using MAFFT #
system('mafft --auto --thread -1 ./reads/root_input.fasta > ./reads/roots_aligned.fasta')

# Construct a tree using IQTree with a general time reversible model with a gamma distribution and invariant site copies #
system('iqtree -s ./reads/roots_aligned.fasta -m GTR+G+I -nt AUTO')

# Read the tree using ape and check that the tip labels match #
root.tre <- read.tree("./reads/roots_aligned.fasta.treefile")
root.tre$tip.label <- sub("^(ASV[0-9]+)_([^_]+)_$", "\\1(\\2)", root.tre$tip.label)

# Combine final phyloseq object for root samples #
decompose_ps(root.ps, 'root')
root$tax$ASV <- rownames(root$tax)
root$tax <- as.matrix(root$tax)
root.ps <- phyloseq(otu_table(root$otu, taxa_are_rows = TRUE),
                    sample_data(root$met),
                    tax_table(root$tax),
                    refseq(root$dna),
                    phy_tree(root.tre))
resave(root.ps, file = './psf_abridged.RData')
save.image("./test.RData")

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
root_nod.colr[ntaxa(root_nod.ps) + 1,] <- "#D4D4D4" 
rownames(root_nod.colr)[ntaxa(root_nod.ps) + 1] <- "Other" 

# Create a tree for ASVs that are found in the nodules and are in the Order Hyphomicrobiales #
root_hyph.ps <- subset_taxa(root_nod.ps, Order == "Hyphomicrobiales")
root_hyph_tre.plot <- plot_tree(root_hyph.ps, label.tips = "taxa_names", ladderize = TRUE, color = "Family")

# save a new phyloseq object without the nodules #
root.ps <- subset_samples(root.ps, Compartment != "Nodule")
root.ps <- subset_taxa(root.ps, taxa_sums(root.ps) > 0)
decompose_ps(root.ps, 'root')
save.image("PSF.RData")

# This just outputs the final time in which the Rscript finishes running #
cat("\n## Script finished at", Sys.time(), "\n")
sink(type = "message")
sink()