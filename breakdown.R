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

load("./soil2.RData")
# Load the metadata #
library(phyloseq); packageVersion("phyloseq")
library(dada2)
library(Biostrings)
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

save.image("./soil2.RData")
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

save.image("./soil2.RData")
# produce final decomposed phyloseq object
decompose_ps(soil.ps, 'soil')

#### Phylogenetic Tree Construction for Soils ####
# Output the reads into a fasta file #
writeXStringSet(soil$dna, "./reads/soil_input.fasta")

# Perform a multiple sequence alignment using MAFFT #
system('mafft --auto --thread -1 ./reads/soil_input.fasta > ./reads/soils_aligned.fasta')

# Construct a tree using IQTree with a general time reversible model with a gamma distribution and invariant site copies #
system('iqtree -s ./reads/soils_aligned.fasta -m GTR+G+I -nt AUTO')

#### Root Primer Removal ####
# Ensure you have the right files #
if(file.exists("./soil2.RData")){
  rm(list=ls())
}
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

save.image("./root2.RData")
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

# track reads through the pipeline #
getN <- function(x) sum(getUniques(x))
root_final.track <- cbind(root_prefilt.track[,1], root_prefilt.track[,2], root_postfilt.track[,2], sapply(root.fdada, getN), sapply(root.rdada, getN), sapply(root.remerged, getN), colSums(root_nochim.st))
colnames(root_final.track) <- c("pre-cutadapt", "post-cutadapt", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(root_final.track) <- root.names
root_final.track <- as.data.frame(root_final.track)

# Assign Taxonomy #
root_rdp.taxa <- assignTaxonomy(rownames(root_nochim.st), refFasta = reference, multithread = TRUE, verbose = TRUE)
root_rdp.taxa <- as.matrix(root_rdp.taxa)

save.image("./root2.RData")
# Load the metadata #
root_raw.met <- read.csv2(root_metadata, sep = ',')
rownames(root_raw.met) <- root_raw.met$Sample
root_raw.met <- root_raw.met[,c('Sample.Name', 'Plant.Species', 'Soil.Origin', 'Compartment')]
rownames(root_raw.met) <- sub("^([^_]+_[^_]+)_.*$", "\\1", rownames(root_raw.met))

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
taxa_names(root.ps) <- paste0('ASV', seq(ntaxa(root.ps)))
decompose_ps(root.ps, 'root')
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

save.image("./root.RData")
#### Phylogenetic Tree Construction for roots ####
# Output the reads into a fasta file #
writeXStringSet(root$dna, "./reads/root_input.fasta")

# Perform a multiple sequence alignment using MAFFT #
system('mafft --auto --thread -1 ./reads/root_input.fasta > ./reads/roots_aligned.fasta')

# Construct a tree using IQTree with a general time reversible model with a gamma distribution and invariant site copies #
system('iqtree -s ./reads/roots_aligned.fasta -m GTR+G+I -nt AUTO')

save.image("./root2.RData")
