# Rscript ~/PSF_MBIOME/analysis.R --raw_soil ~/test_PSF/reads/soil_reads --raw_root ~/test_PSF/reads/endo_reads --soil_metadata ~/PSF_MBIOME/metadata/soil_metadata.csv --pheno ~/PSF_MBIOME/nodnbio.csv --reference ~/PSF_MBIOME/reference/rdp_19_toGenus_trainset.fa.gz | cat > PSF_log.txt #

#### Argument Parsing ####
library(optparse); packageVersion("optparse")
option_list <- list(
  make_option("--raw_soil", type = "character", help = "filepath containing the raw, untrimmed reads for the reads generated from the v4 primers (bulk soil, rhizosphere, and some nodule samples)"),
  make_option("--raw_root", type = "character", help = "filepath containing the raw, untrimmed reads for the reads generated from the v5-v7 primers( root endosphere and some nodule samples)"),
  make_option("--soil_metadata", type = "character", help = "filepath that contains the Comma Separated Values (csv) file of the soil metadata"),
  make_option("--pheno", type = "character", help = "filepath that contains the Comma Separated Values (csv) file of the phenotype data (biomass and nodule counts)"),
  make_option("--reference", type = "character", help = "filepath that contains the reference database to assign taxonomy to the reads"),
  make_option("--root_metadata", type = "character", help = "filepath that contains the Comma Separated Values (csv) file of the root metadata"),
  make_option("--zymo", type = "character", help = "filepath to where the reference fasta files are for the ZYMO mock culture"))

opt <- parse_args(OptionParser(option_list=option_list))

soil.dir <- opt$raw_soil
root.dir <- opt$raw_root
soil.met <- opt$soil_metadata
nodnbio <- opt$pheno
reference <- opt$reference
root.met <- opt$root_metadata
zymo <- opt$zymo

#### Phylogenetic Tree Construction for roots ####
# Output the reads into a fasta file #
library(Biostrings)
writeXStringSet(as.character(root$dna, filepath = "./reads/root_input.fasta")

# Perform a multiple sequence alignment using MAFFT #
system('mafft --auto --thread -1 ./reads/root_input.fasta > ./reads/roots_aligned.fasta')

# Construct a tree using IQTree with a general time reversible model with a gamma distribution and invariant site copies #
system('iqtree -s ./reads/roots_aligned.fasta -m GTR+G+I -nt AUTO')


save.image("./root2.RData")
