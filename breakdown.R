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


# Perform a multiple sequence alignment using MAFFT #
system('mafft --auto --thread -1 ./reads/soil_input.fasta > ./reads/soils_aligned.fasta')

# Construct a tree using IQTree with a general time reversible model with a gamma distribution and invariant site copies #
system('iqtree -s ./reads/soils_aligned.fasta -m GTR+G+I -nt AUTO')