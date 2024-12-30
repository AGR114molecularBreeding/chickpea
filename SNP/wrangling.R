# Install libraries (if needed)
# BiocManager::install("rtracklayer")

# Load library
library(rtracklayer)
library(tibble)
library(dplyr)


# Load GFF file
#  ----------------------------------------------------------------------------
gff_file = "data/GCF_000331145.1_ASM33114v1_genomic.gff"


# Make GFF gene from genome annotation
#  ----------------------------------------------------------------------------
gff_gene = import.gff(con = gff_file,
                      feature.type = "gene",
                      colnames = c("type", "gene"))

gff_gene = unique(gff_gene)




# Make GFF exon from genome annotation
#  ----------------------------------------------------------------------------
gff_exon = import.gff(con = gff_file, 
                      feature.type = "exon", 
                      colnames = c("type", "gene"))

gff_exon = unique(gff_exon)




# Make GFF CDS from genome annotation
#  ----------------------------------------------------------------------------
gff_cds = import.gff(con = gff_file, 
                     feature.type = "CDS", 
                     colnames = c("type", "gene"))

gff_cds = unique(gff_cds)
length(gff_cds)




# Save gff objects into .rda file
#  ----------------------------------------------------------------------------
save(gff_gene, gff_exon, gff_cds, file = "data/gff.rda")

