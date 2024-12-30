# a set of common functions for analysis of SNP in Ca
# Copyright (C) 2023 Jose V. Die  <jose.die@uco.es>
# Distributed under terms of the MIT license. 

# Load GFF objects
#  ----------------------------------------------------------------------------
load("data/gff.rda")


# Load VCF file (from the sequencing service)
#  ----------------------------------------------------------------------------
# Make table of SNPS per chromosome
# system("egrep -v '^#' Chickpea.MCR50.snps.vcf | cut -f 1 | sort | uniq -c | sort -r > data/snp_per_chrs")

snp_chr <- read.table("data/snp_per_chrs.txt")

snp_chr <- snp_chr %>% 
  as_tibble(snp_chr) %>% 
  select(V2, V1) %>% 
  rename("chromosome" = V2, "n" = V1)

snp_chr %>% top_n(8)

# keep data only for the 8 chromosomes
ca_chrs <- snp_chr %>% top_n(8) %>% pull(chromosome)


# Make a GRanges object with the SNP coordinates
#system("egrep -v '^#' data/Chickpea.MCR50.snps.vcf | cut -f 1,2 > data/coordinates_vcf.txt")
#snp_coords <- read.table("data/coordinates_vcf")             #5891 SNP
#snp_coords %>% count(V1) %>% top_n(8)

# Alejandro maneja un archivo final con 1,071 SNP (8 chroms.) en lugar del original 2,564 (8 chroms.)
#snp_coords <- read.delim("data/coordinates_original.txt")    #2564 SNP
snp_coords <- read.delim("data/coordenadas_final.txt")        #1071 SNP

snp_coords <- snp_coords %>% 
  select(2:3) %>% 
  as_tibble(snp_chr) %>% 
  # keep SNP only in the chormosomes
  filter(chrom %in% ca_chrs) %>%
  # unify names : GenBank [.vcf file -> RefSeq (.gff file) ]
  mutate(chrom = case_when(chrom == "CM001764.1" ~ "NC_021160.1", 
                           chrom == "CM001765.1" ~ "NC_021161.1", 
                           chrom == "CM001766.1" ~ "NC_021162.1",
                           chrom == "CM001767.1" ~ "NC_021163.1",
                           chrom == "CM001768.1" ~ "NC_021164.1",
                           chrom == "CM001769.1" ~ "NC_021165.1",
                           chrom == "CM001770.1" ~ "NC_021166.1",
                           chrom == "CM001771.1" ~ "NC_021167.1"))

# SNP total & SNP per Chr
snp_coords %>% nrow()
snp_coords %>% count(chrom)

# Make GRanges with SNP
snp = GRanges(seqnames = snp_coords$chrom, 
              IRanges(start = snp_coords$pos , width = 1), 
              score = 5)

  
# make some tidy 
rm(ca_chrs, snp_chr)
rm(snp_coords)


# GFF gene : find overlaps
#  ----------------------------------------------------------------------------
hits_gene <- findOverlaps(snp, gff_gene)
snp[unique(queryHits(hits_gene))] # coordenada de 395 SNP en genes
by_genes <- gff_gene[unique(subjectHits(hits_gene))] # genes que presentan 395 SNP
true_genes <- unique(elementMetadata(by_genes)[,2])

snp[unique(queryHits(hits_gene))] # coordenada de 395 SNP en genes
# SNP genicos agrupados por chromosomas
seqnames(snp[unique(queryHits(hits_gene))])

snp[-queryHits(hits_gene)]# coordenada de 676 SNP INTERGENICOS
# intergenicos agrupados por chromosomas
seqnames(snp[-queryHits(hits_gene)])


# GFF exon : find overlaps
#  ----------------------------------------------------------------------------
hits_exon <- findOverlaps(snp, gff_exon)
snp[queryHits(hits_exon)] # coordenada de 199 SNP en exones

## OJO : exones pueden ser :
#* genes
#* pseudo-genes

# # Ejemplo pseudo-gene con exones :
# system("egrep -v '^#' data/GCF_000331145.1_ASM33114v1_genomic.gff | grep LOC105851630 | less -S")
# system("egrep -v '^#' data/GCF_000331145.1_ASM33114v1_genomic.gff | grep LOC101488310 | less -S")

#194 SNP en exones de true_genes
gff_exon[subjectHits(hits_exon)][elementMetadata(gff_exon[subjectHits(hits_exon)])[,2] %in% true_genes]

#5 SNP en exones de pseudo_genes
gff_exon[subjectHits(hits_exon)][!elementMetadata(gff_exon[subjectHits(hits_exon)])[,2] %in% true_genes]
pseudo_genes <- unique(gff_exon[subjectHits(hits_exon)][!elementMetadata(gff_exon[subjectHits(hits_exon)])[,2] %in% true_genes]$gene)
pseudo_genes



# GFF CDS : find overlaps
#  ----------------------------------------------------------------------------
hits_cds <- findOverlaps(snp, gff_cds)
snp[queryHits(hits_cds)] # coordenada de 122 SNP

#122 SNP en CDS
gff_cds[subjectHits(hits_cds)]
#122 SNP en CDS de true_gene
gff_cds[subjectHits(hits_cds)][elementMetadata(gff_cds[subjectHits(hits_cds)])[,2] %in% true_genes]

gff_cds[subjectHits(hits_cds)] # CDS que presentan 122 SNP
unique(gff_cds[subjectHits(hits_cds)]) # 78 CDS unicos que presentan 122 SNPs

# agrupados por LOC
sort(table(elementMetadata(unique(gff_cds[subjectHits(hits_cds)]))[,2]), decreasing = T)

by_exon_true_gene <- gff_exon[subjectHits(hits_exon)][elementMetadata(gff_exon[subjectHits(hits_exon)])[,2] %in% true_genes]
by_cds_true_gene <- gff_cds[subjectHits(hits_cds)][elementMetadata(gff_cds[subjectHits(hits_cds)])[,2] %in% true_genes]

#63 SNP en exones de true_gene que no son CDS de true gene : UTR de true_gene
by_exon_true_gene[!elementMetadata(by_exon_true_gene)[,2] %in% elementMetadata(by_cds_true_gene)[,2]]

# ojo !!! : 63 + 122 != 194
