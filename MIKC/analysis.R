# MADS-box.R -- a set of common functions for analysis of MADS-box MIKC-type family in Ca
# Copyright (C) 2023 Jose V. Die  <jose.die@uco.es> 
# Copyright (C) 2023 Alejandro Carmona  <b62cajia@uco.es> 
# Distributed under terms of the MIT license.


# Dependencies
library(dplyr)
library(stringr)
library(refseqR)
library(rentrez)

# Load functions 
source("MIKC/functions.R")

# Set wd
setwd("MIKC/")

## -----------------------------------------------------------------------------
##                    Gene Family Identification
## -----------------------------------------------------------------------------


# --- Get the mRNA sequence ----

# Read object 
xps <- read.csv("dat/ZT55WJYH01N-Alignment-HitTable.csv")

# Make a vector 
xps <- xps$list 
# Fix some typo : 
xps[88] = "XP_004498126.1"

# Get XM ids from XP ids. 
xms <- sapply(xps, function(xp) getXM(xp), USE.NAMES = FALSE)

# Save an object containing the corresponding XM, XP
seqs <- list(unname(xms), xps)
save(seqs, file = "res/sequences.rda")

# Make a multi-fasta file from XM ids 
save_CDSfasta_from_xms(xms, nameFile = "res/mads_cds")



# --- Iterative BLAST  ----
hitFile <- "dat/YN34UYFE013-Alignment-HitTable.csv"
res <- blast_best_homolog(hitFile)
ca_MADSbox <-  res %>% pull(subject)
xps <- unique(xps) # Unique XP ids from the original first list 

# Add 1st XP id 
ca_MADSbox <- c(ca_MADSbox, res$query[1])



# --- Get the mRNA sequence ----

# Read object 
xps <- read.csv("dat/Mads_XP_unique.csv")

# Make a vector 
xps <- xps$list
xps <- unique(xps)


# --- Iterative BLAST  ----
hitFile <- "dat/ZT55WJYH01N-Alignment-HitTable.csv"

blast_res <- blast_homologs(hitFile, ident = 55, cover = 90)


# Gene family is made of : 
#   original xp ids ('xp' object)
#   'subject' column ('blast_res' object)

ca_MADSbox <- unique(c(xps, blast_res$subject))

length(xps)
length(ca_MADSbox)


# Make some plots for checking 
hist(blast_res$identity, main = "iterative BLAST", xlab = "%identity")
hist(blast_res$cov_long, main = "iterative BLAST", xlab = "%coverage")
hist(blast_res$evalue,   main = "iterative BLAST", xlab = "Evalue")



# Get XM ids from XP ids. 
xms <- sapply(ca_MADSbox, function(xp) getXM(xp), USE.NAMES = FALSE)

# Save an object containing the corresponding XM, XP
seqs <- list(unname(xms), ca_MADSbox)


save(seqs, file = "res/sequences.rda")

# Make a multi-fasta file from XM ids 
save_CDSfasta_from_xms(xms, nameFile = "res/mads_cds")





## -----------------------------------------------------------------------------
##                    Gene Family Table
## -----------------------------------------------------------------------------

# --- Make a gene family Table  

targets <- c("XP_004513719.1", "XP_004485955.1", "XP_004498126.1", 
             "XP_004492666.1", "XP_027188084.1")

tdat <- characterizeTable(targets)

# Sort data set by Chr + start coordinate
tdat %>% 
  arrange(Chr, chr_s)

# Data set based on 1 gene model / locus
tdat %>% 
  arrange(Chr, chr_s) %>% 
  distinct(LOC, .keep_all = TRUE)

# Get info on isoforms : n. loci with > 1 protein
tdat %>% 
  arrange(LOC, desc(AA)) %>% 
  count(LOC, sort = T) %>% 
  filter(n > 1)

nisof <- tdat %>% 
  arrange(LOC, desc(AA)) %>% 
  count(LOC) %>% 
  pull(n)

tdat <- tdat %>%  
  arrange(LOC, desc(AA)) %>% 
  distinct(LOC, .keep_all = TRUE) %>% 
  mutate(n_isof = nisof)

write.csv(tdat, file = "res/my_table.csv", row.names=FALSE)


# Make some tidy 
rm(nisof, targets)
rm(ca_MADSbox, hitFile, ids30, ids50, ids70)
rm(tlengths_query)
rm(tlengths_sub)





## -----------------------------------------------------------------------------
##                    GRanges object
## -----------------------------------------------------------------------------

# Negative widths are not allowed by IRanges, so we have to define the true 
# coordinates and build a new data frame, which is the input for the function 
# 'makeGRangesFromDataFrame'

# Install package and load library 
#BiocManager::install("GenomicRanges")
library(GenomicRanges)

gr <- tdat %>% 
  transmute(LOC, Chr, Strand, 
            Coordstart = ifelse(Strand == "-", chr_e, chr_s),
            Coordend = ifelse(Strand == "+", chr_e, chr_s))


# Add new column : seq.length ('bp_cut' = 150) to map the TSS coordinate (promoter analysis) 
# it requires to customize the bp_length when the ATG is in the end of the exon
gr <- gr %>% 
  dplyr::mutate(bp_cut = 150) %>% 
  dplyr::mutate(bp_cut = ifelse(LOC == "LOC101493118", 14, bp_cut), 
                bp_cut = ifelse(LOC == "LOC101504656", 27, bp_cut), 
                bp_cut = ifelse(LOC == "LOC101509413", 3,  bp_cut), 
                bp_cut = ifelse(LOC == "LOC101510075", 39, bp_cut))

gr <- makeGRangesFromDataFrame(gr, start.field = "Coordstart", 
                               end.field = "Coordend", 
                               strand.field = "Strand", ignore.strand = F, 
                               seqnames.field = "Chr", 
                               keep.extra.columns = TRUE)

# Add genome
genome(gr) = "ASM33114v1"

# Filter out 'unplaced' seqs
gr <- gr[seqnames(gr) != "Un"]

# Add seqlengths 
seqinfo(gr)
              
# NCBI, Frontier : https://www.ncbi.nlm.nih.gov/genome/?term=chickpea 
# size = c(48.36, 36.63, 39.99, 49.19, 48.17, 59.46, 48.96, 16.48)
genome(gr)
isCircular(gr)
#seqlengths(gr) = c(NA, 48.36, 39.99, 49.19)

# Order the GR object by chr and then by region location
sort(table(seqnames(gr)), decreasing = TRUE)
# gr[order(gr), ]
gr <- sort(gr)

# Save object 
save(seqs, gr, file = "res/sequences.rda")



## Manipulating GR object
ranges(gr)[1:3]
sort(table(seqnames(gr)), decreasing = TRUE)
gr[seqnames(gr) == "Ca1"]
gr[seqnames(gr) == "Ca3"]



# Save gr into .bed file)
##The only trick is remembering the BED uses 0-based coordinates. 
##If you have a GRangesList rather than a GRanges object, just use unlist(gr) 
##In place of gr (things should still be in the same order).

df <- data.frame(seqnames=seqnames(gr),
                 starts=start(gr)-1,
                 ends=end(gr),
                 names=c(rep(".", length(gr))),
                 scores=c(rep(".", length(gr))),
                 strands=strand(gr))


write.table(df, file="tdat.bed", quote=F, sep="\t", row.names=F, col.names=F)

## Now, opposite direction : read .bed file to build a GRanges object
bed = read.table("tdat.bed",header=F)
colnames(bed) <- c('chr','start','end','id','strand')
head(bed)
gr <- with(bed, GRanges(chr, IRanges(start+1, end), strand=strand, id=id))


# --- Promoter sequence -----
# GR contains the locus (LOC) coordenates 
# We need the transcriptional start site (TSS) coordinate

library(Biostrings)
## {Biostrings} Extract basic information about FASTA files
## without actually loading the sequence data:
length(fasta.seqlengths("res/mads_cds.fasta"))
fasta.seqlengths("res/mads_cds.fasta")
names(fasta.seqlengths("res/mads_cds.fasta"))

# Load mRNA sequences + convert into a DNAStringSet {Biostrings}
mads_cds = readDNAStringSet("res/mads_cds.fasta")

# Convert names to LOCnames
my_names = names(mads_cds)
my_names <- sapply(my_names, function(i) names2LOC(i), USE.NAMES = F)
names(mads_cds) = my_names
mads_cds


library(BSgenome.Carietinum.NCBI.v1)
## {BSgenome.Carietinum.NCBI.v1} Cicer arietinum full genome (ver. NCBI)
genome = BSgenome.Carietinum.NCBI.v1
genome$Ca1 #this loads Chr1 into computer memory
ca1 = genome[[1]]  #loads Chr1 into computer memory


# Parse `TSScoordinates` on Chrs presented in the GR object
levels(seqnames(gr))

gr.tss= GRangesList(sapply(paste0("Ca", 1:8), function(i) TSScoordinates(gr, mads_cds, i)))
gr.tss = unlist(gr.tss)

# Save object 
save(seqs, gr, gr.tss, file = "res/sequences.rda")



# Now we use the TSS coordinates to extract -1500 nt representing the promoter
gr.promoters = promoters(gr.tss, upstream=1500, downstream=9)

prom = sapply(paste0("Ca", 1:8), function(i) {
  get_promoters(gr = gr.promoters, chr = i)})

prom = DNAStringSet(unlist(prom))

names(prom) =  gr.promoters$LOC
names(prom) = paste0(names(prom), "_promoter")

# Check ATG
subseq(prom, start = 1501, end = 1509)


# Write sequences to a fasta file
writeXStringSet(prom, filepath="res/promoters.fasta", format="fasta")  


# Estimate base frequency to create synthetic data 
## {Biostrings} Extract basic information about FASTA files
## without actually loading the sequence data:
length(fasta.seqlengths("res/promoters.fasta"))
names(fasta.seqlengths("res/promoters.fasta"))
fasta.seqlengths("res/promoters.fasta")


# Read StringSet
promoters <- readDNAStringSet("res/promoters.fasta")
alphabetFrequency(promoters)
dinucleotideFrequency(promoters)
dinucleotideFrequency(promoters, as.prob = TRUE)
dinucleotideFrequency(promoters[[1]], as.prob = TRUE, as.matrix = T)
dinucleotideFrequency(promoters[[2]], as.prob = TRUE, as.matrix = T)
sort(dinucleotideFrequency(promoters[[1]], as.prob = TRUE))
sort(dinucleotideFrequency(promoters[[1]], as.prob = TRUE), decreasing = TRUE)
head(sort(trinucleotideFrequency(promoters[[1]], as.prob = TRUE), 
          decreasing = T))


# Get the counts of oligonucleotides overlapping by one nucleotide 
head(sort(oligonucleotideFrequency(promoters[[1]], width = 5, step = 4), decreasing = T))
head(sort(oligonucleotideFrequency(promoters[[1]], width = 6, step = 5), decreasing = T))
head(sort(oligonucleotideFrequency(promoters[[1]], width = 7, step = 6), decreasing = T))
head(sort(oligonucleotideFrequency(promoters[[1]], width = 8, steo = 7), decreasing = T))

letterFrequency(promoters[[1]], "GC", as.prob = T)
matchPattern("TATATATA", promoters[[1]])
countPattern("TATATATA", promoters[[1]])


letterFrequency(promoters, "GC", as.prob = T)


## Get the less and most represented 6-mers:
f6 <- oligonucleotideFrequency(promoters[[1]], 6)
f6[f6 == min(f6)]
f6[f6 == max(f6)]


## With 20-nucleotide forward and reverse probes:
Fprobe <- "AGCTCCGAGTTCCTGCAATA"
Rprobe <- "CGTTGTTCACAAATATGCGG"
matchProbePair(Fprobe, Rprobe, promoters[[1]]) #theoretical amplicon : None


# Base composition of the dataset 
apply(alphabetFrequency(promoters, as.prob = T)[,1:4], 2, mean)
prom_freq <- apply(alphabetFrequency(promoters, as.prob = T)[,1:4], 2, mean)


# Weighted random sample
sample(c(1, 2, 3), size = 10, replace = TRUE, prob = c(0.5, 0.1, 0.4))
sample(c("A", "T", "C", "G"), size = 10, replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25))
table(sample(c("A", "T", "C", "G"), size = 1000, replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25)))

# Generate synthetic promoters

prom_synth <- function(size, frequencies) {
    
  paste0(sample(c("A", "C", "G", "T"), size = 1500, replace = TRUE, prob = frequencies), collapse = "")
}

prom_synth(1500, prom_freq)


# Run Monte Carlo simulation
times = 5
mc_prom_synth <-  replicate(times, prom_synth(1500, prom_freq))

writeXStringSet(DNAStringSet(mc_prom_synth), filepath="res/promoters_synthetic.fasta", format="fasta")  



# --- Outline Prommoter Analyis ----

# Search for number of occurrences in the actual dataset
# Estimate %GC content
# Estimate for number of occurrences in simulated dataset
# Estimate P value for the n. occurrences in the actual dataset
