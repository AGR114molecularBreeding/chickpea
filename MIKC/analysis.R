# MADS-box.R -- a set of common functions for analysis of MADS-box family in Ca
# Copyright (C) 2023 Jose V. Die  <jose.die@uco.es>
# Distributed under terms of the MIT license.


# Dependencies
library(dplyr)
library(stringr)
library(refseqR)
library(rentrez)

# Load functions 
source("mads_box/functions.R")

# Set wd
setwd("mads_box/")

## -----------------------------------------------------------------------------
##                    Gene Family Identification
## -----------------------------------------------------------------------------


# # --- Get the mRNA sequence (Feb 13, 2023) ----
# 
# # Read object 
# xps <- read.csv("dat/HITS_Hom_chickpea_Ilaria.csv")
# 
# # Make a vector 
# xps <- xps$list
# 
# # fix some typo : 
# xps[88] = "XP_004498126.1"
# 
# # Get XM ids from XP ids. 
# xms <- sapply(xps, function(xp) getXM(xp), USE.NAMES = FALSE)
# 
# # Save an object containing the corresponding XM, XP
# seqs <- list(unname(xms), xps)
# save(seqs, file = "res/sequences.rda")
# 
# # Make a multi-fasta file from XM ids 
# save_CDSfasta_from_xms(xms, nameFile = "res/mads_cds")
# 
# 
# 
# # --- Iterative BLAST (Feb 13, 2023) ----
# hitFile <- "dat/YN34UYFE013-Alignment-HitTable.csv"
# res <- blast_best_homolog(hitFile)
# ca_MADSbox <-  res %>% pull(subject)
# 
# xps # Illatia first list 
# xps <- unique(xps) # unique XP ids from Ilaria first list 
# 
# # add 1st XP id 
# ca_MADSbox <- c(ca_MADSbox, res$query[1])
# 
# 
# # --- checks 
# hist(res$evalue)
# 
# hist(res$identity)
# abline(v = 70, col = "red")
# 
# ids70 <- res %>% 
#   filter(identity >=70) %>% 
#   pull(subject)
# 
# ids70 %in% xps
# ids70[!ids70 %in% xps]
# 
# 
# "XP_027187420.1" %in% xps
# 
# 
# ids50 <- res %>% 
#   filter(identity >=50) %>% 
#   pull(subject)
# 
# ids50[!ids50 %in% xps]
# 
# 
# ids30 <- res %>% 
#   filter(identity >=30) %>% 
#   pull(subject)
# 
# ids30[!ids30 %in% xps]





# --- Get the mRNA sequence (Feb 27, 2023) ----

# Read object 
xps <- read.csv("dat/Mads_XP_unici.csv")

# Make a vector 
xps <- xps$list
xps <- unique(xps)


# --- Iterative BLAST (Feb 27, 2023) ----
hitFile <- "dat/ZT55WJYH01N-Alignment-HitTable.csv"

blast_res <- blast_homologs(hitFile, ident = 55, cover = 85)


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


# Make some plots
boxplot(tdat$exon, main = "MADS-box Family", ylab = "Exon number")
stripchart(tdat$exon, vertical = TRUE, method = "jitter", 
           add = TRUE, col = "steelblue", pch = 19, cex = 0.9)

boxplot(tdat$AA, main = "MADS-box Family", ylab = "Protein (aa)")
stripchart(tdat$AA, vertical = TRUE, method = "jitter", 
           add = TRUE, col = "steelblue", pch = 19, cex = 0.9)

boxplot(tdat$mol_wt, main = "MADS-box Family", ylab = "Molecular Weight (kDa)")
stripchart(tdat$mol_wt, vertical = TRUE, method = "jitter", 
           add = TRUE, col = "steelblue", pch = 19, cex = 0.9)

# Make some tidy 
rm(nisof, targets)
rm(ca_MADSbox, hitFile, ids30, ids50, ids70)
rm(tlengths_query)
rm(tlengths_sub)





## -----------------------------------------------------------------------------
##                    GRanges object
## -----------------------------------------------------------------------------

# Negative widths are not allowed by IRanges, so let's define the true 
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

# filter out 'unplaced' seqs (9 seqs! ; filter out??)
gr <- gr[seqnames(gr) != "Un"]

# add seqlengths 
seqinfo(gr)
# NCBI, Frontier : https://www.ncbi.nlm.nih.gov/genome/?term=chickpea 
# size = c(48.36, 36.63, 39.99, 49.19, 48.17, 59.46, 48.96, 16.48)
genome(gr)
isCircular(gr)
#seqlengths(gr) = c(NA, 48.36, 39.99, 49.19)

# order the GR object by chr and then by region location
sort(table(seqnames(gr)), decreasing = TRUE)
# gr[order(gr), ]
gr <- sort(gr)

# Save object 
save(seqs, gr, file = "res/sequences.rda")


## Puedo a??adir metadatos GENEID to match gene information across different Genbank databases) 
# Introduction to Bioconductor / Week 2/ Genes as GRanges / min 1.21


## Manipulating GR object
ranges(gr)[1:3]
sort(table(seqnames(gr)), decreasing = TRUE)
gr[seqnames(gr) == "Ca1"]
gr[seqnames(gr) == "Ca3"]

# si hago plotGRanges el plot es similar a plotIRanges pero sale el
# nombre del chromosoma. Y admite color!
# Introduction to Bioconductor / Week 3/ Operating with GRanges / min 6. 21

# Grabo gr (into .bed file) para no tener que hacer todo el recorrido
##The only trick is remembering the BED uses 0-based coordinates. 
##If you have a GRangesList rather than a GRanges object, just use unlist(gr) 
##in place of gr (things should still be in the same order).

# df <- data.frame(seqnames=seqnames(gr),
#                  starts=start(gr)-1,
#                  ends=end(gr),
#                  names=c(rep(".", length(gr))),
#                  scores=c(rep(".", length(gr))),
#                  strands=strand(gr))
# 
# 
# write.table(df, file="tdat.bed", quote=F, sep="\t", row.names=F, col.names=F)

## Now, opposite direction : read .bed file to build a GRanges object
# bed = read.table("tdat.bed",header=F)
# colnames(bed) <- c('chr','start','end','id','strand')
# head(bed)
# gr <- with(bed, GRanges(chr, IRanges(start+1, end), strand=strand, id=id))


# --- Promoter sequence -----
# GR contains the locus (LOC) coordenates 
# We need the transcriptional start site (TSS) coordinate

library(Biostrings)
## {Biostrings} Extract basic information about FASTA files
## without actually loading the sequence data:
length(fasta.seqlengths("res/mads_cds.fasta"))
fasta.seqlengths("res/mads_cds.fasta")
names(fasta.seqlengths("res/mads_cds.fasta"))

# Load mRNA seqs + convert into a DNAStringSet {Biostrings}
mads_cds = readDNAStringSet("res/mads_cds.fasta")

# Convert names to LOCname
my_names = names(mads_cds)
my_names <- sapply(my_names, function(i) names2LOC(i), USE.NAMES = F)
names(mads_cds) = my_names
mads_cds


library(BSgenome.Carietinum.NCBI.v1)
## {BSgenome.Carietinum.NCBI.v1} Cicer arietinum full genome (ver. NCBI)
genome = BSgenome.Carietinum.NCBI.v1
#genome$Ca1 #this loads Chr1 into computer memory
#ca1 = genome[[1]]  #loads Chr1 into computer memory


# Parse `TSScoordinates` on Chrs presented in the GR object
# [esto es mas eficiente que hacer loops sobre el BSGenome]
levels(seqnames(gr))

# gr.tss = GRangesList(TSScoordinates(gr, mads_cds, "Ca1"),
#                      TSScoordinates(gr, mads_cds, "Ca3"), 
#                      TSScoordinates(gr, mads_cds, "Ca4"))
# 

gr.tss= GRangesList(sapply(paste0("Ca", 1:8), function(i) TSScoordinates(gr, mads_cds, i)))
gr.tss = unlist(gr.tss)

# Save object 
save(seqs, gr, gr.tss, file = "res/sequences.rda")



# Now we use the TSS coordinates to extract -1500 nt representing the promoter
gr.promoters = promoters(gr.tss, upstream=1500, downstream=9)

# Con las coordenadas del promotor puedo extraer la seq del promotor 
# Hago TSScoordinates con los Chr que tenga en el GR
# esto es mas eficiente que hacer loops sobre el BSGenome 
# prom = DNAStringSet(c(get_promoters(gr = gr.promoters, chr = "Ca1"),
#                       get_promoters(gr = gr.promoters, chr = "Ca3"), 
#                       get_promoters(gr = gr.promoters, chr = "Ca4")))
# 

prom = sapply(paste0("Ca", 1:8), function(i) {
  get_promoters(gr = gr.promoters, chr = i)})

prom = DNAStringSet(unlist(prom))

names(prom) =  gr.promoters$LOC #ojo, debo haberlos corrido 'get_promoter' por orden de Chr
names(prom) = paste0(names(prom), "_promoter")

#check ATG
subseq(prom, start = 1501, end = 1509)


# Write sequences to a fasta file
writeXStringSet(prom, filepath="res/promoters.fasta", format="fasta")  


# hacer estima de frecuencia de bases para crear datos sinteticos 
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






letterFrequency(promoters, c("A", "C"), as.prob = T)
alphabetFrequency(promoters, as.prob = T)[,1:4]

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



# --- Outline ----

# Search for number of occurrences in the actual dataset
# Estimate %GC content
# Estimate for number of occurrences in simulated dataset
# Estimate P value for the n. occurrences in the actual dataset
