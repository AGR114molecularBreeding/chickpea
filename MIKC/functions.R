# MADS-box.R -- a set of common functions for analysis of MADS-box family in Ca
# Copyright (C) 2024 Jose V. Die  <jose.die@uco.es>
# Distributed under terms of the MIT license.


get_xp_length <- function(xp_id){
  ### ''' Estimate protein length (aa)
  
  xpinfo <- entrez_summary(db = "protein", id = xp_id)
  xpinfo$slen
  
}

blast_homologs <- function(hitFile, ident, cover) {
  
  ### ''' Take a csv table with Blast hits and return a 
  ### subset filteres by %Identity & %Coverage (longer sequence)
  
  # read downloaded Hit file
  hit = read.csv(hitFile, header = FALSE, stringsAsFactors = FALSE)
  
  # change colnames
  colnames(hit) = c("query", "subject", "identity", "align_length", "mismatches",
                    "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue",
                    "bit_score", "%positives")
  
  # make a table_length for queries & subjects 
  queries <- hit %>% 
    pull(query) %>% 
    unique()
  
  tlengths_query = tibble(query = queries, len.query = sapply(query, function(i) get_xp_length(i)))
  
  subjects <- hit %>% 
    pull(subject) %>% 
    unique()
  
  tlengths_sub = tibble(subject = subjects, len.sub = sapply(subject, function(i) get_xp_length(i)))
  
  
  # combine hit_table with queries & subjects lengths 
  # estimate coverage (longer sequence)
  # apply filters : identity * coverage
  hit %>% 
    tibble() %>% 
    select(-c(align_length, bit_score, "%positives", mismatches, gap_opens)) %>% 
    right_join(tlengths_query, by = "query") %>% 
    right_join(tlengths_sub,   by = "subject") %>% 
    mutate(cov.query = (q.end - q.start + 1)/len.query, 
           cov.sub = (s.end - s.start + 1)/len.sub) %>% 
    mutate(cov_long = case_when(len.query > len.sub ~ cov.query*100, 
                                len.query < len.sub ~ cov.sub*100, 
                                TRUE ~ cov.query*100)) %>% 
    select(-c(q.start, q.end, s.start, s.end, len.query, len.sub, cov.query, cov.sub)) %>% 
    filter(identity >= ident, cov_long >= cover) %>% 
    select(query, subject, identity, cov_long, evalue) %>% 
    arrange(identity)
  
  
}
 


 
blast_best_homolog <- function(hitFile) {
  
  ### ''' Take a csv table with Blast hits and return a 
  ### subset containing only the best hit for each query '''
  
  # Read Downloaded Hit file
  hit = read.csv(hitFile, header = FALSE, stringsAsFactors = FALSE)
  
  # Change colnames
  colnames(hit) = c("query", "subject", "identity", "align_length", "mismatches",
                    "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue",
                    "bit_score", "%positives")
  
  res <- as_tibble(hit) %>% 
    filter(query != subject) %>% 
    group_by(query) %>% 
    slice(1) %>% 
    ungroup()
  
  res
  
}

## Usage
## Ex. multiple blastp searches (with 2 sequences)
## best_homolog("B1J2TDTZ01R-Alignment-HitTable.csv"





characterizeTable <- function(targets) {
  
  # targets, character vector with XM protein accession ids
  
  ## Initialize the vectors containing the features that we need for each sequence
  LOC    <-  vector("character") 
  Chr    <-  vector("character")
  chr_s  <-  vector("integer")
  chr_e  <-  vector("integer")
  exon   <-  vector("integer") 
  AA     <-  vector("integer")
  mol_wt <-  vector("integer")
  xm     <-  vector("integer")
  xp     <-  vector("integer")
  
  ## Extract feature info from NCBI
  for(xi in as.character(targets)) {
    xp = c(xp, xi)
    xm = c(xm, getXM(xi))
    
    xpinfo <- entrez_summary(db = "protein", id = xi)
    AA = c(AA, xpinfo$slen)
    
    protein_gb <- entrez_fetch(db = "protein", id = xi, rettype = "gp")
    mol_wt = c(mol_wt, extract_mol.wt_from_xp(protein_gb)/1000) # kDa
    
    xplink <- entrez_link(dbfrom = "protein", id = xi, db = "gene")
    genesummary = entrez_summary(db = "gene", id = xplink$links[1])
    LOC = c(LOC, paste0("LOC",genesummary$uid))
    Chr = c(Chr, genesummary$chromosome)
    
    # defensive programming
    # some XP may not have chrstart, chrstop or exon count
    if(length(genesummary$genomicinfo$chrstart) != 0) {
      chr_s = c(chr_s, genesummary$genomicinfo$chrstart)
      chr_e = c(chr_e, genesummary$genomicinfo$chrstop)
      exon = c(exon, genesummary$genomicinfo$exoncount)  
      
    } else {
      chr_s = c(chr_s, 0)
      chr_e = c(chr_e, 0)
      exon = c(exon, 0)  
      
    }
  }
  
  ## Build the dataset 
  t1 <- data.frame(LOC, XP = xp, Chr, chr_s, chr_e, AA, mol_wt, XM = xm, exon)
  
  ## Add strand
  t1 <- t1 %>% mutate(Strand = ifelse(chr_e > chr_s, "+", "-"))
  
  ## Sort data set by Chr + start coordinate
  t1 = t1 %>% arrange(Chr, chr_s)
  
  t1
  
}

## Usage
##targets <- c("XP_004498126.1")
##targets <- c("XP_004513719.1", "XP_004485955.1", "XP_004498126.1", "XP_004498126.1")
##tdat <- characterizeTable(targets)



TSScoordinates <- function(gr, CDSstringset, chr) {
  
  CDSstringset = unique(CDSstringset)
  gr.tss = gr[seqnames(gr) == chr]
  
  for(i in seq_along(gr.tss)) {  
    text = genome[[which(names(genome) == chr)]]
    locid = gr.tss$LOC[i]
    
    locid = ifelse(locid=="LOC101488241", "CAAP1", locid) #customize for this study
    pattern = CDSstringset[which(names(CDSstringset) == locid)][[1]][1:gr.tss$bp_cut[i]]
    
    if(gr.tss[i] %in% gr.tss[strand(gr.tss) == "+"]) {
      
      s = start(matchPattern(pattern, text, max.mismatch = 0))
      if(gr.tss$LOC[i] == "LOC101493118") {s = s[1]} #customize for this study
      start(gr.tss[i]) = s
      
    }
    
    if(gr.tss[i] %in% gr.tss[strand(gr.tss) == "-"]) {#do not use 'else' bc some genes might be unmapped (strand = *)
      
      pattern = reverseComplement(pattern)
      e = end(matchPattern(pattern, text, max.mismatch = 0))
      if(gr.tss$LOC[i] == "LOC101509413") {e = e[655749]} #customize for this study
      end(gr.tss[i]) = e
      
    }
    
  }
  
  
  gr.tss  
  
  
}


## Usage
## gr.tss = GRangesList(TSScoordinates(gr, mads_cds, "Ca1"),
##                      TSScoordinates(gr, mads_cds, "Ca3"), 
##                      TSScoordinates(gr, mads_cds, "Ca4"))



names2LOC <- function(str) {
  # Take a vector string (header from multifasta file)
  # Return the LOC id in a tidy format 
  str_split(str_split(str, "\\(")[[1]][2], "\\)")[[1]][1]
}

## Usage 
## mads_cds = readDNAStringSet("mads_cds.fasta")
## my_names = names(mads_cds)
## my_names <- sapply(my_names, function(i) names2LOC(i), USE.NAMES = F)


names2XM <- function(str) {
  # Take a vector string (header from multifasta file) 
  # Return the XM id in a tidy format
  
  # str, a string vector 
  
  str_split(str, "\\:")[[1]][1]  
}

get_before_period <- function(str) {
  
  # Return a string vector without the characters 
  # before a period (excluding the period)
  
  # str, a string vector 
  
  str_split(str[1], "\\.")[[1]][1]
  
}


get_promoters <- function(gr, chr = "Ca3") {
  
  inputGR = gr[seqnames(gr) == chr]  
  text = genome[[which(seqnames(genome) == chr)]]
  promoter.seq = c()
  
  for(i in seq_along(inputGR)) {
    seq = text[start(inputGR[i]):end(inputGR[i])]
    
    if(inputGR[i] %in% inputGR[strand(inputGR) == "+"]) {
      promoter.seq = c(promoter.seq,seq)
    }
    if(inputGR[i] %in% inputGR[strand(inputGR) == "-"]) {
      seq = reverseComplement(seq)
      promoter.seq = c(promoter.seq, seq)
    }
    
  }
  
  promoter.seq
  
}
