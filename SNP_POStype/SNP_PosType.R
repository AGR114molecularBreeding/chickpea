#Dependencies
library(dplyr)
library(stringr)
library(tidyr)
library(GenomicFeatures)
library(rtracklayer)

#Import GFF
unzip("dat/C.arietinum_GenomeGFF.zip", exdir = "dat/")
gff <- import("dat/C.arietinum_GenomeGFF.gff") #Your directory path

#Import SNP information file
Variantes <- read.delim("dat/SNP_poschr.txt") #Your directory path

#Create Dataframe with GFF
gff_df <- data.frame(CHROM = seqnames(gff),
                     start = start(gff),
                     end = end(gff),
                     strand = strand(gff),
                     score = mcols(gff)$score,
                     type = mcols(gff)$type,
                     ID = mcols(gff)$ID)
#Number of different region types
gff_df %>% count(type) 

# Dataframe with SNP positions
snp_df <- data.frame(POS = (Variantes$POS), CHROM = (Variantes$CHROM))

# Unique chromosomes in our data.frame
snp_df <- snp_df %>% filter(grepl("NC_",CHROM))
snp_df <- unique(snp_df)

unique_chromosomes <- unique(snp_df$CHROM)


###--------------------------------------------------------------------------------------------###
################################## Genes ID affected by SNPs #####################################
###--------------------------------------------------------------------------------------------###
# Create a new data.frame to store the results
results_LOCid_df <- data.frame(CHROM = character(),
                               POS = numeric(),
                               ID = character(),
                               STRAND = character(),
                               stringsAsFactors = FALSE)

# Iterate over each chromosome
for (current_chromosome in unique_chromosomes) {
  cat("Verifying positions for the chromosome:", current_chromosome, "\n")
  
  # Filter the data.frame for the current chromosome and types "gene" or "pseudogene"
  gff_chromosome <- gff_df[gff_df$CHROM == current_chromosome & (gff_df$type == "gene" | gff_df$type == "pseudogene"), ]
  
  # Iterate over positions of the current chromosome
  for (POS in snp_df$POS[snp_df$CHROM == current_chromosome]) {
    overlap_plus <- gff_chromosome$start <= POS & gff_chromosome$end >= POS & gff_chromosome$strand == "+"
    overlap_minus <- gff_chromosome$start <= POS & gff_chromosome$end >= POS & gff_chromosome$strand == "-"
    
    if (any(overlap_plus) || any(overlap_minus)) {
      # Get the IDs and strands of the genes if there is overlap
      gene_ids_plus <- gff_chromosome$ID[overlap_plus]
      gene_ids_minus <- gff_chromosome$ID[overlap_minus]
      strands <- character()  # Inicializar un vector vacío para strands
      
      if (any(overlap_plus)) {
        strands <- c(strands, "+")
      }
      
      if (any(overlap_minus)) {
        strands <- c(strands, "-")
      }
      
      # Aggregate IDs based on strand
      gene_ids <- c(gene_ids_plus, gene_ids_minus)
      
      # Add the information to the new data.frame
      # Add the information to the new data.frame
      result_row <- data.frame(CHROM = as.character(current_chromosome), POS = as.numeric(POS), ID = as.character(gene_ids), STRAND = as.character(strands), stringsAsFactors = FALSE)
      results_LOCid_df <- rbind(results_LOCid_df, result_row)
    }
  }
}
mios_result <-results_LOCid_df
print(mios_result)


###--------------------------------------------------------------------------------------------###
################################## SNP_in_two_genes ###################################
###--------------------------------------------------------------------------------------------###
# From the results_LOCid_df file
# Filter the SNPs that are in two different genes (overlapping independently of the strand)
# Remove rows that have the same coordinates (POS) 
SNP_varios_genees <- results_LOCid_df[duplicated(results_LOCid_df$POS) | duplicated(results_LOCid_df$POS, fromLast = TRUE), ]
print(SNP_varios_genees)

# Now, check if there are any SNPs in two or more genes overlapping on the negative strand
# First, select rows with duplicated coordinates on the negative strand
línea_POS_neg <- results_LOCid_df[duplicated(results_LOCid_df$POS) & results_LOCid_df$STRAND == '-' | duplicated(results_LOCid_df$POS, fromLast = TRUE) & results_LOCid_df$STRAND == '-' , ]
print(línea_POS_neg)
SNP_solape_neg <- línea_POS_neg[duplicated(línea_POS_neg$POS) & línea_POS_neg$STRAND == '-' | duplicated(línea_POS_neg$POS, fromLast = TRUE) & línea_POS_neg$STRAND == '-' , ]
print(SNP_solape_neg)

# Same process for the positive strand
línea_POS_posi <- results_LOCid_df[duplicated(results_LOCid_df$POS) & results_LOCid_df$STRAND == '+' | duplicated(results_LOCid_df$POS, fromLast = TRUE) & results_LOCid_df$STRAND == '+' , ]
print(línea_POS_posi)
SNP_solape_posi <- línea_POS_posi[duplicated(línea_POS_posi$POS) & línea_POS_posi$STRAND == '+' | duplicated(línea_POS_posi$POS, fromLast = TRUE) & línea_POS_posi$STRAND == '+' , ]
print(SNP_solape_posi)

# MERGE THE SNP OVERLAP DATAFRAMES IN TWO GENES THAT OVERLAP AND ARE ON THE SAME STRAND

# Merge the SNP_overlap_pos and SNP_overlap_neg dataframes
SNP_solape_total <- rbind(SNP_solape_posi, SNP_solape_neg)

# Print the result
print(SNP_solape_total)


## Create a list to group lines with the same genes into groups
# Create a list of dataframes grouped by ID
grupos <- split(SNP_solape_total, gsub("gene-", "", SNP_solape_total$ID))

# Create different dataframes with descriptive names (without "gene-")
for (i in seq_along(grupos)) {
  assign(paste0("snp_", names(grupos)[i]), grupos[[i]], envir = .GlobalEnv)
}

print(grupos)


# SNP_overlap_total will be used for the filters #

###--------------------------------------------------------------------------------------------###
### Extract the 2 genes that overlap on the same strand and are in the previous SNPs ###
###--------------------------------------------------------------------------------------------###
# Create a vector of unique genes from SNP_overlap_total
vector_interest_genes <- SNP_solape_total %>% distinct(ID)
print(vector_interest_genes)

# Create a list to store the dataframes
genes_list <- list()

# Create a loop to process each gene in the file
for (i in 1:nrow(vector_interest_genes)) {
  
  # Get information about the current gene
  current_gene_id <- vector_interest_genes$ID[i]
  
  cat("Processing gene:", current_gene_id, "\n")
  
  # Find the index of the gene of interest in the GFF
  index_interest_gene <- which(gff_df$ID == current_gene_id)
  
  # Find the index of the first line containing "gene" after the gene of interest
  index_next_gene <- which(grepl("gene", gff_df$ID[(index_interest_gene + 1):nrow(gff_df)]))[1]
  
  # Check if a valid index was found
  if (length(index_next_gene) == 0) {
    cat("No line containing 'gene' was found after the gene of interest.\n")
  } else {
    # Calculate the absolute index of the line with "gene" in the ID column
    index_next_gene_abs <- index_interest_gene + index_next_gene
    
    # Filter lines from the gene of interest to the line with "gene" in the ID column (excluding)
    current_gene_gff <- as.data.frame(gff_df[index_interest_gene:(index_next_gene_abs - 1), ])    
    # Print information about the current gene
    print(current_gene_gff)
    
    # Store the gene in the list with a descriptive name
    id_without_prefix <- gsub("gene-", "", current_gene_id)
    genes_list[[id_without_prefix]] <- current_gene_gff
    
    cat("Gene saved in the list:", id_without_prefix, "\n")
  }
  
  cat("\n")
}

# Print the names of the created dataframes
print(names(genes_list))
# Print specific dataframe(s)
print(genes_list$LOC101504401)
print(genes_list$LOC101504943)

## WITH THE GENES ALREADY IN DIFFERENT DATAFRAMES, PROCEED TO EXAMINE THE GENETIC REGION OF THE SNPS IN PARTS ##
## FIRST, THE SPECIAL SNPS ##


###--------------------------------------------------------------------------------------------###
############################### SNPs Region Type Overlapep Genes #################################
##### mRNA (CDS, Exon or Intron); RNA_type: lnc_RNA, rRNA, snoRNA, snRNA, tRNA; pseudogene #######
###--------------------------------------------------------------------------------------------###
# Create a new data.frame to store the results
results_df <- data.frame(CHROM = character(),
                         POS = numeric(),
                         TYPE = character(),
                         ID = character(),
                         stringsAsFactors = FALSE)

# Iterate over each grupo in grupos
for (grupo in names(grupos)) {
  cat("Checking positions for gene:", grupo, "\n")
  
  # Extract relevant data from the lists
  snp_data <- grupos[[grupo]]
  gff_data <- genes_list[[grupo]]
  
  # Iterate over positions in the current grupo
  for (POS in snp_data$POS) {
    overlaps_plus <- gff_data$start <= POS & gff_data$end >= POS & gff_data$strand == "+"
    overlaps_minus <- gff_data$start <= POS & gff_data$end >= POS & gff_data$strand == "-"
    
    tipo_plus <- if (any(overlaps_plus)) {
      if ("CDS" %in% gff_data$type[overlaps_plus]) {
        "CDS_plus"
      } else if ("mRNA" %in% gff_data$type[overlaps_plus]&& "exon" %in%gff_data$type[overlaps_plus]) {
        "exon_plus"
      } else if ("intron" %in% gff_data$type[overlaps_plus]) {
        "intron_plus"
      } else if ("lncRNA" %in% gff_data$type[overlaps_plus]) {
        "lncRNA_plus"
      } else if ("tRNA" %in% gff_data$type[overlaps_plus]) {
        "tRNA_plus"
      } else if ("snoRNA" %in% gff_data$type[overlaps_plus]) {
        "snoRNA_plus"
      } else if ("snRNA" %in% gff_data$type[overlaps_plus]) {
        "snRNA_plus"
      } else if ("rRNA" %in% gff_data$type[overlaps_plus]) {
        "rRNA_plus"
      } else if ("pseudogene" %in% gff_data$type[overlaps_plus]) {
        "pseudogene_plus"
      } else if ("transcript" %in% gff_data$type[overlaps_plus]) {
        "transcript_plus"
      } else if ("gene" %in% gff_data$type[overlaps_plus]) {
        "gene_plus"
      } else {
        "NA_plus"  # Does not meet any condition on the plus strand = INTERGENIC
      }
    } else {
      tipo_plus <- "NA_plus"
    }
    
    tipo_minus <- if (any(overlaps_minus)) {
      if ("CDS" %in% gff_data$type[overlaps_minus]) {
        "CDS_minus"
      } else if ("mRNA" %in% gff_data$type[overlaps_minus]&& "exon" %in%gff_data$type[overlaps_minus]) {
        "exon_minus"
      } else if ("intron" %in% gff_data$type[overlaps_minus]) {
        "intron_minus"
      } else if ("lncRNA" %in% gff_data$type[overlaps_minus]) {
        "lncRNA_minus"
      } else if ("tRNA" %in% gff_data$type[overlaps_minus]) {
        "tRNA_minus"
      } else if ("snoRNA" %in% gff_data$type[overlaps_minus]) {
        "snoRNA_minus"
      } else if ("snRNA" %in% gff_data$type[overlaps_minus]) {
        "snRNA_minus"
      } else if ("rRNA" %in% gff_data$type[overlaps_minus]) {
        "rRNA_minus"
      } else if ("pseudogene" %in% gff_data$type[overlaps_minus]) {
        "pseudogene_minus"
      } else if ("transcript" %in% gff_data$type[overlaps_minus]) {
        "transcript_minus"
      } else if ("gene" %in% gff_data$type[overlaps_minus]) {
        "gene_minus"
      } else {
        "NA_minus"  # Does not meet any condition on the minus strand = INTERGENIC
      }
    } else {
      tipo_minus <- "NA_minus"
    }
    
    # Process according to the specific conditions of each strand
    tipo <- if (tipo_plus != "NA_plus" && tipo_minus != "NA_minus") {
      paste(tipo_plus, tipo_minus, sep="_and_")
    } else if (tipo_plus != "NA_plus") {
      tipo_plus
    } else {
      tipo_minus
    }
    
    # Get the ID corresponding to the current position
    current_chrom <- snp_data$CHROM[snp_data$POS == POS]
    current_id <- snp_data$ID[snp_data$POS == POS]
    
    # Add information to the new data.frame
    results_df <- rbind(results_df, data.frame(CHROM = current_chrom, POS = POS, TYPE = tipo, ID = current_id))
  }
}

# Print the results
SNP_raritos <- results_df
print(SNP_raritos)
SNP_raritos <- SNP_raritos %>% 
  mutate(strand = ifelse(str_detect(TYPE, "_plus"), "+", ifelse(str_detect(TYPE, "_minus"), "-", NA)))
print(SNP_raritos)

## FROM THE SNP FILE, REMOVE THE SPECIAL SNPS (USED BEFORE) AND SEARCH FOR REGIONS ##
############################### Remove the peculiar SNPs and find regions ###########################

# Remove peculiar SNPs from the snp file
snp_sinraros <- snp_df[!(snp_df$CHROM %in% SNP_solape_total$CHROM & snp_df$POS %in% SNP_solape_total$POS), ]


# Create a new data.frame to store the results
results_df <- data.frame(CHROM = character(),
                         POS = numeric(),
                         TYPE = character(),
                         stringsAsFactors = FALSE)

# Iterate over each chromosome
for (current_chromosome in unique_chromosomes) {
  cat("Checking positions for chromosome:", current_chromosome, "\n")
  
  # Filter the data.frame for the current chromosome
  gff_chromosome <- gff_df[gff_df$CHROM == current_chromosome, ]
  
  # Iterate over positions of the current chromosome
  for (POS in snp_sinraros$POS[snp_sinraros$CHROM == current_chromosome]) {
    overlaps_plus <- gff_chromosome$start <= POS & gff_chromosome$end >= POS & gff_chromosome$strand == "+"
    overlaps_minus <- gff_chromosome$start <= POS & gff_chromosome$end >= POS & gff_chromosome$strand == "-"
    
    tipo_plus <- if (any(overlaps_plus)) {
      if ("CDS" %in% gff_chromosome$type[overlaps_plus]) {
        "CDS_plus"
      } else if ("mRNA" %in% gff_chromosome$type[overlaps_plus] && "exon" %in% gff_chromosome$type[overlaps_plus]) {
        "exon_plus"
      } else if ("mRNA" %in% gff_chromosome$type[overlaps_plus]) {
        "intron_plus"
      } else if ("lnc_RNA" %in% gff_chromosome$type[overlaps_plus]) {
        "lncRNA_plus"
      } else if ("tRNA" %in% gff_chromosome$type[overlaps_plus]) {
        "tRNA_plus"
      } else if ("snoRNA" %in% gff_chromosome$type[overlaps_plus]) {
        "snoRNA_plus"
      } else if ("snRNA" %in% gff_chromosome$type[overlaps_plus]) {
        "snRNA_plus"
      } else if ("rRNA" %in% gff_chromosome$type[overlaps_plus]) {
        "rRNA_plus"
      } else if ("pseudogene" %in% gff_chromosome$type[overlaps_plus]) {
        "pseudogene_plus"
      } else if ("transcript" %in% gff_chromosome$type[overlaps_plus]) {
        "transcript_plus"
      } else if ("gene" %in% gff_chromosome$type[overlaps_plus]) {
        "gene_plus"
      } else {
        "NA_plus"  # Does not meet any condition on the plus strand = INTERGENIC
      }
    } else {
      tipo_plus <- "NA_plus"
    }
    
    tipo_minus <- if (any(overlaps_minus)) {
      if ("CDS" %in% gff_chromosome$type[overlaps_minus]) {
        "CDS_minus"
      } else if ("mRNA" %in% gff_chromosome$type[overlaps_minus] && "exon" %in% gff_chromosome$type[overlaps_minus]) {
        "exon_minus"
      } else if ("mRNA" %in% gff_chromosome$type[overlaps_minus]) {
        "intron_minus"
      } else if ("lnc_RNA" %in% gff_chromosome$type[overlaps_minus]) {
        "lncRNA_minus"
      } else if ("tRNA" %in% gff_chromosome$type[overlaps_minus]) {
        "tRNA_minus"
      } else if ("snoRNA" %in% gff_chromosome$type[overlaps_minus]) {
        "snoRNA_minus"
      } else if ("snRNA" %in% gff_chromosome$type[overlaps_minus]) {
        "snRNA_minus"
      } else if ("rRNA" %in% gff_chromosome$type[overlaps_minus]) {
        "rRNA_minus"
      } else if ("pseudogene" %in% gff_chromosome$type[overlaps_minus]) {
        "pseudogene_minus"
      } else if ("transcript" %in% gff_chromosome$type[overlaps_minus]) {
        "transcript_minus"
      } else if ("gene" %in% gff_chromosome$type[overlaps_minus]) {
        "gene_minus"
      } else {
        "NA_minus"  # Does not meet any condition on the minus strand = INTERGENIC
      }
    } else {
      tipo_minus <- "NA_minus"
    }
    
    # Process according to the specific conditions of each strand
    tipo <- if (tipo_plus != "NA_plus" && tipo_minus != "NA_minus") {
      paste(tipo_plus, tipo_minus, sep="_and_")
    } else if (tipo_plus != "NA_plus") {
      tipo_plus
    } else {
      tipo_minus
    }
    
    # Add information to the new data.frame
    results_df <- rbind(results_df, data.frame(CHROM = current_chromosome, POS = POS, TYPE = tipo))
  }
}
print(results_df)

## RESULTS THAT HAVE OVERLAP WITH GENES ON BOTH POSITIVE AND NEGATIVE STRANDS HAVE "_and_"
## PROCEED TO SEPARATE BOTH RESULTS SO THAT THEY ARE IN DIFFERENT ROWS
results_df <- results_df %>% 
  separate_rows(TYPE, sep = "_and_")
results_df <- results_df %>% 
  mutate(strand = ifelse(str_detect(TYPE, "_plus"), "+", ifelse(str_detect(TYPE, "_minus"), "-", NA)))
print(results_df)

## MERGE THE DATAFRAME WITH RARE SNPS AND THE ONE WITH "NORMAL" SNPS

# To have the genes of rare snps in the dataframe, merge results_df and mine
# mine has the genes in which the snps are found
# results_df, on the other hand, has the genetic regions of the non-rare snps

# Perform the combination using merge
SNP_normal_F <- merge(results_df, mios_result,  by.x = c("POS", "CHROM","strand"), by.y = c("POS", "CHROM","STRAND" ))
# Obtained a dataframe with the genes of the normal snps

# Print the result
print(SNP_normal_F)

#################Unir ambas ###################
#SNP_rarito_F y SNP_normal_F

SNP_combined <- rbind(SNP_normal_F, SNP_raritos)
print(SNP_combined)


SNP_combined <- SNP_combined %>%
  mutate(TYPE = str_remove(TYPE, "_.*")) %>% 
  mutate(ID = str_remove(ID, "gene-"))

SNP_combined <- SNP_combined %>% 
  arrange(CHROM, POS)

###--------------------------------------------------------------------------------------------###
################################## CLEAN final data.frame  #######################################
###--------------------------------------------------------------------------------------------###

##---Show FINAL_results
SNP_combined %>% count(TYPE)
SNP_combined %>% count(ID)

#SNP with "NA" -- intergenic SNP
SNP_combined %>% filter(TYPE %in% "NA")

#SNPs in genes
SNP_combined %>% filter(!TYPE %in% "NA") %>% count(POS)

#Non-coding genes
SNP_combined %>% filter(TYPE %in% c("tRNA", "lncRNA", "snoRNA", "snRNA", "rRNA")) %>% count(TYPE, sort = TRUE)
SNP_combined %>% filter(TYPE %in% c("tRNA", "lncRNA", "snoRNA", "snRNA", "rRNA")) %>% group_by(TYPE) %>% count(ID, sort = TRUE)                                                                             

#Coding genes
SNP_combined %>% filter(TYPE %in% c("CDS","exon")) %>% count(CHROM)
SNP_combined %>% filter(TYPE %in% c("CDS","exon")) %>% distinct(POS, .keep_all = TRUE) %>% count(CHROM)

SNP_combined %>% filter(TYPE %in% c("CDS","exon")) %>% filter(duplicated(POS) | duplicated(POS, fromLast = TRUE))

SNP_combined %>% filter(TYPE %in% c("CDS","exon")) %>% count(TYPE)

#pseudogenes and misc_RNA
SNP_combined %>% filter(TYPE %in% c("transcript", "pseudogene")) %>% count(TYPE)

#SNPs in CDS or exon
SNP_combined %>% filter(TYPE %in% c("CDS","exon"))



# Save final data
write.table(SNP_combined, "SNP_SpecificRegionType.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

