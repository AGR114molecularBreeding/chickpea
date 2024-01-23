# Dependencies
library(readxl)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)


# Define a function
fpkm_to_tpm <- function(fpkm_dat){
  
  # ''' take a dataset of FPKM values
  
  # ''' First column with Gene ids. 
  # ''' next columns with samples
  
  
  
  pivot_longer(-names(fpkm_dat)[1], names_to = "tissue", values_to = "fpkm") %>% # reshape the data 
  
  group_by(tissue) %>% 
  
  mutate(total_fpkm_per_sample = sum(fpkm),                     # sum of FPKM values per sample
           scaling_factor  = total_fpkm_per_sample/ 1e6,        # scaling factor per sample
           tpm_values = fpkm / scaling_factor) %>%              # calculate TPM values
  
  select(names(fpkm_dat)[1], tissue, tpm_values) %>% 
  pivot_wider(names_from = "tissue", values_from = tpm_values)   # return tpm values in the original data shape
} 


# Load data 
load("dat/edata.rda")

# Pre-processing : convert from FPKM to TPM 

## select samples (in columns) 
edata <- as_tibble(edata) %>% 
  select('Transcript id', 'Genomic Coordinates', mRNA, FB1:FB4, FL1:FL5, YL, ML)

## prepare dataset with gene id in 1st column + samples in the rest of columns 
fpkm_dat <- edata %>% 
  select(-c(2,3))

# Calculate TPM values from FPKM using the function `fpkm_to_tpm`
tpm_dat <- fpkm_to_tpm(fpkm_dat)
tpm_dat


# Calculate mean in BUDS, FLOWERS, LEAVES for each gene
tpm_long <- tpm_dat %>% 
  pivot_longer(-1, names_to = "sample", values_to = "tpm") %>% 
  mutate(tissue = dplyr::case_match(sample, 
                                    paste0('FB',1:4) ~ 'FB',
                                    paste0('FL',1:5) ~ 'FL', 
                                    'YL' ~ 'leaves', 
                                    'ML' ~ 'leaves'))
tpm_bytissue <- tpm_long%>% 
  group_by(`Transcript id`, tissue) %>% 
  summarize(mean = mean(tpm), 
            sd = sd(tpm))
  
  
# Pair-wise comparison : Flowers vs flower buds ------------------------------------------------------

# Flowers vs flower buds - specific from flowers (buds = 0)
transcript_0bud = tpm_bytissue %>% 
  filter(tissue == 'FB' & mean == 0) %>% 
  pull(`Transcript id`)

top_20 <- tpm_bytissue %>% 
  filter(! `Transcript id` %in% transcript_0bud) %>% 
  pivot_wider(names_from = tissue, values_from = c(mean, sd)) %>% 
  select(`Transcript id`, mean_FL, sd_FL) %>% 
  arrange(desc(mean_FL)) %>% 
  head(20)

top_20

# Flowers vs flower buds - specific from buds (flowers = 0)
transcript_0flower = tpm_bytissue %>% 
  filter(tissue == 'FL' & mean == 0) %>% 
  pull(`Transcript id`)

top_20 <- tpm_bytissue %>% 
  filter(! `Transcript id` %in% transcript_0flower) %>% 
  pivot_wider(names_from = tissue, values_from = c(mean, sd)) %>% 
  select(`Transcript id`, mean_FL, sd_FL) %>% 
  arrange(desc(mean_FL)) %>% 
  head(20)

top_20


# Flowers vs flower buds - upregulated either in flowers or buds
flower_bud <- tpm_bytissue %>% 
  filter(! `Transcript id` %in% transcript_0flower) %>% 
  filter(! `Transcript id` %in% transcript_0bud) %>% 
  pivot_wider(names_from = tissue, values_from = c(mean, sd))

# set a minimum threshold of meanTPM = 5
flower_bud  %>% 
  select(`Transcript id`, mean_FL, sd_FL, mean_FB, sd_FB) %>% 
  filter(mean_FL >= 5, mean_FB >=5) %>% 
  mutate(ratio = mean_FL/mean_FB) %>% 
  arrange(desc(ratio)) %>% 
  head(20) # higher expression in FLOWERS

flower_bud  %>% 
  select(`Transcript id`, mean_FL, sd_FL, mean_FB, sd_FB) %>% 
  filter(mean_FL >= 5, mean_FB >=5) %>% 
  mutate(ratio = mean_FL/mean_FB) %>% 
  arrange(ratio) %>%  
  head(20) # higher expression in BUDS

# make plots of representative genes 
tpm_long %>% 
  filter(`Transcript id` %in% c("Ca_TC13690", "Ca_TC23853")) %>% 
  filter(tissue %in% c('FB', 'FL')) %>% 
  ggplot(aes(x = tissue, y = tpm, col = sample)) + 
  geom_jitter(width = 0.1, alpha = 0.8) + 
  ylab("TPM") + xlab("") + 
  facet_wrap(~ `Transcript id`, scales = "free_y")

# We have observed that the TPM values of numerous genes in FL5 closely resemble those in FB values. Notably, it is intriguing that in late flower stages, 
# the expression experiences a decline. These particular observations will be excluded from the analysis for the general plot but are deemed noteworthy for 
# a time-course plot within the FL category.

# Filter out : FL5
tpm_long %>% 
  filter(`Transcript id` %in% c("Ca_TC13690", "Ca_TC23853")) %>% 
  filter(tissue %in% c('FB', 'FL')) %>% 
  filter(!sample == "FL5") %>%
  ggplot(aes(x = tissue, y = tpm)) + 
  geom_jitter(width = 0.1, alpha = 0.8, col = "steelblue") + 
  ylab("TPM") + xlab("") + 
  facet_wrap(~ `Transcript id`, scales = "free_y")



# Pair-wise comparison : Flower buds vs leaves -------------------------------------------------------
# Flower buds vs leaves - specific from leaves (buds = 0)
transcript_0bud = tpm_bytissue %>% 
  filter(tissue == 'FB' & mean == 0) %>% 
  pull(`Transcript id`)

top_20 <- tpm_bytissue %>% 
  filter(! `Transcript id` %in% transcript_0bud) %>% 
  pivot_wider(names_from = tissue, values_from = c(mean, sd)) %>% 
  select(`Transcript id`, mean_leaves, sd_leaves) %>% 
  arrange(desc(mean_leaves)) %>% 
  head(20)

top_20


# Flowers vs flower buds - specific from buds (leaves = 0)
transcript_0leaves = tpm_bytissue %>% 
  filter(tissue == 'leaves' & mean == 0) %>% 
  pull(`Transcript id`)

top_20 <- tpm_bytissue %>% 
  filter(! `Transcript id` %in% transcript_0leaves) %>% 
  pivot_wider(names_from = tissue, values_from = c(mean, sd)) %>% 
  select(`Transcript id`, mean_FB, sd_FB) %>% 
  arrange(desc(mean_FB)) %>% 
  head(20)

top_20

# Flowers vs flower buds - upregulated in buds / leaves
buds_leaves <- tpm_bytissue %>% 
  filter(! `Transcript id` %in% transcript_0leaves) %>% 
  filter(! `Transcript id` %in% transcript_0bud) %>% 
  pivot_wider(names_from = tissue, values_from = c(mean, sd))

# minimum threshold of meanTPM = 5
buds_leaves  %>% 
  select(`Transcript id`, mean_FB, sd_FB, mean_leaves, sd_leaves) %>% 
  filter(mean_leaves >= 5, mean_FB >=5) %>% 
  mutate(ratio = mean_FB/mean_leaves) %>% 
  arrange(desc(ratio)) %>% 
  head(20)

buds_leaves  %>% 
  select(`Transcript id`, mean_FB, sd_FB, mean_leaves, sd_leaves) %>% 
  filter(mean_leaves >= 5, mean_FB >=5) %>% 
  mutate(ratio = mean_FB/mean_leaves) %>% 
  arrange(ratio) %>% 
  head(20)

# make plots of representative genes 
tpm_long %>% 
  filter(`Transcript id` %in% c("Ca_TC05559", "Ca_TC18116")) %>% 
  filter(tissue %in% c('FB', 'leaves')) %>% 
  ggplot(aes(x = tissue, y = tpm, col = sample)) + 
  geom_jitter(width = 0.1, alpha = 0.8) + 
  ylab("TPM") + xlab("") + 
  facet_wrap(~ `Transcript id`, scales = "free_y")
