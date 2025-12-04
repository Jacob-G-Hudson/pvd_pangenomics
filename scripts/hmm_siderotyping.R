# script to take final hmms, and siderotype strains using them
## TODO I have redone the pangnome typing since I have made the hmms. This measn that the hmm name
## Does not necessarily match the loci name, especially for the HP_genes
## I need to either translate the name or remake the hmm
# dependencies ---
library(tidyverse)  # data wrangling
library(rhmmer)     # utilities for hmmer files, from https://github.com/arendsee/rhmmer

# setup ----
annotations_file            <- "./data/summaries/annotations_v2.rds"
hmmer_dir                   <- "./data/jackhmmer/20251112_hmms_5"
thresholds_file             <- "./data/jackhmmer/20251112_hmms_5/hmm_thresholds.csv"
pangenome_siderotypes_file  <- "./data/summaries/siderotype_strain_predictions.csv"

# load files ---- 
all_loci              <- read_rds(annotations_file)
all_thresholds        <- read_csv(thresholds_file)
pangenome_siderotypes <- read_csv(pangenome_siderotypes_file)


hmmer_output_files <- list.files(path       = hmmer_dir,
                                 pattern    = "^results.out$",
                                 full.names = TRUE,
                                 recursive  = TRUE)

all_hmm_results <- lapply(X   = as.list(hmmer_output_files),
                          FUN = read_tblout) %>%
  do.call(rbind, .) %>% # here we have 276701
  transmute(Locus_Tag = domain_name,
            query_name= query_name,
            evalue    = sequence_evalue,
            bitscore  = sequence_score) %>%
  group_by(Locus_Tag) %>%
  # for loci with multiple annotations, take the best match
  filter(evalue   == min(evalue))   %>% # here 87378
  filter(bitscore == max(bitscore)) %>% # here 75268
  # filter based on thresholds %>%
  left_join(all_thresholds,
            by = c("query_name" = "hmm")) %>%
  filter(bitscore >= threshold) # here 2312
  
# do we have loci with multiple annotations
sum(duplicated(all_hmm_results$Locus_Tag))

# identify pvd region -----
hmm_pvd_region <- all_hmm_results %>%
  # add in pangenome-based annotation
  # TODO what if I have not at all annotated something in the region?
  left_join(all_loci %>%
              select(Locus_Tag, HOG_annotation_4),
            by = "Locus_Tag") %>%
  # get to one annotation per locus
  group_by(Locus_Tag) %>%
  summarise(HOG_annotation_4 = unique(HOG_annotation_4),
            hmm_annotation   = paste(query_name, collapse = "/")) %>%
  # filter to only adjacent loci (should be in pvd region)
  separate_wider_delim(cols      = Locus_Tag,
                       delim     = "_",
                       names     = c("bacteria_code", "locus"),
                       too_many  = "merge",
                       cols_remove = FALSE) %>%
  mutate(locus = as.numeric(locus)) %>%
  ungroup() %>%
  group_by(bacteria_code) %>%
  arrange(Locus_Tag) %>%
  mutate(locus_lead_dif      = lead(locus) - locus,
         locus_lag_dif       = locus - lag(locus)) %>%
  ungroup() %>%
  rowwise() %>%
  filter(any(locus_lead_dif == 5,
             locus_lag_dif == 5,
             na.rm = TRUE)) %>%
  select(-c(locus_lead_dif, locus_lag_dif))

# Compare annotation accuracy with pangenome annotations ----
hmm_annotation_accuacty <- hmm_pvd_region %>%
  # remove cip strains as we don't have pangenome annotation for them
  filter(!grepl("CIP", bacteria_code)) %>%  
  # to do work out a summary statistic for this stage such as senstivity/specificity/accuracy
  mutate(agreement = case_when(HOG_annotation_4 == hmm_annotation ~ "Agree",
                               .default = "Disagree")) %>%
  group_by(agreement) %>%
  summarise(n = n(),
            per = n/(2269) * 100)
  
# hmm-based pvd_predictions ----
## Use same pvd_typing function as for pangenome
hmm_pvd_prediction <- hmm_pvd_region %>%
  group_by(bacteria_code) %>% 
  summarise(pvd_region = paste(hmm_annotation, collapse = ", ")) %>%
  rowwise() %>%
  mutate(region_stats = pvd_typing(region    = pvd_region,
                                   core      = NULL,
                                   branch    = list(pvd1 = c("pvdE_1", "fpvA_1", "pvdD_1", "pvdJ_1", "pvdI_1"),
                                                    pvd2 = c("HP_586", "HP_585", "aes_1", "solG_1", "dltA_3"),
                                                    pvd3 = c("yojI_1", "fpvA_2", "dltA_1", "HP_1", "dltA_2", "sip_1")))) %>%
  separate_wider_delim(region_stats, delim = ", ", names_sep = "_") %>%
  pivot_longer(cols = starts_with("region_stats"),
               values_to = "pair") %>%
  select(-name) %>%
  separate_wider_delim(pair,
                       delim = "=",
                       names = c("key", "value")) %>%
  pivot_wider(names_from = "key",
              values_from = "value") %>%
  transmute(bacteria_code          = bacteria_code,
            hmm_siderophore_region = pvd_region,
            hmm_pvd_type           = pvd_type)


pangenome_prediction <- pangenome_siderotypes %>%
  mutate(bacteria_code = sub("^.*_(\\d+)$", "\\1", bacteria)) %>%
  # remove cip strains as we don't have pangenome annotation for them
  filter(!grepl("CIP", bacteria_code))
  
  

final_comparison <- hmm_pvd_prediction %>%
  right_join(pangenome_prediction,
            by = "bacteria_code") %>%
  mutate(agreement = case_when(is.na(hmm_pvd_type) ~ "Disagree",    # some hmms couldn't predict the branch at all
                               pvd_type == hmm_pvd_type ~ "Agree",
                               .default = "Disagree")) %>%
  group_by(agreement) %>%
  summarise(n = n(),
            per = n/477 * 100)

# Type the cip strains ----
## get previous siderotypes
in_vivo_pvd_types <- data.frame(bacteria_code = c("CIP6092",
                                        "CIP5940",
                                        "CIP5937",
                                        "CIP5938",
                                        "107",
                                        "109"),
                        pvd_type    = c(3,
                                        2,
                                        3,
                                        2,
                                        1,
                                        1))
regions_to_check <- c("CIP6092",
                      "CIP5940",
                      "CIP5937",
                      "CIP5938",
                      "107",
                      "109")

cip_prediciton <- hmm_pvd_prediction %>%
  filter(bacteria_code %in% regions_to_check) %>%
  left_join(in_vivo_pvd_types,
            by = "bacteria_code")




