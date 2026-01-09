# script to take final hmms, and siderotype strains using them
# dependencies ---
library(tidyverse)  # data wrangling
library(rhmmer)     # utilities for hmmer files, from https://github.com/arendsee/rhmmer
library(pROC)       # roc curve plotting

# setup ----
annotations_file            <- "./data/reannotated_genomes.csv"
hmmer_dir                   <- "./data/hmm_building/"
true_positives_file         <- "./data/hmm_building/true_positives.txt"
thresholds_file             <- "./data/jackhmmer/hmm_building/hmm_thresholds.csv"
pangenome_siderotypes_file  <- "./data/siderotype_strain_predictions.csv"
in_vitro_siderotypes_file   <- "./raw_data/iats_siderotypes.csv"

figure_dir            <- "./data/figures/"

# Functions ----
pvd_typing <- function(region, core, branch){
  
  # Split region into a vector
  region <- str_split_1(string = region, pattern = ", ")
  
  # get region length
  region_length <- paste0("region_length=", length(region))
  
  # Check core presence
  # find genes in region from catagory
  ## NB this will pull them out even if they are not adjacent, and if they're in different blocks!
  present_core <- region[region %in% core]
  
  # count how many genes are in the region, and compare to expected
  n_core <- paste0(length(present_core), "/", length(core))
  
  # find if any genes are missing
  missing_core <- !all(core %in% present_core)
  
  ## Which genes are missing?
  if(missing_core){
    
    missing_core <- paste0(core[!(core %in% present_core)], collapse = "~")
    
  }
  
  ## are any genes duplicated?
  multiple_core <- any(duplicated(present_core))
  
  if(multiple_core){
    
    multiple_core <- paste0(present_core[duplicated(present_core)], collapse = "~")
    
  }
  
  core_result <-  c(paste0("core_length=", n_core),
                    paste0("core_missing=", missing_core),
                    paste0("core_multple=", multiple_core))
  
  # Check optional branch presence 
  branch_results <- lapply(names(branch), function(branch_name){
    
    # Get genes of interest for the category
    catagory_genes <- branch[[branch_name]]
    
    # find genes in region from catagory
    ## NB this will pull them out even if they are not adjacent, and if they're in different blocks!
    present_genes <- region[region %in% catagory_genes]
    
    # count how many genes are in the region, and compare to expected
    n_present <- paste0(length(present_genes), "/", length(catagory_genes))
    
    # find if any genes are missing
    missing_genes <- !all(catagory_genes %in% present_genes)
    
    ## Which genes are missing?
    if(missing_genes){
      
      missing_genes <- paste0(catagory_genes[!(catagory_genes %in% present_genes)], collapse = "~")
      
    }
    
    
    ## are any genes duplicated?
    multiple_genes <- any(duplicated(present_genes))
    
    if(multiple_genes){
      
      multiple_genes <- paste0(present_genes[duplicated(present_genes)], collapse = "~")
      
    }
    
    # Make output
    out <- c(paste0(branch_name, "_length=", n_present),
             paste0(branch_name, "_missing=", missing_genes),
             paste0(branch_name, "_multple=", multiple_genes))
    
    return(paste(out, collapse = ", "))
  })
  
  # Work out pvd_type
  ## This fails for low qualty predictions 
  branch_data <- do.call(rbind, branch_results) %>%
    data.frame() %>%
    rename(data = ".") %>%
    filter(grepl("_missing=FALSE", data)) %>%
    mutate(branch = sub("^([^_]+)_.*", "\\1", data),
           missing = sub(".*missing=([^,]+).*", "\\1", data))
  
  # filter to branches with components
  branch_data <- do.call(rbind, branch_results) %>%
    data.frame() %>%
    rename(data = ".") %>%
    filter(!grepl("^.*_length=0/.*$", data)) %>%
    mutate(branch = sub("^([^_]+)_.*", "\\1", data),
           missing = sub(".*missing=([^,]+).*", "\\1", data))
  
  possible_branches <- branch_data %>%
    pull(branch) %>%
    unique()
  
  branch_missing <-  branch_data %>%
    pull(missing) %>%
    unique()
  
  # if bits are missing then the missing_ is the name of the genes
  missing_core   <- !(missing_core == "FALSE")
  branch_missing <- !(branch_missing == "FALSE") 
  
  if(length(possible_branches) == 1){
    
    pvd_type <- possible_branches
    
    # work out prediction quality
    ## High quality has core + and branch +
    if(!missing_core & !branch_missing){
      
      quality <- "high"
      
      ## Medium quality has core - and branch + 
    } else if(missing_core & !branch_missing){
      
      
      quality <- "medium"
      
      ## Low quality has branch - and core + or - 
    } else {
      
      quality <- "low"
      
    }  
    
    
    ## Unknown type if we have no brnach components or multiple components  
  } else {
    
    pvd_type <- "Unknown"
    quality  <- NA
    
  }
  
  # Make final output
  out <- c(core_result,
           branch_results,
           paste0("pvd_type=", pvd_type),
           paste0("prediction_quality=", quality))
  
  return(paste(out, collapse = ", "))
  
}



# load files ---- 
all_loci              <- read.csv(annotations_file)
pangenome_siderotypes <- read.csv(pangenome_siderotypes_file) %>%
  mutate(bacteria_code = as.character(bacteria_code))
all_true_positives <- read.csv(true_positives_file) %>%
  mutate(prior_classification = 1)

in_vitro_pvd_types <- read.csv(in_vitro_siderotypes_file,
                               colClasses = "character") %>%
  filter(!is.na(siderotype))

## load data for 4-fold training hmms
training_output_files <- list.files(path       = hmmer_dir,
                                    pattern    = "*._results.out",
                                    full.names = TRUE,
                                    recursive  = TRUE)

training_data <- lapply(X   = as.list(training_output_files),
                        FUN = read_tblout) %>%
  do.call(rbind, .) %>% # here we have 276701
  transmute(Locus_Tag = domain_name,
            query_name= query_name,
            evalue    = sequence_evalue,
            bitscore  = sequence_score)

## Load data for the full hmms
final_hmm_files <- list.files(path       = hmmer_dir,
                              pattern    = "^results.out$",
                              full.names = TRUE,
                              recursive  = TRUE)

final_hmm_results <- lapply(X   = as.list(final_hmm_files),
                            FUN = read_tblout) %>%
  do.call(rbind, .) %>% # here we have 276701
  transmute(Locus_Tag = domain_name,
            query_name= query_name,
            evalue    = sequence_evalue,
            bitscore  = sequence_score)


# Make ROC curves --------------------------------------------------------------
training_files <- data.frame(output_file = training_output_files) %>%
  mutate(gene         = sub("^./data/hmm_building/(.*)/.*$", "\\1", output_file),
         hmm_version  = sub("^.*_(\\d+)_r.*$", "\\1", output_file),
         training_loci= paste0("./data/hmm_building/", gene, "/train_loci_", hmm_version, ".txt"))

hmm_load <- function(files, gene){
  
  # Filter to relevant files
  gene_files <- files %>%
    filter(gene == UQ(gene))
  
  # Load and process files
  hmm_results <- lapply(X           = as.list(gene_files$output_file),
                        FUN         = read_tblout) %>%
    do.call(rbind, .)  %>%
    rename(Locus_Tag = domain_name) %>%
    # get highes bitscore for each locus
    group_by(Locus_Tag) %>%
    summarise(bitscore    = max(sequence_score))  %>% # some true negatives are repeated so lets get the worst possible score
    mutate(hmm_annotation = gene)
  return(hmm_results)
  
}

## Load all hmm results
test_genes <- unique(training_files$gene)
train_hmm_results <- lapply(X     = as.list(test_genes),
                            FUN   = hmm_load,
                            files = training_files) %>%
  do.call(rbind, .) %>%
  # add in true positives
  left_join(all_true_positives,
            by = c("hmm_annotation" = "HOG_annotation",
                   "Locus_Tag"      = "Locus_Tag")) %>%
  replace_na(list(prior_classification = 0))

## Fit ROC curves
roc_fits <- train_hmm_results %>%
  mutate(hmm_annotation = factor(hmm_annotation, levels = test_genes)) %>%
  group_split(hmm_annotation) %>%
  map(~{roc(.x$prior_classification, .x$bitscore)})

names(roc_fits) <- test_genes

## Plot ROC curves
roc_plot <- ggroc(roc_fits, aes = c("color"),
                  legacy.axes = FALSE) +
  geom_line(color = "blue",
            linewidth = 1) +
  scale_x_reverse(name   = "\nSpecificity",
                  expand = c(0.01, 0.01)) +
  scale_y_continuous(name = "Sensitivity\n",
                     expand = c(0.01, 0.01)) +
  facet_wrap(~name) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.title = element_text(size = 9),
        axis.text  = element_text(size = 7),
        panel.spacing.x = unit(0.5, "cm"),
        panel.spacing.y = unit(0, "cm"))

roc_plot

ggsave(filename = "supplementary/supplementary_figure_3_roc_curves.png",
       path     = figure_dir,
       width    = 15,
       height   = 15,
       units    = "cm",
       device   = "png",
       dpi      = 300)

# Find bitscore thresholds ------------------------------------------------
hmm_thresholds <- map_df(roc_fits,
                         ~ tibble(threshold   = .x$threshold,
                                  sensitivity = .x$sensitivities,
                                  specificity = .x$specificities),
                         .id = "hmm")  %>%
  group_by(hmm) %>%
  filter(is.finite(threshold)) %>%
  filter(sensitivity == 1) %>%
  filter(specificity == max(specificity)) 

## Get per annotation accuracy for each hmm
#pre_hmm_accuracy <- all_loci %>%
#  select(Locus_Tag)
#  left_join(final_hmm_results,
#            by = "Locus_Tag") %>%
#  left_join(hmm_thresholds %>%
#              select(hmm, threshold),
#            by = c("query_name" = "hmm")) %>%
#  filter(bitscore > threshold) %>%
#  left_join()

# Save out thresholds for supplementary table 1
write.csv(x = hmm_thresholds,
          file = paste0(figure_dir, "supplementary/s_table_1_hmm_thresholds.csv"),
          row.names = FALSE)


# Siderotype strains with hmms -------------------------------------------------
# Get final gene annotations 
## Filter to one annotation per locus
hmm_annotation <- final_hmm_results %>%
  group_by(Locus_Tag) %>%
  # for loci with multiple annotations, take the best match
  filter(evalue   == min(evalue))   %>%
  filter(bitscore == max(bitscore)) %>% 
  # filter based on thresholds %>%
  left_join(hmm_thresholds %>%
              select(-c(sensitivity, specificity)),
            by = c("query_name" = "hmm")) %>% 
  filter(bitscore >= threshold) %>%
  rename(hmm_annotation = query_name)

# get correct rate 
per_loci_correct_rate <- all_loci %>%
  select(Locus_Tag, HOG_annotation) %>%
  left_join(hmm_annotation,
            by = "Locus_Tag") %>%
  mutate(class = case_when(HOG_annotation == hmm_annotation |
                           (!(HOG_annotation %in% unique(hmm_annotation[[, "hmm_annotation"]])) & is.na(hmm_annotation)) ~ "correct",
                           .default = "incorrect")) %>%
  group_by(class) %>%
  summarise(n = n())

# identify pvd region 
hmm_pvd_region <- hmm_annotation %>%
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
  select(-c(locus_lead_dif, locus_lag_dif)) %>%
  ungroup()

# Compare annotation accuracy with pangenome annotations
in_region_correct_rate <- hmm_pvd_region %>%
  # remove cip strains as we don't have pangenome annotation for them
  filter(!grepl("CIP", bacteria_code)) %>%  
  # to do work out a summary statistic for this stage such as senstivity/specificity/accuracy
  mutate(agreement = case_when(HOG_annotation == hmm_annotation ~ "Agree",
                               .default = "Disagree")) %>%
  group_by(agreement, hmm_annotation) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(hmm_annotation) %>%
  reframe(agreement = agreement,
          n = n,
          per = n/sum(n) * 100) %>%
  arrange(hmm_annotation) %>%
  filter(agreement == "Agree") %>%
  transmute(hmm_annotation = hmm_annotation,
            in_region_accuracy = per)

in_region_correct_rate <- all_loci %>%
  select(Locus_Tag, HOG_annotation) %>%
  left_join(hmm_pvd_region,
            by = "Locus_Tag") %>%
  mutate(class = case_when(HOG_annotation == query_name |
                             (!(HOG_annotation %in% hmm_annotation$query_name) & is.na(query_name)) ~ "correct",
                           .default = "incorrect")) %>%
  group_by(class) %>%
  summarise(n = n())

# hmm-based pvd_predictions ----
## Use same pvd_typing function as for pangenome
pvd1_branch <- c("pvdE_1", "fpvA_1", "pvdD_1", "pvdJ_1", "pvdI_1")
pvd2_branch <- c("HP_3", "HP_589", "aes_1", "solG_1", "dltA_3")
pvd3_branch <- c("yojI_1", "fpvA_2", "dltA_1", "HP_1", "dltA_2", "sip_1")


hmm_pvd_prediction <- hmm_pvd_region %>%
  group_by(bacteria_code) %>% 
  summarise(pvd_region = paste(hmm_annotation, collapse = ", ")) %>%
  rowwise() %>%
  mutate(region_stats = pvd_typing(region    = pvd_region,
                                   core      = NULL,
                                   branch    = list(pvd1 = pvd1_branch,
                                                    pvd2 = pvd2_branch,
                                                    pvd3 = pvd3_branch))) %>%
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

hmm_pvd_prediction_stats <- hmm_pvd_prediction %>%
  group_by(hmm_pvd_type) %>%
  summarise(n = n())

final_comparison <- hmm_pvd_prediction %>%
  right_join(pangenome_siderotypes %>%
               mutate(bacteria_code = as.character(bacteria_code)),
            by = "bacteria_code") %>%
  mutate(agreement = case_when(pvd_type == hmm_pvd_type ~ "Agree",
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

cip_prediciton <- in_vitro_pvd_types %>%
  mutate(long_read_id = case_when(CIP.number == "" ~ NA,
                                  .default = paste0("CIP", sub("\\.", "", CIP.number)))) %>%
  left_join(hmm_pvd_prediction %>%
              select(bacteria_code, hmm_pvd_type),
            by = c("long_read_id" = "bacteria_code")) %>%
  rename(long_read_hmm_type = hmm_pvd_type) %>%
  left_join(hmm_pvd_prediction %>%
              select(bacteria_code, hmm_pvd_type), 
            by = c("Strain.ID" = "bacteria_code")) %>%
  rename(short_read_hmm_type = hmm_pvd_type) %>%
  left_join(pangenome_siderotypes %>%
              select(bacteria_code, pvd_type),
            by = c("Strain.ID" = "bacteria_code"))  %>%
  rename(pangenome_type = pvd_type) %>%
 # filter(!is.na(hmm_pvd_type)) %>%
  arrange(siderotype) %>%
  transmute(Strain.ID                  = Strain.ID,
            in_vitro_pvd               = siderotype,
            pangenome_prediction       = sub("pvd", "", pangenome_type),
            hmm_prediction             = sub("pvd", "", short_read_hmm_type),
            resequenced_hmm_prediction = sub("pvd", "", long_read_hmm_type))




