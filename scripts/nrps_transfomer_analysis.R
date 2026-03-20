# Script to take outputs of NRPStransformer and and see how NPRS domains link to
# PVD type 
# Dependences ----
library(tidyverse)    # data wrangling
library(gggenes)      # gene plots
library(RColorBrewer) # nicer pallets


# setup ----
annotations_file       <- "./data/reannotated_genomes.csv"
NRPStransformer_dir    <- "./data/NRPS_transformer/results"
PVD_predictions_file   <- "./data/siderotype_strain_predictions.csv"
  
# Functions ----
csv_read <- function(file, skip = 0, ...){
  
  # Extract information from the file name based on provided names in "..."
  extracted_info   <- list(...)
  extracted_values <- lapply(extracted_info, function(pattern){
    
    if(!is.null(pattern)){
      
      gsub(pattern, replacement = "\\1", x = file)
      
    }else{
      
      NULL
      
    }
  })
  
  # Read the CSV file
  temp_file <- read.csv(file = file,
                        skip = skip)
  
  # Add extracted values as new columns
  for(i in seq_along(extracted_info)){
    
    if(!is.null(extracted_values[[i]])){
      
      col_name              <- names(extracted_info)[i]
      temp_file[, col_name] <- extracted_values[[i]]
      
    }
  }
  
  return(temp_file)
}

# Load data ----
all_loci <- read.csv(annotations_file)

PVD_predictions <- read.csv(PVD_predictions_file)
  
raw_nrps <- lapply(X   = list.files(path      = NRPStransformer_dir,
                                    full.name = TRUE),
                   FUN = csv_read,
                   bacteria_code = "^.*results/(.*).csv$") %>%
  do.call(rbind, .)

nrps_data <- raw_nrps %>%
  transmute(bacteria_code  = as.numeric(bacteria_code),
            Locus_Tag = ID,
            Domain    = Domain,
            Top_1     = Top.1.score.,
            Top_2     = Top.2.score.,
            Top_3     = Top.3.score.) %>%
  pivot_longer(cols      = Top_1:Top_3,
               values_to = "annotation",
               names_to  = "domain_rank") %>%
  separate_wider_delim(cols = annotation,
                       delim = "(",
                       names = c("amino_acid", "score")) %>%
  mutate(score = as.numeric(sub(")", "", score)),
         domain_rank = as.numeric(sub("Top_", "", domain_rank))) %>%
  left_join(all_loci %>%
              select(Locus_Tag, HOG_annotation),
            by = c("Locus_Tag" = "Locus_Tag")) %>%
  left_join(PVD_predictions,
            by = c("bacteria_code" = "bacteria_code"))

pvd1_branch <- c("pvdE_1", "fpvA_1", "pvdD_1", "pvdJ_1", "pvdI_1")
pvd2_branch <- c("HP_3", "HP_589", "aes_1", "solG_1", "dltA_3")
pvd3_branch <- c("yojI_1", "fpvA_2", "dltA_1", "HP_1", "dltA_2", "sip_1")

# expected sequences from smith et al
## common bit from pvdL
pvdL_1_aa <- "glu-tyr-24dab"

## PVDI
pvdI_1_aa <- "ser-arg-ser-hxorn"
pvdJ_1_aa <- "lys-hxorn"
pvdD_1_aa <- "thr-thr"

pvd1_sequence <- paste(c(pvdL_1_aa, pvdI_1_aa, pvdJ_1_aa, pvdD_1_aa),
                       collapse = "*")

## PVDII
solG_1_aa <- "ser-hxorn-orn"
dltA_3_aa <- "gly-thr-ser-hxorn"

pvd2_sequence <- paste(c(pvdL_1_aa, solG_1_aa, dltA_3_aa),
                       collapse = "*")

## PVDIII
dltA_1_aa <- "ser-24dab-hxorn"
HP_1_aa   <- "gln"
dltaA_2_aa<- "gln-hxorn-gly"

pvd3_sequence <- paste(c(pvdL_1_aa, dltA_1_aa, HP_1_aa, dltaA_2_aa),
                       collapse = "*")

# Find any missing bacteria, and what they're missing 
all_bacteria <- unique(all_loci$bacteria_code)
nrps_missing_bacteria <- all_loci %>%
  filter(!(bacteria_code %in% nrps_data$bacteria_code)) %>%
  select(bacteria_code) %>%
  distinct() %>%
  left_join(PVD_predictions,
            by = "bacteria_code")

# Get amino acid sequences for each enzyme ----
## NB is this keeping the domain order the same? I imagine NRPS translocator keeps it preserved?
nrps_sequences <- nrps_data %>% 
  # remove low confidence strains, as they are missing core genes by definition
  filter(prediction_quality %in% c("medium", "high")) %>% 
  # filter to just branch specific NRPS genes
  # get best score
  filter(domain_rank == 1) %>%
  group_by(bacteria_code, HOG_annotation, Locus_Tag) %>%
  summarise(aa_sequence = paste(amino_acid,
                                collapse = "-")) %>%
  ungroup()

unique_nrps_sequences <- nrps_sequences %>%
  group_by(HOG_annotation, aa_sequence) %>%
  summarise(n = n(),
            Locus_Tags = paste(Locus_Tag, collapse = ", "))

# Get amino acid sequences for each pvd ----
## Get the sequences of the pvd made by each bacteria
pvd_sequences <- nrps_sequences %>%
  left_join(all_loci %>%
              select(Locus_Tag, Locus, Strand),
            by = "Locus_Tag") %>%
  group_by(bacteria_code) %>%
  arrange(Locus,
          by.groups = TRUE) %>%
  summarise(Strand = unique(Strand),
            # Enzyme sequence includes info on what enzyme is adding the amino acid
            enzyme_sequence = case_when(Strand == "+" ~ paste(aa_sequence, collapse = "*"),
                                        Strand == "-" ~ paste(rev(aa_sequence), collapse = "*")),
            # PVD sequences is only the sequence of the final PVD moleucle
            pvd_sequence    = case_when(Strand == "+" ~ paste(aa_sequence, collapse = "-"),
                                        Strand == "-" ~ paste(rev(aa_sequence), collapse = "-"))) %>%
  left_join(PVD_predictions,
            by = "bacteria_code") 
  
## Get just distinct sequences for analysis
### From enzymes
distinct_enzyme_sequences <- pvd_sequences %>%
  group_by(enzyme_sequence, pvd_type, prediction_quality) %>%
  summarise(n = n(),
            bacteria_code = paste(bacteria_code,
                             collapse = ", ")) %>%
  ungroup() %>%
  mutate(cannonical = case_when(pvd_type == "pvd1" & enzyme_sequence == pvd1_sequence ~ TRUE,
                                pvd_type == "pvd2" & enzyme_sequence == pvd2_sequence ~ TRUE,
                                pvd_type == "pvd3" & enzyme_sequence == pvd3_sequence ~ TRUE,
                                .default = FALSE),
         prediction_quality = factor(prediction_quality, levels = c("high", "medium"))) %>%
  arrange(pvd_type, prediction_quality)

### From pvd pathway
distinct_pvd_sequences <- pvd_sequences %>%
  group_by(pvd_sequence, pvd_type, prediction_quality) %>%
  summarise(n = n(),
            bacteria_code = paste(bacteria_code,
                             collapse = ", ")) %>%
  ungroup() %>%
  mutate(cannonical = case_when(pvd_type == "pvd1" & pvd_sequence == gsub("\\*", "-", pvd1_sequence) ~ TRUE,
                                pvd_type == "pvd2" & pvd_sequence == gsub("\\*", "-", pvd2_sequence) ~ TRUE,
                                pvd_type == "pvd3" & pvd_sequence == gsub("\\*", "-", pvd3_sequence) ~ TRUE,
                                .default = FALSE),
         prediction_quality = factor(prediction_quality, levels = c("high", "medium", "low"))) %>%
  arrange(pvd_type, prediction_quality)

cannonical_v_not <- distinct_pvd_sequences %>%
  group_by(cannonical) %>%
  summarise(n = sum(n),
            per = n/408 * 100)

cannonical_sequences <- distinct_pvd_sequences %>%
  filter(cannonical) %>%
  group_by(pvd_sequence, pvd_type)  %>%
  summarise(n = sum(n)) 

non_cannonical_sequences <- distinct_pvd_sequences %>%
  filter(!cannonical) %>%
  group_by(pvd_sequence, pvd_type, prediction_quality)  %>%
  summarise(n = sum(n),
            bacteria_codes = paste(bacteria_code, collapse = ", ")) %>%
  arrange(pvd_type, prediction_quality)

example_sequence_table <- distinct_pvd_sequences %>%
  group_by(pvd_sequence, pvd_type, cannonical) %>%
  summarise(n = sum(n)) %>%
  mutate(cannonical = case_when(cannonical ~ "Y",
                                .default = "N"),
         cannonical = factor(cannonical, levels = c("Y", "N"))) %>%
  arrange(pvd_type, cannonical)

write.csv(x         = example_sequence_table,
          file      = "./data/figures/figure_3_PVD_aa_table.csv",
          row.names = FALSE)
# Plot an example region of each mis-typed genome-------------------------------
# Choose bacteria for my example
example_bacteria <- c(107,    #PVDI cannonical
                      15445,
                      172,     #PVDII cannonical
                      5865,
                      496,
                      5532,
                      16389,
                      489,    #PVDIII cannonical 
                      15674,
                      11967,
                      11152)

# Get NRPS enzymes
NRPS_enzymes <- nrps_sequences %>%
  arrange(HOG_annotation) %>%
  pull(HOG_annotation) %>%
  unique()

# Get the boundaries of the branch
# NB we include NRPS_enzymes as some genes have NRPS domains but are not branch specific
example_boundaries <- all_loci %>% # Pa_F092021_1635 is here
  filter(bacteria_code %in% example_bacteria) %>%
  filter(HOG_annotation %in% unique(c(pvd1_branch, pvd2_branch, pvd3_branch, NRPS_enzymes))) %>%
  group_by(bacteria_code) %>%
  summarise(branch_start    = min(Locus),
            branch_end      = max(Locus),
            branch_start_bp = min(Start))

# Filter to boundaires, plus 2 extra loci 
example_regions <- all_loci %>%
  select(bacteria_code, Locus, HOG_annotation, Start, Stop, Strand) %>%
  filter(bacteria_code %in% example_bacteria) %>%
  left_join(example_boundaries,
              by = "bacteria_code") %>%
  rowwise() %>%
  filter(Locus >= branch_start - 5  &
         Locus <= branch_end   + 5) %>%
  mutate(bacteria = factor(bacteria_code, levels = example_bacteria))

# Set colours
genes_of_interest <- example_regions %>%
  filter(HOG_annotation %in% NRPS_enzymes) %>%
  pull(HOG_annotation) %>%
  unique() %>%
  sort()

gene_colors <- data.frame(annotation = genes_of_interest,
                          color      = c("#F28500", # dltA_1,
                                         "#FF7F50", # dltA_2
                                         "#b7410e", # dltA_3
                                         "#928885", # HP_1
                                         "#d896ff", # pvdD_1
                                         "#7D8491", # pvdI_1
                                         "#FF6978", # pvdJ_1
                                         "#458b74", # pvdL_1
                                         "#CD143C")) # solG_1


# Align the regions so they all start at 0
## Get start co-ordinate
example_region_starts <- example_regions %>%
  group_by(bacteria_code) %>%
  summarise(start_bp = min(Start))

## adjust co-ordinates, and take into account the strand of the gene
example_aligned_regions <- example_regions %>%
  left_join(example_region_starts,
            by = "bacteria_code") %>%
  # align genes to common start point
  mutate(align_start  = Start - start_bp,
         aligned_stop = Stop  - start_bp) %>%
  rowwise() %>%
  mutate(forward    = Strand == "+",
         Gene       = case_when(HOG_annotation %in% NRPS_enzymes ~ HOG_annotation,
                                grepl("^HP_\\d+$", HOG_annotation) ~ sub("_\\d+", "", HOG_annotation),
                                .default = NA),
         Label      = case_when(HOG_annotation %in% c(pvd1_branch, pvd2_branch, pvd3_branch) ~ HOG_annotation,
                                .default = NA)) %>%
  left_join(nrps_sequences %>%
              mutate(Locus = as.numeric(sub("^\\d+_(\\d+)$", "\\1", Locus_Tag))),
            by = c("bacteria_code" = "bacteria_code",
                   "Locus"    = "Locus",
                   "HOG_annotation" = "HOG_annotation")) %>%
  left_join(gene_colors,
            by = c("Gene" = "annotation")) %>%
  mutate(Gene  = factor(Gene, levels  = gene_colors$annotation),
         color = factor(color, levels = gene_colors$color),
         bacteria_code = factor(bacteria_code, levels = rev(example_bacteria)))

# Plot regions
arrowhead_height <- 0.5

example_plot <- ggplot(data = example_aligned_regions,
                       mapping = aes(xmin    = align_start,
                                     xmax    = aligned_stop,
                                     forward = forward,
                                     y       = bacteria_code,
                                     label   = aa_sequence,
                                     fill    = color)) +
  geom_gene_arrow(show.legend = TRUE,
                  arrow_body_height = unit(arrowhead_height * 0.75, "cm"),
                  arrowhead_height  = unit(arrowhead_height, "cm")) +
  geom_gene_label(align    = "centre",
                  min.size = 0,
                  height   = unit(arrowhead_height * 0.75, "cm")) + # this should be 9pt font
  scale_fill_identity(guide  = "legend",
                      name   = "HOG",
                      labels = gene_colors$annotation) +
  scale_y_discrete(name = NULL) +
  theme_genes() + 
  theme_grey(base_size = 9) +
  theme(panel.background   = element_blank(), 
        panel.grid.major.y = element_line(colour = "grey", 
                                          linewidth = 1), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour    = "grey20", 
                                   linewidth = 0.5),
        axis.ticks.x = element_line(colour = "grey20",
                                    linewidth = 0.5),
        strip.text.y = element_text(angle = 0),
        strip.background = element_blank())

example_plot

# Save out
ggsave(filename = "./data/figures/figure_3_erronous_regions_plot.svg",
       plot = example_plot,
       device = "svg",
       width  = 21 - 2.54 - 2.54,
       height = 10,
       units  = "cm")





















