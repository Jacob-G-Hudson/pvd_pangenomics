# Script to take outputs of NRPStransformer and and see how NPRS domains link to
# PVD type 
# TODO work out why the number of unique pvd sequences does not equal the number of typed strains
# Dependences ----
library(tidyverse)    # data wrangling
library(gggenes)      # gene plots
library(RColorBrewer) # nicer pallets


# setup ----
annotations_file       <- "./data/summaries/reannotated_genomes.csv"
NRPStransformer_dir    <- "./data/NRPS_transformer/results"
PVD_predictions_file   <- "./data/summaries/siderotype_strain_predictions.csv"
  
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
                   bacteria = "^.*results/(.*).csv$") %>%
  do.call(rbind, .)


nrps_data <- raw_nrps %>%
  transmute(bacteria  = bacteria,
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
              select(Locus_Tag, HOG_annotation_4),
            by = c("Locus_Tag" = "Locus_Tag")) %>%
  left_join(PVD_predictions,
            by = c("bacteria" = "bacteria"))

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

pvd1_sequence <- paste(c(pvdI_1_aa, pvdJ_1_aa, pvdD_1_aa),
                       collapse = "*")

## PVDII
solG_1_aa <- "ser-hxorn-orn"
dltA_3_aa <- "gly-thr-ser-hxorn"

pvd2_sequence <- paste(c(solG_1_aa, dltA_3_aa),
                       collapse = "*")

## PVDIII
dltA_1_aa <- "ser-24dab-hxorn"
HP_1_aa   <- "gln"
dltaA_2_aa<- "gln-hxorn-gly"

pvd3_sequence <- paste(c(dltA_1_aa, HP_1_aa, dltaA_2_aa),
                       collapse = "*")

# In smith et al we have fOHOrn (formyl-N5-hydroxyornithine)
# In NRPS translator we have hxorn (N(5)-hydroxy-N(5)-acyl ornithine)
# I assume that these are equivalent

# Look at gene sizes to better compare with smith et al ----
gene_sizes <- all_loci %>%
  filter(HOG_annotation_4 %in% c(pvd1_path, pvd2_path, pvd3_path)) %>%
  mutate(length = abs(Stop - Start))

ggplot(data = gene_sizes,
       mapping = aes(x = HOG_annotation_4,
                     y = length)) +
  geom_boxplot()

# Get amino acid sequences for each enzyme ----
PVD_predictions %>% # 87 in here are NA!
  group_by(siderotype_prediction) %>%
  summarise(n = n())

## NB is this keeping the domain order the same? I imagine NRPS translocator keeps it preserved?
nrps_sequences <- nrps_data %>%                 #446 bacteria
  # remove strains for which we couldn't type 
  filter(!is.na(pvd_type)) %>%     #389! 
  # filter to just branch specific genes
  filter(HOG_annotation_4 != "pvdL_1") %>%
  # get best score
  filter(domain_rank == 1) %>%
  group_by(bacteria, HOG_annotation_4, Locus_Tag) %>%
  summarise(aa_sequence = paste(amino_acid,
                                collapse = "-")) %>%
  ungroup()

unique_nrps_sequences <- nrps_sequences %>%
  group_by(HOG_annotation_4, aa_sequence) %>%
  summarise(n = n(),
            Locus_Tags = paste(Locus_Tag, collapse = ", "))

# Get amino acid sequences for each pvd ----
## Get the sequences of the pvd made by each bacteria
pvd_sequences <- nrps_sequences %>%
  left_join(all_loci %>%
              select(Locus_Tag, Locus, Strand),
            by = "Locus_Tag") %>%
  group_by(bacteria) %>%
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
            by = "bacteria") 
  
## Get just distinct sequences for analysis
### From enzymes
distinct_enzyme_sequences <- pvd_sequences %>%
  group_by(enzyme_sequence, pvd_type, prediction_quality) %>%
  summarise(n = n(),
            bacteria = paste(bacteria,
                             collapse = ", ")) %>%
  ungroup() %>%
  mutate(cannonical = case_when(pvd_type == "pvd1" & enzyme_sequence == pvd1_sequence ~ TRUE,
                                pvd_type == "pvd2" & enzyme_sequence == pvd2_sequence ~ TRUE,
                                pvd_type == "pvd3" & enzyme_sequence == pvd3_sequence ~ TRUE,
                                .default = FALSE),
         prediction_quality = factor(prediction_quality, levels = c("high", "medium", "low"))) %>%
  arrange(pvd_type, prediction_quality)

### From pvd pathway
distinct_pvd_sequences <- pvd_sequences %>%
  group_by(pvd_sequence, pvd_type, prediction_quality) %>%
  summarise(n = n(),
            bacteria = paste(bacteria,
                             collapse = ", ")) %>%
  ungroup() %>%
  mutate(cannonical = case_when(pvd_type == "pvd1" & pvd_sequence == gsub("\\*", "-", pvd1_sequence) ~ TRUE,
                                pvd_type == "pvd2" & pvd_sequence == gsub("\\*", "-", pvd2_sequence) ~ TRUE,
                                pvd_type == "pvd3" & pvd_sequence == gsub("\\*", "-", pvd3_sequence) ~ TRUE,
                                .default = FALSE),
         prediction_quality = factor(prediction_quality, levels = c("high", "medium", "low"))) %>%
  arrange(pvd_sequence, prediction_quality)

cannonical_sequences <- distinct_pvd_sequences %>%
  filter(cannonical) %>%
  group_by(pvd_sequence, pvd_type)  %>%
  summarise(n = sum(n)) 

non_cannonical_sequences <- distinct_pvd_sequences %>%
  filter(!cannonical) %>%
  group_by(pvd_sequence, pvd_type)  %>%
  summarise(n = sum(n)) 

non_cannonical_sequences <- distinct_pvd_sequences %>%
  filter(!cannonical) %>%
  group_by(prediction_quality) %>%
  summarise(n = sum(n))

write.csv(x         = distinct_sequences,
          file      = "./data/NRPS_transformer/run_2/aa_sequences.csv",
          row.names = FALSE)
# Lets work out whats going on with Pa_i_M37351_3639 ----
# TODO make an example plot for every possible region to put alongside the table of NRPS sequences
## its region is pvdQ_1, pvdA_1, fpvI_1, fpvR_1, pvdR_1, pvdT_1, opmQ_1, pvdP_1, pvdM_1, pvdN_1, pvdO_1, pvdF_1, pvdE_1, fpvA_1, pvdD_1, pvdJ_1, pvdI_1, pvdI_1, fpvG_1, fpvH_1, fpvJ_1, fpvK_1, fpvC_1, fpvD_1, fpvE_1, fpvF_1, pA2411_1, pA2412_1, pvdH_1, yliI_1, HP_584, treA_1, lysR_49, yhaK_2, pncA_7, oprD_6, sphA_3, HP_2056, cas1f_1, cas3f_1, csy1_1, csy2_1, csy3_1, cas6f_1, HP_582, pvdL_1, pvdG_1, pvdS_1, pvdY_1, pvdX_1
## relevant bit is pvdE_1, fpvA_1, pvdD_1, pvdJ_1, pvdI_1, pvdI_1
## ~ should be NA but its giving 1 for some reason
if(FALSE){
sequences_investivation <- data.frame(bacteria = c("Pa_i_M37351_3639",
                                                   "Pa_PA1R_496",
                                                   "Pa_MRSN12280_6456",
                                                   "Pa_X78812_6522"),
                                      reason   = c("duplicate pvdI_1",
                                                   "duplicate dltA_3",
                                                   "duplicated dltA_3",
                                                   "duplicated pvdI_1"))

## lets plot the regions of these above strains to show becky
pvd1_to_plot <- distinct_sequences %>%
  filter(siderotype_prediction == 1) %>%
  filter(!cannonical) %>%
  pull(bacteria)
  
  
regions_of_interest <- all_loci %>%
  filter(bacteria %in% c(pvd1_to_plot, "Pa_DVT729_15615"))  %>%
  left_join(region_boundaries,
            by  = c("bacteria" = "bacteria",
                    "Sequence" = "Sequence")) %>%
  filter(Locus >= region_locus_start - 5 &     # Expand region to include genes with 1 ORFs either end
           Locus <= region_locus_end + 5) %>%
  filter(HOG_annotation_4 %in% c(pvd1_path, pvd2_path, pvd3_path)) %>%
  group_by(bacteria, Sequence) %>%
  # todo this isn't quite right for flipped regions. 
  mutate(aligned_start = Start - min(Start, Stop),
         aligned_stop  = Stop  - min(Start, Stop)) %>%
  select(bacteria, bacteria_code, Sequence, Strand, HOG_annotation_4, Locus_Tag,
         Start, Stop, aligned_start, aligned_stop, siderotype)

## Plot 
split_regions <- ggplot(data = regions_of_interest ,
                                   mapping = aes(xmin  = aligned_start,
                                                 xmax  = aligned_stop,
                                                 y     = bacteria,
                                                 label = HOG_annotation_4,
                                                 color = HOG_annotation_4)) +
  geom_gene_arrow(show.legend = FALSE) +
  geom_gene_label(grow = TRUE,
                  min.size = 0) +
  scale_fill_discrete(name = "Region") +
  scale_y_discrete(name = NULL) +
  theme_genes() + 
  facet_grid(cols = vars(HOG_annotation_4),
             scales = "free") +
  theme_grey() +
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

split_regions
}

# Choose bacteria for my example
# TODO why is Pa_F092021_1635 missing from the final plot?
example_bacteria <- c("Pa_X78812_6522",
                      "Pa_F092021_16351",
                      "Pa_PAO1_107",
                      "Pa_DVT729_15615",
                      "Pa_PA1R_496",
                      "Pa_MRSN12280_6456",
                      "Pa_VRFPA04_230",
                      "Pa_BA7823_9659",
                      "Pa_12_4_4_59_3618",
                      "Pa_152962_15963")

example_bacteria <- distinct_pvd_sequences %>%
  filter(!cannonical) %>%
  filter(prediction_quality != "low") %>%
  pull(bacteria)

## some others just for troubleshooting
example_bacteria <- c("Pa_PAO1_107",
                      "Pa_UCBPP_PA14_109")

# Get NRPS enzymes
NRPS_enzymes <- nrps_sequences %>%
  arrange(HOG_annotation_4) %>%
  pull(HOG_annotation_4) %>%
  unique()

# Set colours
NRPS_colors <-  data.frame(annotation = c(NRPS_enzymes, "HP"),
                           color      = c(brewer.pal(length(NRPS_enzymes), "Set3"),
                                          "#676795"))

# Get the boundaries of the branch
example_boundaries <- all_loci %>% # Pa_F092021_1635 is here
  filter(bacteria %in% example_bacteria) %>%
  filter(HOG_annotation_4 %in% c(pvd1_path, pvd2_path, pvd3_path)) %>%
  group_by(bacteria) %>%
  summarise(branch_start    = min(Locus),
            branch_end      = max(Locus),
            branch_start_bp = min(Start))
  
# Filter to boundaires, plus 2 extra loci 
example_regions <- all_loci %>%
  select(bacteria, Locus, HOG_annotation_4, Start, Stop, Strand) %>%
  filter(bacteria %in% example_bacteria) %>%
  left_join(example_boundaries,
              by = "bacteria") %>%
  rowwise() %>%
  filter(Locus >= branch_start - 5  &
         Locus <= branch_end   + 5) %>%
  mutate(bacteria = factor(bacteria, levels = example_bacteria))

# Align the regions so they all start at 0
## Get start co-ordinate
example_region_starts <- example_regions %>%
  group_by(bacteria) %>%
  summarise(start_bp = min(Start))

## adjust co-ordinates, and take into account the strand of the gene
example_aligned_regions <- example_regions %>%
  left_join(example_region_starts,
            by = "bacteria") %>%
  mutate(align_start  = Start - start_bp,
         aligned_stop = Stop  - start_bp) %>%
  rowwise() %>%
  mutate(forward    = Strand == "+",
         Gene       = case_when(HOG_annotation_4 %in% NRPS_enzymes ~ HOG_annotation_4,
                                grepl("^HP_\\d+$", HOG_annotation_4) ~ sub("_\\d+", "", HOG_annotation_4),
                                .default = NA),
         Label      = case_when(HOG_annotation_4 %in% c(pvd1_path, pvd2_path, pvd3_path) ~ HOG_annotation_4,
                                .default = NA)) %>%
  left_join(nrps_sequences %>%
              mutate(Locus = as.numeric(sub("^\\d+_(\\d+)$", "\\1", Locus_Tag))),
            by = c("bacteria" = "bacteria",
                   "Locus"    = "Locus",
                   "HOG_annotation_4" = "HOG_annotation_4")) %>%
  left_join(NRPS_colors,
            by = c("Gene" = "annotation")) %>%
  mutate(Gene  = factor(Gene, levels  = NRPS_colors$annotation),
         color = factor(color, levels = NRPS_colors$color))

# Plot regions
arrowhead_height <- 0.5

example_plot <- ggplot(data = example_aligned_regions,
                       mapping = aes(xmin    = align_start,
                                     xmax    = aligned_stop,
                                     forward = forward,
                                     y       = bacteria,
                                     label   = aa_sequence,
                                     fill    = color)) +
  geom_gene_arrow(show.legend = TRUE,
                  arrow_body_height = unit(arrowhead_height * 0.75, "cm"),
                  arrowhead_height  = unit(arrowhead_height, "cm")) +
  geom_gene_label(align    = "centre",
                  min.size = 0,
                  height   = unit(arrowhead_height * 0.75, "cm")) + # this should be 9pt font
  scale_fill_identity(guide  = "legend",
                      name   = "Gene",
                      labels = NRPS_colors$annotation) +
  scale_y_discrete(name = NULL) +
  theme_genes() + 
#  facet_grid(rows = vars(bacteria),
#             scales = "free") +
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
ggsave(filename = paste0("./data/NRPS_transformer/run_2/", arrowhead_height, "_erronous_regions_plot.png"),
       plot = example_plot,
       device = "png",
       width  = 21 - 2.54 - 2.54,
       height = 10,
       units  = "cm")

# Examine Pa_TL1285_16256 
# ~ it has a glu being added but I can't see it in the plot
test_bacteria <- "Pa_TL1285_16256"

test_region <- nrps_sequences %>%
  filter(bacteria %in% test_bacteria) %>%
  left_join(all_loci %>%
              select(Locus_Tag, Locus, Strand),
            by = "Locus_Tag") %>%
  arrange(Locus,
          by.groups = TRUE)
























