# Dependencies -----------------------------------------------------------------
library(tidyverse)      # Data Wrangling
library(treedataverse)  # phylogenetic tree plotting
library(aplot)          # helpers for phylogenetic plotting
library(scales)         # Easy log10 axis
library(tidygraph)      # Network plot
library(ggraph)         # Network plotting, tidygraph companion



# Options ----------------------------------------------------------------------
## Inputs
HOG_file               <- "./data/orthofinder/N0.tsv"
bacteria_metadata_file <- "./raw_data/pgdb_strain_summary.tsv"
query_metadata_file    <- "./raw_data/20240711_smith_et_al_Genes.csv"
bakta_annotations      <- "./data/bakta/compiled_tsvs.csv"
iats_siderotypes_file  <- "./raw_data/iats_siderotypes.csv"


#iats_siderotypes_file  <- "./raw_data/iats_siderotypes.csv"

## Outputs
figure_dir       <- "./data/figures/"
annotations_file <- "./data/reannotated_genomes"

# Load data --------------------------------------------------------------------
## Load metadata from previously known genes
query_metadata <- read.csv(file = query_metadata_file)

## Load gff files
gffs <- read.csv(bakta_annotations) %>%
  rename(bacteria_code = Organism) %>%
  filter(Type == "cds") %>%
  rename(bakta_gene    = Gene,
         bakta_product = Product) %>%
  mutate(Pseudogene    = grepl("\\(pseudo\\) ", bakta_product),
         bakta_gene    = gsub("Unknown", "hypothetical protein", bakta_gene),
         bakta_product = gsub("\\(pseudo\\) ", "", bakta_product),
         bakta_product = gsub("conserved, ", "", bakta_product),
         bakta_product = gsub("conserved hypothetical protein", "hypothetical protein", bakta_product),
         bakta_product = gsub("Conserved ", "", bakta_product),
         bakta_product = gsub("hypothetical protein CHP02391 domain-containing protein", "CHP02391 domain-containing protein", bakta_product),
         bakta_product = gsub("Uncharacterized conserved protein, ", "", bakta_product),
         bakta_product = gsub("Helicase conserved C-terminal domain-containing protein", "Helicase conserved C-terminal domain", bakta_product),
         bakta_product = gsub("Uncharacterized conserved protein", "hypothetical protein", bakta_product),
         bakta_product = gsub("Uncharacterized conserved small protein", "hypothetical protein", bakta_product),
         bakta_product = gsub("Uncharacterized protein conserved in bacteria", "hypothetical protein", bakta_product),
         bakta_product = gsub("Uncharacterized protein ", "", bakta_product),
         bakta_product = gsub("Uncharacterized protein", "hypothetical protein", bakta_product),
         bakta_product = gsub("Uncharacterized ", "", bakta_product),
         bakta_product = gsub("Uncharacterized, putative ", "", bakta_product),
         bakta_annotation = ifelse(bakta_gene == "Unknown",
                                   bakta_product,
                                   bakta_gene)) 

## Load bacteria metadata
bacteria_metadata <- read.delim(bacteria_metadata_file) %>%
  select(!(Number.of.genes:GTF.format.annotations...file.name)) %>%
  select(!c(Species, Database.version)) %>%
  rename(bacteria = Full.organism.name)


## Load orthofinder annotations
HOGs <- read.delim(file = HOG_file,
                   header = TRUE,
                   sep = "\t") %>%
  pivot_longer(cols = -c(HOG, OG, Gene.Tree.Parent.Clade),
               names_to = "bacteria",
               values_to = "Loci") %>%
  select(-bacteria) %>%
  mutate(Locus_Tag = strsplit(Loci, ", ")) %>%
  unnest(Locus_Tag) %>%
  select(-Loci)

## Merge assignments into one table
## NB I have deviated here by using bacteria_metadata instead of codes
## This changes the "bacteria" field to be their full names
## This could be better or it could break something
all_loci <- gffs %>%
  left_join(HOGs,
            by = c("Locus_Tag" = "Locus_Tag")) %>%
  left_join(bacteria_metadata,
            by = c("bacteria_code" = "Strain.ID"))

## genes in an OG of 1 member are lost from N0.tsv, as you can't make a HOG from one gene
## I need to add them as HOGs, with a number 
## TODO some OGs are also NA, not sure I need them but if so add here
max_HOG <- max(sub("N0.HOG", "", all_loci$HOG),
               na.rm = TRUE)

new_HOG_numbers <- all_loci %>%
  select(Locus_Tag, HOG) %>%
  filter(is.na(HOG)) %>%
  rownames_to_column(var = "HOG_number") %>%
  rowwise() %>%
  mutate(manual_HOG = sprintf("N0.HOG%07d", as.numeric(max_HOG) + as.numeric(HOG_number))) %>%
  select(Locus_Tag, manual_HOG)

all_loci <- all_loci %>%
  left_join(new_HOG_numbers,
            by = c("Locus_Tag" = "Locus_Tag")) %>%
  mutate(HOG = case_when(is.na(manual_HOG) ~ HOG,
                         is.na(HOG)        ~ manual_HOG)) %>%
  select(-c(manual_HOG))

## Load info on IATS strains (previously siderotyped)
known_siderotypes <- read.csv(iats_siderotypes_file) %>%
  filter(!is.na(siderotype))

iats_codes <- known_siderotypes %>%
  filter(IATS) %>%
  pull(Strain.ID)

# Annotate HOGs ----------------------------------------------------------------
HOG_hypothetcial_percentages <- all_loci %>%
  mutate(bakta_gene = gsub("Unknown", "hypothetical protein", bakta_gene)) %>%
  group_by(HOG, bakta_gene) %>%
  summarise(annotation_freq = n(),
            .groups = "drop") %>%
  group_by(HOG) %>%
  reframe(HOG = HOG,
          bakta_gene = bakta_gene,
          total_proteins = sum(annotation_freq),
          annotation_percent = annotation_freq/total_proteins * 100) %>%
  filter(bakta_gene == "hypothetical protein") %>%
  transmute(HOG           = HOG,
            hypothetical_percent = annotation_percent)

HOG_annotations_without_hypothetical <- all_loci %>%
  mutate(bakta_gene = gsub("Unknown", "hypothetical protein", bakta_gene)) %>%
  filter(bakta_gene != "hypothetical protein") %>%
  group_by(HOG, bakta_gene) %>%
  summarise(annotation_freq = n(),
            .groups = "drop") %>%
  group_by(HOG) %>%
  reframe(bakta_gene = bakta_gene,
          n_annotations = n(),
          annotations   = paste0(bakta_gene, collapse = ", "),
          annotation_freq = annotation_freq,
          total_proteins = sum(annotation_freq),
          annotation_percent = annotation_freq/total_proteins * 100) %>%
  group_by(HOG) %>%
  filter(annotation_percent >= 90) %>%
  transmute(HOG = HOG,
            annotation_without_hypothetical = paste0(bakta_gene, "_", gsub("N0.HOG0*(?=0$)|N0.HOG0*", "", HOG, perl = TRUE)))

HOG_annotation <- all_loci %>%
  mutate(bakta_gene = gsub("Unknown", "hypothetical protein", bakta_gene)) %>%
  group_by(HOG, bakta_gene) %>%
  summarise(annotation_freq = n(),
            .groups = "drop") %>%
  group_by(HOG) %>%
  reframe(bakta_gene = bakta_gene,
          n_annotations = n(),
          annotations   = paste0(bakta_gene, collapse = ", "),
          annotation_freq = annotation_freq,
          total_proteins = sum(annotation_freq),
          annotation_percent = annotation_freq/total_proteins * 100) %>%
  group_by(HOG) %>%
  left_join(HOG_hypothetcial_percentages,
            by = c("HOG" = "HOG")) %>%
  mutate(percent_with_hypothetical = annotation_percent + hypothetical_percent,
         HOG_annotation = gsub(", hypothetical protein", "", annotations),
         HOG_annotation = gsub("hypothetical protein, ", "", HOG_annotation),
         HOG_annotation   = ifelse(test = hypothetical_percent == 100,
                                     yes  = paste0("HP_", gsub("N0.HOG0*", "", HOG)),
                                     no   = ifelse(test = str_count(HOG_annotation, ",") == 0,
                                                   yes  = paste0(HOG_annotation, "_", gsub("N0.HOG0*", "", HOG)),
                                                   no   = paste0("Unclear_", gsub("N0.HOG0*(?=0$)|N0.HOG0*", "", HOG, perl = TRUE))))) %>%
  left_join(HOG_annotations_without_hypothetical,
            by = c("HOG" = "HOG")) %>%
  ungroup() %>%
  mutate(HOG_annotation = ifelse(test = !is.na(annotation_without_hypothetical),
                                   yes  = annotation_without_hypothetical,
                                   no   = HOG_annotation)) %>%
  select(HOG, HOG_annotation) %>%
  distinct()

all_loci <- all_loci %>%
  left_join(HOG_annotation,
            by = c("HOG" = "HOG"))


# Section 1 summary stats ------------------------------------------------------
n_ORFs <- nrow(all_loci)
n_HOGs <- all_loci %>%
  pull(HOG) %>%
  unique() %>%
  length()

HOG_classes <- all_loci %>%
  select(HOG_annotation) %>%
  distinct() %>%
  mutate(HOG_group = case_when(grepl("^HP_\\d+$", HOG_annotation) ~ "Hypothetical",
                               grepl("Unclear_\\d+$", HOG_annotation) ~ "Unclear",
                               .default = "Clear")) %>%
  group_by(HOG_group) %>%
  summarise(n = n()) %>%
  mutate(percentage = n/n_HOGs * 100)

# Manually resolve annotations for some HOGs -----------------------------------
# Find unclear orthogroups which contain our genes of interest
OGs_to_examine <- all_loci %>%
  filter(bakta_annotation %in% query_metadata$locus_name_in_source) %>%
  filter(grepl("Unclear", HOG_annotation)) %>%
  select(OG, HOG, HOG_annotation, Gene.Tree.Parent.Clade) %>%
  distinct()

# Make colour scheme for the genes
unclear_genes <- all_loci %>%
  filter(HOG %in% OGs_to_examine$HOG) %>%
  select(HOG, bakta_annotation) %>%
  pull(bakta_annotation) %>%
  unique()

gene_pal <- c("#470024", # solG
              "#040F0F", # Hypothetical
              "#0FA3B1", # viuB
              "#EDDEA4", # sip
              "#FF9B42", # dltA
              "#BBBE64", # yojI
              "#7D8491", # pvdI
              "#FF6978", # pvdJ
              "#340068", # pvdD
              "#2BA84A", # pvdE
              "#357DED", # entF
              "#483D03") # pvdY 

gene_colors <- data.frame(gene  = c(unclear_genes, "pvdY"),
                          color = gene_pal)

## Unclear 59 ------------------------------------------------------------------
### Load raw tree
raw_tree_59 <- read.newick(file = "./data/orthofinder/Resolved_Gene_Trees/raw_treefiles/OG0000005_tree.txt")

### Add in metadata
treedata_59 <- as_tibble(raw_tree_59) %>%
  mutate(Locus_Tag = gsub(".*_(?=[^_]*_[^_]*$)", "", label, perl = TRUE)) %>%
  left_join(all_loci %>%
              select(Locus_Tag, bakta_gene, HOG),
            by = c("Locus_Tag" = "Locus_Tag")) %>%
  left_join(gene_colors,
            by = c("bakta_gene" = "gene"))

### Filter to "N0.HOG0000059"
HOG_56_node    <- 1549 # determined from viewing entire tree
tree_59        <- tree_subset(tree = as.treedata(treedata_59),
                              node = 1549,
                              levels_back = 0)

### View tree to work out what nodes I should subset from
exploratory_59 <- ggtree(tr     =  tree_59,
                         layout = "rect",
                         x      = length) +
  geom_tippoint(mapping = aes(color = color),
                show.legend = T) +
  geom_nodelab(aes(label = node),
               size = 2) +
  scale_color_identity(guide = "legend",
                       name  = "Annotation",
                       labels = c("Hypothetical\nprotein",
                                  "pvdD",
                                  "entF",
                                  "solG",
                                  "pvdI",
                                  "pvdJ",
                                  "dltA")) +
  coord_flip() +
  scale_x_reverse()

dltA_bottom_node <- 1868
dltA_middle_node <- 1855
dltA_top_node    <- 1845
entF_node        <- NA
pvdD_node        <- 1478
pvdI_node        <- 1282
pvdJ_node        <- 1674
solG_node        <- 1062

node_assigments_59 <- data.frame(node = c(dltA_bottom_node, dltA_middle_node, dltA_top_node,
                                          entF_node, pvdD_node,
                                          pvdI_node, pvdJ_node, solG_node),
                                 gene = c("dltA", "dltA", "hypothetical protein", "entF", "pvdD",
                                          "pvdI", "pvdJ", "solG")) %>%
  left_join(gene_colors,
            by = c("gene" = "gene"))

### Assign split out into new HOGs
for(gene in c("entF", "solG",
              "pvdD", "pvdI", "pvdJ")){
  
  node        <- paste0(gene,"_node")
  node_number <- get(node)
  
  tmp_genes <- as_tibble(tree_59) %>%
    filter(node %in% offspring(tree_59, node_number, type = "tips")) %>%
    transmute(gene                   = gene,
              Locus_Tag              = Locus_Tag,
              Gene.Tree.Parent.Clade = paste0("n", node_number))
  
  assign(x = paste0(gene, "_loci"),
         value = tmp_genes)
  
}

### dltA is a special case
dltA_bottom_loci <-  as_tibble(tree_59) %>%
  filter(node %in% offspring(tree_59, dltA_bottom_node, type = "tips")) %>%
  transmute(gene                   = "HP",
            Locus_Tag              = Locus_Tag,
            Gene.Tree.Parent.Clade = paste0("n", dltA_bottom_node))

dltA_middle_loci <-  as_tibble(tree_59) %>%
  filter(node %in% offspring(tree_59, dltA_middle_node, type = "tips")) %>%
  filter(!(Locus_Tag %in% dltA_bottom_loci$Locus_Tag)) %>%
  transmute(gene                   = "dltA",
            Locus_Tag              = Locus_Tag,
            Gene.Tree.Parent.Clade = paste0("n", dltA_middle_node))

dltA_top_loci <-  as_tibble(tree_59) %>%
  filter(node %in% offspring(tree_59, dltA_top_node, type = "tips")) %>%
  filter(!(Locus_Tag %in% dltA_bottom_loci$Locus_Tag)) %>%
  filter(!(Locus_Tag %in% dltA_middle_loci$Locus_Tag)) %>%
  transmute(gene                   = "dltA",
            Locus_Tag              = Locus_Tag,
            Gene.Tree.Parent.Clade = paste0("n", dltA_top_node))


hp_loci <- as_tibble(tree_59) %>%
  filter(!(Locus_Tag %in% c(dltA_bottom_loci$Locus_Tag, dltA_middle_loci$Locus_Tag,
                            dltA_top_loci$Locus_Tag, entF_loci$Locus_Tag,
                            solG_loci$Locus_Tag, pvdD_loci$Locus_Tag,
                            pvdI_loci$Locus_Tag, pvdJ_loci$Locus_Tag))) %>%
  na.omit() %>%
  transmute(gene                   = "HP",
            Locus_Tag              = Locus_Tag,
            Gene.Tree.Parent.Clade = NA)

Unclear_59_reassignments <- do.call(rbind, list(dltA_bottom_loci, dltA_middle_loci,
                                                dltA_top_loci, entF_loci, solG_loci,
                                                pvdD_loci, pvdI_loci, pvdJ_loci,
                                                hp_loci))

## Unclear_235 -----------------------------------------------------------------
### Load raw tree
raw_tree_235 <- read.newick(file = "./data/orthofinder/Resolved_Gene_Trees/raw_treefiles/OG0000111_tree.txt")

### Add in metadata
treedata_235 <- as_tibble(raw_tree_235) %>%
  mutate(Locus_Tag = gsub(".*_(?=[^_]*_[^_]*$)", "", label, perl = TRUE)) %>%
  left_join(all_loci %>%
              select(Locus_Tag, bakta_gene, HOG),
            by = c("Locus_Tag" = "Locus_Tag")) %>%
  left_join(gene_colors,
            by = c("bakta_gene" = "gene"))

### Plot to find which node to split at
exploratory_235 <- ggtree(tr      =  as.treedata(treedata_235 %>%
                                                   mutate(branch.length = branch.length + 0.5)),
                          layout = "rect",
                          x      = length) +
  geom_tippoint(mapping = aes(color = color),
                show.legend = T) +
  geom_nodelab(aes(label = node),
                 size = 2) +
  scale_color_identity(guide = "legend",
                       name  = "Annotation",
                       labels = c("Hypothetical\nprotein",
                                  "viuB",
                                  "sip")) +
  coord_flip() +
  scale_x_reverse()

exploratory_235

sip_node <- 841

### Assign split out into new HOGs
tree_235 <- as.treedata(treedata_235)

sip_genes <- as_tibble(tree_235) %>%
  filter(node %in% offspring(tree_235, sip_node, type = "tips")) %>%
  transmute(gene                   = "sip",
            Locus_Tag              = Locus_Tag,
            Gene.Tree.Parent.Clade = paste0("n", sip_node))

viuB_genes <- as_tibble(tree_235) %>%
  filter(!(node %in% offspring(tree_235, sip_node, type = "tips"))) %>%
  na.omit() %>%
  transmute(gene                   = "viuB",
            Locus_Tag              = Locus_Tag,
            Gene.Tree.Parent.Clade = NA)

Unclear_235_reassignments <- do.call(rbind, list(sip_genes,
                                                 viuB_genes))

## Unclear_1075 ----------------------------------------------------------------
### Get raw tree
raw_tree_1075 <- read.newick(file = "./data/orthofinder/Resolved_Gene_Trees/raw_treefiles/OG0000933_tree.txt")

### Add in metadata
treedata_1075 <- as_tibble(raw_tree_1075) %>%
  mutate(Locus_Tag = gsub(".*_(?=[^_]*_[^_]*$)", "", label, perl = TRUE)) %>%
  left_join(all_loci %>%
              select(Locus_Tag, bakta_gene, HOG),
            by = c("Locus_Tag" = "Locus_Tag")) %>%
  left_join(gene_colors,
            by = c("bakta_gene" = "gene"))

### Plot to find which node to split at
exploratory_1075 <- ggtree(tr      =  as.treedata(treedata_1075 %>%
                                                    mutate(branch.length = branch.length + 0.5)),
                           layout = "rect",
                           x      = length) +
  geom_tippoint(mapping = aes(color = color),
                show.legend = T) +
  geom_nodelab(aes(label = node),
                 size = 1.5,
                 vjust = 0,
                 hjust = 0) +
  scale_color_identity(guide = "legend",
                       name  = "Annotation",
                       labels = c("Hypothetical\nprotein",
                                  "pvdE",
                                  "yojI")) +
  coord_flip() +
  scale_x_reverse()

exploratory_1075

pvdE_node <- 746
yoji_node <- 651

## Assign split out into new HOGs
tree_1075 <- as.treedata(treedata_1075)

pvdE_genes <- as_tibble(tree_1075) %>%
  filter(node %in% offspring(tree_1075, pvdE_node, type = "tips")) %>%
  transmute(gene = "pvdE",
            Locus_Tag = Locus_Tag,
            Gene.Tree.Parent.Clade = paste0("n", pvdE_node))

yojI_genes <- as_tibble(tree_1075) %>%
  filter((node %in% offspring(tree_1075, yoji_node, type = "tips"))) %>%
  transmute(gene = "yojI",
            Locus_Tag              = Locus_Tag,
            Gene.Tree.Parent.Clade = paste0("n", yoji_node))

hp_1075_genes <-  as_tibble(tree_1075) %>%
  filter(!(node %in% offspring(tree_1075, yoji_node, type = "tips"))) %>%
  filter(!(node %in% offspring(tree_1075, pvdE_node, type = "tips"))) %>%
  na.omit() %>%
  transmute(gene = "HP",
            Locus_Tag              = Locus_Tag,
            Gene.Tree.Parent.Clade = NA)


Unclear_1075_reassignments <- do.call(rbind, list(pvdE_genes,
                                                  yojI_genes,
                                                  hp_1075_genes))
## pvdY  -----------------------------------------------------------------------
### Find OG pvdY is in
pvdY_1_OGs <- all_loci %>%
  filter(grepl("pvdY", HOG_annotation)) %>%
  select(OG, HOG, HOG_annotation) %>%
  distinct()

pvdY_OG <- "OG0004708" 

### Load raw tree
raw_tree_4708 <- read.newick(file = paste0("./data/orthofinder/Resolved_Gene_Trees/raw_treefiles/", pvdY_OG, "_tree.txt"))

### Add in metadata
treedata_4708 <- as_tibble(raw_tree_4708) %>%
  mutate(Locus_Tag = gsub(".*_(?=[^_]*_[^_]*$)", "", label, perl = TRUE)) %>%
  left_join(all_loci %>%
              select(Locus_Tag, bakta_gene, HOG),
            by = c("Locus_Tag" = "Locus_Tag")) %>%
  left_join(gene_colors,
            by = c("bakta_gene" = "gene"))


### View tree to work out what nodes I should subset from
exploratory_4708 <- ggtree(tr      =  as.treedata(treedata_4708),
                           layout = "rect",
                           x      = length) +
  geom_tippoint(mapping = aes(color = color),
                show.legend = T) +
  geom_nodelab(aes(label = node),
               size = 2) +
  scale_color_identity(guide = "legend",
                       name  = "Annotation",
                       labels = c("Hypothetical\nprotein",
                                  "pvdY")) +
  coord_flip() +
  scale_x_reverse()

pvdY_bottom_node <- 521
pvdY_top_node    <- 472 

### Assign split out into new HOGs
pvdY_bottom_genes <- as_tibble(treedata_4708) %>%
  filter(node %in% offspring(as.treedata(treedata_4708), pvdY_bottom_node, type = "tips")) %>%
  transmute(gene                   = "pvdY",
            Locus_Tag              = Locus_Tag,
            Gene.Tree.Parent.Clade = paste0("n", pvdY_bottom_node))

pvdY_top_genes <- as_tibble(treedata_4708) %>% # 941 rows
  filter(node %in% offspring(as.treedata(treedata_4708), pvdY_top_node, type = "tips")) %>%
  filter(!(Locus_Tag %in% pvdY_bottom_genes$Locus_Tag)) %>%
  transmute(gene                   = "HP",
            Locus_Tag              = Locus_Tag,
            Gene.Tree.Parent.Clade = paste0("n", pvdY_top_node))


pvdY_4708_reassignments <- do.call(rbind, list(pvdY_bottom_genes,
                                               pvdY_top_genes))

## Combine all manual reassigment data together
reassignment_data <- do.call(what = rbind,
                             args = list(Unclear_59_reassignments,
                                         Unclear_235_reassignments, 
                                         Unclear_1075_reassignments,
                                         pvdY_4708_reassignments)) %>%
  rename(new_annotation = gene,
         new_Gene.Tree.Parent.Clade = Gene.Tree.Parent.Clade) %>%
  # add the HOGs back in
  left_join(all_loci,
            by = c("Locus_Tag" = "Locus_Tag"))

# re-annotate HOGs off of unclear reassignments --------------------------------
unclear_reassignemtns <- reassignment_data %>%
  select(new_annotation, HOG, new_Gene.Tree.Parent.Clade) %>%
  distinct() %>%
  group_by(HOG) %>%
  mutate(new_HOG = paste0(HOG, "_", row_number())) %>%
  left_join(reassignment_data,
            by = c("HOG" = "HOG",
                   "new_annotation" = "new_annotation",
                   "new_Gene.Tree.Parent.Clade" = "new_Gene.Tree.Parent.Clade")) %>%
  ungroup() %>%
  select(new_annotation, Locus_Tag, new_HOG) %>%
  right_join(all_loci,
             by = c("Locus_Tag" = "Locus_Tag")) %>%
  mutate(HOG = case_when(is.na(new_HOG) ~ HOG,
                         .default = new_HOG),
         HOG_annotation_2 = case_when(is.na(new_annotation) ~ HOG_annotation,
                                      .default = paste0(new_annotation, "_", gsub("N0.HOG0*", "", HOG))))

## redo numbering
new_numbering <- unclear_reassignemtns %>%
  mutate(short_HOG_annotation = sub("_\\d+.*$", "", HOG_annotation_2)) %>%
  select(HOG, short_HOG_annotation) %>%
  distinct() %>%
  group_by(short_HOG_annotation) %>%
  mutate(HOG_annotation_2 = paste0(sub("_\\d+", "", short_HOG_annotation), "_", row_number())) %>%
  ungroup() %>%
  select(-short_HOG_annotation) %>%
  right_join(unclear_reassignemtns,
             by = c("HOG" = "HOG")) %>%
  transmute(Locus_Tag        = Locus_Tag,
            HOG_annotation_2 = HOG_annotation_2.x,
            new_HOG          = HOG) %>%
  distinct()

## merge into all_loci
all_loci <- all_loci %>%
  left_join(new_numbering,
            by = c("Locus_Tag" = "Locus_Tag")) %>%
  mutate(HOG = case_when(is.na(new_HOG) ~ HOG,
                         .default = new_HOG)) %>%
  left_join(reassignment_data %>%
              select(Locus_Tag, new_Gene.Tree.Parent.Clade),
            by = c("Locus_Tag" = "Locus_Tag")) %>%
  mutate(Gene.Tree.Parent.Clade = case_when(is.na(Gene.Tree.Parent.Clade) ~ Gene.Tree.Parent.Clade,
                                            .default = new_Gene.Tree.Parent.Clade),
         HOG_annotation = HOG_annotation_2) %>%
  select(-HOG_annotation_2)

## Make new list of potential genes off of the reassignments
potential_siderophore_HOGs <- all_loci %>%
  filter(gsub("_\\d+", "", HOG_annotation) %in% query_metadata$locus_name_in_source) %>%
  pull(HOG_annotation) %>%
  unique()

potential_siderophore_HOGs <- c(potential_siderophore_HOGs, as.character(gene_colors_2$gene))


## Save out annotations for later use
write.csv(x         = all_loci,
          file      = annotations_file,
          row.names = FALSE)

# Plot phylogenetic trees of reassignemnts -------------------------------------
tree_labels <- do.call(rbind, list(treedata_235,
                                   treedata_1075,
                                   treedata_4708,
                                   treedata_59)) %>%
  select(label, Locus_Tag) %>%
  filter(!grepl("^n.*$", label))

## Make new colour pallet 
## TODO  move this up to where the pallet is made at the start, so it is included in my pvdY tree
genes <- all_loci %>%
  filter(Locus_Tag %in% tree_labels$Locus_Tag) %>%
  pull(HOG_annotation) %>%
  unique() %>%
  sort()

gene_colors_2 <- data.frame(gene = genes,
                            color = c("#F28500", # dltA_1,
                                      "#FF7F50", # dltA_2
                                      "#b7410e" , # dltA_3
                                      "#E8E9EB", # HP_1
                                      "#D3D3D3", # HP_2,
                                      "#BEBEBE", # HP_3,
                                      "#A9A9A9", # HP_4
                                      "#928885", # HP_5
                                      "#340068", # pvdD_1
                                      "#2BA84A", # pvdE_1
                                      "#7D8491", # pvdI_1
                                      "#FF6978", # pvdJ_1
                                      "#483D03", # pvdY_1
                                      "#EDDEA4", # sip_1
                                      "#CD143C", # solG_1
                                      "#801818", # solG_2,
                                      "#0FA3B1", # viuB
                                      "#BBBE64")) %>%# yojI_1 
  mutate(gene = factor(gene, levels = genes),
         color= factor(color, levels = color))

## Make heatmapt of new assignments
new_assigments_data <- all_loci_2 %>%
  filter(Locus_Tag %in% tree_labels$Locus_Tag) %>%
  select(Locus_Tag, HOG_annotation) %>%
  left_join(tree_labels,
            by = "Locus_Tag") %>%
  select(label, HOG_annotation) %>%
  pivot_longer(cols = HOG_annotation,
               names_to = "name",
               values_to = "value") %>%
  left_join(gene_colors_2,
            by = c("value" = "gene"))

gene_annotation_heatmap <- ggplot(data = new_assigments_data,
                                  mapping = aes(x = name,
                                                y = label,
                                                fill = color)) +
  geom_tile() +
  scale_fill_identity(guide  = "legend",
                      name   = "Annotation",
                      labels = gene_colors_2$gene) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())

## lets do treedata_59 to start as it was broken
tree_plot <- ggtree(tr = as.treedata(treedata_59),
                    layout = "rect",
                    x      = length) +
  geom_tippoint(mapping = aes(color = color),
                show.legend = T) +
  scale_color_identity(guide = "legend",
                       name = "Bakta\nAnnotation",
                       labels = c("HP",
                                  "pvdD",
                                  "entF",
                                  "solG",
                                  "pvdI",
                                  "pvdJ",
                                  "dltA")) +
  theme_tree2()

full_tree_plot <- gene_annotation_heatmap %>%
  insert_left(tree_plot, width = 10)

for(tree in c("treedata_235",
              "treedata_1075",
              "treedata_4708",
              "treedata_59")){
  
  tmp <- get(tree)
  
  tree_plot <- ggtree(tr = as.treedata(tmp),
                      layout = "rect",
                      x      = length) +
    geom_tippoint(mapping = aes(color = color),
                  show.legend = F) +
    scale_color_identity() +
    theme_tree2()
  
  full_tree_plot <- gene_annotation_heatmap %>%
    insert_left(tree_plot, width = 10)
  
  assign(x = paste0(tree, "_plot"),
         value = full_tree_plot)
  
}

all_trees <- plot_list(treedata_59_plot, treedata_235_plot,
                       treedata_1075_plot, treedata_4708_plot,
                       design = "AABB\nAACC\nAADD",
                       widths = c(1, 0.5))
all_trees

ggsave(path     = figure_dir,
       filename = "figure_1_phylogenetic_trees.svg",
       plot     = all_trees,
       device   = "svg",
       width    = 15,
       height   = 15,
       dpi      = 300)

# Identify PVD region in complete genomes -------------------------------------
## Are genes in a consistent position?
window_size_bp  <- 50000 

sid_HOG_density <- ggplot(data    = all_loci %>%
                                    filter(HOG_annotation %in% potential_siderophore_HOGs) %>%
                                    filter(!(bacteria_code %in% iats_codes)) %>%
                                    mutate(start_kb = Start/1000,
                                           start_mb = start_kb/1000),
                                  mapping = aes(x    = start_kb)) +
  geom_histogram(binwidth    = 50,
                 show.legend = TRUE) +
  scale_y_continuous(name = "Count (n)",
                     expand = c(0, 0)) +
  scale_x_continuous(name = "\n Gene Start (kb)",
                     expand = c(0, 0)) +
  scale_fill_discrete(name = "Gene") +
  scale_color_discrete(name = "Gene") +
  coord_cartesian(xlim = c(0, 8000),
                  ylim = c(0, 2000)) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.6))

sid_HOG_density

## Identify region size using all potential genes 
region_sizes <- all_loci %>%
  filter(HOG_annotation %in% potential_siderophore_HOGs) %>%
  filter(!(bacteria_code %in% iats_codes)) %>%
  select(bacteria, Locus, Start, Stop, HOG_annotation) %>%
  group_by(bacteria) %>%
  filter(Locus == min(Locus) |
           Locus == max(Locus)) %>%
  summarise(region_size_bp = max(Start, Stop) - min(Start, Stop),
            region_size_kb = region_size_bp/1000)

all_gene_size_plot <- ggplot(data = region_sizes,
                             mapping = aes(x = region_size_kb)) +
  geom_histogram(binwidth = 50) +
  scale_x_continuous(name = "Region size (Kb)",
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Count (n)",
                     expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 200),
                  xlim = c(0, 4050)) +
  theme_bw() 

all_gene_size_plot

## Identify genes always in the region
sid_genes_in_region <- all_loci %>%
  filter(HOG_annotation %in% potential_siderophore_HOGs) %>%
  select(bacteria, Sequence, HOG_annotation, HOG, Locus, Start, Stop) %>%
  group_by(bacteria, Sequence) %>%
  arrange(Locus) %>%
  mutate(distance_to_next_sid = (lead(Locus, 1) - Locus)/5 - 1,
         distance_to_prev_sid = (Locus - lag(Locus, 1))/5 - 1) %>%
  rowwise() %>%
  filter(any(!is.na(distance_to_next_sid), !is.na(distance_to_prev_sid))) %>%
  mutate(distance_to_nearest_sid = min(c(distance_to_next_sid, distance_to_prev_sid), na.rm = TRUE))

intergene_distance_plot <- ggplot(data = sid_genes_in_region %>%
                                    mutate(distance_to_nearest_sid_plus_0.1 = distance_to_nearest_sid + 0.1), 
                                  mapping = aes(y     = HOG_annotation,
                                                x     = distance_to_nearest_sid_plus_0.1,
                                                color = HOG_annotation,
                                                fill  = HOG_annotation)) +
  geom_point(show.legend = FALSE,
             alpha       = 0.2) +
  scale_y_discrete(name = "",
                   limits = rev) +
  geom_vline(xintercept = 3,
             linetype   = "dashed") +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),
                     name   = "Distance to nearest siderophore-related ORF",
                     expand = c(0.01, 0)) +
  theme_bw() +
  theme(strip.text.y.left = element_text(angle = 0),
        panel.spacing = unit(0,'lines'),
        panel.border = element_blank())

intergene_distance_plot

confident_sid_orthogroups <- sid_genes_in_region %>%
  group_by(HOG_annotation) %>%
  summarise(max_ORF_distance = max(distance_to_nearest_sid)) %>%
  filter(max_ORF_distance < 3) %>%
  pull(HOG_annotation)

## Define siderophore region boundaries
region_boundaries <- all_loci %>%                         
  filter(HOG_annotation %in% confident_sid_orthogroups) %>%
  group_by(bacteria, Sequence) %>%
  mutate(distance_to_next_sid = (lead(Locus, 1) - Locus)/5,
         distance_to_prev_sid = (Locus - lag(Locus, 1))/5) %>%
  rowwise() %>%
  mutate(distance_to_furthest_sid = max(distance_to_next_sid, distance_to_prev_sid, na.rm = TRUE), # this returns Inf for isolated genes
         distance_to_nearest_sid  = min(distance_to_next_sid, distance_to_prev_sid, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(distance_to_nearest_sid < 3) %>%
  group_by(bacteria, Sequence) %>%
  summarise(region_locus_start = min(Locus),
            region_locus_end   = max(Locus))

## define full siderophore region, containing all genes
full_sid_regions <- all_loci %>%
  left_join(region_boundaries,
            by = c("bacteria" = "bacteria",
                   "Sequence" = "Sequence")) %>%
  group_by(bacteria, Sequence) %>%
  filter(Locus >= region_locus_start - 5 &
         Locus <= region_locus_end   + 5)

## Plot size of siderophore regions
high_confidence_region_sizes <- full_sid_regions %>%
  filter(!(bacteria_code %in% iats_codes)) %>%
  select(bacteria, Locus, Start, Stop, HOG_annotation) %>%
  group_by(bacteria) %>%
  filter(Locus == min(Locus) |
           Locus == max(Locus)) %>%
  summarise(region_size_bp = max(Start, Stop) - min(Start, Stop),
            region_size_kb = region_size_bp/1000)

high_conf_gene_size_plot <- ggplot(data = high_confidence_region_sizes,
                             mapping = aes(x = region_size_kb)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(name = "Region size (Kb)",
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Count (n)",
                     expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100),
                  xlim = c(0, 200)) +
  theme_bw() 

high_conf_gene_size_plot

left_peak_median <- high_confidence_region_sizes %>%
  filter(region_size_kb < 93) %>%
  pull(region_size_bp) %>%
  median()

right_peak_median <- high_confidence_region_sizes %>%
  filter(region_size_kb > 93) %>%
  pull(region_size_bp) %>%
  median()

## Save out plots
region_sizes_plot <- plot_list(sid_HOG_density, all_gene_size_plot, high_conf_gene_size_plot)

ggsave(path     = figure_dir,
       filename = "figure_2_region_sizes.png",
       plot     = region_sizes_plot,
       device   = "png",
       width    = 15,
       height   = 5,
       dpi      = 300)

ggsave(path     = figure_dir,
       filename = "./supplementary/s_figure_1_intergene_distance.png",
       plot     = intergene_distance_plot,
       device   = "png",
       width    = 15,
       height   = 10,
       dpi      = 300)

# Make pangenome ---------------------------------------------------------------
region_invert <- function(region){
  
  vector_region <- str_split_1(string = region,
                               pattern = ", ")
  
  vector_region_no_numbers <- gsub(pattern     = ":\\d+",
                                   replacement = "",
                                   x           = vector_region)
  
  flipped_region <- rev(vector_region_no_numbers)
  
  flipped_region_numbers <- data.frame(x = flipped_region) %>%
    mutate(gene_with_locus = paste0(x, ":", (row_number() - 1) * 5)) %>%
    pull(gene_with_locus) %>%
    paste(collapse = ", ")
  
  return(flipped_region_numbers)
} 

## Get all possible edges 
gene_pairs <- full_sid_regions %>%
  group_by(bacteria, Sequence) %>%
  mutate(position_in_region = Locus - min(Locus),
         annotation_with_position = paste0(HOG_annotation, ":", Locus - min(Locus))) %>%
  summarise(region_contents = paste0(annotation_with_position, collapse = ", "),
            .groups = "drop") %>%
  select(region_contents) %>%
  distinct() %>%
  rowid_to_column(var = "sequence_number") %>%
  rowwise() %>%
  # check region orientation WRT PAO1
  mutate(orientation = case_when(grepl("^pvdX", region_contents) | grepl("pvdQ_\\d:\\d+$", region_contents) ~ "-",
                                 grepl("^pvdQ", region_contents) | grepl("pvdX_\\d:\\d+$", region_contents) ~ "+",
                                 .default = NA)) %>%
  # Flip inverted regions 
  mutate(region_contents = case_when(orientation == "-" ~ region_invert(region_contents),
                                     .default = region_contents)) %>%
  # Split regions back out into vectors
  mutate(potential_gene = strsplit(region_contents, ", ")) %>%
  unnest(potential_gene) %>%
  select(-region_contents) %>%
  separate(col  = potential_gene,
           into = c("potential_gene", "relative_locus"),
           sep  = ":") %>%
  # Calculate region pairs
  group_by(sequence_number) %>%
  mutate(next_gene = lead(potential_gene, 1)) %>%
  # remove genes which do not have a proceeding gene
  na.omit() %>%
  # Add our edges together 
  rowwise() %>%
  transmute(sequence_number = sequence_number,
            edge        = paste(potential_gene, next_gene, sep = "~")) %>%
  rownames_to_column(var = "edge_number")

## Get number of sequences each edge occurs in
n_regions_with_edge <- gene_pairs %>%
  group_by(edge) %>%
  summarise(n_regions_edge_in = n())

## filter edges based on threshold
edges_to_use <- n_regions_with_edge %>%
  filter(n_regions_edge_in >= length(unique(gene_pairs$sequence_number)) * 0.1) 

## Calcualte co-occurence 
co_occurence <- gene_pairs %>%
  # filter to chosen edges
  filter(edge %in% unique(edges_to_use$edge)) %>%
  separate(edge,
           into = c("gene_1", "gene_2"),
           sep  = "~") %>%
  group_by(gene_1, gene_2) %>%
  # calcualte number of times an edge pair occurs
  summarise(co_occurence = n()) %>%
  # calcualte number of times a gene_1 starts an edge
  group_by(gene_1) %>%
  reframe(gene_1 = gene_1,
          gene_2 = gene_2,
          co_occurence = co_occurence,
          n_edges_gene_1_starts = sum(co_occurence)) %>%
  # Calcaulte co-occurence frequencey
  mutate(gene1_ass_gene2 = co_occurence/n_edges_gene_1_starts)

## Plot fancy network
gene_network_graph <- as_tbl_graph(co_occurence %>%
                                     select(gene_1, gene_2, gene1_ass_gene2),
                                   directed = TRUE) 

## Gene labels
pvd1_branch <- c("pvdE_1", "fpvA_1", "pvdD_1", "pvdJ_1", "pvdI_1")
pvd2_branch <- c("HP_3", "HP_589", "aes_1", "solG_1", "dltA_3")
pvd3_branch <- c("yojI_1", "fpvA_2", "dltA_1", "HP_1", "dltA_2", "sip_1")

pvd_core    <- c("pvdQ_1", "pvdA_1", "fpvI_1", "fpvR_1", "pvdR_1", "pvdT_1", "opmQ_1",
                 "pvdP_1", "pvdM_1", "pvdN_1", "pvdO_1", "pvdF_1", "fpvG_1", "fpvH_1", 
                 "fpvJ_1", "fpvK_1", "fpvC_1", "fpvD_1", "fpvE_1","fpvF_1", "pA2411_1",
                 "pA2412_1", "pvdH_1", "pvdL_1", "pvdG_1", "pvdS_1", "pvdX_1")
  
gene_labels <- data.frame(gene = c(pvd1_branch,
                                   pvd2_branch,
                                   pvd3_branch),
                          class = c("PVD Export", "PVD-Fe Import", "NRPS", "NRPS", "NRPS",
                                    "PVD Export", "PVD-Fe Import", "Hydrolase", "NRPS", "NRPS",
                                    "PVD Export", "PVD-Fe Import", "NRPS", "NRPS", "Reductase"))

fancy_network_graph <- gene_network_graph %>%
  activate(edges) %>%
  mutate(similar   = gene1_ass_gene2,
         distance  = 1 - gene1_ass_gene2) %>%
  activate(nodes) %>%
  left_join(gene_labels,
            by = c("name" = "gene")) %>%
  ggraph(layout = "stress",
         weight = similar) +
  geom_edge_link(arrow = arrow(length = unit(1, 'mm')), 
                 end_cap = circle(1, 'mm')) +
  geom_edge_loop(arrow = arrow(length = unit(1, 'mm')), 
                 end_cap = circle(1, 'mm')) +
  geom_node_point(show.legend = TRUE,
                  size = 11,
                  mapping = aes(color = class)) +
  scale_color_discrete(name = "Function",
                       labels = c("Hydrolysis",
                                  "Synthesis",
                                  "Import",
                                  "Export",
                                  "Reduction")) +
  geom_node_text(aes(label = name,
                     color = class),
                 size = 10,
                 repel = TRUE,
                 point.padding = unit(1, "cm"),
                 show.legend = FALSE) +
#  scale_color_manual(name   = "Gene\nClassification",
#                     labels = c("Not Previously Associated",
#                                "Hypothetical Protein",
#                                "Previously Assoicated",
 #                               "High Confidence Orthogroup",
#                                "Unclear Function"),
#                     values = c("black",
#                                "blue",
#                                "green",
#                                "orange",
#                                "red")) +
  coord_flip() +  
  scale_x_reverse() +
  scale_y_reverse() +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text  = element_blank(),
        axis.line  = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.1, 0.85),
        legend.title = element_text(size = 50),
        legend.text = element_text(size = 40)) +
  guides(color = guide_legend(ncol = 1))

fancy_network_graph

ggsave(path     = figure_dir,
       filename = "figure_3_pangenome_portrait.png",
       plot     = fancy_network_graph,
       device   = "png",
       width    = 15,
       height   = 20,
       dpi      = 300)


# Look at hypothetical genes in branches ---------------------------------------
HP_HOGs <- c("HP_1", "HP_3", "HP_589")

HP_anotations <- full_sid_regions %>%
  filter(HOG_annotation %in% HP_HOGs) %>%
  group_by(HOG_annotation, DbXrefs) %>%
  summarise(n = n())

HP_counts <- full_sid_regions %>%
  filter(HOG_annotation %in% HP_HOGs) %>%
  group_by(HOG_annotation) %>%
  summarise(n = n())

HP_dbxrefs <- full_sid_regions %>%
  filter(HOG_annotation %in% HP_HOGs) %>%
  select(HOG_annotation, Locus_Tag, DbXrefs, bakta_product) %>%
  mutate(RefSeq = case_when(grepl("RefSeq:", DbXrefs) ~ sub(".*RefSeq:([^, ]+).*", "\\1", DbXrefs),
                            .default = NA),
         SO     = case_when(grepl("SO:", DbXrefs) ~ sub(".*SO:([^, ]+).*", "\\1", DbXrefs),
                            .default = NA),
         UniParc = case_when(grepl("UniParc", DbXrefs) ~ sub(".*UniParc:([^, ]+).*", "\\1", DbXrefs),
                             .default = NA),
         UniRef_100 = case_when(grepl("UniRef:UniRef100_", DbXrefs) ~ sub(".*UniRef:UniRef100_([^, ]+).*", "\\1", DbXrefs),
                                .default = NA),
         UniRef_50  = case_when(grepl("UniRef:UniRef50_", DbXrefs) ~ sub(".*UniRef:UniRef50_([^, ]+).*", "\\1", DbXrefs),
                                .default = NA),
         UniRef_90 = case_when(grepl("UniRef:UniRef90_", DbXrefs) ~ sub(".*UniRef:UniRef90_([^, ]+).*", "\\1", DbXrefs),
                               .default = NA)) %>%
  select(-DbXrefs) %>%
  pivot_longer(cols = -c(HOG_annotation, Locus_Tag),
               values_to = "annotation",
               names_to  = "database") %>%
  group_by(HOG_annotation, database, annotation) %>%
  summarise(n_annotations = n()) %>%
  left_join(HP_counts,
            by = "HOG_annotation") %>%
  rename(n_occurences = n)

top_HP_annotations <- HP_dbxrefs %>%
  group_by(HOG_annotation, database) %>%
  filter(n_annotations == max(n_annotations))


# Predict PVD type -------------------------------------------------------------
# Get expected lengths of all the genes
expected_lengths <- full_sid_regions %>% 
  group_by(HOG_annotation) %>%
  summarise(median_size = median(Stop - Start))

# Repair regions
## TODO why does this even need a genes of interest?
## ~ remake so all genes are repaired? 
## ~ then it could happen before pangenome is drawn
region_repair <- function(split_all_loci, region_boundaries, region_pad, Locus_col,
                          annotation_col, expected_lengths, genes_of_interest){
  
  used_loci         <- NULL
  repaired_region   <- NULL
  
  region_of_interest <- split_all_loci %>%
    left_join(region_boundaries,
              by = c("bacteria" = "bacteria",
                     "Sequence" = "Sequence")) %>%
    filter(Locus >= region_locus_start - region_pad &
             Locus <= region_locus_end   + region_pad)
  
  # First repair the region considering adjacent genes etc
  for(Locus in region_of_interest[, Locus_col]){
    
    # Add start/end of string to ensure exact regex matches
    Locus <- paste0("^", Locus, "$")
    
    #writeLines(as.character(Locus))
    #writeLines(as.character(repaired_region))
    
    # Get info on the gene I want
    gene_info     <- region_of_interest[grep(Locus, region_of_interest[, "Locus"]), ]
    gene          <- gene_info[, annotation_col]
    gene_length   <- gene_info[, "Stop"] - gene_info[, "Start"]
    target_length <- unlist(expected_lengths[expected_lengths$"HOG_annotation" == gene, "median_size"])
    length_padding<- 0.2
    gene_strand   <- gene_info[, "Strand"]
    
    # Check if already assesed
    if(Locus %in% used_loci){
      
      next
      #writeLines("Locus already used")
      
    }
    
    #writeLines("Locus not used")
    # Is the gene one I care about?
    #if(!(gene %in% genes_of_interest)){
      
    #  next
      # writeLines("Gene not of interest")
      
      
    #}
    
    #writeLines("Gene of interest")
    
    # Is the length correct?
    if(gene_length >= target_length * (1 - length_padding) &
       gene_length <= target_length * (1 + length_padding)){
      
      repaired_region <- c(repaired_region, gene)
      used_loci       <- c(used_loci, Locus)
      #writeLines("Length Correct")
      
      # if incorrect consider adjacent genes  
      # TODO include special case where we are at the boundaries of a region and therefore
      # we have no next locus. 
      # Or we should consider all_loci? 
    } else {
      
      #writeLines("Length incorrect")
      
      if(gene_strand == "+"){
        
        next_gene_info <- split_all_loci[grep(Locus, split_all_loci[, "Locus"]) + 1, ]
        #writeLines("Strand +")
        
      } else {
        
        next_gene_info <- split_all_loci[grep(Locus, split_all_loci[, "Locus"]) - 1, ]
        #writeLines("Strand -")
        
      }
      
      # Get info for next locus
      next_gene_locus  <- paste0("^", next_gene_info[, "Locus"], "$")
      next_gene        <- next_gene_info[, annotation_col]
      next_gene_strand <- next_gene_info[, "Strand"]
      
      
      # Check locus hasn't already been used
      if(next_gene_locus %in% used_loci){
        
        next
        
      }
      
      
      
      # Get new length if they were merged
      new_start <- min(gene_info[, "Start"], next_gene_info[, "Start"])
      new_stop  <- max(gene_info[, "Stop"], next_gene_info[, "Stop"])
      
      combined_length <- new_stop - new_start
      
      # Check all atributes match
      strand_match <- next_gene_strand == gene_strand
      name_match   <- next_gene == gene | grepl("^HP_\\d+$", next_gene)
      length_match <- combined_length >= target_length * (1 - length_padding) &
        combined_length <= target_length * (1 + length_padding)
      
      if(strand_match & name_match & length_match){
        
        repaired_region <- c(repaired_region, gene)
        used_loci       <- c(used_loci, Locus, next_gene_locus)
        #writeLines("Repair with next gene")
        
      }
      
      
      
    }
    
    
    
  }
  
  return(paste0(repaired_region, collapse = ", "))
  
  
}

repaired_regions <- lapply(X                 = split(x = all_loci,
                                                     f = ~ bacteria),
                           FUN               = region_repair,
                           Locus_col         = "Locus",
                           region_boundaries = region_boundaries,
                           region_pad        = 5,
                           annotation_col    = "HOG_annotation",
                           expected_lengths  = expected_lengths) %>%
  do.call(rbind, .) %>%
  data.frame() %>%
  rownames_to_column(var = "bacteria") %>%
  rename(siderophore_region = ".")

# Type strains
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
    
    pvd_type <- "unknown"
    quality  <- NA
    
  }
  
  # Make final output
  out <- c(core_result,
           branch_results,
           paste0("pvd_type=", pvd_type),
           paste0("prediction_quality=", quality))
  
  return(paste(out, collapse = ", "))
  
}

siderotype_strain_predictions <- repaired_regions %>%
  rowwise() %>%
  mutate(region_stats = pvd_typing(region    = siderophore_region,
                                   core      = pvd_core,
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
              values_from = "value") 

## get summary stats
pvd_strain_stats <- siderotype_strain_predictions %>%
  group_by(pvd_type, prediction_quality) %>%
  summarise(n   = n()) %>%
  mutate(prediction_quality = factor(prediction_quality, levels = c("high", "medium", "low")),
         per = n/477 * 100) %>%
  arrange(pvd_type, prediction_quality,
          by.groups = TRUE)

## stats combining all prediction types together
overall_strain_stats <- pvd_strain_stats %>%
  group_by(pvd_type) %>%
  summarise(n = sum(n),
            per = n/477 * 100)

# Type unique regions 
siderotype_region_predictions <- repaired_regions %>%
  select(siderophore_region) %>%
  distinct() %>%
  rowwise() %>%
  mutate(region_stats = pvd_typing(region    = siderophore_region,
                                   core      = pvd_core,
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
              values_from = "value") 

## get summary stats
pvd_region_stats <- siderotype_region_predictions %>%
  group_by(pvd_type, prediction_quality) %>%
  summarise(n   = n(),
            per = n/276 * 100) %>%
  mutate(prediction_quality = factor(prediction_quality, levels = c("high", "medium", "low"))) %>%
  arrange(pvd_type, prediction_quality,
          by.groups = TRUE)

## stats combining all prediction types together
overall_region_stats <- pvd_region_stats %>%
  group_by(pvd_type) %>%
  summarise(n = sum(n),
            per = n/sum(pvd_region_stats$n) * 100)

# Check predictions of IATS strains
known_predictions <- siderotype_strain_predictions %>%
  left_join(all_loci %>%
              select(bacteria, bacteria_code) %>%
              distinct(),
            by = c("bacteria" = "bacteria")) %>%
  filter(bacteria_code %in% known_siderotypes$Strain.ID) %>%
  left_join(known_siderotypes,
            by = c("bacteria_code" = "Strain.ID")) %>%
  transmute(bacteria            = bacteria,
            known_siderotype    = siderotype,
            predicted_siderotpe = as.numeric(sub("pvd", "", pvd_type)), 
            prediction_quality  = prediction_quality) %>%
  arrange(known_siderotype)

known_accuracy <- known_predictions %>%
  mutate(correct = known_siderotype == predicted_siderotpe) %>%
  group_by(correct) %>%
  summarise(n = n())





