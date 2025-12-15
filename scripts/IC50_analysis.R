# IC50 analysis 
# Dependencies -----------------------------------------------------------------
library(tidyverse)
library(ggpubr)           # statistics for graphs
library(drc)              # DRC curve fit functions
library(rstatix)          # Statistical helper functions
library(patchwork)        # Multiplot functions

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

# Setup ------------------------------------------------------------------------
## ggplot
my_theme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


## Exclude some repeats
repeats_to_exclude <- c(3, # Bacterial contamination across plate as seen in image
                        4, # dose response curves are crap. PAO1 had no inhibition without iron
                        8) # 59.37 control well showed growth, as did its supernatant

# Load Data --------------------------------------------------------------------
well_key_r1 <- read.csv(file = "./raw_data/ic50/r1_plate_plan.csv") %>%
  mutate(Repeat = "1")
well_key_r2 <- read.csv(file = "./raw_data/ic50/r2_plate_plan.csv") %>%
  mutate(Repeat = "2")
well_key_r3 <- read.csv(file = "./raw_data/ic50/r1_plate_plan.csv") %>%
  mutate(Repeat = "3")
well_key_r4 <- read.csv(file = "./raw_data/ic50/r1_plate_plan.csv") %>%
  mutate(Repeat = "4")
well_key_r7 <- read.csv(file = "./raw_data/ic50/r2_plate_plan.csv") %>%
  mutate(Repeat = "7")
well_key_r8 <- read.csv(file = "./raw_data/ic50/r2_plate_plan.csv") %>%
  mutate(Repeat = "8")

well_key <- do.call(rbind, list(well_key_r1, well_key_r2, well_key_r3, well_key_r4,
                                well_key_r7, well_key_r8))

all_files <- list.files(path        = "./raw_data/ic50/",
                        pattern     = ".*_cip_ic50_.*.csv$",
                        full.names  = TRUE,
                        ignore.case = TRUE) 

raw_OD_data <- lapply(X = as.list(list.files(path        = "./raw_data/ic50/",
                                             pattern     = ".*_cip_ic50_.*.csv$",
                                             full.names  = TRUE,
                                             ignore.case = TRUE)),
                      FUN            = csv_read,
                      skip           = 5,
                      file_timestamp = "^.*/(\\d{8}).*",
                      Repeat         = "^.*_r(\\d+)_.*",
                      iron           = "^.*_Fe(.).*",
                      timepoint      = "^.*_t(\\d+).*") 

pvd_types <- data.frame(supernatant = c("PAO1",
                                        "PA14_O",
                                        "60.92",
                                        "59.40",
                                        "59.37",
                                        "59.38"),
                        bacteria_code   = c(107,
                                        109,
                                        3330,
                                        3320,
                                        3314,
                                        3319),
                        pvd_type    = c(1,
                                        1,
                                        3,
                                        2,
                                        3,
                                        2))

# Merge and process OD600 data -------------------------------------------------
## unify column names
col_names <- names(raw_OD_data[[1]])

all_OD_raw <- lapply(raw_OD_data,
                     function(x){
                       names(x) <- col_names
                       return(x)
                     }) %>%
  do.call(rbind, .) %>%
  dplyr::select(-Content, -X) %>%
  right_join(well_key,
             by = c("Well.Col" = "col",
                    "Repeat"   = "Repeat",
                    "Well.Row" = "row")) %>%
  left_join(pvd_types,
            by = "supernatant") %>%
  mutate(well   = ifelse(nchar(Well.Col) == 1,
                         paste0(Well.Row, 0, Well.Col),
                         paste0(Well.Row, Well.Col)),
         time_h = factor(timepoint),
         iron     = ifelse(iron == "-",
                           "No Iron",
                           "Iron")) %>%
  rename(OD600  = Raw.Data..600.,
         row    = Well.Row,
         col    = Well.Col) %>%
  filter(!(Repeat %in% repeats_to_exclude))

## Calculate blank values
blanks <- all_OD_raw %>%
  filter(time_h == 0) %>%
  group_by(well, Repeat, iron) %>%
  summarise(blank_OD600 = mean(OD600))

OD600_plot_data <- all_OD_raw %>%
  left_join(blanks,
            by = c("well" = "well",
                   "iron" = "iron",
                   "Repeat" = "Repeat")) %>%
  mutate(blank_corrected_OD600 = ifelse(OD600 - blank_OD600 < 0,
                                        0,
                                        OD600 - blank_OD600),
         Repeat                = factor(Repeat),
         supernatant           = factor(supernatant))

## Percentage inhibition
### Get controls
control_OD600 <- OD600_plot_data %>%
  filter(grepl("fungi in 50:50 SDB:LB control", notes)) %>%
  group_by(Repeat, time_h, iron) %>%                                    
  summarise(raw_control_OD600       = mean(OD600),
            corrected_control_OD600 = mean(blank_corrected_OD600))

## Calculate % inhibition
per_inhibition <- OD600_plot_data %>%
  filter(!grepl("control", notes)) %>%
  left_join(control_OD600,
            by = c("Repeat" = "Repeat",
                   "iron"   = "iron",
                   "time_h" = "time_h")) %>%
  mutate(raw_per_growth       = (OD600/raw_control_OD600) * 100,
         corrected_per_growth = (blank_corrected_OD600/corrected_control_OD600) * 100,
         iron                 = factor(iron, levels = c("No Iron", "Iron")),
         bacteria_code        = factor(bacteria_code),
         grouping             = paste(bacteria_code, iron, sep = ","),
         row                  = factor(row, levels = rev(c("A", "B", "C", "D",
                                                       "E", "F", "G", "H")))) %>%
  filter(time_h == 48)



## Fit models
supernatant_models <-  drm(formula = corrected_per_growth ~ concentration,
                           curveid = c(iron, bacteria_code),
                           data    = per_inhibition,
                           fct     = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "IC50"),
                                          fixed = c(NA, NA, NA, NA)))

## Try fitting model per figure panel to see if iron has a signficant impact
test_model <- drm(formula = corrected_per_growth ~ concentration ,
                  curveid = iron,
                  data    = per_inhibition %>%
                    filter(bacteria_code == 109),
                  fct     = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "IC50"),
                                 fixed = c(NA, NA, NA, NA)))
test_model_null <- drm(formula = corrected_per_growth ~ concentration ,
                  data    = per_inhibition %>%
                    filter(bacteria_code == 109),
                  fct     = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "IC50"),
                                 fixed = c(NA, NA, NA, NA)))

modelFit(test_model)
EDcomp(test_model, c(50))
## Extract IC50s
supernatant_coefficients <- as.data.frame(summary(supernatant_models)[["coefficients"]]) %>%
  rownames_to_column(var = "rowname") %>%
  separate(col  = rowname,
           into = c("parameter", "group"),
           sep  = ":",
           extra = "merge") %>%
  separate(col   = "group",
           into  = c("bacteria_code", "iron"),
           sep   = ",") %>%
  mutate(lower_estimate = Estimate - `Std. Error`,
         upper_estimate = Estimate + `Std. Error`,
         bacteria       = factor(bacteria_code,
                                 levels = pvd_types$bacteria_code),
         iron           = factor(iron, levels = c("No Iron", "Iron"))) %>%
  filter(parameter == "IC50")


## Calculate p_values
supernatant_p_values <- as.data.frame(compParm(object   = supernatant_models,
                                               strVal   = "IC50",
                                               display  = FALSE,
                                               operator = "-")) %>%
  rownames_to_column(var = "comparison") %>%
  separate(col  = comparison,
           into = c("group1", "group2"),
           sep  = "-") %>%
  separate(col  = group1,
           into = c("bacteria1", "iron1"),
           sep  = ",",
           remove = FALSE) %>%
  separate(col  = group2,
           into = c("bacteria2", "iron2"),
           sep  = ",",
           remove = FALSE) %>%
  filter(bacteria1 == bacteria2) %>%
  transmute(group1 = group1,
            group2 = group2,
            p      = `p-value`,
            p.adj  = p.adjust(`p-value`,
                              method = "BH"),
            y.position = c(1, 1, 1, 1, 1, 1)) %>%
  add_significance("p") %>%
  add_significance("p.adj")

## Work out AUC
library(pracma) # AUC functions

### test
x = c(0, 1, 2)
y = c(2, 1, 0)

trapz(x, y)

auc_data <- per_inhibition %>%
  group_by(iron, bacteria_code, Repeat) %>%
  arrange(concentration,
          by.groups = TRUE) %>%
  summarise(per_AUC = trapz(x = concentration,
                            y = corrected_per_growth),
            raw_AUC = trapz(x = concentration,
                            y = blank_corrected_OD600))

# Plot OD600 data --------------------------------------------------------------
OD600_supernatant_plot <- ggplot(data    = per_inhibition %>%
                                   mutate(supernatant = factor(bacteria_code,
                                                               levels = c("107",
                                                                          "109",
                                                                          "3319",
                                                                          "3320",
                                                                          "3314",
                                                                          "3330"))),
                                 mapping = aes(x     = concentration,
                                               y     = corrected_per_growth,
                                               color = iron)) +
  geom_point() +
  stat_smooth(method      = drm,
              method.args = list(fct    = drc::LL.4(fixed = c(NA, NA, NA, NA)),
                                 lowerl = c(NULL, NULL, NULL, NULL)),
              mapping     = aes(color = iron),
              formula     = y ~ x,
              se          = FALSE,
              fullrange   = TRUE,
              show.legend = FALSE) +
  stat_summary(geom        = "errorbar",
               position    = position_dodge(width = 0.5),
               width       = 1,
               fun.data    = "mean_cl_boot",
               fun.args    = list(conf.int = 0.95),
               show.legend = FALSE) +
  stat_summary(geom        = "point",
               size        = 2,
               position    = position_dodge(width = 0.5),
               fun         = "mean",
               show.legend = TRUE) +
  stat_compare_means(method = "t.test",
                     show.legend = FALSE,
                     hide.ns = TRUE,
                     label.y = 125,
                     mapping = aes(label = after_stat(p.signif))) +
  labs(title = "") +
  scale_x_continuous(name = "\nSupernatant Concentration (%)",
                     limits = c(0, 100),
                     breaks = seq(0, 100, 10),
                     expand = c(0, 0)) +
  scale_y_continuous(name   = "Percentage of\nControl OD600",
                     breaks = seq(0, 150, 25),
                     expand = c(0, 0)) +
  scale_color_discrete(name = "Iron (μM)",
                       labels = c("0",
                                  "100")) +
  coord_cartesian(ylim  = c(0, 150),
                  xlim  = c(0, 55)) +
  facet_wrap(vars(bacteria_code),
             ncol = 2,
             nrow = 3) +
  my_theme +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_blank(),
        axis.title.y = element_text(angle = 0,
                                    vjust = 0.5))

OD600_supernatant_plot 

ggsave(filename = "figure_5_suparnatant_v_OD600_curves.png",
       plot   = OD600_supernatant_plot,
       path   = "./data/figures",
       width  = 20,
       height = 20,
       units = "cm")


# Plot AUC data ----
auc_plot <- ggplot(data = auc_data,
                   mapping = aes(x = bacteria_code,
                                 y = raw_AUC,
                                 color = iron)) +
  geom_point(position = position_dodge(width = 1)) +
  stat_summary(geom        = "errorbar",
               position    = position_dodge(width = 1),
               width       = 1,
               fun.data    = "mean_cl_boot",
               fun.args    = list(conf.int = 0.95),
               show.legend = FALSE) +
  stat_summary(geom        = "point",
               size        = 2,
               position    = position_dodge(width = 1),
               fun         = "mean",
               show.legend = TRUE) +
  stat_compare_means(show.legend = FALSE,
                     method      = "t.test",
                     paired      = TRUE,
                     mapping = aes(label = round(after_stat(p.adj), 2))) +
  coord_cartesian(ylim = c(0, 5100)) +
  my_theme

auc_plot
# Plot IC50 values ----
supernatant_IC50_plot <- ggplot(data = supernatant_coefficients,
                                mapping = aes(x = bacteria,
                                              y = Estimate,
                                              ymin = Estimate - `Std. Error`,
                                              ymax = Estimate + `Std. Error`,
                                              color = iron)) +
  geom_errorbar(position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_y_continuous(name = "IC50") +
  labs(caption = "Points are estimate \u00b1 standard error") +
  coord_cartesian(ylim = c(0, 100)) +
  #facet_grid(cols = vars(siderotype)) + # TODO add this
  #coord_flip(ylim = c(0, 100)) +
  my_theme

supernatant_IC50_plot


ggsave(filename = "supernatant_ic50s.png",
       plot   = supernatant_IC50_plot,
       path   = "./data/graphs",
       width  = 15,
       height = 15,
       units = "cm")

# Plot pvalues ----
p_value_data <- supernatant_p_values %>%
  mutate(comparison = factor(paste(group1, " v ", group2),
                             levels = c("PAO1,No Iron  v  PAO1,Iron",
                                        "PA14_O,No Iron  v  PA14_O,Iron",
                                        "59.38,No Iron  v  59.38,Iron",
                                        "59.40,No Iron  v  59.40,Iron",
                                        "59.37,No Iron  v  59.37,Iron",
                                        "60.92,No Iron  v  60.92,Iron")))


super_OD_p <- ggplot(data    = p_value_data, 
                     mapping = aes(x     = comparison,
                                   y     = Estimate)) +
  stat_pvalue_manual(data = p_value_data,
                     label = "p.adj.signif",
                     hide.ns = F)
super_OD_p

super_OD_heatmap_with_p <- OD600_heatmap_plot + super_OD_p +
  plot_layout(design = "12",
              widths = c(0.8, 0.2))

super_OD_heatmap_with_p

# Plot just OD600 at 50% supernatant
supernatant_bar_plot <- ggplot(data = per_inhibition %>%
                                 filter(concentration == 50),
                               mapping = aes(x = supernatant,
                                             y = corrected_per_growth,
                                             color = iron)) +
  stat_summary(geom        = "errorbar",
               position    = position_dodge(width = 1),
               width       = 1,
               fun.data    = "mean_cl_boot",
               fun.args    = list(conf.int = 0.95),
               show.legend = FALSE) +
  stat_summary(geom        = "point",
               size        = 2,
               position    = position_dodge(width = 1),
               fun         = "mean",
               show.legend = TRUE) +
  geom_point(position    = position_dodge(width = 1)) +
  stat_compare_means(method = "t.test",
                     label  = "p.signif") +
  my_theme

supernatant_bar_plot

ggsave(filename = "supernatant_50_per.png",
       plot   = supernatant_bar_plot,
       path   = "./data/graphs",
       width  = 15,
       height = 15,
       units = "cm")


# Plot model fits on a heatmap to make it look nicer ---------------------------
prediction_data <- per_inhibition %>%
  select(grouping, concentration)


prediction_data$prediction <- predict(object  = supernatant_models,
                                      newdata = prediction_data)

prediction_plot_data <- prediction_data %>%
  separate_wider_delim(cols = grouping,
                       delim = ",",
                       names = c("Bacteria", "Iron"))

pred_heatmap_plot <- ggplot(data    = prediction_plot_data, 
                             mapping = aes(x     = concentration,
                                           y     = Iron,
                                           fill  = prediction,
                                           label = round(prediction))) +
  geom_tile() + 
#  geom_text(color = "white") +
  scale_fill_gradient2(high      = "black",
                       mid       = "purple",
                       low       = "white",
                       midpoint  = 100,
                       name      = "Growth (%)",
                       transform = "reverse") +
  scale_y_discrete(name = NULL,
                   expand = c(0,0)) +
  scale_x_continuous(name   = "\nSupernatant Concentration",
                     expand = c(0,0)) +
  labs(title = "\n") +
  facet_grid(rows = vars(Bacteria),
             switch   = "y") +
  my_theme +
  theme(legend.position      = c(0, 1),
        legend.justification = c(0, 0),
        legend.direction     = "horizontal",
        legend.key.width     = unit(25, "pt"),
        legend.key.height    = unit(10, "pt"),
        legend.spacing.x     = unit(2, "pt"),
        legend.text          = element_text(margin = margin(0, 0, 0, 0, "pt"), size = 10),
        panel.spacing.y      = unit(0, "cm"))
 
pred_heatmap_plot

# fir moedels for each individual one 







