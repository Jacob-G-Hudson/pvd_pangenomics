#  IC50 analysis
# Dependencies -----------------------------------------------------------------
library(tidyverse)
library(ggpubr)           # statistics for graphs
library(drc)              # DRC curve fit functions
library(rstatix)          # Statistical helper functions
library(patchwork)        # Multiplot functions

# Setup ------------------------------------------------------------------------
## ggplot
my_theme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Load data ----
per_inhibition <- read.csv(file = "./raw_data/ic50/compiled_results.csv")

# Compute f-test p-values ----
f_test <- function(data, x, y, group){
  
  # Make models
  reduced_model <- drm(formula = as.formula(paste(x, "~", y)),
                       data    = data,
                       curveid = data[[group]],
                       fct     = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "IC50"),
                                      fixed = c(NA, 0, NA, NA)),
                       pmodels = list( ~ 1, ~ 1, ~ 1, ~data[[group]]-1),  # model is reduced because these params are contrained
                       upperl =  c(+Inf, NULL, 225, +Inf))
  
  full_model   <- drm(formula = as.formula(paste(x, "~", y)),
                      data    = data,
                      curveid = data[[group]],
                      fct     = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "IC50"),
                                     fixed = c(NA, 0, NA, NA)),
                      upperl =  c(+Inf, NULL, 225, +Inf))
  
  
  # F-test
  invisible(capture.output(anova_results <- anova(reduced_model, full_model)))
  
  # Get IC50s
  invisible(capture.output(reduced_IC50s <- as.data.frame(ED(reduced_model, 50))))
  invisible(capture.output(full_IC50s <- as.data.frame(ED(full_model, 50))))
  
  iron_IC50s <- c(reduced_IC50s[1, 1],
                  full_IC50s[1, 1])
  
  iron_IC50s_SEM <- c(reduced_IC50s[1, 2],
                      full_IC50s[1, 2])
  
  no_iron_IC50s <- c(reduced_IC50s[2, 1],
                     full_IC50s[2, 1])
  
  no_iron_IC50s_SEM <- c(reduced_IC50s[2, 2],
                         full_IC50s[2, 2])
  
  
  # Save for output
  anova_results_df <- data.frame(models   = c("reduced", "full"),
                                 model_df = anova_results$ModelDf,
                                 RSS      = anova_results$RSS,
                                 Df       = anova_results$Df,
                                 F_value  = anova_results$`F value`,
                                 p_value  = anova_results$`p value`,
                                 iron_IC50 = iron_IC50s,
                                 iron_IC50_SEM = iron_IC50s_SEM,
                                 no_iron_IC50 = no_iron_IC50s,
                                 no_iron_IC50_SEM = no_iron_IC50s_SEM)
  
  return(anova_results_df)
  
}

f_test_results <- lapply(X = split(per_inhibition ,
                                   ~ bacteria_code),
                         FUN = f_test,
                         x    = "corrected_per_growth",
                         y    = "concentration",
                         group = "iron") %>%
  bind_rows(.id = "bacteria_code") %>%
  mutate(p.adj = p.adjust(p = p_value,
                          method = "BH",
                          n = 6))

IC50_range <- f_test_results %>%
  filter(models == "full")

# Plot curves where curve is calcalted from the significant model ----
for(bacteria_code in unique(per_inhibition$bacteria_code)){
  
  writeLines(paste(bacteria_code))
  
  # Filter to bacteria we are plotting
  tmp_data <- per_inhibition %>%
    filter(bacteria_code == UQ(bacteria_code))
  
  # Get p-value 
  tmp_p_value <- f_test_results %>%
    filter(bacteria_code == UQ(bacteria_code)) %>%
    na.omit() %>%
    transmute(p.adj   = p.adj,
              p.code = case_when(p.adj <  0.0001 ~ "****",
                                 p.adj <= 0.001 ~ "***",
                                 p.adj <= 0.01  ~ "**",
                                 p.adj <= 0.05  ~ "*",
                                 .default = "ns"),
              x     = 3,
              y     = 225)
  
  
  # Re-calculate the best model 
  tmp_model <- drm(formula = corrected_per_growth ~ concentration,
                   data    = tmp_data,
                   curveid = iron,
                   fct     = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "IC50"),
                                  fixed = c(NA, 0, NA, NA)),
                   upperl =  c(+Inf, NULL, 225, +Inf))
  
  
  # fit curve
  tmp_new_data <- expand.grid(concentration = seq(0, 60, length.out = 200),
                              iron          = unique(tmp_data$iron))
  
  tmp_prediction <- predict(object  = tmp_model, 
                            newdata = tmp_new_data,
                            interval = "confidence")
  
  tmp_new_data$fit_lower <- tmp_prediction[, "Lower"]
  tmp_new_data$fit       <- tmp_prediction[, "Prediction"]
  tmp_new_data$fit_upper <- tmp_prediction[, "Upper"]
  
  # plot curve
  tmp_plot <- ggplot(data    = tmp_data,
                     mapping = aes(x = concentration,
                                   y = corrected_per_growth,
                                   color = iron)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_line(data    = tmp_new_data,
              mapping = aes(x     = concentration,
                            y     = fit,
                            color = iron),
              linewidth = 0.6) +
    geom_point(size = 0.75)  +
    stat_summary(geom        = "errorbar",
                 width       = 2,
                 linewidth   = 0.3,
                 fun.data    = mean_se,
                 show.legend = FALSE) +
    stat_summary(geom        = "point",
                 size        = 2,
                 fun         = "mean",
                 show.legend = legend) +
    geom_text(data    = tmp_p_value,
              mapping = aes(x = x, 
                            y = y,
                            label = p.code),
              inherit.aes = FALSE,
              size        = 8/.pt) +
    scale_x_continuous(name  = x_axis_title,
                       limits = c(0, 100),
                       breaks = seq(0, 100, 10),
                       minor_breaks = NULL,
                       expand = c(0, 0)) +
    scale_y_continuous(name   = NULL,
                       breaks = seq(0, 250, 50),
                       minor_breaks = NULL) +
    scale_color_manual(name   = "Iron (μM)",
                       labels = c("0",
                                  "100"),
                       values = c("#E66100",
                                  "#5D3A9B")) +
    facet_grid(cols = vars(bacteria_code),
               rows = vars(pvd_type)) +
    coord_cartesian(ylim  = c(0, 250),
                    xlim  = c(0, 55)) +
    theme_minimal() +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_blank(),
          axis.title.y = element_text(angle = 0,
                                      vjust = 0.5),
          text = element_text(size = 8))
  
  tmp_plot
  
  assign(x = paste0("bacteria_", bacteria_code, "_plot"),
         value = tmp_plot)
  
}

dose_response_curves <- (bacteria_107_plot  + bacteria_109_plot)  /
  (bacteria_3319_plot + bacteria_3320_plot) /
  (bacteria_3314_plot + bacteria_3330_plot)

ggsave(filename = "figure_4_suparnatant_v_OD600_curves.svg",
       plot   = dose_response_curves,
       path   = "./data/figures",
       width  = 21- 2.54 - 2.54,
       height = 20,
       units = "cm")

# Get some stats ----
mean_inhibitions <- per_inhibition %>%
  group_by(bacteria_code, iron, concentration) %>%
  summarise(mean_growth = mean(corrected_per_growth))

