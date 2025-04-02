#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

# Read input arguments
args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
plot_dir <- args[2]
trait <- args[3]
my_functions_R <- args[4] 

fdr_threshold <- 0.05
value_threshold <- 1.0E-99 # Define the threshold for any value
# Read the GWAS results file
# Read the GWAS results file with blanks as missing data
gwas_data <- read.csv(input_file, na.strings = c(""," ","0"))
# Replace any value below the threshold with NA
# Function to create Manhattan plot
plot_manhattan <- function(data, pval_column, title, threshold.SIG=NULL, threshold.SUG=NULL) {
  # Select relevant columns
  gwas_data <- data %>% select(Chromosome, Position, !!sym(pval_column)) %>%
    rename(log_p = !!sym(pval_column))
  gwas_data$log_p[gwas_data$log_p  < value_threshold] <- NA
  gwas_data <- na.omit(gwas_data)  
  gwas_data$log_p <- -log10(gwas_data$log_p)

  # Calculate the cumulative BP position for all chromosomes
  data_cum <- gwas_data %>% 
    group_by(Chromosome) %>% 
    summarise(max_Position = max(Position)) %>% 
    mutate(Position_add = lag(cumsum(as.numeric(max_Position)), default = 0)) %>% 
    select(Chromosome, Position_add)

  gwas_data <- gwas_data %>% 
    inner_join(data_cum, by = "Chromosome") %>% 
    mutate(Position_cum = Position + Position_add)

  axis_set <- gwas_data %>% 
    group_by(Chromosome) %>% 
    summarize(center = mean(Position_cum))

  # Set Y-axis limits and threshold
  ylims.min <- 0 
  if (is.null(threshold.SIG)) {
    threshold.SIG <- -log10( 0.05/length(unique(gwas_data$log_p)))

  }
    if (is.null(threshold.SUG)) {
     threshold.SUG  <- -log10(1/length(unique(gwas_data$log_p)))

  }  
  max.effect <- max(abs(max(gwas_data$log_p)),  threshold.SIG) + 1
  ylims.max <- round(max.effect)
  print(ylims.max)
  # Create Manhattan plot
  colors_UF <- c("#FA4616", rep(c("#4169E1", "#0021A5", "#FA4616"), 9), "#4169E1")
  manhattan_plot <- ggplot(gwas_data, aes(x = Position_cum, y = log_p)) +
    geom_point(aes(color = as.factor(Chromosome)), size = 1, alpha = 0.75) +
    scale_x_continuous(label = axis_set$Chromosome, breaks = axis_set$center) + 
    scale_color_manual(values = colors_UF) +
    scale_y_continuous(expand = c(0,0), breaks = seq(from=ylims.min, to=ylims.max, by = 2),
                       limits = c(ylims.min, ylims.max)) +
    geom_hline(yintercept = threshold.SIG, color ="red", linetype = "dashed") +
    geom_hline(yintercept = threshold.SUG, color = "green", linetype = "dashed") + 
    theme_bw() +
    labs(title = title,
         x = "Chromosome",
         y = "-log(p-value)") +
    theme(legend.position = "none",
          panel.spacing.y = unit(1, "lines"),
          axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 9, hjust = 1),
          axis.text.y = element_text(size = 10))
  
  return(manhattan_plot)
}

# Create Manhattan plots for each p-value type
pval_types <- c("Reduced_SNP_P_value", "Reduced_BOA_P_value", "Reduced_Combined_P_value", 
                 "SNP_vs_Null_P_value", "BOA_vs_Null_P_value")
manhattan_plots <- list()

for (pval_type in pval_types) {
  title <- paste(trait, pval_type, "GWAS Results")
  manhattan_plots[[pval_type]] <- plot_manhattan(gwas_data, pval_type, title)
}
Ind_pval  <- c("Reduced_Combined_P_value",
                 "SNP_vs_Null_P_value", "BOA_vs_Null_P_value")
# Step 3: Save individual Manhattan plots
for (pval_type in Ind_pval) {
  plot_file <- file.path(plot_dir, paste0(pval_type, "_manhattan_plot.png"))
  ggsave(plot_file, manhattan_plots[[pval_type]], width=10, height=6)
}

plot_manhattan_mixed <- function(data, pval_columns, title, threshold.SIG=NULL, threshold.SUG=NULL) {
  # Select relevant columns
  # Manually combine columns using bind_rows
  gwas_data_list <- list()
  for (col in pval_columns) {
    temp_data <- data[, c("SNP_ID", "Chromosome", "Position", col)]
    temp_data$Trait <- col
    names(temp_data)[4] <- "log_p"
    gwas_data_list[[col]] <- temp_data
  }
  gwas_data <- do.call(rbind, gwas_data_list)
  gwas_data$log_p[gwas_data$log_p  < value_threshold ] <- NA
  gwas_data <- gwas_data[!is.na(gwas_data$log_p), ]
  gwas_data$log_p <- -log10(gwas_data$log_p)

  # Calculate the cumulative BP position for all chromosomes
  data_cum <- gwas_data %>% 
    group_by(Chromosome) %>% 
    summarise(max_Position = max(Position)) %>% 
    mutate(Position_add = lag(cumsum(as.numeric(max_Position)), default = 0)) %>% 
    select(Chromosome, Position_add)

  gwas_data <- gwas_data %>% 
    inner_join(data_cum, by = "Chromosome") %>% 
    mutate(Position_cum = Position + Position_add)

  axis_set <- gwas_data %>% 
    group_by(Chromosome) %>% 
    summarize(center = mean(Position_cum))

  # Set Y-axis limits and threshold
  ylims.min <- 0 
  if (is.null(threshold.SIG)) {
    threshold.SIG <- -log10(0.05/nrow( data))

  }
    if (is.null(threshold.SUG)) {
     threshold.SUG  <- -log10(1/nrow( data))

  }  
  max.effect <- max(abs(max(gwas_data$log_p)),  threshold.SIG) + 1
  ylims.max <- round(max.effect)
  print(ylims.max)
  # Create Manhattan plot
  colors_UF <- c("#FA4616", rep(c("#4169E1", "#0021A5", "#FA4616"), 9), "#4169E1")
  manhattan_plot <- ggplot(gwas_data, aes(x = Position_cum, y = log_p)) +
    geom_point(aes(color = as.factor(Chromosome)), size = 1, alpha = 0.75) +
    scale_x_continuous(label = axis_set$Chromosome, breaks = axis_set$center) + 
    scale_color_manual(values = colors_UF) +
    scale_y_continuous(expand = c(0,0), breaks = seq(from=ylims.min, to=ylims.max, by = 2),
                       limits = c(ylims.min, ylims.max)) +
    geom_hline(yintercept = threshold.SIG, color ="red", linetype = "dashed") +
    geom_hline(yintercept = threshold.SUG, color = "green", linetype = "dashed") + 
    theme_bw() +
    labs(title = title,
         x = "Chromosome",
         y = "-log(p-value)") +
    theme(legend.position = "none",
          panel.spacing.y = unit(1, "lines"),
          axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 9, hjust = 1),
          axis.text.y = element_text(size = 10))
  
  return(manhattan_plot)
}


manhattan_plot_model1 <- plot_manhattan(gwas_data, c("SNP_vs_Null_P_value"),
                                            paste(trait," Model 1 GWAS Results"))

manhattan_plot_model2 <- plot_manhattan(gwas_data, c("BOA_vs_Null_P_value"),
                                            paste(trait," Model 2 GWAS Results"))

manhattan_plot_model3 <-plot_manhattan_mixed(gwas_data, c("Reduced_SNP_P_value", "Reduced_BOA_P_value", "Reduced_Combined_P_value"),
                                            paste(trait," Model 3 GWAS Results"))

# Step 4: Stack the plots vertically
stacked_plot <- grid.arrange(manhattan_plot_model1,
                             manhattan_plot_model2,
                             manhattan_plot_model3, 
                             ncol=1)

# Save the stacked plot
stacked_plot_file <- file.path(plot_dir, paste0(trait, "_Main_stacked_manhattan_plot.png"))
ggsave(stacked_plot_file, stacked_plot, width=8, height=8)

# Step 4: Stack the plots vertically
stacked_plot_Reduced_Model <- grid.arrange(manhattan_plots[["Reduced_Combined_P_value"]],
                              manhattan_plots[["Reduced_SNP_P_value"]], 
                              manhattan_plots[["Reduced_BOA_P_value"]], 
                             ncol=1)

# Save the stacked plot
stacked_plot_file_Reduced_Model  <- file.path(plot_dir, paste0(trait, "_Reduced_Model_manhattan_plot.png"))
ggsave(stacked_plot_file_Reduced_Model , stacked_plot_Reduced_Model  , width=9, height=16)
