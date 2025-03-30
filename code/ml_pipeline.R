library(purrr)
library(rstatix)
library(openxlsx)
library(tidyr)
library(stringr)
library(dplyr)
library(tune)
library(tidymodels)
library(vip)
library(pROC)
library(ggpubr)
library(gridExtra)
library(doParallel)
library(glmnet)
library(ranger)


# Import functions
source("code/pipeline.utils.R")

# Set paths
result_path = "result/" 
figure_path = "figures/"

# set up parallel backend for grid search
cl <- makePSOCKcluster(parallel::detectCores(logical = FALSE))
registerDoParallel(cl)


# Load data
immunogenicity <-  read.xlsx("data/immunogenicity_raw_data.xlsx")
endpoints_raw_data <- read.xlsx("data/endpoints_raw_data.xlsx")



# Filter and transform predictor data
immunogenicity_long <- immunogenicity %>%
  #melt(., id.vars = c("SUBJID", "Treatment")) %>%
  pivot_longer(
    cols = -c(SUBJID, Treatment),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(Day = str_split(variable, pattern = "_", simplify = T)[,2]) %>%
  mutate(variable = str_split(variable, pattern = "_", simplify = T)[,1]) %>%
  
  filter(Day %in% c("D28", "D8")) %>%
  mutate(value = case_when(variable == "ASCRaw" ~ value + 1, # add pseudocount to ASCs
                           TRUE ~ value)) %>%
  mutate(value = case_when(variable %in% c("ASCRaw", "FecalIgAVP1onTotal", "NBAA", "SIgA", "SIgG", "NIgAGI1", "SALIgAGI1") ~log10(value), # log10 transform predictors
                           TRUE ~ value)) %>%
  mutate(variable = paste0(variable, "_", Day) ) %>%
  select(-Day) %>%
  pivot_wider(
    names_from = variable,
    values_from = value
  ) %>%
  # join with response
  left_join(., endpoints_raw_data[,c(1,4)]) %>%
  dplyr::rename(response = `qPCR+`) %>%
  select(1:2, response, 3:ncol(.)) %>%
  mutate(response = factor(response, levels = c(0,1))) 


# Run the analysis pipeline
immunogenicity_long %>%
  # Analyze per group
  split(.$Treatment) %>%
  purrr::map(., function(x) {
    Treatment_i <- unique(x$Treatment)
    result <- list()
    i<-1
    for (i in 1:100) {
      print(i)
      result[[i]] <- evaluateModels(x, iteration = i, prop = 0.64, weighted = T)
    }
    file <- paste0(result_path, Sys.Date(), "_", Treatment_i, "_", "NVinfection", "_predictions.rds")
    saveRDS(result, file = file)
  })

# Load results
resultVXA <- readRDS("result/2024-08-25_VXA_NVinfection_predictions.rds")
resultPlacebo <- readRDS("result/2024-08-25_Placebo_NVinfection_predictions.rds")


# Analyze and generate plots
NVinfection_resultVXA <- plotModels(resultVXA)
NVinfection_resultPlacebo <- plotModels(resultPlacebo)


# Lasso
p1 <- NVinfection_resultVXA$multivariable_plots$lasso_reg %>% 
  annotate_figure(fig.lab = "A (VXA)", fig.lab.face = "bold", fig.lab.size = 8, fig.lab.pos = c("top.left"))

p2 <- NVinfection_resultPlacebo$multivariable_plots$lasso_reg %>% 
  annotate_figure(fig.lab = "B (Placebo)", fig.lab.face = "bold", fig.lab.size = 8, fig.lab.pos = c("top.left"))

plot_tmp <- grid.arrange(p1, p2, ncol = 1)
ggsave(plot = plot_tmp,
       file = paste("NVInfection_Lasso", ".pdf", sep=""),
       path = figure_path,
       width = 11,
       height = 9,
       units = "cm",
       device = "pdf",
       useDingbats = FALSE)

# Random Forest
p1 <- NVinfection_resultVXA$multivariable_plots$rf %>% 
  annotate_figure(fig.lab = "A (VXA)", fig.lab.face = "bold", fig.lab.size = 8, fig.lab.pos = c("top.left"))

p2 <- NVinfection_resultPlacebo$multivariable_plots$rf %>% 
  annotate_figure(fig.lab = "B (Placebo)", fig.lab.face = "bold", fig.lab.size = 8, fig.lab.pos = c("top.left"))

plot_tmp <- grid.arrange(p1, p2, ncol = 1)
ggsave(plot = plot_tmp,
       file = paste("NVInfection_RF", ".pdf", sep=""),
       path = figure_path,
       width = 11,
       height = 9,
       units = "cm",
       device = "pdf",
       useDingbats = FALSE)


