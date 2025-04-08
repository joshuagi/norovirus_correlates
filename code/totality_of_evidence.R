library(dplyr)
library(openxlsx)
library(tidyr)
library(ggplot2)
library(lemon)
path = "figures/" # Store results

# Load data
immunogenicity <-  read.xlsx("data/Data file S1.xlsx", sheet = 3)
endpoints_raw_data <- read.xlsx("data/Data file S1.xlsx", sheet = 2)

endpoints_raw_data_Primary <- endpoints_raw_data %>%
  select(SUBJID, Treatment, `qPCR+`, `qPCR+AGE+`) %>%
  mutate(Treatment = case_when(Treatment == "VXA-G1.1-NN" ~ "VXA",  
                               TRUE ~ Treatment))
immunogenicity_raw_data_Primary <- immunogenicity %>%
  select(SUBJID, Treatment, ASCRaw_D8, NBAA_D28, SIgA_D28, SIgG_D28) %>%
  mutate_at(vars(NBAA_D28, SIgA_D28, SIgG_D28), funs(log10(.)))

primary_endpoints <- full_join(endpoints_raw_data_Primary, immunogenicity_raw_data_Primary, by = c("SUBJID", "Treatment")) %>%
  mutate(Treatment = factor(Treatment, levels = c("VXA", "Placebo")))


# Define function to calculate test statistics
estimateDifference <- function(x, iteration = 1) {
  result <- x %>%
    pivot_longer(
      cols = -c(SUBJID, Treatment),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(type_test = case_when(variable %in% c("qPCR+", "qPCR+AGE+") ~ "proportion",
                                 TRUE ~ "continuous")) %>%
    split(.$variable) %>%
    purrr::map(., function(df) {
      type_test <- unique(df$type_test)
      variable <- unique(df$variable)
      
      if(type_test == "proportion") {
        prop.table <- df %>%
          group_by(Treatment) %>%
          dplyr::summarise(counts = sum(value),
                           total = n())
        result <- prop.test(x = prop.table$counts, n = prop.table$total, correct = F, alternative = "two.sided")
        result <- tibble(estimate = result$estimate[2] - result$estimate[1], # Change sign so that improvement in VXA is a positive quantity
                         lower = -result$conf.int[2], 
                         upper = -result$conf.int[1], 
                         p.value = result$p.value,
                         score = sqrt(result$statistic)) %>%
          mutate(variable = variable)
        
      } else if (type_test == "continuous") {
        result <- stats::t.test(df$value ~ df$Treatment, alternative = "two.sided", var.equal = FALSE)
        result <- tibble(estimate = result$estimate[1] - result$estimate[2], # Change sign so that improvement in VXA is a positive quantity
                         lower = result$conf.int[1], 
                         upper = result$conf.int[2], 
                         p.value = result$p.value,
                         score = result$statistic)%>%
          mutate(variable = variable)
        
      }
      
    }) %>%
    bind_rows() %>%
    mutate(iteration = iteration)
  return(result)
}

# Calculate test statistics
True_estimates <- estimateDifference(primary_endpoints) # Material for table and forest plot
True_estimates <- True_estimates %>%
  mutate(group = case_when(variable %in% c("qPCR+", "qPCR+AGE+") ~ "efficacy",
                           variable %in% c("SIgA_D28", "SIgG_D28", "NBAA_D28") ~ "antibody",
                           TRUE ~variable)) %>%
  mutate(group = factor(group, levels = c("efficacy", "antibody", "NBAA_D28", "ASCRaw_D8"))) %>%
  mutate(variable = factor(variable, levels = rev(c("qPCR+", "qPCR+AGE+", "SIgA_D28", "SIgG_D28", "NBAA_D28", "ASCRaw_D8")))) %>%
  mutate(estimate = case_when(group == "efficacy" ~ estimate*100,
                              group %in% c("antibody") ~ 10^estimate,
                              TRUE ~ estimate)) %>%
  mutate(lower = case_when(group == "efficacy" ~ lower*100,
                           group %in% c("antibody") ~ 10^lower,
                           TRUE ~ lower)) %>%
  mutate(upper = case_when(group == "efficacy" ~ upper*100,
                           group %in% c("antibody") ~ 10^upper,
                           TRUE ~ upper)) %>%
  mutate(`Estimate (95% CI)` = paste0(round(estimate,2), 
                                      " (",
                                      round(lower,2),
                                      " to ",
                                      round(upper,2),
                                      ")"))


# Generate null distribution
set.seed(2024)
iter <- 10000
perm_result <- list()
primary_endpoints_perm <- primary_endpoints
for(i in 1:iter) {
  primary_endpoints_perm$Treatment <- sample(primary_endpoints$Treatment, size = nrow(primary_endpoints), replace = FALSE)
  perm_result[[i]] <- estimateDifference(primary_endpoints_perm, iteration = i)
}

perm_result_scores <- perm_result %>%
  bind_rows() %>%
  group_by(iteration) %>%
  dplyr::summarise(mean = mean(score)) 


# Calculate p.value for plotting
test_value <- mean(True_estimates$score)
p_value <- (sum(perm_result_scores$mean >= test_value) + 1) / (nrow(perm_result_scores) + 1) # One-sided p.value, with continuity correction / pseudocount adjustment
p_value < 0.0001
plot_Z <- data.frame(label = paste0("Mean Z = ", round(mean(True_estimates$score), 2)),
                     x = 2,
                     y = 300)
plot_p <- data.frame(label = paste0("p.value = <0.0001"),
                     x = 2,
                     y = 200)

totality_histogram <- perm_result_scores %>%
  ggplot(aes(x = mean)) +
  geom_histogram(fill = "darkgrey", colour = "white", bins = 50, size = 0.1) +
  geom_vline(xintercept = mean(True_estimates$score), linetype = "dashed", size = 0.25) +
  geom_text(data = plot_Z,
            aes(x = x, y = y, label = label), hjust = 0, size=5 / (14/5))+
  geom_text(data = plot_p,
            aes(x = x, y = y, label = label), hjust = 0, size=5 / (14/5)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        strip.text = element_text(size = 5),
        axis.title.x = element_text(size = 5),
        axis.title.y = element_text(size = 5), 
        axis.ticks = element_line(colour = "black", size = 0.2),
        panel.border = element_rect(fill=NA, colour = "black", size=0.2),
        panel.spacing = unit(0.05, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.2, "cm"), 
        strip.background = element_blank(), 
        strip.clip = "off",
        panel.grid = element_blank(),
        strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) + # Space at the top but not at the bottom of the y-axis
  labs(y = "Frequency", x = "Simulated mean z-score")


ggsave(plot = totality_histogram,
       file = paste("totality_histogram", ".pdf", sep=""),
       path = path,
       width = 5,
       height = 4,
       units = "cm",
       device = "pdf",
       useDingbats = FALSE)



