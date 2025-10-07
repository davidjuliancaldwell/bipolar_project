# Load libraries

library(easystats)
library(readxl)
library(dplyr)
library(ggplot2)
library(here)

# Load data

base_dir = here()

data <- read_excel(file.path(base_dir, "LastScripts", "band_power_output.xlsx"))
data <- data %>% mutate(patient = as.factor(patient),
                       freq = as.factor(freq),
                       bpd = as.factor(bpd),
                       power = as.double(power))

data <- data %>% group_by(freq, patient) %>%
mutate(baseline = power[bpd == 0]) %>%
mutate(power_ref = power - baseline) 

data_no_baseline <- data %>%
  filter(bpd != 0)  

# Run stats


data_summary <- data_no_baseline %>%
  group_by(freq, bpd) %>%
  summarise(median_power_ref = median(power_ref),
            IQR_power_ref = IQR(power_ref),
            p_value = wilcox.test(power_ref, mu=0)$p.value,
            effect_size = rank_biserial(power_ref, mu = 0)$r)
# Adjust p-values for multiple comparisons

data_summary <- data_summary %>%
  mutate(p_value_correct = p.adjust(p_value, method="BH"))

data_summary <- data_summary %>%
  mutate(significant = ifelse(p_value_correct < 0.05, "Yes", "No"))

# plot
    ggplot(data_no_baseline, aes(x=bpd, y=power_ref)) +
    geom_violin(position=position_dodge(0.8), width=0.7,trim=FALSE) +
    geom_boxplot(width=0.1) +
    scale_y_continuous() +
    theme_minimal() +
    labs(title="Band Power by Frequency and BPD Status",
        x="Reference Band Distance",
        y="Band Power") +
    theme(legend.position="top") +
    facet_wrap(~ freq) 
    
    # plot
    ggplot(data_summary, aes(x=bpd, y=effect_size)) +
    geom_point() +
    scale_y_continuous() +
    theme_minimal() +
    labs(title="Effect Size by Frequency and BPD Status",
        x="Frequency Band",
        y="Effect Size") +
    theme(legend.position="top") +
    facet_wrap(~ freq) 