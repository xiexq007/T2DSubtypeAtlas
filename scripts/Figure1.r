
# Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(transplantr)
library(cliot)

# Load pre-processed clinical data for all donors
GROUP_LEVELS <- c("ND", "MARD", "MOD", "SIDD", "SIRD")
GROUP_COLS <- c('#999898', '#87a7d5', '#4b5988', '#f1c5a5', '#c03f3f')
STAT_COMPARISONS <- list(c("ND", "MARD"), c("ND", "MOD"), c("ND", "SIDD"), c("ND", "SIRD"))

all_lab <- read.csv("./ND_T2DAll_infor.csv")
all_lab$kclusterGroup <- factor(all_lab$kclusterGroup, levels = GROUP_LEVELS)

plot_clinical_boxplot <- function(df, y_var, y_label) {
  p <- ggplot(df, aes(x = kclusterGroup, y = .data[[y_var]], fill = kclusterGroup)) +
    geom_boxplot(color = "black", size = 0.6, width = 0.7, outlier.shape = 16, outlier.alpha = 0.5) +
    scale_fill_manual(values = GROUP_COLS) +
    theme_classic() +
    labs(x = "", y = y_label) +
    stat_compare_means(
      method = 'wilcox.test', 
      comparisons = STAT_COMPARISONS, 
      label = "p.signif", 
      method.args = list(exact = FALSE)
    ) +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black", size = 12),
      axis.title.y = element_text(face = "bold")
    )
  return(p)
}

# Figure 1G-I ----------

p_gsis <- plot_clinical_boxplot(all_lab, "gsir_si", "GSIS")
p_wbc  <- plot_clinical_boxplot(all_lab, "WBC.THO.uL.", "WBC (THO/uL)")
p_cre  <- plot_clinical_boxplot(all_lab, "Creatinine.mg.dL.", "Creatinine (mg/dL)")

# Figure 1J ----------
# Kidney Function (eGFR and CKD Stages)

all_lab$ethnicity <- ifelse(all_lab$race %in% "African American","black","non-black")

all_lab <- all_lab %>%
  mutate(
    eGFR_epi = ckd_epi(creat = Creatinine.mg.dL., age = age_years, sex = gender, eth = ethnicity, units = "US")
  ) %>%
  mutate(
    stage_epi = case_when(
      eGFR_epi >= 90 ~ "Stage 1",
      eGFR_epi >= 60 & eGFR_epi < 90 ~ "Stage 2",
      eGFR_epi < 60 ~ "Stage 3",
      TRUE ~ NA_character_
    )
  )

allDonor_CKD <- all_lab[,c(1,6,11,25:28)]

plot_data <- allDonor_CKD %>%
  filter(!is.na(stage_epi), !is.na(kclusterGroup)) %>%
  dplyr::group_by(kclusterGroup, stage_epi) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  dplyr::group_by(kclusterGroup) %>%
  dplyr::mutate(percentage = count / sum(count) * 100) %>%
  dplyr::ungroup() %>%
  mutate(stage_epi = factor(stage_epi))

ggplot(plot_data, aes(x = kclusterGroup, y = percentage, fill = stage_epi)) +
  geom_bar(stat = "identity", position = "stack", width = 0.80) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_brewer(palette = "YlOrRd") + 
  theme_bw() +
  labs(y = "Percentage (%)",fill = "CKD Stage",x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = "black", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = "black", linewidth = 1.7)) 

# Figure 1K-L ----------
# Liver Fibrosis Indices (FIB-4 and APRI)

calc_fib4 <- function(age, ast, alt, platelets) {
  score = ifelse(is.na(age) | is.na(ast) | is.na(alt) | is.na(platelets),
                 NA_real_,
                 (age * ast) / (platelets * sqrt(alt)))
  return(score)
}

all_lab <- all_lab %>%
  mutate(
    APRI = transplantr::apri(ast = SGOT.AST..u.L., plt = Platelets.THO.uL., ast_uln = 40),
    De_Ritis = (SGOT.AST..u.L. / SGPT..ALT..u.L.),
    FIB4 = calc_fib4(age = age_years, 
                     ast = SGOT.AST..u.L., 
                     alt = SGPT..ALT..u.L., 
                     platelets = Platelets.THO.uL.))

apri_summary <- all_lab %>%
  dplyr::group_by(kclusterGroup) %>%
  dplyr::summarise( 
    total = dplyr::n(), 
    over = sum(APRI > 2, na.rm = TRUE),
    freq = over / total,
    .groups = "drop")%>%
  mutate(kclusterGroup = factor(kclusterGroup, 
                                levels = c("ND", "MARD", "MOD", "SIDD", "SIRD")))

fib4_summary <- all_lab %>%
  dplyr::group_by(kclusterGroup) %>%
  dplyr::summarise(
    total = dplyr::n(), 
    over = sum(FIB4 > 3.25, na.rm = TRUE),
    freq = over / total,
    .groups = "drop")%>%
  mutate(kclusterGroup = factor(kclusterGroup, 
                                levels = c("ND", "MARD", "MOD", "SIDD", "SIRD")))

cols = c('#999898','#87a7d5','#4b5988','#f1c5a5','#c03f3f')

ggplot(apri_summary, aes(x = kclusterGroup, y = freq, fill = kclusterGroup)) +
  geom_bar(stat = "identity", width = 0.80) +
  labs(x = "",y = "Percentage of APRI > 2") +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 12),
        legend.position = "none",
        panel.border = element_rect(color = "black", linewidth = 1.7))

ggplot(fib4_summary, aes(x = kclusterGroup, y = freq, fill = kclusterGroup)) +
  geom_bar(stat = "identity", width = 0.80) +
  labs(x = "",y = "Percentage of FIB-4 > 3.25") +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 12),
        legend.position = "none",
        panel.border = element_rect(color = "black", linewidth = 1.7))