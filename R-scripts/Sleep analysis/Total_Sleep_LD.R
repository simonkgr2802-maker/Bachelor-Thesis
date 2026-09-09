# ============================================================
# TOTAL SLEEP BY PHASE — IMPORT AND DATA PREPARATION
# ============================================================


# ============================================================
# 1) PACKAGES
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# Models / diagnostics
library(glmmTMB)
library(lme4)
library(lmerTest)
library(DHARMa)
library(emmeans)
library(performance)


# ============================================================
# 2) SETTINGS
# ============================================================

FRAGRANCE_NAME <- "Ethanol"
CONCENTRATION_FILE <- "150uL"

DATA_FILE <- paste0(
  FRAGRANCE_NAME,
  "_",
  CONCENTRATION_FILE,
  "_total_sleep_by_phase_for_R.csv"
)

# Expected factor order
CONDITION_LEVELS <- c(
  "before",
  "after"
)

PHASE_LEVELS <- c(
  "Light",
  "Dark"
)

REPLICATE_LEVELS <- c(
  "R1",
  "R2"
)


# ============================================================
# 3) IMPORT DATA
# ============================================================

df_sleep_phase_raw <- read.csv(
  DATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

df_sleep_phase_raw <- tibble::as_tibble(
  df_sleep_phase_raw,
  .name_repair = "unique"
)


cat(
  "\n============================================\n",
  "RAW DATA IMPORTED\n",
  "============================================\n"
)

print(
  head(df_sleep_phase_raw)
)

cat(
  "\nDimensions:",
  nrow(df_sleep_phase_raw),
  "rows x",
  ncol(df_sleep_phase_raw),
  "columns\n"
)

cat(
  "\nColumn names:\n"
)

print(
  names(df_sleep_phase_raw)
)


# ============================================================
# 4) REMOVE ACCIDENTAL INDEX COLUMNS
# ============================================================
#
# Python exports with index = FALSE should not contain these,
# but this makes the import robust if an index column was
# accidentally written at some point.
# ============================================================

index_columns <- intersect(
  names(df_sleep_phase_raw),
  c(
    "...1",
    "X",
    "Unnamed: 0"
  )
)

if (
  length(index_columns) > 0
) {
  
  df_sleep_phase_raw <- df_sleep_phase_raw %>%
    select(
      -all_of(index_columns)
    )
}


# ============================================================
# 5) CHECK REQUIRED COLUMNS
# ============================================================

required_columns <- c(
  "ID",
  "Replicate",
  "condition",
  "Phase",
  "total_sleep_min"
)

missing_columns <- setdiff(
  required_columns,
  names(df_sleep_phase_raw)
)

if (
  length(missing_columns) > 0
) {
  
  stop(
    paste0(
      "Missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


cat(
  "\nAll required columns are present.\n"
)


# ============================================================
# 6) CLEAN VARIABLES
# ============================================================

df_sleep_phase <- df_sleep_phase_raw %>%
  
  mutate(
    
    # -----------------------------------------
    # ID
    # -----------------------------------------
    
    ID = trimws(
      as.character(ID)
    ),
    
    
    # -----------------------------------------
    # Replicate
    # -----------------------------------------
    
    Replicate = trimws(
      as.character(Replicate)
    ),
    
    
    # -----------------------------------------
    # Condition
    # -----------------------------------------
    
    condition = trimws(
      tolower(
        as.character(condition)
      )
    ),
    
    
    # -----------------------------------------
    # Phase
    # -----------------------------------------
    
    Phase = trimws(
      as.character(Phase)
    ),
    
    Phase = case_when(
      
      tolower(Phase) == "light" ~
        "Light",
      
      tolower(Phase) == "dark" ~
        "Dark",
      
      TRUE ~
        Phase
    ),
    
    
    # -----------------------------------------
    # Response variable
    # -----------------------------------------
    
    total_sleep_min = as.numeric(
      total_sleep_min
    )
  )


# ============================================================
# 7) FACTOR LEVELS
# ============================================================

df_sleep_phase <- df_sleep_phase %>%
  
  mutate(
    
    ID = factor(
      ID
    ),
    
    Replicate = factor(
      Replicate,
      levels = REPLICATE_LEVELS
    ),
    
    condition = factor(
      condition,
      levels = CONDITION_LEVELS
    ),
    
    Phase = factor(
      Phase,
      levels = PHASE_LEVELS
    )
  )


# ============================================================
# 8) REMOVE GENUINELY MISSING OBSERVATIONS
# ============================================================

cat(
  "\n============================================\n",
  "MISSING VALUES\n",
  "============================================\n"
)

print(
  colSums(
    is.na(
      df_sleep_phase[
        required_columns
      ]
    )
  )
)


n_before_na_removal <- nrow(
  df_sleep_phase
)


df_sleep_phase <- df_sleep_phase %>%
  
  filter(
    !is.na(ID),
    !is.na(Replicate),
    !is.na(condition),
    !is.na(Phase),
    !is.na(total_sleep_min)
  )


cat(
  "\nRows before NA removal:",
  n_before_na_removal,
  "\n"
)

cat(
  "Rows after NA removal:",
  nrow(df_sleep_phase),
  "\n"
)


# ============================================================
# 9) BIOLOGICAL RANGE CHECK
# ============================================================
#
# One Light or Dark phase lasts 12 h:
#
# 12 h x 60 min = 720 min
#
# Therefore:
#
# 0 <= total_sleep_min <= 720
#
# ============================================================

invalid_sleep_values <- df_sleep_phase %>%
  
  filter(
    total_sleep_min < 0 |
      total_sleep_min > 720
  )


cat(
  "\n============================================\n",
  "BIOLOGICAL RANGE CHECK\n",
  "============================================\n"
)

cat(
  "Number of observations outside 0–720 min:",
  nrow(invalid_sleep_values),
  "\n"
)


if (
  nrow(invalid_sleep_values) > 0
) {
  
  print(
    invalid_sleep_values
  )
}


# IMPORTANT:
# Do NOT automatically remove these values.
# If any exist, first investigate why.


# ============================================================
# 10) CHECK 5-MINUTE BIN STRUCTURE
# ============================================================
#
# Sleep was defined from 5-min bins.
# Therefore total_sleep_min should normally be divisible by 5.
# ============================================================

non_5min_values <- df_sleep_phase %>%
  
  filter(
    total_sleep_min %% 5 != 0
  )


cat(
  "\n============================================\n",
  "5-MINUTE BIN CHECK\n",
  "============================================\n"
)

cat(
  "Values not divisible by 5:",
  nrow(non_5min_values),
  "\n"
)


if (
  nrow(non_5min_values) > 0
) {
  
  print(
    non_5min_values
  )
}


# ============================================================
# 11) CHECK DUPLICATES
# ============================================================
#
# There should be exactly one observation for each:
#
# ID x condition x Phase
# ============================================================

duplicate_rows <- df_sleep_phase %>%
  
  count(
    ID,
    condition,
    Phase,
    name = "n"
  ) %>%
  
  filter(
    n > 1
  )


cat(
  "\n============================================\n",
  "DUPLICATE CHECK\n",
  "============================================\n"
)

cat(
  "Duplicate ID x condition x Phase combinations:",
  nrow(duplicate_rows),
  "\n"
)


if (
  nrow(duplicate_rows) > 0
) {
  
  print(
    duplicate_rows
  )
}


# ============================================================
# 12) CHECK OBSERVATIONS PER FLY
# ============================================================
#
# Normally each fly should have:
#
# before-Light
# before-Dark
# after-Light
# after-Dark
#
# = 4 observations per fly
# ============================================================

observations_per_fly <- df_sleep_phase %>%
  
  count(
    ID,
    name = "n_observations"
  )


cat(
  "\n============================================\n",
  "OBSERVATIONS PER FLY\n",
  "============================================\n"
)

print(
  table(
    observations_per_fly$n_observations
  )
)


incomplete_flies <- observations_per_fly %>%
  
  filter(
    n_observations != 4
  )


cat(
  "\nFlies without exactly 4 observations:",
  nrow(incomplete_flies),
  "\n"
)


if (
  nrow(incomplete_flies) > 0
) {
  
  print(
    incomplete_flies
  )
}


# ============================================================
# 13) CHECK EXPERIMENTAL DESIGN
# ============================================================

cat(
  "\n============================================\n",
  "EXPERIMENTAL DESIGN\n",
  "============================================\n"
)


cat(
  "\nCondition x Phase:\n"
)

print(
  table(
    df_sleep_phase$condition,
    df_sleep_phase$Phase
  )
)


cat(
  "\nReplicate x Condition x Phase:\n"
)

print(
  xtabs(
    ~ Replicate + condition + Phase,
    data = df_sleep_phase
  )
)


cat(
  "\nNumber of unique flies:",
  n_distinct(
    df_sleep_phase$ID
  ),
  "\n"
)


# ============================================================
# 14) OPTIONAL DERIVED VARIABLES
# ============================================================
#
# Because each phase contains:
#
# 12 hours x 12 five-minute bins/hour
# = 144 possible bins
#
# total_sleep_min / 5 = number of sleep bins
#
# These variables may be useful later when deciding
# which statistical distribution is most appropriate.
# ============================================================

df_sleep_phase <- df_sleep_phase %>%
  
  mutate(
    
    n_sleep_bins = total_sleep_min / 5,
    
    n_awake_bins = 144 - n_sleep_bins,
    
    sleep_fraction = total_sleep_min / 720
  )


# ============================================================
# 15) DESCRIPTIVE SUMMARY
# ============================================================

sleep_phase_summary <- df_sleep_phase %>%
  
  group_by(
    condition,
    Phase
  ) %>%
  
  summarise(
    
    n = n(),
    
    mean_sleep_min = mean(
      total_sleep_min,
      na.rm = TRUE
    ),
    
    sd_sleep_min = sd(
      total_sleep_min,
      na.rm = TRUE
    ),
    
    median_sleep_min = median(
      total_sleep_min,
      na.rm = TRUE
    ),
    
    min_sleep_min = min(
      total_sleep_min,
      na.rm = TRUE
    ),
    
    max_sleep_min = max(
      total_sleep_min,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


cat(
  "\n============================================\n",
  "DESCRIPTIVE SUMMARY\n",
  "============================================\n"
)

print(
  sleep_phase_summary
)


# ============================================================
# 16) FINAL DATAFRAME FOR MODELS
# ============================================================

df_total_sleep_phase <- df_sleep_phase %>%
  
  select(
    ID,
    Replicate,
    condition,
    Phase,
    total_sleep_min,
    n_sleep_bins,
    n_awake_bins,
    sleep_fraction,
    everything()
  )


cat(
  "\n============================================\n",
  "FINAL DATAFRAME FOR MODELLING\n",
  "============================================\n"
)


print(
  head(
    df_total_sleep_phase,
    8
  )
)


str(
  df_total_sleep_phase
)


# ============================================================
# 17) TOTAL SLEEP — LINEAR MIXED MODEL
# ============================================================

# Output folder
OUTPUT_DIR <- file.path(
  paste0("sleep_analysis_", FRAGRANCE_NAME),
  CONCENTRATION_FILE
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Final model — REML
# ------------------------------------------------------------

model_sleep_REML <- lmer(
  total_sleep_min ~ condition * Phase + Replicate + (1 | ID),
  data = df_total_sleep_phase,
  REML = TRUE
)

summary(model_sleep_REML)



# ============================================================
# 18) MODEL DIAGNOSTICS
# ============================================================

# Numerical checks
sleep_singular <- isSingular(
  model_sleep_REML,
  tol = 1e-4
)

sleep_heteroscedasticity <- performance::check_heteroscedasticity(
  model_sleep_REML
)

sleep_outliers <- performance::check_outliers(
  model_sleep_REML
)

cat("\nSingular fit:\n")
print(sleep_singular)

cat("\nHeteroscedasticity:\n")
print(sleep_heteroscedasticity)

cat("\nPotentially influential observations:\n")
print(sleep_outliers)


# ------------------------------------------------------------
# Diagnostic plots
# ------------------------------------------------------------

sleep_resid <- residuals(model_sleep_REML)
sleep_fitted <- fitted(model_sleep_REML)

sleep_random_intercepts <- ranef(
  model_sleep_REML
)$ID[, 1]


# Save all main diagnostics in one PDF
pdf(
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME, "_",
      CONCENTRATION_FILE,
      "_total_sleep_LMM_diagnostics.pdf"
    )
  ),
  width = 8,
  height = 8
)

par(mfrow = c(2, 2))


# 1) Residuals vs fitted
plot(
  sleep_fitted,
  sleep_resid,
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "Residuals vs fitted"
)

abline(h = 0, lty = 2)

lines(
  lowess(sleep_fitted, sleep_resid),
  lwd = 2
)


# 2) Q-Q residuals
qqnorm(
  sleep_resid,
  main = "Q-Q plot of residuals"
)

qqline(sleep_resid)


# 3) Q-Q random intercepts
qqnorm(
  sleep_random_intercepts,
  main = "Q-Q plot of random intercepts"
)

qqline(sleep_random_intercepts)


# 4) Residual spread by group
boxplot(
  sleep_resid ~ interaction(
    df_total_sleep_phase$condition,
    df_total_sleep_phase$Phase
  ),
  xlab = "",
  ylab = "Residuals",
  main = "Residuals by group",
  names = c(
    "Before\nLight",
    "After\nLight",
    "Before\nDark",
    "After\nDark"
  )
)

abline(h = 0, lty = 2)

dev.off()



# ============================================================
# 19) GLOBAL CONDITION × PHASE TEST — LRT
# ============================================================
#
# Fixed-effect comparisons require ML rather than REML.
# ============================================================

model_sleep_full_ML <- update(
  model_sleep_REML,
  REML = FALSE
)

model_sleep_additive_ML <- update(
  model_sleep_full_ML,
  . ~ . - condition:Phase
)


sleep_interaction_LRT <- anova(
  model_sleep_additive_ML,
  model_sleep_full_ML
)

cat(
  "\n============================================\n",
  "GLOBAL CONDITION × PHASE LRT\n",
  "============================================\n"
)

print(sleep_interaction_LRT)



# ============================================================
# 20) PLANNED SIMPLE LRTs
# ============================================================
#
# Scientific questions:
#
# Light: before vs after
# Dark:  before vs after
#
# The full four-cell model is equivalent to condition * Phase.
# ============================================================

df_sleep_LRT <- df_total_sleep_phase %>%
  mutate(
    
    PhaseCondition = factor(
      paste(Phase, condition, sep = "_")
    ),
    
    # H0 Light:
    # Light_before = Light_after
    LRT_Light = factor(
      case_when(
        Phase == "Light" ~ "Light_common",
        Phase == "Dark" & condition == "before" ~ "Dark_before",
        Phase == "Dark" & condition == "after"  ~ "Dark_after"
      )
    ),
    
    # H0 Dark:
    # Dark_before = Dark_after
    LRT_Dark = factor(
      case_when(
        Phase == "Dark" ~ "Dark_common",
        Phase == "Light" & condition == "before" ~ "Light_before",
        Phase == "Light" & condition == "after"  ~ "Light_after"
      )
    )
  )


# Full four-cell model
model_sleep_cells_ML <- lmer(
  total_sleep_min ~ PhaseCondition + Replicate + (1 | ID),
  data = df_sleep_LRT,
  REML = FALSE
)


# Check that both parameterizations are equivalent
sleep_logLik_check <- tibble(
  model = c(
    "condition * Phase",
    "four-cell model"
  ),
  logLik = c(
    as.numeric(logLik(model_sleep_full_ML)),
    as.numeric(logLik(model_sleep_cells_ML))
  )
)

cat("\nParameterization check:\n")
print(sleep_logLik_check)


# ------------------------------------------------------------
# Light null model
# ------------------------------------------------------------

model_sleep_Light_null_ML <- lmer(
  total_sleep_min ~ LRT_Light + Replicate + (1 | ID),
  data = df_sleep_LRT,
  REML = FALSE
)

sleep_LRT_Light <- anova(
  model_sleep_Light_null_ML,
  model_sleep_cells_ML
)


# ------------------------------------------------------------
# Dark null model
# ------------------------------------------------------------

model_sleep_Dark_null_ML <- lmer(
  total_sleep_min ~ LRT_Dark + Replicate + (1 | ID),
  data = df_sleep_LRT,
  REML = FALSE
)

sleep_LRT_Dark <- anova(
  model_sleep_Dark_null_ML,
  model_sleep_cells_ML
)


cat("\nLight — before vs after:\n")
print(sleep_LRT_Light)

cat("\nDark — before vs after:\n")
print(sleep_LRT_Dark)



# ============================================================
# 21) EXTRACT LRT P-VALUES + HOLM CORRECTION
# ============================================================

get_LRT_result <- function(x) {
  
  x <- as.data.frame(x)
  
  tibble(
    Chisq = x$Chisq[2],
    df = x$Df[2],
    p_LRT = x$`Pr(>Chisq)`[2]
  )
}


contrast_test_holm <- bind_rows(
  Light = get_LRT_result(sleep_LRT_Light),
  Dark  = get_LRT_result(sleep_LRT_Dark),
  .id = "Phase"
) %>%
  
  mutate(
    Phase = factor(
      Phase,
      levels = PHASE_LEVELS
    ),
    
    p.value = p.adjust(
      p_LRT,
      method = "holm"
    )
  )


cat(
  "\n============================================\n",
  "PLANNED LRTs + HOLM\n",
  "============================================\n"
)

print(contrast_test_holm)



# ============================================================
# 22) ESTIMATED MARGINAL MEANS + 95% CI
# ============================================================
#
# EMMs and CIs come from the final REML model.
# Primary p-values come from the LRTs above.
# ============================================================

emm_sleep <- emmeans(
  model_sleep_REML,
  ~ condition | Phase,
  lmer.df = "satterthwaite"
)


emm_sleep_df <- summary(
  emm_sleep,
  infer = c(TRUE, FALSE)
) %>%
  as.data.frame() %>%
  as_tibble()


cat(
  "\n============================================\n",
  "ESTIMATED MARGINAL MEANS\n",
  "============================================\n"
)

print(emm_sleep_df)



# ============================================================
# 23) AFTER − BEFORE EFFECT IN MINUTES
# ============================================================

sleep_difference_df <- contrast(
  emm_sleep,
  method = list(
    "after - before" = c(-1, 1)
  ),
  by = "Phase"
) %>%
  
  summary(
    infer = c(TRUE, FALSE),
    adjust = "none"
  ) %>%
  
  as.data.frame() %>%
  as_tibble()


# Combine effect estimate + CI with LRT p-values
sleep_results_final <- sleep_difference_df %>%
  
  select(
    Phase,
    contrast,
    estimate,
    SE,
    df,
    lower.CL,
    upper.CL
  ) %>%
  
  left_join(
    contrast_test_holm %>%
      select(
        Phase,
        Chisq,
        df_LRT = df,
        p_LRT,
        p_LRT_Holm = p.value
      ),
    by = "Phase"
  )


cat(
  "\n============================================\n",
  "FINAL TOTAL SLEEP RESULTS\n",
  "============================================\n"
)

print(sleep_results_final)



# ============================================================
# 24) RESULTS PLOT — TOTAL SLEEP
# ============================================================
#
# Same visual style as the no-LD Total Sleep plot.
# Total sleep remains on the original linear minute scale.
# ============================================================

condition_colors <- c(
  "before" = "#B79F00",
  "after"  = "#2F7F78"
)


# P-values for subtitle
p_light <- contrast_test_holm$p.value[
  contrast_test_holm$Phase == "Light"
]

p_dark <- contrast_test_holm$p.value[
  contrast_test_holm$Phase == "Dark"
]


format_p <- function(p) {
  ifelse(
    p < 0.001,
    "p < 0.001",
    paste0(
      "p = ",
      formatC(p, format = "f", digits = 3)
    )
  )
}


p_light_text <- format_p(p_light)
p_dark_text  <- format_p(p_dark)


# EMM dataframe for plotting
emm_plot_df <- emm_sleep_df %>%
  dplyr::select(
    condition,
    Phase,
    emmean,
    lower.CL,
    upper.CL
  )


nudge_amt <- 0.25


total_sleep_plot <- ggplot(
  df_total_sleep_phase,
  aes(
    x = condition,
    y = total_sleep_min
  )
) +
  
  # Individual flies
  geom_jitter(
    aes(
      color = condition,
      shape = "Individual flies"
    ),
    width = 0.08,
    height = 0,
    alpha = 0.60,
    size = 1.6
  ) +
  
  # White halo behind CI
  geom_errorbar(
    data = emm_plot_df,
    aes(
      x = condition,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    inherit.aes = FALSE,
    position = position_nudge(x = nudge_amt),
    width = 0.14,
    linewidth = 2.4,
    color = "white"
  ) +
  
  # 95% CI
  geom_errorbar(
    data = emm_plot_df,
    aes(
      x = condition,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    inherit.aes = FALSE,
    position = position_nudge(x = nudge_amt),
    width = 0.14,
    linewidth = 1,
    color = "black"
  ) +
  
  # Model-estimated means
  geom_point(
    data = emm_plot_df,
    aes(
      x = condition,
      y = emmean,
      fill = condition,
      shape = "Model-estimated mean ± 95% CI"
    ),
    inherit.aes = FALSE,
    position = position_nudge(x = nudge_amt),
    size = 4,
    stroke = 1.1,
    color = "black"
  ) +
  
  facet_wrap(
    ~ Phase
  ) +
  
  scale_color_manual(
    values = condition_colors,
    guide = "none"
  ) +
  
  scale_fill_manual(
    values = condition_colors,
    guide = "none"
  ) +
  
  scale_shape_manual(
    name = NULL,
    values = c(
      "Individual flies" = 16,
      "Model-estimated mean ± 95% CI" = 21
    )
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 720, by = 120),
    expand = expansion(
      mult = c(0.02, 0.05)
    )
  ) +
  
  labs(
    x = NULL,
    y = "Total sleep time (min)",
    
    title = paste0(
      FRAGRANCE_NAME,
      " ",
      CONCENTRATION_FILE,
      ": total sleep time"
    ),
    
    subtitle = paste0(
      "Before vs. after:   Light ",
      p_light_text,
      "   |   Dark ",
      p_dark_text
    )
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.7
    ),
    
    strip.text = element_text(
      size = 12
    ),
    
    plot.title = element_text(
      size = 14
    ),
    
    plot.subtitle = element_text(
      size = 10.5,
      margin = margin(b = 6)
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    panel.spacing = unit(
      0.3,
      "cm"
    ),
    
    legend.position = "bottom",
    
    legend.text = element_text(
      size = 9.5
    )
  ) +
  
  guides(
    shape = guide_legend(
      override.aes = list(
        size = c(2, 4),
        alpha = c(0.75, 1),
        color = c("grey30", "black"),
        fill = c(NA, "white")
      )
    )
  )


print(total_sleep_plot)



# ============================================================
# 25) SAVE RESULTS
# ============================================================

file_prefix <- paste0(
  FRAGRANCE_NAME,
  "_",
  CONCENTRATION_FILE,
  "_total_sleep"
)


# ------------------------------------------------------------
# Result plot
# ------------------------------------------------------------

ggsave(
  file.path(
    OUTPUT_DIR,
    paste0(file_prefix, "_results.png")
  ),
  total_sleep_plot,
  width = 8,
  height = 5.3,
  units = "in",
  dpi = 300,
  bg = "white"
)


ggsave(
  file.path(
    OUTPUT_DIR,
    paste0(file_prefix, "_results.pdf")
  ),
  total_sleep_plot,
  width = 8,
  height = 5.3,
  units = "in",
  bg = "white"
)


# ------------------------------------------------------------
# Statistical tables
# ------------------------------------------------------------

write.csv(
  as.data.frame(sleep_interaction_LRT),
  file.path(
    OUTPUT_DIR,
    paste0(file_prefix, "_interaction_LRT.csv")
  ),
  row.names = TRUE
)


write.csv(
  emm_sleep_df,
  file.path(
    OUTPUT_DIR,
    paste0(file_prefix, "_EMMeans.csv")
  ),
  row.names = FALSE
)


write.csv(
  sleep_results_final,
  file.path(
    OUTPUT_DIR,
    paste0(file_prefix, "_planned_LRT_Holm.csv")
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Save fitted models
# ------------------------------------------------------------

saveRDS(
  model_sleep_REML,
  file.path(
    OUTPUT_DIR,
    paste0(file_prefix, "_LMM_REML.rds")
  )
)


saveRDS(
  model_sleep_full_ML,
  file.path(
    OUTPUT_DIR,
    paste0(file_prefix, "_LMM_ML.rds")
  )
)



# ============================================================
# 26) SAVE COMPLETE STATISTICAL OUTPUT
# ============================================================

statistics_file <- file.path(
  OUTPUT_DIR,
  paste0(
    file_prefix,
    "_statistics.txt"
  )
)


capture.output(
  
  {
    
    cat(
      "TOTAL SLEEP — LINEAR MIXED MODEL\n\n"
    )
    
    cat(
      "Model:\n",
      "total_sleep_min ~ condition * Phase + ",
      "Replicate + (1 | ID)\n\n"
    )
    
    
    cat(
      "====================\n",
      "REML MODEL\n",
      "====================\n"
    )
    
    print(
      summary(model_sleep_REML)
    )
    
    
    cat(
      "\n\n====================\n",
      "DIAGNOSTICS\n",
      "====================\n"
    )
    
    cat(
      "\nSingular fit:",
      sleep_singular,
      "\n\n"
    )
    
    print(
      sleep_heteroscedasticity
    )
    
    cat("\n")
    
    print(
      sleep_outliers
    )
    
    
    cat(
      "\n\n====================\n",
      "GLOBAL INTERACTION LRT\n",
      "====================\n"
    )
    
    print(
      sleep_interaction_LRT
    )
    
    
    cat(
      "\n\n====================\n",
      "PARAMETERIZATION CHECK\n",
      "====================\n"
    )
    
    print(
      sleep_logLik_check
    )
    
    
    cat(
      "\n\n====================\n",
      "LIGHT LRT\n",
      "====================\n"
    )
    
    print(
      sleep_LRT_Light
    )
    
    
    cat(
      "\n\n====================\n",
      "DARK LRT\n",
      "====================\n"
    )
    
    print(
      sleep_LRT_Dark
    )
    
    
    cat(
      "\n\n====================\n",
      "PLANNED TESTS + HOLM\n",
      "====================\n"
    )
    
    print(
      contrast_test_holm
    )
    
    
    cat(
      "\n\n====================\n",
      "EMMeans + 95% CI\n",
      "====================\n"
    )
    
    print(
      emm_sleep_df
    )
    
    
    cat(
      "\n\n====================\n",
      "FINAL EFFECT TABLE\n",
      "====================\n"
    )
    
    print(
      sleep_results_final
    )
    
    
    cat(
      "\n\n====================\n",
      "SESSION INFO\n",
      "====================\n"
    )
    
    print(
      sessionInfo()
    )
  },
  
  file = statistics_file
)



# ============================================================
# 27) DONE
# ============================================================

cat(
  "\n============================================\n",
  "TOTAL SLEEP ANALYSIS COMPLETE\n",
  "============================================\n",
  "\nResults saved in:\n",
  normalizePath(OUTPUT_DIR),
  "\n"
)