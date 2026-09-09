# ============================================================
# TOTAL SLEEP WITHOUT LIGHT/DARK PHASE — IMPORT AND DATA PREPARATION
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

# New output folder for the analysis without Light/Dark phase
OUTPUT_DIR <- paste0(
  FRAGRANCE_NAME,
  "_Total_sleep_no_LD"
)

if (
  !dir.exists(OUTPUT_DIR)
) {
  dir.create(
    OUTPUT_DIR,
    recursive = TRUE
  )
}

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
    dplyr::select(
      -all_of(index_columns)
    )
}


# ============================================================
# 5) CHECK REQUIRED COLUMNS
# ============================================================
#
# The phase-specific export is used as input and Light + Dark
# are summed later to obtain total sleep across the complete
# 24-h baseline and exposure periods.
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
# 9) CHECK PHASE-SPECIFIC INPUT BEFORE 24-H POOLING
# ============================================================
#
# The input file contains one observation for each:
#
# ID x condition x Phase
#
# Each fly should therefore have:
#
# before-Light
# before-Dark
# after-Light
# after-Dark
#
# = 4 observations per fly
#
# This check is performed before Light and Dark are summed.
# ============================================================

duplicate_phase_rows <- df_sleep_phase %>%
  
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
  "PHASE-SPECIFIC INPUT CHECK\n",
  "============================================\n"
)

cat(
  "Duplicate ID x condition x Phase combinations:",
  nrow(duplicate_phase_rows),
  "\n"
)


if (
  nrow(duplicate_phase_rows) > 0
) {
  
  print(
    duplicate_phase_rows
  )
}


observations_per_fly_phase <- df_sleep_phase %>%
  
  count(
    ID,
    name = "n_observations"
  )


print(
  table(
    observations_per_fly_phase$n_observations
  )
)


incomplete_flies_phase <- observations_per_fly_phase %>%
  
  filter(
    n_observations != 4
  )


cat(
  "\nFlies without exactly 4 phase-specific observations:",
  nrow(incomplete_flies_phase),
  "\n"
)


if (
  nrow(incomplete_flies_phase) > 0
) {
  
  print(
    incomplete_flies_phase
  )
}


# ============================================================
# 10) POOL LIGHT + DARK TO TOTAL 24-H SLEEP
# ============================================================
#
# For the analysis without LD phase, total sleep is calculated
# across the complete 24-h baseline or exposure period:
#
# total_sleep_min =
#   Light total sleep + Dark total sleep
#
# This gives one observation per:
#
# ID x Replicate x condition
# ============================================================

df_total_sleep <- df_sleep_phase %>%
  
  group_by(
    ID,
    Replicate,
    condition
  ) %>%
  
  summarise(
    
    total_sleep_min = sum(
      total_sleep_min,
      na.rm = FALSE
    ),
    
    n_phases = n_distinct(
      Phase
    ),
    
    .groups = "drop"
  )


# Check that both Light and Dark contributed to every 24-h total
incomplete_24h_rows <- df_total_sleep %>%
  
  filter(
    n_phases != 2
  )


cat(
  "\n============================================\n",
  "24-H POOLING CHECK\n",
  "============================================\n"
)

cat(
  "24-h observations without both Light and Dark phases:",
  nrow(incomplete_24h_rows),
  "\n"
)


if (
  nrow(incomplete_24h_rows) > 0
) {
  
  print(
    incomplete_24h_rows
  )
}


# Keep only variables required for the pooled analysis
df_total_sleep <- df_total_sleep %>%
  
  dplyr::select(
    ID,
    Replicate,
    condition,
    total_sleep_min
  )


# ============================================================
# 11) BIOLOGICAL RANGE CHECK
# ============================================================
#
# One complete baseline or exposure period lasts 24 h:
#
# 24 h x 60 min = 1440 min
#
# Therefore:
#
# 0 <= total_sleep_min <= 1440
#
# ============================================================

invalid_sleep_values <- df_total_sleep %>%
  
  filter(
    total_sleep_min < 0 |
      total_sleep_min > 1440
  )


cat(
  "\n============================================\n",
  "BIOLOGICAL RANGE CHECK\n",
  "============================================\n"
)

cat(
  "Number of observations outside 0–1440 min:",
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
# 12) CHECK 5-MINUTE BIN STRUCTURE
# ============================================================
#
# Sleep was defined from 5-min bins.
# Therefore total_sleep_min should normally be divisible by 5.
# ============================================================

non_5min_values <- df_total_sleep %>%
  
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
# 13) CHECK DUPLICATES
# ============================================================
#
# After pooling Light + Dark, there should be exactly one
# observation for each:
#
# ID x condition
# ============================================================

duplicate_rows <- df_total_sleep %>%
  
  count(
    ID,
    condition,
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
  "Duplicate ID x condition combinations:",
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
# 14) CHECK OBSERVATIONS PER FLY
# ============================================================
#
# Normally each fly should now have:
#
# before
# after
#
# = 2 observations per fly
# ============================================================

observations_per_fly <- df_total_sleep %>%
  
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
    n_observations != 2
  )


cat(
  "\nFlies without exactly 2 observations:",
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
# 15) CHECK EXPERIMENTAL DESIGN
# ============================================================

cat(
  "\n============================================\n",
  "EXPERIMENTAL DESIGN\n",
  "============================================\n"
)


cat(
  "\nCondition:\n"
)

print(
  table(
    df_total_sleep$condition
  )
)


cat(
  "\nReplicate x Condition:\n"
)

print(
  xtabs(
    ~ Replicate + condition,
    data = df_total_sleep
  )
)


cat(
  "\nNumber of unique flies:",
  n_distinct(
    df_total_sleep$ID
  ),
  "\n"
)


# ============================================================
# 16) OPTIONAL DERIVED VARIABLES
# ============================================================
#
# Each complete 24-h period contains:
#
# 24 hours x 12 five-minute bins/hour
# = 288 possible bins
#
# total_sleep_min / 5 = number of sleep bins
#
# These variables may be useful later when deciding
# which statistical distribution is most appropriate.
# ============================================================

df_total_sleep <- df_total_sleep %>%
  
  mutate(
    
    n_sleep_bins = total_sleep_min / 5,
    
    n_awake_bins = 288 - n_sleep_bins,
    
    sleep_fraction = total_sleep_min / 1440
  )


# ============================================================
# 17) DESCRIPTIVE SUMMARY
# ============================================================

sleep_summary <- df_total_sleep %>%
  
  group_by(
    condition
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
  "DESCRIPTIVE SUMMARY — TOTAL 24-H SLEEP\n",
  "============================================\n"
)

print(
  sleep_summary
)


# ============================================================
# 18) FINAL DATAFRAME FOR MODELS
# ============================================================

df_total_sleep_no_LD <- df_total_sleep %>%
  
  dplyr::select(
    ID,
    Replicate,
    condition,
    total_sleep_min,
    n_sleep_bins,
    n_awake_bins,
    sleep_fraction,
    everything()
  )


cat(
  "\n============================================\n",
  "FINAL DATAFRAME FOR MODELLING — NO LD PHASE\n",
  "============================================\n"
)


print(
  head(
    df_total_sleep_no_LD,
    8
  )
)


str(
  df_total_sleep_no_LD
)

# ============================================================
# MODELS START BELOW THIS POINT
# ============================================================

# ============================================================
# 19) TOTAL SLEEP — LINEAR MIXED MODEL WITHOUT LD PHASE
# ============================================================

# Output folder
OUTPUT_DIR <- file.path(
  paste0(
    FRAGRANCE_NAME,
    "_Total_sleep_no_LD"
  ),
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
  total_sleep_min ~ condition + Replicate + (1 | ID),
  data = df_total_sleep_no_LD,
  REML = TRUE
)

summary(model_sleep_REML)



# ============================================================
# 20) MODEL DIAGNOSTICS
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

sleep_resid <- residuals(
  model_sleep_REML
)

sleep_fitted <- fitted(
  model_sleep_REML
)

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
      "_total_sleep_no_LD_LMM_diagnostics.pdf"
    )
  ),
  width = 8,
  height = 8
)

par(
  mfrow = c(2, 2)
)


# 1) Residuals vs fitted
plot(
  sleep_fitted,
  sleep_resid,
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "Residuals vs fitted"
)

abline(
  h = 0,
  lty = 2
)

lines(
  lowess(
    sleep_fitted,
    sleep_resid
  ),
  lwd = 2
)


# 2) Q-Q residuals
qqnorm(
  sleep_resid,
  main = "Q-Q plot of residuals"
)

qqline(
  sleep_resid
)


# 3) Q-Q random intercepts
qqnorm(
  sleep_random_intercepts,
  main = "Q-Q plot of random intercepts"
)

qqline(
  sleep_random_intercepts
)


# 4) Residual spread by condition
boxplot(
  sleep_resid ~ df_total_sleep_no_LD$condition,
  xlab = "",
  ylab = "Residuals",
  main = "Residuals by condition",
  names = c(
    "Before",
    "After"
  )
)

abline(
  h = 0,
  lty = 2
)

dev.off()



# ============================================================
# 21) GLOBAL BEFORE–AFTER TEST — LRT
# ============================================================
#
# Primary question:
# Does total sleep across the complete 24-h period differ
# between baseline and exposure?
#
# Fixed-effect comparisons require ML rather than REML.
# ============================================================

model_sleep_full_ML <- update(
  model_sleep_REML,
  REML = FALSE
)

model_sleep_no_condition_ML <- update(
  model_sleep_full_ML,
  . ~ . - condition
)


sleep_condition_LRT <- anova(
  model_sleep_no_condition_ML,
  model_sleep_full_ML
)


cat(
  "\n============================================\n",
  "GLOBAL 24-H BEFORE–AFTER LRT\n",
  "============================================\n"
)

print(
  sleep_condition_LRT
)

# ============================================================
# 22) EXTRACT GLOBAL 24-H LRT RESULT
# ============================================================

get_LRT_result <- function(x) {
  
  x <- as.data.frame(x)
  
  tibble(
    Chisq = x$Chisq[2],
    df_LRT = x$Df[2],
    p_LRT = x$`Pr(>Chisq)`[2]
  )
}


sleep_condition_test <- get_LRT_result(
  sleep_condition_LRT
)


cat(
  "\n============================================\n",
  "GLOBAL 24-H BEFORE–AFTER LRT RESULT\n",
  "============================================\n"
)

print(
  sleep_condition_test
)



# ============================================================
# 23) ESTIMATED MARGINAL MEANS + 95% CI
# ============================================================
#
# EMMs and CIs come from the final REML model.
# Since Light/Dark phase is not included in this analysis,
# EMMs are estimated only for condition (before vs after).
#
# The primary p-value comes from the global 24-h LRT above.
# ============================================================

emm_sleep <- emmeans(
  model_sleep_REML,
  ~ condition,
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
  "ESTIMATED MARGINAL MEANS — TOTAL 24-H SLEEP\n",
  "============================================\n"
)

print(
  emm_sleep_df
)



# ============================================================
# 24) AFTER − BEFORE EFFECT IN MINUTES
# ============================================================
#
# Positive estimate:
# more total sleep during exposure
#
# Negative estimate:
# less total sleep during exposure
#
# Only one predefined before–after comparison is performed,
# therefore no multiplicity adjustment is required.
# ============================================================

sleep_difference_df <- contrast(
  emm_sleep,
  method = list(
    "after - before" = c(-1, 1)
  )
) %>%
  
  summary(
    infer = c(TRUE, FALSE),
    adjust = "none"
  ) %>%
  
  as.data.frame() %>%
  as_tibble()


# Combine effect estimate + CI with primary LRT result
sleep_results_final <- sleep_difference_df %>%
  
  dplyr::select(
    contrast,
    estimate,
    SE,
    df,
    lower.CL,
    upper.CL
  ) %>%
  
  mutate(
    Chisq = sleep_condition_test$Chisq,
    df_LRT = sleep_condition_test$df_LRT,
    p_LRT = sleep_condition_test$p_LRT
  )


cat(
  "\n============================================\n",
  "FINAL TOTAL SLEEP RESULTS — NO LD PHASE\n",
  "============================================\n"
)

print(
  sleep_results_final
)



# ============================================================
# 25) RESULTS PLOT — TOTAL 24-H SLEEP
# ============================================================
#
# Same visual style as the previous Total Sleep plot.
# Total sleep remains on the original linear minute scale.
#
# Colors are slightly darker to improve visibility of
# individual observations on a white background.
# ============================================================

condition_colors <- c(
  "before" = "#B79F00",
  "after"  = "#2F7F78"
)


# ------------------------------------------------------------
# P-value for subtitle
# ------------------------------------------------------------

format_p <- function(p) {
  
  ifelse(
    p < 0.001,
    "p < 0.001",
    paste0(
      "p = ",
      formatC(
        p,
        format = "f",
        digits = 3
      )
    )
  )
}


p_condition <- sleep_condition_test$p_LRT
p_condition_text <- format_p(
  p_condition
)


# ------------------------------------------------------------
# EMM dataframe for plotting
# ------------------------------------------------------------

emm_plot_df <- emm_sleep_df %>%
  
  dplyr::select(
    condition,
    emmean,
    lower.CL,
    upper.CL
  )


nudge_amt <- 0.25


total_sleep_plot <- ggplot(
  df_total_sleep_no_LD,
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
    position = position_nudge(
      x = nudge_amt
    ),
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
    position = position_nudge(
      x = nudge_amt
    ),
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
    position = position_nudge(
      x = nudge_amt
    ),
    size = 4,
    stroke = 1.1,
    color = "black"
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
    breaks = seq(
      0,
      1440,
      by = 240
    ),
    limits = c(
      0,
      1440
    ),
    expand = expansion(
      mult = c(
        0.02,
        0.05
      )
    )
  ) +
  
  labs(
    x = NULL,
    y = "Total sleep time (min)",
    
    title = paste0(
      FRAGRANCE_NAME,
      " ",
      CONCENTRATION_FILE,
      ": total 24-h sleep time"
    ),
    
    subtitle = paste0(
      "Before vs. after: ",
      p_condition_text
    )
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    plot.title = element_text(
      size = 14
    ),
    
    plot.subtitle = element_text(
      size = 10.5,
      margin = margin(
        b = 6
      )
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    legend.position = "bottom",
    
    legend.text = element_text(
      size = 9.5
    )
  ) +
  
  guides(
    shape = guide_legend(
      override.aes = list(
        size = c(
          2,
          4
        ),
        alpha = c(
          0.75,
          1
        ),
        color = c(
          "grey30",
          "black"
        ),
        fill = c(
          NA,
          "white"
        )
      )
    )
  )


print(
  total_sleep_plot
)



# ============================================================
# 26) SAVE RESULTS
# ============================================================

file_prefix <- paste0(
  FRAGRANCE_NAME,
  "_",
  CONCENTRATION_FILE,
  "_total_sleep_no_LD"
)


# ------------------------------------------------------------
# Result plot
# ------------------------------------------------------------

ggsave(
  file.path(
    OUTPUT_DIR,
    paste0(
      file_prefix,
      "_results.png"
    )
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
    paste0(
      file_prefix,
      "_results.pdf"
    )
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
  as.data.frame(
    sleep_condition_LRT
  ),
  file.path(
    OUTPUT_DIR,
    paste0(
      file_prefix,
      "_condition_LRT.csv"
    )
  ),
  row.names = TRUE
)


write.csv(
  emm_sleep_df,
  file.path(
    OUTPUT_DIR,
    paste0(
      file_prefix,
      "_EMMeans.csv"
    )
  ),
  row.names = FALSE
)


write.csv(
  sleep_results_final,
  file.path(
    OUTPUT_DIR,
    paste0(
      file_prefix,
      "_effect_and_LRT.csv"
    )
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
    paste0(
      file_prefix,
      "_LMM_REML.rds"
    )
  )
)


saveRDS(
  model_sleep_full_ML,
  file.path(
    OUTPUT_DIR,
    paste0(
      file_prefix,
      "_LMM_ML.rds"
    )
  )
)


saveRDS(
  model_sleep_no_condition_ML,
  file.path(
    OUTPUT_DIR,
    paste0(
      file_prefix,
      "_LMM_null_ML.rds"
    )
  )
)



# ============================================================
# 27) SAVE COMPLETE STATISTICAL OUTPUT
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
      "TOTAL 24-H SLEEP — LINEAR MIXED MODEL WITHOUT LD PHASE\n\n"
    )
    
    cat(
      "Model:\n",
      "total_sleep_min ~ condition + ",
      "Replicate + (1 | ID)\n\n"
    )
    
    
    cat(
      "====================\n",
      "REML MODEL\n",
      "====================\n"
    )
    
    print(
      summary(
        model_sleep_REML
      )
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
    
    cat(
      "\n"
    )
    
    print(
      sleep_outliers
    )
    
    
    cat(
      "\n\n====================\n",
      "GLOBAL 24-H BEFORE–AFTER LRT\n",
      "====================\n"
    )
    
    print(
      sleep_condition_LRT
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
      "AFTER − BEFORE EFFECT + LRT\n",
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
# 28) DONE
# ============================================================

cat(
  "\n============================================\n",
  "TOTAL 24-H SLEEP ANALYSIS COMPLETE\n",
  "============================================\n",
  "\nResults saved in:\n",
  normalizePath(
    OUTPUT_DIR
  ),
  "\n"
)


