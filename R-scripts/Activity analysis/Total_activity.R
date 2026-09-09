# ============================================================
# DAM2 TOTAL ACTIVITY — 24-H ANALYSIS WITHOUT LIGHT/DARK PHASE
# ============================================================
#
# MAIN MODEL:
#   Negative-binomial GLMM (nbinom2)
#
# INPUT:
#   one row = one fly × one condition × one phase
#
# The Light and Dark observations are summed first so that the
# final analysis contains one total 24-h activity value for
# each fly in each condition:
#
#   before = complete 24-h baseline period
#   after  = complete 24-h exposure period
#
# PRIMARY QUESTION:
#   Does total locomotor activity across the complete 24-h
#   period differ between baseline and exposure?
#
# before is the reference condition.
# R1 is the reference replicate.
#
# For a new fragrance/concentration, normally only edit SECTION 2.
# ============================================================


# ============================================================
# 1) PACKAGES
# ============================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(DHARMa)
library(emmeans)


# ============================================================
# 2) SETTINGS — EDIT HERE
# ============================================================

FRAGRANCE_NAME <- "Ethanol"
CONCENTRATION_FILE <- "150uL"

DATA_FILE <- paste0(
  FRAGRANCE_NAME,
  "_",
  CONCENTRATION_FILE,
  "_activity_by_phase_for_R.csv"
)

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

# Number of DHARMa simulations
N_SIMULATIONS <- 1000


# Output folder — separate from the LD analysis
OUTPUT_DIR <- file.path(
  paste0(
    FRAGRANCE_NAME,
    "_Total_activity_no_LD"
  ),
  CONCENTRATION_FILE
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3) DATA IMPORT
# ============================================================

df_total_activity_phase <- read.csv(
  DATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

df_total_activity_phase <- tibble::as_tibble(
  df_total_activity_phase,
  .name_repair = "unique"
)


# Remove accidental CSV index columns if present
index_columns <- intersect(
  names(df_total_activity_phase),
  c(
    "...1",
    "X",
    "Unnamed: 0"
  )
)

if (
  length(index_columns) > 0
) {

  df_total_activity_phase <- df_total_activity_phase %>%
    dplyr::select(
      -all_of(index_columns)
    )
}


cat(
  "\nImported file:",
  DATA_FILE,
  "\nRows:",
  nrow(df_total_activity_phase),
  "| Columns:",
  ncol(df_total_activity_phase),
  "\nOutput folder:",
  normalizePath(OUTPUT_DIR),
  "\n\n"
)

print(
  names(df_total_activity_phase)
)

print(
  head(df_total_activity_phase)
)


# ============================================================
# 4) CHECK REQUIRED COLUMNS
# ============================================================
#
# Phase is required only because the existing export contains
# separate Light and Dark totals. The two phases are pooled
# before modelling.
# ============================================================

required_columns <- c(
  "ID",
  "Replicate",
  "condition",
  "Phase",
  "total_activity"
)

missing_columns <- setdiff(
  required_columns,
  names(df_total_activity_phase)
)

if (
  length(missing_columns) > 0
) {

  stop(
    paste(
      "Missing required columns:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 5) PREPARE VARIABLES
# ============================================================

df_total_activity_phase <- df_total_activity_phase %>%

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
    ),

    total_activity = as.integer(
      total_activity
    )
  )


# ============================================================
# 6) BASIC CHECKS OF PHASE-SPECIFIC INPUT
# ============================================================

cat(
  "\n================ PHASE-SPECIFIC INPUT CHECKS ================\n"
)

cat(
  "\nNumber of flies:\n"
)

print(
  n_distinct(
    df_total_activity_phase$ID
  )
)

cat(
  "\nCondition × Phase counts:\n"
)

print(
  table(
    df_total_activity_phase$condition,
    df_total_activity_phase$Phase
  )
)

cat(
  "\nReplicate × Condition × Phase counts:\n"
)

print(
  xtabs(
    ~ Replicate + condition + Phase,
    data = df_total_activity_phase
  )
)

cat(
  "\nObservations per fly:\n"
)

print(
  table(
    table(
      df_total_activity_phase$ID
    )
  )
)

cat(
  "\nNegative total activity values:\n"
)

print(
  any(
    df_total_activity_phase$total_activity < 0,
    na.rm = TRUE
  )
)

cat(
  "\nMissing total activity values:\n"
)

print(
  anyNA(
    df_total_activity_phase$total_activity
  )
)

cat(
  "\nInteger counts only:\n"
)

print(
  all(
    df_total_activity_phase$total_activity ==
      floor(
        df_total_activity_phase$total_activity
      ),
    na.rm = TRUE
  )
)


# ------------------------------------------------------------
# Duplicate check before pooling
# ------------------------------------------------------------

duplicate_phase_rows <- df_total_activity_phase %>%

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
  "\nDuplicated ID × condition × Phase combinations:",
  nrow(duplicate_phase_rows),
  "\n"
)

if (
  nrow(duplicate_phase_rows) > 0
) {

  print(
    duplicate_phase_rows
  )

  warning(
    "Duplicated ID × condition × Phase observations were found."
  )
}


# ------------------------------------------------------------
# Every fly should normally have four phase-specific rows:
#
# before-Light
# before-Dark
# after-Light
# after-Dark
# ------------------------------------------------------------

observations_per_fly_phase <- df_total_activity_phase %>%

  count(
    ID,
    name = "n_observations"
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
# 7) POOL LIGHT + DARK TO TOTAL 24-H ACTIVITY
# ============================================================
#
# total_activity =
#   Light activity + Dark activity
#
# This produces one observation for each:
#
# ID × Replicate × condition
# ============================================================

df_total_activity <- df_total_activity_phase %>%

  group_by(
    ID,
    Replicate,
    condition
  ) %>%

  summarise(

    total_activity = sum(
      total_activity,
      na.rm = FALSE
    ),

    n_phases = n_distinct(
      Phase
    ),

    .groups = "drop"
  )


# Check that both Light and Dark contributed to every 24-h total
incomplete_24h_rows <- df_total_activity %>%

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


# Keep only variables required for the no-LD analysis
df_total_activity <- df_total_activity %>%

  dplyr::select(
    ID,
    Replicate,
    condition,
    total_activity
  ) %>%

  mutate(

    total_activity = as.integer(
      total_activity
    ),

    # Optional sensitivity-analysis response
    log_total_activity = log1p(
      total_activity
    )
  )


# ============================================================
# 8) CHECK FINAL 24-H DATASET
# ============================================================

cat(
  "\n================ 24-H DATA CHECKS ================\n"
)

cat(
  "\nNumber of flies:\n"
)

print(
  n_distinct(
    df_total_activity$ID
  )
)

cat(
  "\nCondition counts:\n"
)

print(
  table(
    df_total_activity$condition
  )
)

cat(
  "\nReplicate × Condition counts:\n"
)

print(
  xtabs(
    ~ Replicate + condition,
    data = df_total_activity
  )
)

cat(
  "\nTotal activity summary:\n"
)

print(
  summary(
    df_total_activity$total_activity
  )
)

cat(
  "\nNegative total activity values:\n"
)

print(
  any(
    df_total_activity$total_activity < 0,
    na.rm = TRUE
  )
)

cat(
  "\nMissing total activity values:\n"
)

print(
  anyNA(
    df_total_activity$total_activity
  )
)

cat(
  "\nInteger counts only:\n"
)

print(
  all(
    df_total_activity$total_activity ==
      floor(
        df_total_activity$total_activity
      ),
    na.rm = TRUE
  )
)


# ------------------------------------------------------------
# Duplicate check after pooling
# ------------------------------------------------------------

duplicate_rows <- df_total_activity %>%

  count(
    ID,
    condition,
    name = "n"
  ) %>%

  filter(
    n > 1
  )


cat(
  "\nDuplicated ID × condition combinations:",
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


# ------------------------------------------------------------
# Each fly should now have exactly two observations:
#
# before
# after
# ------------------------------------------------------------

observations_per_fly <- df_total_activity %>%

  count(
    ID,
    name = "n_observations"
  )


cat(
  "\nObservations per fly:\n"
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
# 9) DESCRIPTIVE STATISTICS
# ============================================================

total_activity_summary <- df_total_activity %>%

  group_by(
    Replicate,
    condition
  ) %>%

  summarise(

    n = n(),

    n_flies = n_distinct(
      ID
    ),

    mean_total_activity = mean(
      total_activity,
      na.rm = TRUE
    ),

    sd_total_activity = sd(
      total_activity,
      na.rm = TRUE
    ),

    median_total_activity = median(
      total_activity,
      na.rm = TRUE
    ),

    minimum = min(
      total_activity,
      na.rm = TRUE
    ),

    maximum = max(
      total_activity,
      na.rm = TRUE
    ),

    se_total_activity =
      sd_total_activity / sqrt(
        n_flies
      ),

    .groups = "drop"
  )


cat(
  "\n============================================\n",
  "DESCRIPTIVE STATISTICS — TOTAL 24-H ACTIVITY\n",
  "============================================\n"
)

print(
  total_activity_summary
)


write.csv(
  total_activity_summary,
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_total_activity_no_LD_descriptive_statistics.csv"
    )
  ),
  row.names = FALSE
)


# ============================================================
# 10) PRIMARY MODEL — NEGATIVE-BINOMIAL GLMM
# ============================================================
#
# Primary 24-h model:
#
# total_activity ~ condition + Replicate + (1 | ID)
# ============================================================

model_nb <- glmmTMB(

  total_activity ~
    condition +
    Replicate +
    (1 | ID),

  family = nbinom2(
    link = "log"
  ),

  data = df_total_activity
)


cat(
  "\n\n============================================\n",
  "PRIMARY MODEL: 24-H NEGATIVE-BINOMIAL GLMM\n",
  "============================================\n"
)

print(
  summary(
    model_nb
  )
)


# Save model summary
capture.output(
  summary(
    model_nb
  ),
  file = file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_total_activity_no_LD_NB_GLMM_summary.txt"
    )
  )
)


# Save fixed effects
nb_fixed_effects <- as.data.frame(
  summary(
    model_nb
  )$coefficients$cond
)

nb_fixed_effects$term <- rownames(
  nb_fixed_effects
)

rownames(
  nb_fixed_effects
) <- NULL

nb_fixed_effects <- nb_fixed_effects %>%
  relocate(
    term
  )


write.csv(
  nb_fixed_effects,
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_total_activity_no_LD_NB_GLMM_fixed_effects.csv"
    )
  ),
  row.names = FALSE
)


# ============================================================
# 11) PRIMARY MODEL — CONVERGENCE CHECKS
# ============================================================

cat(
  "\nPositive-definite Hessian:\n"
)

print(
  model_nb$sdr$pdHess
)

cat(
  "\nglmmTMB diagnose():\n"
)

diagnose(
  model_nb
)


# ============================================================
# 12) PRIMARY MODEL — DHARMa DIAGNOSTICS
# ============================================================

set.seed(
  123
)

sim_nb <- simulateResiduals(
  fittedModel = model_nb,
  n = N_SIMULATIONS
)


# Save diagnostic plots
pdf(
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_total_activity_no_LD_NB_GLMM_DHARMa_diagnostics.pdf"
    )
  ),
  width = 10,
  height = 6
)

plot(
  sim_nb
)

plotResiduals(
  sim_nb,
  form = df_total_activity$condition
)

plotResiduals(
  sim_nb,
  form = df_total_activity$Replicate
)

dev.off()


# Show main plot interactively
plot(
  sim_nb
)


# Formal tests
dh_uniformity <- testUniformity(
  sim_nb
)

dh_dispersion <- testDispersion(
  sim_nb
)

dh_zero_inflation <- testZeroInflation(
  sim_nb
)

dh_outliers <- testOutliers(
  sim_nb
)


cat(
  "\nDHARMa uniformity test:\n"
)

print(
  dh_uniformity
)

cat(
  "\nDHARMa dispersion test:\n"
)

print(
  dh_dispersion
)

cat(
  "\nDHARMa zero-inflation test:\n"
)

print(
  dh_zero_inflation
)

cat(
  "\nDHARMa outlier test:\n"
)

print(
  dh_outliers
)


capture.output(

  {

    cat(
      "DHARMa uniformity test\n"
    )

    print(
      dh_uniformity
    )


    cat(
      "\nDHARMa dispersion test\n"
    )

    print(
      dh_dispersion
    )


    cat(
      "\nDHARMa zero-inflation test\n"
    )

    print(
      dh_zero_inflation
    )


    cat(
      "\nDHARMa outlier test\n"
    )

    print(
      dh_outliers
    )
  },

  file = file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_total_activity_no_LD_NB_GLMM_DHARMa_tests.txt"
    )
  )
)


# ============================================================
# 13) GLOBAL 24-H BEFORE–AFTER TEST — LRT
# ============================================================
#
# Primary question:
#
# Does total locomotor activity across the complete 24-h
# period differ between baseline and exposure?
#
# Full model:
#   total_activity ~ condition + Replicate + (1 | ID)
#
# Null model:
#   total_activity ~ Replicate + (1 | ID)
# ============================================================

model_nb_no_condition <- glmmTMB(

  total_activity ~
    Replicate +
    (1 | ID),

  family = nbinom2(
    link = "log"
  ),

  data = df_total_activity
)


activity_condition_LRT <- anova(
  model_nb_no_condition,
  model_nb
)


cat(
  "\n\n============================================\n",
  "GLOBAL 24-H BEFORE–AFTER LRT\n",
  "============================================\n"
)

print(
  activity_condition_LRT
)


# ------------------------------------------------------------
# AIC comparison
# ------------------------------------------------------------

activity_aic_comparison <- AIC(
  model_nb_no_condition,
  model_nb
)


cat(
  "\nAIC comparison:\n"
)

print(
  activity_aic_comparison
)


# ------------------------------------------------------------
# Extract LRT result
# ------------------------------------------------------------

activity_condition_LRT_df <- as.data.frame(
  activity_condition_LRT
)


activity_condition_test <- tibble(

  Chisq = activity_condition_LRT_df[
    2,
    "Chisq"
  ],

  df_LRT = activity_condition_LRT_df[
    2,
    "Chi Df"
  ],

  p_LRT = activity_condition_LRT_df[
    2,
    "Pr(>Chisq)"
  ]
)


cat(
  "\nExtracted LRT result:\n"
)

print(
  activity_condition_test
)


# Save LRT comparison
write.csv(
  activity_condition_LRT_df,
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_total_activity_no_LD_condition_LRT.csv"
    )
  ),
  row.names = TRUE
)


write.csv(
  as.data.frame(
    activity_aic_comparison
  ),
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_total_activity_no_LD_AIC_comparison.csv"
    )
  ),
  row.names = TRUE
)


# ============================================================
# 14) EMMEANS — EFFECT ESTIMATES AND CONFIDENCE INTERVALS
# ============================================================
#
# emmeans is used for:
#
#   - model-estimated marginal means
#   - before / after ratio
#   - 95% confidence intervals
#
# It is NOT used for the primary p-value.
#
# Primary inference comes from the global LRT above.
# ============================================================


# ------------------------------------------------------------
# 14.1 Estimated marginal means
# ------------------------------------------------------------

emm_condition <- emmeans(
  model_nb,
  ~ condition,
  component = "cond"
)


cat(
  "\n\n============================================\n",
  "EMMEANS: ORIGINAL MODEL-SCALE OUTPUT\n",
  "============================================\n"
)

print(
  emm_condition
)


# ------------------------------------------------------------
# 14.2 Estimated marginal means on response scale
# ------------------------------------------------------------

emm_condition_response <- summary(
  emm_condition,
  type = "response",
  infer = c(
    TRUE,
    FALSE
  )
)


cat(
  "\n\n============================================\n",
  "EMMEANS: MODEL-ESTIMATED TOTAL 24-H ACTIVITY\n",
  "============================================\n"
)

print(
  emm_condition_response
)


# ------------------------------------------------------------
# 14.3 Before-vs-after contrast
# ------------------------------------------------------------
#
# Because before is the reference / first factor level,
# pairs() returns:
#
#   before / after
#
# on the response scale for the log-link model.
# ============================================================

condition_pair <- pairs(
  emm_condition
)


cat(
  "\n\n============================================\n",
  "BEFORE–AFTER CONTRAST: ORIGINAL EMMEANS OUTPUT\n",
  "============================================\n"
)

print(
  condition_pair
)


# ------------------------------------------------------------
# 14.4 Response-scale effect estimate
# ------------------------------------------------------------
#
# before / after = 1
#   -> no difference
#
# before / after < 1
#   -> before lower than after
#
# before / after > 1
#   -> before higher than after
#
# Confidence intervals come from emmeans.
# The p-value reported for the primary hypothesis comes from
# the LRT rather than the Wald test.
# ============================================================

condition_pair_response <- summary(
  condition_pair,
  type = "response",
  infer = c(
    TRUE,
    FALSE
  ),
  adjust = "none"
)


cat(
  "\n\n============================================\n",
  "BEFORE–AFTER EFFECT ESTIMATE\n",
  "Ratio + 95% confidence interval\n",
  "Primary p-value comes from LRT\n",
  "============================================\n"
)

print(
  condition_pair_response
)


# ============================================================
# 15) FINAL REPORTING TABLE
# ============================================================

condition_effects <- as.data.frame(
  condition_pair_response
)


condition_contrast_final <- condition_effects %>%

  mutate(

    Chisq = activity_condition_test$Chisq,

    df_LRT = activity_condition_test$df_LRT,

    p_LRT = activity_condition_test$p_LRT
  )


cat(
  "\n\n============================================\n",
  "FINAL CONTRAST TABLE FOR REPORTING\n",
  "EMMEANS EFFECT ESTIMATE + LRT P-VALUE\n",
  "============================================\n"
)

print(
  condition_contrast_final
)


# ============================================================
# 16) EXPORT EMMEANS AND FINAL CONTRAST
# ============================================================

write.csv(
  as.data.frame(
    emm_condition_response
  ),
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_total_activity_no_LD_NB_GLMM_emmeans_before_after.csv"
    )
  ),
  row.names = FALSE
)


write.csv(
  condition_contrast_final,
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_total_activity_no_LD_NB_GLMM_before_after_LRT_final.csv"
    )
  ),
  row.names = FALSE
)


# ============================================================
# 17) RESULTS PLOT — TOTAL 24-H ACTIVITY
# ============================================================
#
# Same visual structure and colors as the no-LD Total Sleep
# figure:
#
# - individual flies
# - model-estimated means
# - model-based 95% confidence intervals
# - global 24-h LRT p-value in subtitle
#
# The sqrt y-axis is retained for activity because the
# count data are strongly right-skewed. No observations are
# removed by the transformation.
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


p_condition <- activity_condition_test$p_LRT

p_condition_text <- format_p(
  p_condition
)


# ------------------------------------------------------------
# Prepare EMM dataframe
# ------------------------------------------------------------

emm_plot_df <- as.data.frame(
  emm_condition_response
) %>%

  mutate(

    condition = factor(
      condition,
      levels = CONDITION_LEVELS
    )
  )


# ------------------------------------------------------------
# Layout helper
# ------------------------------------------------------------

nudge_amt <- 0.25


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

total_activity_plot <- ggplot(
  df_total_activity,
  aes(
    x = condition,
    y = total_activity
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
      ymin = asymp.LCL,
      ymax = asymp.UCL
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
      ymin = asymp.LCL,
      ymax = asymp.UCL
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
      y = response,
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

  # sqrt y-axis:
  # keeps all observations visible while reducing compression
  # caused by strongly right-skewed activity counts
  scale_y_continuous(
    trans = "sqrt",
    breaks = scales::breaks_pretty(
      n = 7
    ),
    expand = expansion(
      mult = c(
        0.02,
        0.06
      )
    )
  ) +

  labs(

    x = NULL,

    y = "Total activity (beam-crossing counts)",

    title = paste0(
      FRAGRANCE_NAME,
      " ",
      CONCENTRATION_FILE,
      ": total 24-h locomotor activity"
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
    ),

    legend.margin = margin(
      t = 2,
      b = 0
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
  total_activity_plot
)


# ============================================================
# 18) SAVE PLOT
# ============================================================

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_total_activity_no_LD_results_sqrt.png"
    )
  ),
  plot = total_activity_plot,
  width = 8,
  height = 5.3,
  units = "in",
  dpi = 300,
  bg = "white"
)


ggsave(
  filename = file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_total_activity_no_LD_results_sqrt.pdf"
    )
  ),
  plot = total_activity_plot,
  width = 8,
  height = 5.3,
  units = "in",
  bg = "white"
)


# ============================================================
# 19) SAVE FITTED MODELS
# ============================================================

saveRDS(
  model_nb,
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_total_activity_no_LD_NB_GLMM_full.rds"
    )
  )
)


saveRDS(
  model_nb_no_condition,
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_total_activity_no_LD_NB_GLMM_null.rds"
    )
  )
)


# ============================================================
# 20) SAVE COMPLETE STATISTICAL OUTPUT
# ============================================================

statistics_file <- file.path(
  OUTPUT_DIR,
  paste0(
    FRAGRANCE_NAME,
    "_",
    CONCENTRATION_FILE,
    "_total_activity_no_LD_statistics.txt"
  )
)


capture.output(

  {

    cat(
      "TOTAL 24-H LOCOMOTOR ACTIVITY — ",
      "NEGATIVE-BINOMIAL GLMM WITHOUT LD PHASE\n\n",
      sep = ""
    )


    cat(
      "Model:\n",
      "total_activity ~ condition + ",
      "Replicate + (1 | ID)\n\n"
    )


    cat(
      "====================\n",
      "PRIMARY MODEL\n",
      "====================\n"
    )

    print(
      summary(
        model_nb
      )
    )


    cat(
      "\n\n====================\n",
      "CONVERGENCE\n",
      "====================\n"
    )

    cat(
      "\nPositive-definite Hessian:\n"
    )

    print(
      model_nb$sdr$pdHess
    )


    cat(
      "\n\n====================\n",
      "DHARMa DIAGNOSTICS\n",
      "====================\n"
    )

    print(
      dh_uniformity
    )

    print(
      dh_dispersion
    )

    print(
      dh_zero_inflation
    )

    print(
      dh_outliers
    )


    cat(
      "\n\n====================\n",
      "GLOBAL 24-H BEFORE–AFTER LRT\n",
      "====================\n"
    )

    print(
      activity_condition_LRT
    )


    cat(
      "\n\n====================\n",
      "EMMEANS — RESPONSE SCALE\n",
      "====================\n"
    )

    print(
      emm_condition_response
    )


    cat(
      "\n\n====================\n",
      "BEFORE / AFTER RATIO + 95% CI + LRT\n",
      "====================\n"
    )

    print(
      condition_contrast_final
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
# 21) DONE
# ============================================================

cat(
  "\n============================================\n",
  "TOTAL 24-H ACTIVITY ANALYSIS COMPLETE\n",
  "============================================\n",
  "\nResults saved in:\n",
  normalizePath(
    OUTPUT_DIR
  ),
  "\n"
)

