# ============================================================
# DAM2 TOTAL ACTIVITY — FINAL ANALYSIS WORKFLOW
#
# MAIN MODEL:
#   Negative-binomial GLMM (nbinom2)
#
#
# Input:
# one row = one fly × one condition × one phase
#
# Expected columns:
# Fragrance, Concentration, ID, Replicate, condition,
# Phase, total_activity
#
# IMPORTANT:
# - Total activity is analysed with condition × Phase.
# - Primary planned contrasts compare before vs after within Light and within Dark.
# - before is the reference condition.
# - Light is the reference phase.
# - R1 is the reference replicate.
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

CONDITION_LEVELS <- c("before", "after")
PHASE_LEVELS <- c("Light", "Dark")
REPLICATE_LEVELS <- c("R1", "R2")

# Number of DHARMa simulations
N_SIMULATIONS <- 1000

# Output folder
OUTPUT_DIR <- file.path(
  paste0("total_activity_statistics_", FRAGRANCE_NAME),
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

# Remove accidental CSV index column if present
if ("...1" %in% names(df_total_activity_phase)) {
  df_total_activity_phase <- df_total_activity_phase %>%
    select(-`...1`)
}

cat(
  "\nImported file:", DATA_FILE,
  "\nRows:", nrow(df_total_activity_phase),
  "| Columns:", ncol(df_total_activity_phase),
  "\nOutput folder:", normalizePath(OUTPUT_DIR),
  "\n\n"
)

print(names(df_total_activity_phase))
print(head(df_total_activity_phase))


# ============================================================
# 4) CHECK REQUIRED COLUMNS
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

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}


# ============================================================
# 5) PREPARE VARIABLES
# ============================================================

df_total_activity_phase <- df_total_activity_phase %>%
  mutate(
    ID = factor(ID),
    
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
    
    total_activity = as.integer(total_activity),
    
    # Sensitivity-analysis response
    log_total_activity = log1p(total_activity)
  )


# ============================================================
# 6) BASIC DATA CHECKS
# ============================================================

cat("\n================ DATA CHECKS ================\n")

cat("\nNumber of flies:\n")
print(n_distinct(df_total_activity_phase$ID))

cat("\nCondition counts:\n")
print(table(df_total_activity_phase$condition))

cat("\nPhase counts:\n")
print(table(df_total_activity_phase$Phase))

cat("\nReplicate counts:\n")
print(table(df_total_activity_phase$Replicate))

cat("\nCondition × Phase counts:\n")
print(
  table(
    df_total_activity_phase$condition,
    df_total_activity_phase$Phase
  )
)

cat("\nObservations per fly:\n")
print(
  table(
    table(df_total_activity_phase$ID)
  )
)

cat("\nTotal activity summary:\n")
print(summary(df_total_activity_phase$total_activity))

cat("\nNegative total activity values:\n")
print(
  any(
    df_total_activity_phase$total_activity < 0,
    na.rm = TRUE
  )
)

cat("\nMissing total activity values:\n")
print(anyNA(df_total_activity_phase$total_activity))

cat("\nInteger counts only:\n")
print(
  all(
    df_total_activity_phase$total_activity ==
      floor(df_total_activity_phase$total_activity),
    na.rm = TRUE
  )
)

# Every ID × condition × Phase combination should normally occur once
duplicate_mask <- duplicated(
  df_total_activity_phase[
    c(
      "ID",
      "condition",
      "Phase"
    )
  ]
)

cat(
  "\nDuplicated ID × condition × Phase rows:",
  sum(duplicate_mask),
  "\n"
)

if (sum(duplicate_mask) > 0) {
  warning(
    "Duplicated ID × condition × Phase observations were found."
  )
}


# ============================================================
# 7) DESCRIPTIVE STATISTICS
# ============================================================

total_activity_summary <- df_total_activity_phase %>%
  group_by(
    Replicate,
    condition,
    Phase
  ) %>%
  summarise(
    n = n(),
    
    n_flies = n_distinct(ID),
    
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
      sd_total_activity / sqrt(n_flies),
    
    .groups = "drop"
  )

print(total_activity_summary)

write.csv(
  total_activity_summary,
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_total_activity_descriptive_statistics.csv"
    )
  ),
  row.names = FALSE
)


# ============================================================
# 8) PRIMARY MODEL — NEGATIVE-BINOMIAL GLMM
# ============================================================

model_nb <- glmmTMB(
  total_activity ~
    condition * Phase +
    Replicate +
    (1 | ID),
  
  family = nbinom2(
    link = "log"
  ),
  
  data = df_total_activity_phase
)

cat(
  "\n\n============================================\n",
  "PRIMARY MODEL: NEGATIVE-BINOMIAL GLMM\n",
  "============================================\n"
)

print(summary(model_nb))


# Save model summary
capture.output(
  summary(model_nb),
  file = file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_NB_GLMM_summary.txt"
    )
  )
)


# Save fixed effects
nb_fixed_effects <- as.data.frame(
  summary(model_nb)$coefficients$cond
)

nb_fixed_effects$term <- rownames(
  nb_fixed_effects
)

rownames(nb_fixed_effects) <- NULL

nb_fixed_effects <- nb_fixed_effects %>%
  relocate(term)

write.csv(
  nb_fixed_effects,
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_NB_GLMM_fixed_effects.csv"
    )
  ),
  row.names = FALSE
)


# ============================================================
# 9) PRIMARY MODEL — CONVERGENCE CHECKS
# ============================================================

cat("\nPositive-definite Hessian:\n")
print(model_nb$sdr$pdHess)

cat("\nglmmTMB diagnose():\n")
diagnose(model_nb)


# ============================================================
# 10) PRIMARY MODEL — DHARMa DIAGNOSTICS
# ============================================================

set.seed(123)

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
      "_NB_GLMM_DHARMa_diagnostics.pdf"
    )
  ),
  width = 10,
  height = 6
)

plot(sim_nb)

plotResiduals(
  sim_nb,
  form = df_total_activity_phase$condition
)

plotResiduals(
  sim_nb,
  form = df_total_activity_phase$Phase
)

plotResiduals(
  sim_nb,
  form = df_total_activity_phase$Replicate
)

dev.off()

# Show main plot interactively
plot(sim_nb)


# Formal tests
dh_uniformity <- testUniformity(sim_nb)
dh_dispersion <- testDispersion(sim_nb)
dh_zero_inflation <- testZeroInflation(sim_nb)
dh_outliers <- testOutliers(sim_nb)

cat("\nDHARMa uniformity test:\n")
print(dh_uniformity)

cat("\nDHARMa dispersion test:\n")
print(dh_dispersion)

cat("\nDHARMa zero-inflation test:\n")
print(dh_zero_inflation)

cat("\nDHARMa outlier test:\n")
print(dh_outliers)

capture.output(
  {
    cat("DHARMa uniformity test\n")
    print(dh_uniformity)
    
    cat("\nDHARMa dispersion test\n")
    print(dh_dispersion)
    
    cat("\nDHARMa zero-inflation test\n")
    print(dh_zero_inflation)
    
    cat("\nDHARMa outlier test\n")
    print(dh_outliers)
  },
  file = file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_NB_GLMM_DHARMa_tests.txt"
    )
  )
)

# ============================================================
# 11) PRIMARY MODEL — TEST CONDITION × PHASE INTERACTION
# ============================================================

model_nb_no_interaction <- glmmTMB(
  total_activity ~
    condition +
    Phase +
    Replicate +
    (1 | ID),
  
  family = nbinom2(
    link = "log"
  ),
  
  data = df_total_activity_phase
)


cat(
  "\n\n============================================\n",
  "NB-GLMM MODEL COMPARISON\n",
  "============================================\n"
)


# ------------------------------------------------------------
# AIC comparison
# ------------------------------------------------------------

nb_aic_comparison <- AIC(
  model_nb,
  model_nb_no_interaction
)

print(nb_aic_comparison)


# ------------------------------------------------------------
# Likelihood-ratio test:
# Does the before/after effect differ between Light and Dark?
# ------------------------------------------------------------

nb_interaction_test <- anova(
  model_nb_no_interaction,
  model_nb
)

cat(
  "\nLikelihood-ratio test for condition x Phase interaction:\n"
)

print(nb_interaction_test)


# ------------------------------------------------------------
# Export interaction model comparison
# ------------------------------------------------------------

write.csv(
  as.data.frame(nb_aic_comparison),
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_NB_GLMM_AIC_comparison.csv"
    )
  ),
  row.names = FALSE
)


write.csv(
  as.data.frame(nb_interaction_test),
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_NB_GLMM_interaction_LRT.csv"
    )
  ),
  row.names = FALSE
)



# ============================================================
# 11B) PLANNED SIMPLE CONTRASTS — LIKELIHOOD-RATIO TESTS
# ============================================================
#
# Scientific questions:
#
#   Light: before vs after
#   Dark:  before vs after
#
# These are the two primary planned comparisons.
#
# P-values are obtained using likelihood-ratio tests (LRT),
# NOT Wald z-tests.
#
# For each test:
#
# Full model:
#   Light-before, Light-after,
#   Dark-before, Dark-after
#   are estimated separately.
#
# Null model:
#   before and after are forced to be equal only within the
#   phase being tested.
#
# The other phase remains unconstrained.
# ============================================================


# ------------------------------------------------------------
# 11B.1 Create cell factors
# ------------------------------------------------------------

df_total_activity_phase <- df_total_activity_phase %>%
  mutate(
    
    # Full 2 x 2 representation
    PhaseCondition = factor(
      paste(
        Phase,
        condition,
        sep = "_"
      )
    ),
    
    
    # --------------------------------------------------------
    # Null hypothesis for Light:
    # Light-before = Light-after
    #
    # Dark-before and Dark-after remain separate.
    # --------------------------------------------------------
    
    LRT_Light = factor(
      case_when(
        
        Phase == "Light" ~
          "Light_common",
        
        Phase == "Dark" &
          condition == "before" ~
          "Dark_before",
        
        Phase == "Dark" &
          condition == "after" ~
          "Dark_after"
      )
    ),
    
    
    # --------------------------------------------------------
    # Null hypothesis for Dark:
    # Dark-before = Dark-after
    #
    # Light-before and Light-after remain separate.
    # --------------------------------------------------------
    
    LRT_Dark = factor(
      case_when(
        
        Phase == "Dark" ~
          "Dark_common",
        
        Phase == "Light" &
          condition == "before" ~
          "Light_before",
        
        Phase == "Light" &
          condition == "after" ~
          "Light_after"
      )
    )
  )



# ------------------------------------------------------------
# 11B.2 Full cell-means model
# ------------------------------------------------------------
#
# This model is statistically equivalent to:
#
# condition * Phase
#
# It is only parameterized differently so that the
# phase-specific null hypotheses can be imposed easily.
# ------------------------------------------------------------

model_nb_cells <- glmmTMB(
  total_activity ~
    PhaseCondition +
    Replicate +
    (1 | ID),
  
  family = nbinom2(
    link = "log"
  ),
  
  data = df_total_activity_phase
)


# ------------------------------------------------------------
# Check equivalence with original model
# ------------------------------------------------------------

cat(
  "\n\n============================================\n",
  "FULL MODEL PARAMETERIZATION CHECK\n",
  "============================================\n"
)

logLik_check <- c(
  
  original_model = as.numeric(
    logLik(model_nb)
  ),
  
  cell_model = as.numeric(
    logLik(model_nb_cells)
  )
)

print(logLik_check)



# ============================================================
# 11B.3 LIGHT — before vs after
# ============================================================

model_nb_Light_null <- glmmTMB(
  total_activity ~
    LRT_Light +
    Replicate +
    (1 | ID),
  
  family = nbinom2(
    link = "log"
  ),
  
  data = df_total_activity_phase
)


lrt_Light <- anova(
  model_nb_Light_null,
  model_nb_cells
)


cat(
  "\n\n============================================\n",
  "LRT: LIGHT — BEFORE VS AFTER\n",
  "============================================\n"
)

print(lrt_Light)



# ============================================================
# 11B.4 DARK — before vs after
# ============================================================

model_nb_Dark_null <- glmmTMB(
  total_activity ~
    LRT_Dark +
    Replicate +
    (1 | ID),
  
  family = nbinom2(
    link = "log"
  ),
  
  data = df_total_activity_phase
)


lrt_Dark <- anova(
  model_nb_Dark_null,
  model_nb_cells
)


cat(
  "\n\n============================================\n",
  "LRT: DARK — BEFORE VS AFTER\n",
  "============================================\n"
)

print(lrt_Dark)



# ============================================================
# 11B.5 Extract LRT p-values
# ============================================================

lrt_Light_df <- as.data.frame(
  lrt_Light
)

lrt_Dark_df <- as.data.frame(
  lrt_Dark
)


# Second row = comparison with the full model
p_Light_LRT <- lrt_Light_df[
  2,
  "Pr(>Chisq)"
]

p_Dark_LRT <- lrt_Dark_df[
  2,
  "Pr(>Chisq)"
]



# ============================================================
# 11B.6 Holm correction
# ============================================================
#
# Light and Dark are treated as one family of
# two planned hypothesis tests.
# ============================================================

p_LRT_Holm <- p.adjust(
  c(
    p_Light_LRT,
    p_Dark_LRT
  ),
  method = "holm"
)



# ============================================================
# 11B.7 Final LRT result table
# ============================================================
#
# IMPORTANT:
#
# p_LRT   = raw likelihood-ratio-test p-value
#
# p.value = Holm-adjusted LRT p-value
#
# The column is deliberately called p.value so that the
# existing plotting code can continue using:
#
# contrast_test_holm$p.value
# ============================================================

contrast_test_holm <- tibble(
  
  Phase = factor(
    c(
      "Light",
      "Dark"
    ),
    levels = PHASE_LEVELS
  ),
  
  contrast = c(
    "before / after",
    "before / after"
  ),
  
  Chisq = c(
    lrt_Light_df[2, "Chisq"],
    lrt_Dark_df[2, "Chisq"]
  ),
  
  df = c(
    lrt_Light_df[2, "Chi Df"],
    lrt_Dark_df[2, "Chi Df"]
  ),
  
  p_LRT = c(
    p_Light_LRT,
    p_Dark_LRT
  ),
  
  p.value = p_LRT_Holm
)


cat(
  "\n\n============================================\n",
  "FINAL PLANNED CONTRASTS — LRT + HOLM\n",
  "============================================\n"
)

print(
  contrast_test_holm
)



# ------------------------------------------------------------
# Export LRT results
# ------------------------------------------------------------

write.csv(
  as.data.frame(
    contrast_test_holm
  ),
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_NB_GLMM_planned_contrasts_LRT_Holm.csv"
    )
  ),
  row.names = FALSE
)



# ============================================================
# 12) EMMEANS — EFFECT ESTIMATES AND CONFIDENCE INTERVALS
# ============================================================
#
# emmeans is now used for:
#
#   - model-estimated marginal means
#   - before / after ratios
#   - 95% confidence intervals
#
# It is NOT used for the primary p-values.
#
# Primary p-values come from the LRTs above.
# ============================================================



# ------------------------------------------------------------
# 12.1 Estimated marginal means
# ------------------------------------------------------------

emm_condition_by_phase <- emmeans(
  model_nb,
  ~ condition | Phase,
  component = "cond"
)


cat(
  "\n\n============================================\n",
  "EMMEANS: ORIGINAL MODEL-SCALE OUTPUT\n",
  "============================================\n"
)


# Inspect the original emmeans object first
print(
  emm_condition_by_phase
)



# ------------------------------------------------------------
# 12.2 Estimated marginal means on response scale
# ------------------------------------------------------------
#
# type = "response":
#
# Back-transforms the estimates from the log scale to
# the original total-activity count scale.
#
# infer = c(TRUE, FALSE):
#
# TRUE  -> calculate 95% confidence intervals
#
# FALSE -> do NOT perform Wald hypothesis tests
# ------------------------------------------------------------

emm_condition_response <- summary(
  emm_condition_by_phase,
  type = "response",
  infer = c(
    TRUE,
    FALSE
  )
)


cat(
  "\n\n============================================\n",
  "EMMEANS: MODEL-ESTIMATED TOTAL ACTIVITY\n",
  "============================================\n"
)


print(
  emm_condition_response
)



# ------------------------------------------------------------
# 12.3 Planned before-vs-after contrasts within each Phase
# ------------------------------------------------------------

condition_pairs <- pairs(
  emm_condition_by_phase
)


cat(
  "\n\n============================================\n",
  "PLANNED CONTRASTS: ORIGINAL EMMEANS OUTPUT\n",
  "============================================\n"
)


# Model-scale contrast
print(
  condition_pairs
)



# ------------------------------------------------------------
# 12.4 Response-scale effect estimates
# ------------------------------------------------------------
#
# The model uses a log link.
#
# Therefore:
#
# difference on log scale
#
# becomes:
#
# ratio on response scale.
#
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
#
# IMPORTANT:
#
# infer = c(TRUE, FALSE)
#
# We calculate confidence intervals,
# but NO Wald p-values.
# ------------------------------------------------------------

condition_pairs_response <- summary(
  condition_pairs,
  type = "response",
  infer = c(
    TRUE,
    FALSE
  ),
  adjust = "none"
)


cat(
  "\n\n============================================\n",
  "PLANNED CONTRASTS: EFFECT ESTIMATES\n",
  "Ratios + 95% confidence intervals\n",
  "NO Wald p-values\n",
  "============================================\n"
)


print(
  condition_pairs_response
)



# ============================================================
# 12.5 Combine EMM effect estimates with LRT inference
# ============================================================
#
# Final reporting table contains:
#
#   ratio
#   95% CI
#   LRT chi-square
#   raw LRT p-value
#   Holm-adjusted LRT p-value
#
# ============================================================

condition_effects <- as.data.frame(
  condition_pairs_response
)


condition_contrasts_final <- condition_effects %>%
  
  left_join(
    
    contrast_test_holm %>%
      select(
        Phase,
        contrast,
        Chisq,
        df,
        p_LRT,
        p.value
      ),
    
    by = c(
      "Phase",
      "contrast"
    )
  ) %>%
  
  rename(
    p_LRT_Holm = p.value
  )



cat(
  "\n\n============================================\n",
  "FINAL CONTRAST TABLE FOR REPORTING\n",
  "EMMEANS EFFECT ESTIMATES + LRT P-VALUES\n",
  "============================================\n"
)


print(
  condition_contrasts_final
)



# ============================================================
# 12.6 Export EMMs and final planned contrasts
# ============================================================


# ------------------------------------------------------------
# Estimated marginal means
# ------------------------------------------------------------

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
      "_NB_GLMM_emmeans_before_after_by_phase.csv"
    )
  ),
  row.names = FALSE
)



# ------------------------------------------------------------
# Final contrast table:
# ratios + CI + LRT + Holm
# ------------------------------------------------------------

write.csv(
  condition_contrasts_final,
  file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_NB_GLMM_planned_contrasts_LRT_Holm_final.csv"
    )
  ),
  row.names = FALSE
)



# ============================================================
# 12B) RESULTS FIGURE
# ============================================================
#
# Figure contains:
#
# - raw observations
# - model-estimated marginal means
# - model-based 95% confidence intervals
# - Holm-adjusted LRT p-values
#
# IMPORTANT:
#
# EMMs + confidence intervals:
#   obtained from emmeans
#
# P-values:
#   obtained from likelihood-ratio tests
#   and adjusted using Holm
#
# ============================================================



# ------------------------------------------------------------
# Prepare EMM dataframe
# ------------------------------------------------------------

emm_plot_df <- as.data.frame(
  emm_condition_response
)


emm_plot_df <- emm_plot_df %>%
  mutate(
    
    condition = factor(
      condition,
      levels = CONDITION_LEVELS
    ),
    
    Phase = factor(
      Phase,
      levels = PHASE_LEVELS
    )
  )



# ------------------------------------------------------------
# Prepare Holm-adjusted LRT p-value labels
# ------------------------------------------------------------

pvalue_plot_df <- contrast_test_holm %>%
  
  mutate(
    
    Phase = factor(
      Phase,
      levels = PHASE_LEVELS
    ),
    
    p_label = paste0(
      "Holm-adjusted LRT p = ",
      format.pval(
        p.value,
        digits = 3,
        eps = 0.001
      )
    )
  ) %>%
  
  select(
    Phase,
    p.value,
    p_label
  )



# ------------------------------------------------------------
# Position labels above highest raw observation
# ------------------------------------------------------------

phase_y_positions <- df_total_activity_phase %>%
  
  group_by(
    Phase
  ) %>%
  
  summarise(
    
    y_position = max(
      total_activity,
      na.rm = TRUE
    ) * 1.08,
    
    .groups = "drop"
  )


pvalue_plot_df <- pvalue_plot_df %>%
  
  left_join(
    phase_y_positions,
    by = "Phase"
  )




# ============================================================
# RESULTS PLOT — TOTAL ACTIVITY BY LIGHT/DARK PHASE
# ============================================================
#
# Same visual style as Total_activity_no_LD:
#
# - same colors
# - individual flies
# - model-estimated means
# - model-based 95% confidence intervals
# - Holm-adjusted LRT p-values
#
# y-axis remains sqrt-transformed.
# ============================================================


condition_colors <- c(
  "before" = "#B79F00",
  "after"  = "#2F7F78"
)


# ------------------------------------------------------------
# P-values for subtitle
# ------------------------------------------------------------

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
      formatC(
        p,
        format = "f",
        digits = 3
      )
    )
  )
}


p_light_text <- format_p(
  p_light
)

p_dark_text <- format_p(
  p_dark
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
    ),
    
    Phase = factor(
      Phase,
      levels = PHASE_LEVELS
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
  df_total_activity_phase,
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
  
  # Light / Dark panels
  facet_wrap(
    ~ Phase
  ) +
  
  # Colors
  scale_color_manual(
    values = condition_colors,
    guide = "none"
  ) +
  
  scale_fill_manual(
    values = condition_colors,
    guide = "none"
  ) +
  
  # Legend symbols
  scale_shape_manual(
    name = NULL,
    values = c(
      "Individual flies" = 16,
      "Model-estimated mean ± 95% CI" = 21
    )
  ) +
  
  # sqrt y-axis
  scale_y_continuous(
    trans = "sqrt",
    breaks = c(
      0,
      100,
      250,
      500,
      750,
      1000,
      1250,
      1500
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
      ": total locomotor activity"
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
    
    panel.spacing = unit(
      0.3,
      "cm"
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
# SAVE PLOT
# ============================================================

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    paste0(
      FRAGRANCE_NAME,
      "_",
      CONCENTRATION_FILE,
      "_total_activity_results_sqrt.png"
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
      "_total_activity_results_sqrt.pdf"
    )
  ),
  plot = total_activity_plot,
  width = 8,
  height = 5.3,
  units = "in",
  bg = "white"
)