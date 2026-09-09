# ============================================================
# SLEEP LATENCY — IMPORT AND DATA PREPARATION
# ============================================================


# ============================================================
# 1) PACKAGES
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# Packages for models later
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(performance)
library(lme4)


# ============================================================
# 2) SETTINGS
# ============================================================

FRAGRANCE_NAME <- "F2"
CONCENTRATION_FILE <- "72uL"

DATA_FILE <- paste0(
  FRAGRANCE_NAME,
  "_",
  CONCENTRATION_FILE,
  "_sleep_parameters_for_R.csv"
)

CONDITION_LEVELS <- c(
  "before",
  "after"
)

REPLICATE_LEVELS <- c(
  "R1",
  "R2"
)


OUTPUT_DIR <- file.path(
  paste0("sleep_latency_analysis_", FRAGRANCE_NAME),
  CONCENTRATION_FILE
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

cat(
  "Output folder:",
  normalizePath(OUTPUT_DIR),
  "\n"
)


# ============================================================
# 3) IMPORT DATA
# ============================================================

df_latency_raw <- read.csv(
  DATA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

df_latency_raw <- tibble::as_tibble(
  df_latency_raw,
  .name_repair = "unique"
)

cat("\n============================================\n")
cat("RAW SLEEP LATENCY DATA\n")
cat("============================================\n")

print(head(df_latency_raw))

cat(
  "\nDimensions:",
  nrow(df_latency_raw),
  "rows x",
  ncol(df_latency_raw),
  "columns\n"
)

cat("\nColumns:\n")
print(names(df_latency_raw))


# ============================================================
# 4) REMOVE ACCIDENTAL INDEX COLUMNS
# ============================================================

index_columns <- intersect(
  names(df_latency_raw),
  c(
    "...1",
    "X",
    "Unnamed: 0"
  )
)

if (length(index_columns) > 0) {
  
  df_latency_raw <- df_latency_raw %>%
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
  "sleep_latency_min"
)

missing_columns <- setdiff(
  required_columns,
  names(df_latency_raw)
)

if (length(missing_columns) > 0) {
  
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

cat("\nAll required columns are present.\n")


# ============================================================
# 6) CLEAN VARIABLES
# ============================================================

df_latency <- df_latency_raw %>%
  mutate(
    
    ID = trimws(
      as.character(ID)
    ),
    
    Replicate = trimws(
      as.character(Replicate)
    ),
    
    condition = trimws(
      tolower(
        as.character(condition)
      )
    ),
    
    sleep_latency_min = as.numeric(
      sleep_latency_min
    )
  )


# ============================================================
# 7) FACTOR LEVELS
# ============================================================

df_latency <- df_latency %>%
  mutate(
    
    ID = factor(ID),
    
    Replicate = factor(
      Replicate,
      levels = REPLICATE_LEVELS
    ),
    
    condition = factor(
      condition,
      levels = CONDITION_LEVELS
    )
  )


# ============================================================
# 8) MISSING VALUES
# ============================================================

cat("\n============================================\n")
cat("MISSING VALUES\n")
cat("============================================\n")

print(
  colSums(
    is.na(
      df_latency[
        required_columns
      ]
    )
  )
)


n_before_na_removal <- nrow(
  df_latency
)

df_latency <- df_latency %>%
  filter(
    !is.na(ID),
    !is.na(Replicate),
    !is.na(condition),
    !is.na(sleep_latency_min)
  )

cat(
  "\nRows before NA removal:",
  n_before_na_removal,
  "\n"
)

cat(
  "Rows after NA removal:",
  nrow(df_latency),
  "\n"
)


# ============================================================
# 9) BIOLOGICAL / NUMERICAL RANGE CHECK
# ============================================================
#
# Sleep latency cannot be negative.
#
# 0 min is possible:
# the fly is already classified as asleep in the first
# 5-min bin following lights-off.
#
# ============================================================

invalid_latency <- df_latency %>%
  filter(
    sleep_latency_min < 0
  )

cat("\n============================================\n")
cat("LATENCY RANGE CHECK\n")
cat("============================================\n")

cat(
  "Negative latency values:",
  nrow(invalid_latency),
  "\n"
)

if (nrow(invalid_latency) > 0) {
  print(invalid_latency)
}


# ============================================================
# 10) CHECK 5-MINUTE BIN STRUCTURE
# ============================================================
#
# Sleep latency originates from 5-min DAM bins.
#
# Therefore valid values should be:
#
# 0, 5, 10, 15, 20, ...
#
# ============================================================

non_5min_latency <- df_latency %>%
  filter(
    sleep_latency_min %% 5 != 0
  )

cat("\n============================================\n")
cat("5-MINUTE BIN CHECK\n")
cat("============================================\n")

cat(
  "Values not divisible by 5:",
  nrow(non_5min_latency),
  "\n"
)

if (nrow(non_5min_latency) > 0) {
  print(non_5min_latency)
}


# ============================================================
# 11) CREATE LATENCY COUNT VARIABLE
# ============================================================
#
# Example:
#
#  0 min  -> 0 bins
#  5 min  -> 1 bin
# 10 min  -> 2 bins
# 20 min  -> 4 bins
#
# This variable can later be used for count models.
#
# ============================================================

if (nrow(non_5min_latency) > 0) {
  
  stop(
    "Some sleep-latency values are not divisible by 5. ",
    "Check the data before creating latency_bins."
  )
}


df_latency <- df_latency %>%
  mutate(
    latency_bins = as.integer(
      sleep_latency_min / 5
    )
  )


# ============================================================
# 12) CHECK DUPLICATES
# ============================================================
#
# There should normally be one observation per:
#
# ID x condition
#
# ============================================================

duplicate_latency <- df_latency %>%
  count(
    ID,
    condition,
    name = "n"
  ) %>%
  filter(
    n > 1
  )

cat("\n============================================\n")
cat("DUPLICATE CHECK\n")
cat("============================================\n")

cat(
  "Duplicate ID x condition combinations:",
  nrow(duplicate_latency),
  "\n"
)

if (nrow(duplicate_latency) > 0) {
  print(duplicate_latency)
}


# ============================================================
# 13) CHECK OBSERVATIONS PER FLY
# ============================================================
#
# Expected:
#
# before
# after
#
# = 2 observations per fly
#
# ============================================================

observations_per_fly <- df_latency %>%
  count(
    ID,
    name = "n_observations"
  )

cat("\n============================================\n")
cat("OBSERVATIONS PER FLY\n")
cat("============================================\n")

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

if (nrow(incomplete_flies) > 0) {
  print(incomplete_flies)
}


# ============================================================
# 14) CHECK EXPERIMENTAL DESIGN
# ============================================================

cat("\n============================================\n")
cat("EXPERIMENTAL DESIGN\n")
cat("============================================\n")

cat("\nCondition:\n")
print(
  table(
    df_latency$condition
  )
)

cat("\nReplicate x condition:\n")
print(
  xtabs(
    ~ Replicate + condition,
    data = df_latency
  )
)

cat(
  "\nNumber of unique flies:",
  n_distinct(
    df_latency$ID
  ),
  "\n"
)


# ============================================================
# 15) ZERO LATENCIES
# ============================================================
#
# Zero latency is biologically possible and is NOT
# automatically treated as an outlier.
#
# ============================================================

cat("\n============================================\n")
cat("ZERO LATENCIES\n")
cat("============================================\n")

cat(
  "Total number of zero-latency observations:",
  sum(
    df_latency$sleep_latency_min == 0
  ),
  "\n"
)

print(
  table(
    df_latency$condition,
    df_latency$sleep_latency_min == 0
  )
)


# ============================================================
# 16) DESCRIPTIVE SUMMARY
# ============================================================

latency_summary <- df_latency %>%
  group_by(
    condition
  ) %>%
  summarise(
    
    n = n(),
    
    mean_min = mean(
      sleep_latency_min
    ),
    
    sd_min = sd(
      sleep_latency_min
    ),
    
    median_min = median(
      sleep_latency_min
    ),
    
    IQR_min = IQR(
      sleep_latency_min
    ),
    
    min_min = min(
      sleep_latency_min
    ),
    
    max_min = max(
      sleep_latency_min
    ),
    
    zero_latency_n = sum(
      sleep_latency_min == 0
    ),
    
    .groups = "drop"
  )

cat("\n============================================\n")
cat("DESCRIPTIVE SUMMARY\n")
cat("============================================\n")

print(latency_summary)


# ============================================================
# 17) DISTRIBUTION — MINUTES
# ============================================================

hist(
  df_latency$sleep_latency_min,
  breaks = 20,
  main = "Distribution of sleep latency",
  xlab = "Sleep latency (min)"
)


# ============================================================
# 18) DISTRIBUTION — 5-MINUTE BINS
# ============================================================

hist(
  df_latency$latency_bins,
  breaks = 20,
  main = "Distribution of sleep latency bins",
  xlab = "Sleep latency (5-min bins)"
)


# ============================================================
# 19) FINAL DATAFRAME FOR MODELLING
# ============================================================

df_sleep_latency <- df_latency %>%
  select(
    ID,
    Replicate,
    condition,
    sleep_latency_min,
    latency_bins,
    everything()
  )


cat("\n============================================\n")
cat("FINAL DATAFRAME FOR MODELLING\n")
cat("============================================\n")

print(
  head(
    df_sleep_latency,
    8
  )
)

str(
  df_sleep_latency
)


# ============================================================
# MODEL SELECTION / FITTING NEXT
# ============================================================


mixed_model <- lmer(data = df_sleep_latency,
                    sleep_latency_min ~ condition +
                      Replicate +
                      (1|ID))
summary(mixed_model)

performance(mixed_model)
check_model(mixed_model)

mixed_model_counts <- lmer(data = df_sleep_latency,
                           latency_bins ~ condition +
                             Replicate +
                             (1|ID))

summary(mixed_model_counts)

# Singularity
isSingular(mixed_model_counts, tol = 1e-4)

# Residuals vs fitted
plot(
  fitted(mixed_model_counts),
  residuals(mixed_model_counts),
  xlab = "Fitted values",
  ylab = "Residuals"
)
abline(h = 0, lty = 2)
lines(lowess(fitted(mixed_model_counts), residuals(mixed_model_counts)), lwd = 2)

# Q-Q plot of residuals
qqnorm(residuals(mixed_model_counts))
qqline(residuals(mixed_model_counts))

# Q-Q plot of random intercepts
random_intercepts <- ranef(mixed_model_counts)$ID[, 1]
qqnorm(random_intercepts)
qqline(random_intercepts)

# Heteroscedasticity
performance::check_heteroscedasticity(mixed_model_counts)

# Influential observations
performance::check_outliers(mixed_model_counts)

# Optional overall diagnostic panel
performance::check_model(mixed_model_counts)


#==========================
# GLMM Models #
#==========================

# GLMM formulas
formula_latency <- latency_bins ~ condition + Replicate + (1 | ID)

# Poisson
model_pois <- glmmTMB(
  formula_latency,
  family = poisson(link = "log"),
  data = df_sleep_latency
)

# Negative binomial
model_nb <- glmmTMB(
  formula_latency,
  family = nbinom2(link = "log"),
  data = df_sleep_latency
)

# Zero-inflated Poisson
model_zip <- glmmTMB(
  formula_latency,
  ziformula = ~1,
  family = poisson(link = "log"),
  data = df_sleep_latency
)

# Zero-inflated negative binomial
model_zinb <- glmmTMB(
  formula_latency,
  ziformula = ~1,
  family = nbinom2(link = "log"),
  data = df_sleep_latency
)

# ============================================================
# MODEL SELECTION + FINAL INFERENCE
# MANUAL SELECTION OF BEST MODEL
# ============================================================


# ============================================================
# 1) MODEL COMPARISON
# ============================================================

model_comparison <- data.frame(
  Model = c(
    "poisson",
    "nb",
    "zip",
    "zinb"
  ),
  
  AIC = c(
    AIC(model_pois),
    AIC(model_nb),
    AIC(model_zip),
    AIC(model_zinb)
  ),
  
  BIC = c(
    BIC(model_pois),
    BIC(model_nb),
    BIC(model_zip),
    BIC(model_zinb)
  )
)

print(model_comparison)


# ============================================================
# 2) SELECT FINAL MODEL MANUALLY
# ============================================================
# Choose after checking:
# - AIC / BIC
# - convergence
# - DHARMa diagnostics
#
# Options:
# "poisson"
# "nb"
# "zip"
# "zinb"

BEST_MODEL_NAME <- "nb"


model_list <- list(
  poisson = model_pois,
  nb      = model_nb,
  zip     = model_zip,
  zinb    = model_zinb
)


best_model <- model_list[[BEST_MODEL_NAME]]


cat(
  "\nSelected final model:",
  BEST_MODEL_NAME,
  "\n"
)

summary(best_model)


# ============================================================
# 3) FINAL MODEL CONVERGENCE
# ============================================================

cat(
  "\nOptimizer convergence:\n"
)

print(
  best_model$fit$convergence
)


cat(
  "\nPositive definite Hessian:\n"
)

print(
  best_model$sdr$pdHess
)


glmmTMB::diagnose(
  best_model
)


# ============================================================
# 4) FINAL MODEL DHARMa DIAGNOSTICS
# ============================================================

set.seed(123)

res_final <- simulateResiduals(
  fittedModel = best_model,
  n = 2000
)


plot(
  res_final
)


testUniformity(
  res_final
)


testDispersion(
  res_final
)


testZeroInflation(
  res_final
)


testOutliers(
  res_final,
  type = "bootstrap",
  nBoot = 1000
)


# ============================================================
# 5) CONDITION LRT
# ============================================================
# Remove condition only from the conditional count component.
# Family and zero-inflation structure remain unchanged.

model_no_condition <- update(
  best_model,
  formula. =
    latency_bins ~
    Replicate +
    (1 | ID)
)


stopifnot(
  nobs(best_model) ==
    nobs(model_no_condition)
)


lrt_condition <- anova(
  model_no_condition,
  best_model
)


print(
  lrt_condition
)


chi_square_condition <- lrt_condition$Chisq[2]

df_condition <- lrt_condition$`Chi Df`[2]

p_lrt_condition <- lrt_condition$`Pr(>Chisq)`[2]


# ============================================================
# 6) EFFECT SIZE
# ============================================================

beta_condition <- fixef(
  best_model
)$cond[
  "conditionafter"
]


se_condition <- sqrt(
  vcov(best_model)$cond[
    "conditionafter",
    "conditionafter"
  ]
)


rate_ratio <- exp(
  beta_condition
)


ci_rate_ratio <- exp(
  beta_condition +
    c(-1, 1) *
    1.96 *
    se_condition
)


percent_change <- (
  rate_ratio - 1
) * 100


ci_percent_change <- (
  ci_rate_ratio - 1
) * 100


cat(
  "\nCondition effect: after vs before\n",
  
  "Final model = ",
  BEST_MODEL_NAME,
  "\n",
  
  "LRT: chi-square(",
  df_condition,
  ") = ",
  round(
    chi_square_condition,
    3
  ),
  
  ", p = ",
  signif(
    p_lrt_condition,
    4
  ),
  
  "\nRate ratio = ",
  round(
    rate_ratio,
    3
  ),
  
  "\n95% CI = ",
  round(
    ci_rate_ratio[1],
    3
  ),
  
  " to ",
  round(
    ci_rate_ratio[2],
    3
  ),
  
  "\nPercent change = ",
  round(
    percent_change,
    1
  ),
  "%",
  
  "\n95% CI percent change = ",
  round(
    ci_percent_change[1],
    1
  ),
  
  "% to ",
  round(
    ci_percent_change[2],
    1
  ),
  "%\n",
  
  sep = ""
)


# ============================================================
# 7) ESTIMATED MARGINAL MEANS
# ============================================================

emm_latency <- emmeans(
  best_model,
  ~ condition,
  type = "response"
)


emm_latency_df <- as.data.frame(
  confint(
    emm_latency
  )
)


print(
  emm_latency_df
)


# Identify emmeans columns robustly

estimate_col <- intersect(
  c(
    "rate",
    "response",
    "emmean"
  ),
  names(emm_latency_df)
)[1]


lower_col <- grep(
  "lower.CL|asymp.LCL",
  names(emm_latency_df),
  value = TRUE
)[1]


upper_col <- grep(
  "upper.CL|asymp.UCL",
  names(emm_latency_df),
  value = TRUE
)[1]


emm_latency_df <- emm_latency_df %>%
  mutate(
    
    predicted_latency_min =
      .data[[estimate_col]] * 5,
    
    lower_latency_min =
      .data[[lower_col]] * 5,
    
    upper_latency_min =
      .data[[upper_col]] * 5
  )


print(
  emm_latency_df
)


# ============================================================
# 8) RESULT TABLE
# ============================================================

sleep_latency_results <- tibble(
  
  model = BEST_MODEL_NAME,
  
  comparison = "after vs before",
  
  rate_ratio = rate_ratio,
  
  CI_low = ci_rate_ratio[1],
  
  CI_high = ci_rate_ratio[2],
  
  percent_change = percent_change,
  
  percent_CI_low = ci_percent_change[1],
  
  percent_CI_high = ci_percent_change[2],
  
  chi_square = chi_square_condition,
  
  df = df_condition,
  
  p_LRT = p_lrt_condition
)


print(
  sleep_latency_results
)


# ============================================================
# 10) SLEEP LATENCY PLOT
# ============================================================


# ------------------------------------------------------------
# Colors
# ------------------------------------------------------------

condition_colors <- c(
  "before" = "#B79F00",
  "after"  = "#2F7F78"
)


# ------------------------------------------------------------
# Number of observations above displayed range
# ------------------------------------------------------------

n_not_shown <- sum(
  df_sleep_latency$sleep_latency_min > 120,
  na.rm = TRUE
)


not_shown_text <- paste0(
  n_not_shown,
  " observation",
  ifelse(
    n_not_shown == 1,
    "",
    "s"
  ),
  " >120 min\nnot shown"
)


# ------------------------------------------------------------
# LRT p-value for plot subtitle
# ------------------------------------------------------------

lrt_condition_df <- as.data.frame(
  lrt_condition
)


p_condition <- lrt_condition_df[
  2,
  "Pr(>Chisq)"
]


p_text <- ifelse(
  p_condition < 0.001,
  "p < 0.001",
  paste0(
    "p = ",
    formatC(
      p_condition,
      format = "f",
      digits = 3
    )
  )
)


# ------------------------------------------------------------
# Layout helper
# ------------------------------------------------------------

nudge_amt <- 0.25


# ------------------------------------------------------------
# Main plot
# ------------------------------------------------------------

plot_sleep_latency <- ggplot(
  df_sleep_latency,
  aes(
    x = condition,
    y = sleep_latency_min
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
    data = emm_latency_df,
    aes(
      x = condition,
      ymin = lower_latency_min,
      ymax = upper_latency_min
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
    data = emm_latency_df,
    aes(
      x = condition,
      ymin = lower_latency_min,
      ymax = upper_latency_min
    ),
    inherit.aes = FALSE,
    position = position_nudge(
      x = nudge_amt
    ),
    width = 0.14,
    linewidth = 1,
    color = "black"
  ) +
  
  # Model-estimated mean
  geom_point(
    data = emm_latency_df,
    aes(
      x = condition,
      y = predicted_latency_min,
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
  
  # Y-axis
  scale_y_continuous(
    breaks = seq(
      0,
      120,
      by = 20
    )
  ) +
  
  # Only limits visible range.
  # Values >120 min remain included in statistical analyses.
  coord_cartesian(
    ylim = c(
      0,
      120
    )
  ) +
  
  # Labels
  labs(
    x = NULL,
    
    y = "Sleep latency (min)",
    
    title = paste0(
      FRAGRANCE_NAME,
      " ",
      CONCENTRATION_FILE,
      ": sleep latency"
    ),
    
    subtitle = paste0(
      "Before vs. after: ",
      p_text,
      "   |   Rate ratio = ",
      sprintf(
        "%.2f",
        rate_ratio
      )
    )
  ) +
  
  # Theme
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


# ------------------------------------------------------------
# Add note ONLY if observations >120 min exist
# ------------------------------------------------------------

if (n_not_shown > 0) {
  
  plot_sleep_latency <- plot_sleep_latency +
    
    annotate(
      "text",
      x = 2.15,
      y = 105,
      label = not_shown_text,
      size = 3.2,
      color = "grey40",
      lineheight = 1.3,
      hjust = 0.5
    )
}


# Show plot
print(
  plot_sleep_latency
)



# ============================================================
# 11) FILE PREFIX
# ============================================================

file_prefix <- paste0(
  FRAGRANCE_NAME,
  "_",
  CONCENTRATION_FILE,
  "_sleep_latency_",
  BEST_MODEL_NAME
)



# ============================================================
# 12) SAVE PLOTS
# ============================================================


# ------------------------------------------------------------
# PNG
# ------------------------------------------------------------

ggsave(
  file.path(
    OUTPUT_DIR,
    paste0(
      file_prefix,
      "_results.png"
    )
  ),
  
  plot_sleep_latency,
  
  width = 7,
  height = 5.3,
  units = "in",
  dpi = 300,
  bg = "white"
)


# ------------------------------------------------------------
# PDF
# ------------------------------------------------------------

ggsave(
  file.path(
    OUTPUT_DIR,
    paste0(
      file_prefix,
      "_results.pdf"
    )
  ),
  
  plot_sleep_latency,
  
  width = 7,
  height = 5.3,
  units = "in",
  bg = "white"
)