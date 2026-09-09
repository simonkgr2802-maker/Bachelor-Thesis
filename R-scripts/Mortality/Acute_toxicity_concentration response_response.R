# Pipeline to create automated dose-response curves from toxicity data
# corresponding test system: fly inhalation mortality
# data is not processed or adapted before it is uploaded into the pipeline

rm(list=ls()) #remove ALL objects 
cat("\014") # clear console window prior to new run
Sys.setenv(LANG = "en") #Let's keep stuff in English
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
Sys.setlocale("LC_ALL")
options(scipen = 999)
options(dplyr.summarise.inform = FALSE)

packages<-c("drc", "tidyverse", "car", "ggplot2", "openxlsx", "multcomp") #packages needed in this script 
lapply(packages, library, character.only=TRUE)


#### ============================ Unit: applied volume [µL] ============================####

#### SETTINGS ####
# Only change this value when analysing another fragrance.
# The corresponding CSV file must have the same name, e.g.:
# FRAGRANCE <- "F1"  ->  F1.csv
# FRAGRANCE <- "F2"  ->  F2.csv
# FRAGRANCE <- "F3"  ->  F3.csv

FRAGRANCE <- "F1"
FILE <- paste0(FRAGRANCE, ".csv")


#### 1) import data set ####
# Expected columns:
# Konzentration_µL, BioRep, Replikat, Fliegen_total, Fliegen_tot

df_raw <- read.csv(
  FILE,
  sep = ";",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


#### 2) Clean and tidy data set ####

required_columns <- c(
  "Konzentration_µL",
  "BioRep",
  "Fliegen_total",
  "Fliegen_tot"
)

missing_columns <- setdiff(required_columns, names(df_raw))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Pool technical replicates A and B belonging to the same
# biological replicate and exposure level. Biological replicates
# remain separate observations.
df_prepared <- df_raw %>%
  group_by(Konzentration_µL, BioRep) %>%
  summarise(
    Fliegen_total = sum(Fliegen_total),
    Fliegen_tot = sum(Fliegen_tot),
    .groups = "drop"
  ) %>%
  mutate(
    Concentration = as.numeric(Konzentration_µL),
    Response = (Fliegen_tot / Fliegen_total) * 100,
    Time = 24,
    Substance = FRAGRANCE,
    Fly = "D. melanogaster"
  ) %>%
  select(
    Time,
    Substance,
    Fly,
    BioRep,
    Concentration,
    Response,
    Fliegen_total,
    Fliegen_tot
  )

# Keep the object name expected by the unchanged downstream pipeline.
df.fuels <- df_prepared

#### 3) Fit dose-response model ####

# Set ID to data according to substance, time and fly 
df_grouped <- df.fuels %>%
  group_by(Time, Substance, Fly) %>%
  mutate(ID = cur_group_id()) %>%
  ungroup()

ID.highest <- max(df_grouped$ID) #set highest ID 

# Fit models to data sets
model.results <- list() # Initialize results list

fit_best_model <- function(subset_i) {
  
  # Skip empty or invalid subsets
  if (nrow(subset_i) < 3 || length(unique(subset_i$Concentration)) < 3) {
    warning("Subset too small or insufficient unique concentrations.")
    return(list(best_model = NULL, best_name = NA, aics = NA))
  }
  
  # Try each model safely
  candidate_models <- list(
    weibull_2 = try(drm(Response ~ Concentration, data = subset_i, fct = W2.3()), silent= TRUE), # models fixed to 0 and 100% response caused problems therefore only LL.4 is fixed 
    loglogistic_2 = try(drm(Response ~ Concentration, data = subset_i, fct = LL.2()), silent = TRUE), # log logistic models with different nr of parameters 
    loglogistic_3 = try(drm(Response ~ Concentration, data = subset_i, fct = LL.3()), silent = TRUE),
    #loglogistic_4 = try(drm(Response ~ Concentration, data = subset_i, fct = LL.4(fixed = c(NA, 0, 100, NA))), silent = TRUE),  #4-parametric fit excluded as design of tests sets max. and min. response, see Daniels et al., 2025
    probit_2 = try(drm(Response ~ Concentration, data = subset_i, fct = LN.2 ()), silent = TRUE), # probit models especially suitable for binomial response (dead/alive)
    probit_3 = try(drm(Response ~ Concentration, data = subset_i, fct = LN.3 ()), silent = TRUE)
  )
  
  # Keep only successful fits
  candidate_models <- candidate_models[!sapply(candidate_models, inherits, "try-error")]
  
  if (length(candidate_models) == 0) {
    return(list(best_model = NULL, best_name = NA, aics = NA))
  }
  
  # Calculate AICs safely
  aic_values <- sapply(candidate_models, function(m) {
    tryCatch(AIC(m), error = function(e) NA)
  })
  
  aic_values <- aic_values[!is.na(aic_values)]
  if (length(aic_values) == 0) {
    return(list(best_model = NULL, best_name = NA, aics = NA))
  }
  
  # Select best model based on AIC
  best_name <- names(which.min(aic_values))
  best_model <- candidate_models[[best_name]]
  
  return(list(best_model = best_model, best_name = best_name, aics = aic_values))
}

# loop over all subsets 
results_summary <- data.frame(    #open empty data frame to save results for each ID
  ID = integer(),
  Time = numeric(),
  Substance = character(),
  Fly = character(),
  Best_Model = character(),
  Best_AIC = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:ID.highest) { #go through all possible scenarios (IDs) 
  subset_i <- df_grouped %>% filter(ID == i) #select the current set of data based on ID 
  
  model_fit <- fit_best_model(subset_i)
  
  model.results[[i]] <- list( #fill data frame with the information for each set of data 
    ID = i,
    Time = unique(subset_i$Time),
    Substance = unique(subset_i$Substance),
    Fly = unique(subset_i$Fly),
    Best_Model_Name = model_fit$best_name,
    Best_Model = model_fit$best_model,
    AICs = model_fit$aics
  )
  
  # add summary row
  results_summary <- rbind(results_summary, data.frame(
    ID = i,
    Time = unique(subset_i$Time),
    Substance = unique(subset_i$Substance),
    Fly = unique(subset_i$Fly),
    Best_Model = model_fit$best_name,
    AICs = model_fit$aics
  ))
  
}

# check results: which model fitted best? 
sapply(model.results, function(x) x$Best_Model_Name)


#### 4) Calculate LC values ####

for (i in 1:length(model.results)) {
  model <- model.results[[i]]$Best_Model #calculate only for chosen model 
  
  if (!is.null(model)) {
    LC50_full <- tryCatch({
      ED(model, 50, interval = "delta")  # EC50 + CI
    }, error = function(e) NA)
    
    if (!is.na(LC50_full[1,1])) {
      model.results[[i]]$LC50 <- LC50_full[1, "Estimate"]
      model.results[[i]]$LC50_lower <- LC50_full[1, "Lower"]
      model.results[[i]]$LC50_upper <- LC50_full[1, "Upper"]
    } else {
      model.results[[i]]$LC50 <- NA
      model.results[[i]]$LC50_lower <- NA
      model.results[[i]]$LC50_upper <- NA
    }
  } else {
    model.results[[i]]$LC50 <- NA
    model.results[[i]]$LC50_lower <- NA
    model.results[[i]]$LC50_upper <- NA
  }
  
  if (!is.null(model)) {
    LC10_full <- tryCatch({
      ED(model, 10, interval = "delta")  # EC10 + CI
    }, error = function(e) NA)
    
    if (!is.na(LC10_full[1,1])) {
      model.results[[i]]$LC10 <- LC10_full[1, "Estimate"]
      model.results[[i]]$LC10_lower <- LC10_full[1, "Lower"]
      model.results[[i]]$LC10_upper <- LC10_full[1, "Upper"]
    } else {
      model.results[[i]]$LC10 <- NA
      model.results[[i]]$LC10_lower <- NA
      model.results[[i]]$LC10_upper <- NA
    }
  } else {
    model.results[[i]]$LC10 <- NA
    model.results[[i]]$LC10_lower <- NA
    model.results[[i]]$LC10_upper <- NA
  }

  if (!is.null(model)) {
    LC5_full <- tryCatch({
      ED(model, 5, interval = "delta")  # EC5 + CI
    }, error = function(e) NA)

    if (!is.na(LC5_full[1,1])) {
      model.results[[i]]$LC5 <- LC5_full[1, "Estimate"]
      model.results[[i]]$LC5_lower <- LC5_full[1, "Lower"]
      model.results[[i]]$LC5_upper <- LC5_full[1, "Upper"]
    } else {
      model.results[[i]]$LC5 <- NA
      model.results[[i]]$LC5_lower <- NA
      model.results[[i]]$LC5_upper <- NA
    }
  } else {
    model.results[[i]]$LC5 <- NA
    model.results[[i]]$LC5_lower <- NA
    model.results[[i]]$LC5_upper <- NA
  }

}
# table with all results
LC50_table <- do.call(rbind, lapply(model.results, function(x) {
  data.frame(
    Time = x$Time,
    Substance = x$Substance,
    Fly = x$Fly,
    Best_Model = x$Best_Model_Name,
    LC50_uL = x$LC50,
    LC50_upper_uL = x$LC50_upper,
    LC50_lower_uL = x$LC50_lower
  )
}))

LC10_table <- do.call(rbind, lapply(model.results, function(x) {
  data.frame(
    Time = x$Time,
    Substance = x$Substance,
    Fly = x$Fly,
    Best_Model = x$Best_Model_Name,
    LC10_uL = x$LC10
  )
}))


LC5_table <- do.call(rbind, lapply(model.results, function(x) {
  data.frame(
    Time = x$Time,
    Substance = x$Substance,
    Fly = x$Fly,
    Best_Model = x$Best_Model_Name,
    LC5_uL = x$LC5,
    LC5_upper_uL = x$LC5_upper,
    LC5_lower_uL = x$LC5_lower
  )
}))



#### 5) Calculate Confidence Interval (CI) bands ####

# open empty list for results
CI.results <- list()

for (i in 1:length(model.results)) {    # whole CI calculation unnecessary, is done again in the plot code
  
  result <- model.results[[i]]
  model <- result$Best_Model
  
  # skip if model failed
  if (is.null(model) || inherits(model, "try-error")) next
  
  # get the matching dataset
  subset_i <- df_grouped %>% filter(ID == i)
  subset_i <- select(subset_i, Concentration)
  subset_i <- as.data.frame (subset_i)
  
  # get predictions + SE
  pred <- predict(model, subset_i, interval = "confidence") #calculates Prediction + lower and upper 95% CI 
  
}



#### =========================================================
#### 6) Plot all models — applied volume [µL]
#### =========================================================

# IMPORTANT:
# The complete analysis remains on the original applied-volume scale.
# "Concentration" contains the applied fragrance volume in µL.
# No conversion to µL/L is performed anywhere in this section.

plots <- list()


for (i in seq_along(model.results)) {

  result <- model.results[[i]]
  model <- result$Best_Model

  if (is.null(model)) {
    message(
      paste(
        "No model for ID",
        i
      )
    )
    next
  }


  # ==========================================================
  # Get LC values
  # ==========================================================
  #
  # ED() returns values on the same x-scale used to fit the model.
  # Since the model was fitted with Concentration in µL,
  # LC5, LC10 and LC50 are already in applied volume [µL].
  # ==========================================================

  LC50 <- result$LC50
  LC10 <- result$LC10
  LC5  <- result$LC5


  # ==========================================================
  # Data for this ID
  # ==========================================================

  data_subset <- df_grouped %>%
    dplyr::filter(
      ID == i
    )


  # Keep positive applied volumes for plotting
  data_subset_pos <- data_subset %>%
    dplyr::filter(
      Concentration > 0
    )


  if (nrow(data_subset_pos) == 0) {
    message(
      paste(
        "No positive applied volumes for ID",
        i
      )
    )
    next
  }


  # ==========================================================
  # Plot range
  # ==========================================================

  x_min_plot <- 0
  x_max_plot <- 95


  # ==========================================================
  # Model predictions
  # ==========================================================

  pred_data <- data.frame(
    Concentration = seq(
      min(data_subset_pos$Concentration) * 0.80,
      x_max_plot,
      length.out = 300
    )
  )


  pred_data$Fitted <- suppressWarnings(
    predict(
      model,
      newdata = pred_data,
      interval = "confidence"
    )
  )


  # ==========================================================
  # Legend
  # ==========================================================

  model_label <- paste0(
    "Selected model: ",
    result$Best_Model_Name
  )


  # ==========================================================
  # Base plot
  # ==========================================================

  p <- ggplot(
    data_subset_pos,
    aes(
      x = Concentration
    )
  ) +

    # 95% confidence interval
    geom_ribbon(
      data = pred_data,
      aes(
        x = Concentration,
        ymin = Fitted[, "Lower"],
        ymax = Fitted[, "Upper"],
        fill = "95% CI"
      ),
      inherit.aes = FALSE,
      alpha = 0.4
    ) +

    # Biological replicates
    geom_point(
      data = data_subset_pos,
      aes(
        x = Concentration,
        y = Response,
        shape = "Biological replicates"
      ),
      color = "black",
      size = 2
    ) +

    # Selected concentration-response model
    geom_line(
      data = pred_data,
      aes(
        x = Concentration,
        y = Fitted[, "Prediction"],
        color = model_label
      ),
      inherit.aes = FALSE,
      linewidth = 1
    ) +

    # Linear x-axis in applied volume [µL]
    scale_x_continuous(
      limits = c(
        x_min_plot,
        x_max_plot
      ),
      breaks = scales::breaks_pretty(
        n = 7
      )
    ) +

    coord_cartesian(
      ylim = c(
        0,
        105
      )
    ) +

    scale_y_continuous(
      breaks = c(
        0,
        25,
        50,
        75,
        100
      )
    ) +

    scale_fill_manual(
      values = c(
        "95% CI" = "grey80"
      )
    ) +

    scale_shape_manual(
      values = c(
        "Biological replicates" = 16
      )
    ) +

    scale_color_manual(
      values = setNames(
        "#0066CC",
        model_label
      )
    )


  # ==========================================================
  # LC lines and labels
  # ==========================================================

  extra_layers <- list()


  if (!is.null(LC50) && !is.na(LC50)) {
    extra_layers <- append(
      extra_layers,
      list(
        geom_vline(
          xintercept = LC50,
          linetype = "dashed",
          color = "#CC0000"
        ),
        annotate(
          "text",
          x = LC50 * 0.97,
          y = 94,
          label = paste0(
            "LC50 = ",
            round(LC50, 2),
            " µL"
          ),
          color = "#CC0000",
          hjust = 1,
          size = 4
        )
      )
    )
  }


  if (!is.null(LC10) && !is.na(LC10)) {
    extra_layers <- append(
      extra_layers,
      list(
        geom_vline(
          xintercept = LC10,
          linetype = "dashed",
          color = "#FF8C00"
        ),
        annotate(
          "text",
          x = LC10 * 0.97,
          y = 78,
          label = paste0(
            "LC10 = ",
            round(LC10, 2),
            " µL"
          ),
          color = "#FF8C00",
          hjust = 1,
          size = 4
        )
      )
    )
  }


  if (!is.null(LC5) && !is.na(LC5)) {
    extra_layers <- append(
      extra_layers,
      list(
        geom_vline(
          xintercept = LC5,
          linetype = "dashed",
          color = "#009900"
        ),
        annotate(
          "text",
          x = LC5 * 0.97,
          y = 64,
          label = paste0(
            "LC5 = ",
            round(LC5, 2),
            " µL"
          ),
          color = "#009900",
          hjust = 1,
          size = 4
        )
      )
    )
  }


  # ==========================================================
  # Finalize plot
  # ==========================================================

  p <- p +
    extra_layers +

    labs(
      title = paste(
        "Substance:",
        result$Substance,
        "| Model:",
        result$Best_Model_Name
      ),

      x = "Concentration [µL]",
      y = "Mortality [%]",

      color = NULL,
      fill = NULL,
      shape = NULL
    ) +

    theme_bw() +

    theme(
      axis.text = element_text(
        size = 12
      ),

      axis.title = element_text(
        size = 14,
        face = "bold"
      ),

      plot.title = element_text(
        size = 14
      ),

      legend.position = "top",
      legend.direction = "horizontal",

      legend.text = element_text(
        size = 9
      ),

      legend.key.width = unit(
        1.2,
        "cm"
      ),

      legend.spacing.x = unit(
        0.15,
        "cm"
      )
    )


  plots[[i]] <- p

  # Explicit printing is needed when sourcing the script
  # or plotting inside a loop.
  print(p)
}


plots
