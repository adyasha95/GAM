# =============================================================================
# Survival Analysis
# Author: Adyasha Tejaswi Khuntia
# Date:   20 Sep 2025
# -----------------------------------------------------------------------------
# What this script does
# 1) Reads cdf.xlsx and makes status monotone per subject (cummax, NA->0)
# 2) Builds two endpoints:
#      - MCI/Dementia (event: first status >= 2)
#      - Dementia-only (event: first status >= 3)
#    -> subject-level: time_final = event_time (if any) else last_followup
#                     status_final = 1 if event else 0
# 3) Cox PH with class as STRATA; exports HRs, summaries, and PH tests
# 4) Per-class parametric survival (one flexsurv model per class), fixed DIST
#    -> anchors: mean(age_at_base, yr_educ), mode(gender, ap4/factor) or mean(ap4/numeric)
#    -> solid curves + 95% CI ribbons; dashed KM overlays for QC
# 5) Risk tables (numbers at risk) at 0,2,4,6,8,10 (count time_final >= tick)
# 6) Plots: Panel A (MCI/Dementia), Panel C (Dementia), Panel B/D risk tables
# 7) Combined 2×2 figure; all outputs saved with specific filenames
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(survival)
  library(survminer)
  library(flexsurv)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(purrr)
})

set.seed(123)

# ------------------------- I/O & styling -------------------------------------
setwd('/Users/adyashatejaswikhuntia/work/personal/job/upwork/CANVAdatascience/survival_analysis/')
infile   <- "valle+APS2cut.csv"
out_dir  <- "valle+APS2cut"
out_figs <- file.path(out_dir, "figs");   dir.create(out_figs, TRUE, TRUE)
out_tabs <- file.path(out_dir, "tables"); dir.create(out_tabs, TRUE, TRUE)

CLASS_COLS <- c(LikelyPos="#418979", Neg="#556219", Pos="#e5a335")
x_breaks   <- seq(0, 12, 2)
risk_ticks <- seq(0, 10, 2)

# Choose ONE distribution for all classes (default per spec)
# Valid: "llogis","lnorm","weibull","gompertz","gengamma","genf","gamma"
DIST <- "weibull"
add_km_dashed <- TRUE    # toggle KM dashed overlays for QC (default ON)

# ------------------------- Load & monotone status ----------------------------
#df_raw <- read_excel(infile) |> as.data.frame()
df_raw <- readr::read_csv(infile) |> as.data.frame()
# Keep your column names; only coerce types where necessary (no renaming)
df_raw <- df_raw |>
  dplyr::mutate(
    patient       = as.character(patient),
    class         = as.factor(class),
    yrs_from_base = suppressWarnings(as.numeric(yrs_from_base)),
    status        = suppressWarnings(as.numeric(status)),
    gender        = if (is.factor(gender) || is.character(gender)) factor(gender) else gender,
    ap4           = if (length(unique(na.omit(ap4))) <= 3) factor(ap4) else suppressWarnings(as.numeric(ap4))
  ) |>
  dplyr::filter(!is.na(class)) |>
  droplevels()

# ------------------------- Helpers -------------------------------------------
Mode <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA); ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

# Collapse to subject-level endpoint at threshold (>= thresh)
make_endpoint <- function(df, thresh) {
  df |>
    arrange(patient, yrs_from_base) |>
    group_by(patient) |>
    summarise(
      event_time = {
        idx <- which(status >= thresh)
        if (length(idx)) yrs_from_base[min(idx)] else NA_real_
      },
      has_event  = any(status >= thresh, na.rm = TRUE),
      last_time  = max(yrs_from_base, na.rm = TRUE),
      .groups="drop_last"
    ) |>
    ungroup() |>
    mutate(
      time_final   = ifelse(has_event, event_time, last_time),
      status_final = as.integer(has_event)
    ) |>
    left_join(
      df |>
        arrange(patient, yrs_from_base) |>
        group_by(patient) |>
        summarise(
          class        = class[which(!is.na(class))][1],
          age_at_base  = age_at_base[which(!is.na(age_at_base))][1],
          gender       = gender[which(!is.na(gender))][1],
          yr_educ      = yr_educ[which(!is.na(yr_educ))][1],
          ap4          = ap4[which(!is.na(ap4))][1],
          .groups="drop"
        ),
      by="patient"
    ) |>
    filter(!is.na(class))
}

# Risk counts at fixed ticks (n at risk defined as time_final >= tick)
risk_counts <- function(df) {
  map_dfr(risk_ticks, function(tt) {
    df |>
      group_by(class) |>
      summarise(n = sum(time_final >= tt, na.rm = TRUE), .groups="drop") |>
      mutate(time = tt)
  }) |>
    mutate(class = factor(class, levels = levels(df$class))) |>
    arrange(class, time)
}

# Risk table as a ggplot grid (numbers in cells)
# Replace your existing risk_grid_plot() with this version
risk_grid_plot <- function(riskdf,
                           class_levels = c("LikelyPos", "Neg", "Pos"),
                           risk_ticks   = NULL,                 # <- no self-reference
                           class_cols   = CLASS_COLS,
                           row_spacing  = 0.60,
                           number_size  = 2.8,
                           base_size    = 11) {
  
  # default ticks if not supplied
  if (is.null(risk_ticks)) risk_ticks <- seq(0, 10, 2)
  
  # enforce order (Low top … High bottom)
  riskdf <- riskdf %>%
    dplyr::mutate(class = factor(class, levels = class_levels))
  
  # y positions (squeezed)
  y_pos    <- function(cls) row_spacing * as.numeric(factor(cls, levels = class_levels))
  y_breaks <- row_spacing * seq_along(class_levels)
  y_labels <- class_levels
  
  ggplot(riskdf, aes(x = time, y = y_pos(class))) +
    # white halo then colored numbers
    geom_text(aes(label = n), size = number_size + 0.3, color = "white") +
    geom_text(aes(label = n, color = class), size = number_size) +
    coord_cartesian(clip = "off") +
    scale_color_manual(values = class_cols, drop = FALSE) +
    scale_x_continuous(breaks = risk_ticks,
                       limits = c(-1, 11),                    # <- extra left room
                       expand = expansion(add = c(0, 0))) +
    scale_y_continuous(breaks = y_breaks, labels = y_labels) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      axis.title.y    = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.line.x     = element_blank(),        # <- remove x axis line
      axis.line.y     = element_blank(),
      axis.text.y     = element_text(color = "black", face = "bold"),
      plot.margin     = margin(0, 8, 6, 10)
    )
}
# Fit one class with fixed distribution; z-center numerics for stability (internal columns)
safe_flexsurv <- function(formula, data, dists) {
  # try a sequence of distributions; return the first that converges
  for (d in dists) {
    fit <- try(
      flexsurv::flexsurvreg(formula, data = data, dist = d,
                            control = list(maxit = 2000)),
      silent = TRUE
    )
    if (!inherits(fit, "try-error")) return(list(fit = fit, used = d))
  }
  stop("All distributions failed: ", paste(dists, collapse = ", "))
}

fit_one_class <- function(dcl, dist) {
  # keep only rows usable for modeling
  d <- dcl |>
    dplyr::filter(is.finite(time_final), !is.na(status_final)) |>
    dplyr::mutate(
      time_final = pmax(time_final, .Machine$double.eps)
    )
  
  # quick exit: no rows
  if (nrow(d) < 2) {
    return(list(fit = NULL, used_dist = NA_character_,
                mu_age = NA_real_, sd_age = NA_real_,
                mu_ed = NA_real_,  sd_ed = NA_real_,
                zero_events = TRUE))
  }
  
  # center/scale continuous for stability
  mu_age <- mean(d$age_at_base, na.rm = TRUE); sd_age <- sd(d$age_at_base, na.rm = TRUE)
  mu_ed  <- mean(d$yr_educ,     na.rm = TRUE); sd_ed  <- sd(d$yr_educ,     na.rm = TRUE)
  
  d$age_c  <- if (is.finite(sd_age) && sd_age > 0) (d$age_at_base - mu_age)/sd_age else 0
  d$educ_c <- if (is.finite(sd_ed)  && sd_ed  > 0) (d$yr_educ     - mu_ed )/sd_ed  else 0
  
  # drop rows with NA in covariates used
  d <- d |> dplyr::filter(!is.na(age_c), !is.na(educ_c))
  
  # factors may be single-level within class → drop from formula to avoid singularities
  use_gender <- is.factor(d$gender) && nlevels(droplevels(d$gender)) > 1
  use_ap4    <- (is.factor(d$ap4)   && nlevels(droplevels(d$ap4))   > 1) ||
    (is.numeric(d$ap4)  && sd(d$ap4, na.rm = TRUE) > 0)
  
  # build formula dynamically
  rhs <- c("age_c", "educ_c", if (use_gender) "gender" else NULL, if (use_ap4) "ap4" else NULL)
  fml <- reformulate(rhs, response = "survival::Surv(time_final, status_final)")
  
  # if ZERO events in this class, skip parametric fit and signal fallback
  zero_events <- sum(d$status_final == 1, na.rm = TRUE) == 0
  if (zero_events) {
    return(list(fit = NULL, used_dist = NA_character_,
                mu_age = mu_age, sd_age = sd_age,
                mu_ed = mu_ed,   sd_ed = sd_ed,
                zero_events = TRUE,
                data = d))
  }
  
  # try main dist, then lnorm → weibull → exponential
  tried <- c(dist, "lnorm", "weibull", "exponential")
  tried <- unique(tried)
  res <- safe_flexsurv(fml, d, tried)
  
  list(fit = res$fit, used_dist = res$used,
       mu_age = mu_age, sd_age = sd_age, mu_ed = mu_ed, sd_ed = sd_ed,
       zero_events = FALSE, data = d)
}

# Predict per class to observed horizon; anchors: mean(age, educ), mode(gender, ap4/factor) or mean(ap4/numeric)
parametric_by_class <- function(df_surv, dist, tmax_global = 12, by = 0.1) {
  cl_levels <- levels(df_surv$class)
  preds <- lapply(cl_levels, function(cl) {
    dcl <- df_surv[df_surv$class == cl, , drop = FALSE]
    if (!nrow(dcl)) return(NULL)
    
    # observed horizon
    tmax_cl <- min(tmax_global, max(dcl$time_final, na.rm = TRUE))
    if (!is.finite(tmax_cl)) tmax_cl <- tmax_global
    times <- seq(0, tmax_cl, by = by)
    
    pack <- fit_one_class(dcl, dist)
    
    # anchors within class
    age_ref <- mean(dcl$age_at_base, na.rm = TRUE)
    ed_ref  <- mean(dcl$yr_educ,     na.rm = TRUE)
    g_ref   <- if (is.factor(dcl$gender)) {
      x <- dcl$gender; x <- x[!is.na(x)]; if (length(x)) names(sort(table(x), dec=TRUE))[1] else levels(dcl$gender)[1]
    } else mean(dcl$gender, na.rm = TRUE)
    ap4_ref <- if (is.factor(dcl$ap4)) {
      x <- dcl$ap4; x <- x[!is.na(x)]; if (length(x)) names(sort(table(x), dec=TRUE))[1] else levels(dcl$ap4)[1]
    } else mean(dcl$ap4, na.rm = TRUE)
    
    # zero-events: return flat S=1 with tiny ribbon so plot renders
    if (isTRUE(pack$zero_events) || is.null(pack$fit)) {
      return(tibble::tibble(
        t = times,
        surv  = rep(1, length(times)),
        lower = pmax(0, 1 - 1e-6),
        upper = pmin(1, 1 + 1e-6),
        class = factor(cl, levels = cl_levels)
      ))
    }
    
    # build newdata (respect factor levels)
    age_c  <- if (is.finite(pack$sd_age) && pack$sd_age > 0) (age_ref - pack$mu_age)/pack$sd_age else 0
    educ_c <- if (is.finite(pack$sd_ed)  && pack$sd_ed  > 0) (ed_ref  - pack$mu_ed )/pack$sd_ed  else 0
    nd <- data.frame(
      age_c  = age_c,
      educ_c = educ_c
    )
    if ("gender" %in% names(model.frame(pack$fit))) {
      nd$gender <- if (is.factor(dcl$gender)) factor(g_ref, levels = levels(dcl$gender)) else g_ref
    }
    if ("ap4" %in% names(model.frame(pack$fit))) {
      nd$ap4 <- if (is.factor(dcl$ap4)) factor(ap4_ref, levels = levels(dcl$ap4)) else as.numeric(ap4_ref)
    }
    
    sm  <- summary(pack$fit, newdata = nd, type = "survival", t = times, ci = TRUE)
    dfp <- as.data.frame(sm)
    names(dfp)[1:4] <- c("time","est","lcl","ucl")
    tibble::tibble(
      t     = as.numeric(dfp$time),
      surv  = pmin(pmax(as.numeric(dfp$est), 0), 1),
      lower = pmin(pmax(as.numeric(dfp$lcl), 0), 1),
      upper = pmin(pmax(as.numeric(dfp$ucl), 0), 1),
      class = factor(cl, levels = cl_levels)
    )
  }) |> dplyr::bind_rows()
  
  if (nrow(preds)) preds$class <- factor(preds$class, levels = cl_levels)
  preds
}

# Plot parametric solid + ribbons; optional dashed KM overlays
plot_param <- function(preds, df_surv, title_lab, add_km = TRUE) {
  p <- ggplot(preds, aes(x = t, y = surv, color = class, fill = class)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20, linewidth = 0) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = CLASS_COLS, drop = FALSE) +
    scale_fill_manual(values = CLASS_COLS, drop = FALSE) +
    scale_x_continuous(limits = c(0,12), breaks = x_breaks,
                       expand = expansion(mult = c(0,0))) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
    labs(x = "Time from Baseline (years)", y = title_lab,
         color = NULL, fill = NULL) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "bottom",
      axis.line = element_line(color = "black"),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.text.x  = element_text(face = "plain", size = 12),
      axis.text.y  = element_text(face = "plain", size = 12)
    )
  
  if (isTRUE(add_km)) {
    fk <- survfit(Surv(time_final, status_final) ~ class, data = df_surv)
    km <- survminer::surv_summary(fk, data = df_surv) |>
      mutate(class = factor(gsub("^class=", "", strata),
                            levels = levels(df_surv$class))) |>
      filter(time <= 12)
    
    p <- p + geom_step(
      data = km,
      aes(x = time, y = surv, color = class, group = class),
      inherit.aes = FALSE,
      linetype = "dashed", linewidth = 0.9
    )
  }
  p
}

# ------------------------- Build endpoints -----------------------------------
df_mci_dem <- make_endpoint(df_raw, 2)   # Endpoint 1: MCI/Dementia
df_dem     <- make_endpoint(df_raw, 3)   # Endpoint 2: Dementia-only

# Save subject-level datasets
write_csv(df_mci_dem, file.path(out_tabs, "subject_level_MCIplusDementia.csv"))
write_csv(df_dem,     file.path(out_tabs, "subject_level_Dementia.csv"))

# ------------------------- Cox PH (strata = class) ---------------------------
cox_formula <- as.formula(Surv(time_final, status_final) ~ strata(class) + age_at_base + gender + yr_educ + ap4)

# Put this right after you define cox_formula (so it can see that object),
# and before you call cox_export(...).

cox_export <- function(df, label, out_tabs) {
  # 1) Harden the input to a plain base data.frame (not tibble / grouped)
  df <- df %>%
    dplyr::ungroup() %>%
    base::as.data.frame(stringsAsFactors = FALSE)
  
  # 2) Coerce types safely (no renaming of your columns)
  if (!is.factor(df$class))  df$class  <- factor(df$class, levels = c("LikelyPos", "Neg", "Pos"))
  if (is.character(df$gender)) df$gender <- factor(df$gender)
  # ap4 stays numeric if numeric; if character turn to factor
  if (is.character(df$ap4)) df$ap4 <- factor(df$ap4)
  
  # 3) Fit inside the data's environment to bypass model.frame(data=...) issues
  fit <- with(
    df,
    survival::coxph(
      Surv(time_final, status_final) ~ strata(class) + age_at_base + gender + yr_educ + ap4,
      ties = "efron", x = TRUE, model = TRUE
    )
  )
  
  # 4) Export HRs (covariates only; class is stratified)
  est <- stats::coef(fit)
  if (length(est)) {
    se  <- sqrt(diag(stats::vcov(fit)))
    hr_tbl <- tibble::tibble(
      term = names(est),
      HR   = exp(est),
      LCL  = exp(est - 1.96*se),
      UCL  = exp(est + 1.96*se)
    )
  } else {
    hr_tbl <- tibble::tibble(term = character(), HR = numeric(), LCL = numeric(), UCL = numeric())
  }
  readr::write_csv(hr_tbl, file.path(out_tabs, paste0("coxph_HR_covariates_", label, ".csv")))
  
  # 5) Save full summary + PH test
  sink(file.path(out_tabs, paste0("coxph_summary_", label, ".txt"))); print(summary(fit)); sink()
  sink(file.path(out_tabs, paste0("coxph_PH_",      label, ".txt"))); print(cox.zph(fit));  sink()
  
  return(fit)
}

cox_mci <- cox_export(df_mci_dem, "MCIplusDementia", out_tabs)
cox_dem <- cox_export(df_dem,     "Dementia",        out_tabs)

# ------------------------- Parametric predictions per class -------------------
pred_mci <- parametric_by_class(df_mci_dem, DIST, tmax_global = 12, by = 0.1)
pred_dem <- parametric_by_class(df_dem,     DIST, tmax_global = 12, by = 0.1)

# Save parametric prediction grids
write_csv(pred_mci, file.path(out_tabs, paste0("param_preds_MCIplusDementia_", DIST, ".csv")))
write_csv(pred_dem, file.path(out_tabs, paste0("param_preds_Dementia_", DIST, ".csv")))

# ------------------------- Risk tables ---------------------------------------
risk_mci <- risk_counts(df_mci_dem)
risk_dem <- risk_counts(df_dem)
# Also save wide versions (columns 0..10) 
write_csv(tidyr::pivot_wider(risk_mci, names_from = time, values_from = n),
          file.path(out_tabs, "risk_table_MCIplusDementia_0to10_by2.csv"))
write_csv(tidyr::pivot_wider(risk_dem, names_from = time, values_from = n),
          file.path(out_tabs, "risk_table_Dementia_0to10_by2.csv"))

# ------------------------- Plots (A, C) --------------------------------------
plot_A <- plot_param(pred_mci, df_mci_dem, "Survival Probability (MCI/Dementia)", add_km = add_km_dashed)
plot_C <- plot_param(pred_dem, df_dem,     "Survival Probability (Dementia)",      add_km = add_km_dashed)

# ------------------------- Risk table plots (B, D) ---------------------------
plot_B <- risk_grid_plot(risk_mci)
plot_D <- risk_grid_plot(risk_dem)

# Save panels
ggsave(file.path(out_figs, paste0("panel_A_plot_MCIplusDEM_", DIST, ".png")), plot_A, width = 7, height = 5, dpi = 300)
ggsave(file.path(out_figs, paste0("panel_C_plot_DEM_",       DIST, ".png")), plot_C, width = 7, height = 5, dpi = 300)
ggsave(file.path(out_figs, "panel_B_risktable_MCIplusDEM.png"), plot_B, width = 7, height = 2.2, dpi = 300)
ggsave(file.path(out_figs, "panel_D_risktable_DEM.png"),        plot_D, width = 7, height = 2.2, dpi = 300)

# ------------------------- Combined 2×2 figure -------------------------------
combined <- ((plot_A | plot_C) / (plot_B | plot_D)) + plot_layout(heights = c(3, 1))
print(combined)
ggsave(file.path(out_figs, paste0("KM_Parametric_combined_AD_", DIST, ".png")),
       combined, width = 12, height = 8, dpi = 300)
ggsave(file.path(out_figs, paste0("KM_Parametric_combined_AD_", DIST, ".pdf")),
       combined, width = 12, height = 8)
