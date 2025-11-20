# =============================================================================
# GAM Trajectories: Fit, Summaries, Predictions, and Plots
# =============================================================================
# Author: Adyasha Khuntia
# Date:   18 Sep 2025
#
# What this script does:
#   1) Pre-process data (no renaming of your columns; only adds `time` & `m_score_std`)
#   2) Fit 4 GAMs (mPACC=zscore, MMSE=mlm (binomial logit scaled to 0–30), MEM=m_score_std, EXEC=func_score)
#   3) Export model summaries (parametric mains/interactions + smooths) to CSV and XLSX
#   4) Export prediction data (new_data grid with fit + 95% CI) per response
#   5) Export 4 individual plots + combined 2×2 panel with fixed axes matching the original figure
#
# Notes:
#   - Modeling choices preserved exactly (bs='cr', k=3; your families & predict types).
#   - We use coord_cartesian for y-limits to avoid ribbon "fill length 0" errors.
#   - “Population-level predictions”: we exclude random effects in predict() to show
#     the typical trend by class at reference covariates (not a specific patient).
# =============================================================================

# ---- Libraries --------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
  library(patchwork)
  library(writexl)    # for XLSX output
  library(readr)      # for CSV output
})

# ---- Output folders ---------------------------------------------------------
out_dir  <- "gam_outputs"
fig_dir  <- file.path(out_dir, "figs")
pred_dir <- file.path(out_dir, "predictions")
sum_dir  <- file.path(out_dir, "summaries")
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(pred_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(sum_dir,  showWarnings = FALSE, recursive = TRUE)

# =========================
# 1) PRE-PROCESSING
# =========================

# --- Load data (keep original column names) ---
df <- read_excel("cdf.xlsx") %>% as.data.frame()

# Time as numeric (adds column; does not alter your originals).
# mgcv smooths must operate on numeric time; factors/characters will break abs(), etc.
df$time <- as.numeric(df$yrs_from_base)

# Convert to appropriate types for modeling.
df$patient <- factor(df$patient)   # random effect ID
df$class   <- factor(df$class, levels = c("Low","Intermediate","Elevated","High"))  # ensures palette/legend order
df$gender  <- factor(df$gender)

# ap4: treat as factor if small discrete set (e.g., 0/1/2), else numeric (continuous dosage).
if (length(unique(na.omit(df$ap4))) <= 3) {
  df$ap4 <- factor(df$ap4)
} else {
  df$ap4 <- as.numeric(df$ap4)
}

# --- MEM standardization at baseline (per-patient min time) ------------------
# Why: The original MEM panel is displayed as relative change around baseline.
# We standardize MEM so that baseline ~0 and trajectories are negative afterwards.
# Baseline definition is robust: each patient's minimum time (not strictly time==0).
baseline_rows <- df %>%
  group_by(patient) %>%
  mutate(tmin = min(time, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(time == tmin & !is.na(m_score)) %>%
  select(m_score)

mem_mu <- mean(baseline_rows$m_score, na.rm = TRUE)
mem_sd <- sd(baseline_rows$m_score,  na.rm = TRUE)
df$m_score_std <- (df$m_score - mem_mu) / mem_sd

# --- Config & helpers --------------------------------------------------------
RESPONSES  <- c("zscore","mlm","m_score_std","func_score")
CLASS_COLS <- c(Low="#418979", Intermediate="#556219", Elevated="#e5a335", High="#96410e")
YLABS      <- c(zscore="mPACC", mlm="MMSE", m_score_std="MEM", func_score="EXEC")

# Statistical mode (for reference factor levels)
Mode <- function(x) { ux <- na.omit(unique(x)); if (!length(ux)) return(NA); ux[which.max(tabulate(match(x, ux)))] }

MAX_X      <- 12                               # x-axis ends at 12 (original style)
N_GRID     <- 121
time_seq   <- seq(0, MAX_X, length.out = N_GRID)
max_follow <- max(df$time, na.rm = TRUE)       # curves truncate here (~11–11.5)

# Reference covariates (population-level predictions)
# We fix non-focal covariates to typical values to isolate class differences.
gender_ref <- Mode(df$gender)
ap4_ref    <- Mode(df$ap4)
age_ref    <- mean(df$age_at_base, na.rm=TRUE)
edu_ref    <- mean(df$yr_educ,     na.rm=TRUE)
pat_ref    <- Mode(df$patient)  # used only to keep factor levels valid; RE is excluded at predict-time

# =========================
# 2) MODEL + PREDICT (FUNCTIONS)
# =========================

# Fit ONE GAM and return model + prediction tibble (no plotting here)
fit_one <- function(resp) {
  # ----------------------------
  # Family & LHS (kept as in your code)
  # ----------------------------
  # For MMSE (mlm), we model as binomial with logit link on 30 trials:
  #   LHS is cbind(mlm, 30 - mlm)
  #   predict(..., type="link") returns η (log-odds). We transform via invlink(η) to 0..30 scale.
  # For Gaussian outcomes, we use identity link on the response scale and predict(type="response").
  if (resp == "mlm") {
    fam <- binomial("logit")
    lhs <- sprintf("cbind(%s, 30 - %s)", resp, resp)
    pred_type <- "link"   # predict on the link-scale (η) to compute SEs correctly
    # invlink: transforms linear predictor η back to the response scale.
    # For binomial-logit: inverse link is logistic: p = 1/(1+exp(-η)), scaled to 0..30 by multiplying by 30.
    invlink   <- function(eta) 1/(1+exp(-eta)) * 30
  } else {
    fam <- gaussian()
    lhs <- resp
    pred_type <- "response"  # direct predictions on the response scale for Gaussian
    # For Gaussian-identity, the inverse link is the identity function (returns the same value).
    invlink   <- identity
  }
  
  # ----------------------------
  # Model formula (kept: bs='cr', k=3)
  # ----------------------------
  # - Main effects for class, age_at_base, gender, ap4, yr_educ
  # - Interactions with time (age_at_base:time, gender:time, ap4:time, yr_educ:time)
  #   let slopes vary by baseline covariates
  # - s(time) and s(time, by=class): class-specific temporal smooths
  # - s(patient, bs="re"): random intercept per patient
  FORM <- as.formula(paste0(
    lhs, " ~ ",
    "class + age_at_base + gender + ap4 + yr_educ + ",
    "age_at_base:time + gender:time + ap4:time + yr_educ:time + ",
    "s(time, bs='cr', k=3) + s(time, by=class, bs='cr', k=3) + ",
    "s(patient, bs='re')"
  ))
  
  # Fit the model (drop rows containing NA in used variables)
  need <- all.vars(FORM)
  dfit <- df %>% tidyr::drop_na(any_of(need))
  m <- mgcv::bam(FORM, data=dfit, family=fam, method="fREML",
                 discrete=TRUE, select=TRUE, control=list(maxit=600))
  
  # ----------------------------
  # Prediction grid (population-level)
  # ----------------------------
  # We hold covariates at reference values to depict a "typical" participant.
  # We include 'patient' with a single reference level to keep factor structure,
  # but we EXCLUDE random-effect terms in predict(), so predictions are marginal.
  newdat <- tidyr::expand_grid(
    time  = time_seq,
    class = levels(df$class)
  ) %>%
    mutate(
      age_at_base = age_ref,
      yr_educ     = edu_ref,
      gender      = factor(gender_ref, levels = levels(df$gender)),
      ap4         = if (is.factor(df$ap4)) factor(ap4_ref, levels = levels(df$ap4)) else as.numeric(ap4_ref),
      patient     = factor(pat_ref, levels = levels(df$patient))
    )
  
  # Exclude random effects (population-level predictions).
  # mgcv assigns names like "s(patient)" and possibly "s(time,patient)".
  s_names <- rownames(summary(m)$s.table)
  excl <- s_names[grepl("^s\\(patient", s_names)]
  if (length(excl) == 0) excl <- c("s(patient)", "s(patient,time)")
  
  # Predict on chosen scale:
  #  - For MLMM: type="link" (η), then invlink(η) -> 0..30
  #  - For Gaussian: type="response" (already on data scale)
  pr <- predict(m, newdata=newdat, se.fit=TRUE, type=pred_type, exclude=excl)
  
  # Assemble prediction tibble.
  # We compute 95% CIs as invlink(η ± 1.96*SE) where η is the linear predictor.
  # For Gaussian, 'fit_lin' is on the response scale (identity link), so invlink = identity.
  preds <- newdat %>%
    mutate(
      fit_lin = as.numeric(pr$fit),
      se_lin  = as.numeric(pr$se.fit),
      fit = invlink(fit_lin),
      lwr = invlink(fit_lin - 1.96*se_lin),
      upr = invlink(fit_lin + 1.96*se_lin),
      response = resp
    ) %>%
    filter(!is.na(class), is.finite(fit), is.finite(lwr), is.finite(upr)) %>%
    filter(time <= max_follow)   # draw curves only where data exist
  
  list(model = m, preds = preds)
}

# ---- Summaries: parametric + smooths to a tidy data frame -------------------
# This function inspects mgcv's 'summary' output and extracts:
#   - Parametric terms (main effects and interactions): estimates, SEs, t/z stats, p-values
#   - Smooth terms (s(...)): F or Chi-square stats and p-values
# It detects the correct column names across families.
tidy_gam_summary <- function(m) {
  s <- summary(m)
  
  ## Parametric terms
  pt_raw <- as.data.frame(s$p.table)
  if (nrow(pt_raw)) {
    stat_col <- intersect(c("t value", "z value"), colnames(pt_raw))[1]
    p_col    <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(pt_raw))[1]
    if (is.na(stat_col)) {
      cand <- grep(" value$", colnames(pt_raw), value = TRUE)
      stat_col <- cand[1]
    }
    if (is.na(p_col)) {
      cand <- grep("^Pr", colnames(pt_raw), value = TRUE)
      p_col <- cand[1]
    }
    pt <- tibble::rownames_to_column(pt_raw, "term") |>
      dplyr::transmute(
        effect_type = "parametric",
        term,
        estimate    = .data[["Estimate"]],
        se          = .data[["Std. Error"]],
        stat        = .data[[stat_col]],
        p_value     = .data[[p_col]]
      )
  } else {
    pt <- tibble::tibble(effect_type=character(), term=character(),
                         estimate=double(), se=double(), stat=double(), p_value=double())
  }
  
  ## Smooth terms (s(...) and random effects)
  st_raw <- as.data.frame(s$s.table)
  if (nrow(st_raw)) {
    stat_col_s <- intersect(c("F", "Chi.sq"), colnames(st_raw))[1]
    if (is.na(stat_col_s)) {
      cand <- intersect(c("F value","Statistic"), colnames(st_raw))
      stat_col_s <- if (length(cand)) cand[1] else colnames(st_raw)[max(1, match("edf", colnames(st_raw))+2)]
    }
    p_col_s <- intersect(c("p-value", "p value", "p.value"), colnames(st_raw))[1]
    if (is.na(p_col_s)) {
      cand <- grep("p", colnames(st_raw), ignore.case = TRUE, value = TRUE)
      p_col_s <- cand[1]
    }
    st <- tibble::rownames_to_column(st_raw, "term") |>
      dplyr::transmute(
        effect_type = "smooth",
        term,
        estimate = NA_real_,
        se       = NA_real_,
        stat     = .data[[stat_col_s]],
        p_value  = .data[[p_col_s]]
      )
  } else {
    st <- tibble::tibble(effect_type=character(), term=character(),
                         estimate=double(), se=double(), stat=double(), p_value=double())
  }
  
  dplyr::bind_rows(pt, st)
}

# ---- Fit all responses ------------------------------------------------------
res <- lapply(RESPONSES, fit_one)
names(res) <- RESPONSES

# =========================
# 3) EXPORTS (summaries + predictions)
# =========================

# 3a) Write per-response prediction data (new_data + fit + CI)
for (resp in RESPONSES) {
  preds <- res[[resp]]$preds
  out_csv  <- file.path(pred_dir, paste0("pred_", resp, ".csv"))
  readr::write_csv(preds, out_csv)
}

# 3b) Summaries: parametric + smooths → CSV & XLSX
sum_list <- list()
for (resp in RESPONSES) {
  m  <- res[[resp]]$model
  td <- tidy_gam_summary(m) |> dplyr::mutate(response = resp, .before = 1)
  readr::write_csv(td, file.path(sum_dir, paste0("summary_", resp, ".csv")))
  sum_list[[resp]] <- td
}
writexl::write_xlsx(sum_list, path = file.path(sum_dir, "model_summaries.xlsx"))

# =========================
# 4) VISUALISATION (outside the loop)
# =========================

# Helper: build a plot from a response key
plot_one <- function(resp_key) {
  dfp <- res[[resp_key]]$preds
  
  # Palette aligned to actual factor levels (prevents empty fill errors)
  lv  <- levels(df$class)
  PAL <- CLASS_COLS[lv]
  if (any(is.na(PAL))) stop("Palette missing for class levels: ", paste(lv[is.na(PAL)], collapse=", "))
  
  p <- ggplot(dfp, aes(time, fit, color=class, fill=class)) +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.20, color=NA, na.rm=TRUE) +
    geom_line(linewidth=1.2, alpha=0.95, na.rm=TRUE) +
    scale_color_manual(values=PAL, limits=lv, drop=FALSE, na.translate=FALSE) +
    scale_fill_manual(values =PAL, limits=lv, drop=FALSE, na.translate=FALSE) +
    scale_x_continuous(limits=c(0,12), breaks=seq(0,12,2)) +
    labs(x="Time from Baseline (years)", y=YLABS[[resp_key]], color=NULL, fill=NULL) +
    theme_classic(base_size=13) +
    theme(legend.position="bottom", axis.line=element_line(color="black"))
  
  # Y-axes to match original using coord_cartesian (limits) + scale_y (breaks)
  if (resp_key == "zscore") {
    p <- p + coord_cartesian(ylim = c(-6, 0)) +
      scale_y_continuous(breaks = c(0, -2, -4, -6))
  }
  if (resp_key == "mlm") {
    p <- p + coord_cartesian(ylim = c(10, 30)) +
      scale_y_continuous(breaks = seq(10, 30, 5))
  }
  if (resp_key == "m_score_std") {
    p <- p + coord_cartesian(ylim = c(-6, 0)) +
      scale_y_continuous(breaks = c(0, -2, -4, -6))
  }
  if (resp_key == "func_score") {
    p <- p + coord_cartesian(ylim = c(-1, 0)) +
      scale_y_continuous(breaks = c(0, -0.25, -0.5, -0.75, -1))
  }
  p
}

# Build plots
pA <- plot_one("zscore")
pB <- plot_one("mlm")
pC <- plot_one("m_score_std")
pD <- plot_one("func_score")

# Save individual panels
ggsave(file.path(fig_dir, "panel_A_mPACC.png"), pA, width=6, height=5, dpi=300)
ggsave(file.path(fig_dir, "panel_B_MMSE.png"),  pB, width=6, height=5, dpi=300)
ggsave(file.path(fig_dir, "panel_C_MEM.png"),   pC, width=6, height=5, dpi=300)
ggsave(file.path(fig_dir, "panel_D_EXEC.png"),  pD, width=6, height=5, dpi=300)

# Combined 2×2 with shared legend at bottom
final <- (pA + pB) / (pC + pD) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(final)
ggsave(file.path(fig_dir, "combined_4panels.png"), final, width=12, height=10, dpi=300)
ggsave(file.path(fig_dir, "combined_4panels.pdf"), final,  width=12, height=10)
