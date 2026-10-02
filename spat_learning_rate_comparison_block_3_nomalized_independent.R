#!/usr/bin/env Rscript
# =====================================================================
# SPAT — Learning rates et comparaison Séquence A vs B (BLOC 3, normalisé
#         par rapport au BLOC 1)  —  GROUPES INDÉPENDANTS
#         (traduction R de learningRateComparisonIndependent.py, version
#         corrigée et étendue)
# =====================================================================
#
# Lit les fichiers Excel par participant produits par
# spat_export_scores.py (SPAT_<participant>_scores.xlsx), ajuste pour
# CHAQUE participant x séquence x score x bloc les modèles
# d'apprentissage (SSM, √-model, exponentiel), en extrait des
# indicateurs de "learning rate", NORMALISE le bloc 3 par rapport au
# bloc 1 (variation relative), puis compare la séquence A à la
# séquence B sur ce bloc 3 normalisé avec les méthodes suivantes :
#
#     1. ANOVA inter-groupes + t-test de Welch (groupes indépendants)
#     2. Bootstrap indépendant (moyenne / médiane)
#     3. Test de permutation (non-paramétrique)
#     4. Mixed ANOVA (facteur inter = sequence, facteur intra = block,
#        comparaison intra-sujet bloc 1 vs bloc 3)
#     5. Linear Mixed Effects Model (LME) sur les scores bruts bloc 1 vs 3
#     6. Comparaison des taux d'apprentissage (modèle exponentiel)
#     7. Aire sous la courbe d'apprentissage (AUC)
#     8. Tailles d'effet (Cohen's d, Hedges' g)
#     9. Test d'équivalence (TOST)
#    10. Comparaison bayésienne indépendante (via BayesFactor::ttestBF)
#    11. Comparaison des critères de modèle (AIC, BIC, log-vraisemblance)
#        du modèle exponentiel, sur le bloc 3 BRUT (non normalisé — voir
#        note ci-dessous)
#
# NOTE IMPORTANTE SUR LA NORMALISATION
# -------------------------------------
# Seuls les indicateurs de "learning rate" au sens propre (lr_ssm_1090,
# lr_sqrt_slope, lr_ssm_slope, lr_ssm_A, exp_rate, exp_asymptote, auc)
# sont normalisés bloc3/bloc1 (variation relative (B3-B1)/|B1|).
# Les critères de qualité d'ajustement du modèle exponentiel (AIC, BIC,
# logLik) ne sont PAS normalisés par rapport au bloc 1 : leur échelle
# dépend du nombre de points et de la variance résiduelle de chaque
# bloc séparément, une différence "B3-B1" ou un ratio n'aurait pas
# d'interprétation statistique claire. Ils sont comparés A vs B
# directement sur leurs valeurs brutes au bloc 3 (section 11).
#
# USAGE
# -----
#   Rscript spat_learning_rate_comparison_independent.R scores_folder out_dir
#
# Packages requis : readxl, tools
# Packages optionnels (repli automatique si absents) :
#   nlme (LME + fallback mixed ANOVA), afex (mixed ANOVA préférée),
#   BayesFactor (comparaison bayésienne)
# =====================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(tools)
})

has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# ---------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------

LR_META <- list(
  lr_ssm_1090 = list(label = "10\u201390 % learning rate (SSM)",
                     ylabel = "Learning rate (SSM 10\u201390)", filename = "1090"),
  lr_sqrt_slope = list(label = "Max slope (\u221a-model)",
                       ylabel = "Learning rate (\u221a max slope)", filename = "sqrt_slope"),
  lr_ssm_slope = list(label = "Max slope (SSM)",
                      ylabel = "Learning rate (SSM max slope)", filename = "ssm_slope"),
  lr_ssm_A = list(label = "SSM retention parameter A",
                  ylabel = "SSM parameter A", filename = "ssm_A")
)

SCORE_LABELS <- c(
  score_original      = "Original (inverted time\u00d7DTW\u00b2)",
  score_dtw_rmse       = "DTW-RMSE (composantes normalis\u00e9es)",
  score_harmonic       = "Harmonique  1/\u221a(d\u00e9viation\u00d7temps)",
  score_exp_mult       = "Exponentiel multiplicatif",
  score_hyperbolic     = "Hyperbolique",
  score_harmonic_td    = "Moyenne harmonique temps/DTW",
  score_linear_penalty = "P\u00e9nalit\u00e9 lin\u00e9aire"
)

# ---------------------------------------------------------------------
# 1. CHARGEMENT DES DONNEES
# ---------------------------------------------------------------------

combine_frames <- function(list_of_df) {
  all_cols <- unique(unlist(lapply(list_of_df, names)))
  list_of_df <- lapply(list_of_df, function(d) {
    missing <- setdiff(all_cols, names(d))
    for (m in missing) d[[m]] <- NA
    d[all_cols]
  })
  do.call(rbind, list_of_df)
}

load_participant_scores <- function(scores_folder) {
  files <- list.files(scores_folder, pattern = "^SPAT_.*_scores\\.xlsx$",
                      full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) {
    cat(sprintf("ERREUR : aucun fichier SPAT_*_scores.xlsx trouv\u00e9 dans %s\n", scores_folder))
    return(data.frame())
  }
  frames <- list()
  for (f in files) {
    stem        <- tools::file_path_sans_ext(basename(f))
    participant <- sub("^SPAT_", "", sub("_scores$", "", stem))
    df          <- as.data.frame(readxl::read_excel(f))
    df$participant <- participant
    frames[[length(frames) + 1]] <- df
    cat(sprintf("  Charg\u00e9 : %s  (%d lignes)\n", participant, nrow(df)))
  }
  combine_frames(frames)
}

# ---------------------------------------------------------------------
# 2. BLOCS
# ---------------------------------------------------------------------

get_block <- function(df, block_num) df[df$block == block_num & !is.na(df$block), ]

# ---------------------------------------------------------------------
# 3. MODELES D'APPRENTISSAGE (\u221a-model, SSM)
# ---------------------------------------------------------------------

sqrt_model <- function(x, a, b, c) a - b * sqrt(x) + c

# Le mod\u00e8le a - b*sqrt(x) + c est LIN\u00c9AIRE en sqrt(x) (a et c sont
# confondus : seule leur somme (a+c) est identifiable). On l'ajuste donc
# par une r\u00e9gression lin\u00e9aire classique (lm), qui converge toujours,
# plut\u00f4t que par nls (non-lin\u00e9aire, sensible aux valeurs initiales et
# sujet \u00e0 des \u00e9checs de convergence pour un mod\u00e8le qui est en fait
# lin\u00e9aire).
fit_sqrt_model <- function(x, y) {
  fit <- stats::lm(y ~ sqrt(x))
  co  <- stats::coef(fit)
  # on choisit c = 0 (a absorbe l'ordonn\u00e9e \u00e0 l'origine) pour rester
  # compatible avec la signature de sqrt_model(x, a, b, c)
  c(a = unname(co[1]), b = -unname(co[2]), c = 0)
}

fit_state_space_model <- function(x, y) {
  xi  <- y[-length(y)]
  xi1 <- y[-1]
  X   <- cbind(xi, 1)
  c_  <- qr.solve(X, xi1)
  list(A = as.numeric(c_[1]), B = as.numeric(c_[2]))
}

simulate_ssm <- function(y0, A, B, n) {
  out <- numeric(n); out[1] <- y0
  if (n > 1) for (i in 2:n) out[i] <- A * out[i - 1] + B
  out
}

learning_rate_10_90 <- function(curve, trials) {
  start <- curve[1]; end <- curve[length(curve)]
  span  <- end - start
  if (!is.finite(span) || abs(span) < 1e-12) return(NA_real_)
  t10 <- start + 0.10 * span
  t90 <- start + 0.90 * span
  idx10 <- which.min(abs(curve - t10))
  idx90 <- which.min(abs(curve - t90))
  trials[idx90] - trials[idx10]
}

learning_rate_max_slope <- function(curve) {
  if (length(curve) < 2) return(NA_real_)
  max(diff(curve))
}

# ---------------------------------------------------------------------
# 4. MODELE EXPONENTIEL  y = a + b*exp(-c*t)
#    (\u00e9tendu : renvoie aussi loglik / aic / bic pour la section 11)
# ---------------------------------------------------------------------

fit_exp_model <- function(x, y) {
  a0 <- y[length(y)]; b0 <- y[1] - y[length(y)]; c0 <- 0.05
  out <- tryCatch({
    fit <- suppressWarnings(stats::nls(y ~ a + b * exp(-c * x),
                                       start = list(a = a0, b = b0, c = c0),
                                       lower = c(-Inf, -Inf, 1e-6), algorithm = "port",
                                       control = stats::nls.control(maxiter = 500, warnOnly = TRUE)))
    co <- stats::coef(fit)
    y_pred <- co["a"] + co["b"] * exp(-co["c"] * x)
    ss_res <- sum((y - y_pred)^2)
    ss_tot <- sum((y - mean(y))^2)
    r2 <- if (ss_tot > 1e-12) 1 - ss_res / ss_tot else NA_real_
    list(a = unname(co["a"]), b = unname(co["b"]), c = unname(co["c"]), r2 = r2,
         loglik = as.numeric(stats::logLik(fit)),
         aic    = stats::AIC(fit),
         bic    = stats::BIC(fit))
  }, error = function(e) list(a = NA_real_, b = NA_real_, c = NA_real_, r2 = NA_real_,
                              loglik = NA_real_, aic = NA_real_, bic = NA_real_))
  out
}

# ---------------------------------------------------------------------
# 5. AIRE SOUS LA COURBE (AUC, trap\u00e8zes)
# ---------------------------------------------------------------------

compute_auc <- function(x, y) {
  if (length(x) < 2) return(NA_real_)
  o <- order(x); x <- x[o]; y <- y[o]
  sum(diff(x) * (y[-1] + y[-length(y)]) / 2)
}

# ---------------------------------------------------------------------
# 6. TAILLES D'EFFET : Cohen's d et Hedges' g
# ---------------------------------------------------------------------

cohens_d <- function(a_vals, b_vals) {
  n1 <- length(a_vals); n2 <- length(b_vals)
  if (n1 < 2 || n2 < 2) return(NA_real_)
  s_pooled <- sqrt(((n1 - 1) * stats::var(a_vals) + (n2 - 1) * stats::var(b_vals)) / (n1 + n2 - 2))
  if (!is.finite(s_pooled) || s_pooled < 1e-12) return(NA_real_)
  (mean(a_vals) - mean(b_vals)) / s_pooled
}

hedges_g <- function(a_vals, b_vals) {
  d <- cohens_d(a_vals, b_vals)
  if (is.na(d)) return(NA_real_)
  n1 <- length(a_vals); n2 <- length(b_vals)
  df <- n1 + n2 - 2
  J  <- if (df > 1) 1 - 3 / (4 * df - 1) else NA_real_
  if (is.na(J)) NA_real_ else d * J
}

# ---------------------------------------------------------------------
# 7. TEST D'EQUIVALENCE (TOST)
# ---------------------------------------------------------------------

tost_test <- function(a_vals, b_vals, delta) {
  diff_ <- mean(a_vals) - mean(b_vals)
  n1 <- length(a_vals); n2 <- length(b_vals)
  s1 <- stats::var(a_vals); s2 <- stats::var(b_vals)
  se <- sqrt(s1 / n1 + s2 / n2)
  if (!is.finite(se) || se < 1e-12) return(c(p_tost = NA_real_, p_lower = NA_real_, p_upper = NA_real_))
  
  t_lower <- (diff_ + delta) / se
  t_upper <- (diff_ - delta) / se
  df_welch <- (s1 / n1 + s2 / n2)^2 / ((s1 / n1)^2 / (n1 - 1) + (s2 / n2)^2 / (n2 - 1))
  
  p_lower <- 1 - stats::pt(t_lower, df = df_welch)
  p_upper <- stats::pt(t_upper, df = df_welch)
  p_tost  <- max(p_lower, p_upper)
  c(p_tost = p_tost, p_lower = p_lower, p_upper = p_upper)
}

# ---------------------------------------------------------------------
# 8. MIXED ANOVA (facteur inter = sequence, facteur intra = block)
#    Comparaison intra-sujet BLOC 1 vs BLOC 3 (baseline vs bloc cible),
#    facteur inter = sequence
# ---------------------------------------------------------------------

mixed_anova <- function(lr_df, score_col, lr_metric) {
  sub <- lr_df[lr_df$score_method == score_col & lr_df$block %in% c(1, 3),
               c("participant", "sequence", "block", lr_metric)]
  sub <- sub[stats::complete.cases(sub), ]
  if (nrow(sub) == 0) return(list())
  
  sub$block_f       <- factor(sub$block)
  sub$sequence_f    <- factor(sub$sequence)
  sub$participant_f <- factor(sub$participant)
  names(sub)[names(sub) == lr_metric] <- "dv"
  
  if (has_pkg("afex")) {
    res <- tryCatch({
      aov_fit <- afex::aov_ez(id = "participant_f", dv = "dv", data = sub,
                              between = "sequence_f", within = "block_f",
                              anova_table = list(es = "pes"))
      tab <- aov_fit$anova_table
      row <- tab[grepl(":", rownames(tab)), , drop = FALSE]
      if (nrow(row) == 0) list() else list(
        mixed_anova_F      = as.numeric(row[["F"]][1]),
        mixed_anova_p      = as.numeric(row[["Pr(>F)"]][1]),
        mixed_anova_eta2   = as.numeric(row[["pes"]][1]),
        mixed_anova_method = "afex"
      )
    }, error = function(e) list())
    if (length(res) > 0) return(res)
  }
  
  # Repli : ANOVA \u00e0 mesures r\u00e9p\u00e9t\u00e9es intra-bloc seulement (pas d'interaction propre)
  tryCatch({
    sub_agg <- stats::aggregate(dv ~ participant_f + block_f, data = sub, FUN = mean)
    piv <- stats::reshape(sub_agg[, c("participant_f", "block_f", "dv")],
                          idvar = "participant_f", timevar = "block_f", direction = "wide")
    piv <- piv[stats::complete.cases(piv), ]
    if (nrow(piv) < 3 || ncol(piv) != 3) return(list())
    
    tt <- stats::t.test(piv[[2]], piv[[3]], paired = TRUE)
    list(mixed_anova_F = unname(tt$statistic)^2,
         mixed_anova_p = tt$p.value,
         mixed_anova_eta2 = NA_real_,
         mixed_anova_method = "paired_ttest_2level_fallback")
  }, error = function(e) list())
}

# ---------------------------------------------------------------------
# 9. LINEAR MIXED EFFECTS MODEL (LME)
#    sur les scores bruts, BLOC 1 vs BLOC 3
# ---------------------------------------------------------------------

linear_mixed_effects <- function(all_df, score_col) {
  sub <- all_df[all_df$block %in% c(1, 3), c(score_col, "sequence", "block", "participant")]
  names(sub)[names(sub) == score_col] <- "score"
  sub$sequence_code <- as.integer(sub$sequence == "B")
  sub$block_code    <- as.integer(sub$block == 3)
  sub <- sub[stats::complete.cases(sub), ]
  if (nrow(sub) < 10) return(list())
  
  if (!has_pkg("nlme")) {
    cat("      LME \u00e9chou\u00e9 : le paquet 'nlme' n'est pas install\u00e9\n")
    return(list())
  }
  
  res <- tryCatch({
    model <- nlme::lme(score ~ sequence_code * block_code, random = ~1 | participant,
                       data = sub, method = "REML")
    smry  <- summary(model)$tTable
    out   <- list()
    for (term in c("sequence_code", "block_code", "sequence_code:block_code")) {
      if (term %in% rownames(smry)) {
        key <- paste0("lme_", term)
        out[[paste0(key, "_coef")]] <- unname(smry[term, "Value"])
        out[[paste0(key, "_p")]]    <- unname(smry[term, "p-value"])
      }
    }
    out$lme_aic <- stats::AIC(model)
    out$lme_n   <- nrow(sub)
    out
  }, error = function(e) {
    cat(sprintf("      LME \u00e9chou\u00e9 : %s\n", conditionMessage(e)))
    list()
  })
  res
}

# ---------------------------------------------------------------------
# 10. CALCUL DES LEARNING RATES PAR BLOC
# ---------------------------------------------------------------------
moving_avg3 <- function(y) {
  n <- length(y)
  if (n < 3) return(y)
  out <- y
  for (i in 2:(n - 1)) out[i] <- mean(y[(i - 1):(i + 1)])
  out
}

compute_learning_rates_for_block <- function(block_df, score_col) {
  block <- block_df[order(block_df$trial), ]
  if (nrow(block) < 4) return(NULL)
  
  x <- block$trial
  y <- moving_avg3(block[[score_col]])
  
  lr_ssm_1090 <- lr_ssm_slope <- lr_ssm_A <- NA_real_
  ssm_fit <- tryCatch(fit_state_space_model(x, y), error = function(e) NULL)
  if (!is.null(ssm_fit)) {
    y_ssm        <- simulate_ssm(y[1], ssm_fit$A, ssm_fit$B, length(y))
    lr_ssm_1090  <- learning_rate_10_90(y_ssm, x)
    lr_ssm_slope <- learning_rate_max_slope(y_ssm)
    lr_ssm_A     <- ssm_fit$A
  }
  
  # --- sqrt-model slope: analytic derivative instead of max(diff()) ---
  lr_sqrt_slope <- NA_real_
  sqrt_fit <- tryCatch(fit_sqrt_model(x, y), error = function(e) NULL)
  if (!is.null(sqrt_fit)) {
    x_steep <- min(x)  # dy/dx is largest in magnitude at the smallest x
    if (is.finite(x_steep) && x_steep > 0) {
      # y = a - b*sqrt(x) + c  =>  dy/dx = -b / (2*sqrt(x))
      lr_sqrt_slope <- abs(unname(sqrt_fit["b"]) / (2 * sqrt(x_steep)))
    }
  }
  
  exp_fit <- fit_exp_model(x, y)
  auc     <- compute_auc(x, y)
  
  list(lr_ssm_1090 = lr_ssm_1090, lr_sqrt_slope = lr_sqrt_slope,
       lr_ssm_slope = lr_ssm_slope, lr_ssm_A = lr_ssm_A,
       exp_asymptote = exp_fit$a, exp_rate = exp_fit$c, exp_r2 = exp_fit$r2,
       exp_loglik = exp_fit$loglik, exp_aic = exp_fit$aic, exp_bic = exp_fit$bic,
       auc = auc)
}

compute_all_learning_rates <- function(all_df, score_cols) {
  rows <- list()
  keys <- unique(all_df[, c("participant", "sequence")])
  for (i in seq_len(nrow(keys))) {
    part <- keys$participant[i]; seq_ <- keys$sequence[i]
    pdata <- all_df[all_df$participant == part & all_df$sequence == seq_, ]
    for (score_col in score_cols) {
      if (!(score_col %in% names(pdata))) next
      for (b in c(1, 2, 3, 4)) {
        block <- get_block(pdata, b)
        lr <- tryCatch(compute_learning_rates_for_block(block, score_col), error = function(e) NULL)
        if (is.null(lr)) next
        rows[[length(rows) + 1]] <- data.frame(
          participant = part, sequence = seq_, score_method = score_col, block = b,
          lr_ssm_1090 = lr$lr_ssm_1090, lr_sqrt_slope = lr$lr_sqrt_slope,
          lr_ssm_slope = lr$lr_ssm_slope, lr_ssm_A = lr$lr_ssm_A,
          exp_asymptote = lr$exp_asymptote, exp_rate = lr$exp_rate, exp_r2 = lr$exp_r2,
          exp_loglik = lr$exp_loglik, exp_aic = lr$exp_aic, exp_bic = lr$exp_bic,
          auc = lr$auc, stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------
# 10bis. NORMALISATION DES LEARNING RATES DU BLOC 3 PAR RAPPORT AU BLOC 1
# ---------------------------------------------------------------------
#
# method = "pct_change" : (B3 - B1) / |B1|   (variation relative)
# method = "ratio"       : B3 / B1
# method = "diff"        : B3 - B1           (variation absolue)
#
# /!\ N'inclure dans lr_metrics QUE des indicateurs de learning rate
# (pas AIC/BIC/logLik, voir note en t\u00eate de fichier).
normalize_lr_block3_vs_block1 <- function(lr_df, lr_metrics,
                                          method = c("pct_change", "ratio", "diff")) {
  method <- match.arg(method)
  
  b1 <- lr_df[lr_df$block == 1, c("participant", "sequence", "score_method", lr_metrics)]
  b3 <- lr_df[lr_df$block == 3, c("participant", "sequence", "score_method", lr_metrics)]
  names(b1)[names(b1) %in% lr_metrics] <- paste0(lr_metrics, "_b1")
  names(b3)[names(b3) %in% lr_metrics] <- paste0(lr_metrics, "_b3")
  
  merged <- merge(b1, b3, by = c("participant", "sequence", "score_method"))
  
  for (m in lr_metrics) {
    v1 <- merged[[paste0(m, "_b1")]]
    v3 <- merged[[paste0(m, "_b3")]]
    norm_col <- switch(method,
                       pct_change = (v3 - v1) / abs(v1),
                       ratio      = v3 / v1,
                       diff       = v3 - v1
    )
    norm_col[!is.finite(norm_col)] <- NA_real_
    merged[[m]] <- norm_col
  }
  
  merged$block <- 3
  merged[, c("participant", "sequence", "score_method", "block", lr_metrics)]
}

# ---------------------------------------------------------------------
# 11. ANOVA INTER-GROUPES + T-TEST DE WELCH
# ---------------------------------------------------------------------

independent_anova_ttest <- function(a_vals, b_vals) {
  tt <- stats::t.test(a_vals, b_vals, var.equal = FALSE)
  aov_res <- tryCatch({
    df_long <- data.frame(value = c(a_vals, b_vals),
                          group = factor(c(rep("A", length(a_vals)), rep("B", length(b_vals)))))
    f <- stats::oneway.test(value ~ group, data = df_long, var.equal = TRUE)
    c(F = unname(f$statistic), p = unname(f$p.value))
  }, error = function(e) c(F = NA_real_, p = NA_real_))
  list(t_stat = unname(tt$statistic), t_p = unname(tt$p.value),
       F = aov_res["F"], p_anova = aov_res["p"])
}

# ---------------------------------------------------------------------
# 12. BOOTSTRAP INDEPENDANT
# ---------------------------------------------------------------------

bootstrap_independent_mean_median <- function(a_vals, b_vals, n_boot = 10000) {
  bm   <- numeric(n_boot)
  bmed <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    sa <- sample(a_vals, length(a_vals), replace = TRUE)
    sb <- sample(b_vals, length(b_vals), replace = TRUE)
    bm[i]   <- mean(sa) - mean(sb)
    bmed[i] <- stats::median(sa) - stats::median(sb)
  }
  ci_m   <- stats::quantile(bm,   c(0.025, 0.975))
  ci_med <- stats::quantile(bmed, c(0.025, 0.975))
  list(bm = bm, bmed = bmed, ci_m = ci_m, ci_med = ci_med,
       p_m = mean(bm <= 0), p_med = mean(bmed <= 0))
}

plot_independent_bootstrap <- function(boot, a_vals, b_vals, title, xlabel, fname) {
  obs_mean_diff   <- mean(a_vals) - mean(b_vals)
  obs_median_diff <- stats::median(a_vals) - stats::median(b_vals)
  
  grDevices::png(fname, width = 1400, height = 550, res = 130)
  par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))
  
  panels <- list(
    list(samp = boot$bm,   ci = boot$ci_m,   p = boot$p_m,   obs = obs_mean_diff,   lbl = "moyenne", col = "#4DBBD5"),
    list(samp = boot$bmed, ci = boot$ci_med, p = boot$p_med, obs = obs_median_diff, lbl = "m\u00e9diane", col = "#00A087")
  )
  for (pn in panels) {
    hist(pn$samp, breaks = 60, freq = FALSE, col = grDevices::adjustcolor(pn$col, 0.4),
         border = NA, main = sprintf("Bootstrap %s\np = %.4f", pn$lbl, pn$p),
         xlab = xlabel, ylab = "Densit\u00e9", cex.main = 0.9)
    abline(v = 0, lty = 2, lwd = 1.5)
    abline(v = pn$ci, col = "#E64B35", lty = 3)
    abline(v = pn$obs, col = "black", lty = 3, lwd = 2)
    legend("topright", legend = sprintf("IC95%%=[%.3f, %.3f]\nobs=%.3f", pn$ci[1], pn$ci[2], pn$obs),
           bty = "n", cex = 0.7)
  }
  mtext(sprintf("%s\nBootstrap (groupes ind\u00e9pendants)", title), outer = TRUE, cex = 1.0, font = 2)
  grDevices::dev.off()
}

# ---------------------------------------------------------------------
# 13. TEST DE PERMUTATION
# ---------------------------------------------------------------------

permutation_test <- function(a_vals, b_vals, n_boot = 10000) {
  obs_diff <- mean(a_vals) - mean(b_vals)
  combined <- c(a_vals, b_vals)
  n_a <- length(a_vals)
  
  perm_diffs <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    perm <- sample(combined)
    perm_diffs[i] <- mean(perm[1:n_a]) - mean(perm[(n_a + 1):length(perm)])
  }
  p_perm <- mean(abs(perm_diffs) >= abs(obs_diff))
  
  boot_props <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    sa <- sample(a_vals, length(a_vals), replace = TRUE)
    sb <- sample(b_vals, length(b_vals), replace = TRUE)
    m  <- min(length(sa), length(sb))
    boot_props[i] <- mean(sample(sa, m) > sample(sb, m))
  }
  ci <- stats::quantile(boot_props, c(0.025, 0.975))
  
  list(obs_diff = obs_diff, p_perm = p_perm, perm_diffs = perm_diffs,
       boot_props = boot_props, ci = ci)
}

plot_permutation_test <- function(result, title, fname) {
  grDevices::png(fname, width = 1400, height = 550, res = 130)
  par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))
  
  hist(result$perm_diffs, breaks = 60, freq = FALSE, col = grDevices::adjustcolor("#4DBBD5", 0.5),
       border = NA, main = sprintf("Permutation\np (bilat\u00e9ral) = %.4f", result$p_perm),
       xlab = "Diff\u00e9rence de moyennes (A \u2212 B)", ylab = "Densit\u00e9", cex.main = 0.9)
  abline(v = result$obs_diff, col = "#E64B35", lwd = 2)
  abline(v = -abs(result$obs_diff), col = "#E64B35", lwd = 2, lty = 2)
  abline(v = 0, lty = 3)
  
  hist(result$boot_props, breaks = 40, freq = FALSE, col = grDevices::adjustcolor("#00A087", 0.5),
       border = NA, main = "Distribution bootstrap de P(A > B)",
       xlab = "P(A > B)", ylab = "Densit\u00e9", cex.main = 0.9)
  abline(v = 0.5, lty = 2)
  abline(v = result$ci, col = "#E64B35", lty = 3)
  
  mtext(sprintf("%s\nTest de permutation (groupes ind\u00e9pendants)", title), outer = TRUE, cex = 1.0, font = 2)
  grDevices::dev.off()
}

# ---------------------------------------------------------------------
# 14. COMPARAISON BAYESIENNE INDEPENDANTE
# ---------------------------------------------------------------------

bayesian_independent_mean_diff <- function(a_vals, b_vals, iterations = 10000) {
  if (!has_pkg("BayesFactor")) {
    cat("      Bay\u00e9sien ignor\u00e9 : le paquet 'BayesFactor' n'est pas install\u00e9\n")
    return(NULL)
  }
  tryCatch({
    bf     <- BayesFactor::ttestBF(x = a_vals, y = b_vals)
    chains <- BayesFactor::posterior(bf, iterations = iterations, progress = FALSE)
    beta   <- as.numeric(chains[, "beta (x - y)"])   # diff\u00e9rence de moyennes standardis\u00e9e * sigma
    sigma  <- as.numeric(chains[, "sig2"])^0.5
    mean_diff_draws <- beta * sigma
    list(mean_diff_draws = mean_diff_draws)
  }, error = function(e) {
    cat(sprintf("      Bay\u00e9sien \u00e9chou\u00e9 : %s\n", conditionMessage(e)))
    NULL
  })
}

plot_independent_bayesian <- function(draws, a_vals, b_vals, title, xlabel, fname) {
  obs_mean_diff <- mean(a_vals) - mean(b_vals)
  hdi <- stats::quantile(draws, c(0.03, 0.97))  # approx. HDI 94% par quantiles
  prob_gt <- mean(draws > 0)
  
  grDevices::png(fname, width = 900, height = 550, res = 130)
  hist(draws, breaks = 80, freq = FALSE, col = grDevices::adjustcolor("#4DBBD5", 0.45),
       border = NA, main = sprintf("Post\u00e9rieur diff\u00e9rence de moyennes\nP(diff>0)=%.3f", prob_gt),
       xlab = xlabel, ylab = "Densit\u00e9")
  abline(v = 0, lty = 2)
  abline(v = hdi, col = "#00A087", lty = 3, lwd = 2)
  abline(v = obs_mean_diff, col = "black", lty = 3, lwd = 2)
  legend("topright", legend = sprintf("HDI94%%~[%.3f, %.3f]\nobs=%.3f", hdi[1], hdi[2], obs_mean_diff),
         bty = "n", cex = 0.7)
  mtext(title, cex = 0.85)
  grDevices::dev.off()
}

# ---------------------------------------------------------------------
# 15. COMPARAISON DES CRITERES DE MODELE (AIC / BIC / logLik)
#     sur le BLOC 3 BRUT (non normalis\u00e9, voir note en t\u00eate de fichier)
# ---------------------------------------------------------------------

compare_model_criterion_A_vs_B <- function(lr_df, score_col, criterion, block_num = 3,
                                           n_boot = 10000, equiv_delta = 2) {
  sub <- lr_df[lr_df$score_method == score_col & lr_df$block == block_num, ]
  a_vals <- stats::na.omit(sub[sub$sequence == "A", criterion])
  b_vals <- stats::na.omit(sub[sub$sequence == "B", criterion])
  
  if (length(a_vals) < 3 || length(b_vals) < 3) {
    return(data.frame(
      score_method = score_col, criterion = criterion, block = block_num,
      n_A = length(a_vals), n_B = length(b_vals),
      note = "n insuffisant (minimum 3 par groupe)",
      stringsAsFactors = FALSE
    ))
  }
  
  tt <- stats::t.test(a_vals, b_vals, var.equal = FALSE)
  
  aov_res <- tryCatch({
    df_long <- data.frame(value = c(a_vals, b_vals),
                          group = factor(c(rep("A", length(a_vals)), rep("B", length(b_vals)))))
    f <- stats::oneway.test(value ~ group, data = df_long, var.equal = TRUE)
    c(F = unname(f$statistic), p = unname(f$p.value))
  }, error = function(e) c(F = NA_real_, p = NA_real_))
  
  d <- cohens_d(a_vals, b_vals)
  g <- hedges_g(a_vals, b_vals)
  
  boot <- bootstrap_independent_mean_median(a_vals, b_vals, n_boot = n_boot)
  perm <- permutation_test(a_vals, b_vals, n_boot = n_boot)
  tost <- tost_test(a_vals, b_vals, delta = equiv_delta)
  
  data.frame(
    score_method = score_col, criterion = criterion, block = block_num,
    n_A = length(a_vals), n_B = length(b_vals),
    mean_A = mean(a_vals), mean_B = mean(b_vals),
    median_A = stats::median(a_vals), median_B = stats::median(b_vals),
    mean_diff = mean(a_vals) - mean(b_vals),
    welch_t = unname(tt$statistic), welch_p = tt$p.value,
    oneway_anova_F = unname(aov_res["F"]), oneway_anova_p = unname(aov_res["p"]),
    cohens_d = d, hedges_g = g,
    boot_mean_ci_low  = unname(boot$ci_m[1]),
    boot_mean_ci_high = unname(boot$ci_m[2]),
    boot_p_mean = boot$p_m,
    boot_median_ci_low  = unname(boot$ci_med[1]),
    boot_median_ci_high = unname(boot$ci_med[2]),
    boot_p_median = boot$p_med,
    perm_obs_diff = perm$obs_diff, perm_p = perm$p_perm,
    tost_delta = equiv_delta, tost_p = unname(tost["p_tost"]),
    note = NA_character_,
    stringsAsFactors = FALSE
  )
}

run_model_criteria_comparisons <- function(lr_df, score_cols,
                                           criteria = c("exp_aic", "exp_bic", "exp_loglik"),
                                           blocks = c(3),
                                           out_dir = NULL,
                                           n_boot = 10000, equiv_delta = 2) {
  rows <- list()
  for (score_col in score_cols) {
    for (crit in criteria) {
      if (!(crit %in% names(lr_df))) {
        cat(sprintf("  Ignor\u00e9 : colonne '%s' absente de lr_df\n", crit))
        next
      }
      for (b in blocks) {
        cat(sprintf("  %s / %s / bloc %d\n", score_col, crit, b))
        res <- tryCatch(
          compare_model_criterion_A_vs_B(lr_df, score_col, crit, b,
                                         n_boot = n_boot, equiv_delta = equiv_delta),
          error = function(e) {
            cat(sprintf("    ERREUR : %s\n", conditionMessage(e)))
            NULL
          })
        if (!is.null(res)) rows[[length(rows) + 1]] <- res
      }
    }
  }
  
  out <- if (length(rows) > 0) combine_frames(rows) else data.frame()
  
  if (!is.null(out_dir) && nrow(out) > 0) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    csv_path <- file.path(out_dir, "model_criteria_AB_comparison.csv")
    utils::write.csv(out, csv_path, row.names = FALSE)
    cat(sprintf("-> R\u00e9sum\u00e9 crit\u00e8res de mod\u00e8le \u00e9crit dans : %s\n", csv_path))
  }
  
  out
}

# ---------------------------------------------------------------------
# 16. ORCHESTRATION PRINCIPALE : COMPARAISON A vs B, BLOC 3 (normalis\u00e9)
# ---------------------------------------------------------------------

run_block3_comparison <- function(lr_df, all_df, score_col, lr_metric, out_dir,
                                  n_boot = 10000, run_bayesian = TRUE,
                                  bayes_iterations = 10000, equiv_delta = 0.5) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # lr_df ici est le data.frame NORMALIS\u00e9 (bloc3 vs bloc1), qui ne
  # contient que des lignes marqu\u00e9es block == 3 (voir
  # normalize_lr_block3_vs_block1()).
  block3 <- lr_df[lr_df$score_method == score_col & lr_df$block == 3, ]
  a_vals <- stats::na.omit(block3[block3$sequence == "A", lr_metric])
  b_vals <- stats::na.omit(block3[block3$sequence == "B", lr_metric])
  
  if (length(a_vals) < 3 || length(b_vals) < 3) {
    cat(sprintf("    Ignor\u00e9 %s/%s : n_A=%d, n_B=%d (minimum 3 par groupe).\n",
                score_col, lr_metric, length(a_vals), length(b_vals)))
    return(NULL)
  }
  
  score_label <- if (score_col %in% names(SCORE_LABELS)) SCORE_LABELS[[score_col]] else score_col
  meta_entry  <- LR_META[[lr_metric]]
  lr_label    <- if (!is.null(meta_entry)) meta_entry$label else lr_metric
  lr_filename <- if (!is.null(meta_entry)) meta_entry$filename else lr_metric
  title       <- sprintf("%s  |  %s  |  Bloc 3 normalis\u00e9/Bloc 1 (A vs B)", lr_label, score_label)
  fname_base  <- file.path(out_dir, sprintf("%s_%s", score_col, lr_filename))
  
  cat(sprintf("    %s / %s : n_A=%d  n_B=%d  mean_A=%.4f  mean_B=%.4f  mean_diff=%.4f\n",
              score_col, lr_metric, length(a_vals), length(b_vals),
              mean(a_vals), mean(b_vals), mean(a_vals) - mean(b_vals)))
  
  # --- 1. Welch t-test + ANOVA inter-groupes ---
  tt <- independent_anova_ttest(a_vals, b_vals)
  cat(sprintf("      [1] t-test Welch : t=%.3f  p=%.4f   ANOVA : F=%.3f  p=%.4f\n",
              tt$t_stat, tt$t_p, tt$F, tt$p_anova))
  
  # --- 2. Bootstrap ---
  boot <- bootstrap_independent_mean_median(a_vals, b_vals, n_boot = n_boot)
  cat(sprintf("      [2] Bootstrap moyenne : IC95%%=[%.4f %.4f]  p=%.4f\n",
              boot$ci_m[1], boot$ci_m[2], boot$p_m))
  plot_independent_bootstrap(boot, a_vals, b_vals, title,
                             sprintf("Diff\u00e9rence A\u2212B (%s)", lr_metric),
                             paste0(fname_base, "_bootstrap.png"))
  
  # --- 3. Test de permutation ---
  perm_res <- permutation_test(a_vals, b_vals, n_boot = n_boot)
  cat(sprintf("      [3] Permutation : diff obs=%.4f  p=%.4f\n", perm_res$obs_diff, perm_res$p_perm))
  plot_permutation_test(perm_res, title, paste0(fname_base, "_permutation.png"))
  
  # --- 4. Mixed ANOVA (bloc 1 vs bloc 3, sur les valeurs BRUTES, pas normalis\u00e9es) ---
  mixed_res <- mixed_anova(all_lr_df_raw, score_col, lr_metric)
  if (length(mixed_res) > 0) {
    cat(sprintf("      [4] Mixed ANOVA interaction : F=%.3f  p=%.4f  \u03b7\u00b2=%s\n",
                mixed_res$mixed_anova_F, mixed_res$mixed_anova_p,
                ifelse(is.na(mixed_res$mixed_anova_eta2), "NA", sprintf("%.3f", mixed_res$mixed_anova_eta2))))
  } else {
    cat("      [4] Mixed ANOVA : donn\u00e9es insuffisantes\n")
  }
  
  # --- 5. Linear Mixed Effects ---
  lme_res <- linear_mixed_effects(all_df, score_col)
  if (length(lme_res) > 0) {
    inter_p <- lme_res[["lme_sequence_code:block_code_p"]]
    cat(sprintf("      [5] LME interaction sequence\u00d7block : p=%s\n",
                ifelse(is.null(inter_p), "NA", sprintf("%.4f", inter_p))))
  } else {
    cat("      [5] LME : donn\u00e9es insuffisantes ou erreur\n")
  }
  
  # --- 6. Taux d'apprentissage exponentiel (normalis\u00e9 bloc3/bloc1) ---
  exp_c_A <- stats::na.omit(block3[block3$sequence == "A", "exp_rate"])
  exp_c_B <- stats::na.omit(block3[block3$sequence == "B", "exp_rate"])
  exp_rate_diff <- NA_real_; exp_rate_p <- NA_real_
  if (length(exp_c_A) >= 3 && length(exp_c_B) >= 3) {
    tt2 <- stats::t.test(exp_c_A, exp_c_B, var.equal = FALSE)
    exp_rate_p    <- tt2$p.value
    exp_rate_diff <- mean(exp_c_A) - mean(exp_c_B)
    cat(sprintf("      [6] Taux exponentiel : mean_A=%.4f  mean_B=%.4f  diff=%.4f  p=%.4f\n",
                mean(exp_c_A), mean(exp_c_B), exp_rate_diff, exp_rate_p))
  } else {
    cat("      [6] Taux exponentiel : donn\u00e9es insuffisantes\n")
  }
  
  # --- 7. AUC (normalis\u00e9e bloc3/bloc1) ---
  auc_A <- stats::na.omit(block3[block3$sequence == "A", "auc"])
  auc_B <- stats::na.omit(block3[block3$sequence == "B", "auc"])
  auc_p <- NA_real_
  if (length(auc_A) >= 3 && length(auc_B) >= 3) {
    tt3 <- stats::t.test(auc_A, auc_B, var.equal = FALSE)
    auc_p <- tt3$p.value
    cat(sprintf("      [7] AUC : mean_A=%.2f  mean_B=%.2f  p=%.4f\n", mean(auc_A), mean(auc_B), auc_p))
  } else {
    cat("      [7] AUC : donn\u00e9es insuffisantes\n")
  }
  
  # --- 8. Tailles d'effet ---
  d <- cohens_d(a_vals, b_vals)
  g <- hedges_g(a_vals, b_vals)
  cat(sprintf("      [8] Cohen's d=%.3f  Hedges' g=%.3f\n", d, g))
  
  # --- 9. TOST ---
  tost <- tost_test(a_vals, b_vals, delta = equiv_delta)
  cat(sprintf("      [9] TOST (delta=%.2f) : p_tost=%.4f  (p_lower=%.4f, p_upper=%.4f)\n",
              equiv_delta, tost["p_tost"], tost["p_lower"], tost["p_upper"]))
  
  # --- Bay\u00e9sien ---
  bayes_row <- list()
  if (run_bayesian) {
    cat("      \u00c9chantillonnage bay\u00e9sien\u2026\n")
    draws_res <- bayesian_independent_mean_diff(a_vals, b_vals, iterations = bayes_iterations)
    if (!is.null(draws_res)) {
      plot_independent_bayesian(draws_res$mean_diff_draws, a_vals, b_vals, title,
                                sprintf("Diff\u00e9rence A\u2212B (%s)", lr_metric),
                                paste0(fname_base, "_bayesian.png"))
      hdi <- stats::quantile(draws_res$mean_diff_draws, c(0.03, 0.97))
      bayes_row <- list(
        bayes_mean_diff     = mean(draws_res$mean_diff_draws),
        bayes_mean_hdi_low  = unname(hdi[1]),
        bayes_mean_hdi_high = unname(hdi[2]),
        bayes_p_mean_gt0    = mean(draws_res$mean_diff_draws > 0)
      )
    }
  }
  
  c(list(
    score_method = score_col, lr_metric = lr_metric,
    n_A = length(a_vals), n_B = length(b_vals),
    mean_A = mean(a_vals), mean_B = mean(b_vals),
    mean_diff = mean(a_vals) - mean(b_vals),
    median_diff = stats::median(a_vals) - stats::median(b_vals),
    welch_ttest_t = tt$t_stat, welch_ttest_p = tt$t_p,
    oneway_anova_F = tt$F, oneway_anova_p = tt$p_anova,
    boot_mean_ci_low = unname(boot$ci_m[1]), boot_mean_ci_high = unname(boot$ci_m[2]), boot_p_mean = boot$p_m,
    boot_median_ci_low = unname(boot$ci_med[1]), boot_median_ci_high = unname(boot$ci_med[2]), boot_p_median = boot$p_med,
    perm_obs_diff = perm_res$obs_diff, perm_p = perm_res$p_perm,
    perm_boot_ci_low = unname(perm_res$ci[1]), perm_boot_ci_high = unname(perm_res$ci[2])
  ), mixed_res, lme_res, list(
    exp_rate_diff_A_minus_B = exp_rate_diff, exp_rate_p = exp_rate_p,
    auc_p = auc_p, cohens_d = d, hedges_g = g,
    tost_delta = equiv_delta, tost_p = unname(tost["p_tost"]),
    tost_p_lower = unname(tost["p_lower"]), tost_p_upper = unname(tost["p_upper"])
  ), bayes_row)
}

# ---------------------------------------------------------------------
# 17. PIPELINE PRINCIPAL
# ---------------------------------------------------------------------

main <- function(scores_folder, out_dir = "./output_learning_rates_normalized/SPAT_unpaired",
                 score_cols = NULL, n_boot = 10000, run_bayesian = TRUE,
                 bayes_iterations = 10000, equiv_delta = 0.5) {
  
  cat(strrep("=", 70), "\n")
  cat("SPAT \u2014 LEARNING RATES & COMPARAISON SEQUENCE A vs B\n")
  cat("         BLOC 3, normalis\u00e9 par rapport au BLOC 1\n")
  cat("         [Groupes ind\u00e9pendants \u2014 non appari\u00e9s]\n")
  cat(strrep("=", 70), "\n")
  
  cat(sprintf("\n[1/4] Chargement des scores depuis : %s\n", scores_folder))
  all_df <- load_participant_scores(scores_folder)
  if (nrow(all_df) == 0) return(invisible(NULL))
  
  available_scores <- names(all_df)[startsWith(names(all_df), "score_")]
  if (is.null(score_cols)) {
    score_cols <- available_scores
  } else {
    score_cols <- intersect(score_cols, available_scores)
  }
  cat(sprintf("  Scores utilis\u00e9s : %s\n", paste(score_cols, collapse = ", ")))
  
  cat("\n[2/4] Ajustement des mod\u00e8les par bloc\u2026\n")
  lr_df <- compute_all_learning_rates(all_df, score_cols)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  lr_csv <- file.path(out_dir, "learning_rates_all_blocks.csv")
  utils::write.csv(lr_df, lr_csv, row.names = FALSE)
  cat(sprintf("  %d lignes -> %s\n", nrow(lr_df), lr_csv))
  
  # rendu disponible globalement pour mixed_anova() dans run_block3_comparison()
  # (donn\u00e9es BRUTES, non normalis\u00e9es, n\u00e9cessaires pour comparer bloc1 vs bloc3)
  all_lr_df_raw <<- lr_df
  
  cat("\n[3/4] Normalisation Bloc 3 / Bloc 1 (indicateurs de learning rate uniquement)\u2026\n")
  # Seuls les indicateurs de learning rate proprement dits sont normalis\u00e9s.
  # AIC / BIC / logLik sont volontairement exclus (voir note en t\u00eate de fichier).
  lr_metrics_all <- c(names(LR_META), "exp_rate", "exp_asymptote", "auc")
  
  lr_df_norm <- normalize_lr_block3_vs_block1(lr_df, lr_metrics_all, method = "pct_change")
  utils::write.csv(lr_df_norm, file.path(out_dir, "learning_rates_block3_normalized.csv"),
                   row.names = FALSE)
  cat(sprintf("  %d lignes normalis\u00e9es -> %s\n", nrow(lr_df_norm),
              file.path(out_dir, "learning_rates_block3_normalized.csv")))
  
  cat("\n[4/4] Comparaison S\u00e9quence A vs B \u2014 Bloc 3 normalis\u00e9\u2026\n")
  stats_dir <- file.path(out_dir, "block3_normalized_comparisons")
  results <- list()
  
  for (score_col in score_cols) {
    for (lr_metric in lr_metrics_all) {
      if (!(lr_metric %in% names(lr_df))) next
      cat(sprintf("\n  --- %s / %s ---\n", score_col, lr_metric))
      res <- tryCatch(
        run_block3_comparison(lr_df_norm, all_df, score_col, lr_metric, stats_dir,
                              n_boot = n_boot, run_bayesian = run_bayesian,
                              bayes_iterations = bayes_iterations, equiv_delta = equiv_delta),
        error = function(e) {
          cat(sprintf("    ERREUR sur %s/%s : %s\n", score_col, lr_metric, conditionMessage(e)))
          NULL
        })
      if (!is.null(res)) results[[length(results) + 1]] <- res
    }
  }
  
  if (length(results) > 0) {
    results_df <- combine_frames(lapply(results, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  } else {
    results_df <- data.frame()
  }
  summary_csv <- file.path(out_dir, "block3_normalized_AB_comparison_summary.csv")
  utils::write.csv(results_df, summary_csv, row.names = FALSE)
  
  cat("\n[Bonus] Comparaison des crit\u00e8res de mod\u00e8le exponentiel (AIC/BIC/logLik) \u2014 Bloc 3 brut\u2026\n")
  criteria_dir <- file.path(out_dir, "model_criteria")
  criteria_results <- run_model_criteria_comparisons(
    lr_df      = lr_df,     # valeurs BRUTES, non normalis\u00e9es (voir note)
    score_cols = score_cols,
    criteria   = c("exp_aic", "exp_bic", "exp_loglik"),
    blocks     = c(3),
    out_dir    = criteria_dir,
    n_boot     = n_boot,
    equiv_delta = equiv_delta
  )
  
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("TERMIN\u00c9.\n")
  cat(sprintf("  Learning rates (tous blocs)         : %s\n", lr_csv))
  cat(sprintf("  Learning rates Bloc3/Bloc1 normalis\u00e9s : %s\n",
              file.path(out_dir, "learning_rates_block3_normalized.csv")))
  cat(sprintf("  R\u00e9sum\u00e9 comparaison Bloc 3 normalis\u00e9   : %s\n", summary_csv))
  cat(sprintf("  R\u00e9sum\u00e9 crit\u00e8res de mod\u00e8le (AIC/BIC/LL) : %s\n",
              file.path(criteria_dir, "model_criteria_AB_comparison.csv")))
  cat(sprintf("  Figures                              : %s/\n", stats_dir))
  cat(strrep("=", 70), "\n")
  invisible(list(lr_df = lr_df, lr_df_norm = lr_df_norm,
                 results_df = results_df, criteria_results_df = criteria_results))
}

# ---------------------------------------------------------------------
# EXECUTION EN LIGNE DE COMMANDE
# ---------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  main(args[1], args[2])
} else {
  main("C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_scores_harmonic/SPAT", "C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_plots/all")
}