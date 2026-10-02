#!/usr/bin/env Rscript
# =====================================================================
# SPAT — Learning rates et comparaison Séquence A vs B (Bloc 3)
#         GROUPES APPARIÉS (même participant en séquence A et B)
# =====================================================================
#
# Corrige la version précédente qui plantait sur presque tous les tests
# (2 à 9) parce que les fonctions attendaient soit des vecteurs alignés
# par participant (a_vals/b_vals appariés) soit un vecteur de
# différences (diffs), mais étaient appelées avec les mauvais
# arguments (participant_ids manquant, a_vals/b_vals passés là où
# diffs était attendu, fonction plot_paired_bootstrap non définie...).
#
# Le principe retenu partout est celui que vous avez donné pour l'AUC :
#
#   df_X_A <- setNames(block3[seq=="A", "X"], block3[seq=="A","participant"])
#   df_X_B <- setNames(block3[seq=="B", "X"], block3[seq=="B","participant"])
#   df_X_A <- df_X_A[!is.na(df_X_A)] ; df_X_B <- df_X_B[!is.na(df_X_B)]
#   common <- intersect(names(df_X_A), names(df_X_B))
#   X_A <- df_X_A[common] ; X_B <- df_X_B[common]
#   t.test(X_A, X_B, paired = TRUE)
#
# -> on ne garde QUE les participants présents dans les deux séquences,
#    on les remet dans le même ordre via l'intersection des noms, et on
#    applique un test apparié (paired = TRUE) sur ce sous-ensemble.
#
# Ce même schéma est maintenant appliqué à : t-test principal (1),
# bootstrap des différences (2), test du signe (3), taux exponentiel
# (6), AUC (7), TOST (9) et comparaison bayésienne appariée (10).
#
# USAGE
# -----
#   Rscript spat_learning_rate_comparison_paired.R scores_folder out_dir
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
  lr_ssm_A = list(label = "SSM speed of learning (A)",
                  ylabel = "SSM parameter A", filename = "ssm_A"),
  lr_ssm_mu = list(label = "SSM asymptotic performance (\u03bc)",
                   ylabel = "SSM parameter mu", filename = "ssm_mu"),
  lr_ssm_Q = list(label = "SSM process noise (Q)",
                  ylabel = "SSM parameter Q", filename = "ssm_Q"),
  lr_ssm_R = list(label = "SSM measurement noise (R)",
                  ylabel = "SSM parameter R", filename = "ssm_R")
)

SCORE_LABELS <- c(
  score_harmonic_norm_participant = "Harmonic score (normalized)"
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
# par une r\u00e9gression lin\u00e9aire classique (lm), qui converge toujours.
fit_sqrt_model <- function(x, y) {
  fit <- stats::lm(y ~ sqrt(x))
  co  <- stats::coef(fit)
  # on choisit c = 0 (a absorbe l'ordonn\u00e9e \u00e0 l'origine) pour rester
  # compatible avec la signature de sqrt_model(x, a, b, c)
  c(a = unname(co[1]), b = -unname(co[2]), c = 0)
}
# ---------------------------------------------------------------------
# 3bis. STATE-SPACE MODEL : x_{t+1} = mu + A*(x_t - mu) + w_t
#       (retour vers l'asymptote mu), observation y_t = x_t + v_t
#       Ajuste (A, mu, Q, R) par maximum de vraisemblance via un
#       filtre de Kalman (plus stable que d'estimer B puis de deriver
#       mu = B/(1-A), qui explose quand A -> 1).
# ---------------------------------------------------------------------

# Log-vraisemblance negative du filtre de Kalman pour une trajectoire y.
# par = c(A, mu, logQ, logR) ; on optimise logQ/logR pour garantir Q,R > 0.
kalman_negloglik <- function(par, y) {
  A <- par[1]; mu <- par[2]
  Q <- exp(par[3]); R <- exp(par[4])
  
  ok <- which(!is.na(y))
  if (length(ok) < 2) return(1e10)
  
  n <- length(y)
  x_pred <- y[ok[1]]          # etat initial = premiere observation
  P_pred <- R + Q             # incertitude initiale
  ll <- 0
  
  for (t in seq_len(n)) {
    yt <- y[t]
    if (!is.na(yt)) {
      innov <- yt - x_pred
      S <- P_pred + R
      if (!is.finite(S) || S <= 0) return(1e10)
      K <- P_pred / S
      x_upd <- x_pred + K * innov
      P_upd <- (1 - K) * P_pred
      ll <- ll - 0.5 * (log(2 * pi * S) + innov^2 / S)
    } else {
      x_upd <- x_pred; P_upd <- P_pred   # pas d'observation : on ne met pas a jour
    }
    x_pred <- mu + A * (x_upd - mu)
    P_pred <- A^2 * P_upd + Q
  }
  
  if (!is.finite(ll)) return(1e10)
  -ll
}

# Rejoue le filtre avec les parametres optimaux pour recuperer la
# trajectoire filtree xhat (utilisee ensuite pour LR 10-90 / max slope).
kalman_filter_states <- function(y, A, mu, Q, R) {
  ok <- which(!is.na(y))
  n <- length(y)
  xhat <- rep(NA_real_, n)
  if (length(ok) < 1) return(xhat)
  
  x_pred <- y[ok[1]]; P_pred <- R + Q
  for (t in seq_len(n)) {
    yt <- y[t]
    if (!is.na(yt)) {
      innov <- yt - x_pred
      S <- P_pred + R
      K <- P_pred / S
      x_upd <- x_pred + K * innov
      P_upd <- (1 - K) * P_pred
    } else {
      x_upd <- x_pred; P_upd <- P_pred
    }
    xhat[t] <- x_upd
    x_pred <- mu + A * (x_upd - mu)
    P_pred <- A^2 * P_upd + Q
  }
  xhat
}

fit_state_space_model <- function(x, y) {
  ok <- which(!is.na(y))
  if (length(ok) < 4) return(NULL)   # au moins 4 points pour identifier 4 parametres
  
  var_y <- stats::var(y[ok])
  if (!is.finite(var_y) || var_y < 1e-12) var_y <- 1e-6
  
  # Point de depart : A modere, mu = dernieres valeurs (asymptote),
  # Q et R se partagent une fraction raisonnable de la variance totale
  mu0 <- mean(utils::tail(y[ok], max(2, floor(length(ok) * 0.3))))
  start <- c(A = 0.8, mu = mu0, logQ = log(var_y * 0.3), logR = log(var_y * 0.3))
  
  opt <- tryCatch(
    stats::optim(par = start, fn = kalman_negloglik, y = y, method = "L-BFGS-B",
                 lower = c(0.001, -Inf, log(1e-8), log(1e-8)),
                 upper = c(0.999,  Inf, log(1e4),  log(1e4)),
                 control = list(maxit = 500)),
    error = function(e) NULL
  )
  if (is.null(opt)) return(NULL)
  
  A <- opt$par[1]; mu <- opt$par[2]
  Q <- exp(opt$par[3]); R <- exp(opt$par[4])
  
  xhat <- kalman_filter_states(y, A, mu, Q, R)
  
  list(A = as.numeric(A), mu = as.numeric(mu), Q = as.numeric(Q), R = as.numeric(R),
       B = as.numeric(mu * (1 - A)),   # gardee pour compatibilite descendante
       xhat = xhat, convergence = opt$convergence)
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
    list(a = unname(co["a"]), b = unname(co["b"]), c = unname(co["c"]), r2 = r2)
  }, error = function(e) list(a = NA_real_, b = NA_real_, c = NA_real_, r2 = NA_real_))
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
# 6. TAILLES D'EFFET : Cohen's d et Hedges' g (version appari\u00e9e "dz")
# ---------------------------------------------------------------------

# Pour des donn\u00e9es appari\u00e9es, la taille d'effet standard est
# dz = mean(diff) / sd(diff). On garde aussi la version "pooled"
# classique pour compatibilit\u00e9 / comparaison.
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

cohens_dz_paired <- function(diffs) {
  diffs <- diffs[is.finite(diffs)]
  n <- length(diffs)
  if (n < 2) return(NA_real_)
  sdd <- stats::sd(diffs)
  if (!is.finite(sdd) || sdd < 1e-12) return(NA_real_)
  mean(diffs) / sdd
}

# ---------------------------------------------------------------------
# 7. TEST D'EQUIVALENCE (TOST) — version appari\u00e9e (test \u00e0 un \u00e9chantillon
#    sur les diff\u00e9rences)
# ---------------------------------------------------------------------

tost_test_paired <- function(diffs, delta) {
  diffs <- diffs[is.finite(diffs)]
  n <- length(diffs)
  if (n < 2) return(c(p_tost = NA_real_, p_lower = NA_real_, p_upper = NA_real_))
  
  m  <- mean(diffs)
  se <- stats::sd(diffs) / sqrt(n)
  if (!is.finite(se) || se < 1e-12) return(c(p_tost = NA_real_, p_lower = NA_real_, p_upper = NA_real_))
  
  df_t <- n - 1
  t_lower <- (m + delta) / se
  t_upper <- (m - delta) / se
  
  p_lower <- 1 - stats::pt(t_lower, df = df_t)
  p_upper <- stats::pt(t_upper, df = df_t)
  p_tost  <- max(p_lower, p_upper)
  c(p_tost = p_tost, p_lower = p_lower, p_upper = p_upper)
}

# ---------------------------------------------------------------------
# 8. MIXED ANOVA (facteur inter = sequence, facteur intra = block)
# ---------------------------------------------------------------------

mixed_anova <- function(lr_df, score_col, lr_metric) {
  sub <- lr_df[lr_df$score_method == score_col & lr_df$block %in% c(3, 4),
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
    if (nrow(piv) < 3) return(list())
    long <- stats::reshape(piv, idvar = "participant_f", varying = names(piv)[-1],
                           v.names = "dv", timevar = "block_f", direction = "long")
    aov_fit <- stats::aov(dv ~ block_f + Error(participant_f / block_f), data = long)
    smry <- summary(aov_fit)
    F_val <- smry[[2]][[1]][["F value"]][1]
    p_val <- smry[[2]][[1]][["Pr(>F)"]][1]
    list(mixed_anova_F = F_val, mixed_anova_p = p_val,
         mixed_anova_eta2 = NA_real_, mixed_anova_method = "statsmodels_rm_only")
  }, error = function(e) list())
}

# ---------------------------------------------------------------------
# 9. LINEAR MIXED EFFECTS MODEL (LME)
# ---------------------------------------------------------------------

linear_mixed_effects <- function(all_df, score_col) {
  sub <- all_df[all_df$block %in% c(3, 4), c(score_col, "sequence", "block", "participant")]
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
  
  lr_ssm_1090 <- lr_ssm_slope <- lr_ssm_A <- lr_ssm_mu <- lr_ssm_Q <- lr_ssm_R <- NA_real_
  ssm_fit <- tryCatch(fit_state_space_model(x, y), error = function(e) NULL)
  if (!is.null(ssm_fit)) {
    y_ssm        <- ssm_fit$xhat
    lr_ssm_1090  <- learning_rate_10_90(y_ssm, x)
    lr_ssm_slope <- learning_rate_max_slope(y_ssm)
    lr_ssm_A     <- ssm_fit$A
    lr_ssm_mu    <- ssm_fit$mu
    lr_ssm_Q     <- ssm_fit$Q
    lr_ssm_R     <- ssm_fit$R
  }
  
  
  # --- sqrt-model slope: analytic derivative instead of max(diff()) ---
  lr_sqrt_slope <- NA_real_
  sqrt_fit <- tryCatch(fit_sqrt_model(x, y), error = function(e) NULL)
  if (!is.null(sqrt_fit)) {
    x_steep <- min(x)  # dy/dx is largest in magnitude at the smallest x
    if (is.finite(x_steep) && x_steep > 0) {
      # y = a - b*sqrt(x) + c  =>  dy/dx = -b / (2*sqrt(x))
      # sign flipped to report "rate of improvement" as positive when
      # the score is increasing with practice (b > 0 case); keep the
      # natural sign otherwise so callers can interpret decreasing scores too
      lr_sqrt_slope <- abs(unname(sqrt_fit["b"]) / (2 * sqrt(x_steep)))
    }
  }
  
  exp_fit <- fit_exp_model(x, y)
  auc     <- compute_auc(x, y)
  
  list(lr_ssm_1090 = lr_ssm_1090, lr_sqrt_slope = lr_sqrt_slope,
       lr_ssm_slope = lr_ssm_slope, lr_ssm_A = lr_ssm_A,
       lr_ssm_mu = lr_ssm_mu, lr_ssm_Q = lr_ssm_Q, lr_ssm_R = lr_ssm_R,
       exp_asymptote = exp_fit$a, exp_rate = exp_fit$c, exp_r2 = exp_fit$r2,
       auc = auc)
}

learning_rate_sqrt_analytic <- function(fs) {
  if (is.null(fs) || length(fs$A_t) == 0) return(NA_real_)
  fs$A_t[1] / 2   # dy/dn = A/(2*sqrt(n)) at n=1
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
          lr_ssm_mu = lr$lr_ssm_mu, lr_ssm_Q = lr$lr_ssm_Q, lr_ssm_R = lr$lr_ssm_R,
          exp_asymptote = lr$exp_asymptote, exp_rate = lr$exp_rate, exp_r2 = lr$exp_r2,
          auc = lr$auc, stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------
# 11. T-TEST APPARIE + ANOVA A MESURES REPETEES (facteur sequence)
# ---------------------------------------------------------------------
# a_vals / b_vals doivent d\u00e9j\u00e0 \u00eatre align\u00e9s (m\u00eame ordre de
# participants) ; participant_ids sert \u00e0 construire le facteur
# participant pour l'ANOVA \u00e0 mesures r\u00e9p\u00e9t\u00e9es.

paired_anova_ttest <- function(a_vals, b_vals, participant_ids) {
  tt <- stats::t.test(a_vals, b_vals, paired = TRUE)
  
  long_df <- data.frame(
    participant = factor(rep(participant_ids, 2)),
    sequence    = factor(c(rep("A", length(participant_ids)), rep("B", length(participant_ids)))),
    value       = c(a_vals, b_vals)
  )
  aov_res <- tryCatch({
    fit  <- stats::aov(value ~ sequence + Error(participant / sequence), data = long_df)
    smry <- summary(fit)
    F_val <- smry[[2]][[1]][["F value"]][1]
    p_val <- smry[[2]][[1]][["Pr(>F)"]][1]
    c(F = F_val, p = p_val)
  }, error = function(e) c(F = NA_real_, p = NA_real_))
  
  list(t_stat = unname(tt$statistic), t_p = unname(tt$p.value),
       F = unname(aov_res["F"]), p_anova = unname(aov_res["p"]))
}

# ---------------------------------------------------------------------
# 12. BOOTSTRAP APPARIE (sur les diff\u00e9rences A-B par participant)
# ---------------------------------------------------------------------

bootstrap_paired_mean_median <- function(diffs, n_boot = 10000) {
  diffs <- diffs[is.finite(diffs)]
  bm   <- numeric(n_boot)
  bmed <- numeric(n_boot)
  n <- length(diffs)
  for (i in seq_len(n_boot)) {
    s       <- sample(diffs, n, replace = TRUE)
    bm[i]   <- mean(s)
    bmed[i] <- stats::median(s)
  }
  ci_m   <- stats::quantile(bm,   c(0.025, 0.975))
  ci_med <- stats::quantile(bmed, c(0.025, 0.975))
  list(bm = bm, bmed = bmed, ci_m = ci_m, ci_med = ci_med,
       p_m = mean(bm <= 0), p_med = mean(bmed <= 0))
}

plot_paired_bootstrap <- function(boot, diffs, title, xlabel, fname) {
  diffs <- diffs[is.finite(diffs)]
  obs_mean_diff   <- mean(diffs)
  obs_median_diff <- stats::median(diffs)
  
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
  mtext(sprintf("%s\nBootstrap pair\u00e9 (diff\u00e9rences par participant)", title), outer = TRUE, cex = 1.0, font = 2)
  grDevices::dev.off()
}

# ---------------------------------------------------------------------
# 13. TEST DU SIGNE (bootstrap binomial sur les diff\u00e9rences)
# ---------------------------------------------------------------------

binomial_bootstrap_sign_test <- function(diffs, n_boot = 10000) {
  diffs <- diffs[is.finite(diffs)]
  sign_ind <- as.integer(diffs > 0)          # 1 si A > B
  n_pos <- sum(sign_ind[diffs != 0])
  n_neg <- sum(diffs < 0)
  n_eff <- n_pos + n_neg
  
  p_obs  <- if (n_eff == 0) NA_real_ else n_pos / n_eff
  p_sign <- if (n_eff == 0) NA_real_ else stats::binom.test(n_pos, n_eff, p = 0.5)$p.value
  
  n <- length(sign_ind)
  boot_props <- numeric(n_boot)
  for (i in seq_len(n_boot)) boot_props[i] <- mean(sample(sign_ind, n, replace = TRUE))
  ci <- stats::quantile(boot_props, c(0.025, 0.975))
  
  list(n_pos = n_pos, n_neg = n_neg, p_obs = p_obs, p_sign_test = p_sign,
       boot_props = boot_props, ci = ci)
}

plot_binomial_bootstrap <- function(result, title, fname) {
  grDevices::png(fname, width = 1400, height = 550, res = 130)
  par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))
  
  barplot(c(result$n_pos, result$n_neg), names.arg = c("A > B", "B > A"),
          col = c("#E64B35", "#4DBBD5"), ylab = "Nombre de participants",
          main = sprintf("Signes observ\u00e9s\np (test du signe) = %s",
                         ifelse(is.na(result$p_sign_test), "N/A", sprintf("%.4f", result$p_sign_test))),
          cex.main = 0.9)
  
  hist(result$boot_props, breaks = 40, freq = FALSE,
       col = grDevices::adjustcolor("#00A087", 0.5), border = NA,
       main = "Distribution bootstrap de la proportion",
       xlab = "Proportion de participants A > B", ylab = "Densit\u00e9", cex.main = 0.9)
  abline(v = 0.5, lty = 2)
  abline(v = result$ci, col = "#E64B35", lty = 3)
  if (!is.na(result$p_obs)) abline(v = result$p_obs, col = "black", lty = 3, lwd = 2)
  legend("topright", legend = sprintf("IC95%%=[%.3f, %.3f]\nobs=%.3f",
                                      result$ci[1], result$ci[2], result$p_obs),
         bty = "n", cex = 0.7)
  
  mtext(sprintf("%s\nBootstrap binomial (test du signe)", title), outer = TRUE, cex = 1.0, font = 2)
  grDevices::dev.off()
}

# ---------------------------------------------------------------------
# 14. COMPARAISON BAYESIENNE APPARIEE (test \u00e0 un \u00e9chantillon sur diffs)
# ---------------------------------------------------------------------

bayesian_paired_mean_diff <- function(diffs, iterations = 10000) {
  diffs <- diffs[is.finite(diffs)]
  if (!has_pkg("BayesFactor")) {
    cat("      Bay\u00e9sien ignor\u00e9 : le paquet 'BayesFactor' n'est pas install\u00e9\n")
    return(NULL)
  }
  tryCatch({
    bf     <- BayesFactor::ttestBF(x = diffs, mu = 0)
    chains <- BayesFactor::posterior(bf, iterations = iterations, progress = FALSE)
    mu     <- as.numeric(chains[, "mu"])
    list(mean_diff_draws = mu)
  }, error = function(e) {
    cat(sprintf("      Bay\u00e9sien \u00e9chou\u00e9 : %s\n", conditionMessage(e)))
    NULL
  })
}

plot_paired_bayesian <- function(draws, diffs, title, xlabel, fname) {
  diffs <- diffs[is.finite(diffs)]
  obs_mean_diff <- mean(diffs)
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
# Utilitaire : construit deux vecteurs align\u00e9s par participant
#              (m\u00eame logique que l'exemple donn\u00e9 pour l'AUC), et
#              renvoie NULL si moins de min_n participants communs.
# ---------------------------------------------------------------------

paired_vectors <- function(df, col, min_n = 3) {
  df_A <- stats::setNames(df[df$sequence == "A", col], df[df$sequence == "A", "participant"])
  df_B <- stats::setNames(df[df$sequence == "B", col], df[df$sequence == "B", "participant"])
  df_A <- df_A[!is.na(df_A)]
  df_B <- df_B[!is.na(df_B)]
  common <- intersect(names(df_A), names(df_B))
  if (length(common) < min_n) return(NULL)
  list(A = df_A[common], B = df_B[common], common = common)
}

# ---------------------------------------------------------------------
# 14bis. VIOLONS BLOC 3 : learning rates Sequence A vs Sequence B
# ---------------------------------------------------------------------
# Violon "maison" (densite miroir + mediane/IQR), sans dependance
# externe (pas de package vioplot/ggplot2 requis).
draw_violin <- function(x, at, col, width = 0.35, border = "black") {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(invisible(NULL))
  if (length(x) == 1 || stats::sd(x) < 1e-8) {
    segments(at - width, x[1], at + width, x[1], col = border, lwd = 2)
    points(at, x[1], pch = 19, col = border, cex = 0.8)
    return(invisible(NULL))
  }
  d <- stats::density(x, n = 256)
  d$y <- d$y / max(d$y) * width
  polygon(c(at - d$y, rev(at + d$y)), c(d$x, rev(d$x)),
          col = grDevices::adjustcolor(col, alpha.f = 0.55), border = border)
  q <- stats::quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  segments(at - width * 0.25, q[2], at + width * 0.25, q[2], lwd = 2.5, col = border)
  segments(at, q[1], at, q[3], lwd = 1.5, col = border)
  points(at, mean(x), pch = 23, bg = "white", col = border, cex = 0.9)
}

# Trace un violon (Seq. A vs Seq. B) pour chaque metrique de `metrics`,
# cote a cote, pour un score x bloc donne. Le n de participants
# apparies est affiche au-dessus de chaque panneau (il peut differer
# d'une metrique a l'autre si certains ajustements ont echoue pour
# certains participants).
plot_violin_lr_panel <- function(lr_df, score_col, metrics, out_file, title_suffix,
                                 block_num = 3, panel_width = 350, panel_height =500) {
  block_data <- lr_df[lr_df$score_method == score_col & lr_df$block == block_num, ]
  if (nrow(block_data) == 0) return(NULL)
  
  score_label <- if (score_col %in% names(SCORE_LABELS)) SCORE_LABELS[[score_col]] else score_col
  n_metrics   <- length(metrics)
  
  # --- Largeur/hauteur désormais paramétrables (au lieu de 480*n+150 fixe) ---
  grDevices::png(out_file, width = panel_width * n_metrics + 150, height = panel_height, res = 130)
  par(mfrow = c(1, n_metrics), mar = c(5, 4.5, 4, 1), oma = c(0, 0, 3, 0), bg = "#F8F9FA")
  
  col_A <- "#4C8CC5"; col_B <- "#F1BC18"
  
  # --- Convertit un p-value en étoiles de significativité ---
  sig_stars <- function(p) {
    if (is.na(p)) return("n.s.")
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    return("n.s.")
  }
  
  for (m in metrics) {
    meta  <- LR_META[[m]]
    label <- if (!is.null(meta)) meta$label else m
    ylab_ <- if (!is.null(meta)) meta$ylabel else m
    
    pv <- paired_vectors(block_data, m, min_n = 1)
    if (is.null(pv)) { plot.new(); title(main = label, cex.main = 0.9); next }
    
    a_vals <- pv$A; b_vals <- pv$B; n_paired <- length(pv$common)
    all_vals <- c(a_vals, b_vals); all_vals <- all_vals[is.finite(all_vals)]
    if (length(all_vals) == 0) { plot.new(); title(main = label, cex.main = 0.9); next }
    
    # --- p-value du test apparié (si assez de points), pour la barre de sig. ---
    p_val <- NA_real_
    if (n_paired >= 3) {
      tt <- tryCatch(stats::t.test(a_vals, b_vals, paired = TRUE), error = function(e) NULL)
      if (!is.null(tt)) p_val <- tt$p.value
    }
    
    yr  <- range(all_vals, na.rm = TRUE)
    pad <- diff(yr) * 0.12; if (!is.finite(pad) || pad == 0) pad <- 1
    # espace supplémentaire en haut pour la barre de significativité
    yr  <- c(yr[1] - pad, yr[2] + pad * 2.2)
    
    # --- x labels horizontaux (las = 1 au lieu de las = 2) ---
    plot(NA, xlim = c(0.5, 2.5), ylim = yr, xaxt = "n", xlab = "",
         ylab = ylab_, main = label, cex.main = 0.9)
    axis(1, at = c(1, 2), labels = c("Seq. A", "Seq. B"),
         las = 1, cex.axis = 0.8)
    grid(nx = NA, ny = NULL, col = "gray90")
    
    draw_violin(a_vals, at = 1, col = col_A, width= 0.20)
    draw_violin(b_vals, at = 2, col = col_B,width= 0.20)
    
    # --- Barre de significativité (pont + étoiles) ---
    if (!is.na(p_val)) {
      bar_y  <- max(all_vals, na.rm = TRUE) + pad * 0.9
      tick_h <- pad * 0.15
      segments(1, bar_y, 2, bar_y, lwd = 1.2)
      segments(1, bar_y - tick_h, 1, bar_y, lwd = 1.2)
      segments(2, bar_y - tick_h, 2, bar_y, lwd = 1.2)
      text(1.5, bar_y + pad * 0.45, labels = sig_stars(p_val), cex = 1.1, font = 2)
    }
    
    mtext(sprintf("n=%d", n_paired), side = 3, line = 0.2, cex = 0.78)
  }
  
  mtext(sprintf("Learning rates Bloc %s (A vs B)", block_num),
        outer = TRUE, cex = 1.05, font = 2)
  grDevices::dev.off()
  cat(sprintf("    Sauvegard\u00e9 : %s\n", out_file))
  out_file
}


# ---------------------------------------------------------------------
# 15. ORCHESTRATION PRINCIPALE : COMPARAISON A vs B, BLOC 3 (APPARIE)
# ---------------------------------------------------------------------

run_block3_comparison <- function(lr_df, all_df, score_col, lr_metric, out_dir,
                                  n_boot = 10000, run_bayesian = TRUE,
                                  bayes_iterations = 10000, equiv_delta = 0.5) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  block3 <- lr_df[lr_df$score_method == score_col & lr_df$block == 4, ]
  
  pv <- paired_vectors(block3, lr_metric, min_n = 3)
  if (is.null(pv)) {
    cat(sprintf("    Ignor\u00e9 %s/%s : moins de 3 participant(s) appari\u00e9(s).\n", score_col, lr_metric))
    return(NULL)
  }
  a_vals <- pv$A
  b_vals <- pv$B
  common <- pv$common
  diffs  <- a_vals - b_vals
  
  score_label <- if (score_col %in% names(SCORE_LABELS)) SCORE_LABELS[[score_col]] else score_col
  meta_entry  <- LR_META[[lr_metric]]
  lr_label    <- if (!is.null(meta_entry)) meta_entry$label else lr_metric
  lr_filename <- if (!is.null(meta_entry)) meta_entry$filename else lr_metric
  title       <- sprintf("%s  |  %s  |  Bloc 4 (A vs B, paired)", lr_label, score_label)
  fname_base  <- file.path(out_dir, sprintf("%s_%s", score_col, lr_filename))
  
  cat(sprintf("    %s / %s : n=%d participants appari\u00e9s  mean_A=%.4f  mean_B=%.4f  mean_diff=%.4f\n",
              score_col, lr_metric, length(common), mean(a_vals), mean(b_vals), mean(diffs)))
  
  # --- 1. t-test appari\u00e9 + ANOVA \u00e0 mesures r\u00e9p\u00e9t\u00e9es ---
  tt <- paired_anova_ttest(a_vals, b_vals, common)
  cat(sprintf("      [1] t-test appari\u00e9 : t=%.3f  p=%.4f   ANOVA : F=%s  p=%s\n",
              tt$t_stat, tt$t_p,
              ifelse(is.na(tt$F), "NA", sprintf("%.3f", tt$F)),
              ifelse(is.na(tt$p_anova), "NA", sprintf("%.4f", tt$p_anova))))
  
  # --- 2. Bootstrap appari\u00e9 (sur les diff\u00e9rences) ---
  boot <- bootstrap_paired_mean_median(diffs, n_boot = n_boot)
  cat(sprintf("      [2] Bootstrap moyenne (diff) : IC95%%=[%.4f %.4f]  p=%.4f\n",
              boot$ci_m[1], boot$ci_m[2], boot$p_m))
  plot_paired_bootstrap(boot, diffs, title,
                        sprintf("Diff\u00e9rence A\u2212B (%s)", lr_metric),
                        paste0(fname_base, "_bootstrap.png"))
  
  # --- 3. Test du signe ---
  binom_res <- binomial_bootstrap_sign_test(diffs, n_boot = n_boot)
  cat(sprintf("      [3] Test du signe : A>B=%d  B>A=%d  p=%s\n",
              binom_res$n_pos, binom_res$n_neg,
              ifelse(is.na(binom_res$p_sign_test), "NA", sprintf("%.4f", binom_res$p_sign_test))))
  plot_binomial_bootstrap(binom_res, title, paste0(fname_base, "_binomial.png"))
  
  # --- 4. Mixed ANOVA ---
  mixed_res <- mixed_anova(lr_df, score_col, lr_metric)
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
  
  # --- 6. Taux d'apprentissage exponentiel (appari\u00e9) ---
  exp_rate_diff <- NA_real_; exp_rate_p <- NA_real_
  pv_exp <- paired_vectors(block3, "exp_rate", min_n = 3)
  if (!is.null(pv_exp)) {
    tt2 <- stats::t.test(pv_exp$A, pv_exp$B, paired = TRUE)
    exp_rate_p    <- tt2$p.value
    exp_rate_diff <- mean(pv_exp$A) - mean(pv_exp$B)
    cat(sprintf("      [6] Taux exponentiel : mean_A=%.4f  mean_B=%.4f  diff=%.4f  p=%.4f\n",
                mean(pv_exp$A), mean(pv_exp$B), exp_rate_diff, exp_rate_p))
  } else {
    cat("      [6] Taux exponentiel : donn\u00e9es insuffisantes (appariement)\n")
  }
  
  # --- 7. AUC (appari\u00e9, comme l'exemple fourni) ---
  auc_p <- NA_real_
  df_auc_A <- stats::setNames(block3[block3$sequence == "A", "auc"],
                              block3[block3$sequence == "A", "participant"])
  df_auc_B <- stats::setNames(block3[block3$sequence == "B", "auc"],
                              block3[block3$sequence == "B", "participant"])
  df_auc_A <- df_auc_A[!is.na(df_auc_A)]
  df_auc_B <- df_auc_B[!is.na(df_auc_B)]
  
  common_auc <- intersect(names(df_auc_A), names(df_auc_B))
  
  if (length(common_auc) >= 3) {
    auc_A <- df_auc_A[common_auc]
    auc_B <- df_auc_B[common_auc]
    tt3   <- stats::t.test(auc_A, auc_B, paired = TRUE)
    auc_p <- tt3$p.value
    cat(sprintf("      [7] AUC : mean_A=%.2f  mean_B=%.2f  p=%.4f\n", mean(auc_A), mean(auc_B), auc_p))
  } else {
    cat("      [7] AUC : donn\u00e9es insuffisantes (appariement)\n")
  }
  
  # --- 8. Tailles d'effet ---
  d  <- cohens_d(a_vals, b_vals)          # version "pooled" classique
  g  <- hedges_g(a_vals, b_vals)
  dz <- cohens_dz_paired(diffs)           # version appari\u00e9e recommand\u00e9e
  cat(sprintf("      [8] Cohen's d=%.3f  Hedges' g=%.3f  Cohen's dz (appari\u00e9)=%.3f\n", d, g, dz))
  
  # --- 9. TOST appari\u00e9 (sur les diff\u00e9rences) ---
  tost <- tost_test_paired(diffs, delta = equiv_delta)
  cat(sprintf("      [9] TOST appari\u00e9 (delta=%.2f) : p_tost=%.4f  (p_lower=%.4f, p_upper=%.4f)\n",
              equiv_delta, tost["p_tost"], tost["p_lower"], tost["p_upper"]))
  
  # --- 10. Bay\u00e9sien appari\u00e9 ---
  bayes_row <- list()
  if (run_bayesian) {
    cat("      \u00c9chantillonnage bay\u00e9sien\u2026\n")
    draws_res <- bayesian_paired_mean_diff(diffs, iterations = bayes_iterations)
    if (!is.null(draws_res)) {
      plot_paired_bayesian(draws_res$mean_diff_draws, diffs, title,
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
    n_paired = length(common),
    mean_A = mean(a_vals), mean_B = mean(b_vals),
    mean_diff = mean(diffs),
    median_diff = stats::median(diffs),
    paired_ttest_t = tt$t_stat, paired_ttest_p = tt$t_p,
    rm_anova_F = tt$F, rm_anova_p = tt$p_anova,
    boot_mean_ci_low = unname(boot$ci_m[1]), boot_mean_ci_high = unname(boot$ci_m[2]), boot_p_mean = boot$p_m,
    boot_median_ci_low = unname(boot$ci_med[1]), boot_median_ci_high = unname(boot$ci_med[2]), boot_p_median = boot$p_med,
    binom_n_Awins      = binom_res$n_pos,
    binom_n_Bwins      = binom_res$n_neg,
    binom_p_sign_test  = binom_res$p_sign_test,
    binom_boot_ci_low  = unname(binom_res$ci[1]),
    binom_boot_ci_high = unname(binom_res$ci[2])
  ), mixed_res, lme_res, list(
    exp_rate_diff_A_minus_B = exp_rate_diff, exp_rate_p = exp_rate_p,
    auc_p = auc_p, cohens_d = d, hedges_g = g, cohens_dz = dz,
    tost_delta = equiv_delta, tost_p = unname(tost["p_tost"]),
    tost_p_lower = unname(tost["p_lower"]), tost_p_upper = unname(tost["p_upper"])
  ), bayes_row)
}

# ---------------------------------------------------------------------
# 16. PIPELINE PRINCIPAL
# ---------------------------------------------------------------------

main <- function(scores_folder, out_dir = "./output_learning_rates_normalized/SPAT_paired",
                 score_cols = NULL, n_boot = 10000, run_bayesian = TRUE,
                 bayes_iterations = 10000, equiv_delta = 0.5) {
  
  cat(strrep("=", 70), "\n")
  cat("SPAT \u2014 LEARNING RATES & COMPARAISON SEQUENCE A vs B (BLOC 3)\n")
  cat("         [Groupes appari\u00e9s \u2014 m\u00eame participant en A et B]\n")
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
  
  cat("\n[3/4] Violons Bloc 3 (learning rates, S\u00e9quence A vs B)\u2026\n")
  violin_dir <- file.path(out_dir, "block3_violin_plots")
  dir.create(violin_dir, showWarnings = FALSE, recursive = TRUE)
  for (score_col in score_cols) {
    plot_violin_lr_panel(
      lr_df, score_col,
      metrics    = c("lr_sqrt_slope", "lr_ssm_A", "lr_ssm_mu"),
      out_file   = file.path(violin_dir, sprintf("%s_violin_sqrt_A_mu.png", score_col)),
      title_suffix = "max slope \u221a-model + A + mu (SSM)",
      block_num  = 4
    )
    plot_violin_lr_panel(
      lr_df, score_col,
      metrics    = c("lr_sqrt_slope", "lr_ssm_A"),
      out_file   = file.path(violin_dir, sprintf("%s_violin_sqrt_A.png", score_col)),
      title_suffix = "max slope \u221a-model + A (SSM)",
      block_num  = 4
    )
  }
  
  cat("\n[4/4] Comparaison S\u00e9quence A vs B \u2014 Bloc 3 (appari\u00e9)\u2026\n")
  stats_dir <- file.path(out_dir, "block3_comparisons")
  results <- list()
  lr_metrics_all <- c(names(LR_META), "exp_rate", "exp_asymptote", "auc")
  
  for (score_col in score_cols) {
    for (lr_metric in lr_metrics_all) {
      if (!(lr_metric %in% names(lr_df))) next
      cat(sprintf("\n  --- %s / %s ---\n", score_col, lr_metric))
      res <- tryCatch(
        run_block3_comparison(lr_df, all_df, score_col, lr_metric, stats_dir,
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
  summary_csv <- file.path(out_dir, "block3_AB_comparison_summary_paired.csv")
  utils::write.csv(results_df, summary_csv, row.names = FALSE)
  
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("TERMIN\u00c9.\n")
  cat(sprintf("  Learning rates (tous blocs) : %s\n", lr_csv))
  cat(sprintf("  R\u00e9sum\u00e9 comparaison Bloc 3   : %s\n", summary_csv))
  cat(sprintf("  Violons Bloc 3 (A vs B)     : %s/\n", violin_dir))
  cat(sprintf("  Figures                     : %s/\n", stats_dir))
  cat(strrep("=", 70), "\n")
  invisible(list(lr_df = lr_df, results_df = results_df))
}

# ---------------------------------------------------------------------
# EXECUTION EN LIGNE DE COMMANDE
# ---------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  main(args[1], args[2])
} else {
  main("C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_scores_harmonic",
       "C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_plots_paired/SPAT_block4")
}
