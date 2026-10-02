#!/usr/bin/env Rscript
# =====================================================================
# SPAT — R\u00e9cup\u00e9ration de performance en bloc 4 par rapport au bloc 1,
#         \u00e9valu\u00e9e \u00e0 l'essai n\u00b0 20, compar\u00e9e entre S\u00e9quence A et B
# =====================================================================
#
# QUESTION
# --------
# "\u00c0 l'essai 20, quel pourcentage de la performance du bloc 1 le
#  bloc 4 a-t-il r\u00e9cup\u00e9r\u00e9 ?" — puis compare cette r\u00e9cup\u00e9ration entre
#  la S\u00e9quence A et la S\u00e9quence B.
#
# D\u00c9FINITION
# -----------
# Pour chaque participant x s\u00e9quence x score :
#   v1 = valeur du score \u00e0 l'essai local n\u00b0 TARGET_TRIAL (20 par
#        d\u00e9faut) DANS LE BLOC 1 (interpolation lin\u00e9aire si l'essai
#        20 n'existe pas exactement, born\u00e9e \u00e0 l'intervalle observ\u00e9)
#   v4 = valeur du score \u00e0 l'essai local n\u00b0 TARGET_TRIAL DANS LE
#        BLOC 4 (m\u00eame proc\u00e9d\u00e9)
#
#   recovery_pct_raw      = 100 * v4 / v1
#       ("le bloc 4 repr\u00e9sente X % de la valeur brute du bloc 1")
#
#   recovery_pct_oriented = 100 * v4 / v1        si le score est
#                                                 "plus grand = meilleur"
#                         = 100 * v1 / v4        si le score est
#                                                 "plus petit = meilleur"
#       (orient\u00e9 pour que > 100 % signifie toujours "meilleure
#        performance en bloc 4 qu'en bloc 1 \u00e0 l'essai 20", quel que
#        soit le sens du score — voir SCORE_META$higher_is_better)
#
# Les deux quantit\u00e9s sont calcul\u00e9es \u00e0 la fois sur :
#   - la valeur BRUTE interpol\u00e9e (bruit\u00e9e, m\u00e9thode "raw")
#   - la valeur LISS\u00c9E par un mod\u00e8le \u00e9tat-espace lin\u00e9aire ajust\u00e9 par
#     filtre de Kalman + lisseur RTS SUR CHAQUE BLOC S\u00c9PAR\u00c9MENT (m\u00e9thode
#     "smooth"), plus robuste au bruit essai-par-essai
#
# Puis compar\u00e9 entre S\u00e9quence A et S\u00e9quence B (appari\u00e9 par
# participant si possible, sinon test de groupes ind\u00e9pendants de
# Welch en repli), avec : t-test, Wilcoxon, bootstrap (IC + p),
# permutation, Cohen's dz / Hedges' g, TOST (\u00e9quivalence).
#
# SORTIES (dans out_dir) :
#   recovery_block4_vs_block1_trial20.csv   toutes les valeurs calcul\u00e9es
#   recovery_comparison_stats.csv           comparaison stat. A vs B, par score/m\u00e9thode
#   recovery_points_<score>_<methode>.png   nuage de points appari\u00e9s A vs B
#
# USAGE
# -----
#   Rscript spat_block4_recovery_trial20.R scores_folder out_dir [target_trial]
#
# Packages requis : readxl, tools
# =====================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(tools)
})

# ---------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------

TARGET_TRIAL_DEFAULT <- 20
N_BOOT <- 10000
MIN_POINTS_SSM <- 5     # nb minimum de points pour ajuster le mod\u00e8le d'\u00e9tat-espace sur un bloc
MIN_PAIRED_N   <- 4     # nb minimum de participants appari\u00e9s pour un test appari\u00e9

SEQ_COLORS <- c(A = "#4C8CC5", B = "#F1BC18")

SCORE_META <- list(
  dtw_dist = list(label = "Dynamic time warping distance", higher_is_better = FALSE),
  time_spent_raw = list(label = "Time spent", higher_is_better = FALSE),
  time_spent_norm_participant = list(label = "Time spent (norm. participant)", higher_is_better = FALSE),
  time_spent_norm_global = list(label = "Time spent (norm. global)", higher_is_better = FALSE),
  mean_dev_raw = list(label = "Mean distance", higher_is_better = FALSE),
  mean_dev_norm_participant = list(label = "Mean distance (norm. participant)", higher_is_better = FALSE),
  mean_dev_norm_global = list(label = "Mean distance (norm. global)", higher_is_better = FALSE),
  score_harmonic_raw = list(label = "Harmonic score", higher_is_better = TRUE),
  score_harmonic_norm_participant = list(label = "Harmonic score (norm. participant)", higher_is_better = TRUE),
  score_harmonic_norm_global = list(label = "Harmonic score (norm. global)", higher_is_better = TRUE)
)

# ---------------------------------------------------------------------
# 1. CHARGEMENT DES DONNEES (identique aux autres scripts SPAT)
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

read_participant_file <- function(path) {
  sheets <- readxl::excel_sheets(path)
  sheet  <- if ("scores" %in% sheets) "scores" else sheets[1]
  df     <- as.data.frame(readxl::read_excel(path, sheet = sheet))
  stem        <- tools::file_path_sans_ext(basename(path))
  participant <- sub("^SPAT_", "", sub("_scores$", "", stem))
  df$participant <- participant
  df
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
    df <- read_participant_file(f)
    frames[[length(frames) + 1]] <- df
    cat(sprintf("  Charg\u00e9 : %s  (%d lignes)\n", unique(df$participant)[1], nrow(df)))
  }
  all_df <- combine_frames(frames)
  all_df[order(all_df$participant, all_df$sequence, all_df$block, all_df$trial), ]
}

# ---------------------------------------------------------------------
# 2. FILTRE DE KALMAN (mod\u00e8le lin\u00e9aire x[n+1]=A*x[n]+B, MLE + RTS),
#    ajust\u00e9 s\u00e9par\u00e9ment sur chaque bloc — m\u00eame machinerie r\u00e9duite que
#    dans spat_block_transition_drops.R
# ---------------------------------------------------------------------

kf_run <- function(y, Phi, b, Q, R, x0, P0) {
  n <- length(y)
  x_pred <- numeric(n); x_filt <- numeric(n)
  P_pred <- numeric(n); P_filt <- numeric(n)
  loglik <- 0
  x_prev <- x0; P_prev <- P0
  for (t in seq_len(n)) {
    xp <- Phi * x_prev + b
    Pp <- Phi * P_prev * Phi + Q
    x_pred[t] <- xp; P_pred[t] <- Pp
    if (is.na(y[t])) {
      xf <- xp; Pf <- Pp
    } else {
      S <- Pp + R
      if (!is.finite(S) || S <= 0) S <- 1e-6
      K <- Pp / S
      e <- y[t] - xp
      xf <- xp + K * e
      Pf <- (1 - K) * Pp
      loglik <- loglik - 0.5 * (log(2 * pi * S) + e^2 / S)
    }
    x_filt[t] <- xf; P_filt[t] <- Pf
    x_prev <- xf; P_prev <- Pf
  }
  list(x_pred = x_pred, P_pred = P_pred, x_filt = x_filt, P_filt = P_filt, loglik = loglik)
}

kf_smooth <- function(kf_out, Phi) {
  n <- length(kf_out$x_filt)
  x_smooth <- numeric(n); x_smooth[n] <- kf_out$x_filt[n]
  if (n >= 2) {
    for (t in (n - 1):1) {
      J <- kf_out$P_filt[t] * Phi / kf_out$P_pred[t + 1]
      x_smooth[t] <- kf_out$x_filt[t] + J * (x_smooth[t + 1] - kf_out$x_pred[t + 1])
    }
  }
  x_smooth
}

fit_kf_linear_block <- function(y) {
  y <- as.numeric(y)
  n_valid <- sum(!is.na(y))
  if (n_valid < MIN_POINTS_SSM) return(NULL)
  y_obs <- y[!is.na(y)]
  var_y <- stats::var(y_obs)
  if (!is.finite(var_y) || var_y <= 0) var_y <- 1
  init_A <- 0.8
  init_B <- mean(diff(y_obs), na.rm = TRUE) * (1 - init_A)
  if (!is.finite(init_B)) init_B <- 0
  init_logQ <- log(max(var_y * 0.1, 1e-4))
  init_logR <- log(max(var_y * 0.5, 1e-4))
  x0 <- y_obs[1]; P0 <- var_y * 10
  negloglik <- function(par) {
    A <- par[1]; B <- par[2]; Q <- exp(par[3]); R <- exp(par[4])
    out <- tryCatch(kf_run(y, A, B, Q, R, x0, P0), error = function(e) NULL)
    if (is.null(out) || !is.finite(out$loglik)) return(1e10)
    -out$loglik
  }
  lower <- c(-0.999, -Inf, log(var_y * 1e-6), log(var_y * 0.05))
  upper <- c(0.999, Inf, log(var_y * 0.05), log(var_y * 10))
  fit <- tryCatch(
    stats::optim(c(init_A, init_B, init_logQ, init_logR), negloglik,
                method = "L-BFGS-B", lower = lower, upper = upper,
                control = list(maxit = 500)),
    error = function(e) NULL
  )
  if (is.null(fit) || !is.finite(fit$value)) {
    A <- init_A; B <- init_B; Q <- exp(init_logQ); R <- exp(init_logR)
  } else {
    A <- fit$par[1]; B <- fit$par[2]; Q <- exp(fit$par[3]); R <- exp(fit$par[4])
  }
  result <- tryCatch({
    kf_out <- kf_run(y, A, B, Q, R, x0, P0)
    list(x_smooth = kf_smooth(kf_out, A))
  }, error = function(e) NULL)
  result
}

# ---------------------------------------------------------------------
# 3. VALEUR D'UN SCORE A UN ESSAI CIBLE, PAR INTERPOLATION
# ---------------------------------------------------------------------

# Retourne la valeur interpol\u00e9e (lin\u00e9aire) de `values` \u00e0 `target_trial`,
# en fonction des positions `trials`. Pas d'extrapolation : NA si
# target_trial est hors de l'intervalle observ\u00e9 [min(trials), max(trials)].
interp_at_trial <- function(trials, values, target_trial) {
  ok <- is.finite(trials) & is.finite(values)
  trials <- trials[ok]; values <- values[ok]
  if (length(trials) < 2) return(NA_real_)
  o <- order(trials); trials <- trials[o]; values <- values[o]
  trials_u <- unique(trials)
  if (length(trials_u) < 2) return(NA_real_)
  if (target_trial < min(trials) || target_trial > max(trials)) return(NA_real_)
  as.numeric(stats::approx(trials, values, xout = target_trial, ties = mean)$y)
}

# Pour un bloc (data.frame ordonn\u00e9 par trial) : renvoie la valeur brute
# interpol\u00e9e ET la valeur liss\u00e9e (Kalman) interpol\u00e9e, \u00e0 l'essai cible.
value_at_trial_block <- function(block_df, score_col, target_trial) {
  bd <- block_df[order(block_df$trial), ]
  y <- bd[[score_col]]
  trials <- bd$trial

  raw_val <- interp_at_trial(trials, y, target_trial)

  smooth_val <- NA_real_
  fl <- tryCatch(fit_kf_linear_block(y), error = function(e) NULL)
  if (!is.null(fl)) {
    smooth_val <- interp_at_trial(trials, fl$x_smooth, target_trial)
  }

  list(raw = raw_val, smooth = smooth_val, n_points = sum(!is.na(y)),
       trial_range = if (length(trials) > 0) sprintf("%d-%d", min(trials, na.rm = TRUE), max(trials, na.rm = TRUE)) else NA_character_)
}

# ---------------------------------------------------------------------
# 4. CALCUL DE LA RECUPERATION BLOC 4 / BLOC 1 A L'ESSAI CIBLE
# ---------------------------------------------------------------------

compute_recovery_one <- function(pdata, score_col, target_trial, higher_is_better) {
  b1 <- pdata[!is.na(pdata$block) & pdata$block == 1, ]
  b4 <- pdata[!is.na(pdata$block) & pdata$block == 4, ]
  if (nrow(b1) == 0 || nrow(b4) == 0) return(NULL)

  v1 <- value_at_trial_block(b1, score_col, target_trial)
  v4 <- value_at_trial_block(b4, score_col, target_trial)

  mk_pct <- function(val1, val4) {
    if (!is.finite(val1) || !is.finite(val4) || abs(val1) < 1e-12) {
      return(list(raw_pct = NA_real_, oriented_pct = NA_real_))
    }
    raw_pct <- 100 * val4 / val1
    oriented_pct <- if (isTRUE(higher_is_better)) raw_pct else (if (abs(val4) < 1e-12) NA_real_ else 100 * val1 / val4)
    list(raw_pct = raw_pct, oriented_pct = oriented_pct)
  }

  pct_raw    <- mk_pct(v1$raw, v4$raw)
  pct_smooth <- mk_pct(v1$smooth, v4$smooth)

  list(
    v1_raw = v1$raw, v4_raw = v4$raw,
    v1_smooth = v1$smooth, v4_smooth = v4$smooth,
    recovery_pct_raw = pct_raw$raw_pct, recovery_pct_raw_oriented = pct_raw$oriented_pct,
    recovery_pct_smooth = pct_smooth$raw_pct, recovery_pct_smooth_oriented = pct_smooth$oriented_pct,
    n_points_b1 = v1$n_points, n_points_b4 = v4$n_points,
    trial_range_b1 = v1$trial_range, trial_range_b4 = v4$trial_range
  )
}

compute_all_recovery <- function(all_df, score_cols, target_trial) {
  rows <- list()
  keys <- unique(all_df[!is.na(all_df$sequence), c("participant", "sequence")])
  for (i in seq_len(nrow(keys))) {
    part <- keys$participant[i]; seq_ <- keys$sequence[i]
    pdata <- all_df[all_df$participant == part & all_df$sequence == seq_, ]
    for (score_col in score_cols) {
      if (!(score_col %in% names(pdata))) next
      meta <- SCORE_META[[score_col]]
      hib  <- if (!is.null(meta)) meta$higher_is_better else TRUE
      r <- tryCatch(compute_recovery_one(pdata, score_col, target_trial, hib), error = function(e) NULL)
      if (is.null(r)) next
      rows[[length(rows) + 1]] <- data.frame(
        participant = part, sequence = seq_, score = score_col, target_trial = target_trial,
        v1_raw = r$v1_raw, v4_raw = r$v4_raw, v1_smooth = r$v1_smooth, v4_smooth = r$v4_smooth,
        recovery_pct_raw = r$recovery_pct_raw, recovery_pct_raw_oriented = r$recovery_pct_raw_oriented,
        recovery_pct_smooth = r$recovery_pct_smooth, recovery_pct_smooth_oriented = r$recovery_pct_smooth_oriented,
        n_points_b1 = r$n_points_b1, n_points_b4 = r$n_points_b4,
        trial_range_b1 = r$trial_range_b1, trial_range_b4 = r$trial_range_b4,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------
# 5. OUTILS STATISTIQUES — comparaison A vs B (identiques \u00e0
#    spat_block_transition_drops.R : appari\u00e9 par d\u00e9faut, repli sur
#    groupes ind\u00e9pendants de Welch)
# ---------------------------------------------------------------------

paired_stats <- function(x, y, n_boot = N_BOOT, delta = NULL) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < MIN_PAIRED_N) return(NULL)

  d <- x - y
  tt <- stats::t.test(x, y, paired = TRUE)
  wc <- tryCatch(stats::wilcox.test(x, y, paired = TRUE), error = function(e) NULL)

  boot_d <- numeric(n_boot)
  for (i in seq_len(n_boot)) boot_d[i] <- mean(sample(d, n, replace = TRUE))
  ci <- stats::quantile(boot_d, c(0.025, 0.975))
  p_boot <- 2 * min(mean(boot_d <= 0), mean(boot_d >= 0))
  p_boot <- min(p_boot, 1)

  perm_d <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    signs <- sample(c(-1, 1), n, replace = TRUE)
    perm_d[i] <- mean(d * signs)
  }
  p_perm <- mean(abs(perm_d) >= abs(mean(d)))

  sd_d <- stats::sd(d)
  dz <- if (sd_d > 1e-10) mean(d) / sd_d else NA_real_
  J <- if (n > 2) 1 - 3 / (4 * (n - 1) - 1) else NA_real_
  g <- if (is.finite(dz) && is.finite(J)) dz * J else NA_real_

  if (is.null(delta)) delta <- 0.5 * sd_d
  se_d <- sd_d / sqrt(n)
  tost_p <- NA_real_
  if (is.finite(se_d) && se_d > 1e-12) {
    t_lower <- (mean(d) + delta) / se_d
    t_upper <- (mean(d) - delta) / se_d
    p_lower <- 1 - stats::pt(t_lower, df = n - 1)
    p_upper <- stats::pt(t_upper, df = n - 1)
    tost_p <- max(p_lower, p_upper)
  }

  list(design = "paired", n = n,
       mean_group1 = mean(x), mean_group2 = mean(y), mean_diff = mean(d), sd_diff = sd_d,
       t = unname(tt$statistic), t_p = tt$p.value,
       wilcox_p = if (!is.null(wc)) wc$p.value else NA_real_,
       boot_ci_low = unname(ci[1]), boot_ci_high = unname(ci[2]), boot_p = p_boot,
       perm_p = p_perm, cohens_dz = dz, hedges_g = g,
       tost_delta = delta, tost_p = tost_p)
}

cohens_d_indep <- function(a, b) {
  n1 <- length(a); n2 <- length(b)
  if (n1 < 2 || n2 < 2) return(NA_real_)
  sp <- sqrt(((n1 - 1) * stats::var(a) + (n2 - 1) * stats::var(b)) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp < 1e-12) return(NA_real_)
  (mean(a) - mean(b)) / sp
}

independent_stats <- function(a, b, n_boot = N_BOOT) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 3 || length(b) < 3) return(NULL)
  tt <- stats::t.test(a, b, var.equal = FALSE)
  boot_diff <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    boot_diff[i] <- mean(sample(a, length(a), replace = TRUE)) - mean(sample(b, length(b), replace = TRUE))
  }
  ci <- stats::quantile(boot_diff, c(0.025, 0.975))
  combined <- c(a, b); na <- length(a)
  perm_diff <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    perm <- sample(combined)
    perm_diff[i] <- mean(perm[1:na]) - mean(perm[(na + 1):length(perm)])
  }
  obs_diff <- mean(a) - mean(b)
  p_perm <- mean(abs(perm_diff) >= abs(obs_diff))
  d <- cohens_d_indep(a, b)
  list(design = "independent (Welch)", n_group1 = length(a), n_group2 = length(b),
       mean_group1 = mean(a), mean_group2 = mean(b), mean_diff = obs_diff,
       t = unname(tt$statistic), t_p = tt$p.value,
       wilcox_p = NA_real_,
       boot_ci_low = unname(ci[1]), boot_ci_high = unname(ci[2]), boot_p = NA_real_,
       perm_p = p_perm, cohens_dz = d, hedges_g = NA_real_,
       tost_delta = NA_real_, tost_p = NA_real_)
}

compare_recovery_A_vs_B <- function(recovery_df, score_col, metric_col) {
  gA <- recovery_df[recovery_df$score == score_col & recovery_df$sequence == "A", ]
  gB <- recovery_df[recovery_df$score == score_col & recovery_df$sequence == "B", ]

  common <- intersect(gA$participant, gB$participant)
  if (length(common) >= MIN_PAIRED_N) {
    gAo <- gA[match(common, gA$participant), ]
    gBo <- gB[match(common, gB$participant), ]
    res <- paired_stats(gAo[[metric_col]], gBo[[metric_col]])
    if (!is.null(res)) return(list(result = res, values1 = gAo[[metric_col]], values2 = gBo[[metric_col]], ids = common))
  }
  res <- independent_stats(gA[[metric_col]], gB[[metric_col]])
  list(result = res, values1 = gA[[metric_col]], values2 = gB[[metric_col]], ids = NULL)
}

# ---------------------------------------------------------------------
# 6. GRAPHIQUE
# ---------------------------------------------------------------------

plot_recovery_comparison <- function(cmp_out, score_col, metric_name, metric_label, out_dir) {
  v1 <- cmp_out$values1; v2 <- cmp_out$values2
  v1 <- v1[is.finite(v1)]; v2 <- v2[is.finite(v2)]
  if (length(v1) == 0 && length(v2) == 0) return(NULL)

  meta <- SCORE_META[[score_col]]
  score_label <- if (!is.null(meta)) meta$label else score_col

  fname <- file.path(out_dir, sprintf("recovery_points_%s_%s.png", score_col, metric_name))
  grDevices::png(fname, width = 900, height = 600, res = 130)

  all_vals <- c(v1, v2)
  yr <- range(all_vals, na.rm = TRUE)
  pad <- diff(yr) * 0.15
  if (!is.finite(pad) || pad == 0) pad <- 1
  yr <- yr + c(-pad, pad)

  par(mar = c(5, 4.5, 4, 2))
  plot(NA, xlim = c(0.5, 2.5), ylim = yr, xaxt = "n", xlab = "", ylab = metric_label,
       main = sprintf("%s\nR\u00e9cup\u00e9ration Bloc 4 vs Bloc 1 \u00e0 l'essai cible \u2014 %s", score_label, metric_label),
       cex.main = 0.85)
  axis(1, at = c(1, 2), labels = c("S\u00e9quence A", "S\u00e9quence B"))
  grid(nx = NA, ny = NULL, col = "gray90")
  abline(h = 100, lty = 2, col = "gray50")

  paired <- !is.null(cmp_out$ids)
  if (paired) {
    for (i in seq_along(v1)) lines(c(1, 2), c(v1[i], v2[i]), col = "gray70", lwd = 1)
  }
  points(rep(1, length(v1)), v1, pch = 21, bg = SEQ_COLORS[["A"]], col = "black", cex = 1.3)
  points(rep(2, length(v2)), v2, pch = 21, bg = SEQ_COLORS[["B"]], col = "black", cex = 1.3)
  segments(0.8, mean(v1), 1.2, mean(v1), lwd = 3)
  segments(1.8, mean(v2), 2.2, mean(v2), lwd = 3)

  res <- cmp_out$result
  if (!is.null(res)) {
    subtitle <- if (res$design == "paired") {
      sprintf("appari\u00e9, n=%d | t=%.2f, p=%.4f | dz=%.2f | perm p=%.4f", res$n, res$t, res$t_p, res$cohens_dz, res$perm_p)
    } else {
      sprintf("ind\u00e9pendant (Welch), n1=%d n2=%d | t=%.2f, p=%.4f | perm p=%.4f",
              res$n_group1, res$n_group2, res$t, res$t_p, res$perm_p)
    }
    mtext(subtitle, side = 3, line = 0.3, cex = 0.7)
  }

  grDevices::dev.off()
  cat(sprintf("    Sauvegard\u00e9 : %s\n", fname))
  fname
}

# ---------------------------------------------------------------------
# 7. PIPELINE PRINCIPAL
# ---------------------------------------------------------------------

main <- function(scores_folder, out_dir = "./output_recovery_block4",
                 target_trial = TARGET_TRIAL_DEFAULT) {
  cat(strrep("=", 70), "\n")
  cat("SPAT \u2014 R\u00c9CUP\u00c9RATION BLOC 4 vs BLOC 1 \u00c0 L'ESSAI CIBLE\n")
  cat(sprintf("         (essai cible = %d)\n", target_trial))
  cat(strrep("=", 70), "\n")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  cat(sprintf("\n[1/3] Chargement des scores depuis : %s\n", scores_folder))
  all_df <- load_participant_scores(scores_folder)
  if (nrow(all_df) == 0) return(invisible(NULL))

  score_cols <- intersect(names(SCORE_META), colnames(all_df))
  cat(sprintf("  Scores trouv\u00e9s : %s\n", paste(score_cols, collapse = ", ")))

  cat(sprintf("\n[2/3] Calcul de la r\u00e9cup\u00e9ration \u00e0 l'essai %d (brut + liss\u00e9 Kalman)\u2026\n", target_trial))
  recovery_df <- compute_all_recovery(all_df, score_cols, target_trial)
  if (nrow(recovery_df) == 0) {
    cat("  Aucune r\u00e9cup\u00e9ration calculable (donn\u00e9es insuffisantes, ou essai cible hors de port\u00e9e).\n")
    return(invisible(NULL))
  }
  recovery_path <- file.path(out_dir, "recovery_block4_vs_block1_trial20.csv")
  utils::write.csv(recovery_df, recovery_path, row.names = FALSE)
  cat(sprintf("  %d lignes -> %s\n", nrow(recovery_df), recovery_path))

  n_na_raw <- sum(!is.finite(recovery_df$recovery_pct_raw))
  if (n_na_raw > 0) {
    cat(sprintf("  Remarque : %d/%d lignes ont une r\u00e9cup\u00e9ration brute non calculable\n", n_na_raw, nrow(recovery_df)))
    cat(sprintf("             (essai %d hors de l'intervalle observ\u00e9 dans le bloc 1 ou le bloc 4).\n", target_trial))
  }

  cat("\n[3/3] Comparaison S\u00e9quence A vs S\u00e9quence B\u2026\n")
  metrics <- c(
    recovery_pct_raw_oriented    = "R\u00e9cup\u00e9ration % (brute, orient\u00e9e)",
    recovery_pct_smooth_oriented = "R\u00e9cup\u00e9ration % (liss\u00e9e Kalman, orient\u00e9e)"
  )
  stats_rows <- list()

  for (score_col in score_cols) {
    for (metric_name in names(metrics)) {
      metric_label <- metrics[[metric_name]]
      cmp_out <- tryCatch(
        compare_recovery_A_vs_B(recovery_df, score_col, metric_name),
        error = function(e) { cat(sprintf("    ERREUR %s/%s : %s\n", score_col, metric_name, conditionMessage(e))); NULL }
      )
      if (is.null(cmp_out) || is.null(cmp_out$result)) {
        cat(sprintf("  %-32s %-32s : donn\u00e9es insuffisantes\n", score_col, metric_name))
        next
      }
      res <- cmp_out$result
      cat(sprintf("  %-32s %-32s : design=%s  mean_A=%.2f  mean_B=%.2f  p=%.4f\n",
                  score_col, metric_name, res$design, res$mean_group1, res$mean_group2, res$t_p))
      stats_rows[[length(stats_rows) + 1]] <- data.frame(
        score = score_col, metric = metric_name, target_trial = target_trial, res, stringsAsFactors = FALSE
      )
      plot_recovery_comparison(cmp_out, score_col, metric_name, metric_label, out_dir)
    }
  }

  if (length(stats_rows) > 0) {
    stats_df <- combine_frames(stats_rows)
    stats_path <- file.path(out_dir, "recovery_comparison_stats.csv")
    utils::write.csv(stats_df, stats_path, row.names = FALSE)
    cat(sprintf("\n  R\u00e9sum\u00e9 statistique \u00e9crit dans : %s\n", stats_path))
  }

  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(sprintf("TERMIN\u00c9. Sorties dans : %s/\n", out_dir))
  cat(strrep("=", 70), "\n")
  invisible(list(recovery = recovery_df))
}

# ---------------------------------------------------------------------
# EXECUTION EN LIGNE DE COMMANDE
# ---------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  target_trial <- if (length(args) >= 3) as.numeric(args[3]) else TARGET_TRIAL_DEFAULT
  main(args[1], args[2], target_trial)
} else {
  main("C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_scores_harmonic", "C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_recovery_block4")
}
