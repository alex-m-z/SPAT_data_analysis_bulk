#!/usr/bin/env Rscript
# =====================================================================
# SPAT — Comparaison des "drops" (chutes) de performance aux
#         transitions de bloc, crois\u00e9es entre s\u00e9quences
# =====================================================================
#
# CONTEXTE
# --------
# On observe une chute ("drop") de performance \u00e0 certaines transitions
# de bloc. On veut comparer :
#
#   Comparaison 1 :  Drop S\u00e9quence A (bloc 1 -> bloc 2)
#                vs  Drop S\u00e9quence B (bloc 2 -> bloc 3)
#
#   Comparaison 2 :  Drop S\u00e9quence B (bloc 1 -> bloc 2)
#                vs  Drop S\u00e9quence A (bloc 2 -> bloc 3)
#
# ... et ce pour TROIS fa\u00e7ons d'estimer la valeur de "fin de bloc" et
# de "d\u00e9but de bloc" utilis\u00e9es pour calculer le drop :
#
#   (1) "avg"  : moyenne brute des derniers/premiers essais du bloc
#                (fen\u00eatre de EDGE_N essais, robuste au bruit ponctuel)
#   (2) "sqrt" : mod\u00e8le de tendance y[n] \u2248 a - b*sqrt(n), ajust\u00e9
#                S\u00c9PAR\u00c9MENT sur chaque bloc (comme dans
#                spat_learning_rate_comparison.R), \u00e9valu\u00e9 au dernier
#                essai du bloc pr\u00e9c\u00e9dent et au premier essai du bloc
#                suivant
#   (3) "ssm"  : mod\u00e8le \u00e9tat-espace lin\u00e9aire x[n+1] = A*x[n] + B,
#                ajust\u00e9 par filtre de Kalman (MLE) + lisseur RTS,
#                S\u00c9PAR\u00c9MENT sur chaque bloc (comme dans
#                spat_plot_evolution.R), \u00e9valu\u00e9 de la m\u00eame fa\u00e7on
#
# Le drop est d\u00e9fini comme :
#     drop = valeur_estim\u00e9e(d\u00e9but du bloc suivant)
#          - valeur_estim\u00e9e(fin du bloc pr\u00e9c\u00e9dent)
#
# Un drop POSITIF signifie donc une AUGMENTATION du score \u00e0 la
# transition ; pour des scores "co\u00fbt" (dtw_dist, time_spent,
# mean_dev, plus petit = meilleur), un drop positif est une
# D\u00c9GRADATION de performance. Pour un score "b\u00e9n\u00e9fice"
# (score_harmonic, plus grand = meilleur), un drop positif est une
# AM\u00c9LIORATION. Ce script ne r\u00e9interpr\u00e8te pas le signe : il
# rapporte le drop brut pour chaque score, et laisse l'utilisateur
# consulter SCORE_META$higher_is_better pour l'interpr\u00e9tation.
#
# Les deux groupes compar\u00e9s (ex: A:1-2 vs B:2-3) partagent en g\u00e9n\u00e9ral
# les M\u00caMES participants (chaque participant fait les deux s\u00e9quences),
# donc la comparaison est APPARI\u00c9E par d\u00e9faut (appariement sur
# 'participant'). Si les participants ne se recoupent pas suffisamment,
# le script bascule automatiquement sur un test de groupes ind\u00e9pendants
# (Welch) et l'indique dans la sortie.
#
# SORTIES (dans out_dir) :
#   drops_all.csv                     toutes les valeurs de drop calcul\u00e9es
#                                      (participant x s\u00e9quence x score x
#                                      transition x m\u00e9thode)
#   drop_comparisons_stats.csv        r\u00e9sultats statistiques des deux
#                                      comparaisons crois\u00e9es, par score
#                                      et par m\u00e9thode
#   drop_points_<comparaison>_<score>_<methode>.png
#                                      nuage de points appari\u00e9s (1 par
#                                      participant, reli\u00e9 par une ligne)
#
# USAGE
# -----
#   Rscript spat_block_transition_drops.R scores_folder out_dir
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

EDGE_N <- 3          # nb d'essais utilis\u00e9s pour la m\u00e9thode "avg"
N_BOOT <- 10000       # nb de r\u00e9\u00e9chantillonnages bootstrap / permutation
MIN_POINTS_SQRT <- 4  # nb minimum de points pour ajuster le mod\u00e8le sqrt sur un bloc
MIN_POINTS_SSM  <- 5  # nb minimum de points pour ajuster le mod\u00e8le SSM (Kalman) sur un bloc
MIN_PAIRED_N    <- 4  # nb minimum de participants appari\u00e9s pour un test appari\u00e9

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

# Les deux comparaisons crois\u00e9es demand\u00e9es par l'utilisateur :
#   nom            groupe 1 (sequence, transition)   groupe 2 (sequence, transition)
COMPARISONS <- list(
  list(name = "A(1-2)_vs_B(2-3)", seq1 = "A", tr1 = "1-2", seq2 = "B", tr2 = "2-3"),
  list(name = "B(1-2)_vs_A(2-3)", seq1 = "B", tr1 = "1-2", seq2 = "A", tr2 = "2-3")
)

METHODS <- c("avg", "sqrt", "ssm")
METHOD_LABELS <- c(avg = "Moyenne brute (fen\u00eatre de bord)",
                   sqrt = "Mod\u00e8le \u221an (par bloc)",
                   ssm  = "Mod\u00e8le \u00e9tat-espace (Kalman, par bloc)")

# ---------------------------------------------------------------------
# 1. CHARGEMENT DES DONNEES (identique \u00e0 spat_plot_evolution.R)
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
# 2. FILTRE DE KALMAN (mod\u00e8le lin\u00e9aire, m\u00eame machinerie que dans
#    spat_plot_evolution.R, r\u00e9duite au seul mod\u00e8le (1) dont on a
#    besoin ici)
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

# Ajuste x[n+1] = A*x[n] + B par MLE (Kalman), retourne la s\u00e9rie
# LISS\u00c9E (RTS) sur le bloc fourni. Retourne NULL si trop peu de points
# ou en cas d'\u00e9chec num\u00e9rique.
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
# 3. MODELE "RACINE CARREE" PAR BLOC  y[n] \u2248 a - b*sqrt(n)
#    (lin\u00e9aire en sqrt(n) -> ajust\u00e9 par lm(), comme dans
#    spat_learning_rate_comparison.R)
# ---------------------------------------------------------------------

fit_sqrt_block <- function(y) {
  y <- as.numeric(y)
  ok <- !is.na(y)
  if (sum(ok) < MIN_POINTS_SQRT) return(NULL)
  n <- length(y)
  x <- seq_len(n)
  fit <- tryCatch(stats::lm(y[ok] ~ sqrt(x[ok])), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  co <- stats::coef(fit)
  a <- unname(co[1]); b <- unname(co[2])
  if (!is.finite(a) || !is.finite(b)) return(NULL)
  list(predict = function(xx) a + b * sqrt(xx), a = a, b = b)
}

# ---------------------------------------------------------------------
# 4. ESTIMATEURS DE VALEUR DE BORD DE BLOC ("fin de bloc pr\u00e9c\u00e9dent" /
#    "d\u00e9but de bloc suivant")
# ---------------------------------------------------------------------

# y doit \u00eatre tri\u00e9 par ordre chronologique (trial croissant) DANS le bloc.

boundary_avg <- function(y, edge = "end", edge_n = EDGE_N) {
  y <- as.numeric(y); y <- y[!is.na(y)]
  if (length(y) == 0) return(NA_real_)
  k <- min(edge_n, length(y))
  if (edge == "end") mean(utils::tail(y, k)) else mean(utils::head(y, k))
}

boundary_sqrt <- function(y, edge = "end") {
  fs <- fit_sqrt_block(y)
  if (is.null(fs)) return(NA_real_)
  n <- length(y)
  if (edge == "end") fs$predict(n) else fs$predict(1)
}

boundary_ssm <- function(y, edge = "end") {
  fl <- fit_kf_linear_block(y)
  if (is.null(fl)) return(NA_real_)
  if (edge == "end") utils::tail(fl$x_smooth, 1) else utils::head(fl$x_smooth, 1)
}

# ---------------------------------------------------------------------
# 5. CALCUL DES DROPS
# ---------------------------------------------------------------------

# Calcule le drop (fin bloc_start -> d\u00e9but bloc_end) pour un
# participant x s\u00e9quence x score donn\u00e9, pour les 3 m\u00e9thodes.
compute_drop_one <- function(pdata, score_col, block_start, block_end) {
  bs <- pdata[!is.na(pdata$block) & pdata$block == block_start, ]
  be <- pdata[!is.na(pdata$block) & pdata$block == block_end, ]
  bs <- bs[order(bs$trial), ]
  be <- be[order(be$trial), ]
  y_bs <- bs[[score_col]]
  y_be <- be[[score_col]]

  avg_end   <- boundary_avg(y_bs, "end")
  avg_start <- boundary_avg(y_be, "start")
  sqrt_end   <- boundary_sqrt(y_bs, "end")
  sqrt_start <- boundary_sqrt(y_be, "start")
  ssm_end   <- boundary_ssm(y_bs, "end")
  ssm_start <- boundary_ssm(y_be, "start")

  list(
    drop_avg  = avg_start - avg_end,
    drop_sqrt = sqrt_start - sqrt_end,
    drop_ssm  = ssm_start - ssm_end,
    n_bs = sum(!is.na(y_bs)), n_be = sum(!is.na(y_be))
  )
}

# Calcule tous les drops pour toutes les transitions n\u00e9cessaires
# ("1-2" et "2-3"), pour tous les scores, participants, s\u00e9quences.
compute_all_drops <- function(all_df, score_cols) {
  transitions <- list("1-2" = c(1, 2), "2-3" = c(2, 3))
  rows <- list()
  keys <- unique(all_df[!is.na(all_df$sequence), c("participant", "sequence")])
  for (i in seq_len(nrow(keys))) {
    part <- keys$participant[i]; seq_ <- keys$sequence[i]
    pdata <- all_df[all_df$participant == part & all_df$sequence == seq_, ]
    for (score_col in score_cols) {
      if (!(score_col %in% names(pdata))) next
      for (tr_name in names(transitions)) {
        bp <- transitions[[tr_name]]
        d <- tryCatch(compute_drop_one(pdata, score_col, bp[1], bp[2]), error = function(e) NULL)
        if (is.null(d)) next
        rows[[length(rows) + 1]] <- data.frame(
          participant = part, sequence = seq_, score = score_col, transition = tr_name,
          drop_avg = d$drop_avg, drop_sqrt = d$drop_sqrt, drop_ssm = d$drop_ssm,
          n_block_start = d$n_bs, n_block_end = d$n_be, stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------
# 6. OUTILS STATISTIQUES (compare deux vecteurs de drops)
#    Appari\u00e9 par d\u00e9faut (sur 'participant') ; repli automatique sur
#    un test de groupes ind\u00e9pendants (Welch) si l'appariement \u00e9choue.
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

# Compare le groupe 1 (seq1, transition tr1) au groupe 2 (seq2, tr2)
# pour un score et une m\u00e9thode donn\u00e9s. Essaie l'appariement par
# participant en premier ; repli sur Welch si l'appariement n'a pas
# assez de sujets communs.
compare_drop_groups <- function(drops_df, score_col, method_col, seq1, tr1, seq2, tr2) {
  g1 <- drops_df[drops_df$score == score_col & drops_df$sequence == seq1 & drops_df$transition == tr1, ]
  g2 <- drops_df[drops_df$score == score_col & drops_df$sequence == seq2 & drops_df$transition == tr2, ]

  common <- intersect(g1$participant, g2$participant)
  if (length(common) >= MIN_PAIRED_N) {
    g1o <- g1[match(common, g1$participant), ]
    g2o <- g2[match(common, g2$participant), ]
    res <- paired_stats(g1o[[method_col]], g2o[[method_col]])
    if (!is.null(res)) return(list(result = res, values1 = g1o[[method_col]], values2 = g2o[[method_col]], ids = common))
  }
  # Repli : groupes ind\u00e9pendants
  res <- independent_stats(g1[[method_col]], g2[[method_col]])
  list(result = res, values1 = g1[[method_col]], values2 = g2[[method_col]], ids = NULL)
}

# ---------------------------------------------------------------------
# 7. GRAPHIQUE — nuage de points appari\u00e9s (ou groupes ind\u00e9pendants)
# ---------------------------------------------------------------------

plot_drop_comparison <- function(cmp_out, comparison_name, score_col, method_name,
                                 label1, label2, col1, col2, out_dir) {
  v1 <- cmp_out$values1; v2 <- cmp_out$values2
  v1 <- v1[is.finite(v1)]; v2 <- v2[is.finite(v2)]
  if (length(v1) == 0 && length(v2) == 0) return(NULL)

  meta <- SCORE_META[[score_col]]
  score_label <- if (!is.null(meta)) meta$label else score_col

  fname <- file.path(out_dir, sprintf("drop_points_%s_%s_%s.png", comparison_name, score_col, method_name))
  grDevices::png(fname, width = 900, height = 600, res = 130)

  all_vals <- c(v1, v2)
  yr <- range(all_vals, na.rm = TRUE)
  pad <- diff(yr) * 0.15
  if (!is.finite(pad) || pad == 0) pad <- 1
  yr <- yr + c(-pad, pad)

  par(mar = c(5, 4.5, 4, 2))
  plot(NA, xlim = c(0.5, 2.5), ylim = yr, xaxt = "n", xlab = "",
       ylab = sprintf("Drop (%s)", METHOD_LABELS[[method_name]]),
       main = sprintf("%s\n%s \u2014 M\u00e9thode : %s", score_label,
                      gsub("_", " ", comparison_name), METHOD_LABELS[[method_name]]),
       cex.main = 0.85)
  axis(1, at = c(1, 2), labels = c(label1, label2))
  grid(nx = NA, ny = NULL, col = "gray90")
  abline(h = 0, lty = 2, col = "gray50")

  paired <- !is.null(cmp_out$ids)
  if (paired) {
    for (i in seq_along(v1)) lines(c(1, 2), c(v1[i], v2[i]), col = "gray70", lwd = 1)
  }
  points(rep(1, length(v1)), v1, pch = 21, bg = col1, col = "black", cex = 1.3)
  points(rep(2, length(v2)), v2, pch = 21, bg = col2, col = "black", cex = 1.3)
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

draw_violin <- function(x, at, col, width = 0.20, border = "black") {
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

sig_stars <- function(p) {
  if (is.na(p)) return("n.s.")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("n.s.")
}
plot_drop_violin <- function(cmp_out, comparison_name, score_col, method_name,
                             label1, label2, col1, col2, out_dir) {
  v1 <- cmp_out$values1; v2 <- cmp_out$values2
  v1 <- v1[is.finite(v1)]; v2 <- v2[is.finite(v2)]
  if (length(v1) == 0 && length(v2) == 0) return(NULL)
  
  meta <- SCORE_META[[score_col]]
  score_label <- if (!is.null(meta)) meta$label else score_col
  
  fname <- file.path(out_dir, sprintf("drop_violin_%s_%s_%s.png", comparison_name, score_col, method_name))
  grDevices::png(fname, width = 900, height = 600, res = 130)
  
  all_vals <- c(v1, v2)
  yr <- range(all_vals, na.rm = TRUE)
  pad <- diff(yr) * 0.15
  if (!is.finite(pad) || pad == 0) pad <- 1
  yr <- yr + c(-pad, pad)
  
  par(mar = c(5, 4.5, 4, 2))
  plot(NA, xlim = c(0.5, 2.5), ylim = yr, xaxt = "n", xlab = "",
       ylab = sprintf("Drop (%s)", METHOD_LABELS[[method_name]]),
       main = sprintf("%s\n%s \u2014 M\u00e9thode : %s", score_label,
                      gsub("_", " ", comparison_name), METHOD_LABELS[[method_name]]),
       cex.main = 0.85)
  axis(1, at = c(1, 2), labels = c(label1, label2))
  grid(nx = NA, ny = NULL, col = "gray90")
  abline(h = 0, lty = 2, col = "gray50")
  
  draw_violin(v1, at = 1, col = col1)
  draw_violin(v2, at = 2, col = col2)
  
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
# 8. PIPELINE PRINCIPAL
# ---------------------------------------------------------------------

main <- function(scores_folder, out_dir = "./output_drop_comparisons") {
  cat(strrep("=", 70), "\n")
  cat("SPAT \u2014 COMPARAISON DES DROPS DE TRANSITION DE BLOC (CROIS\u00c9E)\n")
  cat(strrep("=", 70), "\n")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  cat(sprintf("\n[1/3] Chargement des scores depuis : %s\n", scores_folder))
  all_df <- load_participant_scores(scores_folder)
  if (nrow(all_df) == 0) return(invisible(NULL))

  score_cols <- intersect("score_harmonic_norm_participant", colnames(all_df))
  cat(sprintf("  Scores trouv\u00e9s : %s\n", paste(score_cols, collapse = ", ")))

  cat("\n[2/3] Calcul des drops (bloc1->bloc2 et bloc2->bloc3), 3 m\u00e9thodes\u2026\n")
  drops_df <- compute_all_drops(all_df, score_cols)
  if (nrow(drops_df) == 0) {
    cat("  Aucun drop calculable (donn\u00e9es insuffisantes).\n")
    return(invisible(NULL))
  }
  drops_path <- file.path(out_dir, "drops_all.csv")
  utils::write.csv(drops_df, drops_path, row.names = FALSE)
  cat(sprintf("  %d lignes -> %s\n", nrow(drops_df), drops_path))

  cat("\n[3/3] Comparaisons crois\u00e9es (A:1-2 vs B:2-3, et B:1-2 vs A:2-3)\u2026\n")
  method_cols <- c(avg = "drop_avg", sqrt = "drop_sqrt", ssm = "drop_ssm")
  stats_rows <- list()

  for (cmp in COMPARISONS) {
    label1 <- sprintf("S\u00e9q. %s, bloc %s\u2192%s", cmp$seq1, substr(cmp$tr1, 1, 1), substr(cmp$tr1, 3, 3))
    label2 <- sprintf("S\u00e9q. %s, bloc %s\u2192%s", cmp$seq2, substr(cmp$tr2, 1, 1), substr(cmp$tr2, 3, 3))
    col1 <- SEQ_COLORS[[cmp$seq1]]; col2 <- SEQ_COLORS[[cmp$seq2]]

    for (score_col in score_cols) {
      for (method_name in METHODS) {
        method_col <- method_cols[[method_name]]
        cmp_out <- tryCatch(
          compare_drop_groups(drops_df, score_col, method_col, cmp$seq1, cmp$tr1, cmp$seq2, cmp$tr2),
          error = function(e) { cat(sprintf("    ERREUR %s/%s/%s : %s\n", cmp$name, score_col, method_name, conditionMessage(e))); NULL }
        )
        if (is.null(cmp_out) || is.null(cmp_out$result)) {
          cat(sprintf("  %-18s %-32s %-5s : donn\u00e9es insuffisantes\n", cmp$name, score_col, method_name))
          next
        }
        res <- cmp_out$result
        cat(sprintf("  %-18s %-32s %-5s : design=%s  mean_diff=%.4f  p=%.4f\n",
                    cmp$name, score_col, method_name, res$design, res$mean_diff, res$t_p))
        stats_rows[[length(stats_rows) + 1]] <- data.frame(
          comparison = cmp$name, score = score_col, method = method_name, c(res), stringsAsFactors = FALSE
        )
        plot_drop_comparison(cmp_out, cmp$name, score_col, method_name, label1, label2, col1, col2, out_dir)
      }
    }
  }

  if (length(stats_rows) > 0) {
    stats_df <- combine_frames(stats_rows)
    stats_path <- file.path(out_dir, "drop_comparisons_stats.csv")
    utils::write.csv(stats_df, stats_path, row.names = FALSE)
    cat(sprintf("\n  R\u00e9sum\u00e9 statistique \u00e9crit dans : %s\n", stats_path))
  }
  
  plot_drop_comparison(cmp_out, cmp$name, score_col, method_name, label1, label2, col1, col2, out_dir)
  plot_drop_violin(cmp_out, cmp$name, score_col, method_name, label1, label2, col1, col2, out_dir)

  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(sprintf("TERMIN\u00c9. Sorties dans : %s/\n", out_dir))
  cat(strrep("=", 70), "\n")
  invisible(list(drops = drops_df))
}

# ---------------------------------------------------------------------
# EXECUTION EN LIGNE DE COMMANDE
# ---------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  main(args[1], args[2])
} else {
  main("C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_scores_harmonic", "C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_drop_comparisons")
}
