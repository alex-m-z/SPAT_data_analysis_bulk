#!/usr/bin/env Rscript
# =====================================================================
# SPAT — Évolution (2 modèles superposés) + violons Bloc 3 (Séq. A vs B)
#         + violons Bloc 4 normalisé par Bloc 1 (appariés par participant)
# =====================================================================
# Score utilisé : score_harmonic_norm_participant (uniquement).
#
# Ce script extrait/assemble 3 figures :
#
#   1) evolution_comparison_<score>_AvsB.png
#        Reprise EXACTE de plot_evolution_group_comparison() du script
#        complet : les 2 modèles (état-espace linéaire ET tendance
#        √trial, ajustés par filtre de Kalman) superposés sur le MÊME
#        graphique, avec bandes de confiance, pour Séquence A et B.
#
#   2) violin_learning_rates_<score>_block3.png
#        Reprise EXACTE de plot_violin_learning_rates() du script
#        complet, restreinte au Bloc 3 : violons "Max slope (√-model)"
#        et "A parameter (retention)", Séquence A vs Séquence B.
#
#   3) violin_block4_over_block1_<score>.png   [NOUVEAU]
#        En suivant la logique APPARIÉE PAR PARTICIPANT du script de
#        comparaison Bloc 4 (spat_learning_rate_comparison_paired.R) :
#        pour chaque participant x séquence, on calcule le ratio
#        Bloc4/Bloc1 de lr_sqrt_slope et lr_ssm_A (participant non
#        présent dans les deux blocs -> NA ; ratio instable si la
#        valeur du Bloc 1 est proche de 0 -> NA, même logique que
#        compute_learning_rate_ratios() pour le ratio A/B). Le
#        résultat est tracé en violon, Séquence A vs B, avec une
#        ligne de référence à 1 (= pas de changement Bloc1->Bloc4).
#
# SORTIES (dans out_dir) :
#   evolution_comparison_score_harmonic_norm_participant_AvsB.png
#   violin_learning_rates_score_harmonic_norm_participant_block3.png
#   violin_block4_over_block1_score_harmonic_norm_participant.png
#   learning_rates_all_blocks.csv
#   block4_over_block1_ratios.csv
#
# USAGE
# -----
#   Rscript spat_evolution_and_block_ratio_violins.R scores_folder out_dir
#
# Packages requis : readxl
# =====================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(tools)
})

# ---------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------

SEQ_COLORS <- c(A = "#4C8CC5", B = "#F1BC18")

SCORE_META <- list(
  score_harmonic_norm_participant = list(
    label  = "Harmonic score (normalized by participant)",
    ylabel = "Harmonic score (normalized by participant)"
  )
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

  all_df <- all_df[order(all_df$participant, all_df$sequence, all_df$block, all_df$trial), ]
  grp    <- interaction(all_df$participant, all_df$sequence, drop = TRUE)
  all_df$trial_global <- ave(seq_len(nrow(all_df)), grp, FUN = seq_along)

  all_df
}

get_block_boundaries <- function(sub) {
  blocks <- sort(unique(sub$block))
  blocks <- blocks[!is.na(blocks)]
  if (length(blocks) <= 1) return(numeric(0))

  boundaries <- c()
  for (b in blocks[-length(blocks)]) {
    bd  <- sub[sub$block == b, ]
    mx  <- tapply(bd$trial_global, bd$participant, max)
    boundaries <- c(boundaries, stats::median(mx, na.rm = TRUE))
  }
  boundaries
}

# ---------------------------------------------------------------------
# 2. FILTRE DE KALMAN GENERIQUE + LISSEUR RTS
# ---------------------------------------------------------------------

kf_run <- function(y, Phi, b, Q, H_fun, R, x0, P0) {
  n <- length(y)
  d <- length(x0)

  x_pred <- matrix(NA_real_, n, d)
  x_filt <- matrix(NA_real_, n, d)
  P_pred <- vector("list", n)
  P_filt <- vector("list", n)
  loglik <- 0

  x_prev <- x0
  P_prev <- P0

  for (t in seq_len(n)) {
    xp <- as.numeric(Phi %*% x_prev + b)
    Pp <- Phi %*% P_prev %*% t(Phi) + Q
    x_pred[t, ] <- xp
    P_pred[[t]] <- Pp

    if (is.na(y[t])) {
      xf <- xp
      Pf <- Pp
    } else {
      Ht   <- H_fun(t)
      yhat <- as.numeric(Ht %*% xp)
      S    <- as.numeric(Ht %*% Pp %*% t(Ht)) + R
      if (!is.finite(S) || S <= 0) S <- 1e-6

      K <- (Pp %*% t(Ht)) / S
      e <- y[t] - yhat

      xf <- xp + as.numeric(K * e)
      Pf <- (diag(d) - K %*% Ht) %*% Pp

      loglik <- loglik - 0.5 * (log(2 * pi * S) + e^2 / S)
    }
    x_filt[t, ] <- xf
    P_filt[[t]] <- Pf

    x_prev <- xf
    P_prev <- Pf
  }

  list(x_pred = x_pred, P_pred = P_pred, x_filt = x_filt, P_filt = P_filt, loglik = loglik)
}

kf_smooth <- function(kf_out, Phi) {
  n <- nrow(kf_out$x_filt)
  d <- ncol(kf_out$x_filt)

  x_smooth <- matrix(NA_real_, n, d)
  P_smooth <- vector("list", n)

  x_smooth[n, ] <- kf_out$x_filt[n, ]
  P_smooth[[n]] <- kf_out$P_filt[[n]]

  if (n >= 2) {
    for (t in (n - 1):1) {
      Pf      <- kf_out$P_filt[[t]]
      Pp_next <- kf_out$P_pred[[t + 1]]
      J       <- Pf %*% t(Phi) %*% solve(Pp_next)

      x_smooth[t, ] <- kf_out$x_filt[t, ] +
        as.numeric(J %*% (x_smooth[t + 1, ] - kf_out$x_pred[t + 1, ]))
      P_smooth[[t]] <- Pf + J %*% (P_smooth[[t + 1]] - Pp_next) %*% t(J)
    }
  }

  list(x_smooth = x_smooth, P_smooth = P_smooth)
}

# ---------------------------------------------------------------------
# 3. MODELE (1) — ETAT-ESPACE LINEAIRE  x[n+1] = A*x[n] + B
# ---------------------------------------------------------------------

fit_kf_linear <- function(y) {
  y <- as.numeric(y)
  n_valid <- sum(!is.na(y))
  if (n_valid < 5) return(NULL)

  y_obs <- y[!is.na(y)]
  var_y <- stats::var(y_obs)
  if (!is.finite(var_y) || var_y <= 0) var_y <- 1

  init_A     <- 0.8
  init_B     <- mean(diff(y_obs), na.rm = TRUE) * (1 - init_A)
  if (!is.finite(init_B)) init_B <- 0
  init_logQ  <- log(max(var_y * 0.1, 1e-4))
  init_logR  <- log(max(var_y * 0.5, 1e-4))

  x0 <- c(y_obs[1])
  P0 <- matrix(var_y * 10, 1, 1)
  H_fun <- function(t) matrix(1, 1, 1)

  negloglik <- function(par) {
    A <- par[1]; B <- par[2]; Q <- exp(par[3]); R <- exp(par[4])
    Phi <- matrix(A, 1, 1)
    out <- tryCatch(kf_run(y, Phi, B, matrix(Q, 1, 1), H_fun, R, x0, P0),
                    error = function(e) NULL)
    if (is.null(out) || !is.finite(out$loglik)) return(1e10)
    -out$loglik
  }

  lower <- c(-0.999, -Inf, log(var_y * 1e-6), log(var_y * 0.05))
  upper <- c( 0.999,  Inf, log(var_y * 0.05), log(var_y * 10))

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

  Phi    <- matrix(A, 1, 1)
  result <- tryCatch({
    kf_out <- kf_run(y, Phi, B, matrix(Q, 1, 1), H_fun, R, x0, P0)
    sm     <- kf_smooth(kf_out, Phi)
    se_vec <- sqrt(vapply(sm$P_smooth, function(P) P[1, 1], numeric(1)))
    list(A = A, B = B, Q = Q, R = R,
         x_smooth = as.numeric(sm$x_smooth[, 1]),
         se       = se_vec)
  }, error = function(e) NULL)

  result
}

# ---------------------------------------------------------------------
# 4. MODELE (2) — TENDANCE "RACINE CARREE"  y[n] ~= A_n*sqrt(n) + B_n
# ---------------------------------------------------------------------

fit_kf_sqrt <- function(y, trial_n) {
  y <- as.numeric(y)
  trial_n <- as.numeric(trial_n)
  ok <- !is.na(y)
  n_valid <- sum(ok)
  if (n_valid < 5) return(NULL)

  lm_fit <- tryCatch(stats::lm(y[ok] ~ sqrt(trial_n[ok])), error = function(e) NULL)
  if (is.null(lm_fit)) {
    B0 <- mean(y[ok]); A0 <- 0
  } else {
    B0 <- unname(stats::coef(lm_fit)[1])
    A0 <- unname(stats::coef(lm_fit)[2])
    if (!is.finite(A0)) A0 <- 0
    if (!is.finite(B0)) B0 <- mean(y[ok])
  }

  var_y <- stats::var(y[ok])
  if (!is.finite(var_y) || var_y <= 0) var_y <- 1

  init_logQa <- log(max(var_y * 0.001, 1e-6))
  init_logQb <- log(max(var_y * 0.001, 1e-6))
  init_logR  <- log(max(var_y * 0.5, 1e-4))

  Phi   <- diag(2)
  b     <- c(0, 0)
  x0    <- c(A0, B0)
  P0    <- diag(c(var_y * 10, var_y * 10))
  H_fun <- function(t) matrix(c(sqrt(trial_n[t]), 1), nrow = 1)

  negloglik <- function(par) {
    Qa <- exp(par[1]); Qb <- exp(par[2]); R <- exp(par[3])
    Qm <- diag(c(Qa, Qb))
    out <- tryCatch(kf_run(y, Phi, b, Qm, H_fun, R, x0, P0), error = function(e) NULL)
    if (is.null(out) || !is.finite(out$loglik)) return(1e10)
    -out$loglik
  }

  lower <- c(log(var_y * 1e-10), log(var_y * 1e-10), log(var_y * 0.05))
  upper <- c(log(var_y * 1e-4),  log(var_y * 1e-4),  log(var_y * 10))
  fit <- tryCatch(
    stats::optim(c(init_logQa, init_logQb, init_logR), negloglik,
                 method = "L-BFGS-B", lower = lower, upper = upper,
                 control = list(maxit = 500)),
    error = function(e) NULL
  )

  if (is.null(fit) || !is.finite(fit$value)) {
    Qa <- exp(init_logQa); Qb <- exp(init_logQb); R <- exp(init_logR)
  } else {
    Qa <- exp(fit$par[1]); Qb <- exp(fit$par[2]); R <- exp(fit$par[3])
  }

  Qm     <- diag(c(Qa, Qb))
  result <- tryCatch({
    kf_out <- kf_run(y, Phi, b, Qm, H_fun, R, x0, P0)
    sm     <- kf_smooth(kf_out, Phi)
    n <- length(y)
    fitted_mean <- vapply(seq_len(n), function(t) as.numeric(H_fun(t) %*% sm$x_smooth[t, ]), numeric(1))
    se_vec <- vapply(seq_len(n), function(t) {
      Ht <- H_fun(t)
      sqrt(as.numeric(Ht %*% sm$P_smooth[[t]] %*% t(Ht)))
    }, numeric(1))
    list(Qa = Qa, Qb = Qb, R = R,
         x_smooth = fitted_mean,
         se       = se_vec)
  }, error = function(e) NULL)

  result
}

# ---------------------------------------------------------------------
# Variantes "par bloc" : réinitialisent l'état au début de CHAQUE bloc.
# ---------------------------------------------------------------------

fit_kf_linear_by_block <- function(y, block) {
  n <- length(y)
  x_smooth <- rep(NA_real_, n)
  se       <- rep(NA_real_, n)
  for (b in sort(unique(block))) {
    idx <- which(block == b)
    fl  <- fit_kf_linear(y[idx])
    if (is.null(fl)) next
    x_smooth[idx] <- fl$x_smooth
    se[idx]       <- fl$se
  }
  list(x_smooth = x_smooth, se = se)
}

fit_kf_sqrt_by_block <- function(y, block) {
  n <- length(y)
  x_smooth <- rep(NA_real_, n)
  se       <- rep(NA_real_, n)
  for (b in sort(unique(block))) {
    idx     <- which(block == b)
    n_local <- seq_along(idx)
    fs <- fit_kf_sqrt(y[idx], n_local)
    if (is.null(fs)) next
    x_smooth[idx] <- fs$x_smooth
    se[idx]       <- fs$se
  }
  list(x_smooth = x_smooth, se = se)
}

# ---------------------------------------------------------------------
# 5. VIOLON "MAISON" (densité miroir + médiane/IQR)
# ---------------------------------------------------------------------

draw_violin <- function(x, at, col, width = 0.35, border = "black") {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(invisible(NULL))
  if (stats::sd(x) < 1e-8) {
    segments(at - width, x[1], at + width, x[1], col = border, lwd = 2)
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

# ---------------------------------------------------------------------
# 6. LEARNING RATES PAR BLOC (score x participant x sequence x block)
# ---------------------------------------------------------------------

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

simulate_ssm_r <- function(y0, A, B, n) {
  out <- numeric(n)
  out[1] <- y0
  if (n > 1) for (i in 2:n) out[i] <- A * out[i - 1] + B
  out
}

# Une ligne par (score, participant, sequence, block) : les 3 learning
# rates (dont lr_sqrt_slope et lr_ssm_A) pour CHAQUE bloc présent dans
# les données (pas seulement le Bloc 3), afin de pouvoir ensuite
# calculer le ratio Bloc4/Bloc1.
compute_learning_rates_all <- function(df, score_cols) {
  rows <- list()
  for (score_col in score_cols) {
    sub_all <- df[!is.na(df[[score_col]]) & !is.na(df$sequence) & !is.na(df$block), ]
    for (part in sort(unique(sub_all$participant))) {
      for (seq_ in sort(unique(sub_all$sequence))) {
        pdata <- sub_all[sub_all$participant == part & sub_all$sequence == seq_, ]
        for (b in sort(unique(pdata$block))) {
          bd <- pdata[pdata$block == b, ]
          bd <- bd[order(bd$trial), ]
          if (nrow(bd) < 5) next

          y  <- bd[[score_col]]
          fl <- fit_kf_linear(y)
          if (is.null(fl)) next
          y_ssm <- simulate_ssm_r(y[1], fl$A, fl$B, length(y))

          fs <- fit_kf_sqrt(y, seq_along(y))
          lr_sqrt_slope <- if (!is.null(fs)) learning_rate_max_slope(fs$x_smooth) else NA_real_

          rows[[length(rows) + 1]] <- data.frame(
            score = score_col, participant = part, sequence = seq_, block = b,
            lr_ssm_1090   = learning_rate_10_90(y_ssm, bd$trial),
            lr_sqrt_slope = lr_sqrt_slope,
            lr_ssm_slope  = learning_rate_max_slope(y_ssm),
            lr_ssm_A      = fl$A,
            lr_ssm_B      = fl$B,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------
# 7. MOYENNE DE GROUPE + SEM PAR POSITION CHRONOLOGIQUE
# ---------------------------------------------------------------------

compute_group_mean_sem <- function(sub, score_col) {
  agg_mean <- stats::aggregate(sub[[score_col]], by = list(trial_global = sub$trial_global),
                               FUN = mean, na.rm = TRUE)
  names(agg_mean)[2] <- "mean"
  agg_sd <- stats::aggregate(sub[[score_col]], by = list(trial_global = sub$trial_global),
                             FUN = stats::sd, na.rm = TRUE)
  names(agg_sd)[2] <- "sd"
  agg_n <- stats::aggregate(sub[[score_col]], by = list(trial_global = sub$trial_global),
                            FUN = function(v) sum(!is.na(v)))
  names(agg_n)[2] <- "n"
  gm_block <- stats::aggregate(block ~ trial_global, data = sub,
                               FUN = function(s) as.numeric(names(sort(table(s), decreasing = TRUE))[1]))

  out <- merge(agg_mean, agg_sd, by = "trial_global")
  out <- merge(out, agg_n,  by = "trial_global")
  out <- merge(out, gm_block, by = "trial_global")
  out$sem <- out$sd / sqrt(pmax(out$n, 1))
  out[order(out$trial_global), ]
}

# =======================================================================
# FIGURE 1 — Comparaison A vs B : moyenne de groupe + LES 2 MODELES
#            (état-espace linéaire + tendance √) SUPERPOSES + IC
#            (reprise exacte de plot_evolution_group_comparison())
# =======================================================================

plot_evolution_group_comparison <- function(df, score_col, out_dir, ci_z = 1.96) {
  meta <- SCORE_META[[score_col]]
  if (is.null(meta)) meta <- list(label = score_col, ylabel = score_col)

  sequences <- intersect(sort(unique(stats::na.omit(df$sequence))), names(SEQ_COLORS))
  if (length(sequences) < 2) return(NULL)

  gm_list <- list(); fit_lin_list <- list(); fit_sqrt_list <- list(); boundaries_list <- list()
  npart_list <- list()

  for (seq_ in sequences) {
    sub <- df[!is.na(df$sequence) & df$sequence == seq_ & !is.na(df[[score_col]]), ]
    if (nrow(sub) == 0) next
    gm <- compute_group_mean_sem(sub, score_col)
    gm_list[[seq_]]        <- gm
    fit_lin_list[[seq_]]   <- fit_kf_linear_by_block(gm$mean, gm$block)
    fit_sqrt_list[[seq_]]  <- fit_kf_sqrt_by_block(gm$mean, gm$block)
    boundaries_list[[seq_]] <- get_block_boundaries(sub)
    npart_list[[seq_]]     <- length(unique(sub$participant))
  }
  if (length(gm_list) < 2) return(NULL)

  all_y <- unlist(lapply(sequences, function(seq_) {
    fl <- fit_lin_list[[seq_]]; fs <- fit_sqrt_list[[seq_]]
    vals <- c()
    if (!is.null(fl)) vals <- c(vals, fl$x_smooth - ci_z * fl$se, fl$x_smooth + ci_z * fl$se)
    if (!is.null(fs)) vals <- c(vals, fs$x_smooth - ci_z * fs$se, fs$x_smooth + ci_z * fs$se)
    vals
  }))
  yr  <- range(all_y, na.rm = TRUE)
  pad <- diff(yr) * 0.05; if (!is.finite(pad) || pad == 0) pad <- 1
  yr  <- yr + c(-pad, pad)
  xr  <- range(unlist(lapply(gm_list, function(gm) gm$trial_global)), na.rm = TRUE)

  n_summary <- paste(
    sprintf("%s (n=%d)", sequences, vapply(sequences, function(s) npart_list[[s]], numeric(1))),
    collapse = " vs "
  )

  fname <- file.path(out_dir, sprintf("evolution_comparison_%s_AvsB.png", score_col))
  grDevices::png(fname, width = 1400, height = 700, res = 130)
  par(mar = c(3, 3, 3, 3), bg = "#F8F9FA")

  plot(NA, xlim = xr, ylim = yr, xlab = "Trial", ylab = meta$ylabel,
       main = sprintf("%s \u2014 %s (only models)", meta$label, n_summary),
       cex.main = 0.95, font.main = 2)
  grid(col = "gray90", lty = 1)
  for (seq_ in sequences) for (bd in boundaries_list[[seq_]]) abline(v = bd, col = "dimgray", lty = 2, lwd = 0.8)

  leg_labels <- c(); leg_colors <- c(); leg_lty <- c(); leg_lwd <- c()

  for (seq_ in sequences) {
    gm      <- gm_list[[seq_]]
    col_seq <- SEQ_COLORS[[seq_]]

    # --- Modèle linéaire x[n+1]=Ax[n]+B, avec IC ---
    fl <- fit_lin_list[[seq_]]
    if (!is.null(fl) && any(!is.na(fl$x_smooth))) {
      okl <- is.finite(fl$x_smooth) & is.finite(fl$se)
      if (any(okl)) {
        xx <- gm$trial_global[okl]
        lo <- fl$x_smooth[okl] - ci_z * fl$se[okl]
        hi <- fl$x_smooth[okl] + ci_z * fl$se[okl]
        polygon(c(xx, rev(xx)), c(lo, rev(hi)),
                col = grDevices::adjustcolor(col_seq, alpha.f = 0.10), border = NA)
      }
      lines(gm$trial_global, fl$x_smooth, col = col_seq, lwd = 2, lty = 2)
    }

    # --- Modèle racine carrée y[n]~=A*sqrt(n)+B, avec IC ---
    fs <- fit_sqrt_list[[seq_]]
    if (!is.null(fs) && any(!is.na(fs$x_smooth))) {
      oks <- is.finite(fs$x_smooth) & is.finite(fs$se)
      if (any(oks)) {
        xx <- gm$trial_global[oks]
        lo <- fs$x_smooth[oks] - ci_z * fs$se[oks]
        hi <- fs$x_smooth[oks] + ci_z * fs$se[oks]
        polygon(c(xx, rev(xx)), c(lo, rev(hi)),
                col = grDevices::adjustcolor(col_seq, alpha.f = 0.10), border = NA)
      }
      lines(gm$trial_global, fs$x_smooth, col = col_seq, lwd = 2, lty = 3)
    }

    leg_labels <- c(leg_labels,
                    sprintf("State space model (%s, n=%d)", seq_, npart_list[[seq_]]),
                    sprintf("Tendency \u221atrial (%s, n=%d)", seq_, npart_list[[seq_]]))
    leg_colors <- c(leg_colors, col_seq, col_seq)
    leg_lty    <- c(leg_lty, 2, 3)
    leg_lwd    <- c(leg_lwd, 2, 2)
  }

  legend(
    "bottomleft",
    legend = leg_labels,
    col = leg_colors,
    lty = leg_lty,
    lwd = leg_lwd,
    pch = 16,
    cex = 0.8,
    bty = "o",
    bg = "#F8F9FA",
    inset = 0.02
  )

  grDevices::dev.off()
  cat(sprintf("    Sauvegard\u00e9 : %s\n", fname))
  fname
}

# =======================================================================
# FIGURE 2 — Violons Séquence A vs B, Bloc 3 : lr_sqrt_slope + lr_ssm_A
#            (reprise exacte de plot_violin_learning_rates())
# =======================================================================

plot_violin_learning_rates <- function(lr_df, score_col, block_num, out_dir) {
  sub <- lr_df[lr_df$score == score_col & lr_df$block == block_num, ]
  if (nrow(sub) == 0) return(NULL)

  metrics <- c("lr_sqrt_slope", "lr_ssm_A")
  metric_labels <- c(
    lr_sqrt_slope = "Maximum slope (\u221a-model)",
    lr_ssm_A      = "A parameter (retention)"
  )

  meta  <- SCORE_META[[score_col]]
  label <- if (!is.null(meta)) meta$label else score_col

  nA <- length(unique(sub$participant[sub$sequence == "A"]))
  nB <- length(unique(sub$participant[sub$sequence == "B"]))

  fname <- file.path(out_dir, sprintf("violin_learning_rates_%s_block%s.png", score_col, block_num))
  grDevices::png(fname, width = 1500, height = 550, res = 130)
  par(mfrow = c(1, length(metrics)), mar = c(5, 4, 3, 1), oma = c(0, 0, 3, 0))

  for (m in metrics) {
    a_vals <- sub[sub$sequence == "A", m]
    b_vals <- sub[sub$sequence == "B", m]
    all_vals <- c(a_vals, b_vals)
    all_vals <- all_vals[is.finite(all_vals)]
    if (length(all_vals) == 0) { plot.new(); next }

    yr  <- range(all_vals, na.rm = TRUE)
    pad <- diff(yr) * 0.1
    if (!is.finite(pad) || pad == 0) pad <- 1
    yr <- yr + c(-pad, pad)

    plot(NA, xlim = c(0.5, 2.5), ylim = yr, xaxt = "n", xlab = "",
         ylab = metric_labels[[m]], main = metric_labels[[m]], cex.main = 0.85)
    axis(1, at = c(1, 2), labels = c("S\u00e9quence A", "S\u00e9quence B"),
         las = 2, cex.axis = 0.75)
    grid(nx = NA, ny = NULL, col = "gray90")
    draw_violin(a_vals, at = 1, col = SEQ_COLORS[["A"]])
    draw_violin(b_vals, at = 2, col = SEQ_COLORS[["B"]])
  }

  mtext(sprintf("Learning rates \u2014 %s (Bloc %s) \u2014 S\u00e9q. A (n=%d) vs S\u00e9q. B (n=%d)",
                label, block_num, nA, nB),
        outer = TRUE, cex = 1.1, font = 2)
  grDevices::dev.off()
  cat(sprintf("    Sauvegard\u00e9 : %s\n", fname))
  fname
}

# =======================================================================
# FIGURE 3 [NOUVEAU] — Violons Bloc 4 normalisé par Bloc 1, Séq. A vs B
# =======================================================================
# Suit la logique APPARIÉE PAR PARTICIPANT de
# spat_learning_rate_comparison_paired.R : seuls les participants
# présents dans les DEUX blocs (1 et 4) sont conservés (appariement par
# nom de participant), pour chaque séquence séparément. Le ratio
# Bloc4/Bloc1 est mis à NA si l'une des deux valeurs manque, ou si la
# valeur au Bloc 1 (dénominateur) est trop proche de 0 (ratio non
# interprétable) — même logique que compute_learning_rate_ratios()
# (ratio A/B) du script complet, appliquée ici aux blocs plutôt qu'aux
# séquences.

compute_block_ratio <- function(lr_df, score_col, block_num = 4, block_denom = 1) {
  metrics <- c("lr_sqrt_slope", "lr_ssm_A")
  sub_score <- lr_df[lr_df$score == score_col, ]
  if (nrow(sub_score) == 0) return(data.frame())

  rows <- list()
  for (seq_ in sort(unique(sub_score$sequence))) {
    sub_seq      <- sub_score[sub_score$sequence == seq_, ]
    participants <- sort(unique(sub_seq$participant))
    if (length(participants) == 0) next

    row_df <- data.frame(score = score_col, sequence = seq_, participant = participants,
                         stringsAsFactors = FALSE)

    for (m in metrics) {
      num_val <- vapply(participants, function(p) {
        v <- sub_seq[sub_seq$participant == p & sub_seq$block == block_num, m]
        if (length(v) == 0) NA_real_ else v[1]
      }, numeric(1))
      denom_val <- vapply(participants, function(p) {
        v <- sub_seq[sub_seq$participant == p & sub_seq$block == block_denom, m]
        if (length(v) == 0) NA_real_ else v[1]
      }, numeric(1))

      ratio <- num_val / denom_val
      unstable <- !is.finite(denom_val) | abs(denom_val) < 1e-8
      ratio[unstable] <- NA_real_

      row_df[[paste0("ratio_", m)]] <- ratio
    }

    rows[[length(rows) + 1]] <- row_df
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

plot_violin_block_ratio <- function(ratio_df, score_col, block_num = 4, block_denom = 1, out_dir) {
  sub <- ratio_df[ratio_df$score == score_col, ]
  if (nrow(sub) == 0) return(NULL)

  metrics <- c("ratio_lr_sqrt_slope", "ratio_lr_ssm_A")
  metric_labels <- c(
    ratio_lr_sqrt_slope = sprintf("Max slope ratio (\u221a-model)\nBloc%d / Bloc%d", block_num, block_denom),
    ratio_lr_ssm_A      = sprintf("A parameter ratio\nBloc%d / Bloc%d", block_num, block_denom)
  )

  meta  <- SCORE_META[[score_col]]
  label <- if (!is.null(meta)) meta$label else score_col

  nA <- length(unique(sub$participant[sub$sequence == "A" & is.finite(sub$ratio_lr_sqrt_slope)]))
  nB <- length(unique(sub$participant[sub$sequence == "B" & is.finite(sub$ratio_lr_sqrt_slope)]))

  fname <- file.path(out_dir, sprintf("violin_block%d_over_block%d_%s.png", block_num, block_denom, score_col))
  grDevices::png(fname, width = 1500, height = 550, res = 130)
  par(mfrow = c(1, length(metrics)), mar = c(5, 4, 3, 1), oma = c(0, 0, 3, 0))

  for (m in metrics) {
    a_vals <- sub[sub$sequence == "A", m]
    b_vals <- sub[sub$sequence == "B", m]
    all_vals <- c(a_vals, b_vals)
    all_vals <- all_vals[is.finite(all_vals)]
    if (length(all_vals) == 0) { plot.new(); next }

    yr  <- range(all_vals, na.rm = TRUE)
    yr  <- range(c(yr, 1))                     # inclut la ligne de référence ratio=1
    pad <- diff(yr) * 0.15
    if (!is.finite(pad) || pad == 0) pad <- 1
    yr <- yr + c(-pad, pad)

    plot(NA, xlim = c(0.5, 2.5), ylim = yr, xaxt = "n", xlab = "",
         ylab = metric_labels[[m]], main = metric_labels[[m]], cex.main = 0.8)
    axis(1, at = c(1, 2), labels = c("S\u00e9quence A", "S\u00e9quence B"),
         las = 2, cex.axis = 0.75)
    grid(nx = NA, ny = NULL, col = "gray90")
    abline(h = 1, lty = 2, lwd = 1.2, col = "gray40")  # ratio=1 -> pas de changement Bloc1->Bloc4
    draw_violin(a_vals, at = 1, col = SEQ_COLORS[["A"]])
    draw_violin(b_vals, at = 2, col = SEQ_COLORS[["B"]])
  }

  mtext(sprintf("Bloc %d normalis\u00e9 par Bloc %d (ratio, appari\u00e9 par participant) \u2014 %s\nS\u00e9q. A (n=%d) vs S\u00e9q. B (n=%d)",
                block_num, block_denom, label, nA, nB),
        outer = TRUE, cex = 1.0, font = 2)
  grDevices::dev.off()
  cat(sprintf("    Sauvegard\u00e9 : %s\n", fname))
  fname
}

# ---------------------------------------------------------------------
# PIPELINE PRINCIPAL
# ---------------------------------------------------------------------

main <- function(scores_folder, out_dir = "./output_evolution_and_ratio") {
  cat(strrep("=", 65), "\n")
  cat("SPAT \u2014 EVOLUTION (2 MODELES) + VIOLONS BLOC3 + RATIO BLOC4/BLOC1\n")
  cat(strrep("=", 65), "\n")

  cat(sprintf("\nChargement des scores depuis : %s\n", scores_folder))
  all_df <- load_participant_scores(scores_folder)
  if (nrow(all_df) == 0) return(invisible(NULL))

  score_cols <- intersect(names(SCORE_META), colnames(all_df))
  if (length(score_cols) == 0) {
    cat("ERREUR : la colonne score_harmonic_norm_participant est absente des donn\u00e9es.\n")
    return(invisible(NULL))
  }
  score_col <- score_cols[1]
  cat(sprintf("  Score utilis\u00e9 : %s\n", score_col))

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  cat("\n[1/3] Calcul des learning rates (tous blocs)...\n")
  lr_df <- compute_learning_rates_all(all_df, score_cols)
  if (nrow(lr_df) > 0) {
    lr_path <- file.path(out_dir, "learning_rates_all_blocks.csv")
    utils::write.csv(lr_df, lr_path, row.names = FALSE)
    cat(sprintf("  Learning rates \u00e9crits dans : %s\n", lr_path))
  }

  cat("\n[2/3] Figures...\n")
  cat("  Figure 1 : \u00e9volution (2 mod\u00e8les superpos\u00e9s), S\u00e9q. A vs B...\n")
  plot_evolution_group_comparison(all_df, score_col, out_dir)

  cat("  Figure 2 : violons learning rates, Bloc 3, S\u00e9q. A vs B...\n")
  plot_violin_learning_rates(lr_df, score_col, block_num = 3, out_dir)
  #!/usr/bin/env Rscript
  # =====================================================================
  # SPAT — Évolution (2 modèles superposés) + violons Bloc 3 (Séq. A vs B)
  #         + violons Bloc 4 normalisé par Bloc 1 (appariés par participant)
  # =====================================================================
  # Score utilisé : score_harmonic_norm_participant (uniquement).
  #
  # Ce script extrait/assemble 3 figures :
  #
  #   1) evolution_comparison_<score>_AvsB.png
  #        Reprise EXACTE de plot_evolution_group_comparison() du script
  #        complet : les 2 modèles (état-espace linéaire ET tendance
  #        √trial, ajustés par filtre de Kalman) superposés sur le MÊME
  #        graphique, avec bandes de confiance, pour Séquence A et B.
  #
  #   2) violin_learning_rates_<score>_block3.png
  #        Reprise de plot_violin_learning_rates() du script complet,
  #        restreinte au Bloc 3 : violons "Max slope (√-model)",
  #        "A parameter (retention)" ET "B parameter (intercept)",
  #        Séquence A vs Séquence B.
  #
  #   3) violin_block4_over_block1_<score>.png   [NOUVEAU]
  #        En suivant la logique APPARIÉE PAR PARTICIPANT du script de
  #        comparaison Bloc 4 (spat_learning_rate_comparison_paired.R) :
  #        pour chaque participant x séquence, on calcule le ratio
  #        Bloc4/Bloc1 de lr_sqrt_slope, lr_ssm_A ET lr_ssm_B
  #        (participant non présent dans les deux blocs -> NA ; ratio
  #        instable si la valeur du Bloc 1 est proche de 0 -> NA, même
  #        logique que compute_learning_rate_ratios() pour le ratio A/B).
  #        Le résultat est tracé en violon, Séquence A vs B, avec une
  #        ligne de référence à 1 (= pas de changement Bloc1->Bloc4).
  #
  # SORTIES (dans out_dir) :
  #   evolution_comparison_score_harmonic_norm_participant_AvsB.png
  #   violin_learning_rates_score_harmonic_norm_participant_block3.png
  #   violin_block4_over_block1_score_harmonic_norm_participant.png
  #   learning_rates_all_blocks.csv
  #   block4_over_block1_ratios.csv
  #
  # USAGE
  # -----
  #   Rscript spat_evolution_and_block_ratio_violins.R scores_folder out_dir
  #
  # Packages requis : readxl
  # =====================================================================
  
  suppressPackageStartupMessages({
    library(readxl)
    library(tools)
  })
  
  # ---------------------------------------------------------------------
  # CONFIGURATION
  # ---------------------------------------------------------------------
  
  SEQ_COLORS <- c(A = "#4C8CC5", B = "#F1BC18")
  
  SCORE_META <- list(
    score_harmonic_norm_participant = list(
      label  = "Harmonic score (normalized by participant)",
      ylabel = "Harmonic score (normalized by participant)"
    )
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
    
    all_df <- all_df[order(all_df$participant, all_df$sequence, all_df$block, all_df$trial), ]
    grp    <- interaction(all_df$participant, all_df$sequence, drop = TRUE)
    all_df$trial_global <- ave(seq_len(nrow(all_df)), grp, FUN = seq_along)
    
    all_df
  }
  
  get_block_boundaries <- function(sub) {
    blocks <- sort(unique(sub$block))
    blocks <- blocks[!is.na(blocks)]
    if (length(blocks) <= 1) return(numeric(0))
    
    boundaries <- c()
    for (b in blocks[-length(blocks)]) {
      bd  <- sub[sub$block == b, ]
      mx  <- tapply(bd$trial_global, bd$participant, max)
      boundaries <- c(boundaries, stats::median(mx, na.rm = TRUE))
    }
    boundaries
  }
  
  # ---------------------------------------------------------------------
  # 2. FILTRE DE KALMAN GENERIQUE + LISSEUR RTS
  # ---------------------------------------------------------------------
  
  kf_run <- function(y, Phi, b, Q, H_fun, R, x0, P0) {
    n <- length(y)
    d <- length(x0)
    
    x_pred <- matrix(NA_real_, n, d)
    x_filt <- matrix(NA_real_, n, d)
    P_pred <- vector("list", n)
    P_filt <- vector("list", n)
    loglik <- 0
    
    x_prev <- x0
    P_prev <- P0
    
    for (t in seq_len(n)) {
      xp <- as.numeric(Phi %*% x_prev + b)
      Pp <- Phi %*% P_prev %*% t(Phi) + Q
      x_pred[t, ] <- xp
      P_pred[[t]] <- Pp
      
      if (is.na(y[t])) {
        xf <- xp
        Pf <- Pp
      } else {
        Ht   <- H_fun(t)
        yhat <- as.numeric(Ht %*% xp)
        S    <- as.numeric(Ht %*% Pp %*% t(Ht)) + R
        if (!is.finite(S) || S <= 0) S <- 1e-6
        
        K <- (Pp %*% t(Ht)) / S
        e <- y[t] - yhat
        
        xf <- xp + as.numeric(K * e)
        Pf <- (diag(d) - K %*% Ht) %*% Pp
        
        loglik <- loglik - 0.5 * (log(2 * pi * S) + e^2 / S)
      }
      x_filt[t, ] <- xf
      P_filt[[t]] <- Pf
      
      x_prev <- xf
      P_prev <- Pf
    }
    
    list(x_pred = x_pred, P_pred = P_pred, x_filt = x_filt, P_filt = P_filt, loglik = loglik)
  }
  
  kf_smooth <- function(kf_out, Phi) {
    n <- nrow(kf_out$x_filt)
    d <- ncol(kf_out$x_filt)
    
    x_smooth <- matrix(NA_real_, n, d)
    P_smooth <- vector("list", n)
    
    x_smooth[n, ] <- kf_out$x_filt[n, ]
    P_smooth[[n]] <- kf_out$P_filt[[n]]
    
    if (n >= 2) {
      for (t in (n - 1):1) {
        Pf      <- kf_out$P_filt[[t]]
        Pp_next <- kf_out$P_pred[[t + 1]]
        J       <- Pf %*% t(Phi) %*% solve(Pp_next)
        
        x_smooth[t, ] <- kf_out$x_filt[t, ] +
          as.numeric(J %*% (x_smooth[t + 1, ] - kf_out$x_pred[t + 1, ]))
        P_smooth[[t]] <- Pf + J %*% (P_smooth[[t + 1]] - Pp_next) %*% t(J)
      }
    }
    
    list(x_smooth = x_smooth, P_smooth = P_smooth)
  }
  
  # ---------------------------------------------------------------------
  # 3. MODELE (1) — ETAT-ESPACE LINEAIRE  x[n+1] = A*x[n] + B
  # ---------------------------------------------------------------------
  
  fit_kf_linear <- function(y) {
    y <- as.numeric(y)
    n_valid <- sum(!is.na(y))
    if (n_valid < 5) return(NULL)
    
    y_obs <- y[!is.na(y)]
    var_y <- stats::var(y_obs)
    if (!is.finite(var_y) || var_y <= 0) var_y <- 1
    
    init_A     <- 0.8
    init_B     <- mean(diff(y_obs), na.rm = TRUE) * (1 - init_A)
    if (!is.finite(init_B)) init_B <- 0
    init_logQ  <- log(max(var_y * 0.1, 1e-4))
    init_logR  <- log(max(var_y * 0.5, 1e-4))
    
    x0 <- c(y_obs[1])
    P0 <- matrix(var_y * 10, 1, 1)
    H_fun <- function(t) matrix(1, 1, 1)
    
    negloglik <- function(par) {
      A <- par[1]; B <- par[2]; Q <- exp(par[3]); R <- exp(par[4])
      Phi <- matrix(A, 1, 1)
      out <- tryCatch(kf_run(y, Phi, B, matrix(Q, 1, 1), H_fun, R, x0, P0),
                      error = function(e) NULL)
      if (is.null(out) || !is.finite(out$loglik)) return(1e10)
      -out$loglik
    }
    
    lower <- c(-0.999, -Inf, log(var_y * 1e-6), log(var_y * 0.05))
    upper <- c( 0.999,  Inf, log(var_y * 0.05), log(var_y * 10))
    
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
    
    Phi    <- matrix(A, 1, 1)
    result <- tryCatch({
      kf_out <- kf_run(y, Phi, B, matrix(Q, 1, 1), H_fun, R, x0, P0)
      sm     <- kf_smooth(kf_out, Phi)
      se_vec <- sqrt(vapply(sm$P_smooth, function(P) P[1, 1], numeric(1)))
      list(A = A, B = B, Q = Q, R = R,
           x_smooth = as.numeric(sm$x_smooth[, 1]),
           se       = se_vec)
    }, error = function(e) NULL)
    
    result
  }
  
  # ---------------------------------------------------------------------
  # 4. MODELE (2) — TENDANCE "RACINE CARREE"  y[n] ~= A_n*sqrt(n) + B_n
  # ---------------------------------------------------------------------
  
  fit_kf_sqrt <- function(y, trial_n) {
    y <- as.numeric(y)
    trial_n <- as.numeric(trial_n)
    ok <- !is.na(y)
    n_valid <- sum(ok)
    if (n_valid < 5) return(NULL)
    
    lm_fit <- tryCatch(stats::lm(y[ok] ~ sqrt(trial_n[ok])), error = function(e) NULL)
    if (is.null(lm_fit)) {
      B0 <- mean(y[ok]); A0 <- 0
    } else {
      B0 <- unname(stats::coef(lm_fit)[1])
      A0 <- unname(stats::coef(lm_fit)[2])
      if (!is.finite(A0)) A0 <- 0
      if (!is.finite(B0)) B0 <- mean(y[ok])
    }
    
    var_y <- stats::var(y[ok])
    if (!is.finite(var_y) || var_y <= 0) var_y <- 1
    
    init_logQa <- log(max(var_y * 0.001, 1e-6))
    init_logQb <- log(max(var_y * 0.001, 1e-6))
    init_logR  <- log(max(var_y * 0.5, 1e-4))
    
    Phi   <- diag(2)
    b     <- c(0, 0)
    x0    <- c(A0, B0)
    P0    <- diag(c(var_y * 10, var_y * 10))
    H_fun <- function(t) matrix(c(sqrt(trial_n[t]), 1), nrow = 1)
    
    negloglik <- function(par) {
      Qa <- exp(par[1]); Qb <- exp(par[2]); R <- exp(par[3])
      Qm <- diag(c(Qa, Qb))
      out <- tryCatch(kf_run(y, Phi, b, Qm, H_fun, R, x0, P0), error = function(e) NULL)
      if (is.null(out) || !is.finite(out$loglik)) return(1e10)
      -out$loglik
    }
    
    lower <- c(log(var_y * 1e-10), log(var_y * 1e-10), log(var_y * 0.05))
    upper <- c(log(var_y * 1e-4),  log(var_y * 1e-4),  log(var_y * 10))
    fit <- tryCatch(
      stats::optim(c(init_logQa, init_logQb, init_logR), negloglik,
                   method = "L-BFGS-B", lower = lower, upper = upper,
                   control = list(maxit = 500)),
      error = function(e) NULL
    )
    
    if (is.null(fit) || !is.finite(fit$value)) {
      Qa <- exp(init_logQa); Qb <- exp(init_logQb); R <- exp(init_logR)
    } else {
      Qa <- exp(fit$par[1]); Qb <- exp(fit$par[2]); R <- exp(fit$par[3])
    }
    
    Qm     <- diag(c(Qa, Qb))
    result <- tryCatch({
      kf_out <- kf_run(y, Phi, b, Qm, H_fun, R, x0, P0)
      sm     <- kf_smooth(kf_out, Phi)
      n <- length(y)
      fitted_mean <- vapply(seq_len(n), function(t) as.numeric(H_fun(t) %*% sm$x_smooth[t, ]), numeric(1))
      se_vec <- vapply(seq_len(n), function(t) {
        Ht <- H_fun(t)
        sqrt(as.numeric(Ht %*% sm$P_smooth[[t]] %*% t(Ht)))
      }, numeric(1))
      list(Qa = Qa, Qb = Qb, R = R,
           x_smooth = fitted_mean,
           se       = se_vec)
    }, error = function(e) NULL)
    
    result
  }
  
  # ---------------------------------------------------------------------
  # Variantes "par bloc" : réinitialisent l'état au début de CHAQUE bloc.
  # ---------------------------------------------------------------------
  
  fit_kf_linear_by_block <- function(y, block) {
    n <- length(y)
    x_smooth <- rep(NA_real_, n)
    se       <- rep(NA_real_, n)
    for (b in sort(unique(block))) {
      idx <- which(block == b)
      fl  <- fit_kf_linear(y[idx])
      if (is.null(fl)) next
      x_smooth[idx] <- fl$x_smooth
      se[idx]       <- fl$se
    }
    list(x_smooth = x_smooth, se = se)
  }
  
  fit_kf_sqrt_by_block <- function(y, block) {
    n <- length(y)
    x_smooth <- rep(NA_real_, n)
    se       <- rep(NA_real_, n)
    for (b in sort(unique(block))) {
      idx     <- which(block == b)
      n_local <- seq_along(idx)
      fs <- fit_kf_sqrt(y[idx], n_local)
      if (is.null(fs)) next
      x_smooth[idx] <- fs$x_smooth
      se[idx]       <- fs$se
    }
    list(x_smooth = x_smooth, se = se)
  }
  
  # ---------------------------------------------------------------------
  # 5. VIOLON "MAISON" (densité miroir + médiane/IQR)
  # ---------------------------------------------------------------------
  
  draw_violin <- function(x, at, col, width = 0.35, border = "black") {
    x <- x[is.finite(x)]
    if (length(x) < 2) return(invisible(NULL))
    if (stats::sd(x) < 1e-8) {
      segments(at - width, x[1], at + width, x[1], col = border, lwd = 2)
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
  
  # ---------------------------------------------------------------------
  # 6. LEARNING RATES PAR BLOC (score x participant x sequence x block)
  # ---------------------------------------------------------------------
  
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
  
  simulate_ssm_r <- function(y0, A, B, n) {
    out <- numeric(n)
    out[1] <- y0
    if (n > 1) for (i in 2:n) out[i] <- A * out[i - 1] + B
    out
  }
  
  # Une ligne par (score, participant, sequence, block) : les 3 learning
  # rates (dont lr_sqrt_slope et lr_ssm_A) pour CHAQUE bloc présent dans
  # les données (pas seulement le Bloc 3), afin de pouvoir ensuite
  # calculer le ratio Bloc4/Bloc1.
  compute_learning_rates_all <- function(df, score_cols) {
    rows <- list()
    for (score_col in score_cols) {
      sub_all <- df[!is.na(df[[score_col]]) & !is.na(df$sequence) & !is.na(df$block), ]
      for (part in sort(unique(sub_all$participant))) {
        for (seq_ in sort(unique(sub_all$sequence))) {
          pdata <- sub_all[sub_all$participant == part & sub_all$sequence == seq_, ]
          for (b in sort(unique(pdata$block))) {
            bd <- pdata[pdata$block == b, ]
            bd <- bd[order(bd$trial), ]
            if (nrow(bd) < 5) next
            
            y  <- bd[[score_col]]
            fl <- fit_kf_linear(y)
            if (is.null(fl)) next
            y_ssm <- simulate_ssm_r(y[1], fl$A, fl$B, length(y))
            
            fs <- fit_kf_sqrt(y, seq_along(y))
            lr_sqrt_slope <- if (!is.null(fs)) learning_rate_max_slope(fs$x_smooth) else NA_real_
            
            rows[[length(rows) + 1]] <- data.frame(
              score = score_col, participant = part, sequence = seq_, block = b,
              lr_ssm_1090   = learning_rate_10_90(y_ssm, bd$trial),
              lr_sqrt_slope = lr_sqrt_slope,
              lr_ssm_slope  = learning_rate_max_slope(y_ssm),
              lr_ssm_A      = fl$A,
              lr_ssm_B      = fl$B,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    if (length(rows) == 0) return(data.frame())
    do.call(rbind, rows)
  }
  
  # ---------------------------------------------------------------------
  # 7. MOYENNE DE GROUPE + SEM PAR POSITION CHRONOLOGIQUE
  # ---------------------------------------------------------------------
  
  compute_group_mean_sem <- function(sub, score_col) {
    agg_mean <- stats::aggregate(sub[[score_col]], by = list(trial_global = sub$trial_global),
                                 FUN = mean, na.rm = TRUE)
    names(agg_mean)[2] <- "mean"
    agg_sd <- stats::aggregate(sub[[score_col]], by = list(trial_global = sub$trial_global),
                               FUN = stats::sd, na.rm = TRUE)
    names(agg_sd)[2] <- "sd"
    agg_n <- stats::aggregate(sub[[score_col]], by = list(trial_global = sub$trial_global),
                              FUN = function(v) sum(!is.na(v)))
    names(agg_n)[2] <- "n"
    gm_block <- stats::aggregate(block ~ trial_global, data = sub,
                                 FUN = function(s) as.numeric(names(sort(table(s), decreasing = TRUE))[1]))
    
    out <- merge(agg_mean, agg_sd, by = "trial_global")
    out <- merge(out, agg_n,  by = "trial_global")
    out <- merge(out, gm_block, by = "trial_global")
    out$sem <- out$sd / sqrt(pmax(out$n, 1))
    out[order(out$trial_global), ]
  }
  
  # =======================================================================
  # FIGURE 1 — Comparaison A vs B : moyenne de groupe + LES 2 MODELES
  #            (état-espace linéaire + tendance √) SUPERPOSES + IC
  #            (reprise exacte de plot_evolution_group_comparison())
  # =======================================================================
  
  plot_evolution_group_comparison <- function(df, score_col, out_dir, ci_z = 1.96) {
    meta <- SCORE_META[[score_col]]
    if (is.null(meta)) meta <- list(label = score_col, ylabel = score_col)
    
    sequences <- intersect(sort(unique(stats::na.omit(df$sequence))), names(SEQ_COLORS))
    if (length(sequences) < 2) return(NULL)
    
    gm_list <- list(); fit_lin_list <- list(); fit_sqrt_list <- list(); boundaries_list <- list()
    npart_list <- list()
    
    for (seq_ in sequences) {
      sub <- df[!is.na(df$sequence) & df$sequence == seq_ & !is.na(df[[score_col]]), ]
      if (nrow(sub) == 0) next
      gm <- compute_group_mean_sem(sub, score_col)
      gm_list[[seq_]]        <- gm
      fit_lin_list[[seq_]]   <- fit_kf_linear_by_block(gm$mean, gm$block)
      fit_sqrt_list[[seq_]]  <- fit_kf_sqrt_by_block(gm$mean, gm$block)
      boundaries_list[[seq_]] <- get_block_boundaries(sub)
      npart_list[[seq_]]     <- length(unique(sub$participant))
    }
    if (length(gm_list) < 2) return(NULL)
    
    all_y <- unlist(lapply(sequences, function(seq_) {
      fl <- fit_lin_list[[seq_]]; fs <- fit_sqrt_list[[seq_]]
      vals <- c()
      if (!is.null(fl)) vals <- c(vals, fl$x_smooth - ci_z * fl$se, fl$x_smooth + ci_z * fl$se)
      if (!is.null(fs)) vals <- c(vals, fs$x_smooth - ci_z * fs$se, fs$x_smooth + ci_z * fs$se)
      vals
    }))
    yr  <- range(all_y, na.rm = TRUE)
    pad <- diff(yr) * 0.05; if (!is.finite(pad) || pad == 0) pad <- 1
    yr  <- yr + c(-pad, pad)
    xr  <- range(unlist(lapply(gm_list, function(gm) gm$trial_global)), na.rm = TRUE)
    
    n_summary <- paste(
      sprintf("%s (n=%d)", sequences, vapply(sequences, function(s) npart_list[[s]], numeric(1))),
      collapse = " vs "
    )
    
    fname <- file.path(out_dir, sprintf("evolution_comparison_%s_AvsB.png", score_col))
    grDevices::png(fname, width = 1400, height = 700, res = 130)
    par(mar = c(3, 3, 3, 3), bg = "#F8F9FA")
    
    plot(NA, xlim = xr, ylim = yr, xlab = "Trial", ylab = meta$ylabel,
         main = sprintf("%s \u2014 %s (only models)", meta$label, n_summary),
         cex.main = 0.95, font.main = 2)
    grid(col = "gray90", lty = 1)
    for (seq_ in sequences) for (bd in boundaries_list[[seq_]]) abline(v = bd, col = "dimgray", lty = 2, lwd = 0.8)
    
    leg_labels <- c(); leg_colors <- c(); leg_lty <- c(); leg_lwd <- c()
    
    for (seq_ in sequences) {
      gm      <- gm_list[[seq_]]
      col_seq <- SEQ_COLORS[[seq_]]
      
      # --- Modèle linéaire x[n+1]=Ax[n]+B, avec IC ---
      fl <- fit_lin_list[[seq_]]
      if (!is.null(fl) && any(!is.na(fl$x_smooth))) {
        okl <- is.finite(fl$x_smooth) & is.finite(fl$se)
        if (any(okl)) {
          xx <- gm$trial_global[okl]
          lo <- fl$x_smooth[okl] - ci_z * fl$se[okl]
          hi <- fl$x_smooth[okl] + ci_z * fl$se[okl]
          polygon(c(xx, rev(xx)), c(lo, rev(hi)),
                  col = grDevices::adjustcolor(col_seq, alpha.f = 0.10), border = NA)
        }
        lines(gm$trial_global, fl$x_smooth, col = col_seq, lwd = 2, lty = 2)
      }
      
      # --- Modèle racine carrée y[n]~=A*sqrt(n)+B, avec IC ---
      fs <- fit_sqrt_list[[seq_]]
      if (!is.null(fs) && any(!is.na(fs$x_smooth))) {
        oks <- is.finite(fs$x_smooth) & is.finite(fs$se)
        if (any(oks)) {
          xx <- gm$trial_global[oks]
          lo <- fs$x_smooth[oks] - ci_z * fs$se[oks]
          hi <- fs$x_smooth[oks] + ci_z * fs$se[oks]
          polygon(c(xx, rev(xx)), c(lo, rev(hi)),
                  col = grDevices::adjustcolor(col_seq, alpha.f = 0.10), border = NA)
        }
        lines(gm$trial_global, fs$x_smooth, col = col_seq, lwd = 2, lty = 3)
      }
      
      leg_labels <- c(leg_labels,
                      sprintf("State space model (%s, n=%d)", seq_, npart_list[[seq_]]),
                      sprintf("Tendency \u221atrial (%s, n=%d)", seq_, npart_list[[seq_]]))
      leg_colors <- c(leg_colors, col_seq, col_seq)
      leg_lty    <- c(leg_lty, 2, 3)
      leg_lwd    <- c(leg_lwd, 2, 2)
    }
    
    legend(
      "bottomleft",
      legend = leg_labels,
      col = leg_colors,
      lty = leg_lty,
      lwd = leg_lwd,
      pch = 16,
      cex = 0.8,
      bty = "o",
      bg = "#F8F9FA",
      inset = 0.02
    )
    
    grDevices::dev.off()
    cat(sprintf("    Sauvegard\u00e9 : %s\n", fname))
    fname
  }
  
  # =======================================================================
  # FIGURE 2 — Violons Séquence A vs B, Bloc 3 : lr_sqrt_slope + lr_ssm_A
  #            (reprise exacte de plot_violin_learning_rates())
  # =======================================================================
  
  plot_violin_learning_rates <- function(lr_df, score_col, block_num, out_dir) {
    sub <- lr_df[lr_df$score == score_col & lr_df$block == block_num, ]
    if (nrow(sub) == 0) return(NULL)
    
    metrics <- c("lr_sqrt_slope", "lr_ssm_A", "lr_ssm_B")
    metric_labels <- c(
      lr_sqrt_slope = "Maximum slope (\u221a-model)",
      lr_ssm_A      = "A parameter (retention)",
      lr_ssm_B      = "B parameter (intercept)"
    )
    
    meta  <- SCORE_META[[score_col]]
    label <- if (!is.null(meta)) meta$label else score_col
    
    nA <- length(unique(sub$participant[sub$sequence == "A"]))
    nB <- length(unique(sub$participant[sub$sequence == "B"]))
    
    fname <- file.path(out_dir, sprintf("violin_learning_rates_%s_block%s.png", score_col, block_num))
    grDevices::png(fname, width = 2100, height = 550, res = 130)
    par(mfrow = c(1, length(metrics)), mar = c(5, 4, 3, 1), oma = c(0, 0, 3, 0))
    
    for (m in metrics) {
      a_vals <- sub[sub$sequence == "A", m]
      b_vals <- sub[sub$sequence == "B", m]
      all_vals <- c(a_vals, b_vals)
      all_vals <- all_vals[is.finite(all_vals)]
      if (length(all_vals) == 0) { plot.new(); next }
      
      yr  <- range(all_vals, na.rm = TRUE)
      pad <- diff(yr) * 0.1
      if (!is.finite(pad) || pad == 0) pad <- 1
      yr <- yr + c(-pad, pad)
      
      plot(NA, xlim = c(0.5, 2.5), ylim = yr, xaxt = "n", xlab = "",
           ylab = metric_labels[[m]], main = metric_labels[[m]], cex.main = 0.85)
      axis(1, at = c(1, 2), labels = c("S\u00e9quence A", "S\u00e9quence B"),
           las = 2, cex.axis = 0.75)
      grid(nx = NA, ny = NULL, col = "gray90")
      draw_violin(a_vals, at = 1, col = SEQ_COLORS[["A"]])
      draw_violin(b_vals, at = 2, col = SEQ_COLORS[["B"]])
    }
    
    mtext(sprintf("Learning rates \u2014 %s (Bloc %s) \u2014 S\u00e9q. A (n=%d) vs S\u00e9q. B (n=%d)",
                  label, block_num, nA, nB),
          outer = TRUE, cex = 1.1, font = 2)
    grDevices::dev.off()
    cat(sprintf("    Sauvegard\u00e9 : %s\n", fname))
    fname
  }
  
  # =======================================================================
  # FIGURE 3 [NOUVEAU] — Violons Bloc 4 normalisé par Bloc 1, Séq. A vs B
  # =======================================================================
  # Suit la logique APPARIÉE PAR PARTICIPANT de
  # spat_learning_rate_comparison_paired.R : seuls les participants
  # présents dans les DEUX blocs (1 et 4) sont conservés (appariement par
  # nom de participant), pour chaque séquence séparément. Le ratio
  # Bloc4/Bloc1 est mis à NA si l'une des deux valeurs manque, ou si la
  # valeur au Bloc 1 (dénominateur) est trop proche de 0 (ratio non
  # interprétable) — même logique que compute_learning_rate_ratios()
  # (ratio A/B) du script complet, appliquée ici aux blocs plutôt qu'aux
  # séquences.
  
  compute_block_ratio <- function(lr_df, score_col, block_num = 4, block_denom = 1) {
    metrics <- c("lr_sqrt_slope", "lr_ssm_A", "lr_ssm_B")
    sub_score <- lr_df[lr_df$score == score_col, ]
    if (nrow(sub_score) == 0) return(data.frame())
    
    rows <- list()
    for (seq_ in sort(unique(sub_score$sequence))) {
      sub_seq      <- sub_score[sub_score$sequence == seq_, ]
      participants <- sort(unique(sub_seq$participant))
      if (length(participants) == 0) next
      
      row_df <- data.frame(score = score_col, sequence = seq_, participant = participants,
                           stringsAsFactors = FALSE)
      
      for (m in metrics) {
        num_val <- vapply(participants, function(p) {
          v <- sub_seq[sub_seq$participant == p & sub_seq$block == block_num, m]
          if (length(v) == 0) NA_real_ else v[1]
        }, numeric(1))
        denom_val <- vapply(participants, function(p) {
          v <- sub_seq[sub_seq$participant == p & sub_seq$block == block_denom, m]
          if (length(v) == 0) NA_real_ else v[1]
        }, numeric(1))
        
        ratio <- num_val / denom_val
        unstable <- !is.finite(denom_val) | abs(denom_val) < 1e-8
        ratio[unstable] <- NA_real_
        
        row_df[[paste0("ratio_", m)]] <- ratio
      }
      
      rows[[length(rows) + 1]] <- row_df
    }
    if (length(rows) == 0) return(data.frame())
    do.call(rbind, rows)
  }
  
  plot_violin_block_ratio <- function(ratio_df, score_col, block_num = 4, block_denom = 1, out_dir) {
    sub <- ratio_df[ratio_df$score == score_col, ]
    if (nrow(sub) == 0) return(NULL)
    
    metrics <- c("ratio_lr_sqrt_slope", "ratio_lr_ssm_A", "ratio_lr_ssm_B")
    metric_labels <- c(
      ratio_lr_sqrt_slope = sprintf("Max slope ratio (\u221a-model)\nBloc%d / Bloc%d", block_num, block_denom),
      ratio_lr_ssm_A      = sprintf("A parameter ratio\nBloc%d / Bloc%d", block_num, block_denom),
      ratio_lr_ssm_B      = sprintf("B parameter ratio\nBloc%d / Bloc%d", block_num, block_denom)
    )
    
    meta  <- SCORE_META[[score_col]]
    label <- if (!is.null(meta)) meta$label else score_col
    
    nA <- length(unique(sub$participant[sub$sequence == "A" & is.finite(sub$ratio_lr_sqrt_slope)]))
    nB <- length(unique(sub$participant[sub$sequence == "B" & is.finite(sub$ratio_lr_sqrt_slope)]))
    
    fname <- file.path(out_dir, sprintf("violin_block%d_over_block%d_%s.png", block_num, block_denom, score_col))
    grDevices::png(fname, width = 2100, height = 550, res = 130)
    par(mfrow = c(1, length(metrics)), mar = c(5, 4, 3, 1), oma = c(0, 0, 3, 0))
    
    for (m in metrics) {
      a_vals <- sub[sub$sequence == "A", m]
      b_vals <- sub[sub$sequence == "B", m]
      all_vals <- c(a_vals, b_vals)
      all_vals <- all_vals[is.finite(all_vals)]
      if (length(all_vals) == 0) { plot.new(); next }
      
      yr  <- range(all_vals, na.rm = TRUE)
      yr  <- range(c(yr, 1))                     # inclut la ligne de référence ratio=1
      pad <- diff(yr) * 0.15
      if (!is.finite(pad) || pad == 0) pad <- 1
      yr <- yr + c(-pad, pad)
      
      plot(NA, xlim = c(0.5, 2.5), ylim = yr, xaxt = "n", xlab = "",
           ylab = metric_labels[[m]], main = metric_labels[[m]], cex.main = 0.8)
      axis(1, at = c(1, 2), labels = c("S\u00e9quence A", "S\u00e9quence B"),
           las = 2, cex.axis = 0.75)
      grid(nx = NA, ny = NULL, col = "gray90")
      abline(h = 1, lty = 2, lwd = 1.2, col = "gray40")  # ratio=1 -> pas de changement Bloc1->Bloc4
      draw_violin(a_vals, at = 1, col = SEQ_COLORS[["A"]])
      draw_violin(b_vals, at = 2, col = SEQ_COLORS[["B"]])
    }
    
    mtext(sprintf("Bloc %d normalis\u00e9 par Bloc %d (ratio, appari\u00e9 par participant) \u2014 %s\nS\u00e9q. A (n=%d) vs S\u00e9q. B (n=%d)",
                  block_num, block_denom, label, nA, nB),
          outer = TRUE, cex = 1.0, font = 2)
    grDevices::dev.off()
    cat(sprintf("    Sauvegard\u00e9 : %s\n", fname))
    fname
  }
  
  # ---------------------------------------------------------------------
  # PIPELINE PRINCIPAL
  # ---------------------------------------------------------------------
  
  main <- function(scores_folder, out_dir = "./output_evolution_and_ratio") {
    cat(strrep("=", 65), "\n")
    cat("SPAT \u2014 EVOLUTION (2 MODELES) + VIOLONS BLOC3 + RATIO BLOC4/BLOC1\n")
    cat(strrep("=", 65), "\n")
    
    cat(sprintf("\nChargement des scores depuis : %s\n", scores_folder))
    all_df <- load_participant_scores(scores_folder)
    if (nrow(all_df) == 0) return(invisible(NULL))
    
    score_cols <- intersect(names(SCORE_META), colnames(all_df))
    if (length(score_cols) == 0) {
      cat("ERREUR : la colonne score_harmonic_norm_participant est absente des donn\u00e9es.\n")
      return(invisible(NULL))
    }
    score_col <- score_cols[1]
    cat(sprintf("  Score utilis\u00e9 : %s\n", score_col))
    
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    cat("\n[1/3] Calcul des learning rates (tous blocs)...\n")
    lr_df <- compute_learning_rates_all(all_df, score_cols)
    if (nrow(lr_df) > 0) {
      lr_path <- file.path(out_dir, "learning_rates_all_blocks.csv")
      utils::write.csv(lr_df, lr_path, row.names = FALSE)
      cat(sprintf("  Learning rates \u00e9crits dans : %s\n", lr_path))
    }
    
    cat("\n[2/3] Figures...\n")
    cat("  Figure 1 : \u00e9volution (2 mod\u00e8les superpos\u00e9s), S\u00e9q. A vs B...\n")
    plot_evolution_group_comparison(all_df, score_col, out_dir)
    
    cat("  Figure 2 : violons learning rates, Bloc 3, S\u00e9q. A vs B...\n")
    plot_violin_learning_rates(lr_df, score_col, block_num = 3, out_dir)
    
    cat("\n[3/3] Ratio Bloc4/Bloc1 (appari\u00e9 par participant)...\n")
    ratio_df <- compute_block_ratio(lr_df, score_col, block_num = 4, block_denom = 1)
    if (nrow(ratio_df) > 0) {
      ratio_path <- file.path(out_dir, "block4_over_block1_ratios.csv")
      utils::write.csv(ratio_df, ratio_path, row.names = FALSE)
      cat(sprintf("  Ratios \u00e9crits dans : %s\n", ratio_path))
      cat("  Figure 3 : violons ratio Bloc4/Bloc1, S\u00e9q. A vs B...\n")
      plot_violin_block_ratio(ratio_df, score_col, block_num = 4, block_denom = 1, out_dir)
    } else {
      cat("  Aucun ratio calculable (Bloc 1 et/ou Bloc 4 absents ou insuffisants).\n")
    }
    
    cat("\n", strrep("=", 65), "\n", sep = "")
    cat(sprintf("TERMIN\u00c9. Figures dans : %s/\n", out_dir))
    cat(strrep("=", 65), "\n")
    invisible(NULL)
  }
  
  # ---------------------------------------------------------------------
  # EXECUTION EN LIGNE DE COMMANDE
  # ---------------------------------------------------------------------
  
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) {
    main(args[1], args[2])
  } else {
    main("C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_scores_harmonic",
         "C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_evolution_and_ratio")
  }
  cat("\n[3/3] Ratio Bloc4/Bloc1 (appari\u00e9 par participant)...\n")
  ratio_df <- compute_block_ratio(lr_df, score_col, block_num = 4, block_denom = 1)
  if (nrow(ratio_df) > 0) {
    ratio_path <- file.path(out_dir, "block4_over_block1_ratios.csv")
    utils::write.csv(ratio_df, ratio_path, row.names = FALSE)
    cat(sprintf("  Ratios \u00e9crits dans : %s\n", ratio_path))
    cat("  Figure 3 : violons ratio Bloc4/Bloc1, S\u00e9q. A vs B...\n")
    plot_violin_block_ratio(ratio_df, score_col, block_num = 4, block_denom = 1, out_dir)
  } else {
    cat("  Aucun ratio calculable (Bloc 1 et/ou Bloc 4 absents ou insuffisants).\n")
  }

  cat("\n", strrep("=", 65), "\n", sep = "")
  cat(sprintf("TERMIN\u00c9. Figures dans : %s/\n", out_dir))
  cat(strrep("=", 65), "\n")
  invisible(NULL)
}

# ---------------------------------------------------------------------
# EXECUTION EN LIGNE DE COMMANDE
# ---------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  main(args[1], args[2])
} else {
  main("C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_scores_harmonic",
       "C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_evolution_and_ratio")
}
