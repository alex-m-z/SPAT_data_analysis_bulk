# =====================================================================
# SPAT — Figures harmonic raw/norm + modèles Kalman, PAR PARTICIPANT
# =====================================================================
# Reprend la logique de spat_plot_evolution.R, restreinte à un
# participant à la fois, et produit 2 figures :
#
#   1) harmonic_raw_<id>.png
#        Points/trajectoire du score harmonique BRUT
#        (score_harmonic_raw), par séquence, avec frontières de bloc.
#
#   2) harmonic_norm_<id>.png
#        Points/trajectoire du score harmonique NORMALISÉ
#        (score_harmonic_norm_participant par défaut), par séquence,
#        + superposition des 2 modèles état-espace ajustés par filtre
#        de Kalman :
#          - modèle d'ADAPTATION SENSORIMOTRICE (modèle standard) :
#                x[n+1] = A * x[n] + B * e[n]
#                y[n]   = x[n] + eps[n]
#            où x[n] = état d'adaptation interne, A = facteur de
#            rétention, B = sensibilité à l'erreur (taux
#            d'apprentissage), e[n] = erreur de prédiction sensorielle
#            à l'essai n. CONVENTION RETENUE ICI : e[n] = y[n],
#            c'est-à-dire que l'observation (le score) sert elle-même
#            de signal d'erreur qui pilote la mise à jour de l'état
#            interne (faute de mesure séparée de l'erreur de
#            prédiction sensorielle dans les données). Comme le filtre
#            prédit x[t] à partir de x[t-1], c'est e[t-1] (l'erreur de
#            l'essai PRÉCÉDENT) qui est utilisée ; avant le premier
#            essai d'un bloc, il n'y a pas d'erreur antérieure
#            (e_prev = 0).
#          - modèle "racine carrée" (tendance d'apprentissage) :
#              y[n] ~= A_n*sqrt(n) + B_n
#        (les deux réinitialisés à chaque bloc, comme dans le script
#        original, et ajustés séparément par séquence.)
#
# USAGE
# -----
#   Rscript spat_plot_participant37.R scores_folder out_dir [participant_id]
#
#   - Avec un participant_id fourni : ne traite QUE ce participant,
#     figures directement dans out_dir.
#       Rscript spat_plot_participant37.R scores_folder out_dir SPAT009
#
#   - SANS participant_id : détecte automatiquement TOUS les fichiers
#     SPAT_*_scores.xlsx présents dans scores_folder et génère les
#     figures pour chacun (dans out_dir/participant_<id>/), en
#     sautant silencieusement (avec message) tout participant dont
#     le traitement échoue (fichier illisible, colonnes manquantes,
#     pas assez de points, etc.) plutôt que d'arrêter tout le run.
#       Rscript spat_plot_participant37.R scores_folder out_dir
#
# Packages requis : readxl
# =====================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(tools)
  library(svglite)
})

# ---------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------

SEQ_COLORS <- c(A = "#006880", B = "#D52B1E")
SEQ_LABELS <- c(A = "Space-first", B = "Time-first")

# Colonne utilisée pour la version "normalisée". Changez pour
# "score_harmonic_norm_global" si vous voulez la normalisation globale
# plutôt que la normalisation par participant.
NORM_COL <- "score_harmonic_norm_participant"
RAW_COL  <- "score_harmonic_raw"

SCORE_META <- list(
  score_harmonic_raw = list(
    label  = "Harmonic score",
    ylabel = "Harmonic score (raw)"
  ),
  score_harmonic_norm_participant = list(
    label  = "Harmonic score (normalized)",
    ylabel = "Harmonic score (normalized)"
  ),
  score_harmonic_norm_global = list(
    label  = "Harmonic score (normalized globally)",
    ylabel = "Harmonic score (normalized)"
  )
)

# ---------------------------------------------------------------------
# 1. CHARGEMENT DES DONNEES (un seul participant)
# ---------------------------------------------------------------------

read_participant_file <- function(path) {
  sheets <- readxl::excel_sheets(path)
  sheet  <- if ("scores" %in% sheets) "scores" else sheets[1]
  df     <- as.data.frame(readxl::read_excel(path, sheet = sheet))
  
  stem        <- tools::file_path_sans_ext(basename(path))
  participant <- sub("^SPAT_", "", sub("_scores$", "", stem))
  df$participant <- participant
  df
}

# Liste tous les participants disponibles dans le dossier, à partir
# des noms de fichiers SPAT_<id>_scores.xlsx (insensible à la casse).
list_all_participants <- function(scores_folder) {
  files <- list.files(scores_folder, pattern = "^SPAT_.*_scores\\.xlsx$",
                      full.names = FALSE, ignore.case = TRUE)
  if (length(files) == 0) return(character(0))
  ids <- sub("^SPAT_", "", files, ignore.case = TRUE)
  ids <- sub("_scores\\.xlsx$", "", ids, ignore.case = TRUE)
  sort(unique(ids))
}

load_one_participant <- function(scores_folder, participant_id) {
  fname <- file.path(scores_folder, sprintf("SPAT_%s_scores.xlsx", participant_id))
  if (!file.exists(fname)) {
    # recherche insensible à la casse, au cas où
    candidates <- list.files(scores_folder,
                             pattern = sprintf("^SPAT_%s_scores\\.xlsx$", participant_id),
                             full.names = TRUE, ignore.case = TRUE)
    if (length(candidates) == 0) {
      stop(sprintf("Fichier introuvable pour le participant %s dans %s", participant_id, scores_folder))
    }
    fname <- candidates[1]
  }
  
  df <- read_participant_file(fname)
  df <- df[order(df$sequence, df$block, df$trial), ]
  grp <- interaction(df$sequence, drop = TRUE)
  df$trial_global <- ave(seq_len(nrow(df)), grp, FUN = seq_along)
  df
}

get_block_boundaries <- function(sub) {
  blocks <- sort(unique(sub$block))
  blocks <- blocks[!is.na(blocks)]
  if (length(blocks) <= 1) return(numeric(0))
  
  boundaries <- c()
  for (b in blocks[-length(blocks)]) {
    bd <- sub[sub$block == b, ]
    boundaries <- c(boundaries, max(bd$trial_global, na.rm = TRUE))
  }
  boundaries
}

# ---------------------------------------------------------------------
# 2. FILTRE DE KALMAN GENERIQUE + LISSEUR RTS
# ---------------------------------------------------------------------
# NOTE : b_fun est maintenant une FONCTION de l'indice t (comme H_fun),
# et non plus un vecteur constant. Cela permet un terme d'entrée
# VARIABLE DANS LE TEMPS (nécessaire pour le modèle d'adaptation
# sensorimotrice x[n+1] = A*x[n] + B*e[n], où e[n] change à chaque
# essai). Pour un modèle sans entrée (ex. modèle "racine carrée"), on
# passe simplement b_fun <- function(t) c(0, 0).

kf_run <- function(y, Phi, b_fun, Q, H_fun, R, x0, P0) {
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
    xp <- as.numeric(Phi %*% x_prev + b_fun(t))
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
# 3. MODELE (1) — RETOUR VERS L'ASYMPTOTE  x[n+1] = mu + A*(x[n]-mu) + w[n]
# ---------------------------------------------------------------------
# x[n]  = état interne (latent), qui converge vers l'asymptote mu
# A     = vitesse d'apprentissage (proche de 1 = convergence lente)
# mu    = performance asymptotique attendue
# w[n]  = bruit de processus, variance Q
# y[n]  = x[n] + eps[n]   (observation), eps[n] variance R
#
# Reformulé sous forme affine standard x[n+1] = Phi*x[n] + b(t) avec
# Phi = A et b(t) = (1-A)*mu (CONSTANT dans le temps), ce qui permet de
# réutiliser tel quel le filtre/lisseur générique kf_run / kf_smooth.
# A, mu, Q, R estimés par maximum de vraisemblance, puis état lissé
# (RTS). Retourne NULL si trop peu de points ou échec numérique.

fit_kf_linear <- function(y) {
  y <- as.numeric(y)
  n_valid <- sum(!is.na(y))
  if (n_valid < 5) return(NULL)
  
  y_obs <- y[!is.na(y)]
  var_y <- stats::var(y_obs)
  if (!is.finite(var_y) || var_y <= 0) var_y <- 1
  
  mu0       <- mean(utils::tail(y_obs, max(2, floor(length(y_obs) * 0.3))))
  init_A    <- 0.8
  init_mu   <- mu0
  init_logQ <- log(max(var_y * 0.1, 1e-4))
  init_logR <- log(max(var_y * 0.5, 1e-4))
  
  x0 <- c(y_obs[1])
  P0 <- matrix(var_y * 10, 1, 1)
  H_fun <- function(t) matrix(1, 1, 1)
  
  # b_fun(t) = (1-A)*mu, constant pour un couple (A, mu) donné
  make_b_fun <- function(A, mu) {
    b_val <- (1 - A) * mu
    function(t) matrix(b_val, 1, 1)
  }
  
  negloglik <- function(par) {
    A <- par[1]; mu <- par[2]; Q <- exp(par[3]); R <- exp(par[4])
    Phi <- matrix(A, 1, 1)
    out <- tryCatch(kf_run(y, Phi, make_b_fun(A, mu), matrix(Q, 1, 1), H_fun, R, x0, P0),
                    error = function(e) NULL)
    if (is.null(out) || !is.finite(out$loglik)) return(1e10)
    -out$loglik
  }
  
  # A borné dans (0,1) comme Smith et al. (2004). mu libre (asymptote
  # potentiellement hors de la plage observée). Q/R bornés en fonction
  # de la variance des données, comme dans le script de comparaison.
  lower <- c(0.001, -Inf, log(var_y * 1e-6), log(var_y * 0.05))
  upper <- c(0.999,  Inf, log(var_y * 0.05), log(var_y * 10))
  
  fit <- tryCatch(
    stats::optim(c(init_A, init_mu, init_logQ, init_logR), negloglik,
                 method = "L-BFGS-B", lower = lower, upper = upper,
                 control = list(maxit = 500)),
    error = function(e) NULL
  )
  
  if (is.null(fit) || !is.finite(fit$value)) {
    A <- init_A; mu <- init_mu; Q <- exp(init_logQ); R <- exp(init_logR)
  } else {
    A <- fit$par[1]; mu <- fit$par[2]; Q <- exp(fit$par[3]); R <- exp(fit$par[4])
  }
  
  Phi    <- matrix(A, 1, 1)
  result <- tryCatch({
    kf_out <- kf_run(y, Phi, make_b_fun(A, mu), matrix(Q, 1, 1), H_fun, R, x0, P0)
    sm     <- kf_smooth(kf_out, Phi)
    se_vec <- sqrt(vapply(sm$P_smooth, function(P) P[1, 1], numeric(1)))
    list(A = A, mu = mu, Q = Q, R = R,
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
  b_fun <- function(t) c(0, 0)   # pas d'entrée exogène pour ce modèle
  x0    <- c(A0, B0)
  P0    <- diag(c(var_y * 10, var_y * 10))
  H_fun <- function(t) matrix(c(sqrt(trial_n[t]), 1), nrow = 1)
  
  negloglik <- function(par) {
    Qa <- exp(par[1]); Qb <- exp(par[2]); R <- exp(par[3])
    Qm <- diag(c(Qa, Qb))
    out <- tryCatch(kf_run(y, Phi, b_fun, Qm, H_fun, R, x0, P0), error = function(e) NULL)
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
    kf_out <- kf_run(y, Phi, b_fun, Qm, H_fun, R, x0, P0)
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
# 5. FIGURE 1 — Score harmonique BRUT, participant unique
# ---------------------------------------------------------------------

plot_raw_participant <- function(df, participant_id, out_dir) {
  sub <- df[!is.na(df[[RAW_COL]]), ]
  if (nrow(sub) == 0) {
    cat(sprintf("  Aucune donn\u00e9e pour %s, figure ignor\u00e9e.\n", RAW_COL))
    return(NULL)
  }
  
  meta <- SCORE_META[[RAW_COL]]
  sequences <- sort(unique(stats::na.omit(sub$sequence)))
  
  yr  <- range(sub[[RAW_COL]], na.rm = TRUE)
  pad <- diff(yr) * 0.08
  if (!is.finite(pad) || pad == 0) pad <- 1
  yr  <- yr + c(-pad, pad)
  xr  <- range(sub$trial_global, na.rm = TRUE)
  
  fname <- file.path(out_dir, sprintf("harmonic_raw_%s.pdf", participant_id))
  svglite::svglite(fname, width = 15, height = 7)
  par(mar = c(5.5, 5.5, 4, 4),
      mgp = c(3.4, 1.3, 0)
  )
  
  plot(NA, xlim = xr, ylim = yr, xlab = "Trial", ylab = meta$ylabel,
       main = "Raw harmonic score", 
       cex.lab = 1.8,   # x/y labels
       cex.axis = 1.5,  # tick labels
       cex.main = 2.5,   # title
       font.main = 2)
  grid(col = "gray90", lty = 1)
  
  leg_labels <- c(); leg_colors <- c()
  
  ref_seq <- sequences[1]
  gm_ref  <- sub[sub$sequence == ref_seq, ]
  if (nrow(gm_ref) > 0 && "block" %in% names(gm_ref)) {
    blocks_ref    <- sort(unique(stats::na.omit(gm_ref$block)))
    y_label_block <- yr[2] - diff(yr) * 0.03
    for (b in blocks_ref) {
      idx <- gm_ref$block == b & !is.na(gm_ref$block)
      if (!any(idx)) next
      xmid <- mean(range(gm_ref$trial_global[idx], na.rm = TRUE))
      text(xmid, y_label_block, labels = sprintf("Block %d", b),
           col = "gray30", cex = 1.5, font = 2)
    }
  }
  
  for (seq_ in sequences) {
    sd_ <- sub[sub$sequence == seq_, ]
    sd_ <- sd_[order(sd_$trial_global), ]
    col_seq <- if (seq_ %in% names(SEQ_COLORS)) SEQ_COLORS[[seq_]] else "gray30"
    
    for (bd in get_block_boundaries(sd_)) abline(v = bd, col = "dimgray", lty = 2, lwd = 1.0)
    
    lines(sd_$trial_global, sd_[[RAW_COL]], col = col_seq, lwd = 1.8)
    points(sd_$trial_global, sd_[[RAW_COL]], col = col_seq, pch = 16, cex = 0.7)
    
    leg_labels <- c(leg_labels, SEQ_LABELS[[seq_]])
    leg_colors <- c(leg_colors, col_seq)
  }

  
  legend(
    "bottomleft",
    legend = leg_labels,
    col = leg_colors,
    lty = 1,
    lwd = 1.8,
    pch = 16,
    cex = 1,
    bty = "o",
    inset = 0.02
  )
  
  grDevices::dev.off()
  cat(sprintf("  Sauvegard\u00e9 : %s\n", fname))
  fname
}

# ---------------------------------------------------------------------
# 6. FIGURE 2 — Score harmonique NORMALISE + mod\u00e8les Kalman
# ---------------------------------------------------------------------

plot_norm_with_models_participant <- function(df, participant_id, out_dir, ci_z = 1.96) {
  sub <- df[!is.na(df[[NORM_COL]]), ]
  if (nrow(sub) == 0) {
    cat(sprintf("  Aucune donn\u00e9e pour %s, figure ignor\u00e9e.\n", NORM_COL))
    return(NULL)
  }
  
  meta <- SCORE_META[[NORM_COL]]
  sequences <- sort(unique(stats::na.omit(sub$sequence)))
  
  # --- Ajuste les deux mod\u00e8les par s\u00e9quence (r\u00e9initialis\u00e9s par bloc) ---
  fits <- list()
  for (seq_ in sequences) {
    sd_ <- sub[sub$sequence == seq_, ]
    sd_ <- sd_[order(sd_$trial_global), ]
    fits[[seq_]] <- list(
      data = sd_,
      lin  = fit_kf_linear_by_block(sd_[[NORM_COL]], sd_$block),
      sqrt = fit_kf_sqrt_by_block(sd_[[NORM_COL]], sd_$block)
    )
  }
  
  # --- Bornes des axes (donn\u00e9es brutes + IC des mod\u00e8les) ---
  all_y <- sub[[NORM_COL]]
  for (seq_ in sequences) {
    fl <- fits[[seq_]]$lin; fs <- fits[[seq_]]$sqrt
    if (any(!is.na(fl$x_smooth))) all_y <- c(all_y, fl$x_smooth - ci_z * fl$se, fl$x_smooth + ci_z * fl$se)
    if (any(!is.na(fs$x_smooth))) all_y <- c(all_y, fs$x_smooth - ci_z * fs$se, fs$x_smooth + ci_z * fs$se)
  }
  yr  <- range(all_y, na.rm = TRUE)
  pad <- diff(yr) * 0.1
  if (!is.finite(pad) || pad == 0) pad <- 1
  yr  <- yr + c(-pad, pad)
  xr  <- range(sub$trial_global, na.rm = TRUE)
  
  fname <- file.path(out_dir, sprintf("harmonic_norm_%s.pdf", participant_id))
  svglite::svglite(fname, width = 15, height = 7)
  par(mar = c(5.5, 5.5, 4, 4),
      mgp = c(3.4, 1.3, 0)
  )
  
  plot(NA, xlim = xr, ylim = yr, xlab = "Trial", ylab = meta$ylabel,
       main = "Normalized harmonic score",
       cex.lab = 1.8,   # x/y labels
       cex.axis = 1.5,  # tick labels
       cex.main = 2.5,   # title
       font.main = 2)
  grid(col = "gray90", lty = 1)
  
  leg_labels <- c(); leg_colors <- c(); leg_lty <- c(); leg_lwd <- c()
  
  ref_seq <- sequences[1]
  gm_ref  <- fits[[ref_seq]]$data
  if (nrow(gm_ref) > 0 && "block" %in% names(gm_ref)) {
    blocks_ref    <- sort(unique(stats::na.omit(gm_ref$block)))
    y_label_block <- yr[2] - diff(yr) * 0.03
    for (b in blocks_ref) {
      idx <- gm_ref$block == b & !is.na(gm_ref$block)
      if (!any(idx)) next
      xmid <- mean(range(gm_ref$trial_global[idx], na.rm = TRUE))
      text(xmid, y_label_block, labels = sprintf("Block %d", b),
           col = "gray30", cex = 1.5, font = 2)
    }
  }
  
  for (seq_ in sequences) {
    sd_     <- fits[[seq_]]$data
    col_seq <- if (seq_ %in% names(SEQ_COLORS)) SEQ_COLORS[[seq_]] else "gray30"
    
    for (bd in get_block_boundaries(sd_)) abline(v = bd, col = "dimgray", lty = 2, lwd = 1.0)
    
    # --- Numéros de blocs (une seule fois, basés sur la 1re séquence) ---
    
    
    # points/ligne brute normalisée
    lines(sd_$trial_global, sd_[[NORM_COL]], col = grDevices::adjustcolor(col_seq, alpha.f = 0.45), lwd = 1.2)
    points(sd_$trial_global, sd_[[NORM_COL]], col = grDevices::adjustcolor(col_seq, alpha.f = 0.6), pch = 16, cex = 0.6)
    
    # modèle d'adaptation sensorimotrice x[n+1] = A*x[n] + B*e[n], avec IC
    fl <- fits[[seq_]]$lin
    if (any(!is.na(fl$x_smooth))) {
      ok <- is.finite(fl$x_smooth) & is.finite(fl$se)
      xx <- sd_$trial_global[ok]
      lo <- fl$x_smooth[ok] - ci_z * fl$se[ok]
      hi <- fl$x_smooth[ok] + ci_z * fl$se[ok]
      polygon(c(xx, rev(xx)), c(lo, rev(hi)),
              col = grDevices::adjustcolor(col_seq, alpha.f = 0.12), border = NA)
      lines(sd_$trial_global, fl$x_smooth, col = col_seq, lwd = 2.4, lty = 2)
    }
    
    # modèle racine carrée y[n] ~= A*sqrt(n) + B, avec IC
    fs <- fits[[seq_]]$sqrt
    if (any(!is.na(fs$x_smooth))) {
      ok <- is.finite(fs$x_smooth) & is.finite(fs$se)
      xx <- sd_$trial_global[ok]
      lo <- fs$x_smooth[ok] - ci_z * fs$se[ok]
      hi <- fs$x_smooth[ok] + ci_z * fs$se[ok]
      polygon(c(xx, rev(xx)), c(lo, rev(hi)),
              col = grDevices::adjustcolor(col_seq, alpha.f = 0.12), border = NA)
      lines(sd_$trial_global, fs$x_smooth, col = col_seq, lwd = 2.4, lty = 3)
    }
    
    leg_labels <- c(leg_labels,
                    sprintf("%s data", SEQ_LABELS[[seq_]]),
                    sprintf("State-space model (%s)", SEQ_LABELS[[seq_]]),
                    sprintf("Square-root model (%s)", SEQ_LABELS[[seq_]]))
    leg_colors <- c(leg_colors, col_seq, col_seq, col_seq)
    leg_lty    <- c(leg_lty, 1, 2, 3)
    leg_lwd    <- c(leg_lwd, 1.2, 2.4, 2.4)
  }
  
  legend(
    "bottomright",
    legend = leg_labels,
    col = leg_colors,
    lty = leg_lty,
    lwd = leg_lwd,
    ncol = 2,
    cex = 1,
    bty = "o",
  )
  
  grDevices::dev.off()
  cat(sprintf("  Sauvegard\u00e9 : %s\n", fname))
  fname
}

# ---------------------------------------------------------------------
# 7. PIPELINE PRINCIPAL
# ---------------------------------------------------------------------

main <- function(scores_folder, out_dir = "./output_plots_37", participant_id = "37") {
  cat(strrep("=", 65), "\n")
  cat(sprintf("SPAT \u2014 FIGURES PARTICIPANT %s\n", participant_id))
  cat(strrep("=", 65), "\n")
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat(sprintf("\nChargement des donn\u00e9es du participant %s depuis : %s\n", participant_id, scores_folder))
  df <- load_one_participant(scores_folder, participant_id)
  cat(sprintf("  %d lignes charg\u00e9es.\n", nrow(df)))
  
  cat("\nFigure 1 : score harmonique brut...\n")
  f1 <- plot_raw_participant(df, participant_id, out_dir)
  
  cat("\nFigure 2 : score harmonique normalis\u00e9 + mod\u00e8les Kalman...\n")
  f2 <- plot_norm_with_models_participant(df, participant_id, out_dir)
  
  cat("\n", strrep("=", 65), "\n", sep = "")
  cat("TERMIN\u00c9.\n")
  cat(strrep("=", 65), "\n")
  invisible(c(f1, f2))
}

# ---------------------------------------------------------------------
# 8. TRAITEMENT AUTOMATIQUE DE TOUS LES PARTICIPANTS (avec skip)
# ---------------------------------------------------------------------
# Détecte tous les fichiers SPAT_<id>_scores.xlsx du dossier et lance
# main() pour chacun. Chaque participant est isolé dans son propre
# tryCatch : si un participant échoue (fichier corrompu, colonnes
# manquantes, pas assez de trials valides, etc.), un message est
# affiché et on passe au suivant SANS interrompre le run global.
# Les figures de chaque participant vont dans out_dir_base/participant_<id>/.

run_all_participants <- function(scores_folder, out_dir_base) {
  ids <- list_all_participants(scores_folder)
  if (length(ids) == 0) {
    cat(sprintf("Aucun fichier SPAT_*_scores.xlsx trouv\u00e9 dans %s\n", scores_folder))
    return(invisible(NULL))
  }
  
  cat(sprintf("%d participant(s) d\u00e9tect\u00e9(s) : %s\n\n", length(ids), paste(ids, collapse = ", ")))
  
  ok_ids   <- c()
  skip_ids <- c()
  
  for (pid in ids) {
    out_dir <- file.path(out_dir_base, sprintf("participant_%s", pid))
    cat(strrep("-", 65), "\n")
    result <- tryCatch({
      main(scores_folder, out_dir, pid)
      "ok"
    }, error = function(e) {
      cat(sprintf("  [SKIP] Participant %s ignor\u00e9 : %s\n", pid, conditionMessage(e)))
      "skip"
    }, warning = function(w) {
      # Une warning ne doit pas arrêter le participant ; on la note et on continue.
      cat(sprintf("  [WARN] Participant %s : %s\n", pid, conditionMessage(w)))
      "ok"
    })
    if (identical(result, "ok")) ok_ids <- c(ok_ids, pid) else skip_ids <- c(skip_ids, pid)
  }
  
  cat("\n", strrep("=", 65), "\n", sep = "")
  cat(sprintf("BILAN : %d participant(s) trait\u00e9(s), %d ignor\u00e9(s).\n", length(ok_ids), length(skip_ids)))
  if (length(ok_ids)   > 0) cat(sprintf("  Trait\u00e9s : %s\n", paste(ok_ids, collapse = ", ")))
  if (length(skip_ids) > 0) cat(sprintf("  Ignor\u00e9s : %s\n", paste(skip_ids, collapse = ", ")))
  cat(strrep("=", 65), "\n")
  
  invisible(list(ok = ok_ids, skipped = skip_ids))
}

# ---------------------------------------------------------------------
# EXECUTION EN LIGNE DE COMMANDE
# ---------------------------------------------------------------------
#   Rscript spat_plot_participant37.R scores_folder out_dir            -> TOUS les participants (skip auto)
#   Rscript spat_plot_participant37.R scores_folder out_dir SPAT009    -> UN seul participant

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  main(args[1], args[2], args[3])
} else if (length(args) == 2) {
  run_all_participants(args[1], args[2])
} else {
  main("C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_scores_harmonic",
       "C:/Users/Alex/Documents/Stage_IRCN/SPAT/plots_poster",
       "SPAT001")
}