#!/usr/bin/env Rscript
# =====================================================================
# SPAT — Visualisation des courbes d'évolution (R)
# =====================================================================
# Reproduit spat_plot_evolution.py :
#   - lit les fichiers SPAT_<participant>_scores.xlsx d'un dossier
#   - pour CHAQUE score x CHAQUE séquence, trace :
#       * la trajectoire de chaque participant (ligne fine, transparente)
#       * la trajectoire moyenne du groupe (ligne épaisse)
#       * des lignes verticales aux frontières de bloc (calculées à
#         partir des données, pas de seuils fixes)
#
# EN PLUS, pour chaque score x séquence, ajuste par FILTRE DE KALMAN
# deux modèles état-espace linéaires-gaussiens :
#
#   (1) Modèle linéaire autorégressif :
#         x[n+1] = A * x[n] + B + w[n]     w[n] ~ N(0, Q)
#         y[n]   = x[n] + v[n]             v[n] ~ N(0, R)
#       -> "lisse" la trajectoire d'un score par filtrage/lissage de
#          Kalman (filtre + lisseur RTS). A, B, Q, R estimés par
#          maximum de vraisemblance (décomposition en erreurs de
#          prédiction du filtre de Kalman).
#
#   (2) Modèle de tendance "racine carrée" (courbe d'apprentissage) :
#         y[n] ≈ A_n * sqrt(n) + B_n
#       représenté comme un modèle état-espace où l'état est le couple
#       de coefficients [A_n, B_n], qui évolue en marche aléatoire lente
#       (théta[n] = théta[n-1] + w[n]) et où la matrice d'observation
#       H[n] = [sqrt(n), 1] est VARIABLE DANS LE TEMPS. C'est le filtre
#       de Kalman utilisé comme régression récursive à coefficients
#       variables (même machinerie filtre+lisseur que le modèle (1),
#       mais état de dimension 2 et H[n] dépendant de n).
#
# Les deux modèles sont ajustés à la fois sur la moyenne du groupe et
# sur chaque participant. Les courbes lissées du GROUPE sont superposées
# à la figure d'évolution principale ; les courbes par participant (+
# groupe, pour comparaison) sont tracées dans deux figures dédiées.
#
# SORTIES (dans out_dir), pour chaque score x séquence :
#   evolution_<score>_seq<seq>.png          (figure principale + overlay)
#   state_space_linear_<score>_seq<seq>.png (modèle (1), par participant)
#   state_space_sqrt_<score>_seq<seq>.png   (modèle (2), par participant)
# + state_space_params.csv : A,B,Q,R (ou Qa,Qb,R) par score/séquence/sujet
#
# USAGE
# -----
#   Rscript spat_plot_evolution.R scores_folder out_dir
#
# Packages requis : readxl
# =====================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(tools)
})

# ---------------------------------------------------------------------
# CONFIGURATION (identique à spat_plot_evolution.py)
# ---------------------------------------------------------------------

SEQ_COLORS <- c(A = "#4C8CC5", B = "#F1BC18")

SCORE_META <- list(
  
  score_harmonic_norm_participant = list(
    label  = "Harmonic score (normailzed by partipant)",
    ylabel = "Harmonic score (normailzed by partipant)"
  )
)

# ---------------------------------------------------------------------
# 1. CHARGEMENT DES DONNEES
# ---------------------------------------------------------------------

# Concatène une liste de data.frames aux colonnes potentiellement
# différentes (remplit les colonnes manquantes avec NA), pour éviter
# une dépendance à dplyr::bind_rows.
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

# Charge tous les SPAT_<participant>_scores.xlsx d'un dossier, les
# concatène, et ajoute 'trial_global' : la position chronologique
# réelle (1, 2, 3, ...) de chaque essai pour chaque participant x
# séquence, triée sur (block, trial).
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

# Frontières de bloc : pour chaque bloc (sauf le dernier), la médiane
# (sur les participants) de la dernière position chronologique
# ('trial_global') atteinte dans ce bloc. Calculé à partir des données
# réelles plutôt que de seuils fixes.
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
# État de dimension d quelconque (d=1 pour le modèle linéaire,
# d=2 pour le modèle "racine carrée" à coefficients variables).
#
#   x[t] = Phi %*% x[t-1] + b + w[t]      w[t] ~ N(0, Q)
#   y[t] = H(t) %*% x[t] + v[t]           v[t] ~ N(0, R)   (obs. scalaire)
#
# H_fun : fonction de l'indice t -> matrice ligne (1 x d)
# Les observations manquantes (NA) sont gérées : seule l'étape de
# prédiction est appliquée (pas de mise à jour).

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
    xp <- as.numeric(Phi %*% x_prev + b_fun)
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
      
      K <- (Pp %*% t(Ht)) / S            # d x 1
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

# Lisseur de Rauch-Tung-Striebel (RTS)
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
# A, B, Q, R estimés par maximum de vraisemblance (via kf_run), puis
# la série est lissée (RTS). Retourne NULL si trop peu de points ou en
# cas d'échec numérique (le graphique ignore alors ce sujet).
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
  
  # b constant (pas une fonction de t) : b = (1-A)*mu, comme attendu par kf_run
  b_vec <- function(A, mu) matrix((1 - A) * mu, 1, 1)
  
  negloglik <- function(par) {
    A <- par[1]; mu <- par[2]; Q <- exp(par[3]); R <- exp(par[4])
    Phi <- matrix(A, 1, 1)
    out <- tryCatch(kf_run(y, Phi, b_vec(A, mu), matrix(Q, 1, 1), H_fun, R, x0, P0),
                    error = function(e) NULL)
    if (is.null(out) || !is.finite(out$loglik)) return(1e10)
    -out$loglik
  }
  
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
    kf_out <- kf_run(y, Phi, b_vec(A, mu), matrix(Q, 1, 1), H_fun, R, x0, P0)
    sm     <- kf_smooth(kf_out, Phi)
    se_vec <- sqrt(vapply(sm$P_smooth, function(P) P[1, 1], numeric(1)))
    list(A = A, mu = mu, Q = Q, R = R,
         x_filt   = as.numeric(kf_out$x_filt[, 1]),
         x_smooth = as.numeric(sm$x_smooth[, 1]),
         se       = se_vec,
         loglik   = kf_out$loglik,
         n_params = 4)
  }, error = function(e) NULL)
  
  result
}

# ---------------------------------------------------------------------
# 4. MODELE (2) — TENDANCE "RACINE CARREE"  y[n] ~= A_n*sqrt(n) + B_n
# ---------------------------------------------------------------------
# État théta[n] = [A_n, B_n], marche aléatoire (Phi = I2), observation
# H[n] = [sqrt(n), 1] (VARIABLE dans le temps). Qa, Qb (bruit de
# process sur A_n et B_n) et R estimés par maximum de vraisemblance ;
# état initial théta0 par régression y ~ sqrt(n) sur les points valides.

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
  b_fun <- c(0, 0)
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
  
  # --- inside fit_kf_sqrt, replace the optim() call ---
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
         A_t = sm$x_smooth[, 1], B_t = sm$x_smooth[, 2],
         x_smooth = fitted_mean,
         se       = se_vec,
         loglik   = kf_out$loglik,     # <-- idem
         n_params = 3)                 # Qa, Qb, R
  }, error = function(e) NULL)
  
  result
}

# ---------------------------------------------------------------------
# Variantes "par bloc" : réinitialisent l'état (x0, P0) au début de
# CHAQUE bloc plutôt que de laisser le filtre accumuler l'historique
# sur toute la séquence. Les segments ajustés sont recollés dans
# l'ordre d'origine du vecteur y.
# ---------------------------------------------------------------------

fit_kf_linear_by_block <- function(y, block) {
  n <- length(y)
  x_smooth <- rep(NA_real_, n)
  se       <- rep(NA_real_, n)
  params   <- list()
  for (b in sort(unique(block))) {
    idx <- which(block == b)
    fl  <- fit_kf_linear(y[idx])
    if (is.null(fl)) next
    x_smooth[idx] <- fl$x_smooth
    se[idx]       <- fl$se
    params[[length(params) + 1]] <- data.frame(block = b, A = fl$A, mu = fl$mu, Q = fl$Q, R = fl$R)
  }
  list(x_smooth = x_smooth, se = se,
       params   = if (length(params)) do.call(rbind, params) else NULL)
}

fit_kf_sqrt_by_block <- function(y, block) {
  n <- length(y)
  x_smooth <- rep(NA_real_, n)
  se       <- rep(NA_real_, n)                            # <-- NOUVEAU
  params   <- list()
  for (b in sort(unique(block))) {
    idx     <- which(block == b)
    n_local <- seq_along(idx)
    fs <- fit_kf_sqrt(y[idx], n_local)
    if (is.null(fs)) next
    x_smooth[idx] <- fs$x_smooth
    se[idx]       <- fs$se                                # <-- NOUVEAU
    params[[length(params) + 1]] <- data.frame(block = b, Qa = fs$Qa, Qb = fs$Qb, R = fs$R)
  }
  list(x_smooth = x_smooth, se = se,                        # <-- se ajouté
       params   = if (length(params)) do.call(rbind, params) else NULL)
}


# ---------------------------------------------------------------------
# 4bis. LEARNING RATES (façon spat_learning_rate_comparison.py) +
#       VIOLONS Séquence A vs Séquence B
# ---------------------------------------------------------------------

# Violon "maison" : densité miroir + médiane/IQR, sans dépendance externe.
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

# --- Indicateurs de learning rate, identiques à la définition Python ---
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

# Une ligne par (score, participant, sequence, block) : les 4 LR +
# les paramètres A, B du modèle linéaire (utilisés ensuite pour
# les violons). Réutilise fit_kf_linear / fit_kf_sqrt (bornés, cf.
# la correction précédente) plutôt qu'une simple régression moindres
# carrés, pour rester cohérent avec le reste du script R.
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
          y_ssm <- fl$x_smoot
          
          fs <- fit_kf_sqrt(y, seq_along(y))
          lr_sqrt_slope <- if (!is.null(fs)) learning_rate_max_slope(fs$x_smooth) else NA_real_
          
          rows[[length(rows) + 1]] <- data.frame(
            score = score_col, participant = part, sequence = seq_, block = b,
            lr_ssm_1090   = learning_rate_10_90(y_ssm, bd$trial),
            lr_sqrt_slope = lr_sqrt_slope,
            lr_ssm_slope  = learning_rate_max_slope(y_ssm),
            lr_ssm_A      = fl$A,
            lr_ssm_mu     = fl$mu,
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
# 4bis-bis. RATIO Séquence A / Séquence B, pour chaque définition de
#           learning rate (une ligne par score x bloc x participant)
# ---------------------------------------------------------------------
# Pour chaque métrique de lr_df (lr_ssm_1090, lr_sqrt_slope, lr_ssm_slope,
# lr_ssm_A, lr_ssm_B), calcule ratio_<metrique> = valeur(A) / valeur(B),
# appariée par participant (même score, même bloc). Le ratio est mis à
# NA si l'une des deux valeurs est manquante, ou si la valeur B est trop
# proche de 0 (ratio numériquement instable/non interprétable) — ceci
# est surtout pertinent pour lr_ssm_B et lr_ssm_A, qui peuvent être
# négatifs ou proches de 0 ; dans ce cas le ratio A/B n'a pas vraiment
# de sens et un ratio brut serait trompeur (valeurs énormes ou signe
# qui s'inverse). Envisager un log-ratio si cela pose problème en aval.
compute_learning_rate_ratios <- function(lr_df) {
  metrics <- c("lr_ssm_1090", "lr_sqrt_slope", "lr_ssm_slope", "lr_ssm_A", "lr_ssm_mu")
  if (nrow(lr_df) == 0) return(data.frame())
  
  rows <- list()
  for (score_col in sort(unique(lr_df$score))) {
    sub_score <- lr_df[lr_df$score == score_col, ]
    for (b in sort(unique(sub_score$block))) {
      sub <- sub_score[sub_score$block == b, ]
      participants <- sort(unique(sub$participant))
      if (length(participants) == 0) next
      
      row_df <- data.frame(score = score_col, block = b, participant = participants,
                           stringsAsFactors = FALSE)
      
      for (m in metrics) {
        a_val <- vapply(participants, function(p) {
          v <- sub[sub$participant == p & sub$sequence == "A", m]
          if (length(v) == 0) NA_real_ else v[1]
        }, numeric(1))
        b_val <- vapply(participants, function(p) {
          v <- sub[sub$participant == p & sub$sequence == "B", m]
          if (length(v) == 0) NA_real_ else v[1]
        }, numeric(1))
        
        ratio <- a_val / b_val
        unstable <- !is.finite(b_val) | abs(b_val) < 1e-8
        ratio[unstable] <- NA_real_
        
        row_df[[paste0("ratio_", m)]] <- ratio
      }
      
      rows[[length(rows) + 1]] <- row_df
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# Violons Séquence A vs Séquence B pour un score x bloc donné :
# 4 learning rates + A + B côte à côte.
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
  
  # Nombre de participants ayant au moins une valeur dans chaque
  # séquence, pour ce score x bloc (affiché dans le titre global).
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
    # las = 2 (labels verticaux) + cex.axis réduit : évite que R ne
    # supprime silencieusement un des deux labels par manque de place
    # lorsque les panneaux sont étroits (5 panneaux côte à côte).
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

# ---------------------------------------------------------------------
# 4ter. LEARNING RATES — nuage de points labellisé (A vs B, par participant)
# ---------------------------------------------------------------------
# Pour chaque métrique : un point par participant x séquence, coloré
# par séquence, étiqueté avec l'ID du participant, et relié par une
# ligne fine à sa valeur dans l'autre séquence (comparaison appariée
# visuelle). Alternative plus lisible que les violons quand le nombre
# de participants est petit (n < ~10).

plot_points_learning_rates <- function(lr_df, score_col, block_num, out_dir) {
  sub <- lr_df[lr_df$score == score_col & lr_df$block == block_num, ]
  if (nrow(sub) == 0) return(NULL)
  
  metrics <- c("lr_sqrt_slope", "lr_ssm_A")
  metric_labels <- c(
    lr_sqrt_slope = "Maximum slope (\u221a-model)",
    lr_ssm_A      = "A parameter (retention parameter)"
  )
  
  meta  <- SCORE_META[[score_col]]
  label <- if (!is.null(meta)) meta$label else score_col
  
  # Nombre de participants ayant au moins une valeur dans chaque
  # séquence, pour ce score x bloc (affiché dans le titre global).
  nA <- length(unique(sub$participant[sub$sequence == "A"]))
  nB <- length(unique(sub$participant[sub$sequence == "B"]))
  
  fname <- file.path(out_dir, sprintf("points_learning_rates_%s_block%s.png", score_col, block_num))
  grDevices::png(fname, width = 1600, height = 550, res = 130)
  par(mfrow = c(1, length(metrics)), mar = c(5, 4, 3, 4.5), oma = c(0, 0, 3, 0))
  
  for (m in metrics) {
    wide <- stats::reshape(
      sub[, c("participant", "sequence", m)],
      idvar = "participant", timevar = "sequence", direction = "wide"
    )
    names(wide) <- sub("^.*\\.", "", names(wide))  # -> participant, A, B
    wide <- wide[stats::complete.cases(wide[, c("A", "B")]) |
                   !is.na(wide$A) | !is.na(wide$B), ]
    if (nrow(wide) == 0) { plot.new(); next }
    
    all_vals <- c(wide$A, wide$B)
    all_vals <- all_vals[is.finite(all_vals)]
    yr  <- range(all_vals, na.rm = TRUE)
    pad <- diff(yr) * 0.15
    if (!is.finite(pad) || pad == 0) pad <- 1
    yr <- yr + c(-pad, pad)
    
    plot(NA, xlim = c(0.5, 2.5), ylim = yr, xaxt = "n", xlab = "",
         ylab = metric_labels[[m]], main = metric_labels[[m]], cex.main = 0.85)
    # las = 2 (labels verticaux) + cex.axis réduit : évite que R ne
    # supprime silencieusement un des deux labels par manque de place
    # lorsque les panneaux sont étroits (5 panneaux côte à côte). C'était
    # la cause du "Séquence B" invisible malgré le code correct.
    axis(1, at = c(1, 2), labels = c("S\u00e9quence A", "S\u00e9quence B"),
         las = 2, cex.axis = 0.75)
    grid(nx = NA, ny = NULL, col = "gray90")
    
    for (i in seq_len(nrow(wide))) {
      a_val <- wide$A[i]; b_val <- wide$B[i]
      if (is.finite(a_val) && is.finite(b_val)) {
        lines(c(1, 2), c(a_val, b_val), col = "gray70", lwd = 1)
      }
      if (is.finite(a_val)) {
        points(1, a_val, pch = 21, bg = SEQ_COLORS[["A"]], col = "black", cex = 1.3)
        text(1, a_val, labels = wide$participant[i], pos = 2, cex = 0.55, offset = 0.35)
      }
      if (is.finite(b_val)) {
        points(2, b_val, pch = 21, bg = SEQ_COLORS[["B"]], col = "black", cex = 1.3)
        text(2, b_val, labels = wide$participant[i], pos = 4, cex = 0.55, offset = 0.35)
      }
    }
    
    # Moyennes de groupe superposées
    mean_a <- mean(wide$A, na.rm = TRUE)
    mean_b <- mean(wide$B, na.rm = TRUE)
    segments(0.85, mean_a, 1.15, mean_a, lwd = 3, col = "black")
    segments(1.85, mean_b, 2.15, mean_b, lwd = 3, col = "black")
  }
  
  mtext(sprintf("Learning rates par participant \u2014 %s (Bloc %s) \u2014 S\u00e9q. A (n=%d) vs S\u00e9q. B (n=%d)",
                label, block_num, nA, nB),
        outer = TRUE, cex = 1.1, font = 2)
  grDevices::dev.off()
  cat(sprintf("    Sauvegard\u00e9 : %s\n", fname))
  fname
}
# ---------------------------------------------------------------------
# 4quater. MOYENNE DE GROUPE + ERREUR STANDARD PAR POSITION CHRONOLOGIQUE
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

# ---------------------------------------------------------------------
# 4quinquies. COMPARAISON D'AJUSTEMENT ENTRE LES DEUX MODELES
#              (AIC, BIC, log-vraisemblance, RMSE, correlation)
# ---------------------------------------------------------------------

model_fit_metrics <- function(y, x_smooth, loglik, n_params) {
  ok <- is.finite(y) & is.finite(x_smooth)
  n_obs <- sum(ok)
  if (n_obs < 2 || is.null(loglik) || !is.finite(loglik)) {
    return(list(loglik = NA_real_, aic = NA_real_, bic = NA_real_,
                rmse = NA_real_, cor = NA_real_, n_obs = n_obs))
  }
  resid <- y[ok] - x_smooth[ok]
  rmse  <- sqrt(mean(resid^2))
  corr  <- if (stats::sd(y[ok]) > 1e-10 && stats::sd(x_smooth[ok]) > 1e-10) {
    stats::cor(y[ok], x_smooth[ok])
  } else NA_real_
  aic <- 2 * n_params - 2 * loglik
  bic <- n_params * log(n_obs) - 2 * loglik
  list(loglik = loglik, loglik_per_obs = loglik / n_obs,   # <-- NOUVEAU
       aic = aic, bic = bic, rmse = rmse, cor = corr, n_obs = n_obs)
}

# Compare les 2 modèles sur UNE série y (un bloc, un participant/séquence)
compare_models_on_series <- function(y, label_cols = list()) {
  fl <- fit_kf_linear(y)
  fs <- fit_kf_sqrt(y, seq_along(y))
  
  rows <- list()
  if (!is.null(fl)) {
    m <- model_fit_metrics(y, fl$x_smooth, fl$loglik, fl$n_params)
    rows[[length(rows) + 1]] <- c(label_cols, model = "ssm_x=Ax+B", m)
  }
  if (!is.null(fs)) {
    m <- model_fit_metrics(y, fs$x_smooth, fs$loglik, fs$n_params)
    rows[[length(rows) + 1]] <- c(label_cols, model = "sqrt_trend", m)
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
}

# Applique la comparaison à tout le jeu de données : score x participant
# x séquence x bloc (même granularité que les learning rates).
compute_model_comparison_all <- function(df, score_cols) {
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
          res <- compare_models_on_series(
            bd[[score_col]],
            label_cols = list(score = score_col, participant = part, sequence = seq_, block = b)
          )
          if (nrow(res) > 0) rows[[length(rows) + 1]] <- res
        }
      }
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# Violons Linéaire vs √-model pour AIC/BIC/RMSE/corrélation (un score donné)
plot_model_comparison_metrics <- function(cmp_df, score_col, out_dir) {
  sub <- cmp_df[cmp_df$score == score_col, ]
  if (nrow(sub) == 0) return(NULL)
  
  metrics <- c("aic", "bic", "rmse", "cor")
  metric_labels <- c(aic = "AIC", bic = "BIC", rmse = "RMSE",
                     cor = "Correlation")
  
  meta  <- SCORE_META[[score_col]]
  label <- if (!is.null(meta)) meta$label else score_col
  
  fname <- file.path(out_dir, sprintf("model_comparison_%s.png", score_col))
  grDevices::png(fname, width = 1400, height = 500, res = 130)
  par(mfrow = c(1, length(metrics)), mar = c(3, 4, 3, 1), oma = c(0, 0, 3, 0))
  
  for (m in metrics) {
    lin_vals  <- as.numeric(sub[sub$model == "ssm_x=Ax+B", m])
    sqrt_vals <- as.numeric(sub[sub$model == "sqrt_trend", m])
    all_vals  <- c(lin_vals, sqrt_vals); all_vals <- all_vals[is.finite(all_vals)]
    if (length(all_vals) == 0) { plot.new(); next }
    
    yr  <- range(all_vals, na.rm = TRUE)
    pad <- diff(yr) * 0.1; if (!is.finite(pad) || pad == 0) pad <- 1
    yr  <- yr + c(-pad, pad)
    
    plot(NA, xlim = c(0.5, 2.5), ylim = yr, xaxt = "n", xlab = "",
         ylab = metric_labels[[m]], main = metric_labels[[m]], cex.main = 0.85)
    axis(1, at = c(1, 2), labels = c("SSM", "\u221a-model"))
    grid(nx = NA, ny = NULL, col = "gray90")
    draw_violin(lin_vals,  at = 1, col = "#4C8CC5")
    draw_violin(sqrt_vals, at = 2, col = "#F1BC18")
  }
  
  mtext(sprintf("Model comparison \u2014 %s", label), outer = TRUE, cex = 1.1, font = 2)
  grDevices::dev.off()
  cat(sprintf("    Sauvegard\u00e9 : %s\n", fname))
  fname
}

# Violons Linéaire vs √-model, BLOC 3 uniquement
plot_model_comparison_metrics_block3 <- function(cmp_df, score_col, out_dir, block_num = 3) {
  sub <- cmp_df[cmp_df$score == score_col & cmp_df$block == block_num, ]
  if (nrow(sub) == 0) return(NULL)
  
  metrics <- c("aic", "bic", "rmse", "cor")
  metric_labels <- c(aic = "AIC", bic = "BIC", rmse = "RMSE",
                     cor = "Correlation")
  
  meta  <- SCORE_META[[score_col]]
  label <- if (!is.null(meta)) meta$label else score_col
  
  fname <- file.path(out_dir, sprintf("model_comparison_%s_block%s.png", score_col, block_num))
  grDevices::png(fname, width = 1400, height = 500, res = 130)
  par(mfrow = c(1, length(metrics)), mar = c(3, 4, 3, 1), oma = c(0, 0, 3, 0))
  
  for (m in metrics) {
    lin_vals  <- as.numeric(sub[sub$model == "ssm_x=Ax+B", m])
    sqrt_vals <- as.numeric(sub[sub$model == "sqrt_trend", m])
    all_vals  <- c(lin_vals, sqrt_vals); all_vals <- all_vals[is.finite(all_vals)]
    if (length(all_vals) == 0) { plot.new(); next }
    
    yr  <- range(all_vals, na.rm = TRUE)
    pad <- diff(yr) * 0.1; if (!is.finite(pad) || pad == 0) pad <- 1
    yr  <- yr + c(-pad, pad)
    
    plot(NA, xlim = c(0.5, 2.5), ylim = yr, xaxt = "n", xlab = "",
         ylab = metric_labels[[m]], main = metric_labels[[m]], cex.main = 0.85)
    axis(1, at = c(1, 2), labels = c("SSM", "\u221a-model"))
    grid(nx = NA, ny = NULL, col = "gray90")
    draw_violin(lin_vals,  at = 1, col = "#4C8CC5")
    draw_violin(sqrt_vals, at = 2, col = "#F1BC18")
  }
  
  mtext(sprintf("Model comparison \u2014 %s (Bloc %s)", label, block_num),
        outer = TRUE, cex = 1.1, font = 2)
  grDevices::dev.off()
  cat(sprintf("    Sauvegard\u00e9 : %s\n", fname))
  fname
}

# ---------------------------------------------------------------------
# 4sexies. TABLEAU RECAPITULATIF (mean/median) — BLOC 3 UNIQUEMENT
# ---------------------------------------------------------------------

summarize_model_comparison_block3 <- function(cmp_df, block_num = 3) {
  sub <- cmp_df[cmp_df$block == block_num, ]
  if (nrow(sub) == 0) return(data.frame())
  
  metrics <- c("loglik", "aic", "bic", "rmse", "cor")
  rows <- list()
  
  for (score_col in sort(unique(sub$score))) {
    for (model_name in sort(unique(sub$model))) {
      part <- sub[sub$score == score_col & sub$model == model_name, ]
      row  <- list(score = score_col, model = model_name, block = block_num,
                   n = sum(stats::complete.cases(part[, metrics, drop = FALSE])))
      for (m in metrics) {
        vals <- as.numeric(part[[m]])
        vals <- vals[is.finite(vals)]
        row[[paste0(m, "_mean")]]   <- if (length(vals) > 0) mean(vals)             else NA_real_
        row[[paste0(m, "_median")]] <- if (length(vals) > 0) stats::median(vals)    else NA_real_
        row[[paste0(m, "_sd")]]     <- if (length(vals) > 1) stats::sd(vals)        else NA_real_
      }
      rows[[length(rows) + 1]] <- as.data.frame(row, stringsAsFactors = FALSE)
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------
# 5. FIGURE D'EVOLUTION (un score x une séquence) + modèles état-espace
# ---------------------------------------------------------------------

plot_evolution_for_score_sequence <- function(df, score_col, sequence, out_dir) {
  sub <- df[!is.na(df$sequence) & df$sequence == sequence, ]
  sub <- sub[!is.na(sub[[score_col]]), ]
  if (nrow(sub) == 0) return(NULL)
  
  meta <- SCORE_META[[score_col]]
  if (is.null(meta)) meta <- list(label = score_col, ylabel = score_col)
  seq_col <- if (sequence %in% names(SEQ_COLORS)) SEQ_COLORS[[sequence]] else "gray30"
  
  participants <- sort(unique(sub$participant))
  n_part       <- length(participants)
  part_colors  <- grDevices::hcl.colors(max(n_part, 3), palette = "Dark 3")[seq_len(n_part)]
  names(part_colors) <- participants
  
  boundaries <- get_block_boundaries(sub)
  yr <- range(sub[[score_col]], na.rm = TRUE)
  pad <- diff(yr) * 0.05
  if (!is.finite(pad) || pad == 0) pad <- 1
  yr <- c(yr[1] - pad, yr[2] + pad)
  
  # ---- Moyenne du groupe par position chronologique --------------
  group_mean <- stats::aggregate(sub[[score_col]], by = list(trial_global = sub$trial_global),
                                 FUN = mean, na.rm = TRUE)
  names(group_mean)[2] <- "score"
  
  gm_block <- stats::aggregate(block ~ trial_global, data = sub,
                               FUN = function(s) as.numeric(names(sort(table(s), decreasing = TRUE))[1]))
  group_mean <- merge(group_mean, gm_block, by = "trial_global")
  group_mean <- group_mean[order(group_mean$trial_global), ]
  
  # ---- Ajustement des modèles état-espace (Kalman) ----------------
  fit_group_linear <- fit_kf_linear_by_block(group_mean$score, group_mean$block)
  fit_group_sqrt   <- fit_kf_sqrt_by_block(group_mean$score, group_mean$block)
  
  fits_part_linear <- list()
  fits_part_sqrt   <- list()
  part_data        <- list()
  for (p in participants) {
    pdata <- sub[sub$participant == p, ]
    pdata <- pdata[order(pdata$trial_global), ]
    part_data[[p]]        <- pdata
    fits_part_linear[[p]] <- fit_kf_linear_by_block(pdata[[score_col]], pdata$block)
    fits_part_sqrt[[p]]   <- fit_kf_sqrt_by_block(pdata[[score_col]], pdata$block)
  }
  
  # ---- Table récapitulative des paramètres ajustés ----------------
  params_rows <- list()
  add_row <- function(subject, model, block = NA, A = NA, mu = NA, Q = NA, Qa = NA, Qb = NA, R = NA) {
    data.frame(score = score_col, sequence = sequence, subject = subject, model = model,
               block = block, A = A, mu = mu, Q = Q, Qa = Qa, Qb = Qb, R = R, stringsAsFactors = FALSE)
  }
  
  append_block_rows <- function(subject, model, fit_res) {
    if (is.null(fit_res$params)) return(NULL)
    lapply(seq_len(nrow(fit_res$params)), function(i) {
      rr <- fit_res$params[i, ]
      if (model == "ssm_x=Ax+B") add_row(subject, model, rr$block, A = rr$A, mu = rr$mu, Q = rr$Q, R = rr$R)
      else                          add_row(subject, model, rr$block, Qa = rr$Qa, Qb = rr$Qb, R = rr$R)
    })
  }
  
  params_rows <- c(params_rows,
                   append_block_rows("GROUP", "ssm_x=Ax+B", fit_group_linear),
                   append_block_rows("GROUP", "sqrt_trend",    fit_group_sqrt))
  for (p in participants) {
    params_rows <- c(params_rows,
                     append_block_rows(p, "ssm_x=Ax+B", fits_part_linear[[p]]),
                     append_block_rows(p, "sqrt_trend",    fits_part_sqrt[[p]]))
  }
  params_df <- if (length(params_rows) > 0) do.call(rbind, params_rows) else NULL
  
  files <- c()
  xr <- range(sub$trial_global, na.rm = TRUE)
  
  draw_common <- function(main_title) {
    plot(NA, xlim = xr, ylim = yr, xlab = "Trial",
         ylab = meta$ylabel, main = main_title, cex.main = 0.95, font.main = 2)
    grid(col = "gray90", lty = 1)
    for (bd in boundaries) abline(v = bd, col = "dimgray", lty = 2, lwd = 1.2)
  }
  
  draw_legend <- function(labels, colors, lty, lwd) {
    par(xpd = TRUE)
    xpos <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
    legend(xpos, par("usr")[4], legend = labels, col = colors, lty = lty, lwd = lwd,
           bty = "n", cex = 0.62, xjust = 0, seg.len = 1.6)
    par(xpd = FALSE)
  }
  
  # =========================================================
  # (a) FIGURE PRINCIPALE (comme spat_plot_evolution.py) + overlay
  #     des 2 courbes état-espace du GROUPE
  # =========================================================
  fname_main <- file.path(out_dir, sprintf("evolution_%s_seq%s.png", score_col, sequence))
  grDevices::png(fname_main, width = 1300, height = 650, res = 130)
  par(mar = c(4, 4.5, 3, 13))
  
  draw_common(sprintf("%s  \u2014  S\u00e9quence %s (n=%d)", meta$label, sequence, n_part))
  
  for (p in participants) {
    pdata <- part_data[[p]]
    col_a <- grDevices::adjustcolor(part_colors[[p]], alpha.f = 0.35)
    lines(pdata$trial_global, pdata[[score_col]], col = col_a, lwd = 1.2)
    points(pdata$trial_global, pdata[[score_col]], col = col_a, pch = 16, cex = 0.4)
  }
  
  lines(group_mean$trial_global, group_mean$score, col = seq_col, lwd = 3.2)
  points(group_mean$trial_global, group_mean$score, col = seq_col, pch = 16, cex = 0.6)
  
  if (any(!is.na(fit_group_linear$x_smooth)))
    lines(group_mean$trial_global, fit_group_linear$x_smooth, col = "black", lwd = 2.3, lty = 2)
  if (any(!is.na(fit_group_sqrt$x_smooth)))
    lines(group_mean$trial_global, fit_group_sqrt$x_smooth, col = "black", lwd = 2.3, lty = 3)
  
  leg_labels <- c(participants, "Moyenne du groupe")
  leg_colors <- c(unname(part_colors), seq_col)
  leg_lty    <- c(rep(1, n_part), 1)
  leg_lwd    <- c(rep(1.2, n_part), 3.2)
  if (any(!is.na(fit_group_linear$x_smooth))) {
    leg_labels <- c(leg_labels, "\u00c9tat-espace lin\u00e9aire (groupe)")
    leg_colors <- c(leg_colors, "black"); leg_lty <- c(leg_lty, 2); leg_lwd <- c(leg_lwd, 2.3)
  }
  if (any(!is.na(fit_group_sqrt$x_smooth))) {
    leg_labels <- c(leg_labels, "Tendance \u221atrial (groupe, Kalman)")
    leg_colors <- c(leg_colors, "black"); leg_lty <- c(leg_lty, 3); leg_lwd <- c(leg_lwd, 2.3)
  }
  draw_legend(leg_labels, leg_colors, leg_lty, leg_lwd)
  
  grDevices::dev.off()
  cat(sprintf("    Sauvegard\u00e9 : %s\n", fname_main))
  files <- c(files, fname_main)
  
  # =========================================================
  # (b) FIGURE DEDIEE — modèle état-espace linéaire x[n+1]=A*x[n]+B
  # =========================================================
  fname_lin <- file.path(out_dir, sprintf("state_space_linear_%s_seq%s.png", score_col, sequence))
  grDevices::png(fname_lin, width = 1300, height = 650, res = 130)
  par(mar = c(4, 4.5, 3, 13))
  
  # dans (b) FIGURE DEDIEE — section plot_evolution_for_score_sequence
  draw_common(sprintf("State-space model  x[n+1] = \u03bc + A\u00b7(x[n]-\u03bc)  \u2014  %s (Seq. %s, n=%d)",
                      meta$label, sequence, n_part))
  
  for (p in participants) {
    fl <- fits_part_linear[[p]]
    if (is.null(fl)) next
    lines(part_data[[p]]$trial_global, fl$x_smooth, col = part_colors[[p]], lwd = 1.5)
  }
  if (any(!is.na(fit_group_linear$x_smooth)))
    lines(group_mean$trial_global, fit_group_linear$x_smooth, col = "black", lwd = 3)
  
  leg_labels <- c(participants, "Groupe (liss\u00e9, Kalman)")
  leg_colors <- c(unname(part_colors), "black")
  draw_legend(leg_labels, leg_colors, rep(1, length(leg_labels)), c(rep(1.5, n_part), 3))
  
  grDevices::dev.off()
  cat(sprintf("    Sauvegard\u00e9 : %s\n", fname_lin))
  files <- c(files, fname_lin)
  
  # =========================================================
  # (c) FIGURE DEDIEE — modèle "racine carrée" y[n] ~= A*sqrt(n)+B
  #     (coefficients A_n, B_n estimés par filtre de Kalman)
  # =========================================================
  fname_sqrt <- file.path(out_dir, sprintf("state_space_sqrt_%s_seq%s.png", score_col, sequence))
  grDevices::png(fname_sqrt, width = 1300, height = 650, res = 130)
  par(mar = c(4, 4.5, 3, 13))
  
  draw_common(sprintf("Square root model  y[n] \u2248 A\u00b7\u221an + B (Kalman)  \u2014  %s (Seq. %s, n=%d)",
                      meta$label, sequence, n_part))
  
  for (p in participants) {
    fs <- fits_part_sqrt[[p]]
    if (is.null(fs)) next
    lines(part_data[[p]]$trial_global, fs$x_smooth, col = part_colors[[p]], lwd = 1.5)
  }
  if (any(!is.na(fit_group_sqrt$x_smooth)))
    lines(group_mean$trial_global, fit_group_sqrt$x_smooth, col = "black", lwd = 3)
  
  leg_labels <- c(participants, "Group (smoothed, Kalman)")
  leg_colors <- c(unname(part_colors), "black")
  draw_legend(leg_labels, leg_colors, rep(1, length(leg_labels)), c(rep(1.5, n_part), 3))
  
  grDevices::dev.off()
  cat(sprintf("    Sauvegard\u00e9 : %s\n", fname_sqrt))
  files <- c(files, fname_sqrt)
  
  list(files = files, params = params_df)
}

# ---------------------------------------------------------------------
# 5bis. COMPARAISON A vs B — MOYENNE DE GROUPE + MODELES + IC (un score)
#       (seuls les 2 modèles sont représentés ; la moyenne brute du
#       groupe n'est plus tracée — cf. section précédente)
# ---------------------------------------------------------------------

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
  
  # --- Bornes Y : basées uniquement sur les modèles (lin + sqrt), plus leurs IC ---
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
  
  # Résumé "A (n=12) vs B (n=10)" pour le titre, générique sur les séquences trouvées.
  n_summary <- paste(
    sprintf("%s (n=%d)", sequences, vapply(sequences, function(s) npart_list[[s]], numeric(1))),
    collapse = " vs "
  )
  
  fname <- file.path(out_dir, sprintf("evolution_comparison_%s_AvsB.pdf", score_col))
  grDevices::pdf(fname, width = 17, height = 7, bg="transparent")
  # Marges agrandies (surtout à gauche) pour que les titres d'axes
  # (xlab/ylab) ne soient jamais coupés.
  par(mar = c(4, 4.5, 3, 2), mgp = c(2.4, 0.7, 0), bg = "#F8F9FA")
  
  plot(NA, xlim = xr, ylim = yr, xlab = "Trial", ylab = "Harmonic score normalized",
       main = sprintf("Modeled overview"),
       cex.main = 2.5, cex.lab = 1.5, cex.axis = 1.2, font.main = 2)
  grid(col = "gray90", lty = 1)
  for (seq_ in sequences) for (bd in boundaries_list[[seq_]]) abline(v = bd, col = "dimgray", lty = 2, lwd = 0.8)
  
  # --- Numéros de blocs, affichés en haut du graphique ---
  # Basés sur la structure de blocs de la première séquence disponible
  # (en pratique identique pour toutes les séquences d'un même groupe).
  ref_seq <- sequences[1]
  gm_ref  <- gm_list[[ref_seq]]
  if (!is.null(gm_ref) && "block" %in% names(gm_ref)) {
    blocks_ref    <- sort(unique(stats::na.omit(gm_ref$block)))
    y_label_block <- yr[2] - diff(yr) * 0.03
    for (b in blocks_ref) {
      idx <- gm_ref$block == b & !is.na(gm_ref$block)
      if (!any(idx)) next
      xmid <- mean(range(gm_ref$trial_global[idx], na.rm = TRUE))
      text(xmid, y_label_block, labels = sprintf("Bloc %d", b),
           col = "gray30", cex = 0.85, font = 2)
    }
  }
  
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
                    sprintf("State space model"),
                    sprintf("Tendency \u221a model"))
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
    cex = 1,
    bty = "o",
    bg = "#F8F9FA",
    inset = 0.02
  )
  
  grDevices::dev.off()
  cat(sprintf("    Sauvegard\u00e9 : %s\n", fname))
  fname
}
# ---------------------------------------------------------------------
# 6. PIPELINE PRINCIPAL
# ---------------------------------------------------------------------

main <- function(scores_folder, out_dir = "./output_plots") {
  cat(strrep("=", 65), "\n")
  cat("SPAT \u2014 VISUALISATION DES COURBES D'\u00c9VOLUTION (R)\n")
  cat(strrep("=", 65), "\n")
  
  cat(sprintf("\n[1/2] Chargement des scores depuis : %s\n", scores_folder))
  all_df <- load_participant_scores(scores_folder)
  if (nrow(all_df) == 0) return(invisible(NULL))
  
  cat("\n  Blocs pr\u00e9sents par participant x s\u00e9quence :\n")
  cov <- stats::aggregate(block ~ participant + sequence, data = all_df,
                          FUN = function(s) paste(sort(unique(s)), collapse = ","))
  print(cov)
  
  score_cols <- intersect(names(SCORE_META), colnames(all_df))
  sequences  <- sort(unique(stats::na.omit(all_df$sequence)))
  cat(sprintf("  Scores trouv\u00e9s     : %s\n", paste(score_cols, collapse = ", ")))
  cat(sprintf("  S\u00e9quences trouv\u00e9es : %s\n", paste(sequences, collapse = ", ")))
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("\n  Calcul des learning rates + violons Bloc 3 (A vs B)...\n")
  lr_df <- compute_learning_rates_all(all_df, score_cols)
  if (nrow(lr_df) > 0) {
    lr_path <- file.path(out_dir, "learning_rates_all_blocks.csv")
    utils::write.csv(lr_df, lr_path, row.names = FALSE)
    cat(sprintf("  Learning rates \u00e9crits dans : %s\n", lr_path))
    
    cat("\n  Calcul des ratios Learning Rate A/B (toutes d\u00e9finitions)...\n")
    lr_ratio_df <- compute_learning_rate_ratios(lr_df)
    if (nrow(lr_ratio_df) > 0) {
      lr_ratio_path <- file.path(out_dir, "learning_rate_ratios_AoverB.csv")
      utils::write.csv(lr_ratio_df, lr_ratio_path, row.names = FALSE)
      cat(sprintf("  Ratios Learning Rate A/B \u00e9crits dans : %s\n", lr_ratio_path))
    }
    
    for (score_col in score_cols) {
      plot_violin_learning_rates(lr_df, score_col, block_num = 3, out_dir)
    }
  }
  
  for (score_col in score_cols) {
    plot_violin_learning_rates(lr_df, score_col, block_num = 3, out_dir)
  }
  
  cat("\n  Comparaison d'ajustement des mod\u00e8les (AIC/BIC/loglik/RMSE/corr)...\n")
  model_cmp_df <- compute_model_comparison_all(all_df, score_cols)
  if (nrow(model_cmp_df) > 0) {
    cmp_path <- file.path(out_dir, "model_comparison_linear_vs_sqrt.csv")
    utils::write.csv(model_cmp_df, cmp_path, row.names = FALSE)
    cat(sprintf("  Comparaison des mod\u00e8les \u00e9crite dans : %s\n", cmp_path))
    
    for (score_col in score_cols) {
      plot_model_comparison_metrics(model_cmp_df, score_col, out_dir)
    }
    
    # --- NOUVEAU : tableau récapitulatif + violon, Bloc 3 uniquement ---
    cat("\n  Tableau r\u00e9capitulatif (mean/median) \u2014 Bloc 3 uniquement...\n")
    summary_block3 <- summarize_model_comparison_block3(model_cmp_df, block_num = 3)
    if (nrow(summary_block3) > 0) {
      summary_path <- file.path(out_dir, "model_comparison_summary_block3.csv")
      utils::write.csv(summary_block3, summary_path, row.names = FALSE)
      cat(sprintf("  Tableau r\u00e9capitulatif \u00e9crit dans : %s\n", summary_path))
      print(summary_block3)
    }
    
    for (score_col in score_cols) {
      plot_model_comparison_metrics_block3(model_cmp_df, score_col, out_dir, block_num = 3)
    }
  }
  
  cat(sprintf("\n[2/2] G\u00e9n\u00e9ration des figures dans : %s\n", out_dir))
  written     <- c()
  params_list <- list()
  
  for (score_col in score_cols) {
    for (sequence in sequences) {
      res <- plot_evolution_for_score_sequence(all_df, score_col, sequence, out_dir)
      if (!is.null(res)) {
        written <- c(written, res$files)
        if (!is.null(res$params)) params_list[[length(params_list) + 1]] <- res$params
      }
    }
  }
  
  cat("\n  G\u00e9n\u00e9ration des figures de comparaison A vs B (moyenne + mod\u00e8les + IC)...\n")
  for (score_col in score_cols) {
    cf <- plot_evolution_group_comparison(all_df, score_col, out_dir)
    if (!is.null(cf)) written <- c(written, cf)
  }
  
  if (length(params_list) > 0) {
    params_all <- do.call(rbind, params_list)
    params_path <- file.path(out_dir, "state_space_params.csv")
    utils::write.csv(params_all, params_path, row.names = FALSE)
    cat(sprintf("  Param\u00e8tres des mod\u00e8les \u00e9tat-espace \u00e9crits dans : %s\n", params_path))
  }
  
  cat("\n", strrep("=", 65), "\n", sep = "")
  cat(sprintf("TERMIN\u00c9. %d figure(s) g\u00e9n\u00e9r\u00e9e(s) dans : %s/\n", length(written), out_dir))
  cat(strrep("=", 65), "\n")
  invisible(written)
}

# ---------------------------------------------------------------------
# EXECUTION EN LIGNE DE COMMANDE
# ---------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  main(args[1], args[2])
} else {
  main("C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_scores_harmonic", "C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_evolution_plot/all4")
}