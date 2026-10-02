#!/usr/bin/env Rscript
# =====================================================================
# SPAT — Calcul du score harmonique (R)
# =====================================================================
# Reproduit la logique de score_harmonic de spat_export_scores.py, y
# compris le REGROUPEMENT DES FICHIERS PAR PARTICIPANT :
#
#   1. Reconstruction du contour lisse de la forme fixe (courbes de
#      Bezier quadratiques, identique à la boucle quadraticCurveTo du
#      script JS de l'expérience)
#   2. Chargement des fichiers SPAT_<participant>_sequence_<seq>_<ts>.xlsx
#      d'un dossier, regroupement par (participant, sequence), tri par
#      timestamp croissant pour assigner le numéro de bloc (comme
#      load_files() en Python), puis concaténation des blocs
#   3. Dérotation de chaque essai (angle inféré à partir de la position
#      du tout premier point du stylet, comparée au point de départ
#      canonique)
#   4. Pour chaque essai :
#        - dtw_dist : distance DTW (sous-échantillonnée à <=200 points)
#                     entre la trajectoire et le segment cible idéal
#        - mean_dev : distance moyenne au plus proche voisin du contour
#                     complet de la forme
#        - time_spent : durée de l'essai (dernière valeur de Time)
#   5. Pour mean_dev, time_spent et score_harmonic (= 1 / sqrt(mean_dev
#      * time_spent)), TROIS versions sont calculées :
#        - _raw               : valeur brute
#        - _norm_participant  : mise à l'échelle linéaire sur [0, 100]
#                                par rapport au min/max observés parmi
#                                TOUS les essais DU PARTICIPANT (toutes
#                                séquences/blocs confondus), comme
#                                add_scores() en Python
#        - _norm_global       : mise à l'échelle linéaire sur [0, 100]
#                                par rapport au min/max observés parmi
#                                TOUS les essais de TOUS LES
#                                PARTICIPANTS
#   6. UN FICHIER EXCEL EST ÉCRIT PAR PARTICIPANT (toutes séquences et
#      tous blocs confondus), comme export_participant_excel() en
#      Python — au lieu d'un fichier par fichier d'entrée.
#
# ENTREE
# ------
# Un dossier contenant des fichiers nommés :
#   SPAT_<participant>_sequence_<seq>_<timestamp>.xlsx
# Chaque fichier doit contenir une feuille avec les colonnes (au
# minimum, insensible à la casse) : Trial, Time, PenX, PenY. La feuille
# "PenAndMouseData" est utilisée si présente, sinon la première feuille
# du classeur.
#
# SORTIE
# ------
# Un fichier Excel par participant : SPAT_<participant>_scores.xlsx,
# feuille "scores", une ligne par essai (toutes séquences/blocs
# confondus), colonnes :
#   block | trial | sequence | dtw_dist |
#   time_spent_raw      | time_spent_norm_participant      | time_spent_norm_global |
#   mean_dev_raw         | mean_dev_norm_participant         | mean_dev_norm_global |
#   score_harmonic_raw   | score_harmonic_norm_participant   | score_harmonic_norm_global
#
# USAGE
# -----
#   Rscript spat_score_harmonic.R input_folder output_dir
#
# Packages requis :
#   install.packages(c("readxl", "writexl", "RANN"))
# =====================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(RANN)
  library(tools)
})

# ---------------------------------------------------------------------
# 1. FORME FIXE (identique à FIXED_MIDDLE_POINTS du script Python)
# ---------------------------------------------------------------------

FIXED_MIDDLE_POINTS <- matrix(c(
  395.091, 270.000,
  378.397, 293.040,
  369.752, 314.412,
  357.767, 333.766,
  367.450, 378.229,
  324.543, 364.472,
  313.689, 404.461,
  281.808, 382.348,
  258.742, 377.109,
  230.016, 393.057,
  200.614, 390.179,
  190.547, 358.242,
  154.773, 353.717,
  169.555, 314.721,
  134.792, 298.739,
  157.273, 270.000,
  155.101, 245.577,
  149.477, 216.340,
  154.720, 186.244,
  187.689, 178.585,
  215.982, 176.437,
  225.221, 132.185,
  257.355, 149.687,
  284.837, 128.835,
  314.638, 132.618,
  323.966, 176.529,
  367.462, 161.758,
  372.795, 195.315,
  371.717, 224.713,
  402.549, 241.826
), ncol = 2, byrow = TRUE)

CANVAS_CENTER             <- c(270.0, 270.0)
FIXED_START_XY_CANONICAL  <- c(195.09, 186.81)
FIXED_TARGET_XY_CANONICAL <- c(165.51, 223.48)

# ---------------------------------------------------------------------
# 2. RECONSTRUCTION DU CONTOUR (boucle quadraticCurveTo)
# ---------------------------------------------------------------------

reconstruct_outline_middle <- function(raw = FIXED_MIDDLE_POINTS,
                                       num_points_per_segment = 100) {
  n     <- nrow(raw)
  pts   <- vector("list", n * num_points_per_segment)
  k     <- 1
  t_seq <- (0:(num_points_per_segment - 1)) / num_points_per_segment
  
  for (i in 0:(n - 1)) {
    p_prev <- raw[((i - 1) %% n) + 1, ]
    p_curr <- raw[i + 1, ]
    p_next <- raw[((i + 1) %% n) + 1, ]
    
    start_pt <- (p_prev + p_curr) / 2
    ctrl_pt  <- p_curr
    end_pt   <- (p_curr + p_next) / 2
    
    for (tv in t_seq) {
      pts[[k]] <- (1 - tv)^2 * start_pt + 2 * (1 - tv) * tv * ctrl_pt + tv^2 * end_pt
      k <- k + 1
    }
  }
  do.call(rbind, pts)
}

# ---------------------------------------------------------------------
# 3. STRUCTURE DU CONTOUR (longueur d'arc cumulée) + SEGMENT CIBLE
# ---------------------------------------------------------------------

precompute_outline_structure <- function(shape_outline) {
  diffs       <- diff(shape_outline)
  seg_lengths <- sqrt(rowSums(diffs^2))
  arc_lengths <- c(0, cumsum(seg_lengths))
  closing     <- sqrt(sum((shape_outline[1, ] - shape_outline[nrow(shape_outline), ])^2))
  total       <- arc_lengths[length(arc_lengths)] + closing
  list(arc_lengths = arc_lengths, total = total)
}

compute_shape_path_params <- function(shape_outline, arc_lengths, total_length) {
  idx_start <- RANN::nn2(shape_outline, matrix(FIXED_START_XY_CANONICAL, nrow = 1), k = 1)$nn.idx[1]
  idx_end   <- RANN::nn2(shape_outline, matrix(FIXED_TARGET_XY_CANONICAL, nrow = 1), k = 1)$nn.idx[1]
  
  s_start <- arc_lengths[idx_start]
  s_end   <- arc_lengths[idx_end]
  
  path_length <- (s_end - s_start) %% total_length
  if (path_length < 1e-6) path_length <- total_length
  
  list(s_start = s_start, path_length = path_length)
}

get_target_segment <- function(shape_outline, arc_lengths, total_length,
                               s_start, path_length, num_points = 100) {
  s_vals <- seq(0, path_length, length.out = num_points)
  M      <- nrow(shape_outline)
  pts    <- matrix(NA_real_, nrow = num_points, ncol = 2)
  
  for (k in seq_along(s_vals)) {
    actual_s <- (s_start + s_vals[k]) %% total_length
    
    # Equivalent de np.searchsorted(arc_lengths, actual_s, side='right') - 1
    idx0 <- findInterval(actual_s, arc_lengths) - 1
    idx0 <- min(max(idx0, 0), M - 2)
    
    seg <- arc_lengths[idx0 + 2] - arc_lengths[idx0 + 1]
    if (seg < 1e-10) {
      pts[k, ] <- shape_outline[idx0 + 1, ]
    } else {
      frac <- (actual_s - arc_lengths[idx0 + 1]) / seg
      pts[k, ] <- shape_outline[idx0 + 1, ] + frac * (shape_outline[idx0 + 2, ] - shape_outline[idx0 + 1, ])
    }
  }
  pts
}

# ---------------------------------------------------------------------
# 4. DEROTATION PAR ESSAI (utilise la colonne trial_id, unique au sein
#    du data.frame passé en argument)
# ---------------------------------------------------------------------

derotate_trials <- function(df) {
  cx <- CANVAS_CENTER[1]; cy <- CANVAS_CENTER[2]
  theta0 <- atan2(FIXED_START_XY_CANONICAL[2] - cy, FIXED_START_XY_CANONICAL[1] - cx)
  
  df$PenX_derot <- NA_real_
  df$PenY_derot <- NA_real_
  
  for (tr in unique(df$trial_id)) {
    idx   <- which(df$trial_id == tr)
    first <- df[idx[1], ]
    alpha <- atan2(first$PenY - cy, first$PenX - cx) - theta0
    
    cos_a <- cos(-alpha); sin_a <- sin(-alpha)
    dx <- df$PenX[idx] - cx
    dy <- df$PenY[idx] - cy
    df$PenX_derot[idx] <- cx + dx * cos_a - dy * sin_a
    df$PenY_derot[idx] <- cy + dx * sin_a + dy * cos_a
  }
  df
}

# ---------------------------------------------------------------------
# 5. DTW (distance normalisée par la somme des longueurs des chemins)
# ---------------------------------------------------------------------

dtw_distance <- function(traj, target) {
  n <- nrow(traj); m <- nrow(target)
  if (n == 0 || m == 0) return(0)
  
  M <- matrix(Inf, nrow = n + 1, ncol = m + 1)
  M[1, 1] <- 0
  for (i in 1:n) {
    for (j in 1:m) {
      cost <- sqrt(sum((traj[i, ] - target[j, ])^2))
      M[i + 1, j + 1] <- cost + min(M[i, j + 1], M[i + 1, j], M[i, j])
    }
  }
  M[n + 1, m + 1] / (n + m)
}

# ---------------------------------------------------------------------
# 6. MESURES BRUTES PAR ESSAI (dtw_dist, mean_dev, time_spent)
# ---------------------------------------------------------------------

compute_raw_trial_measurements <- function(df, shape_outline, target_segment) {
  rows <- list()
  
  for (tr in unique(df$trial_id)) {
    tdf <- df[df$trial_id == tr, ]
    if (nrow(tdf) == 0) next
    
    t_max <- max(tdf$Time, na.rm = TRUE)
    pts   <- as.matrix(tdf[, c("PenX_derot", "PenY_derot")])
    
    # DTW sous-échantillonné à <= 200 points (comme en Python)
    if (nrow(pts) > 200) {
      idx <- floor(seq(0, nrow(pts) - 1, length.out = 200)) + 1
      sub <- pts[idx, , drop = FALSE]
    } else {
      sub <- pts
    }
    dtw <- dtw_distance(sub, target_segment)
    
    # Distance moyenne au plus proche voisin du contour complet
    nn       <- RANN::nn2(shape_outline, pts, k = 1)
    mean_dev <- mean(nn$nn.dists[, 1])
    
    rows[[length(rows) + 1]] <- data.frame(
      trial      = tr,
      time_spent = t_max,
      dtw_dist   = dtw,
      mean_dev   = mean_dev
    )
  }
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------
# 7. NORMALISATION ET SCORE HARMONIQUE
# ---------------------------------------------------------------------
#
# normalize_minmax() met à l'échelle linéaire un vecteur sur [0, 100]
# par rapport à son propre min/max. Utilisée à la fois :
#   - au sein d'un participant (min/max calculés sur ses essais)
#   - sur l'ensemble des participants (min/max calculés sur tous les
#     essais de tous les participants)
#
# compute_score_harmonic_raw() reproduit score_harmonic =
# 1 / sqrt(mean_dev * time_spent), sur les valeurs BRUTES.
# ---------------------------------------------------------------------

normalize_minmax <- function(x) {
  lo <- min(x, na.rm = TRUE)
  hi <- max(x, na.rm = TRUE)
  if (hi - lo < 1e-12) {
    return(rep(50.0, length(x)))
  }
  100 * (x - lo) / (hi - lo)
}

compute_score_harmonic_raw <- function(mean_dev, time_spent) {
  1 / sqrt(pmax(mean_dev, 1e-6) * pmax(time_spent, 1e-6))
}

# Ajoute, sur un data.frame contenant TOUS les essais de TOUS les
# participants (colonne `participant` requise), les colonnes :
#   time_spent_raw / _norm_participant / _norm_global
#   mean_dev_raw   / _norm_participant / _norm_global
#   score_harmonic_raw / _norm_participant / _norm_global
add_all_score_columns <- function(df) {
  df$time_spent_raw    <- df$time_spent
  df$mean_dev_raw      <- df$mean_dev
  df$score_harmonic_raw <- compute_score_harmonic_raw(df$mean_dev_raw, df$time_spent_raw)
  
  # --- Normalisation sur l'ensemble des participants (global) ---
  df$time_spent_norm_global     <- normalize_minmax(df$time_spent_raw)
  df$mean_dev_norm_global       <- normalize_minmax(df$mean_dev_raw)
  df$score_harmonic_norm_global <- normalize_minmax(df$score_harmonic_raw)
  
  cat(sprintf("  [GLOBAL] time_spent      — lo=%.6f  hi=%.6f\n", min(df$time_spent_raw), max(df$time_spent_raw)))
  cat(sprintf("  [GLOBAL] mean_dev        — lo=%.6f  hi=%.6f\n", min(df$mean_dev_raw), max(df$mean_dev_raw)))
  cat(sprintf("  [GLOBAL] score_harmonic  — lo=%.6f  hi=%.6f\n", min(df$score_harmonic_raw), max(df$score_harmonic_raw)))
  
  # --- Normalisation au sein de chaque participant ---
  df$time_spent_norm_participant     <- ave(df$time_spent_raw, df$participant, FUN = normalize_minmax)
  df$mean_dev_norm_participant       <- ave(df$mean_dev_raw, df$participant, FUN = normalize_minmax)
  df$score_harmonic_norm_participant <- ave(df$score_harmonic_raw, df$participant, FUN = normalize_minmax)
  
  for (p in unique(df$participant)) {
    sub <- df[df$participant == p, ]
    cat(sprintf("  [%s] time_spent      — lo=%.6f  hi=%.6f\n", p, min(sub$time_spent_raw), max(sub$time_spent_raw)))
    cat(sprintf("  [%s] mean_dev        — lo=%.6f  hi=%.6f\n", p, min(sub$mean_dev_raw), max(sub$mean_dev_raw)))
    cat(sprintf("  [%s] score_harmonic  — lo=%.6f  hi=%.6f\n", p, min(sub$score_harmonic_raw), max(sub$score_harmonic_raw)))
  }
  
  df
}

# ---------------------------------------------------------------------
# 8. CHARGEMENT ET REGROUPEMENT DES FICHIERS PAR (PARTICIPANT, SEQUENCE)
#    — équivalent de load_files() dans spat_export_scores.py
# ---------------------------------------------------------------------

LOAD_RE <- "^SPAT_(.+?)_sequence_([A-Za-z])_([0-9]+)\\.xlsx$"

read_one_file <- function(path) {
  sheets <- readxl::excel_sheets(path)
  sheet  <- if ("PenAndMouseData" %in% sheets) "PenAndMouseData" else sheets[1]
  df <- as.data.frame(readxl::read_excel(path, sheet = sheet))
  
  # Résolution des colonnes requises, insensible à la casse
  find_col <- function(target) {
    hit <- names(df)[tolower(trimws(names(df))) == target]
    if (length(hit) == 0) {
      stop(sprintf("Colonne '%s' introuvable dans '%s'.", target, basename(path)))
    }
    hit[1]
  }
  names(df)[names(df) == find_col("trial")] <- "Trial"
  names(df)[names(df) == find_col("time")]  <- "Time"
  names(df)[names(df) == find_col("penx")]  <- "PenX"
  names(df)[names(df) == find_col("peny")]  <- "PenY"
  df
}

load_files <- function(folder) {
  if (!dir.exists(folder)) {
    stop(sprintf("Dossier d'entrée introuvable : %s", folder))
  }
  files <- list.files(folder, pattern = "^SPAT_.*\\.xlsx$",
                      full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) {
    stop(sprintf("Aucun fichier SPAT_*.xlsx trouvé dans : %s", folder))
  }
  
  # Analyse des noms de fichiers -> participant / séquence / timestamp
  info <- lapply(files, function(f) {
    fname <- basename(f)
    m <- regmatches(fname, regexec(LOAD_RE, fname, ignore.case = TRUE))[[1]]
    if (length(m) == 0) {
      cat(sprintf("  [SKIP] '%s' — ne correspond pas au format attendu\n", fname))
      return(NULL)
    }
    list(path = f, participant = m[2], sequence = toupper(m[3]), timestamp = m[4])
  })
  info <- Filter(Negate(is.null), info)
  if (length(info) == 0) {
    stop("Aucun fichier ne correspond au format 'SPAT_<participant>_sequence_<seq>_<timestamp>.xlsx'.")
  }
  
  # Regroupement par (participant, sequence)
  keys   <- vapply(info, function(x) paste(x$participant, x$sequence, sep = "|||"), character(1))
  groups <- split(info, keys)
  
  data <- list()
  for (key in names(groups)) {
    entries <- groups[[key]]
    # Tri par timestamp croissant = ordre des blocs (comme en Python)
    entries <- entries[order(vapply(entries, function(x) x$timestamp, character(1)))]
    
    frames <- vector("list", length(entries))
    for (block_num in seq_along(entries)) {
      e  <- entries[[block_num]]
      df <- read_one_file(e$path)
      df$Block <- block_num  # numéro de bloc = rang temporel du fichier
      frames[[block_num]] <- df
      cat(sprintf("  Chargé : %s  seq %s  bloc %d  (%d lignes)  [%s]\n",
                  e$participant, e$sequence, block_num, nrow(df), basename(e$path)))
    }
    
    combined <- do.call(rbind, frames)
    # Trial unique à travers les blocs, pour éviter les collisions
    combined$Trial <- combined$Block * 10000 + combined$Trial
    
    parts <- strsplit(key, "\\|\\|\\|")[[1]]
    data[[key]] <- list(participant = parts[1], sequence = parts[2], df = combined)
    cat(sprintf("  -> %s seq %s : %d lignes totales (%d bloc(s))\n\n",
                parts[1], parts[2], nrow(combined), length(frames)))
  }
  
  data
}

# ---------------------------------------------------------------------
# 9. TRAITEMENT D'UN GROUPE (PARTICIPANT, SEQUENCE) : dérotation +
#    mesures brutes. Les scores (normalisations) ne sont PAS calculés
#    ici : ils sont calculés une fois toutes les séquences/blocs de
#    TOUS les participants rassemblés (cf. main()), pour pouvoir
#    normaliser à la fois par participant et sur l'ensemble.
# ---------------------------------------------------------------------

process_group_measurements <- function(df, shape_outline, target_segment) {
  df$trial_id <- df$Trial  # déjà unique (Block*10000+Trial)
  df <- derotate_trials(df)
  raw_meas <- compute_raw_trial_measurements(df, shape_outline, target_segment)
  
  block_map <- unique(df[, c("trial_id", "Block")])
  raw_meas  <- merge(raw_meas, block_map, by.x = "trial", by.y = "trial_id")
  raw_meas$block <- raw_meas$Block
  raw_meas$Block <- NULL
  raw_meas$trial <- raw_meas$trial %% 10000  # restaure le numéro de trial d'origine
  raw_meas
}

# ---------------------------------------------------------------------
# 10. EXPORT EXCEL — UN FICHIER PAR PARTICIPANT (toutes séquences et
#     blocs confondus), équivalent de export_participant_excel() en
#     Python
# ---------------------------------------------------------------------

export_participant_excel <- function(participant, participant_df, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  df <- participant_df[order(participant_df$sequence, participant_df$trial), ]
  
  cols <- c(
    "block", "trial", "sequence", "dtw_dist",
    "time_spent_raw", "time_spent_norm_participant", "time_spent_norm_global",
    "mean_dev_raw", "mean_dev_norm_participant", "mean_dev_norm_global",
    "score_harmonic_raw", "score_harmonic_norm_participant", "score_harmonic_norm_global"
  )
  df <- df[, cols]
  rownames(df) <- NULL
  
  fname <- file.path(out_dir, sprintf("SPAT_%s_scores.xlsx", participant))
  writexl::write_xlsx(list(scores = df), fname)
  cat(sprintf("  Écrit : %s  (%d lignes)\n", fname, nrow(df)))
  fname
}

# ---------------------------------------------------------------------
# 11. PIPELINE PRINCIPAL
# ---------------------------------------------------------------------

main <- function(input_folder, output_dir) {
  cat("======================================================================\n")
  cat("SPAT — EXPORT DES SCORES PAR PARTICIPANT (R)\n")
  cat("======================================================================\n")
  
  cat(sprintf("\n[1/3] Chargement des données depuis : %s\n", input_folder))
  raw_data <- load_files(input_folder)
  
  cat("\n[2/3] Reconstruction de la forme et du chemin cible idéal…\n")
  shape_outline  <- reconstruct_outline_middle()
  struct         <- precompute_outline_structure(shape_outline)
  path_params    <- compute_shape_path_params(shape_outline, struct$arc_lengths, struct$total)
  target_segment <- get_target_segment(shape_outline, struct$arc_lengths, struct$total,
                                       path_params$s_start, path_params$path_length)
  cat(sprintf("  Points du contour : %d   s_start=%.1f px   path_length=%.1f px\n",
              nrow(shape_outline), path_params$s_start, path_params$path_length))
  
  cat("\n[3/3] Dérotation, mesures brutes, normalisation (participant + globale) et export…\n")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Mesures brutes pour chaque (participant, sequence), rassemblées
  # dans un seul data.frame couvrant TOUS les participants
  all_meas <- list()
  for (key in names(raw_data)) {
    grp      <- raw_data[[key]]
    raw_meas <- process_group_measurements(grp$df, shape_outline, target_segment)
    raw_meas$sequence    <- grp$sequence
    raw_meas$participant <- grp$participant
    all_meas[[length(all_meas) + 1]] <- raw_meas
    cat(sprintf("  %s Seq %s : %d essai(s) mesuré(s)\n", grp$participant, grp$sequence, nrow(raw_meas)))
  }
  all_meas <- do.call(rbind, all_meas)
  
  cat("\n  Calcul des colonnes brutes / normalisées par participant / normalisées globalement…\n")
  all_scored <- add_all_score_columns(all_meas)
  
  written_files <- character(0)
  for (participant in unique(all_scored$participant)) {
    cat(sprintf("\n  --- %s ---\n", participant))
    scored <- all_scored[all_scored$participant == participant, ]
    fpath  <- export_participant_excel(participant, scored, output_dir)
    written_files <- c(written_files, fpath)
  }
  
  cat("\n======================================================================\n")
  cat(sprintf("TERMINÉ. %d fichier(s) Excel écrit(s) dans : %s\n", length(written_files), output_dir))
  cat("======================================================================\n")
  invisible(written_files)
}

# ---------------------------------------------------------------------
# EXECUTION EN LIGNE DE COMMANDE
# ---------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  main(args[1], args[2])
} else {
  main("C:/Users/Alex/Documents/Stage_IRCN/SPAT/data_pen_separated/SPAT",
       "C:/Users/Alex/Documents/Stage_IRCN/SPAT/output_scores_harmonic")
}
