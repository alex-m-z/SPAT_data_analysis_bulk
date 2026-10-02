#!/usr/bin/env Rscript
# =====================================================================
# SPAT — Figure PDF PAR ESSAI, pour UN SEUL PARTICIPANT
# =====================================================================
# Reprend :
#   - la reconstruction du contour lisse de la forme fixe (Bezier
#     quadratiques) de spat_score_harmonic.R
#   - la logique de chargement/regroupement par (participant,
#     séquence), tri par timestamp -> numéro de bloc, et dérotation
#     de spat_score_harmonic.R
#   - l'esthétique / l'organisation "un participant à la fois, sortie
#     en pdf" de spat_plot_participant37.R
#
# Pour CHAQUE essai (trial) du participant demandé, une figure est
# produite (en PDF ET en PNG), superposant :
#   - le contour complet de la forme fixe (gris, fermé)
#   - la trajectoire RÉELLEMENT PRODUITE par le participant pour cet
#     essai, DEROTÉE selon la même logique que spat_score_harmonic.R
#     (angle du tout premier point du stylet aligné sur le point de
#     départ canonique de la forme), pour permettre une comparaison
#     directe avec le contour de la forme
#   - les marqueurs de DÉPART et de FIN placés sur le PREMIER et le
#     DERNIER point RÉELLEMENT tracés par le participant pour cet
#     essai (et non les points canoniques de la forme, qui ne
#     coïncident pas nécessairement avec eux après dérotation)
#
# ENTREE
# ------
# Un dossier contenant les fichiers bruts :
#   SPAT_<participant>_sequence_<seq>_<timestamp>.xlsx
# avec (entre autres) les colonnes :
#   Trial, Block, Time, MouseX, MouseY, MovementX, MovementY,
#   PenX, PenY, Allowed, MouseUp, Error, Score, Source, Angle, Lag
# Seules Trial, Time, PenX, PenY sont utilisées ici (la colonne Block
# du fichier brut, si présente, est ignorée : le numéro de bloc est
# recalculé à partir du rang temporel des fichiers, comme dans
# spat_score_harmonic.R, pour rester cohérent avec lui).
#
# SORTIE
# ------
# Deux fichiers par essai, dans out_dir :
#   trial_seq<seq>_block<bloc>_trial<trial>.pdf
#   trial_seq<seq>_block<bloc>_trial<trial>.png
#
# USAGE
# -----
#   Rscript spat_plot_trials_participant.R input_folder out_dir participant_id
#
# Packages requis :
#   install.packages(c("readxl"))
# =====================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(tools)
  library(svglite)

})

# ---------------------------------------------------------------------
# 1. FORME FIXE (identique à spat_score_harmonic.R)
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
# 2. RECONSTRUCTION DU CONTOUR (boucle quadraticCurveTo, identique)
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
# 3. DEROTATION PAR ESSAI (identique à spat_score_harmonic.R)
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
# 4. CHARGEMENT DES FICHIERS D'UN SEUL PARTICIPANT (regroupement par
#    séquence, tri par timestamp -> bloc), équivalent restreint de
#    load_files() dans spat_score_harmonic.R
# ---------------------------------------------------------------------

LOAD_RE <- "^SPAT_(.+?)_sequence_([A-Za-z])_([0-9]+)\\.xlsx$"

find_col <- function(df, target, fname) {
  hit <- names(df)[tolower(trimws(names(df))) == target]
  if (length(hit) == 0) {
    stop(sprintf("Colonne '%s' introuvable dans '%s'.", target, basename(fname)))
  }
  hit[1]
}

read_one_file <- function(path) {
  sheets <- readxl::excel_sheets(path)
  sheet  <- if ("PenAndMouseData" %in% sheets) "PenAndMouseData" else sheets[1]
  df <- as.data.frame(readxl::read_excel(path, sheet = sheet))

  names(df)[names(df) == find_col(df, "trial", path)] <- "Trial"
  names(df)[names(df) == find_col(df, "time",  path)] <- "Time"
  names(df)[names(df) == find_col(df, "penx",  path)] <- "PenX"
  names(df)[names(df) == find_col(df, "peny",  path)] <- "PenY"
  df
}

load_participant_files <- function(folder, participant_id) {
  if (!dir.exists(folder)) {
    stop(sprintf("Dossier d'entrée introuvable : %s", folder))
  }
  files <- list.files(folder, pattern = "^SPAT_.*\\.xlsx$",
                      full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) {
    stop(sprintf("Aucun fichier SPAT_*.xlsx trouvé dans : %s", folder))
  }

  info <- lapply(files, function(f) {
    fname <- basename(f)
    m <- regmatches(fname, regexec(LOAD_RE, fname, ignore.case = TRUE))[[1]]
    if (length(m) == 0) return(NULL)
    list(path = f, participant = m[2], sequence = toupper(m[3]), timestamp = m[4])
  })
  info <- Filter(Negate(is.null), info)
  info <- Filter(function(x) tolower(x$participant) == tolower(participant_id), info)

  if (length(info) == 0) {
    stop(sprintf("Aucun fichier trouvé pour le participant '%s' dans : %s", participant_id, folder))
  }

  keys   <- vapply(info, function(x) x$sequence, character(1))
  groups <- split(info, keys)

  data <- list()
  for (seq_key in names(groups)) {
    entries <- groups[[seq_key]]
    entries <- entries[order(vapply(entries, function(x) x$timestamp, character(1)))]

    frames <- vector("list", length(entries))
    for (block_num in seq_along(entries)) {
      e  <- entries[[block_num]]
      df <- read_one_file(e$path)
      df$Block <- block_num
      frames[[block_num]] <- df
      cat(sprintf("  Chargé : %s  seq %s  bloc %d  (%d lignes)  [%s]\n",
                  e$participant, e$sequence, block_num, nrow(df), basename(e$path)))
    }

    combined <- do.call(rbind, frames)
    combined$Trial <- combined$Block * 10000 + combined$Trial

    data[[seq_key]] <- list(participant = entries[[1]]$participant,
                            sequence    = seq_key,
                            df          = combined)
    cat(sprintf("  -> seq %s : %d lignes totales (%d bloc(s))\n\n", seq_key, nrow(combined), length(frames)))
  }

  data
}

# ---------------------------------------------------------------------
# 5. FIGURE PAR ESSAI : contour de la forme + trajectoire dérotée
# ---------------------------------------------------------------------

plot_one_trial <- function(shape_outline, trial_df, participant_id, sequence, block, trial, out_dir,
                           save_pdf = TRUE, save_png = TRUE, png_res = 150) {
  traj <- trial_df[order(trial_df$Time), c("PenX_derot", "PenY_derot")]
  
  # Départ / fin = premier et dernier point RÉELLEMENT tracés par le
  # participant pour cet essai (et non les points canoniques de la
  # forme, qui ne coïncident pas forcément avec eux après dérotation).
  start_xy <- c(traj$PenX_derot[1],        traj$PenY_derot[1])
  end_xy   <- c(traj$PenX_derot[nrow(traj)], traj$PenY_derot[nrow(traj)])
  
  all_x <- c(shape_outline[, 1], traj$PenX_derot)
  all_y <- c(shape_outline[, 2], traj$PenY_derot)
  
  xr <- range(all_x, na.rm = TRUE)
  yr <- range(all_y, na.rm = TRUE)
  pad_x <- diff(xr) * 0.06; if (!is.finite(pad_x) || pad_x == 0) pad_x <- 10
  pad_y <- diff(yr) * 0.06; if (!is.finite(pad_y) || pad_y == 0) pad_y <- 10
  xr <- xr + c(-pad_x, pad_x)
  yr <- yr + c(-pad_y, pad_y)
  
  # --- Fonction interne : dessine le graphique sur le périphérique
  #     graphique actuellement ouvert (SVG ou PNG), afin de ne pas
  #     dupliquer le code de tracé entre les deux formats de sortie. ---
  draw_plot <- function() {
    par(mar = c(2, 2, 4, 2), mgp = c(2.3, 0.7, 0), bg = "transparent")
    
    plot(NA, xlim = xr, ylim = rev(yr), xlab = "", ylab = "",
         main = sprintf(""),
         asp = 1, axes = FALSE, cex.main = 1.3, font.main = 2)
    
    # Contour de la forme fixe (fermé)
    polygon(shape_outline[, 1], shape_outline[, 2],
            border = "#082A39", col = NA, lwd = 4)
    
    # Trajectoire réellement produite (dérotée)
    lines(traj$PenX_derot, traj$PenY_derot, col = "#4E8CC3", lwd = 4)
    
    # Départ / fin = premier et dernier point de la trajectoire tracée
    points(start_xy[1], start_xy[2], pch = 17, col = "#6E7F36", cex = 1.4)
    points(end_xy[1],   end_xy[2],   pch = 15, col = "#CC5A71", cex = 1.4)
  }

  base_name <- sprintf("trial_seq%s_block%d_trial%d", sequence, block, trial)
  fnames <- character(0)

  if (isTRUE(save_pdf)) {
    fname_pdf <- file.path(out_dir, paste0(base_name, ".svg"))
    svglite::svglite(fname_pdf, width = 7, height = 7, bg = "transparent")
    draw_plot()
    grDevices::dev.off()
    fnames <- c(fnames, fname_pdf)
  }

  if (isTRUE(save_png)) {
    fname_png <- file.path(out_dir, paste0(base_name, ".png"))
    grDevices::png(fname_png, width = 7, height = 7, units = "in", res = png_res, bg = "transparent")
    draw_plot()
    grDevices::dev.off()
    fnames <- c(fnames, fname_png)
  }

  fnames
}

# ---------------------------------------------------------------------
# 6. PIPELINE PRINCIPAL
# ---------------------------------------------------------------------

TRIALS_TO_PLOT <- data.frame(
  sequence = c("A", "A", "A", "A", "A", "B", "B"),
  block    = c(1,   2,   2,   3,   4,   2,   3),
  trial    = c(18,  1,   14,  1,   2,   1,   1)
)

is_trial_wanted <- function(sequence, block, trial, trials_to_plot) {
  if (is.null(trials_to_plot)) return(TRUE)
  any(toupper(trials_to_plot$sequence) == toupper(sequence) &
        trials_to_plot$block == block &
        trials_to_plot$trial == trial)
}

main <- function(input_folder, out_dir, participant_id, trials_to_plot = TRIALS_TO_PLOT) {
  cat(strrep("=", 70), "\n")
  cat(sprintf("SPAT — FIGURES PAR ESSAI — PARTICIPANT %s\n", participant_id))
  cat(strrep("=", 70), "\n")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  cat(sprintf("\n[1/3] Chargement des fichiers du participant depuis : %s\n", input_folder))
  raw_data <- load_participant_files(input_folder, participant_id)

  cat("\n[2/3] Reconstruction du contour de la forme fixe…\n")
  shape_outline <- reconstruct_outline_middle()
  cat(sprintf("  Points du contour : %d\n", nrow(shape_outline)))

  cat("\n[3/3] Dérotation par essai et génération des figures…\n")
  written <- character(0)
  n_trials <- 0

  for (seq_key in names(raw_data)) {
    grp <- raw_data[[seq_key]]
    df  <- grp$df
    df$trial_id <- df$Trial
    df <- derotate_trials(df)

    block_map <- unique(df[, c("trial_id", "Block")])

    for (tr_id in sort(unique(df$trial_id))) {
      block       <- block_map$Block[block_map$trial_id == tr_id][1]
      trial_local <- tr_id %% 10000
      
      # <-- AJOUT : on ignore tout essai qui n'est pas dans la liste
      if (!is_trial_wanted(grp$sequence, block, trial_local, trials_to_plot)) next
      
      tdf <- df[df$trial_id == tr_id, ]
      if (nrow(tdf) == 0) next
      
      fpaths <- plot_one_trial(shape_outline, tdf, participant_id, grp$sequence,
                               block, trial_local, out_dir)
      written <- c(written, fpaths)
      n_trials <- n_trials + 1
      for (fp in fpaths) cat(sprintf("  Sauvegardé : %s\n", fp))
    }
  }

  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(sprintf("TERMINÉ. %d essai(s) traité(s), %d fichier(s) écrit(s) (PDF + PNG par essai) dans : %s\n",
              n_trials, length(written), out_dir))
  cat(strrep("=", 70), "\n")
  invisible(written)
}

# ---------------------------------------------------------------------
# EXECUTION EN LIGNE DE COMMANDE
# ---------------------------------------------------------------------
#   Rscript spat_plot_trials_participant.R input_folder out_dir participant_id



args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  main(args[1], args[2], args[3])
} else {
  main("C:/Users/Alex/Documents/Stage_IRCN/SPAT/data_pen_separated/SPAT",
       "C:/Users/Alex/Documents/Stage_IRCN/SPAT/plots_par_essai",
       "SPAT001")
}
