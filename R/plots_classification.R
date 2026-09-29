# =============================================================================
# plots_classification.R -- ROC, précision-rappel, calibration, confusion
# (family = "binomial"), une seule implémentation pour toutes les méthodes
# =============================================================================
# Construction manuelle (aucune dépendance à {pROC}).

#' Coordonnées d'une courbe ROC
#' @noRd
.roc_df <- function(prob, actual) {
  ord <- order(-prob)
  a <- actual[ord]
  tibble::tibble(FPR = c(0, cumsum(1 - a) / sum(1 - a)), TPR = c(0, cumsum(a) / sum(a)))
}

#' Points de calibration par tranche de probabilité
#' @noRd
.calib_df <- function(df, bins) {
  df %>%
    dplyr::mutate(bin = cut(Prediction, breaks = seq(0, 1, length.out = bins + 1), include.lowest = TRUE)) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(prob_moyenne = mean(Prediction), freq_observee = mean(Reel), n = dplyr::n(), .groups = "drop") %>%
    tidyr::drop_na(bin)
}

.lab_roc <- function() {
  list(x = "Taux de faux positifs (1 - sp\u00e9cificit\u00e9)", y = "Taux de vrais positifs (sensibilit\u00e9)")
}

#' Courbe ROC + AUC
#'
#' @param fit objet "cv_fit" (family = "binomial")
#' @param mode "pooled" (défaut, CV honnête), "final" (modèle final sur ses
#'   données d'entraînement, optimiste) ou "overlay" (les deux)
#' @export
plot_roc <- function(fit, mode = c("pooled", "final", "overlay")) {
  .stop_if_not_cv_fit(fit)
  .stop_if_not_binomial(fit, "plot_roc")
  mode <- match.arg(mode)
  lbl <- fit$method_label
  fp  <- fit$fold_predictions
  auc_pooled <- .compute_auc(fp$Prediction, fp$Reel)
  base <- function(df, colour) {
    ggplot2::ggplot(df, ggplot2::aes(x = FPR, y = TPR)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
      ggplot2::geom_path(colour = colour, linewidth = 0.9) + ggplot2::coord_equal()
  }
  if (mode == "pooled") {
    return(base(.roc_df(fp$Prediction, fp$Reel), "#D95F02") +
             ggplot2::labs(title = sprintf("Courbe ROC (%s, pr\u00e9dictions pool\u00e9es, CV honn\u00eate)", lbl),
                           subtitle = sprintf("AUC = %.3f \nPointill\u00e9 = classifieur al\u00e9atoire (AUC = 0.5)", auc_pooled),
                           x = .lab_roc()$x, y = .lab_roc()$y) + .theme_cv(13))
  }
  fin <- .final_predictions(fit)
  auc_final <- .compute_auc(fin$pred, fin$Y)
  if (mode == "final") {
    return(base(.roc_df(fin$pred, fin$Y), "#D95F02") +
             ggplot2::labs(title = sprintf("Courbe ROC (%s) -- mod\u00e8le final sur ses donn\u00e9es d'entra\u00eenement", lbl),
                           subtitle = sprintf("\u26a0 AUC = %.3f, OPTIMISTE (entra\u00eenement) \ncomparer \u00e0 l'AUC pool\u00e9e (%.3f) via mode = 'overlay'", auc_final, auc_pooled),
                           x = .lab_roc()$x, y = .lab_roc()$y) +
             .theme_cv(13, 9, "#D95F02"))
  }
  combined <- dplyr::bind_rows(
    .roc_df(fp$Prediction, fp$Reel) %>% dplyr::mutate(Source = sprintf("CV pool\u00e9e (AUC=%.3f)", auc_pooled)),
    .roc_df(fin$pred, fin$Y) %>% dplyr::mutate(Source = sprintf("Mod\u00e8le final, entra\u00eenement (AUC=%.3f)", auc_final))
  )
  ggplot2::ggplot(combined, ggplot2::aes(x = FPR, y = TPR, colour = Source)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_path(linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = c("#2C7FB8", "#D95F02")) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = sprintf("Courbe ROC (%s) : CV honn\u00eate vs mod\u00e8le final (optimiste)", lbl),
                  subtitle = "L'\u00e9cart entre les 2 courbes = biais d'\u00e9valuation intra-\u00e9chantillon",
                  x = .lab_roc()$x, y = .lab_roc()$y, colour = NULL) +
    .theme_cv(13, 8.5) + ggplot2::theme(legend.position = "bottom")
}

#' Courbe précision-rappel (prédictions poolées)
#' @param fit objet "cv_fit" (family = "binomial")
#' @export
plot_precision_recall <- function(fit) {
  .stop_if_not_cv_fit(fit)
  .stop_if_not_binomial(fit, "plot_precision_recall")
  fp <- fit$fold_predictions
  a  <- fp$Reel[order(-fp$Prediction)]
  pr <- tibble::tibble(Recall = cumsum(a) / sum(a), Precision = cumsum(a) / seq_along(a))
  prevalence <- mean(fp$Reel)
  ggplot2::ggplot(pr, ggplot2::aes(x = Recall, y = Precision)) +
    ggplot2::geom_hline(yintercept = prevalence, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_path(colour = "#D95F02", linewidth = 0.9) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) + ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = sprintf("Courbe Pr\u00e9cision-Rappel (%s, pr\u00e9dictions pool\u00e9es)", fit$method_label),
                  subtitle = sprintf("Pointill\u00e9 = classifieur al\u00e9atoire (pr\u00e9valence = %.3f)", prevalence),
                  x = "Rappel (sensibilit\u00e9)", y = "Pr\u00e9cision (VPP)") +
    .theme_cv(13)
}

#' Calibration : probabilité prédite (moyenne par tranche) vs fréquence observée
#'
#' @param fit objet "cv_fit" (family = "binomial")
#' @param bins nombre de tranches (défaut 10)
#' @param mode "pooled" (défaut), "final" ou "overlay"
#' @export
plot_calibration <- function(fit, bins = 10, mode = c("pooled", "final", "overlay")) {
  .stop_if_not_cv_fit(fit)
  .stop_if_not_binomial(fit, "plot_calibration")
  mode <- match.arg(mode)
  lbl  <- fit$method_label
  note <- method_notes(fit, "calibration")
  sous_titre <- paste0("Pointill\u00e9 = calibration parfaite \u00b7 Taille = nb de pr\u00e9dictions dans le bin",
                       if (!is.null(note)) paste0("\n", note) else "")
  base <- function(df, colour) {
    ggplot2::ggplot(df, ggplot2::aes(x = prob_moyenne, y = freq_observee)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
      ggplot2::geom_line(colour = colour) +
      ggplot2::geom_point(ggplot2::aes(size = n), colour = colour, alpha = 0.8)
  }
  common <- list(ggplot2::scale_x_continuous(limits = c(0, 1)), ggplot2::scale_y_continuous(limits = c(0, 1)),
                 ggplot2::coord_equal())
  labs_xy <- list(x = "Probabilit\u00e9 pr\u00e9dite (moyenne du bin)", y = "Fr\u00e9quence observ\u00e9e de l'\u00e9v\u00e9nement")

  if (mode == "pooled") {
    return(base(.calib_df(fit$fold_predictions, bins), "#2C7FB8") + common +
             ggplot2::labs(title = sprintf("Calibration (%s, pr\u00e9dictions pool\u00e9es, CV honn\u00eate)", lbl),
                           subtitle = sous_titre, x = labs_xy$x, y = labs_xy$y, size = "N") +
             .theme_cv(13, 8.5))
  }
  fin <- .final_predictions(fit)
  df_final <- tibble::tibble(Reel = fin$Y, Prediction = fin$pred)
  if (mode == "final") {
    return(base(.calib_df(df_final, bins), "#D95F02") + common +
             ggplot2::labs(title = sprintf("Calibration (%s) -- mod\u00e8le final sur ses donn\u00e9es d'entra\u00eenement", lbl),
                           subtitle = "\u26a0 OPTIMISTE (entra\u00eenement) : ne garantit rien sur de nouvelles donn\u00e9es \ncomparer \u00e0 la version CV via mode = 'overlay'",
                           x = labs_xy$x, y = labs_xy$y, size = "N") +
             .theme_cv(13, 8.5, "#D95F02"))
  }
  niv <- c("CV pool\u00e9e (honn\u00eate)", "Mod\u00e8le final (entra\u00eenement, optimiste)")
  combined <- dplyr::bind_rows(
    .calib_df(fit$fold_predictions, bins) %>% dplyr::mutate(Source = niv[1]),
    .calib_df(df_final, bins) %>% dplyr::mutate(Source = niv[2])
  ) %>% dplyr::mutate(Source = factor(Source, levels = niv))
  ggplot2::ggplot(combined, ggplot2::aes(x = prob_moyenne, y = freq_observee, colour = Source)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_line() + ggplot2::geom_point(ggplot2::aes(size = n), alpha = 0.8) +
    ggplot2::scale_colour_manual(values = c("#2C7FB8", "#D95F02")) + common +
    ggplot2::labs(title = sprintf("Calibration (%s) : CV honn\u00eate vs mod\u00e8le final (optimiste)", lbl),
                  subtitle = sous_titre, x = labs_xy$x, y = labs_xy$y, size = "N", colour = NULL) +
    .theme_cv(13, 8.5) + ggplot2::theme(legend.position = "bottom")
}

#' Matrice de confusion (prédictions poolées)
#' @param fit objet "cv_fit" (family = "binomial")
#' @param threshold seuil ; NULL = \code{fit$classification_threshold}
#' @export
plot_confusion_matrix <- function(fit, threshold = NULL) {
  .stop_if_not_cv_fit(fit)
  .stop_if_not_binomial(fit, "plot_confusion_matrix")
  if (is.null(threshold)) threshold <- fit$classification_threshold
  conf <- fit$fold_predictions %>%
    dplyr::mutate(Pred_class = as.numeric(Prediction >= threshold)) %>%
    dplyr::count(Reel, Pred_class) %>%
    tidyr::complete(Reel = c(0, 1), Pred_class = c(0, 1), fill = list(n = 0L)) %>%
    dplyr::mutate(Reel_lab = factor(Reel, levels = c(0, 1), labels = c("R\u00e9el : 0", "R\u00e9el : 1")),
                  Pred_lab = factor(Pred_class, levels = c(0, 1), labels = c("Pr\u00e9dit : 0", "Pr\u00e9dit : 1")))
  ggplot2::ggplot(conf, ggplot2::aes(x = Pred_lab, y = Reel_lab, fill = n)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = n), size = 6, colour = "white", fontface = "bold") +
    ggplot2::scale_fill_gradient(low = "#2C7FB8", high = "#08306B", guide = "none") +
    ggplot2::labs(title = sprintf("Matrice de confusion (%s, pr\u00e9dictions pool\u00e9es)", fit$method_label),
                  subtitle = sprintf("Seuil de classification = %.2f \n(somme = %d pr\u00e9dictions)", threshold, sum(conf$n)),
                  x = NULL, y = NULL) +
    .theme_cv(13)
}
