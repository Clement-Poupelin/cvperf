# =============================================================================
# check_model.R -- le grand panel de diagnostics, commun à toutes les méthodes
# =============================================================================

#' Panel complet de diagnostics (à la manière de performance::check_model())
#'
#' Gaussien : qualité de prédiction, résidus vs prédits, Q-Q plot,
#' scale-location, panneau "variables", outliers. Binomial : qualité de
#' prédiction, ROC, calibration, matrice de confusion, panneau "variables",
#' outliers. Le panneau "variables" dépend de la méthode : VIF (Lasso),
#' fréquence de sélection (MARS), stabilité de l'importance (autres).
#'
#' @param fit objet "cv_fit"
#' @param labels libellés optionnels des variables
#' @return patchwork
#' @export
check_model <- function(fit, labels = NULL) {
  .stop_if_not_cv_fit(fit)
  p5 <- .importance_panel(fit, labels = labels)
  p6 <- plot(check_outliers(fit))
  p1 <- plot_prediction_quality(fit)
  if (identical(fit$family, "binomial")) {
    return((p1 | plot_roc(fit)) / (plot_calibration(fit) | plot_confusion_matrix(fit)) / (p5 | p6))
  }
  (p1 | plot_residuals_vs_fitted(fit)) /
    (plot(check_normality(fit)) | plot(check_heteroscedasticity(fit))) /
    (p5 | p6)
}
