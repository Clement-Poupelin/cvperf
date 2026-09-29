# =============================================================================
# performance.R -- model_performance(), commun à toutes les méthodes
# =============================================================================

#' Indices de performance en validation croisée répétée
#'
#' Deux lectures complémentaires :
#' \itemize{
#'   \item \code{pooled} : toutes les prédictions de test poolées avant de
#'     calculer un seul jeu de métriques (le plus stable) ;
#'   \item \code{par_fold} / \code{par_fold_summary} : métriques fold par fold,
#'     puis résumées (montre la variabilité d'un tirage à l'autre).
#' }
#' Gaussien : R², RMSE, MAE, MSE. Binomial : accuracy, sensibilité,
#' spécificité, précision, F1, AUC, Brier, log-loss.
#'
#' @param model objet "cv_fit"
#' @param ... inutilisé
#' @return liste de classe "performance_cv"
#' @export
model_performance <- function(model, ...) UseMethod("model_performance")

#' Métriques fold par fold (partagé avec compare_methods())
#' @noRd
.per_fold_metrics <- function(fp, family, threshold = 0.5) {
  if (identical(family, "binomial")) {
    fp %>%
      dplyr::group_by(Fold) %>%
      dplyr::group_modify(~ .compute_classification_metrics(.x$Prediction, .x$Reel, threshold)) %>%
      dplyr::ungroup()
  } else {
    fp %>%
      dplyr::group_by(Fold) %>%
      dplyr::group_modify(~ .compute_regression_metrics(.x$Reel, .x$Prediction)) %>%
      dplyr::ungroup()
  }
}

#' @rdname model_performance
#' @export
model_performance.cv_fit <- function(model, ...) {
  fp <- model$fold_predictions
  par_fold <- .per_fold_metrics(fp, model$family, model$classification_threshold)

  if (identical(model$family, "binomial")) {
    pooled <- .compute_classification_metrics(fp$Prediction, fp$Reel, model$classification_threshold) %>%
      dplyr::mutate(N = nrow(fp)) %>%
      dplyr::relocate("N")
    par_fold_summary <- par_fold %>%
      dplyr::summarise(dplyr::across(
        c(accuracy, sensitivity, specificity, precision, f1, auc, brier, logloss),
        list(moyen = ~ mean(.x, na.rm = TRUE), sd = ~ stats::sd(.x, na.rm = TRUE)),
        .names = "{.col}_{.fn}"
      ))
  } else {
    sst_pool <- sum((fp$Reel - mean(fp$Reel))^2)
    pooled <- tibble::tibble(
      R2   = 1 - sum(fp$Residu^2) / sst_pool,
      RMSE = sqrt(mean(fp$Residu^2)),
      MAE  = mean(abs(fp$Residu)),
      MSE  = mean(fp$Residu^2),
      N    = nrow(fp)
    )
    par_fold_summary <- par_fold %>%
      dplyr::summarise(
        R2_moyen = mean(rsq, na.rm = TRUE), R2_median = stats::median(rsq, na.rm = TRUE),
        R2_sd = stats::sd(rsq, na.rm = TRUE),
        RMSE_moyen = mean(rmse), RMSE_sd = stats::sd(rmse),
        MAE_moyen = mean(mae),
        MSE_moyen = mean(mse), MSE_sd = stats::sd(mse)
      )
  }
  structure(
    list(pooled = pooled, par_fold = par_fold, par_fold_summary = par_fold_summary,
         family = model$family, method_label = model$method_label),
    class = "performance_cv"
  )
}

#' @export
print.performance_cv <- function(x, ...) {
  cat(sprintf("Performance poolee -- %s (%d folds, un seul jeu de metriques global) :\n",
              x$method_label %||% "", nrow(x$par_fold)))
  print(knitr::kable(x$pooled, digits = 3))
  cat("\nPerformance par fold, resumee (variabilite entre folds) :\n")
  print(knitr::kable(x$par_fold_summary, digits = 3))
  invisible(x)
}
