# =============================================================================
# utils.R -- utilitaires partagés par toutes les méthodes
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Étoiles de significativité conventionnelles
#'
#' @param p vecteur de p-values
#' @return vecteur character ("***", "**", "*", ".", "ns")
#' @export
p_to_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.10  ~ ".",
    TRUE      ~ "ns"
  )
}

#' AUC par la formule de Mann-Whitney / Wilcoxon rank-sum
#' @noRd
.compute_auc <- function(prob, actual) {
  actual <- as.numeric(actual)
  n1 <- sum(actual == 1); n0 <- sum(actual == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(prob)
  (sum(r[actual == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Métriques de classification (probabilités vs réponse 0/1)
#' @noRd
.compute_classification_metrics <- function(prob, actual, threshold = 0.5) {
  actual     <- as.numeric(actual)
  pred_class <- as.numeric(prob >= threshold)
  tp <- sum(pred_class == 1 & actual == 1); tn <- sum(pred_class == 0 & actual == 0)
  fp <- sum(pred_class == 1 & actual == 0); fn <- sum(pred_class == 0 & actual == 1)
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  precision   <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  f1 <- if (!is.na(precision) && !is.na(sensitivity) && (precision + sensitivity) > 0) {
    2 * precision * sensitivity / (precision + sensitivity)
  } else NA_real_
  eps    <- 1e-15
  prob_c <- pmin(pmax(prob, eps), 1 - eps)
  tibble::tibble(
    accuracy    = (tp + tn) / length(actual),
    sensitivity = sensitivity,
    specificity = specificity,
    precision   = precision,
    f1          = f1,
    auc         = .compute_auc(prob, actual),
    brier       = mean((prob - actual)^2),
    logloss     = -mean(actual * log(prob_c) + (1 - actual) * log(1 - prob_c))
  )
}

#' Métriques de régression (un fold)
#' @noRd
.compute_regression_metrics <- function(reel, prediction) {
  residu <- reel - prediction
  sst <- sum((reel - mean(reel))^2)
  sse <- sum(residu^2)
  tibble::tibble(
    sst  = sst,
    sse  = sse,
    rmse = sqrt(mean(residu^2)),
    mae  = mean(abs(residu)),
    rsq  = if (sst > 0) 1 - sse / sst else NA_real_,
    mse  = mean(residu^2)
  )
}

#' Nombre de coeurs réellement utilisés
#'
#' NULL = 2/3 des coeurs logiques détectés (au moins 1), plafonné au nombre de
#' tâches. Toujours 1 si \code{parallelisation = FALSE}.
#' @noRd
.resolve_n_cores <- function(parallelisation, n_cores, n_tasks) {
  if (!isTRUE(parallelisation)) return(1L)
  n_dispo <- parallel::detectCores(logical = TRUE)
  if (is.na(n_dispo) || n_dispo < 1) n_dispo <- 1L
  if (is.null(n_cores)) {
    n <- max(1L, floor(2 * n_dispo / 3))
  } else {
    stopifnot("n_cores doit etre un entier positif" = n_cores >= 1)
    if (n_cores > n_dispo) {
      warning(sprintf("n_cores = %d demande, mais seulement %d coeur(s) detecte(s) -- utilisation de %d.",
                      as.integer(n_cores), n_dispo, n_dispo))
    }
    n <- min(as.integer(n_cores), n_dispo)
  }
  as.integer(max(1L, min(n, max(1L, n_tasks))))
}

#' Résumé d'une ligne de `$computational_info`
#' @noRd
.format_computational_info <- function(ci) {
  sprintf(
    "Temps de calcul : %.1fs (folds) + %.1fs (modele final) = %.1fs au total | %s | %.2fs/fold en moyenne",
    ci$time_folds_sec, ci$time_final_model_sec, ci$time_total_sec,
    if (isTRUE(ci$parallelisation)) sprintf("parallelise sur %d coeurs", ci$n_cores_used) else "sequentiel (1 coeur)",
    ci$time_folds_sec / max(1, ci$n_folds)
  )
}

#' Vérifie qu'un package optionnel (Suggests) est installé
#' @noRd
.require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Le package {%s} est requis pour cette fonction : install.packages('%s').", pkg, pkg),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Contrôle de classe commun
#' @noRd
.stop_if_not_cv_fit <- function(fit, class = "cv_fit") {
  if (!inherits(fit, class)) {
    stop(sprintf(
      "`fit` doit etre un objet '%s' (sortie d'une fonction fit_*_repeated_cv()). Classe recue : %s. Verifier que le fit n'a pas ete ecrase par le resultat de model_performance().",
      class, paste(class(fit), collapse = "/")), call. = FALSE)
  }
  invisible(TRUE)
}

#' Contrôle family = "binomial"
#' @noRd
.stop_if_not_binomial <- function(fit, fun) {
  if (!identical(fit$family, "binomial")) {
    stop(sprintf("%s() requiert un modele family = 'binomial'.", fun), call. = FALSE)
  }
}

#' Libellés lisibles de variables (vecteur nommé optionnel)
#' @noRd
.label_lookup <- function(v, labels = NULL) {
  if (is.null(labels)) return(v)
  dplyr::coalesce(unname(labels[v]), v)
}

#' Thème commun des graphiques du package
#' @noRd
.theme_cv <- function(base_size = 12, subtitle_size = 9, subtitle_colour = "grey45") {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(colour = subtitle_colour, size = subtitle_size))
}

#' Affiche un data.frame comme table interactive exportable (CSV / Excel)
#'
#' Pour une sortie HTML (Quarto / R Markdown) uniquement.
#' @param df data.frame / tibble
#' @param caption titre optionnel
#' @param page_length lignes par page (défaut 10)
#' @return objet \code{DT::datatable}
#' @export
datatable_export <- function(df, caption = NULL, page_length = 10) {
  .require_pkg("DT")
  DT::datatable(
    df, caption = caption, filter = "top", rownames = FALSE, extensions = "Buttons",
    options = list(
      dom = "Bfrtip",
      buttons = list(list(extend = "csv", filename = "export"),
                     list(extend = "excel", filename = "export")),
      pageLength = page_length, scrollX = TRUE
    )
  )
}
