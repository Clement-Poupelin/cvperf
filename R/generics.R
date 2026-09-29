# =============================================================================
# generics.R -- les "crochets" S3 qu'une méthode implémente
# =============================================================================
# Toute la couche commune (diagnostics, plots, SHAP, comparaison) ne s'appuie
# que sur les champs communs d'un "cv_fit" et sur ces génériques :
#
#   model_design(fit)       [OBLIGATOIRE] X, Y, predict(newdata), var_map du
#                           modèle FINAL -> mode "final" des plots,
#                           plot_effects(), check_shap()
#   importance_table(fit)   table d'importance du modèle
#   method_notes(fit, topic)  avertissements d'interprétation propres à la méthode
#   check_hyperparam_stability(), plot_hyperparam_path(), plot_hyperparam_stability()
#   plot_variable_importance(), variables_summary()
#   .print_details(), .print_body(), .importance_panel(), .effect_variables()

#' Design du modèle final (crochet principal d'une méthode)
#'
#' @param fit objet "cv_fit"
#' @param ... inutilisé
#' @return liste : \code{X} (matrice ou data.frame des prédicteurs, tel que le
#'   modèle les consomme), \code{Y} (réponse numérique, 0/1 si binomial),
#'   \code{predict} (fonction(newdata) -> prédictions sur l'échelle de la
#'   réponse, P(Y=1) si binomial), \code{var_map} (vecteur nommé colonne de X ->
#'   variable d'origine)
#' @export
model_design <- function(fit, ...) UseMethod("model_design")

#' @export
model_design.default <- function(fit, ...) {
  stop("model_design() n'est pas implemente pour la classe ", paste(class(fit), collapse = "/"),
       ". Toute nouvelle methode doit definir model_design.<classe>().", call. = FALSE)
}

#' Prédictions du modèle FINAL sur ses propres données d'entraînement
#' (optimistes par construction)
#' @noRd
.final_predictions <- function(fit) {
  d <- model_design(fit)
  list(Y = d$Y, pred = as.numeric(d$predict(d$X)))
}

#' Table d'importance des variables
#'
#' Méthode par défaut : importance moyenne à travers les folds de CV
#' (\code{fit$fold_importance}). Les méthodes peuvent la spécialiser
#' (coefficients pour le Lasso, gain natif pour XGBoost, evimp pour MARS...).
#' @param fit objet "cv_fit"
#' @param ... arguments propres à la méthode
#' @export
importance_table <- function(fit, ...) UseMethod("importance_table")

#' @export
importance_table.cv_fit <- function(fit, ...) {
  fit$fold_importance %>%
    dplyr::group_by(Variable) %>%
    dplyr::summarise(
      mean_importance    = mean(Importance),
      sd_importance      = stats::sd(Importance),
      pct_folds_positive = 100 * mean(Importance > 0),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_importance))
}

#' Avertissements d'interprétation propres à une méthode
#'
#' @param fit objet "cv_fit"
#' @param topic "normality", "heteroscedasticity", "calibration", "final",
#'   "residuals"
#' @return chaîne de caractères ou NULL
#' @export
method_notes <- function(fit, topic) UseMethod("method_notes")

#' @export
method_notes.cv_fit <- function(fit, topic) {
  switch(topic,
    normality = sprintf(
      "DESCRIPTIF uniquement : %s ne suppose pas la normalite des residus -- ce test sert a comparer la FORME des residus entre methodes, pas a diagnostiquer une mauvaise specification.",
      fit$method_label),
    NULL
  )
}

#' Stabilité des hyperparamètres à travers les folds
#' @param fit objet "cv_fit"
#' @export
check_hyperparam_stability <- function(fit) UseMethod("check_hyperparam_stability")

#' Chemin de réglage des hyperparamètres (modèle final)
#' @param fit objet "cv_fit"
#' @export
plot_hyperparam_path <- function(fit) UseMethod("plot_hyperparam_path")

#' Distribution des hyperparamètres retenus à travers les folds
#' @param fit objet "cv_fit"
#' @export
plot_hyperparam_stability <- function(fit) UseMethod("plot_hyperparam_stability")

#' @export
check_hyperparam_stability.default <- function(fit) stop("Non implemente pour cette classe.", call. = FALSE)
#' @export
plot_hyperparam_path.default <- function(fit) stop("Non implemente pour cette classe.", call. = FALSE)
#' @export
plot_hyperparam_stability.default <- function(fit) stop("Non implemente pour cette classe.", call. = FALSE)

#' Résumé "variables en entrée / en sortie" d'un fit (utilisé par compare_models())
#' @param fit objet "cv_fit"
#' @export
variables_summary <- function(fit) UseMethod("variables_summary")

#' @export
variables_summary.cv_fit <- function(fit) {
  imp <- importance_table(fit)
  tibble::tibble(
    N_variables_entree = length(.vars_candidates(fit)),
    Variables_top5     = paste(utils::head(imp$Variable, 5), collapse = ", ")
  )
}

# ── Crochets internes ---------------------------------------------------------

#' Lignes spécifiques (réglage, modèle final) affichées par print() avant la performance
#' @noRd
.print_details <- function(x) UseMethod(".print_details")
#' @export
.print_details.cv_fit <- function(x) invisible(NULL)

#' Bloc spécifique affiché par print() après la performance
#' @noRd
.print_body <- function(x) UseMethod(".print_body")
#' @export
.print_body.cv_fit <- function(x) {
  imp <- importance_table(x)
  n <- min(10, nrow(imp))
  cat(sprintf("\nImportance des variables (%s) -- top %d :\n", x$importance_type, n))
  print(utils::head(imp, n), n = Inf)
  if (!is.null(x$fold_group_importance) && nrow(x$fold_group_importance) > 0) {
    cat("\nImportance des groupes de variables (permutation d'un bloc) :\n")
    print(group_importance(x), n = Inf)
  }
  invisible(NULL)
}

#' Panneau "variables" du panel check_model()
#' @noRd
.importance_panel <- function(fit, labels = NULL) UseMethod(".importance_panel")
#' @export
.importance_panel.cv_fit <- function(fit, labels = NULL) {
  plot_importance(check_importance(fit), labels = labels, fit = fit)
}

#' Variables tracées par défaut par plot_effects()
#' @noRd
.effect_variables <- function(fit, top_n) UseMethod(".effect_variables")
#' @export
.effect_variables.cv_fit <- function(fit, top_n) {
  utils::head(importance_table.cv_fit(fit)$Variable, top_n)
}

# ── print() commun ------------------------------------------------------------

#' @export
print.cv_fit <- function(x, ...) {
  titre <- sprintf("\u2500\u2500 %s en validation crois\u00e9e r\u00e9p\u00e9t\u00e9e ", x$method_label)
  cat(titre, strrep("\u2500", max(3, 68 - nchar(titre))), "\n", sep = "")
  cat(sprintf("Formule       : %s\n", paste(deparse(x$formula), collapse = " ")))
  cat(sprintf("Famille       : %s%s\n", x$family,
              if (identical(x$family, "binomial")) sprintf(" (classe positive = \"%s\", seuil = %.2f)",
                                                           x$positive_class, x$classification_threshold) else ""))
  cat(sprintf("N individus   : %d | CV : %d folds x %d repetitions | folds valides : %d / %d\n",
              nrow(x$data), x$v, x$repeats, x$n_folds_valides, x$n_folds_total))
  .print_details(x)
  if (!is.null(x$computational_info)) cat(.format_computational_info(x$computational_info), "\n")

  perf <- model_performance(x)$pooled
  if (identical(x$family, "binomial")) {
    cat(sprintf("Performance (poolee, %d predictions) : Accuracy = %.3f | AUC = %.3f | Brier = %.3f\n",
                perf$N, perf$accuracy, perf$auc, perf$brier))
  } else {
    cat(sprintf("Performance (poolee, %d predictions) : R\u00b2 = %.3f | RMSE = %.3g | MAE = %.3g\n",
                perf$N, perf$R2, perf$RMSE, perf$MAE))
  }
  .print_body(x)
  cat(strrep("\u2500", 68), "\n", sep = "")
  cat("Fonctions associees : model_performance(), check_model(), plot_prediction_quality(),\n")
  cat("  plot_variable_importance(), plot_effects(), check_hyperparam_stability(), check_shap()\n")
  invisible(x)
}
