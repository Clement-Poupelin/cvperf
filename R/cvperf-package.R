#' cvperf : performance de modèles prédictifs en validation croisée répétée
#'
#' Architecture en deux couches :
#' \itemize{
#'   \item une couche COMMUNE (classe \code{"cv_fit"}) : moteur de CV répétée,
#'     métriques, diagnostics de résidus, courbes ROC/calibration, importance,
#'     dépendance partielle, SHAP, comparaison de méthodes ;
#'   \item une couche SPÉCIFIQUE par méthode (\code{fit_lasso_repeated_cv()},
#'     \code{fit_rf_repeated_cv()}, \code{fit_xgb_repeated_cv()},
#'     \code{fit_svm_repeated_cv()}, \code{fit_mars_repeated_cv()}) qui
#'     n'implémente que l'ajustement, le réglage des hyperparamètres et
#'     quelques crochets S3 (\code{model_design()}, \code{importance_table()},
#'     \code{method_notes()}, hyperparamètres).
#' }
#' Ajouter une méthode = écrire un fichier \code{R/method_xxx.R} qui appelle
#' \code{.run_repeated_cv()} et implémente ces crochets.
#'
#' @keywords internal
#' @importFrom dplyr %>%
#' @importFrom rlang .data :=
#' @importFrom patchwork plot_annotation plot_layout wrap_plots
"_PACKAGE"

utils::globalVariables(c(
  ".", "Fold", "Repeat", "ID", "Reel", "Prediction", "Residu", "Residu_lag",
  "Variable", "Importance", "Rank", "Used", "Coefficient", "Label", "Direction",
  "mean_importance", "sd_importance", "rank_sd", "mean_rank", "in_topk",
  "pct_in_topk", "n_in_topk", "sqrt_abs_resid", "residu_moyen", "residu_sd",
  "is_outlier", "FPR", "TPR", "Recall", "Precision", "Source", "bin",
  "prob_moyenne", "freq_observee", "n", "Pred_class", "Pred_lab", "Reel_lab",
  "pred_min", "pred_max", "pred_moyenne", "pred_final_moyenne", "x", "y", "m",
  "Modele", "Metric", "Rang", "rang_moyen", "Score", "mean_abs_shap", "shap_value",
  "feature_value", "feature_value_scaled", ".id", "best_rank", "Composante",
  "temps_sec", "temps_grand_total_sec", "temps_shap_sec", "Term", "VIF",
  "Category", "n_selected", "pct_selected", "mean_coef_sel", "log_lambda", "cvm",
  "cvlo", "cvup", "lambda_min", "mtry", "oob_error", "mtry_used", "nrounds",
  "error", "depth_f", "max_depth_used", "eta_used", "nrounds_used", "cost",
  "cost_used", "gamma_used", "pct_sv", "cost_lab", "gamma_lab", "degree_f",
  "degree_used", "n_terms", "n_used", "pct_used", "mean_importance_when_used",
  "pct_folds_used", "Gain", "Groupe", "selectionnee_ce_fold",
  "Lasso_pct_folds_selectionnee", "Lasso_selectionnee", "n_methodes",
  "sst", "sse", "rmse", "mae", "rsq", "mse", "p_value", "p_adj",
  "accuracy", "sensitivity", "specificity", "precision", "f1", "auc", "brier",
  "logloss", "rank", "nprune", "density", "mean_shap", "eta_lab", "Tolerance"
))
