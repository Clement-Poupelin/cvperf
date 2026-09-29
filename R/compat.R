# =============================================================================
# compat.R -- anciens noms des toolkits (GENERE, ne pas editer a la main)
# =============================================================================
# Les tutoriels et scripts ecrits avec les toolkits sources (plot_roc_rf(),
# check_model_xgb(), plot.check_normality_svm()...) continuent de fonctionner :
# chaque alias verifie la classe puis appelle la fonction commune.
# Pour du nouveau code, preferer les noms generiques (plot_roc(), check_model()...).

#' Anciens noms des toolkits (compatibilite)
#'
#' Alias vers les fonctions communes, pour que les scripts ecrits avec les
#' toolkits \code{*_performance_toolkit.R} fonctionnent sans modification.
#' @param fit,x,fits,model objets d'entree
#' @param ... arguments transmis
#' @name cvperf-compat
#' @rawNamespace export(plot.check_normality_lasso)
#' @rawNamespace export(plot.check_heteroscedasticity_lasso)
#' @rawNamespace export(plot.check_outliers_lasso)
#' @rawNamespace export(plot.check_independence_lasso)
#' @rawNamespace export(plot.check_normality_rf)
#' @rawNamespace export(plot.check_heteroscedasticity_rf)
#' @rawNamespace export(plot.check_outliers_rf)
#' @rawNamespace export(plot.check_independence_rf)
#' @rawNamespace export(plot.check_normality_xgb)
#' @rawNamespace export(plot.check_heteroscedasticity_xgb)
#' @rawNamespace export(plot.check_outliers_xgb)
#' @rawNamespace export(plot.check_independence_xgb)
#' @rawNamespace export(plot.check_normality_svm)
#' @rawNamespace export(plot.check_heteroscedasticity_svm)
#' @rawNamespace export(plot.check_outliers_svm)
#' @rawNamespace export(plot.check_independence_svm)
#' @rawNamespace export(plot.check_normality_mars)
#' @rawNamespace export(plot.check_heteroscedasticity_mars)
#' @rawNamespace export(plot.check_outliers_mars)
#' @rawNamespace export(plot.check_independence_mars)
NULL

# ---- lasso ----
#' @export
model_performance.lasso_cv_fit <- function(model, ...) model_performance.cv_fit(model, ...)

#' @rdname cvperf-compat
#' @export
check_normality_lasso <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "lasso_cv_fit")
  check_normality(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_heteroscedasticity_lasso <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "lasso_cv_fit")
  check_heteroscedasticity(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_outliers_lasso <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "lasso_cv_fit")
  check_outliers(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_independence_lasso <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "lasso_cv_fit")
  check_independence(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_roc_lasso <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "lasso_cv_fit")
  plot_roc(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_precision_recall_lasso <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "lasso_cv_fit")
  plot_precision_recall(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_calibration_lasso <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "lasso_cv_fit")
  plot_calibration(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_confusion_matrix_lasso <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "lasso_cv_fit")
  plot_confusion_matrix(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_model_lasso <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "lasso_cv_fit")
  check_model(fit, ...)
}

plot.check_normality_lasso <- function(x, ...) {
  plot(x, ...)
}

plot.check_heteroscedasticity_lasso <- function(x, ...) {
  plot(x, ...)
}

plot.check_outliers_lasso <- function(x, ...) {
  plot(x, ...)
}

plot.check_independence_lasso <- function(x, ...) {
  plot(x, ...)
}

#' @rdname cvperf-compat
#' @export
compare_models_lasso <- function(fits, ...) {
  compare_models(fits, ...)
}

# ---- rf ----
#' @export
model_performance.rf_cv_fit <- function(model, ...) model_performance.cv_fit(model, ...)

#' @rdname cvperf-compat
#' @export
check_normality_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  check_normality(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_heteroscedasticity_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  check_heteroscedasticity(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_outliers_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  check_outliers(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_independence_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  check_independence(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_roc_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  plot_roc(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_precision_recall_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  plot_precision_recall(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_calibration_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  plot_calibration(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_confusion_matrix_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  plot_confusion_matrix(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_model_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  check_model(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_prediction_quality_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  plot_prediction_quality(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_residuals_vs_fitted_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  plot_residuals_vs_fitted(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_residuals_distribution_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  plot_residuals_distribution(fit, ...)
}

plot.check_normality_rf <- function(x, ...) {
  plot(x, ...)
}

plot.check_heteroscedasticity_rf <- function(x, ...) {
  plot(x, ...)
}

plot.check_outliers_rf <- function(x, ...) {
  plot(x, ...)
}

plot.check_independence_rf <- function(x, ...) {
  plot(x, ...)
}

#' @rdname cvperf-compat
#' @export
check_hyperparam_stability_rf <- function(fit) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  check_hyperparam_stability(fit)
}

#' @rdname cvperf-compat
#' @export
plot_hyperparam_path_rf <- function(fit) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  plot_hyperparam_path(fit)
}

#' @rdname cvperf-compat
#' @export
plot_hyperparam_stability_rf <- function(fit) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  plot_hyperparam_stability(fit)
}

#' @rdname cvperf-compat
#' @export
plot_variable_importance_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  plot_variable_importance(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_effects_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  plot_effects(fit, ...)
}

#' @rdname cvperf-compat
#' @export
importance_table_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  importance_table(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_importance_rf <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "rf_cv_fit")
  check_importance(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_importance_rf <- function(x, ...) {
  plot_importance(x, ...)
}

#' @rdname cvperf-compat
#' @export
compare_models_rf <- function(fits, ...) {
  compare_models(fits, ...)
}

# ---- xgb ----
#' @export
model_performance.xgb_cv_fit <- function(model, ...) model_performance.cv_fit(model, ...)

#' @rdname cvperf-compat
#' @export
check_normality_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  check_normality(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_heteroscedasticity_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  check_heteroscedasticity(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_outliers_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  check_outliers(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_independence_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  check_independence(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_roc_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  plot_roc(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_precision_recall_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  plot_precision_recall(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_calibration_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  plot_calibration(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_confusion_matrix_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  plot_confusion_matrix(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_model_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  check_model(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_prediction_quality_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  plot_prediction_quality(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_residuals_vs_fitted_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  plot_residuals_vs_fitted(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_residuals_distribution_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  plot_residuals_distribution(fit, ...)
}

plot.check_normality_xgb <- function(x, ...) {
  plot(x, ...)
}

plot.check_heteroscedasticity_xgb <- function(x, ...) {
  plot(x, ...)
}

plot.check_outliers_xgb <- function(x, ...) {
  plot(x, ...)
}

plot.check_independence_xgb <- function(x, ...) {
  plot(x, ...)
}

#' @rdname cvperf-compat
#' @export
check_hyperparam_stability_xgb <- function(fit) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  check_hyperparam_stability(fit)
}

#' @rdname cvperf-compat
#' @export
plot_hyperparam_path_xgb <- function(fit) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  plot_hyperparam_path(fit)
}

#' @rdname cvperf-compat
#' @export
plot_hyperparam_stability_xgb <- function(fit) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  plot_hyperparam_stability(fit)
}

#' @rdname cvperf-compat
#' @export
plot_variable_importance_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  plot_variable_importance(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_effects_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  plot_effects(fit, ...)
}

#' @rdname cvperf-compat
#' @export
importance_table_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  importance_table(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_importance_xgb <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  check_importance(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_importance_xgb <- function(x, ...) {
  plot_importance(x, ...)
}

#' @rdname cvperf-compat
#' @export
group_importance_xgb <- function(fit) {
  .stop_if_not_cv_fit(fit, "xgb_cv_fit")
  group_importance(fit)
}

# ---- svm ----
#' @export
model_performance.svm_cv_fit <- function(model, ...) model_performance.cv_fit(model, ...)

#' @rdname cvperf-compat
#' @export
check_normality_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  check_normality(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_heteroscedasticity_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  check_heteroscedasticity(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_outliers_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  check_outliers(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_independence_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  check_independence(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_roc_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  plot_roc(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_precision_recall_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  plot_precision_recall(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_calibration_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  plot_calibration(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_confusion_matrix_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  plot_confusion_matrix(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_model_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  check_model(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_prediction_quality_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  plot_prediction_quality(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_residuals_vs_fitted_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  plot_residuals_vs_fitted(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_residuals_distribution_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  plot_residuals_distribution(fit, ...)
}

plot.check_normality_svm <- function(x, ...) {
  plot(x, ...)
}

plot.check_heteroscedasticity_svm <- function(x, ...) {
  plot(x, ...)
}

plot.check_outliers_svm <- function(x, ...) {
  plot(x, ...)
}

plot.check_independence_svm <- function(x, ...) {
  plot(x, ...)
}

#' @rdname cvperf-compat
#' @export
check_hyperparam_stability_svm <- function(fit) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  check_hyperparam_stability(fit)
}

#' @rdname cvperf-compat
#' @export
plot_hyperparam_path_svm <- function(fit) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  plot_hyperparam_path(fit)
}

#' @rdname cvperf-compat
#' @export
plot_hyperparam_stability_svm <- function(fit) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  plot_hyperparam_stability(fit)
}

#' @rdname cvperf-compat
#' @export
plot_variable_importance_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  plot_variable_importance(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_effects_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  plot_effects(fit, ...)
}

#' @rdname cvperf-compat
#' @export
importance_table_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  importance_table(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_importance_svm <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  check_importance(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_importance_svm <- function(x, ...) {
  plot_importance(x, ...)
}

#' @rdname cvperf-compat
#' @export
group_importance_svm <- function(fit) {
  .stop_if_not_cv_fit(fit, "svm_cv_fit")
  group_importance(fit)
}

# ---- mars ----
#' @export
model_performance.mars_cv_fit <- function(model, ...) model_performance.cv_fit(model, ...)

#' @rdname cvperf-compat
#' @export
check_normality_mars <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  check_normality(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_heteroscedasticity_mars <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  check_heteroscedasticity(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_outliers_mars <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  check_outliers(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_independence_mars <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  check_independence(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_roc_mars <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  plot_roc(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_precision_recall_mars <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  plot_precision_recall(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_calibration_mars <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  plot_calibration(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_confusion_matrix_mars <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  plot_confusion_matrix(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_model_mars <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  check_model(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_prediction_quality_mars <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  plot_prediction_quality(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_residuals_vs_fitted_mars <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  plot_residuals_vs_fitted(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_residuals_distribution_mars <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  plot_residuals_distribution(fit, ...)
}

plot.check_normality_mars <- function(x, ...) {
  plot(x, ...)
}

plot.check_heteroscedasticity_mars <- function(x, ...) {
  plot(x, ...)
}

plot.check_outliers_mars <- function(x, ...) {
  plot(x, ...)
}

plot.check_independence_mars <- function(x, ...) {
  plot(x, ...)
}

#' @rdname cvperf-compat
#' @export
check_hyperparam_stability_mars <- function(fit) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  check_hyperparam_stability(fit)
}

#' @rdname cvperf-compat
#' @export
plot_hyperparam_path_mars <- function(fit) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  plot_hyperparam_path(fit)
}

#' @rdname cvperf-compat
#' @export
plot_hyperparam_stability_mars <- function(fit) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  plot_hyperparam_stability(fit)
}

#' @rdname cvperf-compat
#' @export
plot_variable_importance_mars <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  plot_variable_importance(fit, ...)
}

#' @rdname cvperf-compat
#' @export
plot_effects_mars <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  plot_effects(fit, ...)
}

#' @rdname cvperf-compat
#' @export
importance_table_mars <- function(fit, ...) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  importance_table(fit, ...)
}

#' @rdname cvperf-compat
#' @export
check_hyperparam_stability_lasso <- function(fit) check_lambda_range_lasso(fit)
