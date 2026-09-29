# =============================================================================
# method_rf.R -- forêt aléatoire ({ranger})
# =============================================================================
# Spécifique : mtry choisi par erreur OOB (sous-produit gratuit du bagging,
# rôle de cv.glmnet() pour lambda) ; importance native de ranger
# (permutation par défaut) ; interface `dependent.variable.name` (les noms de
# colonnes "exotiques" -- "Oligo-SNCA", "Aβ38", "4HNE" -- sont rejetés par
# l'interface formule de ranger).

.default_mtry_grid <- function(p, family) {
  centre <- if (identical(family, "binomial")) sqrt(p) else p / 3
  sort(unique(pmin(pmax(round(centre * c(0.5, 0.75, 1, 1.5, 2, 3)), 1), p)))
}

.ranger_fit <- function(df_fit, response_var, family, num.trees, mtry, min.node.size, importance) {
  ranger::ranger(dependent.variable.name = response_var, data = df_fit, num.trees = num.trees, mtry = mtry,
                 min.node.size = min.node.size, probability = identical(family, "binomial"),
                 importance = importance, oob.error = TRUE)
}

.tune_mtry_oob <- function(df_fit, response_var, family, num.trees, min.node.size, mtry_grid, seed) {
  p <- ncol(df_fit) - 1
  mtry_grid <- if (is.null(mtry_grid)) .default_mtry_grid(p, family) else unique(pmin(pmax(round(mtry_grid), 1), p))
  res <- purrr::map_dfr(mtry_grid, function(m) {
    set.seed(seed)
    f <- .ranger_fit(df_fit, response_var, family, num.trees, m, min.node.size, "none")
    tibble::tibble(mtry = m, oob_error = f$prediction.error)
  })
  best <- res$mtry[which.min(res$oob_error)]
  list(best_mtry = best, grid_results = res, at_edge = best == min(mtry_grid) || best == max(mtry_grid))
}

#' Forêt aléatoire (ranger) en validation croisée répétée
#'
#' @inheritParams fit_lasso_repeated_cv
#' @param num.trees nombre d'arbres (500)
#' @param mtry_grid grille de mtry testée par OOB (NULL = autour de sqrt(p) / p/3)
#' @param min.node.size taille minimale de noeud (NULL = défaut ranger)
#' @param importance mesure d'importance ranger ("permutation" par défaut, non
#'   biaisée par le nombre de modalités ; ou "impurity")
#' @return objet c("rf_cv_fit", "cv_fit")
#' @export
fit_rf_repeated_cv <- function(data, formula, id_col = "NUM_PAT", v = 5, repeats = 20, num.trees = 500,
                               mtry_grid = NULL, min.node.size = NULL, family = c("gaussian", "binomial"),
                               classification_threshold = 0.5, seed = 42, cv_repeats = NULL,
                               importance = "permutation", parallelisation = FALSE, n_cores = NULL) {
  .require_pkg("ranger")
  family <- match.arg(family)
  force(num.trees); force(mtry_grid); force(min.node.size); force(importance)
  response_var <- .response_var(formula)
  prep_full <- .prepare_df_data(data, formula, id_col, family)
  if (is.null(prep_full)) stop("Jeu de donnees inutilisable (facteur mono-niveau ou reponse binaire mal formee).", call. = FALSE)

  fold_fun <- function(train, test, fold_num) {
    tr <- .prepare_df_data(train, formula, id_col, family)
    te <- .prepare_df_data(test,  formula, id_col, family)
    if (is.null(tr) || is.null(te)) return(NULL)
    tun <- .tune_mtry_oob(tr$df_fit, response_var, family, num.trees, min.node.size, mtry_grid, seed + fold_num)
    set.seed(seed + fold_num)
    rf <- .ranger_fit(tr$df_fit, response_var, family, num.trees, tun$best_mtry, min.node.size, importance)
    pr <- stats::predict(rf, data = te$df_fit)$predictions
    imp <- ranger::importance(rf)
    list(Y = te$Y, pred = if (identical(family, "binomial")) as.numeric(pr[, "1"]) else as.numeric(pr),
         extras = list(
           importance = tibble::tibble(Variable = names(imp), Importance = as.numeric(imp)),
           hyperparams = tibble::tibble(mtry_used = tun$best_mtry, mtry_grid_min = min(tun$grid_results$mtry),
                                        mtry_grid_max = max(tun$grid_results$mtry),
                                        oob_error = min(tun$grid_results$oob_error), at_edge = tun$at_edge)))
  }

  cv <- run_repeated_cv(data, formula, id_col, family, v, repeats, seed, cv_repeats, fold_fun,
                        parallelisation, n_cores, pkgs = "ranger")

  fin <- .time_it({
    tun <- .tune_mtry_oob(prep_full$df_fit, response_var, family, num.trees, min.node.size, mtry_grid, seed)
    set.seed(seed)
    list(tuning = tun, model = .ranger_fit(prep_full$df_fit, response_var, family, num.trees, tun$best_mtry, min.node.size, importance))
  })

  new_cv_fit(
    cv, data, formula, id_col, family, classification_threshold, v, repeats,
    fold_importance = cv$extras$importance, fold_hyperparams = cv$extras$hyperparams,
    final_model = fin$value$model, final_time_sec = fin$time_sec,
    method = "rf", method_label = "For\u00eat al\u00e9atoire", importance_type = importance,
    final_tuning = fin$value$tuning, final_mtry = fin$value$tuning$best_mtry,
    num_trees = num.trees, min_node_size = min.node.size,
    mtry_grid = mtry_grid %||% .default_mtry_grid(ncol(prep_full$df_fit) - 1, family),
    class = "rf_cv_fit"
  )
}

#' @export
model_design.rf_cv_fit <- function(fit, ...) {
  p <- .prepare_df_data(fit$data, fit$formula, fit$id_col, fit$family)
  X <- p$df_fit[, setdiff(names(p$df_fit), fit$response_var), drop = FALSE]
  list(X = X, Y = p$Y, var_map = stats::setNames(names(X), names(X)),
       predict = function(newdata) {
         pr <- stats::predict(fit$final_model, data = as.data.frame(newdata))$predictions
         if (identical(fit$family, "binomial")) as.numeric(pr[, "1"]) else as.numeric(pr)
       })
}

#' @export
importance_table.rf_cv_fit <- function(fit, ...) {
  imp <- ranger::importance(fit$final_model)
  tibble::tibble(Variable = names(imp), Importance = as.numeric(imp)) %>% dplyr::arrange(dplyr::desc(Importance))
}

#' @export
method_notes.rf_cv_fit <- function(fit, topic) {
  switch(topic,
    heteroscedasticity = "Rappel for\u00eat : le moyennage d'arbres tire les pr\u00e9dictions vers la moyenne (shrinkage), ce qui peut \u00e9largir le r\u00e9sidu aux extr\u00eames de Y sans d\u00e9faut de sp\u00e9cification.",
    calibration = "Une for\u00eat de probabilit\u00e9 tend \u00e0 \u00ab tasser \u00bb les probabilit\u00e9s vers 0.5 (moyennage d'arbres) : surveiller un aplatissement syst\u00e9matique.",
    final = sprintf("Le mod\u00e8le final dispose aussi de son erreur OOB (%.3g), plus honn\u00eate que cette lecture d'entra\u00eenement.",
                    fit$final_model$prediction.error),
    residuals = "Une courbe qui s'\u00e9carte de 0 peut signaler une zone de la plage de Y mal couverte par les arbres",
    NextMethod())
}

#' @export
.print_details.rf_cv_fit <- function(x) {
  cat(sprintf("mtry final    : %d (grille : %s)%s\n", x$final_mtry, paste(x$mtry_grid, collapse = ", "),
              if (isTRUE(x$final_tuning$at_edge)) "  \u26a0 AU BORD -> elargir `mtry_grid`" else ""))
  cat(sprintf("num.trees = %d | min.node.size = %s | erreur OOB (modele final) = %.4g\n", x$num_trees,
              if (is.null(x$min_node_size)) "defaut ranger" else x$min_node_size, x$final_model$prediction.error))
}

#' @export
check_hyperparam_stability.rf_cv_fit <- function(fit) {
  structure(list(final_at_edge = isTRUE(fit$final_tuning$at_edge),
                 pct_folds_edge = 100 * mean(fit$fold_hyperparams$at_edge, na.rm = TRUE)),
            class = "check_hyperparam_stability_rf")
}

#' @export
print.check_hyperparam_stability_rf <- function(x, ...) {
  cat(sprintf("Modele final : mtry retenu %s au bord de sa grille.\n", if (x$final_at_edge) "EST" else "n'est PAS"))
  cat(sprintf("%.1f%% des folds ont un mtry au bord de LEUR grille.\n", x$pct_folds_edge))
  if (x$final_at_edge) cat("\u26a0 Grille mal calee : elargir `mtry_grid`.\n")
  else if (x$pct_folds_edge > 30) cat("\u2139 Proportion elevee de folds au bord (modele final OK).\n")
  else cat("\u2713 Grille correctement calee.\n")
  invisible(x)
}

#' @export
plot_hyperparam_path.rf_cv_fit <- function(fit) {
  df <- fit$final_tuning$grid_results
  bord <- isTRUE(fit$final_tuning$at_edge)
  ggplot2::ggplot(df, ggplot2::aes(x = mtry, y = oob_error)) +
    ggplot2::geom_line(colour = "grey60") + ggplot2::geom_point(colour = "#D95F02", size = 2.2) +
    ggplot2::geom_vline(xintercept = fit$final_mtry, colour = "#2C3E50", linewidth = 0.8) +
    ggplot2::annotate("text", x = fit$final_mtry, y = max(df$oob_error), label = sprintf("mtry = %d", fit$final_mtry),
                      colour = "#2C3E50", hjust = -0.1, size = 3.2) +
    ggplot2::labs(title = "Erreur OOB en fonction de mtry (mod\u00e8le final)",
                  subtitle = if (bord) "\u26a0 mtry retenu au bord de la grille : \u00e9largir `mtry_grid`" else "Trait vertical = mtry retenu (OOB minimale)",
                  x = "mtry (variables tir\u00e9es \u00e0 chaque split)",
                  y = if (identical(fit$family, "binomial")) "Erreur OOB (Brier)" else "Erreur OOB (MSE)") +
    .theme_cv(12, 9, if (bord) "#D95F02" else "grey45")
}

#' @export
plot_hyperparam_stability.rf_cv_fit <- function(fit) {
  hp <- fit$fold_hyperparams
  ggplot2::ggplot(hp, ggplot2::aes(x = mtry_used)) +
    ggplot2::geom_bar(fill = "#2C7FB8", alpha = 0.7, colour = "white") +
    ggplot2::geom_vline(xintercept = stats::median(hp$mtry_used), colour = "#D95F02", linetype = "dashed", linewidth = 0.8) +
    ggplot2::labs(title = "Stabilit\u00e9 de mtry \u00e0 travers les folds",
                  subtitle = sprintf("M\u00e9diane = %.0f \u00b7 distribution large = complexit\u00e9 optimale instable", stats::median(hp$mtry_used)),
                  x = "mtry retenu (par fold)", y = "Nombre de folds") +
    .theme_cv(12)
}

#' @export
variables_summary.rf_cv_fit <- function(fit) {
  dplyr::bind_cols(NextMethod(), tibble::tibble(mtry_final = fit$final_mtry))
}
