# =============================================================================
# method_lasso.R -- Lasso / Elastic Net / Ridge ({glmnet})
# =============================================================================
# Spécifique : lambda choisi par cv.glmnet() dans chaque fold ; coefficients
# (signe + amplitude, effet PARTIEL) ; VIF (check_collinearity_lasso()) ;
# chemin et stabilité de lambda. Tout le reste vient de la couche commune.
#
# fold_importance : |coefficient| max parmi les colonnes d'une même variable
# d'origine (facteur étalé par model.matrix()), Used = coefficient != 0.

#' Lasso (glmnet) en validation croisée répétée
#'
#' @param data data.frame complet (AVANT tout split)
#' @param formula formule, ex. \code{Y ~ .} (glmnet ajuste son propre intercept)
#' @param id_col colonne identifiant, exclue des prédicteurs
#' @param v folds par répétition (5)
#' @param repeats répétitions (20)
#' @param alpha 1 = Lasso, 0 = Ridge, entre les deux = Elastic Net
#' @param lambda grille de lambda (NULL = grille automatique de glmnet)
#' @param family "gaussian" ou "binomial"
#' @param classification_threshold seuil de décision (métriques uniquement)
#' @param seed graine (découpage + folds internes : seed + numéro de fold)
#' @param cv_repeats objet rsample à PARTAGER entre méthodes pour les comparer
#' @param parallelisation,n_cores parallélisation des folds (résultat identique)
#' @return objet de classe c("lasso_cv_fit", "cv_fit")
#' @export
fit_lasso_repeated_cv <- function(data, formula, id_col = "NUM_PAT", v = 5, repeats = 20,
                                  alpha = 1, lambda = NULL, family = c("gaussian", "binomial"),
                                  classification_threshold = 0.5, seed = 42, cv_repeats = NULL,
                                  parallelisation = FALSE, n_cores = NULL) {
  .require_pkg("glmnet")
  family <- match.arg(family)
  force(alpha); force(lambda)

  prep_full <- .prepare_matrix_data(data, formula, id_col, family, pkg = "glmnet")
  if (is.null(prep_full)) stop("Jeu de donnees inutilisable (facteur mono-niveau ou reponse binaire mal formee).", call. = FALSE)
  vc <- .vars_candidates(list(formula = formula, data = data, id_col = id_col))
  var_map <- stats::setNames(.map_to_original_variable(colnames(prep_full$X), vc), colnames(prep_full$X))

  at_edge <- function(cv) {
    isTRUE(all.equal(cv$lambda.min, min(cv$lambda))) || isTRUE(all.equal(cv$lambda.min, max(cv$lambda)))
  }

  fold_fun <- function(train, test, fold_num) {
    tr <- .prepare_matrix_data(train, formula, id_col, family, pkg = "glmnet")
    te <- .prepare_matrix_data(test,  formula, id_col, family, pkg = "glmnet")
    if (is.null(tr) || is.null(te) || !identical(colnames(tr$X), colnames(te$X))) return(NULL)
    cv_fit <- glmnet::cv.glmnet(tr$X, tr$Y, alpha = alpha, lambda = lambda, family = family)
    pred <- as.numeric(stats::predict(cv_fit, s = "lambda.min", newx = te$X, type = "response"))
    cm <- as.matrix(stats::coef(cv_fit, s = "lambda.min"))
    list(Y = te$Y, pred = pred, extras = list(
      coefs = tibble::tibble(Variable = rownames(cm), Coefficient = cm[, 1]) %>% dplyr::filter(Variable != "(Intercept)"),
      hyperparams = tibble::tibble(lambda_min = cv_fit$lambda.min, lambda_range_min = min(cv_fit$lambda),
                                   lambda_range_max = max(cv_fit$lambda), at_edge = at_edge(cv_fit))
    ))
  }

  cv <- run_repeated_cv(data, formula, id_col, family, v, repeats, seed, cv_repeats, fold_fun,
                        parallelisation, n_cores, pkgs = "glmnet")

  fin <- .time_it({
    set.seed(seed)
    final_cv <- glmnet::cv.glmnet(prep_full$X, prep_full$Y, alpha = alpha, lambda = lambda, family = family)
    list(cv = final_cv,
         model = glmnet::glmnet(prep_full$X, prep_full$Y, alpha = alpha, lambda = final_cv$lambda.min, family = family))
  })

  fold_coefs <- cv$extras$coefs
  fold_importance <- fold_coefs %>%
    dplyr::mutate(Variable = dplyr::coalesce(unname(var_map[Variable]), Variable)) %>%
    dplyr::group_by(Variable, Fold) %>%
    dplyr::summarise(Importance = max(abs(Coefficient)), Used = any(Coefficient != 0), .groups = "drop")

  new_cv_fit(
    cv, data, formula, id_col, family, classification_threshold, v, repeats,
    fold_importance = fold_importance, fold_hyperparams = cv$extras$hyperparams,
    final_model = fin$value$model, final_time_sec = fin$time_sec,
    method = "lasso",
    method_label = if (alpha == 1) "Lasso" else if (alpha == 0) "Ridge" else sprintf("Elastic Net (alpha=%.2g)", alpha),
    importance_type = "|coefficient|",
    fold_coefs = fold_coefs, lambda_traj = cv$extras$hyperparams,
    final_cv = fin$value$cv, final_lambda = fin$value$cv$lambda.min,
    alpha = alpha, x_names = colnames(prep_full$X), var_map = var_map,
    class = "lasso_cv_fit"
  )
}

#' @export
model_design.lasso_cv_fit <- function(fit, ...) {
  p <- .prepare_matrix_data(fit$data, fit$formula, fit$id_col, fit$family, pkg = "glmnet")
  list(X = p$X, Y = p$Y, var_map = fit$var_map,
       predict = function(newdata) as.numeric(stats::predict(fit$final_model, newx = as.matrix(newdata),
                                                             s = fit$final_lambda, type = "response")))
}

#' @export
method_notes.lasso_cv_fit <- function(fit, topic) {
  switch(topic,
    normality = NULL,  # hypothèse "classique" d'un modèle linéaire : pas de mise en garde
    NextMethod())
}

#' Coefficients non nuls du modèle final (jamais tronqués)
#'
#' @param fit objet "lasso_cv_fit"
#' @return tibble Variable / Coefficient / VariableOriginale, par |coefficient| décroissant
#' @export
coef_table_lasso <- function(fit) {
  .stop_if_not_cv_fit(fit, "lasso_cv_fit")
  cf <- as.matrix(stats::coef(fit$final_model, s = fit$final_lambda))
  tibble::tibble(Variable = rownames(cf), Coefficient = cf[, 1]) %>%
    dplyr::filter(Variable != "(Intercept)", Coefficient != 0) %>%
    dplyr::mutate(VariableOriginale = .map_to_original_variable(Variable, .vars_candidates(fit))) %>%
    dplyr::arrange(dplyr::desc(abs(Coefficient)))
}

#' @export
importance_table.lasso_cv_fit <- function(fit, ...) coef_table_lasso(fit)

#' @export
variables_summary.lasso_cv_fit <- function(fit) {
  retenues <- unique(coef_table_lasso(fit)$VariableOriginale)
  tibble::tibble(N_variables_entree = length(.vars_candidates(fit)), N_variables_sortie = length(retenues),
                 Variables_sortie = if (length(retenues) > 0) paste(retenues, collapse = ", ") else "(aucune)")
}

#' @export
.effect_variables.lasso_cv_fit <- function(fit, top_n) utils::head(unique(coef_table_lasso(fit)$VariableOriginale), top_n)

#' @export
.print_details.lasso_cv_fit <- function(x) {
  cv <- x$final_cv
  bord <- isTRUE(all.equal(cv$lambda.min, min(cv$lambda))) || isTRUE(all.equal(cv$lambda.min, max(cv$lambda)))
  cat(sprintf("Lambda final  : %.4g  (grille testee : %.4g a %.4g)%s\n", x$final_lambda, min(cv$lambda), max(cv$lambda),
              if (bord) "  \u26a0 AU BORD de la grille -> l'elargir (argument `lambda`)" else ""))
}

#' @export
.print_body.lasso_cv_fit <- function(x) {
  cd <- coef_table_lasso(x)
  cat(sprintf("\nModele final : %d / %d variable(s) retenue(s) (coefficient \u2260 0)%s\n",
              nrow(cd), length(x$x_names), if (identical(x$family, "binomial")) " -- log-odds" else ""))
  if (nrow(cd) == 0) cat("Aucune variable retenue -- modele reduit a l'intercept (cf. plot_hyperparam_path()).\n")
  else print(dplyr::select(cd, Variable, Coefficient), n = 20)
  invisible(NULL)
}

#' @export
.importance_panel.lasso_cv_fit <- function(fit, labels = NULL) {
  plot_collinearity(check_collinearity_lasso(fit$data, fit$formula, fit$id_col, fit = fit, selected_only = TRUE),
                    labels = labels)
}

# ── Lambda ---------------------------------------------------------------------

#' @export
check_hyperparam_stability.lasso_cv_fit <- function(fit) {
  cv <- fit$final_cv
  structure(list(final_at_edge = isTRUE(all.equal(cv$lambda.min, min(cv$lambda))) || isTRUE(all.equal(cv$lambda.min, max(cv$lambda))),
                 pct_folds_edge = 100 * mean(fit$lambda_traj$at_edge, na.rm = TRUE)),
            class = "check_lambda_range_lasso")
}

#' Lambda.min au bord de sa grille ? (modèle final et folds)
#' @param fit objet "lasso_cv_fit"
#' @export
check_lambda_range_lasso <- function(fit) check_hyperparam_stability(fit)

#' @export
print.check_lambda_range_lasso <- function(x, ...) {
  cat(sprintf("Modele final : lambda.min %s au bord de sa grille testee.\n", if (x$final_at_edge) "EST" else "n'est PAS"))
  cat(sprintf("%.1f%% des folds ont un lambda.min au bord de LEUR grille.\n", x$pct_folds_edge))
  if (x$final_at_edge) {
    cat("\u26a0 MODELE FINAL au bord : elargir `lambda` (ex. lambda = exp(seq(-10, 6, length.out = 100))).\n")
  } else if (x$pct_folds_edge > 30) {
    cat("\u2139 Proportion elevee de folds au bord (modele final OK) : marge a elargir, sans urgence.\n")
  } else {
    cat("\u2713 Grille correctement calee.\n")
  }
  invisible(x)
}

#' @export
plot_hyperparam_path.lasso_cv_fit <- function(fit) {
  cv <- fit$final_cv
  bord <- check_hyperparam_stability(fit)$final_at_edge
  df <- tibble::tibble(log_lambda = log(cv$lambda), cvm = cv$cvm, cvlo = cv$cvlo, cvup = cv$cvup)
  yr <- range(df$cvup, df$cvlo); y_top <- yr[2] - 0.04 * diff(yr)
  ggplot2::ggplot(df, ggplot2::aes(x = log_lambda, y = cvm)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = cvlo, ymax = cvup), colour = "grey80", width = 0) +
    ggplot2::geom_point(colour = "#D95F02", size = 1.6) +
    ggplot2::geom_vline(xintercept = log(cv$lambda.min), colour = "#2C3E50", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = log(cv$lambda.1se), colour = "#7B3294", linewidth = 0.8, linetype = "dotdash") +
    ggplot2::annotate("text", x = log(cv$lambda.min), y = y_top, hjust = -0.05, size = 3, colour = "#2C3E50",
                      label = sprintf("lambda = %.4g\nlog(lambda) = %.2f", cv$lambda.min, log(cv$lambda.min))) +
    ggplot2::labs(title = "Chemin de validation crois\u00e9e du lambda (mod\u00e8le final)",
                  subtitle = if (bord) "\u26a0 lambda.min au bord de la grille : \u00e9largir `lambda`" else
                    "Trait plein = lambda.min (pr\u00e9dictions) \u00b7 Tiret-point violet = lambda.1se",
                  x = "log(lambda)", y = sprintf("%s (CV interne)", cv$name %||% "Erreur")) +
    .theme_cv(12, 9, if (bord) "#D95F02" else "grey45")
}

#' @rdname check_lambda_range_lasso
#' @export
plot_lambda_path_lasso <- function(fit) plot_hyperparam_path(fit)

#' @export
plot_hyperparam_stability.lasso_cv_fit <- function(fit) {
  lt <- fit$lambda_traj
  ggplot2::ggplot(lt, ggplot2::aes(x = log(lambda_min))) +
    ggplot2::geom_histogram(bins = 20, fill = "#2C7FB8", alpha = 0.7, colour = "white") +
    ggplot2::geom_vline(xintercept = log(stats::median(lt$lambda_min)), colour = "#D95F02", linetype = "dashed", linewidth = 0.8) +
    ggplot2::labs(title = "Stabilit\u00e9 du lambda optimal",
                  subtitle = sprintf("M\u00e9diane = %.3g \u00b7 distribution large = p\u00e9nalisation instable selon l'\u00e9chantillon", stats::median(lt$lambda_min)),
                  x = "log(lambda.min)", y = "Nombre de folds") +
    .theme_cv(12)
}

#' @rdname check_lambda_range_lasso
#' @export
plot_lambda_stability <- function(fit) plot_hyperparam_stability(fit)

# ── Coefficients ----------------------------------------------------------------

#' @rdname plot_variable_importance
#' @param colors couleurs c(positif =, negatif =) (Lasso)
#' @export
plot_variable_importance.lasso_cv_fit <- function(fit, labels = NULL, top_k = 15, top_n = 30,
                                                  colors = c(positif = "#D95F02", negatif = "#2C7FB8"), ...) {
  full <- fit$fold_coefs %>%
    dplyr::group_by(Variable) %>%
    dplyr::summarise(n_selected = sum(Coefficient != 0), pct_selected = 100 * sum(Coefficient != 0) / fit$n_folds_valides,
                     mean_coef_sel = if (any(Coefficient != 0)) mean(Coefficient[Coefficient != 0]) else NA_real_, .groups = "drop") %>%
    dplyr::filter(pct_selected > 0) %>%
    dplyr::mutate(Label = .label_lookup(Variable, labels)) %>%
    dplyr::arrange(dplyr::desc(pct_selected))
  if (nrow(full) == 0) stop("Aucune variable selectionnee dans aucun fold.", call. = FALSE)
  n_jamais <- dplyr::n_distinct(fit$fold_coefs$Variable) - nrow(full)
  shown <- dplyr::slice_head(full, n = top_n); n_tronque <- nrow(full) - nrow(shown)
  ord <- shown %>% dplyr::arrange(pct_selected) %>% dplyr::pull(Label)

  p_freq <- ggplot2::ggplot(shown, ggplot2::aes(x = factor(Label, levels = ord), y = pct_selected, fill = mean_coef_sel > 0)) +
    ggplot2::geom_col(width = 0.7, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 50, colour = "grey40", linetype = "dashed") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(`TRUE` = unname(colors["positif"]), `FALSE` = unname(colors["negatif"])), na.value = "grey70", guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, 100)) +
    ggplot2::labs(title = sprintf("Fr\u00e9quence de s\u00e9lection (%s)", fit$method_label),
                  subtitle = sprintf("%d variable(s) s\u00e9lectionn\u00e9e(s) \u2265 1 fois sur %d folds (%d jamais)%s\nOrange = effet positif en moyenne \u00b7 Bleu = n\u00e9gatif",
                                     nrow(full), fit$n_folds_valides, n_jamais,
                                     if (n_tronque > 0) sprintf(" \u00b7 top %d affich\u00e9es", top_n) else ""),
                  x = NULL, y = "% des folds o\u00f9 la variable est retenue") +
    .theme_cv(12, 8.5)
  p_coef <- fit$fold_coefs %>%
    dplyr::filter(Coefficient != 0, Variable %in% shown$Variable) %>%
    dplyr::mutate(Label = factor(.label_lookup(Variable, labels), levels = ord)) %>%
    ggplot2::ggplot(ggplot2::aes(x = Label, y = Coefficient)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
    ggplot2::geom_boxplot(fill = unname(colors["negatif"]), alpha = 0.4, width = 0.5, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, height = 0, alpha = 0.35, size = 1) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Amplitude du coefficient (folds o\u00f9 s\u00e9lectionn\u00e9e)",
                  subtitle = "Dispersion large = effet instable",
                  x = NULL, y = if (identical(fit$family, "binomial")) "Coefficient (log-odds)" else "Coefficient") +
    .theme_cv(12) + ggplot2::theme(axis.text.y = ggplot2::element_blank())
  p_freq | p_coef
}

#' @rdname plot_variable_importance
#' @export
plot_coefficient_importance <- function(fit, labels = NULL, top_n = 30, colors = c(positif = "#D95F02", negatif = "#2C7FB8")) {
  .stop_if_not_cv_fit(fit, "lasso_cv_fit")
  plot_variable_importance.lasso_cv_fit(fit, labels = labels, top_n = top_n, colors = colors)
}

# ── VIF ------------------------------------------------------------------------

#' Multicolinéarité des prédicteurs (VIF), à la façon de performance::check_collinearity()
#'
#' Calculé sur un lm() avec une réponse arbitraire (le VIF ne dépend que de X).
#' Non calculable si p >= n (précisément le régime où le Lasso est nécessaire),
#' sauf à restreindre aux variables retenues (\code{selected_only = TRUE}).
#' @param data,formula,id_col données et formule
#' @param fit objet "lasso_cv_fit" (requis si selected_only = TRUE)
#' @param selected_only restreindre aux variables du modèle final
#' @return tibble de classe "check_collinearity_lasso"
#' @export
check_collinearity_lasso <- function(data, formula, id_col = "NUM_PAT", fit = NULL, selected_only = FALSE) {
  .require_pkg("car")
  predictors <- .vars_candidates(list(formula = formula, data = data, id_col = id_col))
  if (isTRUE(selected_only)) {
    stopifnot("selected_only = TRUE necessite `fit`" = !is.null(fit))
    predictors <- intersect(predictors, unique(coef_table_lasso(fit)$VariableOriginale))
  }
  df <- tidyr::drop_na(as.data.frame(data)[, predictors, drop = FALSE])
  vide <- function(reason) structure(
    tibble::tibble(Term = character(0), VIF = numeric(0), Tolerance = numeric(0), Category = character(0)),
    class = c("check_collinearity_lasso", "tbl_df", "tbl", "data.frame"),
    not_computable = TRUE, reason = reason, selected_only = selected_only)
  if (length(predictors) < 2) return(vide(sprintf("%d variable(s) -- le VIF en necessite au moins 2", length(predictors))))
  if (length(predictors) >= nrow(df)) return(vide(sprintf("%d predicteurs >= %d observations", length(predictors), nrow(df))))
  set.seed(1)
  df$.vif_dummy_y <- stats::rnorm(nrow(df))
  vif <- car::vif(stats::lm(stats::reformulate(paste0("`", predictors, "`"), response = ".vif_dummy_y"), data = df))
  if (is.matrix(vif)) vif <- vif[, "GVIF^(1/(2*Df))"]^2
  out <- tibble::tibble(Term = names(vif), VIF = as.numeric(vif), Tolerance = 1 / as.numeric(vif)) %>%
    dplyr::mutate(Category = dplyr::case_when(VIF < 5 ~ "low", VIF < 10 ~ "moderate", TRUE ~ "high")) %>%
    dplyr::arrange(VIF)
  structure(out, class = c("check_collinearity_lasso", class(out)), not_computable = FALSE, selected_only = selected_only)
}

#' @export
print.check_collinearity_lasso <- function(x, ...) {
  if (isTRUE(attr(x, "not_computable"))) {
    cat(sprintf("VIF non calculable : %s.\n", attr(x, "reason")))
    return(invisible(x))
  }
  print(tibble::as_tibble(unclass(x)[c("Term", "VIF", "Tolerance", "Category")]))
  invisible(x)
}

#' Lollipop du VIF (vert < 5, bleu 5-10, rouge >= 10)
#' @param x objet "check_collinearity_lasso"
#' @param labels libellés optionnels
#' @param log_scale axe log10 (TRUE)
#' @export
plot_collinearity <- function(x, labels = NULL, log_scale = TRUE) {
  stopifnot(inherits(x, "check_collinearity_lasso"))
  if (isTRUE(attr(x, "not_computable"))) {
    return(ggplot2::ggplot() + ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
             ggplot2::annotate("text", x = 0.5, y = 0.5, size = 4.2, colour = "grey30",
                               label = sprintf("VIF non calculable\n(%s)", attr(x, "reason"))) +
             ggplot2::theme_void() + ggplot2::labs(title = "Multicolin\u00e9arit\u00e9 entre pr\u00e9dicteurs (VIF)"))
  }
  restreint <- isTRUE(attr(x, "selected_only"))
  df <- tibble::as_tibble(unclass(x)[c("Term", "VIF", "Tolerance", "Category")]) %>%
    dplyr::mutate(Label = .label_lookup(Term, labels), Category = factor(Category, levels = c("low", "moderate", "high"))) %>%
    dplyr::arrange(VIF) %>% dplyr::mutate(Label = factor(Label, levels = unique(Label)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Label, y = VIF, colour = Category)) +
    ggplot2::geom_segment(ggplot2::aes(xend = Label, y = 1, yend = VIF), linewidth = 0.9) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::geom_hline(yintercept = 5, colour = "#1b6ca8", linetype = "dashed", linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 10, colour = "#cd201f", linetype = "dashed", linewidth = 0.6) +
    ggplot2::scale_colour_manual(values = c(low = "#3aaf85", moderate = "#1b6ca8", high = "#cd201f"), drop = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Multicolin\u00e9arit\u00e9 entre pr\u00e9dicteurs (VIF)",
                  subtitle = paste0("Vert < 5 \u00b7 Bleu 5-10 \u00b7 Rouge \u2265 10",
                                    if (restreint) "\nRestreint aux variables retenues par le mod\u00e8le final" else ""),
                  x = NULL, y = if (log_scale) "VIF (\u00e9chelle log10)" else "VIF", colour = "Corr\u00e9lation") +
    .theme_cv(12)
  if (log_scale) p <- p + ggplot2::scale_y_log10()
  p
}
