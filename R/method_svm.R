# =============================================================================
# method_svm.R -- SVM à noyau ({e1071})
# =============================================================================
# Spécifique : standardisation interne (X, et Y en régression) mémorisée dans
# le modèle ; réglage (C, gamma) par CV interne ; probabilités de Platt en
# binomial ; importance par permutation sur le fold de test (+ groupes).

.svm_fit_scaler <- function(X) {
  sdv <- apply(X, 2, stats::sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  list(center = colMeans(X), scale = sdv)
}
.svm_apply_scaler <- function(X, sc) sweep(sweep(X, 2, sc$center, "-"), 2, sc$scale, "/")

.fit_svm <- function(X, Y, family, cost, gamma, epsilon, kernel, svm_args = list()) {
  X <- as.matrix(X); sc <- .svm_fit_scaler(X)
  args <- list(x = .svm_apply_scaler(X, sc), kernel = kernel, cost = cost, scale = FALSE)
  if (!identical(kernel, "linear")) args$gamma <- gamma
  y_center <- 0; y_scale <- 1
  if (identical(family, "binomial")) {
    args$y <- factor(Y, levels = c(0, 1)); args$type <- "C-classification"; args$probability <- TRUE
  } else {
    y_center <- mean(Y); y_scale <- stats::sd(Y); if (!is.finite(y_scale) || y_scale == 0) y_scale <- 1
    args$y <- (Y - y_center) / y_scale; args$type <- "eps-regression"; args$epsilon <- epsilon
  }
  m <- tryCatch(suppressWarnings(do.call(e1071::svm, utils::modifyList(args, svm_args))), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  m$x_scaler <- sc; m$y_center <- y_center; m$y_scale <- y_scale; m$x_names <- colnames(X)
  m
}

.predict_svm_matrix <- function(model, X, family) {
  X <- as.matrix(X)
  if (!is.null(model$x_names)) X <- X[, model$x_names, drop = FALSE]
  Xs <- .svm_apply_scaler(X, model$x_scaler)
  if (identical(family, "binomial")) {
    as.numeric(attr(stats::predict(model, Xs, probability = TRUE), "probabilities")[, "1"])
  } else {
    as.numeric(stats::predict(model, Xs)) * model$y_scale + model$y_center
  }
}

.tune_svm <- function(X, Y, family, kernel, cost_grid, gamma_grid, epsilon, svm_args, inner_v, seed) {
  cost_grid <- sort(unique(cost_grid %||% c(0.25, 1, 4, 16)))
  gamma_grid <- if (identical(kernel, "linear")) NA_real_ else sort(unique(gamma_grid %||% (c(0.1, 0.3, 1, 3) / max(1, ncol(X)))))
  grid <- expand.grid(cost = cost_grid, gamma = gamma_grid, KEEP.OUT.ATTRS = FALSE)
  folds <- .make_inner_folds(Y, inner_v, family, seed)
  errs <- vapply(seq_len(nrow(grid)), function(g) {
    mean(vapply(folds, function(it) {
      set.seed(seed)
      m <- .fit_svm(X[-it, , drop = FALSE], Y[-it], family, grid$cost[g], grid$gamma[g], epsilon, kernel, svm_args)
      if (is.null(m)) NA_real_ else mean((.predict_svm_matrix(m, X[it, , drop = FALSE], family) - Y[it])^2)
    }, numeric(1)), na.rm = TRUE)
  }, numeric(1))
  errs[is.nan(errs)] <- NA_real_
  res <- tibble::as_tibble(grid) %>% dplyr::mutate(error = errs)
  b <- if (all(is.na(errs))) 1L else which.min(errs)
  list(best_cost = res$cost[b], best_gamma = res$gamma[b], grid_results = res, cost_grid = cost_grid,
       gamma_grid = gamma_grid, inner_error = res$error[b])
}

#' SVM à noyau (e1071) en validation croisée répétée
#'
#' @inheritParams fit_xgb_repeated_cv
#' @param kernel "radial" (défaut) ou "linear"
#' @param cost_grid grille de C (NULL = 0.25, 1, 4, 16)
#' @param gamma_grid grille de gamma (NULL = c(0.1, 0.3, 1, 3) / p)
#' @param epsilon epsilon de la régression (0.1, sur Y standardisée)
#' @param svm_args arguments supplémentaires pour e1071::svm()
#' @return objet c("svm_cv_fit", "cv_fit")
#' @export
fit_svm_repeated_cv <- function(data, formula, id_col = "NUM_PAT", v = 5, repeats = 20,
                                kernel = c("radial", "linear"), cost_grid = NULL, gamma_grid = NULL,
                                epsilon = 0.1, svm_args = list(), inner_v = 5, n_perm = 20, perm_groups = NULL,
                                family = c("gaussian", "binomial"), classification_threshold = 0.5, seed = 42,
                                cv_repeats = NULL, parallelisation = FALSE, n_cores = NULL) {
  .require_pkg("e1071")
  family <- match.arg(family); kernel <- match.arg(kernel)
  force(cost_grid); force(gamma_grid); force(epsilon); force(svm_args); force(inner_v); force(n_perm); force(perm_groups)

  prep_full <- .prepare_matrix_data(data, formula, id_col, family, pkg = "e1071")
  if (is.null(prep_full)) stop("Jeu de donnees inutilisable (facteur mono-niveau ou reponse binaire mal formee).", call. = FALSE)
  x_names <- colnames(prep_full$X)
  vc <- .vars_candidates(list(formula = formula, data = data, id_col = id_col))
  var_map <- stats::setNames(.map_to_original_variable(x_names, vc), x_names)
  groups_cols <- .resolve_perm_groups(perm_groups, var_map)

  fold_fun <- function(train, test, fold_num) {
    tr <- .prepare_matrix_data(train, formula, id_col, family, pkg = "e1071")
    te <- .prepare_matrix_data(test,  formula, id_col, family, pkg = "e1071")
    if (is.null(tr) || is.null(te) || !identical(colnames(tr$X), colnames(te$X))) return(NULL)
    tun <- .tune_svm(tr$X, tr$Y, family, kernel, cost_grid, gamma_grid, epsilon, svm_args, inner_v, seed + fold_num)
    set.seed(seed + fold_num)
    model <- .fit_svm(tr$X, tr$Y, family, tun$best_cost, tun$best_gamma, epsilon, kernel, svm_args)
    if (is.null(model)) return(NULL)
    perm <- .perm_importance(function(Xn) .predict_svm_matrix(model, Xn, family), te$X, te$Y, family,
                             var_map, groups_cols, n_perm, seed + fold_num)
    lin <- identical(kernel, "linear")
    list(Y = te$Y, pred = .predict_svm_matrix(model, te$X, family), extras = list(
      importance = perm$vars, group_importance = perm$groups,
      hyperparams = tibble::tibble(cost_used = tun$best_cost, gamma_used = tun$best_gamma, inner_error = tun$inner_error,
                                   n_sv = as.integer(model$tot.nSV), pct_sv = 100 * as.integer(model$tot.nSV) / nrow(tr$X),
                                   at_edge_cost = tun$best_cost >= max(tun$cost_grid),
                                   at_edge_gamma = if (lin) FALSE else tun$best_gamma >= max(tun$gamma_grid) || tun$best_gamma <= min(tun$gamma_grid))))
  }

  cv <- run_repeated_cv(data, formula, id_col, family, v, repeats, seed, cv_repeats, fold_fun,
                        parallelisation, n_cores, pkgs = "e1071")

  fin <- .time_it({
    tun <- .tune_svm(prep_full$X, prep_full$Y, family, kernel, cost_grid, gamma_grid, epsilon, svm_args, inner_v, seed)
    set.seed(seed)
    m <- .fit_svm(prep_full$X, prep_full$Y, family, tun$best_cost, tun$best_gamma, epsilon, kernel, svm_args)
    if (is.null(m)) stop("Echec de l'ajustement du modele final (e1071::svm()).", call. = FALSE)
    list(tuning = tun, model = m)
  })
  tun <- fin$value$tuning

  new_cv_fit(
    cv, data, formula, id_col, family, classification_threshold, v, repeats,
    fold_importance = cv$extras$importance, fold_hyperparams = cv$extras$hyperparams,
    final_model = fin$value$model, final_time_sec = fin$time_sec,
    method = "svm", method_label = if (identical(kernel, "linear")) "SVM lin\u00e9aire" else "SVM",
    importance_type = "permutation (fold de test)",
    fold_group_importance = cv$extras$group_importance, final_tuning = tun,
    final_cost = tun$best_cost, final_gamma = tun$best_gamma, x_names = x_names, var_map = var_map,
    kernel = kernel, epsilon = epsilon, n_perm = n_perm, perm_groups = perm_groups,
    cost_grid = tun$cost_grid, gamma_grid = tun$gamma_grid,
    class = "svm_cv_fit"
  )
}

#' @export
model_design.svm_cv_fit <- function(fit, ...) {
  p <- .prepare_matrix_data(fit$data, fit$formula, fit$id_col, fit$family, pkg = "e1071")
  list(X = p$X, Y = p$Y, var_map = fit$var_map,
       predict = function(newdata) .predict_svm_matrix(fit$final_model, newdata, fit$family))
}

#' @export
.print_details.svm_cv_fit <- function(x) {
  lin <- identical(x$kernel, "linear")
  cat(sprintf("Noyau         : %s | C testes {%s}%s\n", x$kernel, paste(format(x$cost_grid), collapse = ", "),
              if (lin) "" else sprintf(" | gamma testes {%s}", paste(signif(x$gamma_grid, 3), collapse = ", "))))
  cat(sprintf("Modele final  : C = %s%s | %d vecteur(s) de support / %d individus\n", format(x$final_cost),
              if (lin) "" else sprintf(" | gamma = %s", signif(x$final_gamma, 3)),
              as.integer(x$final_model$tot.nSV), nrow(x$data)))
}

#' @export
check_hyperparam_stability.svm_cv_fit <- function(fit) {
  hp <- fit$fold_hyperparams
  structure(list(pct_cost_edge = 100 * mean(hp$at_edge_cost, na.rm = TRUE),
                 pct_gamma_edge = 100 * mean(hp$at_edge_gamma, na.rm = TRUE),
                 cost_table = table(C = hp$cost_used), gamma_table = table(gamma = signif(hp$gamma_used, 3)),
                 pct_sv_median = stats::median(hp$pct_sv), kernel = fit$kernel),
            class = "check_hyperparam_stability_svm")
}

#' @export
print.check_hyperparam_stability_svm <- function(x, ...) {
  cat(sprintf("%.1f%% des folds retiennent le plus grand C de la grille.\n", x$pct_cost_edge))
  if (!identical(x$kernel, "linear")) cat(sprintf("%.1f%% des folds retiennent un gamma en bordure de grille.\n", x$pct_gamma_edge))
  cat(sprintf("Part de vecteurs de support (mediane) : %.0f%%\n", x$pct_sv_median))
  cat("C retenu :\n"); print(x$cost_table)
  if (!identical(x$kernel, "linear")) { cat("gamma retenu :\n"); print(x$gamma_table) }
  if (x$pct_cost_edge > 30 || x$pct_gamma_edge > 30) cat("\u26a0 Proportion elevee de folds en bordure : elargir `cost_grid` / `gamma_grid`.\n")
  else cat("\u2713 Grille correctement centree.\n")
  if (x$pct_sv_median > 90) cat("\u26a0 Presque tous les individus sont vecteurs de support : sur-apprentissage probable.\n")
  invisible(x)
}

#' @export
plot_hyperparam_path.svm_cv_fit <- function(fit) {
  df <- fit$final_tuning$grid_results
  ylab <- if (identical(fit$family, "binomial")) "Erreur CV interne (Brier)" else "Erreur CV interne (MSE)"
  if (identical(fit$kernel, "linear")) {
    return(ggplot2::ggplot(df, ggplot2::aes(x = cost, y = error)) +
             ggplot2::geom_line(linewidth = 0.8, colour = "#2C7FB8") + ggplot2::geom_point(size = 2.2, colour = "#2C7FB8") +
             ggplot2::geom_point(data = df[df$cost == fit$final_cost, ], shape = 4, size = 5, stroke = 1.6, colour = "#2C3E50") +
             ggplot2::scale_x_log10() +
             ggplot2::labs(title = "R\u00e9glage du SVM lin\u00e9aire (mod\u00e8le final)", x = "C (\u00e9chelle log)", y = ylab) + .theme_cv(12))
  }
  df$gamma_lab <- factor(signif(df$gamma, 3), levels = sort(unique(signif(df$gamma, 3))))
  df$cost_lab <- factor(df$cost, levels = sort(unique(df$cost)))
  best <- df[df$cost == fit$final_cost & abs(df$gamma - fit$final_gamma) < 1e-12, ]
  ggplot2::ggplot(df, ggplot2::aes(x = cost_lab, y = gamma_lab, fill = error)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", error)), size = 3) +
    ggplot2::geom_point(data = best, shape = 4, size = 7, stroke = 1.8, colour = "#2C3E50") +
    ggplot2::scale_fill_gradient(low = "#e8f5e9", high = "#d95f02", name = ylab, na.value = "grey85") +
    ggplot2::labs(title = "R\u00e9glage du SVM (mod\u00e8le final)",
                  subtitle = "Croix = (C, gamma) retenu \u00b7 une vall\u00e9e large de cases vertes = r\u00e9glage peu sensible, donc fiable",
                  x = "C (co\u00fbt)", y = "gamma") +
    .theme_cv(12) + ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' @export
plot_hyperparam_stability.svm_cv_fit <- function(fit) {
  hp <- fit$fold_hyperparams
  p_c <- ggplot2::ggplot(hp, ggplot2::aes(x = factor(cost_used))) +
    ggplot2::geom_bar(fill = "#2C7FB8", alpha = 0.75, colour = "white") +
    ggplot2::labs(title = "C retenu", x = "C", y = "Nombre de folds") + .theme_cv(12)
  p_sv <- ggplot2::ggplot(hp, ggplot2::aes(x = pct_sv)) +
    ggplot2::geom_histogram(bins = 20, fill = "#2C7FB8", alpha = 0.75, colour = "white") +
    ggplot2::geom_vline(xintercept = stats::median(hp$pct_sv), colour = "#D95F02", linetype = "dashed", linewidth = 0.8) +
    ggplot2::labs(title = "Part de vecteurs de support",
                  subtitle = sprintf("M\u00e9diane = %.0f%% \u00b7 proche de 100%% = le mod\u00e8le m\u00e9morise", stats::median(hp$pct_sv)),
                  x = "% d'individus d'entra\u00eenement", y = "Nombre de folds") + .theme_cv(12)
  if (identical(fit$kernel, "linear")) return(p_c | p_sv)
  p_g <- ggplot2::ggplot(hp, ggplot2::aes(x = factor(signif(gamma_used, 3)))) +
    ggplot2::geom_bar(fill = "#2C7FB8", alpha = 0.75, colour = "white") +
    ggplot2::labs(title = "gamma retenu", x = "gamma", y = "Nombre de folds") + .theme_cv(12)
  p_c | p_g | p_sv
}
