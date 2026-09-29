# =============================================================================
# method_xgb.R -- gradient boosting ({xgboost})
# =============================================================================
# Spécifique : réglage (max_depth, eta, nrounds) par CV interne (on entraîne une
# fois avec le max de tours puis on prédit avec les k premiers arbres) ; deux
# importances : permutation sur le fold de test (défaut, honnête) et gain natif
# (entraînement) ; NA gérés nativement dans les prédicteurs ; nthread = 1 (le
# parallélisme se fait sur les folds). Toute l'API {xgboost} est isolée ici.

.xgb_env <- new.env(parent = emptyenv())

#' Convention de `iterationrange` de la version installée (détectée une fois)
#' @noRd
.xgb_iter_offset <- function() {
  if (!is.null(.xgb_env$offset)) return(.xgb_env$offset)
  set.seed(1)
  Xp <- matrix(stats::rnorm(60), 30, 2, dimnames = list(NULL, c("a", "b")))
  m <- xgboost::xgb.train(params = list(objective = "reg:squarederror", eta = 0.5, max_depth = 2, nthread = 1),
                          data = xgboost::xgb.DMatrix(data = Xp, label = Xp[, 1] + stats::rnorm(30, 0, 0.1)),
                          nrounds = 5, verbose = 0)
  full <- as.numeric(stats::predict(m, Xp))
  essai <- function(off) tryCatch(as.numeric(stats::predict(m, Xp, iterationrange = c(1L, 5L + off))), error = function(e) NULL)
  a <- essai(0L); b <- essai(1L)
  off <- if (!is.null(a) && isTRUE(all.equal(a, full))) 0L else if (!is.null(b) && isTRUE(all.equal(b, full))) 1L else {
    warning("Convention de `iterationrange` non detectee : intervalle ferme suppose.", call. = FALSE); 0L
  }
  .xgb_env$offset <- off
  off
}

.default_xgb_params <- function(xgb_params = list()) {
  utils::modifyList(list(min_child_weight = 1, subsample = 0.8, colsample_bytree = 1, lambda = 1, gamma = 0), xgb_params)
}

.fit_xgb <- function(X, Y, family, params, nrounds, seed = 1) {
  p <- params
  p$objective <- if (identical(family, "binomial")) "binary:logistic" else "reg:squarederror"
  p$seed <- as.integer(seed); p$nthread <- 1L
  tryCatch(xgboost::xgb.train(params = p, data = xgboost::xgb.DMatrix(data = as.matrix(X), label = Y),
                              nrounds = as.integer(nrounds), verbose = 0),
           error = function(e) { .xgb_env$last_error <- conditionMessage(e); NULL })
}

.predict_xgb_matrix <- function(model, X, family, nrounds = NULL) {
  X <- as.matrix(X)
  as.numeric(if (is.null(nrounds)) stats::predict(model, X) else
    stats::predict(model, X, iterationrange = c(1L, as.integer(nrounds) + .xgb_iter_offset())))
}

.xgb_gain_by_variable <- function(model, x_names, var_map) {
  gain <- stats::setNames(rep(0, length(x_names)), x_names)
  imp <- tryCatch(as.data.frame(xgboost::xgb.importance(feature_names = x_names, model = model)), error = function(e) NULL)
  if (!is.null(imp) && nrow(imp) > 0 && all(c("Feature", "Gain") %in% names(imp))) {
    ok <- imp$Feature %in% x_names
    gain[imp$Feature[ok]] <- imp$Gain[ok]
  }
  tibble::tibble(Gain = unname(gain), Variable = unname(var_map[x_names])) %>%
    dplyr::group_by(Variable) %>% dplyr::summarise(Gain = sum(Gain), .groups = "drop")
}

.tune_xgb <- function(X, Y, family, depth_grid, eta_grid, nrounds_grid, params, inner_v, seed) {
  depth_grid <- sort(unique(as.integer(depth_grid))); eta_grid <- sort(unique(eta_grid))
  nrounds_grid <- sort(unique(as.integer(nrounds_grid))); nmax <- max(nrounds_grid)
  folds <- .make_inner_folds(Y, inner_v, family, seed)
  combos <- expand.grid(max_depth = depth_grid, eta = eta_grid, KEEP.OUT.ATTRS = FALSE)
  res <- purrr::map_dfr(seq_len(nrow(combos)), function(g) {
    pr <- params; pr$max_depth <- combos$max_depth[g]; pr$eta <- combos$eta[g]
    err <- do.call(cbind, lapply(folds, function(it) {
      m <- .fit_xgb(X[-it, , drop = FALSE], Y[-it], family, pr, nmax, seed)
      if (is.null(m)) return(rep(NA_real_, length(nrounds_grid)))
      vapply(nrounds_grid, function(k) mean((.predict_xgb_matrix(m, X[it, , drop = FALSE], family, k) - Y[it])^2), numeric(1))
    }))
    e <- rowMeans(err, na.rm = TRUE); e[is.nan(e)] <- NA_real_
    tibble::tibble(max_depth = combos$max_depth[g], eta = combos$eta[g], nrounds = nrounds_grid, error = e)
  })
  b <- if (all(is.na(res$error))) 1L else which.min(res$error)
  list(best_depth = res$max_depth[b], best_eta = res$eta[b], best_nrounds = res$nrounds[b], grid_results = res,
       depth_grid = depth_grid, eta_grid = eta_grid, nrounds_grid = nrounds_grid, inner_error = res$error[b])
}

#' XGBoost en validation croisée répétée
#'
#' @inheritParams fit_lasso_repeated_cv
#' @param max_depth_grid profondeurs testées (1, 2, 3 ; 1 = modèle additif)
#' @param eta_grid taux d'apprentissage testés (0.03, 0.1)
#' @param nrounds_grid nombres d'arbres testés (NULL = 10 à 300)
#' @param xgb_params paramètres xgboost supplémentaires (défaut prudent petit n)
#' @param inner_v folds de la CV interne de réglage (5)
#' @param importance "permutation" (fold de test, défaut) ou "gain" (natif)
#' @param n_perm permutations par variable (20)
#' @param perm_groups liste NOMMÉE de groupes de variables permutés d'un bloc
#' @return objet c("xgb_cv_fit", "cv_fit")
#' @export
fit_xgb_repeated_cv <- function(data, formula, id_col = "NUM_PAT", v = 5, repeats = 20,
                                max_depth_grid = c(1, 2, 3), eta_grid = c(0.03, 0.1), nrounds_grid = NULL,
                                xgb_params = list(), inner_v = 5, importance = c("permutation", "gain"),
                                n_perm = 20, perm_groups = NULL, family = c("gaussian", "binomial"),
                                classification_threshold = 0.5, seed = 42, cv_repeats = NULL,
                                parallelisation = FALSE, n_cores = NULL) {
  .require_pkg("xgboost")
  family <- match.arg(family); importance <- match.arg(importance)
  force(max_depth_grid); force(eta_grid); force(xgb_params); force(inner_v); force(n_perm); force(perm_groups)
  if (is.null(nrounds_grid)) nrounds_grid <- c(10L, 20L, 40L, 80L, 120L, 200L, 300L)
  params <- .default_xgb_params(xgb_params)

  prep_full <- .prepare_matrix_data(data, formula, id_col, family, allow_na_x = TRUE, pkg = "xgboost")
  if (is.null(prep_full)) stop("Jeu de donnees inutilisable (facteur mono-niveau ou reponse binaire mal formee).", call. = FALSE)
  x_names <- colnames(prep_full$X)
  vc <- .vars_candidates(list(formula = formula, data = data, id_col = id_col))
  var_map <- stats::setNames(.map_to_original_variable(x_names, vc), x_names)
  groups_cols <- .resolve_perm_groups(perm_groups, var_map)

  fold_fun <- function(train, test, fold_num) {
    tr <- .prepare_matrix_data(train, formula, id_col, family, allow_na_x = TRUE, pkg = "xgboost")
    te <- .prepare_matrix_data(test,  formula, id_col, family, allow_na_x = TRUE, pkg = "xgboost")
    if (is.null(tr) || is.null(te) || !identical(colnames(tr$X), colnames(te$X))) return(NULL)
    tun <- .tune_xgb(tr$X, tr$Y, family, max_depth_grid, eta_grid, nrounds_grid, params, inner_v, seed + fold_num)
    pb <- params; pb$max_depth <- tun$best_depth; pb$eta <- tun$best_eta
    model <- .fit_xgb(tr$X, tr$Y, family, pb, tun$best_nrounds, seed + fold_num)
    if (is.null(model)) return(NULL)
    perm <- .perm_importance(function(Xn) .predict_xgb_matrix(model, Xn, family), te$X, te$Y, family,
                             var_map, groups_cols, n_perm, seed + fold_num)
    gain <- .xgb_gain_by_variable(model, x_names, var_map)
    list(Y = te$Y, pred = .predict_xgb_matrix(model, te$X, family), extras = list(
      perm = perm$vars, gain = dplyr::transmute(gain, Variable, Importance = Gain), group_importance = perm$groups,
      hyperparams = tibble::tibble(max_depth_used = tun$best_depth, eta_used = tun$best_eta, nrounds_used = tun$best_nrounds,
                                   inner_error = tun$inner_error, n_vars_used = sum(gain$Gain > 0),
                                   at_edge_nrounds = tun$best_nrounds >= max(tun$nrounds_grid),
                                   at_edge_depth = length(tun$depth_grid) > 1 && tun$best_depth >= max(tun$depth_grid))))
  }

  cv <- run_repeated_cv(data, formula, id_col, family, v, repeats, seed, cv_repeats, fold_fun,
                        parallelisation, n_cores, pkgs = "xgboost")

  fin <- .time_it({
    tun <- .tune_xgb(prep_full$X, prep_full$Y, family, max_depth_grid, eta_grid, nrounds_grid, params, inner_v, seed)
    fp <- params; fp$max_depth <- tun$best_depth; fp$eta <- tun$best_eta
    m <- .fit_xgb(prep_full$X, prep_full$Y, family, fp, tun$best_nrounds, seed)
    if (is.null(m)) stop("Echec de l'ajustement du modele final (xgboost) : ", .xgb_env$last_error %||% "", call. = FALSE)
    list(tuning = tun, params = fp, model = m)
  })
  tun <- fin$value$tuning

  new_cv_fit(
    cv, data, formula, id_col, family, classification_threshold, v, repeats,
    fold_importance = if (identical(importance, "gain")) cv$extras$gain else cv$extras$perm,
    fold_hyperparams = cv$extras$hyperparams, final_model = fin$value$model, final_time_sec = fin$time_sec,
    method = "xgb", method_label = "XGBoost",
    importance_type = if (identical(importance, "gain")) "gain natif (entra\u00eenement)" else "permutation (fold de test)",
    fold_perm = cv$extras$perm, fold_gain = cv$extras$gain, fold_group_importance = cv$extras$group_importance,
    final_tuning = tun, final_params = fin$value$params, final_depth = tun$best_depth, final_eta = tun$best_eta,
    final_nrounds = tun$best_nrounds, x_names = x_names, var_map = var_map, importance = importance,
    n_perm = n_perm, perm_groups = perm_groups,
    depth_grid = tun$depth_grid, eta_grid = tun$eta_grid, nrounds_grid = tun$nrounds_grid,
    class = "xgb_cv_fit"
  )
}

#' @export
model_design.xgb_cv_fit <- function(fit, ...) {
  p <- .prepare_matrix_data(fit$data, fit$formula, fit$id_col, fit$family, allow_na_x = TRUE, pkg = "xgboost")
  list(X = p$X, Y = p$Y, var_map = fit$var_map,
       predict = function(newdata) .predict_xgb_matrix(fit$final_model, newdata, fit$family))
}

#' @export
importance_table.xgb_cv_fit <- function(fit, ...) {
  perm <- fit$fold_perm %>% dplyr::group_by(Variable) %>%
    dplyr::summarise(mean_importance = mean(Importance), sd_importance = stats::sd(Importance),
                     pct_folds_positive = 100 * mean(Importance > 0), .groups = "drop")
  gain <- .xgb_gain_by_variable(fit$final_model, fit$x_names, fit$var_map) %>% dplyr::rename(gain_final = Gain)
  dplyr::left_join(perm, gain, by = "Variable") %>% dplyr::arrange(dplyr::desc(mean_importance))
}

#' @export
method_notes.xgb_cv_fit <- function(fit, topic) {
  switch(topic,
    calibration = "Le boosting peut produire des probabilit\u00e9s sur-confiantes (trop proches de 0/1) quand beaucoup de tours sont retenus.",
    NextMethod())
}

#' @export
.print_details.xgb_cv_fit <- function(x) {
  cat(sprintf("Reglage       : CV interne | max_depth {%s} | eta {%s} | nrounds {%s}\n",
              paste(x$depth_grid, collapse = ", "), paste(x$eta_grid, collapse = ", "), paste(x$nrounds_grid, collapse = ", ")))
  cat(sprintf("Modele final  : max_depth = %d | eta = %s | %d arbre(s)\n", x$final_depth, format(x$final_eta), x$final_nrounds))
}

#' @export
check_hyperparam_stability.xgb_cv_fit <- function(fit) {
  hp <- fit$fold_hyperparams
  structure(list(pct_nrounds_edge = 100 * mean(hp$at_edge_nrounds, na.rm = TRUE),
                 pct_depth_edge = 100 * mean(hp$at_edge_depth, na.rm = TRUE),
                 depth_table = table(max_depth = hp$max_depth_used), eta_table = table(eta = hp$eta_used),
                 nrounds_median = stats::median(hp$nrounds_used), n_vars_median = stats::median(hp$n_vars_used)),
            class = "check_hyperparam_stability_xgb")
}

#' @export
print.check_hyperparam_stability_xgb <- function(x, ...) {
  cat(sprintf("%.1f%% des folds retiennent le plus grand nombre de tours ; %.1f%% la plus grande profondeur.\n",
              x$pct_nrounds_edge, x$pct_depth_edge))
  cat(sprintf("Tours (mediane) : %.0f | variables utilisees (mediane) : %.0f\n", x$nrounds_median, x$n_vars_median))
  cat("Profondeur retenue :\n"); print(x$depth_table)
  cat("eta retenu :\n"); print(x$eta_table)
  if (x$pct_nrounds_edge > 30) cat("\u26a0 Souvent au maximum de tours : elargir `nrounds_grid` ou augmenter `eta_grid`.\n")
  else cat("\u2713 Nombre de tours correctement cale.\n")
  if (x$pct_depth_edge > 30) cat("\u2139 Profondeur maximale souvent retenue : envisager `max_depth_grid` plus large.\n")
  invisible(x)
}

#' @export
plot_hyperparam_path.xgb_cv_fit <- function(fit) {
  df <- fit$final_tuning$grid_results
  df$depth_f <- factor(df$max_depth)
  df$eta_lab <- factor(paste0("eta = ", df$eta), levels = paste0("eta = ", sort(unique(df$eta))))
  best <- df[df$max_depth == fit$final_depth & df$eta == fit$final_eta & df$nrounds == fit$final_nrounds, ]
  ggplot2::ggplot(df, ggplot2::aes(x = nrounds, y = error, colour = depth_f, group = depth_f)) +
    ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 2) +
    ggplot2::geom_point(data = best, shape = 4, size = 5, stroke = 1.6, colour = "#2C3E50") +
    ggplot2::facet_wrap(~ eta_lab) +
    ggplot2::labs(title = "R\u00e9glage de XGBoost (mod\u00e8le final)",
                  subtitle = "Croix = combinaison retenue \u00b7 courbe qui remonte = sur-apprentissage quand on ajoute des arbres",
                  x = "Nombre d'arbres", y = if (identical(fit$family, "binomial")) "Erreur CV interne (Brier)" else "Erreur CV interne (MSE)",
                  colour = "max_depth") +
    .theme_cv(12)
}

#' @export
plot_hyperparam_stability.xgb_cv_fit <- function(fit) {
  hp <- fit$fold_hyperparams
  bar <- function(var, title, xlab) {
    ggplot2::ggplot(hp, ggplot2::aes(x = factor(.data[[var]], levels = sort(unique(.data[[var]]))))) +
      ggplot2::geom_bar(fill = "#2C7FB8", alpha = 0.75, colour = "white") +
      ggplot2::labs(title = title, x = xlab, y = "Nombre de folds") + .theme_cv(12)
  }
  bar("max_depth_used", "Profondeur retenue", "max_depth") | bar("eta_used", "eta retenu", "eta") |
    bar("nrounds_used", "Nombre d'arbres retenu", "nrounds")
}
