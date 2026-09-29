# =============================================================================
# method_mars.R -- MARS ({earth})
# =============================================================================
# Spécifique : fonctions charnière + élagage ; réglage (degree, nprune) par CV
# interne (ou GCV) ; sélection de variables RÉELLE (méthode creuse, comme le
# Lasso) ; lecture explicite des termes (basis_table_mars()) ; diagnostic de
# complexité (la non-linéarité est-elle utilisée ?).
#
# fold_importance : à l'échelle de la variable d'ORIGINE (max evimp sur ses
# colonnes, Used = au moins une colonne dans le modèle élagué).
# fold_importance_columns : même chose au niveau des colonnes de X.

.fit_earth <- function(X, Y, family, degree, nprune, penalty, nk, earth_args = list()) {
  args <- utils::modifyList(list(x = X, y = Y, degree = as.integer(degree), nprune = as.integer(nprune),
                                 penalty = penalty, nk = as.integer(max(nk, nprune + 1L))), earth_args)
  if (identical(family, "binomial")) args$glm <- list(family = stats::binomial)
  tryCatch(suppressWarnings(do.call(earth::earth, args)), error = function(e) NULL)
}

.predict_mars_matrix <- function(model, X, family) {
  X <- as.matrix(X)
  as.numeric(if (identical(family, "binomial")) stats::predict(model, newdata = X, type = "response") else stats::predict(model, newdata = X))
}

.earth_used_variables <- function(model, var_names) {
  used <- stats::setNames(rep(FALSE, length(var_names)), var_names)
  dirs <- model$dirs; sel <- model$selected.terms
  if (is.null(dirs) || length(sel) == 0) return(used)
  u <- colSums(abs(dirs[sel, , drop = FALSE])) > 0
  if (length(u) == length(var_names)) used[] <- as.logical(u)
  else if (!is.null(colnames(dirs))) used[intersect(colnames(dirs)[u], var_names)] <- TRUE
  used
}

.earth_term_counts <- function(model) {
  sel <- model$selected.terms; dirs <- model$dirs
  if (is.null(dirs) || length(sel) == 0) return(list(n_terms = length(sel), n_hinge_terms = 0L, n_interaction_terms = 0L))
  d <- dirs[sel, , drop = FALSE]
  list(n_terms = length(sel), n_hinge_terms = as.integer(sum(rowSums(abs(d) == 1) > 0)),
       n_interaction_terms = as.integer(sum(rowSums(d != 0) >= 2)))
}

.earth_importance_table <- function(model, var_names, criterion = "gcv") {
  used <- .earth_used_variables(model, var_names)
  score <- stats::setNames(rep(0, length(var_names)), var_names)
  if (any(used)) {
    imp <- tryCatch(unclass(earth::evimp(model, trim = FALSE)), error = function(e) NULL)
    if (is.matrix(imp) && nrow(imp) > 0 && criterion %in% colnames(imp) && !is.null(rownames(imp))) {
      rn <- sub("-unused$", "", rownames(imp)); vals <- as.numeric(imp[, criterion])
      ok <- rn %in% var_names & !is.na(vals)
      score[rn[ok]] <- vals[ok]
    }
    score[used & score <= 0] <- 1e-6
  }
  score[!used] <- 0
  tibble::tibble(Variable = var_names, Importance = unname(score), Used = unname(used))
}

.default_nprune_grid <- function(n_train) {
  mx <- max(5, min(30, floor(n_train / 4)))
  sort(unique(pmax(as.integer(round(mx * c(0.2, 0.35, 0.5, 0.75, 1))), 3L)))
}

.tune_mars <- function(X, Y, family, degree_grid, nprune_grid, penalty, nk, earth_args, tuning, inner_v, seed) {
  nprune_grid <- sort(unique(as.integer(nprune_grid %||% .default_nprune_grid(nrow(X)))))
  degree_grid <- sort(unique(as.integer(degree_grid)))
  nk_eff <- if (is.null(nk)) max(21L, 2L * max(nprune_grid) + 1L) else as.integer(nk)
  cap <- max(nprune_grid)
  res <- if (identical(tuning, "gcv")) {
    purrr::map_dfr(degree_grid, function(d) {
      m <- .fit_earth(X, Y, family, d, cap, penalty, nk_eff, earth_args)
      tibble::tibble(degree = d, nprune = cap, error = if (is.null(m)) NA_real_ else as.numeric(m$gcv))
    })
  } else {
    folds <- .make_inner_folds(Y, inner_v, family, seed)
    purrr::map_dfr(degree_grid, function(d) {
      err <- do.call(cbind, lapply(folds, function(it) {
        Xtr <- X[-it, , drop = FALSE]; Ytr <- Y[-it]; Xte <- X[it, , drop = FALSE]; Yte <- Y[it]
        m_max <- .fit_earth(Xtr, Ytr, family, d, cap, penalty, nk_eff, earth_args)
        if (is.null(m_max)) return(rep(NA_real_, length(nprune_grid)))
        s_max <- length(m_max$selected.terms)
        e_max <- mean((.predict_mars_matrix(m_max, Xte, family) - Yte)^2)
        vapply(nprune_grid, function(k) {
          if (k >= s_max) return(e_max)
          m <- .fit_earth(Xtr, Ytr, family, d, k, penalty, nk_eff, earth_args)
          if (is.null(m)) NA_real_ else mean((.predict_mars_matrix(m, Xte, family) - Yte)^2)
        }, numeric(1))
      }))
      e <- rowMeans(err, na.rm = TRUE); e[is.nan(e)] <- NA_real_
      tibble::tibble(degree = d, nprune = nprune_grid, error = e)
    })
  }
  best <- if (all(is.na(res$error))) res[1, ] else res[which.min(res$error), ]
  list(best_degree = best$degree, best_nprune = best$nprune, grid_results = res, nprune_grid = nprune_grid,
       nk = nk_eff, mode = tuning, inner_error = if (all(is.na(res$error))) NA_real_ else best$error)
}

#' MARS (earth) en validation croisée répétée
#'
#' @inheritParams fit_lasso_repeated_cv
#' @param degree_grid degrés d'interaction testés (1, 2)
#' @param nprune_grid plafonds de termes testés (NULL = selon n)
#' @param penalty pénalité GCV par noeud (3, conservateur)
#' @param nk taille max de la passe avant (NULL = 2 x max(nprune_grid) + 1, >= 21)
#' @param tuning "cv" (CV interne, défaut) ou "gcv" (plus rapide, plus optimiste)
#' @param inner_v folds de la CV interne (5)
#' @param earth_args arguments supplémentaires pour earth::earth()
#' @param importance critère evimp : "gcv", "rss" ou "nsubsets"
#' @return objet c("mars_cv_fit", "cv_fit")
#' @export
fit_mars_repeated_cv <- function(data, formula, id_col = "NUM_PAT", v = 5, repeats = 20, degree_grid = c(1, 2),
                                 nprune_grid = NULL, penalty = 3, nk = NULL, tuning = c("cv", "gcv"), inner_v = 5,
                                 earth_args = list(), family = c("gaussian", "binomial"), classification_threshold = 0.5,
                                 seed = 42, cv_repeats = NULL, importance = c("gcv", "rss", "nsubsets"),
                                 parallelisation = FALSE, n_cores = NULL) {
  .require_pkg("earth")
  family <- match.arg(family); tuning <- match.arg(tuning); importance <- match.arg(importance)
  force(degree_grid); force(nprune_grid); force(penalty); force(nk); force(earth_args); force(inner_v)

  prep_full <- .prepare_matrix_data(data, formula, id_col, family, pkg = "earth")
  if (is.null(prep_full)) stop("Jeu de donnees inutilisable (facteur mono-niveau ou reponse binaire mal formee).", call. = FALSE)
  x_names <- colnames(prep_full$X)
  vc <- .vars_candidates(list(formula = formula, data = data, id_col = id_col))
  var_map <- stats::setNames(.map_to_original_variable(x_names, vc), x_names)

  fold_fun <- function(train, test, fold_num) {
    tr <- .prepare_matrix_data(train, formula, id_col, family, pkg = "earth")
    te <- .prepare_matrix_data(test,  formula, id_col, family, pkg = "earth")
    if (is.null(tr) || is.null(te) || !identical(colnames(tr$X), colnames(te$X))) return(NULL)
    tun <- .tune_mars(tr$X, tr$Y, family, degree_grid, nprune_grid, penalty, nk, earth_args, tuning, inner_v, seed + fold_num)
    model <- .fit_earth(tr$X, tr$Y, family, tun$best_degree, tun$best_nprune, penalty, tun$nk, earth_args)
    if (is.null(model)) return(NULL)
    imp <- .earth_importance_table(model, colnames(tr$X), importance)
    ct <- .earth_term_counts(model)
    list(Y = te$Y, pred = .predict_mars_matrix(model, te$X, family), extras = list(
      importance = imp,
      hyperparams = tibble::tibble(degree_used = tun$best_degree, nprune_cap = tun$best_nprune,
                                   nprune_grid_max = max(tun$nprune_grid), n_terms = ct$n_terms,
                                   n_hinge_terms = ct$n_hinge_terms, n_interaction_terms = ct$n_interaction_terms,
                                   n_vars_used = sum(imp$Used), inner_error = tun$inner_error,
                                   at_edge = tun$best_nprune >= max(tun$nprune_grid) && ct$n_terms >= tun$best_nprune)))
  }

  cv <- run_repeated_cv(data, formula, id_col, family, v, repeats, seed, cv_repeats, fold_fun,
                        parallelisation, n_cores, pkgs = "earth")

  fin <- .time_it({
    set.seed(seed)
    tun <- .tune_mars(prep_full$X, prep_full$Y, family, degree_grid, nprune_grid, penalty, nk, earth_args, tuning, inner_v, seed)
    m <- .fit_earth(prep_full$X, prep_full$Y, family, tun$best_degree, tun$best_nprune, penalty, tun$nk, earth_args)
    if (is.null(m)) stop("Echec de l'ajustement du modele final (earth::earth()).", call. = FALSE)
    list(tuning = tun, model = m)
  })
  tun <- fin$value$tuning
  ct <- .earth_term_counts(fin$value$model)
  tun$at_edge <- tun$best_nprune >= max(tun$nprune_grid) && ct$n_terms >= tun$best_nprune

  cols_imp <- cv$extras$importance
  fold_importance <- cols_imp %>%
    dplyr::mutate(Variable = dplyr::coalesce(unname(var_map[Variable]), Variable)) %>%
    dplyr::group_by(Variable, Fold) %>%
    dplyr::summarise(Importance = max(Importance), Used = any(Used), .groups = "drop")

  new_cv_fit(
    cv, data, formula, id_col, family, classification_threshold, v, repeats,
    fold_importance = fold_importance, fold_hyperparams = cv$extras$hyperparams,
    final_model = fin$value$model, final_time_sec = fin$time_sec,
    method = "mars", method_label = "MARS", importance_type = sprintf("evimp (%s)", importance),
    fold_importance_columns = cols_imp, final_tuning = tun, final_degree = tun$best_degree,
    final_nprune = tun$best_nprune, final_n_terms = ct$n_terms, x_names = x_names, var_map = var_map,
    tuning = tuning, penalty = penalty, importance_criterion = importance,
    degree_grid = sort(unique(as.integer(degree_grid))), nprune_grid = tun$nprune_grid,
    class = "mars_cv_fit"
  )
}

#' @export
model_design.mars_cv_fit <- function(fit, ...) {
  p <- .prepare_matrix_data(fit$data, fit$formula, fit$id_col, fit$family, pkg = "earth")
  list(X = p$X, Y = p$Y, var_map = fit$var_map,
       predict = function(newdata) .predict_mars_matrix(fit$final_model, newdata, fit$family))
}

#' Termes (fonctions charnière) du modèle final et leurs coefficients
#'
#' Coefficient exprimé dans l'unité de SA variable (non comparable entre
#' variables) ; log-odds en binomial.
#' @param fit objet "mars_cv_fit"
#' @export
basis_table_mars <- function(fit) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  m <- fit$final_model
  cf <- as.matrix(if (identical(fit$family, "binomial") && !is.null(m$glm.coefficients)) m$glm.coefficients else m$coefficients)
  rn <- rownames(cf) %||% paste0("terme_", seq_len(nrow(cf)))
  tibble::tibble(Term = rn, Coefficient = as.numeric(cf[, 1])) %>%
    dplyr::filter(Term != "(Intercept)") %>% dplyr::arrange(dplyr::desc(abs(Coefficient)))
}

#' @export
importance_table.mars_cv_fit <- function(fit, only_used = TRUE, ...) {
  out <- .earth_importance_table(fit$final_model, fit$x_names, fit$importance_criterion) %>%
    dplyr::mutate(VariableOriginale = unname(fit$var_map[Variable])) %>%
    dplyr::arrange(dplyr::desc(Importance))
  if (isTRUE(only_used)) out <- dplyr::filter(out, Used)
  out
}

#' @export
variables_summary.mars_cv_fit <- function(fit) {
  retenues <- unique(importance_table(fit)$VariableOriginale)
  tibble::tibble(N_variables_entree = length(.vars_candidates(fit)), N_variables_sortie = length(retenues),
                 Variables_sortie = if (length(retenues) > 0) paste(retenues, collapse = ", ") else "(aucune)",
                 n_termes_final = fit$final_n_terms)
}

#' @export
.effect_variables.mars_cv_fit <- function(fit, top_n) utils::head(unique(importance_table(fit)$VariableOriginale), top_n)

#' @export
.importance_panel.mars_cv_fit <- function(fit, labels = NULL) plot_selection_mars(check_selection_mars(fit), labels = labels)

#' @export
.print_details.mars_cv_fit <- function(x) {
  cat(sprintf("Reglage       : %s | degree {%s} | nprune {%s} | penalty = %s\n",
              if (identical(x$tuning, "cv")) "CV interne" else "GCV", paste(x$degree_grid, collapse = ", "),
              paste(x$nprune_grid, collapse = ", "), format(x$penalty)))
  cat(sprintf("Modele final  : degree = %d | plafond nprune = %d | %d terme(s) (intercept compris)%s\n",
              x$final_degree, x$final_nprune, x$final_n_terms,
              if (isTRUE(x$final_tuning$at_edge)) "  \u26a0 PLAFOND ATTEINT -> elargir `nprune_grid`" else ""))
}

#' @export
.print_body.mars_cv_fit <- function(x) {
  bt <- basis_table_mars(x)
  cat(sprintf("\nModele final : %d terme(s) hors intercept, par |coefficient| decroissant%s :\n",
              nrow(bt), if (identical(x$family, "binomial")) " (log-odds)" else ""))
  if (nrow(bt) == 0) cat("Aucun terme retenu -- modele reduit a l'intercept.\n") else print(bt, n = 15)
  invisible(NULL)
}

#' @export
check_hyperparam_stability.mars_cv_fit <- function(fit) {
  hp <- fit$fold_hyperparams
  structure(list(final_at_edge = isTRUE(fit$final_tuning$at_edge), pct_folds_edge = 100 * mean(hp$at_edge, na.rm = TRUE),
                 degree_table = table(degree = hp$degree_used), n_terms_median = stats::median(hp$n_terms)),
            class = "check_hyperparam_stability_mars")
}

#' @export
print.check_hyperparam_stability_mars <- function(x, ...) {
  cat(sprintf("Modele final : le plafond nprune %s sature.\n", if (x$final_at_edge) "EST atteint et" else "n'est PAS"))
  cat(sprintf("%.1f%% des folds saturent le plafond maximal de leur grille | termes (mediane) : %.0f\n",
              x$pct_folds_edge, x$n_terms_median))
  cat("Degre retenu :\n"); print(x$degree_table)
  if (x$final_at_edge) cat("\u26a0 Modele final sature : elargir `nprune_grid` (et `nk`).\n")
  else if (x$pct_folds_edge > 30) cat("\u2139 Proportion elevee de folds satures (modele final OK).\n")
  else cat("\u2713 Plafond correctement cale.\n")
  invisible(x)
}

#' @export
plot_hyperparam_path.mars_cv_fit <- function(fit) {
  df <- fit$final_tuning$grid_results; df$degree_f <- factor(df$degree)
  bord <- isTRUE(fit$final_tuning$at_edge); gcv <- identical(fit$tuning, "gcv")
  ylab <- if (gcv) "GCV" else if (identical(fit$family, "binomial")) "Erreur CV interne (Brier)" else "Erreur CV interne (MSE)"
  p <- if (gcv) {
    ggplot2::ggplot(df, ggplot2::aes(x = degree_f, y = error, fill = degree_f)) +
      ggplot2::geom_col(alpha = 0.8, show.legend = FALSE) + ggplot2::labs(x = "degree", y = ylab)
  } else {
    ggplot2::ggplot(df, ggplot2::aes(x = nprune, y = error, colour = degree_f, group = degree_f)) +
      ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 2.2) +
      ggplot2::geom_point(data = df[df$degree == fit$final_degree & df$nprune == fit$final_nprune, ],
                          shape = 4, size = 5, stroke = 1.6, colour = "#2C3E50") +
      ggplot2::labs(x = "nprune (plafond de termes)", y = ylab, colour = "degree")
  }
  p + ggplot2::labs(title = "R\u00e9glage de MARS (mod\u00e8le final)",
                    subtitle = if (bord) "\u26a0 plafond nprune au max de la grille ET satur\u00e9 : \u00e9largir `nprune_grid`" else
                      "Croix = (degree, nprune) retenu \u00b7 plateau = le GCV a \u00e9lagu\u00e9 sous le plafond") +
    .theme_cv(12, 9, if (bord) "#D95F02" else "grey45")
}

#' @export
plot_hyperparam_stability.mars_cv_fit <- function(fit) {
  hp <- fit$fold_hyperparams
  p1 <- ggplot2::ggplot(hp, ggplot2::aes(x = factor(degree_used))) +
    ggplot2::geom_bar(fill = "#2C7FB8", alpha = 0.75, colour = "white") +
    ggplot2::labs(title = "Degr\u00e9 d'interaction retenu", x = "degree", y = "Nombre de folds") + .theme_cv(12)
  p2 <- ggplot2::ggplot(hp, ggplot2::aes(x = n_terms)) +
    ggplot2::geom_bar(fill = "#2C7FB8", alpha = 0.75, colour = "white") +
    ggplot2::geom_vline(xintercept = stats::median(hp$n_terms), colour = "#D95F02", linetype = "dashed", linewidth = 0.8) +
    ggplot2::labs(title = "Nombre de termes du mod\u00e8le \u00e9lagu\u00e9",
                  subtitle = sprintf("M\u00e9diane = %.0f", stats::median(hp$n_terms)),
                  x = "termes (intercept compris)", y = "Nombre de folds") + .theme_cv(12)
  p1 | p2
}

#' La non-linéarité est-elle réellement utilisée par MARS ?
#'
#' Part des folds avec au moins une charnière / une interaction. Presque
#' aucune charnière = MARS se comporte en modèle linéaire creux.
#' @param fit objet "mars_cv_fit"
#' @export
check_complexity_mars <- function(fit) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  hp <- fit$fold_hyperparams
  s <- tibble::tibble(
    n_folds = nrow(hp), pct_folds_avec_charniere = 100 * mean(hp$n_hinge_terms > 0),
    pct_folds_avec_interaction = 100 * mean(hp$n_interaction_terms > 0),
    termes_mediane = stats::median(hp$n_terms), termes_min = min(hp$n_terms), termes_max = max(hp$n_terms),
    variables_mediane = stats::median(hp$n_vars_used))
  structure(list(summary = s, per_fold = hp), class = "check_complexity_mars")
}

#' @export
print.check_complexity_mars <- function(x, ...) {
  s <- x$summary
  cat(sprintf("Complexite MARS sur %d folds : termes mediane %.0f (min %.0f, max %.0f) | variables mediane %.0f\n",
              s$n_folds, s$termes_mediane, s$termes_min, s$termes_max, s$variables_mediane))
  cat(sprintf("  - %.1f%% des folds avec au moins une CHARNIERE ; %.1f%% avec une INTERACTION\n",
              s$pct_folds_avec_charniere, s$pct_folds_avec_interaction))
  if (s$pct_folds_avec_charniere < 20) cat("\u2139 Tres peu de charnieres : MARS ~ modele lineaire creux (comparer au Lasso).\n")
  else if (s$pct_folds_avec_charniere > 80) cat("\u2139 Charnieres frequentes : effets de seuil/courbure probables (cf. plot_effects()).\n")
  invisible(x)
}

#' Stabilité de la sélection des variables par MARS
#' @param fit objet "mars_cv_fit"
#' @param top_n variables conservées (30)
#' @export
check_selection_mars <- function(fit, top_n = 30) {
  .stop_if_not_cv_fit(fit, "mars_cv_fit")
  fi <- fit$fold_importance; nf <- dplyr::n_distinct(fi$Fold)
  out <- fi %>% dplyr::group_by(Fold) %>% dplyr::mutate(Rank = rank(-Importance, ties.method = "average")) %>% dplyr::ungroup() %>%
    dplyr::group_by(Variable) %>%
    dplyr::summarise(pct_folds_used = 100 * sum(Used) / nf,
                     mean_importance_when_used = if (any(Used)) mean(Importance[Used]) else NA_real_,
                     mean_rank = mean(Rank), rank_sd = stats::sd(Rank), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(pct_folds_used), dplyr::desc(mean_importance_when_used)) %>%
    dplyr::slice_head(n = top_n)
  structure(out, class = c("check_selection_mars", class(out)), n_folds = nf, top_n = top_n)
}

#' @export
print.check_selection_mars <- function(x, ...) {
  cat(sprintf("Stabilite de la selection sur %d folds, top %d variable(s) :\n", attr(x, "n_folds"), attr(x, "top_n")))
  print(knitr::kable(tibble::as_tibble(x), digits = 2))
  cat("\npct_folds_used faible = selection opportuniste (souvent une variable parmi plusieurs correlees).\n")
  invisible(x)
}

#' Lollipop de la fréquence de sélection MARS
#' @param x objet "check_selection_mars"
#' @param labels libellés optionnels
#' @export
plot_selection_mars <- function(x, labels = NULL) {
  stopifnot(inherits(x, "check_selection_mars"))
  df <- tibble::as_tibble(x) %>% dplyr::mutate(Label = .label_lookup(Variable, labels)) %>%
    dplyr::arrange(pct_folds_used, mean_importance_when_used) %>% dplyr::mutate(Label = factor(Label, levels = unique(Label)))
  ggplot2::ggplot(df, ggplot2::aes(x = Label, y = pct_folds_used, colour = mean_importance_when_used)) +
    ggplot2::geom_hline(yintercept = 50, colour = "grey60", linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(xend = Label, y = 0, yend = pct_folds_used), linewidth = 0.9) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::scale_colour_gradient(low = "#9ecae1", high = "#08519c", na.value = "grey70", name = "Importance moyenne\n(quand retenue)") +
    ggplot2::scale_y_continuous(limits = c(0, 100)) + ggplot2::coord_flip() +
    ggplot2::labs(title = "Fr\u00e9quence de s\u00e9lection des variables par MARS",
                  subtitle = "% des folds o\u00f9 la variable figure dans le mod\u00e8le \u00e9lagu\u00e9", x = NULL, y = "% des folds") +
    .theme_cv(12, 8.5)
}

#' @rdname plot_variable_importance
#' @export
plot_variable_importance.mars_cv_fit <- function(fit, labels = NULL, top_k = 15, top_n = 30, ...) {
  fi <- fit$fold_importance; nf <- dplyr::n_distinct(fi$Fold)
  full <- fi %>% dplyr::group_by(Variable) %>%
    dplyr::summarise(n_used = sum(Used), pct_used = 100 * sum(Used) / nf,
                     mean_importance_when_used = if (any(Used)) mean(Importance[Used]) else NA_real_, .groups = "drop") %>%
    dplyr::filter(n_used > 0) %>% dplyr::mutate(Label = .label_lookup(Variable, labels)) %>%
    dplyr::arrange(dplyr::desc(pct_used), dplyr::desc(mean_importance_when_used))
  if (nrow(full) == 0) stop("Aucune variable retenue dans aucun fold (cf. check_complexity_mars()).", call. = FALSE)
  shown <- dplyr::slice_head(full, n = top_n); ord <- shown %>% dplyr::arrange(pct_used) %>% dplyr::pull(Label)
  p_freq <- ggplot2::ggplot(shown, ggplot2::aes(x = factor(Label, levels = ord), y = pct_used)) +
    ggplot2::geom_col(width = 0.7, fill = "#2C7FB8") +
    ggplot2::geom_hline(yintercept = 50, colour = "grey40", linetype = "dashed") +
    ggplot2::coord_flip() + ggplot2::scale_y_continuous(limits = c(0, 100)) +
    ggplot2::labs(title = "Fr\u00e9quence de s\u00e9lection par MARS",
                  subtitle = sprintf("%d variable(s) retenue(s) au moins une fois sur %d folds", nrow(full), nf),
                  x = NULL, y = "% des folds o\u00f9 la variable est dans le mod\u00e8le") + .theme_cv(12, 8.5)
  p_imp <- fi %>% dplyr::filter(Used, Variable %in% shown$Variable) %>%
    dplyr::mutate(Label = factor(.label_lookup(Variable, labels), levels = ord)) %>%
    ggplot2::ggplot(ggplot2::aes(x = Label, y = Importance)) +
    ggplot2::geom_boxplot(fill = "#2C7FB8", alpha = 0.4, width = 0.5, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, height = 0, alpha = 0.25, size = 1) + ggplot2::coord_flip() +
    ggplot2::labs(title = "Importance quand la variable est retenue",
                  subtitle = "Normalis\u00e9e \u00e0 100 pour la meilleure variable de chaque fold (comparable en rang)",
                  x = NULL, y = sprintf("Importance evimp (%s)", fit$importance_criterion)) +
    .theme_cv(12, 8.5) + ggplot2::theme(axis.text.y = ggplot2::element_blank())
  p_freq | p_imp
}
