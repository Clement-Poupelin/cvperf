# =============================================================================
# compare.R -- comparaison de modèles ajustés sur les MÊMES folds
# =============================================================================
# PRÉREQUIS : construire rsample::vfold_cv() UNE fois et le passer en
# `cv_repeats` à CHAQUE fit_*_repeated_cv().
#
# p-values INDICATIVES : les folds d'une CV répétée ne sont pas indépendants
# (Nadeau & Bengio, 2003) -> tests anti-conservateurs. Lire l'AMPLEUR
# (diff_mediane, leaderboard) avant les étoiles.
#
# compare_models()  : 2 à 4 fits de la MÊME méthode (ex. 2 formules) + table
#                     "variables en entrée / en sortie" (variables_summary()).
# compare_methods() : 2 à 5 fits de méthodes quelconques + leaderboard.
# Les deux partagent le même moteur (.paired_comparison()).

# ── Adaptateur minimal ---------------------------------------------------------

#' Structure minimale commune pour la comparaison
#'
#' Pour tout "cv_fit", rien à écrire : la méthode \code{cv_fit} lit les champs
#' communs. \code{importance} (Variable/Score/Fold) n'est comparable qu'en RANG
#' entre méthodes. Pour une méthode creuse (Lasso, MARS : colonne \code{Used}
#' dans \code{fold_importance}), Score = 0 si la variable est écartée.
#' @param fit objet "cv_fit" ou "comparable_fit"
#' @param label libellé d'affichage (défaut : \code{fit$method_label})
#' @export
as_comparable_fit <- function(fit, label = NULL) UseMethod("as_comparable_fit")

#' @export
as_comparable_fit.default <- function(fit, label = NULL) {
  stop("as_comparable_fit() ne sait pas traiter la classe '", paste(class(fit), collapse = "/"),
       "'. Tout objet 'cv_fit' est accepte.", call. = FALSE)
}

#' @export
as_comparable_fit.comparable_fit <- function(fit, label = NULL) {
  if (!is.null(label)) fit$label <- label
  fit
}

#' @export
as_comparable_fit.cv_fit <- function(fit, label = NULL) {
  fi <- fit$fold_importance
  sparse <- "Used" %in% names(fi)
  importance <- if (sparse) {
    dplyr::transmute(fi, Variable, Score = dplyr::if_else(Used, pmax(Importance, 1e-6), 0), Fold)
  } else {
    dplyr::transmute(fi, Variable, Score = Importance, Fold)
  }
  structure(
    list(label = label %||% fit$method_label, method = fit$method, sparse = sparse, family = fit$family,
         classification_threshold = fit$classification_threshold, positive_class = fit$positive_class,
         cv_repeats = fit$cv_repeats, fold_predictions = fit$fold_predictions, importance = importance,
         n_variables_candidates = length(.vars_candidates(fit)),
         computational_info = fit$computational_info),
    class = "comparable_fit"
  )
}

#' Liste nommée de "comparable_fit" (labels : noms de liste > `labels` > défaut)
#' @noRd
.standardize_fits <- function(fits, labels = NULL, min_n = 2) {
  stopifnot("fits doit etre une liste de modeles" = is.list(fits) && length(fits) >= min_n)
  if (is.null(labels)) labels <- names(fits)
  cmp <- lapply(seq_along(fits), function(i) {
    lbl <- if (!is.null(labels) && !is.na(labels[i]) && nzchar(labels[i])) labels[i] else NULL
    as_comparable_fit(fits[[i]], label = lbl)
  })
  lab <- make.unique(vapply(cmp, `[[`, "", "label"), sep = " ")
  for (i in seq_along(cmp)) cmp[[i]]$label <- lab[i]
  stats::setNames(cmp, lab)
}

#' Vérifie (best-effort) que les modèles partagent les mêmes folds
#' @noRd
.check_same_folds <- function(fits) {
  n_folds <- vapply(fits, function(f) nrow(f$cv_repeats), integer(1))
  if (length(unique(n_folds)) > 1) {
    warning("Les modeles n'ont pas le meme nombre de folds (", paste(n_folds, collapse = " vs "),
            ") : reconstruire chaque fit avec le MEME objet `cv_repeats`.", call. = FALSE)
    return(invisible(FALSE))
  }
  ref <- lapply(fits[[1]]$cv_repeats$splits, rsample::complement)
  ok <- vapply(fits[-1], function(f) identical(ref, lapply(f$cv_repeats$splits, rsample::complement)), logical(1))
  if (!all(ok)) {
    warning("Les folds ne sont PAS identiques entre les modeles compares : le test apparie suppose que chaque ",
            "'Fold' designe le meme decoupage train/test. Reconstruire les fits avec le MEME `cv_repeats`.",
            call. = FALSE)
  }
  invisible(all(ok))
}

.metric_labels <- function(family) {
  if (identical(family, "binomial")) {
    c(accuracy = "Accuracy", sensitivity = "Sensibilit\u00e9", specificity = "Sp\u00e9cificit\u00e9",
      precision = "Pr\u00e9cision", f1 = "F1", auc = "AUC", brier = "Brier", logloss = "Log-loss")
  } else {
    c(rsq = "R\u00b2", rmse = "RMSE", mae = "MAE", mse = "MSE")
  }
}

.default_metrics <- function(family) {
  if (identical(family, "binomial")) c("accuracy", "auc", "brier", "logloss") else c("rmse", "mse", "rsq", "mae")
}

#' Moteur de comparaison appariée (tests + plots + leaderboard)
#' @noRd
.paired_comparison <- function(comparable, metrics, p_adjust_method, title) {
  familles <- unique(vapply(comparable, `[[`, "", "family"))
  if (length(familles) != 1) stop("Tous les modeles compares doivent partager la meme famille (gaussian ou binomial).", call. = FALSE)
  family <- familles
  labels_final <- names(comparable)
  .check_same_folds(comparable)

  metric_labels <- .metric_labels(family)
  if (is.null(metrics)) metrics <- .default_metrics(family)
  if (!all(metrics %in% names(metric_labels))) {
    stop("Metrique(s) invalide(s) pour family = '", family, "'. Disponibles : ",
         paste(names(metric_labels), collapse = ", "), call. = FALSE)
  }

  perf_list <- lapply(comparable, function(cf) {
    .per_fold_metrics(cf$fold_predictions, family, cf$classification_threshold) %>% dplyr::mutate(Modele = cf$label)
  })
  combined <- dplyr::bind_rows(perf_list) %>% dplyr::mutate(Modele = factor(Modele, levels = labels_final))

  res <- lapply(metrics, function(metric) {
    wide <- purrr::reduce(
      purrr::map2(perf_list, labels_final, ~ dplyr::select(.x, Fold, dplyr::all_of(metric)) %>%
                    dplyr::rename(!!.y := dplyr::all_of(metric))),
      dplyr::inner_join, by = "Fold")
    if (length(labels_final) == 2) {
      gt <- stats::wilcox.test(wide[[labels_final[1]]], wide[[labels_final[2]]], paired = TRUE, exact = FALSE)
      gt_label <- sprintf("Wilcoxon appari\u00e9 (%d folds) : p = %.3f (%s)", nrow(wide), gt$p.value, p_to_stars(gt$p.value))
    } else {
      gt <- stats::friedman.test(as.matrix(wide[, labels_final]))
      gt_label <- sprintf("Friedman (%d folds) : p = %.3f (%s)", nrow(wide), gt$p.value, p_to_stars(gt$p.value))
    }
    pw <- purrr::map_dfr(utils::combn(labels_final, 2, simplify = FALSE), function(pr) {
      x <- wide[[pr[1]]]; y <- wide[[pr[2]]]
      tibble::tibble(Metric = metric, Modele_A = pr[1], Modele_B = pr[2],
                     diff_mediane = stats::median(x - y, na.rm = TRUE),
                     p_value = stats::wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value)
    }) %>% dplyr::mutate(p_adj = stats::p.adjust(p_value, method = p_adjust_method), stars = p_to_stars(p_adj))

    p <- ggplot2::ggplot(combined, ggplot2::aes(x = Modele, y = .data[[metric]], group = Fold)) +
      ggplot2::geom_line(colour = "grey65", alpha = 0.22, linewidth = 0.3) +
      ggplot2::geom_point(colour = "grey45", alpha = 0.35, size = 1.2) +
      ggplot2::stat_summary(ggplot2::aes(group = 1), fun = stats::median, geom = "line", colour = "#D95F02", linewidth = 1.1) +
      ggplot2::stat_summary(ggplot2::aes(group = 1), geom = "errorbar", width = 0.08, colour = "#D95F02", linewidth = 0.6,
                            fun.min = function(z) stats::quantile(z, 0.25), fun.max = function(z) stats::quantile(z, 0.75)) +
      ggplot2::stat_summary(ggplot2::aes(group = 1), fun = stats::median, geom = "point", colour = "#D95F02", size = 3) +
      ggplot2::labs(title = metric_labels[[metric]], subtitle = gt_label, x = NULL,
                    y = sprintf("%s (fold de test)", metric_labels[[metric]])) +
      .theme_cv(12, 8) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
    if (metric == "rsq") p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40")
    list(global_test = gt, pairwise = pw, plot = p)
  })
  names(res) <- metrics

  panel <- purrr::reduce(lapply(res, `[[`, "plot"), `+`) +
    patchwork::plot_layout(ncol = ceiling(sqrt(length(metrics)))) +
    patchwork::plot_annotation(
      title = title,
      subtitle = "Trait fin gris = 1 fold (reli\u00e9 d'un mod\u00e8le \u00e0 l'autre) \u00b7 Point + trait orange = m\u00e9diane et IQR",
      caption = "p-values INDICATIVES : folds non ind\u00e9pendants (Nadeau & Bengio, 2003) -- lire l'ampleur de l'\u00e9cart avant les \u00e9toiles")

  list(fits = comparable, table_comparaison = combined,
       global_test = lapply(res, `[[`, "global_test"),
       pairwise = dplyr::bind_rows(lapply(res, `[[`, "pairwise")),
       plot = panel,
       leaderboard = leaderboard_table(combined, metrics, metric_labels, labels_final))
}

#' Compare 2 à 4 modèles de la MÊME méthode (ex. 2 formules), mêmes folds
#'
#' @param fits liste (idéalement nommée) de 2 à 4 fits de même classe
#' @param labels libellés optionnels
#' @param p_adjust_method correction multi-tests (défaut "BH")
#' @param metrics métriques ; NULL = 4 par défaut selon la famille
#' @return liste(fits, table_comparaison, global_test, pairwise, plot,
#'   leaderboard, variables_summary)
#' @export
compare_models <- function(fits, labels = NULL, p_adjust_method = "BH", metrics = NULL) {
  stopifnot("fits doit etre une liste de 2 a 4 objets cv_fit" =
              is.list(fits) && length(fits) >= 2 && length(fits) <= 4 && all(vapply(fits, inherits, logical(1), "cv_fit")))
  if (length(unique(vapply(fits, function(f) class(f)[1], ""))) > 1) {
    stop("compare_models() compare des fits de la MEME methode ; utiliser compare_methods() pour des methodes differentes.",
         call. = FALSE)
  }
  if (is.null(labels)) labels <- names(fits)
  if (is.null(labels)) labels <- paste0("Modele ", seq_along(fits))
  comparable <- .standardize_fits(fits, labels)
  out <- .paired_comparison(comparable, metrics, p_adjust_method,
                            sprintf("Comparaison de %d mod\u00e8les %s sur les m\u00eames folds", length(fits), fits[[1]]$method_label))
  out$variables_summary <- purrr::map2_dfr(fits, names(comparable), function(f, l) {
    dplyr::bind_cols(tibble::tibble(Modele = l), variables_summary(f))
  })
  out
}

#' Compare 2 à 5 modèles de méthodes quelconques, mêmes folds
#'
#' @inheritParams compare_models
#' @param fits liste (idéalement nommée) de 2 à 5 "cv_fit" / "comparable_fit"
#' @return liste(fits, table_comparaison, global_test, pairwise, plot, leaderboard)
#' @export
compare_methods <- function(fits, labels = NULL, p_adjust_method = "BH", metrics = NULL) {
  comparable <- .standardize_fits(fits, labels)
  stopifnot("compare_methods() accepte 2 a 5 modeles" = length(comparable) >= 2 && length(comparable) <= 5)
  .paired_comparison(comparable, metrics, p_adjust_method,
                     sprintf("Comparaison de %d m\u00e9thode(s) sur les m\u00eames folds : %s",
                             length(comparable), paste(names(comparable), collapse = " vs ")))
}

#' Tableau de synthèse par modèle (médiane / moyenne / sd + rangs)
#'
#' @param combined tibble long (Fold, Modele, métriques) -- \code{$table_comparaison}
#' @param metrics métriques à résumer
#' @param metric_labels libellés (non utilisés dans le calcul, conservés pour compatibilité)
#' @param labels_final ordre des modèles
#' @export
leaderboard_table <- function(combined, metrics, metric_labels = NULL, labels_final = levels(combined$Modele)) {
  erreurs <- c("rmse", "mae", "mse", "brier", "logloss")
  resume <- combined %>%
    dplyr::group_by(Modele) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(metrics),
                                   list(mediane = ~ stats::median(.x, na.rm = TRUE), moyenne = ~ mean(.x, na.rm = TRUE),
                                        sd = ~ stats::sd(.x, na.rm = TRUE)),
                                   .names = "{.col}_{.fn}"), .groups = "drop")
  for (m in metrics) {
    col <- resume[[paste0(m, "_mediane")]]
    resume[[paste0("rang_", m)]] <- if (m %in% erreurs) rank(col) else rank(-col)
  }
  resume %>%
    dplyr::mutate(rang_moyen = rowMeans(dplyr::across(dplyr::starts_with("rang_")))) %>%
    dplyr::arrange(rang_moyen) %>%
    dplyr::mutate(Modele = factor(Modele, levels = labels_final))
}

#' Heatmap du rang par modèle x métrique
#' @param leaderboard sortie de leaderboard_table() / \code{$leaderboard}
#' @param metrics ignoré (déduit des colonnes \code{rang_*}) ; conservé pour compatibilité
#' @export
plot_metric_ranks <- function(leaderboard, metrics = NULL) {
  rl <- leaderboard %>%
    dplyr::select("Modele", dplyr::starts_with("rang_"), -"rang_moyen") %>%
    tidyr::pivot_longer(-"Modele", names_to = "Metric", values_to = "Rang") %>%
    dplyr::mutate(Metric = sub("^rang_", "", Metric))
  ggplot2::ggplot(rl, ggplot2::aes(x = Metric, y = Modele, fill = Rang)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", Rang)), colour = "white", fontface = "bold") +
    ggplot2::scale_fill_gradient(low = "#3aaf85", high = "#cd201f", name = "Rang\n(1 = meilleur)") +
    ggplot2::labs(title = "Rang par m\u00e9thode et par m\u00e9trique", x = NULL, y = NULL) +
    .theme_cv(12)
}

# ── Importance / sélection ------------------------------------------------------

.rank_correlations <- function(rank_wide, labels) {
  purrr::map_dfr(utils::combn(labels, 2, simplify = FALSE), function(pr) {
    w <- tidyr::drop_na(dplyr::select(rank_wide, "Variable", dplyr::all_of(pr)))
    tibble::tibble(Modele_A = pr[1], Modele_B = pr[2], n_variables_communes = nrow(w),
                   rho = if (nrow(w) < 3) NA_real_ else stats::cor(w[[pr[1]]], w[[pr[2]]], method = "spearman"))
  })
}

.bump_chart <- function(df, rank_col, labels, top_n, title, ylab, shape_col = NULL) {
  best <- df %>% dplyr::group_by(Variable) %>%
    dplyr::summarise(best_rank = min(.data[[rank_col]]), .groups = "drop") %>%
    dplyr::arrange(best_rank) %>% dplyr::slice_head(n = top_n) %>% dplyr::pull(Variable)
  pd <- dplyr::filter(df, Variable %in% best)
  p <- ggplot2::ggplot(pd, ggplot2::aes(x = Modele, y = .data[[rank_col]], group = Variable, colour = Variable)) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.8)
  p <- if (is.null(shape_col)) p + ggplot2::geom_point(size = 2.5) else p + ggplot2::geom_point(ggplot2::aes(shape = .data[[shape_col]]), size = 3)
  p + ggrepel::geom_text_repel(data = dplyr::filter(pd, Modele == labels[length(labels)]),
                               ggplot2::aes(label = Variable), size = 3, direction = "y", show.legend = FALSE) +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(title = title, subtitle = sprintf("Top %d variable(s) (meilleur rang toutes m\u00e9thodes confondues) \u00b7 rang 1 = plus importante", top_n),
                  x = NULL, y = ylab) +
    .theme_cv(12, 8.5) + ggplot2::theme(legend.position = if (is.null(shape_col)) "none" else "right") +
    ggplot2::guides(colour = "none")
}

#' Compare le CLASSEMENT des variables entre méthodes (rang moyen sur les folds)
#'
#' Coefficient, importance de permutation et evimp ne sont pas sur la même
#' échelle : seule la comparaison des RANGS a un sens.
#' @param fits liste de fits
#' @param labels libellés optionnels
#' @param top_n variables du bump chart (15)
#' @return liste(table, rank_correlation, plot)
#' @export
compare_importance <- function(fits, labels = NULL, top_n = 15) {
  cmp <- .standardize_fits(fits, labels)
  lab <- names(cmp)
  rank_table <- purrr::map_dfr(cmp, function(cf) {
    cf$importance %>% dplyr::group_by(Fold) %>% dplyr::mutate(Rank = rank(-Score, ties.method = "average")) %>%
      dplyr::ungroup() %>% dplyr::group_by(Variable) %>%
      dplyr::summarise(mean_rank = mean(Rank), .groups = "drop") %>% dplyr::mutate(Modele = cf$label)
  }) %>% dplyr::mutate(Modele = factor(Modele, levels = lab))
  wide <- tidyr::pivot_wider(rank_table, names_from = "Modele", values_from = "mean_rank")
  list(table = rank_table, rank_correlation = .rank_correlations(wide, lab),
       plot = .bump_chart(rank_table, "mean_rank", lab, top_n,
                          "Classement des variables selon la m\u00e9thode (rang moyen sur les folds)",
                          "Rang moyen (1 = plus importante)"))
}

#' Recouvrement des variables jugées "importantes" (indice de Jaccard)
#'
#' "Importante" : méthodes creuses (Lasso, MARS) = retenue (Score != 0) dans
#' au moins \code{min_pct_folds}\% des folds ; méthodes denses (forêt, SVM,
#' XGBoost) = dans le top-\code{top_k} d'au moins \code{min_pct_folds}\% des
#' folds. Avec peu de variables, prendre \code{top_k} ~ la moitié de p.
#' @param fits liste de fits
#' @param labels libellés optionnels
#' @param min_pct_folds seuil de fréquence (50)
#' @param top_k top-K des méthodes denses (15)
#' @param top_k_rf alias de \code{top_k} (compatibilité)
#' @return liste(table, jaccard, plot)
#' @export
compare_variable_selection <- function(fits, labels = NULL, min_pct_folds = 50, top_k = 15, top_k_rf = NULL) {
  if (!is.null(top_k_rf)) top_k <- top_k_rf
  cmp <- .standardize_fits(fits, labels)
  lab <- names(cmp)
  sel <- lapply(cmp, function(cf) {
    nf <- dplyr::n_distinct(cf$importance$Fold)
    imp <- cf$importance
    if (!isTRUE(cf$sparse)) {
      imp <- imp %>% dplyr::group_by(Fold) %>% dplyr::mutate(Score = as.numeric(rank(-Score, ties.method = "average") <= top_k)) %>% dplyr::ungroup()
    }
    imp %>% dplyr::group_by(Variable) %>% dplyr::summarise(pct = 100 * sum(Score != 0) / nf, .groups = "drop") %>%
      dplyr::filter(.data$pct >= min_pct_folds) %>% dplyr::pull(Variable)
  })
  toutes <- unique(unlist(sel))
  presence <- tibble::tibble(Variable = toutes)
  for (l in lab) presence[[l]] <- toutes %in% sel[[l]]
  presence <- presence %>% dplyr::mutate(n_methodes = rowSums(dplyr::across(dplyr::all_of(lab)))) %>%
    dplyr::arrange(dplyr::desc(n_methodes))
  jaccard <- purrr::map_dfr(utils::combn(lab, 2, simplify = FALSE), function(pr) {
    a <- sel[[pr[1]]]; b <- sel[[pr[2]]]; u <- length(union(a, b))
    tibble::tibble(Modele_A = pr[1], Modele_B = pr[2], n_A = length(a), n_B = length(b),
                   n_intersection = length(intersect(a, b)), n_union = u,
                   jaccard = if (u > 0) length(intersect(a, b)) / u else NA_real_)
  })
  p <- presence %>% dplyr::mutate(Variable = factor(Variable, levels = rev(Variable))) %>%
    ggplot2::ggplot(ggplot2::aes(x = Variable, y = n_methodes, fill = factor(n_methodes))) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::scale_fill_brewer(palette = "Blues", name = "Nb m\u00e9thodes") +
    ggplot2::labs(title = "Variables jug\u00e9es importantes, par nombre de m\u00e9thodes d'accord",
                  subtitle = sprintf("Retenue dans \u2265 %d%% des folds (m\u00e9thodes creuses) ou dans le top-%d (m\u00e9thodes denses)", min_pct_folds, top_k),
                  x = NULL, y = "Nombre de m\u00e9thodes") +
    .theme_cv(12, 8.5)
  list(table = presence, jaccard = jaccard, plot = p)
}

#' Compare l'importance SHAP de plusieurs modèles (même algorithme partout)
#'
#' Comparable en RANG et en MAGNITUDE. Si un Lasso figure parmi les fits,
#' ajoute \code{Lasso_selectionnee} (coefficient non nul dans au moins
#' \code{lasso_selection_threshold}\% des folds).
#' @param fits liste (nommée) de 2 à 5 "cv_fit" bruts
#' @param labels libellés optionnels
#' @param nsim,sample_size,seed cf. check_shap()
#' @param top_n variables du bump chart
#' @param lasso_selection_threshold seuil en \% (50)
#' @param parallelisation,n_cores par modèle si assez de coeurs, sinon par observation
#' @param methods ignoré (compatibilité)
#' @return liste(shap, table, rank_correlation, plot, time_sec)
#' @export
compare_shap <- function(fits, labels = NULL, nsim = 50, sample_size = NULL, seed = 42, top_n = 15,
                         lasso_selection_threshold = 50, parallelisation = FALSE, n_cores = NULL, methods = NULL) {
  stopifnot("compare_shap() attend une liste de 2 a 5 cv_fit" =
              is.list(fits) && length(fits) >= 2 && length(fits) <= 5 && all(vapply(fits, inherits, logical(1), "cv_fit")))
  if (is.null(labels)) labels <- names(fits)
  if (is.null(labels) || any(is.na(labels) | !nzchar(labels))) labels <- vapply(fits, `[[`, "", "method_label")
  labels <- make.unique(labels, sep = " ")
  names(fits) <- labels

  n_cores_used <- .resolve_n_cores(parallelisation, n_cores, length(fits))
  if (.Platform$OS.type == "windows") n_cores_used <- 1L
  t0 <- Sys.time()
  shap_objs <- if (n_cores_used >= length(fits) && n_cores_used > 1L) {
    parallel::mclapply(fits, check_shap, nsim = nsim, sample_size = sample_size, seed = seed, mc.cores = n_cores_used)
  } else {
    lapply(labels, function(l) {
      message(sprintf("compare_shap() : '%s' (nsim = %d)...", l, nsim))
      check_shap(fits[[l]], nsim = nsim, sample_size = sample_size, seed = seed,
                 parallelisation = isTRUE(parallelisation), n_cores = n_cores)
    })
  }
  names(shap_objs) <- labels
  time_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  table <- purrr::map2_dfr(shap_objs, labels, ~ dplyr::mutate(.x$summary, Modele = .y)) %>%
    dplyr::mutate(Modele = factor(Modele, levels = labels))
  lasso_i <- which(vapply(fits, inherits, logical(1), "lasso_cv_fit"))
  shape_col <- NULL
  if (length(lasso_i) > 0) {
    sel <- fits[[lasso_i[1]]]$fold_importance %>% dplyr::group_by(Variable) %>%
      dplyr::summarise(Lasso_pct_folds_selectionnee = 100 * mean(Used), .groups = "drop")
    table <- table %>% dplyr::left_join(sel, by = "Variable") %>%
      dplyr::mutate(Lasso_selectionnee = Lasso_pct_folds_selectionnee >= lasso_selection_threshold)
    shape_col <- "Lasso_selectionnee"
  }
  wide <- tidyr::pivot_wider(dplyr::select(table, "Variable", "Modele", "rank"), names_from = "Modele", values_from = "rank")
  p <- .bump_chart(table, "rank", labels, top_n,
                   "Classement des variables par importance SHAP (m\u00eame algorithme pour toutes les m\u00e9thodes)",
                   "Rang SHAP (1 = plus importante)", shape_col)
  if (!is.null(shape_col)) {
    p <- p + ggplot2::scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), na.value = 4,
                                         name = sprintf("S\u00e9lectionn\u00e9e par\nle Lasso (\u2265%d%%)", lasso_selection_threshold))
  }
  list(shap = shap_objs, table = table, rank_correlation = .rank_correlations(wide, labels), plot = p, time_sec = time_sec)
}

# ── Courbes superposées (binomial) ------------------------------------------------

.legend_inside <- function() {
  ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.95, 0.05),
                 legend.justification = c("right", "bottom"))
}

#' Superpose les courbes ROC de plusieurs modèles binomiaux
#' @param fits liste de fits binomiaux
#' @param labels libellés optionnels
#' @export
compare_roc <- function(fits, labels = NULL) {
  cmp <- .standardize_fits(fits, labels)
  stopifnot("compare_roc() requiert family = 'binomial'" = all(vapply(cmp, `[[`, "", "family") == "binomial"))
  roc <- purrr::map_dfr(cmp, function(cf) {
    fp <- cf$fold_predictions
    .roc_df(fp$Prediction, fp$Reel) %>%
      dplyr::mutate(Modele = sprintf("%s (AUC = %.3f)", cf$label, .compute_auc(fp$Prediction, fp$Reel)))
  })
  ggplot2::ggplot(roc, ggplot2::aes(x = FPR, y = TPR, colour = Modele)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_path(linewidth = 0.9) + ggplot2::coord_equal() +
    ggplot2::labs(title = "Comparaison des courbes ROC (pr\u00e9dictions pool\u00e9es)",
                  subtitle = "Pointill\u00e9 = classifieur al\u00e9atoire (AUC = 0.5)",
                  x = .lab_roc()$x, y = .lab_roc()$y, colour = NULL) +
    .theme_cv(13) + .legend_inside()
}

#' Superpose les courbes de calibration de plusieurs modèles binomiaux
#' @param fits liste de fits binomiaux
#' @param labels libellés optionnels
#' @param bins nombre de tranches (10)
#' @export
compare_calibration <- function(fits, labels = NULL, bins = 10) {
  cmp <- .standardize_fits(fits, labels)
  stopifnot("compare_calibration() requiert family = 'binomial'" = all(vapply(cmp, `[[`, "", "family") == "binomial"))
  cal <- purrr::map_dfr(cmp, function(cf) dplyr::mutate(.calib_df(cf$fold_predictions, bins), Modele = cf$label))
  ggplot2::ggplot(cal, ggplot2::aes(x = prob_moyenne, y = freq_observee, colour = Modele)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_line() + ggplot2::geom_point(ggplot2::aes(size = n), alpha = 0.8) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) + ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "Comparaison des courbes de calibration",
                  subtitle = "Pointill\u00e9 = calibration parfaite \u00b7 Taille = nb de pr\u00e9dictions dans le bin",
                  x = "Probabilit\u00e9 pr\u00e9dite (moyenne du bin)", y = "Fr\u00e9quence observ\u00e9e", colour = NULL, size = "N") +
    .theme_cv(13) + .legend_inside()
}

#' Temps de calcul par modèle (+ SHAP optionnel)
#'
#' @param fits liste d'au moins 1 fit
#' @param labels libellés optionnels
#' @param shap_result sortie de compare_shap() (ou check_shap() si un seul fit)
#' @return liste de classe "compare_computational_cost" (table, plot)
#' @export
compare_computational_cost <- function(fits, labels = NULL, shap_result = NULL) {
  cmp <- .standardize_fits(fits, labels, min_n = 1)
  tab <- purrr::map_dfr(cmp, function(cf) {
    ci <- cf$computational_info
    tibble::tibble(Modele = cf$label, temps_folds_sec = ci$time_folds_sec, temps_modele_final_sec = ci$time_final_model_sec,
                   temps_total_sec = ci$time_total_sec, temps_par_fold_sec = ci$time_folds_sec / max(1, ci$n_folds),
                   parallelisation = ci$parallelisation, n_coeurs = ci$n_cores_used, n_folds = ci$n_folds)
  })
  comp <- c("temps_folds_sec", "temps_modele_final_sec"); noms <- c("Folds (CV r\u00e9p\u00e9t\u00e9e)", "Mod\u00e8le final")
  if (!is.null(shap_result)) {
    ts <- if (!is.null(shap_result$shap)) {
      purrr::imap_dfr(shap_result$shap, ~ tibble::tibble(Modele = .y, temps_shap_sec = .x$time_sec))
    } else if (inherits(shap_result, "check_shap") && length(cmp) == 1) {
      tibble::tibble(Modele = names(cmp)[1], temps_shap_sec = shap_result$time_sec)
    } else NULL
    if (!is.null(ts)) {
      tab <- dplyr::left_join(tab, ts, by = "Modele") %>% dplyr::mutate(temps_shap_sec = tidyr::replace_na(temps_shap_sec, 0))
      comp <- c(comp, "temps_shap_sec"); noms <- c(noms, "SHAP")
    }
  }
  tab$temps_grand_total_sec <- rowSums(tab[, comp])
  tab <- dplyr::arrange(tab, temps_grand_total_sec)
  pd <- tab %>% dplyr::select("Modele", dplyr::all_of(comp)) %>%
    tidyr::pivot_longer(-"Modele", names_to = "Composante", values_to = "temps_sec") %>%
    dplyr::mutate(Composante = factor(Composante, levels = comp, labels = noms),
                  Modele = factor(Modele, levels = tab$Modele))
  p <- ggplot2::ggplot(pd, ggplot2::aes(x = Modele, y = temps_sec, fill = Composante)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(data = dplyr::mutate(tab, Modele = factor(Modele, levels = tab$Modele)),
                       ggplot2::aes(x = Modele, y = temps_grand_total_sec, label = sprintf("%.1fs", temps_grand_total_sec)),
                       inherit.aes = FALSE, vjust = -0.4, size = 3.3) +
    ggplot2::scale_fill_manual(values = stats::setNames(c("#2C7FB8", "#8AB4D6", "#D95F02")[seq_along(noms)], noms)) +
    ggplot2::labs(title = "Temps de calcul par m\u00e9thode", subtitle = "Empil\u00e9 : folds + mod\u00e8le final (+ SHAP)",
                  x = NULL, y = "Temps (secondes)", fill = NULL) +
    .theme_cv(12, 8.5) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  structure(list(table = tab, plot = p), class = c("compare_computational_cost", "list"))
}

#' @export
print.compare_computational_cost <- function(x, ...) { print(x$table); invisible(x) }
