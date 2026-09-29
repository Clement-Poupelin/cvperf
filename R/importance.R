# =============================================================================
# importance.R -- stabilité et lecture de l'importance, communes à toutes les
# méthodes (fit$fold_importance : Variable / Importance / Fold, à l'échelle de
# la variable d'ORIGINE)
# =============================================================================

#' Stabilité de l'importance des variables à travers les folds
#'
#' Pendant du VIF pour les méthodes qui n'estiment pas de coefficients : la
#' question n'est plus "les prédicteurs sont-ils trop corrélés pour être
#' estimés" mais "l'importance attribuée est-elle stable d'un tirage à l'autre".
#' \code{rank_sd} (écart-type du rang) complète le coefficient de variation,
#' peu fiable quand l'importance moyenne est proche de 0.
#'
#' @param fit objet "cv_fit"
#' @param top_n nombre de variables conservées (défaut 30)
#' @return tibble de classe "check_importance_cv"
#' @export
check_importance <- function(fit, top_n = 30) {
  .stop_if_not_cv_fit(fit)
  summary_df <- fit$fold_importance %>%
    dplyr::group_by(Fold) %>%
    dplyr::mutate(Rank = rank(-Importance, ties.method = "average")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Variable) %>%
    dplyr::summarise(
      mean_importance = mean(Importance),
      sd_importance   = stats::sd(Importance),
      cv_importance   = dplyr::if_else(mean(Importance) > 0, stats::sd(Importance) / mean(Importance), NA_real_),
      mean_rank       = mean(Rank),
      rank_sd         = stats::sd(Rank),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_importance)) %>%
    dplyr::slice_head(n = top_n)
  structure(summary_df, class = c("check_importance_cv", class(summary_df)),
            n_folds = dplyr::n_distinct(fit$fold_importance$Fold), top_n = top_n,
            importance_type = fit$importance_type, method_label = fit$method_label)
}

#' @export
print.check_importance_cv <- function(x, ...) {
  cat(sprintf("Stabilite de l'importance (%s) sur %d folds, top %d variable(s) :\n",
              attr(x, "importance_type"), attr(x, "n_folds"), attr(x, "top_n")))
  print(knitr::kable(dplyr::select(tibble::as_tibble(x), Variable, mean_importance, sd_importance, mean_rank, rank_sd), digits = 3))
  cat("\nrank_sd faible = variable classee de facon coherente d'un fold a l'autre (signal stable).\n")
  cat("rank_sd eleve  = le classement de cette variable varie beaucoup selon l'echantillon.\n")
  invisible(x)
}

#' Lollipop de l'importance moyenne, coloré par instabilité du rang
#'
#' @param x objet "check_importance_cv"
#' @param labels vecteur nommé optionnel de libellés
#' @param fit objet "cv_fit" optionnel : ajoute le signe de la corrélation de
#'   Spearman MARGINALE entre chaque variable numérique et la réponse (repère de
#'   lecture, PAS un effet partiel)
#' @export
plot_importance <- function(x, labels = NULL, fit = NULL) {
  stopifnot(inherits(x, "check_importance_cv"))
  df <- tibble::as_tibble(x) %>%
    dplyr::mutate(Label = .label_lookup(Variable, labels)) %>%
    dplyr::arrange(mean_importance) %>%
    dplyr::mutate(Label = factor(Label, levels = unique(Label)))
  if (!is.null(fit)) {
    direction <- purrr::map_dbl(as.character(df$Variable), function(v) {
      if (!v %in% names(fit$data) || !is.numeric(fit$data[[v]])) return(NA_real_)
      suppressWarnings(stats::cor(fit$data[[v]], as.numeric(fit$data[[fit$response_var]]),
                                  method = "spearman", use = "pairwise.complete.obs"))
    })
    df$Direction <- dplyr::case_when(is.na(direction) ~ "n/d (facteur ou N/A)",
                                     direction > 0 ~ "corr\u00e9lation marginale +",
                                     direction < 0 ~ "corr\u00e9lation marginale -",
                                     TRUE ~ "\u2248 0")
  } else {
    df$Direction <- "non calcul\u00e9e"
  }
  ggplot2::ggplot(df, ggplot2::aes(x = Label, y = mean_importance, colour = rank_sd)) +
    ggplot2::geom_segment(ggplot2::aes(xend = Label, y = 0, yend = mean_importance), linewidth = 0.9) +
    ggplot2::geom_point(ggplot2::aes(shape = Direction), size = 3.5) +
    ggplot2::scale_colour_gradient(low = "#3aaf85", high = "#cd201f", name = "Instabilit\u00e9\n(\u00e9cart-type du rang)") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = sprintf("Importance des variables (%s, %s) et stabilit\u00e9 de leur rang",
                      attr(x, "method_label"), attr(x, "importance_type")),
      subtitle = paste0("Vert = classement stable d'un fold \u00e0 l'autre \u00b7 Rouge = classement instable",
                        if (!is.null(fit)) "\nForme du point = signe de la corr\u00e9lation marginale (brute) avec la r\u00e9ponse -- PAS un effet partiel" else ""),
      x = NULL, y = "Importance moyenne", shape = NULL) +
    .theme_cv(12, 8.5)
}

#' Fréquence de "sélection" et amplitude de l'importance à travers les folds
#'
#' Méthode par défaut (forêt, SVM, XGBoost) : une variable est comptée
#' "retenue" sur un fold si elle est dans le top-\code{top_k} de ce fold.
#' Le Lasso (fréquence de coefficient non nul + amplitude signée) et MARS
#' (fréquence de présence dans le modèle élagué) ont leur propre méthode.
#'
#' @param fit objet "cv_fit"
#' @param labels libellés optionnels
#' @param top_k seuil de rang pour être "retenue" sur un fold (défaut 15)
#' @param top_n nombre maximal de variables affichées (défaut 30)
#' @param ... arguments propres à la méthode
#' @return patchwork (2 panneaux)
#' @export
plot_variable_importance <- function(fit, labels = NULL, top_k = 15, top_n = 30, ...) {
  UseMethod("plot_variable_importance")
}

#' @export
plot_variable_importance.cv_fit <- function(fit, labels = NULL, top_k = 15, top_n = 30, ...) {
  ranks <- fit$fold_importance %>%
    dplyr::group_by(Fold) %>%
    dplyr::mutate(Rank = rank(-Importance, ties.method = "average"), in_topk = Rank <= top_k) %>%
    dplyr::ungroup()
  full <- ranks %>%
    dplyr::group_by(Variable) %>%
    dplyr::summarise(n_in_topk = sum(in_topk), pct_in_topk = 100 * sum(in_topk) / fit$n_folds_valides,
                     .groups = "drop") %>%
    dplyr::mutate(Label = .label_lookup(Variable, labels)) %>%
    dplyr::arrange(dplyr::desc(pct_in_topk))
  shown <- dplyr::slice_head(full, n = top_n)
  n_tronque <- nrow(full) - nrow(shown)
  label_order <- shown %>% dplyr::arrange(pct_in_topk) %>% dplyr::pull(Label)

  p_freq <- ggplot2::ggplot(shown, ggplot2::aes(x = factor(Label, levels = label_order), y = pct_in_topk)) +
    ggplot2::geom_col(width = 0.7, fill = "#2C7FB8", na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 50, colour = "grey40", linetype = "dashed") +
    ggplot2::coord_flip() + ggplot2::scale_y_continuous(limits = c(0, 100)) +
    ggplot2::labs(title = sprintf("Fr\u00e9quence dans le top-%d par importance (%s)", top_k, fit$method_label),
                  subtitle = sprintf("%d variable(s), %d folds%s", nrow(full), fit$n_folds_valides,
                                     if (n_tronque > 0) sprintf(" \ntop %d affich\u00e9es, %d masqu\u00e9es", top_n, n_tronque) else ""),
                  x = NULL, y = sprintf("%% des folds o\u00f9 la variable est dans le top-%d", top_k)) +
    .theme_cv(12, 8.5)

  p_imp <- fit$fold_importance %>%
    dplyr::filter(Variable %in% shown$Variable) %>%
    dplyr::mutate(Label = factor(.label_lookup(Variable, labels), levels = label_order)) %>%
    ggplot2::ggplot(ggplot2::aes(x = Label, y = Importance)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
    ggplot2::geom_boxplot(fill = "#2C7FB8", alpha = 0.4, width = 0.5, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, height = 0, alpha = 0.25, size = 1) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Amplitude de l'importance (tous folds)",
                  subtitle = "Dispersion large = importance instable selon l'\u00e9chantillon",
                  x = NULL, y = sprintf("Importance (%s)", fit$importance_type)) +
    .theme_cv(12) + ggplot2::theme(axis.text.y = ggplot2::element_blank())
  p_freq | p_imp
}

#' Importance de GROUPES de variables (permutation d'un bloc)
#'
#' Disponible pour les méthodes à importance par permutation (SVM, XGBoost),
#' si \code{perm_groups} a été passé à l'ajustement. Un groupe est "utile" si
#' son importance est positive dans la quasi-totalité des folds et si
#' [q025 ; q975] exclut 0 (variabilité entre folds, pas une incertitude
#' d'échantillonnage sur de nouveaux individus).
#' @param fit objet "cv_fit"
#' @return tibble Groupe / mean_importance / sd_importance / q025 / q975 / pct_folds_positive
#' @export
group_importance <- function(fit) {
  .stop_if_not_cv_fit(fit)
  gi <- fit$fold_group_importance
  if (is.null(gi) || nrow(gi) == 0) {
    stop("Aucun groupe calcule : passer `perm_groups = list(Nom = c(\"var1\", \"var2\"))` a l'ajustement (SVM, XGBoost).",
         call. = FALSE)
  }
  gi %>%
    dplyr::group_by(Groupe = Variable) %>%
    dplyr::summarise(
      mean_importance = mean(Importance), sd_importance = stats::sd(Importance),
      q025 = unname(stats::quantile(Importance, 0.025)), q975 = unname(stats::quantile(Importance, 0.975)),
      pct_folds_positive = 100 * mean(Importance > 0), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_importance))
}
