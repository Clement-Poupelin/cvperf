# =============================================================================
# shap.R -- valeurs de Shapley, UN SEUL algorithme pour toutes les méthodes
# =============================================================================
# Approximation Monte Carlo de Štrumbelj & Kononenko (2014), via {fastshap},
# appliquée au MODÈLE FINAL. Plus aucune branche par méthode : tout passe par
# model_design(fit) (X, fonction de prédiction, correspondance colonne ->
# variable d'origine). Les valeurs de Shapley des colonnes indicatrices d'un
# même facteur sont SOMMÉES sous le nom de la variable d'origine, pour que
# toutes les méthodes désignent une variable par le même nom.
#
# Coût : ~ 2 x nsim x p appels à predict(). Leviers : parallelisation (sans
# effet sur la précision), sample_size, nsim.

#' Valeurs de Shapley (Monte Carlo, model-agnostic) du modèle final
#'
#' @param fit objet "cv_fit"
#' @param nsim répétitions Monte Carlo (défaut 50)
#' @param sample_size nombre d'observations expliquées (NULL = toutes)
#' @param seed graine
#' @param parallelisation TRUE = répartit les observations sur plusieurs coeurs
#'   (\code{parallel::mclapply()} : Unix/macOS ; séquentiel sous Windows)
#' @param n_cores NULL = 2/3 des coeurs détectés
#' @param method ignoré (conservé pour compatibilité ; la méthode est déduite de la classe)
#' @return liste de classe "check_shap"
#' @export
check_shap <- function(fit, nsim = 50, sample_size = NULL, seed = 42,
                       parallelisation = FALSE, n_cores = NULL, method = NULL) {
  .stop_if_not_cv_fit(fit)
  .require_pkg("fastshap")
  d <- model_design(fit)
  X <- as.data.frame(d$X)
  if (anyNA(X)) warning("Valeurs manquantes dans les predicteurs : calcul SHAP non teste dans ce cas.", call. = FALSE)
  pred_wrapper <- function(object, newdata) as.numeric(d$predict(newdata))

  n_total <- nrow(X)
  idx <- if (!is.null(sample_size) && sample_size < n_total) { set.seed(seed); sort(sample(seq_len(n_total), sample_size)) } else seq_len(n_total)
  X_explain <- X[idx, , drop = FALSE]

  n_cores_used <- .resolve_n_cores(parallelisation, n_cores, nrow(X_explain))
  if (n_cores_used > 1L && .Platform$OS.type == "windows") {
    message("check_shap() : parallelisation par fork indisponible sous Windows -- execution sequentielle.")
    n_cores_used <- 1L
  }
  t0 <- Sys.time()
  if (n_cores_used <= 1L) {
    set.seed(seed)
    shap_raw <- fastshap::explain(object = fit$final_model, X = X, newdata = X_explain,
                                  pred_wrapper = pred_wrapper, nsim = nsim, adjust = TRUE)
  } else {
    morceaux <- split(seq_len(nrow(X_explain)), cut(seq_len(nrow(X_explain)), n_cores_used, labels = FALSE))
    message(sprintf("check_shap() : %d observation(s) reparties sur %d coeurs...", nrow(X_explain), n_cores_used))
    res <- parallel::mclapply(seq_along(morceaux), function(i) {
      set.seed(seed + i)
      fastshap::explain(object = fit$final_model, X = X, newdata = X_explain[morceaux[[i]], , drop = FALSE],
                        pred_wrapper = pred_wrapper, nsim = nsim, adjust = TRUE)
    }, mc.cores = n_cores_used)
    shap_raw <- do.call(rbind, lapply(res, as.data.frame))[order(unlist(morceaux)), , drop = FALSE]
  }
  time_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  ids <- if (!is.null(fit$id_col) && fit$id_col %in% names(fit$data)) as.character(fit$data[[fit$id_col]][idx]) else as.character(idx)
  shap_long <- as.data.frame(shap_raw)
  names(shap_long) <- colnames(X)
  shap_long <- shap_long %>%
    dplyr::mutate(.id = ids) %>%
    tidyr::pivot_longer(-".id", names_to = "Variable", values_to = "shap_value") %>%
    dplyr::mutate(Variable = unname(d$var_map[Variable])) %>%
    dplyr::group_by(.id, Variable) %>%
    dplyr::summarise(shap_value = sum(shap_value), .groups = "drop")

  predictor_cols <- intersect(unique(unname(d$var_map)), names(fit$data))
  feature_long <- as.data.frame(fit$data)[idx, predictor_cols, drop = FALSE] %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
    dplyr::mutate(.id = ids) %>%
    tidyr::pivot_longer(-".id", names_to = "Variable", values_to = "feature_value")

  shap_long <- dplyr::left_join(shap_long, feature_long, by = c(".id", "Variable")) %>%
    dplyr::group_by(Variable) %>%
    dplyr::mutate(feature_value_scaled = {
      rng <- suppressWarnings(range(feature_value, na.rm = TRUE))
      if (all(is.finite(rng)) && diff(rng) > 0) (feature_value - rng[1]) / diff(rng) else NA_real_
    }) %>%
    dplyr::ungroup()

  summary_df <- shap_long %>%
    dplyr::group_by(Variable) %>%
    dplyr::summarise(mean_abs_shap = mean(abs(shap_value)), mean_shap = mean(shap_value),
                     sd_shap = stats::sd(shap_value), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_abs_shap)) %>%
    dplyr::mutate(rank = dplyr::row_number())

  structure(list(method = fit$method, method_label = fit$method_label, shap_long = shap_long,
                 summary = summary_df, n_explained = length(idx), n_total = n_total, nsim = nsim,
                 fit_family = fit$family, time_sec = time_sec,
                 parallelisation = n_cores_used > 1L, n_cores_used = n_cores_used),
            class = "check_shap")
}

#' @export
print.check_shap <- function(x, top_n = 15, ...) {
  cat(sprintf("Valeurs de Shapley (Monte Carlo, {fastshap}, nsim = %d) -- %s\n", x$nsim, x$method_label))
  cat(sprintf("%d observation(s) expliquee(s) sur %d | %.1fs | %s\n\n", x$n_explained, x$n_total, x$time_sec,
              if (isTRUE(x$parallelisation)) sprintf("parallelise sur %d coeurs", x$n_cores_used) else "sequentiel"))
  cat(sprintf("Top %d variable(s) par |SHAP| moyen :\n", min(top_n, nrow(x$summary))))
  print(knitr::kable(utils::head(dplyr::select(x$summary, Variable, mean_abs_shap, mean_shap, rank), top_n), digits = 4))
  invisible(x)
}

#' Représentation des valeurs de Shapley
#'
#' @param x objet "check_shap"
#' @param labels libellés optionnels
#' @param top_n variables affichées (défaut 20)
#' @param type "beeswarm" (1 point = 1 observation x 1 variable, couleur =
#'   valeur du prédicteur) ou "bar" (|SHAP| moyen)
#' @export
plot_shap <- function(x, labels = NULL, top_n = 20, type = c("beeswarm", "bar")) {
  type <- match.arg(type)
  stopifnot("x doit etre un objet 'check_shap'" = inherits(x, "check_shap"))
  top_vars <- utils::head(x$summary$Variable, top_n)
  order_df <- x$summary %>% dplyr::filter(Variable %in% top_vars) %>% dplyr::arrange(mean_abs_shap) %>%
    dplyr::mutate(Label = .label_lookup(Variable, labels))
  if (identical(type, "bar")) {
    return(
      ggplot2::ggplot(dplyr::mutate(order_df, Label = factor(Label, levels = unique(Label))),
                      ggplot2::aes(x = Label, y = mean_abs_shap)) +
        ggplot2::geom_col(fill = "#2C7FB8", alpha = 0.85) + ggplot2::coord_flip() +
        ggplot2::labs(title = sprintf("Importance SHAP moyenne -- %s", x$method_label),
                      subtitle = sprintf("Moyenne de |SHAP| sur %d observation(s) (nsim = %d)", x$n_explained, x$nsim),
                      x = NULL, y = "Moyenne(|SHAP|)") +
        .theme_cv(12, 8.5))
  }
  df <- x$shap_long %>% dplyr::filter(Variable %in% top_vars) %>%
    dplyr::mutate(Label = factor(.label_lookup(Variable, labels), levels = unique(order_df$Label)))
  ggplot2::ggplot(df, ggplot2::aes(x = shap_value, y = Label, colour = feature_value_scaled)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
    ggplot2::geom_jitter(height = 0.15, width = 0, alpha = 0.7, size = 1.5, na.rm = TRUE) +
    ggplot2::scale_colour_gradient(low = "#2C7FB8", high = "#D95F02", na.value = "grey70",
                                   name = "Valeur du\npr\u00e9dicteur\n(bas \u2192 haut)") +
    ggplot2::labs(title = sprintf("SHAP summary plot -- %s", x$method_label),
                  subtitle = sprintf("%d observation(s) (nsim = %d) \u00b7 top %d variable(s) par |SHAP| moyen",
                                     x$n_explained, x$nsim, length(top_vars)),
                  x = "Valeur de Shapley (impact sur la pr\u00e9diction)", y = NULL) +
    .theme_cv(12, 8.5)
}
