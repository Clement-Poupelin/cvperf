# =============================================================================
# diagnostics_residuals.R -- check_normality / check_heteroscedasticity /
# check_outliers / check_independence, communs à toutes les méthodes
# =============================================================================
# Ils n'opèrent que sur fit$fold_predictions (résidus POOLÉS sur folds répétés :
# non indépendants -> p-values indicatives). Seule l'INTERPRÉTATION change
# selon la méthode : elle vient de method_notes().

.stop_if_binomial_resid <- function(fit, fun, alternative) {
  if (identical(fit$family, "binomial")) {
    stop(sprintf("%s() n'est pas pertinent pour une reponse binaire. Utiliser %s a la place.", fun, alternative),
         call. = FALSE)
  }
}

# ── Normalité -----------------------------------------------------------------

#' Normalité des résidus poolés (Shapiro-Wilk + Anderson-Darling)
#'
#' Pour un modèle linéaire (Lasso), une vraie question de spécification ; pour
#' les méthodes non paramétriques, une lecture DESCRIPTIVE permettant de
#' comparer la forme des résidus entre méthodes.
#' @param fit objet "cv_fit" (family = "gaussian")
#' @return liste de classe "check_normality_cv" (méthodes print() et plot())
#' @export
check_normality <- function(fit) {
  .stop_if_not_cv_fit(fit)
  .stop_if_binomial_resid(fit, "check_normality", "plot_calibration() ou plot_roc()")
  residus <- fit$fold_predictions$Residu
  # shapiro.test() est plafonné à n = 5000 : échantillon aléatoire au-delà
  residus_sw <- if (length(residus) > 5000) sample(residus, 5000) else residus
  sw <- stats::shapiro.test(residus_sw)
  ad <- nortest::ad.test(residus)
  structure(
    list(
      shapiro = list(statistic = unname(sw$statistic), p = sw$p.value, stars = p_to_stars(sw$p.value)),
      anderson_darling = list(statistic = unname(ad$statistic), p = ad$p.value, stars = p_to_stars(ad$p.value)),
      n = length(residus), residus = residus,
      method_label = fit$method_label, note = method_notes(fit, "normality")
    ),
    class = "check_normality_cv"
  )
}

#' @export
print.check_normality_cv <- function(x, ...) {
  cat(sprintf("Shapiro-Wilk       : W = %.3f, p = %.3f (%s)\nAnderson-Darling   : A = %.3f, p = %.3f (%s)\n",
              x$shapiro$statistic, x$shapiro$p, x$shapiro$stars,
              x$anderson_darling$statistic, x$anderson_darling$p, x$anderson_darling$stars))
  if (!is.null(x$note)) cat(strwrap(paste("\u2139", x$note), width = 90, exdent = 2), sep = "\n")
  cat("\u26a0 Residus POOLES sur folds repetes : non-independance -> p-values indicatives.\n")
  invisible(x)
}

#' @export
plot.check_normality_cv <- function(x, ...) {
  ggplot2::ggplot(data.frame(Residu = x$residus), ggplot2::aes(sample = Residu)) +
    ggplot2::stat_qq(colour = "#2C7FB8", alpha = 0.35, size = 1.4) +
    ggplot2::stat_qq_line(colour = "#D95F02", linewidth = 0.8) +
    ggplot2::labs(
      title    = sprintf("Q-Q plot des r\u00e9sidus (%s)", x$method_label),
      subtitle = if (is.null(x$note)) sprintf("n = %d r\u00e9sidus pool\u00e9s", x$n)
                 else "Lecture descriptive -- pas une hypoth\u00e8se du mod\u00e8le",
      x = "Quantiles th\u00e9oriques", y = "Quantiles observ\u00e9s"
    ) +
    .theme_cv(13)
}

# ── Hétéroscédasticité --------------------------------------------------------

#' Homogénéité de la variance des résidus (Spearman |résidu| vs prédiction)
#' @param fit objet "cv_fit" (family = "gaussian")
#' @return liste de classe "check_heteroscedasticity_cv"
#' @export
check_heteroscedasticity <- function(fit) {
  .stop_if_not_cv_fit(fit)
  .stop_if_binomial_resid(fit, "check_heteroscedasticity", "plot_calibration()")
  fp <- fit$fold_predictions
  test <- suppressWarnings(stats::cor.test(abs(fp$Residu), fp$Prediction, method = "spearman"))
  structure(
    list(rho = unname(test$estimate), p = test$p.value, stars = p_to_stars(test$p.value), data = fp,
         method_label = fit$method_label, note = method_notes(fit, "heteroscedasticity")),
    class = "check_heteroscedasticity_cv"
  )
}

#' @export
print.check_heteroscedasticity_cv <- function(x, ...) {
  verdict <- if (x$p < 0.05) "Tendance detectee : la variance de l'erreur semble varier avec le niveau de prediction."
             else "Pas de tendance forte detectee entre variance de l'erreur et niveau de prediction."
  cat(sprintf("Spearman(|residu|, prediction) : rho = %.3f, p = %.3f (%s)\n%s\n", x$rho, x$p, x$stars, verdict))
  if (!is.null(x$note)) cat(strwrap(paste("\u2139", x$note), width = 90, exdent = 2), sep = "\n")
  invisible(x)
}

#' @export
plot.check_heteroscedasticity_cv <- function(x, ...) {
  df <- x$data %>% dplyr::mutate(sqrt_abs_resid = sqrt(abs(scale(Residu)[, 1])))
  ggplot2::ggplot(df, ggplot2::aes(x = Prediction, y = sqrt_abs_resid)) +
    ggplot2::geom_point(size = 1.6, alpha = 0.25, colour = "#2C7FB8") +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE, colour = "#D95F02",
                         linewidth = 0.8, fill = "#D95F02", alpha = 0.15) +
    ggplot2::labs(
      title    = sprintf("Scale-Location (%s)", x$method_label),
      subtitle = paste0(sprintf("Spearman(|r\u00e9sidu|, pr\u00e9diction) rho = %.3f (%s)", x$rho, x$stars),
                        if (!is.null(x$note)) paste0("\n", x$note) else ""),
      x = "Pr\u00e9diction", y = expression(sqrt("|r\u00e9sidu standardis\u00e9|"))
    ) +
    .theme_cv(13, 8.5)
}

# ── Outliers ------------------------------------------------------------------

#' Individus à erreur systématique (résidu moyen à travers les répétitions)
#'
#' Propriété de la PROCÉDURE de CV répétée : identique pour toutes les méthodes.
#' @param fit objet "cv_fit"
#' @param threshold nombre d'écarts-types au-delà duquel un résidu moyen est signalé
#' @return tibble de classe "check_outliers_cv"
#' @export
check_outliers <- function(fit, threshold = 2) {
  .stop_if_not_cv_fit(fit)
  pe <- fit$fold_predictions %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(n_apparitions = dplyr::n(), residu_moyen = mean(Residu),
                     residu_sd = stats::sd(Residu), .groups = "drop")
  if (all(pe$n_apparitions <= 1)) {
    warning("check_outliers() : chaque individu n'apparait qu'une fois (repeats = 1) -- residu_sd est NA. Augmenter `repeats`.",
            call. = FALSE)
  }
  seuil <- threshold * stats::sd(pe$residu_moyen)
  out <- pe %>%
    dplyr::mutate(is_outlier = abs(residu_moyen) > seuil) %>%
    dplyr::arrange(dplyr::desc(abs(residu_moyen)))
  structure(out, class = c("check_outliers_cv", class(out)),
            threshold = threshold, seuil = seuil, method_label = fit$method_label)
}

#' @export
print.check_outliers_cv <- function(x, ...) {
  n_out <- sum(x$is_outlier)
  cat(sprintf("%d individu(s) signale(s) comme outlier (|residu moyen| > %.1f ecarts-types, seuil = %.2f)\n",
              n_out, attr(x, "threshold"), attr(x, "seuil")))
  if (n_out > 0) print(knitr::kable(dplyr::filter(tibble::as_tibble(x), is_outlier), digits = 2))
  invisible(x)
}

#' @export
plot.check_outliers_cv <- function(x, ...) {
  lbl <- attr(x, "method_label")
  df <- tibble::as_tibble(x)
  ggplot2::ggplot(df, ggplot2::aes(x = residu_moyen, y = residu_sd)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(ggplot2::aes(colour = is_outlier), size = 3, alpha = 0.7) +
    ggrepel::geom_text_repel(data = dplyr::filter(df, is_outlier), ggplot2::aes(label = ID),
                             size = 3, colour = "#D95F02") +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#D95F02", `FALSE` = "#2C7FB8"), guide = "none") +
    ggplot2::labs(
      title    = sprintf("Individus \u00e0 erreur syst\u00e9matique (%s)", lbl),
      subtitle = "x = biais moyen de l'individu \ny = variabilit\u00e9 de son erreur selon le fold \nOrange = signal\u00e9",
      x = "R\u00e9sidu moyen (r\u00e9el - pr\u00e9dit)", y = "\u00c9cart-type du r\u00e9sidu"
    ) +
    .theme_cv(13)
}

# ── Indépendance --------------------------------------------------------------

#' Cohérence du résidu d'un même individu entre répétitions
#' @param fit objet "cv_fit"
#' @return liste de classe "check_independence_cv"
#' @export
check_independence <- function(fit) {
  .stop_if_not_cv_fit(fit)
  rl <- fit$fold_predictions %>%
    dplyr::group_by(ID) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::arrange(Repeat, .by_group = TRUE) %>%
    dplyr::mutate(Residu_lag = dplyr::lag(Residu)) %>%
    dplyr::ungroup() %>%
    tidyr::drop_na(Residu_lag)
  if (nrow(rl) == 0) {
    stop("check_independence() : aucun individu n'apparait plus d'une fois (repeats = 1). Augmenter `repeats` (>= 2, 20 recommande).",
         call. = FALSE)
  }
  structure(list(correlation = stats::cor(rl$Residu, rl$Residu_lag), data = rl,
                 method_label = fit$method_label),
            class = "check_independence_cv")
}

#' @export
print.check_independence_cv <- function(x, ...) {
  cat(sprintf("Correlation du residu entre repetitions (meme individu) : %.3f\n", x$correlation))
  cat("Elevee = le residu depend surtout de l'INDIVIDU, peu du tirage train/test.\n")
  invisible(x)
}

#' @export
plot.check_independence_cv <- function(x, ...) {
  ggplot2::ggplot(x$data, ggplot2::aes(x = Residu_lag, y = Residu)) +
    ggplot2::geom_point(alpha = 0.3, size = 1.4, colour = "#2C7FB8") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#D95F02", linewidth = 0.8) +
    ggplot2::labs(
      title    = sprintf("Coh\u00e9rence du r\u00e9sidu d'un individu entre 2 r\u00e9p\u00e9titions (%s)", x$method_label),
      subtitle = sprintf("Corr\u00e9lation = %.3f (r\u00e9p\u00e9titions interchangeables)", x$correlation),
      x = "R\u00e9sidu (r\u00e9p\u00e9tition A)", y = "R\u00e9sidu (r\u00e9p\u00e9tition B, suivante)"
    ) +
    .theme_cv(13)
}
