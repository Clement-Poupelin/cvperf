# =============================================================================
# effects.R -- dépendance partielle générique (toutes méthodes)
# =============================================================================

#' Forme de l'effet des variables (dépendance partielle du modèle final)
#'
#' Pour chaque variable v et chaque valeur g d'une grille : prédiction moyenne
#' du modèle final quand v est fixée à g pour TOUS les individus. Tirets gris =
#' meilleure droite ajustée sur la courbe (écart = non-linéarité). Seules les
#' variables représentées par UNE colonne numérique sont tracées (variables
#' numériques, facteurs à 2 niveaux encodés en indicatrice) ; les autres
#' facteurs sont ignorés avec un message.
#'
#' Limite : si deux variables sont très corrélées, la grille évalue le modèle
#' sur des combinaisons irréalistes.
#'
#' @param fit objet "cv_fit"
#' @param variables variables d'ORIGINE à tracer ; NULL = les \code{top_n} plus
#'   importantes (MARS : parmi celles retenues par le modèle final ; Lasso :
#'   parmi les coefficients non nuls)
#' @param top_n nombre de variables par défaut (6)
#' @param n_grid nombre de points de grille (40)
#' @param labels libellés optionnels
#' @param ncol colonnes de facettes
#' @return ggplot
#' @export
plot_effects <- function(fit, variables = NULL, top_n = 6, n_grid = 40, labels = NULL, ncol = 3) {
  .stop_if_not_cv_fit(fit)
  d <- model_design(fit)
  X <- d$X
  if (is.null(variables)) variables <- .effect_variables(fit, top_n)

  cols <- vapply(variables, function(v) {
    cc <- names(d$var_map)[unname(d$var_map) == v]
    if (length(cc) == 1 && is.numeric(X[, cc])) cc else NA_character_
  }, character(1))
  ignorees <- variables[is.na(cols)]
  if (length(ignorees) > 0) {
    message("plot_effects() : variable(s) ignoree(s) (facteur a plus de 2 niveaux, ou introuvable) : ",
            paste(ignorees, collapse = ", "))
  }
  keep <- !is.na(cols)
  variables <- variables[keep]; cols <- cols[keep]
  if (length(variables) == 0) stop("Aucune variable numerique a tracer.", call. = FALSE)

  pdp <- purrr::map2_dfr(variables, cols, function(v, cl) {
    x_ok <- X[, cl]; x_ok <- x_ok[!is.na(x_ok)]
    grille <- if (length(unique(x_ok)) <= 2) sort(unique(x_ok)) else
      seq(stats::quantile(x_ok, 0.02), stats::quantile(x_ok, 0.98), length.out = n_grid)
    tibble::tibble(Variable = v, x = grille, y = vapply(grille, function(g) {
      Xg <- X; Xg[, cl] <- g
      mean(d$predict(Xg), na.rm = TRUE)
    }, numeric(1)))
  })
  pdp$Label <- factor(.label_lookup(pdp$Variable, labels), levels = unique(.label_lookup(variables, labels)))
  rug_df <- purrr::map2_dfr(variables, cols, function(v, cl) tibble::tibble(Variable = v, x = X[, cl]))
  rug_df <- rug_df[!is.na(rug_df$x), ]
  rug_df$Label <- factor(.label_lookup(rug_df$Variable, labels), levels = levels(pdp$Label))

  ggplot2::ggplot(pdp, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey55", linetype = "dashed", linewidth = 0.6) +
    ggplot2::geom_line(colour = "#D95F02", linewidth = 1) +
    ggplot2::geom_point(colour = "#D95F02", size = 1.2) +
    ggplot2::geom_rug(data = rug_df, ggplot2::aes(x = x), inherit.aes = FALSE, sides = "b",
                      alpha = 0.3, length = ggplot2::unit(0.04, "npc")) +
    ggplot2::facet_wrap(~ Label, scales = "free", ncol = ncol) +
    ggplot2::labs(
      title = sprintf("Forme de l'effet des variables les plus importantes (%s, mod\u00e8le final)", fit$method_label),
      subtitle = "Trait orange = d\u00e9pendance partielle \u00b7 Tirets gris = meilleure droite (\u00e9cart = non-lin\u00e9arit\u00e9) \u00b7 Rug = valeurs observ\u00e9es",
      x = "Valeur de la variable",
      y = if (identical(fit$family, "binomial")) "Probabilit\u00e9 pr\u00e9dite moyenne" else "Pr\u00e9diction moyenne") +
    .theme_cv(12, 8.5) + ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
}
