# =============================================================================
# plots_prediction.R -- plot_prediction_quality(), plot_residuals_vs_fitted(),
# plot_residuals_distribution() : une seule implémentation pour toutes les méthodes
# =============================================================================

#' Qualité de prédiction : 5 lectures complémentaires
#'
#' @param fit objet "cv_fit"
#' @param mode
#'   \itemize{
#'   \item "pooled" (défaut) : trait bleu = étendue min-max des prédictions de CV
#'     par valeur réelle ; point bleu = moyenne CV ; point orange = modèle FINAL ;
#'   \item "pooled_mean" : étendue + moyenne CV uniquement (lecture honnête) ;
#'   \item "final" : modèle final sur ses données d'entraînement (optimiste) ;
#'   \item "overlay" : toutes les prédictions brutes de CV + modèle final ;
#'   \item "panel" : les 3 lectures de base assemblées.
#'   }
#'   Pour family = "binomial", densités de probabilité par classe réelle
#'   ("pooled_mean" retombe sur "pooled").
#' @param base_size taille de police de base
#' @return ggplot ou patchwork
#' @export
plot_prediction_quality <- function(fit, mode = c("pooled", "pooled_mean", "final", "overlay", "panel"),
                                    base_size = 13) {
  .stop_if_not_cv_fit(fit)
  mode <- match.arg(mode)
  lbl  <- fit$method_label
  fp   <- fit$fold_predictions
  perf <- model_performance(fit)$pooled
  fin  <- .final_predictions(fit)
  df_final <- tibble::tibble(Reel = fin$Y, Prediction = fin$pred)
  note_final <- method_notes(fit, "final")
  n_folds <- dplyr::n_distinct(fp$Fold)

  # ── Réponse binaire : densités ---------------------------------------------
  if (identical(fit$family, "binomial")) {
    if (mode == "pooled_mean") {
      warning("mode = 'pooled_mean' n'a pas d'equivalent pour family = 'binomial' -- utilisation de 'pooled'.",
              call. = FALSE)
      mode <- "pooled"
    }
    lab_reel <- function(d) dplyr::mutate(d, Reel_lab = factor(Reel, levels = c(0, 1), labels = c("R\u00e9el : 0", "R\u00e9el : 1")))
    fp_lab <- lab_reel(fp)
    df_final_lab <- lab_reel(df_final)
    perf_final <- .compute_classification_metrics(df_final$Prediction, df_final$Reel, fit$classification_threshold)
    pal <- c("R\u00e9el : 0" = "#2C7FB8", "R\u00e9el : 1" = "#D95F02")

    build_density <- function(df, pv, title, subtitle, compact = FALSE) {
      ggplot2::ggplot(df, ggplot2::aes(x = Prediction, fill = Reel_lab, colour = Reel_lab)) +
        ggplot2::geom_density(alpha = 0.45, linewidth = 0.6) +
        ggplot2::geom_rug(sides = "b", alpha = 0.25, length = ggplot2::unit(0.04, "npc"), show.legend = FALSE) +
        ggplot2::geom_vline(xintercept = fit$classification_threshold, linetype = "dashed", colour = "grey40") +
        ggplot2::scale_fill_manual(values = pal, name = NULL) +
        ggplot2::scale_colour_manual(values = pal, guide = "none") +
        ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
                          label = sprintf(" Accuracy=%.3f\nAUC=%.3f\nBrier=%.3f", pv$accuracy, pv$auc, pv$brier),
                          size = if (compact) 3 else 4, fontface = "italic") +
        ggplot2::scale_x_continuous(limits = c(0, 1)) +
        ggplot2::labs(title = title, subtitle = subtitle, x = "Probabilit\u00e9 pr\u00e9dite", y = "Densit\u00e9") +
        .theme_cv(if (compact) base_size * 0.9 else base_size, if (compact) base_size * 0.6 else base_size * 0.8) +
        ggplot2::theme(legend.position = "top")
    }

    if (mode == "pooled") {
      return(build_density(fp_lab, perf,
        sprintf("Distribution des probabilit\u00e9s pr\u00e9dites (%s) \u2014 %d folds pool\u00e9s", lbl, n_folds),
        sprintf("Pointill\u00e9 = seuil (%.2f) \u00b7 %d pr\u00e9dictions (classe positive = \"%s\")",
                fit$classification_threshold, nrow(fp), fit$positive_class)))
    }
    if (mode == "final") {
      return(build_density(df_final_lab, perf_final,
        sprintf("Distribution des probabilit\u00e9s pr\u00e9dites (%s) \u2014 mod\u00e8le final uniquement", lbl),
        sprintf("\u26a0 Pr\u00e9dictions sur les donn\u00e9es d'ENTRA\u00ceNEMENT (%d obs.) \noptimiste par construction, ne pas reporter comme performance",
                nrow(df_final_lab))))
    }
    if (mode == "overlay") {
      combined <- dplyr::bind_rows(
        dplyr::mutate(fp_lab, Source = "CV (honn\u00eate)"),
        dplyr::mutate(df_final_lab, Source = "Mod\u00e8le final (entra\u00eenement)")
      ) %>% dplyr::mutate(Source = factor(Source, levels = c("CV (honn\u00eate)", "Mod\u00e8le final (entra\u00eenement)")))
      return(
        ggplot2::ggplot(combined, ggplot2::aes(x = Prediction, fill = Reel_lab, colour = Reel_lab, linetype = Source)) +
          ggplot2::geom_density(alpha = 0.2, linewidth = 0.7) +
          ggplot2::geom_vline(xintercept = fit$classification_threshold, linetype = "dotted", colour = "grey40") +
          ggplot2::scale_fill_manual(values = pal, name = NULL) +
          ggplot2::scale_colour_manual(values = pal, guide = "none") +
          ggplot2::scale_linetype_manual(values = c("CV (honn\u00eate)" = "solid", "Mod\u00e8le final (entra\u00eenement)" = "dashed"), name = NULL) +
          ggplot2::scale_x_continuous(limits = c(0, 1)) +
          ggplot2::labs(
            title = sprintf("CV (honn\u00eate) vs mod\u00e8le final (entra\u00eenement, optimiste) -- %s", lbl),
            subtitle = sprintf("Trait plein = CV pool\u00e9e (AUC=%.3f) \u00b7 \nTirets = mod\u00e8le final sur son propre entra\u00eenement (AUC=%.3f, optimiste)",
                               perf$auc, perf_final$auc),
            x = "Probabilit\u00e9 pr\u00e9dite", y = "Densit\u00e9") +
          .theme_cv(base_size, base_size * 0.8) + ggplot2::theme(legend.position = "top")
      )
    }
    p1 <- build_density(fp_lab, perf, "CV (honn\u00eate)", sprintf("AUC = %.3f", perf$auc), compact = TRUE)
    p2 <- build_density(df_final_lab, perf_final, "Mod\u00e8le final (entra\u00eenement)",
                        sprintf("\u26a0 optimiste \u00b7 AUC = %.3f", perf_final$auc), compact = TRUE)
    return((p1 | p2) + patchwork::plot_annotation(
      title = sprintf("Distribution des probabilit\u00e9s (%s) : CV honn\u00eate vs mod\u00e8le final", lbl),
      subtitle = "Gauche = performance CV honn\u00eate \u00b7 Droite = mod\u00e8le final sur ses propres donn\u00e9es d'entra\u00eenement (optimiste)",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = base_size * 1.6),
                             plot.subtitle = ggplot2::element_text(size = base_size * 0.9))))
  }

  # ── Réponse continue ---------------------------------------------------------
  pred_summary <- fp %>%
    dplyr::group_by(Reel) %>%
    dplyr::summarise(pred_min = min(Prediction), pred_max = max(Prediction),
                     pred_moyenne = mean(Prediction), n = dplyr::n(), .groups = "drop") %>%
    dplyr::left_join(df_final %>% dplyr::group_by(Reel) %>%
                       dplyr::summarise(pred_final_moyenne = mean(Prediction), .groups = "drop"),
                     by = "Reel")
  lims <- range(c(fp$Reel, fp$Prediction, df_final$Reel, df_final$Prediction), na.rm = TRUE)
  axes <- list(ggplot2::scale_x_continuous(limits = lims), ggplot2::scale_y_continuous(limits = lims),
               ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50"))
  th <- function(compact) .theme_cv(if (compact) base_size * 0.9 else base_size,
                                    if (compact) base_size * 0.6 else base_size * 0.8)
  lab_perf <- sprintf(" R\u00b2=%.3f\nRMSE=%.3g\nMAE=%.3g", perf$R2, perf$RMSE, perf$MAE)

  build_final <- function(compact = FALSE) {
    r2_final <- 1 - sum((df_final$Reel - df_final$Prediction)^2) / sum((df_final$Reel - mean(df_final$Reel))^2)
    ggplot2::ggplot(df_final, ggplot2::aes(x = Reel, y = Prediction)) + axes +
      ggplot2::geom_point(colour = "#D95F02", size = 1.6, alpha = 0.7) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#D95F02", linewidth = 0.8) +
      ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
                        label = sprintf(" R\u00b2 (entra\u00eenement) = %.3f", r2_final),
                        size = if (compact) 3 else 4, fontface = "italic") +
      ggplot2::labs(
        title = if (compact) "Entra\u00eenement (mod\u00e8le final)" else sprintf("Qualit\u00e9 de pr\u00e9diction (%s) \u2014 mod\u00e8le final uniquement", lbl),
        subtitle = if (compact) "\u26a0 optimiste" else paste0(
          sprintf("\u26a0 Pr\u00e9dictions sur les donn\u00e9es d'ENTRA\u00ceNEMENT (%d obs.) : optimiste par construction", nrow(df_final)),
          if (!is.null(note_final)) paste0("\n", note_final) else ""),
        x = "Valeur r\u00e9elle", y = "Valeur pr\u00e9dite") +
      th(compact) + ggplot2::theme(plot.subtitle = ggplot2::element_text(colour = "#D95F02"))
  }

  build_pooled_mean <- function(compact = FALSE) {
    ggplot2::ggplot(pred_summary, ggplot2::aes(x = Reel, y = pred_moyenne)) + axes +
      ggplot2::geom_linerange(ggplot2::aes(ymin = pred_min, ymax = pred_max), colour = "#2C7FB8", alpha = 0.45, linewidth = 0.4) +
      ggplot2::geom_point(colour = "#2C7FB8", size = if (compact) 1.2 else 1.4, alpha = 0.85) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#2C7FB8", linewidth = 0.8) +
      ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, label = lab_perf,
                        size = if (compact) 3 else 4, fontface = "italic") +
      ggplot2::labs(
        title = if (compact) "Performance CV (moyenne)" else sprintf("Qualit\u00e9 de pr\u00e9diction (%s) : %d folds pool\u00e9s (moyenne)", lbl, n_folds),
        subtitle = if (compact) "Point = moyenne CV (honn\u00eate)" else
          sprintf("Point BLEU = pr\u00e9diction moyenne de la CV par valeur r\u00e9elle \nTrait = \u00e9tendue min-max (%d valeurs distinctes)", nrow(pred_summary)),
        x = "Valeur r\u00e9elle", y = "Valeur pr\u00e9dite") +
      th(compact)
  }

  build_pooled_final <- function(compact = FALSE) {
    ggplot2::ggplot(pred_summary, ggplot2::aes(x = Reel, y = pred_final_moyenne)) + axes +
      ggplot2::geom_linerange(ggplot2::aes(ymin = pred_min, ymax = pred_max), colour = "#2C7FB8", alpha = 0.45, linewidth = 0.4) +
      ggplot2::geom_point(ggplot2::aes(y = pred_moyenne), colour = "#2C7FB8", size = if (compact) 1.2 else 1.4, alpha = 0.85) +
      ggplot2::geom_smooth(ggplot2::aes(y = pred_moyenne), method = "lm", formula = y ~ x, se = FALSE, colour = "#2C7FB8", linewidth = 0.8) +
      ggplot2::geom_point(colour = "#D95F02", size = if (compact) 1.4 else 1.6, alpha = 0.85) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#D95F02", linewidth = 0.8) +
      ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, label = lab_perf,
                        size = if (compact) 3 else 4, fontface = "italic") +
      ggplot2::labs(
        title = if (compact) "Mod\u00e8le final + \u00e9tendue CV" else sprintf("Qualit\u00e9 de pr\u00e9diction (%s) : %d folds pool\u00e9s + mod\u00e8le final", lbl, n_folds),
        subtitle = if (compact) "Orange = mod\u00e8le final \u00b7 bleu = moyenne CV" else paste0(
          "Point orange = pr\u00e9diction du mod\u00e8le FINAL \u00b7 Point bleu = pr\u00e9diction MOYENNE de la CV",
          sprintf("\nTrait bleu = \u00e9tendue min-max des pr\u00e9dictions de CV (%d valeurs distinctes)", nrow(pred_summary)),
          "\n\u26a0 R\u00b2/RMSE/MAE affich\u00e9s = estimation honn\u00eate par CV (pool\u00e9e) ; les points orange sont optimistes"),
        x = "Valeur r\u00e9elle", y = "Valeur pr\u00e9dite") +
      th(compact)
  }

  if (mode == "final")       return(build_final())
  if (mode == "pooled_mean") return(build_pooled_mean())
  if (mode == "pooled")      return(build_pooled_final())
  if (mode == "overlay") {
    return(
      ggplot2::ggplot() + axes +
        ggplot2::geom_jitter(data = fp, ggplot2::aes(x = Reel, y = Prediction), colour = "#2C7FB8",
                             alpha = 0.12, size = 1, width = diff(lims) * 0.004, height = 0) +
        ggplot2::geom_point(data = df_final, ggplot2::aes(x = Reel, y = Prediction), colour = "#D95F02", size = 1.8, alpha = 0.9) +
        ggplot2::labs(
          title = sprintf("Qualit\u00e9 de pr\u00e9diction (%s) \u2014 tous les folds + mod\u00e8le final", lbl),
          subtitle = sprintf("Bleu = %d pr\u00e9dictions des %d folds (CV) \nOrange = mod\u00e8le final (entra\u00eenement)", nrow(fp), n_folds),
          x = "Valeur r\u00e9elle", y = "Valeur pr\u00e9dite") +
        .theme_cv(base_size, base_size * 0.7)
    )
  }
  (build_pooled_mean(TRUE) / (build_final(TRUE) | build_pooled_final(TRUE))) +
    patchwork::plot_annotation(
      title = sprintf("Qualit\u00e9 de pr\u00e9diction (%s) : 3 lectures compl\u00e9mentaires", lbl),
      subtitle = "Haut : performance CV honn\u00eate \u00b7 Bas gauche : ajustement sur toutes les donn\u00e9es (optimiste) \u00b7 Bas droite : mod\u00e8le final avec l'\u00e9tendue de la CV",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = base_size * 1.6),
                             plot.subtitle = ggplot2::element_text(size = base_size * 0.9)))
}

#' Résidus vs valeurs prédites (résidus poolés de CV)
#' @param fit objet "cv_fit" (family = "gaussian")
#' @export
plot_residuals_vs_fitted <- function(fit) {
  .stop_if_not_cv_fit(fit)
  ggplot2::ggplot(fit$fold_predictions, ggplot2::aes(x = Prediction, y = Residu)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(size = 1.6, alpha = 0.25, colour = "#2C7FB8") +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE, colour = "#D95F02",
                         linewidth = 0.8, fill = "#D95F02", alpha = 0.15) +
    ggplot2::labs(
      title = sprintf("R\u00e9sidus vs valeurs pr\u00e9dites (%s)", fit$method_label),
      subtitle = method_notes(fit, "residuals") %||%
        "Une courbe qui s'\u00e9carte franchement de 0 signale une zone de la plage de Y mal mod\u00e9lis\u00e9e",
      x = "Pr\u00e9diction", y = "R\u00e9sidu (r\u00e9el - pr\u00e9dit)") +
    .theme_cv(13, 8.5)
}

#' Distribution des résidus poolés
#' @param fit objet "cv_fit" (family = "gaussian")
#' @export
plot_residuals_distribution <- function(fit) {
  .stop_if_not_cv_fit(fit)
  ggplot2::ggplot(fit$fold_predictions, ggplot2::aes(x = Residu)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 25,
                            fill = "#2C7FB8", alpha = 0.7, colour = "white") +
    ggplot2::geom_density(colour = "#D95F02", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::labs(
      title = sprintf("Distribution des r\u00e9sidus pool\u00e9s (%s)", fit$method_label),
      subtitle = sprintf("n = %d r\u00e9sidus (%d folds)", nrow(fit$fold_predictions), fit$n_folds_valides),
      x = "R\u00e9sidu", y = "Densit\u00e9") +
    .theme_cv(13)
}
