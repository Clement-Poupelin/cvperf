# =============================================================================
# cv_engine.R -- LE moteur de validation croisée répétée, commun à toutes les
# méthodes (remplace les 5 copies de fit_one_fold()/.run_folds_*()).
# =============================================================================
# Chaque méthode fournit seulement une fonction
#     fold_fun(train, test, fold_num) -> NULL (fold ignoré) ou
#                                        list(Y, pred, extras = list(<tibbles>))
# Le moteur se charge du découpage, des graines, de la parallélisation, des
# identifiants, des résidus, de l'assemblage et du chronométrage.
#
# REPRODUCTIBILITÉ : set.seed(seed + fold_num) est fixé AVANT chaque fold, pour
# toutes les méthodes. Le résultat ne dépend donc ni de l'ordre d'exécution ni
# du nombre de coeurs (y compris pour le Lasso, dont cv.glmnet() tire ses folds
# internes au hasard).

#' Moteur de CV répétée (interne, utilisé par tous les fit_*_repeated_cv())
#'
#' @param data,formula,id_col,family,v,repeats,seed,cv_repeats cf. fit_*_repeated_cv()
#' @param fold_fun fonction(train, test, fold_num)
#' @param parallelisation,n_cores cf. fit_*_repeated_cv()
#' @param pkgs packages à charger sur les workers parallèles
#' @return liste : cv_repeats, fold_predictions, extras, n_folds_valides,
#'   n_folds_total, n_cores_used, parallelisation, time_folds_sec
#' @keywords internal
#' @export
run_repeated_cv <- function(data, formula, id_col, family, v, repeats, seed, cv_repeats,
                             fold_fun, parallelisation = FALSE, n_cores = NULL,
                             pkgs = character()) {

  if (is.null(cv_repeats)) {
    set.seed(seed)
    cv_repeats <- rsample::vfold_cv(data, v = v, repeats = repeats)
  }
  n_total <- nrow(cv_repeats)
  repeat_labels <- if ("id2" %in% names(cv_repeats)) cv_repeats$id else rep("Repeat1", n_total)

  fold_args <- lapply(seq_len(n_total), function(i) {
    list(split = cv_repeats$splits[[i]], fold_num = i, repeat_label = repeat_labels[i])
  })

  run_one <- function(args) {
    set.seed(seed + args$fold_num)
    train <- rsample::analysis(args$split)
    test  <- rsample::assessment(args$split)
    res   <- fold_fun(train, test, args$fold_num)
    if (is.null(res)) return(NULL)

    # Identifiant : colonne id_col si présente, sinon NUMÉRO DE LIGNE dans
    # `data` (unique d'un fold à l'autre -- indispensable à check_outliers()).
    ids <- if (!is.null(id_col) && id_col %in% names(test)) {
      test[[id_col]]
    } else {
      as.character(rsample::complement(args$split))
    }
    preds <- tibble::tibble(
      ID = ids, Reel = res$Y, Prediction = res$pred, Residu = res$Y - res$pred,
      Fold = args$fold_num, Repeat = args$repeat_label
    )
    extras <- lapply(res$extras, function(tb) {
      tb <- tibble::as_tibble(tb)
      tb$Fold <- args$fold_num
      tb$Repeat <- args$repeat_label
      dplyr::relocate(tb, "Fold", "Repeat")
    })
    list(preds = preds, extras = extras)
  }

  n_cores_used <- .resolve_n_cores(parallelisation, n_cores, n_total)
  t0 <- Sys.time()
  results <- purrr::compact(.run_folds(fold_args, run_one, n_cores_used, pkgs))
  t1 <- Sys.time()

  n_valides <- length(results)
  if (n_valides == 0) {
    stop("Aucun fold n'a pu etre ajuste (facteur mono-niveau, colonnes incoherentes entre train/test, ",
         "ou echec de l'ajustement). Relancer avec parallelisation = FALSE pour voir les messages d'erreur.",
         call. = FALSE)
  }
  if (n_valides < n_total) {
    message(sprintf("%d fold(s) sur %d ignore(s) (facteur a un seul niveau, colonnes incoherentes ou echec d'ajustement).",
                    n_total - n_valides, n_total))
  }

  extra_names <- unique(unlist(lapply(results, function(r) names(r$extras))))
  extras <- stats::setNames(
    lapply(extra_names, function(nm) dplyr::bind_rows(lapply(results, function(r) r$extras[[nm]]))),
    extra_names
  )

  list(
    cv_repeats       = cv_repeats,
    fold_predictions = dplyr::bind_rows(lapply(results, `[[`, "preds")),
    extras           = extras,
    n_folds_valides  = n_valides,
    n_folds_total    = n_total,
    n_cores_used     = n_cores_used,
    parallelisation  = n_cores_used > 1L,
    time_folds_sec   = as.numeric(difftime(t1, t0, units = "secs"))
  )
}

#' Exécute les folds, séquentiellement ou sur un cluster PSOCK
#'
#' Cluster PSOCK ({parallel}, base R) : fonctionne à l'identique sous
#' Windows/macOS/Linux. Les workers chargent le namespace de cvperf (plus aucun
#' clusterExport() manuel). Si cvperf n'est pas INSTALLÉ (ex. devtools::load_all()),
#' retombe en séquentiel avec un avertissement.
#' @noRd
.run_folds <- function(fold_args, fun, n_cores_used = 1L, pkgs = character()) {
  if (n_cores_used <= 1L) return(lapply(fold_args, fun))

  message(sprintf("Parallelisation activee : %d coeur(s).", n_cores_used))
  cl <- parallel::makeCluster(n_cores_used)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  ok <- tryCatch({
    for (pk in unique(c("cvperf", pkgs))) parallel::clusterCall(cl, base::loadNamespace, pk)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    warning("Impossible de charger {cvperf} (ou une dependance) sur les workers paralleles : ",
            "le package doit etre INSTALLE (devtools::install(), pas seulement load_all()). ",
            "Execution sequentielle.", call. = FALSE)
    return(lapply(fold_args, fun))
  }

  tryCatch(
    parallel::parLapply(cl, fold_args, fun),
    error = function(e) {
      stop("Echec de l'execution parallele des folds : ", conditionMessage(e),
           "\nRelancer avec parallelisation = FALSE pour isoler le probleme.", call. = FALSE)
    }
  )
}

#' Chronomètre l'ajustement du modèle final
#' @noRd
.time_it <- function(expr) {
  t0 <- Sys.time()
  value <- force(expr)
  list(value = value, time_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

#' Constructeur commun des objets "*_cv_fit"
#'
#' Champs COMMUNS garantis pour toutes les méthodes (utilisés par toutes les
#' fonctions génériques) : data, formula, id_col, response_var, family,
#' classification_threshold, positive_class, cv_repeats, fold_predictions
#' (ID/Reel/Prediction/Residu/Fold/Repeat), fold_importance (Variable /
#' Importance / Fold [/ Used pour les méthodes creuses], à l'échelle de la
#' variable D'ORIGINE), fold_hyperparams, final_model, method, method_label,
#' importance_type, n_folds_valides, n_folds_total, v, repeats,
#' computational_info.
#'
#' @param cv sortie de run_repeated_cv()
#' @param ... champs spécifiques à la méthode
#' @param class classe spécifique (ex. "rf_cv_fit")
#' @keywords internal
#' @export
new_cv_fit <- function(cv, data, formula, id_col, family, classification_threshold, v, repeats,
                        fold_importance, fold_hyperparams, final_model, final_time_sec,
                        method, method_label, importance_type, ..., class) {
  response_var <- .response_var(formula)
  structure(
    c(
      list(
        data             = data,
        formula          = formula,
        id_col           = id_col,
        response_var     = response_var,
        family           = family,
        classification_threshold = classification_threshold,
        positive_class   = .positive_class(data[[response_var]], family),
        cv_repeats       = cv$cv_repeats,
        fold_predictions = cv$fold_predictions,
        fold_importance  = fold_importance,
        fold_hyperparams = fold_hyperparams,
        final_model      = final_model,
        method           = method,
        method_label     = method_label,
        importance_type  = importance_type,
        n_folds_valides  = cv$n_folds_valides,
        n_folds_total    = cv$n_folds_total,
        v                = v,
        repeats          = repeats,
        computational_info = list(
          time_folds_sec       = cv$time_folds_sec,
          time_final_model_sec = final_time_sec,
          time_total_sec       = cv$time_folds_sec + final_time_sec,
          parallelisation      = cv$parallelisation,
          n_cores_used         = cv$n_cores_used,
          n_folds              = cv$n_folds_valides
        )
      ),
      list(...)
    ),
    class = c(class, "cv_fit")
  )
}
