# =============================================================================
# data_prep.R -- préparation des données, COMMUNE à toutes les méthodes
# =============================================================================
# Remplace .make_model_matrix() (Lasso), .prepare_xgb_data(), .prepare_svm_data(),
# .prepare_mars_data() -- qui ne différaient que par leur tolérance aux NA --
# et .prepare_rf_data() (data.frame natif pour {ranger}).
#
# Convention de "classe positive" (réponse binaire), identique partout :
#   facteur   -> le 2e niveau (ordre de levels()) devient 1 ;
#   numérique -> la plus grande des 2 valeurs devient 1.

#' Nom de la variable réponse d'une formule
#' @noRd
.response_var <- function(formula) all.vars(formula)[1]

#' Variables CANDIDATES (colonnes du data.frame d'origine proposées comme
#' prédicteurs), que la formule soit `Y ~ .` ou explicite
#' @param fit objet (ou liste) avec $formula, $data, $response_var, $id_col
#' @noRd
.vars_candidates <- function(fit) {
  vc <- all.vars(fit$formula)[-1]
  if ("." %in% vc) {
    vc <- setdiff(names(fit$data), c(fit$response_var %||% .response_var(fit$formula), fit$id_col))
  }
  vc
}

#' Ramène un nom de colonne "étalé" par model.matrix() (ex. "SEXEM") à sa
#' variable d'origine (ex. "SEXE"). Candidats les plus longs testés d'abord.
#' @noRd
.map_to_original_variable <- function(variable_names, vars_candidates) {
  vt <- vars_candidates[order(-nchar(vars_candidates))]
  purrr::map_chr(variable_names, function(v) {
    direct <- vt[vt == v]
    if (length(direct) > 0) return(direct[1])
    prefixe <- vt[startsWith(v, vt)]
    if (length(prefixe) > 0) return(prefixe[1])
    v
  })
}

#' Encode la réponse : numérique (gaussian) ou 0/1 (binomial) ; NULL si la
#' réponse binaire n'a pas exactement 2 valeurs sur ce sous-ensemble
#' @noRd
.encode_response <- function(y_raw, family) {
  if (identical(family, "binomial")) {
    if (is.factor(y_raw)) {
      y_raw <- droplevels(y_raw)
      if (nlevels(y_raw) != 2) return(NULL)
      return(as.numeric(y_raw) - 1)
    }
    vu <- unique(stats::na.omit(y_raw))
    if (length(vu) != 2) return(NULL)
    return(as.numeric(y_raw == max(vu)))
  }
  as.numeric(y_raw)
}

#' Libellé de la classe traitée comme "1"
#' @noRd
.positive_class <- function(y_raw, family) {
  if (!identical(family, "binomial")) return(NA_character_)
  if (is.factor(y_raw)) levels(droplevels(y_raw))[2] else as.character(max(unique(stats::na.omit(y_raw))))
}

#' Retire l'identifiant et restreint aux variables de la formule (sauf `Y ~ .`)
#' @noRd
.restrict_df <- function(df, formula, id_col) {
  response_var <- .response_var(formula)
  df <- as.data.frame(df)
  if (!is.null(id_col) && id_col %in% names(df)) df <- df[, setdiff(names(df), id_col), drop = FALSE]
  rhs <- all.vars(formula)[-1]
  if (!("." %in% rhs)) df <- df[, unique(c(response_var, rhs)), drop = FALSE]
  df
}

#' TRUE si un prédicteur facteur/caractère n'a qu'un niveau sur ce sous-ensemble
#' @noRd
.has_single_level_factor <- function(df, response_var) {
  preds <- setdiff(names(df), response_var)
  any(vapply(preds, function(nm) {
    x <- df[[nm]]
    (is.factor(x) || is.character(x)) && length(unique(stats::na.omit(as.character(x)))) < 2
  }, logical(1)))
}

#' Matrice de design X + réponse numérique Y (Lasso, XGBoost, SVM, MARS)
#'
#' Facteurs encodés par model.matrix() (référence simple), colonne
#' "(Intercept)" retirée, backticks retirés des noms de colonnes. Aucune
#' ligne n'est supprimée silencieusement (na.pass) : une réponse manquante
#' fait échouer, un prédicteur manquant aussi sauf si `allow_na_x = TRUE`.
#'
#' @param allow_na_x TRUE pour les méthodes qui gèrent nativement les NA (xgboost)
#' @param pkg nom du package d'ajustement (pour le message d'erreur)
#' @return list(X, Y) ou NULL (facteur mono-niveau, réponse binaire mal formée)
#' @noRd
.prepare_matrix_data <- function(df, formula, id_col, family = "gaussian",
                                 allow_na_x = FALSE, pkg = "ce modele") {
  response_var <- .response_var(formula)
  df_model <- .restrict_df(df, formula, id_col)
  if (.has_single_level_factor(df_model, response_var)) return(NULL)

  mf     <- stats::model.frame(formula, data = df_model, na.action = stats::na.pass)
  X_full <- stats::model.matrix(formula, data = mf)
  X <- X_full[, colnames(X_full) != "(Intercept)", drop = FALSE]
  colnames(X) <- gsub("`", "", colnames(X), fixed = TRUE)
  rownames(X) <- NULL
  storage.mode(X) <- "double"

  Y <- .encode_response(df_model[[response_var]], family)
  if (is.null(Y)) return(NULL)
  if (anyNA(Y)) {
    stop("Valeurs manquantes dans la reponse : retirer les lignes concernees avant l'ajustement.", call. = FALSE)
  }
  if (!isTRUE(allow_na_x) && anyNA(X)) {
    stop(sprintf("Valeurs manquantes dans les predicteurs : {%s} n'en tolere aucune. Imputer ou retirer les lignes concernees.", pkg),
         call. = FALSE)
  }
  list(X = X, Y = Y)
}

#' data.frame prêt pour {ranger} (facteurs natifs, réponse factor "0"/"1" en binomial)
#' @return list(df_fit, Y) ou NULL
#' @noRd
.prepare_df_data <- function(df, formula, id_col, family = "gaussian") {
  response_var <- .response_var(formula)
  df_model <- .restrict_df(df, formula, id_col)
  if (.has_single_level_factor(df_model, response_var)) return(NULL)
  Y <- .encode_response(df_model[[response_var]], family)
  if (is.null(Y)) return(NULL)
  if (anyNA(Y)) {
    stop("Valeurs manquantes dans la reponse : retirer les lignes concernees avant l'ajustement.", call. = FALSE)
  }
  if (identical(family, "binomial")) {
    df_model[[response_var]] <- factor(Y, levels = c(0, 1), labels = c("0", "1"))
  } else {
    df_model[[response_var]] <- Y
  }
  list(df_fit = df_model, Y = Y)
}

#' Construit une formule "réponse ~ variables" sûre (backticks partout)
#'
#' @param variables vecteur de noms (ex. \code{coef_table_lasso(fit)$VariableOriginale})
#' @param response nom de la réponse
#' @return formule
#' @export
build_reduced_formula <- function(variables, response) {
  variables <- unique(variables[!is.na(variables) & nzchar(variables)])
  if (length(variables) == 0) {
    stop("Aucune variable valide fournie a build_reduced_formula() -- verifier que le modele a retenu au moins une variable.",
         call. = FALSE)
  }
  stats::reformulate(paste0("`", variables, "`"), response = response)
}
