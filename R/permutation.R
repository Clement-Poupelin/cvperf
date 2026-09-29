# =============================================================================
# permutation.R -- briques partagées par les méthodes réglées par CV INTERNE
# (SVM, XGBoost, MARS) : folds internes, importance par permutation, groupes
# =============================================================================

#' Folds internes (stratifiés sur la classe si binomial), dans le jeu
#' d'ENTRAÎNEMENT du fold externe -- le fold de test externe n'intervient jamais
#' dans le réglage
#' @return liste d'indices TEST (un vecteur par fold interne)
#' @noRd
.make_inner_folds <- function(Y, k, family, seed) {
  n <- length(Y)
  k <- max(2L, min(as.integer(k), n %/% 2L))
  set.seed(seed)
  if (identical(family, "binomial")) {
    foldid <- integer(n)
    for (cl in unique(Y)) {
      idx <- which(Y == cl)
      foldid[idx] <- sample(rep_len(seq_len(k), length(idx)))
    }
  } else {
    foldid <- sample(rep_len(seq_len(k), n))
  }
  lapply(seq_len(k), function(f) which(foldid == f))
}

#' Résout perm_groups (variables d'origine) en groupes de colonnes de X
#' @noRd
.resolve_perm_groups <- function(perm_groups, var_map) {
  if (is.null(perm_groups) || length(perm_groups) == 0) return(list())
  stopifnot("perm_groups doit etre une liste NOMMEE de vecteurs de variables" =
              is.list(perm_groups) && !is.null(names(perm_groups)) && all(nzchar(names(perm_groups))))
  vars_ok <- unique(unname(var_map))
  purrr::imap(perm_groups, function(vars, nom) {
    inconnues <- setdiff(vars, vars_ok)
    if (length(inconnues) > 0) {
      stop(sprintf("perm_groups[['%s']] : variable(s) inconnue(s) : %s", nom, paste(inconnues, collapse = ", ")),
           call. = FALSE)
    }
    names(var_map)[unname(var_map) %in% vars]
  })
}

#' Importance par permutation sur des données NON VUES (fold de test)
#'
#' Importance(v) = perte(v permutée) - perte(intacte), moyennée sur n_perm
#' permutations ; perte = RMSE (gaussien) ou Brier (binomial). Toutes les
#' colonnes indicatrices d'un facteur sont permutées ENSEMBLE. Les groupes
#' sont permutés d'un bloc.
#' @return list(vars = tibble(Variable, Importance), groups = tibble(Variable, Importance))
#' @noRd
.perm_importance <- function(predict_fun, X, Y, family, var_map, groups_cols, n_perm = 20, seed = 1) {
  set.seed(seed)
  loss <- if (identical(family, "gaussian")) function(p) sqrt(mean((p - Y)^2)) else function(p) mean((p - Y)^2)
  base <- loss(predict_fun(X))
  n <- nrow(X)
  delta <- function(cols) {
    mean(vapply(seq_len(n_perm), function(r) {
      Xp <- X
      Xp[, cols] <- X[sample.int(n), cols]
      loss(predict_fun(Xp)) - base
    }, numeric(1)))
  }
  vars <- unique(unname(var_map))
  list(
    vars = tibble::tibble(Variable = vars,
                          Importance = vapply(vars, function(v) delta(names(var_map)[unname(var_map) == v]), numeric(1))),
    groups = if (length(groups_cols) > 0) {
      tibble::tibble(Variable = names(groups_cols), Importance = vapply(groups_cols, delta, numeric(1)))
    } else {
      tibble::tibble(Variable = character(), Importance = numeric())
    }
  )
}
