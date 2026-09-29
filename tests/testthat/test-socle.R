test_that("preparation des donnees : encodage binaire et correspondance des variables", {
  df <- data.frame(Y = factor(c("a", "b", "a", "b")), X1 = 1:4, SEXE = factor(c("F", "M", "M", "F")))
  p <- cvperf:::.prepare_matrix_data(df, Y ~ ., id_col = "ID", family = "binomial")
  expect_equal(p$Y, c(0, 1, 0, 1))
  expect_equal(colnames(p$X), c("X1", "SEXEM"))
  expect_equal(cvperf:::.map_to_original_variable(colnames(p$X), c("X1", "SEXE")), c("X1", "SEXE"))
})

test_that("moteur de CV : structure commune et identifiants uniques", {
  set.seed(1)
  df <- data.frame(Y = rnorm(40), X1 = rnorm(40), X2 = rnorm(40))
  fold_fun <- function(train, test, fold_num) {
    m <- lm(Y ~ ., data = train)
    list(Y = test$Y, pred = unname(predict(m, test)),
         extras = list(importance = tibble::tibble(Variable = c("X1", "X2"), Importance = abs(coef(m)[-1]))))
  }
  cv <- run_repeated_cv(df, Y ~ ., id_col = "ID", family = "gaussian", v = 4, repeats = 2, seed = 1,
                        cv_repeats = NULL, fold_fun = fold_fun)
  expect_equal(cv$n_folds_valides, 8)
  expect_named(cv$fold_predictions, c("ID", "Reel", "Prediction", "Residu", "Fold", "Repeat"))
  expect_equal(length(unique(cv$fold_predictions$ID)), 40)   # une ligne = un identifiant
})

test_that("Lasso de bout en bout (si glmnet est installe)", {
  skip_if_not_installed("glmnet")
  set.seed(1)
  df <- data.frame(Y = rnorm(60), X1 = rnorm(60), X2 = rnorm(60))
  df$Y <- df$Y + 2 * df$X1
  fit <- fit_lasso_repeated_cv(df, Y ~ ., id_col = "ID", v = 3, repeats = 2)
  expect_s3_class(fit, c("lasso_cv_fit", "cv_fit"))
  expect_true(model_performance(fit)$pooled$R2 > 0.5)
  expect_s3_class(check_model(fit), "patchwork")
})
