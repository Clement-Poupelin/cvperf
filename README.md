# cvperf

**Performance de modèles prédictifs en validation croisée répétée** : Lasso / Elastic Net, forêt aléatoire, XGBoost, SVM et MARS, avec **une seule interface commune** pour les métriques, les diagnostics, les graphiques, l'importance des variables, les valeurs de Shapley et la comparaison de méthodes.

Le package regroupe les toolkits `*_performance_toolkit.R`, `shap_toolkit.R` et `model_comparison_toolkit.R`. Les fonctions communes, auparavant dupliquées dans chaque toolkit (≈ 60 % du code), n'existent plus qu'en un exemplaire.

## Installation

```r
# depuis la racine du package (le dossier qui contient DESCRIPTION)
devtools::document()
devtools::install()
# ou depuis le dossier parent : devtools::install("cvperf")
# méthodes utilisées (installer celles dont on a besoin)
install.packages(c("glmnet", "ranger", "xgboost", "e1071", "earth", "fastshap", "car", "DT"))
```

> ⚠ **Installer** le package (et pas seulement `devtools::load_all()`) si vous utilisez `parallelisation = TRUE` : les processus parallèles chargent `cvperf` depuis la bibliothèque installée. Avec `load_all()`, l'exécution retombe en séquentiel avec un avertissement.

## Démarrage rapide

```r
library(cvperf)
data(diabetes, package = "lars")
df <- data.frame(Y = diabetes$y, unclass(diabetes$x))

# Mêmes folds pour toutes les méthodes -> comparaison appariée possible
set.seed(42)
cv <- rsample::vfold_cv(df, v = 5, repeats = 20)

fit_lasso <- fit_lasso_repeated_cv(df, Y ~ ., id_col = "ID", cv_repeats = cv)
fit_rf    <- fit_rf_repeated_cv(df,    Y ~ ., id_col = "ID", cv_repeats = cv, parallelisation = TRUE)
fit_mars  <- fit_mars_repeated_cv(df,  Y ~ ., id_col = "ID", cv_repeats = cv)

fit_rf                          # résumé
model_performance(fit_rf)       # R², RMSE, MAE poolés et par fold
check_model(fit_rf)             # panel complet de diagnostics
plot_variable_importance(fit_rf)
plot_effects(fit_mars)
check_hyperparam_stability(fit_mars)
plot_shap(check_shap(fit_rf, nsim = 30))

cmp <- compare_methods(list(Lasso = fit_lasso, RF = fit_rf, MARS = fit_mars))
cmp$plot; cmp$leaderboard
compare_importance(list(Lasso = fit_lasso, RF = fit_rf, MARS = fit_mars))$plot
```

## Architecture

```
R/
├── cv_engine.R            moteur de CV répétée unique : run_repeated_cv(), new_cv_fit()
├── data_prep.R            préparation des données (une seule fonction, tolérance aux NA en option)
├── generics.R             crochets S3 que chaque méthode implémente + print() commun
├── performance.R          model_performance()
├── diagnostics_residuals.R check_normality/heteroscedasticity/outliers/independence
├── plots_prediction.R     plot_prediction_quality(), résidus
├── plots_classification.R ROC, précision-rappel, calibration, confusion
├── importance.R           check_importance(), plot_importance(), plot_variable_importance(), group_importance()
├── permutation.R          folds internes + importance par permutation (SVM, XGBoost)
├── effects.R              plot_effects() : dépendance partielle, toutes méthodes
├── check_model.R          panel check_model()
├── shap.R                 check_shap(), plot_shap() : un seul algorithme pour toutes les méthodes
├── compare.R              compare_models(), compare_methods(), compare_* ...
├── method_lasso.R  method_rf.R  method_xgb.R  method_svm.R  method_mars.R
└── compat.R               anciens noms (plot_roc_rf(), check_model_xgb()...) -- généré
```

Chaque `fit_*_repeated_cv()` renvoie un objet de classe `c("<méthode>_cv_fit", "cv_fit")`. Ses champs communs sont garantis : `fold_predictions` (ID / Reel / Prediction / Residu / Fold / Repeat), `fold_importance` (Variable / Importance / Fold, à l'échelle de la variable d'origine), `fold_hyperparams`, `final_model`, `cv_repeats`, `computational_info`, etc. Toute la couche commune ne lit que ces champs, plus quelques génériques S3 :

| Crochet | Rôle | Obligatoire |
|---|---|---|
| `model_design(fit)` | X, Y, `predict(newdata)` et correspondance colonne → variable du modèle final (modes « final », `plot_effects()`, SHAP) | **oui** |
| `importance_table(fit)` | table d'importance affichée | non (défaut : moyenne sur les folds) |
| `method_notes(fit, topic)` | avertissements d'interprétation propres à la méthode | non |
| `check_hyperparam_stability()`, `plot_hyperparam_path()`, `plot_hyperparam_stability()` | réglage des hyperparamètres | non |
| `plot_variable_importance()`, `variables_summary()` | lecture spécifique (sélection Lasso/MARS) | non (défaut : top-K) |

**Ajouter une méthode** revient à écrire `R/method_xxx.R` : une fonction `fold_fun(train, test, fold_num)` passée à `run_repeated_cv()`, l'ajustement du modèle final, un appel à `new_cv_fit()` et une méthode `model_design.xxx_cv_fit()`. Tout le reste (diagnostics, SHAP, comparaison) fonctionne immédiatement. Un exemple complet figure dans `tutoriels/tutoriel_model_comparison_toolkit.qmd`, section « Étendre à une méthode supplémentaire ».

## Anciens noms → noms génériques

Les anciens noms restent disponibles (`R/compat.R`), les tutoriels n'ont donc eu besoin que de remplacer `source(...)` par `library(cvperf)`. Pour du nouveau code :

| Toolkit | Package |
|---|---|
| `check_normality_rf(fit)`, `plot.check_normality_rf(x)` | `check_normality(fit)`, `plot(x)` |
| `plot_roc_xgb()`, `plot_calibration_svm()`, `check_model_mars()`... | `plot_roc()`, `plot_calibration()`, `check_model()` |
| `plot_prediction_quality_rf()`... | `plot_prediction_quality()` |
| `check_hyperparam_stability_rf()`, `check_hyperparam_stability_xgb()`... | `check_hyperparam_stability()` |
| `plot_lambda_path_lasso()`, `plot_hyperparam_path_svm()`... | `plot_hyperparam_path()` |
| `plot_variable_importance_rf()`, `plot_coefficient_importance()` | `plot_variable_importance()` |
| `plot_effects_mars()`, `plot_effects_svm()`... | `plot_effects()` |
| `importance_table_xgb()`, `coef_table_lasso()` | `importance_table()` |
| `check_importance_rf()`, `group_importance_svm()` | `check_importance()`, `group_importance()` |
| `compare_models_lasso()`, `compare_models_rf()` | `compare_models()` |
| `check_shap(fit, method = "rf")` | `check_shap(fit)` (argument `method` ignoré) |

## Changements de comportement par rapport aux toolkits

- **Reproductibilité** : chaque fold est précédé de `set.seed(seed + numéro de fold)` pour toutes les méthodes. Le Lasso donne donc désormais le même résultat en séquentiel et en parallèle ; ses résultats peuvent différer légèrement des anciennes sorties.
- **Identifiant par défaut** : sans colonne `id_col`, l'identifiant est le numéro de ligne dans `data` (et non plus dans le fold de test, qui produisait des doublons faussant `check_outliers()`).
- **Forêt aléatoire** : `fit$positive_class` contient désormais le vrai libellé de la classe positive (comme les autres méthodes), et plus `"1"`.
- **Réponse manquante** : une valeur manquante dans la réponse arrête l'ajustement avec un message clair (le Lasso supprimait auparavant ces lignes en silence).
- **Lasso** : les prédicteurs sont restreints aux variables de la formule, comme pour les autres méthodes.
- **MARS** : `plot_effects()` attend des noms de variables d'origine (et non plus de colonnes de `model.matrix()`).
- **Comparaison** : `compare_variable_selection()` détecte les méthodes creuses par la colonne `Used` de `fold_importance` ; une méthode inconnue est traitée comme dense (top-K).
- `compare_roc()` / `compare_calibration()` requièrent **ggplot2 ≥ 3.5.0**.
- `model_performance()`, `check_model()`, `check_normality()`, `check_heteroscedasticity()` et `check_outliers()` portent les mêmes noms que dans {performance} : si les deux packages sont chargés, le dernier chargé masque l'autre (utiliser `cvperf::` au besoin).

## Tutoriels

Les six tutoriels Quarto sont dans `tutoriels/` (exclus du build du package), prêts pour le site : ils chargent `library(cvperf)` au lieu de sourcer les toolkits.

## Statut

Version 0.1.0 : le code a été restructuré sans pouvoir être exécuté. Avant publication :

```r
devtools::document()   # régénère NAMESPACE et man/ depuis les balises roxygen
devtools::test()
devtools::check()
quarto::quarto_render("tutoriels/")
```
