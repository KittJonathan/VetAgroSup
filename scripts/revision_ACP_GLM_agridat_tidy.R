# =============================================================================
#  MODULE "DÉMARCHE SCIENTIFIQUE" - VetAgro Sup, 3e année
#  SCRIPT COMPAGNON : équivalents tidyverse / easystats / rstatix / tidymodels
#  du script "revision_ACP_GLM_agridat.R" (R base)
# -----------------------------------------------------------------------------
#  La numérotation des sections est IDENTIQUE à celle du script R base, pour
#  pouvoir les projeter côte à côte. Les rappels théoriques ne sont pas répétés :
#  se reporter au script de base. Ici on se concentre sur la traduction du code.
#
#  Rôle de chaque famille de packages (à expliquer aux étudiants) :
#   - tidyverse (dplyr, tidyr, purrr, ggplot2, tibble) : manipuler, itérer, tracer
#   - broom         : transformer n'importe quel modèle en tableaux (tidy/glance/augment)
#   - easystats     : performance (qualité, diagnostics), parameters (coefficients),
#                     effectsize, modelbased (moyennes/prédictions marginales),
#                     correlation, datawizard, report, see (graphiques)
#   - rstatix       : tests statistiques "pipe-friendly" (cor_test, shapiro_test,
#                     levene_test, anova_test, mahalanobis_distance, ...)
#   - tidymodels    : cadre unifié MODÉLISATION PRÉDICTIVE (recipes, parsnip,
#                     workflows, rsample, yardstick, tune)
#
#  PLAN (identique au script de base)
#   0. Préparation
#   1. Variance-covariance et régression linéaire (matricielle)
#   2. Analyses multivariées : ACP, CAH, interaction G x E
#   3. GLM : Poisson, IRLS, surdispersion, binomial, Gamma, diagnostics
#   4. Fonction de diagnostic réutilisable + exercices
#
#  REMARQUE PÉDAGOGIQUE IMPORTANTE (démarche scientifique)
#   Inférence (comprendre, quantifier un effet, tester une hypothèse) :
#     -> lm/glm + broom + easystats + rstatix
#   Prédiction (prédire de nouvelles données, comparer des modèles par leur
#   performance hors-échantillon) : -> tidymodels (rééchantillonnage, métriques)
#   Les deux démarches utilisent les mêmes modèles mais répondent à des questions
#   différentes. Les blocs "TIDYMODELS" ci-dessous servent à illustrer cela.
#
#  Nécessite R >= 4.1 (pipe natif |> et fonctions anonymes \(x)).
#  Les sorties des packages easystats évoluent vite : en cas d'erreur d'argument,
#  consulter la documentation de la version installée (packageVersion()).
# =============================================================================


# =============================================================================
# 0. PRÉPARATION
# =============================================================================
pkgs <- c("agridat", "tidyverse", "broom", "easystats", "rstatix",
          "tidymodels", "poissonreg", "MASS")
a_installer <- setdiff(pkgs, rownames(installed.packages()))
if (length(a_installer) > 0) install.packages(a_installer)

# Packages optionnels
dispo <- function(p) requireNamespace(p, quietly = TRUE)

library(agridat)
library(tidyverse)     # dplyr, tidyr, purrr, ggplot2, tibble, readr, stringr, forcats
library(broom)         # tidy(), glance(), augment()
library(easystats)     # performance, parameters, effectsize, modelbased, correlation, ...
library(rstatix)
library(tidymodels)    # recipes, parsnip, workflows, rsample, yardstick, tune, ...
library(poissonreg)    # parsnip : poisson_reg()

# ATTENTION AUX CONFLITS DE NOMS
#  - MASS::select masque dplyr::select -> on N'ATTACHE PAS MASS, on appelle MASS::glm.nb
#  - tidymodels_prefer() donne priorité aux fonctions tidymodels en cas de conflit
tidymodels_prefer()

theme_set(theme_minimal(base_size = 12))
options(digits = 5, width = 100)
set.seed(2026)


# =============================================================================
# 1. VARIANCE-COVARIANCE ET RÉGRESSION LINÉAIRE MATRICIELLE
# =============================================================================
# Équivalent de la section 1 du script de base (lasrosas.corn, année 2001).

# --- Préparation des données (tibble + dplyr) --------------------------------
maize <- lasrosas.corn |>
  as_tibble() |>
  filter(year == 2001) |>
  select(yield, nitro, topo, bv) |>
  drop_na() |>
  mutate(topo = factor(topo))

glimpse(maize)                                   # équivalent moderne de str()
maize |> get_summary_stats(type = "common")      # rstatix : n, min, max, médiane, moyenne, sd, se, IC

# --- 1.1 Variance, covariance, corrélation ----------------------------------
# dplyr::summarise : tout dans un tableau
maize |>
  summarise(var_nitro = var(nitro),
            cov_nitro_yield = cov(nitro, yield),
            cor_nitro_yield = cor(nitro, yield))

# rstatix : test de corrélation avec IC et p-value dans un tibble
maize |> rstatix::cor_test(yield, nitro, method = "pearson")

# Matrice de corrélation : rstatix::cor_mat() + p-values
mat_cor <- maize |> select(where(is.numeric)) |> cor_mat()
mat_cor
mat_cor |> cor_get_pval()                        # matrice des p-values

# easystats : correlation::correlation() (tableau de paires, avec IC, p ajustées)
cor_es <- maize |> select(where(is.numeric)) |> correlation::correlation()
cor_es
cor_es |> summary() |> plot()                    # matrice de corrélation (package see)

# --- 1.2 Matrice de variance-covariance --------------------------------------
# L'algèbre matricielle reste la même : on convertit le tibble en matrice.
X_vars <- maize |> select(yield, nitro, bv) |> as.matrix()
Xc <- scale(X_vars, center = TRUE, scale = FALSE)
S_manuelle <- crossprod(Xc) / (nrow(Xc) - 1)     # crossprod(A) = t(A) %*% A
S_manuelle
all.equal(S_manuelle, cov(X_vars), check.attributes = FALSE)

# --- 1.3 Régression simple : formules scalaires avec dplyr -------------------
maize |>
  summarise(b1 = cov(nitro, yield) / var(nitro),
            b0 = mean(yield) - b1 * mean(nitro))

# --- 1.4 Résolution matricielle ----------------------------------------------
X <- model.matrix(~ nitro, data = maize)
y <- maize$yield
n <- nrow(X); p <- ncol(X)
beta_hat <- solve(crossprod(X), crossprod(X, y))  # (X'X)^-1 X'y
sigma2   <- sum((y - X %*% beta_hat)^2) / (n - p)
V_beta   <- sigma2 * solve(crossprod(X))

# Résultats rangés dans un tibble (même format que broom::tidy)
res_matriciel <- tibble(
  term      = colnames(X),
  estimate  = as.numeric(beta_hat),
  std.error = sqrt(diag(V_beta))
) |>
  mutate(statistic = estimate / std.error,
         p.value   = 2 * pt(abs(statistic), df = n - p, lower.tail = FALSE))
res_matriciel

# Comparaison avec lm() via broom : même structure, on peut empiler les deux
mod_lm <- lm(yield ~ nitro, data = maize)
bind_rows(matriciel = res_matriciel,
          lm        = tidy(mod_lm),
          .id = "methode") |>
  arrange(term, methode)

# broom::tidy  -> coefficients (avec IC)     broom::glance -> statistiques globales
# broom::augment -> données + ajustées, résidus, leviers (.hat), Cook (.cooksd)
tidy(mod_lm, conf.int = TRUE)
glance(mod_lm)                                    # r.squared, sigma, statistic (F), AIC, BIC, ...
aug <- augment(mod_lm)
aug

# easystats : tableaux de paramètres et de performance
model_parameters(mod_lm)                         # coefficients + IC + t + p + ddl
model_performance(mod_lm)                        # AIC, BIC, R2, R2 ajusté, RMSE, sigma
effectsize::effectsize(mod_lm)                   # coefficients standardisés
# report : génère une phrase de résultat (en anglais) -> à discuter, pas à recopier !
report::report(mod_lm)

# --- 1.5 Matrice chapeau, décomposition de la variance -----------------------
H <- X %*% solve(crossprod(X)) %*% t(X)
all.equal(diag(H), aug$.hat, check.attributes = FALSE)   # leviers = diag(H)

aug |>
  summarise(SCT = sum((yield - mean(yield))^2),
            SCE = sum((.fitted - mean(yield))^2),
            SCR = sum(.resid^2)) |>
  mutate(SCE_plus_SCR = SCE + SCR, R2 = SCE / SCT)

# --- 1.7 Régression multiple -------------------------------------------------
mod_mult <- lm(yield ~ nitro + bv + topo, data = maize)
model_parameters(mod_mult)
Xm <- model.matrix(yield ~ nitro + bv + topo, data = maize)
tibble(term = colnames(Xm),
       manuel = as.numeric(solve(crossprod(Xm), crossprod(Xm, maize$yield))),
       lm     = coef(mod_mult))

# --- 1.8 Hypothèses du modèle linéaire ---------------------------------------
# easystats : un seul appel pour tous les graphiques de diagnostic
check_model(mod_mult)
# ... ou test par test
check_normality(mod_mult)                        # normalité des résidus (Shapiro)
check_heteroscedasticity(mod_mult)               # test de Breusch-Pagan
check_collinearity(mod_mult)                     # VIF (colinéarité)
check_outliers(mod_mult)                         # observations influentes
# rstatix : tests "à la main" sur les résidus ou sur les groupes
augment(mod_mult) |> shapiro_test(.resid)
maize |> levene_test(yield ~ topo)               # homogénéité des variances entre groupes
maize |> anova_test(yield ~ topo)                # ANOVA à 1 facteur
maize |> tukey_hsd(yield ~ topo)                 # comparaisons multiples

# --- TIDYMODELS : inférence vs prédiction -----------------------------------
# Question prédictive : "quel modèle prédit le mieux le rendement d'une NOUVELLE
# parcelle ?" -> on ne juge plus sur le R² d'ajustement, mais sur des données
# non utilisées pour l'ajuster (validation croisée).
lm_spec <- linear_reg() |> set_engine("lm")      # parsnip : spécification du modèle

# (a) Un workflow = formule (ou recipe) + modèle
lm_wf  <- workflow() |> add_formula(yield ~ nitro + bv + topo) |> add_model(lm_spec)
lm_fit <- fit(lm_wf, data = maize)
tidy(lm_fit)                                     # mêmes coefficients que lm()
extract_fit_engine(lm_fit) |> coef()             # récupérer l'objet lm sous-jacent

# (b) Découpage apprentissage / test
split <- initial_split(maize, prop = 0.75, strata = topo)
fit_train <- fit(lm_wf, data = training(split))
predict(fit_train, new_data = testing(split)) |>
  bind_cols(testing(split)) |>
  metrics(truth = yield, estimate = .pred)       # yardstick : rmse, rsq, mae

# (c) Validation croisée + comparaison de plusieurs modèles
folds <- vfold_cv(maize, v = 5, repeats = 5, strata = topo)
mes_metriques <- metric_set(rmse, mae, rsq)
jeu_modeles <- workflow_set(
  preproc = list(nitro_seul = yield ~ nitro,
                 nitro_bv   = yield ~ nitro + bv,
                 complet    = yield ~ nitro + bv + topo),
  models  = list(lm = lm_spec)
) |>
  workflow_map("fit_resamples", resamples = folds, metrics = mes_metriques,
               seed = 2026, verbose = FALSE)
rank_results(jeu_modeles, rank_metric = "rmse") |>
  filter(.metric == "rmse") |>
  select(wflow_id, mean, std_err, rank)
autoplot(jeu_modeles, metric = "rmse")
# Discussion : un modèle plus complexe a TOUJOURS un meilleur R² d'ajustement ;
# la validation croisée pénalise le sur-ajustement (comparer avec AIC/BIC).


# =============================================================================
# 2. ANALYSES MULTIVARIÉES
# =============================================================================
# Jeu de données : australia.soybean, moyenne par génotype.
traits <- c("yield", "height", "lodging", "size", "protein", "oil")

soja_tbl <- australia.soybean |>
  as_tibble() |>
  group_by(gen) |>
  summarise(across(all_of(traits), \(v) mean(v, na.rm = TRUE))) |>
  drop_na() |>
  mutate(gen = as.character(gen))

# Version "data.frame avec noms de lignes" pour les fonctions de base/FactoMineR
soja <- soja_tbl |> as.data.frame() |> column_to_rownames("gen")
dim(soja)

# --- 2.1 Exploration préalable -----------------------------------------------
soja_tbl |> get_summary_stats(all_of(traits), type = "mean_sd")

# Distributions : pivot_longer + ggplot (tout en une chaîne)
soja_tbl |>
  pivot_longer(all_of(traits), names_to = "variable", values_to = "valeur") |>
  ggplot(aes(valeur)) +
  geom_histogram(bins = 15, fill = "steelblue", colour = "white") +
  facet_wrap(~ variable, scales = "free")

# Corrélations
soja |> cor_mat() |> cor_reorder() |> pull_lower_triangle() |> cor_plot()
if (dispo("GGally")) GGally::ggpairs(soja)

# Détection multivariée des valeurs atypiques (distance de Mahalanobis, rstatix)
soja |> mahalanobis_distance() |> filter(is.outlier)
# Normalité multivariée (test de Shapiro multivarié)
soja |> mshapiro_test()

# Standardisation (3 écritures équivalentes)
Z1 <- scale(soja)                                                  # base
Z2 <- soja |> mutate(across(everything(), \(v) as.numeric(scale(v))))  # tidyverse
Z3 <- datawizard::standardize(soja)                                # easystats
all.equal(as.numeric(Z1), as.numeric(as.matrix(Z2)))
all.equal(as.numeric(Z1), as.numeric(as.matrix(Z3)))

# --- 2.2 ACP "à la main" (version tibble) -----------------------------------
R <- cor(soja)
ed <- eigen(R)
p  <- ncol(soja); n <- nrow(soja)
lambda <- ed$values

tab_vp <- tibble(axe = paste0("PC", 1:p), valeur_propre = lambda) |>
  mutate(pct_inertie = 100 * valeur_propre / sum(valeur_propre),
         pct_cumule  = cumsum(pct_inertie))
tab_vp

# Éboulis avec ggplot2
ggplot(tab_vp, aes(axe, valeur_propre)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 1, linetype = 2, colour = "red") +   # critère de Kaiser
  labs(title = "Éboulis des valeurs propres", x = NULL, y = "Valeur propre")

# Corrélations variables-axes (tibble)
cor_var_tbl <- (ed$vectors %*% diag(sqrt(lambda))) |>
  `colnames<-`(paste0("PC", 1:p)) |>
  as_tibble() |>
  mutate(variable = colnames(soja), .before = 1)
cor_var_tbl

# Contributions (en %) des variables aux axes = vecteurs propres au carré
contr_var_tbl <- (100 * ed$vectors^2) |>
  `colnames<-`(paste0("PC", 1:p)) |>
  as_tibble() |>
  mutate(variable = colnames(soja), .before = 1)
contr_var_tbl |> select(variable, PC1, PC2) |> arrange(desc(PC1))

# --- 2.3 ACP : broom, easystats, tidymodels ---------------------------------
# (a) prcomp + broom : tidy() / augment() sur objet prcomp
acp <- prcomp(soja, center = TRUE, scale. = TRUE)
tidy(acp, matrix = "eigenvalues")                # std.dev, percent, cumulative
tidy(acp, matrix = "rotation")                   # vecteurs propres au format long
tidy(acp, matrix = "scores")
augment(acp, data = soja_tbl)                    # données + coordonnées .fittedPC1, ...

# Cercle des corrélations avec ggplot2
cercle <- tibble(t = seq(0, 2 * pi, length.out = 200), x = cos(t), y = sin(t))
ggplot() +
  geom_path(data = cercle, aes(x, y)) +
  geom_hline(yintercept = 0, linetype = 3) + geom_vline(xintercept = 0, linetype = 3) +
  geom_segment(data = cor_var_tbl,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), colour = "red") +
  geom_text(data = cor_var_tbl, aes(PC1 * 1.1, PC2 * 1.1, label = variable), colour = "red") +
  coord_fixed() + labs(title = "Cercle des corrélations", x = "PC1", y = "PC2")

# Nuage des individus
scores_tbl <- augment(acp, data = soja_tbl)
ggplot(scores_tbl, aes(.fittedPC1, .fittedPC2, label = gen)) +
  geom_point(colour = "darkgreen") +
  geom_text(vjust = -0.6, size = 2.5, check_overlap = TRUE) +
  geom_hline(yintercept = 0, linetype = 3) + geom_vline(xintercept = 0, linetype = 3) +
  labs(x = sprintf("PC1 (%.1f %%)", tab_vp$pct_inertie[1]),
       y = sprintf("PC2 (%.1f %%)", tab_vp$pct_inertie[2]),
       title = "Plan factoriel des génotypes")

# Version factoextra (ggplot2 sous le capot) si disponible
if (dispo("factoextra")) {
  print(factoextra::fviz_eig(acp, addlabels = TRUE))
  print(factoextra::fviz_pca_var(acp, col.var = "contrib", repel = TRUE))
  print(factoextra::fviz_pca_biplot(acp, repel = TRUE))
}

# (b) easystats : parameters::principal_components + choix du nombre d'axes
pca_es <- parameters::principal_components(soja, n = "auto")
pca_es                                          # chargements (loadings) avec "complexité"
summary(pca_es)                                 # variance expliquée par axe
predict(pca_es) |> head()                       # scores des individus
plot(parameters::principal_components(soja))    # nécessite le package see
# Nombre d'axes : consensus de plusieurs méthodes (Kaiser, éboulis, parallel analysis...)
if (dispo("psych") && dispo("nFactors")) {
  nf <- parameters::n_factors(soja)
  print(nf)
  plot(nf)
}

# (c) tidymodels : recipes::step_pca  (ACP = étape de prétraitement)
# C'est l'ACP "en tant que transformation" réutilisable sur de nouvelles données.
pca_rec <- recipe(~ ., data = soja_tbl) |>
  update_role(gen, new_role = "id") |>           # gen n'entre pas dans l'ACP
  step_normalize(all_numeric_predictors()) |>    # centrage-réduction (étape 1)
  step_pca(all_numeric_predictors(), num_comp = 6)  # ACP (étape 2)
pca_prep <- prep(pca_rec)                         # "apprend" moyennes, sd, axes

tidy(pca_prep, number = 2, type = "coef")        # vecteurs propres (loadings)
tidy(pca_prep, number = 2, type = "variance") |>  # valeurs propres, % inertie, cumul
  filter(terms %in% c("variance", "percent variance", "cumulative percent variance")) |>
  pivot_wider(names_from = terms, values_from = value)
scores_rec <- bake(pca_prep, new_data = NULL)    # scores des individus
scores_rec |> ggplot(aes(PC1, PC2, label = gen)) +
  geom_point() + geom_text(vjust = -0.6, size = 2.5, check_overlap = TRUE)
# Projection de NOUVEAUX individus avec les paramètres appris (= éléments suppl.)
nouveau <- soja_tbl |> slice(1) |> mutate(across(all_of(traits), \(v) v * 1.05))
bake(pca_prep, new_data = nouveau) |> select(gen, PC1, PC2)
# Variante : garder le nombre d'axes nécessaire pour 80 % de variance
#   step_pca(all_numeric_predictors(), threshold = 0.8)
# Les signes des axes peuvent différer entre prcomp, recipes et FactoMineR : sans importance.

# --- 2.4 Classification (CAH) -----------------------------------------------
nb_axes <- sum(lambda > 1)
coord <- scores_tbl |> select(gen, num_range(".fittedPC", 1:nb_axes)) |>
  column_to_rownames("gen")
hc <- coord |> dist() |> hclust(method = "ward.D2")
if (dispo("factoextra")) print(factoextra::fviz_dend(hc, k = 3, cex = 0.5))

groupes <- tibble(gen = hc$labels[hc$order], dummy = 1) |> select(gen)  # ordre de l'arbre
groupes <- tibble(gen = names(cutree(hc, 3)), classe = factor(cutree(hc, 3)))

# Caractérisation des classes : jointure + group_by/summarise
soja_tbl |>
  left_join(groupes, by = "gen") |>
  group_by(classe) |>
  summarise(n = n(), across(all_of(traits), \(v) round(mean(v), 2)))

# Visualiser les classes sur le plan factoriel
scores_tbl |>
  left_join(groupes, by = "gen") |>
  ggplot(aes(.fittedPC1, .fittedPC2, colour = classe)) +
  geom_point(size = 2) + labs(x = "PC1", y = "PC2", colour = "Classe")

# Test des différences entre classes pour chaque variable (rstatix)
soja_tbl |>
  left_join(groupes, by = "gen") |>
  pivot_longer(all_of(traits), names_to = "variable", values_to = "valeur") |>
  group_by(variable) |>
  anova_test(valeur ~ classe) |>
  adjust_pvalue(method = "holm")
# ATTENTION (point d'esprit critique) : les classes ont été construites pour
# maximiser les différences, donc ces tests sont circulaires : p-values
# descriptives uniquement, pas de véritable inférence.

# Approche tidymodels de la classification non supervisée : package {tidyclust}
if (dispo("tidyclust")) {
  library(tidyclust)
  km_fit <- k_means(num_clusters = 3) |> set_engine("stats") |>
    fit(~ ., data = soja)
  print(extract_cluster_assignment(km_fit) |> count(.cluster))
  hc_fit <- hier_clust(num_clusters = 3, linkage_method = "ward.D2") |>
    fit(~ ., data = as_tibble(scale(soja)))
  print(extract_cluster_assignment(hc_fit) |> count(.cluster))
}

# --- 2.5 Interaction génotype x milieu (version dplyr/tidyr) ----------------
wheat <- yan.winterwheat |> as_tibble() |> select(gen, env, yield)

# Double centrage en dplyr : interaction = y - moy_gen - moy_env + moy_générale
wheat_int <- wheat |>
  mutate(mu = mean(yield)) |>
  group_by(gen) |> mutate(m_gen = mean(yield)) |> ungroup() |>
  group_by(env) |> mutate(m_env = mean(yield)) |> ungroup() |>
  mutate(interaction = yield - m_gen - m_env + mu)

# Lien avec le modèle linéaire : l'interaction est le RÉSIDU du modèle additif !
mod_add <- lm(yield ~ gen + env, data = wheat)
all.equal(wheat_int$interaction, unname(resid(mod_add)))   # TRUE si tableau complet

# Passage au format large et SVD
M_int <- wheat_int |>
  pivot_wider(id_cols = gen, names_from = env, values_from = interaction) |>
  as.data.frame() |> column_to_rownames("gen") |> as.matrix()
sv <- svd(M_int)
tibble(axe = paste0("IPC", seq_along(sv$d)), pct = 100 * sv$d^2 / sum(sv$d^2)) |> head(4)

sc_gen <- tibble(nom = rownames(M_int), type = "génotype",
                 IPC1 = sv$u[, 1] * sqrt(sv$d[1]), IPC2 = sv$u[, 2] * sqrt(sv$d[2]))
sc_env <- tibble(nom = colnames(M_int), type = "environnement",
                 IPC1 = sv$v[, 1] * sqrt(sv$d[1]), IPC2 = sv$v[, 2] * sqrt(sv$d[2]))
bind_rows(sc_gen, sc_env) |>
  ggplot(aes(IPC1, IPC2, label = nom, colour = type)) +
  geom_hline(yintercept = 0, linetype = 3) + geom_vline(xintercept = 0, linetype = 3) +
  geom_text(size = 3) + labs(title = "Biplot de l'interaction G x E (AMMI)")


# =============================================================================
# 3. MODÈLES LINÉAIRES GÉNÉRALISÉS
# =============================================================================
# --- 3.1 Visualisation des lois avec expand_grid + ggplot -------------------
expand_grid(lambda = c(1, 4, 10), y = 0:25) |>
  mutate(proba = dpois(y, lambda)) |>
  ggplot(aes(y, proba)) +
  geom_col(width = 0.4, fill = "steelblue") +
  facet_wrap(~ lambda, labeller = label_both) +
  labs(title = "Loi de Poisson : moyenne = variance = lambda")

expand_grid(p = c(0.1, 0.5, 0.8), y = 0:20) |>
  mutate(proba = dbinom(y, 20, p)) |>
  ggplot(aes(y, proba)) + geom_col(width = 0.4, fill = "darkorange") +
  facet_wrap(~ p, labeller = label_both) +
  labs(title = "Loi binomiale (n = 20)")

expand_grid(shape = c(2, 8, 20), x = seq(0.5, 100, by = 0.5)) |>
  mutate(rate = shape / 20,                       # moyenne = shape/rate = 20 fixée
         densite = dgamma(x, shape = shape, rate = rate)) |>
  ggplot(aes(x, densite, colour = factor(shape))) + geom_line() +
  labs(title = "Lois Gamma de même moyenne (20)", colour = "forme")

# --- 3.2 GLM de Poisson (beall.webworms) -------------------------------------
ww <- beall.webworms |> as_tibble() |>
  mutate(across(c(block, spray, lead), factor))

# Moyenne vs variance par traitement : la variance croît-elle avec la moyenne ?
ww |>
  group_by(spray, lead) |>
  summarise(moyenne = mean(y), variance = var(y), rapport = variance / moyenne,
            .groups = "drop")
ww |> ggplot(aes(interaction(spray, lead), y)) + geom_boxplot() +
  labs(x = "spray.lead", y = "Nb de larves")

mod_pois <- glm(y ~ block + spray * lead, family = poisson, data = ww)

# broom : coefficients en rapports multiplicatifs avec IC
tidy(mod_pois, exponentiate = TRUE, conf.int = TRUE) |>
  filter(!str_detect(term, "^block"))
glance(mod_pois)                                  # null.deviance, deviance, AIC, df.residual...

# easystats : même information, mise en forme plus riche
model_parameters(mod_pois, exponentiate = TRUE)
model_performance(mod_pois)                       # AIC, BIC, R2 de Nagelkerke/Tjur, RMSE...
r2_mcfadden(mod_pois)

# Prédictions avec IC (échelle du lien puis retransformation) : broom::augment
nd <- expand_grid(block = factor(levels(ww$block)[1], levels = levels(ww$block)),
                  spray = levels(ww$spray), lead = levels(ww$lead))
augment(mod_pois, newdata = nd, type.predict = "link", se_fit = TRUE) |>
  mutate(mu = exp(.fitted),
         ic_inf = exp(.fitted - 1.96 * .se.fit),
         ic_sup = exp(.fitted + 1.96 * .se.fit))
# modelbased : moyennes marginales sur l'échelle de la réponse
# (selon la version : argument `by` (récent) ou `at` (ancien))
modelbased::estimate_means(mod_pois, by = c("spray", "lead"))

# Tests : analyse de déviance en tibble + test du rapport de vraisemblance
anova(mod_pois, test = "Chisq") |> tidy()
mod_pois_add <- update(mod_pois, . ~ . - spray:lead)
test_likelihoodratio(mod_pois_add, mod_pois)      # performance : LRT / comparaison

# Surdispersion et excès de zéros (performance)
check_overdispersion(mod_pois)
check_zeroinflation(mod_pois)
check_model(mod_pois)                             # panneau complet de diagnostics

# Résidus simulés (équivalent DHARMa) intégrés à performance >= 0.12
if (packageVersion("performance") >= "0.12.0" && dispo("DHARMa")) {
  res_sim <- simulate_residuals(mod_pois)
  print(check_residuals(res_sim))
  print(check_overdispersion(res_sim))
}

# --- TIDYMODELS : Poisson via parsnip + validation croisée ------------------
pois_spec <- poisson_reg() |> set_engine("glm")
pois_wf   <- workflow() |> add_formula(y ~ block + spray * lead) |> add_model(pois_spec)
pois_fit  <- fit(pois_wf, data = ww)
tidy(pois_fit) |> filter(!str_detect(term, "^block"))   # = coef(mod_pois)

# Validation croisée ALÉATOIRE vs GROUPÉE par bloc :
# les parcelles d'un même bloc ne sont pas indépendantes. Une validation croisée
# aléatoire est trop optimiste ; la validation par blocs répond à la question
# "le modèle prédit-il un NOUVEAU bloc ?".
metriques_comptage <- metric_set(rmse, mae)
folds_alea  <- vfold_cv(ww, v = 5, repeats = 3)
folds_bloc  <- group_vfold_cv(ww, group = block, v = 5)
res_alea <- fit_resamples(pois_wf, folds_alea, metrics = metriques_comptage)
res_bloc <- fit_resamples(pois_wf, folds_bloc, metrics = metriques_comptage)
bind_rows(aleatoire = collect_metrics(res_alea),
          par_bloc  = collect_metrics(res_bloc), .id = "validation")
# Variante : ajouter la perte log-Poisson (yardstick >= 1.2) :
#   metric_set(rmse, mae, poisson_log_loss)
# Remarque : la binomiale négative n'est pas dans tous les cadres parsnip ;
# vérifier la doc de poissonreg (negative_binomial()) selon votre version, sinon
# conserver MASS::glm.nb avec broom/easystats (cf. 3.4).

# --- 3.3 IRLS en style fonctionnel (purrr::accumulate) -----------------------
# Chaque itération est une fonction pure b_{k+1} = f(b_k) ; accumulate() garde
# toute la trajectoire, ce qui permet de VISUALISER la convergence.
Xww <- model.matrix(y ~ block + spray * lead, data = ww)
yww <- ww$y

irls_step <- function(b, X, y) {
  eta <- as.numeric(X %*% b); mu <- exp(eta)
  z <- eta + (y - mu) / mu                        # variable de travail
  w <- mu                                         # poids (Poisson / log)
  as.numeric(solve(crossprod(X, X * w), crossprod(X, w * z)))
}
dev_pois <- function(y, mu) {
  mu <- as.numeric(mu)
  2 * sum(ifelse(y == 0, 0, y * log(y / mu)) - (y - mu))
}
b0 <- as.numeric(solve(crossprod(Xww), crossprod(Xww, log(yww + 0.5))))  # départ : MCO sur log(y+0.5)
trajet <- accumulate(1:8, \(b, i) irls_step(b, Xww, yww), .init = b0)

convergence <- tibble(iteration = 0:8, beta = trajet) |>
  mutate(deviance  = map_dbl(beta, \(b) dev_pois(yww, exp(Xww %*% b))),
         ecart_glm = map_dbl(beta, \(b) max(abs(b - coef(mod_pois)))))
convergence |> select(-beta)
ggplot(convergence, aes(iteration, deviance)) + geom_line() + geom_point() +
  labs(title = "Convergence de l'IRLS", y = "Déviance")
deviance(mod_pois)                                # valeur cible

# --- 3.4 Surdispersion : quasi-Poisson et binomiale négative ----------------
mod_qp <- update(mod_pois, family = quasipoisson)
mod_nb <- MASS::glm.nb(y ~ block + spray * lead, data = ww)

# Comparaison des erreurs standard des trois modèles dans un seul tableau
bind_rows(Poisson = tidy(mod_pois), QuasiPoisson = tidy(mod_qp), NegBin = tidy(mod_nb),
          .id = "modele") |>
  filter(term %in% c("sprayY", "(Intercept)") | str_detect(term, "^lead")) |>
  select(modele, term, estimate, std.error) |>
  pivot_wider(names_from = modele, values_from = c(estimate, std.error))

# Comparaison globale AIC / BIC / R2 / RMSE (quasi-Poisson n'a pas d'AIC)
compare_performance(mod_pois, mod_nb, rank = TRUE)
check_overdispersion(mod_nb)
model_parameters(mod_nb, exponentiate = TRUE) |> as_tibble() |> slice(1:5)

# --- 3.5 GLM binomial (crowder.seeds) ---------------------------------------
cs <- crowder.seeds |> as_tibble() |>
  mutate(across(c(gen, extract), factor), prop = germ / n)
cs |> ggplot(aes(interaction(gen, extract), prop)) + geom_boxplot() + geom_jitter(width = .1) +
  labs(x = "gen.extract", y = "Proportion germée")

mod_bin <- glm(cbind(germ, n - germ) ~ gen * extract, family = binomial, data = cs)
tidy(mod_bin, exponentiate = TRUE, conf.int = TRUE)      # odds ratios
model_parameters(mod_bin, exponentiate = TRUE)
check_overdispersion(mod_bin)                             # phi > 1 : surdispersion
binned_residuals(mod_bin) |> plot()                       # résidus "binnés" (spécifique binomial)
modelbased::estimate_means(mod_bin, by = c("gen", "extract"))   # probabilités prédites

# Comparaison avec l'approche "classique" rstatix sur les proportions (ANOVA) :
# ne tient pas compte des effectifs n_i ni de la variance p(1-p)/n
cs |> anova_test(prop ~ gen * extract)
# vs GLM quasi-binomial (corrige la surdispersion)
mod_qb <- update(mod_bin, family = quasibinomial)
anova(mod_qb, test = "F") |> tidy()

# GLMM avec effet aléatoire de la boîte de Petri (ouverture)
if (dispo("lme4")) {
  mod_glmm <- lme4::glmer(cbind(germ, n - germ) ~ gen * extract + (1 | plate),
                          family = binomial, data = cs)
  print(model_parameters(mod_glmm, exponentiate = TRUE))
  print(compare_performance(mod_bin, mod_glmm))
}

# --- TIDYMODELS : régression logistique ------------------------------------
# parsnip::logistic_reg attend une réponse binaire (un facteur) : on "désagrège"
# les données (1 ligne = 1 graine) avec tidyr::uncount().
cs_bin <- cs |>
  mutate(echec = n - germ) |>
  pivot_longer(c(germ, echec), names_to = "issue", values_to = "k") |>
  uncount(k) |>
  mutate(issue = factor(issue, levels = c("germ", "echec")))   # 1er niveau = événement pour yardstick
count(cs_bin, gen, extract, issue)

logit_spec <- logistic_reg() |> set_engine("glm")
logit_wf   <- workflow() |> add_formula(issue ~ gen * extract) |> add_model(logit_spec)
logit_fit  <- fit(logit_wf, cs_bin)
tidy(logit_fit, exponentiate = TRUE)
# Point pédagogique : parsnip considère le 1er niveau du facteur comme l'événement ;
# l'objet glm sous-jacent modélise lui le 2e niveau. Comparer avec mod_bin :
# les coefficients peuvent être de signe opposé selon l'ordre des niveaux.
extract_fit_engine(logit_fit) |> coef()
coef(mod_bin)

# Validation croisée GROUPÉE par boîte de Petri (unité expérimentale) :
folds_plate <- group_vfold_cv(cs_bin, group = plate, v = 5)
res_logit <- fit_resamples(logit_wf, folds_plate,
                           metrics = metric_set(roc_auc, mn_log_loss, accuracy))
collect_metrics(res_logit)
# Attention : sur données désagrégées en 0/1, la surdispersion n'est plus
# détectable (check_overdispersion n'a plus de sens). La structure en boîtes de
# Petri doit être traitée explicitement (GLMM) ou par un rééchantillonnage groupé.

# --- 3.6 GLM Gamma (lasrosas.corn) ------------------------------------------
maize2 <- maize   # rendement, nitro, topo déjà prêts (cf. section 1)
mod_gauss <- glm(yield ~ nitro + topo, family = gaussian, data = maize2)
mod_gamma <- glm(yield ~ nitro + topo, family = Gamma(link = "log"), data = maize2)

model_parameters(mod_gamma, exponentiate = TRUE)          # effets multiplicatifs
compare_performance(mod_gauss, mod_gamma, rank = TRUE)
check_model(mod_gamma)

# tidymodels : le moteur glm accepte une famille personnalisée
gamma_spec <- linear_reg() |> set_engine("glm", family = Gamma(link = "log"))
gamma_wf   <- workflow() |> add_formula(yield ~ nitro + topo) |> add_model(gamma_spec)
gauss_wf   <- workflow() |> add_formula(yield ~ nitro + topo) |>
  add_model(linear_reg() |> set_engine("lm"))
folds_m <- vfold_cv(maize2, v = 5, repeats = 3)
bind_rows(gaussien = fit_resamples(gauss_wf, folds_m, metrics = metric_set(rmse, mae)) |> collect_metrics(),
          gamma_log = fit_resamples(gamma_wf, folds_m, metrics = metric_set(rmse, mae)) |> collect_metrics(),
          .id = "modele") |>
  select(modele, .metric, mean, std_err)
# Discussion : un AIC meilleur (inférence) n'implique pas forcément un meilleur
# RMSE en validation croisée (prédiction) ; les critères répondent à des questions différentes.

# --- 3.7 Diagnostics et sélection de modèle ---------------------------------
mod_ww0 <- MASS::glm.nb(y ~ block, data = ww)
mod_ww1 <- MASS::glm.nb(y ~ block + spray + lead, data = ww)
mod_ww2 <- MASS::glm.nb(y ~ block + spray * lead, data = ww)
compare_performance(mod_ww0, mod_ww1, mod_ww2, rank = TRUE)  # AIC, BIC, R2, RMSE + score global
test_likelihoodratio(mod_ww0, mod_ww1, mod_ww2)
# Observations influentes : broom::augment fournit .cooksd et .hat
augment(mod_ww2) |> slice_max(.cooksd, n = 5) |> select(y, spray, lead, .fitted, .hat, .cooksd)
check_outliers(mod_ww2)
# Résidus de déviance vs valeurs ajustées avec ggplot2
augment(mod_ww2) |>
  mutate(res_dev = residuals(mod_ww2, type = "deviance")) |>
  ggplot(aes(.fitted, res_dev)) + geom_point(alpha = .5) + geom_hline(yintercept = 0) +
  geom_smooth(se = FALSE) + labs(x = "Valeurs ajustées (échelle du lien)", y = "Résidus de déviance")


# =============================================================================
# 4. FONCTION DE DIAGNOSTIC RÉUTILISABLE + EXERCICES
# =============================================================================
# Fonction à donner (ou à faire écrire) aux étudiants : résume un GLM en un seul
# appel avec les outils easystats + broom. Montre l'intérêt de la programmation
# fonctionnelle en R pour la reproductibilité.
diagnostic_glm <- function(modele, exponentiate = FALSE) {
  list(
    parametres  = model_parameters(modele, exponentiate = exponentiate),
    performance = model_performance(modele),
    surdisp     = tryCatch(check_overdispersion(modele), error = \(e) NA),
    zeros       = tryCatch(check_zeroinflation(modele), error = \(e) NA),
    influents   = augment(modele) |> slice_max(.cooksd, n = 3)
  )
}
diagnostic_glm(mod_pois, exponentiate = TRUE)

# Appliquer la même fonction à plusieurs modèles avec purrr
liste_modeles <- list(poisson = mod_pois, quasi = mod_qp, negbin = mod_nb)
map(liste_modeles, \(m) tidy(m)$std.error[1:3])

# --- Exercices (en miroir de ceux du script R base) --------------------------
#  1. Refaire l'exercice A (régression matricielle) avec un autre jeu agridat,
#     comparer (X'X)^-1 X'y, lm + broom::tidy et parsnip::linear_reg.
#  2. ACP : reproduire l'ACP de soja avec prcomp + broom, puis avec
#     recipes::step_pca ; vérifier que les valeurs propres coïncident.
#     Discuter : dans quel cas préférer la version "recipe" ? (réutiliser les
#     paramètres de normalisation et d'ACP sur de nouvelles données)
#  3. GLM : sur webworms, comparer Poisson et NB (a) par AIC (easystats),
#     (b) par validation croisée groupée par bloc (tidymodels). Les deux
#     critères conduisent-ils au même choix ? Pourquoi ?
#  4. Écrire une fonction qui, pour un modèle donné, renvoie un tibble
#     « terme, estimation, IC, p » avec rapports multiplicatifs si lien log/logit.
#  5. Rédiger 5 lignes distinguant démarche INFÉRENTIELLE et PRÉDICTIVE sur un
#     des exemples (quelle question ? quels outils ? quels critères ?).
# =============================================================================

sessionInfo()
