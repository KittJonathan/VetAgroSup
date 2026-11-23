# Tutoriel régression linéaire
# 2026-11-12

# Table des matières ----

# 0. Préparation : packages, données, conventions
# 1. Exploration des données
# 2. Covariance
# 3. Corrélation (Pearson, Spearman, Kendall, partielle, paradoxe de Simpson)
# 4. Régression linéaire simple - x quantitatif (calcul "à la main", puis lm)
# 5. Régression linéaire simple - x qualitatif (codages, contrastes, ANOVA, t-test)
# 6. Régression linéaire multiple (additif, interaction, colinéarité, sélection)
# 7. Présupposés du modèle linéaire (diagnostics, violations simulées, solutions)
# 8. Aperçu tidymodels (approche "prédictive" : train/test, validation croisée)
# 9. Modèles linéaires généralisés (GLM)
#     9.1 Principe
#     9.2 Aparté sur les distributions et les liens
#     9.3 Logistique
#     9.4 Poisson / sur-dispersion / zéros
#     9.5 Gamma
#     9.6 Beta
# 10. Modèles mixtes (LMM puis GLMM)
# 11. Check-list finale

# 0. Préparation ----

## 0.1 Packages ----

# Installer les packages nécessaires

pkgs <- c("tidyverse", "palmerpenguins", "easystats", 
          "tidymodels", "broom", "lme4", "lmerTest", 
          "glmmTMB", "MASS", "patchwork", "emmeans",
          "sandwich", "lmtest", "car", "nlme")

to_install <- setdiff(pkgs, rownames(installed.packages()))

if (length(to_install) > 0) install.packages(to_install)

rm(pkgs, to_install)

# Charger les packages

library(tidyverse)       # dplyr, ggplot2, tidyr, purrr, stringr, lubridate...
library(tidymodels)      # rsample, recipes, parsnip, workflows, yardstick...
library(easystats)       # performance, parameters, effectsize, modelbased,
                         # correlation, datawizard, report, see, insight...
library(lme4)            # modèles mixtes (lmer, glmer)
library(lmerTest)        # p-values (Satterthwaite) pour lmer ; masque lme4::lmer
library(glmmTMB)         # GLMM très flexibles (NB, Beta, zero-inflation...)
library(patchwork)       # assembler des ggplot : p1 + p2, p1 / p2

# ATTENTION aux conflits de noms : MASS::select masque dplyr::select.
# On n'attache donc JAMAIS MASS, on l'appelle toujours avec MASS::

# Thème ggplot
theme_set(theme_bw())

# Palette de couleurs
my_cols <- c(Adelie = "#FF8C00", Chinstrap = "#A034F0", Gentoo = "#159090")

# Nombre de chiffres à afficher
options(digits = 5)

## 0.2 Données ----

# ATTENTION : depuis R 4.5, un jeu `penguins` existe aussi dans {datasets}, 
# avec d'autres noms de colonnes (bill_len, ...). 
# On appelle donc explicitement le package palmerpenguins.

penguins_raw <- palmerpenguins::penguins

# Les données brutes comportent : 
#   - 344 individus
#   - 3 espèces
#   - 3 îles
#   - 3 années

# Variables : 
#   - species
#   - island
#   - bill_length_mm
#   - bill_depth_mm
#   - flipper_length_mm
#   - body_mass_g
#   - sex
#   - year

# Données manquantes : 2 individus sans aucune mesure, 11 sans sexe
colSums(is.na(penguins_raw))

# Pour ne pas compliquer les premières parties, on travaille 
# sur les cas complets.
# (En vrai projet : documenter les NA, comprendre POURQUOI ils manquent. 
# Si c'est lié à la variable réponse, supprimer les lignes biaise l'analyse).

penguins <- penguins_raw |> 
  drop_na()

nrow(penguins)

# peng <- peng_all |>
#   drop_na() |>
#   mutate(mass_100g = body_mass_g / 100,                       # unité plus parlante
#          mass_c    = mass_100g - mean(mass_100g))               # centrée
# nrow(peng)   # 333 individus

# 1. Exploration ----

glimpse(penguins)
skimr::skim(penguins)

datawizard::describe_distribution(penguins)
datawizard::data_codebook(penguins)

penguins |> count(species, island)  # plan d'échantillonnage déséquilibré !
penguins |> count(species, sex)     # sex ratio ~ 1:1 dans chaque espèce

# Toujours regarder les données AVANT de modéliser !

p_hist <- ggplot(penguins, aes(x = body_mass_g, fill = species)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity",
                 show.legend = FALSE) +
  scale_fill_manual(values = my_cols) +
  labs(title = "Distribution de la masse corporelle",
       x = "Masse (g)")

p_scatter <- ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g)) +
  geom_point(aes(color = species), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  scale_color_manual(values = my_cols) +
  labs(title = "Masse ~ longueur de nageoire",
       x = "Nageoire (mm)",
       y = "Masse (g)")

p_hist + p_scatter

# 2. Covariance ----

# La covariance mesure comment deux variables varient ENSEMBLE.

# cov(X, Y) = Σ (xi - x̄)(yi - ȳ) / (n - 1)
#   > 0 : quand X est au-dessus de sa moyenne, Y l'est en général aussi
#   < 0 : relation inverse
#   ≈ 0 : pas de relation LINÉAIRE

# La variance est un cas particulier : var(X) = cov(X, X)

x <- penguins$flipper_length_mm
y <- penguins$body_mass_g
n <- length(x)

## 2.1 Calcul "à la main" ----

cov_manuel <- sum((x - mean(x)) * (y - mean(y))) / (n - 1)
cov_manuel

cov(x, y)   

all.equal(cov_manuel, cov(x, y))

# La variance est la covariance d'une variable avec elle-même

all.equal(var(x), cov(x, x))

## 2.2 Visualisation ----

# Les quatre quadrants autour du point moyen :
# Chaque point contribue par le signe du produit (xi - x̄)(yi - ȳ).

ggplot() +
  geom_segment(aes(x = -2, xend = 2, y = 0, yend = 0)) +
  geom_segment(aes(x = 0, xend = 0, y = -2, yend = 2)) +
  annotate(geom = "text", x = 1.5, y = 1.5, label = "I") +
  annotate(geom = "text", x = 1, y = 1, label = "(+, +)") +
  annotate(geom = "text", x = -1.5, y = 1.5, label = "II") +
  annotate(geom = "text", x = -1, y = 1, label = "(-, +)") +
  annotate(geom = "text", x = -1.5, y = -1.5, label = "III") +
  annotate(geom = "text", x = -1, y = -1, label = "(-, -)") +
  annotate(geom = "text", x = 1.5, y = -1.5, label = "IV") +
  annotate(geom = "text", x = 1, y = -1, label = "(+, -)") +
  labs(x = "x", y = "y")


penguins |>
  mutate(produit = (flipper_length_mm - mean(flipper_length_mm)) *
                   (body_mass_g - mean(body_mass_g)),
         contribution = if_else(produit > 0, "positive", "négative")) |>
  ggplot(aes(flipper_length_mm, body_mass_g, color = contribution)) +
  geom_point(alpha = .7) +
  geom_vline(xintercept = mean(x), linetype = "dashed") +
  geom_hline(yintercept = mean(y), linetype = "dashed") +
  scale_color_manual(values = c(positive = "#2c7fb8", négative = "#d95f0e")) +
  labs(title = "Chaque point contribue + ou - à la covariance",
       subtitle = "Majorité de points 'positifs' => covariance positive")

## 2.3 Limite : la covariance dépend des unités ----

cov(x, y)         # mm * g  -> 9852.2
cov(x / 10, y)    # cm * g  -> 985.22
cov(x, y / 1000)  # mm * kg -> 9.8522

# => impossible de comparer des covariances entre paires de variables.
# D'où la corrélation, qui est sans unité.

## 2.4 Matrice de variance-covariance ----

penguins |> 
  select(where(is.numeric), -year) |> 
  cov()

# sur la diagonale : variances

# Cette matrice est l'objet central de la statistique multivariée ;
# on la retrouvera dans les modèles mixtes (matrice de covariance 
# des observations).

# 3. Corrélation ----

# coefficient de corrélation de Pearson : 
# r = cov(X, Y) / (sd(X) * sd(Y)) : covariance des variables standardisées.

# Toujours dans [-1, 1], sans unité, symétrique.

## 3.1 Pearson à la main ----

r_manuel <- cov(x, y) / (sd(x) * sd(y))
r_manuel
cor(x, y, method = "pearson")

# Corrélation = covariance des variables centrées-réduites
cov(scale(x), scale(y))

# Test d'hypothèse du coefficient de corrélation
# H0 : rho = 0
# test t à n-2 degrés de liberté

t_stat <- r_manuel * sqrt((n - 2) / (1 - r_manuel^2))
p_val <- 2 * pt(-abs(t_stat), df = n - 2)
c(t = t_stat, p = p_val)

cor_test <- cor.test(x, y)
cor_test$statistic
cor_test$p.value

# Intervalle de confiance à 95% par la transformation de Fisher
# z = atanh(r), se = 1 / sqrt(n - 3)

ic_fisher <- tanh(atanh(r_manuel) + c(-1, 1) * qnorm(0.975) / sqrt(n - 3))
ic_fisher

cor_test$conf.int


# r² = part de variance de Y "expliquée" linéairement par X (cf. section 4)
r_manuel^2

## 3.2 Spearman et Kendall : corrélations de rangs ----

# Spearman = Pearson sur les rangs ; détecte toute relation MONOTONE
# et est robuste aux valeurs extrêmes.

# Kendall = concordance de paires.

cor(rank(x), rank(y))
cor(x, y, method = "spearman")
cor(x, y, method = "kendall")

## 3.3 Les pièges de la corrélation ----

# (a) Quartet d'Anscombe : mêmes moyennes, variances ET corrélation,
#     mais structures très différentes => TOUJOURS visualiser les données.

anscombe_long <- datasets::anscombe |> 
  pivot_longer(everything(), 
               names_to = c(".value", "jeu"),
               names_pattern = "(.)(.)") 

anscombe_long |> 
  summarise(moy_x = mean(x), moy_y = mean(y), r = cor(x, y),
            .by = jeu)

anscombe_long |> 
  ggplot(aes(x, y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~jeu)

# (b) Relation parfaite mais non linéaire : r = 0

set.seed(1)
nl <- tibble(x = runif(200, -3, 3), y = x^2)
cor(nl$x, nl$y)                # ~ 0 alors que y est fonction de x !

# (c) Paradoxe de Simpson : relation inversée quand on ignore un facteur

cor(penguins$bill_length_mm, penguins$bill_depth_mm)  # globalement négative ...

penguins |> 
  summarise(r = cor(bill_length_mm, bill_depth_mm),   # ... positive par espèce
            .by = species)                        

penguins |> 
  ggplot(aes(x = bill_length_mm, bill_depth_mm)) +
  geom_point(aes(color = species)) +
  geom_smooth(method = "lm", color = "black") +                   # pente globale
  geom_smooth(aes(color = species), method = "lm", se = FALSE) +  # pentes intra-espèces
  scale_color_manual(values = my_cols) +
  labs(title = "Paradoxe de Simpson",
       subtitle = "Négatif globalement, positif au sein de chaque espèce")

# Leçon : une corrélation brute peut être produite par une variable de confusion.
# La régression multiple est l'outil pour "contrôler" ce genre d'effet.

## 3.4 Avec easystats::correlation() ----

penguins_num <- penguins |> 
  select(bill_length_mm:body_mass_g)

cor_tab <- correlation(penguins_num)               # Pearson par défaut, avec IC et p
cor_tab                                            # p-values ajustées (Holm) par défaut
summary(cor_tab)                                   # format matriciel

correlation(penguins_num, method = "spearman")
correlation(penguins_num, partial = TRUE)          # corrélations PARTIELLES (autres variables contrôlées)


## 3.4 Avec easystats::correlation ----
cor_tab <- correlation(peng_num)                 # Pearson par défaut, avec IC et p
cor_tab                                          # p-values ajustées (Holm) par défaut
summary(cor_tab)                                 # format matriciel
correlation(peng_num, method = "spearman")
correlation(peng_num, partial = TRUE)            # corrélations PARTIELLES (autres var. contrôlées)

# Corrélation par groupe

penguins |> 
  group_by(species) |> 
  correlation(select = "bill_length_mm", select2 = "bill_depth_mm")

# Heatmap maison

cor(penguins_num) |> 
  as_tibble(rownames = "v1") |> 
  pivot_longer(-v1, names_to = "v2", values_to = "r") |> 
  ggplot(aes(v1, v2, fill = r)) +
  geom_tile() +
  geom_text(aes(label = round(r, 2))) +
  scale_fill_gradient2(limits = c(-1, 1), low = "#d95f0e", mid = "white", high = "#2c7fb8") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# 4. Régression linéaire simple, x quantitatif ----

# Modèle :  y_i = β0 + β1 * x_i + ε_i ,  ε_i ~ N(0, σ²) indépendants

# Question : de combien la masse varie-t-elle, en moyenne, pour +1 mm de nageoire ?

## 4.1 Estimation par moindres carrés ordinaires (MCO), à la main ----

# Principe : choisir β0, β1 qui minimisent la somme des carrés des résidus SCE.

# Solution analytique :   
#     β1 = cov(x, y) / var(x)      
#     β0 = ȳ - β1 * x̄

beta1 <- cov(x, y) / var(x)
beta0 <- mean(y) - beta1 * mean(x)
c(beta0 = beta0, beta1 = beta1)

# Lien avec la corrélation : β1 = r * sd(y)/sd(x)

cor(x, y) * sd(y) / sd(x)

# Valeurs ajustées, résidus, décomposition de la variance

y_hat <- beta0 + beta1 * x
res <- y - y_hat

SCE <- sum(res^2)                 # somme des carrés résiduelle (non expliquée)
SCT <- sum((y - mean(y))^2)       # somme des carrés totale
SCR <- sum((y_hat - mean(y))^2)   # somme des carrés expliquée par le modèle

all.equal(SCT, SCR + SCE)         # SCT = SCR + SCE

R2 <- SCR / SCT                   # = r² en régression simple

all.equal(R2, cor(x, y)^2)

sigma2 <- SCE / (n - 2)           # estimation de σ² (2 paramètres estimés)
sqrt(sigma2)                      # "residual standard error" (en g)

# Erreur-type de β1 : σ / sqrt(Σ(xi - x̄)²)

Sxx   <- sum((x - mean(x))^2)
se_b1 <- sqrt(sigma2 / Sxx)
t_b1  <- beta1 / se_b1
p_b1  <- 2 * pt(-abs(t_b1), df = n - 2)
ic_b1 <- beta1 + c(-1, 1) * qt(.975, n - 2) * se_b1
c(se = se_b1, t = t_b1, p = p_b1); ic_b1

# F global = (SCR/1) / (SCE/(n-2)) ; en simple, F = t²

F_stat <- (SCR / 1) / (SCE / (n - 2))
c(F = F_stat, t_carre = t_b1^2)

## 4.2 Avec lm() ----

m1 <- lm(body_mass_g ~ flipper_length_mm, data = peng)
summary(m1)                       # retrouver chaque nombre calculé ci-dessus
# Interprétation : +1 mm de nageoire <=> +49,7 g de masse en moyenne.
# L'intercept (-5780 g) = masse prédite pour une nageoire de 0 mm : sans sens
# biologique (extrapolation hors du domaine). On centre donc x (cf. 4.5).

## 4.3 Version matricielle (la vraie formule générale) ----
# β̂ = (X'X)⁻¹ X'y ; Var(β̂) = σ² (X'X)⁻¹ . Valable pour TOUS les lm (multiple aussi).
X <- model.matrix(m1)             # colonne de 1 + flipper
head(X)
XtX_inv  <- solve(crossprod(X))
beta_mat <- XtX_inv %*% crossprod(X, y)
beta_mat
se_mat   <- sqrt(diag(sigma2 * XtX_inv))
se_mat

# Matrice chapeau H = X (X'X)⁻¹ X' ; ses diagonales = leviers (influence de chaque
# observation sur sa propre valeur ajustée).
H <- X %*% XtX_inv %*% t(X)
all.equal(unname(diag(H)), unname(hatvalues(m1)))

## 4.4 Sorties "propres" : broom et easystats ----
tidy(m1, conf.int = TRUE)         # coefficients (tibble)
glance(m1)                        # R², AIC, sigma, F... (1 ligne)
augment(m1) |> head()             # données + .fitted, .resid, .hat, .cooksd...

model_parameters(m1)              # coefficients + IC + t + p
model_performance(m1)             # AIC, BIC, R², R² ajusté, RMSE, sigma
report(m1)                        # texte automatique pour un rapport (à relire !)

## 4.5 Centrer / standardiser pour interpréter ----
# Centrer x : l'intercept devient la masse prédite pour une nageoire MOYENNE.
m1c <- lm(body_mass_g ~ I(flipper_length_mm - mean(flipper_length_mm)), data = peng)
coef(m1c)                         # intercept ≈ mean(y), pente inchangée
mean(y)

# Standardiser x ET y : la pente devient la corrélation (en régression simple).
standardize_parameters(m1, method = "refit")
cor(x, y)

## 4.6 Prédiction : intervalle de confiance vs intervalle de prédiction ----
# IC (confidence)   : incertitude sur la DROITE (moyenne conditionnelle).
# IP (prediction)   : incertitude sur une NOUVELLE observation (inclut σ), plus large.
new <- tibble(flipper_length_mm = seq(170, 235, by = 1))
conf <- predict(m1, new, interval = "confidence") |> as_tibble()
pred <- predict(m1, new, interval = "prediction") |> as_tibble()

ggplot() +
  geom_point(data = peng, aes(flipper_length_mm, body_mass_g), alpha = .4) +
  geom_ribbon(aes(x = new$flipper_length_mm, ymin = pred$lwr, ymax = pred$upr),
              alpha = .15, fill = "grey40") +
  geom_ribbon(aes(x = new$flipper_length_mm, ymin = conf$lwr, ymax = conf$upr),
              alpha = .4, fill = "#2c7fb8") +
  geom_line(aes(x = new$flipper_length_mm, y = conf$fit), linewidth = 1) +
  labs(title = "Droite de régression, IC 95 % (bleu) et intervalle de prédiction (gris)")

# Avec modelbased
estimate_relation(m1, by = "flipper_length_mm", length = 30) |> plot()

# Les 3 sources de dispersion à bien distinguer : variabilité des données (σ),
# incertitude sur les coefficients (SE), incertitude sur une prédiction (les deux).


# 5. RÉGRESSION LINÉAIRE SIMPLE, X QUALITATIF -----------------------------
# Une variable qualitative est transformée en variables indicatrices (dummy).
# Par défaut dans R : contraste "treatment" (la 1re modalité = référence).

## 5.1 Que fait lm() avec un facteur ? ----
levels(peng$species)
contrasts(peng$species)           # matrice de codage
m2 <- lm(body_mass_g ~ species, data = peng)
head(model.matrix(m2))            # 1 | Chinstrap | Gentoo
summary(m2)

# Interprétation :
#   (Intercept)      = moyenne du groupe de référence (Adelie)
#   speciesChinstrap = écart Chinstrap - Adelie
#   speciesGentoo    = écart Gentoo - Adelie
peng |> group_by(species) |> summarise(moyenne = mean(body_mass_g)) |>
  mutate(ecart_vs_Adelie = moyenne - moyenne[1])

# Le modèle ne fait que comparer des MOYENNES, et suppose (cf. section 7) :
# résidus normaux, variance identique dans tous les groupes.

## 5.2 Changer le codage change l'interprétation, pas le modèle ----
# (a) Autre référence
peng_rel <- peng |> mutate(species = relevel(species, ref = "Gentoo"))
coef(lm(body_mass_g ~ species, data = peng_rel))

# (b) Sans intercept : un coefficient = une moyenne de groupe
coef(lm(body_mass_g ~ 0 + species, data = peng))

# (c) Contrastes "somme" : intercept = moyenne des moyennes (grand mean),
#     coefficients = écart de chaque groupe à cette moyenne générale.
m2_sum <- lm(body_mass_g ~ species, data = peng,
             contrasts = list(species = contr.sum))
coef(m2_sum)
# IMPORTANT en présence d'interactions : le codage influence l'interprétation
# des effets "principaux" et les tests de type III. Le piège classique.

# Quel que soit le codage : mêmes prédictions, même R², même F global.
all.equal(fitted(m2), fitted(m2_sum))

## 5.3 Test global du facteur : ANOVA ----
anova(m2)                         # F de l'effet espèce (1 seul test, 2 df)
# ANOVA à un facteur = régression linéaire sur variable qualitative.
aov(body_mass_g ~ species, data = peng) |> summary()
eta_squared(m2)                   # taille d'effet : part de variance expliquée

## 5.4 Comparaisons deux à deux (post-hoc) ----
# Après un F significatif, quelles espèces diffèrent ? Il faut corriger la
# multiplicité des tests (Tukey, Holm...).
estimate_means(m2, by = "species")                              # moyennes estimées + IC
estimate_contrasts(m2, contrast = "species", p_adjust = "tukey")
emmeans::emmeans(m2, pairwise ~ species)                        # alternative classique

## 5.5 Cas à 2 modalités : t-test = régression ----
m_sex <- lm(body_mass_g ~ sex, data = peng)
summary(m_sex)$coefficients
t.test(body_mass_g ~ sex, data = peng, var.equal = TRUE)       # même t (signe près), même p
# F = t² :
c(F = anova(m_sex)$`F value`[1], t2 = summary(m_sex)$coefficients[2, "t value"]^2)
# Welch (variances inégales) = t.test par défaut, sans var.equal = TRUE.


# 6. RÉGRESSION LINÉAIRE MULTIPLE ----------------------------------------
# y = β0 + β1 x1 + β2 x2 + ... + ε. Chaque βj est un effet "toutes choses égales
# par ailleurs" : variation de y pour +1 de xj, les AUTRES variables étant FIXÉES.

## 6.1 Modèle additif ----
m_a <- lm(body_mass_g ~ flipper_length_mm, data = peng)
m_b <- lm(body_mass_g ~ flipper_length_mm + bill_length_mm, data = peng)
m_c <- lm(body_mass_g ~ flipper_length_mm + bill_length_mm + species, data = peng)
model_parameters(m_c)

# Notez comme le coefficient d'une variable change selon les covariables incluses
compare_parameters(m_a, m_b, m_c)

## 6.2 Même calcul matriciel qu'en simple ----
X <- model.matrix(m_c)
y <- peng$body_mass_g
beta <- solve(crossprod(X)) %*% crossprod(X, y)
cbind(manuel = beta[, 1], lm = coef(m_c))

## 6.3 Interprétation de "toutes choses égales par ailleurs" : Frisch-Waugh-Lovell ----
# Le coefficient de x1 en régression multiple = pente de la régression de la part
# de y NON expliquée par les autres variables sur la part de x1 NON expliquée
# par les autres variables.
e_y <- resid(lm(body_mass_g ~ bill_length_mm + species, data = peng))
e_x <- resid(lm(flipper_length_mm ~ bill_length_mm + species, data = peng))
coef(lm(e_y ~ e_x))["e_x"]
coef(m_c)["flipper_length_mm"]     # identique

## 6.4 Retour sur le paradoxe de Simpson ----
coef(lm(bill_depth_mm ~ bill_length_mm, data = peng))                     # pente < 0
coef(lm(bill_depth_mm ~ bill_length_mm + species, data = peng))           # pente > 0
# Ajouter l'espèce "contrôle" la confusion : on compare des individus de même espèce.

## 6.5 Colinéarité ----
# Quand deux prédicteurs sont très corrélés, leurs coefficients sont instables
# (écarts-types gonflés) même si le modèle prédit bien.
cor(peng$flipper_length_mm, peng$body_mass_g)
check_collinearity(m_c)            # VIF : < 5 OK ; 5-10 modéré ; > 10 problème
# VIF_j = 1 / (1 - R²_j) où R²_j = R² de la régression de xj sur les autres
r2_j <- summary(lm(flipper_length_mm ~ bill_length_mm + species, data = peng))$r.squared
1 / (1 - r2_j)
# Pistes : retirer/regrouper des variables, ACP, régularisation (ridge/lasso)...
plot(check_collinearity(m_c))

## 6.6 Interaction : l'effet de x dépend d'un facteur ----
m_d <- lm(body_mass_g ~ flipper_length_mm * species, data = peng)
model_parameters(m_d)
# flipper_length_mm            = pente pour la référence (Adelie)
# flipper_length_mm:speciesX   = DIFFÉRENCE de pente entre X et Adelie
# speciesX                     = différence d'ordonnée quand flipper = 0 (hors domaine !)
# => centrer la variable continue rend ces "effets principaux" interprétables.

# Pentes par espèce :
estimate_slopes(m_d, trend = "flipper_length_mm", by = "species")
estimate_relation(m_d, by = c("flipper_length_mm", "species")) |> plot()

# Test de l'interaction (modèles emboîtés) :
anova(m_c, m_d)                    # F-test de comparaison
# Types de sommes de carrés (matière à débat !) :
car::Anova(m_d, type = 2)          # type II : effets principaux sans interaction
# type III : exige des contrastes "somme" (options(contrasts = c("contr.sum","contr.poly")))

## 6.7 Coefficients standardisés et taille d'effet ----
# Pour comparer l'importance relative de variables d'unités différentes.
standardize_parameters(m_c, method = "refit")
# (Ne pas sur-interpréter : "importance" ≠ effet causal.)

## 6.8 Comparer des modèles ----
# R² ajusté, AIC/BIC (parcimonie), test F pour modèles emboîtés, validation croisée.
compare_performance(m_a, m_b, m_c, m_d, rank = TRUE)
plot(compare_performance(m_a, m_b, m_c, m_d, rank = TRUE))
anova(m_a, m_b, m_c, m_d)          # séquence de modèles emboîtés
# AIC = -2 logLik + 2k : pénalise la complexité ; ΔAIC > 2 : différence notable.
# Règle : éviter de "pêcher" des variables par sélection automatique (step) sans
# cadre ; formuler les modèles à partir des hypothèses biologiques.


# 7. PRÉSUPPOSÉS DU MODÈLE LINÉAIRE ---------------------------------------
# Moyen mnémotechnique "LINE" :
#   L inéarité      : E[y|x] est linéaire en les coefficients (et la forme de x correcte)
#   I ndépendance   : les erreurs sont indépendantes (plan d'échantillonnage !)
#   N ormalité      : les ERREURS (donc les résidus) suivent une loi normale
#   E galité des variances (homoscédasticité) : Var(ε) constante
# + absence de colinéarité forte, pas de points aberrants/influents excessifs.
#
# Importance pratique : I (la plus grave, non détectable par les graphiques,
# vient du protocole) > L > E > N. La normalité porte sur les RÉSIDUS (pas sur y,
# ni sur les x) et concerne surtout les tests/IC en petits échantillons (TCL).

m_diag <- lm(body_mass_g ~ flipper_length_mm + species, data = peng)

## 7.1 Les 4 graphiques classiques, construits à la main ----
aug <- augment(m_diag)

g1 <- ggplot(aug, aes(.fitted, .resid)) +
  geom_point(alpha = .5) + geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
  labs(title = "Résidus vs ajustées", subtitle = "Linéarité (courbe plate) et homoscédasticité (bande régulière)")

g2 <- ggplot(aug, aes(sample = .std.resid)) +
  stat_qq() + stat_qq_line() +
  labs(title = "QQ-plot des résidus standardisés", subtitle = "Normalité : points sur la droite")

g3 <- ggplot(aug, aes(.fitted, sqrt(abs(.std.resid)))) +
  geom_point(alpha = .5) + geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
  labs(title = "Scale-location", subtitle = "Homoscédasticité : courbe plate")

g4 <- ggplot(aug, aes(.hat, .std.resid, size = .cooksd)) +
  geom_point(alpha = .5) +
  labs(title = "Résidus vs levier", subtitle = "Points influents : levier ET résidu forts (grande distance de Cook)")

(g1 + g2) / (g3 + g4)

## 7.2 Le tout-en-un d'easystats ----
check_model(m_diag)                # panneau complet avec aide à la lecture

## 7.3 Vérifications ciblées (tests + graphiques) ----
# ATTENTION : les tests de présupposés sont sensibles à n. Pour n grand, tout
# écart minime est "significatif" ; pour n petit, ils manquent de puissance.
# => privilégier la lecture des graphiques, utiliser les tests comme aide.

# Normalité
check_normality(m_diag)            # Shapiro-Wilk sur résidus
plot(check_normality(m_diag), type = "qq")

# Homoscédasticité
check_heteroscedasticity(m_diag)   # test de Breusch-Pagan
plot(check_heteroscedasticity(m_diag))

# Indépendance : Durbin-Watson (utile pour séries temporelles/ordre des données)
check_autocorrelation(m_diag)      # ici l'ordre des lignes n'a pas de sens biologique
# Pour de l'indépendance spatiale/temporelle/groupée : voir variogramme, ACF,
# ou mieux, raisonner sur le plan d'échantillonnage => modèles mixtes (section 10).

# Points aberrants / influents
check_outliers(m_diag, method = "cook")
cooks <- cooks.distance(m_diag)
which(cooks > 4 / nrow(peng))      # règle empirique : D > 4/n
plot(check_outliers(m_diag))
# Que faire d'un point influent ? Vérifier une erreur de saisie ; refaire le modèle
# avec/sans et COMPARER ; ne jamais supprimer sans justification biologique.

# Colinéarité : voir section 6.5
check_collinearity(m_diag)

# Vérification prédictive : le modèle génère-t-il des données ressemblant aux observées ?
check_predictions(m_diag)

## 7.4 Voir ce que "violer" un présupposé veut dire (simulations) ----
# On connaît le vrai modèle, donc on SAIT ce qui cloche.
set.seed(2026)
n_sim <- 200
sim <- tibble(x = runif(n_sim, 1, 10)) |>
  mutate(y_ok     = 2 + 3 * x + rnorm(n_sim, 0, 2),               # tout va bien
         y_nonlin = 2 + 0.8 * x^2 + rnorm(n_sim, 0, 2),           # non-linéarité
         y_hetero = 2 + 3 * x + rnorm(n_sim, 0, 0.6 * x),         # variance croissante
         y_skew   = 2 + 3 * x + (rexp(n_sim, rate = .5) - 2))     # erreurs asymétriques

check_model(lm(y_ok ~ x, data = sim))        # référence
check_model(lm(y_nonlin ~ x, data = sim))    # courbe en U sur résidus vs ajustées
check_model(lm(y_hetero ~ x, data = sim))    # "entonnoir"
check_model(lm(y_skew ~ x, data = sim))      # QQ-plot courbé

## 7.5 Remèdes ----
# (a) Non-linéarité : changer la forme fonctionnelle
m_poly <- lm(y_nonlin ~ poly(x, 2), data = sim)     # ou  y ~ x + I(x^2)
check_model(m_poly)
# (alternative flexible : GAM, mgcv::gam(y ~ s(x)))

# (b) Hétéroscédasticité :
#   - transformer y (log, racine) ; ou modéliser la variance avec un GLM adapté
#   - erreurs-types robustes (HC3) : corrige les SE sans changer les coefficients
m_het <- lm(y_hetero ~ x, data = sim)
model_parameters(m_het)                             # SE classiques
model_parameters(m_het, vcov = "HC3")               # SE robustes
lmtest::coeftest(m_het, vcov = sandwich::vcovHC(m_het, type = "HC3"))  # idem
#   - moindres carrés pondérés quand la structure de variance est connue :
m_wls <- lm(y_hetero ~ x, data = sim, weights = 1 / x^2)
model_parameters(m_wls)

# (c) Non-normalité :
#   - Si n est grand et seulement l'inférence sur les coefficients importe, le TCL
#     protège (mais les intervalles de PRÉDICTION restent sensibles à la forme).
#   - Bootstrap des coefficients :
model_parameters(lm(y_skew ~ x, data = sim), bootstrap = TRUE, iterations = 500)
#   - Transformation de Box-Cox (λ proche de 0 = log, 0.5 = racine, 1 = rien)
MASS::boxcox(body_mass_g ~ flipper_length_mm + species, data = peng)
#   - Ou un GLM avec la bonne distribution (section 9) : souvent la meilleure option.

# (d) Non-indépendance : modèles mixtes (section 10), pas de rustine possible.


# 8. APERÇU TIDYMODELS : AJUSTER ET ÉVALUER LA PRÉDICTION -------------------
# Autre philosophie : l'objectif est de PRÉDIRE hors échantillon (on évalue sur des
# données non vues) plutôt que d'interpréter des coefficients et des p-values.
# Les deux approches sont complémentaires.

set.seed(123)
split <- initial_split(peng, prop = 0.8, strata = species)
train <- training(split)
test  <- testing(split)

# Recipe : prétraitement reproductible (appris sur train, appliqué à test)
rec <- recipe(body_mass_g ~ flipper_length_mm + bill_length_mm + species, data = train) |>
  step_dummy(all_nominal_predictors()) |>
  step_normalize(all_numeric_predictors())

lm_spec <- linear_reg() |> set_engine("lm")     # moteur interchangeable (glmnet, stan...)

wf <- workflow() |> add_recipe(rec) |> add_model(lm_spec)
fit_wf <- fit(wf, data = train)
tidy(fit_wf)                                    # coefficients (variables normalisées)

# Évaluation sur l'échantillon test
preds <- predict(fit_wf, new_data = test) |> bind_cols(test |> select(body_mass_g))
metric_set(rmse, rsq, mae)(preds, truth = body_mass_g, estimate = .pred)

# Validation croisée 10-fold sur le train : estimation plus stable de l'erreur
folds <- vfold_cv(train, v = 10, strata = species)
cv_res <- fit_resamples(wf, resamples = folds, metrics = metric_set(rmse, rsq))
collect_metrics(cv_res)

# Récupérer le lm sous-jacent pour utiliser easystats dessus :
lm_engine <- extract_fit_engine(fit_wf)
model_parameters(lm_engine)


# 9. MODÈLES LINÉAIRES GÉNÉRALISÉS (GLM) ----------------------------------

## 9.1 Principe ----
# Le LM impose : y ~ Normale, variance constante, y non borné. Souvent faux :
#   - y binaire (sexe, survie, présence)            -> Binomiale/Bernoulli
#   - y comptage (nb de parasites, d'œufs)          -> Poisson / Binomiale négative
#   - y positif asymétrique (durées, masses)        -> Gamma / log-normale
#   - y proportion continue dans ]0, 1[             -> Beta
#
# Un GLM a TROIS composantes :
#   1. Composante aléatoire : y_i ~ loi de la famille exponentielle (moyenne μ_i)
#   2. Prédicteur linéaire :  η_i = β0 + β1 x1i + ... 
#   3. Fonction de lien g :   g(μ_i) = η_i   (μ_i = g⁻¹(η_i) reste dans le bon domaine)
# Le LM est le cas particulier : loi Normale + lien identité.
#
# Estimation : maximum de vraisemblance (algorithme IRLS), et non plus MCO.
# Mesure d'écart du modèle : la DÉVIANCE (généralise la somme des carrés),
# D = 2 [logLik(modèle saturé) - logLik(modèle)].
# Comparaison de modèles emboîtés : test du rapport de vraisemblance (LRT, χ²).

# Le LM est un GLM gaussien :
m_glm_gauss <- glm(body_mass_g ~ flipper_length_mm, family = gaussian(link = "identity"),
                   data = peng)
all.equal(coef(m_glm_gauss), coef(m1))

## 9.2 APARTÉ : LES DISTRIBUTIONS ET LES LIENS ----

### 9.2.1 Tableau récapitulatif (à garder sous la main) ----
# | Loi            | Support          | Moyenne | Variance            | Lien canonique | Usage typique            |
# |----------------|------------------|---------|---------------------|----------------|--------------------------|
# | Normale        | ℝ                | μ       | σ² (constante)      | identité       | mesures continues        |
# | Bernoulli      | {0,1}            | p       | p(1-p)              | logit          | présence/absence         |
# | Binomiale      | {0,…,N}          | Np      | Np(1-p)             | logit          | nb succès sur N essais   |
# | Poisson        | {0,1,2,…}        | λ       | λ  (= moyenne)      | log            | comptages (sans surdisp.)|
# | Binom. négative| {0,1,2,…}        | μ       | μ + μ²/θ            | log            | comptages sur-dispersés  |
# | Gamma          | ]0, +∞[          | μ       | μ²/shape  (φ·μ²)    | inverse (log)  | positifs asymétriques    |
# | Log-normale    | ]0, +∞[          | e^(μ+σ²/2) | ...              | (lm sur log y) | positifs, effets multipl.|
# | Beta           | ]0, 1[           | μ       | μ(1-μ)/(1+φ)        | logit          | proportions continues    |
#
# À retenir : ce sont les RELATIONS MOYENNE-VARIANCE qui distinguent les lois.
#   Poisson  : var = moyenne ;  NB : var > moyenne ;  Gamma : var ∝ moyenne² ;
#   Binomiale/Beta : variance maximale à p = 0,5 et nulle aux bornes.
# Autres : quasi-Poisson/quasi-binomiale (dispersion libre), zero-inflated/hurdle
# (excès de zéros), Tweedie (positif avec zéros exacts), ordinale/multinomiale
# (ordered/nnet::multinom ; ordinal::clm).

### 9.2.2 Voir les distributions ----
x_cont <- tibble(x = seq(-4, 4, length.out = 300))
d_cont <- bind_rows(
  tibble(x = seq(-4, 12, length.out = 400)) |>
    mutate(densite = dnorm(x, 4, 1.5), loi = "Normale (μ=4, σ=1,5)"),
  tibble(x = seq(0.01, 15, length.out = 400)) |>
    mutate(densite = dgamma(x, shape = 2, rate = .5), loi = "Gamma (forme=2, taux=0,5)"),
  tibble(x = seq(0.01, 15, length.out = 400)) |>
    mutate(densite = dgamma(x, shape = 8, rate = 2), loi = "Gamma (forme=8, taux=2)"),
  tibble(x = seq(0.001, .999, length.out = 400)) |>
    mutate(densite = dbeta(x, 2, 5), loi = "Beta (2, 5)"),
  tibble(x = seq(0.001, .999, length.out = 400)) |>
    mutate(densite = dbeta(x, .5, .5), loi = "Beta (0,5 ; 0,5)")
)
g_cont <- ggplot(d_cont, aes(x, densite)) + geom_line(linewidth = 1, color = "#2c7fb8") +
  facet_wrap(~loi, scales = "free") + labs(title = "Lois continues")

d_disc <- bind_rows(
  tibble(k = 0:10)  |> mutate(prob = dbinom(k, 10, .3),  loi = "Binomiale (N=10, p=0,3)"),
  tibble(k = 0:20)  |> mutate(prob = dpois(k, 3),        loi = "Poisson (λ=3)"),
  tibble(k = 0:20)  |> mutate(prob = dpois(k, 10),       loi = "Poisson (λ=10)"),
  tibble(k = 0:20)  |> mutate(prob = dnbinom(k, mu = 3, size = 1), loi = "Bin. négative (μ=3, θ=1)"),
  tibble(k = 0:20)  |> mutate(prob = dnbinom(k, mu = 10, size = 2), loi = "Bin. négative (μ=10, θ=2)")
)
g_disc <- ggplot(d_disc, aes(k, prob)) + geom_col(fill = "#d95f0e") +
  facet_wrap(~loi, scales = "free") + labs(title = "Lois discrètes")
g_cont / g_disc

### 9.2.3 Relations moyenne-variance, par simulation ----
set.seed(1); N <- 1e5
sims <- list(Normale   = rnorm(N, 10, 2),
             Binomiale = rbinom(N, 20, .5),
             Poisson   = rpois(N, 10),
             `Bin. négative` = rnbinom(N, mu = 10, size = 2),
             Gamma     = rgamma(N, shape = 5, rate = .5))
map_dfr(sims, ~ tibble(moyenne = mean(.x), variance = var(.x)), .id = "loi") |>
  mutate(variance_sur_moyenne = variance / moyenne)
# Poisson : ratio ≈ 1 ; NB : ≈ 1 + μ/θ = 6 ; Gamma : ≈ μ/shape = 2 ; Normale : ratio sans signification.

### 9.2.4 Les fonctions de lien ----
# Elles transforment le prédicteur linéaire η (dans ℝ) vers l'échelle de la moyenne μ.
eta <- seq(-6, 6, length.out = 300)
tibble(eta,
       logit   = plogis(eta),                 # g⁻¹ = 1/(1+exp(-η))
       probit  = pnorm(eta),                  # fonction de répartition normale
       cloglog = 1 - exp(-exp(eta))) |>       # asymétrique
  pivot_longer(-eta, names_to = "lien", values_to = "probabilite") |>
  ggplot(aes(eta, probabilite, color = lien)) + geom_line(linewidth = 1) +
  labs(title = "Liens pour données binaires : η (échelle du lien) -> p",
       x = "prédicteur linéaire η")
# Toutes les familles sont des objets R : lien direct et inverse accessibles
binomial()$linkfun(0.7)        # logit(0.7) = log(0.7/0.3)
binomial()$linkinv(0.847)      # retour à la probabilité
poisson()$linkfun(5)           # log(5)
poisson()$linkinv(2)           # exp(2)
Gamma(link = "log")$linkinv(1)
# INTERPRÉTATION selon le lien :
#   identité : effet ADDITIF sur la moyenne (LM)
#   log      : effet MULTIPLICATIF : exp(β) = rapport de moyennes (ou de taux)
#   logit    : effet additif sur le log-odds : exp(β) = odds ratio (OR)

## 9.3 RÉGRESSION LOGISTIQUE (réponse binaire) ----
# Question : peut-on prédire le sexe d'un pingouin à partir de sa masse ?
# Modèle : sexe_i ~ Bernoulli(p_i) ;  logit(p_i) = log(p_i / (1-p_i)) = β0 + β1 x_i
# Ici R modélise P(sexe = "male") car "male" est le 2e niveau du facteur.
levels(peng$sex)

### 9.3.1 Intuition : odds et log-odds ----
peng |>
  count(species, sex) |>
  pivot_wider(names_from = sex, values_from = n) |>
  mutate(p_male = male / (male + female),
         odds   = male / female,            # côte
         log_odds = log(odds))
# odds = p/(1-p) ; log-odds ∈ ℝ => modélisable linéairement.

### 9.3.2 Ajustement ----
# mass_c = (masse en 100 g) centrée : l'intercept correspond à la masse moyenne,
# et la pente à +100 g.
m_logit0 <- glm(sex ~ mass_c, data = peng, family = binomial(link = "logit"))
summary(m_logit0)
model_parameters(m_logit0)                      # échelle log-odds
model_parameters(m_logit0, exponentiate = TRUE) # odds ratios
# Lecture : exp(β1) = facteur multiplicatif des ODDS de "male" pour +100 g.
# ATTENTION : odds ratio ≠ risque relatif, sauf si l'événement est rare.

### 9.3.3 Estimation par maximum de vraisemblance "à la main" ----
# logL(β) = Σ [ y_i log(p_i) + (1-y_i) log(1-p_i) ], avec p_i = plogis(β0 + β1 x_i)
y01 <- as.numeric(peng$sex == "male")
nll <- function(par, x, y) {                    # -log vraisemblance
  p <- plogis(par[1] + par[2] * x)
  -sum(dbinom(y, size = 1, prob = p, log = TRUE))
}
opt <- optim(c(0, 0), nll, x = peng$mass_c, y = y01, method = "BFGS", hessian = TRUE)
rbind(optim = opt$par, glm = coef(m_logit0))
# Erreurs-types = racine de la diagonale de l'inverse de la hessienne (info de Fisher)
rbind(optim = sqrt(diag(solve(opt$hessian))), glm = summary(m_logit0)$coefficients[, 2])
c(optim = -opt$value, glm = as.numeric(logLik(m_logit0)))

### 9.3.4 Déviance, tests, intervalles ----
m_logit0$null.deviance        # déviance du modèle nul (intercept seul)
m_logit0$deviance             # déviance résiduelle
# LRT : le modèle avec masse est-il meilleur que le nul ?
anova(m_logit0, test = "Chisq")
# Intervalles : Wald (rapide, approximatif) vs profil de vraisemblance (plus fiable)
confint.default(m_logit0)     # Wald
confint(m_logit0)             # profil

### 9.3.5 Ajouter l'espèce (la masse seule est confondue avec l'espèce) ----
m_logit1 <- glm(sex ~ mass_c + species, data = peng, family = binomial)
model_parameters(m_logit1, exponentiate = TRUE)
anova(m_logit0, m_logit1, test = "Chisq")
compare_performance(m_logit0, m_logit1)
# Variantes de lien : probit / cloglog (même structure, interprétation du lien différente)
m_probit <- glm(sex ~ mass_c + species, data = peng, family = binomial("probit"))
compare_performance(m_logit1, m_probit)

### 9.3.6 Représenter le modèle sur l'échelle de la probabilité ----
estimate_relation(m_logit1, by = c("mass_c", "species"), length = 40) |> plot()
# Courbe en S. Les effets sont additifs sur le logit mais PAS sur la probabilité :
# l'écart en probabilité dépend du point de départ (non-linéarité du lien).
estimate_means(m_logit1, by = "species")               # probabilités marginales

### 9.3.7 Qualité de l'ajustement et pouvoir prédictif ----
r2_tjur(m_logit1)                                      # pseudo-R² de Tjur (le R² des GLM n'existe pas !)
model_performance(m_logit1)                            # AIC, BIC, Tjur, RMSE, log-loss, score...
plot(binned_residuals(m_logit1))                       # résidus binés : ~ 95 % dans les bandes
check_model(m_logit1)

prob <- predict(m_logit1, type = "response")           # type="response" => probabilités
res_cls <- tibble(truth = peng$sex,
                  prob  = prob,
                  pred  = factor(if_else(prob > .5, "male", "female"), levels = levels(peng$sex)))
conf_mat(res_cls, truth, pred)                         # matrice de confusion
accuracy(res_cls, truth, pred)
roc_auc(res_cls, truth, prob, event_level = "second")  # AUC (le 2e niveau = "male")
roc_curve(res_cls, truth, prob, event_level = "second") |> autoplot()
# Ici l'évaluation est faite sur les données d'ajustement (optimiste) : pour une vraie
# évaluation prédictive, utiliser train/test ou validation croisée (cf. section 8).

### 9.3.8 Binomiale avec effectifs agrégés ----
# Quand on dispose du nombre de succès et d'échecs par groupe : cbind(succès, échecs)
agg <- peng |>
  count(species, island, year, sex) |>
  pivot_wider(names_from = sex, values_from = n, values_fill = 0)
head(agg)
m_agg <- glm(cbind(male, female) ~ species, data = agg, family = binomial)
model_parameters(m_agg, exponentiate = TRUE)
# Équivalent au modèle sur données individuelles ; l'intercept teste (en log-odds) un
# sex-ratio différent de 1:1. Penser à vérifier la sur-dispersion (check_overdispersion)
# quand les effectifs par groupe sont grands.

## 9.4 RÉGRESSION DE POISSON ET COMPTAGES ----
# Le jeu palmerpenguins ne contient pas de vrais comptages. On SIMULE donc un nombre de
# parasites par individu, avec une vérité connue :
#   log(λ_i) = log(jours_i) + 0.5 + 0.4 × masse_z + effet_espèce
#   (Chinstrap : +0.3 ; Gentoo : -0.4 ; Adelie : référence)
# "jours" = durée d'observation différente par individu => OFFSET.
set.seed(42)
peng_sim <- peng |>
  mutate(mass_z = as.numeric(scale(body_mass_g)),
         jours  = sample(5:15, n(), replace = TRUE),
         lambda = jours * exp(0.5 + 0.4 * mass_z +
                              case_when(species == "Chinstrap" ~ 0.3,
                                        species == "Gentoo"    ~ -0.4,
                                        .default = 0)),
         parasites_pois = rpois(n(), lambda),                       # Poisson pur
         parasites_nb   = rnbinom(n(), mu = lambda, size = 1.2))   # sur-dispersé

ggplot(peng_sim, aes(parasites_pois)) + geom_histogram(bins = 30) +
  labs(title = "Comptages (Poisson simulé)")

### 9.4.1 Ajustement et lecture ----
# Offset : on modélise un TAUX (parasites par jour) ; le coefficient de l'offset est fixé à 1.
m_pois <- glm(parasites_pois ~ mass_z + species + offset(log(jours)),
              data = peng_sim, family = poisson(link = "log"))
model_parameters(m_pois)                        # échelle log
model_parameters(m_pois, exponentiate = TRUE)   # IRR = rapports de taux
# Les estimations retrouvent-elles la vérité ?
tibble(parametre = names(coef(m_pois)),
       vrai = c(0.5, 0.4, 0.3, -0.4),
       estime = coef(m_pois))
# exp(β_mass_z) ≈ exp(0.4) = 1,49 : +1 écart-type de masse multiplie le taux par ~1,5.

### 9.4.2 Hypothèse clé : Var = moyenne. Sur-dispersion ? ----
# Sur-dispersion = variance supérieure à celle prévue => SE sous-estimés, p-values trop optimistes.
# Dispersion de Pearson à la main : doit être ≈ 1.
sum(residuals(m_pois, type = "pearson")^2) / df.residual(m_pois)
check_overdispersion(m_pois)                    # OK sur données Poisson

# Même modèle sur des données sur-dispersées :
m_pois_bad <- glm(parasites_nb ~ mass_z + species + offset(log(jours)),
                  data = peng_sim, family = poisson)
check_overdispersion(m_pois_bad)                # sur-dispersion détectée
plot(check_overdispersion(m_pois_bad))

### 9.4.3 Remèdes ----
# (1) quasi-Poisson : même coefficients, SE corrigés par la dispersion estimée.
#     Pas de vraisemblance => pas d'AIC.
m_quasi <- glm(parasites_nb ~ mass_z + species + offset(log(jours)),
               data = peng_sim, family = quasipoisson)
compare_parameters(m_pois_bad, m_quasi)         # coefficients identiques, SE différents

# (2) Binomiale négative : var = μ + μ²/θ. Vraie vraisemblance : AIC, LRT possibles.
m_nb <- MASS::glm.nb(parasites_nb ~ mass_z + species + offset(log(jours)), data = peng_sim)
# Équivalent glmmTMB (famille nbinom2, la plus courante) :
m_nb_tmb <- glmmTMB(parasites_nb ~ mass_z + species + offset(log(jours)),
                    data = peng_sim, family = nbinom2)
check_overdispersion(m_nb_tmb)
compare_performance(m_pois_bad, m_nb, m_nb_tmb)
# (θ vrai = 1,2 : retrouvé ?)  
m_nb$theta

# Contrôle par résidus simulés (DHARMa), optionnel mais très recommandé :
if (requireNamespace("DHARMa", quietly = TRUE)) {
  plot(DHARMa::simulateResiduals(m_nb_tmb))
}

### 9.4.4 Excès de zéros ----
# Si trop de zéros pour la loi choisie : modèle à inflation de zéros (ZIP/ZINB)
# ou "hurdle". ZI : un processus produit des zéros "structurels" (ex. individu non
# exposé) en plus des zéros d'échantillonnage.
set.seed(7)
peng_sim <- peng_sim |>
  mutate(parasites_zi = if_else(rbinom(n(), 1, .3) == 1, 0L, rpois(n(), lambda)))
m_zi_bad <- glm(parasites_zi ~ mass_z + species + offset(log(jours)),
                data = peng_sim, family = poisson)
check_zeroinflation(m_zi_bad)                   # trop de zéros observés vs attendus
m_zip <- glmmTMB(parasites_zi ~ mass_z + species + offset(log(jours)),
                 ziformula = ~ 1,               # probabilité de zéro structurel (constante)
                 data = peng_sim, family = poisson)
model_parameters(m_zip)                         # blocs "conditional" (comptage) et "zero_inflated"
plogis(fixef(m_zip)$zi)                         # proportion de zéros structurels (vrai : 0,30)

## 9.5 GLM GAMMA (positif continu asymétrique) ----
# La masse est strictement positive ; si la variabilité croît avec la moyenne
# (écart-type ∝ moyenne = coefficient de variation constant), Gamma + lien log convient.
m_lm    <- lm(body_mass_g ~ flipper_length_mm + species, data = peng)
m_loglm <- lm(log(body_mass_g) ~ flipper_length_mm + species, data = peng)   # alternative classique
m_gam   <- glm(body_mass_g ~ flipper_length_mm + species,
               data = peng, family = Gamma(link = "log"))
model_parameters(m_gam, exponentiate = TRUE)   # rapports multiplicatifs de masse moyenne
# +1 mm de nageoire <=> masse moyenne × exp(β) (≈ +1 %). Gentoo vs Adelie : × exp(β).
summary(m_gam)$dispersion                       # φ = 1/shape ; CV² ≈ φ
# Différence avec lm(log y) : ici on modélise log(E[y]) et non E[log(y)] (plus direct
# à interpréter, pas de retour-transformation biaisé).
compare_performance(m_lm, m_gam)                # AIC comparable entre ces deux (même y) ; pas avec m_loglm
check_model(m_gam)

## 9.6 GLM BETA (proportions continues dans ]0,1[) ----
# Exemple : rapport épaisseur/longueur du bec, naturellement compris entre 0 et 1.
peng <- peng |> mutate(ratio_bec = bill_depth_mm / bill_length_mm)
range(peng$ratio_bec)
m_beta <- glmmTMB(ratio_bec ~ species + sex, data = peng, family = beta_family(link = "logit"))
model_parameters(m_beta)
# Interprétation sur l'échelle logit, comme pour la logistique (mais y n'est PAS binaire).
# Si présence de 0 et/ou 1 exacts : modèle zero-one-inflated beta (glmmTMB/brms).
estimate_means(m_beta, by = "species")
check_model(m_beta)

## 9.7 Présupposés d'un GLM : ce qui change ----
# - Plus d'hypothèse de normalité des résidus (ni d'homoscédasticité) ;
#   MAIS la relation moyenne-variance de la famille doit être correcte.
# - Linéarité sur l'échelle du LIEN (η = Xβ).
# - Indépendance des observations (toujours !).
# - Pas de sur-dispersion non modélisée (Poisson/binomiale), pas d'excès de zéros.
# - Pas de séparation complète (logistique : un prédicteur sépare parfaitement 0/1
#   => coefficients qui "explosent", message "fitted probabilities 0 or 1").
# - Taille d'échantillon suffisante (≥ 10 événements par paramètre, règle empirique).
# - Les résidus bruts sont peu lisibles : utiliser résidus de déviance, binés, ou
#   résidus quantiles simulés (DHARMa).


# 10. MODÈLES MIXTES -------------------------------------------------------

## 10.1 Pourquoi ? ----
# Observations NON indépendantes : individus mesurés plusieurs fois, animaux d'un
# même nid / élevage / troupeau / parcelle, sites groupés...
# Ignorer ce regroupement => pseudo-réplication : SE trop petits, faux positifs.
#
# Effets FIXES : ce qu'on veut estimer (effet du sexe, de l'espèce, du traitement) ;
#                niveaux d'intérêt en soi, reproductibles.
# Effets ALÉATOIRES : variabilité entre unités d'un échantillon d'une population
#                de groupes (nids, individus, blocs...). On estime leur VARIANCE,
#                pas un coefficient par groupe : "partial pooling" (retrait vers la moyenne).
#                Nécessite un nombre suffisant de niveaux (≥ 5-6, idéalement 10+).
#
# Syntaxe lme4 (formule) :
#   (1 | g)            intercept aléatoire par groupe g
#   (1 + x | g)        intercept ET pente de x aléatoires, corrélés
#   (1 + x || g)       idem mais non corrélés
#   (1 | a) + (1 | b)  deux effets aléatoires croisés
#   (1 | a / b)        b emboîté dans a  (= (1|a) + (1|a:b))
# Modèle : y_ij = β0 + β1 x_ij + u_j + ε_ij ,  u_j ~ N(0, σ²_u), ε_ij ~ N(0, σ²)

## 10.2 Données réelles avec regroupement : les nids ----
# Dans penguins_raw, `Individual ID` (ex. "N12A1") code le NID (N12) et l'adulte (A1/A2).
# Les deux partenaires d'un nid partagent le même environnement et ne sont pas
# indépendants. Chaque nid est mesuré (au plus) 2 fois : cas minimal mais réel.
nests <- palmerpenguins::penguins_raw |>
  transmute(study = studyName,
            species = factor(word(Species, 1)),
            island  = Island,
            ind_id  = `Individual ID`,
            nest    = str_extract(`Individual ID`, "^N\\d+"),
            sex     = factor(str_to_lower(Sex)),
            body_mass_g       = `Body Mass (g)`,
            flipper_length_mm = `Flipper Length (mm)`,
            year    = year(`Date Egg`)) |>
  drop_na(nest, sex, body_mass_g, flipper_length_mm) |>
  # Un nid = étude + espèce + île + n° de nid (les numéros sont réutilisés d'une étude à l'autre)
  mutate(nest_id = factor(paste(study, species, island, nest, sep = "_")))

nests |> count(nest_id) |> count(n, name = "nb_nids")   # taille des groupes
nlevels(nests$nest_id)

## 10.3 Modèle linéaire mixte (LMM) à intercept aléatoire ----
m_naif  <- lm(body_mass_g ~ sex + species, data = nests)                 # ignore les nids
m_mixte <- lmer(body_mass_g ~ sex + species + (1 | nest_id), data = nests)  # REML par défaut

summary(m_mixte)
# Lire : "Random effects" (variance entre nids, variance résiduelle) puis "Fixed effects".
VarCorr(m_mixte)
model_parameters(m_mixte)                       # effets fixes (df de Satterthwaite) + effets aléatoires

# Comparaison avec le modèle naïf : coefficients proches, mais les SE changent.
# Le sexe varie À L'INTÉRIEUR des nids (1 mâle + 1 femelle) => comparaison intra-nid,
# plus précise : le SE de l'effet sexe DIMINUE ; à l'inverse, un prédicteur variant
# ENTRE nids (ex. l'espèce) verrait son SE augmenter ou rester proche.
compare_parameters(m_naif, m_mixte)

### 10.3.1 Variance, ICC, R² ----
# ICC = σ²_u / (σ²_u + σ²) : part de variance due aux groupes (ressemblance intra-groupe).
icc(m_mixte)
# R² marginal : variance expliquée par les effets fixes ; conditionnel : fixes + aléatoires.
r2_nakagawa(m_mixte)

### 10.3.2 Effets aléatoires prédits (BLUPs) et retrait (shrinkage) ----
blup <- as_tibble(ranef(m_mixte)$nest_id, rownames = "nest_id") |>
  rename(blup = `(Intercept)`)
brut <- nests |>
  mutate(res_fixe = body_mass_g - predict(m_mixte, re.form = NA)) |>   # résidu vs effets fixes
  group_by(nest_id) |>
  summarise(moyenne_brute = mean(res_fixe), n_obs = n())
left_join(brut, blup, by = "nest_id") |>
  ggplot(aes(moyenne_brute, blup)) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Retrait vers 0 des effets aléatoires (shrinkage)",
       subtitle = "Les nids à valeurs extrêmes sont ramenés vers la moyenne ; plus n_obs est petit, plus le retrait est fort",
       x = "Moyenne brute des résidus du nid", y = "BLUP (prédiction de l'effet aléatoire)")

### 10.3.3 Tests ----
# Effets fixes : tests de Satterthwaite (lmerTest)
anova(m_mixte)
# Effet aléatoire : LRT (valeur-p à interpréter avec prudence : test à la frontière)
ranova(m_mixte)
# Comparer des structures de EFFETS FIXES => ajuster en ML (REML=FALSE).
# Comparer des structures d'EFFETS ALÉATOIRES (mêmes fixes) => REML.
m_ml0 <- update(m_mixte, . ~ . - sex, REML = FALSE)
m_ml1 <- update(m_mixte, REML = FALSE)
anova(m_ml0, m_ml1)

### 10.3.4 Diagnostics ----
# Présupposés : résidus ~ N(0, σ²) ; effets aléatoires ~ N(0, σ²_u) ; indépendance
# CONDITIONNELLE aux effets aléatoires ; linéarité ; homoscédasticité.
check_model(m_mixte)                            # inclut le QQ-plot des effets aléatoires
check_singularity(m_mixte)                      # variance aléatoire estimée à 0 / corrélation ±1 ?
# Piège : (1|island) avec 3 îles seulement : trop peu de niveaux pour estimer une variance ;
# mieux vaut mettre island en effet FIXE.

# Rapport automatique
report(m_mixte)

## 10.4 Pentes aléatoires : données longitudinales (simulation) ----
# Suivi hebdomadaire de 40 individus (8 pesées) pendant la couvaison (jeûne) : la masse
# diminue, mais le rythme de perte VARIE selon l'individu => intercept ET pente aléatoires.
# On simule à partir de vraies masses de départ.
set.seed(7)
n_ind <- 40
ind <- peng |>
  slice_sample(n = n_ind) |>
  transmute(id = factor(row_number()), species, mass0 = body_mass_g) |>
  mutate(pente = -35 - 0.03 * (mass0 - mean(mass0)) + rnorm(n_ind, 0, 10))   # g/semaine

suivi <- ind |>
  crossing(semaine = 0:7) |>
  mutate(masse = mass0 + pente * semaine + rnorm(n(), 0, 60))

ggplot(suivi, aes(semaine, masse, group = id, color = species)) +
  geom_line(alpha = .6) + scale_color_manual(values = pal) +
  labs(title = "Évolution individuelle de la masse (simulation)")

# Trois structures aléatoires, de la plus simple à la plus riche :
m_ri   <- lmer(masse ~ semaine + species + (1 | id),             data = suivi)   # intercept seul
m_rs0  <- lmer(masse ~ semaine + species + (1 + semaine || id),  data = suivi)   # pentes, non corrélées
m_rs   <- lmer(masse ~ semaine + species + (1 + semaine | id),   data = suivi)   # pentes, corrélées
VarCorr(m_rs)                    # écarts-types intercept, pente, corrélation intercept-pente
# Vrai : sd pente = 10 ; corr négative (les gros perdent plus vite).

anova(m_ri, m_rs0, m_rs)         # LRT + AIC/BIC
compare_performance(m_ri, m_rs0, m_rs, rank = TRUE)

model_parameters(m_rs)           # effet moyen de la semaine : ≈ -35 g/sem (vrai)

# Visualiser : population (effets fixes) vs individus (effets fixes + aléatoires)
suivi <- suivi |>
  mutate(pred_pop = predict(m_rs, re.form = NA),
         pred_ind = predict(m_rs))
ggplot(suivi, aes(semaine)) +
  geom_point(aes(y = masse), alpha = .25) +
  geom_line(aes(y = pred_ind, group = id), alpha = .4) +
  geom_line(aes(y = pred_pop, group = species), color = "red", linewidth = 1) +
  facet_wrap(~species) +
  labs(title = "Prédictions individuelles (gris) et de population (rouge)")

# Bonnes pratiques :
# - Le modèle "maximal" justifié par le plan (Barr et al.) est un point de départ, mais
#   avec peu de données il devient singulier : simplifier (|| puis supprimer la pente).
# - Centrer le temps réduit la corrélation intercept-pente (intercept = masse à mi-suivi).
# - Si les résidus intra-individu sont autocorrélés (mesures rapprochées),
#   nlme::lme(..., correlation = corAR1(form = ~ semaine | id)) ou glmmTMB (ar1()).
check_model(m_rs)

## 10.5 GLMM : modèles mixtes généralisés ----
# = GLM + effets aléatoires. Même logique, y non gaussien (comptage, binaire...).
# On reprend les parasites, avec cette fois un EFFET NID simulé : certains nids sont
# plus infestés que d'autres, quelles que soient la masse ou l'espèce.
set.seed(99)
effet_nid <- tibble(nest_id = levels(nests$nest_id),
                    u = rnorm(nlevels(nests$nest_id), 0, 0.5))       # sd vraie = 0.5
nests_sim <- nests |>
  left_join(effet_nid, by = "nest_id") |>
  mutate(mass_z = as.numeric(scale(body_mass_g)),
         parasites = rpois(n(), exp(1 + 0.3 * mass_z + u)))         # vérité : β0=1, β_masse=0.3

# (a) GLM Poisson qui ignore les nids
glm_p <- glm(parasites ~ mass_z + species, data = nests_sim, family = poisson)
check_overdispersion(glm_p)       # variance excessive : due à l'effet nid non modélisé

# (b) GLMM Poisson avec intercept aléatoire de nid (approximation de Laplace)
glmm_p <- glmer(parasites ~ mass_z + species + (1 | nest_id),
                data = nests_sim, family = poisson)
# Équivalent glmmTMB :
glmm_p_tmb <- glmmTMB(parasites ~ mass_z + species + (1 | nest_id),
                      data = nests_sim, family = poisson)
check_overdispersion(glmm_p_tmb)  # corrigée
compare_parameters(glm_p, glmm_p)
VarCorr(glmm_p)                   # sd du nid ≈ 0.5 ?
compare_performance(glm_p, glmm_p_tmb)
icc(glmm_p)                       # ICC sur l'échelle latente (log)

# Interprétation d'un GLMM :
#  - Les coefficients sont conditionnels aux effets aléatoires (individu "typique", u=0) :
#    ils diffèrent des effets marginaux (moyenne de population) dès que le lien est non
#    linéaire (log, logit) — particularité propre aux GLMM, pas aux LMM.
#  - Pour un résultat binaire : glmer(y ~ x + (1 | g), family = binomial).
#  - Convergence difficile ? Essayer glmmTMB, ou lme4 avec control = glmerControl(optimizer = "bobyqa"),
#    standardiser les variables continues, simplifier la structure aléatoire.
#  - Résidus : DHARMa::simulateResiduals() est l'outil de référence pour les GLMM.
estimate_means(glmm_p_tmb, by = "species")        # moyennes marginales (échelle de la réponse)


# 11. CHECK-LIST FINALE POUR LES ÉTUDIANTS ---------------------------------
# 1. Question et PLAN D'ÉCHANTILLONNAGE d'abord : unités statistiques, répétitions,
#    regroupements (=> effets aléatoires ?), covariables, hypothèses à tester.
# 2. Explorer : graphiques, NA, valeurs aberrantes, distribution de y, corrélations.
# 3. Choisir la distribution selon la NATURE de y :
#      continu symétrique -> Normale | binaire -> Bernoulli | succès/N -> Binomiale |
#      comptage -> Poisson / NB (offset si exposition) | positif asymétrique -> Gamma /
#      log-normale | proportion ]0,1[ -> Beta | trop de zéros -> ZI / hurdle.
# 4. Données non indépendantes ? -> modèle mixte (≥ 5-6 niveaux pour un effet aléatoire).
# 5. Ajuster, puis VÉRIFIER les présupposés (check_model, check_overdispersion,
#    check_zeroinflation, DHARMa, check_collinearity, check_outliers).
# 6. Comparer des modèles (AIC / LRT / validation croisée) sur la base d'hypothèses,
#    pas de pêche aux p-values. Ne jamais comparer AIC de modèles sur des y différents.
# 7. Interpréter sur la bonne échelle (lien !), avec IC et tailles d'effet, pas
#    seulement des p-values : model_parameters(), estimate_means(), estimate_relation().
# 8. Rapporter : modèle complet (formule, famille, lien), N, structure aléatoire,
#    diagnostics, reproductibilité (set.seed, sessionInfo).

session_info()    # easystats::session_info() : versions des packages utilisés

# RESSOURCES UTILES
#  - easystats : https://easystats.github.io/easystats/
#  - Zuur et al. (2009) "Mixed effects models and extensions in ecology with R"
#  - Bolker et al. (2009) "GLMMs: a practical guide for ecology and evolution" (TREE)
#  - Harrison et al. (2018) "A brief introduction to mixed effects modelling..." (PeerJ)
#  - Barr et al. (2013) "Random effects structure for confirmatory hypothesis testing"
#  - GLMM FAQ de Ben Bolker : https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html
###############################################################################
