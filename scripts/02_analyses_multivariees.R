# Analyses multivariées

# Table des matières ----

# 0. Comment utiliser ce script
# 1. Installer et charger les packages
# 2. Importer les données

# Plan :
#   0. Packages et données
#   1. Exploration : structure, valeurs manquantes, corrélations
#   2. Valeurs atypiques multivariées (distance de Mahalanobis)
#   3. ACP (Analyse en Composantes Principales)
#   4. Distances mixtes et PCoA (Gower)
#   5. Classification non supervisée (CAH et k-means)
#   6. Tests multivariés : MANOVA et PERMANOVA
#   7. Analyse discriminante linéaire (LDA)
#   8. AFDM : variables quantitatives ET qualitatives
#   9. Pièges classiques et check-list


# 0. Comment utiliser ce script ----

# Pour exécuter une commande : se placer sur la ligne de la commande et utiliser 
# la combinaison `Ctrl` + `Entrée`, ou cliquer sur `Run`

# Pour exécuter plusieurs lignes : sélectionner les lignes à exécuter et utiliser
# la combinaison `Ctrl` + `Entrée`, ou cliquer sur `Run`

# Pour exécuter l'ensemble du script : utiliser la combinaison 
# `Ctrl` + `Shift` + `Entrée`, ou cliquer sur `Source`

# Pour afficher les différentes sections : cliquer sur l'icône 📄 à droite du script
# (à côté de `Source`)

# Pour replier les sections : `Alt` + `O`

# Pour développer les sections : `Alt` + `Shift` + `O`

# 1. Installer et charger les packages ----

# Décommenter les lignes ci-dessous si besoin et les exécuter pour 
# télécharger les packages nécessaires

# install.packages("tidyverse")
# install.packages("skimr")
# install.packages("palmerpenguins")
# install.packages("naniar")
# install.packages("broom")
# install.packages("GGally")
# install.packages("FactoMineR")
# install.packages("factoextra")
# install.packages("cluster")
# install.packages("vegan")
# install.packages("naniar")
# install.packages("rstatix")

# ⌨️ Raccourci clavier pour insérer ou retirer un commentaire : 
#     `Ctrl` + `Shift` + `C`

# Charger les packages :

library(tidyverse)          # manipulation et visualisation des données
library(naniar)             # exploration des données manquantes
library(rstatix)            # calculs et tests statistiques
library(GGally)             # ggpairs : matrice de nuages de points
library(correlation)        # calcul et visualisation des corrélations
library(see)                # visualisation de résultats de tests statistiques
library(FactoMineR)         # ACP

# library(skimr)              # exploration des données
# library(palmerpenguins)     # jeu de données
# library(broom)              # mise en forme "tidy" des sorties de modèles
# library(FactoMineR)         # PCA, FAMD avec variables supplémentaires
# library(factoextra)         # visualisations ggplot2 pour les analyses factorielles
# library(cluster)            # daisy (Gower), silhouette
# library(vegan)              # adonis2 (PERMANOVA), betadisper
# MASS est appelé avec MASS:: car MASS::select masque dplyr::select !

# Définir le thème ggplot
theme_set(theme_bw())

# 2. Importer les données ----

# Nous importons les données nettoyées du script `00_exploration_des_donnees.R`
penguins <- readRDS("data/penguins_clean.rds")

# Si besoin, exécuter le code ci-dessous pour nettoyer les données
penguins <- palmerpenguins::penguins |> 
  # transformer la variable `year` en facteur
  mutate(year = factor(year)) |> 
  # supprimer les données manquantes
  drop_na()

# Vérifier les données importées
glimpse(penguins)

# ⌨️ Raccourci clavier pour insérer une flèche d'assignation : # `Alt` + `-`

# 3. Exploration des données ----

# Les méthodes multivariées classiques ne gèrent pas les données 
# manquantes : il faut soit supprimer celles-ci, soit effectuer leur imputation.
# Nous avons déjà supprimé les données manquantes :

miss_var_summary(penguins)

# Nous calculons la moyenne et l'écart-type pour chaque variable quantitative,
# en fonction de l'espèce : 

penguins |> 
  summarise(
    # nombre d'observations
    n = n(),
    # appliquer des calculs à une sélection de colonnes
    across(
      bill_length_mm:body_mass_g,
      # liste des statistiques à calculer
      list(mean = mean, sd = sd),
      # appliquer le format "variable_fonction" aux noms des nouvelles colonnes
      .names = "{.col}_{.fn}"),
    # grouper les calculs par la variable `species`
    .by = species)

# Le package {rstatix} propose des fonctions utiles pour des calculs et 
# tests statistiques, ces fonctions pouvant être utilisées dans des pipes :

penguins |> 
  group_by(species) |> 
  rstatix::get_summary_stats(type = "mean_sd")

# La fonction get_summary_stats() s'applique par défaut aux variables
# quantitatives d'un tableau.

# Les écarts-types et les unités (g vs. mm) diffèrent entre les variables : il faudra
# standardiser les variables (centrer-réduire) pour l'ACP.

# 4. Covariance et variance ----

## 4.1 Calcul de la covariance ----

# La covariance mesure comment deux variables varient l'une par rapport à l'autre.

# cov(X, Y) = Σ (xi - x̄)(yi - ȳ) / (n - 1)
#   > 0 : quand X est au-dessus de sa moyenne, Y l'est en général aussi
#   < 0 : relation inverse
#   ≈ 0 : pas de relation LINÉAIRE

# Exemple : covariance de `flipper_length_mm` et `body_mass_g`
x <- penguins$flipper_length_mm
y <- penguins$body_mass_g  
n <- length(x)

cov_manuel <- sum((x - mean(x)) * (y - mean(y))) / (n - 1)
cov_manuel

# Calcul de la covariance à l'aide de la fonction cov()
cov(x, y)

# Les deux valeurs sont-elles strictement égales ?
all.equal(cov_manuel, cov(x, y))

# La variance est un cas particulier : var(X) = cov(X, X)
var(x)
cov(x, x)

## 4.2 Visualisation ----

# A partir d'un scatter plot, imaginer quatre quadrants, les axes étant 
# définis par les moyennes des deux variables.

# Chaque point du scatter plot contribue par le signe du produit
# (xi - x̄)(yi - ȳ) :

#   - covariance positive quand x et y varient dans la même direction
#   - covariance négative quand x et y varient dans des directions opposées

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

# Visualisation de la covariance de `flipper_length_mm` et `body_mass_g`

penguins |> 
  mutate(flipper_length_mm_diff = flipper_length_mm - mean(flipper_length_mm),
         body_mass_g_diff = body_mass_g - mean(body_mass_g),
         sign = case_when(flipper_length_mm_diff * body_mass_g_diff > 0 ~ "positive",
                          .default = "negative")
  ) |> 
  ggplot(aes(x = flipper_length_mm, y = body_mass_g,
             color = sign)) +
  geom_point(show.legend = FALSE) +
  geom_hline(yintercept = mean(penguins$body_mass_g), linewidth = 0.5, linetype = "dotted") +
  geom_vline(xintercept = mean(penguins$flipper_length_mm), linewidth = 0.5, linetype = "dotted") +
  labs(title = "Covariance entre la longueur des nageoires et la masse",
       subtitle = "Contribution de chaque observation à la covariance",
       caption = "Les moyennes sont indiquées par les traits pointillés")

## 4.3 Matrice de variance-covariance ----

# Appliquée à un data.frame, la fonction cov() génère une matrice de
# variance-covariance : 

#   - variance d'une variable sur la diagonale
#   - covariance entre deux variables dans les autres cellules

penguins |> 
  # Sélectionner les variables numériques
  select(where(is.numeric)) |> 
  # Matrice de variance-covariance
  cov()

# La valeur de la covariance dépend des unités des variables : 

cov(x, y)       # longueur des nageoires en mm ; masse en g
cov(x/10, y)    # longueur des nageoires en cm ; masse en g  
cov(x, y/1000)  # longueur des nageoires en mm ; masse en kg

# La corrélation permet de s'affranchir des unités des variables.

# 5. Corrélation ----

## 5.1 Calcul de la corrélation ----

# Coefficient de corrélation de Pearson : 
#   - se situe entre -1 et 1
#   - le signe indique la direction de la corrélation
#   - une valeur absolue proche de 1 indique une très forte corrélation
#   - une valeur proche de 0 indique l'absence de corrélation

# r = cov(X, Y) / (sd(X) * sd(Y))

r_manuel <- cov(x, y) / (sd(x) * sd(y))

# Calcul de la corrélation à l'aide de la fonction cor()
cor(x, y, method = "pearson")  # Pearson par défaut

# Les deux valeurs sont-elles strictement égales ?
all.equal(r_manuel, cor(x, y))

# La corrélation est la covariance des variables centrées-réduites : 
cov(scale(x), scale(y))

## 5.2 Visualisation ----

penguins |> 
  ggplot() +
  geom_point(aes(x = scale(flipper_length_mm),
                 y = scale(body_mass_g)),
             color = "grey90") +
  geom_abline(slope = r_manuel) +
  annotate(geom = "segment", x = 0, xend = 1, y = 0, yend = 0,
           color = "red") +
  annotate(geom = "text", x = 0.5, y = -0.2, label = "1",
           color = "red", size = 6) +
  annotate(geom = "segment", x = 1, xend = 1, y = 0, yend = r_manuel,
           color = "red") +
  annotate(geom = "text", x = 1, y = 0.5, label = paste0("r = ", round(r_manuel, 2)),
           color = "red", size = 6, hjust = -0.1) +
  labs(title = "Corrélation (variables centrées-réduites)",
       x = "Flipper length (mm)",
       y = "Body mass (g)")

## 5.3 Test d'hypothèse du coefficient de corrélation ----

### 5.3.1 Formulation du test ----

# On teste le coefficient de corrélation de la population, noté rho.
# On estime rho par le coefficient observé r sur l'échantillon

# H0 : rho = 0 (pas de relation linéaire entre les deux variables)
# H1 : rho != 0 (test bilatéral), rho < 0 ou rho > 0 (test unilatéral)

# Sous l'hypothèse nulle H0, t suit une loi de Student à n-2 degrés de liberté

### 5.3.2 Calcul de la statistique de test ----

# La statistique de test est définie par :
# t = r * sqrt(n - 2) / sqrt(1 - r^2)

t_manuel <- r_manuel * sqrt(n - 2) / sqrt(1 - r_manuel^2)

### 5.3.3 Calcul de la p-value ----

# Calcul de la p-value associée (probabilité sous H0 d'observer une valeur
# au moins aussi extrême) : 

# q = -abs(t_manuel) ..... calcul de l'aire sous la courbe à gauche
# df = n - 2 ............. degrés de libertés de la distribution t
# on multiplie la valeur obtenue par 2 pour compter les deux queues (test bilatéral)

p_manuel <- 2 * pt(-abs(t_manuel), df = n - 2)

### 5.3.4 Calcul de l'intervalle de confiance de r ----

# La distribution de r n'est pas symétrique (elle est bornée entre -1 et 1).
# On ne peut donc pas utiliser la formule r +/- t * erreur type.

# On doit passer par la transformation de Fisher, qui rend la distribution 
# approximativement normale :

# 1) Transformation de Fisher
# z = 0.5 * log((1 + r) / (1 - r)) = arctanh(r)

z <- 0.5 * log((1 + r_manuel) / (1 - r_manuel))

all.equal(z, atanh(r_manuel))

# 2) Erreur-type de z (indépendante de rho)
se <- 1 / sqrt(n - 3)

# 3) Intervalle de confiance sur z
z_inf <- z - qnorm(p = 1 - (0.05 / 2)) * se
z_sup <- z + qnorm(p = 1 - (0.05 / 2)) * se

# 4) Retour à l'échelle de r avec la transformation inverse
r_lower <- tanh(z_inf)
r_upper <- tanh(z_sup)

### 5.3.5 Visualisation du résultat ----

ggplot(data = tibble(x = c(-35, 35)), aes(x)) +
  stat_function(fun = dt, args = list(df = n - 2)) +
  geom_rect(aes(xmin = -Inf, xmax = -t_manuel, ymin = 0, ymax = Inf),
            fill = "red", alpha = 0.2)  +
  geom_rect(aes(xmin = t_manuel, xmax = Inf, ymin = 0, ymax = Inf),
            fill = "red", alpha = 0.2)  +
  geom_vline(xintercept = c(-t_manuel, t_manuel),
             color = "red", linewidth = 0.6, linetype = "dashed") +
  labs(title = "Test d'hypothèse du coefficient de corrélation",
       subtitle = paste0("t = ", round(t_manuel, 2), " - pval = ", signif(p_manuel, 2)))

### 5.3.6 cor.test() ----

# Le test d'hypothèse du coefficient de corrélation peut se faire à l'aide
# de la fonction cor.test()

cor_test <- cor.test(x, y)
cor_test

# Les valeurs affichées sont-elles les mêmes que celles que nous avons calculées ?
# Nous utilisons dans certains cas la fonction unname() - quand il s'agit d'un vecteur nommé

all.equal(unname(cor_test$estimate), r_manuel)
all.equal(unname(cor_test$statistic), t_manuel)
all.equal(cor_test$p.value, p_manuel)
all.equal(cor_test$conf.int[1], r_lower)
all.equal(cor_test$conf.int[2], r_upper)

### 5.2.7 Visualisation des corrélations deux à deux ----

# Proposée dans le package {GGally}, la fonction ggpairs() permet 
# de visualiser :
#     - un density plot pour chaque variable quantitative
#     - pour chaque paire de variables quantitatives : 
#         - un scatter plot
#         - la corrélation

penguins |> 
  select(where(is.numeric)) |> 
  ggpairs()

# L'ajout d'une variable qualitative permet d'afficher :
#     - un bar plot pour chaque variable qualitative (facteur)
#     - un histogramme représentant la distribution de chaque variable quantitative
#       avec une séparation par classe présente dans la variable qualitative
#     - un boxplot, un histogramme et un density plot pour chaque variable qualitative
#     - pour chaque combinaison de deux variables : 
#         - un density plot
#         - la corrélation globale
#         - la corrélation par classe (en mappant la couleur à une variable)

penguins |> 
  # sélectionner les variables : `species` + variables quantitatives
  select(species, ends_with("mm"), body_mass_g) |> 
  # colorer la figure par la variable `species`
  ggpairs(aes(color = species)) +
  # appliquer une palette spécifique au contour (y compris les points)
  scale_color_manual(values = c("darkorange", "purple", "cyan4")) +
  # appliquer une palette spécifique à l'intérieur
  scale_fill_manual(values = c("darkorange", "purple", "cyan4"))


# Nous pouvons noter :
#     - une forte corrélation entre `body_mass_g` et `flipper_length_mm`
#     - des cas de figures où la corrélation globale est de signe différent
#       de la corrélation par espèce (on parle de "paradoxe de Simpson") : 
#         - `bill_length_mm` vs. `bill_depth_mm`
#         - `bill_depth_mm` vs. `flipper_length_mm`
#         - `bill_depth_mm` vs. `body_mass_g`



### 5.2.8 Matrice de corrélations ----

# {easystats} propose des packages pour les analyses statistiques : 
#     - {correlation} : analyses de corrélation
#     - {datawizard}  : exploration et préparation des données
#     - {effectsize}  : calcul de taille d'effet
#     - {performance} : évaluation d'un modèle
#     - {report}      : générer un rapport d'analyse
#     - ....

cor_table <- correlation::correlation(penguins)
cor_table

# Comparons avec les valeurs obtenues plus haut : 

cor_verif <- cor_table |> 
  filter(Parameter1 == "flipper_length_mm",
         Parameter2 == "body_mass_g")

cor_verif$r == r_manuel
cor_verif$t == t_manuel
cor_verif$p == p_manuel       # correlation() applique une correction au calcul de la p-value
cor_verif$CI_low == r_lower
cor_verif$CI_high == r_upper

# Par défaut, correlation() applique la correction de Holm au calcul de la p-value, et propose
# également les corrections suivantes : 
#   - "hochberg" ... Hochberg
#   - "hommel" ..... Hommel
#   - "bonferroni" .. Bonferroni
#   - "BH" .......... Benjamini-Hochberg
#   - "BY" .......... Benjamini-Yekutieli
#   - "fdr" ......... False Discovery Rate
#   - "somers" ...... Somers

# Pour retrouver les mêmes valeurs que celles calculées manuellement, nous
# utilisons l'argument `p_adjust = "none"` : 

correlation(penguins, p_adjust = "none")$p[6] == p_manuel

# L'ajustement selon Holm s'effectue de la manière suivante : 
# 1) Trier les p-values par ordre croissant de p_i à p_m
# 2) Multiplier la i-ième p-value par (m-i+1) : la plus petite par m, la suivante par m-1, ...
# 3) Imposer la monotonie : chaque p-value ajustée doit être supérieure ou égale à la précédente (maximum cumulé)
# 4) On plafonne à 1

iter <- tibble(var1 = penguins |> 
                 select(where(is.numeric)) |> 
                 names()
) |>
  mutate(var2 = var1) |> 
  complete(var1, var2) |> 
  filter_out(var1 == var2) |> 
  mutate(key1 = pmin(var1, var2),
         key2 = pmax(var1, var2)) |> 
  distinct(key1, key2, .keep_all = TRUE) |> 
  select(var1, var2)

p_values <- numeric(length = nrow(iter))

for (i in 1:nrow(iter)) {
  p_values[i] <- cor.test(penguins[[iter$var1[i]]], penguins[[iter$var2[i]]])$p.value
}

p_res <- iter |> 
  mutate(p = p_values) |> 
  arrange(p) |> 
  mutate(factor = seq(n(), 1),
         product = factor * p,
         p_holm = cummax(product))

all.equal(sort(cor_table$p), p_res$p_holm)

### 5.2.9 Heatmap des coefficients de corrélation ----

# Le package {see}, de la suite {easystats}, propose des fonctions 
# de visualisation :

cor_table |> 
  summary(redundant = FALSE) |> 
  plot() +
  theme_blackboard()


## 5.4 Les pièges de la corrélation ----

### 5.4.1 Le quartet d'Anscombe ----

# Ce jeu de données contient 4 paires de variables x,y ayant
# des statistiques descriptives et des coefficients de corrélation
# quasi identiques :

anscombe_long <- datasets::anscombe |> 
  pivot_longer(everything(), 
               names_to = c(".value", "jeu"),
               names_pattern = "(.)(.)") 

anscombe_long |> 
  summarise(moy_x = mean(x), 
            moy_y = mean(y), 
            r = cor(x, y),
            .by = jeu)

# La visualisation de ces paires de variables montrent des 
# données bien différentes ! 

anscombe_long |> 
  ggplot(aes(x, y)) +
  geom_point() +
  # geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~jeu)

# Ceci illustre l'importance fondamentale de ne pas se contenter de
# statistiques descriptives, et de visualiser les données.

### 5.4.2 Relation parfaite mais non linéaire ----

set.seed(42)
non_lineaire <- tibble(x = runif(200, -3, 3), y = x^2)

non_lineaire |> 
  ggplot(aes(x, y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

correlation(non_lineaire)  # r ~ 0 alors que y est fonction de x

### 5.4.3 Paradoxe de Simpson ----

# La corrélation entre `bill_length_mm` et `bill_depth_mm` est globalement
# négative ...

penguins |> 
  select(starts_with("bill")) |> 
  correlation()

# Mais devient positive si on la calcule pour chaque espèce !

penguins |> 
  group_by(species) |> 
  correlation(select = "bill_length_mm", select2 = "bill_depth_mm")

# Visualisation de ce paradoxe :

# Palette de couleurs {palmerpenguins}
my_cols <- c("Adelie" = "darkorange",
             "Chinstrap" = "purple",
             "Gentoo" = "cyan4")

penguins |> 
  ggplot(aes(x = bill_length_mm, y = bill_depth_mm)) +
  geom_point(aes(color = species)) +
  scale_color_manual(values = my_cols) +
  # corrélation globale négative
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  # corrélation positive par espèce
  geom_smooth(aes(color = species), method = "lm", se = FALSE,
              show.legend = FALSE) +
  labs(title = "Paradoxe de Simpson",
       subtitle = "Négatif globalement, positif au sein de chaque espèce")

# Une corrélation peut être due à une variable de confusion.
# La régression multiple est l'outil pour "contrôler" ce genre d'effet.

# 6. Analyse en Composantes Principales (ACP) ----

# L'objectif d'une ACP est de résumer p variables corrélées par q (< p) axes
# non corrélés (on parle de composantes principales) qui maximisent la variance
# projetée.

## 6.1 ACP à l'aide du package {FactoMineR} ----

pca_res <- PCA(select(penguins, where(is.numeric)))
# par défaut, PCA() standardise les données

# Valeurs propres et choix du nombre d'axes
pca_res$eig

# eigenvalue ....... valeur propre
# percentage of variance 
# cumulative percentage of variance

factoextra::fviz_screeplot(pca_res, addlabels = TRUE)

# Critères de choix du nombre d'axes :
#  - Kaiser (ACP normée) : garder val_propre > 1
#  - Coude de l'éboulis
#  - Bâton brisé : garder les axes dont pct > baton
#  - Interprétabilité, % cumulé (70-80 % indicatif)

# Ici PC1 ~ 69 % et PC2 ~ 19 % : 2 axes suffisent (~ 88 %).

p <- ncol(X)
eig <- tibble(
  axe        = factor(paste0("PC", 1:p)),
  val_propre = pca$sdev^2,                          # variance de chaque axe
  pct        = val_propre / sum(val_propre) * 100,  # % d'inertie
  pct_cum    = cumsum(pct),
  # modèle du bâton brisé : % attendu sous hypothèse d'absence de structure
  baton      = map_dbl(1:p, ~ sum(1 / (.x:p)) / p * 100)
)
eig


ggplot(eig, aes(axe, pct)) +
  geom_col(fill = "steelblue") +
  geom_line(aes(y = baton, group = 1), colour = "red") +
  geom_point(aes(y = baton), colour = "red") +
  geom_text(aes(label = round(pct, 1)), vjust = -0.5) +
  labs(y = "% d'inertie", title = "Éboulis + bâton brisé (rouge)")
# Ici PC1 ~ 69 % et PC2 ~ 19 % : 2 axes suffisent (~ 88 %).

# 3.3 Interprétation des variables
# - rotation  = vecteurs propres (loadings) : coefficients de la combinaison
# - corrélations variables-axes = loadings * sdev (ACP normée) : plus parlant
pca$rotation |> round(2)

cor_var_axes <- sweep(pca$rotation, 2, pca$sdev, `*`)
cor_var_axes |> round(2)
# PC1 : toutes les variables sauf bill_depth corrélées positivement
#       -> axe "TAILLE" globale (effet taille). Gentoo à droite.
# PC2 : essentiellement bill_depth_mm -> axe de forme du bec.

# Qualité de représentation des variables : cos2 = corrélation^2
(cor_var_axes^2)[, 1:2] |> round(2)

# 3.4 Interprétation des individus
scores <- as_tibble(pca$x) |> bind_cols(peng)

# Contribution des individus à l'axe k : x_ik^2 / ((n-1) * lambda_k)
n <- nrow(X)
contrib_ind <- as_tibble(sweep(pca$x^2, 2, (n - 1) * pca$sdev^2, `/`) * 100) |>
  rename_with(~ paste0("ctr_", .x))

# Qualité de représentation (cos2) d'un individu sur le plan 1-2
cos2_ind <- tibble(cos2_plan12 = (pca$x[, 1]^2 + pca$x[, 2]^2) / rowSums(pca$x^2))

scores <- bind_cols(scores, contrib_ind, cos2_ind)
scores |> arrange(desc(ctr_PC1)) |> select(id, species, ctr_PC1, cos2_plan12) |> head()
# Un individu mal représenté (cos2 faible) ne doit pas être interprété sur le plan.

# 3.5 Représentations graphiques

# Plan des individus, groupes (variable ILLUSTRATIVE : non utilisée dans le calcul)
ggplot(scores, aes(PC1, PC2, colour = species)) +
  geom_point(alpha = 0.7) +
  stat_ellipse(level = 0.9) +
  coord_fixed() +
  labs(x = sprintf("PC1 (%.1f %%)", eig$pct[1]),
       y = sprintf("PC2 (%.1f %%)", eig$pct[2]),
       title = "ACP normée - plan des individus")

# Cercle des corrélations (fait à la main avec ggplot2)
cercle <- tibble(t = seq(0, 2 * pi, length.out = 200),
                 x = cos(t), y = sin(t))
cor_df <- as_tibble(cor_var_axes, rownames = "variable")

ggplot() +
  geom_path(data = cercle, aes(x, y), colour = "grey60") +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
  geom_segment(data = cor_df, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), colour = "firebrick") +
  geom_text(data = cor_df, aes(PC1, PC2, label = variable), vjust = -0.6) +
  coord_fixed() +
  labs(x = "PC1", y = "PC2", title = "Cercle des corrélations")
# Lecture : flèche proche du cercle = bien représentée ; angle aigu = corrélation
# positive ; angle droit = pas de corrélation ; opposées = corrélation négative.

# Version rapide avec factoextra (biplot)
fviz_pca_biplot(pca, geom.ind = "point", col.ind = peng$species,
                addEllipses = TRUE, legend.title = "Espèce")

# 3.6 Variables supplémentaires avec FactoMineR
# Les variables qualitatives (espèce, île, sexe) ne participent pas au calcul
# des axes mais sont projetées (barycentres des modalités).
peng_df <- peng |>
  select(species, island, sex, all_of(vars_quanti)) |>
  as.data.frame()

res_pca <- PCA(peng_df, quali.sup = 1:3, scale.unit = TRUE, graph = FALSE)
plot(res_pca, choix = "ind", invisible = "ind", title = "Modalités supplémentaires")
res_pca$quali.sup$coord[, 1:2] |> round(2)
dimdesc(res_pca, axes = 1:2)   # description automatique des axes

# 3.7 Pourquoi standardiser ? Comparaison ACP non normée
pca_brut <- prcomp(X, center = TRUE, scale. = FALSE)
summary(pca_brut)$importance[2, ] |> round(3)
# body_mass_g (en grammes, variance énorme) écrase tout : PC1 ~ 100 % !
# -> ACP non normée uniquement si les variables ont la même unité/échelle.

# 3.8 Rappel des hypothèses/limites
#  - relations LINEAIRES entre variables
#  - sensible aux valeurs atypiques
#  - ACP = exploration (pas de test) ; ne pas interpréter les axes faibles
#  - ne pas confondre : axes d'ACP (variance max) != axes séparant les groupes
#    (c'est le rôle de la LDA, section 7)


# 4. DISTANCES MIXTES ET PCoA ----------------------------------------------------
# PCoA (ou MDS métrique) : représentation euclidienne à partir de n'importe
# quelle matrice de distances. Utile avec des variables mixtes -> Gower.

dat_mixte <- peng |> select(island, sex, all_of(vars_quanti))

d_gower <- daisy(dat_mixte, metric = "gower")   # gère quanti + quali
pcoa <- cmdscale(d_gower, k = 2, eig = TRUE)

tibble(Dim1 = pcoa$points[, 1], Dim2 = pcoa$points[, 2],
       species = peng$species, sex = peng$sex) |>
  ggplot(aes(Dim1, Dim2, colour = species, shape = sex)) +
  geom_point(alpha = 0.7) +
  labs(title = "PCoA sur distance de Gower")
# Avec une distance euclidienne sur données standardisées, la PCoA donne
# EXACTEMENT les mêmes coordonnées que l'ACP normée (à un signe près).
# Attention aux valeurs propres négatives si la distance n'est pas euclidienne.


# 5. CLASSIFICATION NON SUPERVISEE -------------------------------------------------
# On cherche des groupes "naturels" SANS utiliser l'espèce, puis on compare.

X_sc <- scale(X)                     # standardisation indispensable
d_euc <- dist(X_sc)                  # distance euclidienne

# 5.1 CAH (classification ascendante hiérarchique), critère de Ward
hc <- hclust(d_euc, method = "ward.D2")
# Ward minimise la perte d'inertie intra-classe à chaque fusion.
plot(hc, labels = FALSE, hang = -1, main = "Dendrogramme (Ward)")
rect.hclust(hc, k = 3, border = "red")

# Choix de k : saut d'inertie, silhouette, gap statistic
fviz_nbclust(X_sc, FUNcluster = hcut, method = "silhouette")
fviz_nbclust(X_sc, FUNcluster = hcut, method = "wss")

cluster_cah <- cutree(hc, k = 3)

# 5.2 k-means : toujours fixer la graine ET utiliser plusieurs départs
set.seed(42)
km <- kmeans(X_sc, centers = 3, nstart = 25)
km$tot.withinss / km$totss   # part d'inertie intra (plus petit = mieux)

# 5.3 Validation : silhouette (proche de 1 = bien classé, < 0 = mal classé)
sil <- silhouette(km$cluster, d_euc)
summary(sil)$avg.width
fviz_silhouette(sil)

# 5.4 Confrontation avec l'espèce réelle (table de contingence)
table(espece = peng$species, cluster_kmeans = km$cluster)
table(espece = peng$species, cluster_CAH = cluster_cah)
# Les clusters recoupent les espèces, mais Adelie/Chinstrap se mélangent :
# les variables morphologiques ne séparent pas parfaitement ces deux espèces.
# (Indice de Rand ajusté : mclust::adjustedRandIndex pour quantifier.)

# 5.5 Visualisation des clusters dans le plan de l'ACP
scores |>
  mutate(cluster = factor(km$cluster)) |>
  ggplot(aes(PC1, PC2, colour = cluster, shape = species)) +
  geom_point(size = 2, alpha = 0.8) +
  coord_fixed()


# 6. TESTS MULTIVARIES : MANOVA ET PERMANOVA -----------------------------------
# Question : les vecteurs moyens des groupes diffèrent-ils (sur l'ensemble des
# variables simultanément) ?

# 6.1 MANOVA : H0 : mu_Adelie = mu_Chinstrap = mu_Gentoo (vecteurs)
fit_manova <- manova(cbind(bill_length_mm, bill_depth_mm,
                           flipper_length_mm, body_mass_g) ~ species,
                     data = peng)
summary(fit_manova, test = "Pillai")   # Pillai : le plus robuste
tidy(fit_manova)                       # sortie en tibble (broom)

# Hypothèses : indépendance, normalité multivariée par groupe, homogénéité des
# matrices de covariance (test de Box, package heplots::boxM, très sensible),
# pas de multicolinéarité parfaite.

# 6.2 Tests post-hoc : ANOVA variable par variable (si MANOVA significative)
summary.aov(fit_manova)
# Corriger les tests multiples (Bonferroni/Holm) :
peng |>
  pivot_longer(all_of(vars_quanti), names_to = "variable", values_to = "valeur") |>
  group_by(variable) |>
  group_modify(~ tidy(aov(valeur ~ species, data = .x))) |>
  filter(term == "species") |>
  ungroup() |>
  mutate(p_ajustee = p.adjust(p.value, method = "holm"))

# 6.3 PERMANOVA : alternative sans hypothèse de normalité, sur toute distance
adonis2(d_euc ~ species, data = peng, permutations = 999)
# R2 = part de la variation (inertie) expliquée par l'espèce.

# La PERMANOVA est sensible aux différences de DISPERSION : on la vérifie.
disp <- betadisper(d_euc, peng$species)
permutest(disp, permutations = 999)   # p < 0.05 => dispersions hétérogènes
boxplot(disp)


# 7. ANALYSE DISCRIMINANTE LINEAIRE (LDA) --------------------------------------
# Méthode SUPERVISEE : cherche les combinaisons linéaires des variables qui
# maximisent le rapport variance inter-groupes / variance intra-groupes.
# À la fois descriptive (axes discriminants) et prédictive (classement).

dat_lda <- peng |> select(species, all_of(vars_quanti))

fit_lda <- MASS::lda(species ~ ., data = dat_lda)
fit_lda                      # probabilités a priori, moyennes, coefficients
fit_lda$svd^2 / sum(fit_lda$svd^2)   # part de la séparation portée par LD1, LD2

# Représentation dans le plan discriminant
pred <- predict(fit_lda)
tibble(as_tibble(pred$x), species = peng$species) |>
  ggplot(aes(LD1, LD2, colour = species)) +
  geom_point(alpha = 0.7) +
  stat_ellipse() +
  labs(title = "LDA - plan discriminant")
# Comparer à l'ACP : la LDA sépare nettement mieux les groupes (c'est son rôle).

# Performance : l'erreur de resubstitution est OPTIMISTE -> validation croisée
# leave-one-out (CV = TRUE)
lda_cv <- MASS::lda(species ~ ., data = dat_lda, CV = TRUE)
conf <- table(observé = dat_lda$species, prédit = lda_cv$class)
conf
sum(diag(conf)) / sum(conf)   # taux de bon classement en validation croisée

# Hypothèses : normalité multivariée intra-groupe, matrices de covariance
# égales (sinon : QDA, MASS::qda), n >> p, pas de variables colinéaires.
# Pour un vrai projet : découpage apprentissage/test ou k-fold répété.


# 8. AFDM : QUANTITATIFS ET QUALITATIFS -----------------------------------------
# Analyse factorielle de données mixtes (FAMD) = ACP + ACM combinées.
# Les variables quanti sont standardisées, les quali sont codées en disjonctif
# complet et pondérées de façon à équilibrer leur influence.

famd_dat <- peng |>
  select(island, sex, all_of(vars_quanti), species) |>
  as.data.frame()

res_famd <- FAMD(famd_dat, sup.var = "species", graph = FALSE)
get_eigenvalue(res_famd)
fviz_famd_var(res_famd, repel = TRUE)
fviz_famd_ind(res_famd, habillage = "species", addEllipses = TRUE, geom = "point")

# Pour que des variables qualitatives seules : ACM (FactoMineR::MCA).
# Pour deux variables quali croisées : AFC (FactoMineR::CA) sur le tableau de
# contingence, test du chi2 associé.


# 9. PIEGES CLASSIQUES ET CHECK-LIST -------------------------------------------
# [ ] Données : NA traités ? Variables numériques bien typées ? Unités ?
# [ ] Standardisation justifiée (unités différentes => ACP normée) ?
# [ ] Outliers examinés (et pas supprimés mécaniquement) ?
# [ ] Colinéarité forte entre variables (redondance, double comptage) ?
# [ ] Nombre d'axes retenu justifié (Kaiser, coude, bâton brisé) ?
# [ ] Interprétation des axes : corrélations variables-axes ET cos2 ?
# [ ] Variables illustratives vs actives clairement distinguées ?
# [ ] Corrélation != causalité ; effet taille à repérer (PC1 souvent "taille") ?
# [ ] Attention à la structure en groupes : paradoxe de Simpson, corrélations
#     globales trompeuses (cf. section 1.4)
# [ ] Clustering : k choisi par plusieurs critères, stabilité vérifiée (nstart,
#     graines multiples), pas de "preuve" d'existence de groupes
# [ ] Tests : hypothèses vérifiées, corrections pour tests multiples
# [ ] Prédiction : validation croisée ou jeu test, jamais le taux sur
#     l'échantillon d'apprentissage seul
# [ ] Reproductibilité : set.seed(), sessionInfo(), script commenté et versionné

sessionInfo()
