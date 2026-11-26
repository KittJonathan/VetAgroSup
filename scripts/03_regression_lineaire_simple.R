# Régression linéaire simple

# Table des matières ----

# 0. Comment utiliser ce script
# 1. Installer et charger les packages
# 2. Importer les données
# 3. Covariance et variance
#     3.1 Calcul

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
# install.packages("GGally")

# ⌨️ Raccourci clavier pour insérer ou retirer un commentaire : 
#     `Ctrl` + `Shift` + `C`

# Charger les packages :

library(tidyverse)          # manipulation et visualisation des données
library(skimr)              # exploration des données
library(palmerpenguins)     # jeu de données
library(GGally)             # matrice de comparaisons

theme_set(theme_bw())       # thème ggplot2

# 2. Importer les données ----

# Nous importons les données nettoyées du script `00_exploration_des_donnees.R`
penguins <- readRDS("data/penguins_clean.rds")

# Si besoin, exécuter le code ci-dessous pour nettoyer les données
penguins <- palmerpenguins::penguins |> 
  # transformer la variable `year` en facteur
  mutate(year = factor(year)) |> 
  # supprimer les données manquantes
  drop_na()

# ⌨️ Raccourci clavier pour insérer une flèche d'assignation : # `Alt` + `-`

# Vérifier les données importées
glimpse(penguins)

# 3. Covariance et variance ----

## 3.1 Calcul de la covariance ----

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

## 3.2 Visualisation ----

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

## 3.3 Matrice de variance-covariance ----

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

# 4. Corrélation ----

## 4.1 Calcul de la corrélation ----

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

## 4.2 Visualisation ----

penguins |> 
  ggplot() +
  geom_point(aes(x = scale(flipper_length_mm),
                 y = scale(body_mass_g)),
             color = "grey90") +
  geom_abline(slope = r_manuel) +
  annotate(geom = "segment", x = 0, xend = 1, y = 0, yend = 0,
           color = "red") +
  annotate(geom = "text", x = 0.5, y = -0.2, label = "1",
           color = "red", size = 10) +
  annotate(geom = "segment", x = 1, xend = 1, y = 0, yend = r_manuel,
           color = "red") +
  annotate(geom = "text", x = 1, y = 0.5, label = paste0("r = ", round(r_manuel, 2)),
           color = "red", size = 10, hjust = -0.1) +
  labs(title = "Corrélation (variables centrées-réduites)",
       x = "Flipper length (mm)",
       y = "Body mass (g)")

## 4.3 Test d'hypothèse du coefficient de corrélation ----

### 4.3.1 Formulation du test ----

# On teste le coefficient de corrélation de la population, noté rho.
# On estime rho par le coefficient observé r sur l'échantillon

# H0 : rho = 0 (pas de relation linéaire entre les deux variables)
# H1 : rho != 0 (test bilatéral), rho < 0 ou rho > 0 (test unilatéral)

# Sous l'hypothèse nulle H0, t suit une loi de Student à n-2 degrés de liberté

### 4.3.2 Calcul de la statistique de test ----

# La statistique de test est définie par :
# t = r * sqrt(n - 2) / sqrt(1 - r^2)

t_manuel <- r_manuel * sqrt(n - 2) / sqrt(1 - r_manuel^2)

### 4.3.3 Calcul de la p-value ----

# Calcul de la p-value associée (probabilité sous H0 d'observer une valeur
# au moins aussi extrême) : 

# q = -abs(t_manuel) ..... calcul de l'aire sous la courbe à gauche
# df = n - 2 ............. degrés de libertés de la distribution t
# on multiplie la valeur obtenue par 2 pour compter les deux queues (test bilatéral)

p_manuel <- 2 * pt(-abs(t_stat), df = n - 2)

### 4.3.4 Calcul de l'intervalle de confiance de r ----

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

### 4.3.5 Visualisation du résultat ----

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

### 4.3.6 cor.test() ----

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

### 4.2.7 Visualisation des corrélations deux à deux ----

# Proposée dans le package {GGally}, la fonction ggpairs() permet 
# de visualiser :
#     - un density plot pour chaque variable quantitative
#     - pour chaque paire de variables quantitatives : 
#         - un scatter plot
#         - la corrélation

penguins |> 
  select(where(is.numeric)) |> 
  ggpairs()

# Proposée dans le package {GGally}, la fonction ggpairs() permet 
# de visualiser :
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
  ggpairs()

# L'ajout d'une variable qualitative permet d'afficher des figures et des
# coefficients de corrélation regroupés par classe pour cette variable :
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



### 4.2.8 Matrice de corrélations ----

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

### 4.2.9 Heatmap des coefficients de corrélation ----

# Le package {see}, de la suite {easystats}, propose des fonctions 
# de visualisation :

cor_table |> 
  summary(redundant = FALSE) |> 
  plot() +
  theme_blackboard()


## 4.4 Les pièges de la corrélation ----

### 4.4.1 Le quartet d'Anscombe ----

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

# Cependant, la visualisation de ces paires de variables montrent des 
# données bien différentes ! 

anscombe_long |> 
  ggplot(aes(x, y)) +
  geom_point() +
  # geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~jeu)

# Ceci illustre l'importance fondamentale de ne pas se contenter de
# statistiques descriptives, et de visualiser les données.

### 4.4.2 Relation parfaite mais non linéaire ----

set.seed(42)
non_lineaire <- tibble(x = runif(200, -3, 3), y = x^2)

correlation(non_lineaire)  # r ~ 0 alors que y est fonction de x

non_lineaire |> 
  ggplot(aes(x, y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

### 4.4.3 Paradoxe de Simpson ----

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

# Leçon : une corrélation brute peut être produite par une variable de confusion.
# La régression multiple est l'outil pour "contrôler" ce genre d'effet.

# Cela illustre le fait qu'une corrélation globale peut être due à une variable
# de confusion.
# La régression multiple est l'outil de choix pour "contrôler" ce genre d'effet.

# 5. Régression linéaire simple ----

## 5.1 Variable explicative quantitative ----

# Modèle :  y = b0 + b1 * x + e 

#   - y ..... variable à expliquer
#   - b0 .... intercept (valeur de y quand x = 0)
#   - b1 .... pente (changement dans y pour une augmentation d'une unité dans x)
#   - x ..... variable explicative
#   - e ..... erreurs (résidus)

# Question : de combien la masse varie-t-elle, en moyenne, pour +1 mm de nageoire ?

### 5.1.1 Estimation par moindres carrés ordinaires ----

# Moindres carrés ordinaires (MCO) ou Ordinary Least Squares (OLS) : 
# On cherche b0 et b1 qui minimisent la somme des carrés des erreurs

# Solution analytique : 
# b1 = cov(x, y) / var(x) 
# b0 = mean(y) - b1 * mean(x)

x <- penguins$flipper_length_mm
y <- penguins$body_mass_g

b1_manuel <- cov(x, y) / var(x)
b0_manuel <- mean(y) - b1_manuel * mean(x)

# Lien avec la corrélation : b1 = r * sd(y) / sd(x)
all.equal(b1_manuel, r_manuel * sd(y) /sd(x))

# Calcul des valeurs théoriques
y_pred <- b0_manuel + b1_manuel * x

# Calcul des résidus
residus <- y - y_pred

# Somme des carrés totale : écart entre les valeurs observées et la moyenne
SCT <- sum((y - mean(y))^2)

# Somme des carrés expliquée par le modèle : écart entre les valeurs théoriques et la moyenne
SCR <- sum((y_pred - mean(y))^2)

# Somme des carrés des résidus (non expliqué par le modèle)
SCE <- sum(residus^2)

# Vérification que SCT = SCR + SCE
all.equal(SCT, SCR + SCE)

# Le coefficient r² est la part de la variabilité expliquée par le modèle.
r2 <- SCR / SCT

# Ce coefficient est égal au carré du coefficient de corrélation
all.equal(r2, r_manuel^2)

# Body mass vs. flipper length ----

penguins |> 
  ggplot(aes(x = flipper_length_mm, y = body_mass_g)) +
  geom_point()

# Résolution matricielle

y <- penguins$body_mass_g

X <- matrix(
  cbind(
    rep(1, length(y)),
    penguins$flipper_length_mm),
  ncol = 2, nrow = length(y)
  )

head(X)
head(y)

n <- nrow(X)  # nombre d'observations
p <- ncol(X)  # nombre de paramètres à estimer

solve(t(X) %*% X) %*% t(X) %*% y

coefficients(lm(body_mass_g ~ flipper_length_mm, data = penguins))
