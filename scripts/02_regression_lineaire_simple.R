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

# Nous créons une palette de couleurs personnalisée, reprenant les couleurs
# du package {palmerpenguins} : un vecteur nommé, construit en associant
# chaque espèce à une couleur.

my_cols <- c(
  "Adelie" = "darkorange",
  "Chinstrap" = "purple",
  "Gentoo" = "cyan4"
)

# ⌨️ Raccourci clavier pour insérer une flèche d'assignation : # `Alt` + `-`

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

# 3. Covariance et variance ----

## 3.1 Calcul ----

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

cov(x, y)

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

# 3. Exploration des données ----

penguins |> 
  # sélectionner les variables : `species` + variables quantitatives
  select(species, ends_with("mm"), body_mass_g) |> 
  # colorer la figure par la variable `species`
  ggpairs(aes(color = species)) +
  # appliquer une palette spécifique au contour (y compris les points)
  scale_color_manual(values = c("darkorange", "purple", "cyan4")) +
  # appliquer une palette spécifique à l'intérieur
  scale_fill_manual(values = c("darkorange", "purple", "cyan4"))


# La fonction ggpairs() du package {GGally} affiche : 
#     - un bar plot pour la variable `species`
#     - un boxplot, un histogramme et un density plot par variable 
#     - pour chaque combinaison de deux variables : 
#         - un density plot
#         - la corrélation globale
#         - la corrélation par espèce

# Nous pouvons noter :
#     - une forte corrélation entre `body_mass_g` et `flipper_length_mm`
#     - des cas de figures où la corrélation globale est de signe différent
#       de la corrélation par espèce (on parle de "paradoxe de Simpson") : 
#         - `bill_length_mm` vs. `bill_depth_mm`
#         - `bill_depth_mm` vs. `flipper_length_mm`
#         - `bill_depth_mm` vs. `body_mass_g`

penguins |> 
  group_by(species) |> 
  correlation::correlation()


correlation::correlation(penguins)

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
