# Analyses univariées

# Table des matières ----

# 0. Comment utiliser ce script
# 1. Installer et charger les packages
# 2. Importer les données


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

# ⌨️ Raccourci clavier pour insérer ou retirer un commentaire : 
#     `Ctrl` + `Shift` + `C`

# Charger les packages :

library(tidyverse)          # manipulation et visualisation des données
library(skimr)              # exploration des données
library(naniar)             # exploration des données manquantes
library(arrow)              # exporter des données en conservant les métadonnées
library(feather)            # exporter des données en conservant les métadonnées
library(palmerpenguins)     # jeu de données

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

# 3. Tests ----

datawizard::describe_distribution(penguins)

ggplot(data = tibble(x = c(4207.06 - 4 * 805.22, 4207.06 + 4 * 805.22)), aes(x)) +
  stat_function(fun = dnorm, args = list(mean = 4207.06, sd = 805.22),
                linewidth = 1, color = "steelblue") +
  geom_density(data = penguins, aes(x = body_mass_g))


penguins |> 
  filter(species == "Gentoo") |> 
  datawizard::describe_distribution()

ggplot(data = tibble(x = c(47.57 - 4 * 3.11 , 47.57 + 4 * 3.11)), aes(x)) +
  stat_function(fun = dnorm, args = list(mean = 47.57, sd = 3.11),
                linewidth = 1, color = "steelblue") +
  geom_density(data = filter(penguins, species == "Gentoo"), aes(x = bill_length_mm))
