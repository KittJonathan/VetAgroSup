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
# install.packages("datawizard")
# install.packages("rstatix")
# install.packages("palmerpenguins")

# ⌨️ Raccourci clavier pour insérer ou retirer un commentaire : 
#     `Ctrl` + `Shift` + `C`

# Charger les packages :

library(tidyverse)          # manipulation et visualisation des données
library(skimr)              # exploration des données
library(datawizard)         # description des variables
library(rstatix)            # calculs et tests statistiques
# library(naniar)             # exploration des données manquantes
# library(arrow)              # exporter des données en conservant les métadonnées
# library(feather)            # exporter des données en conservant les métadonnées
# library(palmerpenguins)     # jeu de données

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

# 3. Exploration des variables ----

# La fonction skim() du package {skimr} affiche le détail des variables : 

skim(penguins)

## 3.1 Variables qualitatives ----

### 3.1.1 Synthèse ----

# Le package {datawizard} (de la suite {easystats}) propose des fonctions permettant
# la description des variables.

data_codebook(penguins |> select(where(is.factor)))

### 3.1.2 Fréquence des observations ----

# La fonction count() du package {dplyr} permet de compter les observations :
penguins |> 
  count(species)

# L'argument `sort = TRUE` permet de trier les données par fréquence décroissante 
penguins |> 
  count(species, sort = TRUE)

# Cette fonction est similaire à la fonction table()
table(penguins$species)

# Il est possible de calculer les fréquences pour des combinaisons de variables 
penguins |> 
  count(species, island)

# Par défaut, les combinaisons de variables n'ayant aucune observation sont 
# retirées de la table de contingence. L'argument `.drop = FALSE` permet 
# d'afficher toutes les combinaisons possibles
penguins |> 
  count(species, island, .drop = FALSE)

# Les données sont triées par ordre alphabétique

# Calculer manuellement les proportions par groupe
penguins |> 
  count(species, island, .drop = FALSE) |> 
  mutate(prop = n / sum(n),
         .by = species)

# Note : l'argument `.by` permet de grouper les données pour la fonction mutate()
#        Il est également possible d'utiliser la fonction group_by(), mais il 
#        faudra dans ce cas ajouter la fonction ungroup() à la suite pour éviter
#        de conserver une structure interne qui pourrait fausser les calculs suivants

penguins |> 
  count(species, island, .drop = FALSE) |> 
  group_by(species) |> 
  mutate(prop = n / sum(n)) |> 
  ungroup()

### 3.1.3 Visualisation des fréquences ----

# Le graphique de référence pour visualiser les fréquences des observations de
# variables qualitatives est le bar plot.
# (éviter le pie chart, la lecture des différences reposant sur des angles, 
# ce qui est plus compliqué que la comparaison de hauteurs de barres).

penguins |> 
  ggplot(aes(x = species)) +
  geom_bar()

# Par défaut, les classes sont triées par ordre alphabétique. Une bonne pratique
# de visualisation graphique consiste à ordonner les classes par ordre croissant
# (ou décroissant). On parle alors de diagramme de Pareto.

# Le package {forcats} (de la suite {tidyverse}) propose des fonctions permettant
# la manipulation des facteurs. Les noms de ces fonctions commencent par 'fct_' : 

penguins |> 
  ggplot(aes(x = fct_infreq(species))) +
  geom_bar()

penguins |> 
  ggplot(aes(x = fct_rev(fct_infreq(species)))) +
  geom_bar()

# Visualiser les fréquences de combinaisons de variables : 

# 1) Fréquences absolues
penguins |> 
  ggplot(aes(x = species, fill = island)) +
  geom_bar()

# 2) Fréquences relatives
penguins |> 
  ggplot(aes(x = species, fill = island)) +
  geom_bar(position = "fill")

## 3.2 Variables quantitatives ----

### 3.2.1 Synthèse ----

# La fonction get_summary_stats() du package {rstatix} permet de calculer
# les statistiques descriptives de variables quantitatives :

get_summary_stats(penguins)

# L'argument `type` permet de sélectionner les statistiques descriptives
# à afficher :

get_summary_stats(penguins, type = "common")
get_summary_stats(penguins, type = "robust")
get_summary_stats(penguins, type = "mean_sd") 
get_summary_stats(penguins, type = "five_number")

### 3.1.2 Description des variables ----

# La fonction describe_distribution() du package {datawizard} permet d'étudier
# la distribution des variables quantitatives :
describe_distribution(penguins)

#### 3.1.2.1 Symétrie de la distribution ----

# La skewness décrit la symétrie de la distribution :
# - skewness ~ 0 ...... distribution symétrique
# - skewness < 0 ....... distribution biaisée négativement
# - skewness > 0 ....... distribution biaisée positivement

# Une variable avec un fort biais négatif aura une moyenne inférieure à la médiane
# Une variable avec un fort biais positif aura une moyenne supérieure à la médiane

# Exemple : distribution de la variable `body_mass_g`

# Le package {datawizard} propose des fonctions permet de calculer la skewness et la kurtosis
datawizard::skewness(penguins$body_mass_g)

# Visualisation : 
penguins |> 
  ggplot(aes(x = body_mass_g)) +
  geom_density() +
  geom_vline(xintercept = mean(penguins$body_mass_g),
             color = "red") +
  annotate(geom = "text", x = mean(penguins$body_mass_g), y = 0, 
           label = "Moyenne", color = "red", hjust = -0.1) +
  geom_vline(xintercept = median(penguins$body_mass_g),
             color = "blue") +
  annotate(geom = "text", x = median(penguins$body_mass_g), y = 0, 
           label = "Médiane", color = "blue", hjust = 1.1) +
  labs(title = "Distribution de la masse",
       subtitle = paste0("Skew = ", round(skewness(penguins$body_mass_g), 2)))

#### 3.1.2.1 Forme de la distribution ----

# La kurtosis décrit l'épaisseur des queues la distribution (et pas la forme plus ou moins applatie du pic).

# Une distribution normale a une valeur de kurtosis ~ 3. 
# On compare généralement une distribution a cette valeur : on parlera alors d'excès de 
# kurtosis (~ 0 dans le cas d'une distribution normale)

# - excès de kurtosis ~ 0 ...... distribution mésokurtique (similaire à une distribution suivant la loi normale)
# - excès de kurtosis < 0 ....... distribution leptokurtique (queues plus épaisses qu'une distribution normale)
# - excès de kurtosis > 0 ....... distribution platykurtique (queues moins épaisses qu'une distribution normale)

# Visualisation d'un excès de kurtosis < 0

penguins |> 
  ggplot(aes(x = flipper_length_mm)) +
  geom_density() 
  

# 3. Tests ----

datawizard::describe_distribution(penguins)

ggplot(data = tibble(x = c(200.97 - 4 * 14.02, 200.97 + 4 * 14.02)), aes(x)) +
  stat_function(fun = dnorm, args = list(mean = 200.97, sd = 14.02),
                linewidth = 1, color = "steelblue") +
  geom_density(data = penguins, aes(x = flipper_length_mm))


penguins |> 
  filter(species == "Gentoo") |> 
  datawizard::describe_distribution()

ggplot(data = tibble(x = c(47.57 - 4 * 3.11 , 47.57 + 4 * 3.11)), aes(x)) +
  stat_function(fun = dnorm, args = list(mean = 47.57, sd = 3.11),
                linewidth = 1, color = "steelblue") +
  geom_density(data = filter(penguins, species == "Gentoo"), aes(x = bill_length_mm))
