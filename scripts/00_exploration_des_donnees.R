# Exploration des données

# Table des matières ----

# 0. Comment utiliser ce script
# 1. Installer et charger les packages
# 2. Importer les données
# 3. Explorer les données
#     3.1 Afficher les données
#     3.2 Type et classe
#     3.3 Structure des données
#     3.4 Synthèse des données 
#     3.5 Explorer les données manquantes
#     4.1 Données manquantes par variable
#     4.2 Données manquantes par observation
#     4.3 Répartition des données manquantes
# 5. Nettoyer les données

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
# install.packages("naniar")
# install.packages("arrow")
# install.packages("feather")
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

# Nous allons utiliser le jeu de données `penguins`, disponible dans le 
# packages {palmerpenguins}.

# ⚠️ Depuis la version 4.5.0 de R, un jeu de données `penguins` est disponible
# sans faire appel à un package externe. Cependant, les noms de certaines variables
# diffèrent de ceux utilisés dans le package {palmerpenguins}. Nous privilégions donc
# le jeu de données proposé dans ce package plutôt que le jeu de données 
# proposé par défaut dans R.

penguins <- palmerpenguins::penguins

# ⌨️ Raccourci clavier pour insérer une flèche d'assignation : # `Alt` + `-`

# 3. Explorer les données ----

## 3.1 Afficher les données ----

# Afficher les données dans un onglet : 
View(penguins)

# Ou cliquer sur l'icône à droite de `penguins` dans l'onglet "Environment".

# Dans cet onglet, il est possible de trier et filtrer les données, d'effectuer
# une recherche.

# Dans les versions récentes de RStudio, une synthèse des variables s'affiche : 
# histogramme pour les variables quantitatives, bar plot pour les variables
# qualitatives.

# Afficher les six premières lignes du tableau
head(penguins)

# Afficher les quatre dernières lignes du tableau
tail(penguins, n = 4)

# Afficher une sélection aléatoire de lignes du tableau

# Dans R, les tirages aléatoires sont en réalité des processus "pseudo-aléatoires",
# faisant appel à des listes prédéterminées de valeurs.
# Pour retrouver les mêmes résultats à chaque tirage "aléatoire", on fixe une valeur
# de référence dans la fonction `set.seed()` (ex : 123, 42, 2026, ...)
# Cela permet de "reproduire" des tirages aléatoires (ce qui peut sembler 
# contre-intuitif), afin de trouver les mêmes valeurs lorsqu'on exécute le code
# à plusieurs reprises.

set.seed(42)                    # pour la reproductibilité
slice_sample(penguins, n = 5)

# ℹ️ Utilisation du pipe `|>`
# Le pipe permet d'enchaîner les commandes : la sortie d'une commmande est
# passée à la commande suivante. 
# Cela rend le code plus lisible quand les commandes s'enchaînent, 
# et permmet de construire un workflow de façon logique (notamment pour 
# la manipulation de données en enchaînant par exemple la sélection de colonnes,
# un filtre conditionnel sur une variable, le calcul d'une nouvelle variable, ...)

set.seed(42)

penguins |> 
  slice_sample(n = 5)

# On peut trouver dans les ouvrages et articles le pipe sous la forme `%>%`.
# Il s'agit d'un pipe proposé dans le package {magrittr}.
# Depuis la version 4.1, R propose le pipe dir "natif" `|>`.

# ⌨️ Raccourci clavier pour insérer un pipe : `Ctrl` + `Shift` + `M`
# Pour définir le pipe à utiliser (`%>%` ou `|>`) : 
#     Tools > Global Options > Code > Editing > Use native pipe operator

## 3.2 Type et classe ----

typeof(penguins)            
# liste (un data.frame est une liste particulière)

class(penguins)             
# data.frame = tableau
# tbl = `tibble`, le data.frame spécifique au tidyverse
# tbl_df = `tibble` et `data.frame`

## 3.3 Structure des données ----

# Afficher la structure des données : 
str(penguins)

# Cette commande retourne dans la console : 
#     - les dimensions du tableau : nb de lignes (observations) x nb de colonnes (variables) 
#     - les noms des variables précédés du symbole `$` (permet d'accéder à une variable)
#     - les types des variables (Factor, num, int, ...)
#         - pour une variable qualitative : les niveaux (classes)
#         - pour une variable quantitative : le nombre d'observations
#     - les premiers éléments des variables

# ℹ️ Dans R, les variables qualitatives peuvent être recodées en facteurs.
# La commande str() retourne dans le cas d'une variable de type facteur les
# niveaux (qui sont par défaut ordonnés alphabétiquement) et affiche les éléments
# d'un facteur sous forme de valeur numérique correspondant à l'index du niveau.

head(penguins$species)               # affiche les 6 premiers éléments de la variable et les 3 niveaux
typeof(penguins$species)             # integer : les éléments sont stockés sous forme de nombre
class(penguins$species)              # facteur
as.integer(head(penguins$species))   # affiche les valeurs numériques correspondantes

# Equivalent de la fonction str(), mais spécifique au tidyverse
glimpse(penguins)

## 3.4 Synthèse des données ----

# Afficher une synthèse détaillée des données
skim(penguins)

# Cette commande affiche : 
#   - les dimensions du tableau
#   - la répartition des variables entre quantitatif/qualitatif
#   - une synthèse des variables qualitatives (facteurs) : données manquantes, classes, fréquences
#   - une synthèse des variables quantitatives : statistiques descriptives, histogramme

# La fonction summary() affiche une synthèse du tableau : 
summary(penguins)

## 3.5 Explorer les données manquantes ----

# Compter le nombre de données manquantes par variable
colSums(is.na(penguins))

### 3.5.1 Données manquantes par variable ----

# Afficher le nombre et la proportion de NA par variable
miss_var_summary(penguins) |> 
  # conserver les variable ayant au moins une donnée manquante
  filter(n_miss != 0)

# pct_miss : 
# 11 NA / 344 observations = 3.20%
# 2 NA / 344 observations = 0.581%

# Table de contingence
miss_var_table(penguins)
# 3 variables n'ont aucune donnée manquante (3/8 = 37.5%)
# 4 variables ont 2 données manquantes (4/8 = 50%)
# 1 variable a 11 données manquantes (1/8 = 12.5%)

# Visualiser le nombre de NA par variable
gg_miss_var(penguins)

### 3.5.2 Données manquantes par observation ----

# Afficher le nombre et la proportion de NA par observation
miss_case_summary(penguins) |> 
  # conserver les observations ayant au moins une donnée manquante
  filter(n_miss != 0)
# pct_miss : 
# 5 NA / 8 variables = 62.5%
# 1 NA / 8 variables = 12.5%

penguins |> 
  # ajouter une colonne "index"
  rowid_to_column() |>
  # conserver les observations (lignes) ayant au moins 1 NA
  filter(if_any(everything(), is.na))

# Table de contingence
miss_case_table(penguins)
# 333 observations n'ont aucune donnée manquante (333/344 = 96.8%)
# 9 observations ont 1 donnée manquante (9/344 = 2.62%)
# 2 observations ont 2 données manquantes (2/344 = 0.581%)

# Visualiser le nombre de données manquantes par observation
gg_miss_case(penguins)

### 3.5.3 Répartition des données manquantes ----

# Visualiser la répartition des données manquantes à l'aide
# d'un "upset plot" :

gg_miss_upset(penguins)

# Lecture : 
#   - Barplot horizontal à gauche du plot :
#       - 2 observations ont 1 NA dans les variables `bill_length_mm`, `bill_depth_mm`, 
#         `flipper_length_mm` et `body_mass_g`
#       - 11 observations ont 1 NA dans la variable `sex`
#   - Barplot vertical et liens entre variables au centre du plot :
#       - 9 observations ont 1 NA dans la variable `sex`
#       - 2 observations ont 1 NA dans les variables `bill_length_mm`, `bill_depth_mm`, 
#         `flipper_length_mm`, `body_mass_g` et `sex` 


# 5. Nettoyer les données ----

# Suite à l'exploration des données, nous décidons de : 
# - transformer la variable `year` en facteur
# - supprimer les observations ayant des NA

penguins_clean <- penguins |> 
  # transformer la variable `year` en facteur
  mutate(year = factor(year)) |> 
  # supprimer les données manquantes
  drop_na()

# Nous souhaitons exporter ce tableau de données nettoyées pour la suite
# des analyses.
# En exportant dans un format "classique" (.csv, .txt, ...), l'information
# des types de variables est perdue, et les variables concernées devront
# à nouveau être modifiées

# Il existe plusieurs solutions pour exporter le tableau tout en conservant
# les métadonnées : 

# Exporter au format RDS
saveRDS(penguins_clean, "data/penguins_clean.rds")
readRDS("data/penguins_clean.rds") # lit les données

# Exporter au format RData
save(penguins_clean, file = "data/penguins_clean.RData")
load("data/penguins_clean.RData")  # charge les données dans l'environnement

# Exporter au format `parquet` (utile pour les données volumineuses)
arrow::write_parquet(penguins_clean, "data/penguins_clean.parquet")
arrow::read_parquet("data/penguins_clean.parquet")

# Exporter au format `feather`
feather::write_feather(penguins_clean, "data/penguins_clean.feather")
feather::read_feather("data/penguins_clean.feather")
