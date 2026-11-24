# Digital Soil Mapping with R
# https://zia207.quarto.pub/

# Table des matières ----

# 0. Comment utiliser ce script
# 1. Installer et charger les packages
# 2. Importer les données
# 3. Explorer les données
#     3.1 Afficher les données
#     3.2 Type et classe
#     3.3 Structure des données
#     3.4 Synthèse des données 
# 4. Explorer les données manquantes
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

# Pour importer des données, il y a deux façons de faire : 
#   - en "clique-bouton" : dans l'onglet "Environment" -> "Import Dataset"
#   - en ligne de commande :

soil <- read_csv("data/gp_soil_data_na.csv")

# ⌨️ Raccourci clavier pour insérer une flèche d'assignation : # `Alt` + `-`

# 3. Explorer les données ----

## 3.1 Afficher les données ----

# Afficher les données dans un onglet : 
View(soil)

# Ou cliquer sur l'icône à droite de `penguins` dans l'onglet "Environment".

# Dans cet onglet, il est possible de trier et filtrer les données, d'effectuer
# une recherche.

# Dans les versions récentes de RStudio, une synthèse des variables s'affiche : 
# histogramme pour les variables quantitatives, bar plot pour les variables
# qualitatives.

# Afficher les six premières lignes du tableau
head(soil)

# Afficher les quatre dernières lignes du tableau
tail(soil, n = 4)

# Afficher une sélection aléatoire de lignes du tableau
# Dans R, les tirages aléatoires sont en réalité des processus "pseudo-aléatoires",
# faisant appel à des listes prédéterminées de valeurs.
# Pour retrouver les mêmes résultats à chaque tirage "aléatoire", on fixe une valeur
# de référence dans la fonction `set.seed()` (ex : 123, 42, 2026, ...)
# Cela permet de "reproduire" des tirages aléatoires (ce qui peut sembler 
# contre-intuitif), afin de trouver les mêmes valeurs lorsqu'on exécute le code
# à plusieurs reprises.

set.seed(42)                    
slice_sample(soil, n = 5)

# ℹ️ Utilisation du pipe `|>`
# Le pipe permet d'enchaîner les commandes : la sortie d'une commmande est
# passée à la commande suivante. 
# Cela rend le code plus lisible quand les commandes s'enchaînent, 
# et permmet de construire un workflow de façon logique (notamment pour 
# la manipulation de données en enchaînant par exemple la sélection de colonnes,
# un filtre conditionnel sur une variable, le calcul d'une nouvelle variable, ...)

set.seed(42)
soil |> 
  slice_sample(n = 5)

# On peut trouver dans les ouvrages et articles le pipe sous la forme `%>%`.
# Il s'agit d'un pipe proposé dans le package {magrittr}.
# Depuis la version 4.1, R propose le pipe dir "natif" `|>`.

# ⌨️ Raccourci clavier pour insérer un pipe : `Ctrl` + `Shift` + `M`
# Pour définir le pipe à utiliser (`%>%` ou `|>`) : 
#     Tools > Global Options > Code > Editing > Use native pipe operator

## 3.2 Type et classe ----

typeof(soil)            
# liste (un data.frame est une liste particulière)

class(soil)             
# `data.frame` ...... tableau
# `tbl` ............. `tibble`, le data.frame spécifique au tidyverse
# `tbl_df` .......... `tibble` et `data.frame`
# `spec_tbl_df` ..... sous-classe de `tibble` créée lors de l'importation avec read_csv()

## 3.3 Structure des données ----

# Afficher la structure des données : 
str(soil)

# Cette commande retourne dans la console : 
#     - les dimensions du tableau (nb lignes * nb colonnes)
#     - les noms des variables précédés du symbole `$` (permet d'accéder à une variable)
#     - les types des variables (num, chr, ...)
#         - pour une variable qualitative : les niveaux (classes)
#         - pour une variable quantitative : le nombre d'observations
#     - les premiers éléments des variables

# ℹ️ Dans R, les variables qualitatives peuvent être recodées en facteurs.

as.factor(soil$STATE)

# Les niveaux (classes) d'un facteur sont par défaut ordonnés alphabétiquement :
levels(as.factor(soil$STATE))

# Les éléments d'un facteur sont stockés sous forme de valeurs numériques entières,
# chaque valeur correspondant à l'un des niveaux : 
as.integer(as.factor(soil$STATE))

# La fonction glimpse() est l'équivalent de la fonction str(), propre au tidyverse :
glimpse(soil)

## 3.4 Synthèse des données ----

# Afficher une synthèse détaillée des données
skim(soil)

# Cette commande affiche : 
#   - les dimensions du tableau
#   - la répartition des variables entre quantitatif/qualitatif
#   - une synthèse des variables qualitatives : données manquantes, longueur des chaînes de caractères, nombre de classes
#   - une synthèse des variables quantitatives : statistiques descriptives, histogramme

# La fonction summary() affiche une synthèse du tableau : 
summary(soil)

## 3.5 Explorer les données manquantes ----

# Compter le nombre de données manquantes par variable
colSums(is.na(soil))

### 3.5.1 Données manquantes par variable ----

# Afficher le nombre et la proportion de NA par variable
miss_var_summary(soil) |> 
  # conserver les variable ayant au moins une donnée manquante
  filter(n_miss != 0)

# pct_miss : 4 NA / 471 observations = 0.849%

# Table de contingence
miss_var_table(soil)

# Visualiser le nombre de NA par variable
gg_miss_var(soil)

### 3.5.2 Données manquantes par observation ----

# Afficher le nombre et la proportion de NA par observation
miss_case_summary(soil) |> 
  # conserver les observations ayant au moins une donnée manquante
  filter(n_miss != 0)

# pct_miss : 1 NA / 19 variables = 5.26%

soil |> 
  # ajouter une colonne "index"
  rowid_to_column() |>
  # conserver les observations (lignes) ayant au moins 1 NA
  filter(if_any(everything(), is.na))

# Table de contingence
miss_case_table(soil)

# pct_cases : 
# 467 observations complètes / 471 observations = 99.2%
# 4 observations avec 1 NA / 471 observations   = 0.849%

# Visualiser le nombre de données manquantes par observation
gg_miss_case(soil)

### 3.5.3 Répartition des données manquantes ----

# Nous ajoutons des données manquantes : 
#   - NA dans la colonne `Slope` quand NA dans la colonne `SOC`
#   - NA dans la colonne `DEM` quand ID est 95 ou 153

# Nous visualisons ensuite la répartition des données manquantes à l'aide
# d'un "upset plot" :

soil |> 
  mutate(
    Slope = replace_when(Slope, is.na(SOC) ~ NA),
    DEM = replace_when(DEM, ID %in% c(95, 153) ~ NA)
  ) |> 
  gg_miss_upset()

# Lecture : 
#   - Barplot horizontal à gauche du plot :
#       - 2 observations ont 1 NA dans la variable `DEM`
#       - 4 observations ont 1 NA dans la variable `SOC`
#   - Barplot vertical et liens entre variables au centre du plot :
#       - 2 observations ont 1 NA dans les variables `SOC` et `Slope`
#       - 2 observations ont 1 NA dans les variables `DEM`, SOC` et `Slope`

# 4. Nettoyer et exporter les données ----

# Suite à l'exploration des données, nous décidons de : 
# - uniformiser les noms de variables (minuscules)
# - supprimer les observations ayant des NA
# - transformer en facteurs les variables `ID`, `FIPS`, `STATE_ID`, `STATE`, 
# `COUNTY`, `NLCD` et `FRG`

soil_clean <- soil |> 
  # nettoyer les noms de variables (par défaut, la fonction janitor::clean_names() utilise des minuscules)
  janitor::clean_names() |> 
  # supprimer les données manquantes
  drop_na() |> 
  # transformer les variables qualitatives en facteurs
  mutate(across(.cols = c(id:county, nlcd, frg), 
                .fns = as.factor))

# Nous souhaitons exporter ce tableau de données nettoyées pour la suite
# des analyses.
# En exportant dans un format "classique" (.csv, .txt, ...), l'information
# des types de variables est perdue.

# Il existe plusieurs solutions pour exporter le tableau tout en conservant
# les métadonnées : 

# Exporter au format RDS
saveRDS(soil_clean, "data/soil_clean.rds")
readRDS("data/soil_clean.rds") # lit les données

# Exporter au format RData
save(soil_clean, file = "data/soil_clean.RData")
load("data/soil_clean.RData")  # charge les données dans l'environnement

# Exporter au format `parquet` (utile pour les données volumineuses)
arrow::write_parquet(soil_clean, "data/soil_clean.parquet")
arrow::read_parquet("data/soil_clean.parquet")

# Exporter au format `feather`
feather::write_feather(soil_clean, "data/soil_clean.feather")
feather::read_feather("data/soil_clean.feather")

# 5. Régression linéaire simple ----

# Un modèle de régression linéaire simple est décrit par la formule suivante :
# y = b0 + b1*x + e

# y ....... variable à expliquer
# x ....... variable explicative
# b0 ...... intercept
# b1 ...... pente
# e ....... erreur

# La pente b1 représente le changement dans la variable y pour une augmentation
# d'une unité de la variable x

# L'intercept b0 représente la valeur de la variable y quand x = 0

## 5.1 Matrice de variance-covariance ----

# La covariance mesure comment deux variables varient l'une par rapport à l'autre.

# cov(X, Y) = Σ (xi - x̄)(yi - ȳ) / (n - 1)
#   > 0 : quand X est au-dessus de sa moyenne, Y l'est en général aussi
#   < 0 : relation inverse
#   ≈ 0 : pas de relation LINÉAIRE

# Calcul manuel de la covariance de NDVI et SOC :

x <- soil_clean$ndvi
y <- soil_clean$soc  
n <- length(x)

cov_manuel <- sum((x - mean(x)) * (y - mean(y))) / (n - 1)
cov_manuel

# Calcul à l'aide de la fonction cov() :
cov(x, y)

# Les deux valeurs sont-elles strictement égales ?
all.equal(cov_manuel, cov(x, y))

# La variance est un cas particulier : var(X) = cov(X, X)
var(x)
cov(x, x)
