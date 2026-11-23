# Analyses multivariées

# Table des matières ----

# 0. Packages et données
# 1. Exploration : structure, valeurs manquantes, corrélations 
# 2. Valeurs atypiques multivariées (distance de Mahalanobis)
# 3. ACP (Analyse en Composantes Principales)
# 4. Distances mixtes et PCoA (Gower)
# 5. Classification non supervisée (CAH et k-means)
# 6. Tests multivariés : MANOVA et PERMANOVA
# 7. Analyse discriminante linéaire (LDA)
# 8. AFDM : variables quantitatives ET qualitatives
# 9. Pièges classiques et check-list

# 0. Packages et données ----

library(tidyverse)       # dplyr, ggplot2, tidyr, purrr, tibble...
library(palmerpenguins)  # jeu de données `penguins`
library(broom)           # mise en forme "tidy" des sorties de modèles
library(GGally)          # ggpairs : matrice de nuages de points
library(FactoMineR)      # PCA, FAMD avec variables supplémentaires
library(factoextra)      # visualisations ggplot2 pour les analyses factorielles
library(cluster)         # daisy (Gower), silhouette
library(vegan)           # adonis2 (PERMANOVA), betadisper
# MASS est appelé avec MASS:: car MASS::select masque dplyr::select !