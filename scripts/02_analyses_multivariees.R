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
# library(skimr)              # exploration des données
# library(naniar)             # exploration des données manquantes
# library(palmerpenguins)     # jeu de données
# library(broom)              # mise en forme "tidy" des sorties de modèles
# library(GGally)             # ggpairs : matrice de nuages de points
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

## 3.1 Synthèse des variables quantitatives ----

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

# Les écarts-types et les unités diffèrent entre les variables : il faudra
# standardiser les variables (centrer-réduire) pour l'ACP.


# 1.3 Matrice de nuages de points + corrélations, colorée par espèce
peng |>
  ggpairs(columns = vars_quanti, mapping = aes(colour = species, alpha = 0.6))

# 1.4 Matrice de corrélation globale
X <- peng |> select(all_of(vars_quanti))
cor(X) |> round(2)
# Observation : forte corrélation positive masse-nageoire (~0.87) ;
# la corrélation longueur-profondeur du bec est NEGATIVE globalement
# alors qu'elle est POSITIVE au sein de chaque espèce
# -> paradoxe de Simpson ! Illustration :
ggplot(peng, aes(bill_length_mm, bill_depth_mm)) +
  geom_point(aes(colour = species)) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +            # global
  geom_smooth(aes(colour = species), method = "lm", se = FALSE) +       # par groupe
  labs(title = "Paradoxe de Simpson : tendance globale vs intra-espèce")


# 2. VALEURS ATYPIQUES MULTIVARIEES --------------------------------------------

# Un individu peut être normal sur chaque variable prise isolément mais
# atypique dans la combinaison. Distance de Mahalanobis :
# D2 = (x - mu)' S^-1 (x - mu), ~ chi2(p) si normalité multivariée.
peng <- peng |>
  mutate(maha = mahalanobis(X, center = colMeans(X), cov = cov(X)),
         atypique = maha > qchisq(0.975, df = length(vars_quanti)))

peng |> count(atypique)
peng |> filter(atypique) |> select(id, species, sex, all_of(vars_quanti), maha)
# Ne pas supprimer automatiquement : comprendre (erreur de saisie ? vrai
# individu extrême ? effet d'une structure en groupes ?).
# Ici, les "atypiques" globaux reflètent surtout la structure en 3 espèces.


# 3. ACP -----------------------------------------------------------------------
# Objectif : résumer p variables corrélées par q < p axes non corrélés
# (composantes principales) qui maximisent la variance projetée.

# 3.1 Réalisation : prcomp (SVD), center = TRUE, scale. = TRUE => ACP NORMEE
pca <- prcomp(X, center = TRUE, scale. = TRUE)
summary(pca)

# 3.2 Valeurs propres et choix du nombre d'axes
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

# Critères de choix du nombre d'axes :
#  - Kaiser (ACP normée) : garder val_propre > 1
#  - Coude de l'éboulis
#  - Bâton brisé : garder les axes dont pct > baton
#  - Interprétabilité, % cumulé (70-80 % indicatif)
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
