# =============================================================================
#  MODULE "DÉMARCHE SCIENTIFIQUE" - VetAgro Sup, 3e année
#  Révisions enseignant : analyses multivariées, ACP, modèles linéaires généralisés
#  Jeux de données : package {agridat} (données agronomiques réelles)
# -----------------------------------------------------------------------------
#  PLAN
#   0. Préparation (packages, options)
#   1. Rappels : variance-covariance et régression linéaire (résolution matricielle)
#        jeu de données : lasrosas.corn (rendement de maïs vs azote)
#   2. Analyses multivariées
#        2.1 Principe, standardisation, distances
#        2.2 ACP "à la main" (diagonalisation de la matrice de corrélation)
#        2.3 ACP avec prcomp / FactoMineR, aides à l'interprétation
#        2.4 Classification ascendante hiérarchique (CAH) sur les axes de l'ACP
#        2.5 Interaction génotype x milieu : lien ACP / AMMI (SVD)
#   3. Modèles linéaires généralisés
#        3.1 Rappels : lois de probabilité, famille exponentielle, fonctions de lien
#        3.2 GLM de Poisson (comptages)     : beall.webworms
#        3.3 Algorithme IRLS programmé à la main
#        3.4 Surdispersion : quasi-Poisson, binomiale négative
#        3.5 GLM binomial (proportions)     : crowder.seeds
#        3.6 GLM Gamma (variable positive asymétrique) : lasrosas.corn
#        3.7 Diagnostics et sélection de modèle
#   4. Pistes d'exercices pour les étudiants
#
#  CONSEIL D'USAGE : exécuter section par section (Ctrl+Entrée) et
#  regarder systématiquement str(), summary() et les graphiques.
# =============================================================================


# =============================================================================
# 0. PRÉPARATION
# =============================================================================

# Packages indispensables
pkgs_requis <- c("agridat", "MASS")
# Packages optionnels (le script fonctionne sans, avec moins de graphiques)
pkgs_optionnels <- c("FactoMineR", "factoextra", "DHARMa", "lme4")

for (p in pkgs_requis) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(agridat)
library(MASS)

# Teste la disponibilité d'un package optionnel
dispo <- function(p) requireNamespace(p, quietly = TRUE)
for (p in pkgs_optionnels) {
  if (!dispo(p)) message("Package optionnel non installé : ", p,
                         "  -> install.packages('", p, "')")
}

options(digits = 5, width = 100)
set.seed(2026)  # reproductibilité (simulations pédagogiques)

# Fonction utilitaire : vérifie que les colonnes attendues existent
verifie_colonnes <- function(df, cols) {
  manquantes <- setdiff(cols, names(df))
  if (length(manquantes) > 0)
    stop("Colonnes absentes : ", paste(manquantes, collapse = ", "),
         "\nColonnes disponibles : ", paste(names(df), collapse = ", "))
  invisible(TRUE)
}


# =============================================================================
# 1. RAPPELS : VARIANCE-COVARIANCE ET RÉGRESSION LINÉAIRE MATRICIELLE
# =============================================================================
# Jeu de données : lasrosas.corn (Las Rosas, Argentine) - agriculture de
# précision. Rendement en maïs (yield, qq/ha), dose d'azote (nitro, kg/ha),
# topographie (topo), valeur de brillance du sol (bv), etc. Années 1999 et 2001.

data(lasrosas.corn)
str(lasrosas.corn)
verifie_colonnes(lasrosas.corn, c("year", "yield", "nitro", "topo", "bv"))

maize <- subset(lasrosas.corn, year == 2001)   # une seule année
maize <- na.omit(maize[, c("yield", "nitro", "topo", "bv")])
cat("n =", nrow(maize), "parcelles\n")

# --- 1.1 Variance, covariance, corrélation : définitions --------------------

# Pour deux variables x et y observées sur n individus :
#   Var(x)      = 1/(n-1) * sum( (x_i - xbar)^2 )
#   Cov(x, y)   = 1/(n-1) * sum( (x_i - xbar)(y_i - ybar) )
#   Cor(x, y)   = Cov(x, y) / ( sd(x) * sd(y) )      ∈ [-1, 1]
# Cov > 0 : x et y varient dans le même sens ; la corrélation, elle, est
# sans unité (la covariance dépend des unités de mesure).

x <- maize$nitro
y <- maize$yield
n <- length(x)

var_x  <- sum((x - mean(x))^2) / (n - 1)
cov_xy <- sum((x - mean(x)) * (y - mean(y))) / (n - 1)
c(var_manuelle = var_x, var_R = var(x))
c(cov_manuelle = cov_xy, cov_R = cov(x, y))
c(cor_manuelle = cov_xy / (sd(x) * sd(y)), cor_R = cor(x, y))

# --- 1.2 Matrice de variance-covariance (cas multivarié) --------------------
# Pour p variables, Σ est la matrice p x p symétrique dont l'élément (j,k)
# vaut Cov(X_j, X_k) ; la diagonale contient les variances.
# En notation matricielle, avec Xc la matrice centrée (colonnes de moyenne 0) :
#       S = Xc' Xc / (n - 1)
X_vars <- as.matrix(maize[, c("yield", "nitro", "bv")])
Xc <- scale(X_vars, center = TRUE, scale = FALSE)   # centrage
S_manuelle <- t(Xc) %*% Xc / (n - 1)
S_manuelle
all.equal(S_manuelle, cov(X_vars), check.attributes = FALSE)  # doit être TRUE

# La matrice de corrélation R est la matrice de variance-covariance des
# variables CENTRÉES-RÉDUITES (c'est le point clé pour l'ACP normée, section 2).
R <- cov(scale(X_vars))
all.equal(R, cor(X_vars), check.attributes = FALSE)

# --- 1.3 Régression linéaire simple : formules "scalaires" -------------------
# Modèle : y_i = b0 + b1 * x_i + e_i,  e_i ~ N(0, sigma^2) indépendants
# Estimateurs des moindres carrés ordinaires (MCO) :
#       b1 = Cov(x, y) / Var(x)          b0 = ybar - b1 * xbar
b1_scal <- cov_xy / var_x
b0_scal <- mean(y) - b1_scal * mean(x)
c(b0 = b0_scal, b1 = b1_scal)

# --- 1.4 Régression linéaire : résolution MATRICIELLE -----------------------
# Écriture : y = X b + e, avec X la matrice du plan d'expérience (n x (p+1)),
# 1ère colonne = 1 (intercept).
# On minimise ||y - X b||^2 ; en annulant le gradient on obtient les
# ÉQUATIONS NORMALES :     (X'X) b = X'y     =>   b_hat = (X'X)^{-1} X'y
# Variance des estimateurs : Var(b_hat) = sigma^2 (X'X)^{-1}
# avec sigma2_hat = SCR / (n - p - 1)   (SCR = somme des carrés résiduelle)
X <- cbind(1, x)                          # matrice du plan (n x 2)
colnames(X) <- c("Intercept", "nitro")
XtX <- t(X) %*% X
Xty <- t(X) %*% y
beta_hat <- solve(XtX, Xty)               # préférable à solve(XtX) %*% Xty
beta_hat

# Résidus et variance résiduelle
y_hat  <- X %*% beta_hat
resid  <- y - y_hat
p      <- ncol(X)                         # nb de paramètres (ici 2)
sigma2 <- sum(resid^2) / (n - p)
V_beta <- sigma2 * solve(XtX)             # matrice de var-cov des estimateurs
se_beta <- sqrt(diag(V_beta))
t_stat  <- beta_hat / se_beta
p_val   <- 2 * pt(abs(t_stat), df = n - p, lower.tail = FALSE)
cbind(Estimation = beta_hat, ErreurStd = se_beta, t = t_stat, p = p_val)

# Comparaison avec lm : doit donner exactement les mêmes valeurs
mod_lm <- lm(yield ~ nitro, data = maize)
summary(mod_lm)$coefficients

# --- 1.5 Matrice chapeau, décomposition de la variance, R² ------------------
# Matrice chapeau H = X (X'X)^{-1} X' : projection orthogonale de y sur
# l'espace engendré par les colonnes de X (y_hat = H y).
# Ses éléments diagonaux h_ii = leviers (leverage) des observations.
H <- X %*% solve(XtX) %*% t(X)
all.equal(as.numeric(H %*% y), as.numeric(y_hat))
all.equal(diag(H), unname(hatvalues(mod_lm)))
sum(diag(H))                               # = nombre de paramètres (trace de H)

# Décomposition de la variabilité : SCT = SCE + SCR
SCT <- sum((y - mean(y))^2)                # totale
SCE <- sum((y_hat - mean(y))^2)            # expliquée par le modèle
SCR <- sum(resid^2)                        # résiduelle
c(SCT = SCT, SCE_plus_SCR = SCE + SCR)
R2 <- SCE / SCT
c(R2_manuel = R2, R2_lm = summary(mod_lm)$r.squared,
  r_au_carre = cor(x, y)^2)                # en régression simple R² = r²

# Statistique de Fisher : F = (SCE/(p-1)) / (SCR/(n-p))
F_stat <- (SCE / (p - 1)) / (SCR / (n - p))
c(F_manuel = F_stat, F_lm = unname(summary(mod_lm)$fstatistic[1]))

# --- 1.6 Lien avec la variance-covariance -----------------------------------
# Si on CENTRE x et y, l'intercept disparaît et b1 = (xc' yc) / (xc' xc),
# c'est-à-dire exactement Cov(x,y) / Var(x). Les deux écritures sont identiques.
xc <- x - mean(x); yc <- y - mean(y)
as.numeric(crossprod(xc, yc) / crossprod(xc))   # = b1

# --- 1.7 Régression multiple : même formule, X plus large -------------------
# Variables qualitatives : model.matrix() fabrique les indicatrices (codage
# "traitement" par défaut : une modalité de référence).
maize$topo <- factor(maize$topo)
Xm <- model.matrix(yield ~ nitro + bv + topo, data = maize)
head(Xm)
beta_m <- solve(t(Xm) %*% Xm, t(Xm) %*% maize$yield)
cbind(manuel = beta_m, lm = coef(lm(yield ~ nitro + bv + topo, data = maize)))

# --- 1.8 Vérification graphique des hypothèses du modèle linéaire -----------
# Hypothèses : linéarité, indépendance, homoscédasticité, normalité des résidus.
par(mfrow = c(2, 2)); plot(mod_lm); par(mfrow = c(1, 1))
# Ici la réponse "rendement" dépend aussi de la topographie : un modèle sans
# topo laisse de la structure dans les résidus. À discuter avec les étudiants.
boxplot(resid(mod_lm) ~ maize$topo, xlab = "Topographie", ylab = "Résidus")


# =============================================================================
# 2. ANALYSES MULTIVARIÉES
# =============================================================================
# 2.1 PRINCIPE
# -----------------------------------------------------------------------------
# Objectif : décrire un tableau individus x variables quantitatives (n x p)
# sans variable "à expliquer" (méthodes EXPLORATOIRES / non supervisées).
# Familles de méthodes :
#   - ACP  : variables quantitatives -> axes factoriels
#   - AFC  : deux variables qualitatives (tableau de contingence)
#   - ACM  : plusieurs variables qualitatives
#   - AFM  : plusieurs groupes de variables
#   - CAH, k-means : classification (regroupement d'individus)
# Jeu de données : australia.soybean (essai multi-environnements de soja,
# Australie) - on regroupe par génotype pour obtenir un profil de variables
# agronomiques par variété.

data(australia.soybean)
str(australia.soybean)
traits <- c("yield", "height", "lodging", "size", "protein", "oil")
verifie_colonnes(australia.soybean, c("gen", traits))

# Moyenne par génotype pour chaque caractère (on ignore les valeurs manquantes)
soja <- aggregate(australia.soybean[, traits],
                  by = list(gen = australia.soybean$gen),
                  FUN = mean, na.rm = TRUE)
rownames(soja) <- soja$gen
soja$gen <- NULL
soja <- soja[complete.cases(soja), ]
dim(soja)                                  # n génotypes x 6 caractères
summary(soja)

# Exploration préalable : TOUJOURS regarder les distributions et les
# corrélations avant de faire une ACP.
round(cor(soja), 2)
pairs(soja, pch = 19, col = "darkgreen", main = "Caractères moyens par génotype")

# --- 2.1b Pourquoi standardiser ? -------------------------------------------
# Les variables ont des unités et variances très différentes. Sans réduction,
# les variables à forte variance dominent l'ACP (ACP "non normée" sur la
# matrice de covariance). Avec réduction (ACP "normée" sur la matrice de
# corrélation), chaque variable pèse 1 dans l'inertie totale.
round(apply(soja, 2, var), 3)
# Inertie totale en ACP normée = p (nombre de variables)

# --- 2.1c Distance entre individus -------------------------------------------
# L'ACP normée repose sur la distance euclidienne entre individus calculée
# sur les variables centrées-réduites : d²(i,k) = sum_j ((z_ij - z_kj))²
Z <- scale(soja)                           # centrage-réduction
D <- dist(Z)                               # distances euclidiennes
round(as.matrix(D)[1:5, 1:5], 2)


# -----------------------------------------------------------------------------
# 2.2 ACP "À LA MAIN"
# -----------------------------------------------------------------------------
# L'ACP cherche des combinaisons linéaires u_k des variables (axes
# principaux), orthogonales, qui maximisent la variance projetée.
#   - On maximise Var(Z u) = u' R u   sous la contrainte u'u = 1
#   - Solution : u = vecteur propre de R associé à la plus grande valeur propre
#     lambda_1 ; puis lambda_2 pour l'axe 2 (orthogonal au 1er), etc.
#   - lambda_k = variance (inertie) portée par l'axe k.
#   - Somme des lambda_k = trace(R) = p (en ACP normée).
#   - Coordonnées des individus (composantes principales) : C = Z U
#   - Corrélation variable j - composante k : cor(z_j, c_k) = u_jk * sqrt(lambda_k)
n <- nrow(Z); p <- ncol(Z)
R <- cor(soja)                              # = cov(Z)
ed <- eigen(R)                              # diagonalisation d'une matrice symétrique
lambda <- ed$values                         # valeurs propres
U      <- ed$vectors                        # vecteurs propres (colonnes)
rownames(U) <- colnames(soja); colnames(U) <- paste0("PC", 1:p)

lambda
sum(lambda)                                 # = p = 6 : inertie totale
pct_inertie <- 100 * lambda / sum(lambda)
round(rbind(valeur_propre = lambda,
            pct_inertie   = pct_inertie,
            pct_cumule    = cumsum(pct_inertie)), 2)

# Orthogonalité des axes : U'U = I
round(t(U) %*% U, 10)

# Composantes principales (coordonnées des individus)
C <- Z %*% U
colnames(C) <- paste0("PC", 1:p)
# La variance de chaque composante vaut lambda_k, et elles sont non corrélées
round(apply(C, 2, var), 4)
round(cor(C), 10)

# Corrélations variables - composantes ("cercle des corrélations")
cor_var <- U %*% diag(sqrt(lambda))
colnames(cor_var) <- paste0("PC", 1:p)
all.equal(cor_var, cor(Z, C), check.attributes = FALSE)
round(cor_var[, 1:3], 2)

# --- 2.2b Décomposition en valeurs singulières (SVD) -------------------------
# Autre écriture : Z/sqrt(n-1) = P D Q' (SVD). Les colonnes de Q sont les
# vecteurs propres de R et les d_k² les valeurs propres : même résultat que
# eigen(), mais numériquement plus stable (c'est ce qu'utilise prcomp).
sv <- svd(Z / sqrt(n - 1))
all.equal(sv$d^2, lambda)


# -----------------------------------------------------------------------------
# 2.3 ACP AVEC prcomp ET AIDES À L'INTERPRÉTATION
# -----------------------------------------------------------------------------
acp <- prcomp(soja, center = TRUE, scale. = TRUE)
summary(acp)
# Attention : prcomp renvoie des écarts-types (sdev) -> valeurs propres = sdev^2
all.equal(acp$sdev^2, lambda)
# Les signes des axes sont arbitraires (u et -u sont solutions) : ne pas
# s'étonner si les graphiques sont "en miroir" par rapport à la version manuelle.

# --- 2.3a Combien d'axes retenir ? -------------------------------------------
# Critères usuels (à combiner, pas à appliquer mécaniquement) :
#   - Kaiser : valeur propre > 1 (axe qui "résume plus qu'une variable")
#   - Éboulis : chercher le "coude" dans le diagramme des valeurs propres
#   - Pourcentage cumulé d'inertie (ex. > 70-80 %)
#   - Interprétabilité des axes
barplot(lambda, names.arg = paste0("PC", 1:p), col = "steelblue",
        ylab = "Valeur propre", main = "Éboulis des valeurs propres")
abline(h = 1, lty = 2, col = "red")        # critère de Kaiser

# --- 2.3b Qualité de représentation et contributions ------------------------
# cos² (individu i, axe k) = c_ik² / ||z_i||²  : qualité de représentation
# Contribution (individu i, axe k) = c_ik² / (n * lambda_k) * 100
cos2_ind  <- C^2 / rowSums(C^2)
contr_ind <- 100 * sweep(C^2, 2, n * lambda, "/")
# Pour les variables : cos² = corrélation², contribution = u_jk² * 100
cos2_var  <- cor_var^2
contr_var <- 100 * U^2
round(contr_var[, 1:2], 1)                  # quelles variables construisent PC1/PC2 ?

# Individus les mieux représentés / qui contribuent le plus à l'axe 1
round(head(sort(contr_ind[, 1], decreasing = TRUE), 5), 1)

# --- 2.3c Représentations graphiques ----------------------------------------
# Cercle des corrélations (variables)
theta <- seq(0, 2 * pi, length.out = 200)
plot(cos(theta), sin(theta), type = "l", asp = 1, xlab = "PC1", ylab = "PC2",
     main = "Cercle des corrélations")
abline(h = 0, v = 0, lty = 3)
arrows(0, 0, cor_var[, 1], cor_var[, 2], length = 0.1, col = "red")
text(cor_var[, 1] * 1.1, cor_var[, 2] * 1.1, rownames(cor_var), col = "red")
# Lecture : flèches proches = variables corrélées ; opposées = anti-corrélées ;
# orthogonales = non corrélées ; flèche longue (proche du cercle) = bien représentée.

# Nuage des individus
plot(C[, 1], C[, 2], pch = 19, col = "darkgreen",
     xlab = sprintf("PC1 (%.1f %%)", pct_inertie[1]),
     ylab = sprintf("PC2 (%.1f %%)", pct_inertie[2]),
     main = "Plan factoriel des génotypes de soja")
text(C[, 1], C[, 2], rownames(C), pos = 3, cex = 0.5)
abline(h = 0, v = 0, lty = 3)

# Biplot (individus + variables)
biplot(acp, cex = c(0.5, 0.9), main = "Biplot ACP")

# Version FactoMineR / factoextra (si installés) : sorties plus riches
if (dispo("FactoMineR") && dispo("factoextra")) {
  library(FactoMineR); library(factoextra)
  res_pca <- PCA(soja, scale.unit = TRUE, graph = FALSE)
  print(res_pca$eig)
  print(fviz_pca_var(res_pca, col.var = "contrib"))
  print(fviz_pca_ind(res_pca, col.ind = "cos2", repel = TRUE))
}

# --- 2.3d Individus / variables SUPPLÉMENTAIRES ------------------------------
# Elles ne participent pas au calcul des axes mais sont projetées a posteriori.
# Ex. : projeter un nouveau génotype (profil moyen hypothétique)
nouveau <- soja[1, ] * 1.05
predict(acp, newdata = nouveau)[, 1:2]

# --- 2.3e Pièges classiques à discuter --------------------------------------
# 1. ACP non normée vs normée : comparer les résultats
acp_nn <- prcomp(soja, center = TRUE, scale. = FALSE)
round(rbind(normee     = summary(acp)$importance[2, 1:3],
            non_normee = summary(acp_nn)$importance[2, 1:3]), 3)
round(acp_nn$rotation[, 1], 2)              # dominé par la variable à forte variance ?
# 2. Une corrélation visible sur le plan 1-2 n'a de sens que si les variables
#    sont bien représentées (cos² élevé).
# 3. L'ACP décrit, elle ne teste pas : pas de p-value, pas de causalité.
# 4. L'ACP ne capte que les relations LINÉAIRES.


# -----------------------------------------------------------------------------
# 2.4 CLASSIFICATION ASCENDANTE HIÉRARCHIQUE (CAH)
# -----------------------------------------------------------------------------
# Principe : partir de n classes (1 individu chacune) et fusionner à chaque
# étape les deux classes les plus proches. Critère de Ward : minimiser
# l'augmentation de l'inertie intra-classe (-> classes compactes).
# Bonne pratique : classer sur les premiers axes de l'ACP (bruit éliminé).
nb_axes <- sum(lambda > 1)                  # critère de Kaiser
cat("Axes retenus :", nb_axes, "\n")
coord <- acp$x[, 1:nb_axes]
hc <- hclust(dist(coord), method = "ward.D2")
plot(hc, cex = 0.5, main = "CAH (Ward) sur les composantes principales")
# Le "saut" d'inertie suggère le nombre de classes
barplot(rev(hc$height)[1:10], names.arg = 1:10,
        xlab = "Nb de classes", ylab = "Hauteur de fusion")
k <- 3
groupes <- factor(cutree(hc, k = k))
table(groupes)

# Visualisation des classes sur le plan factoriel
plot(acp$x[, 1], acp$x[, 2], col = groupes, pch = 19,
     xlab = "PC1", ylab = "PC2", main = "Classes de génotypes")
legend("topright", legend = levels(groupes), col = 1:k, pch = 19, title = "Classe")
# Caractérisation des classes : moyennes par classe des variables d'origine
round(aggregate(soja, by = list(classe = groupes), FUN = mean), 2)
# Alternative non hiérarchique : kmeans(coord, centers = k, nstart = 25)


# -----------------------------------------------------------------------------
# 2.5 INTERACTION GÉNOTYPE x MILIEU : LIEN ACP / AMMI
# -----------------------------------------------------------------------------
# Jeu de données : yan.winterwheat (blé d'hiver, Ontario) - rendement de
# génotypes (gen) dans des environnements (env). On forme le tableau
# gen x env et on étudie l'interaction.
# Modèle additif : y_ge = mu + a_g + b_e + (ab)_ge
# Le terme d'interaction est la partie résiduelle du tableau DOUBLEMENT CENTRÉ,
# qu'on résume par une SVD (c'est le principe du modèle AMMI).
data(yan.winterwheat)
str(yan.winterwheat)
verifie_colonnes(yan.winterwheat, c("gen", "env", "yield"))
M <- xtabs(yield ~ gen + env, data = yan.winterwheat)   # tableau gen x env
M <- unclass(M)
mu <- mean(M)
Mc <- sweep(sweep(M, 1, rowMeans(M)), 2, colMeans(M)) + mu   # double centrage
Mc <- Mc - mu                                                # = interaction pure
sv_ge <- svd(Mc)
pct_gxe <- 100 * sv_ge$d^2 / sum(sv_ge$d^2)
round(pct_gxe, 1)
# Scores des génotypes et des environnements sur l'axe 1 d'interaction
sc_gen <- sv_ge$u[, 1] * sqrt(sv_ge$d[1])
sc_env <- sv_ge$v[, 1] * sqrt(sv_ge$d[1])
names(sc_gen) <- rownames(M); names(sc_env) <- colnames(M)
# Biplot de l'interaction
plot(range(c(sc_gen, sc_env)), range(c(sv_ge$u[, 2] * sqrt(sv_ge$d[2]),
                                       sv_ge$v[, 2] * sqrt(sv_ge$d[2]))),
     type = "n", xlab = "IPC1", ylab = "IPC2", main = "Biplot de l'interaction G x E")
text(sc_gen, sv_ge$u[, 2] * sqrt(sv_ge$d[2]), names(sc_gen), col = "darkgreen", cex = 0.7)
text(sc_env, sv_ge$v[, 2] * sqrt(sv_ge$d[2]), names(sc_env), col = "red", cex = 0.7)
abline(h = 0, v = 0, lty = 3)
# Lecture : un génotype proche d'un environnement (même direction) y est
# favorisé ; un génotype proche de l'origine est stable.


# =============================================================================
# 3. MODÈLES LINÉAIRES GÉNÉRALISÉS (GLM)
# =============================================================================
# 3.1 RAPPELS THÉORIQUES
# -----------------------------------------------------------------------------
# Limites du modèle linéaire gaussien : il suppose y continue, de variance
# constante, et peut prédire des valeurs impossibles (comptages négatifs,
# probabilités > 1). Le GLM généralise le modèle en 3 composantes :
#
#  (1) COMPOSANTE ALÉATOIRE : y_i suit une loi de la famille exponentielle
#         f(y; theta, phi) = exp( (y*theta - b(theta)) / phi + c(y, phi) )
#      avec  E[y] = mu = b'(theta)   et   Var(y) = phi * b''(theta) = phi * V(mu)
#  (2) COMPOSANTE SYSTÉMATIQUE : prédicteur linéaire  eta = X b
#  (3) FONCTION DE LIEN g :  g(mu) = eta       (mu = g^{-1}(eta))
#
#  Loi        | Données                     | V(mu)       | Lien canonique | phi
#  -----------+-----------------------------+-------------+----------------+-----
#  Normale    | continues                   | 1           | identité       | sigma²
#  Poisson    | comptages (0,1,2,...)       | mu          | log            | 1
#  Binomiale  | succès sur n essais         | mu(1-mu)    | logit          | 1
#  Gamma      | continues > 0, asymétriques | mu²         | inverse (log   | à estimer
#             |                             |             | souvent utilisé)|
#
# Estimation : maximum de vraisemblance, par l'algorithme IRLS (section 3.3).
# Qualité d'ajustement : DÉVIANCE  D = 2 * [ logLik(modèle saturé) - logLik(modèle) ]
# Comparaison de modèles emboîtés : test du rapport de vraisemblance (LRT),
#   différence de déviances ~ chi² (phi connu) ou F (phi estimé).

# --- Visualisation des lois ---
par(mfrow = c(1, 3))
# Poisson : moyenne = variance = lambda
for (lam in c(1, 4, 10))
  plot(0:25, dpois(0:25, lam), type = "h", main = paste("Poisson, lambda =", lam),
       xlab = "y", ylab = "P(Y = y)", ylim = c(0, 0.4))
par(mfrow = c(1, 3))
# Binomiale : n essais, probabilité p. Moyenne np, variance np(1-p)
for (pp in c(0.1, 0.5, 0.8))
  plot(0:20, dbinom(0:20, 20, pp), type = "h", main = paste("Binomiale(20, p =", pp, ")"),
       xlab = "y", ylab = "P(Y = y)")
par(mfrow = c(1, 1))
# Gamma : forme a, taux b. Moyenne a/b, variance a/b² (= mu²/a : variance
# proportionnelle au carré de la moyenne -> coefficient de variation constant)
curve(dgamma(x, shape = 2, rate = 0.1), 0, 100, ylab = "densité", col = 1,
      main = "Lois Gamma de même moyenne ? (a/b = 20)")
curve(dgamma(x, shape = 8, rate = 0.4), 0, 100, add = TRUE, col = 2)
curve(dgamma(x, shape = 20, rate = 1), 0, 100, add = TRUE, col = 3)

# Fonctions de lien : logit et son inverse (la sigmoïde)
#  logit(p) = log( p / (1 - p) ) ∈ ]-inf, +inf[  ;  p = 1 / (1 + exp(-eta))
curve(plogis(x), -6, 6, xlab = "eta (prédicteur linéaire)", ylab = "p = g^-1(eta)",
      main = "Inverse du lien logit")

# Rappel d'interprétation des coefficients :
#  - lien log (Poisson)  : exp(b) = rapport multiplicatif de la moyenne
#  - lien logit (binom.) : exp(b) = odds ratio (rapport de cotes)

# Illustration : pourquoi ne pas faire de lm() sur des comptages ?
#   -> variance croissante avec la moyenne, résidus non normaux, prédictions < 0.


# -----------------------------------------------------------------------------
# 3.2 GLM DE POISSON : COMPTAGES (beall.webworms)
# -----------------------------------------------------------------------------
# Essai sur des vers de la betterave (webworms) : nombre de larves par parcelle
# (y), en fonction d'un traitement insecticide (spray) et d'un apport de plomb
# (lead), en blocs.
data(beall.webworms)
str(beall.webworms)
verifie_colonnes(beall.webworms, c("y", "block", "spray", "lead"))
ww <- beall.webworms
ww$block <- factor(ww$block)
ww$spray <- factor(ww$spray)
ww$lead  <- factor(ww$lead)

table(ww$y)                                 # distribution des comptages
c(moyenne = mean(ww$y), variance = var(ww$y))   # Poisson : moyenne ≈ variance ?
boxplot(y ~ spray:lead, data = ww, las = 2, ylab = "Nb de larves")
# Moyenne et variance par traitement : la variance augmente-t-elle avec la moyenne ?
aggregate(y ~ spray + lead, data = ww, FUN = function(v) c(m = mean(v), s2 = var(v)))

# Modèle de Poisson, lien log (canonique) :
#   y_i ~ Poisson(mu_i),   log(mu_i) = b0 + effet_bloc + effet_spray + effet_lead
mod_pois <- glm(y ~ block + spray * lead, family = poisson(link = "log"), data = ww)
summary(mod_pois)
# Lecture du summary :
#  - Coefficients sur l'échelle du LIEN (log)
#  - Null deviance : modèle avec intercept seul ; Residual deviance : modèle ajusté
#  - AIC : critère pour comparer des modèles (plus bas = meilleur compromis)
exp(coef(mod_pois))                         # effets multiplicatifs
# Prédictions : type = "link" (échelle de eta) ou "response" (échelle de mu)
nd <- expand.grid(block = levels(ww$block)[1], spray = levels(ww$spray),
                  lead = levels(ww$lead))
pr <- predict(mod_pois, newdata = nd, type = "link", se.fit = TRUE)
nd$mu    <- exp(pr$fit)
nd$ic_inf <- exp(pr$fit - 1.96 * pr$se.fit)       # IC construit sur l'échelle du lien
nd$ic_sup <- exp(pr$fit + 1.96 * pr$se.fit)       # puis retransformé
nd

# Tests des effets : analyse de déviance (séquentielle, type I)
anova(mod_pois, test = "Chisq")
# Test d'un effet précis par comparaison de modèles emboîtés :
mod_pois_add <- update(mod_pois, . ~ . - spray:lead)
anova(mod_pois_add, mod_pois, test = "Chisq")      # interaction significative ?

# Test de SURDISPERSION : rapport déviance résiduelle (ou chi² de Pearson) / ddl
# Valeur >> 1 => la variance dépasse celle de Poisson => erreurs standard trop
# petites, p-values trop optimistes.
phi_pearson <- sum(residuals(mod_pois, type = "pearson")^2) / df.residual(mod_pois)
c(phi_Pearson = phi_pearson,
  deviance_sur_ddl = deviance(mod_pois) / df.residual(mod_pois))


# -----------------------------------------------------------------------------
# 3.3 ALGORITHME IRLS PROGRAMMÉ À LA MAIN (Poisson, lien log)
# -----------------------------------------------------------------------------
# IRLS = Iteratively Reweighted Least Squares (équivalent à Newton-Raphson /
# scoring de Fisher). À chaque itération on résout une régression linéaire
# pondérée sur une "variable de travail" z :
#     eta = X b             mu = g^{-1}(eta)
#     z   = eta + (y - mu) * g'(mu)          (g' = dérivée du lien)
#     W   = diag( 1 / ( V(mu) * g'(mu)^2 ) )
#     b_nouveau = (X' W X)^{-1} X' W z
# Pour Poisson/log : g'(mu) = 1/mu, V(mu) = mu, donc W = diag(mu) et
#     z = eta + (y - mu)/mu
# Variance asymptotique de b : phi * (X' W X)^{-1}
irls_poisson <- function(X, y, tol = 1e-10, maxit = 50) {
  mu  <- y + 0.1                            # initialisation (éviter log(0))
  eta <- log(mu)
  b_old <- rep(0, ncol(X))
  for (it in 1:maxit) {
    z <- eta + (y - mu) / mu
    w <- mu
    XtWX <- t(X) %*% (X * w)
    b <- solve(XtWX, t(X) %*% (w * z))
    eta <- as.numeric(X %*% b)
    mu  <- exp(eta)
    if (max(abs(b - b_old)) < tol) break
    b_old <- b
  }
  list(coef = drop(b), vcov = solve(XtWX), iterations = it, mu = mu)
}
Xww <- model.matrix(y ~ block + spray * lead, data = ww)
res_irls <- irls_poisson(Xww, ww$y)
res_irls$iterations                          # convergence rapide (4-6 itérations)
# Comparaison avec glm
all.equal(unname(res_irls$coef), unname(coef(mod_pois)), tolerance = 1e-6)
all.equal(unname(sqrt(diag(res_irls$vcov))),
          unname(summary(mod_pois)$coefficients[, 2]), tolerance = 1e-6)

# Déviance de Poisson : D = 2 * sum( y*log(y/mu) - (y - mu) ), avec 0*log(0) = 0
dev_pois <- function(y, mu) 2 * sum(ifelse(y == 0, 0, y * log(y / mu)) - (y - mu))
c(manuelle = dev_pois(ww$y, res_irls$mu), glm = deviance(mod_pois))


# -----------------------------------------------------------------------------
# 3.4 SURDISPERSION : QUASI-POISSON ET BINOMIALE NÉGATIVE
# -----------------------------------------------------------------------------
# Causes fréquentes : hétérogénéité non modélisée (variables manquantes),
# agrégation spatiale des insectes, excès de zéros, non-indépendance.
# Solution 1 - quasi-Poisson : Var(y) = phi * mu, phi estimé (Pearson).
#   Les coefficients sont identiques ; les erreurs standard sont multipliées
#   par sqrt(phi) ; on utilise des tests F. Pas de vraisemblance -> pas d'AIC.
mod_qp <- update(mod_pois, family = quasipoisson)
summary(mod_qp)$dispersion
anova(mod_qp, test = "F")
# Solution 2 - binomiale négative : Var(y) = mu + mu²/theta (vraie vraisemblance
# donc AIC disponible ; variance quadratique en la moyenne).
mod_nb <- glm.nb(y ~ block + spray * lead, data = ww)
summary(mod_nb)
AIC(mod_pois, mod_nb)
# Comparer les erreurs standard : Poisson vs quasi-Poisson vs NB
round(cbind(Poisson = summary(mod_pois)$coefficients[, 2],
            QuasiPoisson = summary(mod_qp)$coefficients[, 2],
            NegBin = summary(mod_nb)$coefficients[, 2])[1:6, ], 3)

# Diagnostic par simulation de résidus (DHARMa), très utile en pédagogie
if (dispo("DHARMa")) {
  library(DHARMa)
  sim <- simulateResiduals(mod_pois)
  plot(sim)
  testDispersion(sim)
  testZeroInflation(sim)
}


# -----------------------------------------------------------------------------
# 3.5 GLM BINOMIAL : PROPORTIONS (crowder.seeds)
# -----------------------------------------------------------------------------
# Taux de germination de graines (Crowder, 1978) : 21 boîtes de Petri (plate),
# 2 variétés (gen) de graines de pois, 2 extraits racinaires (extract).
# germ = nb de graines germées sur n graines semées.
data(crowder.seeds)
str(crowder.seeds)
verifie_colonnes(crowder.seeds, c("plate", "gen", "extract", "germ", "n"))
cs <- crowder.seeds
cs$gen <- factor(cs$gen); cs$extract <- factor(cs$extract)
cs$prop <- cs$germ / cs$n
boxplot(prop ~ gen:extract, data = cs, ylab = "Proportion de graines germées")

# Modèle : germ_i ~ Binomiale(n_i, p_i),  logit(p_i) = b0 + gen + extract + gen:extract
# Réponse = matrice à 2 colonnes (succès, échecs).
# Alternative équivalente : glm(prop ~ ..., weights = n, family = binomial)
mod_bin <- glm(cbind(germ, n - germ) ~ gen * extract, family = binomial, data = cs)
summary(mod_bin)

# Interprétation : exp(coef) = odds ratios (rapports de cotes)
#  odds = p / (1 - p) ; si exp(b) = 1.5, la cote de germination est multipliée
#  par 1.5 (ce n'est PAS 1.5 fois la probabilité).
round(exp(cbind(OR = coef(mod_bin), confint.default(mod_bin))), 2)

# Probabilités prédites dans chaque condition
nd2 <- expand.grid(gen = levels(cs$gen), extract = levels(cs$extract))
nd2$p_hat <- predict(mod_bin, nd2, type = "response")
nd2

# Surdispersion en binomial : mêmes principes (phi > 1 ?)
phi_b <- sum(residuals(mod_bin, type = "pearson")^2) / df.residual(mod_bin)
phi_b
# Ici phi > 1 : variabilité entre boîtes de Petri supérieure à la variabilité
# binomiale pure (c'est le cas d'école de Crowder).
mod_qb <- update(mod_bin, family = quasibinomial)
anova(mod_qb, test = "F")
# Approche alternative, plus fidèle au dispositif : effet aléatoire boîte (GLMM)
if (dispo("lme4")) {
  library(lme4)
  mod_glmm <- glmer(cbind(germ, n - germ) ~ gen * extract + (1 | plate),
                    family = binomial, data = cs)
  print(summary(mod_glmm))
}
# (À présenter comme ouverture : GLM -> GLMM quand les données sont groupées.)

# Pourquoi pas un lm() sur la proportion ? Faire la comparaison : variance non
# constante (p(1-p)/n), pas de prise en compte des effectifs n_i, prédictions
# hors de [0,1] possibles.
mod_lm_prop <- lm(prop ~ gen * extract, data = cs)
plot(fitted(mod_lm_prop), resid(mod_lm_prop), main = "lm sur proportions : résidus")


# -----------------------------------------------------------------------------
# 3.6 GLM GAMMA : VARIABLE CONTINUE POSITIVE ASYMÉTRIQUE (lasrosas.corn)
# -----------------------------------------------------------------------------
# Quand la variabilité croît avec la moyenne (écart-type ∝ moyenne), un modèle
# Gamma avec lien log est souvent plus adapté que le gaussien. Ici : rendement
# du maïs en fonction de l'azote et de la topographie.
maize2 <- na.omit(subset(lasrosas.corn, year == 2001)[, c("yield", "nitro", "topo")])
maize2$topo <- factor(maize2$topo)
mod_gauss <- glm(yield ~ nitro + topo, family = gaussian, data = maize2)
mod_gamma <- glm(yield ~ nitro + topo, family = Gamma(link = "log"), data = maize2)
summary(mod_gamma)
# Dans un GLM Gamma, phi = 1/shape est estimé (c'est le "dispersion parameter").
# exp(b) = effet multiplicatif sur le rendement moyen (ex. +x% par unité d'azote).
round(100 * (exp(coef(mod_gamma)) - 1), 2)
AIC(mod_gauss, mod_gamma)
par(mfrow = c(1, 2))
plot(fitted(mod_gauss), residuals(mod_gauss, "pearson"), main = "Gaussien")
plot(fitted(mod_gamma), residuals(mod_gamma, "pearson"), main = "Gamma / log")
par(mfrow = c(1, 1))


# -----------------------------------------------------------------------------
# 3.7 DIAGNOSTICS ET SÉLECTION DE MODÈLE
# -----------------------------------------------------------------------------
# Types de résidus :
#  - "response"  : y - mu_hat (peu utile : variance non constante)
#  - "pearson"   : (y - mu_hat) / sqrt(V(mu_hat))
#  - "deviance"  : signe(y - mu) * sqrt(contribution à la déviance) (défaut)
# Pour un modèle bien spécifié, les résidus de déviance sont approx. centrés,
# sans structure vs valeurs ajustées.
par(mfrow = c(2, 2))
plot(mod_nb)                                # résidus vs ajustés, QQ, levier, Cook
par(mfrow = c(1, 1))
# Distance de Cook et leviers : identifier des observations influentes
cd <- cooks.distance(mod_nb)
head(sort(cd, decreasing = TRUE), 5)

# Sélection de modèle :
#  - Modèles emboîtés : test de déviance (anova ... test = "Chisq" ou "F")
#  - Modèles non emboîtés ou de familles différentes : AIC = -2 logL + 2 k
#    (BIC pénalise plus : -2 logL + k log n). Différence d'AIC < 2 : modèles
#    équivalents.
mod_ww0 <- glm.nb(y ~ block, data = ww)
mod_ww1 <- glm.nb(y ~ block + spray + lead, data = ww)
mod_ww2 <- glm.nb(y ~ block + spray * lead, data = ww)
AIC(mod_ww0, mod_ww1, mod_ww2)
BIC(mod_ww0, mod_ww1, mod_ww2)
anova(mod_ww0, mod_ww1, mod_ww2)

# Pseudo-R² de McFadden / de déviance (il n'existe pas de R² "vrai" en GLM)
1 - deviance(mod_pois) / mod_pois$null.deviance


# =============================================================================
# 4. PISTES D'EXERCICES POUR LES ÉTUDIANTS (groupe de 8)
# =============================================================================
# Chaque exercice se travaille en binôme ou en groupe de 4 ; chaque groupe
# présente démarche, résultats, limites. Insister sur la DÉMARCHE :
# question -> exploration -> modèle -> diagnostic -> interprétation -> limites.
#
# A. Régression linéaire (matrice)
#    A1. Sur un autre jeu agridat (ex. data(mead.cauliflower), data(gomez.fractionalfactorial),
#        ou lasrosas.corn 1999), recalculer (X'X)^{-1} X'y et comparer à lm().
#    A2. Montrer que Var(b1) = sigma2 / sum((x - xbar)^2) à partir de V_beta.
#    A3. Que se passe-t-il quand deux colonnes de X sont presque colinéaires ?
#        Simuler et regarder cond. de X'X (kappa(XtX)) et variances des b.
#
# B. ACP
#    B1. Refaire l'ACP de soja en excluant 'yield' puis projeter 'yield' en
#        variable supplémentaire : que conclure ?
#    B2. Comparer ACP normée / non normée sur un jeu où les unités diffèrent.
#    B3. Appliquer l'ACP à lasrosas.corn (lat, long, yield, nitro, bv) ou à un
#        tableau multi-caractères (data(cochran.lattice), data(jenkyn.mildew)
#        selon le cas) ; proposer une typologie des parcelles par CAH.
#    B4. Rédiger un paragraphe d'interprétation d'un axe (3-4 phrases) en
#        s'appuyant sur corrélations ET contributions.
#
# C. GLM
#    C1. Poisson sur beall.webworms : modèle sans/avec interaction, test de
#        surdispersion, choix entre quasi-Poisson et binomiale négative.
#    C2. crowder.seeds : calculer à la main un odds ratio et son IC de Wald
#        (exp(b ± 1.96 se)), vérifier avec confint.default().
#    C3. Simuler des comptages surdispersés (rnbinom), ajuster un Poisson,
#        et mesurer le taux de faux positifs sur un effet nul (simulation de
#        1000 jeux de données) : montrer l'inflation du risque de 1re espèce.
#    C4. Choisir une famille de GLM adaptée pour : (i) nombre de
#        plantes malades sur 20, (ii) masse de grains par épi, (iii) nombre de
#        pucerons par feuille, (iv) présence/absence d'une adventice. Justifier.
#
# D. Rédaction scientifique (fil rouge du module)
#    Pour un des modèles : Matériel et méthodes (modèle, famille, lien,
#    logiciel), Résultats (estimations + IC, pas seulement p-values),
#    Discussion des hypothèses vérifiées et limites.
# =============================================================================

sessionInfo()   # à joindre aux scripts rendus pour la reproductibilité
