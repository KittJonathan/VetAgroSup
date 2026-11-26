library(tidyverse)
library(palmerpenguins)

# ------------------------------------------------------------------
# 1. Données : on retire les NA et on standardise x ET y
#    (surface de loss bien "ronde" -> convergence facile et comparable)
# ------------------------------------------------------------------
df <- penguins |>
  drop_na(flipper_length_mm, body_mass_g)

x_raw <- df$flipper_length_mm
y_raw <- df$body_mass_g

x <- as.numeric(scale(x_raw))
y <- as.numeric(scale(y_raw))
n <- length(x)

# Solution analytique (référence) dans l'espace standardisé
w_ols <- cov(x, y) / var(x)           # = cor(x, y)
b_ols <- mean(y) - w_ols * mean(x)    # = 0

# ------------------------------------------------------------------
# 2. Fonctions : MSE et descente de gradient générique
#    batch_size = n   -> descente de gradient "batch"
#    batch_size = 1   -> SGD (stochastique)
#    1 < batch_size < n -> mini-batch
# ------------------------------------------------------------------
mse <- function(w, b, x, y) mean((y - (b + w * x))^2)

gradient_descent <- function(x, y, w, b, lr, epochs, batch_size, seed = 1) {
  set.seed(seed)
  n <- length(x)
  n_batches <- ceiling(n / batch_size)
  n_steps <- epochs * n_batches
  
  # On enregistre les paramètres après chaque mise à jour
  w_hist <- numeric(n_steps + 1)
  b_hist <- numeric(n_steps + 1)
  epoch_hist <- integer(n_steps + 1)
  w_hist[1] <- w
  b_hist[1] <- b
  
  step <- 1
  for (epoch in 1:epochs) {
    idx <- sample(n)                                   # mélange des données
    batches <- split(idx, ceiling(seq_along(idx) / batch_size))
    
    for (batch in batches) {
      xb <- x[batch]
      yb <- y[batch]
      
      error <- yb - (b + w * xb)
      w_grad <- -2 * mean(error * xb)
      b_grad <- -2 * mean(error)
      
      w <- w - lr * w_grad
      b <- b - lr * b_grad
      
      step <- step + 1
      w_hist[step] <- w
      b_hist[step] <- b
      epoch_hist[step] <- epoch
    }
  }
  
  tibble(step = 0:n_steps, epoch = epoch_hist, w = w_hist, b = b_hist) |>
    # loss calculée sur TOUTES les données pour comparer les méthodes équitablement
    mutate(loss = map2_dbl(w, b, mse, x = x, y = y))
}

# ------------------------------------------------------------------
# 3. Configurations à comparer
# ------------------------------------------------------------------
w_init <- -1.5   # point de départ choisi loin de l'optimum (lisibilité du graphique)
b_init <- 1.5
epochs <- 30

configs <- tribble(
  ~method,                  ~batch_size, ~lr,
  "Batch (n = 342)",         n,          0.10,
  "Mini-batch (32)",         32,         0.05,
  "SGD (1 observation)",     1,          0.01
)

results <- pmap(configs, function(method, batch_size, lr) {
  gradient_descent(x, y, w_init, b_init, lr, epochs, batch_size) |>
    mutate(method = method)
}) |>
  bind_rows() |>
  mutate(method = factor(method, levels = configs$method))

# ------------------------------------------------------------------
# 4. Surface de loss L(w, b) : un paraboloïde (cuvette), pas une vraie cloche
# ------------------------------------------------------------------
w_seq <- seq(-2, 2, length.out = 80)
b_seq <- seq(-2, 2, length.out = 80)

loss_fun <- Vectorize(function(w, b) mse(w, b, x, y))
z <- outer(w_seq, b_seq, loss_fun)

grid <- expand_grid(w = w_seq, b = b_seq) |>
  mutate(loss = as.numeric(z))   # outer: w varie en premier -> ordre cohérent

# 4a. Surface 3D (base R) avec la trajectoire du batch GD
batch_path <- results |> filter(method == levels(method)[1])

pmat <- persp(
  w_seq, b_seq, z,
  theta = 35, phi = 30, expand = 0.6, col = "lightblue", border = NA,
  shade = 0.4, ticktype = "detailed",
  xlab = "w (pente)", ylab = "b (biais)", zlab = "MSE",
  main = "Surface de la loss L(w, b) et trajectoire (batch)"
)
lines(trans3d(batch_path$w, batch_path$b, batch_path$loss, pmat),
      col = "red", lwd = 2)
points(trans3d(w_init, b_init, batch_path$loss[1], pmat),
       pch = 19, col = "red", cex = 1.2)
points(trans3d(w_ols, b_ols, mse(w_ols, b_ols, x, y), pmat),
       pch = 8, col = "black", cex = 1.5)

# 4b. Courbes de niveau + trajectoires des trois méthodes
p_contour <- ggplot() +
  geom_contour_filled(data = grid, aes(w, b, z = loss), bins = 15, alpha = 0.8) +
  geom_path(data = results, aes(w, b), colour = "white", linewidth = 0.3) +
  geom_point(data = tibble(w = w_init, b = b_init), aes(w, b),
             colour = "red", size = 3) +
  geom_point(data = tibble(w = w_ols, b = b_ols), aes(w, b),
             shape = 8, colour = "black", size = 3, stroke = 1.2) +
  facet_wrap(~method) +
  labs(
    title = "Trajectoires dans l'espace des paramètres",
    subtitle = "Point rouge : initialisation | étoile : solution des moindres carrés",
    x = "w (pente)", y = "b (biais)", fill = "MSE"
  ) +
  theme_minimal()

print(p_contour)

# ------------------------------------------------------------------
# 5. Évolution de la loss : par epoch et par mise à jour des paramètres
# ------------------------------------------------------------------
loss_epoch <- results |>
  group_by(method, epoch) |>
  slice_tail(n = 1) |>      # loss en fin d'epoch (epoch 0 = initialisation)
  ungroup()

p_epoch <- ggplot(loss_epoch, aes(epoch, loss, colour = method)) +
  geom_line(linewidth = 1) +
  scale_y_log10() +
  labs(
    title = "Loss en fonction des epochs",
    subtitle = "Le batch fait 1 mise à jour/epoch, le SGD en fait 342",
    x = "Epoch", y = "MSE (échelle log)", colour = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

p_steps <- ggplot(results, aes(step + 1, loss, colour = method)) +
  geom_line(alpha = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Loss en fonction du nombre de mises à jour",
    subtitle = "Le SGD est bruité : la loss oscille autour du minimum",
    x = "Nombre de mises à jour (échelle log)", y = "MSE (échelle log)", colour = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_epoch)
print(p_steps)

# ------------------------------------------------------------------
# 6. Résultats finaux et retour à l'échelle d'origine
# ------------------------------------------------------------------
final <- results |>
  group_by(method) |>
  slice_tail(n = 1) |>
  ungroup() |>
  mutate(
    # y_raw = mean(y) + sd(y) * (b + w * (x_raw - mean(x)) / sd(x))
    w_orig = w * sd(y_raw) / sd(x_raw),
    b_orig = mean(y_raw) + sd(y_raw) * b - w_orig * mean(x_raw)
  ) |>
  select(method, w, b, loss, w_orig, b_orig)

print(final)

cat("\nSolution analytique (échelle d'origine) :\n")
print(coef(lm(body_mass_g ~ flipper_length_mm, data = df)))
