# Régression linéaire simple

# Packages ----

library(tidyverse)
library(palmerpenguins)
library(skimr)
library(GGally)

theme_set(theme_bw())

my_cols <- c(
  "Adelie" = "darkorange",
  "Chinstrap" = "purple",
  "Gentoo" = "cyan4"
)

# Data ----

penguins <- palmerpenguins::penguins

# Exploratory Data Analysis ----

head(penguins)
glimpse(penguins)
skim(penguins)

# Clean data ----

penguins <- penguins |> 
  mutate(year = factor(year)) |> 
  drop_na()

penguins |> 
  select(where(is.numeric)) |> 
  ggpairs()

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
