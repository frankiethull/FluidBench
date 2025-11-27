# Gray-Scott Reaction-Diffusion – “Turing Pattern”
# -------------------------------------------------
# Grid 100×100  |  Parameters: Du = 0.16, Dv = 0.08, F = 0.035, k = 0.060
# 2 000 iterations  |  Output: hi-res ggplot2 heat-map of chemical V

library(tidyverse)

# -------------------- parameters --------------------
N <- 100
Du <- 0.16
Dv <- 0.08
F <- 0.035
k <- 0.060
dt <- 1.0 # stable for these parameters
steps <- 2000

# -------------------- initial state --------------------
U <- matrix(1, N, N)
V <- matrix(0, N, N)

# central 10×10 noise patch
mid <- N %/% 2
V[(mid - 4):(mid + 5), (mid - 4):(mid + 5)] <- runif(100)

# -------------------- Laplacian (5-point stencil) --------------------
lap <- function(x) {
  rbind(x[N, ], x[-N, ]) +
    rbind(x[-1, ], x[1, ]) +
    cbind(x[, N], x[, -N]) +
    cbind(x[, -1], x[, 1]) -
    4 * x
}

# -------------------- time loop --------------------
for (i in 1:steps) {
  Lu <- lap(U)
  Lv <- lap(V)

  uv2 <- U * V^2

  U <- U + dt * (Du * Lu - uv2 + F * (1 - U))
  V <- V + dt * (Dv * Lv + uv2 - (F + k) * V)
}

# -------------------- tidy for ggplot --------------------
df <- expand.grid(x = 1:N, y = 1:N) %>%
  mutate(V = c(V))

# -------------------- plot --------------------
ggplot(df, aes(x, y, fill = V)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", begin = 0, end = 1) +
  coord_fixed(expand = FALSE) +
  theme_void() +
  labs(
    title = "Gray-Scott Turing Pattern  (Coral preset)",
    subtitle = "Chemical V after 2 000 iterations",
    fill = "V"
  ) +
  theme(legend.position = "right")
