# 1-D Linear Advection – “Space-Time Waterfall”
# ---------------------------------------------
# PDE:  ∂u/∂t + c ∂u/∂x = 0   with periodic BC
# Domain: x ∈ [0,2]   c = 1
# Output: ggplot2 heat-map (x-axis = space, y-axis = time, colour = u)

library(tidyverse)

# ----------------- user controls -----------------
Nx <- 200 # grid points
x_max <- 2 # domain length
c <- 1 # advection speed
CFL <- 0.9 # Courant number
t_max <- 2 # final physical time
# -------------------------------------------------

dx <- x_max / Nx
dt <- CFL * dx / c
Nt <- ceiling(t_max / dt) # number of time steps
dt <- t_max / Nt # exact dt that divides t_max

x <- seq(0, x_max, length.out = Nx + 1)[-1] # drop 0 to have periodic x_max

# initial condition: Gaussian pulse
u0 <- exp(-(x - 0.5)^2 / 0.05)

# storage matrix: each row = one time level
U <- matrix(0, nrow = Nt + 1, ncol = Nx)
U[1, ] <- u0

# ----------------- solver (1st-order upwind) -----------------
u <- u0
for (n in 1:Nt) {
  u <- u - c * dt / dx * (u - c(c(u[-1], u[1]))) # periodic wrap-around
  U[n + 1, ] <- u
}

# ----------------- reshape for ggplot -----------------
space_time <- expand.grid(x = x, t = seq(0, t_max, length.out = Nt + 1)) %>%
  mutate(u = c(U))

# ----------------- plot -----------------
ggplot(space_time, aes(x, t, fill = u)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C") +
  labs(
    x = "Space  x",
    y = "Time  t",
    fill = "Amplitude  u",
    title = "1-D Advection: Space-Time Waterfall"
  ) +
  theme_minimal(base_size = 14) +
  coord_fixed(ratio = 1) # keeps the 45° diagonal stripes
