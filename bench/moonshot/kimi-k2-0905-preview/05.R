# Rayleigh–Bénard convection – Boussinesq approximation
# ----------------------------------------------------
# 2-D incompressible Navier–Stokes + energy equation
# Aspect ratio 2:1  |  Bottom T=1  |  Top T=0
# Output: interactive plotly overlay
#   - heat-map: temperature
#   - quiver: velocity vectors

library(tidyverse)
library(plotly)
library(pracma) # for fast 2-D finite differences

# -------------------- domain & numerics --------------------
nx <- 128
ny <- 64
Lx <- 2
Ly <- 1
dx <- Lx / (nx - 1)
dy <- Ly / (ny - 1)
Ra <- 1e5 # Rayleigh number
Pr <- 0.71 # Prandtl number
g <- 9.81
alpha <- 1 # thermal expansion (absorbed into Ra)
nu <- sqrt(Ra / Pr) * sqrt(g * alpha * Ly^3) / Ra # viscosity
kappa <- nu / Pr # thermal diffusivity
dt <- 0.2 * min(dx, dy)^2 / nu # viscous CFL
t_end <- 50 # run until rolls develop

# grid
x <- seq(0, Lx, length.out = nx)
y <- seq(0, Ly, length.out = ny)
xx <- outer(x, y, function(x, y) x)
yy <- outer(x, y, function(x, y) y)

# fields
u <- matrix(0, nx, ny) # x-velocity
v <- matrix(0, nx, ny) # y-velocity
T <- outer(x, y, function(x, y) 1 - y) # linear conduction profile
p <- matrix(0, nx, ny)

# -------------------- finite-difference helpers --------------------
lap <- function(f) del2(f, dx, dy)
grad <- function(f) {
  list(
    x = gradient(f, dx, dy)[[2]], # d/dx
    y = gradient(f, dx, dy)[[1]]
  )
} # d/dy
div <- function(u, v) grad(u)$x + grad(v)$y

# -------------------- projection method step --------------------
project <- function(u, v) {
  # compute divergence of tentative field
  div_uv <- div(u, v)
  # solve Poisson ∇²p = ∇·u / dt  (Neumann BC)
  p <<- poisson2D(div_uv / dx^2, dx, dy, maxit = 100, tol = 1e-4)
  # correct velocities
  grad_p <- grad(p)
  u <- u - dt * grad_p$x
  v <- v - dt * grad_p$y
  list(u = u, v = v)
}

# -------------------- boundary conditions --------------------
bc <- function() {
  # no-slip walls
  u[, 1] <- 0
  u[, ny] <- 0
  u[1, ] <- 0
  u[nx, ] <- 0
  v[, 1] <- 0
  v[, ny] <- 0
  v[1, ] <- 0
  v[nx, ] <- 0
  # temperature
  T[, 1] <- 1 # bottom hot
  T[, ny] <- 0 # top cold
  # insulating side walls (Neumann)
  T[1, ] <- T[2, ]
  T[nx, ] <- T[nx - 1, ]
}

# -------------------- time loop --------------------
t <- 0
while (t < t_end) {
  bc()

  # tentative momentum (buoyancy in y)
  lap_u <- lap(u)
  lap_v <- lap(v)
  buoy <- Ra / Pr * T # Boussinesq term Fy = alpha*g*T (scaled)

  u_star <- u + dt * (nu * lap_u - grad(u * u)$x - grad(u * v)$y)
  v_star <- v + dt * (nu * lap_v - grad(u * v)$x - grad(v * v)$y + buoy)

  # projection to enforce ∇·u = 0
  sol <- project(u_star, v_star)
  u <- sol$u
  v <- sol$v

  # temperature advection–diffusion
  lap_T <- lap(T)
  adv_T <- -grad(u * T)$x - grad(v * T)$y
  T <- T + dt * (kappa * lap_T + adv_T)

  t <- t + dt
}

# -------------------- tidy for plotly --------------------
df <- expand.grid(x = x, y = y) %>%
  mutate(T = c(T), u = c(u), v = c(v))

# down-sample arrows for clarity
arrow_df <- df %>%
  filter(x %% 4 == 0, y %% 4 == 0)

# -------------------- interactive overlay --------------------
p <- plot_ly() %>%
  add_trace(
    data = df,
    x = ~x,
    y = ~y,
    z = ~T,
    type = "heatmap",
    colors = "RdBu",
    zmin = 0,
    zmax = 1,
    colorbar = list(title = "T")
  ) %>%
  add_trace(
    data = arrow_df,
    x = ~x,
    y = ~y,
    u = ~u,
    v = ~v,
    type = "scatter",
    mode = "lines",
    line = list(color = "black", width = 1),
    hoverinfo = "skip"
  ) %>%
  layout(
    title = "Rayleigh–Bénard convection rolls",
    xaxis = list(title = "x"),
    yaxis = list(title = "y", scaleanchor = "x")
  )

p
# (Optional: export via plotly::export(p, "convection.html"))
