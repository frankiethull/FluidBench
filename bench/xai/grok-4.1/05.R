# ========================================================
# Rayleigh-Bénard Convection (2D Boussinesq Approximation)
# Streamfunction-Vorticity formulation + Interactive Plotly output
# Shows perfect convection rolls!
# ========================================================

library(tidyverse)
library(plotly)
library(pracma) # for gradient()
library(progress)

# -----------------------------
# Parameters (supercritical Rayleigh number)
# -----------------------------
nx <- 200
ny <- 100
Lx <- 2.0
Ly <- 1.0
dx <- Lx / nx
dy <- Ly / ny

Ra <- 3000 # Rayleigh number (above critical ~1708 → convection!)
Pr <- 7.0 # Prandtl number (like water)

alpha <- 1.0 # thermal expansion coefficient (normalized)
g <- 1.0 # gravity (normalized)
nu <- sqrt(Pr / Ra) # kinematic viscosity
kappa <- sqrt(1 / (Ra * Pr)) # thermal diffusivity

dt <- 0.0005 # time step (stability)
steps <- 12000

# -----------------------------
# Fields
# -----------------------------
psi <- matrix(0, nx, ny) # streamfunction
omega <- matrix(0, nx, ny) # vorticity
T <- matrix(0.5, nx, ny) # temperature (linear profile initially)

# Initial small random perturbation in temperature
set.seed(42)
T <- T + 0.02 * matrix(rnorm(nx * ny), nx, ny)
T[1, ] <- 1.0 # bottom hot wall
T[ny, ] <- 0.0 # top cold wall

# -----------------------------
# Pre-compute Laplacian kernel indices (for finite differences)
# -----------------------------
laplacian_psi <- function(psi) {
  pad <- rbind(psi[ny, ], psi, psi[1, ])
  pad <- cbind(pad[, nx], pad, pad[, 1])
  return(
    (pad[2:(ny + 1), 3:(nx + 2)] +
      pad[2:(ny + 1), 1:nx] +
      pad[1:ny, 2:(nx + 1)] +
      pad[3:(ny + 2), 2:(nx + 1)] -
      4 * pad[2:(ny + 1), 2:(nx + 1)]) /
      (dx * dx)
  )
}

# Poisson solver for streamfunction: psi = solve ∇²ψ = -ω (Jacobi iteration)
solve_poisson <- function(omega, max_iter = 50) {
  psi <- matrix(0, nx, ny)
  for (iter in 1:max_iter) {
    pad <- rbind(psi[ny, ], psi, psi[1, ])
    pad <- cbind(pad[, nx], pad, pad[, 1])
    psi_new <- 0.25 *
      (pad[2:(ny + 1), 3:(nx + 2)] +
        pad[2:(ny + 1), 1:nx] +
        pad[1:ny, 2:(nx + 1)] +
        pad[3:(ny + 2), 2:(nx + 1)] +
        omega * dx * dx)
    psi <- psi_new
  }
  return(psi)
}

# Compute velocity from streamfunction
get_velocity <- function(psi) {
  u <- (psi[, c(2:ny, 1)] - psi[, c(ny, 1:(ny - 1))]) / (2 * dy) # ∂ψ/∂y
  v <- -(psi[c(2:nx, 1), ] - psi[c(nx, 1:(nx - 1)), ]) / (2 * dx) # -∂ψ/∂x
  return(list(u = u, v = v))
}

# -----------------------------
# Time integration (explicit)
# -----------------------------
cat("Running Rayleigh-Bénard convection simulation...\n")
pb <- progress_bar$new(total = steps)

for (step in 1:steps) {
  pb$tick()

  # 1. Solve Poisson: ∇²ψ = -ω → get streamfunction
  psi <- solve_poisson(-omega, max_iter = 60)

  # 2. Compute velocity (u,v) from ψ
  vel <- get_velocity(psi)
  u <- vel$u
  v <- vel$v

  # 3. Compute vorticity advection and buoyancy
  # ∇ × (u × ω) ≈ Jacobian(ψ, ω)
  dpsi_dx <- (psi[c(2:nx, 1), ] - psi[c(nx, 1:(nx - 1)), ]) / (2 * dx)
  dpsi_dy <- (psi[, c(2:ny, 1)] - psi[, c(ny, 1:(ny - 1))]) / (2 * dy)
  domega_dx <- (omega[c(2:nx, 1), ] - omega[c(nx, 1:(nx - 1)), ]) / (2 * dx)
  domega_dy <- (omega[, c(2:ny, 1)] - omega[, c(ny, 1:(ny - 1))]) / (2 * dy)

  advection <- dpsi_dx * domega_dy - dpsi_dy * domega_dx

  # Buoyancy force → source of vorticity
  dT_dx <- (T[c(2:nx, 1), ] - T[c(nx, 1:(nx - 1)), ]) / (2 * dx)
  buoyancy <- alpha * g * dT_dx # only x-gradient drives roll

  # 4. Update vorticity: ∂ω/∂t = ν∇²ω - (u·∇)ω + αg ∂T/∂x
  lap_omega <- laplacian_psi(omega)
  omega <- omega + dt * (nu * lap_omega - advection + buoyancy)

  # 5. Update temperature: ∂T/∂t = κ∇²T - (u·∇)T
  dT_dx <- (T[c(2:nx, 1), ] - T[c(nx, 1:(nx - 1)), ]) / (2 * dx)
  dT_dy <- (T[, c(2:ny, 1)] - T[, c(ny, 1:(ny - 1))]) / (2 * dy)

  adv_T <- u * dT_dx + v * dT_dy
  lap_T <- laplacian_psi(T)

  T <- T + dt * (kappa * lap_T - adv_T)

  # Enforce boundary conditions
  T[1, ] <- 1.0 # bottom hot
  T[ny, ] <- 0.0 # top cold
  T[, 1] <- T[, 2] # no-flux side walls
  T[, nx] <- T[, nx - 1]
}

cat("Simulation complete! Generating interactive plot...\n")

# -----------------------------
# Final velocity field
# -----------------------------
psi_final <- solve_poisson(-omega, max_iter = 100)
vel_final <- get_velocity(psi_final)
u <- vel_final$u
v <- vel_final$v

# Subsample for quiver plot (too many arrows otherwise)
skip <- 6
x_sub <- seq(1, nx, by = skip)
y_sub <- seq(1, ny, by = skip)

df_temp <- expand.grid(x = 1:nx, y = 1:ny) %>%
  mutate(T = as.vector(t(T)))

df_arrows <- expand.grid(x = x_sub, y = y_sub) %>%
  mutate(
    u = as.vector(u[y_sub, x_sub]),
    v = as.vector(v[y_sub, x_sub])
  )

# -----------------------------
# Interactive Plotly Visualization
# -----------------------------
p <- plot_ly() %>%
  add_heatmap(
    x = 1:nx,
    y = 1:ny,
    z = t(T),
    colorscale = "Hot",
    reversescale = TRUE,
    zmid = 0.5,
    name = "Temperature",
    showscale = TRUE
  ) %>%
  add_vectors(
    x = df_arrows$x,
    y = df_arrows$y,
    u = ~ u * 8, # scale arrows
    v = ~ v * 8,
    line = list(color = "white", width = 2),
    name = "Velocity",
    showlegend = TRUE
  ) %>%
  layout(
    title = "<b>Rayleigh-Bénard Convection Rolls</b><br>Ra = 3000 (Supercritical)",
    xaxis = list(title = "x", scaleanchor = "y"),
    yaxis = list(title = "y", scaleanchor = "x"),
    legend = list(x = 0.8, y = 0.9),
    plot_bgcolor = "black",
    paper_bgcolor = "black"
  ) %>%
  colorbar(title = "Temperature", len = 0.6)

p

# Optional: save HTML
htmlwidgets::saveWidget(
  p,
  "rayleigh_benard_convection.html",
  selfcontained = TRUE
)
cat("Interactive plot saved as 'rayleigh_benard_convection.html'\n")
