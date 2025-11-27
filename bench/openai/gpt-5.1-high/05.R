############################################################
# Rayleigh–Bénard Convection (2D, Boussinesq approximation)
# - 2:1 aspect ratio domain
# - Heated from below (T = 1), cooled from above (T = 0)
# - Streamfunction–vorticity formulation (incompressible)
# - Plot: Temperature heatmap + velocity vector field
############################################################

# You will need:
# install.packages(c("ggplot2", "viridis"))

library(ggplot2)
library(viridis)

####------------------------------------------------------####
#### 1. Parameters and grid                              ####
####------------------------------------------------------####

# Dimensionless parameters
Pr <- 1.0 # Prandtl number
Ra <- 1e4 # Rayleigh number (>= 1708 for convection onset)

# Grid (aspect ratio 2:1)
Nx <- 64 # number of points in x
Ny <- 32 # number of points in y
Lx <- 2.0 # domain length in x
Ly <- 1.0 # domain length in y

dx <- Lx / (Nx - 1)
dy <- Ly / (Ny - 1)

dt <- 1e-3 # time step (explicit scheme – small for stability)
nsteps <- 1500 # number of time steps (increase for more developed rolls)
psi_iter <- 60 # Poisson iterations per step (streamfunction solve)

dx2 <- dx^2
dy2 <- dy^2
inv_denom <- 1 / (2 * (1 / dx2 + 1 / dy2)) # for Poisson Gauss–Seidel

# Coordinates (for plotting)
x <- seq(0, Lx, length.out = Nx)
y <- seq(0, Ly, length.out = Ny)

####------------------------------------------------------####
#### 2. Fields and initial conditions                     ####
####------------------------------------------------------####

# Temperature T, vorticity omega, streamfunction psi
T <- matrix(0, Nx, Ny)
omega <- matrix(0, Nx, Ny)
psi <- matrix(0, Nx, Ny)

# Initial temperature: linear conduction profile + small noise
for (j in 1:Ny) {
  yj <- (j - 1) * dy / Ly
  T[, j] <- 1 - yj # 1 at bottom, 0 at top (linear)
}
set.seed(1)
T[, 2:(Ny - 1)] <- T[, 2:(Ny - 1)] +
  0.01 * matrix(rnorm(Nx * (Ny - 2)), Nx, Ny - 2)

# Enforce boundary temperatures
T[, 1] <- 1.0 # bottom hot
T[, Ny] <- 0.0 # top cold

# Velocities (computed from streamfunction later)
u <- matrix(0, Nx, Ny)
v <- matrix(0, Nx, Ny)

####------------------------------------------------------####
#### 3. Helper: Poisson solver ∇²ψ = -ω (psi from omega)  ####
####    - periodic in x, ψ = 0 at top/bottom (no-slip)   ####
####------------------------------------------------------####

solve_streamfunction <- function(psi, omega, niter, dx2, dy2, inv_denom) {
  Nx <- nrow(psi)
  Ny <- ncol(psi)

  for (k in 1:niter) {
    # Copy previous psi to update interior
    psi_old <- psi

    # Interior points 2..(Ny-1) in y, all x (1..Nx) with periodic x
    for (j in 2:(Ny - 1)) {
      for (i in 1:Nx) {
        iL <- if (i == 1) Nx else i - 1
        iR <- if (i == Nx) 1 else i + 1

        psi[i, j] <- ((psi_old[iL, j] + psi_old[iR, j]) /
          dx2 +
          (psi_old[i, j - 1] + psi_old[i, j + 1]) / dy2 +
          omega[i, j]) *
          inv_denom
      }
    }

    # Boundary conditions for psi: no-slip walls at y=0 and y=1 => ψ = 0
    psi[, 1] <- 0
    psi[, Ny] <- 0
    # Periodic in x is naturally enforced via iL, iR wrap
  }

  psi
}

####------------------------------------------------------####
#### 4. Helper: velocities from streamfunction            ####
####    u = ∂ψ/∂y, v = -∂ψ/∂x                             ####
####------------------------------------------------------####

compute_velocity <- function(psi, dx, dy) {
  Nx <- nrow(psi)
  Ny <- ncol(psi)
  u <- matrix(0, Nx, Ny)
  v <- matrix(0, Nx, Ny)

  for (j in 2:(Ny - 1)) {
    for (i in 1:Nx) {
      iL <- if (i == 1) Nx else i - 1
      iR <- if (i == Nx) 1 else i + 1
      jB <- j - 1
      jT <- j + 1

      u[i, j] <- (psi[i, jT] - psi[i, jB]) / (2 * dy) # u = dψ/dy
      v[i, j] <- -(psi[iR, j] - psi[iL, j]) / (2 * dx) # v = -dψ/dx
    }
  }

  # No-slip at top/bottom
  u[, 1] <- 0
  u[, Ny] <- 0
  v[, 1] <- 0
  v[, Ny] <- 0

  list(u = u, v = v)
}

####------------------------------------------------------####
#### 5. Time stepping                                      ####
####    Vorticity & temperature equations:                 ####
####    dω/dt + u·∇ω = Pr ∇²ω + Ra·Pr ∂T/∂x                ####
####    dT/dt + u·∇T = ∇²T                                 ####
####------------------------------------------------------####

for (n in 1:nsteps) {
  # 1) Solve for streamfunction from vorticity
  psi <- solve_streamfunction(psi, omega, psi_iter, dx2, dy2, inv_denom)

  # 2) Velocities from streamfunction
  vel <- compute_velocity(psi, dx, dy)
  u <- vel$u
  v <- vel$v

  # 3) Update temperature and vorticity (explicit Euler, upwind advection)
  T_new <- T
  omega_new <- omega

  for (j in 2:(Ny - 1)) {
    for (i in 1:Nx) {
      iL <- if (i == 1) Nx else i - 1
      iR <- if (i == Nx) 1 else i + 1
      jB <- j - 1
      jT <- j + 1

      # ----- Temperature -----
      # Laplacian (central)
      lapT <- (T[iL, j] - 2 * T[i, j] + T[iR, j]) /
        dx2 +
        (T[i, jB] - 2 * T[i, j] + T[i, jT]) / dy2

      # Upwind derivatives for advection
      dTdx <- if (u[i, j] > 0) {
        (T[i, j] - T[iL, j]) / dx
      } else {
        (T[iR, j] - T[i, j]) / dx
      }

      dTdy <- if (v[i, j] > 0) {
        (T[i, j] - T[i, jB]) / dy
      } else {
        (T[i, jT] - T[i, j]) / dy
      }

      advT <- u[i, j] * dTdx + v[i, j] * dTdy

      T_new[i, j] <- T[i, j] + dt * (lapT - advT)

      # ----- Vorticity -----
      # Laplacian (central)
      lapW <- (omega[iL, j] - 2 * omega[i, j] + omega[iR, j]) /
        dx2 +
        (omega[i, jB] - 2 * omega[i, j] + omega[i, jT]) / dy2

      # Upwind derivatives for advection
      dWdx <- if (u[i, j] > 0) {
        (omega[i, j] - omega[iL, j]) / dx
      } else {
        (omega[iR, j] - omega[i, j]) / dx
      }

      dWdy <- if (v[i, j] > 0) {
        (omega[i, j] - omega[i, jB]) / dy
      } else {
        (omega[i, jT] - omega[i, j]) / dy
      }

      advW <- u[i, j] * dWdx + v[i, j] * dWdy

      # Buoyancy contribution Ra*Pr*∂T/∂x (central in x)
      dTdx_c <- (T[iR, j] - T[iL, j]) / (2 * dx)

      omega_new[i, j] <- omega[i, j] +
        dt * (Pr * lapW - advW + Ra * Pr * dTdx_c)
    }
  }

  # 4) Enforce boundary conditions

  # Temperature: fixed top/bottom, periodic in x automatically via stencil
  T_new[, 1] <- 1.0
  T_new[, Ny] <- 0.0

  # Streamfunction walls already enforced in solver; now impose wall vorticity
  # No-slip walls: ω_wall ≈ -2*(ψ_interior - ψ_wall)/dy², ψ_wall = 0
  for (i in 1:Nx) {
    omega_new[i, 1] <- -2 * psi[i, 2] / dy2 # bottom
    omega_new[i, Ny] <- -2 * psi[i, Ny - 1] / dy2 # top
  }

  # Update fields
  T <- T_new
  omega <- omega_new

  if (n %% 100 == 0) {
    cat("Step", n, "of", nsteps, "\n")
  }
}

####------------------------------------------------------####
#### 6. Final velocity (solve psi once more)              ####
####------------------------------------------------------####

psi <- solve_streamfunction(psi, omega, psi_iter * 2, dx2, dy2, inv_denom)
vel <- compute_velocity(psi, dx, dy)
u <- vel$u
v <- vel$v

####------------------------------------------------------####
#### 7. Build data for ggplot                             ####
####------------------------------------------------------####

# Temperature field (for heatmap)
df_T <- expand.grid(
  ix = 1:Nx,
  jy = 1:Ny
)
df_T$x <- x[df_T$ix]
df_T$y <- y[df_T$jy]
df_T$T <- as.vector(T)

# Velocity arrows: thin grid to avoid clutter
skip <- 2
ix_vec <- seq(1, Nx, by = skip)
jy_vec <- seq(1, Ny, by = skip)

arrow_grid <- expand.grid(
  ix = ix_vec,
  jy = jy_vec
)

arrow_grid$x <- x[arrow_grid$ix]
arrow_grid$y <- y[arrow_grid$jy]
arrow_grid$u <- u[cbind(arrow_grid$ix, arrow_grid$jy)]
arrow_grid$v <- v[cbind(arrow_grid$ix, arrow_grid$jy)]

# Scale arrows visually
arrow_scale <- 0.05
arrow_grid$xend <- arrow_grid$x + arrow_scale * arrow_grid$u
arrow_grid$yend <- arrow_grid$y + arrow_scale * arrow_grid$v

####------------------------------------------------------####
#### 8. Plot: Temperature heatmap + velocity vectors      ####
####------------------------------------------------------####

ggplot() +
  geom_raster(data = df_T, aes(x = x, y = y, fill = T)) +
  scale_fill_viridis_c(
    option = "inferno",
    direction = -1,
    name = "Temperature"
  ) +
  geom_segment(
    data = arrow_grid,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(type = "closed", length = unit(0.12, "cm")),
    color = "white",
    size = 0.3,
    alpha = 0.8
  ) +
  coord_fixed(ratio = 1) +
  labs(
    title = "Rayleigh–Bénard Convection (Ra = 1e4, Pr = 1)",
    x = "x",
    y = "y"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black")
  )
