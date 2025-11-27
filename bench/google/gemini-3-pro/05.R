# Load necessary libraries
library(ggplot2)
library(reshape2) # For data wrangling
library(scales) # For vector scaling

# ==========================================
# 1. CONFIGURATION & PHYSICS
# ==========================================

# Domain
Nx <- 100 # Grid width
Ny <- 50 # Grid height (Aspect ratio 2:1)
Lx <- 2.0
Ly <- 1.0
dx <- Lx / Nx
dy <- Ly / Ny

# Physical Parameters (Rayleigh-Bénard)
# Pr (Prandtl) = Viscosity / Thermal Diffusivity
# Ra (Rayleigh) = Drives the convection strength
Pr <- 10.0
Ra <- 20000.0 # High enough to cause instability (Critical Ra ~ 1708)

# Time stepping
dt <- 0.0001 # Very small time step required for stability in explicit methods
steps <- 4000 # Number of iterations

# ==========================================
# 2. INITIALIZATION
# ==========================================

# Allocate matrices
T_mat <- matrix(0, Nx, Ny) # Temperature
W <- matrix(0, Nx, Ny) # Vorticity (Omega)
P <- matrix(0, Nx, Ny) # Streamfunction (Psi)

# Initial Condition: Linear Temperature gradient + Noise
# T=1 at bottom (y=1), T=0 at top (y=Ny)
for (j in 1:Ny) {
  T_mat[, j] <- 1 - (j - 1) / (Ny - 1)
}

# Add random noise to kickstart the instability
set.seed(123)
noise <- matrix(runif(Nx * Ny, -0.05, 0.05), Nx, Ny)
T_mat <- T_mat + noise

# Enforce Boundary Conditions on T immediately
T_mat[, 1] <- 1.0 # Bottom Hot
T_mat[, Ny] <- 0.0 # Top Cold

# ==========================================
# 3. VECTORIZED HELPER FUNCTIONS
# ==========================================
# Instead of loops, we shift matrices to access neighbors
# Periodic boundaries on X, Rigid (No-Slip) on Y

shift_L <- function(m) cbind(m[, 2:ncol(m)], m[, ncol(m)]) # Neumann approx at walls for internal calc
shift_R <- function(m) cbind(m[, 1], m[, 1:(ncol(m) - 1)])
shift_U <- function(m) rbind(m[2:nrow(m), ], m[1, ]) # Periodic X
shift_D <- function(m) rbind(m[nrow(m), ], m[1:(nrow(m) - 1)]) # Periodic X

# Laplacian (Finite Difference)
laplacian <- function(M) {
  (shift_U(M) + shift_D(M) - 2 * M) /
    (dx^2) +
    (shift_L(M) + shift_R(M) - 2 * M) / (dy^2)
}

# ==========================================
# 4. SIMULATION LOOP
# ==========================================
cat(
  "Simulating Convection... (This uses an iterative Poisson solver, please wait)\n"
)

for (step in 1:steps) {
  # --- A. SOLVE STREAMFUNCTION (Poisson Eq: Lap(Psi) = -Omega) ---
  # We use Successive Over-Relaxation (SOR) or simple Jacobi iteration.
  # Since the fluid changes slowly, we iterate P a few times per time step.
  for (iter in 1:20) {
    P_new <- 0.25 *
      (shift_U(P) + shift_D(P) + shift_L(P) + shift_R(P) + dx * dy * W)
    # Enforce Boundary Conditions on Psi (Rigid Walls at Top/Bottom: Psi=0)
    P_new[, 1] <- 0
    P_new[, Ny] <- 0
    P <- P_new
  }

  # --- B. CALCULATE VELOCITY (u = dP/dy, v = -dP/dx) ---
  # Central difference
  u <- (shift_L(P) - shift_R(P)) / (2 * dy)
  v <- -(shift_U(P) - shift_D(P)) / (2 * dx)

  # --- C. CALCULATE DERIVATIVES FOR TRANSPORT ---
  # 1. Temperature Diffusion (Laplacian)
  lapT <- laplacian(T_mat)

  # 2. Vorticity Diffusion (Laplacian)
  lapW <- laplacian(W)

  # 3. Buoyancy Term (Ra * Pr * dT/dx)
  # Gravity acts in Y, so horizontal T gradients cause rotation
  dT_dx <- (shift_U(T_mat) - shift_D(T_mat)) / (2 * dx)
  buoyancy <- Ra * Pr * dT_dx

  # 4. Advection (u * d/dx + v * d/dy)
  # Using simple central difference for gradients (Note: Upwind is more stable but complex)
  # dT
  dT_dx_adv <- (shift_U(T_mat) - shift_D(T_mat)) / (2 * dx)
  dT_dy_adv <- (shift_L(T_mat) - shift_R(T_mat)) / (2 * dy)
  adv_T <- u * dT_dx_adv + v * dT_dy_adv

  # dW
  dW_dx_adv <- (shift_U(W) - shift_D(W)) / (2 * dx)
  dW_dy_adv <- (shift_L(W) - shift_R(W)) / (2 * dy)
  adv_W <- u * dW_dx_adv + v * dW_dy_adv

  # --- D. TIME INTEGRATION (Euler) ---
  # dT/dt = Laplacian(T) - Advection(T)
  T_mat <- T_mat + dt * (lapT - adv_T)

  # dW/dt = Pr * Laplacian(W) - Advection(W) + Buoyancy
  W <- W + dt * (Pr * lapW - adv_W + buoyancy)

  # --- E. BOUNDARY CONDITIONS ---
  # Temperature: Fixed Top/Bottom
  T_mat[, 1] <- 1.0
  T_mat[, Ny] <- 0.0

  # Vorticity: No-Slip Walls (Thom's Formula)
  # W_wall = -2 * Psi_neighbor / dy^2
  W[, 1] <- -2 * P[, 2] / dy^2 # Bottom Wall
  W[, Ny] <- -2 * P[, Ny - 1] / dy^2 # Top Wall

  if (step %% 500 == 0) cat(sprintf("Step: %d\n", step))
}

# ==========================================
# 5. DATA PREPARATION FOR PLOTTING
# ==========================================

# 1. Prepare Temperature Heatmap Data
df_temp <- melt(T_mat)
colnames(df_temp) <- c("x_idx", "y_idx", "Temp")
df_temp$x <- (df_temp$x_idx - 1) * dx
df_temp$y <- (df_temp$y_idx - 1) * dy

# 2. Prepare Vector Field Data (Downsampling)
# We can't plot 5000 arrows, it will be a mess. Let's take every 4th point.
skip <- 4
ix <- seq(1, Nx, by = skip)
iy <- seq(1, Ny, by = skip)

# Extract u and v at these indices
u_sub <- u[ix, iy]
v_sub <- v[ix, iy]
x_sub <- (ix - 1) * dx
y_sub <- (iy - 1) * dy

df_vec <- data.frame(
  x = as.vector(outer(x_sub, rep(1, length(iy)))),
  y = as.vector(outer(rep(1, length(ix)), y_sub)),
  u = as.vector(u_sub),
  v = as.vector(v_sub)
)

# Normalize vectors for visualization (so arrows aren't huge)
df_vec$len <- sqrt(df_vec$u^2 + df_vec$v^2)
scale_factor <- 0.05 / max(df_vec$len)
df_vec$u_plot <- df_vec$u * scale_factor
df_vec$v_plot <- df_vec$v * scale_factor


# ==========================================
# 6. VISUALIZATION
# ==========================================

ggplot() +
  # A. The Temperature Heatmap
  geom_raster(
    data = df_temp,
    aes(x = x, y = y, fill = Temp),
    interpolate = TRUE
  ) +
  scale_fill_distiller(palette = "RdYlBu", name = "Temp") +

  # B. The Vector Field (Arrows)
  geom_segment(
    data = df_vec,
    aes(x = x, y = y, xend = x + u_plot, yend = y + v_plot),
    arrow = arrow(length = unit(0.15, "cm")),
    color = "black",
    alpha = 0.6,
    size = 0.4
  ) +

  # C. Formatting
  coord_fixed(ratio = 1) +
  labs(
    title = "Rayleigh-Bénard Convection",
    subtitle = "Hot Bottom / Cold Top (Ra=20,000). Color=Temp, Arrows=Velocity",
    x = "X Position",
    y = "Y Position"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
