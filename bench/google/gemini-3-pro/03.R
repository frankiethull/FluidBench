# Load libraries
library(ggplot2)
library(reshape2) # For melting matrix to data frame

# ==========================================
# 1. PARAMETERS
# ==========================================

# Grid size
N <- 100

# Physics Parameters ("Coral" preset)
Du <- 0.16 # Diffusion rate of U
Dv <- 0.08 # Diffusion rate of V
F <- 0.035 # Feed rate
k <- 0.060 # Kill rate

# Simulation settings
dt <- 1.0 # Time step
steps <- 2500 # Number of iterations (2000-3000 yields good patterns)

# ==========================================
# 2. INITIALIZATION
# ==========================================

# Create grids
# U is completely filled with 1.0 initially
U <- matrix(1.0, nrow = N, ncol = N)

# V is completely empty (0.0) initially
V <- matrix(0.0, nrow = N, ncol = N)

# Seeding the reaction:
# Create a 10x10 square in the center where we disturb the equilibrium
# We lower U and raise V to start the autocatalytic process
mid <- N / 2
range <- (mid - 5):(mid + 5)

U[range, range] <- 0.50
V[range, range] <- 0.25

# Add a tiny bit of random noise to the square to break perfect symmetry
# allowing more organic growth
noise <- matrix(runif(length(range)^2, 0, 0.02), nrow = length(range))
V[range, range] <- V[range, range] + noise

# ==========================================
# 3. SIMULATION LOOP (Vectorized)
# ==========================================
# We use integer vectors for periodic boundary shifting
# This avoids creating new index vectors inside the loop for speed
idx_curr <- 1:N
idx_up <- c(2:N, 1) # Shift Up: Row 2 becomes 1, Row 1 wraps to N
idx_down <- c(N, 1:(N - 1)) # Shift Down: Row N becomes 1
idx_left <- c(2:N, 1) # Shift Left
idx_right <- c(N, 1:(N - 1)) # Shift Right

cat("Running simulation... Please wait.\n")

for (t in 1:steps) {
  # --- 1. CALCULATE LAPLACIAN (Diffusion) ---
  # Laplacian = Neighbor_Sum - 4 * Center
  # We calculate neighbors by shifting the matrices

  # Laplacian of U
  Lu <- (U[idx_up, ] + U[idx_down, ] + U[, idx_left] + U[, idx_right]) - 4 * U

  # Laplacian of V
  Lv <- (V[idx_up, ] + V[idx_down, ] + V[, idx_left] + V[, idx_right]) - 4 * V

  # --- 2. REACTION TERM ---
  # Reaction: uv^2
  uvv <- U * V * V

  # --- 3. UPDATE STATE (Euler Integration) ---
  # dU/dt = Du*LapU - uv^2 + F*(1-U)
  U <- U + (Du * Lu - uvv + F * (1 - U)) * dt

  # dV/dt = Dv*LapV + uvv - (F+k)*V
  V <- V + (Dv * Lv + uvv - (F + k) * V) * dt

  # Progress bar every 500 steps
  if (t %% 500 == 0) cat(sprintf("Step %d / %d\n", t, steps))
}

# ==========================================
# 4. VISUALIZATION
# ==========================================

# Reshape V-matrix for ggplot
# melt converts the matrix into: Var1 (Row), Var2 (Col), value (Concentration)
df_v <- melt(V)
colnames(df_v) <- c("x", "y", "concentration")

# Plot
ggplot(df_v, aes(x = x, y = y, fill = concentration)) +
  geom_raster(interpolate = TRUE) + # Interpolate smooths the pixels slightly
  scale_fill_viridis_c(option = "magma") + # "magma" gives a nice coral/biological look
  coord_fixed() + # Ensures the grid remains square
  labs(
    title = "Gray-Scott Reaction-Diffusion",
    subtitle = paste0("F=", F, ", k=", k, " (Coral Pattern)"),
    x = NULL,
    y = NULL,
    fill = "[V]"
  ) +
  theme_void() + # Removes axes and background for a clean art look
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      margin = margin(b = 10)
    ),
    plot.subtitle = element_text(hjust = 0.5)
  )
