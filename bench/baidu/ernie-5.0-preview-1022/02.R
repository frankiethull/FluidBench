# --------------------------------------------------------------------------
# Problem 2: Von Kármán Vortex Street (Lattice Boltzmann)
# --------------------------------------------------------------------------
# Load necessary libraries
library(ggplot2)
library(gganimate)
library(dplyr)
library(tidyr)
library(png) # For writing the GIF

# 1. SETUP: Define simulation parameters and LBM constants
# ---------------------------------------------------------
# Grid parameters
nx <- 80
ny <- 40
n_steps <- 2000 # Number of time steps to simulate

# Fluid and LBM parameters
tau <- 0.6 # Relaxation time (must be > 0.5)
omega <- 1 / tau # Relaxation frequency
rho0 <- 1.0 # Reference density

# Obstacle parameters
obs_x <- 20
obs_y <- 20
obs_r <- 3

# Inlet/Outlet parameters
u_in <- 0.05 # Inlet velocity (low for stability)

# D2Q9 Lattice weights and velocity vectors
# e_x, e_y are the discrete velocity vectors for each direction 'i'
weights <- c(4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36)
ex <- c(0, 1, 0, -1, 0, 1, -1, -1, 1)
ey <- c(0, 0, 1, 0, -1, 1, 1, -1, -1)

# Pre-calculate terms for equilibrium distribution
# This is a common optimization
e_sq <- ex^2 + ey^2
u_eq_term1 <- 3 * (ex + ey)
u_eq_term2 <- 9 / 2 * (ex^2 + ey^2)
u_eq_term3 <- 3 / 2 * (ex + ey)

# 2. INITIALIZATION: Create grids and set initial state
# -----------------------------------------------------
# Create the distribution function array: f[i, x, y]
# We use a 3D array: 9 directions x nx columns x ny rows
f <- array(0, dim = c(9, nx, ny))
f_eq <- array(0, dim = c(9, nx, ny))

# Initialize with equilibrium distribution at zero velocity
for (i in 1:9) {
  f_eq[i, , ] <- weights[i] * rho0
}
f <- f_eq

# Create obstacle mask (1 where obstacle is, 0 otherwise)
obstacle <- array(0, dim = c(nx, ny))
for (x in 1:nx) {
  for (y in 1:ny) {
    if ((x - obs_x)^2 + (y - obs_y)^2 < obs_r^2) {
      obstacle[x, y] <- 1
    }
  }
}

# Storage for animation frames (store every 5th step)
frame_interval <- 5
vorticity_history <- list()
frame_count <- 0

# 3. SIMULATION LOOP
# ------------------
cat("Starting LBM simulation...\n")
for (step in 1:n_steps) {
  # --- MACROSCOPIC PROPERTIES ---
  # Calculate density (rho) and velocity (u, v) from f
  rho <- apply(f, c(2, 3), sum)

  # Momentum calculation
  momentum_x <- apply(f * ex, c(2, 3), sum)
  momentum_y <- apply(f * ey, c(2, 3), sum)

  # Velocity calculation (with zero-velocity inside obstacle)
  u <- momentum_x / rho
  v <- momentum_y / rho
  u[obstacle == 1] <- 0
  v[obstacle == 1] <- 0

  # --- BOUNDARY CONDITIONS ---
  # Inlet (Left): Zou/He BC for constant velocity u_in
  # This enforces a density and velocity at the boundary
  u[1, ] <- u_in
  v[1, ] <- 0
  # Density at inlet is calculated to satisfy the BCs
  rho[1, ] <- (1 / (1 - u_in)) *
    (sum(f[c(2, 5, 8), 1, ]) + 2 * sum(f[c(4, 3, 7), 2, ]) + u_in * rho[1, ])
  # Re-calculate f for incoming directions at inlet
  f[1, 1, ] <- weights[1] *
    rho[1, ] +
    weights[3] * rho[1, ] * (-u_in) +
    weights[5] * rho[1, ] * u_in
  f[3, 1, ] <- weights[3] * rho[1, ] + weights[3] * rho[1, ] * (-u_in)
  f[5, 1, ] <- weights[5] * rho[1, ] + weights[5] * rho[1, ] * u_in

  # Outlet (Right): Zero-gradient (outflow)
  u[nx, ] <- u[nx - 1, ]
  v[nx, ] <- v[nx - 1, ]
  rho[nx, ] <- rho[nx - 1, ]

  # --- COLLISION STEP ---
  # Calculate equilibrium distribution f_eq
  u_sq <- u^2 + v^2
  for (i in 1:9) {
    # Dot product of velocity and lattice vector
    u_dot_e <- u * ex[i] + v * ey[i]
    f_eq[i, , ] <- weights[i] *
      rho *
      (1 + 3 * u_dot_e + 4.5 * u_dot_e^2 - 1.5 * u_sq)
  }

  # BGK Collision: f_new = f - omega * (f - f_eq)
  f <- f - omega * (f - f_eq)

  # --- BOUNCE-BACK (OBSTACLE) ---
  # Simple bounce-back: particles hitting the obstacle are reflected back
  # We swap the distribution functions for opposite directions
  # Note: R is 1-based, so direction 2 (ex=1) is opposite to 4 (ex=-1)
  opposite_dirs <- c(4, 3, 2, 1, 6, 8, 7, 5) # Opposite of i is opposite_dirs[i]
  for (i in 1:9) {
    f[i, , ][obstacle == 1] <- f[opposite_dirs[i], , ][obstacle == 1]
  }

  # --- STREAMING STEP (VECTORIZED) ---
  # This is the key optimization. We shift entire matrices at once.
  f_new <- array(0, dim = c(9, nx, ny))
  for (i in 1:9) {
    # Use circshift from the 'png' package or a custom one
    # A simple base R implementation of circular shift:
    shift_x <- ex[i]
    shift_y <- ey[i]
    f_new[i, , ] <- t(apply(f[i, , ], 1, function(row) {
      c(tail(row, shift_x), head(row, -shift_x))
    }))
    f_new[i, , ] <- apply(f_new[i, , ], 2, function(col) {
      c(tail(col, shift_y), head(col, -shift_y))
    })
  }
  f <- f_new

  # --- VORTICITY CALCULATION & STORAGE ---
  # Vorticity = dv/dx - du/dy
  # We use finite differences on the velocity field
  du_dy <- (u[, -1, ] - u[, -ny, ]) / 1 # dy=1
  dv_dx <- (v[-1, ] - v[-nx, ]) / 1 # dx=1

  # Pad to match original grid size
  du_dy_padded <- cbind(du_dy, du_dy[, ny])
  dv_dx_padded <- rbind(dv_dx, dv_dx[nx, ])

  vorticity <- dv_dx_padded - du_dy_padded
  vorticity[obstacle == 1] <- NA # Hide vorticity inside obstacle

  # Store frame for animation
  if (step %% frame_interval == 0) {
    frame_count <- frame_count + 1
    vort_df <- as.data.frame(vorticity)
    names(vort_df) <- 1:ny
    vort_df$x <- 1:nx
    vort_long <- pivot_longer(
      vort_df,
      cols = -x,
      names_to = "y",
      values_to = "vorticity"
    )
    vort_long$y <- as.numeric(vort_long$y)
    vort_long$frame <- frame_count
    vorticity_history[[frame_count]] <- vort_long
    cat(sprintf("Step %d/%d | Frame %d\n", step, n_steps, frame_count))
  }
}

cat("Simulation finished. Preparing animation...\n")

# 4. VISUALIZATION: Create the animated GIF
# -----------------------------------------
# Combine all stored frames into one large dataframe
anim_df <- bind_rows(vorticity_history)

# Create the plot
anim <- ggplot(anim_df, aes(x = x, y = y, fill = vorticity)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    name = "Vorticity"
  ) +
  coord_fixed(ratio = 1) + # Essential for correct aspect ratio
  theme_void() +
  labs(
    title = "Von Kármán Vortex Street",
    subtitle = "Time Step: {frame_time}"
  ) +
  transition_time(frame) +
  shadow_wake(wake_length = 0.1, alpha = FALSE) # Adds a nice trailing effect

# Render and save the animation
anim_save(
  "von_karman_vortex_street.gif",
  animation = anim,
  width = 800,
  height = 400,
  nframes = frame_count,
  fps = 20
)

cat("Animation saved to von_karman_vortex_street.gif\n")
