# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(gganimate)

# ==========================================
# 1. CONFIGURATION & LATTICE SETUP
# ==========================================

# Grid Dimensions
nx <- 80
ny <- 40

# Physical Parameters
u_inlet <- 0.1 # Horizontal inflow velocity
tau <- 0.6 # Relaxation time (controls viscosity). Viscosity = (tau - 0.5)/3
omega <- 1 / tau # Relaxation parameter

# Simulation Steps
nt <- 2000 # Total time steps
snapshot_freq <- 40 # Save frame every N steps

# Lattice Weights (D2Q9)
# 1 Center, 4 Cardinal, 4 Diagonal
w <- c(4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36)

# Lattice Directions (D2Q9)
# indices: 1:C, 2:E, 3:N, 4:W, 5:S, 6:NE, 7:NW, 8:SW, 9:SE
cx <- c(0, 1, 0, -1, 0, 1, -1, -1, 1)
cy <- c(0, 0, 1, 0, -1, 1, 1, -1, -1)

# Opposite directions (for bounce-back)
# e.g., The opposite of East (2) is West (4)
opposite <- c(1, 4, 5, 2, 3, 8, 9, 6, 7)

# ==========================================
# 2. INITIALIZATION
# ==========================================

# Create Obstacle (Cylinder)
# x matrix and y matrix for coordinate calculations
x_grid <- matrix(rep(1:nx, ny), nrow = nx, ncol = ny)
y_grid <- t(matrix(rep(1:ny, nx), nrow = ny, ncol = nx))

# Circle at (20, 20) radius 3
obstacle <- ((x_grid - 20)^2 + (y_grid - 20)^2) <= 3^2

# Initial Flow State (Start with uniform flow to speed up stabilization)
rho <- matrix(1, nx, ny)
ux <- matrix(u_inlet, nx, ny)
uy <- matrix(0, nx, ny)

# Initialize Distribution Function (f) to Equilibrium
f <- array(0, dim = c(nx, ny, 9))

# Function to calculate Equilibrium (Vectorized)
calc_eq <- function(rho, ux, uy) {
  usq <- ux^2 + uy^2
  feq <- array(0, dim = c(nx, ny, 9))
  for (i in 1:9) {
    cu <- cx[i] * ux + cy[i] * uy
    feq[,, i] <- rho * w[i] * (1 + 3 * cu + 4.5 * cu^2 - 1.5 * usq)
  }
  return(feq)
}

f <- calc_eq(rho, ux, uy)

# Storage for animation
animation_data <- list()
frame_counter <- 1

# ==========================================
# 3. MAIN LBM LOOP
# ==========================================
cat("Starting Simulation... this may take a moment.\n")

for (t in 1:nt) {
  # --- A. MACROSCOPIC VARIABLES ---
  # rho = sum(f), u = sum(f*c) / rho
  rho <- apply(f, c(1, 2), sum)

  # Calculate Momentum
  mx <- matrix(0, nx, ny)
  my <- matrix(0, nx, ny)
  for (i in 1:9) {
    mx <- mx + f[,, i] * cx[i]
    my <- my + f[,, i] * cy[i]
  }
  ux <- mx / rho
  uy <- my / rho

  # --- B. INLET/OUTLET BOUNDARIES (Macroscopic level) ---
  # Inlet (Left, x=1): constant velocity
  ux[1, ] <- u_inlet
  uy[1, ] <- 0
  rho[1, ] <- 1 # Approx density at inlet

  # Outlet (Right, x=nx): Zero gradient (copy neighbor)
  ux[nx, ] <- ux[nx - 1, ]
  uy[nx, ] <- uy[nx - 1, ]
  rho[nx, ] <- rho[nx - 1, ]

  # Re-calculate Equilibrium with enforced BCs
  feq <- calc_eq(rho, ux, uy)

  # --- C. COLLISION (BGK) ---
  # f_out = f_in - (f_in - f_eq) / tau
  # Note: We apply collision *before* streaming in this implementation structure
  f_out <- f - (f - feq) * omega

  # --- D. OBSTACLE BOUNCE-BACK ---
  # For solid nodes, outgoing particles bounce back to opposite direction.
  # We do this by extracting the PRE-collision values.
  for (i in 1:9) {
    # If a cell is an obstacle, replace its post-collision value
    # with the value from the opposite direction
    op_idx <- opposite[i]
    f_out[,, i][obstacle] <- f[,, op_idx][obstacle]
  }

  # --- E. STREAMING (Vectorized Shift) ---
  # Move particles to neighboring cells
  f_streamed <- f_out # Copy structure

  for (i in 1:9) {
    # Current matrix layer
    mat <- f_out[,, i]

    # Shift in X
    dx_shift <- cx[i]
    if (dx_shift == 1) {
      mat <- rbind(mat[nx, ], mat[1:(nx - 1), ]) # Periodic (safe for now)
      # Enforce non-periodic inlet for strictness:
      # But standard LBM practice often uses index shifting:
      mat[2:nx, ] <- f_out[1:(nx - 1), , i]
    } else if (dx_shift == -1) {
      mat[1:(nx - 1), ] <- f_out[2:nx, , i]
    }

    # Shift in Y
    dy_shift <- cy[i]
    if (dy_shift == 1) {
      mat[, 2:ny] <- mat[, 1:(ny - 1)]
    } else if (dy_shift == -1) {
      mat[, 1:(ny - 1)] <- mat[, 2:ny]
    }

    f_streamed[,, i] <- mat
  }

  # Update main array
  f <- f_streamed

  # --- F. VORTICITY & STORAGE ---
  if (t %% snapshot_freq == 0) {
    cat(sprintf("Step: %d / %d\n", t, nt))

    # Calculate Curl (Vorticity): dv/dx - du/dy
    # Simple central difference for internal nodes
    # curl = (uy(x+1) - uy(x-1))/2dx - (ux(y+1) - ux(y-1))/2dy

    # X-gradients of Uy
    duy_dx <- (rbind(uy[2:nx, ], uy[1, ]) - rbind(uy[nx, ], uy[1:(nx - 1), ])) /
      2
    # Y-gradients of Ux
    dux_dy <- (cbind(ux[, 2:ny], ux[, 1]) - cbind(ux[, ny], ux[, 1:(ny - 1)])) /
      2

    curl <- duy_dx - dux_dy

    # Mask the obstacle in visualization
    curl[obstacle] <- NA

    # Convert to Data Frame
    df_frame <- expand.grid(x = 1:nx, y = 1:ny)
    df_frame$vorticity <- as.vector(curl)
    df_frame$frame <- t

    animation_data[[frame_counter]] <- df_frame
    frame_counter <- frame_counter + 1
  }
}

# Combine all frames
plot_data <- bind_rows(animation_data)

# ==========================================
# 4. VISUALIZATION
# ==========================================

# Circle polygon for plotting
circle_dat <- data.frame(
  x = 20 + 3 * cos(seq(0, 2 * pi, length.out = 50)),
  y = 20 + 3 * sin(seq(0, 2 * pi, length.out = 50))
)

# Create Animation
p <- ggplot(plot_data, aes(x = x, y = y, fill = vorticity)) +
  geom_raster(interpolate = TRUE) + # Smoother look
  geom_polygon(
    data = circle_dat,
    aes(x = x, y = y),
    fill = "black",
    inherit.aes = FALSE
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-0.02, 0.02),
    oob = scales::squish,
    name = "Vorticity\n(Curl)"
  ) +
  labs(
    title = "Von Kármán Vortex Street (Lattice Boltzmann)",
    subtitle = "Frame: {frame_time}",
    x = NULL,
    y = NULL
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  transition_time(frame) +
  ease_aes('linear')

# Render GIF
animate(p, fps = 10, width = 800, height = 400, renderer = gifski_renderer())
