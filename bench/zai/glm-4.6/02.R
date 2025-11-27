# --------------------------------------------------------------
# Von Kármán Vortex Street Simulation with Lattice Boltzmann (D2Q9)
# --------------------------------------------------------------

# 1. LOAD LIBRARIES -------------------------------------------------
library(ggplot2)
library(gganimate)
library(dplyr)
library(tidyr)

# Optional: For a progress bar
# install.packages("pbapply")
library(pbapply)

# 2. SIMULATION PARAMETERS -------------------------------------------

# Lattice parameters
Nx <- 80 # Lattice width (x-direction)
Ny <- 40 # Lattice height (y-direction)

# Obstacle parameters
obstacle_x <- 20
obstacle_y <- 20
obstacle_r <- 3

# Flow parameters (Reynolds Number)
Re <- 150 # Reynolds number (40-200 is good for vortex street)
u_max <- 0.1 # Maximum inlet velocity (keep < 0.1 for stability)
nu <- u_max * 2 * obstacle_r / Re # Kinematic viscosity
tau <- 3 * nu + 0.5 # Relaxation time for BGK operator

# Time parameters
n_steps <- 4000 # Total simulation steps
plot_every <- 10 # Store a frame for animation every N steps

# 3. D2Q9 LATTICE STRUCTURE ------------------------------------------

# Lattice velocities (c_i)
cx <- c(0, 1, 0, -1, 0, 1, -1, -1, 1)
cy <- c(0, 0, 1, 0, -1, 1, 1, -1, -1)

# Lattice weights (w_i)
w <- c(4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36)

# Opposite directions for bounce-back (1-based indexing)
# [1,2,3,4,5,6,7,8,9] -> [1,4,5,2,3,8,9,6,7]
opp <- c(1, 4, 5, 2, 3, 8, 9, 6, 7)

# 4. INITIALIZATION --------------------------------------------------

# Create coordinate grids
x_coords <- 1:Nx
y_coords <- 1:Ny

# Create obstacle mask (TRUE for obstacle, FALSE for fluid)
is_obstacle <- outer(x_coords, y_coords, function(x, y) {
  ((x - obstacle_x)^2 + (y - obstacle_y)^2) < obstacle_r^2
})

# Initialize distribution functions f(x, y, i)
# f is a 3D array: (x, y, direction)
f <- array(0, dim = c(Nx, Ny, 9))

# Initialize with equilibrium distribution (rho=1, u=v=0)
rho_init <- 1
u_init <- rep(0, Nx * Ny)
v_init <- rep(0, Nx * Ny)
for (i in 1:9) {
  cu <- cx[i] * u_init + cy[i] * v_init
  u_sq <- u_init^2 + v_init^2
  f[,, i] <- rho_init * w[i] * (1 + 3 * cu + 4.5 * cu^2 - 1.5 * u_sq)
}

# Create top/bottom wall mask
is_wall <- array(FALSE, dim = c(Nx, Ny))
is_wall[, 1] <- TRUE # Bottom wall
is_wall[, Ny] <- TRUE # Top wall

# 5. HELPER FUNCTIONS ------------------------------------------------

# Function to calculate equilibrium distribution
calc_feq <- function(rho, u, v) {
  feq <- array(0, dim = c(Nx, Ny, 9))
  for (i in 1:9) {
    cu <- cx[i] * u + cy[i] * v
    u_sq <- u^2 + v^2
    feq[,, i] <- rho * w[i] * (1 + 3 * cu + 4.5 * cu^2 - 1.5 * u_sq)
  }
  return(feq)
}

# 6. MAIN SIMULATION LOOP --------------------------------------------

# Storage for animation frames
vorticity_frames <- list()

cat("Starting LBM simulation...\n")
for (step in 1:n_steps) {
  # --- A. Collision Step (BGK) ---
  # Calculate macroscopic variables (rho, u, v)
  rho <- apply(f, c(1, 2), sum)
  u <- (f[,, 2] - f[,, 4] + f[,, 6] - f[,, 7] + f[,, 9] - f[,, 8]) / rho
  v <- (f[,, 3] - f[,, 5] + f[,, 6] + f[,, 7] - f[,, 8] - f[,, 9]) / rho

  # Apply boundary conditions for macroscopic variables
  u[is_obstacle] <- 0
  v[is_obstacle] <- 0
  u[is_wall] <- 0
  v[is_wall] <- 0

  # Calculate equilibrium distribution
  f_eq <- calc_feq(rho, u, v)

  # Perform collision
  f <- f - (1 / tau) * (f - f_eq)

  # --- B. Streaming Step ---
  # This is the vectorized part. We shift the distribution functions.
  f_new <- array(0, dim = dim(f))
  for (i in 1:9) {
    # Shift the matrix for direction i
    f_new[,, i] <- f[,, i]
    if (cx[i] == 1) {
      f_new[2:Nx, , i] <- f[1:(Nx - 1), , i]
    }
    if (cx[i] == -1) {
      f_new[1:(Nx - 1), , i] <- f[2:Nx, , i]
    }
    if (cy[i] == 1) {
      f_new[, 2:Ny, i] <- f[, 1:(Ny - 1), i]
    }
    if (cy[i] == -1) f_new[, 1:(Ny - 1), i] <- f[, 2:Ny, i]
  }
  f <- f_new

  # --- C. Boundary Conditions ---

  # Bounce-back on obstacle and walls
  # Find all nodes that are solid (obstacle or wall)
  is_solid <- is_obstacle | is_wall
  for (i in 1:9) {
    f[is_solid, i] <- f[is_solid, opp[i]]
  }

  # Inlet (Left boundary): Zou-He velocity boundary
  x_inlet <- 1
  u_inlet <- u_max
  v_inlet <- 0

  # Known distributions at inlet
  f_in_2 <- f[x_inlet, , 2]
  f_in_5 <- f[x_inlet, , 5]
  f_in_6 <- f[x_inlet, , 6]

  # Unknown distributions to be calculated
  rho_inlet <- (f_in_2 +
    f_in_5 +
    f_in_6 +
    2 * (f[x_inlet, , 3] + f[x_inlet, , 7] + f[x_inlet, , 8])) /
    (1 - u_inlet)

  f[x_inlet, , 1] <- f_in_2 - 2 / 3 * rho_inlet * u_inlet
  f[x_inlet, , 3] <- f_in_5 -
    1 / 6 * rho_inlet * u_inlet +
    0.5 * (f[x_inlet, , 2] - f[x_inlet, , 4])
  f[x_inlet, , 9] <- f_in_6 -
    1 / 6 * rho_inlet * u_inlet -
    0.5 * (f[x_inlet, , 2] - f[x_inlet, , 4])

  # Outlet (Right boundary): Simple extrapolation
  f[Nx, , ] <- f[Nx - 1, , ]

  # --- D. Calculate and Store Vorticity ---
  if (step %% plot_every == 0) {
    # Calculate vorticity (curl of velocity)
    # omega = dv/dx - du/dy
    # Use central differences for interior points
    du_dy <- (u[3:Nx, ] - u[1:(Nx - 2), ]) / 2
    dv_dx <- (v[, 3:Ny] - v[, 1:(Ny - 2)]) / 2

    # Vorticity is calculated on a staggered grid
    omega <- dv_dx[2:(Nx - 1), ] - du_dy[, 2:(Ny - 1)]

    # Pad with NA to match original grid size for plotting
    omega_padded <- array(NA, dim = c(Nx, Ny))
    omega_padded[2:(Nx - 1), 2:(Ny - 1)] <- omega

    vorticity_frames[[length(vorticity_frames) + 1]] <- omega_padded
  }

  # Progress update
  # if (step %% 100 == 0) cat(sprintf("Step %d / %d\n", step, n_steps))
}
cat("Simulation finished.\n")


# 7. CREATE ANIMATION ------------------------------------------------

cat("Preparing animation data...\n")
# Convert list of matrices to a long data frame for ggplot
animation_df <- pblapply(seq_along(vorticity_frames), function(i) {
  frame_data <- as.data.frame.table(vorticity_frames[[i]])
  colnames(frame_data) <- c("x", "y", "vorticity")
  frame_data$frame <- i
  frame_data$time <- i * plot_every
  return(frame_data)
}) %>%
  bind_rows()

# Filter out NA values (boundaries)
animation_df <- na.omit(animation_df)

# Create the plot
p <- ggplot(animation_df, aes(x = x, y = y, fill = vorticity)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Vorticity"
  ) +
  coord_fixed() + # Ensure pixels are square
  labs(
    title = "Von Kármán Vortex Street (LBM Simulation)",
    subtitle = "Time step: {frame_time}",
    x = "X (lattice units)",
    y = "Y (lattice units)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right", panel.grid = element_blank())

# Animate the plot
cat("Generating GIF... (this may take a moment)\n")
anim <- animate(
  p,
  renderer = gifski_renderer(),
  width = 800,
  height = 400,
  duration = 15, # Duration of the GIF in seconds
  fps = 20
)

# 8. SAVE THE ANIMATION -----------------------------------------------
anim_save("von_karman_vortex_street.gif", animation = anim)
cat("Animation saved to 'von_karman_vortex_street.gif'\n")

# Display the animation in the RStudio viewer (if available)
print(anim)
