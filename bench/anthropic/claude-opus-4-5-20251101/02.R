# Von Kármán Vortex Street Simulation
# 2D Lattice Boltzmann Method (D2Q9) in R
# Simulates fluid flow past a circular obstacle

# Load required libraries
library(ggplot2)
library(reshape2)
library(gifski) # For creating GIFs

# Check for gganimate, use manual approach if not available
use_gganimate <- requireNamespace("gganimate", quietly = TRUE)
if (use_gganimate) {
  library(gganimate)
  cat("Using gganimate for animation\n")
} else {
  cat("gganimate not found - will create GIF manually with gifski\n")
}

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================

cat("============================================\n")
cat("VON KÁRMÁN VORTEX STREET SIMULATION\n")
cat("Lattice Boltzmann Method (D2Q9)\n")
cat("============================================\n\n")

# Grid dimensions
Nx <- 200 # Extended width for better vortex development
Ny <- 80 # Height
Re <- 100 # Reynolds number (controls vortex shedding)

# Obstacle parameters
cx <- 40 # Obstacle center x
cy <- Ny / 2 # Obstacle center y (middle)
r <- 8 # Obstacle radius

# Flow parameters
u_inlet <- 0.08 # Inlet velocity (keep small for stability)
nu <- u_inlet * 2 * r / Re # Kinematic viscosity from Reynolds number
omega <- 1 / (3 * nu + 0.5) # Relaxation parameter (BGK)

# Time parameters
n_steps <- 8000 # Total simulation steps
save_interval <- 40 # Save frame every N steps
n_frames <- n_steps / save_interval

cat(sprintf("Grid: %d x %d\n", Nx, Ny))
cat(sprintf("Obstacle: center=(%d, %d), radius=%d\n", cx, cy, r))
cat(sprintf("Reynolds number: %d\n", Re))
cat(sprintf("Inlet velocity: %.3f\n", u_inlet))
cat(sprintf("Viscosity: %.6f\n", nu))
cat(sprintf("Relaxation (omega): %.4f\n", omega))
cat(sprintf("Total steps: %d\n", n_steps))
cat(sprintf("Frames to save: %d\n\n", n_frames))

# ============================================================================
# D2Q9 LATTICE CONFIGURATION
# ============================================================================

# D2Q9 velocity vectors (9 directions)
#   6  2  5
#    \ | /
#   3--0--1
#    / | \
#   7  4  8

# Velocity components (indexed 1-9 for R)
ex <- c(0, 1, 0, -1, 0, 1, -1, -1, 1)
ey <- c(0, 0, 1, 0, -1, 1, 1, -1, -1)

# Weights for equilibrium distribution
w <- c(4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36)

# Opposite direction indices (for bounce-back)
opp <- c(1, 4, 5, 2, 3, 8, 9, 6, 7) # 0<->0, 1<->3, 2<->4, etc.

# ============================================================================
# INITIALIZE ARRAYS
# ============================================================================

cat("Initializing simulation arrays...\n")

# Distribution functions: f[x, y, direction]
f <- array(0, dim = c(Nx, Ny, 9))
f_eq <- array(0, dim = c(Nx, Ny, 9))
f_temp <- array(0, dim = c(Nx, Ny, 9))

# Macroscopic variables
rho <- matrix(1, nrow = Nx, ncol = Ny) # Density
ux <- matrix(u_inlet, nrow = Nx, ncol = Ny) # Velocity x
uy <- matrix(0, nrow = Nx, ncol = Ny) # Velocity y

# Create obstacle mask (TRUE where obstacle exists)
obstacle <- matrix(FALSE, nrow = Nx, ncol = Ny)
for (i in 1:Nx) {
  for (j in 1:Ny) {
    if ((i - cx)^2 + (j - cy)^2 <= r^2) {
      obstacle[i, j] <- TRUE
    }
  }
}

# Set velocity to zero inside obstacle
ux[obstacle] <- 0
uy[obstacle] <- 0

cat(sprintf("Obstacle cells: %d\n", sum(obstacle)))

# ============================================================================
# EQUILIBRIUM DISTRIBUTION FUNCTION
# ============================================================================

compute_equilibrium <- function(rho, ux, uy) {
  # Compute equilibrium distribution for all cells
  f_eq <- array(0, dim = c(Nx, Ny, 9))

  # Precompute common terms
  u_sq <- ux^2 + uy^2

  for (k in 1:9) {
    # Velocity dot product
    eu <- ex[k] * ux + ey[k] * uy

    # Equilibrium formula
    f_eq[,, k] <- w[k] * rho * (1 + 3 * eu + 4.5 * eu^2 - 1.5 * u_sq)
  }

  return(f_eq)
}

# Initialize with equilibrium distribution
f <- compute_equilibrium(rho, ux, uy)

# ============================================================================
# STREAMING STEP (VECTORIZED)
# ============================================================================

streaming_step <- function(f) {
  f_new <- array(0, dim = c(Nx, Ny, 9))

  # Direction 0: no movement
  f_new[,, 1] <- f[,, 1]

  # Direction 1: (+1, 0) - shift right
  f_new[2:Nx, , 2] <- f[1:(Nx - 1), , 2]
  f_new[1, , 2] <- f[Nx, , 2] # Periodic/inlet

  # Direction 2: (0, +1) - shift up
  f_new[, 2:Ny, 3] <- f[, 1:(Ny - 1), 3]
  f_new[, 1, 3] <- f[, Ny, 3] # Periodic top-bottom

  # Direction 3: (-1, 0) - shift left
  f_new[1:(Nx - 1), , 4] <- f[2:Nx, , 4]
  f_new[Nx, , 4] <- f[1, , 4] # Outlet

  # Direction 4: (0, -1) - shift down
  f_new[, 1:(Ny - 1), 5] <- f[, 2:Ny, 5]
  f_new[, Ny, 5] <- f[, 1, 5] # Periodic top-bottom

  # Direction 5: (+1, +1) - shift right-up
  f_new[2:Nx, 2:Ny, 6] <- f[1:(Nx - 1), 1:(Ny - 1), 6]
  f_new[1, 2:Ny, 6] <- f[Nx, 1:(Ny - 1), 6]
  f_new[2:Nx, 1, 6] <- f[1:(Nx - 1), Ny, 6]
  f_new[1, 1, 6] <- f[Nx, Ny, 6]

  # Direction 6: (-1, +1) - shift left-up
  f_new[1:(Nx - 1), 2:Ny, 7] <- f[2:Nx, 1:(Ny - 1), 7]
  f_new[Nx, 2:Ny, 7] <- f[1, 1:(Ny - 1), 7]
  f_new[1:(Nx - 1), 1, 7] <- f[2:Nx, Ny, 7]
  f_new[Nx, 1, 7] <- f[1, Ny, 7]

  # Direction 7: (-1, -1) - shift left-down
  f_new[1:(Nx - 1), 1:(Ny - 1), 8] <- f[2:Nx, 2:Ny, 8]
  f_new[Nx, 1:(Ny - 1), 8] <- f[1, 2:Ny, 8]
  f_new[1:(Nx - 1), Ny, 8] <- f[2:Nx, 1, 8]
  f_new[Nx, Ny, 8] <- f[1, 1, 8]

  # Direction 8: (+1, -1) - shift right-down
  f_new[2:Nx, 1:(Ny - 1), 9] <- f[1:(Nx - 1), 2:Ny, 9]
  f_new[1, 1:(Ny - 1), 9] <- f[Nx, 2:Ny, 9]
  f_new[2:Nx, Ny, 9] <- f[1:(Nx - 1), 1, 9]
  f_new[1, Ny, 9] <- f[Nx, 1, 9]

  return(f_new)
}

# ============================================================================
# COMPUTE VORTICITY (CURL)
# ============================================================================

compute_vorticity <- function(ux, uy) {
  # Vorticity = duy/dx - dux/dy
  # Using central differences

  vorticity <- matrix(0, nrow = Nx, ncol = Ny)

  # Interior points (central difference)
  for (i in 2:(Nx - 1)) {
    for (j in 2:(Ny - 1)) {
      duy_dx <- (uy[i + 1, j] - uy[i - 1, j]) / 2
      dux_dy <- (ux[i, j + 1] - ux[i, j - 1]) / 2
      vorticity[i, j] <- duy_dx - dux_dy
    }
  }

  return(vorticity)
}

# ============================================================================
# BOUNDARY CONDITIONS
# ============================================================================

apply_boundary_conditions <- function(f, rho, ux, uy) {
  # --- Inlet (left boundary, x = 1): Zou-He velocity BC ---
  # Set velocity to inlet profile with small perturbation
  i <- 1
  for (j in 1:Ny) {
    # Parabolic profile with perturbation
    y_rel <- (j - Ny / 2) / (Ny / 2)
    u_in <- u_inlet * (1 - 0.5 * y_rel^2) # Parabolic profile

    # Small perturbation to trigger instability
    #u_in <- u_in * (1 + 0.01 * sin(2 * pi * j / Ny))

    rho_in <- 1.0

    # Zou-He formulation for inlet
    f[i, j, 2] <- f[i, j, 4] + (2 / 3) * rho_in * u_in
    f[i, j, 6] <- f[i, j, 8] -
      0.5 * (f[i, j, 3] - f[i, j, 5]) +
      (1 / 6) * rho_in * u_in
    f[i, j, 9] <- f[i, j, 7] +
      0.5 * (f[i, j, 3] - f[i, j, 5]) +
      (1 / 6) * rho_in * u_in
  }

  # --- Outlet (right boundary, x = Nx): Zero gradient ---
  f[Nx, , ] <- f[Nx - 1, , ]

  # --- Top and Bottom: Bounce-back (no-slip walls) ---
  # Bottom wall (j = 1)
  f[, 1, 3] <- f[, 1, 5] # 2 <- 4
  f[, 1, 6] <- f[, 1, 8] # 5 <- 7
  f[, 1, 7] <- f[, 1, 9] # 6 <- 8

  # Top wall (j = Ny)
  f[, Ny, 5] <- f[, Ny, 3] # 4 <- 2
  f[, Ny, 8] <- f[, Ny, 6] # 7 <- 5
  f[, Ny, 9] <- f[, Ny, 7] # 8 <- 6

  return(f)
}

# ============================================================================
# OBSTACLE BOUNCE-BACK
# ============================================================================

apply_obstacle_bounce_back <- function(f, obstacle) {
  # Bounce-back at obstacle boundary
  f_copy <- f

  for (i in 1:Nx) {
    for (j in 1:Ny) {
      if (obstacle[i, j]) {
        # Reverse all directions
        for (k in 1:9) {
          f[i, j, k] <- f_copy[i, j, opp[k]]
        }
      }
    }
  }

  return(f)
}

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

cat("\nStarting simulation...\n")
cat("Progress: ")

# Storage for animation frames
frames_data <- list()
frame_count <- 0

# Add small perturbation to trigger vortex shedding
set.seed(42)
uy <- uy + 0.001 * matrix(rnorm(Nx * Ny), nrow = Nx, ncol = Ny)
uy[obstacle] <- 0

# Main time loop
start_time <- Sys.time()

for (step in 1:n_steps) {
  # --- 1. Collision Step (BGK) ---
  # Compute macroscopic quantities
  rho <- apply(f, c(1, 2), sum)
  ux <- (f[,, 2] - f[,, 4] + f[,, 6] - f[,, 7] - f[,, 8] + f[,, 9]) / rho
  uy <- (f[,, 3] - f[,, 5] + f[,, 6] + f[,, 7] - f[,, 8] - f[,, 9]) / rho

  # Set velocity to zero in obstacle
  ux[obstacle] <- 0
  uy[obstacle] <- 0

  # Compute equilibrium
  f_eq <- compute_equilibrium(rho, ux, uy)

  # BGK collision
  f <- f - omega * (f - f_eq)

  # --- 2. Streaming Step ---
  f <- streaming_step(f)

  # --- 3. Boundary Conditions ---
  f <- apply_boundary_conditions(f, rho, ux, uy)

  # --- 4. Obstacle Bounce-back ---
  f <- apply_obstacle_bounce_back(f, obstacle)

  # --- 5. Save Frame ---
  if (step %% save_interval == 0) {
    frame_count <- frame_count + 1

    # Compute vorticity
    vorticity <- compute_vorticity(ux, uy)
    vorticity[obstacle] <- NA

    # Store frame data
    frames_data[[frame_count]] <- list(
      step = step,
      vorticity = vorticity,
      ux = ux,
      uy = uy,
      rho = rho
    )

    # Progress indicator
    if (frame_count %% 20 == 0) {
      cat(sprintf("%d/%d ", frame_count, n_frames))
    }
  }
}

end_time <- Sys.time()
cat("\n\nSimulation completed!\n")
cat(sprintf(
  "Time elapsed: %.1f seconds\n",
  difftime(end_time, start_time, units = "secs")
))

# ============================================================================
# CREATE ANIMATION
# ============================================================================

cat("\nGenerating animation frames...\n")

# Create temporary directory for frames
dir.create("frames", showWarnings = FALSE)

# Color scale limits for vorticity
vort_max <- 0.015 # Adjusted for visibility

# Generate individual frames
for (i in seq_along(frames_data)) {
  frame <- frames_data[[i]]

  # Prepare data for plotting
  vort_df <- melt(frame$vorticity)
  colnames(vort_df) <- c("x", "y", "vorticity")

  # Remove obstacle points
  vort_df$vorticity[is.na(vort_df$vorticity)] <- NA

  # Create the plot
  p <- ggplot(vort_df, aes(x = x, y = y, fill = vorticity)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-vort_max, vort_max),
      oob = scales::squish,
      name = "Vorticity\n(curl)"
    ) +
    # Add obstacle
    annotate(
      "point",
      x = cx,
      y = cy,
      size = r * 1.5,
      shape = 21,
      fill = "gray30",
      color = "black"
    ) +
    labs(
      title = "Von Kármán Vortex Street",
      subtitle = sprintf(
        "Lattice Boltzmann D2Q9 | Re = %d | Step: %d",
        Re,
        frame$step
      ),
      x = "x",
      y = "y"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "right",
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    coord_fixed(ratio = 1, expand = FALSE)

  # Save frame
  ggsave(
    sprintf("frames/frame_%04d.png", i),
    p,
    width = 10,
    height = 4,
    dpi = 100
  )

  if (i %% 20 == 0) cat(sprintf("Frame %d/%d saved\n", i, length(frames_data)))
}

cat("\nCreating GIF animation...\n")

# Create GIF from frames
png_files <- list.files(
  "frames",
  pattern = "frame_.*\\.png$",
  full.names = TRUE
)
png_files <- sort(png_files)

gifski(
  png_files,
  gif_file = "vonkarman_vortex_street.gif",
  width = 1000,
  height = 400,
  delay = 0.1 # 10 fps
)

cat("Animation saved as 'vonkarman_vortex_street.gif'\n")

# Clean up frames
unlink("frames", recursive = TRUE)

# ============================================================================
# STATIC VISUALIZATION: Final State
# ============================================================================

cat("\nCreating static visualizations...\n")

# Final frame vorticity plot
final_frame <- frames_data[[length(frames_data)]]
vort_df <- melt(final_frame$vorticity)
colnames(vort_df) <- c("x", "y", "vorticity")

p_final <- ggplot(vort_df, aes(x = x, y = y, fill = vorticity)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(
    low = "#0000FF",
    mid = "white",
    high = "#FF0000",
    midpoint = 0,
    limits = c(-vort_max, vort_max),
    oob = scales::squish,
    name = "Vorticity"
  ) +
  annotate(
    "point",
    x = cx,
    y = cy,
    size = r * 1.5,
    shape = 21,
    fill = "gray20",
    color = "black",
    stroke = 1
  ) +
  labs(
    title = "Von Kármán Vortex Street - Final State",
    subtitle = sprintf(
      "Lattice Boltzmann D2Q9 | Re = %d | Grid: %d×%d",
      Re,
      Nx,
      Ny
    ),
    x = "Downstream Distance",
    y = "Cross-stream Position",
    caption = "Red: clockwise rotation | Blue: counter-clockwise rotation"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    plot.caption = element_text(size = 9, color = "gray50"),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  coord_fixed(ratio = 1, expand = FALSE)

print(p_final)
ggsave(
  "vonkarman_final_vorticity.png",
  p_final,
  width = 12,
  height = 5,
  dpi = 150
)

# ============================================================================
# VELOCITY MAGNITUDE PLOT
# ============================================================================

vel_mag <- sqrt(final_frame$ux^2 + final_frame$uy^2)
vel_mag[obstacle] <- NA

vel_df <- melt(vel_mag)
colnames(vel_df) <- c("x", "y", "velocity")

p_vel <- ggplot(vel_df, aes(x = x, y = y, fill = velocity)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "Velocity\nMagnitude",
    na.value = "gray30"
  ) +
  annotate(
    "point",
    x = cx,
    y = cy,
    size = r * 1.5,
    shape = 21,
    fill = "gray20",
    color = "white",
    stroke = 1
  ) +
  labs(
    title = "Flow Velocity Magnitude",
    subtitle = sprintf(
      "Inlet velocity: %.3f | Reynolds number: %d",
      u_inlet,
      Re
    ),
    x = "x",
    y = "y"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  coord_fixed(ratio = 1, expand = FALSE)

print(p_vel)
ggsave("vonkarman_velocity.png", p_vel, width = 12, height = 5, dpi = 150)

# ============================================================================
# STREAMLINES PLOT (APPROXIMATE)
# ============================================================================

# Downsample for streamline arrows
skip <- 6
x_arrow <- seq(1, Nx, by = skip)
y_arrow <- seq(1, Ny, by = skip)

arrow_df <- expand.grid(x = x_arrow, y = y_arrow)
arrow_df$ux <- final_frame$ux[cbind(arrow_df$x, arrow_df$y)]
arrow_df$uy <- final_frame$uy[cbind(arrow_df$x, arrow_df$y)]

# Remove arrows in/near obstacle
arrow_df$dist <- sqrt((arrow_df$x - cx)^2 + (arrow_df$y - cy)^2)
arrow_df <- arrow_df[arrow_df$dist > r + 2, ]

# Scale arrows
arrow_scale <- 3

p_stream <- ggplot() +
  geom_raster(
    data = vort_df,
    aes(x = x, y = y, fill = vorticity),
    interpolate = TRUE
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-vort_max, vort_max),
    oob = scales::squish,
    name = "Vorticity"
  ) +
  geom_segment(
    data = arrow_df,
    aes(
      x = x,
      y = y,
      xend = x + arrow_scale * ux / u_inlet,
      yend = y + arrow_scale * uy / u_inlet
    ),
    arrow = arrow(length = unit(0.08, "cm")),
    color = "black",
    alpha = 0.6,
    linewidth = 0.3
  ) +
  annotate(
    "point",
    x = cx,
    y = cy,
    size = r * 1.5,
    shape = 21,
    fill = "gray20",
    color = "white",
    stroke = 1
  ) +
  labs(
    title = "Vorticity Field with Velocity Vectors",
    subtitle = "Alternating vortices shed from the cylinder",
    x = "x",
    y = "y"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  coord_fixed(ratio = 1, expand = FALSE)

print(p_stream)
ggsave("vonkarman_streamlines.png", p_stream, width = 12, height = 5, dpi = 150)

# ============================================================================
# TIME SERIES: Vorticity at a Probe Point
# ============================================================================

# Extract vorticity time series at a point behind the cylinder
probe_x <- cx + 30
probe_y <- cy

vort_timeseries <- sapply(frames_data, function(f) {
  f$vorticity[probe_x, probe_y]
})

time_df <- data.frame(
  step = sapply(frames_data, function(f) f$step),
  vorticity = vort_timeseries
)

p_time <- ggplot(time_df, aes(x = step, y = vorticity)) +
  geom_line(color = "purple", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Vorticity Time Series at Probe Point",
    subtitle = sprintf(
      "Probe location: (%d, %d) - downstream of cylinder",
      probe_x,
      probe_y
    ),
    x = "Time Step",
    y = "Vorticity"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5)
  )

print(p_time)
ggsave("vonkarman_timeseries.png", p_time, width = 10, height = 5, dpi = 150)

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n============================================\n")
cat("SIMULATION SUMMARY\n")
cat("============================================\n")
cat("Physics: 2D Incompressible Flow Past Cylinder\n")
cat("Method: Lattice Boltzmann D2Q9 with BGK collision\n")
cat(sprintf("Reynolds Number: %d\n", Re))
cat(sprintf("Grid: %d x %d\n", Nx, Ny))
cat(sprintf("Obstacle: Cylinder at (%d, %d), radius %d\n", cx, cy, r))
cat("\nOutput Files:\n")
cat("  1. vonkarman_vortex_street.gif - Animation of vortex shedding\n")
cat("  2. vonkarman_final_vorticity.png - Final vorticity field\n")
cat("  3. vonkarman_velocity.png - Velocity magnitude\n")
cat("  4. vonkarman_streamlines.png - Vorticity with velocity vectors\n")
cat("  5. vonkarman_timeseries.png - Vorticity oscillation at probe\n")
cat("\nVortex Street Characteristics:\n")
cat("  - Red regions: Clockwise rotation (negative vorticity)\n")
cat("  - Blue regions: Counter-clockwise rotation (positive vorticity)\n")
cat("  - Alternating vortices form the characteristic 'street' pattern\n")
cat("============================================\n")
