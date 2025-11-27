# 2D Shallow Water Equations Simulation
# "Pond Ripple" - 3D Rayshader Visualization
#
# Equations:
# ∂h/∂t + ∂(hu)/∂x + ∂(hv)/∂y = 0
# ∂(hu)/∂t + ∂(hu² + ½gh²)/∂x + ∂(huv)/∂y = 0
# ∂(hv)/∂t + ∂(huv)/∂x + ∂(hv² + ½gh²)/∂y = 0

# Load required libraries
library(ggplot2)
library(reshape2)

# Install rayshader if not available
if (!requireNamespace("rayshader", quietly = TRUE)) {
  cat("Installing rayshader package...\n")
  install.packages("rayshader")
}
library(rayshader)

# Optional: for better 3D rendering
if (!requireNamespace("rgl", quietly = TRUE)) {
  install.packages("rgl")
}

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================

cat("============================================\n")
cat("2D SHALLOW WATER EQUATIONS\n")
cat("Pond Ripple Simulation\n")
cat("============================================\n\n")

# Grid parameters
Nx <- 200 # Grid points in x
Ny <- 200 # Grid points in y
Lx <- 10.0 # Domain length in x
Ly <- 10.0 # Domain length in y
dx <- Lx / Nx # Spatial step x
dy <- Ly / Ny # Spatial step y

# Physical parameters
g <- 9.81 # Gravitational acceleration
h0 <- 1.0 # Base water depth

# Time parameters
CFL <- 0.4 # CFL number for stability
c_max <- sqrt(g * h0 * 2) # Maximum wave speed estimate
dt <- CFL * min(dx, dy) / c_max
T_final <- 2.5 # Simulation end time
n_steps <- ceiling(T_final / dt)

# Drop parameters
drop_amplitude <- 0.3 # Height of initial drop
drop_sigma <- 0.4 # Width of Gaussian drop
drop_x <- Lx / 2 # Drop center x
drop_y <- Ly / 2 # Drop center y

cat(sprintf("Grid: %d x %d\n", Nx, Ny))
cat(sprintf("Domain: %.1f x %.1f\n", Lx, Ly))
cat(sprintf("dx = dy = %.4f\n", dx))
cat(sprintf("Base depth (h0): %.2f\n", h0))
cat(sprintf("Gravity (g): %.2f\n", g))
cat(sprintf("Max wave speed: %.3f\n", c_max))
cat(sprintf("Time step (dt): %.5f\n", dt))
cat(sprintf("Total steps: %d\n", n_steps))
cat(sprintf("CFL number: %.2f\n\n", CFL))

# ============================================================================
# INITIALIZE FIELDS
# ============================================================================

cat("Initializing water surface...\n")

# Create coordinate grids
x <- seq(dx / 2, Lx - dx / 2, by = dx)
y <- seq(dy / 2, Ly - dy / 2, by = dy)
X <- matrix(rep(x, Ny), nrow = Nx, ncol = Ny)
Y <- matrix(rep(y, each = Nx), nrow = Nx, ncol = Ny)

# Initialize water height with Gaussian drop
h <- h0 +
  drop_amplitude * exp(-((X - drop_x)^2 + (Y - drop_y)^2) / (2 * drop_sigma^2))

# Initialize momentum (hu, hv) - starts at rest
hu <- matrix(0, nrow = Nx, ncol = Ny)
hv <- matrix(0, nrow = Nx, ncol = Ny)

# Store initial state
h_initial <- h

cat(sprintf(
  "Initial drop: amplitude=%.2f, sigma=%.2f\n",
  drop_amplitude,
  drop_sigma
))
cat(sprintf("Initial h: min=%.3f, max=%.3f\n\n", min(h), max(h)))

# ============================================================================
# LAX-WENDROFF SCHEME HELPER FUNCTIONS
# ============================================================================

# Compute fluxes in x-direction
compute_flux_x <- function(h, hu, hv) {
  u <- hu / pmax(h, 1e-10)

  F1 <- hu # Mass flux
  F2 <- hu * u + 0.5 * g * h^2 # x-momentum flux
  F3 <- hu * hv / pmax(h, 1e-10) # y-momentum flux

  return(list(F1 = F1, F2 = F2, F3 = F3))
}

# Compute fluxes in y-direction
compute_flux_y <- function(h, hu, hv) {
  v <- hv / pmax(h, 1e-10)

  G1 <- hv # Mass flux
  G2 <- hu * hv / pmax(h, 1e-10) # x-momentum flux
  G3 <- hv * v + 0.5 * g * h^2 # y-momentum flux

  return(list(G1 = G1, G2 = G2, G3 = G3))
}

# ============================================================================
# TWO-STEP LAX-WENDROFF SCHEME
# ============================================================================

lax_wendroff_step <- function(h, hu, hv, dt, dx, dy) {
  # Get current fluxes
  Fx <- compute_flux_x(h, hu, hv)
  Fy <- compute_flux_y(h, hu, hv)

  # ===== PREDICTOR STEP (half time step, staggered grid) =====

  # Half-step values at cell interfaces (i+1/2, j) and (i, j+1/2)

  # For x-direction (at i+1/2)
  h_half_x <- matrix(0, nrow = Nx - 1, ncol = Ny)
  hu_half_x <- matrix(0, nrow = Nx - 1, ncol = Ny)
  hv_half_x <- matrix(0, nrow = Nx - 1, ncol = Ny)

  for (i in 1:(Nx - 1)) {
    h_half_x[i, ] <- 0.5 *
      (h[i, ] + h[i + 1, ]) -
      0.5 * (dt / dx) * (Fx$F1[i + 1, ] - Fx$F1[i, ])
    hu_half_x[i, ] <- 0.5 *
      (hu[i, ] + hu[i + 1, ]) -
      0.5 * (dt / dx) * (Fx$F2[i + 1, ] - Fx$F2[i, ])
    hv_half_x[i, ] <- 0.5 *
      (hv[i, ] + hv[i + 1, ]) -
      0.5 * (dt / dx) * (Fx$F3[i + 1, ] - Fx$F3[i, ])
  }

  # For y-direction (at j+1/2)
  h_half_y <- matrix(0, nrow = Nx, ncol = Ny - 1)
  hu_half_y <- matrix(0, nrow = Nx, ncol = Ny - 1)
  hv_half_y <- matrix(0, nrow = Nx, ncol = Ny - 1)

  for (j in 1:(Ny - 1)) {
    h_half_y[, j] <- 0.5 *
      (h[, j] + h[, j + 1]) -
      0.5 * (dt / dy) * (Fy$G1[, j + 1] - Fy$G1[, j])
    hu_half_y[, j] <- 0.5 *
      (hu[, j] + hu[, j + 1]) -
      0.5 * (dt / dy) * (Fy$G2[, j + 1] - Fy$G2[, j])
    hv_half_y[, j] <- 0.5 *
      (hv[, j] + hv[, j + 1]) -
      0.5 * (dt / dy) * (Fy$G3[, j + 1] - Fy$G3[, j])
  }

  # Ensure positivity of height
  h_half_x <- pmax(h_half_x, 1e-10)
  h_half_y <- pmax(h_half_y, 1e-10)

  # Compute fluxes at half-step positions
  Fx_half <- compute_flux_x(h_half_x, hu_half_x, hv_half_x)
  Fy_half <- compute_flux_y(h_half_y, hu_half_y, hv_half_y)

  # ===== CORRECTOR STEP (full time step) =====

  h_new <- h
  hu_new <- hu
  hv_new <- hv

  # Interior points update
  for (i in 2:(Nx - 1)) {
    for (j in 2:(Ny - 1)) {
      h_new[i, j] <- h[i, j] -
        (dt / dx) * (Fx_half$F1[i, j] - Fx_half$F1[i - 1, j]) -
        (dt / dy) * (Fy_half$G1[i, j] - Fy_half$G1[i, j - 1])

      hu_new[i, j] <- hu[i, j] -
        (dt / dx) * (Fx_half$F2[i, j] - Fx_half$F2[i - 1, j]) -
        (dt / dy) * (Fy_half$G2[i, j] - Fy_half$G2[i, j - 1])

      hv_new[i, j] <- hv[i, j] -
        (dt / dx) * (Fx_half$F3[i, j] - Fx_half$F3[i - 1, j]) -
        (dt / dy) * (Fy_half$G3[i, j] - Fy_half$G3[i, j - 1])
    }
  }

  # Ensure positivity
  h_new <- pmax(h_new, 1e-10)

  return(list(h = h_new, hu = hu_new, hv = hv_new))
}

# ============================================================================
# BOUNDARY CONDITIONS (Reflective walls)
# ============================================================================

apply_boundary_conditions <- function(h, hu, hv) {
  # Reflective (solid wall) boundary conditions

  # Left boundary (x = 0): reflect hu
  h[1, ] <- h[2, ]
  hu[1, ] <- -hu[2, ]
  hv[1, ] <- hv[2, ]

  # Right boundary (x = Lx): reflect hu
  h[Nx, ] <- h[Nx - 1, ]
  hu[Nx, ] <- -hu[Nx - 1, ]
  hv[Nx, ] <- hv[Nx - 1, ]

  # Bottom boundary (y = 0): reflect hv
  h[, 1] <- h[, 2]
  hu[, 1] <- hu[, 2]
  hv[, 1] <- -hv[, 2]

  # Top boundary (y = Ly): reflect hv
  h[, Ny] <- h[, Ny - 1]
  hu[, Ny] <- hu[, Ny - 1]
  hv[, Ny] <- -hv[, Ny - 1]

  return(list(h = h, hu = hu, hv = hv))
}

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

cat("Starting simulation...\n")
cat("Progress: ")

start_time <- Sys.time()

# Storage for animation frames
save_interval <- max(1, floor(n_steps / 50))
frames <- list()
frame_times <- c()

for (step in 1:n_steps) {
  # Time integration step
  result <- lax_wendroff_step(h, hu, hv, dt, dx, dy)
  h <- result$h
  hu <- result$hu
  hv <- result$hv

  # Apply boundary conditions
  bc_result <- apply_boundary_conditions(h, hu, hv)
  h <- bc_result$h
  hu <- bc_result$hu
  hv <- bc_result$hv

  # Save frame
  if (step %% save_interval == 0 || step == n_steps) {
    frames[[length(frames) + 1]] <- h
    frame_times <- c(frame_times, step * dt)
  }

  # Progress indicator
  if (step %% 100 == 0) {
    progress <- round(100 * step / n_steps)
    cat(sprintf("%d%% ", progress))
  }
}

end_time <- Sys.time()
cat("\n\nSimulation completed!\n")
cat(sprintf(
  "Time elapsed: %.1f seconds\n",
  difftime(end_time, start_time, units = "secs")
))
cat(sprintf("Final h: min=%.4f, max=%.4f\n\n", min(h), max(h)))

# ============================================================================
# 2D VISUALIZATION (Top-down view)
# ============================================================================

cat("Creating 2D visualizations...\n")

# Convert to dataframe
h_df <- melt(h - h0) # Show deviation from mean
colnames(h_df) <- c("x", "y", "elevation")
h_df$x <- x[h_df$x]
h_df$y <- y[h_df$y]

# 2D heatmap of final state
p_2d <- ggplot(h_df, aes(x = x, y = y, fill = elevation)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(
    low = "#00008B",
    mid = "#4169E1",
    high = "#E0FFFF",
    midpoint = 0,
    name = "Surface\nElevation"
  ) +
  labs(
    title = "Pond Ripple - Water Surface Elevation",
    subtitle = sprintf("2D Shallow Water Equations | t = %.2f s", T_final),
    x = "x (m)",
    y = "y (m)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right",
    aspect.ratio = 1
  ) +
  coord_fixed(expand = FALSE)

print(p_2d)
ggsave("pond_ripple_2d.png", p_2d, width = 10, height = 10, dpi = 150)

# ============================================================================
# RAYSHADER 3D VISUALIZATION
# ============================================================================

cat("\nCreating rayshader 3D visualizations...\n")

# Prepare height matrix for rayshader (transpose for correct orientation)
h_matrix <- t(h)

# Scale the height for better visualization
height_scale <- 2.0
h_scaled <- (h_matrix - h0) * height_scale + h0

# ===== Create water color palette =====
water_colors <- colorRampPalette(c(
  "#001f3f", # Deep navy
  "#003366", # Dark blue
  "#0066cc", # Medium blue
  "#3399ff", # Light blue
  "#66ccff", # Sky blue
  "#99ddff", # Very light blue
  "#ccf2ff" # Almost white
))

# Create hillshade texture
h_texture <- h_scaled %>%
  height_shade(texture = water_colors(256)) %>%
  add_shadow(ray_shade(h_scaled, zscale = 0.05, lambert = TRUE), 0.5) %>%
  add_shadow(ambient_shade(h_scaled, zscale = 0.05), 0.3)

# ===== Static 3D render =====
cat("Rendering static 3D image...\n")

h_texture %>%
  plot_3d(
    h_scaled,
    zscale = 0.3,
    fov = 45,
    theta = 45,
    phi = 30,
    windowsize = c(1200, 900),
    zoom = 0.7,
    background = "#87CEEB",
    shadowcolor = "#1a3a5c",
    water = TRUE,
    waterdepth = h0 - 0.05,
    watercolor = "#3399ff",
    wateralpha = 0.6,
    waterlinecolor = "white",
    waterlinealpha = 0.3,
    solid = TRUE,
    soliddepth = "auto",
    solidcolor = "#8B4513",
    shadow = TRUE
  )

# Render high-quality snapshot
render_snapshot(
  "pond_ripple_3d_rayshader.png",
  clear = TRUE,
  title_text = "Pond Ripple - 3D Water Surface",
  title_size = 30,
  title_color = "white",
  title_bar_color = "#003366"
)

cat("3D render saved as 'pond_ripple_3d_rayshader.png'\n")

# ===== Alternative angle render =====
cat("Rendering alternative view...\n")

h_texture %>%
  plot_3d(
    h_scaled,
    zscale = 0.3,
    fov = 60,
    theta = -30,
    phi = 20,
    windowsize = c(1200, 900),
    zoom = 0.65,
    background = "#1a1a2e",
    shadowcolor = "#0a0a1a",
    water = TRUE,
    waterdepth = h0 - 0.08,
    watercolor = "#4169E1",
    wateralpha = 0.5,
    solid = TRUE,
    solidcolor = "#2c2c54",
    shadow = TRUE
  )

render_snapshot("pond_ripple_3d_angle2.png", clear = TRUE)
cat("Alternative view saved as 'pond_ripple_3d_angle2.png'\n")

# ===== Top-down 3D view =====
cat("Rendering top-down 3D view...\n")

h_texture %>%
  plot_3d(
    h_scaled,
    zscale = 0.4,
    fov = 0,
    theta = 0,
    phi = 89,
    windowsize = c(1000, 1000),
    zoom = 0.8,
    background = "#000033",
    water = TRUE,
    waterdepth = h0 - 0.1,
    watercolor = "#0055aa",
    wateralpha = 0.4,
    shadow = TRUE
  )

render_snapshot("pond_ripple_3d_topdown.png", clear = TRUE)
cat("Top-down view saved as 'pond_ripple_3d_topdown.png'\n")

# ============================================================================
# HIGH-QUALITY PHOTOREALISTIC RENDER
# ============================================================================

cat("\nCreating photorealistic render (this may take a moment)...\n")

# Enhanced water texture
h_photo <- h_scaled %>%
  height_shade(texture = water_colors(512)) %>%
  add_shadow(
    ray_shade(
      h_scaled,
      zscale = 0.05,
      sunaltitude = 30,
      sunangle = 315,
      lambert = TRUE,
      multicore = TRUE
    ),
    0.4
  ) %>%
  add_shadow(lamb_shade(h_scaled, zscale = 0.05, sunaltitude = 30), 0.3) %>%
  add_shadow(ambient_shade(h_scaled, zscale = 0.05), 0.2)

h_photo %>%
  plot_3d(
    h_scaled,
    zscale = 0.25,
    fov = 50,
    theta = 35,
    phi = 25,
    windowsize = c(1400, 1000),
    zoom = 0.6,
    background = "#87CEEB",
    shadowcolor = "#1a3a5c",
    water = TRUE,
    waterdepth = h0 - 0.06,
    watercolor = "#2266aa",
    wateralpha = 0.7,
    waterlinecolor = "#aaddff",
    waterlinealpha = 0.4,
    solid = TRUE,
    soliddepth = -0.5,
    solidcolor = "#654321",
    shadow = TRUE,
    shadowdepth = -0.6
  )

# Use render_highquality if available (requires rayrender)
if (requireNamespace("rayrender", quietly = TRUE)) {
  cat("Using rayrender for high-quality output...\n")
  render_highquality(
    filename = "pond_ripple_photorealistic.png",
    samples = 200,
    light = TRUE,
    lightdirection = c(315, 315),
    lightaltitude = c(30, 60),
    lightintensity = c(400, 200),
    lightcolor = c("#ffffee", "#aaccff"),
    interactive = FALSE,
    clear = TRUE,
    width = 1800,
    height = 1200
  )
  cat("Photorealistic render saved as 'pond_ripple_photorealistic.png'\n")
} else {
  render_snapshot("pond_ripple_photorealistic.png", clear = TRUE)
  cat("High-quality snapshot saved as 'pond_ripple_photorealistic.png'\n")
}

# ============================================================================
# EVOLUTION SEQUENCE
# ============================================================================

cat("\nCreating evolution sequence...\n")

# Select frames for evolution
n_display <- min(6, length(frames))
display_indices <- round(seq(1, length(frames), length.out = n_display))

# Create individual plots for each time step
evolution_plots <- list()
for (i in seq_along(display_indices)) {
  idx <- display_indices[i]
  h_frame <- frames[[idx]]

  h_frame_df <- melt(h_frame - h0)
  colnames(h_frame_df) <- c("xi", "yi", "elevation")
  h_frame_df$x <- x[h_frame_df$xi]
  h_frame_df$y <- y[h_frame_df$yi]

  p <- ggplot(h_frame_df, aes(x = x, y = y, fill = elevation)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient2(
      low = "#00008B",
      mid = "#4169E1",
      high = "#E0FFFF",
      midpoint = 0,
      limits = c(-0.2, 0.4),
      name = "η"
    ) +
    labs(title = sprintf("t = %.2f s", frame_times[idx])) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
      aspect.ratio = 1
    ) +
    coord_fixed(expand = FALSE)

  evolution_plots[[i]] <- p
}

# Combine plots
if (requireNamespace("gridExtra", quietly = TRUE)) {
  library(gridExtra)
  p_evolution <- grid.arrange(
    grobs = evolution_plots,
    nrow = 2,
    top = "Pond Ripple Evolution - Expanding Circular Waves"
  )
  ggsave(
    "pond_ripple_evolution.png",
    p_evolution,
    width = 14,
    height = 10,
    dpi = 150
  )
  cat("Evolution sequence saved as 'pond_ripple_evolution.png'\n")
}

# ============================================================================
# CROSS-SECTION PLOT
# ============================================================================

cat("\nCreating cross-section plot...\n")

# Extract cross-section through center
mid_y <- Ny / 2
cross_section_df <- data.frame()

for (i in seq_along(display_indices)) {
  idx <- display_indices[i]
  h_frame <- frames[[idx]]

  temp_df <- data.frame(
    x = x,
    elevation = h_frame[, mid_y] - h0,
    time = sprintf("t = %.2f s", frame_times[idx])
  )
  cross_section_df <- rbind(cross_section_df, temp_df)
}

p_cross <- ggplot(cross_section_df, aes(x = x, y = elevation, color = time)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_viridis_d(option = "plasma", name = "Time") +
  labs(
    title = "Cross-Section Through Center of Pond",
    subtitle = "Wave height profile at y = Ly/2",
    x = "x (m)",
    y = "Surface Elevation (m)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right"
  )

print(p_cross)
ggsave(
  "pond_ripple_crosssection.png",
  p_cross,
  width = 12,
  height = 6,
  dpi = 150
)

# ============================================================================
# VELOCITY FIELD VISUALIZATION
# ============================================================================

cat("\nCreating velocity field visualization...\n")

# Compute velocity components
u <- hu / pmax(h, 1e-10)
v <- hv / pmax(h, 1e-10)
vel_mag <- sqrt(u^2 + v^2)

# Downsample for arrows
skip <- 8
x_idx <- seq(1, Nx, by = skip)
y_idx <- seq(1, Ny, by = skip)

arrow_df <- expand.grid(xi = x_idx, yi = y_idx)
arrow_df$x <- x[arrow_df$xi]
arrow_df$y <- y[arrow_df$yi]
arrow_df$u <- u[cbind(arrow_df$xi, arrow_df$yi)]
arrow_df$v <- v[cbind(arrow_df$xi, arrow_df$yi)]
arrow_df$vel <- vel_mag[cbind(arrow_df$xi, arrow_df$yi)]

# Velocity magnitude background
vel_df <- melt(vel_mag)
colnames(vel_df) <- c("xi", "yi", "velocity")
vel_df$x <- x[vel_df$xi]
vel_df$y <- y[vel_df$yi]

p_vel <- ggplot() +
  geom_raster(
    data = vel_df,
    aes(x = x, y = y, fill = velocity),
    interpolate = TRUE
  ) +
  scale_fill_viridis_c(option = "turbo", name = "Velocity\n(m/s)") +
  geom_segment(
    data = arrow_df,
    aes(x = x, y = y, xend = x + u * 3, yend = y + v * 3),
    arrow = arrow(length = unit(0.1, "cm")),
    color = "white",
    alpha = 0.7,
    linewidth = 0.3
  ) +
  labs(
    title = "Velocity Field at Final Time",
    subtitle = "Arrow direction shows flow, color shows speed",
    x = "x (m)",
    y = "y (m)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    aspect.ratio = 1
  ) +
  coord_fixed(expand = FALSE)

print(p_vel)
ggsave("pond_ripple_velocity.png", p_vel, width = 10, height = 10, dpi = 150)

# ============================================================================
# 3D RAYSHADER WITH VELOCITY OVERLAY
# ============================================================================

cat("\nCreating final artistic 3D render...\n")

# Create a beautiful final render
artistic_colors <- colorRampPalette(c(
  "#0d1b2a", # Dark navy
  "#1b263b", # Navy
  "#415a77", # Steel blue
  "#778da9", # Light steel
  "#e0e1dd" # Off white
))

h_artistic <- h_scaled %>%
  height_shade(texture = artistic_colors(256)) %>%
  add_shadow(
    ray_shade(
      h_scaled,
      zscale = 0.05,
      sunaltitude = 25,
      sunangle = 290,
      lambert = TRUE
    ),
    0.5
  ) %>%
  add_shadow(ambient_shade(h_scaled, zscale = 0.05), 0.3)

h_artistic %>%
  plot_3d(
    h_scaled,
    zscale = 0.35,
    fov = 55,
    theta = 40,
    phi = 28,
    windowsize = c(1400, 1000),
    zoom = 0.55,
    background = "#f0f0f0",
    shadowcolor = "#3a3a5a",
    water = TRUE,
    waterdepth = h0 - 0.04,
    watercolor = "#4a6fa5",
    wateralpha = 0.65,
    waterlinecolor = "#e0e0ff",
    waterlinealpha = 0.5,
    solid = TRUE,
    soliddepth = "auto",
    solidcolor = "#5a4a3a",
    shadow = TRUE
  )

render_snapshot(
  "pond_ripple_artistic.png",
  clear = TRUE,
  title_text = "Shallow Water Waves",
  title_size = 28,
  title_color = "#1b263b",
  title_bar_color = "#e0e1dd",
  title_bar_alpha = 0.9
)

cat("Artistic render saved as 'pond_ripple_artistic.png'\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n============================================\n")
cat("SIMULATION SUMMARY\n")
cat("============================================\n")
cat("Physics: 2D Shallow Water Equations\n")
cat("Method: Two-step Lax-Wendroff scheme\n")
cat("Boundaries: Reflective (solid walls)\n")
cat(sprintf("Grid: %d x %d cells\n", Nx, Ny))
cat(sprintf("Domain: %.1f x %.1f meters\n", Lx, Ly))
cat(sprintf("Simulation time: %.2f seconds\n", T_final))
cat(sprintf("Time steps: %d\n", n_steps))
cat("\nInitial Condition:\n")
cat(sprintf("  Base depth: %.2f m\n", h0))
cat(sprintf("  Drop amplitude: %.2f m\n", drop_amplitude))
cat(sprintf("  Drop width (σ): %.2f m\n", drop_sigma))
cat("\nOutput Files:\n")
cat("  1. pond_ripple_2d.png - Top-down 2D view\n")
cat("  2. pond_ripple_3d_rayshader.png - Main 3D render\n")
cat("  3. pond_ripple_3d_angle2.png - Alternative 3D angle\n")
cat("  4. pond_ripple_3d_topdown.png - Top-down 3D view\n")
cat("  5. pond_ripple_photorealistic.png - High-quality render\n")
cat("  6. pond_ripple_evolution.png - Time evolution sequence\n")
cat("  7. pond_ripple_crosssection.png - Wave profile cross-section\n")
cat("  8. pond_ripple_velocity.png - Velocity field\n")
cat("  9. pond_ripple_artistic.png - Artistic 3D render\n")
cat("\nPhysics Notes:\n")
cat("  - Circular waves expand from initial drop\n")
cat("  - Wave speed ≈ √(gh) (shallow water approximation)\n")
cat("  - Reflections occur at boundaries\n")
cat("  - Energy conserved (no dissipation in model)\n")
cat("============================================\n")
