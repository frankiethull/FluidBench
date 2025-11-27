# Gray-Scott Reaction-Diffusion Model
# Simulating Turing Pattern Formation (Coral/Fingerprint Patterns)
#
# Equations:
# ∂U/∂t = Du∇²U - UV² + F(1-U)
# ∂V/∂t = Dv∇²V + UV² - (F+k)V

# Load required libraries
library(ggplot2)
library(reshape2)

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================

cat("============================================\n")
cat("GRAY-SCOTT REACTION-DIFFUSION MODEL\n")
cat("Turing Pattern Formation\n")
cat("============================================\n\n")

# Grid parameters
Nx <- 200 # Grid width (increased for better resolution)
Ny <- 200 # Grid height
dx <- 1.0 # Spatial step

# Model parameters (Coral preset)
Du <- 0.16 # Diffusion rate of U
Dv <- 0.08 # Diffusion rate of V
F <- 0.035 # Feed rate
k <- 0.060 # Kill rate

# Time stepping
dt <- 1.0 # Time step
n_iterations <- 15000 # Total iterations (increased for pattern development)
save_interval <- 500 # Save visualization every N steps

cat("Grid: ", Nx, " x ", Ny, "\n")
cat("Parameters (Coral preset):\n")
cat(sprintf("  Du = %.3f (U diffusion)\n", Du))
cat(sprintf("  Dv = %.3f (V diffusion)\n", Dv))
cat(sprintf("  F  = %.3f (feed rate)\n", F))
cat(sprintf("  k  = %.3f (kill rate)\n", k))
cat(sprintf("Iterations: %d\n", n_iterations))
cat(sprintf("Time step: %.2f\n\n", dt))

# ============================================================================
# INITIALIZE CONCENTRATION FIELDS
# ============================================================================

cat("Initializing concentration fields...\n")

# Initialize U = 1, V = 0 everywhere
U <- matrix(1.0, nrow = Nx, ncol = Ny)
V <- matrix(0.0, nrow = Nx, ncol = Ny)

# Add perturbation in the center: 20x20 square with noise
set.seed(42)
center_x <- Nx / 2
center_y <- Ny / 2
perturb_size <- 20

# Define perturbation region
x_start <- floor(center_x - perturb_size / 2)
x_end <- floor(center_x + perturb_size / 2)
y_start <- floor(center_y - perturb_size / 2)
y_end <- floor(center_y + perturb_size / 2)

# Set perturbation with random noise
U[x_start:x_end, y_start:y_end] <- 0.5 +
  0.1 *
    matrix(
      runif((x_end - x_start + 1) * (y_end - y_start + 1)),
      nrow = x_end - x_start + 1
    )
V[x_start:x_end, y_start:y_end] <- 0.25 +
  0.1 *
    matrix(
      runif((x_end - x_start + 1) * (y_end - y_start + 1)),
      nrow = x_end - x_start + 1
    )

# Add additional small random seeds to create more varied patterns
n_seeds <- 5
for (i in 1:n_seeds) {
  seed_x <- sample(20:(Nx - 20), 1)
  seed_y <- sample(20:(Ny - 20), 1)
  seed_size <- sample(5:10, 1)

  sx <- max(1, seed_x - seed_size):min(Nx, seed_x + seed_size)
  sy <- max(1, seed_y - seed_size):min(Ny, seed_y + seed_size)

  U[sx, sy] <- 0.5 +
    0.1 * matrix(runif(length(sx) * length(sy)), nrow = length(sx))
  V[sx, sy] <- 0.25 +
    0.1 * matrix(runif(length(sx) * length(sy)), nrow = length(sx))
}

cat("Initial perturbation added in center and", n_seeds, "random seeds\n\n")

# ============================================================================
# LAPLACIAN OPERATOR (5-POINT STENCIL WITH PERIODIC BC)
# ============================================================================

compute_laplacian <- function(C) {
  # 5-point stencil Laplacian with periodic boundary conditions
  # ∇²C = (C[i+1,j] + C[i-1,j] + C[i,j+1] + C[i,j-1] - 4*C[i,j]) / dx²

  # Shift operations with periodic boundaries
  C_right <- cbind(C[, 2:Ny], C[, 1]) # C[i, j+1]
  C_left <- cbind(C[, Ny], C[, 1:(Ny - 1)]) # C[i, j-1]
  C_up <- rbind(C[2:Nx, ], C[1, ]) # C[i+1, j]
  C_down <- rbind(C[Nx, ], C[1:(Nx - 1), ]) # C[i-1, j]

  laplacian <- (C_right + C_left + C_up + C_down - 4 * C) / (dx^2)

  return(laplacian)
}

# ============================================================================
# ALTERNATIVE: 9-POINT STENCIL FOR SMOOTHER PATTERNS
# ============================================================================

compute_laplacian_9pt <- function(C) {
  # 9-point stencil for better isotropy
  # Weights: corners = 0.25, edges = 0.5, center = -3

  # Shift operations
  C_N <- rbind(C[2:Nx, ], C[1, ])
  C_S <- rbind(C[Nx, ], C[1:(Nx - 1), ])
  C_E <- cbind(C[, 2:Ny], C[, 1])
  C_W <- cbind(C[, Ny], C[, 1:(Ny - 1)])

  C_NE <- rbind(cbind(C[2:Nx, 2:Ny], C[2:Nx, 1]), cbind(C[1, 2:Ny], C[1, 1]))
  C_NW <- rbind(
    cbind(C[2:Nx, Ny], C[2:Nx, 1:(Ny - 1)]),
    cbind(C[1, Ny], C[1, 1:(Ny - 1)])
  )
  C_SE <- rbind(
    cbind(C[Nx, 2:Ny], C[Nx, 1]),
    cbind(C[1:(Nx - 1), 2:Ny], C[1:(Nx - 1), 1])
  )
  C_SW <- rbind(
    cbind(C[Nx, Ny], C[Nx, 1:(Ny - 1)]),
    cbind(C[1:(Nx - 1), Ny], C[1:(Nx - 1), 1:(Ny - 1)])
  )

  laplacian <- (0.5 *
    (C_N + C_S + C_E + C_W) +
    0.25 * (C_NE + C_NW + C_SE + C_SW) -
    3 * C) /
    (dx^2)

  return(laplacian)
}

# ============================================================================
# VISUALIZATION FUNCTION
# ============================================================================

create_pattern_plot <- function(V, iteration, title_suffix = "") {
  # Convert matrix to dataframe for ggplot
  V_df <- melt(V)
  colnames(V_df) <- c("x", "y", "concentration")

  # Create the plot
  p <- ggplot(V_df, aes(x = x, y = y, fill = concentration)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_viridis_c(
      option = "viridis",
      name = "V\nConcentration",
      limits = c(0, max(0.4, max(V)))
    ) +
    labs(
      title = paste0(
        "Gray-Scott Reaction-Diffusion: Coral Pattern",
        title_suffix
      ),
      subtitle = sprintf(
        "Iteration: %d | F = %.3f, k = %.3f, Du = %.2f, Dv = %.2f",
        iteration,
        F,
        k,
        Du,
        Dv
      ),
      x = "x",
      y = "y"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
      legend.position = "right",
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      aspect.ratio = 1
    ) +
    coord_fixed(expand = FALSE)

  return(p)
}

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

cat("Starting simulation...\n")
cat("Progress: ")

start_time <- Sys.time()

# Storage for snapshots
snapshots <- list()
snapshot_times <- c()

for (iter in 1:n_iterations) {
  # Compute Laplacians
  lap_U <- compute_laplacian(U)
  lap_V <- compute_laplacian(V)

  # Reaction term
  UVV <- U * V * V

  # Update equations (forward Euler)
  # dU/dt = Du*∇²U - UV² + F(1-U)
  # dV/dt = Dv*∇²V + UV² - (F+k)V

  dU <- Du * lap_U - UVV + F * (1 - U)
  dV <- Dv * lap_V + UVV - (F + k) * V

  U <- U + dt * dU
  V <- V + dt * dV

  # Clamp values to prevent instability
  U <- pmax(0, pmin(1, U))
  V <- pmax(0, pmin(1, V))

  # Save snapshots at intervals
  if (iter %% save_interval == 0 || iter == n_iterations) {
    snapshots[[length(snapshots) + 1]] <- V
    snapshot_times <- c(snapshot_times, iter)
    cat(sprintf("%d ", iter))
  }

  # Progress indicator
  if (iter %% 1000 == 0 && !(iter %% save_interval == 0)) {
    cat(".")
  }
}

end_time <- Sys.time()
cat("\n\nSimulation completed!\n")
cat(sprintf(
  "Time elapsed: %.1f seconds\n\n",
  difftime(end_time, start_time, units = "secs")
))

# ============================================================================
# FINAL VISUALIZATION
# ============================================================================

cat("Creating visualizations...\n")

# Final pattern plot
p_final <- create_pattern_plot(V, n_iterations)
print(p_final)
ggsave("turing_pattern_coral.png", p_final, width = 10, height = 10, dpi = 200)
cat("Final pattern saved as 'turing_pattern_coral.png'\n")

# ============================================================================
# ALTERNATIVE COLOR SCHEMES
# ============================================================================

# Plasma color scheme (more organic look)
V_df <- melt(V)
colnames(V_df) <- c("x", "y", "concentration")

p_plasma <- ggplot(V_df, aes(x = x, y = y, fill = concentration)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "V",
    limits = c(0, max(0.4, max(V)))
  ) +
  labs(
    title = "Turing Pattern - Plasma Colormap",
    subtitle = sprintf("Gray-Scott Model | F = %.3f, k = %.3f", F, k),
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right",
    aspect.ratio = 1
  )

print(p_plasma)
ggsave(
  "turing_pattern_plasma.png",
  p_plasma,
  width = 10,
  height = 10,
  dpi = 200
)

# Magma color scheme (biological tissue look)
p_magma <- ggplot(V_df, aes(x = x, y = y, fill = concentration)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    option = "magma",
    name = "V",
    limits = c(0, max(0.4, max(V)))
  ) +
  labs(
    title = "Turing Pattern - Organic Texture",
    subtitle = "Self-organizing pattern from reaction-diffusion",
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right",
    aspect.ratio = 1
  )

print(p_magma)
ggsave("turing_pattern_magma.png", p_magma, width = 10, height = 10, dpi = 200)

# ============================================================================
# EVOLUTION SEQUENCE PLOT
# ============================================================================

cat("\nCreating evolution sequence...\n")

# Select frames for evolution display
n_display <- min(6, length(snapshots))
display_indices <- round(seq(1, length(snapshots), length.out = n_display))

# Combine snapshots into single dataframe
evolution_df <- data.frame()
for (i in seq_along(display_indices)) {
  idx <- display_indices[i]
  V_temp <- snapshots[[idx]]
  temp_df <- melt(V_temp)
  colnames(temp_df) <- c("x", "y", "concentration")
  temp_df$iteration <- sprintf("t = %d", snapshot_times[idx])
  evolution_df <- rbind(evolution_df, temp_df)
}

# Create faceted evolution plot
evolution_df$iteration <- factor(
  evolution_df$iteration,
  levels = unique(evolution_df$iteration)
)

p_evolution <- ggplot(evolution_df, aes(x = x, y = y, fill = concentration)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(option = "viridis", name = "V") +
  facet_wrap(~iteration, nrow = 2) +
  labs(
    title = "Pattern Evolution: From Initial Seed to Coral Structure",
    subtitle = sprintf(
      "Gray-Scott Model | Du=%.2f, Dv=%.2f, F=%.3f, k=%.3f",
      Du,
      Dv,
      F,
      k
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    aspect.ratio = 1
  ) +
  coord_fixed(expand = FALSE)

print(p_evolution)
ggsave(
  "turing_pattern_evolution.png",
  p_evolution,
  width = 14,
  height = 10,
  dpi = 150
)
cat("Evolution sequence saved as 'turing_pattern_evolution.png'\n")

# ============================================================================
# BOTH CHEMICALS SIDE BY SIDE
# ============================================================================

# Create combined U and V visualization
U_df <- melt(U)
colnames(U_df) <- c("x", "y", "concentration")
U_df$chemical <- "U (Activator substrate)"

V_df <- melt(V)
colnames(V_df) <- c("x", "y", "concentration")
V_df$chemical <- "V (Activator)"

combined_df <- rbind(U_df, V_df)

p_both <- ggplot(combined_df, aes(x = x, y = y, fill = concentration)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(option = "viridis", name = "Concentration") +
  facet_wrap(~chemical) +
  labs(
    title = "Both Chemical Species in Gray-Scott Model",
    subtitle = "U depleted where V is high (competitive dynamics)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    aspect.ratio = 1
  ) +
  coord_fixed(expand = FALSE)

print(p_both)
ggsave(
  "turing_pattern_both_chemicals.png",
  p_both,
  width = 14,
  height = 7,
  dpi = 150
)

# ============================================================================
# DIFFERENT PARAMETER PRESETS
# ============================================================================

cat("\nSimulating additional pattern presets...\n")

# Function to run simulation with different parameters
run_grayscott <- function(F_val, k_val, n_iter = 10000, label = "") {
  # Initialize
  U_new <- matrix(1.0, nrow = Nx, ncol = Ny)
  V_new <- matrix(0.0, nrow = Nx, ncol = Ny)

  # Central perturbation
  U_new[x_start:x_end, y_start:y_end] <- 0.5 +
    0.1 *
      matrix(
        runif((x_end - x_start + 1) * (y_end - y_start + 1)),
        nrow = x_end - x_start + 1
      )
  V_new[x_start:x_end, y_start:y_end] <- 0.25 +
    0.1 *
      matrix(
        runif((x_end - x_start + 1) * (y_end - y_start + 1)),
        nrow = x_end - x_start + 1
      )

  # Run simulation
  for (iter in 1:n_iter) {
    lap_U <- compute_laplacian(U_new)
    lap_V <- compute_laplacian(V_new)

    UVV <- U_new * V_new * V_new

    dU <- Du * lap_U - UVV + F_val * (1 - U_new)
    dV <- Dv * lap_V + UVV - (F_val + k_val) * V_new

    U_new <- U_new + dt * dU
    V_new <- V_new + dt * dV

    U_new <- pmax(0, pmin(1, U_new))
    V_new <- pmax(0, pmin(1, V_new))
  }

  return(V_new)
}

# Different presets
presets <- list(
  list(F = 0.035, k = 0.060, name = "Coral"),
  list(F = 0.030, k = 0.055, name = "Mitosis"),
  list(F = 0.025, k = 0.060, name = "Maze"),
  list(F = 0.040, k = 0.060, name = "Spots")
)

# Run each preset and collect results
preset_results <- list()
for (preset in presets) {
  cat(sprintf(
    "  Running '%s' preset (F=%.3f, k=%.3f)...\n",
    preset$name,
    preset$F,
    preset$k
  ))
  V_result <- run_grayscott(preset$F, preset$k, n_iter = 12000)

  V_temp_df <- melt(V_result)
  colnames(V_temp_df) <- c("x", "y", "concentration")
  V_temp_df$preset <- sprintf(
    "%s\n(F=%.3f, k=%.3f)",
    preset$name,
    preset$F,
    preset$k
  )

  preset_results[[preset$name]] <- V_temp_df
}

# Combine all presets
all_presets_df <- do.call(rbind, preset_results)

# Create comparison plot
p_presets <- ggplot(all_presets_df, aes(x = x, y = y, fill = concentration)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(option = "viridis", name = "V") +
  facet_wrap(~preset, nrow = 2) +
  labs(
    title = "Gray-Scott Parameter Space: Different Turing Patterns",
    subtitle = "Same equations, different parameters → different self-organized structures"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 10),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    aspect.ratio = 1
  ) +
  coord_fixed(expand = FALSE)

print(p_presets)
ggsave(
  "turing_pattern_presets.png",
  p_presets,
  width = 12,
  height = 12,
  dpi = 150
)
cat("Parameter presets saved as 'turing_pattern_presets.png'\n")

# ============================================================================
# HIGH-RESOLUTION ARTISTIC RENDER
# ============================================================================

cat("\nCreating high-resolution artistic render...\n")

# Artistic version with custom color palette
V_df <- melt(V)
colnames(V_df) <- c("x", "y", "concentration")

# Custom biologically-inspired palette
bio_colors <- c(
  "#0a0a23",
  "#1a1a4e",
  "#2d4a6e",
  "#4a8c7f",
  "#7eb87e",
  "#b8d87a",
  "#f0e68c",
  "#ffffff"
)

p_artistic <- ggplot(V_df, aes(x = x, y = y, fill = concentration)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradientn(
    colors = bio_colors,
    name = "Chemical V",
    limits = c(0, max(V))
  ) +
  labs(
    title = "TURING PATTERN",
    subtitle = "Self-Organization in Reaction-Diffusion Systems"
  ) +
  theme_void(base_size = 16) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 24,
      hjust = 0.5,
      color = "white",
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray70"),
    plot.background = element_rect(fill = "black", color = NA),
    legend.position = "none",
    aspect.ratio = 1,
    plot.margin = margin(20, 20, 20, 20)
  )

print(p_artistic)
ggsave(
  "turing_pattern_artistic.png",
  p_artistic,
  width = 12,
  height = 12,
  dpi = 300
)
cat("Artistic render saved as 'turing_pattern_artistic.png'\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n============================================\n")
cat("SIMULATION SUMMARY\n")
cat("============================================\n")
cat("Model: Gray-Scott Reaction-Diffusion\n")
cat("Pattern Type: Coral/Fingerprint (Turing Pattern)\n")
cat(sprintf("Grid Size: %d x %d\n", Nx, Ny))
cat(sprintf("Iterations: %d\n", n_iterations))
cat("\nFinal State Statistics:\n")
cat(sprintf("  U: min=%.4f, max=%.4f, mean=%.4f\n", min(U), max(U), mean(U)))
cat(sprintf("  V: min=%.4f, max=%.4f, mean=%.4f\n", min(V), max(V), mean(V)))
cat("\nOutput Files:\n")
cat("  1. turing_pattern_coral.png - Main viridis colormap\n")
cat("  2. turing_pattern_plasma.png - Plasma colormap\n")
cat("  3. turing_pattern_magma.png - Magma colormap\n")
cat("  4. turing_pattern_evolution.png - Time evolution sequence\n")
cat("  5. turing_pattern_both_chemicals.png - U and V comparison\n")
cat("  6. turing_pattern_presets.png - Different parameter presets\n")
cat("  7. turing_pattern_artistic.png - High-res artistic render\n")
cat("\nPhysics Notes:\n")
cat("  - Turing instability: diffusion-driven pattern formation\n")
cat("  - U acts as substrate, V as autocatalytic activator\n")
cat("  - Different F,k values → spots, stripes, mazes, etc.\n")
cat("  - Biological relevance: animal coat patterns, coral growth\n")
cat("============================================\n")
