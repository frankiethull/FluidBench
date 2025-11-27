# 1D Advection Simulation: Space-Time Waterfall Visualization
# Simulating ∂u/∂t + c(∂u/∂x) = 0 with periodic boundary conditions

# Load required libraries
library(ggplot2)
library(reshape2)

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================

# Spatial domain
L <- 2.0 # Domain length [0, L]
Nx <- 200 # Number of spatial grid points
dx <- L / Nx # Spatial step size
x <- seq(0, L - dx, by = dx) # Spatial grid (periodic, so exclude endpoint)

# Temporal domain
c <- 1.0 # Advection velocity
CFL <- 0.8 # Courant-Friedrichs-Lewy number (stability condition: CFL <= 1)
dt <- CFL * dx / abs(c) # Time step (derived from CFL condition)
T_final <- 4.0 # Total simulation time
Nt <- ceiling(T_final / dt) # Number of time steps
t <- seq(0, (Nt - 1) * dt, by = dt) # Time grid

cat("============================================\n")
cat("1D ADVECTION SIMULATION PARAMETERS\n")
cat("============================================\n")
cat(sprintf("Domain: x ∈ [0, %.1f]\n", L))
cat(sprintf("Grid points (Nx): %d\n", Nx))
cat(sprintf("Spatial step (dx): %.6f\n", dx))
cat(sprintf("Advection velocity (c): %.1f\n", c))
cat(sprintf("CFL number: %.2f\n", CFL))
cat(sprintf("Time step (dt): %.6f\n", dt))
cat(sprintf("Total time: %.2f\n", T_final))
cat(sprintf("Time steps (Nt): %d\n", Nt))
cat("============================================\n\n")

# ============================================================================
# INITIAL CONDITION: Gaussian Pulse
# ============================================================================

# Gaussian pulse parameters
x0 <- 0.5 # Initial center of the pulse
sigma <- 0.1 # Width of the Gaussian

# Initial condition: u(x, 0) = exp(-((x - x0)^2) / (2 * sigma^2))
u_initial <- exp(-((x - x0)^2) / (2 * sigma^2))

# ============================================================================
# STORAGE MATRIX FOR SPACE-TIME DATA
# ============================================================================

# Create matrix to store solution at all time steps
# Rows = time steps, Columns = spatial points
u_history <- matrix(0, nrow = Nt, ncol = Nx)

# Store initial condition
u_history[1, ] <- u_initial

# ============================================================================
# NUMERICAL SCHEME: Upwind Finite Difference Method
# ============================================================================

# For c > 0, use backward difference (upwind scheme):
# u_j^{n+1} = u_j^n - c * (dt/dx) * (u_j^n - u_{j-1}^n)
#
# This is first-order accurate and stable for CFL <= 1

cat("Running simulation using Upwind scheme...\n")

# Current solution vector
u <- u_initial

# Time-stepping loop
for (n in 2:Nt) {
  # Store previous state
  u_old <- u

  # Apply upwind scheme with periodic boundaries
  if (c > 0) {
    # Backward difference for positive velocity
    # j-1 index with periodic boundary (wrap around)
    u_left <- c(u_old[Nx], u_old[1:(Nx - 1)]) # Shifted array for u_{j-1}
    u <- u_old - c * (dt / dx) * (u_old - u_left)
  } else {
    # Forward difference for negative velocity
    # j+1 index with periodic boundary
    u_right <- c(u_old[2:Nx], u_old[1]) # Shifted array for u_{j+1}
    u <- u_old - c * (dt / dx) * (u_right - u_old)
  }

  # Store current state
  u_history[n, ] <- u

  # Progress indicator
  if (n %% 100 == 0) {
    cat(sprintf("  Time step %d / %d (t = %.3f)\n", n, Nt, t[n]))
  }
}

cat("Simulation complete!\n\n")

# ============================================================================
# DATA RESHAPING FOR VISUALIZATION
# ============================================================================

cat("Reshaping data for visualization...\n")

# Convert matrix to dataframe for ggplot2
# We'll downsample if there are too many time steps for smooth plotting
downsample_factor <- max(1, floor(Nt / 200)) # Keep ~200 time slices max

time_indices <- seq(1, Nt, by = downsample_factor)
u_downsampled <- u_history[time_indices, ]
t_downsampled <- t[time_indices]

# Create a dataframe in long format
spacetime_df <- expand.grid(
  x = x,
  t = t_downsampled
)

# Add amplitude values (need to transpose and vectorize properly)
spacetime_df$u <- as.vector(t(u_downsampled))

cat(sprintf("Data reshaped: %d rows for plotting\n", nrow(spacetime_df)))
cat(sprintf(
  "Time slices: %d (downsampled from %d)\n\n",
  length(t_downsampled),
  Nt
))

# ============================================================================
# VISUALIZATION: Space-Time Heatmap
# ============================================================================

cat("Creating Space-Time Heatmap...\n")

# Main space-time plot
p1 <- ggplot(spacetime_df, aes(x = x, y = t, fill = u)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "Amplitude\n(u)",
    limits = c(0, 1)
  ) +
  labs(
    title = "1D Advection: Space-Time Waterfall Diagram",
    subtitle = expression(paste(
      "Equation: ",
      frac(partialdiff * u, partialdiff * t) +
        c * frac(partialdiff * u, partialdiff * x) ==
        0,
      "   |   c = 1.0, Periodic BCs, Gaussian Initial Condition"
    )),
    x = "Space (x)",
    y = "Time (t)",
    caption = "Diagonal stripes show wave propagation | Upwind finite difference scheme"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
    plot.caption = element_text(size = 9, color = "gray50"),
    legend.position = "right",
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold")
  ) +
  coord_cartesian(expand = FALSE)

# Display the plot
print(p1)

# Save the plot
ggsave("advection_spacetime.png", p1, width = 10, height = 8, dpi = 150)
cat("Plot saved as 'advection_spacetime.png'\n\n")

# ============================================================================
# ADDITIONAL VISUALIZATION: Alternative Color Scheme
# ============================================================================

# Create an alternative visualization with different aesthetics
p2 <- ggplot(spacetime_df, aes(x = x, y = t, fill = u)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("navy", "blue", "cyan", "yellow", "red"),
    name = "Amplitude\n(u)",
    limits = c(0, 1)
  ) +
  labs(
    title = "Space-Time Evolution of Gaussian Pulse",
    subtitle = "Wave travels with velocity c = 1 (period = 2 due to periodic boundaries)",
    x = "Position (x)",
    y = "Time (t)"
  ) +
  theme_dark(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5,
      color = "white"
    ),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray80"),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  coord_cartesian(expand = FALSE)

print(p2)

# Save alternative plot
ggsave("advection_spacetime_alt.png", p2, width = 10, height = 8, dpi = 150)
cat("Alternative plot saved as 'advection_spacetime_alt.png'\n\n")

# ============================================================================
# SNAPSHOT VISUALIZATION: Wave at Different Times
# ============================================================================

# Select specific time snapshots
snapshot_times <- c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5)
snapshot_indices <- sapply(snapshot_times, function(ts) which.min(abs(t - ts)))

# Create snapshot dataframe
snapshot_df <- data.frame()
for (i in seq_along(snapshot_times)) {
  idx <- snapshot_indices[i]
  temp_df <- data.frame(
    x = x,
    u = u_history[idx, ],
    time_label = sprintf("t = %.1f", t[idx])
  )
  snapshot_df <- rbind(snapshot_df, temp_df)
}

# Snapshot plot
p3 <- ggplot(snapshot_df, aes(x = x, y = u, color = time_label)) +
  geom_line(linewidth = 1) +
  scale_color_viridis_d(option = "turbo", name = "Time") +
  labs(
    title = "Wave Snapshots at Different Times",
    subtitle = "Gaussian pulse traveling with periodic boundary conditions",
    x = "Position (x)",
    y = "Amplitude (u)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right"
  ) +
  ylim(0, 1.1)

print(p3)

# Save snapshot plot
ggsave("advection_snapshots.png", p3, width = 10, height = 6, dpi = 150)
cat("Snapshot plot saved as 'advection_snapshots.png'\n\n")

# ============================================================================
# VERIFICATION: Analytical Solution Comparison
# ============================================================================

cat("============================================\n")
cat("VERIFICATION: Comparing with Analytical Solution\n")
cat("============================================\n")

# Analytical solution: u(x,t) = u0(x - ct) with periodic boundaries
analytical_solution <- function(x, t, x0, sigma, c, L) {
  # Position shifted by advection (with periodic wrap)
  x_shifted <- (x - c * t) %% L
  # Handle negative modulo
  x_shifted[x_shifted < 0] <- x_shifted[x_shifted < 0] + L
  # Gaussian at shifted position
  exp(-((x_shifted - x0)^2) / (2 * sigma^2))
}

# Compare at final time
u_numerical <- u_history[Nt, ]
u_analytical <- analytical_solution(x, t[Nt], x0, sigma, c, L)

# Calculate error
L2_error <- sqrt(sum((u_numerical - u_analytical)^2) * dx)
Linf_error <- max(abs(u_numerical - u_analytical))

cat(sprintf("Final time: t = %.3f\n", t[Nt]))
cat(sprintf("L2 Error: %.6e\n", L2_error))
cat(sprintf("L∞ Error: %.6e\n", Linf_error))
cat("(Note: Upwind scheme has numerical diffusion - this is expected)\n")
cat("============================================\n\n")

# Comparison plot
comparison_df <- data.frame(
  x = rep(x, 2),
  u = c(u_numerical, u_analytical),
  type = rep(c("Numerical", "Analytical"), each = Nx)
)

p4 <- ggplot(comparison_df, aes(x = x, y = u, color = type, linetype = type)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Analytical" = "blue", "Numerical" = "red")) +
  scale_linetype_manual(
    values = c("Analytical" = "solid", "Numerical" = "dashed")
  ) +
  labs(
    title = sprintf("Numerical vs Analytical Solution at t = %.2f", t[Nt]),
    subtitle = sprintf(
      "L2 Error = %.2e | Note: Upwind scheme introduces numerical diffusion",
      L2_error
    ),
    x = "Position (x)",
    y = "Amplitude (u)",
    color = "Solution",
    linetype = "Solution"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "top"
  )

print(p4)

ggsave("advection_comparison.png", p4, width = 10, height = 6, dpi = 150)
cat("Comparison plot saved as 'advection_comparison.png'\n\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("============================================\n")
cat("SIMULATION SUMMARY\n")
cat("============================================\n")
cat("Physics: 1D Linear Advection Equation\n")
cat("Method: First-order Upwind Finite Difference\n")
cat("Boundary: Periodic (wave wraps around domain)\n")
cat("Initial Condition: Gaussian pulse\n")
cat("\nOutput Files:\n")
cat("  1. advection_spacetime.png - Main space-time heatmap\n")
cat("  2. advection_spacetime_alt.png - Alternative color scheme\n")
cat("  3. advection_snapshots.png - Wave profiles at different times\n")
cat("  4. advection_comparison.png - Numerical vs analytical solution\n")
cat("\nInterpretation:\n")
cat("  - Diagonal stripes in heatmap show wave propagation\n")
cat("  - Slope of stripes = 1/c (inverse of wave speed)\n")
cat("  - Wave completes one period through domain in t = L/c = 2\n")
cat("============================================\n")
