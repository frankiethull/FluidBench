# Load necessary libraries
library(ggplot2)
library(dplyr)

# Parameters
c <- 1 # Advection velocity
L <- 2 # Domain length (x ∈ [0, L])
T <- 2 # Total simulation time
Nx <- 200 # Number of spatial grid points
Nt <- 200 # Number of time steps
dx <- L / (Nx - 1) # Spatial step size
dt <- T / (Nt - 1) # Time step size

# Ensure stability (CFL condition: c*dt/dx <= 1)
if (c * dt / dx > 1) {
  stop("CFL condition violated. Reduce dt or increase dx.")
}

# Grid: x (space) and t (time)
x_grid <- seq(0, L, length.out = Nx)
t_grid <- seq(0, T, length.out = Nt)

# Initialize matrix to store u(x,t) at all times
# Dimension: Nt rows (time) × Nx columns (space)
u_matrix <- matrix(0, nrow = Nt, ncol = Nx)

# Initial condition (Gaussian pulse centered at x = 1)
u_matrix[1, ] <- exp(-((x_grid - 1) / 0.1)^2)

# 1D Advection Solver (Finite Difference Upwind Method)
for (t_idx in 1:(Nt - 1)) {
  for (x_idx in 1:Nx) {
    # Periodic boundary handling
    x_prev <- ifelse(x_idx == 1, Nx, x_idx - 1) # Wrap around for x=0
    u_matrix[t_idx + 1, x_idx] <- u_matrix[t_idx, x_idx] -
      (c * dt / dx) * (u_matrix[t_idx, x_idx] - u_matrix[t_idx, x_prev])
  }
}

# Reshape data for ggplot (long format: columns for x, t, u)
df_plot <- expand.grid(x = x_grid, t = t_grid) %>%
  mutate(u = as.vector(u_matrix))

# Space-Time Heatmap (ggplot2)
ggplot(df_plot, aes(x = x, y = t, fill = u)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(option = "plasma") +
  labs(
    title = "Space-Time Heatmap of 1D Advection (Gaussian Pulse)",
    subtitle = "∂u/∂t + c ∂u/∂x = 0, c = 1, Periodic Boundaries",
    x = "Space (x)",
    y = "Time (t)",
    fill = "Amplitude (u)"
  ) +
  theme_minimal() +
  theme(legend.position = "right")
