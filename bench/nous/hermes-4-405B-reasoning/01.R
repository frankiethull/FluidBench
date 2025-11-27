# 1D Linear Advection Simulation with Space-Time Heatmap
library(ggplot2)

# Physics Parameters
c <- 1 # Wave speed
x_min <- 0 # Domain start
x_max <- 2 # Domain end
t_max <- 2 # Simulation time

# Numerical Parameters
nx <- 200 # Spatial grid points
nt <- 400 # Temporal grid points
dx <- (x_max - x_min) / (nx - 1)
dt <- t_max / nt

# Stability check (CFL condition)
cfl <- c * dt / dx
if (cfl > 1.0) {
  stop("CFL condition violated - decrease dt or increase dx")
}

# Spatial grid
x <- seq(x_min, x_max, length.out = nx)

# Initial condition (Gaussian pulse)
u0 <- exp(-(x - 0.5)^2 / (2 * (0.1)^2))

# Initialize solution matrix
u <- matrix(0, nrow = nx, ncol = nt)
u[, 1] <- u0

# Time-stepping with periodic BCs
for (n in 1:(nt - 1)) {
  # Interior points
  for (i in 2:nx) {
    u[i, n + 1] <- u[i, n] - c * dt / dx * (u[i, n] - u[i - 1, n])
  }

  # Periodic boundary (left boundary)
  u[1, n + 1] <- u[1, n] - c * dt / dx * (u[1, n] - u[nx, n])
}

# Create time vector
time_steps <- seq(0, t_max, length.out = nt)

# Reshape data for ggplot
df <- data.frame(
  x = rep(x, times = nt),
  t = rep(time_steps, each = nx),
  u = as.vector(u)
)

# Create space-time heatmap
ggplot(df, aes(x = x, y = t, fill = u)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(option = "A", direction = -1) +
  labs(
    title = "1D Advection Space-Time Diagram",
    x = "Space (x)",
    y = "Time (t)",
    fill = "Amplitude"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())
