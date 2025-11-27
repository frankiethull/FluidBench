# Space-Time Waterfall: 1D Advection Simulation with ggplot2 Heatmap

library(ggplot2)
library(reshape2) # for melt() if needed, though we'll use tidyr style manually
library(dplyr)

# Physics parameters
L <- 2.0 # Domain length [0, L]
c <- 1.0 # Advection velocity
T_final <- 2.0 # Total simulation time (so wave travels full domain)

# Numerical parameters
nx <- 200 # Number of spatial points
nt <- 300 # Number of time steps
dx <- L / nx # Spatial step
dt <- T_final / nt # Time step

# Stability check (CFL condition for explicit upwind)
cfl <- c * dt / dx
if (cfl > 1.0) {
  warning("CFL condition violated (cfl = ", cfl, "). Results may be unstable.")
}

# Create spatial grid
x <- seq(0, L - dx, by = dx) # periodic: x[1] and x[nx+1] should be same, so we stop at L-dx

# Initialize solution matrix: rows = time steps, columns = space
u <- matrix(0, nrow = nt + 1, ncol = nx)

# Initial condition: Gaussian pulse centered at x = 0.5
sigma <- 0.1 # width of Gaussian
u[1, ] <- exp(-(x - 0.5)^2 / (2 * sigma^2))

# Time stepping using 1st-order upwind scheme (for c > 0)
for (n in 1:nt) {
  for (i in 1:nx) {
    # Upwind: since c > 0, information comes from left
    if (i == 1) {
      # Periodic boundary: left neighbor of first point is last point
      u_left <- u[n, nx]
    } else {
      u_left <- u[n, i - 1]
    }
    # Update using upwind scheme: u^{n+1}_i = u^n_i - c*dt/dx * (u^n_i - u^n_{i-1})
    u[n + 1, i] <- u[n, i] - c * dt / dx * (u[n, i] - u_left)
  }
}

# Prepare data for ggplot2: reshape to long format
# Create vectors for space, time, and amplitude
space_vec <- rep(x, times = nt + 1)
time_vec <- rep(seq(0, T_final, length.out = nt + 1), each = nx)
amp_vec <- as.vector(u)

# Create dataframe
df <- data.frame(
  x = space_vec,
  t = time_vec,
  u = amp_vec
)

# Plot: Space-Time Heatmap
p <- ggplot(df, aes(x = x, y = t, fill = u)) +
  geom_raster() +
  scale_fill_viridis_c(option = "plasma", name = "Amplitude") +
  labs(
    title = "Space-Time Waterfall: 1D Advection of Gaussian Pulse",
    x = "Space (x)",
    y = "Time (t)",
    fill = "Wave Amplitude (u)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p)

# Optional: save plot
# ggsave("advection_space_time_waterfall.png", plot = p, width = 10, height = 8, dpi = 300)
