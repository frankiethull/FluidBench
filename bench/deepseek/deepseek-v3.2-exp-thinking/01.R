library(ggplot2)
library(dplyr)
library(tidyr)

# Parameters
xmin <- 0
xmax <- 2
nx <- 200
dx <- (xmax - xmin) / (nx - 1)
c <- 1
dt <- 0.5 * dx / c # CFL condition for stability
tmax <- 1.5
nt <- floor(tmax / dt)

# Spatial grid
x <- seq(xmin, xmax, length.out = nx)

# Initial condition: Gaussian pulse
u0 <- exp(-20 * (x - 0.5)^2)

# Matrix to store solutions (time steps x spatial points)
u_matrix <- matrix(0, nrow = nt + 1, ncol = nx)
u_matrix[1, ] <- u0

# Time integration using FTCS scheme
u_old <- u0
for (i in 1:nt) {
  u_new <- numeric(nx)
  for (j in 1:nx) {
    left <- ifelse(j == 1, nx - 1, j - 1)
    u_new[j] <- u_old[j] - c * dt / dx * (u_old[j] - u_old[left])
  }
  u_matrix[i + 1, ] <- u_new
  u_old <- u_new
}

# Create time vector
time <- seq(0, tmax, length.out = nt + 1)

# Reshape data into long format
df <- expand.grid(x = x, t = time)
df$u <- as.vector(u_matrix)

# Plot space-time heatmap
ggplot(df, aes(x = x, y = t, fill = u)) +
  geom_raster() +
  scale_fill_viridis_c() +
  labs(x = "Space (x)", y = "Time (t)", fill = "Amplitude") +
  theme_minimal()
