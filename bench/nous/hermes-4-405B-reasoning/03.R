# Gray-Scott Reaction-Diffusion Simulation
library(ggplot2)

# Simulation Parameters
nx <- 100 # Grid size (x)
ny <- 100 # Grid size (y)
t_max <- 2000 # Time steps
dt <- 1.0 # Time step size

# Gray-Scott Parameters (Coral preset)
Du <- 0.16 # Diffusion rate U
Dv <- 0.08 # Diffusion rate V
F <- 0.035 # Feed rate
k <- 0.060 # Kill rate

# Initialize concentrations
U <- matrix(1.0, nrow = nx, ncol = ny)
V <- matrix(0.0, nrow = nx, ncol = ny)

# Add initial noise square (10x10 in center)
center_x <- seq(nx / 2 - 5, nx / 2 + 4)
center_y <- seq(ny / 2 - 5, ny / 2 + 4)
V[center_x, center_y] <- 1.0

# Laplacian Kernel (5-point stencil)
laplacian <- function(X) {
  # Periodic boundary conditions
  X_c <- X[c(1, 1:nx), c(1:ny, 1)] # Center shifted
  X_l <- X[c(1, 1:nx), 1:ny] # Left
  X_r <- X[c(1, 1:nx), c(2:(ny + 1))] # Right
  X_u <- X[1:nx, 1:ny] # Up
  X_d <- X[c(2:(nx + 1)), 1:ny] # Down

  return(X_l + X_r + X_u + X_d - 4 * X_c)
}

# Time-stepping with explicit Euler method
for (t in 1:t_max) {
  # Calculate Laplacians
  lapU <- laplacian(U)
  lapV <- laplacian(V)

  # Calculate reaction terms
  uv2 <- U * V * V

  # Update concentrations
  dU <- Du * lapU - uv2 + F * (1 - U)
  dV <- Dv * lapV + uv2 - (F + k) * V

  U <- U + dU * dt
  V <- V + dV * dt

  # Ensure concentrations stay bounded
  U[U < 0] <- 0
  V[V < 0] <- 0
}

# Prepare data for plotting
df <- expand.grid(x = 1:nx, y = 1:ny)
df$V <- as.vector(V)

# Create high-resolution heatmap
ggplot(df, aes(x = x, y = y, fill = V)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  labs(
    title = "Gray-Scott Reaction-Diffusion Pattern (Coral)",
    x = "X Position",
    y = "Y Position",
    fill = "Concentration V"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )
