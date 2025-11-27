# Load required packages
library(ggplot2)
library(dplyr)

# Simulation parameters
nx <- 100       # Number of spatial points
nt <- 100       # Number of time steps
L <- 2          # Domain length [0, 2]
c <- 1          # Advection speed
dx <- L / nx    # Spatial step
dt <- 0.01      # Time step (must satisfy CFL condition: c*dt/dx <= 1)

# Initialize spatial grid
x <- seq(0, L, length.out = nx)

# Initial condition: Gaussian pulse centered at x=0.5
sigma <- 0.1    # Width of Gaussian
x0 <- 0.5       # Initial position
u0 <- exp(-(x - x0)^2 / (2 * sigma^2))

# Initialize solution matrix (space-time)
u <- matrix(0, nrow = nt, ncol = nx)
u[1, ] <- u0   # Set initial condition

# Time-stepping loop (using upwind scheme)
for (n in 1:(nt-1)) {
  for (i in 1:nx) {
    # Upwind scheme for advection
    if (c > 0) {
      # Periodic boundary condition
      im1 <- ifelse(i == 1, nx, i - 1)
      u[n+1, i] <- u[n, i] - c * dt / dx * (u[n, i] - u[n, im1])
    } else {
      # For c < 0 (not used here but included for completeness)
      ip1 <- ifelse(i == nx, 1, i + 1)
      u[n+1, i] <- u[n, i] - c * dt / dx * (u[n, ip1] - u[n, i])
    }
  }
}

# Create time vector
t <- seq(0, (nt-1)*dt, length.out = nt)

# Reshape data for ggplot2
df <- data.frame(
  x = rep(x, each = nt),
  t = rep(t, times = nx),
  u = as.vector(u)
)

# Create space-time heatmap
ggplot(df, aes(x = x, y = t, fill = u)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", name = "Amplitude") +
  labs(title = "1D Advection of Gaussian Pulse (Space-Time Heatmap)",
       x = "Space (x)",
       y = "Time (t)") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))