# Load required packages
library(ggplot2)
library(viridis)
library(dplyr)

# Simulation parameters
nx <- 100       # Grid size
ny <- 100
Du <- 0.16     # Diffusion rate for U
Dv <- 0.08     # Diffusion rate for V
F <- 0.035     # Feed rate
k <- 0.060     # Kill rate
dt <- 1.0      # Time step
steps <- 2000  # Number of iterations

# Initialize grids
U <- matrix(1, nrow = nx, ncol = ny)
V <- matrix(0, nrow = nx, ncol = ny)

# Add noise to center (10x10 square)
noise_center <- 45:54  # Center 10x10 area
V[noise_center, noise_center] <- runif(100, 0.25, 0.75)

# Laplacian kernel (5-point stencil)
laplacian <- function(X) {
  # Create shifted matrices
  X_left <- cbind(X[, -1], X[, ncol(X)])
  X_right <- cbind(X[, ncol(X)], X[, -ncol(X)])
  X_up <- rbind(X[-1, ], X[nrow(X), ])
  X_down <- rbind(X[nrow(X), ], X[-nrow(X), ])

  # Calculate Laplacian
  (X_left + X_right + X_up + X_down - 4 * X) / (1.0^2)
}

# Main simulation loop
for (step in 1:steps) {
  # Compute Laplacians
  d2U <- laplacian(U)
  d2V <- laplacian(V)

  # Reaction terms
  UV2 <- U * V^2
  reaction_U <- -UV2 + F * (1 - U)
  reaction_V <- UV2 - (F + k) * V

  # Update U and V
  U <- U + dt * (Du * d2U + reaction_U)
  V <- V + dt * (Dv * d2V + reaction_V)

  # Ensure concentrations stay within bounds
  U <- pmax(0, pmin(1, U))
  V <- pmax(0, pmin(1, V))

  # Print progress
  if (step %% 100 == 0) {
    cat("Step:", step, "/", steps, "\n")
  }
}

# Prepare data for plotting
df <- expand.grid(x = 1:nx, y = 1:ny)
df$V <- as.vector(V)

# Create high-resolution heatmap
ggplot(df, aes(x = x, y = y, fill = V)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", direction = -1, name = "V concentration") +
  labs(title = "Gray-Scott Reaction-Diffusion Pattern (Coral Preset)") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "right",
    legend.key.width = unit(1.5, "cm")
  )