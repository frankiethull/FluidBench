library(ggplot2)
library(viridis)
library(reshape2)

# Parameters for "Coral" pattern
Du <- 0.16 # Diffusion rate for U
Dv <- 0.08 # Diffusion rate for V
F <- 0.035 # Feed rate
k <- 0.060 # Kill rate

# Grid parameters
n <- 100 # Grid size (n x n)
dt <- 1.0 # Time step
nsteps <- 2000 # Number of iterations

# Initialize concentration matrices
U <- matrix(1, nrow = n, ncol = n) # Chemical U
V <- matrix(0, nrow = n, ncol = n) # Chemical V

# Add small square of noise in the center for V
center_start <- floor(n / 2) - 5
center_end <- floor(n / 2) + 5
V[center_start:center_end, center_start:center_end] <- 1

# Add some random noise to break symmetry
set.seed(123) # For reproducibility
V <- V + matrix(runif(n * n, -0.01, 0.01), nrow = n, ncol = n)

# Function to compute Laplacian using 5-point stencil
compute_laplacian <- function(M) {
  laplacian <- matrix(0, nrow = n, ncol = n)

  # Interior points (using 5-point stencil)
  for (i in 2:(n - 1)) {
    for (j in 2:(n - 1)) {
      laplacian[i, j] <- (M[i - 1, j] +
        M[i + 1, j] +
        M[i, j - 1] +
        M[i, j + 1] -
        4 * M[i, j])
    }
  }

  # Boundary conditions (Neumann: zero flux)
  # Top and bottom boundaries
  for (j in 2:(n - 1)) {
    laplacian[1, j] <- (M[2, j] + M[1, j - 1] + M[1, j + 1] - 3 * M[1, j])
    laplacian[n, j] <- (M[n - 1, j] + M[n, j - 1] + M[n, j + 1] - 3 * M[n, j])
  }

  # Left and right boundaries
  for (i in 2:(n - 1)) {
    laplacian[i, 1] <- (M[i - 1, 1] + M[i + 1, 1] + M[i, 2] - 3 * M[i, 1])
    laplacian[i, n] <- (M[i - 1, n] + M[i + 1, n] + M[i, n - 1] - 3 * M[i, n])
  }

  # Corners
  laplacian[1, 1] <- (M[2, 1] + M[1, 2] - 2 * M[1, 1])
  laplacian[1, n] <- (M[2, n] + M[1, n - 1] - 2 * M[1, n])
  laplacian[n, 1] <- (M[n - 1, 1] + M[n, 2] - 2 * M[n, 1])
  laplacian[n, n] <- (M[n - 1, n] + M[n, n - 1] - 2 * M[n, n])

  return(laplacian)
}

# Alternative vectorized Laplacian function (faster)
compute_laplacian_fast <- function(M) {
  laplacian <- matrix(0, nrow = n, ncol = n)

  # Shifted matrices for neighbors
  M_up <- rbind(M[2:n, ], M[n, ])
  M_down <- rbind(M[1, ], M[1:(n - 1), ])
  M_left <- cbind(M[, 2:n], M[, n])
  M_right <- cbind(M[, 1], M[, 1:(n - 1)])

  # 5-point stencil
  laplacian <- M_up + M_down + M_left + M_right - 4 * M

  return(laplacian)
}

# Main simulation loop
for (step in 1:nsteps) {
  # Compute Laplacians
  lap_U <- compute_laplacian_fast(U)
  lap_V <- compute_laplacian_fast(V)

  # Reaction terms
  reaction <- U * V * V

  # Update equations
  dU <- Du * lap_U - reaction + F * (1 - U)
  dV <- Dv * lap_V + reaction - (F + k) * V

  # Euler integration
  U <- U + dt * dU
  V <- V + dt * dV

  # Ensure concentrations stay in reasonable bounds
  U[U < 0] <- 0
  U[U > 1] <- 1
  V[V < 0] <- 0
  V[V > 1] <- 1

  # Progress indicator
  if (step %% 100 == 0) {
    cat("Step:", step, "/", nsteps, "\n")
  }
}

# Create data frame for plotting
df <- melt(V)
names(df) <- c("x", "y", "concentration")

# Plot the final pattern
ggplot(df, aes(x = x, y = y, fill = concentration)) +
  geom_raster() +
  scale_fill_viridis_c(option = "plasma", name = "V Concentration") +
  coord_fixed() +
  labs(
    title = "Turing Pattern: Gray-Scott Model (Coral Preset)",
    subtitle = "Reaction-Diffusion System",
    x = "X",
    y = "Y"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save the plot
ggsave("turing_pattern.png", width = 8, height = 8, dpi = 300)
