# --------------------------------------------------------------------------
# Problem 3: The "Turing Pattern" (Reaction-Diffusion)
# --------------------------------------------------------------------------
# Load necessary libraries
# install.packages(c("ggplot2", "viridis", "dplyr", "tidyr"))
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)

# 1. SETUP: Define simulation parameters
# -------------------------------------
# Grid parameters
N <- 100 # Grid size (N x N)
n_steps <- 2000 # Number of iterations to run

# Gray-Scott model parameters (Coral preset)
Du <- 0.16 # Diffusion rate of U
Dv <- 0.08 # Diffusion rate of V
F <- 0.035 # Feed rate
k <- 0.060 # Kill rate

# Time step for the simulation
dt <- 1.0 # Must be small enough for numerical stability

# 2. INITIALIZATION: Set up grids and initial state
# -------------------------------------------------
# Create matrices for the concentrations of U and V
U <- matrix(1.0, nrow = N, ncol = N)
V <- matrix(0.0, nrow = N, ncol = N)

# --- Initial Condition ---
# A small square of noise in the center to kickstart the reaction
noise_size <- 10
center <- N / 2
start_idx <- floor(center - noise_size / 2)
end_idx <- floor(center + noise_size / 2)

# Set V to 1 in the central patch with some random noise
V[start_idx:end_idx, start_idx:end_idx] <- 1.0
V <- V + matrix(runif(N * N, -0.01, 0.01), nrow = N, ncol = N)
V[V < 0] <- 0 # Ensure V doesn't go negative
V[V > 1] <- 1 # Clamp V to 1

# 3. SIMULATION LOOP (Vectorized)
# -------------------------------
cat("Starting Gray-Scott simulation...\n")

# Pre-calculate the Laplacian kernel weights
laplace_U_weights <- c(Du, Du, -4 * Du, Du, Du)
laplace_V_weights <- c(Dv, Dv, -4 * Dv, Dv, Dv)

for (i in 1:n_steps) {
  # --- Calculate Laplacian (Vectorized) ---
  # This is the core optimization. We use shifted matrices instead of loops.
  # The c(NA, ...) and c(..., NA) create the shifted neighbors.
  U_neighbors <- array(0, dim = c(N, N, 5))
  V_neighbors <- array(0, dim = c(N, N, 5))

  U_neighbors[,, 1] <- cbind(NA, U[, -N]) # Left
  U_neighbors[,, 2] <- rbind(NA, U[-N, ]) # Top
  U_neighbors[,, 3] <- U # Center
  U_neighbors[,, 4] <- rbind(U[-1, ], NA) # Bottom
  U_neighbors[,, 5] <- cbind(U[, -1], NA) # Right

  V_neighbors[,, 1] <- cbind(NA, V[, -N])
  V_neighbors[,, 2] <- rbind(NA, V[-N, ])
  V_neighbors[,, 3] <- V
  V_neighbors[,, 4] <- rbind(V[-1, ], NA)
  V_neighbors[,, 5] <- cbind(V[, -1], NA)

  # Sum the neighbors, ignoring NAs
  laplace_U <- rowSums(U_neighbors * laplace_U_weights, na.rm = TRUE, dims = 2)
  laplace_V <- rowSums(V_neighbors * laplace_V_weights, na.rm = TRUE, dims = 2)

  # --- Reaction Terms ---
  reaction_U <- -U * V^2 + F * (1 - U)
  reaction_V <- U * V^2 - (F + k) * V

  # --- Update Concentrations ---
  # We need U_old for the V update, so we do it in two steps
  U_new <- U + dt * (laplace_U + reaction_U)
  V_new <- V + dt * (laplace_V + reaction_V)

  # Clamp values to be physically meaningful [0, 1]
  U <- pmax(0, pmin(1, U_new))
  V <- pmax(0, pmin(1, V_new))

  # Print progress
  if (i %% 200 == 0) {
    cat(sprintf("Step %d/%d\n", i, n_steps))
  }
}

cat("Simulation finished. Preparing plot...\n")

# 4. VISUALIZATION: Create the high-res heatmap
# ---------------------------------------------
# Reshape the final V matrix into a long data frame for ggplot2
vorticity_df <- as.data.frame(V)
names(vorticity_df) <- 1:N
vorticity_df$x <- 1:N
vorticity_long <- pivot_longer(
  vorticity_df,
  cols = -x,
  names_to = "y",
  values_to = "V_conc"
)
vorticity_long$y <- as.numeric(vorticity_long$y)

# Create the plot
turing_plot <- ggplot(vorticity_long, aes(x = x, y = y, fill = V_conc)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(option = "magma", name = "V Concentration") +
  coord_fixed() + # Crucial for a square, undistorted pattern
  labs(
    title = "Turing Pattern (Gray-Scott 'Coral' Preset)",
    subtitle = paste("After", n_steps, "iterations on a", N, "x", N, "grid"),
    x = "",
    y = ""
  ) +
  theme_void() + # Clean background for the pattern
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right"
  )

# Display the plot
print(turing_plot)

# Optional: Save the plot to a file
# ggsave("turing_pattern_coral.png", plot = turing_plot, width = 8, height = 8, dpi = 300)
