# Libraries
library(ggplot2)
library(dplyr)

# Parameters (The 'Coral' preset)
Du <- 0.16 # Diffusion rate of U
Dv <- 0.08 # Diffusion rate of V
F <- 0.035 # Feed rate
k <- 0.060 # Kill rate
N <- 100 # Grid size (N x N)
dt <- 1 # Time step size
steps <- 2000 # Total iterations

# Initialize concentrations: U(x,y,t=0), V(x,y,t=0)
U <- matrix(1, nrow = N, ncol = N)
V <- matrix(0, nrow = N, ncol = N)

# Initial perturbation (10x10 noisy square in the center)
set.seed(123)
V[(N / 2 - 5):(N / 2 + 4), (N / 2 - 5):(N / 2 + 4)] <- runif(
  100,
  min = 0.5,
  max = 1
)

# Laplacian kernel (finite differences)
laplacian <- function(X) {
  X_ext <- rbind(X[N, ], X, X[1, ]) # Periodic boundary (top-bottom)
  X_ext <- cbind(X_ext[, N], X_ext, X_ext[, 1]) # Periodic (left-right)
  lap <- (
    X_ext[1:N, 2:(N + 1)] +
      X_ext[3:(N + 2), 2:(N + 1)] + # Up/Down
      X_ext[2:(N + 1), 1:N] +
      X_ext[2:(N + 1), 3:(N + 2)] - # Left/Right
      4 * X_ext[2:(N + 1), 2:(N + 1)] # Center
  )
  lap
}

# Main Reaction-Diffusion Loop
for (step in 1:steps) {
  # Compute Laplacians ∇²U and ∇²V
  Lu <- laplacian(U)
  Lv <- laplacian(V)

  # Reaction-diffusion equations (Euler method)
  U_new <- U + dt * (Du * Lu - U * V^2 + F * (1 - U))
  V_new <- V + dt * (Dv * Lv + U * V^2 - (F + k) * V)

  # Enforce boundary conditions (0 ≤ U,V ≤ 1)
  U_new[U_new < 0] <- 0
  U_new[U_new > 1] <- 1
  V_new[V_new < 0] <- 0
  V_new[V_new > 1] <- 1

  # Update concentrations
  U <- U_new
  V <- V_new

  # Progress bar
  if (step %% 100 == 0) cat(sprintf("Step %d/%d\n", step, steps))
}

# Convert V (final state) to a data.frame for ggplot
df_plot <- expand.grid(x = 1:N, y = 1:N) %>%
  mutate(V = as.vector(V))

# High-Res Heatmap Visualization
ggplot(df_plot, aes(x, y, fill = V)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    title = "Gray-Scott Reaction-Diffusion: Coral Pattern",
    subtitle = sprintf("Du=%.2f, Dv=%.2f, F=%.3f, k=%.3f", Du, Dv, F, k),
    x = "",
    y = ""
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

# Optional: Save high-res PNG
ggsave("gray_scott_coral.png", width = 8, height = 8, dpi = 300)
