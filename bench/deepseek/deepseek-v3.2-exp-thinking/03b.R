library(ggplot2)
library(viridis)
library(reshape2)

# Parameters
Du <- 0.16
Dv <- 0.08
F <- 0.035
k <- 0.060

n <- 100
dt <- 1.0
nsteps <- 2000

# Initialize
U <- matrix(1, nrow = n, ncol = n)
V <- matrix(0, nrow = n, ncol = n)

# Initial perturbation
center_start <- floor(n / 2) - 5
center_end <- floor(n / 2) + 5
V[center_start:center_end, center_start:center_end] <- 1
set.seed(123)
V <- V + matrix(runif(n * n, -0.01, 0.01), nrow = n, ncol = n)

# Fast Laplacian function
compute_laplacian_fast <- function(M) {
  n <- nrow(M)
  M_up <- rbind(M[2:n, ], M[n, ])
  M_down <- rbind(M[1, ], M[1:(n - 1), ])
  M_left <- cbind(M[, 2:n], M[, n])
  M_right <- cbind(M[, 1], M[, 1:(n - 1)])
  return(M_up + M_down + M_left + M_right - 4 * M)
}

# Store patterns at different time points
patterns <- list()

for (step in 1:nsteps) {
  lap_U <- compute_laplacian_fast(U)
  lap_V <- compute_laplacian_fast(V)

  reaction <- U * V * V

  dU <- Du * lap_U - reaction + F * (1 - U)
  dV <- Dv * lap_V + reaction - (F + k) * V

  U <- U + dt * dU
  V <- V + dt * dV

  # Clamp values
  U[U < 0] <- 0
  U[U > 1] <- 1
  V[V < 0] <- 0
  V[V > 1] <- 1

  # Save patterns at key time points
  if (step %in% c(100, 500, 1000, 1500, 2000)) {
    df <- melt(V)
    names(df) <- c("x", "y", "concentration")
    df$step <- step
    patterns[[as.character(step)]] <- df
  }
}

# Combine all patterns
all_patterns <- do.call(rbind, patterns)

# Create faceted plot showing evolution
ggplot(all_patterns, aes(x = x, y = y, fill = concentration)) +
  geom_raster() +
  facet_wrap(~step, ncol = 3) +
  scale_fill_viridis_c(option = "plasma", name = "V Concentration") +
  coord_fixed() +
  labs(
    title = "Turing Pattern Evolution",
    subtitle = "Gray-Scott Reaction-Diffusion Model"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
