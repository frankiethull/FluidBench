# ========================================================
# Gray-Scott Reaction-Diffusion Model
# "Coral" / "Fingerprint" pattern
# High-resolution 100x100 grid → 2000 iterations
# ========================================================

library(tidyverse)
library(viridis)
library(scales)

# -----------------------------
# Parameters (Classic Coral preset)
# -----------------------------
Du <- 0.16
Dv <- 0.08
F <- 0.035
k <- 0.060

n <- 100 # grid size: 100x100
steps <- 2000 # total iterations

# Spatial step (grid spacing)
dx <- 1.0

# Time step (must satisfy stability: dt < dx²/(4*max(Du,Dv)))
dt <- 0.9 * dx^2 / (4 * max(Du, Dv)) # very safe

# Laplacian kernel (5-point stencil)
laplacian <- function(field) {
  # Padding with periodic boundaries
  pad <- rbind(field[n, ], field, field[1, ])
  pad <- cbind(pad[, n], pad, pad[, 1])

  # 5-point stencil: center - 4 neighbors
  return(
    pad[2:(n + 1), 3:(n + 2)] +
      pad[2:(n + 1), 1:n] +
      pad[1:n, 2:(n + 1)] +
      pad[3:(n + 2), 2:(n + 1)] -
      4 * pad[2:(n + 1), 2:(n + 1)]
  ) /
    (dx^2)
}

# -----------------------------
# Initialize fields
# -----------------------------
U <- matrix(1.0, n, n)
V <- matrix(0.0, n, n)

# Add small square perturbation in the center
center <- n %/% 2
half <- 5 # 10x10 square
set.seed(1234)
V[(center - half):(center + half), (center - half):(center + half)] <-
  runif((2 * half + 1)^2, 0.9, 1.0)

# Add a tiny bit of noise everywhere for symmetry breaking
V <- V + 0.01 * matrix(runif(n * n), n, n)

# -----------------------------
# Time integration loop
# -----------------------------
cat("Running Gray-Scott simulation (2000 steps)...\n")
pb <- txtProgressBar(min = 0, max = steps, style = 3)

for (i in 1:steps) {
  setTxtProgressBar(pb, i)

  Lu <- laplacian(U)
  Lv <- laplacian(V)

  # Reaction terms
  UV2 <- U * V^2
  dU <- Du * Lu - UV2 + F * (1 - U)
  dV <- Dv * Lv + UV2 - (F + k) * V

  # Update (Euler step)
  U <- U + dt * dU
  V <- V + dt * dV

  # Optional: clamp to [0,1] for stability (not strictly needed here)
  U[U < 0] <- 0
  U[U > 1] <- 1
  V[V < 0] <- 0
  V[V > 1] <- 1
}

close(pb)

# -----------------------------
# Prepare data for ggplot2
# -----------------------------
df <- expand.grid(x = 1:n, y = 1:n)
df$V <- as.vector(t(V)) # transpose so y increases upward

# -----------------------------
# Create the final high-resolution coral pattern
# -----------------------------
p <- ggplot(df, aes(x, y, fill = V)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    option = "A", # magma (perfect for coral!)
    name = "V concentration",
    limits = c(0, 1),
    guide = guide_colorbar(barwidth = 20, barheight = 1.5)
  ) +
  coord_fixed(expand = FALSE) +
  labs(
    title = "Gray-Scott Reaction-Diffusion: Coral Pattern",
    subtitle = paste(
      "F =",
      F,
      ", k =",
      k,
      "| Du =",
      Du,
      ", Dv =",
      Dv,
      "| 2000 steps"
    ),
    caption = "Classic Turing pattern (biological morphogenesis simulation)"
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 20,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    plot.caption = element_text(
      size = 10,
      color = "gray60",
      margin = margin(t = 10)
    ),
    legend.position = "bottom",
    plot.background = element_rect(fill = "black"),
    panel.background = element_rect(fill = "black"),
    plot.margin = margin(20, 20, 20, 20)
  )

print(p)

# -----------------------------
# Save high-res image
# -----------------------------
ggsave(
  "gray_scott_coral.png",
  plot = p,
  width = 10,
  height = 10.3,
  dpi = 300,
  bg = "black"
)

cat("\nCoral pattern saved as 'gray_scott_coral.png'\n")
