# Load required packages
library(rayshader)
library(dplyr)
library(ggplot2)

# Simulation parameters
nx <- 100 # Number of grid points in x
ny <- 100 # Number of grid points in y
dx <- 0.1 # Spatial step
dy <- 0.1 # Spatial step
dt <- 0.01 # Time step
g <- 9.81 # Gravitational acceleration
max_steps <- 200 # Maximum number of time steps

# Domain size
Lx <- nx * dx
Ly <- ny * dy

# Create grid
x <- seq(0, Lx, length.out = nx)
y <- seq(0, Ly, length.out = ny)

# Initialize variables
h <- matrix(1, nrow = nx, ncol = ny) # Water height (flat surface)
hu <- matrix(0, nrow = nx, ncol = ny) # x-momentum
hv <- matrix(0, nrow = nx, ncol = ny) # y-momentum

# Add Gaussian drop to center
sigma <- 2
x0 <- Lx / 2
y0 <- Ly / 2
for (i in 1:nx) {
  for (j in 1:ny) {
    h[i, j] <- h[i, j] +
      0.5 * exp(-((x[i] - x0)^2 + (y[j] - y0)^2) / (2 * sigma^2))
  }
}

# Lax-Wendroff scheme for shallow water equations
simulate_shallow_water <- function(h, hu, hv) {
  # Allocate new matrices
  h_new <- matrix(0, nrow = nx, ncol = ny)
  hu_new <- matrix(0, nrow = nx, ncol = ny)
  hv_new <- matrix(0, nrow = nx, ncol = ny)

  # Intermediate step (predictor)
  h_half <- h
  hu_half <- hu
  hv_half <- hv

  for (i in 2:(nx - 1)) {
    for (j in 2:(ny - 1)) {
      # Calculate fluxes
      Fx <- hu[i, j]
      Gx <- hu[i, j]^2 / h[i, j] + 0.5 * g * h[i, j]^2
      Fy <- hv[i, j]
      Gy <- hu[i, j] * hv[i, j] / h[i, j]

      # Update intermediate values
      h_half[i, j] <- h[i, j] -
        0.5 * dt / dx * (Fx[i + 1, j] - Fx[i - 1, j]) -
        0.5 * dt / dy * (Fy[i, j + 1] - Fy[i, j - 1])
      hu_half[i, j] <- hu[i, j] -
        0.5 * dt / dx * (Gx[i + 1, j] - Gx[i - 1, j]) -
        0.5 * dt / dy * (Gy[i, j + 1] - Gy[i, j - 1])
      hv_half[i, j] <- hv[i, j] -
        0.5 *
          dt /
          dx *
          (Fx[i + 1, j] *
            hv[i + 1, j] /
            h[i + 1, j] -
            Fx[i - 1, j] * hv[i - 1, j] / h[i - 1, j]) -
        0.5 * dt / dy * (Gy[i, j + 1] - Gy[i, j - 1])
    }
  }

  # Corrector step
  for (i in 2:(nx - 1)) {
    for (j in 2:(ny - 1)) {
      # Calculate fluxes at half step
      Fx_half <- hu_half[i, j]
      Gx_half <- hu_half[i, j]^2 / h_half[i, j] + 0.5 * g * h_half[i, j]^2
      Fy_half <- hv_half[i, j]
      Gy_half <- hu_half[i, j] * hv_half[i, j] / h_half[i, j]

      # Update final values
      h_new[i, j] <- h[i, j] -
        dt / dx * (Fx_half[i + 1, j] - Fx_half[i - 1, j]) -
        dt / dy * (Fy_half[i, j + 1] - Fy_half[i, j - 1])
      hu_new[i, j] <- hu[i, j] -
        dt / dx * (Gx_half[i + 1, j] - Gx_half[i - 1, j]) -
        dt / dy * (Gy_half[i, j + 1] - Gy_half[i, j - 1])
      hv_new[i, j] <- hv[i, j] -
        dt /
          dx *
          (Fx_half[i + 1, j] *
            hv_half[i + 1, j] /
            h_half[i + 1, j] -
            Fx_half[i - 1, j] * hv_half[i - 1, j] / h_half[i - 1, j]) -
        dt / dy * (Gy_half[i, j + 1] - Gy_half[i, j - 1])
    }
  }

  # Boundary conditions (reflective)
  h_new[, 1] <- h_new[, 2]
  h_new[, ny] <- h_new[, ny - 1]
  h_new[1, ] <- h_new[2, ]
  h_new[nx, ] <- h_new[nx - 1, ]

  hu_new[, 1] <- -hu_new[, 2]
  hu_new[, ny] <- -hu_new[, ny - 1]
  hu_new[1, ] <- -hu_new[2, ]
  hu_new[nx, ] <- -hu_new[nx - 1, ]

  hv_new[, 1] <- -hv_new[, 2]
  hv_new[, ny] <- -hv_new[, ny - 1]
  hv_new[1, ] <- -hv_new[2, ]
  hv_new[nx, ] <- -hv_new[nx - 1, ]

  return(list(h = h_new, hu = hu_new, hv = hv_new))
}

# Run simulation
for (step in 1:max_steps) {
  result <- simulate_shallow_water(h, hu, hv)
  h <- result$h
  hu <- result$hu
  hv <- result$hv

  # Print progress
  if (step %% 20 == 0) {
    cat("Step:", step, "/", max_steps, "\n")
  }
}

# Prepare data for rayshader
water_surface <- h

# Create 3D visualization with rayshader
water_surface %>%
  sphere_shade(texture = "imhof1") %>%
  add_shadow(ray_shade(water_surface, zscale = 10), 0.5) %>%
  add_water(
    detect_water(water_surface, waterlevel = 1),
    color = "lightblue"
  ) %>%
  plot_3d(
    water_surface,
    zscale = 10,
    fov = 60,
    theta = 45,
    phi = 45,
    windowsize = c(800, 800),
    zoom = 0.7
  )

# Alternative version with custom lighting
render_highquality(
  water_surface,
  zscale = 10,
  fov = 60,
  theta = 45,
  phi = 45,
  zoom = 0.7,
  water = TRUE,
  waterlevel = 1,
  wateralpha = 0.7,
  watercolor = "lightblue",
  waterlinecolor = "white",
  waterlinealpha = 0.5
)
