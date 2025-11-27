# 2D Shallow Water Equations with Rayshader Visualization
library(rayshader)

# Simulation Parameters
nx <- 100 # Grid size (x)
ny <- 100 # Grid size (y)
t_max <- 100 # Time steps
dt <- 0.05 # Time step size
dx <- 1.0 # Spatial step size
g <- 9.81 # Gravity

# Initialize height and momentum fields
h <- matrix(1.0, nrow = nx, ncol = ny) # Base water depth
hu <- matrix(0, nrow = nx, ncol = ny) # Momentum x
hv <- matrix(0, nrow = nx, ncol = ny) # Momentum y

# Add Gaussian drop in center
x <- seq(-5, 5, length.out = nx)
y <- seq(-5, 5, length.out = ny)
drop <- 0.5 * exp(-(outer(x, y, function(x, y) x^2 + y^2)) / 2)
h <- h + drop

# Lax-Wendroff Scheme Implementation
lax_wendroff <- function(h, hu, hv, dt, dx, g) {
  # Calculate derivatives
  dhu_dx <- diff(c(hu, hu[1, ]), dim = 1) / dx
  dhv_dy <- diff(c(hv, hv[, 1]), dim = 2) / dx

  # Update momentum
  hu_half <- hu - g * h * dhu_dx * dt / 2
  hv_half <- hv - g * h * dhv_dy * dt / 2

  h_half <- h -
    (diff(c(hu_half, hu_half[1, ]), dim = 1) +
      diff(c(hv_half, hv_half[, 1]), dim = 2)) *
      dt /
      (2 * dx)

  # Full step update
  hu <- hu - g * h_half * diff(c(h_half, h_half[1, ]), dim = 1) * dt / dx
  hv <- hv - g * h_half * diff(c(h_half, h_half[, 1]), dim = 2) * dt / dx
  h <- h -
    (diff(c(hu, hu[1, ]), dim = 1) + diff(c(hv, hv[, 1]), dim = 2)) * dt / dx

  return(list(h = h, hu = hu, hv = hv))
}

# Main Simulation Loop
for (t in 1:t_max) {
  # Apply Lax-Wendroff scheme
  result <- lax_wendroff(h, hu, hv, dt, dx, g)
  h <- result$h
  hu <- result$hu
  hv <- result$hv

  # Open boundary conditions (simple absorption)
  h[1, ] <- 1.0
  h[nx, ] <- 1.0
  h[, 1] <- 1.0
  h[, ny] <- 1.0
  hu[1, ] <- 0
  hu[nx, ] <- 0
  hu[, 1] <- 0
  hu[, ny] <- 0
  hv[1, ] <- 0
  hv[nx, ] <- 0
  hv[, 1] <- 0
  hv[, ny] <- 0
}

# Rayshader 3D Visualization
elevation <- h - 1.0 # Convert to surface elevation

# Create water texture
water_texture <- height_shade(elevation, texture = "water")

# Generate 3D plot
plot_3d <- plot_3d(
  elevation,
  zscale = 5,
  solid = TRUE,
  background = "sky",
  shadowcolor = "dodgerblue4"
)

# Add realistic water rendering
plot_3d |>
  add_water(
    detect_water(elevation < 0.01),
    color = "dodgerblue",
    alpha = 0.8
  ) |>
  add_shadow(ray_shade(elevation), maxshadow = 0.5, darken = 0.7) |>
  add_shadow(ambient_shade(elevation), darken = 0.4) |>
  add_shadow(water_texture, darken = 0.3)

# Render final image
render_snapshot(
  title_text = "Pond Ripple Simulation",
  title_color = "white",
  title_bar_color = "black"
)
