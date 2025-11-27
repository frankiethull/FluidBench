# --------------------------------------------------------------------------
# Problem 4: The "Pond Ripple" (2D Shallow Water Equations)
# --------------------------------------------------------------------------
# install.packages(c("rayshader", "viridis", "ggplot2", "png"))
library(rayshader)
library(viridis)
library(ggplot2)

# 1. SETUP: Define simulation parameters
# -------------------------------------
# Grid parameters
nx <- 101 # Number of grid points in x (odd number is good for a center)
ny <- 101 # Number of grid points in y
L <- 10.0 # Domain size (L x L meters)
dx <- L / (nx - 1)
dy <- L / (ny - 1)

# Physical parameters
g <- 9.81 # Acceleration due to gravity

# Simulation parameters
t_final <- 3.0 # Run simulation until this time (s)
h_dry <- 1e-6 # A "dry" state to prevent division by zero

# Lax-Wendroff Scheme (2D) requires a CFL condition < 1/sqrt(2)
CFL <- 0.45
dt <- CFL * min(dx, dy) / sqrt(g * 1.0) # dt based on max expected wave speed (sqrt(g*h))

# 2. INITIALIZATION: Set up grids and initial state
# -------------------------------------------------
# State variables: h (height), hu (x-momentum), hv (y-momentum)
h <- matrix(1.0, nrow = nx, ncol = ny) # Base water height of 1.0m
hu <- matrix(0.0, nrow = nx, ncol = ny)
hv <- matrix(0.0, nrow = nx, ncol = ny)

# Initial Condition: A Gaussian drop in the center
x <- seq(0, L, length.out = nx)
y <- seq(0, L, length.out = ny)
drop_center_x <- L / 2
drop_center_y <- L / 2
drop_radius <- 1.0
drop_amplitude <- 0.1

for (i in 1:nx) {
  for (j in 1:ny) {
    dist_sq <- (x[i] - drop_center_x)^2 + (y[j] - drop_center_y)^2
    h[i, j] <- h[i, j] + drop_amplitude * exp(-dist_sq / (2 * drop_radius^2))
  }
}

# Store initial state for comparison
h_initial <- h

# 3. SIMULATION LOOP (Lax-Wendroff)
# ---------------------------------
cat("Starting Shallow Water simulation...\n")
n_steps <- ceiling(t_final / dt)

for (n in 1:n_steps) {
  # --- Create padded matrices for stencil operations (avoids boundary checks) ---
  h_pad <- rbind(h_dry, cbind(h_dry, h, h_dry), h_dry)
  hu_pad <- rbind(0, cbind(0, hu, 0), 0)
  hv_pad <- rbind(0, cbind(0, hv, 0), 0)

  # --- Calculate fluxes at half-steps (n+1/2) ---
  # These are the core of the Lax-Wendroff method

  # Velocities at full step n
  u <- hu_pad[2:(nx + 1), 2:(ny + 1)] /
    pmax(h_pad[2:(nx + 1), 2:(ny + 1)], h_dry)
  v <- hv_pad[2:(nx + 1), 2:(ny + 1)] /
    pmax(h_pad[2:(nx + 1), 2:(ny + 1)], h_dry)

  # Fluxes F (x-direction) and G (y-direction) at half-steps
  F_h_half <- 0.5 *
    ((hu_pad[3:(nx + 2), 2:(ny + 1)] + hu_pad[2:(nx + 1), 2:(ny + 1)]) -
      dt /
        dx *
        ((hu_pad[3:(nx + 2), 2:(ny + 1)]^2 /
          pmax(h_pad[3:(nx + 2), 2:(ny + 1)], h_dry) +
          0.5 * g * h_pad[3:(nx + 2), 2:(ny + 1)]^2) -
          (hu_pad[2:(nx + 1), 2:(ny + 1)]^2 /
            pmax(h_pad[2:(nx + 1), 2:(ny + 1)], h_dry) +
            0.5 * g * h_pad[2:(nx + 1), 2:(ny + 1)]^2)))

  G_h_half <- 0.5 *
    ((hv_pad[2:(nx + 1), 3:(ny + 2)] + hv_pad[2:(nx + 1), 2:(ny + 1)]) -
      dt /
        dy *
        ((hu_pad[2:(nx + 1), 3:(ny + 2)] *
          hv_pad[2:(nx + 1), 3:(ny + 2)] /
          pmax(h_pad[2:(nx + 1), 3:(ny + 2)], h_dry)) -
          (hu_pad[2:(nx + 1), 2:(ny + 1)] *
            hv_pad[2:(nx + 1), 2:(ny + 1)] /
            pmax(h_pad[2:(nx + 1), 2:(ny + 1)], h_dry))))

  F_hu_half <- 0.5 *
    ((hu_pad[3:(nx + 2), 2:(ny + 1)]^2 /
      pmax(h_pad[3:(nx + 2), 2:(ny + 1)], h_dry) +
      0.5 * g * h_pad[3:(nx + 2), 2:(ny + 1)]^2 +
      hu_pad[2:(nx + 1), 2:(ny + 1)]^2 /
        pmax(h_pad[2:(nx + 1), 2:(ny + 1)], h_dry) +
      0.5 * g * h_pad[2:(nx + 1), 2:(ny + 1)]^2) -
      dt /
        dx *
        ((hu_pad[3:(nx + 2), 2:(ny + 1)]^3 /
          pmax(h_pad[3:(nx + 2), 2:(ny + 1)], h_dry)^2 +
          g * hu_pad[3:(nx + 2), 2:(ny + 1)] * h_pad[3:(nx + 2), 2:(ny + 1)]) -
          (hu_pad[2:(nx + 1), 2:(ny + 1)]^3 /
            pmax(h_pad[2:(nx + 1), 2:(ny + 1)], h_dry)^2 +
            g *
              hu_pad[2:(nx + 1), 2:(ny + 1)] *
              h_pad[2:(nx + 1), 2:(ny + 1)])))

  G_hu_half <- 0.5 *
    ((hu_pad[2:(nx + 1), 3:(ny + 2)] *
      hv_pad[2:(nx + 1), 3:(ny + 2)] /
      pmax(h_pad[2:(nx + 1), 3:(ny + 2)], h_dry) +
      hu_pad[2:(nx + 1), 2:(ny + 1)] *
        hv_pad[2:(nx + 1), 2:(ny + 1)] /
        pmax(h_pad[2:(nx + 1), 2:(ny + 1)], h_dry)) -
      dt /
        dy *
        ((hu_pad[2:(nx + 1), 3:(ny + 2)]^2 *
          hv_pad[2:(nx + 1), 3:(ny + 2)] /
          pmax(h_pad[2:(nx + 1), 3:(ny + 2)], h_dry)^2 +
          0.5 * g * h_pad[2:(nx + 1), 3:(ny + 2)]^2) -
          (hu_pad[2:(nx + 1), 2:(ny + 1)]^2 *
            hv_pad[2:(nx + 1), 2:(ny + 1)] /
            pmax(h_pad[2:(nx + 1), 2:(ny + 1)], h_dry)^2 +
            0.5 * g * h_pad[2:(nx + 1), 2:(ny + 1)]^2)))

  # (Similarly for hv fluxes, but we only need h and hu for the update)

  # --- Update to full step (n+1) using Lax-Wendroff formula ---
  h_new <- h -
    dt / dx * (F_h_half - F_h_half_shifted_left) -
    dt / dy * (G_h_half - G_h_half_shifted_down)
  hu_new <- hu -
    dt / dx * (F_hu_half - F_hu_half_shifted_left) -
    dt / dy * (G_hu_half - G_hu_half_shifted_down)
  # (And a similar update for hv)

  # Simplified update based on standard textbook formula
  h_new <- h -
    dt /
      (2 * dx) *
      (hu_pad[3:(nx + 2), 2:(ny + 1)] - hu_pad[1:nx, 2:(ny + 1)]) -
    dt / (2 * dy) * (hv_pad[2:(nx + 1), 3:(ny + 2)] - hv_pad[2:(nx + 1), 1:ny])

  hu_new <- hu -
    dt /
      (2 * dx) *
      ((hu_pad[3:(nx + 2), 2:(ny + 1)]^2 /
        pmax(h_pad[3:(nx + 2), 2:(ny + 1)], h_dry) +
        0.5 * g * h_pad[3:(nx + 2), 2:(ny + 1)]^2) -
        (hu_pad[1:nx, 2:(ny + 1)]^2 /
          pmax(h_pad[1:nx, 2:(ny + 1)], h_dry) +
          0.5 * g * h_pad[1:nx, 2:(ny + 1)]^2)) -
    dt /
      (2 * dy) *
      ((hu_pad[2:(nx + 1), 3:(ny + 2)] *
        hv_pad[2:(nx + 1), 3:(ny + 2)] /
        pmax(h_pad[2:(nx + 1), 3:(ny + 2)], h_dry)) -
        (hu_pad[2:(nx + 1), 1:ny] *
          hv_pad[2:(nx + 1), 1:ny] /
          pmax(h_pad[2:(nx + 1), 1:ny], h_dry)))

  hv_new <- hv -
    dt /
      (2 * dx) *
      ((hu_pad[3:(nx + 2), 2:(ny + 1)] *
        hv_pad[3:(nx + 2), 2:(ny + 1)] /
        pmax(h_pad[3:(nx + 2), 2:(ny + 1)], h_dry)) -
        (hu_pad[1:nx, 2:(ny + 1)] *
          hv_pad[1:nx, 2:(ny + 1)] /
          pmax(h_pad[1:nx, 2:(ny + 1)], h_dry))) -
    dt /
      (2 * dy) *
      ((hv_pad[2:(nx + 1), 3:(ny + 2)]^2 /
        pmax(h_pad[2:(nx + 1), 3:(ny + 2)], h_dry) +
        0.5 * g * h_pad[2:(nx + 1), 3:(ny + 2)]^2) -
        (hv_pad[2:(nx + 1), 1:ny]^2 /
          pmax(h_pad[2:(nx + 1), 1:ny], h_dry) +
          0.5 * g * h_pad[2:(nx + 1), 1:ny]^2))

  # --- Apply Boundary Conditions (Reflecting/Wall) ---
  # Zero-gradient (reflecting) by copying edge values
  h_new[1, ] <- h_new[2, ]
  h_new[nx, ] <- h_new[nx - 1, ]
  h_new[, 1] <- h_new[, 2]
  h_new[, ny] <- h_new[, ny - 1]
  hu_new[1, ] <- -hu_new[2, ]
  hu_new[nx, ] <- -hu_new[nx - 1, ] # Momentum flips sign
  hu_new[, 1] <- hu_new[, 2]
  hu_new[, ny] <- hu_new[, ny - 1]
  hv_new[1, ] <- hv_new[2, ]
  hv_new[nx, ] <- hv_new[nx - 1, ]
  hv_new[, 1] <- -hv_new[, 2]
  hv_new[, ny] <- -hv_new[, ny - 1] # Momentum flips sign

  # Update state
  h <- h_new
  hu <- hu_new
  hv <- hv_new

  # Print progress and check if wave has reached the edge
  if (n %% 100 == 0) {
    max_h <- max(h)
    cat(sprintf("Step %d/%d | Max Height: %.3f m\n", n, n_steps, max_h))
    if (max_h > 1.001 && max_h < 1.01) {
      # Check if wave front is at the edge
      # cat("Wave has reached the edge. Stopping early.\n")
      # break
    }
  }
}
cat("Simulation finished.\n")

# 4. VISUALIZATION with rayshader
# -------------------------------
# For rayshader, the matrix needs to be in x-y orientation, so we transpose it.
# We also reverse the rows so that (1,1) is the front-left corner.
heightmap <- h[nrow(h):1, ]

# Set up the 3D scene
window_size <- 1000
win <- make_scene(heightmap, heightmap_scale = 5) %>%
  add_water(color = "lightblue", depth = 0.5) %>%
  add_shadow(ray_shade(heightmap, zscale = 5, sunangle = 225), 0.5) %>%
  add_shadow(ambient_shade(heightmap), 0) %>%
  plot_3d(
    heightmap,
    zscale = 5,
    fov = 0,
    theta = 135,
    phi = 45,
    windowsize = c(window_size, window_size),
    zoom = 0.6
  )

# Render a high-quality snapshot
render_snapshot(
  filename = "pond_ripple_rayshader.png",
  clear = TRUE,
  title = "Pond Ripple (Shallow Water Equations)"
)

cat("Render saved to pond_ripple_rayshader.png\n")
