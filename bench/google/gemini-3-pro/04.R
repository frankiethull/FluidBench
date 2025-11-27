# ==========================================
# 1. SETUP & LIBRARIES
# ==========================================
library(rayshader)

# Simulation Parameters
nx <- 200 # Grid size x
ny <- 200 # Grid size y
nt <- 300 # Time steps (run enough for ripple to expand)
dt <- 0.1 # Time step
dx <- 1.0 # Spatial step
g <- 9.81 # Gravity
H <- 10.0 # Mean water depth
damping <- 0.995 # Energy loss (friction) to keep it stable

# ==========================================
# 2. INITIALIZATION
# ==========================================

# Allocate Matrices
# h = Deviation from mean water level (surface height)
h <- matrix(0, nrow = nx, ncol = ny)
# u = Velocity in x direction
u <- matrix(0, nrow = nx, ncol = ny)
# v = Velocity in y direction
v <- matrix(0, nrow = nx, ncol = ny)

# Initial Condition: Gaussian Drop
# Create a peak in the center
x <- 1:nx
y <- 1:ny
grid <- expand.grid(x = x, y = y)
center_x <- nx / 2
center_y <- ny / 2

# Apply Gaussian to h matrix
# Height 2.0, width determined by denominator
h_init <- 5.0 * exp(-((grid$x - center_x)^2 + (grid$y - center_y)^2) / 20)
h <- matrix(h_init, nrow = nx, ncol = ny)

# ==========================================
# 3. PHYSICS LOOP (Shallow Water)
# ==========================================
# We use a semi-implicit finite difference scheme (Forward-Backward).
# This is very stable for gravity waves.

cat("Simulating Ripples...\n")

for (t in 1:nt) {
  # --- Step 1: Update Velocity (u, v) from Height Gradients ---
  # u(new) = u(old) - g * (dh/dx) * dt

  # Calculate gradients using shifting (Vectorized)
  # h_ip1 = h[i+1], h_im1 = h[i-1]

  # Shift H for X-gradient
  h_x_next <- rbind(h[2:nx, ], h[nx, ]) # Periodic wrap for simplicity
  h_x_prev <- rbind(h[1, ], h[1:(nx - 1), ])
  dh_dx <- (h_x_next - h_x_prev) / (2 * dx)

  # Shift H for Y-gradient
  h_y_next <- cbind(h[, 2:ny], h[, ny])
  h_y_prev <- cbind(h[, 1], h[, 1:(ny - 1)])
  dh_dy <- (h_y_next - h_y_prev) / (2 * dx)

  # Update Velocities
  u <- u - g * dh_dx * dt
  v <- v - g * dh_dy * dt

  # Apply Damping (Friction)
  u <- u * damping
  v <- v * damping

  # --- Step 2: Update Height (h) from Velocity Divergence ---
  # h(new) = h(old) - H * (du/dx + dv/dy) * dt

  # Shift U for X-divergence
  u_x_next <- rbind(u[2:nx, ], u[nx, ])
  u_x_prev <- rbind(u[1, ], u[1:(nx - 1), ])
  du_dx <- (u_x_next - u_x_prev) / (2 * dx)

  # Shift V for Y-divergence
  v_y_next <- cbind(v[, 2:ny], v[, ny])
  v_y_prev <- cbind(v[, 1], v[, 1:(ny - 1)])
  dv_dy <- (v_y_next - v_y_prev) / (2 * dx)

  # Update Height
  h <- h - H * (du_dx + dv_dy) * dt
}

# ==========================================
# 4. 3D RENDERING (Rayshader)
# ==========================================

cat("Rendering 3D Scene...\n")

# Prepare the water palette
# We create a color map that is mostly teal/blue
water_palette <- create_texture(
  "#1a3b6e", # Deep Blue
  "#2a5b8e", # Mid Blue
  "#4a8bbe", # Light Blue
  "#8adbf0", # Highlight/Foam
  "#ffffff" # Specular
)

# Render
# Note: zscale controls vertical exaggeration.
# Since ripples are small, we need a low zscale number (inverse logic in rayshader)
# to make them look tall enough.
h %>%
  height_shade(texture = water_palette) %>%
  add_shadow(
    ray_shade(h, zscale = 0.5, sunaltitude = 45, sunangle = 315),
    0.5
  ) %>%
  add_shadow(ambient_shade(h, zscale = 0.5), 0.3) %>%
  plot_3d(
    h,
    zscale = 0.5, # Exaggerate height for dramatic effect
    fov = 30, # Field of view
    theta = 45, # Camera Rotation
    phi = 40, # Camera Angle
    windowsize = c(800, 600),
    zoom = 0.6,
    water = TRUE, # Enable glossy water effect
    waterdepth = 0, # Water level relative to matrix
    wateralpha = 0.5, # Transparency
    watercolor = "lightblue",
    solid = FALSE, # Do not draw dirt/sides
    background = "black"
  )

# Adjust Camera slightly to get a good "Pond" view
render_camera(theta = 45, phi = 45, zoom = 0.6)

# Save snapshot
# render_snapshot("ripple_3d.png")
