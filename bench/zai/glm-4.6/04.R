# --------------------------------------------------------------
# 1. Load Required Libraries
# --------------------------------------------------------------
# install.packages(c("ggplot2", "rayshader"))
library(ggplot2)
library(rayshader)

# --------------------------------------------------------------
# 2. Define Simulation Parameters
# --------------------------------------------------------------
# Grid parameters
grid_dim <- 150 # Dimension of the square grid (e.g., 150x150)
dx <- 1.0 # Spatial step size

# Physical parameters
g <- 9.81 # Gravitational acceleration

# Time parameters (must satisfy CFL condition for stability)
# Wave speed c = sqrt(g*h). Max h will be ~1.1, so c_max ~ 3.3
# dt < dx / c_max  =>  dt < 1 / 3.3
dt <- 0.04 # Time step
n_steps <- 800 # Number of simulation steps

# Initial Gaussian drop parameters
drop_amplitude <- 0.5 # Height of the initial drop
drop_sigma <- 5.0 # Width of the drop

# --------------------------------------------------------------
# 3. Initialize State Variables
# --------------------------------------------------------------

# Create coordinate grids
x_grid <- (1:grid_dim) - (grid_dim / 2)
y_grid <- (1:grid_dim) - (grid_dim / 2)
coords <- expand.grid(x = x_grid, y = y_grid)

# Initialize state matrices
# h: water height
# hu: momentum in x-direction (h * u)
# hv: momentum in y-direction (h * v)
h <- matrix(1.0, nrow = grid_dim, ncol = grid_dim) # Start with flat water at height 1
hu <- matrix(0.0, nrow = grid_dim, ncol = grid_dim) # No initial momentum
hv <- matrix(0.0, nrow = grid_dim, ncol = grid_dim) # No initial momentum

# Add the Gaussian drop to the center of the height field
r_sq <- coords$x^2 + coords$y^2
h <- h + drop_amplitude * exp(-r_sq / (2 * drop_sigma^2))


# --------------------------------------------------------------
# 4. The Lax-Wendroff Simulation Loop
# --------------------------------------------------------------

cat("Starting Shallow Water simulation...\n")
pb <- txtProgressBar(min = 0, max = n_steps, style = 3)

for (n in 1:n_steps) {
  # --- Store old values ---
  h_old <- h
  hu_old <- hu
  hv_old <- hv

  # Add small epsilon to prevent division by zero
  epsilon <- 1e-6
  h_safe <- pmax(h_old, epsilon)

  # --- Calculate primitive variables (velocities) ---
  u <- hu_old / h_safe
  v <- hv_old / h_safe

  # --- Calculate fluxes F (in x-direction) and G (in y-direction) ---
  # Fluxes are 3D arrays: [grid_x, grid_y, variable_index]
  # Variable indices: 1=h, 2=hu, 3=hv

  # F_x = [hu, hu*u + 0.5*g*h^2, hu*v]
  F_x_h <- hu_old
  F_x_hu <- hu_old * u + 0.5 * g * h_old^2
  F_x_hv <- hu_old * v

  # G_y = [hv, hv*u, hv*v + 0.5*g*h^2]
  G_y_h <- hv_old
  G_y_hu <- hv_old * u
  G_y_hv <- hv_old * v + 0.5 * g * h_old^2

  # --- Predictor Step (Half-step fluxes at cell faces) ---
  # Calculate fluxes at x+1/2 and y+1/2 faces

  # Fluxes in x-direction
  F_x_half_h <- 0.5 *
    (F_x_h[2:grid_dim, ] + F_x_h[1:(grid_dim - 1), ]) -
    (dt / (2 * dx)) * (F_x_h[2:grid_dim, ] - F_x_h[1:(grid_dim - 1), ])
  F_x_half_hu <- 0.5 *
    (F_x_hu[2:grid_dim, ] + F_x_hu[1:(grid_dim - 1), ]) -
    (dt / (2 * dx)) * (F_x_hu[2:grid_dim, ] - F_x_hu[1:(grid_dim - 1), ])
  F_x_half_hv <- 0.5 *
    (F_x_hv[2:grid_dim, ] + F_x_hv[1:(grid_dim - 1), ]) -
    (dt / (2 * dx)) * (F_x_hv[2:grid_dim, ] - F_x_hv[1:(grid_dim - 1), ])

  # Fluxes in y-direction
  G_y_half_h <- 0.5 *
    (G_y_h[, 2:grid_dim] + G_y_h[, 1:(grid_dim - 1)]) -
    (dt / (2 * dx)) * (G_y_h[, 2:grid_dim] - G_y_h[, 1:(grid_dim - 1)])
  G_y_half_hu <- 0.5 *
    (G_y_hu[, 2:grid_dim] + G_y_hu[, 1:(grid_dim - 1)]) -
    (dt / (2 * dx)) * (G_y_hu[, 2:grid_dim] - G_y_hu[, 1:(grid_dim - 1)])
  G_y_half_hv <- 0.5 *
    (G_y_hv[, 2:grid_dim] + G_y_hv[, 1:(grid_dim - 1)]) -
    (dt / (2 * dx)) * (G_y_hv[, 2:grid_dim] - G_y_hv[, 1:(grid_dim - 1)])

  # --- Corrector Step (Full update) ---
  # h_new = h_old - dt/dx * (F_x_half_i+1/2 - F_x_half_i-1/2) - dt/dy * (G_y_half_j+1/2 - G_y_half_j-1/2)

  # Update h
  h[2:(grid_dim - 1), 2:(grid_dim - 1)] <- h_old[
    2:(grid_dim - 1),
    2:(grid_dim - 1)
  ] -
    (dt / dx) *
      (F_x_half_h[2:(grid_dim - 1), ] - F_x_half_h[1:(grid_dim - 2), ]) -
    (dt / dx) *
      (G_y_half_h[, 2:(grid_dim - 1)] - G_y_half_h[, 1:(grid_dim - 2)])

  # Update hu
  hu[2:(grid_dim - 1), 2:(grid_dim - 1)] <- hu_old[
    2:(grid_dim - 1),
    2:(grid_dim - 1)
  ] -
    (dt / dx) *
      (F_x_half_hu[2:(grid_dim - 1), ] - F_x_half_hu[1:(grid_dim - 2), ]) -
    (dt / dx) *
      (G_y_half_hu[, 2:(grid_dim - 1)] - G_y_half_hu[, 1:(grid_dim - 2)])

  # Update hv
  hv[2:(grid_dim - 1), 2:(grid_dim - 1)] <- hv_old[
    2:(grid_dim - 1),
    2:(grid_dim - 1)
  ] -
    (dt / dx) *
      (F_x_half_hv[2:(grid_dim - 1), ] - F_x_half_hv[1:(grid_dim - 2), ]) -
    (dt / dx) *
      (G_y_half_hv[, 2:(grid_dim - 1)] - G_y_half_hv[, 1:(grid_dim - 2)])

  # --- Apply Reflective Boundary Conditions ---
  # No normal flow at boundaries means momentum perpendicular to wall is zero.
  hu[1, ] <- 0
  hu[grid_dim, ] <- 0
  hv[, 1] <- 0
  hv[, grid_dim] <- 0
  # Also, no gradient in height at the wall (mirror values)
  h[1, ] <- h[2, ]
  h[grid_dim, ] <- h[grid_dim - 1, ]
  h[, 1] <- h[, 2]
  h[, grid_dim] <- h[, grid_dim - 1]

  setTxtProgressBar(pb, n)
}
close(pb)
cat("\nSimulation finished.\n")


# --------------------------------------------------------------
# 5. Visualize with rayshader
# --------------------------------------------------------------

cat("Preparing 3D render with rayshader...\n")

# First, convert the final height matrix to a data frame for ggplot
h_df <- expand.grid(x = 1:grid_dim, y = 1:grid_dim) %>%
  mutate(height = as.vector(h))

# Create the base ggplot object
# We use a blue gradient to color the water depth
water_colors <- c(
  "#003f5c",
  "#2f4b7c",
  "#665191",
  "#a05195",
  "#d45087",
  "#f95d6a",
  "#ff7c43",
  "#ffa600"
)
p <- ggplot(h_df, aes(x = x, y = y, fill = height)) +
  geom_raster() +
  scale_fill_gradientn(colors = water_colors, name = "Height") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")

# Now, use rayshader to render the plot
# plot_gg() creates the 3D heightmap from the ggplot object
plot_gg(
  p,
  width = grid_dim,
  height = grid_dim,
  scale = 250, # Height exaggeration factor
  multicore = TRUE,
  windowsize = c(1200, 800),
  zoom = 0.6,
  phi = 30,
  theta = -30
)

# render_water() adds a translucent, animated water surface on top
render_water(water_color = "#0080ff", waterlinez = -0.2, water_alpha = 0.8)

# render_camera() sets up the final camera view and lighting
render_camera(
  theta = 135,
  phi = 25,
  zoom = 0.65,
  lightdirection = c(100, 100, 100),
  lightintensity = 300
)

# render_snapshot() saves the final image
render_snapshot("Pond_Ripple_3D.png", clear = TRUE)

cat("3D render saved to 'Pond_Ripple_3D.png'\n")

# To view the interactive 3D scene in RStudio, run:
# rgl::rgl.open()
# rgl::rgl.clear()
# rgl::rgl.viewpoint(theta = 135, phi = 25, zoom = 0.65)
# plot_gg(p, multicore = TRUE, windowsize = c(1200, 800), scale = 250)
# render_water(water_color = "#0080ff", waterlinez = -0.2, water_alpha = 0.8)
