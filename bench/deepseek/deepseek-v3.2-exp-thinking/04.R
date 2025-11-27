library(rayshader)
library(ggplot2)
library(viridis)
library(plot3D)

# Simulation parameters
nx <- 200
ny <- 200
Lx <- 10
Ly <- 10
dx <- Lx / (nx - 1)
dy <- Ly / (ny - 1)
g <- 9.81 # gravity

# Time parameters
dt <- 0.01
nt <- 500
cfl <- sqrt(g * 1) * dt / min(dx, dy) # Check CFL condition
cat("CFL number:", cfl, "\n")

# Initialize arrays
h <- matrix(1, nrow = nx, ncol = ny) # Water height
hu <- matrix(0, nrow = nx, ncol = ny) # x-momentum
hv <- matrix(0, nrow = nx, ncol = ny) # y-momentum

# Add Gaussian drop in the center
x_center <- nx / 2
y_center <- ny / 2
for (i in 1:nx) {
  for (j in 1:ny) {
    r2 <- ((i - x_center) / 10)^2 + ((j - y_center) / 10)^2
    h[i, j] <- 1 - 0.2 * exp(-r2 / 0.5)
  }
}

# Function to apply reflective boundary conditions
apply_boundary_conditions <- function(M) {
  # Interior remains the same
  # Left and right boundaries
  M[1, ] <- M[2, ]
  M[nx, ] <- M[nx - 1, ]
  # Top and bottom boundaries
  M[, 1] <- M[, 2]
  M[, ny] <- M[, ny - 1]
  # Corners
  M[1, 1] <- M[2, 2]
  M[1, ny] <- M[2, ny - 1]
  M[nx, 1] <- M[nx - 1, 2]
  M[nx, ny] <- M[nx - 1, ny - 1]
  return(M)
}

# Lax-Wendroff scheme for shallow water equations
simulate_shallow_water <- function(h, hu, hv, nt) {
  # Storage for final state
  h_final <- h
  hu_final <- hu
  hv_final <- hv

  for (n in 1:nt) {
    # Create temporary arrays
    h_new <- h
    hu_new <- hu
    hv_new <- hv

    # Compute fluxes
    for (i in 2:(nx - 1)) {
      for (j in 2:(ny - 1)) {
        # x-direction fluxes
        F1_x <- hu[i, j]
        F2_x <- hu[i, j]^2 / h[i, j] + 0.5 * g * h[i, j]^2
        F3_x <- hu[i, j] * hv[i, j] / h[i, j]

        # y-direction fluxes
        G1_y <- hv[i, j]
        G2_y <- hu[i, j] * hv[i, j] / h[i, j]
        G3_y <- hv[i, j]^2 / h[i, j] + 0.5 * g * h[i, j]^2

        # Lax-Wendroff updates
        # x-direction
        h_new[i, j] <- h[i, j] -
          dt / (2 * dx) * (hu[i + 1, j] - hu[i - 1, j]) +
          (dt^2 / (2 * dx^2)) *
            (F1_x[i + 1, j] - 2 * F1_x[i, j] + F1_x[i - 1, j])

        hu_new[i, j] <- hu[i, j] -
          dt / (2 * dx) * (F2_x[i + 1, j] - F2_x[i - 1, j]) +
          (dt^2 / (2 * dx^2)) *
            (F2_x[i + 1, j] - 2 * F2_x[i, j] + F2_x[i - 1, j])

        # y-direction
        h_new[i, j] <- h_new[i, j] -
          dt / (2 * dy) * (hv[i, j + 1] - hv[i, j - 1]) +
          (dt^2 / (2 * dy^2)) *
            (G1_y[i, j + 1] - 2 * G1_y[i, j] + G1_y[i, j - 1])

        hv_new[i, j] <- hv[i, j] -
          dt / (2 * dy) * (G3_y[i, j + 1] - G3_y[i, j - 1]) +
          (dt^2 / (2 * dy^2)) *
            (G3_y[i, j + 1] - 2 * G3_y[i, j] + G3_y[i, j - 1])
      }
    }

    # Update arrays
    h <- apply_boundary_conditions(h_new)
    hu <- apply_boundary_conditions(hu_new)
    hv <- apply_boundary_conditions(hv_new)

    # Progress indicator
    if (n %% 50 == 0) {
      cat("Time step:", n, "/", nt, "\n")
    }
  }

  return(list(h = h, hu = hu, hv = hv))
}

# Run simulation
cat("Starting shallow water simulation...\n")
result <- simulate_shallow_water(h, hu, hv, nt)
h_final <- result$h

# Alternative simpler simulation using wave equation (more stable)
simulate_waves <- function(h, nt) {
  h_prev <- h
  h_curr <- h
  h_next <- h

  c <- 0.5 # Wave speed

  for (n in 1:nt) {
    for (i in 2:(nx - 1)) {
      for (j in 2:(ny - 1)) {
        # 2D wave equation discretization
        laplacian <- (h_curr[i + 1, j] +
          h_curr[i - 1, j] +
          h_curr[i, j + 1] +
          h_curr[i, j - 1] -
          4 * h_curr[i, j])
        h_next[i, j] <- 2 *
          h_curr[i, j] -
          h_prev[i, j] +
          (c * dt / dx)^2 * laplacian
      }
    }

    # Apply boundary conditions
    h_next <- apply_boundary_conditions(h_next)

    # Update arrays
    h_prev <- h_curr
    h_curr <- h_next

    if (n %% 100 == 0) {
      cat("Wave simulation step:", n, "/", nt, "\n")
    }
  }

  return(h_curr)
}

# Run wave simulation (more stable alternative)
cat("Running wave simulation...\n")
h_final <- simulate_waves(h, nt)

# Create 3D visualization with rayshader
cat("Creating 3D visualization...\n")

# Create water surface matrix
water_surface <- h_final

# Create color palette for water
water_palette <- colorRampPalette(c(
  "#000066",
  "#0000FF",
  "#0066FF",
  "#66CCFF",
  "#FFFFFF"
))(256)

# Generate the 3D plot
water_surface %>%
  sphere_shade(texture = "imhof1") %>%
  add_water(detect_water(water_surface, zscale = 0.1), color = "imhof1") %>%
  plot_3d(
    water_surface,
    zscale = 0.1,
    water = TRUE,
    waterdepth = 0,
    wateralpha = 0.5,
    watercolor = "lightblue",
    theta = -45,
    phi = 30,
    windowsize = c(1000, 800)
  )

# Add lighting
render_snapshot(clear = FALSE)

# Alternative: More advanced rendering with multiple light sources
cat("Rendering advanced visualization...\n")

# Clear previous plot
rgl::clear3d()

# Create enhanced visualization
water_surface %>%
  height_shade(
    texture = grDevices::colorRampPalette(c(
      "#000033",
      "#000088",
      "#0000FF",
      "#0088FF",
      "#88CCFF"
    ))(256)
  ) %>%
  add_overlay(
    sphere_shade(
      water_surface,
      texture = "imhof1",
      colorintensity = 5,
      progbar = FALSE
    ),
    alphalayer = 0.8
  ) %>%
  add_shadow(ray_shade(water_surface, zscale = 0.1, lambert = TRUE), 0.5) %>%
  add_shadow(ambient_shade(water_surface, zscale = 0.1), 0) %>%
  plot_3d(
    water_surface,
    zscale = 0.1,
    water = TRUE,
    waterdepth = min(water_surface),
    wateralpha = 0.7,
    watercolor = "#88AAFF",
    theta = -30,
    phi = 35,
    windowsize = c(1200, 900),
    soliddepth = -max(water_surface),
    background = "lightblue"
  )

# Render high quality snapshot
render_highquality(
  filename = "pond_ripples.png",
  light = TRUE,
  lightdirection = c(315, 45),
  lightcolor = c("white", "lightblue"),
  lightintensity = c(800, 200),
  samples = 128,
  width = 1600,
  height = 1200,
  aperture = 32
)

# Create animation of ripple propagation (optional)
create_ripple_animation <- function() {
  # Re-run simulation and save frames
  h_temp <- h
  frames <- list()

  for (frame in 1:20) {
    h_temp <- simulate_waves(h_temp, 25) # 25 steps per frame

    # Create and save frame
    h_temp %>%
      sphere_shade(texture = "imhof1") %>%
      add_water(detect_water(h_temp, zscale = 0.1), color = "imhof1") %>%
      plot_3d(
        h_temp,
        zscale = 0.1,
        water = TRUE,
        waterdepth = 0,
        wateralpha = 0.5,
        watercolor = "lightblue",
        theta = -45,
        phi = 30
      )

    render_snapshot(filename = paste0("frame_", sprintf("%03d", frame), ".png"))
    rgl::clear3d()

    cat("Rendered frame:", frame, "/20\n")
  }
}

# Uncomment to create animation (takes longer)
# create_ripple_animation()

cat("Simulation complete! Check 'pond_ripples.png' for the final render.\n")
