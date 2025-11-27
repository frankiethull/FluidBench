# Libraries
library(rrayshader)
library(ggplot2)
library(dplyr)

# Parameters
Nx <- 200 # Grid points (x)
Ny <- 200 # Grid points (y)
Lx <- 10 # Domain size (x)
Ly <- 10 # Domain size (y)
dx <- Lx / (Nx - 1) # Spatial step
dy <- Ly / (Ny - 1)
dt <- 0.01 # Time step (CFL condition: dt < min(dx,dy)/sqrt(g*h0))
g <- 9.81 # Gravity
steps <- 150 # Simulation time steps
h0 <- 1 # Base water depth

# Initialize: h(x,y,0), hu(x,y,0), hv(x,y,0)
h <- matrix(h0, nrow = Nx, ncol = Ny)
hu <- matrix(0, nrow = Nx, ncol = Ny)
hv <- matrix(0, nrow = Nx, ncol = Ny)

# Gaussian "drop" initial condition (center perturbation)
x_grid <- seq(0, Lx, length.out = Nx)
y_grid <- seq(0, Ly, length.out = Ny)
x0 <- Lx / 2
y0 <- Ly / 2
sigma <- 0.5
h_bump <- 0.2 *
  exp(
    -((outer(x_grid, x0, "-") / sigma)^2 + (outer(y_grid, y0, "-") / sigma)^2)
  )
h <- h + t(h_bump)

# Lax-Wendroff Solver
for (step in 1:steps) {
  # 1. Compute u, v from hu, hv (avoid division by zero)
  u <- ifelse(h > 1e-6, hu / h, 0)
  v <- ifelse(h > 1e-6, hv / h, 0)

  # 2. Half-step fluxes (F_x, F_y) for h, hu, hv
  F_hx <- hu
  F_hy <- hv
  F_hux <- hu * u + 0.5 * g * h^2
  F_huy <- hu * v
  F_hvx <- hv * u
  F_hvy <- hv * v + 0.5 * g * h^2

  # 3. Half-step (predictor) via central differences
  h_half <- h -
    0.5 *
      dt *
      ((cbind(F_hx[, 2:Ny], F_hx[, Ny]) -
        cbind(F_hx[, 1], F_hx[, 1:(Ny - 1)])) /
        dx +
        (rbind(F_hy[2:Nx, ], F_hy[Nx, ]) -
          rbind(F_hy[1, ], F_hy[1:(Nx - 1), ])) /
          dy)
  hu_half <- hu -
    0.5 *
      dt *
      ((cbind(F_hux_hx[, 2:Ny], F_hux[, Ny]) -
        cbind(F_hux[, 1], F_hux[, 1:(Ny - 1)])) /
        dx +
        (rbind(F_huy[2:Nx, ], F_huy[Nx, ]) -
          rbind(F_huy[1, ], F_huy[1:(Nx - 1), ])) /
          dy)
  hv_half <- hv -
    0.5 *
      dt *
      ((cbind(F_hvx[, 2:Ny], F_hvx[, Ny]) -
        cbind(F_hvx[, 1], F_hvx[, 1:(Ny - 1)])) /
        dx +
        (rbind(F_hvy[2:Nx, ], F_hvy[Nx, ]) -
          rbind(F_hvy[1, ], F_hvy[1:(Nx - 1), ])) /
          dy)

  # 4. Update u, v for half-step
  u_half <- ifelse(h_half > 1e-6, hu_half / h_half, 0)
  v_half <- ifelse(h_half > 1e-6, hv_half / h_half, 0)

  # 5. Recompute fluxes at half-step (F^{n+1/2})
  F_hx_half <- hu_half
  F_hy_half <- hv_half
  F_hux_half <- hu_half * u_half + 0.5 * g * h_half^2
  F_huy_half <- hu_half * v_half
  F_hvx_half <- hv_half * u_half
  F_hvy_half <- hv_half * v_half + 0.5 * g * h_half^2

  # 6. Correct (full-step update)
  h_new <- h -
    dt *
      ((cbind(F_hx_half[, 2:Ny], F_hx_half[, Ny]) -
        cbind(F_hx_half[, 1], F_hx_half[, 1:(Ny - 1)])) /
        dx +
        (rbind(F_hy_half[2:Nx, ], F_hy_half[Nx, ]) -
          rbind(F_hy_half[1, ], F_hy_half[1:(Nx - 1), ])) /
          dy)
  hu_new <- hu -
    dt *
      ((cbind(F_hux_half[, 2:Ny], F_hux_half[, Ny]) -
        cbind(F_hux_half[, 1], F_hux_half[, 1:(Ny - 1)])) /
        dx +
        (rbind(F_huy_half[2:Nx, ], F_huy_half[Nx, ]) -
          rbind(F_huy_half[1, ], F_huy_half[1:(Nx - 1), ])) /
          dy)
  hv_new <- hv -
    dt *
      ((cbind(F_hvx_half[, 2:Ny], F_hvx_half[, Ny]) -
        cbind(F_hvx_half[, 1], F_hvx_half[, 1:(Ny - 1)])) /
        dx +
        (rbind(F_hvy_half[2:Nx, ], F_hvy_half[Nx, ]) -
          rbind(F_hvy_half[1, ], F_hvy_half[1:(Nx - 1), ])) /
          dy)

  # 7. Enforce boundary conditions (reflective)
  h_new[1, ] <- h_new[2, ]
  h_new[Nx, ] <- h_new[Nx - 1, ]
  h_new[, 1] <- h_new[, 2]
  h_new[, Ny] <- h_new[, Ny - 1]
  hu_new[1, ] <- -hu_new[2, ] # No-slip
  hu_new[Nx, ] <- -hu_new[Nx - 1, ]
  hu_new[, 1] <- -hu_new[, 2]
  hu_new[, Ny] <- -hu_new[, Ny - 1]
  hv_new[1, ] <- -hv_new[2, ]
  hv_new[Nx, ] <- -hv_new[Nx - 1, ]
  hv_new[, 1] <- -hv_new[, 2]
  hv_new[, Ny] <- -hv_new[, Ny - 1]

  # Update vars
  h <- h_new
  hu <- hu_new
  hv <- hv_new

  # Progress
  if (step %% 10 == 0) cat(sprintf("Step %d/%d\n", step, steps))
}

# rayshader 3D Rendering
matrix_to_raster <- function(matrix) {
  raster <- matrix %>%
    as.raster() %>%
    `attr<-`("class", c("raster", "matrix"))
  return(raster)
}

# Height matrix -> Raster -> 3D Scene
el_mat <- h - h0 # Elevation relative to base height
el_raster <- matrix_to_raster(el_mat)
scene <- height_shade(el_raster, texture = "water")
scene <- add_overlay(
  generate_2d(el_raster, texture = "water", res_factor = 2),
  scene,
  res_factor = 2
)

# Lighting: Sun + Ambient
scene <- add_lamb_shade(scene, sunaltitude = 45, sunangle = 315)
scene <- add_shadow(scene, texture = "water")

# Render 3D
render_scene(
  scene = scene,
  width = 1200,
  height = 800,
  sunangle = 315,
  sunaltitude = 45,
  fov = 70,
  zoom = 0.65,
  phi = 45,
  theta = -45,
  background = "lightblue" # Sky color
)

# Save as high-res PNG
render_snapshot(
  filename = "pond_ripple_3d.png",
  width = 2400,
  height = 1600,
  title_text = "Ripples in a Pond (2D Shallow Water Equations)"
)
