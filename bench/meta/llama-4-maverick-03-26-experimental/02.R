# Compute equilibrium distribution
equilibrium <- function(rho, ux, uy, w_i, c_i, cu_factor, cu2_factor) {
  cu <- (ux %o% c_i[1, ]) + (uy %o% c_i[2, ])
  feq <- array(rep(0, length(rho) * 9), dim = c(dim(rho), 9))
  for (i in 1:9) {
    feq[, , i] <- w_i[i] * rho * (1 + cu_factor * cu[, , i] + 
                                   cu2_factor * cu[, , i]^2 - 
                                   cu2_factor * (ux^2 + uy^2))
  }
  feq
}

# Shift matrix (periodic boundaries)
shift_matrix <- function(mat, dx, dy) {
  Nx <- nrow(mat)
  Ny <- ncol(mat)
  x_new <- (1:Nx + dx) %% Nx
  x_new[x_new == 0] <- Nx
  y_new <- (1:Ny + dy) %% Ny
  y_new[y_new == 0] <- Ny
  mat[x_new, y_new]
}

# Libraries
library(ggplot2)
library(gganimate)
library(dplyr)

# Lattice Boltzmann D2Q9 Parameters
Nx <- 80  # Lattice width
Ny <- 40  # Lattice height
rho0 <- 1 # Background density
u0 <- 0.1 # Inflow velocity (Left to Right)
tau <- 0.6 # Relaxation time (viscosity = (tau - 0.5)/3)
R <- 3    # Obstacle (cylinder) radius
x_c <- 20 # Obstacle center x
y_c <- 20 # Obstacle center y
steps <- 2000 # Simulation time steps

# D2Q9 Lattice Velocities and Weights
c_i <- matrix(c(
  0,  1,  0, -1,  0,  1, -1, -1,  1, # x-components
  0,  0,  1,  0, -1,  1,  1, -1, -1  # y-components
), nrow = 2, byrow = TRUE)
w_i <- c(4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36)

# Initialize distributions: f_i(x,y,t=0)
f <- array(rep(w_i * rho0, Nx * Ny), dim = c(Nx, Ny, 9))

# Precompute equilibrium helper terms
cu2_factor <- 3/2 * (1 / (rho0 * u0)^2)
cu_factor <- 3 / rho0

# Main LBM Loop
vorticity_frames <- list()
for (step in 1:steps) {
  # 0. Compute macroscopic vars: rho(x,y), u(x,y)
  rho <- apply(f, 1:2, sum)
  ux <- (apply(f * c_i[1, ], 1:2, sum)) / rho
  uy <- (apply(f * c_i[2, ], 1:2, sum)) / rho
  
  # Boundary: Inflow (Left) - Fixed velocity u0
  rho[:, 1] <- 1 / (1 - ux[, 1]) * (apply(f[, 1, c(1, 3, 4, 7, 8)], 1, sum) + 
                                    2 * apply(f[, 1, c(2, 5, 6)], 1, sum))
  ux[, 1] <- u0
  uy[, 1] <- 0
  feq <- equilibrium(rho[, 1], ux[, 1], uy[, 1], w_i, c_i, cu_factor, cu2_factor)
  f[, 1, 2] <- feq[, 2] - f[, 1, 4] + feq[, 4]
  f[, 1, 5] <- feq[, 5] - f[, 1, 7] + feq[, 7]
  f[, 1, 6] <- feq[, 6] - f[, 1, 8] + feq[, 8]
  
  # 1. Collision Step (relax toward equilibrium)
  feq <- equilibrium(rho, ux, uy, w_i, c_i, cu_factor, cu2_factor)
  f <- f - (1 / tau) * (f - feq)
  
  # 2. Bounce-back at obstacle (circular cylinder)
  solid_nodes <- (row(f[, , 1]) - x_c)^2 + (col(f[, , 1]) - y_c)^2 <= R^2
  for (i in 1:9) {
    i_opposite <- c(1, 3, 4, 2, 5, 7, 8, 6, 9)[i]
    f[solid_nodes, i] <- f[solid_nodes, i_opposite]
  }
  
  # 3. Streaming Step (vectorized)
  f_new <- f
  for (i in 2:9) { # Skip rest (c_0 = (0,0))
    x_shift <- c_i[1, i]
    y_shift <- c_i[2, i]
    f_new[, , i] <- shift_matrix(f[, , i], x_shift, y_shift)
  }
  f <- f_new
  
  # 4. Outflow (Right) - Zero-gradient
  f[Nx, , ] <- f[Nx - 1, , ]
  
  # 5. Compute vorticity (curl) for visualization
  vorticity <- (c(ux[-1, ] - ux[-Nx, ]) - c(uy[, -1] - uy[, -Ny])) / 2
  vorticity_frames[[step]] <- data.frame(
    x = rep(1:Nx, Ny),
    y = rep(1:Ny, each = Nx),
    vorticity = as.vector(vorticity)
  )
  
  # Progress bar
  if (step %% 100 == 0) cat(sprintf("Step %d/%d\n", step, steps))
}

# Animate vorticity frames using gganimate
vorticity_df <- bind_rows(vorticity_frames, .id = "step")
ggplot(vorticity_df, aes(x, y, fill = vorticity)) +
  geom_raster() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Von Kármán Vortex Street (Vorticity)", x = "X", y = "Y", fill = "Vorticity") +
  theme_void() + theme(legend.position = "none") +
  transition_manual(step) +
  ease_aes('linear') -> vorticity_animation

# Save animation as GIF
anim_save("vkarman_vorticity.gif", vorticity_animation, 
          nframes = steps, fps = 30, width = 800, height = 400)