library(ggplot2)
library(gganimate)
library(viridis)

# D2Q9 Lattice Boltzmann parameters
nx <- 80
ny <- 40
omega <- 1.0 # Relaxation parameter
tau <- 1 / omega
cs2 <- 1 / 3 # Speed of sound squared

# Lattice velocities (D2Q9)
cx <- c(0, 1, 0, -1, 0, 1, -1, -1, 1)
cy <- c(0, 0, 1, 0, -1, 1, 1, -1, -1)
w <- c(4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36) # weights
opposite <- c(1, 4, 5, 2, 3, 8, 9, 6, 7) # Opposite directions

# Initialize distribution functions
f <- array(0, dim = c(nx, ny, 9))
feq <- array(0, dim = c(nx, ny, 9))
rho <- matrix(1, nrow = nx, ncol = ny) # Density
ux <- matrix(0, nrow = nx, ncol = ny) # x-velocity
uy <- matrix(0, nrow = nx, ncol = ny) # y-velocity

# Create circular obstacle
obstacle <- matrix(FALSE, nrow = nx, ncol = ny)
cx_obs <- 20
cy_obs <- 20
radius <- 3

for (i in 1:nx) {
  for (j in 1:ny) {
    if ((i - cx_obs)^2 + (j - cy_obs)^2 <= radius^2) {
      obstacle[i, j] <- TRUE
    }
  }
}

# Inlet velocity profile (parabolic)
inlet_velocity <- 0.1
for (j in 1:ny) {
  ux[1, j] <- inlet_velocity * 4 * (j / ny) * (1 - j / ny)
}

# Initialize equilibrium distribution
for (k in 1:9) {
  for (i in 1:nx) {
    for (j in 1:ny) {
      if (!obstacle[i, j]) {
        u2 <- ux[i, j]^2 + uy[i, j]^2
        cu <- cx[k] * ux[i, j] + cy[k] * uy[i, j]
        feq[i, j, k] <- w[k] * rho[i, j] * (1 + 3 * cu + 4.5 * cu^2 - 1.5 * u2)
        f[i, j, k] <- feq[i, j, k]
      }
    }
  }
}

# Function to compute macroscopic variables
compute_macroscopic <- function(f) {
  rho_new <- apply(f, c(1, 2), sum)
  ux_new <- apply(sweep(f, 3, cx, "*"), c(1, 2), sum) / rho_new
  uy_new <- apply(sweep(f, 3, cy, "*"), c(1, 2), sum) / rho_new
  list(rho = rho_new, ux = ux_new, uy = uy_new)
}

# Function to compute vorticity
compute_vorticity <- function(ux, uy) {
  vort <- matrix(0, nrow = nx, ncol = ny)
  for (i in 2:(nx - 1)) {
    for (j in 2:(ny - 1)) {
      if (!obstacle[i, j]) {
        dudy <- (uy[i + 1, j] - uy[i - 1, j]) / 2
        dvdx <- (ux[i, j + 1] - ux[i, j - 1]) / 2
        vort[i, j] <- dvdx - dudy
      }
    }
  }
  vort
}

# Simulation parameters
nsteps <- 1000
save_every <- 20
frames <- list()
frame_count <- 1

# Main simulation loop
for (step in 1:nsteps) {
  # Collision step
  for (k in 1:9) {
    for (i in 1:nx) {
      for (j in 1:ny) {
        if (!obstacle[i, j]) {
          u2 <- ux[i, j]^2 + uy[i, j]^2
          cu <- cx[k] * ux[i, j] + cy[k] * uy[i, j]
          feq[i, j, k] <- w[k] *
            rho[i, j] *
            (1 + 3 * cu + 4.5 * cu^2 - 1.5 * u2)
          f[i, j, k] <- f[i, j, k] - omega * (f[i, j, k] - feq[i, j, k])
        }
      }
    }
  }

  # Streaming step (vectorized)
  f_stream <- f
  for (k in 1:9) {
    # Create shifted matrices for streaming
    i_shift <- ((1:nx) - cx[k] - 1) %% nx + 1
    j_shift <- ((1:ny) - cy[k] - 1) %% ny + 1

    f_stream[,, k] <- f[i_shift, j_shift, k]
  }
  f <- f_stream

  # Boundary conditions - bounce back at obstacle
  for (k in 1:9) {
    f[,, k][obstacle] <- f[,, opposite[k]][obstacle]
  }

  # Inlet boundary condition (constant velocity)
  for (j in 1:ny) {
    if (!obstacle[1, j]) {
      rho_in <- (f[1, j, 1] +
        f[1, j, 2] +
        f[1, j, 3] +
        f[1, j, 4] +
        f[1, j, 5] +
        f[1, j, 6] +
        f[1, j, 7] +
        f[1, j, 8] +
        f[1, j, 9])
      for (k in 1:9) {
        cu <- cx[k] * ux[1, j] + cy[k] * uy[1, j]
        f[1, j, k] <- feq[1, j, k] + w[k] * rho_in * 3 * cu
      }
    }
  }

  # Update macroscopic variables
  macro <- compute_macroscopic(f)
  rho <- macro$rho
  ux <- macro$ux
  uy <- macro$uy

  # Compute vorticity
  vorticity <- compute_vorticity(ux, uy)

  # Save frame periodically
  if (step %% save_every == 0) {
    # Create data frame for plotting
    plot_data <- expand.grid(x = 1:nx, y = 1:ny)
    plot_data$vorticity <- as.vector(vorticity)
    plot_data$step <- step
    plot_data$obstacle <- as.vector(obstacle)

    frames[[frame_count]] <- plot_data
    frame_count <- frame_count + 1

    cat("Step:", step, "\n")
  }
}

# Combine all frames
all_frames <- do.call(rbind, frames)

# Create animation
p <- ggplot(all_frames, aes(x = x, y = y, fill = vorticity)) +
  geom_tile() +
  geom_tile(
    data = all_frames[all_frames$obstacle, ],
    aes(x = x, y = y),
    fill = "black"
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-0.02, 0.02)
  ) +
  coord_fixed() +
  labs(
    title = "Von Kármán Vortex Street\nStep: {frame_time}",
    x = "X",
    y = "Y",
    fill = "Vorticity"
  ) +
  theme_minimal() +
  transition_time(step)

# Animate and save
animate(p, nframes = length(frames), fps = 10, width = 800, height = 400)
anim_save("vortex_street.gif")
