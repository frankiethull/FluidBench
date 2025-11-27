# Load required packages
library(ggplot2)
library(gganimate)
library(dplyr)
library(purrr)

# Simulation parameters
nx <- 80    # Grid width
ny <- 40    # Grid height
tau <- 0.6  # Relaxation time
max_steps <- 5000  # Maximum number of time steps
save_every <- 20   # Save frame every n steps

# D2Q9 model parameters
weights <- c(4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36)
cx <- c(0, 1, 0, -1,  0, 1, -1, -1,  1)  # x velocities
cy <- c(0, 0, 1,  0, -1, 1,  1, -1, -1)  # y velocities
opposite <- c(1, 4, 5, 2, 3, 8, 9, 6, 7)  # Opposite directions

# Initialize density distributions
f <- array(0, dim = c(nx, ny, 9))
feq <- array(0, dim = c(nx, ny, 9))

# Initialize macroscopic variables
rho <- matrix(1, nrow = nx, ncol = ny)  # Density
ux <- matrix(0, nrow = nx, ncol = ny)   # x-velocity
uy <- matrix(0, nrow = nx, ncol = ny)   # y-velocity

# Circular obstacle parameters
obstacle_x <- 20
obstacle_y <- 20
radius <- 3

# Create obstacle mask
obstacle <- matrix(FALSE, nrow = nx, ncol = ny)
for (i in 1:nx) {
  for (j in 1:ny) {
    if ((i - obstacle_x)^2 + (j - obstacle_y)^2 <= radius^2) {
      obstacle[i, j] <- TRUE
    }
  }
}

# Function to calculate equilibrium distribution
calculate_feq <- function(rho, ux, uy) {
  usqr <- ux^2 + uy^2
  feq <- array(0, dim = c(nx, ny, 9))

  for (k in 1:9) {
    cu <- 3 * (cx[k] * ux + cy[k] * uy)
    feq[,,k] <- weights[k] * rho * (1 + cu + 0.5 * cu^2 - 1.5 * usqr)
  }

  return(feq)
}

# Function to calculate vorticity
calculate_vorticity <- function(ux, uy) {
  # Calculate gradients using central differences
  dudy <- matrix(0, nrow = nx, ncol = ny)
  dvdx <- matrix(0, nrow = nx, ncol = ny)

  # Inner points
  dudy[2:(nx-1), 2:(ny-1)] <- (uy[2:(nx-1), 3:ny] - uy[2:(nx-1), 1:(ny-2)]) / 2
  dvdx[2:(nx-1), 2:(ny-1)] <- (ux[3:nx, 2:(ny-1)] - ux[1:(nx-2), 2:(ny-1)]) / 2

  # Boundary conditions (simple forward/backward differences)
  dudy[1, ] <- uy[1, 2:ny] - uy[1, 1:(ny-1)]
  dudy[nx, ] <- uy[nx, 2:ny] - uy[nx, 1:(ny-1)]
  dudy[, 1] <- uy[2:nx, 1] - uy[1:(nx-1), 1]
  dudy[, ny] <- uy[2:nx, ny] - uy[1:(nx-1), ny]

  dvdx[1, ] <- ux[2, ] - ux[1, ]
  dvdx[nx, ] <- ux[nx, ] - ux[nx-1, ]
  dvdx[, 1] <- ux[, 2] - ux[, 1]
  dvdx[, ny] <- ux[, ny] - ux[, ny-1]

  return(dvdx - dudy)
}

# Main simulation loop
frames <- list()
step <- 0

while (step < max_steps) {
  step <- step + 1

  # Boundary conditions (inflow from left)
  rho[1, ] <- 1
  ux[1, ] <- 0.1  # Inflow velocity
  uy[1, ] <- 0

  # Bounce-back boundary condition for obstacle
  for (i in 1:nx) {
    for (j in 1:ny) {
      if (obstacle[i, j]) {
        for (k in 1:9) {
          f[i, j, k] <- f[i, j, opposite[k]]
        }
      }
    }
  }

  # Calculate macroscopic variables
  rho <- apply(f, c(1, 2), sum)
  ux <- apply(f * cx, c(1, 2), sum) / rho
  uy <- apply(f * cy, c(1, 2), sum) / rho

  # Apply boundary conditions to velocity
  ux[obstacle] <- 0
  uy[obstacle] <- 0

  # Calculate equilibrium distribution
  feq <- calculate_feq(rho, ux, uy)

  # Collision step
  f <- f - (1/tau) * (f - feq)

  # Streaming step (vectorized)
  for (k in 1:9) {
    f[,,k] <- shift(f[,,k], c(cx[k], cy[k]))
  }

  # Calculate vorticity
  vorticity <- calculate_vorticity(ux, uy)

  # Save frame periodically
  if (step %% save_every == 0) {
    # Create data frame for visualization
    df <- expand.grid(x = 1:nx, y = 1:ny)
    df$vorticity <- as.vector(vorticity)
    df$obstacle <- as.vector(obstacle)
    df$step <- step

    frames[[length(frames) + 1]] <- df
    cat("Step:", step, "\n")
  }
}

# Combine all frames
all_frames <- bind_rows(frames)

# Create animation
p <- ggplot(all_frames, aes(x = x, y = y, fill = vorticity)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limit = c(-0.2, 0.2), name = "Vorticity") +
  geom_tile(data = filter(all_frames, obstacle), aes(x = x, y = y), fill = "black") +
  labs(title = 'Von Kármán Vortex Street (Step: {frame_time})') +
  theme_void() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  transition_time(step) +
  ease_aes('linear')

# Render animation
animate(p, fps = 15, width = 800, height = 400, renderer = gifski_renderer())