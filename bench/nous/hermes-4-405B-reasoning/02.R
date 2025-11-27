# Von Kármán Vortex Street Simulation (Lattice Boltzmann Method)
library(ggplot2)
library(gganimate)

# Simulation Parameters
nx <- 80 # Grid width
ny <- 40 # Grid height
radius <- 3 # Obstacle radius
cx <- 20 # Obstacle x-center
cy <- 20 # Obstacle y-center
omega <- 1.0 # Relaxation parameter
u_max <- 0.1 # Inflow velocity
t_max <- 10000 # Total time steps

# D2Q9 Lattice Parameters
w <- c(4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36) # Weights
cx <- c(0, 1, 0, -1, 0, 1, -1, -1, 1) # Velocity x-components
cy <- c(0, 0, 1, 0, -1, 1, 1, -1, -1) # Velocity y-components

# Initialize distribution functions and macroscopic variables
f <- array(1 / 9, dim = c(nx, ny, 9))
rho <- matrix(1, nrow = nx, ncol = ny)
u <- matrix(0, nrow = nx, ncol = ny)
v <- matrix(0, nrow = nx, ncol = ny)

# Create obstacle mask
obstacle <- matrix(FALSE, nrow = nx, ncol = ny)
for (i in seq(nx)) {
  for (j in seq(ny)) {
    if ((i - cx)^2 + (j - cy)^2 <= radius^2) {
      obstacle[i, j] <- TRUE
    }
  }
}

# Lattice Boltzmann Functions
collide <- function(f, rho, u, v) {
  # Compute equilibrium distribution
  feq <- array(0, dim = c(nx, ny, 9))
  for (k in 1:9) {
    u2 <- u^2 + v^2
    cu <- cx[k] * u + cy[k] * v
    feq[,, k] <- w[k] * rho * (1 + 3 * cu + 4.5 * cu^2 - 1.5 * u2)
  }

  # Collision step
  f <- f - omega * (f - feq)
  return(f)
}

stream <- function(f) {
  # Streaming step with bounce-back BCs
  f_new <- array(0, dim = c(nx, ny, 9))

  # Vectorized streaming for interior points
  for (k in 1:9) {
    f_new[,, k] <- shift(f[,, k], cx[k], cy[k], obstacle)
  }

  return(f_new)
}

shift <- function(f_k, dx, dy, obstacle) {
  # Apply periodic BCs with obstacle handling
  f_shifted <- f_k[c((1:nx) - dx, nx)[1:nx], c((1:ny) - dy, ny)[1:ny]]
  f_shifted[obstacle] <- f_k[obstacle] # Bounce-back
  return(f_shifted)
}

# Main Simulation Loop
for (t in 1:t_max) {
  # Inflow boundary condition (left edge)
  f[1, , c(2, 6, 9)] <- w[c(2, 6, 9)] * (1 + 3 * u_max) # Zou-He velocity BC

  # Collision
  f <- collide(f, rho, u, v)

  # Streaming
  f <- stream(f)

  # Calculate macroscopic variables
  rho <- apply(f, c(1, 2), sum)
  u <- (apply(apply(f, c(1, 2, 3), function(x) x * cx), c(1, 2), sum)) / rho
  v <- (apply(apply(f, c(1, 2, 3), function(x) x * cy), c(1, 2), sum)) / rho

  # Calculate vorticity
  vorticity <- vorticity_field(u, v)

  # Animation frame every 50 steps
  if (t %% 50 == 0) {
    plot_frame(vorticity, t)
  }
}

# Vorticity Calculation
vorticity_field <- function(u, v) {
  # Simple central difference approximation
  dvdx <- (v[, -1] - v[, -ncol(v)])[-1, ]
  dudy <- (u[-1, ] - u[-nrow(u), ])[, -1]
  return(dvdx - dudy)
}

# Plotting Function
plot_frame <- function(vort, step) {
  df <- expand.grid(x = 1:nx, y = 1:ny)
  df$vort <- as.vector(vort)

  p <- ggplot(df, aes(x, y, fill = vort)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-0.1, 0.1)
    ) +
    labs(
      title = paste("Von Kármán Vortex Street - Step", step),
      x = "X Position",
      y = "Y Position"
    ) +
    theme_minimal()

  # For gganimate:
  # anim <- p + transition_manual(step)
  # anim_save("vortex.gif", animation = anim)

  # For frame-by-frame saving:
  ggsave(
    paste0("frame_", sprintf("%04d", step), ".png"),
    plot = p,
    width = 8,
    height = 4,
    dpi = 100
  )
}

# Create animation from frames (requires ImageMagick)
# system("convert -delay 5 frame_*.png vortex.gif")
