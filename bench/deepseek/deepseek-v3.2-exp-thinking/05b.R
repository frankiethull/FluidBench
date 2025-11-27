library(ggplot2)
library(plotly)

# Improved Rayleigh-Bénard simulation using streamfunction-vorticity formulation
simulate_rb_convection_improved <- function() {
  # Parameters
  Ra <- 2000
  Pr <- 0.71
  nx <- 100
  ny <- 50
  Lx <- 4
  Ly <- 1
  dx <- Lx / (nx - 1)
  dy <- Ly / (ny - 1)
  dt <- 0.001
  nsteps <- 5000

  # Initialize fields
  psi <- matrix(0, nrow = nx, ncol = ny) # Streamfunction
  omega <- matrix(0, nrow = nx, ncol = ny) # Vorticity
  T <- matrix(0, nrow = nx, ncol = ny) # Temperature

  # Boundary conditions
  T[, 1] <- 1 # Hot bottom
  T[, ny] <- 0 # Cold top

  # Add perturbation
  set.seed(123)
  for (i in 1:nx) {
    for (j in 1:ny) {
      x <- (i - 1) * dx
      y <- (j - 1) * dy
      T[i, j] <- T[i, j] + 0.01 * sin(pi * x / Lx) * sin(pi * y / Ly)
    }
  }

  # Main simulation loop
  for (step in 1:nsteps) {
    # Solve for streamfunction: ∇²ψ = -ω
    for (iter in 1:20) {
      psi_old <- psi
      for (i in 2:(nx - 1)) {
        for (j in 2:(ny - 1)) {
          psi[i, j] <- 0.25 *
            (psi[i + 1, j] +
              psi[i - 1, j] +
              psi[i, j + 1] +
              psi[i, j - 1] +
              dx * dy * omega[i, j])
        }
      }
      # Convergence check
      if (max(abs(psi - psi_old)) < 1e-6) break
    }

    # Compute velocities from streamfunction
    u <- matrix(0, nrow = nx, ncol = ny)
    v <- matrix(0, nrow = nx, ncol = ny)
    for (i in 2:(nx - 1)) {
      for (j in 2:(ny - 1)) {
        u[i, j] <- (psi[i, j + 1] - psi[i, j - 1]) / (2 * dy)
        v[i, j] <- -(psi[i + 1, j] - psi[i - 1, j]) / (2 * dx)
      }
    }

    # Update vorticity equation
    omega_new <- omega
    for (i in 2:(nx - 1)) {
      for (j in 2:(ny - 1)) {
        # Advection
        u_domega_dx <- u[i, j] * (omega[i + 1, j] - omega[i - 1, j]) / (2 * dx)
        v_domega_dy <- v[i, j] * (omega[i, j + 1] - omega[i, j - 1]) / (2 * dy)

        # Diffusion
        laplacian_omega <- (omega[i + 1, j] +
          omega[i - 1, j] -
          2 * omega[i, j]) /
          (dx^2) +
          (omega[i, j + 1] + omega[i, j - 1] - 2 * omega[i, j]) / (dy^2)

        # Buoyancy term
        dT_dx <- (T[i + 1, j] - T[i - 1, j]) / (2 * dx)

        omega_new[i, j] <- omega[i, j] +
          dt *
            (-u_domega_dx -
              v_domega_dy +
              sqrt(Pr / Ra) * laplacian_omega +
              Ra * Pr * dT_dx)
      }
    }
    omega <- omega_new

    # Update temperature equation
    T_new <- T
    for (i in 2:(nx - 1)) {
      for (j in 2:(ny - 1)) {
        # Advection
        u_dT_dx <- u[i, j] * (T[i + 1, j] - T[i - 1, j]) / (2 * dx)
        v_dT_dy <- v[i, j] * (T[i, j + 1] - T[i, j - 1]) / (2 * dy)

        # Diffusion
        laplacian_T <- (T[i + 1, j] + T[i - 1, j] - 2 * T[i, j]) /
          (dx^2) +
          (T[i, j + 1] + T[i, j - 1] - 2 * T[i, j]) / (dy^2)

        T_new[i, j] <- T[i, j] +
          dt * (-u_dT_dx - v_dT_dy + (1 / sqrt(Ra * Pr)) * laplacian_T)
      }
    }
    T <- T_new

    # Reapply boundary conditions
    T[, 1] <- 1
    T[, ny] <- 0
    T[1, ] <- T[2, ]
    T[nx, ] <- T[nx - 1, ] # Insulating sides

    if (step %% 500 == 0) cat("Step:", step, "/", nsteps, "\n")
  }

  return(list(T = T, u = u, v = v, psi = psi))
}

# Run improved simulation
cat("Running improved Rayleigh-Bénard simulation...\n")
result_improved <- simulate_rb_convection_improved()

# Create enhanced visualization
create_enhanced_visualization <- function(T, u, v) {
  nx <- nrow(T)
  ny <- ncol(T)
  Lx <- 4
  Ly <- 1

  # Create data frames
  x_coords <- seq(0, Lx, length.out = nx)
  y_coords <- seq(0, Ly, length.out = ny)

  temp_data <- expand.grid(x = x_coords, y = y_coords)
  temp_data$Temperature <- as.vector(T)

  # Create vector field (subsampled)
  subsample <- 6
  vector_data <- expand.grid(
    x = x_coords[seq(1, nx, by = subsample)],
    y = y_coords[seq(1, ny, by = subsample)]
  )
  vector_data$u <- as.vector(u[
    seq(1, nx, by = subsample),
    seq(1, ny, by = subsample)
  ])
  vector_data$v <- as.vector(v[
    seq(1, nx, by = subsample),
    seq(1, ny, by = subsample)
  ])

  # Filter small vectors
  speed <- sqrt(vector_data$u^2 + vector_data$v^2)
  vector_data <- vector_data[speed > max(speed) * 0.1, ]

  # Create plot
  p <- ggplot(temp_data, aes(x = x, y = y)) +
    geom_raster(aes(fill = Temperature), interpolate = TRUE) +
    geom_segment(
      data = vector_data,
      aes(x = x, y = y, xend = x + u * 0.2, yend = y + v * 0.2),
      arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
      color = "white",
      alpha = 0.8,
      linewidth = 0.4
    ) +
    scale_fill_viridis_c(option = "inferno") +
    coord_fixed(ratio = Ly / Lx) +
    labs(
      title = "Rayleigh-Bénard Convection Rolls",
      subtitle = "Hot fluid rises, cold fluid sinks forming convection cells",
      x = "X",
      y = "Y"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

  return(p)
}

# Create and save enhanced visualization
p_enhanced <- create_enhanced_visualization(
  result_improved$T,
  result_improved$u,
  result_improved$v
)
print(p_enhanced)
ggsave(
  "rayleigh_benard_enhanced.png",
  p_enhanced,
  width = 12,
  height = 4,
  dpi = 300
)
