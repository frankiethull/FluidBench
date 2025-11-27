library(ggplot2)
library(plotly)
library(viridis)
library(reshape2)
library(deSolve)

# Parameters for Rayleigh-Bénard Convection
Pr <- 0.71 # Prandtl number (air)
Ra <- 2000 # Rayleigh number (moderate convection)
aspect_ratio <- 2

# Grid parameters
nx <- 80
ny <- 40
Lx <- 2 * aspect_ratio
Ly <- 2
dx <- Lx / (nx - 1)
dy <- Ly / (ny - 1)

# Precompute coefficients
alpha <- 1.0 # Thermal expansion coefficient
g <- 9.81 # Gravity

# Initialize fields
T <- matrix(0, nrow = nx, ncol = ny) # Temperature
u <- matrix(0, nrow = nx, ncol = ny) # x-velocity
v <- matrix(0, nrow = nx, ncol = ny) # y-velocity
p <- matrix(0, nrow = nx, ncol = ny) # Pressure

# Set boundary conditions
# Temperature: Bottom hot (T=1), Top cold (T=0)
T[, 1] <- 1 # Bottom wall
T[, ny] <- 0 # Top wall

# Velocity: No-slip at all boundaries
u[, 1] <- 0
u[, ny] <- 0
u[1, ] <- 0
u[nx, ] <- 0
v[, 1] <- 0
v[, ny] <- 0
v[1, ] <- 0
v[nx, ] <- 0

# Add small perturbation to break symmetry
set.seed(123)
perturbation <- matrix(rnorm(nx * ny, 0, 0.01), nrow = nx, ncol = ny)
T <- T + perturbation
T[, 1] <- 1 # Reapply BC after perturbation
T[, ny] <- 0

# Function to compute derivatives using finite differences
compute_derivatives <- function(field, dx, dy) {
  dfdx <- matrix(0, nrow = nx, ncol = ny)
  dfdy <- matrix(0, nrow = nx, ncol = ny)
  laplacian <- matrix(0, nrow = nx, ncol = ny)

  # Interior points
  for (i in 2:(nx - 1)) {
    for (j in 2:(ny - 1)) {
      dfdx[i, j] <- (field[i + 1, j] - field[i - 1, j]) / (2 * dx)
      dfdy[i, j] <- (field[i, j + 1] - field[i, j - 1]) / (2 * dy)
      laplacian[i, j] <- (field[i + 1, j] + field[i - 1, j] - 2 * field[i, j]) /
        (dx^2) +
        (field[i, j + 1] + field[i, j - 1] - 2 * field[i, j]) / (dy^2)
    }
  }

  return(list(dfdx = dfdx, dfdy = dfdy, laplacian = laplacian))
}

# Streamfunction-vorticity formulation for Rayleigh-Bénard
simulate_convection <- function(nsteps, dt) {
  # Storage for results
  results <- list()

  for (step in 1:nsteps) {
    # Compute derivatives
    T_derivs <- compute_derivatives(T, dx, dy)
    u_derivs <- compute_derivatives(u, dx, dy)
    v_derivs <- compute_derivatives(v, dx, dy)

    # Temporary arrays for updates
    T_new <- T
    u_new <- u
    v_new <- v

    # Update equations (Boussinesq approximation)
    for (i in 2:(nx - 1)) {
      for (j in 2:(ny - 1)) {
        # Temperature equation (advection-diffusion)
        T_new[i, j] <- T[i, j] +
          dt *
            (-u[i, j] *
              T_derivs$dfdx[i, j] -
              v[i, j] * T_derivs$dfdy[i, j] +
              (1 / sqrt(Pr * Ra)) * T_derivs$laplacian[i, j])

        # x-momentum equation
        u_new[i, j] <- u[i, j] +
          dt *
            (-u[i, j] *
              u_derivs$dfdx[i, j] -
              v[i, j] * u_derivs$dfdy[i, j] +
              sqrt(Pr / Ra) * u_derivs$laplacian[i, j] -
              (p[i + 1, j] - p[i - 1, j]) / (2 * dx))

        # y-momentum equation with Boussinesq buoyancy term
        v_new[i, j] <- v[i, j] +
          dt *
            (
              -u[i, j] *
                v_derivs$dfdx[i, j] -
                v[i, j] * v_derivs$dfdy[i, j] +
                sqrt(Pr / Ra) * v_derivs$laplacian[i, j] -
                (p[i, j + 1] - p[i, j - 1]) / (2 * dy) +
                alpha * g * T[i, j] # Buoyancy term
            )
      }
    }

    # Pressure correction (simplified)
    p_derivs <- compute_derivatives(p, dx, dy)
    for (i in 2:(nx - 1)) {
      for (j in 2:(ny - 1)) {
        divergence <- (u_new[i + 1, j] - u_new[i - 1, j]) /
          (2 * dx) +
          (v_new[i, j + 1] - v_new[i, j - 1]) / (2 * dy)
        p[i, j] <- p[i, j] - 0.1 * divergence # Relaxation
      }
    }

    # Update fields
    T <- T_new
    u <- u_new
    v <- v_new

    # Reapply boundary conditions
    T[, 1] <- 1
    T[, ny] <- 0
    T[1, ] <- (T[2, ] + T[1, ]) / 2
    T[nx, ] <- (T[nx - 1, ] + T[nx, ]) / 2 # Insulating sides

    u[, 1] <- 0
    u[, ny] <- 0
    u[1, ] <- 0
    u[nx, ] <- 0
    v[, 1] <- 0
    v[, ny] <- 0
    v[1, ] <- 0
    v[nx, ] <- 0

    # Store results periodically
    if (step %% 100 == 0) {
      results[[length(results) + 1]] <- list(
        step = step,
        T = T,
        u = u,
        v = v
      )
      cat("Step:", step, "/", nsteps, "\n")
    }
  }

  return(results)
}

# Run simulation
cat("Starting Rayleigh-Bénard convection simulation...\n")
nsteps <- 2000
dt <- 0.01
results <- simulate_convection(nsteps, dt)

# Get final state
final_state <- results[[length(results)]]
T_final <- final_state$T
u_final <- final_state$u
v_final <- final_state$v

# Create data frame for ggplot2 visualization
create_plot_data <- function(T, u, v, subsample = 4) {
  # Create coordinate grids
  x_coords <- seq(0, Lx, length.out = nx)
  y_coords <- seq(0, Ly, length.out = ny)

  # Create temperature data
  temp_data <- expand.grid(x = x_coords, y = y_coords)
  temp_data$Temperature <- as.vector(T)

  # Create vector field data (subsampled for clarity)
  vector_data <- expand.grid(
    x = x_coords[seq(1, nx, by = subsample)],
    y = y_coords[seq(1, ny, by = subsample)]
  )

  # Get corresponding velocity components
  vector_data$u <- as.vector(u[
    seq(1, nx, by = subsample),
    seq(1, ny, by = subsample)
  ])
  vector_data$v <- as.vector(v[
    seq(1, nx, by = subsample),
    seq(1, ny, by = subsample)
  ])

  # Filter out very small vectors for clarity
  speed <- sqrt(vector_data$u^2 + vector_data$v^2)
  vector_data <- vector_data[speed > 0.001, ]

  return(list(temp_data = temp_data, vector_data = vector_data))
}

# Create plot data
plot_data <- create_plot_data(T_final, u_final, v_final, subsample = 3)

# Create static ggplot2 visualization
p_static <- ggplot(plot_data$temp_data, aes(x = x, y = y)) +
  geom_raster(aes(fill = Temperature)) +
  geom_segment(
    data = plot_data$vector_data,
    aes(x = x, y = y, xend = x + u * 0.5, yend = y + v * 0.5),
    arrow = arrow(length = unit(0.1, "cm")),
    color = "white",
    alpha = 0.7,
    linewidth = 0.3
  ) +
  scale_fill_viridis_c(option = "plasma") +
  coord_fixed(ratio = 1) +
  labs(
    title = "Rayleigh-Bénard Convection",
    subtitle = "Temperature (heatmap) + Velocity (vectors)",
    x = "X",
    y = "Y"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_static)

# Save static plot
ggsave(
  "rayleigh_benard_static.png",
  p_static,
  width = 10,
  height = 6,
  dpi = 300
)

# Create interactive plotly visualization
create_interactive_plot <- function(T, u, v, subsample = 6) {
  plot_data <- create_plot_data(T, u, v, subsample)

  # Create the base temperature plot
  p <- plot_ly() %>%
    add_heatmap(
      x = unique(plot_data$temp_data$x),
      y = unique(plot_data$temp_data$y),
      z = matrix(plot_data$temp_data$Temperature, nrow = nx),
      colorscale = "Viridis",
      showscale = TRUE,
      name = "Temperature"
    ) %>%
    add_annotations(
      x = 0.5,
      y = -0.1,
      text = "Hot Bottom",
      xref = "paper",
      yref = "paper",
      showarrow = FALSE
    ) %>%
    add_annotations(
      x = 0.5,
      y = 1.1,
      text = "Cold Top",
      xref = "paper",
      yref = "paper",
      showarrow = FALSE
    )

  # Add vector field
  p <- p %>%
    add_annotations(
      x = plot_data$vector_data$x,
      y = plot_data$vector_data$y,
      ax = plot_data$vector_data$u * 0.3,
      ay = plot_data$vector_data$v * 0.3,
      xref = "x",
      yref = "y",
      axref = "x",
      ayref = "y",
      showarrow = TRUE,
      arrowhead = 2,
      arrowsize = 0.5,
      arrowwidth = 1.5,
      arrowcolor = "white"
    )

  p <- p %>%
    layout(
      title = "Rayleigh-Bénard Convection - Interactive",
      xaxis = list(title = "X", scaleanchor = "y"),
      yaxis = list(title = "Y"),
      showlegend = FALSE
    )

  return(p)
}

# Create interactive plot
p_interactive <- create_interactive_plot(T_final, u_final, v_final)
p_interactive

# Save interactive plot
htmlwidgets::saveWidget(p_interactive, "rayleigh_benard_interactive.html")

# Alternative: Create animation showing evolution
if (length(results) > 5) {
  # Select a subset of time steps for animation
  selected_steps <- round(seq(1, length(results), length.out = 6))

  animation_data <- list()
  for (i in seq_along(selected_steps)) {
    state <- results[[selected_steps[i]]]
    plot_data <- create_plot_data(state$T, state$u, state$v, subsample = 4)
    plot_data$temp_data$step <- selected_steps[i]
    plot_data$vector_data$step <- selected_steps[i]
    animation_data[[i]] <- plot_data
  }

  # Combine all frames
  all_temp_data <- do.call(
    rbind,
    lapply(animation_data, function(x) x$temp_data)
  )
  all_vector_data <- do.call(
    rbind,
    lapply(animation_data, function(x) x$vector_data)
  )

  # Create faceted animation plot
  p_animation <- ggplot(all_temp_data, aes(x = x, y = y)) +
    geom_raster(aes(fill = Temperature)) +
    geom_segment(
      data = all_vector_data,
      aes(x = x, y = y, xend = x + u * 0.4, yend = y + v * 0.4),
      arrow = arrow(length = unit(0.08, "cm")),
      color = "white",
      alpha = 0.7,
      linewidth = 0.2
    ) +
    scale_fill_viridis_c(option = "plasma") +
    facet_wrap(~step, ncol = 3) +
    coord_fixed(ratio = 1) +
    labs(title = "Rayleigh-Bénard Convection Evolution", x = "X", y = "Y") +
    theme_minimal()

  print(p_animation)
  ggsave(
    "rayleigh_benard_evolution.png",
    p_animation,
    width = 12,
    height = 8,
    dpi = 300
  )
}

cat("Simulation complete!\n")
cat("- Static plot: rayleigh_benard_static.png\n")
cat("- Interactive plot: rayleigh_benard_interactive.html\n")
