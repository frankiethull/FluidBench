# Rayleigh-Bénard Convection Simulation
# 2D Incompressible Navier-Stokes + Energy Equation
# Boussinesq Approximation for Buoyancy-Driven Flow
#
# Governing Equations:
# ∇·u = 0 (continuity)
# ∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u + βgT ĵ (momentum + buoyancy)
# ∂T/∂t + (u·∇)T = α∇²T (energy)

# Load required libraries
library(ggplot2)
library(reshape2)

# Check for plotly
use_plotly <- requireNamespace("plotly", quietly = TRUE)
if (use_plotly) {
  library(plotly)
  cat("plotly available - will create interactive visualizations\n")
} else {
  cat("plotly not found - will use ggplot2 for all visualizations\n")
}

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================

cat("============================================\n")
cat("RAYLEIGH-BÉNARD CONVECTION SIMULATION\n")
cat("Boussinesq Approximation\n")
cat("============================================\n\n")

# Domain parameters (aspect ratio 2:1)
Lx <- 2.0 # Domain width
Ly <- 1.0 # Domain height
Nx <- 128 # Grid points in x
Ny <- 64 # Grid points in y
dx <- Lx / Nx
dy <- Ly / Ny

# Physical parameters
Pr <- 0.71 # Prandtl number (air)
Ra <- 1e5 # Rayleigh number (controls convection intensity)

# Derived parameters (non-dimensional formulation)
nu <- sqrt(Pr / Ra) # Kinematic viscosity
alpha <- 1 / sqrt(Pr * Ra) # Thermal diffusivity
g_beta <- 1.0 # Buoyancy coefficient (g*β*ΔT in non-dim form)

# Time parameters
dt <- 0.0002 # Time step
T_final <- 5.0 # Total simulation time
n_steps <- ceiling(T_final / dt)
save_interval <- 500

# Boundary temperatures
T_hot <- 1.0 # Bottom wall temperature
T_cold <- 0.0 # Top wall temperature

cat(sprintf("Domain: %.1f x %.1f (aspect ratio 2:1)\n", Lx, Ly))
cat(sprintf("Grid: %d x %d\n", Nx, Ny))
cat(sprintf("Rayleigh number: %.2e\n", Ra))
cat(sprintf("Prandtl number: %.2f\n", Pr))
cat(sprintf("Kinematic viscosity (ν): %.6f\n", nu))
cat(sprintf("Thermal diffusivity (α): %.6f\n", alpha))
cat(sprintf("Time step: %.5f\n", dt))
cat(sprintf("Total steps: %d\n\n", n_steps))

# ============================================================================
# GRID SETUP
# ============================================================================

# Staggered grid arrangement:
# - Pressure/Temperature at cell centers
# - u-velocity at vertical cell faces
# - v-velocity at horizontal cell faces

x_c <- seq(dx / 2, Lx - dx / 2, by = dx) # Cell centers (x)
y_c <- seq(dy / 2, Ly - dy / 2, by = dy) # Cell centers (y)
x_u <- seq(0, Lx, by = dx) # u-velocity faces
y_v <- seq(0, Ly, by = dy) # v-velocity faces

# ============================================================================
# INITIALIZE FIELDS
# ============================================================================

cat("Initializing fields...\n")

# Velocity components (on staggered grid)
u <- matrix(0, nrow = Nx + 1, ncol = Ny) # u at vertical faces
v <- matrix(0, nrow = Nx, ncol = Ny + 1) # v at horizontal faces

# Pressure and temperature at cell centers
p <- matrix(0, nrow = Nx, ncol = Ny)

# Temperature: linear profile with small perturbation
set.seed(42)
T_field <- matrix(0, nrow = Nx, ncol = Ny)
for (j in 1:Ny) {
  T_field[, j] <- T_hot - (T_hot - T_cold) * y_c[j] / Ly
}

# Add small random perturbation to break symmetry
perturbation <- 0.05 * matrix(rnorm(Nx * Ny), nrow = Nx, ncol = Ny)
# Apply perturbation more strongly in the middle
for (i in 1:Nx) {
  for (j in 1:Ny) {
    y_factor <- sin(pi * y_c[j] / Ly)
    T_field[i, j] <- T_field[i, j] + perturbation[i, j] * y_factor
  }
}

cat("Initial temperature field with perturbation created\n\n")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Interpolate u to cell centers
u_to_center <- function(u) {
  u_c <- matrix(0, nrow = Nx, ncol = Ny)
  for (i in 1:Nx) {
    u_c[i, ] <- 0.5 * (u[i, ] + u[i + 1, ])
  }
  return(u_c)
}

# Interpolate v to cell centers
v_to_center <- function(v) {
  v_c <- matrix(0, nrow = Nx, ncol = Ny)
  for (j in 1:Ny) {
    v_c[, j] <- 0.5 * (v[, j] + v[, j + 1])
  }
  return(v_c)
}

# Compute divergence of velocity field
compute_divergence <- function(u, v, dx, dy) {
  div <- matrix(0, nrow = Nx, ncol = Ny)
  for (i in 1:Nx) {
    for (j in 1:Ny) {
      div[i, j] <- (u[i + 1, j] - u[i, j]) / dx + (v[i, j + 1] - v[i, j]) / dy
    }
  }
  return(div)
}

# ============================================================================
# PRESSURE POISSON SOLVER (Jacobi iteration)
# ============================================================================

solve_pressure_poisson <- function(div, p, dx, dy, max_iter = 500, tol = 1e-5) {
  dx2 <- dx^2
  dy2 <- dy^2
  coef <- 2 * (1 / dx2 + 1 / dy2)

  for (iter in 1:max_iter) {
    p_old <- p

    # Interior points
    for (i in 2:(Nx - 1)) {
      for (j in 2:(Ny - 1)) {
        p[i, j] <- ((p_old[i + 1, j] + p_old[i - 1, j]) /
          dx2 +
          (p_old[i, j + 1] + p_old[i, j - 1]) / dy2 -
          div[i, j]) /
          coef
      }
    }

    # Periodic boundary in x
    for (j in 2:(Ny - 1)) {
      p[1, j] <- ((p_old[2, j] + p_old[Nx, j]) /
        dx2 +
        (p_old[1, j + 1] + p_old[1, j - 1]) / dy2 -
        div[1, j]) /
        coef
      p[Nx, j] <- ((p_old[1, j] + p_old[Nx - 1, j]) /
        dx2 +
        (p_old[Nx, j + 1] + p_old[Nx, j - 1]) / dy2 -
        div[Nx, j]) /
        coef
    }

    # Neumann BC at top and bottom (∂p/∂y = 0)
    p[, 1] <- p[, 2]
    p[, Ny] <- p[, Ny - 1]

    # Check convergence
    residual <- max(abs(p - p_old))
    if (residual < tol) break
  }

  return(p)
}

# ============================================================================
# ADVECTION TERMS (2nd order central + upwind blend)
# ============================================================================

compute_advection_u <- function(u, v, dx, dy) {
  adv <- matrix(0, nrow = Nx + 1, ncol = Ny)

  for (i in 2:Nx) {
    for (j in 1:Ny) {
      # u at this face
      u_here <- u[i, j]

      # ∂(uu)/∂x at u-face
      u_e <- 0.5 * (u[i, j] + u[min(i + 1, Nx + 1), j])
      u_w <- 0.5 * (u[i, j] + u[max(i - 1, 1), j])
      dudx <- (u_e^2 - u_w^2) / dx

      # ∂(uv)/∂y at u-face
      if (j > 1 && j < Ny) {
        v_n <- 0.5 * (v[i - 1, j + 1] + v[min(i, Nx), j + 1])
        v_s <- 0.5 * (v[i - 1, j] + v[min(i, Nx), j])
        u_n <- 0.5 * (u[i, j] + u[i, j + 1])
        u_s <- 0.5 * (u[i, j] + u[i, j - 1])
        dudy <- (v_n * u_n - v_s * u_s) / dy
      } else if (j == 1) {
        v_n <- 0.5 * (v[max(i - 1, 1), j + 1] + v[min(i, Nx), j + 1])
        u_n <- 0.5 * (u[i, j] + u[i, j + 1])
        dudy <- (v_n * u_n) / dy
      } else {
        v_s <- 0.5 * (v[max(i - 1, 1), j] + v[min(i, Nx), j])
        u_s <- 0.5 * (u[i, j] + u[i, j - 1])
        dudy <- (-v_s * u_s) / dy
      }

      adv[i, j] <- dudx + dudy
    }
  }

  return(adv)
}

compute_advection_v <- function(u, v, dx, dy) {
  adv <- matrix(0, nrow = Nx, ncol = Ny + 1)

  for (i in 1:Nx) {
    for (j in 2:Ny) {
      # ∂(uv)/∂x at v-face
      i_e <- ifelse(i < Nx, i + 1, 1)
      i_w <- ifelse(i > 1, i - 1, Nx)

      u_e <- 0.5 * (u[i + 1, j - 1] + u[i + 1, min(j, Ny)])
      u_w <- 0.5 * (u[i, j - 1] + u[i, min(j, Ny)])
      v_e <- 0.5 * (v[i, j] + v[i_e, j])
      v_w <- 0.5 * (v[i, j] + v[i_w, j])
      dvdx <- (u_e * v_e - u_w * v_w) / dx

      # ∂(vv)/∂y at v-face
      v_n <- 0.5 * (v[i, j] + v[i, min(j + 1, Ny + 1)])
      v_s <- 0.5 * (v[i, j] + v[i, max(j - 1, 1)])
      dvdy <- (v_n^2 - v_s^2) / dy

      adv[i, j] <- dvdx + dvdy
    }
  }

  return(adv)
}

# ============================================================================
# DIFFUSION TERMS (Laplacian)
# ============================================================================

compute_laplacian_u <- function(u, dx, dy) {
  lap <- matrix(0, nrow = Nx + 1, ncol = Ny)

  for (i in 2:Nx) {
    for (j in 1:Ny) {
      # Second derivative in x (periodic)
      u_xx <- (u[i + 1, j] - 2 * u[i, j] + u[i - 1, j]) / dx^2

      # Second derivative in y (with wall BC)
      if (j == 1) {
        u_yy <- (u[i, j + 1] - 2 * u[i, j] + 0) / dy^2 # u=0 at wall
      } else if (j == Ny) {
        u_yy <- (0 - 2 * u[i, j] + u[i, j - 1]) / dy^2 # u=0 at wall
      } else {
        u_yy <- (u[i, j + 1] - 2 * u[i, j] + u[i, j - 1]) / dy^2
      }

      lap[i, j] <- u_xx + u_yy
    }
  }

  return(lap)
}

compute_laplacian_v <- function(v, dx, dy) {
  lap <- matrix(0, nrow = Nx, ncol = Ny + 1)

  for (i in 1:Nx) {
    for (j in 2:Ny) {
      # Second derivative in x (periodic)
      i_e <- ifelse(i < Nx, i + 1, 1)
      i_w <- ifelse(i > 1, i - 1, Nx)
      v_xx <- (v[i_e, j] - 2 * v[i, j] + v[i_w, j]) / dx^2

      # Second derivative in y
      v_yy <- (v[i, j + 1] - 2 * v[i, j] + v[i, j - 1]) / dy^2

      lap[i, j] <- v_xx + v_yy
    }
  }

  return(lap)
}

# ============================================================================
# TEMPERATURE ADVECTION AND DIFFUSION
# ============================================================================

compute_T_advection <- function(T_field, u, v, dx, dy) {
  adv <- matrix(0, nrow = Nx, ncol = Ny)

  u_c <- u_to_center(u)
  v_c <- v_to_center(v)

  for (i in 1:Nx) {
    for (j in 1:Ny) {
      # Upwind scheme for stability
      i_e <- ifelse(i < Nx, i + 1, 1)
      i_w <- ifelse(i > 1, i - 1, Nx)

      if (u_c[i, j] > 0) {
        dTdx <- (T_field[i, j] - T_field[i_w, j]) / dx
      } else {
        dTdx <- (T_field[i_e, j] - T_field[i, j]) / dx
      }

      if (j == 1) {
        dTdy <- (T_field[i, j + 1] - T_field[i, j]) / dy
      } else if (j == Ny) {
        dTdy <- (T_field[i, j] - T_field[i, j - 1]) / dy
      } else {
        if (v_c[i, j] > 0) {
          dTdy <- (T_field[i, j] - T_field[i, j - 1]) / dy
        } else {
          dTdy <- (T_field[i, j + 1] - T_field[i, j]) / dy
        }
      }

      adv[i, j] <- u_c[i, j] * dTdx + v_c[i, j] * dTdy
    }
  }

  return(adv)
}

compute_T_laplacian <- function(T_field, dx, dy) {
  lap <- matrix(0, nrow = Nx, ncol = Ny)

  for (i in 1:Nx) {
    for (j in 1:Ny) {
      i_e <- ifelse(i < Nx, i + 1, 1)
      i_w <- ifelse(i > 1, i - 1, Nx)

      T_xx <- (T_field[i_e, j] - 2 * T_field[i, j] + T_field[i_w, j]) / dx^2

      if (j == 1) {
        T_yy <- (T_field[i, j + 1] - 2 * T_field[i, j] + T_hot) / dy^2
      } else if (j == Ny) {
        T_yy <- (T_cold - 2 * T_field[i, j] + T_field[i, j - 1]) / dy^2
      } else {
        T_yy <- (T_field[i, j + 1] - 2 * T_field[i, j] + T_field[i, j - 1]) /
          dy^2
      }

      lap[i, j] <- T_xx + T_yy
    }
  }

  return(lap)
}

# ============================================================================
# BUOYANCY FORCE
# ============================================================================

compute_buoyancy <- function(T_field) {
  # Buoyancy acts on v (vertical velocity)
  # Interpolate T to v-faces
  buoy <- matrix(0, nrow = Nx, ncol = Ny + 1)

  for (i in 1:Nx) {
    for (j in 2:Ny) {
      T_at_face <- 0.5 * (T_field[i, j - 1] + T_field[i, j])
      # Buoyancy: hot fluid rises (positive buoyancy when T > T_ref)
      T_ref <- 0.5 * (T_hot + T_cold)
      buoy[i, j] <- g_beta * (T_at_face - T_ref)
    }
  }

  return(buoy)
}

# ============================================================================
# BOUNDARY CONDITIONS
# ============================================================================

apply_velocity_bc <- function(u, v) {
  # No-slip at top and bottom walls
  # u = 0 at walls (already handled in interior)

  # v = 0 at top and bottom walls
  v[, 1] <- 0 # Bottom
  v[, Ny + 1] <- 0 # Top

  # Periodic in x
  u[1, ] <- u[Nx + 1, ]
  u[Nx + 1, ] <- u[2, ] # Wrap around

  return(list(u = u, v = v))
}

apply_temperature_bc <- function(T_field) {
  # Fixed temperature at walls (Dirichlet BC applied in Laplacian)
  # Periodic in x (already handled in advection/diffusion)
  return(T_field)
}

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

cat("Starting simulation...\n")
cat("Progress: ")

start_time <- Sys.time()

# Storage for animation
frames <- list()
frame_count <- 0

for (step in 1:n_steps) {
  # ===== Step 1: Compute advection terms =====
  adv_u <- compute_advection_u(u, v, dx, dy)
  adv_v <- compute_advection_v(u, v, dx, dy)

  # ===== Step 2: Compute diffusion terms =====
  lap_u <- compute_laplacian_u(u, dx, dy)
  lap_v <- compute_laplacian_v(v, dx, dy)

  # ===== Step 3: Compute buoyancy force =====
  buoy <- compute_buoyancy(T_field)

  # ===== Step 4: Predictor step (without pressure) =====
  u_star <- u + dt * (-adv_u + nu * lap_u)
  v_star <- v + dt * (-adv_v + nu * lap_v + buoy)

  # Apply velocity BC to intermediate velocity
  bc_result <- apply_velocity_bc(u_star, v_star)
  u_star <- bc_result$u
  v_star <- bc_result$v

  # ===== Step 5: Solve pressure Poisson equation =====
  div <- compute_divergence(u_star, v_star, dx, dy)
  p <- solve_pressure_poisson(div / dt, p, dx, dy)

  # ===== Step 6: Correct velocity (projection) =====
  # u = u* - dt * ∂p/∂x
  for (i in 2:Nx) {
    for (j in 1:Ny) {
      u[i, j] <- u_star[i, j] - dt * (p[i, j] - p[i - 1, j]) / dx
    }
  }

  # v = v* - dt * ∂p/∂y
  for (i in 1:Nx) {
    for (j in 2:Ny) {
      v[i, j] <- v_star[i, j] - dt * (p[i, j] - p[i, j - 1]) / dy
    }
  }

  # Apply velocity BC
  bc_result <- apply_velocity_bc(u, v)
  u <- bc_result$u
  v <- bc_result$v

  # ===== Step 7: Update temperature =====
  T_adv <- compute_T_advection(T_field, u, v, dx, dy)
  T_lap <- compute_T_laplacian(T_field, dx, dy)

  T_field <- T_field + dt * (-T_adv + alpha * T_lap)
  T_field <- apply_temperature_bc(T_field)

  # Clamp temperature
  T_field <- pmax(T_cold, pmin(T_hot, T_field))

  # ===== Save frame =====
  if (step %% save_interval == 0 || step == n_steps) {
    frame_count <- frame_count + 1
    frames[[frame_count]] <- list(
      step = step,
      time = step * dt,
      T = T_field,
      u = u,
      v = v,
      u_c = u_to_center(u),
      v_c = v_to_center(v)
    )
  }

  # Progress
  if (step %% (n_steps %/% 20) == 0) {
    progress <- round(100 * step / n_steps)
    cat(sprintf("%d%% ", progress))
  }
}

end_time <- Sys.time()
cat("\n\nSimulation completed!\n")
cat(sprintf(
  "Time elapsed: %.1f seconds\n\n",
  difftime(end_time, start_time, units = "secs")
))

# ============================================================================
# VISUALIZATION: COMBINED TEMPERATURE + VELOCITY FIELD
# ============================================================================

cat("Creating visualizations...\n")

# Get final frame
final <- frames[[length(frames)]]

# Prepare temperature data
T_df <- melt(final$T)
colnames(T_df) <- c("xi", "yi", "Temperature")
T_df$x <- x_c[T_df$xi]
T_df$y <- y_c[T_df$yi]

# Prepare velocity data (downsampled for arrows)
skip_x <- 4
skip_y <- 2
x_idx <- seq(1, Nx, by = skip_x)
y_idx <- seq(1, Ny, by = skip_y)

vel_df <- expand.grid(xi = x_idx, yi = y_idx)
vel_df$x <- x_c[vel_df$xi]
vel_df$y <- y_c[vel_df$yi]
vel_df$u <- final$u_c[cbind(vel_df$xi, vel_df$yi)]
vel_df$v <- final$v_c[cbind(vel_df$xi, vel_df$yi)]
vel_df$vel_mag <- sqrt(vel_df$u^2 + vel_df$v^2)

# Scale arrows for visibility
arrow_scale <- 0.15 / max(vel_df$vel_mag + 1e-10)

# ===== Main combined plot =====
p_combined <- ggplot() +
  # Temperature heatmap
  geom_raster(
    data = T_df,
    aes(x = x, y = y, fill = Temperature),
    interpolate = TRUE
  ) +
  scale_fill_gradientn(
    colors = c(
      "#00008B",
      "#0000FF",
      "#4169E1",
      "#87CEEB",
      "#FFFFE0",
      "#FFD700",
      "#FF8C00",
      "#FF4500",
      "#8B0000"
    ),
    values = scales::rescale(c(0, 0.2, 0.35, 0.45, 0.5, 0.55, 0.65, 0.8, 1)),
    name = "Temperature",
    limits = c(0, 1)
  ) +
  # Velocity arrows
  geom_segment(
    data = vel_df,
    aes(x = x, y = y, xend = x + u * arrow_scale, yend = y + v * arrow_scale),
    arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
    color = "black",
    alpha = 0.7,
    linewidth = 0.4
  ) +
  labs(
    title = "Rayleigh-Bénard Convection",
    subtitle = sprintf(
      "Ra = %.0e, Pr = %.2f | t = %.2f | Convection Rolls",
      Ra,
      Pr,
      final$time
    ),
    x = "x",
    y = "y",
    caption = "Hot fluid rises (red), cold fluid sinks (blue) | Arrows show velocity field"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.caption = element_text(size = 10, color = "gray40"),
    legend.position = "right",
    panel.grid = element_blank(),
    aspect.ratio = Ly / Lx
  ) +
  coord_fixed(xlim = c(0, Lx), ylim = c(0, Ly), expand = FALSE)

print(p_combined)
ggsave(
  "rayleigh_benard_combined.png",
  p_combined,
  width = 14,
  height = 7,
  dpi = 200
)
cat("Combined plot saved as 'rayleigh_benard_combined.png'\n")

# ============================================================================
# PLOTLY INTERACTIVE VISUALIZATION
# ============================================================================

if (use_plotly) {
  cat("Creating interactive plotly visualization...\n")

  # Create plotly heatmap with velocity overlay
  fig <- plot_ly() %>%
    # Temperature heatmap
    add_heatmap(
      x = x_c,
      y = y_c,
      z = t(final$T),
      colorscale = list(
        c(0, "#00008B"),
        c(0.2, "#0000FF"),
        c(0.35, "#4169E1"),
        c(0.45, "#87CEEB"),
        c(0.5, "#FFFFE0"),
        c(0.55, "#FFD700"),
        c(0.65, "#FF8C00"),
        c(0.8, "#FF4500"),
        c(1, "#8B0000")
      ),
      colorbar = list(title = "Temperature"),
      zmin = 0,
      zmax = 1,
      hovertemplate = "x: %{x:.2f}<br>y: %{y:.2f}<br>T: %{z:.3f}<extra></extra>"
    )

  # Add velocity arrows as annotations
  arrows <- list()
  for (i in 1:nrow(vel_df)) {
    if (vel_df$vel_mag[i] > 0.001) {
      arrows[[length(arrows) + 1]] <- list(
        x = vel_df$x[i] + vel_df$u[i] * arrow_scale * 0.8,
        y = vel_df$y[i] + vel_df$v[i] * arrow_scale * 0.8,
        ax = vel_df$x[i],
        ay = vel_df$y[i],
        xref = "x",
        yref = "y",
        axref = "x",
        ayref = "y",
        showarrow = TRUE,
        arrowhead = 2,
        arrowsize = 1,
        arrowwidth = 1.5,
        arrowcolor = "rgba(0,0,0,0.6)"
      )
    }
  }

  fig <- fig %>%
    layout(
      title = list(
        text = sprintf(
          "<b>Rayleigh-Bénard Convection</b><br>Ra = %.0e, Pr = %.2f",
          Ra,
          Pr
        ),
        font = list(size = 18)
      ),
      xaxis = list(title = "x", scaleanchor = "y"),
      yaxis = list(title = "y"),
      annotations = arrows
    )

  # Save as HTML
  htmlwidgets::saveWidget(fig, "rayleigh_benard_interactive.html")
  cat("Interactive plot saved as 'rayleigh_benard_interactive.html'\n")
}

# ============================================================================
# STREAMLINES VISUALIZATION
# ============================================================================

cat("Creating streamline visualization...\n")

# Compute stream function (ψ) where u = ∂ψ/∂y, v = -∂ψ/∂x
psi <- matrix(0, nrow = Nx, ncol = Ny)
for (j in 2:Ny) {
  for (i in 1:Nx) {
    psi[i, j] <- psi[i, j - 1] + final$u_c[i, j] * dy
  }
}

psi_df <- melt(psi)
colnames(psi_df) <- c("xi", "yi", "psi")
psi_df$x <- x_c[psi_df$xi]
psi_df$y <- y_c[psi_df$yi]

p_stream <- ggplot() +
  geom_raster(
    data = T_df,
    aes(x = x, y = y, fill = Temperature),
    interpolate = TRUE,
    alpha = 0.7
  ) +
  scale_fill_gradientn(
    colors = c(
      "#00008B",
      "#4169E1",
      "#87CEEB",
      "#FFFFE0",
      "#FF8C00",
      "#8B0000"
    ),
    name = "Temperature"
  ) +
  geom_contour(
    data = psi_df,
    aes(x = x, y = y, z = psi),
    color = "black",
    linewidth = 0.5,
    bins = 20
  ) +
  labs(
    title = "Streamlines and Temperature Field",
    subtitle = "Closed streamlines indicate convection cells",
    x = "x",
    y = "y"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    aspect.ratio = Ly / Lx
  ) +
  coord_fixed(expand = FALSE)

print(p_stream)
ggsave(
  "rayleigh_benard_streamlines.png",
  p_stream,
  width = 14,
  height = 7,
  dpi = 150
)

# ============================================================================
# VORTICITY FIELD
# ============================================================================

cat("Creating vorticity visualization...\n")

# Compute vorticity ω = ∂v/∂x - ∂u/∂y
vorticity <- matrix(0, nrow = Nx, ncol = Ny)
for (i in 1:Nx) {
  for (j in 1:Ny) {
    i_e <- ifelse(i < Nx, i + 1, 1)
    i_w <- ifelse(i > 1, i - 1, Nx)

    dvdx <- (final$v_c[i_e, j] - final$v_c[i_w, j]) / (2 * dx)

    if (j == 1) {
      dudy <- (final$u_c[i, j + 1] - final$u_c[i, j]) / dy
    } else if (j == Ny) {
      dudy <- (final$u_c[i, j] - final$u_c[i, j - 1]) / dy
    } else {
      dudy <- (final$u_c[i, j + 1] - final$u_c[i, j - 1]) / (2 * dy)
    }

    vorticity[i, j] <- dvdx - dudy
  }
}

vort_df <- melt(vorticity)
colnames(vort_df) <- c("xi", "yi", "vorticity")
vort_df$x <- x_c[vort_df$xi]
vort_df$y <- y_c[vort_df$yi]

vort_max <- quantile(abs(vorticity), 0.98)

p_vort <- ggplot(vort_df, aes(x = x, y = y, fill = vorticity)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(
    low = "#0000FF",
    mid = "white",
    high = "#FF0000",
    midpoint = 0,
    limits = c(-vort_max, vort_max),
    oob = scales::squish,
    name = "Vorticity\n(ω)"
  ) +
  labs(
    title = "Vorticity Field",
    subtitle = "Red: clockwise, Blue: counter-clockwise rotation",
    x = "x",
    y = "y"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    aspect.ratio = Ly / Lx
  ) +
  coord_fixed(expand = FALSE)

print(p_vort)
ggsave(
  "rayleigh_benard_vorticity.png",
  p_vort,
  width = 14,
  height = 7,
  dpi = 150
)

# ============================================================================
# TIME EVOLUTION SEQUENCE
# ============================================================================

cat("Creating evolution sequence...\n")

# Select frames to display
n_display <- min(6, length(frames))
display_indices <- round(seq(1, length(frames), length.out = n_display))

evolution_df <- data.frame()
for (idx in display_indices) {
  frame <- frames[[idx]]

  temp_T <- melt(frame$T)
  colnames(temp_T) <- c("xi", "yi", "Temperature")
  temp_T$x <- x_c[temp_T$xi]
  temp_T$y <- y_c[temp_T$yi]
  temp_T$time_label <- sprintf("t = %.2f", frame$time)

  evolution_df <- rbind(evolution_df, temp_T)
}

evolution_df$time_label <- factor(
  evolution_df$time_label,
  levels = unique(evolution_df$time_label)
)

p_evolution <- ggplot(evolution_df, aes(x = x, y = y, fill = Temperature)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradientn(
    colors = c(
      "#00008B",
      "#4169E1",
      "#87CEEB",
      "#FFFFE0",
      "#FF8C00",
      "#8B0000"
    ),
    name = "T",
    limits = c(0, 1)
  ) +
  facet_wrap(~time_label, nrow = 2) +
  labs(
    title = "Evolution of Rayleigh-Bénard Convection",
    subtitle = "From conductive state to established convection rolls"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    strip.text = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  coord_fixed(expand = FALSE)

print(p_evolution)
ggsave(
  "rayleigh_benard_evolution.png",
  p_evolution,
  width = 14,
  height = 8,
  dpi = 150
)

# ============================================================================
# VERTICAL TEMPERATURE PROFILE
# ============================================================================

cat("Creating temperature profile plot...\n")

# Average temperature profile
T_profile <- colMeans(final$T)

# Conductive profile for comparison
T_conductive <- T_hot - (T_hot - T_cold) * y_c / Ly

profile_df <- data.frame(
  y = rep(y_c, 2),
  T = c(T_profile, T_conductive),
  type = rep(c("Convective (simulated)", "Conductive (linear)"), each = Ny)
)

p_profile <- ggplot(
  profile_df,
  aes(x = T, y = y, color = type, linetype = type)
) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c(
      "Convective (simulated)" = "#FF4500",
      "Conductive (linear)" = "#4169E1"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Convective (simulated)" = "solid",
      "Conductive (linear)" = "dashed"
    )
  ) +
  labs(
    title = "Vertical Temperature Profile",
    subtitle = "Convection enhances heat transfer compared to pure conduction",
    x = "Temperature",
    y = "Height (y)",
    color = "Profile Type",
    linetype = "Profile Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "top"
  )

print(p_profile)
ggsave(
  "rayleigh_benard_profile.png",
  p_profile,
  width = 8,
  height = 8,
  dpi = 150
)

# ============================================================================
# NUSSELT NUMBER CALCULATION
# ============================================================================

# Nusselt number: ratio of total heat transfer to conductive heat transfer
# Nu = (q_total) / (q_conductive) = 1 + convective contribution

# At the bottom wall: Nu = -Ly * (dT/dy)|_wall / (T_hot - T_cold)
dTdy_bottom <- (final$T[, 1] - T_hot) / (dy / 2)
Nu_local <- -Ly * dTdy_bottom / (T_hot - T_cold)
Nu_avg <- mean(Nu_local)

cat(sprintf("\nNusselt number (average): Nu = %.3f\n", Nu_avg))
cat("(Nu > 1 indicates convection enhances heat transfer)\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n============================================\n")
cat("SIMULATION SUMMARY\n")
cat("============================================\n")
cat("Physics: Rayleigh-Bénard Convection\n")
cat("Model: 2D Incompressible Navier-Stokes + Energy Equation\n")
cat("Approximation: Boussinesq (density varies with temperature)\n")
cat(sprintf("\nParameters:\n"))
cat(sprintf("  Rayleigh number (Ra): %.2e\n", Ra))
cat(sprintf("  Prandtl number (Pr): %.2f\n", Pr))
cat(sprintf("  Domain: %.1f x %.1f (aspect ratio 2:1)\n", Lx, Ly))
cat(sprintf("  Grid: %d x %d\n", Nx, Ny))
cat(sprintf("\nResults:\n"))
cat(sprintf("  Nusselt number: %.3f\n", Nu_avg))
cat(sprintf("  Max velocity: %.4f\n", max(sqrt(final$u_c^2 + final$v_c^2))))
cat("\nOutput Files:\n")
cat("  1. rayleigh_benard_combined.png - Temperature + velocity field\n")
if (use_plotly) {
  cat(
    "  2. rayleigh_benard_interactive.html - Interactive plotly visualization\n"
  )
}
cat("  3. rayleigh_benard_streamlines.png - Streamlines overlay\n")
cat("  4. rayleigh_benard_vorticity.png - Vorticity field\n")
cat("  5. rayleigh_benard_evolution.png - Time evolution sequence\n")
cat("  6. rayleigh_benard_profile.png - Temperature profile\n")
cat("\nPhysical Interpretation:\n")
cat("  - Hot fluid at bottom becomes buoyant and rises\n")
cat("  - Cold fluid at top sinks\n")
cat("  - Results in organized convection rolls (cells)\n")
cat("  - Number of rolls depends on Ra and aspect ratio\n")
cat("============================================\n")
