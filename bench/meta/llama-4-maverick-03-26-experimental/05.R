# Libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

# Parameters
Nx <- 100 # Grid points (x)
Ny <- 50 # Grid points (y)
Lx <- 2 # Domain length (x)
Ly <- 1 # Domain height (y)
dx <- Lx / (Nx - 1)
dy <- Ly / (Ny - 1)
dt <- 0.001 # Time step (stability-critical)
nu <- 0.01 # Kinematic viscosity
kappa <- 0.01 # Thermal diffusivity
alpha <- 1 # Thermal expansion
g <- 10 # Gravity
Ra <- alpha * g * (1 - 0) * Ly^3 / (nu * kappa) # Rayleigh number
cat(sprintf("Rayleigh Number = %.2f\n", Ra))

# Initialize staggered grids
u <- matrix(0, nrow = Nx - 1, ncol = Ny) # Horizontal velocity (cell edges)
v <- matrix(0, nrow = Nx, ncol = Ny - 1) # Vertical velocity (cell edges)
T <- matrix(0.5, nrow = Nx, ncol = Ny) # Temperature (cell centers)
p <- matrix(0, nrow = Nx, ncol = Ny) # Pressure (cell centers)

# Boundary conditions
T[, 1] <- 1 # Bottom hot
T[, Ny] <- 0 # Top cold
T[1, ] <- T[2, ] # Left Neumann
T[Nx, ] <- T[Nx - 1, ] # Right Neumann

# Helper functions
laplacian <- function(F, dx, dy) {
  F_ext <- cbind(F[, 1], F, F[, Ny])
  F_ext <- rbind(F_ext[1, ], F_ext, F_ext[Nx, ])
  lap <- ((F_ext[1:Nx, 2:(Ny + 1)] + F_ext[3:(Nx + 2), 2:(Ny + 1)]) /
    dy^2 +
    (F_ext[2:(Nx + 1), 1:Ny] + F_ext[2:(Nx + 1), 3:(Ny + 2)]) / dx^2 -
    2 * (1 / dx^2 + 1 / dy^2) * F_ext[2:(Nx + 1), 2:(Ny + 1)])
  lap
}

grad_x <- function(F, dx) {
  (cbind(F[, 2:Ny], F[, Ny]) - cbind(F[, 1], F[, 1:(Ny - 1)])) / dx
}
grad_y <- function(F, dy) {
  (rbind(F[2:Nx, ], F[Nx, ]) - rbind(F[1, ], F[1:(Nx - 1), ])) / dy
}

# Time-stepping loop
steps <- 5000
for (step in 1:steps) {
  # 1. Temperature advection-diffusion (FTCS scheme)
  T_new <- T +
    dt *
      (-(u %*% t(grad_x(T, dx))) -
        (t(v) %*% grad_y(T, dy)) +
        kappa * laplacian(T, dx, dy))
  T_new[, 1] <- 1
  T_new[, Ny] <- 0
  T_new[1, ] <- T_new[2, ]
  T_new[Nx, ] <- T_new[Nx - 1, ]

  # 2. Momentum (predictor step without pressure)
  u_star <- u +
    dt *
      (-(0.5 *
        (u + cbind(u[, 2:Ny], u[, Ny])) *
        (u - cbind(u[, 1], u[, 1:(Ny - 1)])) /
        dx) -
        (0.5 *
          (v[1:(Nx - 1), ] + v[2:Nx, ]) *
          (u - rbind(u[1, ], u[1:(Nx - 2), ])) /
          dy) +
        nu * (laplacian(cbind(u, u[, Ny]), dx, dy)[, 1:Ny]))
  v_star <- v +
    dt *
      (
        -(0.5 *
          (u[, 1:(Ny - 1)] + u[, 2:Ny]) *
          (v - cbind(v[, 1], v[, 1:(Ny - 2)])) /
          dx) -
          (0.5 *
            (v + rbind(v[2:Nx, ], v[Nx, ])) *
            (v - rbind(v[1, ], v[1:(Nx - 1), ])) /
            dy) +
          nu * (laplacian(rbind(v, v[Nx, ]), dx, dy)[1:Nx, ]) +
          alpha * g * 0.5 * (T[1:Nx, 1:(Ny - 1)] + T[1:Nx, 2:Ny]) # Buoyancy term
      )

  # 3. Pressure Poisson equation (∇²p = -∇·u_star/dt)
  div_u_star <- ((cbind(u_star[, 2:Ny], u_star[, Ny]) -
    cbind(u_star[, 1], u_star[, 1:(Ny - 1)])) /
    dx +
    (rbind(v_star[2:Nx, ], v_star[Nx, ]) -
      rbind(v_star[1, ], v_star[1:(Nx - 1), ])) /
      dy)
  p_new <- p
  for (iter in 1:20) {
    # Jacobi iterations for ∇²p = div_u_star/dt
    p_old <- p_new
    p_new <- 0.25 *
      (cbind(p_old[, 2:Ny], p_old[, Ny]) +
        cbind(p_old[, 1], p_old[, 1:(Ny - 1)]) +
        rbind(p_old[2:Nx, ], p_old[Nx, ]) +
        rbind(p_old[1, ], p_old[1:(Nx - 1), ]) -
        dx^2 * div_u_star / dt)
    p_new[, 1] <- p_new[, 2] # Neumann
    p_new[, Ny] <- p_new[, Ny - 1]
    p_new[1, ] <- p_new[2, ]
    p_new[Nx, ] <- p_new[Nx - 1, ]
  }

  # 4. Correct velocities (u^{n+1} = u_star - dt * ∇p)
  u_new <- u_star - dt * grad_x(p_new, dx)
  v_new <- v_star - dt * grad_y(p_new, dy)

  # Update fields
  T <- T_new
  u <- u_new
  v <- v_new
  p <- p_new

  # Progress bar
  if (step %% 500 == 0) cat(sprintf("Step %d/%d\n", step, steps))
}

# Visualization (ggplot2)
df_T <- expand.grid(
  x = seq(0, Lx, length.out = Nx),
  y = seq(0, Ly, length.out = Ny)
) %>%
  mutate(T = as.vector(T))
df_u <- expand.grid(
  x = seq(0, Lx, length.out = Nx - 1) + dx / 2,
  y = seq(0, Ly, length.out = Ny)
) %>%
  mutate(u = as.vector(u))
df_v <- expand.grid(
  x = seq(0, Lx, length.out = Nx),
  y = seq(0, Ly, length.out = Ny - 1) + dy / 2
) %>%
  mutate(v = as.vector(v))

ggplot() +
  geom_raster(data = df_T, aes(x, y, fill = T)) +
  scale_fill_gradient(low = "blue", high = "red") +
  geom_segment(
    data = df_u,
    aes(x = x, y = y, xend = x + u / 5, yend = y),
    arrow = arrow(length = unit(0.02, "npc")),
    color = "black"
  ) +
  geom_segment(
    data = df_v,
    aes(x = x, y = y, xend = x, yend = y + v / 5),
    arrow = arrow(length = unit(0.02, "npc")),
    color = "black"
  ) +
  labs(
    title = sprintf("Rayleigh-Bénard Convection (Ra = %.2f)", Ra),
    subtitle = "Temperature (color) + Velocity (arrows)",
    x = "X",
    y = "Y",
    fill = "T"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Optional: Interactive plotly version
# library(plotly)
# plot_ly(df_T, x = ~x, y = ~y, z = ~T, type = 'mesh3d') %>%
#   add_trace(data = df_u, x = ~x, y = ~y, z = ~rep(0, nrow(df_u)),
#             u = ~u, v = ~rep(0, nrow(df_u)), w = ~rep(0, nrow(df_u)),
#             type = 'cone', sizemode = "absolute", sizeref = 0.5) %>%
#   add_trace(data = df_v, x = ~x, y = ~y, z = ~rep(0, nrow(df_v)),
#             u = ~rep(0, nrow(df_v)), v = ~v, w = ~rep(0, nrow(df_v)),
#             type = 'cone', sizemode = "absolute", sizeref = 0.5)
