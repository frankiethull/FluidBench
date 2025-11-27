# --------------------------------------------------------------
# 1. Load Required Libraries
# --------------------------------------------------------------
# install.packages(c("ggplot2", "pracma", "dplyr", "viridis"))
library(ggplot2)
library(pracma) # For the fast Poisson solver
library(dplyr)
library(viridis)

# --------------------------------------------------------------
# 2. Define Simulation Parameters
# --------------------------------------------------------------

# Grid parameters
Nx <- 80 # Grid points in x-direction
Ny <- 40 # Grid points in y-direction (aspect ratio 2:1)
Lx <- 2.0
Ly <- 1.0
dx <- Lx / Nx
dy <- Ly / Ny

# Physical parameters (non-dimensional)
Ra <- 1e4 # Rayleigh number (controls buoyancy strength)
Pr <- 0.71 # Prandtl number (ratio of momentum to thermal diffusivity)

# Time parameters
# Ensure numerical stability with a CFL condition
dt <- 0.2 * min(dx^2, dy^2) / (4 * max(1, Pr))
n_steps <- 5000
plot_every <- 250

# --------------------------------------------------------------
# 3. Initialize Fields and Grids
# --------------------------------------------------------------

# Create coordinate grids
x <- seq(dx / 2, Lx - dx / 2, length.out = Nx)
y <- seq(dy / 2, Ly - dy / 2, length.out = Ny)

# Initialize fields
T_field <- matrix(0, nrow = Nx, ncol = Ny) # Temperature
omega <- matrix(0, nrow = Nx, ncol = Ny) # Vorticity
psi <- matrix(0, nrow = Nx, ncol = Ny) # Stream function

# Set initial temperature condition
# Linear profile from hot (1) at bottom to cold (0) at top
for (j in 1:Ny) {
  T_field[, j] <- 1 - y[j]
}
# Add a small random perturbation to break symmetry and initiate convection
T_field <- T_field + matrix(runif(Nx * Ny, -0.01, 0.01), nrow = Nx, ncol = Ny)

# --------------------------------------------------------------
# 4. The Main Simulation Loop
# --------------------------------------------------------------

cat("Starting Rayleigh-Bénard simulation...\n")
pb <- txtProgressBar(min = 0, max = n_steps, style = 3)

for (step in 1:n_steps) {
  # --- A. Solve Poisson Equation for Stream Function ---
  # ∇²ψ = -ω
  # The pracma::poisson solver handles Dirichlet BCs (psi=0 on boundary)
  psi <- poisson(-omega, Lx, Nx, Ly, Ny)

  # --- B. Calculate Velocities from Stream Function ---
  # u = ∂ψ/∂y, v = -∂ψ/∂x
  # Use central differences, velocities are at cell centers
  u <- matrix(0, nrow = Nx, ncol = Ny)
  v <- matrix(0, nrow = Nx, ncol = Ny)

  u[, 2:(Ny - 1)] <- (psi[, 3:Ny] - psi[, 1:(Ny - 2)]) / (2 * dy)
  v[2:(Nx - 1), ] <- -(psi[3:Nx, ] - psi[1:(Nx - 2), ]) / (2 * dx)

  # --- C. Update Temperature (Advection-Diffusion) ---
  # ∂T/∂t + u∂T/∂x + v∂T/∂y = ∇²T
  dTdx <- matrix(0, nrow = Nx, ncol = Ny)
  dTdy <- matrix(0, nrow = Nx, ncol = Ny)
  lap_T <- matrix(0, nrow = Nx, ncol = Ny)

  dTdx[2:(Nx - 1), ] <- (T_field[3:Nx, ] - T_field[1:(Nx - 2), ]) / (2 * dx)
  dTdy[, 2:(Ny - 1)] <- (T_field[, 3:Ny] - T_field[, 1:(Ny - 2)]) / (2 * dy)
  lap_T[2:(Nx - 1), 2:(Ny - 1)] <- (T_field[3:Nx, 2:(Ny - 1)] -
    2 * T_field[2:(Nx - 1), 2:(Ny - 1)] +
    T_field[1:(Nx - 2), 2:(Ny - 1)]) /
    dx^2 +
    (T_field[2:(Nx - 1), 3:Ny] -
      2 * T_field[2:(Nx - 1), 2:(Ny - 1)] +
      T_field[2:(Nx - 1), 1:(Ny - 2)]) /
      dy^2

  T_field <- T_field + dt * (-u * dTdx - v * dTdy + lap_T)

  # --- D. Update Vorticity (Advection-Diffusion + Buoyancy) ---
  # ∂ω/∂t + u∂ω/∂x + v∂ω/∂y = Pr∇²ω + Ra*Pr*∂T/∂x
  domegadx <- matrix(0, nrow = Nx, ncol = Ny)
  domegady <- matrix(0, nrow = Nx, ncol = Ny)
  lap_omega <- matrix(0, nrow = Nx, ncol = Ny)

  domegadx[2:(Nx - 1), ] <- (omega[3:Nx, ] - omega[1:(Nx - 2), ]) / (2 * dx)
  domegady[, 2:(Ny - 1)] <- (omega[, 3:Ny] - omega[, 1:(Ny - 2)]) / (2 * dy)
  lap_omega[2:(Nx - 1), 2:(Ny - 1)] <- (omega[3:Nx, 2:(Ny - 1)] -
    2 * omega[2:(Nx - 1), 2:(Ny - 1)] +
    omega[1:(Nx - 2), 2:(Ny - 1)]) /
    dx^2 +
    (omega[2:(Nx - 1), 3:Ny] -
      2 * omega[2:(Nx - 1), 2:(Ny - 1)] +
      omega[2:(Nx - 1), 1:(Ny - 2)]) /
      dy^2

  omega <- omega +
    dt * (-u * domegadx - v * domegady + Pr * lap_omega + Ra * Pr * dTdx)

  # --- E. Apply Boundary Conditions ---
  # Temperature: Hot bottom, Cold top, Adiabatic sides
  T_field[, 1] <- 1.0 # Bottom wall
  T_field[, Ny] <- 0.0 # Top wall
  T_field[1, ] <- T_field[2, ] # Left wall (adiabatic: dT/dx=0)
  T_field[Nx, ] <- T_field[Nx - 1, ] # Right wall (adiabatic: dT/dx=0)

  # Stream function: psi = 0 on all boundaries (handled by poisson solver)

  # Vorticity: Use Thom's formula for no-slip walls
  # omega_wall = -2 * psi_first_interior / (delta^2)
  omega[1, ] <- -2 * psi[2, ] / dx^2 # Left wall
  omega[Nx, ] <- -2 * psi[Nx - 1, ] / dx^2 # Right wall
  omega[, 1] <- -2 * psi[, 2] / dy^2 # Bottom wall
  omega[, Ny] <- -2 * psi[, Ny - 1] / dy^2 # Top wall

  setTxtProgressBar(pb, step)
}
close(pb)
cat("\nSimulation finished.\n")


# --------------------------------------------------------------
# 5. Prepare Data for Visualization
# --------------------------------------------------------------

# Create a data frame for the temperature heatmap
temp_df <- expand.grid(x = x, y = y) %>%
  mutate(Temperature = as.vector(T_field))

# Create a data frame for the velocity vector field
# Downsample for a cleaner quiver plot
skip <- 4
x_vec <- x[seq(1, Nx, by = skip)]
y_vec <- y[seq(1, Ny, by = skip)]
u_vec <- u[seq(1, Nx, by = skip), seq(1, Ny, by = skip)]
v_vec <- v[seq(1, Nx, by = skip), seq(1, Ny, by = skip)]

vector_df <- expand.grid(x = x_vec, y = y_vec) %>%
  mutate(
    u = as.vector(u_vec),
    v = as.vector(v_vec),
    # Calculate end points for the arrows
    xend = x + u * 0.2, # Scale factor for arrow length
    yend = y + v * 0.2
  )

# --------------------------------------------------------------
# 6. Create the Combined Plot with ggplot2
# --------------------------------------------------------------

p <- ggplot() +
  # Temperature heatmap as the background
  geom_tile(data = temp_df, aes(x = x, y = y, fill = Temperature)) +
  scale_fill_viridis_c(name = "Temperature") +

  # Velocity vector field overlay
  geom_segment(
    data = vector_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "white",
    alpha = 0.7,
    arrow = arrow(length = unit(0.05, "inches"), type = "closed")
  ) +

  # Aesthetics and labels
  coord_fixed() + # Ensures 1:1 aspect ratio for the plot
  labs(
    title = "Rayleigh-Bénard Convection",
    subtitle = sprintf("Ra = %.0f, Pr = %.2f", Ra, Pr),
    x = "Width",
    y = "Height"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )

# Display the plot
print(p)

# --------------------------------------------------------------
# 7. Save the Plot
# --------------------------------------------------------------
ggsave(
  filename = "Rayleigh_Benard_Convection.png",
  plot = p,
  width = 10,
  height = 5,
  dpi = 300
)

cat("Combined plot saved to 'Rayleigh_Benard_Convection.png'\n")
