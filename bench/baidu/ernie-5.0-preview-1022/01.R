# --------------------------------------------------------------------------
# Problem 1: The "Space-Time" Waterfall (1D Advection)
# --------------------------------------------------------------------------

# Load necessary libraries
# install.packages(c("ggplot2", "viridis", "dplyr", "tidyr"))
library(ggplot2)
library(viridis) # For a perceptually uniform color scale
library(dplyr) # For data manipulation
library(tidyr) # For reshaping data

# 1. SETUP: Define simulation parameters
# -------------------------------------
# Physics parameters
L <- 2.0 # Domain length (x in [0, L])
c <- 1.0 # Advection velocity

# Numerical parameters
nx <- 201 # Number of spatial points (grid resolution)
nt <- 401 # Number of time steps to simulate

# Derived parameters
dx <- L / (nx - 1) # Spatial step size
# For stability of the FTCS scheme, the Courant number C = c*dt/dx must be <= 1.
# We'll set it to exactly 1 for a perfect, non-diffusive solution.
C <- 1.0
dt <- C * dx / c # Time step size
T_final <- nt * dt # Total simulation time

# 2. INITIALIZATION: Set up the grid and initial state
# ----------------------------------------------------
# Spatial grid
x <- seq(0, L, length.out = nx)

# Initial condition: A Gaussian pulse centered at x = L/2
u0 <- exp(-200 * (x - L / 2)^2)

# Storage: Create a matrix to store the state at every time step
# Rows are time steps, columns are spatial points
u_history <- matrix(0, nrow = nt, ncol = nx)
u_history[1, ] <- u0 # Store the initial condition

# 3. SIMULATION: Time-stepping loop
# ---------------------------------
# We use the Forward-Time, Central-Space (FTCS) finite difference scheme:
# u(i, n+1) = u(i, n) - 0.5 * C * (u(i+1, n) - u(i-1, n))

# Loop from the second time step to the end
for (n in 1:(nt - 1)) {
  # Get the current state
  u_current <- u_history[n, ]
  u_next <- rep(0, nx)

  # Loop over all spatial points (interior and boundaries)
  for (i in 1:nx) {
    # --- Periodic Boundary Conditions ---
    # Handle the left neighbor (i-1)
    if (i == 1) {
      u_left <- u_current[nx] # Wrap around to the last point
    } else {
      u_left <- u_current[i - 1]
    }

    # Handle the right neighbor (i+1)
    if (i == nx) {
      u_right <- u_current[1] # Wrap around to the first point
    } else {
      u_right <- u_current[i + 1]
    }

    # --- FTCS Update Rule ---
    u_next[i] <- u_current[i] - 0.5 * C * (u_right - u_left)
  }

  # Store the new state in the history matrix
  u_history[n + 1, ] <- u_next
}

# 4. VISUALIZATION: Create the Space-Time Heatmap
# -----------------------------------------------
# The 'u_history' matrix is in the wrong format for ggplot2.
# We need a "long" data frame with columns: x, t, u.

# Create a time vector corresponding to each row in u_history
time_vector <- seq(0, T_final, by = dt)

# Reshape the matrix into a long data frame
# pivot_longer is a modern and clean way to do this
spacetime_df <- u_history %>%
  as.data.frame() %>%
  # Add a time column
  mutate(t = time_vector) %>%
  pivot_longer(
    cols = -t, # Select all columns except 't'
    names_to = "x_idx", # New column for the original column names
    values_to = "u" # New column for the cell values
  ) %>%
  mutate(
    # Convert the 'x_idx' (e.g., "V1", "V2") to a numeric index
    x_idx = as.numeric(sub("V", "", x_idx)),
    # Map the index to the actual x-coordinate
    x = x[x_idx]
  ) %>%
  select(x, t, u) # Keep only the essential columns

# Create the plot
spacetime_plot <- ggplot(spacetime_df, aes(x = x, y = t, fill = u)) +
  geom_raster(interpolate = TRUE) + # Use geom_raster for performance
  scale_fill_viridis_c(option = "plasma", name = "Amplitude (u)") +
  coord_fixed(ratio = 1 / (L / T_final)) + # Make the wave's slope visually correct (dx/dt = c)
  labs(
    title = "1D Advection of a Gaussian Pulse",
    subtitle = "Space-Time Diagram (Waterfall Plot)",
    x = "Space (x)",
    y = "Time (t)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Display the plot
print(spacetime_plot)

# Optional: Save the plot to a file
# ggsave("space_time_waterfall.png", plot = spacetime_plot, width = 8, height = 6)
