# --------------------------------------------------------------
# 1.  Load required packages
# --------------------------------------------------------------
library(ggplot2) # for plotting
library(dplyr) # for data manipulation
library(tidyr) # for reshaping (optional, but nice)

# --------------------------------------------------------------
# 2.  Physical & numerical parameters
# --------------------------------------------------------------
# Domain
x_min <- 0.0
x_max <- 2.0
Nx <- 200 # number of spatial grid points
dx <- (x_max - x_min) / Nx

# Time
t_min <- 0.0
t_max <- 2.0 # long enough to see a few periods
CFL <- 0.8 # Courant number (<=1 for upwind stability)
c <- 1.0 # advection speed
dt <- CFL * dx / c # time step
Nt <- ceiling((t_max - t_min) / dt) # number of time steps

# Grid vectors
x_grid <- seq(x_min, x_max - dx, length.out = Nx) # periodic, so last point omitted
t_grid <- seq(t_min, t_max, length.out = Nt)

# --------------------------------------------------------------
# 3.  Initial condition: Gaussian pulse
# --------------------------------------------------------------
x0 <- 0.5 # centre of the pulse
sigma <- 0.05 # width
u0 <- exp(-((x_grid - x0)^2) / (2 * sigma^2))

# --------------------------------------------------------------
# 4.  Allocate storage matrix
# --------------------------------------------------------------
# Rows = time, Columns = space
U_mat <- matrix(0, nrow = Nt, ncol = Nx)
U_mat[1, ] <- u0 # t = 0

# --------------------------------------------------------------
# 5.  Time integration (first‑order upwind, periodic BC)
# --------------------------------------------------------------
for (n in 1:(Nt - 1)) {
  # current solution
  u_now <- U_mat[n, ]

  # upwind derivative (c>0 => look left)
  # periodic: index 1 gets value from Nx
  u_left <- c(u_now[Nx], u_now[1:(Nx - 1)])

  # spatial derivative ∂u/∂x ≈ (u_now - u_left) / dx
  dudx <- (u_now - u_left) / dx

  # update using forward Euler
  u_next <- u_now - c * dt * dudx

  # store
  U_mat[n + 1, ] <- u_next
}

# --------------------------------------------------------------
# 6.  Convert matrix → long data.frame for ggplot2
# --------------------------------------------------------------
# Create a data.frame with columns: t, x, u
df <- expand.grid(
  t = t_grid,
  x = x_grid
) %>%
  mutate(
    u = as.vector(U_mat) # fills column‑wise (t changes fastest)
  )

# (optional) sanity check
head(df)
tail(df)

# --------------------------------------------------------------
# 7.  Plot the space‑time heatmap
# --------------------------------------------------------------
p <- ggplot(df, aes(x = x, y = t, fill = u)) +
  geom_raster() + # or geom_tile() – both work
  scale_fill_viridis_c(
    name = "Amplitude (u)",
    option = "inferno" # feel free to try "magma", "plasma", etc.
  ) +
  labs(
    title = "Space‑Time Waterfall of a 1‑D Advected Gaussian Pulse",
    subtitle = sprintf(
      "c = %.1f, CFL = %.2f, Nx = %d, Nt = %d",
      c,
      CFL,
      Nx,
      Nt
    ),
    x = "Space (x)",
    y = "Time (t)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )

# Display
print(p)

# --------------------------------------------------------------
# 8.  Save the figure (optional)
# --------------------------------------------------------------
ggsave(
  filename = "SpaceTime_Waterfall.png",
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)

# --------------------------------------------------------------
# End of script
# --------------------------------------------------------------
