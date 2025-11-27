# 1D Linear Advection of a Gaussian Pulse
# ∂u/∂t + c ∂u/∂x = 0  with periodic boundaries

library(ggplot2)

#-----------------------------
# Parameters and grids
#-----------------------------
c <- 1 # advection speed
x_min <- 0
x_max <- 2
dx <- 0.02
x <- seq(x_min, x_max - dx, by = dx) # avoid double-counting x = 2
Nx <- length(x)

dt <- 0.01
t_max <- 4
time <- seq(0, t_max, by = dt)
Nt <- length(time)

# CFL number (for reference)
CFL <- c * dt / dx # should be <= 1 for this explicit scheme

#-----------------------------
# Initial condition: Gaussian pulse
#-----------------------------
x0 <- 0.5 # center of Gaussian
sigma <- 0.1
u <- exp(-(x - x0)^2 / (2 * sigma^2))

#-----------------------------
# Time stepping (upwind scheme, periodic BC)
#-----------------------------
u_store <- matrix(NA_real_, nrow = Nx, ncol = Nt) # rows: x, cols: t
u_store[, 1] <- u

for (n in 2:Nt) {
  u_new <- u

  # Interior points (i = 2..Nx)
  u_new[2:Nx] <- u[2:Nx] - CFL * (u[2:Nx] - u[1:(Nx - 1)])

  # Periodic boundary: point i = 1 uses i = Nx as "left" neighbor
  u_new[1] <- u[1] - CFL * (u[1] - u[Nx])

  u <- u_new
  u_store[, n] <- u
}

#-----------------------------
# Reshape data to long format for ggplot2
#-----------------------------
df <- expand.grid(
  x = x,
  time = time
)
df$u <- as.vector(u_store) # matches expand.grid(x, time) ordering

#-----------------------------
# Space-Time heatmap
#-----------------------------
ggplot(df, aes(x = x, y = time, fill = u)) +
  geom_raster() +
  scale_fill_viridis_c(option = "plasma") +
  labs(
    title = "1D Linear Advection of a Gaussian Pulse (Space–Time Diagram)",
    x = "Space x",
    y = "Time t",
    fill = "u(x, t)"
  ) +
  theme_minimal()
