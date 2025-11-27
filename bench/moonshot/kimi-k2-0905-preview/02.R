# 2-D Von Kármán vortex street – Lattice-Boltzmann (D2Q9) in base R
# ------------------------------------------------------------------
# Grid: 80×40  |  Obstacle: circle r = 3 centred at (20,20)
# Output: animated GIF of vorticity (curl) via gganimate
# NOTE: pure-R, fully vectorised streaming step

library(tidyverse)
library(gganimate)
library(gifski)          # fast GIF writer

# -------------------- LBM parameters --------------------
nx   <- 80
ny   <- 40
r    <- 3                # cylinder radius
cx   <- 20
cy   <- ny/2
Re   <- 80               # Reynolds number
Umax <- 0.05             # max inlet speed
nu   <- Umax*r/Re        # kinematic viscosity
tau  <- 3*nu + 0.5       # relaxation time
omega <- 1/tau

# D2Q9 lattice constants
e <- matrix(c(0,1,0,-1,0,1,-1,-1,1,
              0,0,1,0,-1,1,1,-1,-1), nrow=2)
w <- c(4,1,1,1,1,1,1,1,1)/9
w[1] <- 4/9

# initialise
rho  <- matrix(1, nx, ny)
ux   <- matrix(0, nx, ny)
uy   <- matrix(0, nx, ny)

# helper: equilibrium distributions
feq <- function(r, u){
  eu <- e[1,]*u$x + e[2,]*u$y
  u2 <- u$x^2 + u$y^2
  r*w*(1 + 3*eu + 4.5*eu^2 - 1.5*u2)
}

# obstacle mask (TRUE = solid)
obstacle <- matrix(FALSE, nx, ny)
for(i in 1:nx) for(j in 1:ny)
  obstacle[i,j] <- (i-cx)^2 + (j-cy)^2 <= r^2

# -------------------- time loop --------------------
f <- array(1, dim=c(nx,ny,9))       # populations
f[] <- feq(1, list(x=0,y=0))       # initialise to rest

frames <- 400                       # total frames
stepEvery <- 10                     # only store every 10 steps

storage <- list()

for(step in 1:(frames*stepEvery)){
  # ---------------- collision ----------------
  rho <- apply(f, c(1,2), sum)
  ux  <- apply(f*aperm(array(e[1,], dim=c(1,1,9)), c(3,1,2)), c(1,2), sum)/rho
  uy  <- apply(f*aperm(array(e[2,], dim=c(1,1,9)), c(3,1,2)), c(1,2), sum)/rho

  # bounce-back inside obstacle
  solid <- which(obstacle, arr.ind=TRUE)
  f[solid, 2:9] <- f[solid, c(2:9,2:9)[c(1:8)+8]]   # simple bounce-back

  # equilibrium & collide
  f <- f - omega*(f - feq(rho, list(x=ux, y=uy)))

  # ---------------- streaming (vectorised) ----------------
  for(k in 1:9){
    f[,,k] <- c(f[ c((1:nx)-e[1,k]) %% nx + 1,
                   c((1:ny)-e[2,k]) %% ny + 1 , k])
  }

  # ---------------- store frame ----------------
  if(step %% stepEvery == 0){
    # curl = ∂uy/∂x - ∂ux/∂y  (2nd-order central)
    curl <- ux*NA
    for(i in 1:nx) for(j in 1:ny){
      ip <- i %% nx + 1; im <- (i-2) %% nx + 1
      jp <- j %% ny + 1; jm <- (j-2) %% ny + 1
      curl[i,j] <- (uy[ip,j] - uy[im,j])/(2*dx) -
                    (ux[i,jp] - ux[i,jm])/(2*dy)
    }
    storage[[step/stepEvery]] <-
      expand.grid(x=1:nx, y=1:ny) %>%
      mutate(curl = c(curl),
             frame = step/stepEvery)
  }
}

# -------------------- build animation --------------------
dx <- dy <- 1
full <- bind_rows(storage)

p <- ggplot(full, aes(x, y, fill = curl)) +
  geom_tile() +
  scale_fill_gradient2(low="steelblue", mid="white", high="firebrick",
                       midpoint = 0, limits=quantile(full$curl, c(.02,.98))) +
  coord_fixed() +
  theme_void() +
  labs(title = "Von Kármán vortex street – step {closest_state}",
       subtitle = "Colour = vorticity (curl)",
       fill = "curl") +
  transition_states(frame, transition_length = 1, state_length = 1)

animate(p, renderer = gifski_renderer(), fps = 12, width = 