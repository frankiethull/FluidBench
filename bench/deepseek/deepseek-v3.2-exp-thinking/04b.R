library(rayshader)
library(Rcpp)
library(RcppArmadillo)

# Use Rcpp for faster computation (optional but recommended)
sourceCpp(
  code = '
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::export]]
mat simulate_ripples_cpp(int nx, int ny, int nt, double dt, double dx) {
  mat h = ones<mat>(nx, ny);
  mat h_prev = h;
  mat h_curr = h;
  mat h_next = h;
  
  // Initial Gaussian drop
  int center_x = nx/2;
  int center_y = ny/2;
  for(int i=0; i<nx; i++) {
    for(int j=0; j<ny; j++) {
      double r2 = pow((i-center_x)/10.0, 2) + pow((j-center_y)/10.0, 2);
      h_curr(i,j) = 1.0 - 0.3 * exp(-r2 / 0.3);
    }
  }
  h_prev = h_curr;
  
  double c = 0.8; // Wave speed
  
  for(int n=0; n<nt; n++) {
    for(int i=1; i<nx-1; i++) {
      for(int j=1; j<ny-1; j++) {
        double laplacian = h_curr(i+1,j) + h_curr(i-1,j) + 
                          h_curr(i,j+1) + h_curr(i,j-1) - 4*h_curr(i,j);
        h_next(i,j) = 2*h_curr(i,j) - h_prev(i,j) + pow(c*dt/dx, 2) * laplacian;
      }
    }
    
    // Reflective boundaries
    h_next.row(0) = h_next.row(1);
    h_next.row(nx-1) = h_next.row(nx-2);
    h_next.col(0) = h_next.col(1);
    h_next.col(ny-1) = h_next.col(ny-2);
    
    h_prev = h_curr;
    h_curr = h_next;
    
    if(n % 100 == 0) Rprintf("Step %d of %d\\n", n, nt);
  }
  
  return h_curr;
}'
)

# Run C++ simulation (much faster)
cat("Running optimized C++ simulation...\n")
h_final <- simulate_ripples_cpp(200, 200, 800, 0.01, 0.05)

# Create photorealistic render
cat("Creating photorealistic render...\n")

water_surface <- h_final

# Enhanced rendering with multiple effects
water_surface %>%
  height_shade(
    texture = colorRampPalette(
      c("#000033", "#000066", "#000099", "#0000CC", "#0066FF", "#66CCFF")
    )(256)
  ) %>%
  add_overlay(
    sphere_shade(
      water_surface,
      texture = "imhof1",
      colorintensity = 3,
      progbar = FALSE
    ),
    alphalayer = 0.6
  ) %>%
  add_shadow(
    ray_shade(water_surface, zscale = 0.05, sunangle = 45, lambert = TRUE),
    0.4
  ) %>%
  add_shadow(ambient_shade(water_surface, zscale = 0.05), 0.2) %>%
  plot_3d(
    water_surface,
    zscale = 0.05,
    water = TRUE,
    waterdepth = quantile(water_surface, 0.3)[1],
    wateralpha = 0.8,
    watercolor = "#88CCFF",
    theta = -25,
    phi = 40,
    windowsize = c(1400, 1000),
    background = "#88BBFF",
    soliddepth = -max(water_surface)
  )

# High quality render with depth of field
render_highquality(
  filename = "photorealistic_pond.png",
  samples = 256,
  light = TRUE,
  lightdirection = c(225, 135),
  lightcolor = c("#FFFFFF", "#88CCFF"),
  lightintensity = c(1000, 300),
  aperture = 16,
  focal_distance = 100,
  width = 2000,
  height = 1500,
  tonemap = "reinhold"
)

cat("Photorealistic render complete! Check 'photorealistic_pond.png'\\n")
