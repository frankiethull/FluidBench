Problem 5: Rayleigh-Bénard Convection (Boussinesq approx)

Difficulty: Expert
Visual Tech: plotly (Interactive Vector Fields) or grid arrows.
Concept: Coupled Temperature-Velocity, Buoyancy, Vector Fields.

    Prompt: "Write an R script to simulate Rayleigh-Bénard Convection (a fluid heated from below and cooled from above).
        Physics: Solve the 2D Incompressible Navier-Stokes equations coupled with an Energy equation for Temperature (TT).
        Approximation: Use the Boussinesq approximation, where fluid density depends on Temperature, creating a buoyant force term in the Y-velocity equation (Fy=α⋅g⋅TFy​=α⋅g⋅T).
        Domain: Aspect ratio 2:1. Bottom wall T=1T=1 (Hot), Top wall T=0T=0 (Cold).
        Output: Create a combined plot using ggplot2 or plotly:
            A Heatmap background representing Temperature.
            A Vector Field (quiver plot) overlay representing Velocity arrows.
        Visual Goal: The visualization must show 'convection rolls' (circular currents) where hot fluid rises in plumes and cold fluid sinks."
