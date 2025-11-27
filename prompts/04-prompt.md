Problem 4: The "Pond Ripple" (2D Shallow Water Equations)

Difficulty: Advanced
Visual Tech: rayshader (3D Rendering)
Concept: Hyperbolic Systems, 3D Surface Rendering, Raytracing.

    Prompt: "Write an R script to simulate the 2D Shallow Water Equations (simulating water height hh in a tank).
        Initial Condition: A flat water surface with a Gaussian 'drop' added to the center height.
        Solver: Use the Lax-Wendroff or Leapfrog scheme to solve for height hh and momentum (hu,hv)(hu,hv).
        Visualization: Run the simulation until the ripples reach the edge of the domain.
        Output: Use the rayshader library to render the final water height matrix as a 3D surface.
            Include a water texture or blue color palette.
            Add a light source to create shadows/refraction effects on the waves.
        Visual Goal: A photorealistic (or near-photorealistic) 3D render of a ripple ring expanding in a pool."
