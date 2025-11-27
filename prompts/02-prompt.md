Problem 2: Von Kármán Vortex Street (Lattice Boltzmann)

Difficulty: Intermediate
Visual Tech: gganimate (or a GIF loop)
Concept: Mesoscopic Simulation, Obstacles, and Animation.

    Prompt: "Write an R script to implement a 2D Lattice Boltzmann (D2Q9) simulation of fluid flowing past a circular obstacle.
        Grid: 80×4080×40 lattice.
        Setup: A fixed circular obstacle (radius 3) is placed at (x=20,y=20)(x=20,y=20). Fluid flows from Left to Right.
        Optimization: Use R matrix operations (vectorization) for the streaming step.
        Physics: Calculate the Curl (Vorticity) of the flow field at each step.
        Output: Create an animated GIF using gganimate (or a loop saving .png files).
        Visual Goal: The animation must clearly show alternating vortices shedding off the back of the cylinder (the Von Kármán vortex street), colored by vorticity (red for clockwise, blue for counter-clockwise)."
