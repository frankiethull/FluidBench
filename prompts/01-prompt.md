Problem 1: The "Space-Time" Waterfall (1D Advection)

Difficulty: Basic
Visual Tech: ggplot2 (geom_raster/geom_tile)
Concept: 1D Linear Advection + Data Reshaping.

    Prompt: "Write an R script to simulate 1D Advection (∂u∂t+c∂u∂x=0∂t∂u​+c∂x∂u​=0) of a Gaussian pulse.
        Physics: Domain x∈[0,2]x∈[0,2]. Velocity c=1c=1. Periodic boundaries (xendxend​ connects to xstartxstart​).
        Storage: Do not just plot the final state. You must store the state of the wave at every time step in a matrix or dataframe.
        Output: Create a Space-Time Heatmap using ggplot2.
            X-axis: Space (xx)
            Y-axis: Time (tt)
            Color: Wave Amplitude (uu)
        Visual Goal: The result should look like diagonal stripes, showing the wave traveling through time."
