Problem 3: The "Turing Pattern" (Reaction-Diffusion)

Difficulty: Intermediate
Visual Tech: ggplot2 (High-Res Heatmap)
Concept: Pattern Formation, Gray-Scott Model, Numerical Stability.

    Prompt: "Write an R script to simulate the Gray-Scott Reaction-Diffusion model on a 100×100100×100 2D grid.
        Equations: The system involves two chemicals, UU and VV, diffusing at different rates and reacting.
            dU/dt=Du∇2U−UV2+F(1−U)dU/dt=Du​∇2U−UV2+F(1−U)
            dV/dt=Dv∇2V+UV2−(F+k)VdV/dt=Dv​∇2V+UV2−(F+k)V
        Parameters: Du=0.16,Dv=0.08,F=0.035,k=0.060Du​=0.16,Dv​=0.08,F=0.035,k=0.060 (The 'Coral' preset).
        Initial Condition: U=1,V=0U=1,V=0 everywhere, except for a small 10×1010×10 square of noise in the center where V=1V=1.
        Output: Run for 2,000 iterations. Plot the concentration of chemical VV using ggplot2 with scale_fill_viridis_c().
        Visual Goal: The result must be a biologically organic, coral-like texture or fingerprint pattern."

