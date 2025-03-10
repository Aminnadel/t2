# t2

Overview

This project provides numerical solutions to the ordinary differential equation (ODE):
\[ \frac{dy}{dt} = -50y + \sin(t), \quad y(0) = 1 \]

using various numerical methods, including:
- Forward Euler Method
- Modified Euler Method
- Backward Euler Method
- Runge-Kutta 4th Order (RK4)
- Adams-Bashforth 2-step Method
- Adams-Moulton 2-step Method

Requirements
This project is designed to run in **GNU Octave**. To avoid rendering issues, you **must** set the correct graphics toolkit before running the script:

graphics_toolkit("gnuplot")

Usage
1. Open GNU Octave.
2. Set the required graphics toolkit using:
   graphics_toolkit("gnuplot")
3. Run the script:
   numerical_ode_extended(h)
   where `h` is the step size (e.g., `0.01`).

Troubleshooting
- If you encounter rendering errors, ensure `gnuplot` is installed and selected as the graphics toolkit.
- If performance issues arise, try using a larger step size (e.g., `h = 0.05`).


