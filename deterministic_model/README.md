# Deterministic Repressilator Simulation

This script (`repressilator_deterministic.py`) simulates the idealized, noise-free behavior of the Elowitz & Leibler (2000) Repressilator using Ordinary Differential Equations (ODEs).

## Script Function

*   Implements the 6-ODE model for mRNA and protein concentrations of LacI, TetR, and λ CI.
*   Solves the system using `scipy.integrate.solve_ivp`.
*   Generates and saves a plot of the smooth, periodic protein oscillations to `figures/deterministic_oscillation.png`.
*   Prints calculated period and amplitude statistics for LacI.

## How to Run

1.  Ensure prerequisites are installed (see main `requirements.txt` and the overview README in the parent directory).
2.  Navigate to this `deterministic_model/` directory.
3.  Execute:
    ```bash
    python repressilator_deterministic.py
    ```

## Key Parameters (Defaults in Script)

*   `beta = 0.2`
*   `n = 2`
*   `alpha = 50`
*   `alpha_0 = alpha * 0.001`

*For detailed parameter justifications and a full model explanation, please see the `README_repressilator_overview.md` in the parent `elowitz_leibler_2000_repressilator/` directory.*

## Expected Output

A plot (`figures/deterministic_oscillation.png`) showing regular, phase-shifted oscillations of the three repressor proteins.
