# Stochastic Repressilator Simulation

This script (`repressilator_stochastic.ipyb`) simulates the Elowitz & Leibler (2000) Repressilator incorporating molecular noise, using Stochastic Differential Equations (SDEs).

## Script Function

*   Implements the 6-component model with SDEs for protein dynamics, solved via the Euler-Maruyama method.
*   Noise is state-dependent (proportional to `sqrt(protein_concentration)`).
*   Generates and saves a plot of one trajectory showing irregular protein oscillations to `figures/stochastic_oscillation.png`.
*   Attempts to calculate and print period and amplitude statistics for LacI (which will show variability).

## How to Run

1.  Ensure prerequisites are installed (see main `requirements.txt` and the overview README in the parent directory).
2.  Navigate to this `stochastic_model/` directory.
3.  Execute:
    ```bash
    python repressilator_stochastic.py
    ```

## Key Parameters (Defaults in Script)

*   `beta = 0.2`
*   `n = 2`
*   `alpha = 50`
*   `alpha_0 = alpha * 0.001`
*   `noise_coeff = 0.5`

*For detailed parameter justifications and a full model explanation, please see the `README_repressilator_overview.md` in the parent `elowitz_leibler_2000_repressilator/` directory.*

## Expected Output

A plot (`figures/stochastic_oscillation.png`) showing irregular oscillations with variability in amplitude and period for the three repressor proteins.
