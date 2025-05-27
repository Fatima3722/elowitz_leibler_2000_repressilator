# elowitz_leibler_2000_repressilator

This directory contains Python code and generated figures for the replication of key aspects of the Repressilator, a synthetic genetic oscillator, as originally described in:

**Original Publication:**
> Elowitz, M. B., & Leibler, S. (2000). A synthetic oscillatory network of transcriptional regulators. *Nature*, *403*(6767), 335-338. DOI: [10.1038/35002125](https://doi.org/10.1038/35002125)

## Model Overview

The Repressilator is a pioneering synthetic genetic circuit consisting of three transcriptional repressor genes (LacI, TetR, and λ CI from bacteriophage lambda). These genes are arranged in a cyclic negative feedback loop: LacI represses TetR, TetR represses CI, and CI represses LacI. This design can lead to sustained oscillations in the concentrations of the three repressor proteins. The work by Elowitz and Leibler was a landmark demonstration of the rational design of dynamic behavior in living cells and highlighted the inherent stochasticity in such systems.

## Replication Goals

The Python scripts in this directory aim to:
1.  Simulate the **deterministic behavior** of the Repressilator model using Ordinary Differential Equations (ODEs) based on the parameters and model structure described in Box 1 of the paper.
2.  Simulate a **stochastic version** of the Repressilator, incorporating state-dependent noise via the Euler-Maruyama method, to illustrate the impact of molecular noise on the oscillations, a key experimental finding of the paper.

## Key Parameters Used in Simulations

These parameters are chosen to align with the dimensionless model presented in Box 1 and the illustrative parameters in Figure 1c of the Elowitz & Leibler (2000) paper:
*   `beta` (Ratio of protein degradation rate to mRNA degradation rate): **0.2**
    *   (Derived from protein half-life ≈ 10 min, mRNA half-life ≈ 2 min)
*   `n` (Hill coefficient for repression): **2**
*   `alpha` (Dimensionless maximum transcription rate): **50** (This value influences amplitude and was chosen to give clear oscillations)
*   `alpha_0` (Dimensionless basal/leaky transcription rate): `alpha * 0.001`
*   For the stochastic simulation:
    *   `noise_coeff` (Scales the state-dependent noise term): **0.5** (This value was chosen to show clear noise effects without completely destabilizing oscillations)

Time in the plots is converted to minutes, assuming 1 unit of dimensionless ODE time corresponds to 2 minutes (approximating the mRNA half-life for visualization). Protein concentrations are presented in dimensionless units, relative to K<sub>M</sub> (the repressor concentration for half-maximal repression).

