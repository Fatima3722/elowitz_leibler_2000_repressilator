# elowitz_leibler_2000_repressilator

This directory contains Python code and generated figures for the replication of key aspects of the Repressilator, a synthetic genetic oscillator, as originally described in:

**Original Publication:**
> Elowitz, M. B., & Leibler, S. (2000). A synthetic oscillatory network of transcriptional regulators. *Nature*, *403*(6767), 335-338. DOI: [10.1038/35002125](https://doi.org/10.1038/35002125)

Scientists made an artificial network which is made up of proteins (regulators) controls the gene expression, and its shows oscillatory behaviour.


## Introduction to the paper

* Cells  consists diff. biomolecules which make some sorts of network to work together but we still didnot understand the "design principle" till now.
* The main aim is to implement an oscillation they designed and made an artificial network.
* They used three trancriptional repressor system which are not part of natural biological clock and they created an oscillation network from it in ecoli and named it "repressilator".
* This repressilator network periodically make GFP. GFP is a protein which shown flouroscence in green light.So that scientist can know when network is "on".
* This oscillation is about an hour long which is longer than the divison of bacteria that means it get passed onto their offspring.
* The observed some noise in their too. 
 

## Replication Goals

The Python scripts in this directory aim to:
1.  Simulate the **deterministic behavior** of the Repressilator model using Ordinary Differential Equations (ODEs) based on the parameters and model structure described in Box 1 of the paper.
2.  Simulate a **stochastic version** of the Repressilator, incorporating state-dependent noise via the Euler-Maruyama method, to illustrate the impact of molecular noise on the oscillations, a key experimental finding of the paper.

## Network Design

![rps](https://github.com/user-attachments/assets/ec0972be-6a4f-4ce9-ac96-550dd6391bb2)

* LacI represses tetR represses cl(lambda phage) represses LacI.
* This negative feedback loop with time produces oscillations.

## Oscillations are favoured by

* Strong Promoters - For the promoters they choosed Hybrid promoter
* Efficient ribosome binding site
* Tight transcriptional repression, decreases leakiness
* Cooperative Repression (repressor mol. together "off" gene effectively).
* Protein and decay mrna rate are comparable to overcome this problem as the mrna lifetime is 2 min. to match with protein they put SSRA RNA seq. at the carboxy terminal tag so that proteases recognize and degrade them .

## EXP. SETUP
* Low copy plasmid of (Repressilator) so that it wont interfere normal functioning of ecoli.
* High copy no. of Reporter plasmid for visualisation.
* IPTG interferes LacI, They used pulse of IPTG to synchronize all cells at one start point.
* When they visualised Ecoli **MC4100** they saw damped oscillation when they grew it in IPTG media and transferred back into normal media.
* As population dont have apparent way to sync. So they used microscope to study single cell their flouroscence intensity.

## Network Structure

* The three repressor proteins (Pi) and their corresponding mRNA conc. are continous variable ( i = LacI, tetR, cl ).

  mRNA Equations

  dmi/dt = -mi + alpha/(1+Pj ^^n) + alpha0

  dmi/dt = mrna conc.
  
  -mi = Degradation of mrna
  
  alpha = Max. production rate
  
  Pj = the repressor protein repressing "i" mrna
  
  alpha0 = Leakiness

  n = Hill coeff. ( Cooperative n>1)

  Protein Equation

  dpi/dt = -beta(Pi - Mi)

  dPi/dt = Change in rate of protein concentration
  Beta = Ratio of protein degradation rate to mRNA degradation rate

  * Time is rescaled in the unit of mrna lifetime

## Objective

My objective was to replicate and analyze the Elowitz and Leibler (2000) Repressilator model. To do this, I developed two computational models in Python: a deterministic ODE-based model to understand the idealized system dynamics, and a stochastic SDE-based model to investigate the impact of molecular noise.

1. Deterministic Model:

  * The deterministic model directly implements the system of six coupled Ordinary Differential Equations outlined in Box 1 of the Elowitz & Leibler paper. These equations describe the rate of change for the concentrations of three mRNAs (m1, m2, m3 for LacI, TetR, λ CI respectively) and their corresponding proteins (p1, p2, p3).
  * The core repressive interactions are modeled using Hill functions of the form α / (1 + P_repressor^n), where P_repressor is the dimensionless concentration of the repressing protein. This term is augmented by α_0 to account for basal transcriptional leakiness.
  * Protein synthesis is proportional to mRNA concentration, and both mRNA and protein species undergo first-order degradation. The system was solved numerically using scipy.integrate.solve_ivp with the 'RK45' adaptive step-size algorithm, which is well-suited for such non-linear ODE systems.
   
## Why I Chose These Specific "Settings" (Parameters):

beta = 0.2: This number, β, is very important. It's the protein breakdown speed divided by the mRNA breakdown speed. In their paper (Figure 1c caption), they mention proteins lasting about 10 minutes and mRNA about 2 minutes. If you calculate that ratio, you get 0.2. Using this value helps get the blinking speed (the period) into a range that makes sense with their experiments.

n = 2 (Hill Coefficient): They used n=2 in their own simulation (Figure 1c left). This n tells you how strongly the proteins 'cooperate' to stop the next gene. A value of 2 means there's some cooperation, which makes the 'on-off' switch a bit sharper and helps the system oscillate better than if n was just 1.

alpha_0 = alpha * 0.001 (Tiny "Leak"): Their Figure 1c caption also implies that even when a gene is 'off,' it's still about 1/1000th as active as when it's fully 'on.' So, I set my α_0 (the leakiness) to be 0.001 times whatever α (the max production rate) was.

alpha = 50 (Max Production Speed): α is how fast mRNA gets made when there's no 'stop' signal. Once I had β, n, and the α_0/α ratio figured out, I chose α=50. Their theory (Figure 1b) suggests that with these other settings, a value like this should make it oscillate. It also gives protein amounts in my simulation that seem reasonable if you think about their K<0xE2><0x82><0x82><0x82>M value of about 40 molecules. The main thing α does here is set how 'high' the blinks go.

Time Scale: When I plot it, the time axis is in minutes. I converted it from the simulation's 'math time' by multiplying by 2. This is like saying each 'math time unit' is 2 minutes long, which matches their mRNA half-life and is a common way to show these plots.

Result-

 we can clearly see the three proteins going up and down in these nice, smooth, regular waves.
 And look how they're out of sync – LacI peaks, then TetR, then CI. It's about a 120-degree difference, which is exactly what we expect from that three-part cycle.
 The blinking happens about every [calculated period, e.g., 65-70] minutes. The speed is mostly set by that β value, but n and α play a part too.
 This plot really confirms that their design, on paper, should totally work to create these oscillations, just like their own theoretical plot (Figure 1c, left) showed.

2. Stochastic Model-

  Building Model

  * Next, I wanted to see what happens when you add the randomness that's always there in real cells. So, I made a 'stochastic' version. The basic math for how proteins are made and 
    stopped is the same as before.
  * But this time, I added a 'random nudge' term to the protein equations. This nudge is noise_coeff * sqrt(P_effective) * dW. P_effective is just the amount of protein (making sure 
    it's not negative), dW is a small random number, and noise_coeff controls how big these random nudges are.
  * Making the noise depend on the square root of the protein amount is a common way to mimic 'intrinsic noise' the idea that when you have more molecules, the random fluctuations can 
    be bigger in absolute terms. It's related to something called the Chemical Langevin Equation.
  * To solve these equations with randomness, I used a method called Euler-Maruyama.

    noise_coeff = 0.5: I set this noise_coeff to 0.5. I tried a few values, and this one was good because it made the randomness clearly visible we could see it messing with the 
    blinking but it didn't completely break the oscillations during the simulation. I wanted it to look a bit like the noisy plots in their paper.

    Results-

    * The height of the blinks (amplitude) is all over the place. Some LacI peaks are tall, some are short.
   * The timing between blinks (period) is also really irregular. It's not a steady rhythm anymore.
   * And the lines themselves are all jagged from those constant random effects.
   * This kind of messy but still sort-of-blinking behavior is exactly what Elowitz and Leibler saw in their actual experiments with bacteria (their Figure 2), and it also looks like 
     their own noisy simulation (Figure 1c, right).


## Conclusion

  * First, my simulation proved that the Repressilator design, based on their math, really should create regular oscillations.
  * Then, my  simulation proved that when you add in the kind of molecular noise that's always present in cells, you get those irregular, variable blinks. This was a huge point in 
    their paper: they designed something to work one way, but the cell's own randomness changed how it actually performed.
  * So, doing these simulations helped me see how their design works in theory, and then how noise makes it behave like it did in their real experiments. It really shows why their 
    paper was so important for understanding both how to build new things in cells and how tricky it can be because of that natural randomness.
 




