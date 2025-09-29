## Synthetic Control and Power Calculations Replication

This repository contains the material for our Applied Economics course.
It replicates figures and methodology from Danilo Freire’s article:

> Freire, D. (2018). *Evaluating the Effect of Homicide Prevention Strategies in São Paulo, Brazil: A Synthetic Control Approach*.  
> Latin American Research Review, 53(2), 231–249. https://doi.org/10.25222/larr.334

---

### Overview  

1. **Synthetic Control (São Paulo case)**  
   - We apply the synthetic control method (Abadie & Gardeazabal, 2003; Abadie, Diamond & Hainmueller, 2010) to replicate Freire’s main figures (except Figure 8).  
   - The method builds a “synthetic São Paulo” from a weighted combination of Brazilian states that did not implement the same policies.  
   - By comparing observed homicide rates with this counterfactual, we estimate the aggregate effect of public safety measures adopted after 1999.  
   - Robustness checks include placebo interventions, leave-one-out exercises, and permutation tests.  

2. **Power Calculations**  
   - Monte Carlo simulations illustrate how power depends on sample size, effect size, error variance, and treatment assignment share.  
   - We replicate the tutorial example, then explore alternative settings (higher variance, unbalanced treatment assignment, inclusion of covariates).  
   - Results are presented in static graphs and briefly discussed.  

We used **Stata**  with the package `synth` (for synthetic control estimation).


