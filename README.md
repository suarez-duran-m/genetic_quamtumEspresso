# genetic_quamtumEspresso
Genetic algorithm and Quantum Espresso to Atom Cluster Optimatisation

# Author: Mauricio Suárez-Durán
# Date: July 5, 2024
# Version: 1.0


This python code is an integration of the genetic algorithm and Quantum Espresso libraries to Optimasie Cluster Atoms.

The methodology implemnted here is based on:
"Geometry Optimisation of Aluminium Clusters Using a Genetic Algorithm", 
<https://doi.org/10.1002/1439-7641(20020517)3:5%3C408::AID-CPHC408%3E3.0.CO;2-G>

## The Steps

1. **Creation of the population**: 5 individuals of 5 atoms each. Each atom is assigned coordinates $[0, N^{1/3}]$ scaled by $r_e=2.0$ Å, where $N$ is the number of atoms.
2. **Evaluation**: Using the *espresso code*.
3. **Selection**: Implemented with the roulette wheel method. The new offspring constitute 80% of the total population. They are constructed as follows:
    1. The *fitness* is calculated:
        \[
        F_i = \exp \left( -\rho \left[ \frac{V_i - V_{V_{\text{max}}}}{V_{V_{\text{min}}}} \right] \right)
        \]
    2. A cluster *i* is randomly selected and chosen if $F_i > R[0, 1]$. If it is not chosen, proceed with another one. The idea is to find pairs to create the offspring.
    3. **Mixing**: Using Deaveu and Ho's *cut-and-splice* operator. Rotate each parent randomly, split into two halves along the $x-y$ plane randomly, and construct the child.
4. **Mutate**: Assigning new coordinates to one-third of the atoms.
5. **Relaxation**: The descendants are relaxed.
6. **Arrangement**: The individuals are arranged from lowest to highest energy. This includes $N_{\text{mut}}$ + $N^{\text{old}}_{\text{pop}}$, and the lowest energy ones are chosen.
7. **Iteration**: Repeat the process for $n$ generations.
