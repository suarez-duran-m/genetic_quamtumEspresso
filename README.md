# genetic_quamtumEspresso
Genetic algorithm and Quantum Espresso to Atom Cluster Optimatisation

# Author: Mauricio Suárez-Durán
# Date: July 5, 2024
# Version: 1.0


This python code is an integration of the genetic algorithm and Quantum Espresso libraries to Optimasie Cluster Atoms.

The methodology implemnted here is based on:
"Geometry Optimisation of Aluminium Clusters Using a Genetic Algorithm", 
<https://doi.org/10.1002/1439-7641(20020517)3:5%3C408::AID-CPHC408%3E3.0.CO;2-G>

## The steps

1. Creation of the population: 5 individuals of 5 atoms. Each is assigned coordinates $[0, N^{1/3}]$ scaled by $r_e=2.0$ Å; with $N$ the number of atoms.
2. Evaluation: *espresso code*.
3. Selection: roulette wheel method. The new offspring constitute 80% of the total population. They are constructed as follows:

3.1. The *fitness* is calculated:
 \begin{equation}
     F_i = \exp \left( -rho \left[ \frac{V_i - V_{V_max}} - V_{V_min}} \right})
 \end{equation}
  
  3.2. A cluster *i* is randomly selected and chosen if $F_i > R[0, 1]$. If it is not chosen, proceed with another one. The idea is to find pairs to create the offspring.

  3.3. Mixing: Deaveu and Ho's *cut-and-splice* operator. Rotate each parent, randomly, and split into two halves, $x-y$ plane randomly, and construct the child.

4. Mutate: assigning new coordinates to one third of the atoms.
5. The descendants are relaxed.
6. They are arranged from lowest to highest energy: $N_mathrm{mut}$ + $N^mathrm{old}_mathrm{pop}$ and the lowest energy ones are chosen.
7. Iterate again for $n$ generations. 
~                                             

Translated with DeepL.com (free version)
