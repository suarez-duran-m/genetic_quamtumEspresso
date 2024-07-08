import math
import random
import itertools
import numpy as np

class selection:

  def __init__(self):
    self.a = 0

  def calculate_fitness(self, V_i=1.0, V_min=1.0, V_max=1.0, rho=1.0):
    # Ensure V_max is not equal to V_min to avoid division by zero
    if V_max == V_min:
      fitness = 0
      #raise ValueError("V_max should not be equal to V_min to avoid division by zero.")
    else:
      # Calculate the fitness
      fitness = math.exp(-rho * (V_i - V_min) / (V_max - V_min))

    return fitness


  # Function to select parents for mating
  def select_unique_pairs_for_mating(self, clusters='nan', num_pairs_needed=1):
    selected_pairs = []
    all_pairs = list(itertools.combinations(clusters.keys(), 2))

    while len(selected_pairs) < num_pairs_needed:
      # Pick a pair at random
      pair = random.choice(all_pairs)
      
      # Check if both individuals in the pair have sufficient fitness
      F_i = [clusters[key]['fitness'] for key in pair]

      # Generate random numbers between 0 and 1 for both individuals
      R = [random.random() for _ in pair]
    
      # Accept the pair for mating if both F_i > R
      if all(F > r for F, r in zip(F_i, R)):
        if pair not in selected_pairs:
          selected_pairs.append(pair)
          # Remove the selected pair from all_pairs to avoid duplicates
          all_pairs.remove(pair)

    return selected_pairs

  # Define the crossover function
  def crossover(self, parent1, parent2, num_atoms):
    child = {}
    for atom_num in range(1, num_atoms+1):
      atom_key = f'atom_coord{atom_num}'
      crossover_point = np.random.randint(1, 3)
      child[atom_key] = parent1[atom_key][:crossover_point] \
        + parent2[atom_key][crossover_point:]
    return child


  def mutate(self, individual, mutation_rate):
    for atom_key in individual.keys():
      if np.random.rand() < mutation_rate:
        individual[atom_key] = [coord + np.random.normal(0, 0.1) \
          for coord in individual[atom_key]]
    return individual

    
