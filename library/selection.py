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


  # Function to rotate individual coordinates
  def rotate_individual(self, individual, Rx, Ry):
    # Extract coordinates into a NumPy array
    coords = np.array([individual[key] for key in individual if key.startswith('atom_coord')])

    # Rotate the coordinates: apply Rx first, then Ry
    rotated_coords = np.dot(Rx, coords.T).T  # Transpose to align dimensions, then transpose back
    rotated_coords = np.dot(Ry, rotated_coords.T).T

    # Update the dictionary with the rotated coordinates
    rotated_individual = individual.copy()
    for i, key in enumerate([key for key in individual if key.startswith('atom_coord')]):
        rotated_individual[key] = rotated_coords[i].tolist()

    return rotated_individual


# Define a function to rotate a cluster about two perpendicular axes
  def rotate_cluster(self, cluster, angle_x, angle_y):
    # Convert angles from degrees to radians
    angle_x = np.radians(angle_x)
    angle_y = np.radians(angle_y)
    
    # Create rotation matrices
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(angle_x), -np.sin(angle_x)],
                   [0, np.sin(angle_x), np.cos(angle_x)]])
    
    Ry = np.array([[np.cos(angle_y), 0, np.sin(angle_y)],
                   [0, 1, 0],
                   [-np.sin(angle_y), 0, np.cos(angle_y)]])
    
    # Rotate each individual in the cluster
    rotated_cluster = {key: rotate_individual(cluster[key], Rx, Ry) for key in cluster}
    
    return rotated_cluster


  # Function to extract and sort coordinates by z-coordinate
  def extract_and_sort(cluster):
    coords = []
    for ind in cluster.values():        
        for key in ind:
            if key.startswith('atom_coord'):
                coords.append(ind[key])
    coords = np.array(coords)    
    sorted_coords = coords[np.argsort(-coords[:, 2])]    
    return sorted_coords


  # Function to perform the mating process
  def mate_clusters(parent1, parent2, N, M):
    # Randomly rotate both parent clusters
    angle_x1, angle_y1 = np.random.rand(2) * 2 * np.pi
    angle_x2, angle_y2 = np.random.rand(2) * 2 * np.pi
    
    rotated_parent1 = rotate_cluster(parent1, angle_x1, angle_y1)
    rotated_parent2 = rotate_cluster(parent2, angle_x2, angle_y2)
    
    sorted_parent1 = extract_and_sort(rotated_parent1)
    sorted_parent2 = extract_and_sort(rotated_parent2)
      
    # Select highest N-M atoms from the first parent
    selected_parent1 = sorted_parent1[:N-M]
    
    # Select lowest M atoms from the second parent
    selected_parent2 = sorted_parent2[-M:]
    
    # Combine the selected fragments to form the child cluster
    child_cluster_coords = np.vstack((selected_parent1, selected_parent2))
    
    # Convert back to dictionary format similar to parent clusters
    child_cluster = {}
    for i, coord in enumerate(child_cluster_coords):        
        child_cluster[f'atom_coord{i+1}'] = coord.tolist()
    
    return child_cluster


  # Define the crossover function
  def crossover(self, parent1, parent2, num_atoms):
    child = {}
    return mate_cluster(parent1, parent2, num_atoms, \
      np.random.randint(1, num_atoms))

    #for atom_num in range(1, num_atoms+1):
    #  atom_key = f'atom_coord{atom_num}'
    #  crossover_point = np.random.randint(1, 3)
    #  child[atom_key] = parent1[atom_key][:crossover_point] \
    #    + parent2[atom_key][crossover_point:]
    #return child


  # Mutation
  def mutate(self, individual, mutation_rate):
    for atom_key in individual.keys():
      if np.random.rand() < mutation_rate:
        individual[atom_key] = [coord + np.random.normal(0, 0.1) \
          for coord in individual[atom_key]]
    return individual 
