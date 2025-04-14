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
      if F_i[0] == 0 and F_i[1] == 0:
        selected_pairs.append(pair)
        break
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
    coords = np.array(individual)

    # Rotate the coordinates: apply Rx first, then Ry
    rotated_coords = np.dot(Rx, coords.T).T  # Transpose to align dimensions, then transpose back
    rotated_coords = np.dot(Ry, rotated_coords.T).T

    return rotated_coords


# Define a function to rotate a cluster about two perpendicular axes
  def rotate_cluster(self, cluster_coord, angle_x, angle_y):
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
    return {key: self.rotate_individual(cluster_coord[key], Rx, Ry) for key in cluster_coord if key.startswith('atom_coord')}

  def rotate_cluster_bimetallic(self, cluster_coord, angle_x, angle_y):
    """
    Rotate a bimetallic cluster. Each key starting with "atom_" holds a dictionary:
      {'coords': [x, y, z], 'spice': <element>}
    The rotation is applied to the coordinate part.
    """
    angle_x = np.radians(angle_x)
    angle_y = np.radians(angle_y)
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(angle_x), -np.sin(angle_x)],
                   [0, np.sin(angle_x), np.cos(angle_x)]])
    Ry = np.array([[np.cos(angle_y), 0, np.sin(angle_y)],
                   [0, 1, 0],
                   [-np.sin(angle_y), 0, np.cos(angle_y)]])
    rotated = {}
    for key, value in cluster_coord.items():
      if key.startswith("atom_"):
        # Extract coordinate vector and rotate it
        coords = np.array(value["coords"])
        rotated_coords = Ry @ (Rx @ coords)
        # Preserve the spice information
        rotated[key] = {"coords": rotated_coords.tolist(), "spice": value["spice"]}
      else:
        # Copy any other items (like energy or fitness)
        rotated[key] = value
    return rotated


  # Function to extract and sort coordinates by z-coordinate
  def extract_and_sort(self, cluster):
    coords = []
    for ind in cluster.values():\
      coords.append(ind)
    coords = np.array(coords)
    sorted_coords = coords[np.argsort(-coords[:, 2])] 
    return sorted_coords

  def extract_and_sort_bimetallic(self, cluster):
    """
    Extract atom dictionaries (with 'coords' key) and sort them in descending order by z-coordinate.
    """
    atoms = []
    for key, value in cluster.items():
      if key.startswith("atom_"):
        atoms.append(value)
    sorted_atoms = sorted(atoms, key=lambda atom: -atom["coords"][2])
    return sorted_atoms

  # Function to perform the mating process
  def mate_clusters(self, parent1, parent2, N, issingle, M):
    # Randomly rotate both parent clusters
    child_cluster = {}
    angle_x1, angle_y1 = np.random.rand(2) * 2 * np.pi
    angle_x2, angle_y2 = np.random.rand(2) * 2 * np.pi

    rotated_parent1 = self.rotate_cluster(parent1, angle_x1, angle_y1)
    rotated_parent2 = self.rotate_cluster(parent2, angle_x2, angle_y2)

    sorted_parent1 = self.extract_and_sort(rotated_parent1)
    sorted_parent2 = self.extract_and_sort(rotated_parent2)

    # Select highest N-M atoms from the first parent
    selected_parent1 = sorted_parent1[:N - M]

    # Select lowest M atoms from the second parent
    selected_parent2 = sorted_parent2[-M:]

    # Combine the selected fragments to form the child cluster
    child_cluster_coords = np.vstack((selected_parent1, selected_parent2))

    # Convert back to dictionary format similar to parent clusters
    for i, coord in enumerate(child_cluster_coords):
        child_cluster[f'atom_coord{i + 1}'] = coord.tolist()
    return child_cluster

  def mate_clusters_bimetallic(self, parent1, parent2, N, M):
    """
    For two-spices clusters:
      - Rotate both parents using rotate_cluster_bimetallic.
      - Extract atoms (each a dict with keys 'coords' and 'spice') from both parents.
      - Sort them by z-coordinate.
      - Select the top (N-M) atoms from parent1 and the bottom M atoms from parent2.
      - Combine these atoms into a new child cluster.
    """
    child_cluster = {}
    angle_x1, angle_y1 = np.random.rand(2) * 360
    angle_x2, angle_y2 = np.random.rand(2) * 360

    rotated_parent1 = self.rotate_cluster_bimetallic(parent1, angle_x1, angle_y1)
    rotated_parent2 = self.rotate_cluster_bimetallic(parent2, angle_x2, angle_y2)

    sorted_parent1 = self.extract_and_sort_bimetallic(rotated_parent1)
    sorted_parent2 = self.extract_and_sort_bimetallic(rotated_parent2)

    selected_parent1 = sorted_parent1[:N - M]
    selected_parent2 = sorted_parent2[-M:]
    all_atoms = selected_parent1 + selected_parent2

    for i, atom in enumerate(all_atoms, start=1):
      child_cluster[f'atom_{i}'] = atom  # atom is already a dict with 'coords' and 'spice'
    return child_cluster

  # Define the crossover function
  def crossover(self, parent1, parent2, num_atoms, issingle=1):
    """
    Perform crossover between two parent clusters.
    - For single-spice clusters (issingle == 1): use mate_clusters.
    - For bimetallic clusters (issingle == 0): use mate_clusters_bimetallic.

    A random integer M in the interval [1, num_atoms) is chosen to split the atoms.
    """
    M = np.random.randint(1, num_atoms)
    if issingle:
      return self.mate_clusters(parent1, parent2, num_atoms, True, M)
    else:
      return self.mate_clusters_bimetallic(parent1, parent2, num_atoms, M)

  # Mutation
  def mutate(self, individual, mutation_rate):
    """
    Mutate an individual (child cluster) with a given mutation rate.

    For single-spice clusters, each key (like "atom_coord1") directly holds a coordinate list.
    For bimetallic (two-spice) clusters, each key (like "atom_1") holds a dictionary containing:
       - "coords": [x, y, z]
       - "spice": element label (unchanged)

    A mutation perturbs the coordinates by adding a random number drawn from a normal
    distribution with mean 0 and standard deviation 0.1.
    """
    for atom_key, value in individual.items():
      if np.random.rand() < mutation_rate:
        # For single-spice case, value is assumed to be a coordinate list.
        if isinstance(value, list):
          individual[atom_key] = [coord + np.random.normal(0, 0.1) for coord in value]
        # For two-spice (bimetallic) case, value is assumed to be a dict with a "coords" key.
        elif isinstance(value, dict) and "coords" in value:
          new_coords = [coord + np.random.normal(0, 0.1) for coord in value["coords"]]
          # Preserve the spice information
          value["coords"] = new_coords
          individual[atom_key] = value
    return individual

