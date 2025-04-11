import os
import sys
import time
import json
import argparse
import numpy as np

import library.selection as selection
import library.population as pop
import library.runspresso as rspress
import library.readparams as parameters

start_time = time.time()

def load_parameters(file_path):
  with open(file_path, "r") as f:
    return json.load(f)

def main():
  parser = argparse.ArgumentParser(description="Run the Genetic Algorithm.")
  parser.add_argument("-paramsfile", type=str, required=True,
                      help="Path to the JSON parameter file")
  parser.add_argument("-singlespice", type=int, required=True,
                      help="singlespice: if you are running for a single atomic spice: 1, or more: 0")
  return parser.parse_args()

args = 0

if __name__ == "__main__":
  args = main()

# Reading parameters from JSON file
params = parameters.readparams()
params.fetchparams(args.paramsfile, args.singlespice)
#
# running the algorithm
#

# Loop over generations

clusters = {}
best_individuals = {}
history = {}

popu = pop.population()
runspresso = rspress.runspresso()
select = selection.selection()
total_energies = []

for gen in range(params.n_generations):
  print("\nRunning generation", gen, "\n")
  print("\nRunning generation", gen, "\n", file=sys.stderr)
  
  # ===============================
  # Step 1. Creating the population
  #
  print("\nMSD: on Step 1. Creating the population")

  if gen == 0:
    if params.getZeroPopu is None:
      # Creating the zero generation and writing .in files for espresso
      clusters = popu.write_espresso_file(params, args.singlespice)

    elif params.getZeroPopu is not None:
      clusters = popu.write_espresso_file(params)

  # ============================
  # Step 2. Espresso calculation
  #
    print("\nMSD: on Step 2. Espresso calculation")
    pids = []
    for key, individual in clusters.items():
      pid = runspresso.run_job_in_screen(job_id=key, nNodes=params.nNodes)
      pids.append(pid)
    # Wait for all PIDs to complete
    print("Waiting for screen")
    while any(os.path.exists(f"/proc/{pid}") for pid in pids):
      time.sleep(1)

    # Fetching the total energy for each individual
    tmpener = 0
    for key, individual in clusters.items():
      tmpener = runspresso.get_total_energies(key)
      individual['energ'] = tmpener
      total_energies.append(tmpener)

  else:
    # if not generation 0, work with the best of the previous generation 
    clusters = best_individuals

  # ===============================
  # Step 3. Calculating the Fitness
  #
  print("\nMSD: on Step 3. Calculating the Fitness")

  # =============================
  # Step 3.1. Fitness Calculation
  print("\nMSD: on Step 3.1. Fitness Calculation")
  V_min = np.min(total_energies)
  V_max = np.max(total_energies)
  for i, (key, individual) in enumerate(clusters.items()):
    print(i, key)
    clusters[key][f'fitness'] = select.calculate_fitness(total_energies[i], \
                                                         V_min=V_min, V_max=V_max, rho=params.rho_scaling)

  # =====================================================
  # Step 3.2. Selection by modified roulette wheel method
  print("\nMSD: on Step 3.2. Selection by modified roulette wheel method")
  selected_pairs = select.select_unique_pairs_for_mating(clusters, params.num_children)

  # =======================
  # Step 3.3. Mixing/Mating
  print("\nMSD: on Step 3.3. Mixing/Mating")
  child_clusters = {}
  for j, (p1_key, p2_key) in enumerate(selected_pairs, start=1):
    child = select.crossover(clusters[p1_key], clusters[p2_key], params.num_atoms)
    child_clusters[f'child{j}'] = child


  # ================
  # Step 4. Mutation
  #
  print("\nMSD: on Step 4. Mutation")
  for child_key in child_clusters.keys():
    mutated_child = select.mutate(child_clusters[child_key], params.mutation_rate)
    child_clusters[child_key] = mutated_child


  # ===========================
  # Step 5. Children evaluation
  #
  print("\nMSD: on Step 5. Children evaluation")

  # Running espresso
  popu.write_espresso_file_children(params, child_clusters, args.singlespice)

  pids = []
  for key, individual in child_clusters.items():
    pid = runspresso.run_job_in_screen(key, params.nNodes)
    pids.append(pid)
  # Wait for all PIDs to complete
  print("Waiting for screen for Children")
  while any(os.path.exists(f"/proc/{pid}") for pid in pids):
    time.sleep(1)

  # Fetching total energies for children
  tmpener = 0
  child_total_energies = []
  for key, individual in child_clusters.items():
    tmpener = runspresso.get_total_energies(key)
    individual['energ'] = tmpener
    child_total_energies.append(tmpener)


  # ===============
  # Step 6. Elitism
  #
  print("\nMSD: on Step 6. Elitism")
  combined_list = [(key, value) for key, value in clusters.items()] \
    + [(key, value) for key, value in child_clusters.items()]
  combined_list_sorted = sorted(combined_list, key=lambda x: x[1]['energ'])
  top_best_individuals = combined_list_sorted[:params.num_indvs]

  # Create/reset the best_individuals dictionary
  best_individuals = {}
  for i, (key, value) in enumerate(top_best_individuals):
    new_key = f'ind{i + 1}'
    # Copy the value dictionary and remove 'fitness' if it exists
    best_individuals[new_key] = {k: v for k, v in value.items() if k != 'fitness'}
    total_energies[i] = value['energ']
  
  # Writing the history
  history[f'record{gen}'] = [clusters, child_clusters]
  with open('history.json', 'w') as json_file:
    json.dump(history, json_file, indent=4)

print("\n\n") 
print("Writing final candidates into .in espresso's file")

popu.write_final_espresso_file(params, best_individuals)

end_time = time.time()

with open("execution_time.txt", "w") as file:
  file.write(f"Executaion time: {end_time - start_time}")
print("\nGood luck pichurria\n")
