import sys
import time
import json
import argparse
import numpy as np

import library.selection as selection
import library.population as pop
import library.runspresso as rspress


start_time = time.time()

def main():
  parser = argparse.ArgumentParser(description="Run the Genetic Algorithm.")
  parser.add_argument("-prefix", type=str, required=True,
    help="Prefix name to identify the Espresso's data card")
  parser.add_argument("-getZeroPopu", type=str, required=False,
    help="Argument with filename from zero population will be create")
  parser.add_argument("-num_indvs", type=int, required=True,
    help="Number of individuals for the cluster")
  parser.add_argument("-num_atoms", type=int, required=True, 
    help="Number of atoms per individual/cluster")
  parser.add_argument("-num_children", type=int, required=True, 
    help="Number of children from best individuals/clusters")
  parser.add_argument("-mutation_rate", type=float, required=True,
    help="Probability of mutation")
  parser.add_argument("-n_generations", type=int, required=True, 
    help="Number of generations to run")
  parser.add_argument("-r_scale", type=float, required=True,
    help="Scale factor for atoms distances inside the cluster")
  parser.add_argument("-rho_scaling", type=float, required=True,
    help="Scale factor for atoms distances inside the cluster")
  parser.add_argument("-ele_name", type=str, required=True,
    help="Atom Element")
  parser.add_argument("-atom_weight", type=float, required=True,
    help="Atomic Weight")
  parser.add_argument("-pseudoDir", type=str, required=True,
    help="Path where pseudo potential is located")
  parser.add_argument("-pseudo", type=str, required=True,
    help="Pseudo potential")
  parser.add_argument("-nNodes", type=int, required=True,
    help="Number of Nodes to use")

  return parser.parse_args()

if __name__ == "__main__":
  args = main()

#
# running the algoritm
#

# Loop over generations

clusters = {}
best_individuals = {}
history = {}

popu = pop.population()
runspresso = rspress.runspresso()
select = selection.selection()
total_energies = []

for gen in range(args.n_generations): 
  print("\nRunning generation", gen, "\n")
  print("\nRunning generation", gen, "\n", file=sys.stderr)
  
  # ===============================
  # Step 1. Creating the population  
  #
  if gen == 0:
    if args.getZeroPopu is None:
      # Creating the zero generation and writing .in files for espresso
      clusters = popu.write_espresso_file(args.prefix, args.num_indvs, \
              args.num_atoms, args.r_scale, args.ele_name, \
              args.atom_weight, args.pseudoDir, args.pseudo)
    elif args.getZeroPopu is not None:
      clusters = popu.write_espresso_file_notRandom(args.prefix, args.num_indvs, \
              args.num_atoms, args.ele_name, args.atom_weight, \
              args.pseudoDir, args.pseudo, args.getZeroPopu)

  # ============================
  # Step 2. Espresso calculation    
  #
    ischild = False
    for i in range(1, args.num_indvs + 1):
      runspresso.run_job_in_screen(job_id=i, ischild=ischild, nNodes=args.nNodes)

    # Reading the energies from espresso outputs
    total_energies = runspresso.get_total_energies(args.num_indvs, ischild=ischild)
    
    # Adding total energies to clusters
    indv_ids = [f'ind{i}' for i in range(1, args.num_indvs + 1)]
    for i in range(len(indv_ids)):
      ind_key = f'ind{i+1}'
      clusters[ind_key][f'energ'] = total_energies[i]

  else:
    # if not generation 0, work with the best of the previous generation 
    clusters = best_individuals

  # ===============================
  # Step 3. Calculating the Fitness
  #

  # =============================
  # Step 3.1. Fitness Calculation

  V_min = np.min(total_energies)
  V_max = np.max(total_energies)
  for i in range(args.num_indvs):
    ind_key = f'ind{i+1}'
    clusters[ind_key][f'fitness'] = select.calculate_fitness(total_energies[i], \
      V_min=V_min, V_max=V_max, rho=args.rho_scaling)

  # =====================================================
  # Step 3.2. Selection by modified roulette wheel method
  selected_pairs = select.select_unique_pairs_for_mating(clusters, args.num_children)

  # =======================
  # Step 3.3. Mixing/Mating
  child_clusters = {}

  for j in range(len(selected_pairs)):
    parent1 = clusters[selected_pairs[j][0]]
    parent2 = clusters[selected_pairs[j][1]]
    child_key = f'child{j + 1}'
    child = select.crossover(parent1, parent2, args.num_atoms)
    child_clusters[child_key] = child


  # ================
  # Step 4. Mutation
  #
  for child_key in child_clusters.keys():
    mutated_child = select.mutate(child_clusters[child_key], args.mutation_rate)
    child_clusters[child_key] = mutated_child


  # ===========================
  # Step 5. Children evaluation
  #

  # Running espresso
  ischild = True
  popu.write_espresso_file_children(args.prefix, child_clusters, \
          args.num_children, args.num_atoms, args.ele_name, \
          args.atom_weight, args.pseudoDir, args.pseudo)
 
  for i in range(1, args.num_children + 1):
    runspresso.run_job_in_screen(job_id=i, ischild=ischild, nNodes=args.nNodes)

  # Fetching total energies for children
  child_total_energies = runspresso.get_total_energies(args.num_children, \
    ischild=ischild)

  # Adding total energies to children clusters
  indv_ids = [f'child{i}' for i in range(1, args.num_children + 1)]
  for i in range(len(indv_ids)):
    ind_key = f'child{i+1}'
    child_clusters[ind_key][f'energ'] = child_total_energies[i]


  # ===============
  # Step 6. Elitism
  #
  combined_list = [(key, value) for key, value in clusters.items()] \
    + [(key, value) for key, value in child_clusters.items()]
  combined_list_sorted = sorted(combined_list, key=lambda x: x[1]['energ'])
  top_best_individuals = combined_list_sorted[:args.num_indvs]

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

popu.write_final_espresso_file(args.prefix, best_individuals, \
        args.num_indvs, args.num_atoms, args.ele_name, \
        args.atom_weight, args.pseudoDir, args.pseudo)

end_time = time.time()

with open("execution_time.txt", "w") as file:
  file.write(f"Executaion time: {end_time - start_time}")
print("\nGood luck pichurria\n")
