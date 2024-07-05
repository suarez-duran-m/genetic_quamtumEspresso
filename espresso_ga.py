import time
import argparse
import numpy as np

import library.selection as selection
import library.population as pop
import library.runspresso as rspress

def main():
  parser = argparse.ArgumentParser(description="Run the Genetic Algorithm.")
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

  return parser.parse_args()

if __name__ == "__main__":
  args = main()

#
# running the algoritm
#

# Loop over generations

clusters = {}
next_generation = {}

popu = pop.population()
runspresso = rspress.runspresso()
select = selection.selection()

for gen in range(args.n_generations): 
  print("\nRunning generation", gen, "\n")
  
  # ===============================
  # Step 1. Creating the population  
  #
  if gen == 0:
    # Creating the zero generation and writing .in files for espresso
    clusters = popu.write_espresso_file(args.num_indvs, args.num_atoms, args.r_scale)
  else:
    clusters = next_generation

  # ============================
  # Step 2. Espresso calculation
  # Step 2. Espresso calculation  #

  # Running espresso 
  ischild = False
  #for i in range(1, args.num_indvs + 1):
  #  runspresso.run_job_in_screen(job_id=i, ischild=ischild)
  #  time.sleep(120)

  # Checking if all jobs finished
  runspresso.check_if_all_finished(args.num_indvs)

  # ===============================
  # Step 3. Calculating the Fitness
  #

  # Reading the energies from espresso outputs  
  total_energies = runspresso.get_total_energies(args.num_indvs, ischild=ischild)
  print("\nEnergies calculated")
  print(total_energies)

  # Adding total energies to clusters
  indv_ids = [f'ind{i}' for i in range(1, args.num_indvs + 1)]
  for i in range(len(indv_ids)):
    ind_key = f'ind{i+1}'
    clusters[ind_key][f'energ'] = total_energies[i]

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
  print("\nSelected pairs")
  print(selected_pairs)

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
  
  print("\nChildren clusters")
  print(child_clusters['child1'])

  # ===========================
  # Step 5. Children evaluation
  #

  # Running espresso
  ischild = True
  popu.write_espresso_file_children(child_clusters, args.num_children, args.num_atoms)

  #for i in range(1, args.num_children + 1):
  #  runspresso.run_job_in_screen(job_id=i, ischild=ischild)
  #  time.sleep(120)

  # Checking if all jobs finished
  runspresso.check_if_all_finished(args.num_children)

  # Fetching total energies for children
  child_total_energies = runspresso.get_total_energies(args.num_children, \
    ischild=ischild)
  print("\nEnergies calculated")
  print(child_total_energies, "\n")

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
  top_5_elements = combined_list_sorted[:5]

  # Create the next_generation dictionary
  next_generation = {}
  for i, (key, value) in enumerate(top_5_elements):
    new_key = f'ind{i + 1}'
    # Copy the value dictionary and remove 'fitness' if it exists
    next_generation[new_key] = {k: v for k, v in value.items() if k != 'fitness'}


print("\n\n") 
print("Writing final candidates into .in espresso's file")

popu.write_final_espresso_file(next_generation, args.num_indvs, args.num_atoms)

print("\nGood luck pichurria\n")
