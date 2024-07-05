import random

class population:

  def __init__(self):
    self.a = 0

  def create_content(self):
    content = """&control
calculation = 'scf',
prefix='Cu-1',!prefijo, debe cambiar segun el numero
outdir='/tmp/',
pseudo_dir = './'
/
&system
ibrav= 0,
nat=  4,!!numero de atomos
ntyp= 1,
ecutwfc = 50.0,
occupations='smearing',
smearing='marzari-vanderbilt',
degauss=0.04
/
&ELECTRONS
conv_thr = 1.0d-3
electron_maxstep = 10,
diagonalization  = 'david'
/
&IONS
ion_dynamics  = 'bfgs'
/
K_POINTS {gamma}
ATOMIC_SPECIES
Cu  63.546 Cu.pbe-dn-rrkjus_psl.1.0.0.UPF
CELL_PARAMETERS (angstrom)
        20.0000000000   0.0000000000   0.0000000000
         0.0000000000  20.0000000000   0.0000000000
         0.0000000000   0.0000000000  20.0000000000
ATOMIC_POSITIONS (angstrom)
  """
    return content

  def generate_random_coordinates(self, scale_factor=1.0, num_atoms=4):
    N_cubed_root = num_atoms ** (1/3)
    return f"{random.uniform(0, N_cubed_root) * scale_factor:.10f} {random.uniform(0, N_cubed_root) * scale_factor:.10f} {random.uniform(0, N_cubed_root) * scale_factor:.10f}"


  def write_espresso_file(self, num_indvs=1, num_atoms=4, r_scale=1.0):
    clusters = {}
    for indv in range(1, num_indvs+1):
      ind_key = f'ind{indv}'
      clusters[ind_key] = {}
      content = self.create_content()

      # Append random coordinates for each Cu atom
      for atom_num in range(1, num_atoms+1):
        tmp_coord = self.generate_random_coordinates(r_scale, num_atoms)
        string_components = tmp_coord.split()
        # Convert each component to a float
        clusters[ind_key][f'atom_coord{atom_num}'] = \
          [float(component) for component in string_components]
        content += f"Cu {tmp_coord}\n"

      # Write the content to the output file
      with open(f"{ind_key}.in", "w") as file:
        file.write(content)
    return clusters

  def write_espresso_file_children(self, child_clusters, num_indvs=1, num_atoms=4):  
    for indv in range(1, num_indvs+1):
      ind_key = f'child{indv}'
      content = self.create_content()
      # Append random coordinates for each Cu atom
      for atom_num in range(1, num_atoms+1):
        tmp_coord = child_clusters[ind_key][f'atom_coord{atom_num}']
        tmp_coord_str = " ".join(f"{coord:.10f}" for coord in tmp_coord)
        content += f"Cu {tmp_coord_str}\n"

      # Write the content to the output file
      with open(f"{ind_key}.in", "w") as file:
        file.write(content)

  def write_final_espresso_file(self, clusters, num_indvs=1, num_atoms=4):    
    for indv in range(1, num_indvs+1):
      ind_key = f'ind{indv}'
      content = self.create_content()

      # Append random coordinates for each Cu atom
      for atom_num in range(1, num_atoms+1):
        tmp_coord = clusters[ind_key][f'atom_coord{atom_num}']
        tmp_coord_str = " ".join(f"{coord:.10f}" for coord in tmp_coord)
        content += f"Cu {tmp_coord_str}\n"

      # Write the content to the output file
      with open(f"final_{ind_key}.in", "w") as file:
        file.write(content)
