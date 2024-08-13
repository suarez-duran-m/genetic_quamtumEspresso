import random

class population:

  def __init__(self):
    self.a = 0

  def create_content(self, prefix='nan', elementName='nan', natom=1, \
          totalatoms=1, atomicweight=1,  pseudodir='nan', \
          pseudopotential='nan'):
    content = """&control
calculation = 'scf',
prefix='PREFIXNAME', !prefijo, debe cambiar segun el numero
outdir='/tmp/',
pseudo_dir = 'PATHPSEUDO'
/
&system
ibrav= 0,
nat= totATOM,!!numero de atomos
ntyp= 1,
ecutwfc = 60.0,
occupations='smearing',
smearing='marzari-vanderbilt',
degauss=0.04
/
&ELECTRONS
conv_thr = 1.0d-6
electron_maxstep = 100,
diagonalization  = 'david'
/
K_POINTS {gamma}
ATOMIC_SPECIES
ELENAME ATOWEIGHT PSEUDOPOT
CELL_PARAMETERS (angstrom)
        20.0000000000   0.0000000000   0.0000000000
         0.0000000000  20.0000000000   0.0000000000
         0.0000000000   0.0000000000  20.0000000000
ATOMIC_POSITIONS (angstrom)
  """
    
    content = content.replace('PREFIXNAME', f'{prefix}')
    content = content.replace('ELENAME', f'{elementName}')
    content = content.replace('-nATOM', f'-{natom}')
    content = content.replace('totATOM', f'{totalatoms}')
    content = content.replace('ATOWEIGHT', f'{atomicweight}')
    content = content.replace('PATHPSEUDO', f'{pseudodir}')
    content = content.replace('PSEUDOPOT', f'{pseudopotential}')

    return content

  def generate_random_coordinates(self, scale_factor=1.0, num_atoms=4):
    N_cubed_root = num_atoms ** (1/3)
    return f"{random.uniform(0, N_cubed_root) * scale_factor:.10f} {random.uniform(0, N_cubed_root) * scale_factor:.10f} {random.uniform(0, N_cubed_root) * scale_factor:.10f}"


  def write_espresso_file(self, prefix='nan', num_indvs=1, num_atoms=4, \
          r_scale=1.0, atomName='nan', atomic_weight=1.0, \
          pseudoDir='nan', pseudo_potential='nan'):
    clusters = {}
    for indv in range(1, num_indvs+1):
      ind_key = f'ind{indv}'
      clusters[ind_key] = {}
      content = self.create_content(prefix, elementName=atomName, \
              natom=indv, totalatoms=num_atoms, atomicweight=atomic_weight, \
              pseudodir=pseudoDir, pseudopotential=pseudo_potential)

      # Append random coordinates for each Cu atom
      for atom_num in range(1, num_atoms+1):
        tmp_coord = self.generate_random_coordinates(r_scale, num_atoms)
        string_components = tmp_coord.split()
        # Convert each component to a float
        clusters[ind_key][f'atom_coord{atom_num}'] = \
          [float(component) for component in string_components]
        content += f"{atomName} {tmp_coord}\n"

      # Write the content to the output file
      with open(f"{ind_key}.in", "w") as file:
        file.write(content)
    return clusters

  def write_espresso_file_children(self, prefix='nan', child_clusters={}, \
          num_indvs=1, num_atoms=4, atomName='nan', atomic_weight=1.0, \
          pseudoDir='nan', pseudo_potential='nan'):
    for indv in range(1, num_indvs+1):
      ind_key = f'child{indv}'
      content = self.create_content(prefix, elementName=atomName, \
              natom=indv, totalatoms=num_atoms, \
              atomicweight=atomic_weight, pseudodir=pseudoDir, \
              pseudopotential=pseudo_potential)

      # Append random coordinates for each atom
      for atom_num in range(1, num_atoms+1):
        tmp_coord = child_clusters[ind_key][f'atom_coord{atom_num}']
        tmp_coord_str = " ".join(f"{coord:.10f}" for coord in tmp_coord)
        content += f"{atomName} {tmp_coord_str}\n"

      # Write the content to the output file
      with open(f"{ind_key}.in", "w") as file:
        file.write(content)

  def write_final_espresso_file(self, prefix='nan', last_cluster={} \
          , num_indvs=1, num_atoms=4, atomName='nan', atomic_weight=1.0, \
          pseudoDir='nan', pseudo_potential='nan'):

    for indv in range(1, num_indvs+1):
      ind_key = f'ind{indv}'
      content = self.create_content(prefix, elementName=atomName, \
              natom=indv, totalatoms=num_atoms, \
              atomicweight=atomic_weight, pseudodir=pseudoDir, \
              pseudopotential=pseudo_potential)

      # Append random coordinates for each atom
      for atom_num in range(1, num_atoms+1):
        tmp_coord = last_cluster[ind_key][f'atom_coord{atom_num}']
        tmp_coord_str = " ".join(f"{coord:.10f}" for coord in tmp_coord)
        content += f"{atomName} {tmp_coord_str}\n"

      # Write the content to the output file
      with open(f"final_{ind_key}.in", "w") as file:
        file.write(content)
