import random
import json
import math

class population:

  def __init__(self):
    self.a = 0

  def create_content(self, prefix='nan', elementName='nan', totalatoms=1, \
                     atomicweight=1,  pseudodir='nan', pseudopotential='nan'):
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
    content = content.replace('totATOM', f'{totalatoms}')
    content = content.replace('ATOWEIGHT', f'{atomicweight}')
    content = content.replace('PATHPSEUDO', f'{pseudodir}')
    content = content.replace('PSEUDOPOT', f'{pseudopotential}')

    return content

  def create_content_4twospices(self, prefix='nan', elementNameA='nan', elementNameB='nan', \
                            nspices=2, totalatoms=2, atomicweightA=1, atomicweightB=1, \
                            pseudodir='nan', pseudopotentialA='nan', pseudopotentialB='nan'):
    content = """&control
calculation = 'relax',
prefix='PREFIXNAME'
outdir='/tmp', 
pseudo_dir = 'PATHPSEUDO', 
/
&system    
ibrav= 0, 
nat= totATOM,
ntyp= nSPICES,
ecutwfc = 60.0,    
occupations='smearing', 
smearing='marzari-vanderbilt', 
degauss=0.04
/
&ELECTRONS
electron_maxstep = 100,
conv_thr = 1.0d-3
mixing_mode  = 'plain'
mixing_beta  = 0.3D0
diagonalization  = 'david'
startingwfc  = 'atomic+random' 
/
&IONS
trust_radius_max = 0.1D0
ion_dynamics      = 'bfgs'
ion_positions     = 'default'
pot_extrapolation = 'atomic'
wfc_extrapolation = 'none'
/
&CELL
cell_dynamics = 'bfgs',
cell_dofree = 'all',
/
CELL_PARAMETERS (angstrom)
       20.0000000000       0.0000000000        0.0000000000
        0.0000000000      20.0000000000        0.0000000000
        0.0000000000       0.0000000000       20.0000000000 
ATOMIC_SPECIES
ELENAMEA  ATOWEIGHTA  PSEUDOPOTA 
ELENAMEB  ATOWEIGHTB  PSEUDOPOTB

K_POINTS {gamma}

ATOMIC_POSITIONS (angstrom)
    """
    # Fitting the based datacard for the simulation
    content = content.replace('PREFIXNAME', f'{prefix}')
    content = content.replace('ELENAMEA', f'{elementNameA}')
    content = content.replace('ELENAMEB', f'{elementNameB}')
    content = content.replace('nSPICES', f'{nspices}')
    content = content.replace('totATOM', f'{totalatoms}')
    content = content.replace('ATOWEIGHTA', f'{atomicweightA}')
    content = content.replace('ATOWEIGHTB', f'{atomicweightB}')
    content = content.replace('PATHPSEUDO', f'{pseudodir}')
    content = content.replace('PSEUDOPOTA', f'{pseudopotentialA}')
    content = content.replace('PSEUDOPOTB', f'{pseudopotentialB}')
    return content

  def create_content_magnetic(self, prefix='nan', elementName='nan', totalatoms=1, \
                              atomicweight=1, pseudodir='nan', pseudopotential='nan'):
    content = """&control
  calculation = 'scf',
  prefix='PREFIXNAME', !prefijo, debe cambiar segun el numero
  outdir='/tmp/',
  pseudo_dir = 'PATHPSEUDO'
  /
  &system
  ibrav= 0,
  nat= totATOM,!!numero de atomos
  ntyp= 2,
  ecutwfc = 60.0,
  occupations='smearing',
  smearing='marzari-vanderbilt',
  degauss=0.04
  nsping=2
  starting_magnetization(1)= -0.1,
  starting_magnetization(2)= -0.1,
  /
  &ELECTRONS
  conv_thr = 1.0d-6
  electron_maxstep = 100,
  diagonalization  = 'david'
  /
  K_POINTS {gamma}
  ATOMIC_SPECIES
  ELENAME1 ATOWEIGHT PSEUDOPOT
  ELENAME2 ATOWEIGHT PSEUDOPOT  
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

  def generate_random_coordinates(self, scale_factor=1.0, num_atoms=4, min_distance=0.):
    def distance(coord1, coord2):
      """Calculate the Euclidean distance between two 3D coordinates."""
      return math.sqrt(
        (coord1[0] - coord2[0]) ** 2 +
        (coord1[1] - coord2[1]) ** 2 +
        (coord1[2] - coord2[2]) ** 2
      )

    N_cubed_root = num_atoms ** (1 / 3)  # Calculate cube root of the number of atoms
    coordinates = []  # Store generated coordinates here
    cube_side = N_cubed_root * scale_factor  # Maximum value for x, y, z

    for i in range(num_atoms):
      # Generate random coordinates for the new atom
      x = random.uniform(0, N_cubed_root) * scale_factor
      y = random.uniform(0, N_cubed_root) * scale_factor
      z = random.uniform(0, N_cubed_root) * scale_factor
      new_atom = (x, y, z)

      # If this is the first atom, checking for distances is not needed
      if not coordinates:
        coordinates.append(new_atom)

      while True:
        too_close = False

        # Check the distance from the new atom to all previously placed atoms
        for atom in coordinates:
          if distance(new_atom, atom) < min_distance:
            too_close = True
            new_atom = ((new_atom[0] + random.uniform(0, 0.2)) % cube_side,
                    (new_atom[1] + random.uniform(0, 0.2)) % cube_side,
                    (new_atom[2] + random.uniform(0, 0.2)) % cube_side )
          break

        # If no atoms are too close, accept this atom's position
        if not too_close:
          coordinates.append(new_atom)
          break  # Exit the while loop and proceed with the next atom

    # Return the coordinates formatted as strings
    return coordinates

  def write_espresso_file(self, params='nan', issingle=1):
    return self.write_espresso_file_4single(params) if issingle else self.write_espresso_file_4twospices(params)

  def write_espresso_file_4single(self, params='nan'):
    clusters = {}

    for indv in range(1, params.num_indvs+1):
      ind_key = f'ind{indv}'
      clusters[ind_key] = {}
      if params.ifMagnetic == 0:
        content = self.create_content(params.prefix, params.ele_name, params.num_atoms, \
                                      params.atom_weight, params.pseudoDir, params.pseudo)

      elif params.ifMagnetic == 1:
        content = self.create_content_magnetic(params.prefix, params.ele_name, \
                                               indv, params.num_atoms, params.atom_weight, \
                                               params.pseudoDir, params.pseudo)

      # Append random coordinates for each atom
      coordinates = self.generate_random_coordinates(params.r_scale, params.num_atoms, 2.0)
      for atom_num in range(1, params.num_atoms+1):
        tmp_coord = f"{coordinates[atom_num-1][0]:.10f} {coordinates[atom_num-1][1]:.10f} {coordinates[atom_num-1][2]:.10f}"
        string_components = tmp_coord.split()
        # Convert each component to a float
        clusters[ind_key][f'atom_coord{atom_num}'] = \
          [float(component) for component in string_components]

        if params.ifMagnetic == 0:
          content += f"{params.ele_name} {tmp_coord}\n"
        elif params.ifMagnetic == 1:
          if atom_num % 2 == 0:
            content += f"{params.ele_name}2 {tmp_coord}\n"
          else:
            content += f"{params.ele_name}1 {tmp_coord}\n"

      # Write the content to the output file
      with open(f"{ind_key}.in", "w") as file:
        file.write(content)
    return clusters

  def write_espresso_file_4twospices(self, params='nan'):
    clusters = {}
    nspices = 2
    totalatoms = params.num_atomsA + params.num_atomsB

    for indv in range(1, params.num_indvs+1):
      ind_key = f'ind{indv}'
      clusters[ind_key] = {}
      content = self.create_content_4twospices(params.prefix, params.ele_nameA, params.ele_nameB, nspices, \
                                               totalatoms, params.atom_weightA, params.atom_weightB, \
                                               params.pseudoDir, params.pseudoA, params.pseudoB)
      # Append random coordinates for each atom
      coordinates = self.generate_random_coordinates(params.r_scale, totalatoms, min_distance=2.0)
      for atom_num in range(1, totalatoms+1):
        # Format the coordinate string
        tmp_coord = f"{coordinates[atom_num-1][0]:.10f} {coordinates[atom_num-1][1]:.10f} {coordinates[atom_num-1][2]:.10f}"
        # Decide on species based on the atom number
        if atom_num <= params.num_atomsA:
          species_value = params.ele_nameA
          content += f"{params.ele_nameA} {tmp_coord}\n"
        else:
          species_value = params.ele_nameB
          content += f"{params.ele_nameB} {tmp_coord}\n"
        # Store each atom as a dictionary containing both coordinates and species
        clusters[ind_key] [f'atom_{atom_num}'] = {
          'coords': [float(component) for component in tmp_coord.split()],
          'spice': species_value
        }
      # Write the content to the output file
      with open(f"{ind_key}.in", "w") as file:
        file.write(content)
    return clusters

  def write_espresso_file_notRandom(self, params='nan'):
    clusters = {}
    with open(params.getZeroPopu) as file:
      data = json.load(file)
    last_record_key = sorted(data.keys(), key=lambda x: int(x.replace('record', '')))[-1]
    last_record_value = data[last_record_key]

    for indv in range(1, params.num_indvs+1):
      ind_key = f'ind{indv}'
      clusters[ind_key] = {}
      content = self.create_content(params.prefix, elementName=params.ele_name, \
              natom=indv, totalatoms=params.num_atoms, atomicweight=params.atom_weight, \
              pseudodir=params.pseudoDir, pseudopotential=params.pseudo)

      for atom_num in range(1, params.num_atoms + 1):
        tmp_coord = last_record_value[0][ind_key][f'atom_coord{atom_num}']
        clusters[ind_key][f'atom_coord{atom_num}'] = [float(component) for component in tmp_coord]
        content += f"{params.ele_name} {' '.join(map(str, tmp_coord))}\n"
      # Write the content to the output file
      with open(f"{ind_key}.in", "w") as file:
        file.write(content)

    return clusters

  def write_espresso_file_children(self, params='nan', child_clusters='nan', issingle=1):
    for ind_key, individual in child_clusters.items():
      if issingle:
        content = self.create_content(params.prefix, params.ele_name, params.num_atoms, \
                                    params.atom_weight, params.pseudoDir, params.pseudo)
        # Append random coordinates for each atom
        for atom_num in range(1, params.num_atoms + 1):
          tmp_coord = individual[f'atom_coord{atom_num}']
          tmp_coord_str = " ".join(f"{coord:.10f}" for coord in tmp_coord)
          content += f"{params.ele_name} {tmp_coord_str}\n"
      else:
        nspices = 2
        totalatoms = params.num_atomsA + params.num_atomsB
        content = self.create_content_4twospices(params.prefix, params.ele_nameA, params.ele_nameB, nspices, \
                                                 totalatoms, params.atom_weightA, params.atom_weightB, \
                                                 params.pseudoDir, params.pseudoA, params.pseudoB)
        # Append random coordinates for each atom
        for atom_num in range(1, totalatoms+1):
          tmp_coord = individual[f'atom_{atom_num}']['coords']
          tmp_coord_str = " ".join(f"{coord:.10f}" for coord in tmp_coord)
          spice = individual[f'atom_{atom_num}']['spice']
          content += f"{spice} {tmp_coord_str}\n"

      # Write the content to the output file
      with open(f"{ind_key}.in", "w") as file:
        file.write(content)

  def write_final_espresso_file(self, params='nan', last_cluster={}, issingle=1):
    for ind_key, individual in last_cluster.items():
      if issingle:
        content = self.create_content(params.prefix, params.ele_name, params.num_atoms, \
                                    params.atom_weight, params.pseudoDir, params.pseudo)
        # Append random coordinates for each atom
        for atom_num in range(1, params.num_atoms+1):
          tmp_coord = individual[f'atom_coord{atom_num}']
          tmp_coord_str = " ".join(f"{coord:.10f}" for coord in tmp_coord)
          content += f"{params.ele_name} {tmp_coord_str}\n"
      else:
        nspices = 2
        totalatoms = params.num_atomsA + params.num_atomsB
        content = self.create_content_4twospices(params.prefix, params.ele_nameA, params.ele_nameB, nspices, \
                                                 totalatoms, params.atom_weightA, params.atom_weightB, \
                                                 params.pseudoDir, params.pseudoA, params.pseudoB)
        # Append random coordinates for each atom
        for atom_num in range(1, totalatoms + 1):
          tmp_coord = individual[f'atom_{atom_num}']['coords']
          tmp_coord_str = " ".join(f"{coord:.10f}" for coord in tmp_coord)
          spice = individual[f'atom_{atom_num}']['spice']
          content += f"{spice} {tmp_coord_str}\n"

      # Write the content to the output file
      with open(f"final_{ind_key}.in", "w") as file:
        file.write(content)

  def getTotalAtoms4single(self, params='nan'):
    return params.num_atoms
  def getTotalAtoms4two(self, params='nan'):
    return params.num_atomsA + params.num_atomsB