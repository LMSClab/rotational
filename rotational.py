import numpy as np
import sys

def rotate_dihedral(atom_coords, dihedral_indices, angle_degrees):
    """
    Rotate the dihedral angle of a set of atoms and return the new coordinates.
    
    Args:
        atom_coords (list of numpy arrays): List of atom coordinates in the format [x, y, z].
        dihedral_indices (list of 4 integers): Indices of the four atoms defining the dihedral angle.
        angle_degrees (float): Angle in degrees by which to rotate the dihedral angle.
    
    Returns:
        list of numpy arrays: New atom coordinates after the dihedral angle rotation.
    """
    # Convert angle to radians
    angle_radians = np.radians(angle_degrees)
    
    # Extract atom coordinates and dihedral atoms
    a, b, c, d = [atom_coords[i] for i in dihedral_indices]
    
    # Calculate rotation axis and rotation matrix
    b0 = -1.0 * (b - a)
    b1 = c - b
    b2 = np.cross(b0, b1)
    b2 /= np.linalg.norm(b2)
    
    R = np.array([[np.cos(angle_radians), -np.sin(angle_radians), 0],
                  [np.sin(angle_radians), np.cos(angle_radians), 0],
                  [0, 0, 1]])
    
    # Rotate atoms
    rotated_coords = []
    for atom_coord in atom_coords:
        rotated = np.dot(R, atom_coord - b) + b
        rotated_coords.append(rotated)
    
    return rotated_coords

def rotate_dihedral_angle(xyz_file, atoms): #, start_angle, end_angle, step
  """
  Rotates the dihedral angle between the four atoms in the xyz file between the given range of angles with the given step increment.

  Args:
    xyz_file: The path to the xyz file.
    atoms: A list of the four atoms that you want to rotate the dihedral angle between.
    start_angle: The starting angle of the dihedral angle.
    end_angle: The ending angle of the dihedral angle.
    step: The step increment of the dihedral angle.

  Returns:
    A new xyz file with the rotated dihedral angle.
  """

  l = []

  # Lê o .gjf e retira o cabeçalho.
  with open(xyz_file, 'r') as f:
    read_lines = f.readlines()
    gjf_data = read_lines[8:]

    # Localiza a secção entre coordenadas e a matriz de ligações
    for key, n in enumerate(gjf_data):
       if n == '\n':
          l.append(key)

    # Organiza em listas as coordenadas e a matriz de ligações
    xyz_atoms = gjf_data[0 : l[0]]
    xyz_bonds = gjf_data[l[0] + 1 : l[1]]

    # Faz a organização dos índices recebidos pelo usuário
    atoms = atoms.split(',')
    atoms_int = []

    for n in atoms:
       n = int(n)
       atoms_int.append(n-1)

    # Separa as coordenadas referentes aos índices fornecidos
    select_atoms = []

    for m in atoms_int:
        for key, n in enumerate(xyz_atoms):
           if key == m:
              select_atoms.append(n)

    # Retira os espaços entre atomos e coordenadas e escreve em uma lista

    tratamento = []

    for n in select_atoms:
       b = n.split(' ')
       for m in b:
          if m != '':
             tratamento.append(m)

    # Indice para montagem do .gjf, remoção dos indices das coordenadas e substituição dos valores com /n

    indice = []

    for key, n in enumerate(tratamento):
       if key == 0 or key == 4 or key == 8 or key == 12:
          indice.append(n)
       if key == 0 or key == 3 or key == 6 or key == 9:
          tratamento.pop(key)
       if key == 2 or key == 5 or key == 8 or key == 11:
          a = n.replace('\n','')
          tratamento.pop(key)
          tratamento.insert(key, a)

    # Comprehension para converter str em float

    tratamento = [float(n) for n in tratamento]

    # Montagem dos arrays

    p1 = []
    p2 = []
    p3 = []
    p4 = []
    atom_coords = []

    for key, n in enumerate(tratamento):
       if key == 0 or key == 1 or key == 2:
          p1.append(n)
       if key == 3 or key == 4 or key == 5:
          p2.append(n)
       if key == 6 or key == 7 or key == 8:
          p3.append(n)
       if key == 9 or key == 10 or key == 11:
          p4.append(n)
    p1 = np.array(p1)
    atom_coords.append(p1)
    p2 = np.array(p2)
    atom_coords.append(p2)    
    p3 = np.array(p3)
    atom_coords.append(p3)    
    p4 = np.array(p4)      
    atom_coords.append(p4)

    dihedral_indices = [0, 1, 2, 3]

    angle_degrees = 45.0

    new_coords = rotate_dihedral(atom_coords, dihedral_indices, angle_degrees)

    print(new_coords)


  '''
  # Identify the four atoms that you want to rotate the dihedral angle between.
  atom_indices = [xyz_data[i].split()[0] for i in atoms]

  # Calculate the dihedral angle between the four atoms.
  dihedral_angle = xyz.calculate_dihedral_angle(xyz_data, atom_indices)

  # Create a loop that iterates over the range of angles.
  for angle in np.arange(start_angle, end_angle, step):
    # Rotate the dihedral angle by the step increment.
    dihedral_angle = dihedral_angle + step

    # Calculate the new dihedral angle.
    new_dihedral_angle = xyz.calculate_dihedral_angle(xyz_data, atom_indices, dihedral_angle)

    # Write the new xyz file with the rotated dihedral angle.
    with open(xyz_file, 'w') as f:
      f.writelines(xyz_data)
    '''
  return xyz_file

if __name__ == '__main__':
   if len(sys.argv) < 2:
        print("Usage: python rotational.py <xyz_file> <list_atoms> <start_angles> <end_angle> <step>")
   else:
        #XYZ - Deve ser gerado com o GEOMCHECK
        xyz_file = sys.argv[1]
        atoms = sys.argv[2]
        #start_angle = sys.argv[3]
        #end_angle = sys.argv[4]
        #step = sys.argv[5]
        # Rotate the dihedral angle and write the new xyz file.
        new_xyz_file = rotate_dihedral_angle(xyz_file, atoms) #, start_angle, end_angle, step