"""
Author: Shenheng Yan, Hao Ma
Last update: 2023-01-06
Mail: oaham@xmu.edu.cn
Usage: Calculating electric field
"""
import numpy as np
#For calculating Electric field. A bug here, what sparc used now is a pdbfile without hydrogen,
# however we need hydrogen to calculate electric field. So it may be performed separately.
def extract_keyword(keyword, file_name):
    # extract the word following sprcified keyword in file_name
    with open(file_name, "r") as f:
        for line in f.readlines():
            if (keyword in line):
                return (line.split()[1])

def extract_keyword3(keyword, file_name):
    # extract the three words following sprcified keyword in file_name
    with open(file_name, "r") as f:
        for line in f.readlines():
            if (keyword in line):
                return (line.split()[1:4])

def extract_keyword_line_num(keyword, file_name):
    # extract the line number (starting from 0)where a specified keyword in file_name occurs
    with open(file_name, "r") as f:
        line_number = 0
        for line in f.readlines():
            if (keyword in line):
                return (line_number)
            line_number = line_number + 1

def read_input(input_info, input_file_name):
    # read from input file input informations, and store in the dictionary input_info
    # For compiling this algorthm into sparc, it is not useful,
    # but we want to retain its features and can be invoked separately.

    input_info["atom_num"] = extract_keyword("atom_number", input_file_name)
    input_info["prmtop_name"] = extract_keyword("prmtop_filename", input_file_name)
    input_info["pdb_name"] = extract_keyword("pdb_filename", input_file_name)
    input_info["charge_deletion_pdb_name"] = extract_keyword("charge_deletion_file", input_file_name)
    input_info["center_coord"] = extract_keyword3("center_atom_coord", input_file_name)
    input_info["atom1_coord"] = extract_keyword3("atom1_coord", input_file_name)
    input_info["atom2_coord"] = extract_keyword3("atom2_coord", input_file_name)
    return ()

def read_pdb(pdb_name, coords, residue_ids):
    # read from the pdb file coordinates and residue numbers of atoms
    with open(pdb_name) as f:
        current_atom_idx = 0
        for line in f.readlines():
            if (line.split()[0] == "ATOM"):
                residue_ids[current_atom_idx] = line.split()[4]
                coords[current_atom_idx, :] = np.array(line.split()[5:8])
                current_atom_idx = current_atom_idx + 1
    return current_atom_idx

def read_prmtop(prmtop_name, charges, key):
    # read from the prmtop file atomic charges
    if key == "amber":
        atom_num = len(charges)
        # line number of the lines which starts to record charge in prmtop file (FORMAT(5E16.8))
        starting_line_num = extract_keyword_line_num("%FLAG CHARGE", prmtop_name) + 2
        # number of lines recording charge in prmtop file (FORMAT(5E16.8))
        complete_line_num = atom_num / 5

        with open(prmtop_name, "r") as prmtop_file:
            all_lines = prmtop_file.readlines()
            current_atom = 0

            for current_line_num in range(starting_line_num, starting_line_num + complete_line_num + 1):
                charges[current_atom: current_atom + 5] = all_lines[current_line_num].split()
                current_atom = current_atom + 5

    # charges are multiplied by 18.2223 in prmtop files, need to transform the units back to e.
    # it depends on which software we used. For amber, charge are multiplied with 18.2223. whereas gmx, charge are charge with units of e.
    # Amber version is recommanded, because it is more precise. (0.0000)  GMX: 0.00
            charges[:] = charges[:] / 18.2223

    if key == "gmx":
        with open(prmtop_name) as f:
            current_atom_idx = 0
            for line in f.readlines():
                if (line.split()[0] == "ATOM"):
                    charges[current_atom_idx] = line.split()[8]
                    current_atom_idx = current_atom_idx + 1

    return ()
# reserve this functions, another method is to accumulate the charges. To read another pdb file(third) is too tedious!
#  We know which residue that needs to be setted to zero.
def modify_charges(pdb_name, residue_num, charges):
    atom_num1 = len(charges)
    # delete charges of those atoms in charge_deletion_pdb
    with open(pdb_name, "r") as charge_deletion_pdb:
        for line in charge_deletion_pdb.readlines():
            if ("ATOM" in line):
                if line.split()[4]== residue_num:
                    charges[atom_num1 - 1] = 0
    return ()

def get_directional_vectors(coords, r0):
    # get directional vectors and distances from atoms to r0
    atom_num = np.shape(coords)[0]
    directional_vectors = np.zeros([atom_num, 3], dtype="float")
    for i in range(0, atom_num):
        directional_vectors[i, :] = coords[i, :] - r0[:]
    return (directional_vectors)

def get_distances(directional_vectors):
    # get distances from atoms to r0
    atom_num = np.shape(directional_vectors)[0]
    distances = np.zeros([atom_num], dtype="float")
    for i in range(0, atom_num):
        distances[i] = (np.dot(directional_vectors[i, :], directional_vectors[i, :])) ** 0.5
    return (distances)

def get_unit_directional_vectors(directional_vectors, distances):
    # get UNIT directional vectors and distances from atoms to r0
    atom_num = np.shape(directional_vectors)[0]
    unit_directional_vectors = np.zeros([atom_num, 3], dtype="float")
    for i in range(0, atom_num):
        unit_directional_vectors[i, :] = directional_vectors[i, :] / distances[i]
    return (unit_directional_vectors)

def get_atoms_electric_field_modulus(distances, charges, r0):
    # get individual atom's e field modulus according to |E| = kQ/r^2
    atom_num = np.shape(distances)[0]
    k = 1439.9757
    atoms_e_field_modulus = np.zeros([atom_num], dtype="float")
    for i in range(0, atom_num):
        if (distances[i] == 0):
            atoms_e_field_modulus[i] = 0
        atoms_e_field_modulus[i] = -1 * k * charges[i] / (distances[i] ** 2)
    return atoms_e_field_modulus

def get_atoms_e_field_vectors(atoms_e_field_modulus, unit_directional_vectors):
    # get individual atom's e contribution to field VECTORS by multiply vectors' modulus by directional vectors
    atom_num = np.shape(atoms_e_field_modulus)[0]
    atoms_e_field_vectors = np.zeros([atom_num, 3], dtype="float")
    for i in range(0, atom_num):
        atoms_e_field_vectors[i, :] = atoms_e_field_modulus[i] * unit_directional_vectors[i, :]
    return (atoms_e_field_vectors)

def get_atoms_projected_e_field_vectors(atoms_e_field_vectors, chosen_axis):
    # project individual atom's contribution to e field to the chosen axis
    atom_num = np.shape(atoms_e_field_vectors)[0]
    atoms_projected_e_field_vectors = np.zeros([atom_num], dtype="float")
    unit_vector_along_chosen_axis = chosen_axis / (np.dot(chosen_axis, chosen_axis)) ** 0.5
    for i in range(0, atom_num):
        atoms_projected_e_field_vectors[i] = np.dot(atoms_e_field_vectors[i], unit_vector_along_chosen_axis)
    return (atoms_projected_e_field_vectors)

def get_residues_projected_e_field_vectors(atoms_projected_e_field_vectors, residue_ids):
    # sum individual atoms' contribution to PROJECTED e field to the chosen axis within the same residue to get
    # residues' contribution to PROJECTED e field
    atom_num = np.shape(residue_ids)[0]
    residue_num = residue_ids[-1]
    residues_projected_e_field_vectors = np.zeros([residue_num], dtype="float")
    for i in range(0, atom_num):
        current_residue = residue_ids[i] - 1  # residue_ids[:] starts with 1, not 0
        residues_projected_e_field_vectors[current_residue] = residues_projected_e_field_vectors[current_residue] \
                                                              + atoms_projected_e_field_vectors[i]

    return (residues_projected_e_field_vectors)

# def write_vis_pdb(total_e_field_vector, r0):
#     # write a pdb file for field direction visualization
#     # this pdb contains two dummy atoms, one at the position of r0, the other is 5A away from r0, and the line connecting r0
#     # is along the direction of the total field
#     dummy_atom_coord = r0 + 5 * total_e_field_vector / np.dot(total_e_field_vector, total_e_field_vector)
#     return ()

def calculate_electric_field(filepath, charge_file, atom_num ,r_num, atom1, atom2, r_init):
    # read from the input file informations about atom numbers, pdb name, etc. and store in dictionary input_info
    # input_file_name = "efield_quan.inp"
    # input_info = dict(atom_num="", prmtop_name="", pdb_name="", charge_deletion_pdb_name="", \
    #                   center_coord="", atom1_coord="", atom2_coord="")
    # read_input(input_info, input_file_name)
    # read from the pdb file coordinates and residue numbers of atoms
    atom_num = atom_num
    coords = np.zeros([atom_num, 3], dtype="float")
    residue_ids = np.zeros([atom_num], dtype="int")

    atom_ = read_pdb(filepath, coords, residue_ids)
    assert atom_ == atom_num
    # read from the prmtop file charges of atoms, key is gmx.
    charges = np.zeros([atom_num], dtype="float")
    read_prmtop(charge_file, charges, "gmx")

    # modify charges of those atoms in charge_deletion_pdb to be zero
    modify_charges(filepath, r_num, charges)

    # r0 is the point where we calculate the e field
    r0 = np.array(r_init, dtype="float")

    # get directional vectors from atoms to r0
    directional_vectors = get_directional_vectors(coords, r0)

    # get distances from atoms to r0
    distances = get_distances(directional_vectors)

    # get UNIT directional vectors from atoms to r0
    unit_directional_vectors = get_unit_directional_vectors(directional_vectors, distances)

    # get individual atom's contribution to e field modulus according to |E| = kQ/r^2
    atoms_e_field_modulus = get_atoms_electric_field_modulus(distances, charges, r0)

    # get individual atom's e contribution to field VECTORS by multiply vectors' modulus by directional vectors
    atoms_e_field_vectors = get_atoms_e_field_vectors(atoms_e_field_modulus, unit_directional_vectors)

    # sum atoms' field to get the total field
    total_e_field_vector = np.sum(atoms_e_field_vectors, axis=0)
    modulus_total_e_field_vector = np.dot(total_e_field_vector, total_e_field_vector) ** 0.5

    chosen_axis = np.array(atom2, dtype="float") - np.array(atom1, dtype="float")

    # project individual atom's contribution to e field to the chosen axis
    atoms_projected_e_field_vectors = get_atoms_projected_e_field_vectors(atoms_e_field_vectors, chosen_axis)
    total_projected_e_field_vectors = np.sum(atoms_projected_e_field_vectors, axis=0)

    # sum individual atoms' contribution to PROJECTED e field to the chosen axis within the same residue to get
    # residues' contribution to PROJECTED e field
    residues_projected_e_field_vectors = get_residues_projected_e_field_vectors(atoms_projected_e_field_vectors,
                                                                                residue_ids)

    # write a pdb file for field direction visualization
    # this pdb contains two dummy atoms, one at the position of r0, the other is 5A away from r0, and the line connecting r0
    # is along the direction of the total field

    # write_vis_pdb(total_e_field_vector, r0)

    print("Total electirc field vector: (MV/cm)")
    print("{:10.3f}{:10.3f}{:10.3f}".format(total_e_field_vector[0], total_e_field_vector[1], total_e_field_vector[2]))
    print("Total electirc field modulus: (MV/cm)")
    print("{:10.3f}".format(modulus_total_e_field_vector))
    print("Projection of the total electric field vector along atom_1 and atom_2 (MV/cm)")
    print("{:10.3f}".format(total_projected_e_field_vectors))

    # debug file
    with open("debug_out", "w") as debug_file:
        for i in range(0, atom_num):
            line = "{:10.3f} {:10.3f} {:10.3f} {:10d} {:10.5f}".format(coords[i, 0], coords[i, 1], coords[i, 2],
                                                                       residue_ids[i], charges[i])
            debug_file.write(line + "\n")
    return total_e_field_vector, modulus_total_e_field_vector, total_projected_e_field_vectors
