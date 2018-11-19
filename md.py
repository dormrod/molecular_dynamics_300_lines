# -*- coding: utf-8 -*-
"""
Single Molecule Molecular Dynamics Code 
Created 2018 by David of Theoretically Speaking
Please Modify!
"""
from __future__ import print_function
import os
import sys
import numpy as np


# Global variables for unit conversions
hartree = 4.35974465e-18   # J, atomic unit of energy
emass = 5.486e-4           # kg
dalton = 1.660539040e-27   # kg
avo = 6.02214086e23        # mol^-1
emass = 9.109534e-28       # g, atomic unit of mass
boltz = 1.38064852e-23 / hartree    # E_h K^-1
bohr = 0.52917721067       # Angstroms
hbar = 6.626070040e-34     # Js
atomic_time = hbar / hartree


# Global files to prevent constant opening/closing
xyz_file = open("coordinates.xyz", "w")
energy_file = open("energies.dat", "w")


def display_header():
    """Write opening message to screen"""

    print_dashed_line()
    print("Welcome to the Theoretically Speaking molecular dynamics code")
    print_dashed_line()


def print_dashed_line(length = 65):
    """Write --- line of given length to screen"""

    line = "-" * length
    print(line)


def string_to_boolean(string):
    """Converts input string of True or False to a boolean True or False"""

    string = string.lower().strip()
    true_strings = ["true", "t"]
    false_strings = ["false", "f"]
    if string in true_strings: return True
    elif string in false_strings: return False
    raise ValueError("Bad Boolean Value: " + string)


def get_input_parameters():
    """Ask user for input file name, read input parameters and store in dictionary"""

    # Get list of available input files
    input_files = get_recursive_file_list("inpt")

    # Ask user to select input file from list
    if len(input_files) == 0: # If cannot find any input files close program
        print("No available input files. Exiting.")
        sys.exit()
    else:
        while True:
            print("Select an input file from the list:")
            for i, file in enumerate(input_files):
                print("[{0}]  {1}".format(i, file))
            try:
                user_selection = int(input())
                input_file = input_files[user_selection]
                print("Input file selected: {0}".format(input_file))
                print_dashed_line()
                break
            except: pass

    # Open input file and read parameters into dictionary
    parameters = {}
    with open(input_file, "r") as file:
        print("Reading input file")
        # Skip header
        for i in range(2): file.readline()
        # Simulation parameters
        try:
            for i in range(2): file.readline()
            parameters["time_total"] = float(file.readline().split()[0]) / (atomic_time * 1e12)
            parameters["time_step"] = float(file.readline().split()[0]) / (atomic_time * 1e12)
            parameters["box_size"] = float(file.readline().split()[0]) / bohr
            parameters["write_freq"] = float(file.readline().split()[0]) / (atomic_time * 1e12)
            print("  - Simulation parameters read")
        except:
            print("Error in simulation parameters")
            sys.exit()
        # Atom data
        try:
            for i in range(2): file.readline()
            num_atoms = parameters["num_atoms"] = int(file.readline().split()[0])
            parameters["random_displacement"] = string_to_boolean(file.readline().split()[0])
            parameters["random_displacement_limit"] = float(file.readline().split()[0]) / bohr
            file.readline() # skip comment
            name_to_index = {}  # dictionary to convert atom name to array index
            parameters["atom_names"] = []  # empty list for names
            parameters["atom_masses"] = np.empty(num_atoms)  # empty array for masses
            parameters["atom_crds"] = np.empty([num_atoms, 3])  # empty array for coordinates
            for i in range(num_atoms):
                line = file.readline().split()
                name_to_index[line[0]] = i
                parameters["atom_names"].append(line[0])
                parameters["atom_masses"][i] = float(line[1]) / (avo * emass)
                parameters["atom_crds"][i] = np.array(line[2:5], dtype = float) / bohr
            print("  - Atom data read")
        except:
            print("Error in atom data")
            sys.exit()
        # Bond Data
        try:
            for i in range(2): file.readline()
            num_bonds = parameters["num_bonds"] = int(file.readline().split()[0])
            file.readline() # skip comment
            parameters["bond_pairs"] = np.empty([num_bonds, 2], dtype=int) # empty array for indices of bonded atom pairs
            parameters["bond_params"] = np.empty([num_bonds, 2]) # empty array for harmonic bond r0 and k
            for i in range(num_bonds):
                line = file.readline().split()
                parameters["bond_pairs"][i, 0] = name_to_index[line[0]]
                parameters["bond_pairs"][i, 1] = name_to_index[line[1]]
                parameters["bond_params"][i, 0] = float(line[2]) / bohr
                parameters["bond_params"][i, 1] = float(line[3]) * (bohr * 1e-10)**2 / hartree
            print("  - Bond data read")
        except:
            print("Error in bond data")
            sys.exit()
        print("Read successful")
        print_dashed_line()
    return parameters


def get_recursive_file_list(ext):
    """Get list of files with specifed extension in current directory and all subdirectories"""

    # Search over all files in all subdirectories, add to list if have required extension
    files = []
    for dirpath, dirname, filenames in os.walk("./"):
        for filename in filenames:
            if filename.endswith(ext):
                filepath = os.path.join(dirpath,filename)
                files.append(filepath)
    return files


def apply_periodic_boundary_condition(crds, box_size):
    """Apply periodicity to keep atoms within simulation box"""

    crds[crds < 0] += box_size
    crds[crds > box_size] -= box_size
    return crds


def minimum_image_displacement(crd_0, crd_1, box_size):
    """Find displacement between nearest periodic images of atom pair"""

    displacement = crd_0 - crd_1
    displacement[displacement < -box_size / 2] += box_size
    displacement[displacement > box_size / 2] -= box_size
    return displacement


def initialise_coordinates(crds, box_size, displace, limit):
    """Recentre atoms in simulation box, apply periodic boundary, apply random displacement"""

    crds += box_size / 2
    crds = apply_periodic_boundary_condition(crds, box_size)
    if displace:
        displacements = np.random.uniform(low = -limit, high = limit, size = crds.shape)
        crds += displacements
    return crds


def calculate_energy(masses, crds, velocities, bond_pairs, bond_params, box_size):
    """Calculate kinetic, potential and total energy of system"""

    kinetic_energy = 0.5 * (masses * np.sum(velocities ** 2, axis=1)).sum()  # U=0.5*m*v^2

    # Calculate harmonic potential energy using: U=0.5*k(r-r0)^2
    for i, bond in enumerate(bond_pairs):
        atom_0, atom_1 = bond[0], bond[1]
        displacement = minimum_image_displacement(crds[atom_0, :], crds[atom_1, :], box_size)
        distance = np.linalg.norm(displacement)
        potential_energy = 0.5 * bond_params[i, 1] * (distance - bond_params[i, 0]) ** 2

    total_energy = kinetic_energy + potential_energy  # Total energy as sum of ke and pe
    return np.array([kinetic_energy, potential_energy, total_energy])


def update_accelerations(masses, crds, bond_pairs, bond_params, box_size):
    """Calculate the acceleration on each atom using potential model and Newton's laws of motion"""

    # Calculate forces using Hooke's law: F=-k(r-r0)
    # Convert to acceleration using Newton's laws: F=ma, action has opposite reaction
    accelerations = np.zeros_like(crds)  # x,y,z accelerations for each atom
    for i, bond in enumerate(bond_pairs):
        atom_0, atom_1 = bond[0], bond[1]
        displacement = minimum_image_displacement(crds[atom_0, :], crds[atom_1, :], box_size)
        distance = np.linalg.norm(displacement)
        force_direction = displacement / distance
        force_magnitude = - bond_params[i, 1] * (distance - bond_params[i, 0])
        force = force_magnitude * force_direction
        accelerations[atom_0] += force / masses[atom_0]
        accelerations[atom_1] -= force / masses[atom_1]
    return accelerations


def update_coordinates(crds, accelerations, velocities, time_step, box_size):
    """Update coordinates using: x(t+dt)=x(t)+v(t)*dt+0.5*a(t)*dt**2"""

    crds += velocities * time_step + 0.5 * accelerations * time_step ** 2
    crds = apply_periodic_boundary_condition(crds, box_size)
    return crds


def update_velocities(velocities, accelerations_start, accelerations_end, time_step):
    """Update velocities using: v(t+dt)=v(t)+0.5*dt*(a(t)+a(t+dt))"""

    velocities += 0.5 * time_step * (accelerations_start + accelerations_end)
    return velocities


def write_output_files(time_step, num_atoms, names, crds, energies):
    """Writes coordinates in XYZ file type to 'coordinates.xyz'
    Write kinetic, potential and total energies to 'energies.dat'"""

    # Write XYZ file
    xyz_file.write("{0}  \n\n".format(num_atoms))
    for i, crd in enumerate(crds):
        xyz = crd * bohr
        xyz_file.write("{0}  {1:.6f}  {2:.6f}  {3:.6f}  \n".format(names[i], xyz[0], xyz[1], xyz[2]))

    # Write energies
    energy = energies * hartree * avo * 1e-3
    energy_file.write("{0}  {1}  {2}  {3}  \n".format(time_step, energy[0], energy[1], energy[2]))


def main():
    """Handle input/output and molecular dynamics velocity-verlet algorithm"""

    # Display opening message
    display_header()

    # Read user parameters from input file
    input_parameters = get_input_parameters()

    # Unpack parameters
    time_total = input_parameters["time_total"]
    time_step = input_parameters["time_step"]
    box_size = input_parameters["box_size"]
    write_freq = input_parameters["write_freq"]
    num_atoms = input_parameters["num_atoms"]
    displace_atoms = input_parameters["random_displacement"]
    displacement_limit = input_parameters["random_displacement_limit"]
    atom_names = input_parameters["atom_names"]
    atom_masses = input_parameters["atom_masses"]
    atom_crds = input_parameters["atom_crds"]
    bond_pairs = input_parameters["bond_pairs"]
    bond_params = input_parameters["bond_params"]

    # Recentre coordinates and apply displacements
    atom_crds = initialise_coordinates(atom_crds, box_size, displace_atoms, displacement_limit)

    # Initialise Molecular Dynamics Variables
    num_steps = int(time_total / time_step)  # total number of steps of md
    write_steps = int(write_freq / time_step)  # number of steps to write out results
    atom_vels = np.zeros_like(atom_crds)  # velocities in x,y,z directions for all atoms
    atom_acc_start = atom_acc_end = np.zeros_like(atom_crds)  # accelerations at start and end of time step
    atom_acc_start = update_accelerations(atom_masses, atom_crds, bond_pairs, bond_params, box_size)  # calculate initial accelerations
    system_energy = calculate_energy(atom_masses, atom_crds, atom_vels, bond_pairs, bond_params, box_size)  # calculate initial energies
    write_output_files(0, num_atoms, atom_names, atom_crds, system_energy)

    # Molecular dynamics
    print("Performing molecular dynamics simulation")
    for step in range(1, num_steps+1):
        # Velocity - Verlet algorithm
        atom_crds = update_coordinates(atom_crds, atom_acc_start, atom_vels, time_step, box_size)
        atom_acc_end = update_accelerations(atom_masses, atom_crds, bond_pairs, bond_params, box_size)
        atom_vels = update_velocities(atom_vels, atom_acc_start, atom_acc_end, time_step)
        atom_acc_start = atom_acc_end
        # Write coordinates and energies
        if step % write_steps == 0:
            system_energy = calculate_energy(atom_masses, atom_crds, atom_vels, bond_pairs, bond_params, box_size)
            write_output_files(step, num_atoms, atom_names, atom_crds, system_energy)
            print("Completion: {:.3f}%".format(100 * float(step) / num_steps))
    print_dashed_line()
    print("Simulation complete \nCoordinates written to coordinates.xyz \nEnergies written to energies.dat")
    print_dashed_line()


# Execute code if main file
if __name__ == "__main__":
    main()
