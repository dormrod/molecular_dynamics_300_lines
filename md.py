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
atomic_time = hbar/hartree


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
                print("[{0}]  {1}".format(i,file))
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
            n_atoms = parameters["n_atoms"] = int(file.readline().split()[0])
            parameters["random_displacement"] = string_to_boolean(file.readline().split()[0])
            file.readline() # skip comment
            name_to_index = {}  # dictionary to convert atom name to array index
            parameters["atom_names"] = []  # empty list for names
            parameters["atom_masses"] = np.empty(n_atoms)  # empty array for masses
            parameters["atom_crds"] = np.empty([n_atoms, 3])  # empty array for coordinates
            for i in range(n_atoms):
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
            n_bonds = parameters["n_bonds"] = int(file.readline().split()[0])
            file.readline() # skip comment
            parameters["bond_pairs"] = np.empty([n_bonds, 2], dtype=int) # empty array for indices of bonded atom pairs
            parameters["bond_params"] = np.empty([n_bonds, 2]) # empty array for harmonic bond r0 and k
            for i in range(n_bonds):
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
    crds[crds < 0.0] += box_size
    crds[crds > box_size] -= box_size
    return crds

def main():
    """Handle input/output and molecular dynamics velocity-verlet algorithm"""

    # Display opening message
    display_header()

    # Read user parameters from input file
    input_parameters = get_input_parameters()


# Execute code if main file
if __name__ == "__main__":
    main()

