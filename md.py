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

    print_dotted_line()
    print("Welcome to the Theoretically Speaking molecular dynamics code")
    print_dotted_line()


def print_dotted_line(length = 65):
    """Write --- line of given length to screen"""

    line = "-" * length
    print(line)


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
                input_file=input_files[user_selection]
                print("Input file: {0} selected".format(input_file))
                print_dotted_line()
                break
            except: pass


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


def main():
    """Handle input/output and molecular dynamics velocity-verlet algorithm"""

    # Display opening message
    display_header()

    # Read user parameters from input file
    input_parameters = get_input_parameters()


# Execute code if main file
if __name__ == "__main__":
    main()

