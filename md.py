# -*- coding: utf-8 -*-
"""
Single Molecule Molecular Dynamics Code 
Created 2018 by David of Theoretically Speaking
Please Modify!
"""
from __future__ import print_function
import os
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


def main():
    """Handle input/output and molecular dynamics velocity-verlet algorithm"""

    print("hello")

# Execute code if main file
if __name__ == "__main__":
    main()

