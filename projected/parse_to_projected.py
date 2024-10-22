import sys
import os
import argparse
from typing import Iterable, Tuple

from pysmt.environment import reset_env, get_env
from pysmt.fnode import FNode
from pysmt.shortcuts import read_smtlib

from allsat_cnf.label_cnfizer import LabelCNFizer
from allsat_cnf.polarity_cnfizer import PolarityCNFizer
from allsat_cnf.utils import get_allsat, is_cnf, SolverOptions, check_sat
from allsat_cnf.utils import get_lra_atoms, get_boolean_variables, check_models

def parse_args():
    parser = argparse.ArgumentParser(description='Enumerate models of formulas')
    parser.add_argument('input', help='Folder with .smt2 files')
    parser.add_argument('-o', '--output', action='store',
                        help='Output folder where to save the DIMACS files')
    return parser.parse_args()

def process_file(input_file, output_folder):
    phi = read_smtlib(input_file)

    atoms = get_boolean_variables(phi)

    phi_cnf = PolarityCNFizer(nnf=True, mutex_nnf_labels=True,
                              label_neg_polarity=False).convert_as_formula(phi)
    
    additional_atoms = get_boolean_variables(phi_cnf) - atoms
    
    n_clauses = len(phi_cnf.args())

    print(f"Processing {input_file}")
    print("Number of clauses: ", n_clauses)
    print("CNF formula: ", phi_cnf.serialize())
    print("Atoms: ", atoms)

    counter = 1
    atom_to_dimacs = {}

    for a in atoms:
        atom_to_dimacs[a] = counter
        counter += 1

    for a in additional_atoms:
        atom_to_dimacs[a] = counter
        counter += 1


    # Create the output file path
    output_file = os.path.join(output_folder, os.path.basename(input_file).replace('.smt2', '.cnf'))
    with open(output_file, "w") as f:
        f.write("p cnf " + str(len(atom_to_dimacs)) + " " + str(n_clauses) + "\n")
        f.write("c p show " + " ".join([str(i) for i in range(1, len(atoms) + 1)]) + " 0\n")
        for c in phi_cnf.args():
            if len(c.args()) < 2:
                if c.is_not():
                    atom = c.arg(0)
                    f.write("-" + str(atom_to_dimacs[atom]) + " 0\n")
                else:
                    atom = c
                    f.write(str(atom_to_dimacs[atom]) + " 0\n")
                continue
            for l in c.args():
                if l.is_not():
                    atom = l.arg(0)
                    f.write("-" + str(atom_to_dimacs[atom]) + " ")
                else:
                    atom = l
                    f.write(str(atom_to_dimacs[atom]) + " ")
            f.write("0\n")

def main():
    sys.setrecursionlimit(10000)
    args = parse_args()

    input_folder = args.input
    output_folder = args.output if args.output else input_folder + "_DIMACS_NEW"

    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Process each .smt2 file in the input folder
    for filename in os.listdir(input_folder):
        if filename.endswith('.smt2'):
            input_file = os.path.join(input_folder, filename)
            process_file(input_file, output_folder)

if __name__ == '__main__':
    main()
