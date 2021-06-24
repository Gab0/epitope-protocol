#!/python

import sys

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP


def execute_dssp(fpath):
    p = PDBParser()
    structure = p.get_structure("S", fpath)
    model = structure[0]
    return DSSP(model, fpath)


if __name__ == "__main__":
    residues = execute_dssp(sys.argv[1])
    for res in residues:
        print(res)
