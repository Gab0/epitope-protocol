#!/bin/python

import sys
import numpy as np
from Bio import PDB

import pyrama


def EvaluateRamachandran(pdbfile):
    normals, outliers = pyrama.calc_ramachandran([pdbfile])
    Values = pyrama.RAMA_PREF_VALUES
    scr = []
    for AAType in normals.keys():
        for x, y in zip(normals[AAType]["x"], normals[AAType]["y"]):
            x, y = [round(v) + 179 for v in [x, y]]
            scr.append(Values[AAType][x, y])

    return sum(scr) / len(scr)


def main():
    res = EvaluateRamachandran(sys.argv[1])
    print(res)


if __name__ == "__main__":
    main()
