#!/bin/python

import os
import sys
import subprocess
import datetime
import argparse
from .. import StructureInfo
EXECUTABLE = "./qvina-w"

def is_mutation_strut(fpath):
    EXT = os.path.splitext(fpath)[-1]
    for n in ["original", "mutat"]:
        if n in fpath and "pdb" in EXT:
            return True

    return False


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", dest="Ligand", required=True)

    return parser.parse_args()


def main():
    options = parse_arguments()
    LIGAND = options.Ligand
    for F in os.listdir():
        if is_mutation_strut(F):

            RECEPTOR = os.path.splitext(os.path.split(F)[-1])[0]
            W = datetime.datetime.now().strftime("%H:%M")
            print(f"[{W}]: runnint {F}")
            CONF_PATH = f"run{F}.conf"

            StructureInfo.read_structure(F, 80, CONF_PATH)

            subprocess.call([EXECUTABLE,
                             "--ligand", LIGAND,
                             "--receptor", F,
                             "--config", CONF_PATH,
                             "--out", f"{RECEPTOR}_{LIGAND}_out.pdbqt"
                             ])


if __name__ == "__main__":
    main()
