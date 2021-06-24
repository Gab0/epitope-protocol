#!/bin/python

import warnings
import argparse
from . StructureUtils import BasicStructureOperations as BSO

def generate_vinaconf(center, side):
    AXIS = ['x', 'y', 'z']
    groups = ["center", "size"]

    Parameters = []
    for g, group in enumerate(groups):
        for a, AX in enumerate(AXIS):
            p = "_".join([group, AX])
            if g:
                v = side
            else:
                v = center[a]
            Parameters.append((p, v))

    return Parameters


def render_vina_parameter(V):
    (w, v) = V
    return w + " = " + str(v)



def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', dest="LigandPDB", default="")
    parser.add_argument('-o', dest="OutputPath", default="")
    parser.add_argument('-s', dest="SearchBoxSide", default=25)
    parser.add_argument('-r', dest="ReceptorPDB", required=True)
    return parser.parse_args()


def main():
    options = parse_arguments()

    message = "Reading center of mass of the "
    if options.LigandPDB:
        struct_path = options.LigandPDB
        message += "ligand."
    else:
        struct_path = options.ReceptorPDB

    print(message)
    read_structure(struct_path, options.SearchBoxSide, options.OutputPath)


def read_structure(struct_path, search_box_side, output_path):

    message = f"Molecule at {struct_path}"
    print(message)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        struct = BSO.loadPDB(struct_path)

    CM = BSO.center_of_mass(struct)

    print(CM)

    Parameters = generate_vinaconf(CM, search_box_side)
    output = map(render_vina_parameter, Parameters)
    vinaconf = "\n".join(output)

    print(vinaconf)

    if output_path:
        with open(output_path, 'w') as f:
            f.write(vinaconf)
