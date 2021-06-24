#!/bin/python

import rmsd
import copy

"""

Calculates RMSD distance between

"""
def CalculateRMSD(pdb_path_a: str, pdb_path_b: str):

    p_all_atoms, p_all = rmsd.get_coordinates(pdb_path_a, "pdb")
    q_all_atoms, q_all = rmsd.get_coordinates(pdb_path_b, "pdb")

    p_size = p_all.shape[0]
    q_size = q_all.shape[0]

    p_coord = copy.deepcopy(p_all)
    q_coord = copy.deepcopy(q_all)
    p_atoms = copy.deepcopy(p_all_atoms)
    q_atoms = copy.deepcopy(q_all_atoms)

    if not p_size == q_size:
        print("error: Structures not same size")
        quit()

    p_cent = rmsd.centroid(p_coord)
    q_cent = rmsd.centroid(q_coord)
    p_coord -= p_cent
    q_coord -= q_cent

    rotation_method = rmsd.kabsch_rmsd
    result_rmsd = rotation_method(p_coord, q_coord)
    return result_rmsd
