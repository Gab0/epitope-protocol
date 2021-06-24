#!/bin/python

import freesasa
import re
from Bio import PDB

from ..StructureUtils import BasicStructureOperations
import argparse


class RosettaHydrogenClassifier(freesasa.Classifier):
    def _wclassify(self, residueName, atomName):
        if re.match('\s*H', atomName):
            return 'Hydrogen'
        if re.match('\s*N', atomName):
            return 'Nitrogen'

        return 'Not-nitrogen'

    def radius(self, residueName, atomName):
        if re.match('\s*H', atomName):  # Hydrogen
            return 1.1
        if re.match('\s*N', atomName):  # Nitrogen
            return 1.6
        if re.match('\s*C', atomName):  # Carbon
            return 1.7
        if re.match('\s*O', atomName):  # Oxygen
            return 1.4
        if re.match('\s*S', atomName):  # Sulfur
            return 1.8

        return 0


def loadChain(pdbpath, wantedChainName=None):
    parser = PDB.PDBParser()
    structure = parser.get_structure("struct", pdbpath)

    if wantedChainName:
        wantedChain = [
            chain for chain in structure.get_chains()
            if chain.id == wantedChainName][0]
    else:
        Chains = list(structure.get_chains())
        if len(Chains) == 1:
            wantedChain = Chains[0]
        else:
            raise(Exception("No chain name provided, " +
                            "yet the structure has more than one chain."))

    return wantedChain


def freesasaSurfaceResidues(wantedChain):
    Residues = list(wantedChain.get_residues())

    From = BasicStructureOperations.GetResidueIndex(Residues[0])
    To = BasicStructureOperations.GetResidueIndex(Residues[-1])

    Structure = freesasa.structureFromBioPDB(
        wantedChain,
        RosettaHydrogenClassifier(),
        # freesasa.Classifier(),
    )
    Result = freesasa.calc(Structure)

    ResidueSASA = []
    for w in range(From, To + 1):
        SelectionQuery = "res, chain %s and resi %i" % (wantedChain.id, w)
        area = freesasa.selectArea([SelectionQuery], Structure, Result)["res"]
        ResidueSASA.append(area)

    return ResidueSASA


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", dest="PDBFile", required=True)
    return parser.parse_args()


def main():
    options = parse_arguments()

    wantedChain = loadChain(options.PDBFile, "F")
    freesasaSurfaceResidues(wantedChain)
