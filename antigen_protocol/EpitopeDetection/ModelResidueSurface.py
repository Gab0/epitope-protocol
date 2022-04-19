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
        radii = {
            "H": 1.1,  # Hydrogen
            "N": 1.6,  # Nitrogen
            "C": 1.7,  # Carbon
            "O": 1.4,  # Oxygen
            "S": 1.8   # Sulfur
        }

        for symbol, radius in radii.items():
            if re.match(fr'\s*{symbol}', atomName):
                return radius

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
