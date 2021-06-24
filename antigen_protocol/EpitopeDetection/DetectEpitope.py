#!/bin/python
"""

Detects epitopes in protein structure files (.pdb)

"""
from Bio import PDB

from ..StructureUtils import BasicStructureOperations
from . import ModelResidueSurface

import argparse
import colored

import numpy as np
import matplotlib.pyplot as plt
import Bio.Data.IUPACData as IUPACD


DefaultDistanceThreshold = 7.0


class HydrogenBond():
    def __init__(self, dist, group, degree, chains):
        self.Distance = dist
        self.AtomGroup = group
        self.Degree = degree
        self.InvolvedChains = chains

    def get_identifier(self):
        return "".join(self.AtomGroup)


def getDockedDegree(Distance, distanceThreshold):
    DockedRules = (
        Distance <= distanceThreshold / 1.5,
        Distance <= distanceThreshold,
        Distance <= distanceThreshold * 2
    )

    for r, Rule in enumerate(DockedRules):
        if Rule:
            return len(DockedRules) - r

    return 0


def ValidHydrogenBond(atomids, Verbose=False):
    names = [w[0] for w in atomids]

    hasH = any(n in "H" for n in names)
    hasEneg = any(n in "FNO" for n in names)

    if hasH and hasEneg:
        return True
    else:
        if Verbose:
            print(atomids)
    return False


def ShowData(Data, FieldSizes):
    for data, size in zip(Data, FieldSizes):
        Spacer = " " * (size - len(data))
        print(data + Spacer, end="")
    print()


def GetWantedChains(structure, chains):
    pass


def ShowEpitope(pdbpath,
                target_chains,
                distanceThreshold=DefaultDistanceThreshold,
                validateBond=True):

    parser = PDB.PDBParser()
    structure = parser.get_structure("X", pdbpath)

    # -- Select source and destination protein chains;
    ExcludedChains = [target_chains[0]]
    ExcludedChains += [
        chain.id for chain in structure.get_chains()
        if chain.id not in target_chains[1:]
    ]

    source_chain = [
        chain for chain in structure.get_chains()
        if chain.id == target_chains[0]
    ]

    dest_chains = [
        chain for chain in structure.get_chains()
        if chain.id not in ExcludedChains
    ]

    # source_chain, other_chain = GetWantedChains(structure, target_chain)
    if source_chain and dest_chains:
        source_chain = source_chain[0]
    else:
        print("ERROR! No suitable chain found on the structure file.")
        exit(1)

    # -- Initialize output field names;
    FieldColumns = (
        "IDX     ",
        "RES   ",
        "RES LETTER",
        "DIST",
        "DOCKED TO CHAIN",
        "DOCKED DIST RATING",
        "DOCKED ATOMS",
        "RESIDUE SASA"
    )

    FieldSizes = [len(f) + 5 for f in FieldColumns]
    MinimumResidueDistances = []
    EpitopeSequences = {"ALL": ""}
    WantedChainResidues = list(source_chain.get_residues())

    ResidueExposedAreas =\
        ModelResidueSurface.freesasaSurfaceResidues(source_chain)

    W = len(ResidueExposedAreas)
    K = len(WantedChainResidues)
    print(W - K)
    assert(W == K)

    # -- Iterate wanted chain residues,
    # -- printing data to screen & writing to file.
    OutputFile = open(pdbpath + "_epitope.tsv", 'w')

    Header = "\t".join(FieldColumns)
    ShowData(FieldColumns, FieldSizes)
    OutputFile.write(Header + "\n")

    for r, residue in enumerate(WantedChainResidues):
        ResidueDocked = False
        MinimumResidueDistance = 1e5

        DockedAtoms = []

        for res_atom in residue.get_atoms():
            for dest_chain in dest_chains:
                for other_atom in dest_chain.get_atoms():
                    Distance = BasicStructureOperations.distance(
                        res_atom.coord,
                        other_atom.coord
                    )

                    Updated = False
                    if Distance < MinimumResidueDistance:
                        MinimumResidueDistance = Distance
                        Updated = True
                    if Updated or Distance < 4.0:
                        CurrentAtomGroup = [
                            res_atom.get_name(),
                            other_atom.get_name()
                        ]
                        BondIsValid =\
                            not validateBond or ValidHydrogenBond(CurrentAtomGroup)

                        if BondIsValid:
                            MinimumResidueDistance = Distance
                            DockedDegree = getDockedDegree(
                                Distance,
                                distanceThreshold
                            )
                            if not DockedAtoms or DockedDegree:
                                DockedPair = HydrogenBond(
                                    Distance,
                                    CurrentAtomGroup,
                                    DockedDegree,
                                    [source_chain.id,
                                     dest_chain.id])

                                DockedAtoms.append(DockedPair)


                # -- continue earlier to improve performance
                other_chain_com = BasicStructureOperations.center_of_mass(
                    dest_chain)
                other_chain_dist = BasicStructureOperations.distance(
                    res_atom.coord,
                    other_chain_com
                )

                if other_chain_dist > 15:
                    continue
                if ResidueDocked == 1:
                    break
            if ResidueDocked == 1:
                break
        ResidueName = residue.resname

        ResidueLetter = IUPACD.protein_letters_3to1[ResidueName.title()]

        Colors = {
            1: "green",
            2: "yellow",
            3: "red"
        }

        DockedAtoms = sorted(DockedAtoms, key=lambda d: d.Distance)

        for sb, ShortestBond in enumerate(DockedAtoms):
            TargetChain = ShortestBond.InvolvedChains[-1]

            if ShortestBond.Degree:
                print(colored.fg(Colors[ShortestBond.Degree]), end="")
            # -- Append AA to epitope sequence record.
            if ResidueDocked >= 2:
                if TargetChain not in EpitopeSequences.keys():
                    EpitopeSequences[TargetChain] = ""
                EpitopeSequences[TargetChain] += ResidueLetter
                EpitopeSequences["ALL"] += ResidueLetter

            Fields = (
                "%i/%i" % (r + 1, len(WantedChainResidues)) if not sb else "",
                ResidueName if not sb else "",
                ResidueLetter if not sb else "",
                "%.2f" % ShortestBond.Distance,
                TargetChain,
                "*" * ResidueDocked,
                " ".join(ShortestBond.AtomGroup),
                "%.2f" % ResidueExposedAreas[r]
            )

            ROW = "\t".join(Fields)
            ShowData(Fields, FieldSizes)
            OutputFile.write(ROW + "\n")

            if ShortestBond.Degree:
                print(colored.attr("reset"), end="")

        MinimumResidueDistances.append(MinimumResidueDistance)

    OutputFile.close()
    return MinimumResidueDistances, EpitopeSequences


def ParseChains(chains: str):
    return chains.split(":")


# FIXME: UNDER DEVELOPMENT;
def WritePymolScript(options, MinimumResidueDistances):
    content = []
    content.append("load %s" % options.PdbFile)

    for i, c in enumerate(MinimumResidueDistances):
        if c < 5:
            pass


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", dest="PdbFile", required=True)
    parser.add_argument("-c", dest="WantedChain", required=True)
    parser.add_argument("--dist", dest="ThresholdDistance",
                        default=DefaultDistanceThreshold, type=float)
    parser.add_argument("--nvbond", dest="ValidateBond",
                        action="store_false", default=True)

    parser.add_argument("--plot", dest="Plot", action="store_true")

    return parser.parse_args()


def main():
    options = parse_arguments()

    ds, es = ShowEpitope(options.PdbFile,
                         ParseChains(options.WantedChain),
                         options.ThresholdDistance,
                         options.ValidateBond)

    print("Atom processed.")
    print(es)
    ds = np.array(ds)
    print(ds.shape)
    if options.Plot:
        plt.plot(np.arange(len(ds)), ds)
        plt.show()


if __name__ == "__main__":
    main()
