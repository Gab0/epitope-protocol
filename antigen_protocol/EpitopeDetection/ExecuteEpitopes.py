#!/bin/python

import argparse
from .StructureUtils import ModelResidueSurface, DetectEpitopeInModel
import pandas as pd
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", dest="EpitopeFile", required=True)
    parser.add_argument("-s", dest="StructureFile", required=True)
    parser.add_argument("--chain", dest="WantedChainName", required=True)
    return parser.parse_args()


def main():
    options = parse_arguments()

    epitopes = pd.read_csv(options.EpitopeFile)

    wantedChainStructure = ModelResidueSurface.loadChain(
        options.StructureFile,
        options.WantedChainName
    )

    ResidueExposedSurfaces =\
        ModelResidueSurface.freesasaSurfaceResidues(wantedChainStructure)

    ResidueSequence =\
        ModelResidueSurface.PDBOps.GetStructureFasta(wantedChainStructure)

    for E in range(epitopes.shape[0]):
        Epitope = epitopes.iloc[E]
        EpitopeSequence = Epitope["Description"]

        if Epitope["Object Type"] != "Linear peptide":
            continue

        EpitopeBounds = DetectEpitopeInModel.LocateEpitopeString(
           EpitopeSequence, wantedChainStructure, True)

        # Linear epitopes should always have a single sequence;
        if len(EpitopeBounds) > 1:
            continue

        EpitopeSequenceInModel = "".join([ResidueSequence[b[0]:b[1]]
                                          for b in EpitopeBounds])
        print()
        print("e: %s" % EpitopeSequence)
        print("s: %s" % EpitopeSequenceInModel)
        print("l: %i" % len(EpitopeBounds))
        print("m: %s" % (EpitopeSequence == EpitopeSequenceInModel))

        ExposedAreas = [
            ResidueExposedSurfaces[b[0]:b[1]]
            for b in EpitopeBounds
        ]

        meanExposedAreas = np.mean([np.mean(e) for e in ExposedAreas])

        print(meanExposedAreas)
        print(np.mean(np.array(ResidueExposedSurfaces)))
        print()
