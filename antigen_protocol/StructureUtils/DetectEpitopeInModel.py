#!/bin/python

import pandas as pd

from . import ModelResidueSurface
from ..StructureUtils import BasicStructureOperations

from .ProteinSequenceAlignment import AlignSequences, ProcessAlignmentPath


def ReadStructureRecord(filepath: str):
    df = pd.read_csv(filepath, index_col=False, sep="\t")
    return df


def StructureRecordToStructureString(StructureRecord: pd.DataFrame):
    StructureString = "".join([StructureRecord.iloc[i]["RESIDUE LETTER"]
                               for i in range(StructureRecord.shape[0])])

    return StructureString


def LocateEpitopeString(EpitopeString: str,
                        Structure,
                        Verbose: bool):

    # StructureRecord = ReadStructureRecord(StructureRecordFilepath)
    # print(StructureRecord.to_string())

    StructureString = BasicStructureOperations.GetStructureFasta(Structure)

    Alignment = AlignSequences(StructureString, EpitopeString)

    AlignmentBounds = ProcessAlignmentPath(Alignment.path)

    for Bound in AlignmentBounds:
        Fragment = StructureString[Bound[0]:Bound[1]]
        assert(Fragment in StructureString)

    return AlignmentBounds


if __name__ == "__main__":
    structurefile = "../../Proteins/1ynt.pdb"
    w = LocateEpitopeString("TALLAYIKGD", structurefile, True)
    ModelResidueSurface.freesasaSurfaceResidues(structurefile)
