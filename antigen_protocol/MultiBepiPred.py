#!/bin;python3

"""

Takes a batch of protein sequences,

"""

from io import StringIO
import os
import hashlib
import pandas as pd

from . import StructureMutator
from . import ReadDSSP

from .EpitopeDetection import bepipred_post
from .OutputConstructor import plotEpitopeScore
from .ProteinSequence.Antigens import processFile, read_epitope_file
from .StructureMutator import sort_unique_sequences
from .StructureUtils.BasicStructureOperations import GetStructureFasta, loadPDB

WORD_VAR = "Variação"


def ordered_set(data):
    output = []
    for k in data:
        if k not in output:
            output.append(k)

    return output


def main():
    options = StructureMutator.parse_arguments(__doc__)

    ProteinAlignment = StructureMutator.loadProteinAlignment(
        options.AlignmentFile
    )
    PDBStructure = loadPDB(options.PDBFile)
    PDBSequence = GetStructureFasta(
        PDBStructure
    )

    _, AllModelSequences =\
        StructureMutator.ExtractStructurallyRelevantMutations(
            ProteinAlignment,
            PDBSequence
        )

    _, (_, EpitopeVector, MutationVector) =\
        processFile(AllModelSequences,
                    read_epitope_file(options.EpitopeFile))

    UniqueSequences = ordered_set([
        str(seq.seq)
        for seq in AllModelSequences
    ])

    UniqueSequences = sort_unique_sequences(UniqueSequences, PDBSequence)

    DSSPVector = ReadDSSP.execute_dssp(options.PDBFile)

    # FIXME: These checking methods are deprecated;
    # assert Sequences[0] == AllModelSequences[0]
    # print(AllModelSequences.index(Sequences[0]))
    for Method in bepipred_post.Methods:
        predict_sequences(
            options,
            UniqueSequences,
            Method,
            DSSPVector,
            MutationVector,
            EpitopeVector
        )
        break


def get_code(seq):
    return hashlib.md5(seq.encode("utf-8")).hexdigest()


def execute_prediction(Sequences, PredictionMethod):
    curves = []
    for s, SEQ in enumerate(Sequences):
        print(f"Sequence {s + 1} of {len(Sequences)}")
        SEQ_CODE = get_code(SEQ)
        SEQ_FILE = f"epitope_pred_{SEQ_CODE[:4]}.csv"

        if os.path.isfile(SEQ_FILE):
            print("Loaded from file.")
            data = pd.read_csv(SEQ_FILE)
        else:
            print("Fetching from remote server...")

            content = bepipred_post.get_prediction(
                str(SEQ),
                prediction_method=PredictionMethod
            )

            buffer_content = StringIO(content.text)
            with open(SEQ_FILE, 'w') as f:
                f.write(content.text)
            data = pd.read_csv(buffer_content)

        curves.append(data)

    return curves


def predict_sequences(options,
                      UniqueSequences,
                      PredictionMethod,
                      DSSPVector,
                      MutationVector,
                      EpitopeVector):

    curve_labels = [
        "Original"
    ] + [
        f"{WORD_VAR} #{i + 1}"
        for i in range(len(UniqueSequences[1:]))
    ]

    curves = execute_prediction(UniqueSequences, PredictionMethod)
    name_prefix = os.path.splitext(os.path.split(options.AlignmentFile)[-1])[0]

    output_filepath = \
        f"{name_prefix}_multi_{PredictionMethod}." + \
        f"{options.OutputFigureExtension}"

    plotEpitopeScore.plotEpitopeScores(
        curves,
        curve_labels=curve_labels,
        output_filepath=output_filepath,
        MutationVector=MutationVector,
        DSSPVector=DSSPVector,
        EpitopeVector=EpitopeVector
    )


if __name__ == "__main__":
    main()
