#!/bin/python

"""

Executable to check all observed versions of nucleotide sequences,
translating them to protein sequences, but only the segment present
on the structure file.

"""

import sys
import argparse
import os
import re
from typing import Tuple, List, Optional

import numpy as np
import Bio.SeqIO

from Bio import PDB

from ..StructureUtils import BasicStructureOperations, RMSD
from ..ProteinSequence import ProteinSequenceAlignment
from .Mutators import PymolMutator

from .Types import OutputMutationFile, Mutation


def ExtractMatchingSmallSequenceFromBigSequence(
        Big_str: str,
        Small_str: str):
    try:
        Z = Big_str.index(Small_str[:5])
    except ValueError:
        raise Exception("Sequence not found.")
    K = Z + len(Small_str)
    return Big_str[Z:K]


def CountSNP(A: str, B: str):
    return [
        (i, a, b)
        for i, (a, b) in enumerate(zip(A, B))
        if a != b
    ]


def ShowSNP(attr: Tuple[int, str, str]):
    return "%i: %s -> %s" % attr


def FindMutationHotspots(Sequences):
    Hotspots = []
    for i in range(len(Sequences[0])):
        try:
            Bases = [seq[i] for seq in Sequences]
            if len(list(set(Bases))) > 1:
                Hotspots.append(i + 1)
        except IndexError as e:
            print(f"Hotspots error: {e}")

    return Hotspots


def ValidateSeqRecords(X):
    return all([isinstance(x, Bio.SeqIO.SeqRecord) for x in X])


def buildSeqRecord(sequence: str,
                   base: Bio.SeqIO.SeqRecord) -> Bio.SeqIO.SeqRecord:

    S = Bio.Seq.Seq(sequence)
    return Bio.SeqIO.SeqRecord(S, id=base.id)


def extract_model_versions(ProteinAlignment):
    pass


def extract_structurally_relevant_mutations(
        ProteinAlignment, PDBSequence):

    AllModelSequences = [
        Bio.SeqIO.SeqRecord(Bio.Seq.Seq(PDBSequence), id="PDB Sequence")
    ]

    ModelVersions = [[]]

    for sequence in ProteinAlignment:
        AASequence = str(sequence.seq).replace("*", "")

        CroppedToModelAAAlignment =\
            ProteinSequenceAlignment.AlignSequences(AASequence,
                                                    PDBSequence)

        CroppedToModelAABounds =\
            ProteinSequenceAlignment.ProcessAlignmentPath(
                CroppedToModelAAAlignment.path)

        CroppedToModelAASequence =\
            ProteinSequenceAlignment.RetrieveSequenceFromBounds(
                CroppedToModelAABounds, AASequence)

        CroppedSequence = buildSeqRecord(
            CroppedToModelAASequence,
            base=sequence
        )

        AllModelSequences.append(CroppedSequence)

        snps = CountSNP(CroppedToModelAASequence, PDBSequence)
        if not snps:
            print("matches: %s" % sequence.id)
        else:
            print(" misses: %s(%i snp) %s" % (
                sequence.id + " " * (32 - len(sequence.id)),
                len(snps),
                " ".join([ShowSNP(snp) for snp in snps])
            ))
            if snps not in ModelVersions:
                ModelVersionIndex = len(ModelVersions)
                ModelVersions.append(snps)
                AllModelSequences.append(CroppedSequence)
            else:
                ModelVersionIndex = ModelVersions.index(snps)

            print("%s version: %i" % (sequence.id, ModelVersionIndex))

    if not ValidateSeqRecords(AllModelSequences):
        print("ERROR: Invalid SeqRecords")
        sys.exit(1)

    return ModelVersions, AllModelSequences


def loadStructureSequence(PDBFile) -> List[Tuple[int, str]]:
    parser = PDB.PDBParser()
    structure = parser.get_structure("strut", PDBFile)

    # chain = list(structure.get_chains())[0]
    return BasicStructureOperations.get_structure_fasta_and_indexes(structure)


def loadProteinAlignment(ProteinFastaFile):
    return Bio.SeqIO.parse(
        ProteinFastaFile,
        format="clustal"
    )


def extract_mutations(ReferenceSequence, AlternativeSequence):
    for idx, base in ReferenceSequence:
        if base != AlternativeSequence[idx - 1]:
            yield Mutation(idx, base,
                           AlternativeSequence[idx-1])


def read_mutations(Input: str) -> List[Tuple[int, str, str]]:
    lines = Input.split("\n")

    return list(filter(None, map(read_mutation, lines)))


def read_mutation(mutation: str) -> Optional[Tuple[int, str, str]]:
    Patterns = [
        r"\d+",
        r"^\D+",
        r"\D+$"
    ]
    Matches = [re.findall(p, mutation.strip(" \n")) for p in Patterns]
    try:
        [idx, a, b] = [M[0] for M in Matches]
        return (int(idx), a, b)
    except IndexError:
        if mutation:
            print("Failure to read mutation string.")
        return None


class BadMutationError(Exception):
    pass


def validate_mutation(ReferenceSequence, mutations, output_file):

    PDBSequence = loadStructureSequence(output_file)
    MutationsPositions = [m[0] for m in mutations]

    for b, refbase in ReferenceSequence:
        MatchingBase = [k for i, k in PDBSequence if i == b]
        if not MatchingBase:
            raise BadMutationError(
                "No matching index found between structures."
            )

        expected_mutation = Mutation(
            b,
            refbase,
            MatchingBase[0]
        )
        smut = expected_mutation.show()
        if refbase != MatchingBase[0]:
            if b in MutationsPositions:
                MutationsPositions.pop(MutationsPositions.index(b))
            else:
                raise BadMutationError(
                    f"Mutation {smut} was not predicted!"
                )

    if MutationsPositions:
        for MUT in MutationsPositions:
            print(MUT)
        raise BadMutationError("Intended mutations were not executed.")


def sort_rule(reference, seq) -> Tuple[int, str]:
    delta = 0
    for bref, bseq in zip(reference, seq):
        if bref != bseq:
            delta += 1

    return delta, seq


def sort_unique_sequences(sequences, reference_sequence) -> List[str]:
    """

    Guarantees a list of similar unique sequences will
    always be in the same order,
    while also starting with the reference sequence.

    """
    if all(isinstance(seq, Bio.SeqRecord.SeqRecord) for seq in sequences):
        sequences = [str(s.seq) for s in sequences]
    elif all(isinstance(seq, str) for seq in sequences):
        pass
    else:
        raise Exception("Unknown sequence type.")

    if reference_sequence is None:
        return sequences

    unique_sequences = set(sequences)
    unique_sequences.remove(reference_sequence)

    mutated_sequences = sorted(
        list(unique_sequences),
        key=lambda seq: sort_rule(reference_sequence, seq)
    )

    return [reference_sequence] + mutated_sequences


def get_pdb_name(filepath: str) -> str:
    return os.path.splitext(os.path.split(filepath)[-1])[0]


def RunStructureAgainstSequences(
        ProteinFastaFile: str,
        PDBFile: str,
        WorkingDirectory: str):

    """

    Index 0 of the lists ModelVersions & ModelSequences holds
    the unchanged sequence
    (no mutations: it's the sequence of the original PDB structure file).

    """
    modelName = get_pdb_name(PDBFile)
    ProteinAlignment = loadProteinAlignment(ProteinFastaFile)
    PDBSequence = loadStructureSequence(PDBFile)

    print("Running for model %s." % modelName)
    ModelVersions, AllModelSequences =\
        extract_structurally_relevant_mutations(ProteinAlignment, PDBSequence)

    unique_model_sequences = sort_unique_sequences(
        AllModelSequences,
        PDBSequence
    )

    assert PDBSequence == unique_model_sequences[0]
    RealWorkingDirectory = os.path.abspath(WorkingDirectory)
    OutputStructureFiles = []
    OutputVariations = []
    AllMutations = []

    for idx in range(len(ModelVersions)):
        mutations = list(extract_mutations(
            unique_model_sequences[0],
            unique_model_sequences[idx]
        ))

        AllMutations += mutations

        w = OutputMutationFile(
            mutations,
            GetModelPrefix(idx),
            os.path.join(
                RealWorkingDirectory,
                GetModelFilename(idx, modelName)
            )
        )
        OutputVariations.append(w)
        OutputStructureFiles.append(w)

    AllMutations = list(set(AllMutations))
    for mut in AllMutations:
        mut_str = mut.show()
        mut_prefix = f"mutation_{mut_str}"
        w = OutputMutationFile(
            [mut],
            mut_prefix, os.path.join(
                RealWorkingDirectory,
                mut_prefix + ".pdb")
        )

        OutputStructureFiles.append(w)

    # -- Create one model per mutation group.
    for idx, output_file in enumerate(OutputStructureFiles):
        write_mutated_structure(
            PDBFile,
            PDBSequence,
            output_file
        )

    # -- Locate and write mutation hotspot information;
    Hotspots = FindMutationHotspots(AllModelSequences)
    HotspotFilepath = os.path.join(WorkingDirectory, "Hotspots.txt")

    with open(HotspotFilepath, 'w') as f:
        f.write("\n".join([str(h) for h in Hotspots]))

    print("Building matrix...")
    try:
        Matrix = buildRMSDMatrix(
            [f.output_filepath for f in OutputStructureFiles],
            ModelVersions)
        print(Matrix)
    except Exception as e:
        print(e)

    return OutputVariations


def write_mutated_structure(
        PDBFile: str,
        PDBSequence: List[Tuple[int, str]],
        output_file: OutputMutationFile) -> None:
    print(f"Generating mutations to file {output_file.output_filepath}")

    for mut in output_file.mutations:
        print("\t" + mut.show_spread())

    if False:
        MUTATOR = None
        # MUTATOR = MutateRosetta

    else:
        MUTATOR = PymolMutator.MutatePymol

    MUTATOR(PDBFile, output_file.mutations, output_file.output_filepath)

    validate_mutation(
        PDBSequence,
        output_file.mutations,
        output_file.output_filepath
    )


# FIXME: DEPRECATED?
def buildRMSDMatrix(AllModelOutputs, ModelVersions):
    # -- Build RMSD matrix between all models;
    MatrixShape = len(AllModelOutputs)
    print(AllModelOutputs)
    print(MatrixShape)

    RMSDMatrix = np.zeros(shape=(MatrixShape, MatrixShape))
    for mi in range(MatrixShape):
        for mj in range(MatrixShape):
            if mi == mj:
                V = 0
            else:
                V = RMSD.CalculateRMSD(
                    AllModelOutputs[mi],
                    AllModelOutputs[mj]
                )

            RMSDMatrix[mi, mj] = V

    GlobalVariationCount = len(list(set(AllModelOutputs)))
    ModelVariationCount = len(ModelVersions)

    print("\n\n")
    print("Global Variation Count: %i" % GlobalVariationCount)
    print("Inside Model Variation Count: %i" % ModelVariationCount)
    print("\n")

    return RMSDMatrix


def GetModelPrefix(idx):
    if idx == 0:
        return "original_%s" % idx
    return "mutate_%s" % idx


def GetModelFilename(idx, modelName):
    return "_".join([GetModelPrefix(idx), modelName, "0001.pdb"])


def parse_arguments(description=__doc__, require_pdb=True):

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-a", dest="AlignmentFile",
                        help="Path to fasta file.")

    parser.add_argument("-p", dest="PDBFile",
                        required=require_pdb,
                        help="Path to PDB structure file.")

    parser.add_argument("-d", dest="WorkingDirectory", default="",
                        help="Path to output working directory.")

    parser.add_argument("--ext", dest="OutputFigureExtension",
                        default="eps")

    parser.add_argument("-e", dest="EpitopeFile")

    parser.add_argument(
        "-S",
        dest="show_strains",
        action="store_true",
        help="Show which strains displays each observed sequence version."
    )

    parser.add_argument(
        "-I",
        dest="mutation_recipe",
        help="Mutate protein according to input file." +
        "It should contain one mutation per line like 'N12A'"
    )

    return parser.parse_args()


def main():
    options = parse_arguments()

    if options.WorkingDirectory:
        if not os.path.isdir(options.WorkingDirectory):
            os.mkdir(options.WorkingDirectory)

    if options.mutation_recipe:
        PDBSequence = loadStructureSequence(options.PDBFile)

        mutations = read_mutations(open(options.mutation_recipe).read())

        PDBName = os.path.splitext(os.path.split(options.PDBFile)[-1])[0]

        output_path = os.path.join(
                os.path.abspath(options.WorkingDirectory),
                PDBName + "_" + options.mutation_recipe + ".pdb"
            )

        output_file = OutputMutationFile(
            mutations,
            'W',
            output_path
        )

        write_mutated_structure(options.PDBFile, PDBSequence, output_file)
    else:
        RunStructureAgainstSequences(
            options.AlignmentFile,
            options.PDBFile,
            options.WorkingDirectory
        )
