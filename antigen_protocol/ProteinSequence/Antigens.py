#!/bin/python

"""

This module combines an alignment
for given protein region into a single figure, containing mutation positions.

Additional data can also be considered:
- Epitope positions, when a .csv file containing IEDB.org data is provided.
- Model span length when a .pdb model of the protein is provided.



"""

import sys
import re
import os
import argparse
from typing import List

import colored

import pandas as pd
from Bio import SeqIO

from ..StructureUtils import BasicStructureOperations
from ..EpitopeDetection import ModelResidueSurface, ParseBepipred

from .. OutputConstructor import plotProteinSequence, plotEpitopeScore
from . import AminoAcidSwaps, ProteinSequenceAlignment

from ..Mutation.Types import MutationSummary, Mutation


class EpitopeMatchResult():
    def __init__(self):
        pass

def ParseSequence(FilePath, re_pattern=""):
    Seq = SeqIO.parse(FilePath, format="clustal")

    if re_pattern:
        Result = [s for s in Seq if re.findall(re_pattern, s.description)]
    else:
        Result = list(Seq)

    return Result


def processDirectory(Directory, Epitopes, FileNamePattern):
    Files = [
        F for F in os.listdir(Directory)
        if re.findall(FileNamePattern, F)
    ]

    BinaryScoreCounter = []

    for F in Files:
        FilePath = os.path.join(Directory, F)
        CompareAgainst = ParseSequence(FilePath)

        print(">> Processing file %s" % F)

        BinaryScores, _ = processFile(CompareAgainst, Epitopes)
        BinaryScoreCounter.append(BinaryScores)

        binaryscore =\
            100 * (BinaryScores[0] + BinaryScores[1]) / BinaryScores[2]

        matchscore =\
            100 * BinaryScores[0] / BinaryScores[2]

        if matchscore:
            print()
            print("Full match or zero match percentage for %s: %.2f%% (%i)"
                  % (F, binaryscore, BinaryScores[0]))
            print("Percentage of epitopes with full match for %s: %.2f%%"
                  % (F, matchscore))
            print("\n")
        else:
            print("no match...\n")

    return BinaryScoreCounter


def findEpitopeInSequences(Epitope, Sequences, SequenceMatches, Verbose=False):
    """

    Find a single epitope in a set of sequences.
    SequenceMatches is modified in-place! (good idea?)

    """

    EpitopeSequence = Epitope["Description"]
    EpitopeCount = 0

    EpitopeFound = False
    for s, FastaSequence in enumerate(Sequences):
        CurrentSequence = str(FastaSequence.seq).replace("*", "")
        if EpitopeSequence in CurrentSequence:
            EpitopeFound = True
            break

    if not EpitopeFound:
        return 0

    for s, FastaSequence in enumerate(Sequences):
        CurrentSequence = str(FastaSequence.seq).replace("*", "")

        Alignment = ProteinSequenceAlignment.AlignSequences(
            CurrentSequence,
            EpitopeSequence
        )

        EpitopeBounds =\
            ProteinSequenceAlignment.ProcessAlignmentPath(Alignment.path)

        ValidAlignment = len(EpitopeBounds) == 1

        if not ValidAlignment:
            print("Invalid: %s" % EpitopeSequence)
            print(Alignment)
            print()
            for E in EpitopeBounds:
                print(CurrentSequence[E[0]: E[1]])

            print()
            # raise(Exception("Invalid Alignment."))

        if ValidAlignment:
            MatchedSequence =\
                CurrentSequence[EpitopeBounds[0][0]: EpitopeBounds[0][1]]

            Warning_A = len(MatchedSequence) != len(EpitopeSequence)
            Warning_B = all([MatchedSequence[0] != EpitopeSequence[0],
                             MatchedSequence[-1] != EpitopeSequence[-1]])

            if Warning_A or Warning_B:
                print(Alignment)
                print("Match/Epitope:")
                print(MatchedSequence)
                print(EpitopeSequence)
                # raise(Exception("BAD ALIGNMENT"))
            EpitopeCount += 1
            SequenceMatches[s].append(MatchedSequence)

        if Verbose:
            if len(EpitopeBounds) > 1:
                print("%s (%i)" % (FastaSequence.id, len(EpitopeBounds)))
            else:
                print(FastaSequence.id)

    if EpitopeCount and Verbose:
        print(Epitope["Antigen Name"])
        print("%i/%i" % (EpitopeCount, len(Sequences)))
        print(Epitope["Description"])
        print()

    return EpitopeCount


def processFile(Sequences, EpitopeData: pd.DataFrame, Verbose=0):
    """

    Process an alignment file.
    Return statistics and MutationVectors/EpitopeVectors

    """

    if not Sequences:
        print("No sequences!")
        sys.exit()

    Sequences = list(Sequences)
    print(Sequences)
    print(type(Sequences[0]))
    # Accept both SeqIO and AlignIO records;
    SequencesAsStrings = [seq.seq for seq in Sequences]
    SequenceIDs = [seq.id for seq in Sequences]

    assert SequenceIDs

    FullMatchCount = 0
    ZeroMatchCount = 0

    SequenceMatches: List[List[str]] = [[] for s in Sequences]
    if EpitopeData is not None:
        EpitopeDataSize = EpitopeData.shape[0]
        for i in range(EpitopeData.shape[0]):
            Epitope = EpitopeData.iloc[i]
            EpitopeCount = findEpitopeInSequences(
                Epitope,
                Sequences,
                SequenceMatches,
                Verbose
            )
            if EpitopeCount == 0:
                ZeroMatchCount += 1
            elif EpitopeCount == len(Sequences):
                FullMatchCount += 1

        SequencesAndMatches = zip(SequencesAsStrings, SequenceMatches)

        BestEpitopeMatch = sorted(
            SequencesAndMatches,
            key=lambda x: len(x[1]),
            reverse=True)[0]

        BestSequence, BestMatch = BestEpitopeMatch
        assert len(list({len(s) for s in SequencesAsStrings})) == 1

        AllowedBestMatches = ValidateEpitopeMatches(BestSequence,
                                                    BestMatch)
    else:
        EpitopeDataSize = 0
        BestSequence = SequencesAsStrings[0]
        AllowedBestMatches = []

    K = len(SequencesAsStrings)
    J = len(SequenceIDs)
    if not K == J:
        print(K)
        print(J)
        print("NO SEQS")
        #sys.exit()

    MutationVector = CreateMutationVector(SequencesAsStrings,
                                          BestSequence,
                                          SequenceIDs)

    EpitopeVector = BuildMatches(
        SequencesAsStrings,
        BestSequence,
        SequenceIDs,
        AllowedBestMatches,
        MutationVector
    )

    Statistics = (FullMatchCount, ZeroMatchCount, EpitopeDataSize)
    Results = (BestSequence, EpitopeVector, MutationVector)

    return Statistics, Results


def ViewMatches(Sequence, EpitopeVector, MutationVector):

    EpitopeColors = [
        0,
        "light_red",
        "red",
        "magenta",
        "purple_4a",
        "purple_3",
        "purple_3",
        "purple_3",
        "purple_3",
        "purple_3"
    ]

    MutationColors = [
        0,
        "orange_4b",
        "dark_orange_3a",
        "gold_3a",
        "gold_3a",
        "gold_3a",
        "gold_3a",
        "gold_3a",
        "gold_3a",
    ]

    for i, c in enumerate(Sequence):
        BG = ""
        FG = ""
        if MutationVector[i]:
            BG = colored.bg(MutationColors[MutationVector[i].get_nbvar()])
        if EpitopeVector[i]:
            FG = colored.fg(EpitopeColors[EpitopeVector[i]])

        print(FG + BG + c + colored.attr("reset"), end="")

    print()


def ShowEpitope(BaseHighlightMatch, BaseHighlightPositions):
    BaseHighlightPositions = list(reversed(BaseHighlightPositions))

    color = "red"
    print(">> " + colored.fg(color), end="")
    for _, BaseIndex in enumerate(BaseHighlightPositions):
        # if abs(BaseHighlightPositions[b+1] - BaseIndex) == 1:

        BaseHighlightMatch = BaseHighlightMatch[:BaseIndex] +\
            colored.fg("yellow") +\
            BaseHighlightMatch[BaseIndex] +\
            colored.fg(color) +\
            BaseHighlightMatch[BaseIndex+1:]

    return BaseHighlightMatch


def ValidateEpitopeMatches(ReferenceSequence, EpitopeMatches):

    AllowedMatches = []
    for Match in EpitopeMatches:
        print(Match)
        # -- Is it a relevant or redundant match?
        ProcessedMatch = True
        for OtherMatch in EpitopeMatches:
            if Match in OtherMatch and len(Match) < len(OtherMatch):
                ProcessedMatch = False
                break
        if not ProcessedMatch:
            continue
        StartIndex = ReferenceSequence.index(Match)
        EndIndex = StartIndex + len(Match)

        AllowedMatches.append((Match, StartIndex, EndIndex))

    return AllowedMatches


def ProcessMutation(Index, ReferenceBase, AllSequenceBases, SequenceIDs):
    Mutation = 0

    # FIXME: Debug section
    K = len(AllSequenceBases)
    J = len(SequenceIDs)
    if not K == J:
        print(K)
        print(J)
        sys.exit()

    BaseSet = list(set(AllSequenceBases))
    if len(BaseSet) > 1:
        if "X" in BaseSet:
            if len(BaseSet) == 2:
                return Mutation

        Mutation = MutationSummary(ReferenceBase, Index)

        for idx, Base in enumerate(AllSequenceBases):
            Mutation.addvar(Base, SequenceIDs[idx])

        Mutation.SpecialSwap = AminoAcidSwaps.EvaluateSwap(BaseSet)

    return Mutation


def CreateMutationVector(AllSequences, ReferenceSequence, SequenceIDs):
    MutationVector = [0 for x in ReferenceSequence]

    for Index, _ in enumerate(ReferenceSequence):
        AllSequenceBases = [s[Index] for s in AllSequences]
        ReferenceBase = ReferenceSequence[Index]
        MutationVector[Index] = ProcessMutation(Index,
                                                ReferenceBase,
                                                AllSequenceBases,
                                                SequenceIDs)
    return MutationVector


def BuildMatches(AllSequences,
                 ReferenceSequence,
                 SequenceIDs: List[str],
                 EpitopeMatches,
                 MutationVector,
                 color="red"):
    ShownSequence = str(ReferenceSequence)
    WillShow = False

    EpitopeVector = [0 for x in ReferenceSequence]

    for Match, StartIndex, EndIndex in EpitopeMatches:
        for i in range(StartIndex, EndIndex):
            EpitopeVector[i] += 1

    for Match, StartIndex, EndIndex in EpitopeMatches:

        BaseHighlightPositions = []
        MatchMutations = []

        # FIXME: What is this?
        for b, _ in enumerate(Match):
            for aSequence in AllSequences:
                try:
                    PointIndexInSequence = StartIndex + b
                except ValueError:
                    continue

            # -- Check if mutations are happening;
            MatchMutation = MutationVector[PointIndexInSequence]

            if MatchMutation:
                WillShow = True

                BaseHighlightPositions.append(b)
                MatchMutations.append(MatchMutation)

        BaseHighlightMatch = ShowEpitope(Match, BaseHighlightPositions)

        OverlappedStart = False
        OverlappedEnd = False
        for om, s, e in EpitopeMatches:
            if om != Match:
                if s < EndIndex < e:
                    OverlappedEnd = True
                if s < StartIndex < e:
                    OverlappedStart = True

        MatchReplacement = ""
        if not OverlappedStart:
            MatchReplacement += colored.fg(color)
        MatchReplacement += BaseHighlightMatch
        if not OverlappedEnd:
            MatchReplacement += colored.attr("reset")

        print(MatchReplacement)
        print(colored.attr("reset"))
        ShownSequence = ShownSequence.replace(Match, MatchReplacement)

        for Mutation in MatchMutations:
            print(Mutation)

    if WillShow:
        print(ReferenceSequence)
        print()
        print(ShownSequence)
        print()
        ViewMatches(ReferenceSequence, EpitopeVector, MutationVector)

    return EpitopeVector


def read_epitope_file(fpath):
    if fpath:
        return pd.read_csv(fpath)
    else:
        return None


def create_isolate_record(alignment, mutation_vector, ModelVectorBounds, output_fpath):

    muts = filter(lambda x: x != 0, mutation_vector)
    print(ModelVectorBounds)
    try:
        model_offset = ModelVectorBounds[0]
    except TypeError:
        model_offset = (0, 1e6)

    records = []
    for mut in muts:
        position = mut.position - model_offset[0] + 1
        if position <= 0:
            continue
        if position > model_offset[1]:
            continue

        for var, isolates in mut.variations.items():
            isolates = ["_".join(x.split("_")[:-3]) for x in isolates]

            n = len(isolates)

            if n > 6:
                isolates = f"{n} isolados."
            else:
                isolates = "; ".join(isolates)

            records.append({
                "Mutação": f"{position}{var}",
                "Isolados Portadores": isolates
            })

    df = pd.DataFrame(data=records)

    df.to_csv(output_fpath, index=False)


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", dest="EpitopeFile", required=True)
    parser.add_argument("-d", dest="WorkingDirectory")
    parser.add_argument(
        "--alnf-pattern",
        dest="FileNamePattern",
        default=r"Protein_[\w\d]+.aln")

    parser.add_argument("-f", dest="AlignmentFile")
    parser.add_argument("-p", dest="ProteinModelFile")
    parser.add_argument("-c", dest="WantedChain")
    parser.add_argument("-x", "--xsize", dest="XSIZE", type=int, default=80)
    parser.add_argument(
        "--aln-pattern", dest="AlignmentDescriptionPattern",
        help="A regex pattern to be applied in the description " +
        "of each alignment sequence. Only matching descriptions " +
        "will be considered."
    )
    parser.add_argument("--bep", dest="BepipredFile")

    parser.add_argument("-o", "--output", dest="OutputFilePath")
    parser.add_argument("--ext", dest="OutputFigureExtension", default="eps")

    parser.add_argument("--isolate-record-file")
    return parser.parse_args()


def main():
    options = parse_arguments()

    EpitopeData = read_epitope_file(options.EpitopeFile)

    # processDirectory(options.WorkingDirectory,
    # EpitopeData, options.FileNamePattern)
    Alignment = ParseSequence(options.AlignmentFile,
                              options.AlignmentDescriptionPattern)

    _, (Sequence, EpitopeVector, MutationVector) =\
        processFile(Alignment, EpitopeData)

    ModelSequenceBounds = None
    if options.ProteinModelFile:
        wantedChain = ModelResidueSurface.loadChain(
            options.ProteinModelFile,
            options.WantedChain
        )

        print("\n\n")
        print("Searching pdb structure sequence in region sequence:\n")
        chainSequence = BasicStructureOperations.GetStructureFasta(
            wantedChain)
        print(">%s\n" % chainSequence)

        print(Sequence)

        ModelSequenceAlignment =\
            ProteinSequenceAlignment.AlignSequences(str(Sequence).replace("*", ""),
                                                    chainSequence,
                                                    mode="global")

        print(ModelSequenceAlignment.path)
        print(ModelSequenceAlignment)

        ModelSequenceBounds =\
            ProteinSequenceAlignment.ProcessAlignmentPath(ModelSequenceAlignment.path)

        print(ModelSequenceBounds)
        # ModelResidueSurface.freesasaSurfaceResidues(wantedChain)

        print(chainSequence)

    AlignmentFileName = os.path.split(options.AlignmentFile)[-1]
    Basename = os.path.splitext(AlignmentFileName)[0]
    OutputFilePath = Basename + "." + options.OutputFigureExtension

    plotProteinSequence.ShowSequence(
        Sequence,
        len(Alignment),
        AlignmentFileName,
        MutationVector,
        EpitopeVector,
        ModelSequenceBounds,
        options.XSIZE,
        OutputFilePath
    )

    if options.isolate_record_file:
        create_isolate_record(
            Alignment,
            MutationVector,
            ModelSequenceBounds,
            options.isolate_record_file
        )

    if options.BepipredFile and False:
        bepOutputFilepath = Basename + "_score." + options.OutputFigureExtension
        epitope_data = ParseBepipred.parse_bepipred(options.BepipredFile)
        plotEpitopeScore.plotEpitopeScores(
            [epitope_data],
            MutationVector,
            EpitopeVector,
            bepOutputFilepath
        )


if __name__ == "__main__":
    main()
