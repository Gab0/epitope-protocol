#!/bin/python

from typing import List, Tuple
from Bio import Align


def AlignSequences(seqA: str, seqB: str, mode="local"):
    """

    This method is the core of this module.
    Aligns the epitope sequence with the full structure chain sequence,
    using suitable scores.

    """
    Aligner = Align.PairwiseAligner()

    Aligner.mode = mode
    Aligner.open_gap_score = -10

    Aligner.internal_open_gap_score = -20
    # Aligner.external_open_gap_score = 0

    Aligner.extend_gap_score = -1
    Aligner.match_score = 300
    Aligner.mismatch_score = -0.5
    # Aligner.gap_score = -1e2

    d = Aligner.align(seqA, seqB)

    return d[0]


def ProcessAlignmentPath(AlignmentPath: List[Tuple[int, int]]):
    Bounds = []
    for i, (X, Y) in enumerate(AlignmentPath):
        if i < len(AlignmentPath) - 1:
            nX, nY = AlignmentPath[i + 1]
            if Y != nY:
                Bounds.append((X, nX))

    return Bounds


def RetrieveSequenceFromBounds(Bounds, Sequence):
    Output = ""
    for B in Bounds:
        Output += Sequence[B[0]:B[1]]
    return Output
