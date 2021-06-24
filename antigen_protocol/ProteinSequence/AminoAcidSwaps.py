#!/bin/python

from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP

aromatic_aas = "YWF"
charged_aas = "KRHDECY"
positive_aas = "HKR"
negative_aas = "CDEY"


def EvaluateSwap(Variations: str):
    AA_GROUPS = [aromatic_aas, charged_aas]
    GROUP_DESCRIPTORS = [
        "aromatic into non aromatic",
        "charged into not charged",
    ]

    Violators = [
        DESCRIPTOR
        for GROUP, DESCRIPTOR in zip(AA_GROUPS, GROUP_DESCRIPTORS)
        if len(list(set([V in GROUP for V in Variations]))) > 1
    ]

    return Violators
