#!/bin/python

from lxml import etree
import Bio.Data.IUPACData as IUPACD

"""
<MutateResidue name="(&string)" target="(&string)" new_res="(&string)" preserve_atom_coords="(false &bool)" mutate_self="(false &bool)" update_polymer_bond_dependent="(false &bool)" />

"""


def GenRosettaScript(Relax=False):
    root = etree.Element("ROSETTASCRIPTS")

    RequiredFields = [
        "SCOREFXNS",
        "RESIDUE_SELECTORS",
        "TASKOPERATIONS",
        "SIMPLE_METRICS",
        "FILTERS",
        "MOVERS",
        "PROTOCOLS",
        "OUTPUT"
    ]

    for Field in RequiredFields:
        w = etree.SubElement(root, Field)
        if Field == "MOVERS":
            movers = w
        if Field == "PROTOCOLS" and Relax:
            relax = etree.SubElement(w, "Idealize")

    return root, movers


def WriteRosettaScript(root, outputpath):
    s = etree.tostring(root, pretty_print=True).decode("utf-8")

    with open(outputpath, 'w') as f:
        f.write(s)


def MutationEntry(movers, FromAA, ToAA, Index):

    mut = etree.SubElement(movers, "MutateResidue")

    ToAA3lcode = IUPACD.protein_letters_1to3[ToAA]

    mut.set("name", "%s%i%s" % (FromAA, Index, ToAA))
    mut.set("target", "%s" % Index)
    mut.set("new_res", ToAA3lcode)


def CreateMutationMover(Mutations,
                        moveroutputpath,
                        Workingdirectory=""):

    RosettaScript, movers = GenRosettaScript()
    for (pos, fromAA, toAA) in Mutations:
        MutationEntry(movers, fromAA, toAA, b)

    WriteRosettaScript(RosettaScript, moveroutputpath)
