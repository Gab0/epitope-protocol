#!/bin/python

import pypdb
import os
import random
import time
from Bio import PDB

import argparse

from ..StructureUtils import BasicStructureOperations as BSO
from .. import Evaluator
from ..Wrapper import Rosetta, Gromacs


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omp", dest="OMP", action="store_true")
    parser.add_argument("--dd", dest="DockingDecoys", type=int, default=10000)
    parser.add_argument("-d", dest="WorkingDirectory", default="")

    return parser.parse_args()


def DownloadPdb(code_name):
    content = pypdb.get_pdb_file(code_name, filetype="pdb")
    with open("%s.pdb" % code_name, 'w') as f:
        f.write(content)


def RunRCSBSearch():
    Results = pypdb.Query("Antibody").search()
    Results = ["1ynt", "1kzq"]
    random.shuffle(Results)

    for pdb_code in Results:
        DownloadPdb(pdb_code)
        result = Rosetta.BuildAntibody("%s.pdb" % pdb_code)
        if not result:
            print("Success!")
            break


def TestAntibody(options):
    # Rosetta.RunAntibody("ab_chains/sag1.fasta")
    # Rosetta.PrepareStructure("dock_input/SAG1.pdb")
    InputStructures = [
        "out.pdb"
    ]

    # Rosetta.RunAntibodyH3()
    outputpdb = "ABSAG1_NOWATER.pdb"
    gmx_dir = Gromacs.getDirectoryName()

    Gromacs.RunGromacs(gmx_dir, "H3_modeling/model-0.relaxed_0001.pdb")
    BSO.RemoveSolvent(Gromacs.GetLatestSimulationStep(gmx_dir), outputpdb)

    W = Evaluator.EvaluateRamachandran(outputpdb)
    T = 0.014403962351065048  # Evaluator.EvaluateRamachandran("1emt.pdb")

    print(W)
    print(T)


def InputsToUname(InputStructures):
    parts = [
        os.path.splitext(os.path.split(x)[-1])[0]
        for x in InputStructures]

    return "_DOCK_".join(parts)


def TestDocking(AntigenPath, AntibodyPath, options):
    # RemoveSolvent("GMX_RUN/step40c.pdb")

    InputStructures = [
        AntigenPath,
        AntibodyPath
    ]

    uname = InputsToUname(InputStructures)
    print(uname)
    inputpdb = "%s.pdb" % uname

    BSO.MergePDB(InputStructures, inputpdb)

    chains = BSO.getPDBChains(inputpdb)
    antigen_chains = BSO.GetAntigenChains(chains)

    if antigen_chains is None:
        return

    BSO.MoveChains(inputpdb, antigen_chains[0], "H")

    Rosetta.ScorePdb(inputpdb)

    Rosetta.RunPrepack(inputpdb)
    partners = BSO.GetPartnersFromChains(chains)

    if not partners:
        return

    base_score = ReadScore(inputpdb)

    Rosetta.RunDocking(inputpdb, partners, nbdecoys=options.DockingDecoys)
    dock_kwargs = {
        "inputpdb": inputpdb,
        "partners": partners,
        "nbdecoys": options.DockingDecoys
    }


def TestEval(inputpdb, options):
    for i in range(1, options.DockingDecoys):
        w = os.path.splitext(inputpdb)
        pdbname = "%s_%.4i" % (w[0], i) + w[1]
        print("Scoring %s" % pdbname)
        Rosetta.ScorePdb(pdbname)


def ReadScore(f):
    try:
        time.sleep(1)
        ScoreData = Rosetta.ParseScore(f)
        s = ScoreData["total_score"].min()
        print((f))
        print("T: %.5f" % s)
        if ["I_sc"] in ScoreData.columns:
            Interface = ScoreData["I_sc"].min()
            print("I: %.5f" % Interface)

        return s
    except Exception as e:
        print("Error on %s" % f)
        print(e)


def ReadAllScores():
    for f in sorted(os.listdir()):
        if f.endswith(".scr"):
            ReadScore(f)

# FIXME: Makes it complicated to organize ligands & receptors files;
def DockTwoDirectories(AntigenDirectory, AntibodyDirectory, options):
    for agf in os.listdir(AntigenDirectory):
        for abf in os.listdir(AntibodyDirectory):
            InputFiles = [agf, abf]
            if any([f.startswith("_") for f in InputFiles]):
                continue
            InputPaths = [
                os.path.join(p, f)
                for p, f in
                zip([AntigenDirectory, AntibodyDirectory], InputFiles)
            ]
            if sum([os.path.isfile(f) for f in InputPaths]) < 2:
                continue
            TestDocking(*InputPaths, options)


def main():
    options = parse_arguments()
    # TestAntibody(options)

    AGD = os.path.join(options.WorkingDirectory, "Antigens")
    ABD = os.path.join(options.WorkingDirectory, "Antibodies")
    DockTwoDirectories(AGD, ABD, options)
    ReadAllScores()


if __name__ == "__main__":
    main()
