#!/bin/python

import os
import shutil
import random
import string
import gromacs
from subprocess import call


def getGromacsProtocolFile(filename):
    base = os.path.dirname(os.path.realpath(__file__))
    subdir = "gromacs_protocols"

    full_path = os.path.join(base, subdir, filename)

    print("Gromacs MDP at\n>> %s" % full_path)
    assert(os.path.isfile(full_path))

    return full_path


def GromacsStep(Arguments):
    result = call(Arguments)

    if result == 1:
        print("Failure!")
        exit(1)


def getDirectoryName():
    run_id = ''.join([random.choice(string.ascii_lowercase) for z in range(4)])
    return "GMX_RUN_" + run_id


def RunGromacs(directory, pdbpath):
    os.mkdir(directory)

    parameterspath = getGromacsProtocolFile("minim.mdp")

    pdbfilename = os.path.split(pdbpath)[-1]
    cloned_pdbpath = os.path.join(directory, pdbfilename)
    shutil.copy(pdbpath, cloned_pdbpath)

    original_dir = os.getcwd()
    os.chdir(directory)

    Arguments = [
        "gmx", "pdb2gmx",
        "-f", pdbfilename,
        "-ff", "amber03",
        "-water", "tip3p",
        "-o", "processed.gro",
        "-ignh"
    ]

    GromacsStep(Arguments)

    """
    gmx editconf -f jz4_ini.pdb -o jz4.gro
    """

    BoxArguments = [
        "gmx", "editconf",
        "-f", "processed.gro",
        "-o", "newbox.gro",
        "-bt", "cubic",
        "-d", "1.0"
    ]

    call(BoxArguments)

    gromacs.solvate(
        cp="newbox.gro",
        cs="spc216.gro",
        p="topol.top",
        o="solv.gro",
    )

    MinimizeArguments = [
        "gmx", "grompp",
        "-f", parameterspath,
        "-c", "solv.gro",
        "-r", "solv.gro",
        "-p", "topol.top",
        "-o", "em.tpr",
        "-maxwarn", "2"
    ]

    GromacsStep(MinimizeArguments)

    gromacs.mdrun(v=True, deffnm="em", nsteps=10000)

    gromacs.trjconv(
        s="em.tpr",
        f="em.trr",
        pbc="mol",
        ur="compact"
    )

    """
    gmx mdrun -v -deffnm em
    gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
    gmx mdrun -deffnm nvt
    gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr

    gmx mdrun -deffnm md_0_10

    """

    os.chdir(original_dir)


def GetLatestSimulationStep(Directory):
    all_files = os.listdir(Directory)
    files = [f for f in all_files if f.startswith("step")]
    files = sorted(files, reverse=True)
    return os.path.join(Directory, files[0])

