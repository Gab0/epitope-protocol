#!/bin/python


from subprocess import Popen
import os
import io
import string
import shutil
import random
import pandas as pd
import multiprocessing


def RunRosettaExecutable(progname,
                         RosettaArguments,
                         omp=True,
                         Overwrite=False,
                         ScoreFileName=None,
                         OutputPdbFileName=None,
                         WorkingDirectory=None):
    ROSETTA_BASE =\
        "/home/Database/Rosetta/rosetta_bin_linux_2019.35.60890_bundle/"

    ROSETTA_BIN = ROSETTA_BASE + "main/source/bin"

    BUILD_TYPE = "default" if not omp else "omp"
    RELEASE = "linuxclangrelease"

    RosettaProgramsBase = [
        "antibody_designer",
        "antibody_numbering_converter",
        "snugdock",
        "antibody",
        "antibody_H3",
        "docking_protocol",
        "relax",
        "docking_prepack_protocol",
        "score_jd2",
        "rosetta_scripts"
    ]

    def executable_name(progname):
        p = progname
        p += ".%s.%s" % (BUILD_TYPE, RELEASE)
        return p

    RosettaPrograms = {
        p: os.path.join(ROSETTA_BIN, executable_name(progname))
        for p in RosettaProgramsBase
    }
    executable = RosettaPrograms[progname]

    # -- Parse Common Rosetta Arguments
    AdditionalArguments = []
    if Overwrite:
        AdditionalArguments.append("-overwrite")

    if ScoreFileName:
        AdditionalArguments += ["-out:file:scorefile", ScoreFileName]

    if OutputPdbFileName:
        AdditionalArguments += ["-out:prefix", OutputPdbFileName]

    if not WorkingDirectory:
        WorkingDirectory = os.getcwd()

    CMD = [executable] + RosettaArguments + AdditionalArguments
    Process = Popen(CMD, cwd=WorkingDirectory)
    Process.communicate()
    return Process


def ScorePdb(pdb_path):
    Arguments = [
        "-in:file:s",
        pdb_path,
        "-out:file:scorefile",
        pdb_path + "strut_score.scr"
    ]

    return RunRosettaExecutable("score_jd2", Arguments)


def BuildAntibody(pdb_path=None):
    Arguments = [
        "-s",
        pdb_path,
    ]
    # ExtraArgs = "-primary_cdrs H3 -graft_design_cdrs H3 -seq_design_cdrs H1 H2 -light_chain lambda -nstruct 1"
    # ExtraArgs1 = "-primary_cdrs H3 -graft_design_cdrs H3 -seq_design_cdrs H1 H2 -light_chain lambda -do_dock -nstruct 1"

    # Arguments += ExtraArgs.split(" ")

    return RunRosettaExecutable("antibody", Arguments)


def RunAntibodyGrafting(chainseqfasta):
    Arguments = [
        "-fasta",
        chainseqfasta,
        "-optimal_graft"
    ]

    return RunRosettaExecutable("antibody", Arguments)


def RunRosettaScript(inputpdb, outputpdbprefix, scriptpath, **kwargs):
    Arguments = [
        "-s", inputpdb,
        "-parser:protocol", scriptpath,
    ]

    return RunRosettaExecutable("rosetta_scripts",
                                Arguments,
                                Overwrite=True,
                                ScoreFileName=outputpdbprefix + ".scr",
                                OutputPdbFileName=outputpdbprefix,
                                **kwargs)


def RunAntibodyH3():
    output_directory = "H3_modeling"
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    Arguments = [
        "-s", "grafting/model-0.relaxed.pdb",
        "-nstruct", "1",
        "-multiple_processes_writing_to_one_directory",
        "-antibody:auto_generate_kink_constraint",
        "-antibody:all_atom_mode_kink_constraint",
        "-out:path:pdb", output_directory,
    ]

    return RunRosettaExecutable("antibody_H3", Arguments)


def RunPrepack(pdbfile):
    Arguments = [
        "-in:file:s",
        pdbfile,
        # "-out:pdb",
        # structurepdb,
        # "-overwrite"
    ]
    r = RunRosettaExecutable("docking_prepack_protocol", Arguments)
    prefix = os.path.splitext(pdbfile)[0]
    shutil.move("%s_0001.pdb" % prefix, pdbfile)
    return r


def RunMulticore(nbproc, fn, dock_kwargs):

    pool = multiprocessing.Pool(processes=nbproc)

    dock_kwargs["with_uname"] = True
    pool.starmap(fn, [dock_kwargs for i in range(nbproc)])

    pool.close()
    pool.join()


def RunDocking(structurepdb, partners, with_uname=False, nbdecoys=10):
    if with_uname:
        uname = "".join([random.choice(string.ascii_lowercase)
                         for x in range(4)])
    else:
        uname = ""

    InputArguments = [
        "-in:file:s", structurepdb,
        "-partners", partners
    ]

    DockArguments = [
        "-spin",
        "-randomize1",
        "-randomize2",
        "-dock_pert", "1", "3",
        "-docking_local_refine",
        # "-use_ellipsoidal_randomization", "true",
        "-docking:sc_min",
        "-ex1",
        "-ex2aro",
    ]

    OutputArguments = [
        "-overwrite",
        "-score:docking_interface_score", "1",
        "-nstruct", str(nbdecoys),
        "-out:prefix", uname,
    ]

    Arguments = InputArguments + DockArguments + OutputArguments

    return RunRosettaExecutable(
        "docking_protocol",
        Arguments,
        ScoreFileName=structurepdb + "_dock_score.scr")


def PrepareStructure(pdbpath):
    Arguments = [
        "-in:file:s",
        pdbpath
    ]
    return RunRosettaExecutable("relax", Arguments)


def ParseScore(scrpath):
    w = open(scrpath).read().split("\n")
    output = []

    for line in w[1:]:
        line = line.split(" ")
        line = [x.strip(" ") for x in line if x.strip(" ")]
        line = line[1:]
        output.append(",".join(line))

    Buffer = io.StringIO()

    OutputContent = "\n".join(output)
    print(OutputContent)
    Buffer.write(OutputContent)
    Buffer.seek(0)

    return pd.read_csv(Buffer)


def MutateRosetta(PDBFile, mutations, ExpectedOutputFilepath):

    print("Looking for file %s ..." % ExpectedOutputFilepath)
    if os.path.isfile(ExpectedOutputFilepath):
        print(f"Mutator file exists at \n\t{ExpectedOutputFilepath},\n\t\tskipping mutation...")
        return
    moverfilepath = os.path.join(
            WorkingDirectory,
        BaseName + ".xml"
    )

    RosettaMover.CreateMutationMover(
        mutations,
        moverfilepath
    )

    Rosetta.RunRosettaScript(
        PDBFile,
        BaseName + "_",
        BaseName + ".xml",
        WorkingDirectory=WorkingDirectory
    )

    assert os.path.isfile(ExpectedOutputFilepath)
