#!/bin/python

import os
import argparse
import json

from .StructureCrossSection import Build, Show, MatrixRecord
from .StructureUtils import BasicStructureOperations


def ReadHotspotFile(Filepath):
    with open(Filepath) as f:
        content = f.read().split("\n")
        hotspots = [int(line) for line in content if line]

    return hotspots


def ProcessDirectory(options):
    files = [
        os.path.join(options.WorkingDirectory, f)
        for f in os.listdir(options.WorkingDirectory)
        if f.endswith(".pdb")
    ]

    HotspotFilepath = os.path.join(options.WorkingDirectory, "Hotspots.txt")
    Hotspots = ReadHotspotFile(HotspotFilepath)

    H = Hotspots[0]

    print(H)
    PlaneResidues = [H - 1, H, H + 1]
    for f in files:
        OutputImageFilepath = os.path.splitext(f)[0]
        ProcessEpitopeFile(f, options, PlaneResidues, OutputImageFilepath)


def GetFilename(path):
    return os.path.splitext(os.path.split(path)[-1])[0]


def ProcessEpitopeFile(pdbfilepath,
                       options,
                       PlaneResidues: [int],
                       OutputFilePathNoExt: str):

    pdbname = GetFilename(pdbfilepath)
    Structure = BasicStructureOperations.loadPDB(pdbfilepath, pdbname)

    # Those measures are in Angstroms (or maybe just PDB coordinates);
    PlaneConfig = {
        "Resolution": options.PlaneResolution,
        "SideSize": options.SideSize,
        "StackSize": options.StackSize
    }

    PlaneContent, PlaneInformation =\
        Build.BuildPlaneContents(Structure, PlaneConfig, PlaneResidues, 0)

    CompressedPlaneContent = Build.CompressPlaneContent(PlaneContent)

    hash_file = os.path.join(options.WorkingDirectory, "plane_hashes.txt")
    with open(hash_file, 'a') as f:
        f.write("%s\n" % Build.MakePlaneArrayHash(CompressedPlaneContent))

    plane_file = os.path.join(options.WorkingDirectory, "plane_info.txt")
    with open(plane_file, 'a') as f:
        f.write(json.dumps(PlaneInformation, indent=2) + "\n")

    print("Plotting plane...")

    Show.PlotPlaneImage(CompressedPlaneContent, OutputFilePathNoExt)
    print("Drawing image...")

    Show.PlotPlaneImageColored(
        PlaneContent,
        OutputFilePathNoExt + "_color",
        Format="png")

    MatrixRecord.WriteMatrix(PlaneContent, OutputFilePathNoExt)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", dest="WorkingDirectory")
    parser.add_argument("-c", dest="CompareMatrices", action="store_true")
    parser.add_argument("--resolution", dest="PlaneResolution",
                        type=float, default=0.1)

    parser.add_argument("--depth", dest="StackSize",
                        type=int, default=10)

    parser.add_argument("--size", dest="SideSize",
                        type=int, default=20)

    return parser.parse_args()


def main():
    options = parse_arguments()
    if options.CompareMatrices:
        Matrices = [
            MatrixRecord.LoadMatrix(
                os.path.join(options.WorkingDirectory, f))
            for f in os.listdir(options.WorkingDirectory)
            if f.endswith(".npy")
        ]
        Similarities = MatrixRecord.CompareMatrices(Matrices)
        print(Similarities)
    else:
        ProcessDirectory(options)
