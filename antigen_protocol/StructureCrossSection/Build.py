#!/bin/python

"""

Takes a PDB structure and a list of wanted residues (as indexes).

Makes a series of parallel Stereographic Projections and sums them in order to get
a better 2D representation of the epitope surface area.


"""


import time
import hashlib

import numpy as np

from typing import List
from Bio import PDB

from . import Engine3D, Show
from ..StructureUtils import BasicStructureOperations

WPR = 3


class Timer():
    def __init__(self):
        self.TimeZero = time.time()

    def measure(self):
        return time.time() - self.TimeZero


def FindPlaneWantedAtoms(Structure: PDB.Structure.Structure,
                         EpitopeIndexes: List[int]):
    PlaneResidues = []
    # -- Find Key residues (the plane will intersect those)
    COM = BasicStructureOperations.center_of_mass(Structure)
    for r, residue in enumerate(Structure.get_residues()):
        if r in EpitopeIndexes:
            PlaneResidues.append(residue)
            if len(PlaneResidues) >= WPR:
                break

    assert PlaneResidues

    # -- Find Plane atoms
    # (from the Key residues... the plane will intersect those);

    PlaneAtoms = []

    for residue in PlaneResidues:
        DistancesFromCOM = [
            BasicStructureOperations.distance(COM, atom.get_coord())
            for atom in residue.get_atoms()
        ]
        DistAtoms = zip(residue.get_atoms(), DistancesFromCOM)
        BestAtom = sorted(DistAtoms, key=lambda da: da[1])[0][0]
        PlaneAtoms.append(BestAtom)

    return PlaneAtoms


def BuildAtomCoords(PlaneAtoms, PlaneConfig):
    #
    # -- Build plane;
    Plane = Engine3D.Plane2D([atom.get_coord() for atom in PlaneAtoms],
                             PlaneConfig["Resolution"],
                             PlaneConfig["SideSize"])

    PlaneCoordinates = Plane.GetCoordinates()
    PX, PY, coords = PlaneCoordinates.shape

    # -- Initialize plane content arrays;
    PlaneCoordinates = np.zeros(shape=(PlaneConfig["StackSize"], PX, PY, 3))

    for S in range(PlaneConfig["StackSize"]):
        PlaneCoordinates[S] = Plane.GetCoordinates()
        Plane = Plane.GetNextLayer()

    return PlaneCoordinates


def GetAtomRadius(Atoms):
    # -- Get all atom radius (squared);
    # Atom Radius data took from freesasa sources;
    AtomRadius = {
        "H": 1.1,
        "C": 1.7,
        "O": 1.4,
        "N": 1.6,
        "S": 1.8
    }

    AtomRadius = np.array([
        AtomRadius[BasicStructureOperations.GetAtomName(atom)]
        for atom in Atoms])

    return AtomRadius


def BuildPlaneContents(Structure, PlaneConfig, EpitopeIndexes, DrawMethod):

    # -- INITIALIZE TIMERS;
    TimePA = Timer()
    PlaneAtoms = FindPlaneWantedAtoms(Structure, EpitopeIndexes)
    print("%i seconds to find plane wanted atoms.\n\n" % TimePA.measure())

    TimePC = Timer()
    print("%i seconds to find plane coordinates.\n\n" % TimePC.measure())

    # PlotPlaneDebug(PlaneCoordinates, Plane.PlaneKeyPoints)

    # -- Get all atom coordinates;
    AllNearAtoms = [
        atom for atom in Structure.get_atoms()
        if BasicStructureOperations.distance(
                atom.get_coord(), PlaneAtoms[0].get_coord()) <= 10
    ]
    AllNearAtomCoords = np.array([atom.get_coord() for atom in AllNearAtoms])

    # -- Get all atom radius;
    AllNearAtomRadius = GetAtomRadius(AllNearAtoms)
    AllNearAtomRadiusSQ = AllNearAtomRadius ** 2

    EpitopeIdentifier = "%s@%s" % (Structure.id,
                                   "|".join([str(x) for x in EpitopeIndexes]))

    # -- Build 3D Topology from model;
    print("Building image for %s..." % EpitopeIdentifier)
    TimerIM = Timer()

    if DrawMethod == 0:
        Plane = Engine3D.Plane2D(
            [atom.get_coord() for atom in PlaneAtoms],
            PlaneConfig["Resolution"],
            PlaneConfig["SideSize"])

        PlaneContent =\
            BuildImageDirect(PlaneConfig, Plane,
                             AllNearAtomCoords, AllNearAtomRadius)

    elif DrawMethod == 1:
        PlaneCoordinates = BuildAtomCoords(PlaneAtoms, PlaneConfig)
        PlaneContent = BuildImageBySpace(
            AllNearAtomCoords,
            AllNearAtomRadiusSQ, PlaneCoordinates)

    elif DrawMethod == 2:
        PlaneCoordinates = BuildAtomCoords(PlaneAtoms, PlaneConfig)
        PlaneContent = BuildImageByAtom(PlaneConfig, AllNearAtomCoords,
                                        AllNearAtomRadiusSQ, PlaneCoordinates)

    PlaneDensity = np.sum(PlaneContent) / np.prod(PlaneContent.shape)

    print("Plane density is %.5f" % PlaneDensity)

    PlaneExtraInfo = {"density": PlaneDensity}
    PlaneInformation = {**PlaneConfig, **PlaneExtraInfo}

    print("Took %0f seconds to build the image." % TimerIM.measure())

    np.save("test1.npy", PlaneContent)

    return PlaneContent, PlaneInformation


def BuildImageDirect(PlaneConfig, Plane, Atoms, AtomRadius, Verbose=False):

    Sides = Plane.GetContentDataWidth()
    PlaneContent = np.zeros(shape=(
        Plane.AngstromToDataSize(PlaneConfig["StackSize"]),
        Sides,
        Sides), dtype=np.float32)

    for A, Atom in enumerate(Atoms):
        ProjectionTimer = Timer()
        Projection = Plane.ProjectPoint(Atom)

        ProjectionTimer.measure()
        AtomZVector = Projection - Atom

        ZMOD = int(
            np.mean(AtomZVector / Plane.ZVector) / PlaneConfig["Resolution"])

        VectorInPlane = Projection - Plane.CenterPoint
        PlaneDistance = BasicStructureOperations.distance(
            Plane.CenterPoint,
            Projection)

        DistanceFromPlane = BasicStructureOperations.distance(Projection, Atom)

        if Verbose:
            print(DistanceFromPlane)
            print(ZMOD)

        print(">ATOM %i/%i" % (A + 1, len(Atoms)))

        if DistanceFromPlane > PlaneConfig["StackSize"]:
            print("Not Drawn.")
            continue
        print()

        PlaneDistanceMod = PlaneDistance / PlaneConfig["Resolution"]

        if Verbose:
            print(Sides)
            print(PlaneDistanceMod)
            print()

        Xcos = Plane.FindAngleBetweenVectors(Plane.XVector, VectorInPlane)
        Ycos = Plane.FindAngleBetweenVectors(Plane.YVector, VectorInPlane)

        if Verbose:
            print(">cos")
            print(Xcos)
            print(Ycos)
            print()

        if Verbose:
            w = np.array(Plane.CenterPoint +
                         Xcos * Plane.XVector * PlaneDistance +
                         Ycos * Plane.YVector * PlaneDistance)
            print(w)
            print(Projection)

        XMOD = (PlaneDistanceMod * Xcos + Sides // 2)
        YMOD = (PlaneDistanceMod * Ycos + Sides // 2)
        # math.cos(math.radians(YAngle))) + Sides // 2

        if False:
            print(XMOD)
            print(YMOD)
            print()

        MatrixRadius = AtomRadius[A] / PlaneConfig["Resolution"]

        AtomCoordinateOnPlane = [ZMOD, XMOD, YMOD]

        DrawTimer = Timer()
        DrawAtomOnPlaneMask(AtomCoordinateOnPlane, MatrixRadius, PlaneContent)
        print("drawing in %.2fs" % DrawTimer.measure())

    return PlaneContent


def DrawAtomOnPlane(MatrixFocus, MatrixRadius, PlaneContent):
    Bounds = []

    ADD = np.float32(1.0)
    for i, p in enumerate(MatrixFocus):
        J = max(0, int(p - MatrixRadius))
        K = min(PlaneContent.shape[i] - 1, int(p + MatrixRadius))

        Bounds.append((J, K))

    for z in range(*Bounds[0]):
        for x in range(*Bounds[1]):
            for y in range(*Bounds[2]):
                Dist = BasicStructureOperations.distance(MatrixFocus,
                                                         (z, x, y))
                if Dist < MatrixRadius:
                    PlaneContent[z, x, y] += ADD


def DrawAtomOnPlaneMask(MatrixFocus, MatrixRadius, PlaneContent):
    z, x, y = PlaneContent.shape
    Z, X, Y = np.ogrid[:z, :x, :y]

    NonEmpty = np.float32(1.0)

    # -- MODE 1
    if False:
        distances = np.sqrt(
            (Z - MatrixFocus[0]) ** 2 +
            (X - MatrixFocus[1]) ** 2 +
            (Y - MatrixFocus[2]) ** 2
        )

    # -- MODE 2
    Origin = np.asarray(MatrixFocus)

    CoordinateDifferences = np.subtract(
        np.indices(PlaneContent.shape).T, Origin)
    # distances = np.linalg.norm(SUB, axis=-1)
    distances = np.sum(CoordinateDifferences ** 2, axis=-1)

    Mask = distances.T <= MatrixRadius
    PlaneContent[Mask] = NonEmpty


def BuildImageByAtom(PlaneConfig, AtomCoords, AtomRadiusSQ, PlaneCoordinates):
    PlaneContent = np.zeros(shape=PlaneCoordinates.shape[:-1],
                            dtype=np.float32)

    for A, Atom in enumerate(AtomCoords):
        AtomDistances = np.sum((PlaneCoordinates - Atom) ** 2, axis=-1)
        M = np.min(AtomDistances)

        if M < PlaneConfig["Resolution"] * 2:
            MatrixFocus = [c[0] for c in np.where(AtomDistances == M)]
            print(M)
            MatrixRadius = AtomRadiusSQ[A] / PlaneConfig["Resolution"]
            DrawAtomOnPlane(MatrixFocus, MatrixRadius, PlaneContent)

    return PlaneContent


def BuildImageBySpace(AtomCoords, AtomRadiusSQ, PlaneCoordinates):
    PlaneContent = np.zeros(shape=PlaneCoordinates.shape[:-1],
                            dtype=np.float32)
    PCStack, PCX, PCY = PlaneContent.shape
    for S in range(PCStack):
        if not S % 5:
            print("%.2f%%" % (100 * S / PCStack))
        for i in range(PCX):
            for j in range(PCY):
                w = np.sum((AtomCoords -
                            PlaneCoordinates[S, i, j]) ** 2, axis=-1)

                Inside = int(np.where((w <= AtomRadiusSQ), 1, 0).any())

                PlaneContent[S, i, j] = Inside

    return PlaneContent


def CompressPlaneContent(PlaneContent):
    N_Stack = PlaneContent.shape[0]

    Output = np.zeros(shape=PlaneContent.shape[1:])
    Output[0, 0] = 1
    PlaneContent = np.flip(PlaneContent, axis=0)
    for i, Layer in enumerate(PlaneContent):
        print(i)
        v = 1 - (i / N_Stack)
        v = v ** 2
        for x in range(Output.shape[0]):
            for y in range(Output.shape[1]):
                if Layer[x, y]:
                    if Output[x, y] == 0:
                        Output[x, y] = v

        # PlaneContent[i] = C * np.log(abs(i -N_Stack))
        # Output = np.where(
        # PlaneContent[i] and not Output, PlaneContent[i], Output)

    if PlaneContent.sum():
        print(Output.sum())
        if Output.sum() == 0:
            exit("Cannot compress plane.")
    # return np.mean(PlaneContent, axis=0)
    return Output


def MakePlaneArrayHash(Array):
    return hashlib.sha1(Array.data).hexdigest()
