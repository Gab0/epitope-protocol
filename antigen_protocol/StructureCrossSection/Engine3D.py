#!/bin/python

import sympy
import itertools
import numpy as np
import math
from pyquaternion import Quaternion


def distance(a, b):
    D = np.sqrt(np.sum((a - b) ** 2))
    assert D != np.nan
    return D


class Plane2D():
    d = None

    def __init__(self, points, Resolution, SideSize):
        assert(len(points) == 3)

        assert Resolution
        assert SideSize

        self.Resolution = Resolution
        self.SideSize = SideSize

        # These three points define our plane;
        self.PlaneKeyPoints = np.array(points)

        # Checking every combination of KeyPoints is not needed.
        # But let's let it stay, maybe it becomes useful someday...
        for A, B, C in itertools.permutations(range(3), 3):
            v1 = self.PlaneKeyPoints[C] - self.PlaneKeyPoints[A]
            v2 = self.PlaneKeyPoints[B] - self.PlaneKeyPoints[A]

            v3 = self.PlaneKeyPoints[C] - self.PlaneKeyPoints[B]
            self.Equation = np.cross(v1, v2)
            a, b, c = self.Equation

            self.altXVector = v1
            self.altYVector = np.cross(v2, v3)

            # print(self.altYVector)
            # print(self.Equation)
            d = np.dot(self.Equation, self.PlaneKeyPoints[B])

            # print("D=%.4f" % d)

            if self.d is not None:
                assert(abs(abs(d) - abs(self.d)) <= 0.01)

            cos = self.FindAngleBetweenVectors(v1, v2)
            # print("angle: %.8f" % cos)
            self.d = d
            break

        # -- Calculate Distances between key plane points;
        WPR = 3
        KeyPlanePointsDistances = np.zeros(shape=(WPR, WPR))
        for i in range(WPR):
            for j in range(WPR):
                KeyPlanePointsDistances[i, j] =\
                    distance(self.PlaneKeyPoints[i],
                             self.PlaneKeyPoints[j])

        # -- Calculate the side of the coverage (in coordinate units) square.
        self.SideLengthCoverage = np.max(KeyPlanePointsDistances) * 4
        assert self.SideLengthCoverage != np.nan

        # -- Calculate Central point;
        self.CenterPoint = np.mean(self.PlaneKeyPoints, axis=0)
        assert(self.inPlane(self.CenterPoint))

        # -- Calculate 2D axis vectors;
        self.XVector = (
            self.PlaneKeyPoints[1] -
            self.PlaneKeyPoints[0]) / KeyPlanePointsDistances[1, 0]

        # - Rotate X Vector 90º to obtain Y vector;
        Q = Quaternion(axis=self.GetOrthogonalVector(), angle=math.radians(90))
        self.YVector = Q.rotate(self.XVector)

        self.ZVector = np.array(self.GetOrthogonalVector())

        #cos = self.FindAngleBetweenVectors(self.XVector, self.YVector)
        #print("cosine: %.4f" % cos)

        #PlotPointsDebug(coords, self.PlaneKeyPoints)
        self.GetOrthogonalVector()

    def FindAngleBetweenVectors(self, A1, A2):
        DOT = np.dot(A1, A2)

        def Mag(V):
            return np.sqrt(sum([w ** 2 for w in V]))

        MAG1 = Mag(A1)
        MAG2 = Mag(A2)

        cos = DOT / (MAG1 * MAG2)

        return cos
        return math.degrees(math.acos(cos))

    def GetContentDataWidth(self):
        assert self.Resolution != np.nan
        assert self.SideLengthCoverage != np.nan

        return int(self.SideLengthCoverage / self.Resolution)

    def AngstromToDataSize(self, A):
        return int(A / self.Resolution)

    def GetCoordinates(self):
        Sides = self.GetContentDataWidth()
        HalfSides = Sides // 2
        RES = np.float32(Sides / self.SideLengthCoverage)

        coords = np.zeros(shape=(Sides, Sides, 3), dtype=np.float32)

        def Coord(i, j):
            I = i - HalfSides
            J = j - HalfSides
            IT = self.XVector * I / RES
            JT = self.YVector * J / RES
            return self.CenterPoint + IT + JT

        for i in range(coords.shape[0]):
            for j in range(coords.shape[1]):
                coords[i, j] = Coord(i, j)

        return coords

    def inPlane(self, point):
        point = np.array(point)
        w = np.dot(self.Equation, point)
        #print("@: %.12f" % w)
        return abs(w - self.d) < 1e-3

    # Here I discovered sympy XD (so maybe port everything to sympy code?)
    def GetOrthogonalVector(self):
        Points = [sympy.Point3D(*pkp) for pkp in self.PlaneKeyPoints]
        P = sympy.Plane(*Points)

        return np.array([float(x) for x in P.normal_vector])

    def GetNextLayer(self):
        normal = np.array(self.GetOrthogonalVector())
        next_keypoints = self.PlaneKeyPoints + ((normal / 100) * self.Resolution)

        return Plane2D(next_keypoints, self.Resolution, self.SideSize)

    def ProjectPoint(self, point):
        SQR = sum([w ** 2 for w in self.Equation])
        a, b, c = self.Equation
        d, e, f = self.CenterPoint
        x, y, z = point
        t = a * d - a * x + b * e - b * y + c * f - c * z
        t /= SQR

        Projection = (x + t * a, y + t * b, z + t * c)
        return np.array(Projection)
