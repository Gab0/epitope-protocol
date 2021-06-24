#!/bin/python


import numpy as np


def WriteMatrix(PlaneContent, OutputFilePathNoExt):
    np.save(OutputFilePathNoExt + ".npy", PlaneContent)


def CompareMatrices(Matrices):
    shapes = list(set([m.shape for m in Matrices]))
    print(shapes)
    #assert len(shapes) == 1

    Counter = 0
    L = len(Matrices)
    Similarities = np.zeros(shape=(L, L))
    m0 = Matrices[0].shape
    for i in range(m0[0]):
        for j in range(m0[1]):
            for z in range(m0[2]):
                Values = [m[i, j, z] for m in Matrices]
                if any(Values):
                    Counter += 1
                    for v1, V1 in enumerate(Values):
                        for v2, V2 in enumerate(Values):
                            if V1 == V2:
                                Similarities[v1, v2] += 1

    return Similarities / Counter


def LoadMatrix(filepath):
    return np.load(filepath, allow_pickle=True)
