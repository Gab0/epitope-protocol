#!/bin/python

import matplotlib.pyplot as plt
import matplotlib as mpl

import numpy as np

from PIL import Image


def PlotPlaneImageColored(PlaneContent3D, OutputFilePathNoExt, Format="eps"):

    print("Plotting colored plane.")
    for i in range(PlaneContent3D.shape[0]):
        PlaneContent = np.repeat(PlaneContent3D[i][:,:,np.newaxis], 3, axis=2)
        img = Image.fromarray(PlaneContent, "RGB")
        OutputFilePath = OutputFilePathNoExt + "_%i.%s" % (i, Format)
        img.save(OutputFilePath, "png")
        print("Saving %s" % OutputFilePath)

    OutputFilePath = OutputFilePathNoExt + ".%s" % Format


def PlotPlaneImage(PointContents, OutputFilePath, Format="eps"):

    PlaneContent = np.repeat(PointContents[:, :, np.newaxis], 3, axis=2)
    img = Image.fromarray(PlaneContent, mode="RGB")
    OutputFilePath += ".%s" % Format

    img.save(OutputFilePath, "png")


def PlotPlaneDebug(PlaneStack,  KeyPoints):
    fig, (ax) = plt.subplots(projection='3d')

    for P in PlaneStack:
        ax.plot_surface(*P.T)

    ax.scatter3D(*KeyPoints.T, cmap='Reds')
    plt.show()


def PlotAtomDebug(PlaneContent):

    Points = []
    S, X, Y = PlaneContent.shape
    for s in range(S):
        for x in range(X):
            for y in range(Y):
                if PlaneContent[s, x, y]:
                    Points.append((s, x, y))

    for i in range(S):
        print(PlaneContent[i].shape)
        plt.matshow(PlaneContent[i])
        plt.show()

    ax = plt.axes(projection='3d')
    ax.scatter3D(*np.array(Points).T)
    plt.show()
