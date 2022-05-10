#!/bin/python


from setuptools import setup, find_packages


setup(
    name="antigen_protocol",
    packages=find_packages(),
    version="0.02",
    entry_points={
        'console_scripts': [
            # -- Molecular Docking;
            # GROMACS
            "autoantigens=antigen_protocol.Dock.Gromacs:main",

            # AutoDock 4/Vina
            "pdbdock=antigen_protocol.Dock.DockEngine:main",

            # -- Epitope Detection;
            "depitope=antigen_protocol.EpitopeDetection.DetectEpitope:main",
            "evalepitope=antigen_protocol.StructureUtils.ModelResidueSurface:main",
            "eptarea=antigen_protocol.ExecuteEpitopes:main",
            "autobepipred=antigen_protocol.EpitopeDetection.ParseBepipred:main",
            "multibepipred=antigen_protocol.MultiBepiPred:main",

            # -- Mutate Protein;
            "strutmut=antigen_protocol.Mutation.StructureMutator:main",

            # -- Structure Visualization;
            "strutgraphic=antigen_protocol.StructureMap:main",
            "strutmutview=antigen_protocol.StructureMutationViewer:main",

            # -- Antigen binding detection;
            "protantigen=antigen_protocol.ProteinSequence.Antigens:main",

            # -- Rosetta Operations;
            "rscore=antigen_protocol.RosettaOperation.Score:main",

            # -- Structure Information;
            "pdbinfo=antigen_protocol.StructureInfo:main",
            "evalstruct=antigen_protocol.Evaluator:main"
        ]
    }
)
