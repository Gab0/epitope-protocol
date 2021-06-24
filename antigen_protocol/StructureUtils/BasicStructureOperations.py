#!/bin/python

from Bio import PDB
import math
from typing import List, Tuple

import Bio.Data.IUPACData as IUPACD

ExpectedChains = ["L", "H"]


def GetAntigenChains(chains: List[str]):
    for echain in ExpectedChains:
        if echain not in chains:
            return None
    AntigenChains = []
    for chain in chains:
        if chain not in ExpectedChains:
            AntigenChains.append(chain)

    return AntigenChains


def GetStructureFasta(structure):
    Residues = structure.get_residues()

    ResidueMap = IUPACD.protein_letters_3to1

    ResidueIdentifiers = [
        residue.resname.title()
        for residue in Residues
    ]

    Letters = [
        ResidueMap[resid]
        for resid in ResidueIdentifiers
        if resid in ResidueMap.keys()
    ]

    return "".join(Letters)


def get_structure_fasta_and_indexes(structure) -> List[Tuple[int, str]]:
    Residues = list(structure.get_residues())

    ResidueMap = IUPACD.protein_letters_3to1

    ResidueIdentifiers = [
        residue.resname.title()
        for residue in Residues
    ]

    ResidueIndexes = [
        residue.get_id()[1]
        for residue in Residues
    ]

    Letters = [
        ResidueMap[resid]
        for resid in ResidueIdentifiers
        if resid in ResidueMap.keys()
    ]

    return list(zip(ResidueIndexes, Letters))


def GetPartnersFromChains(chains):

    AntigenChains = GetAntigenChains(chains)

    gs = ["".join(group) for group in [ExpectedChains, AntigenChains]]
    return "_".join(gs)


def loadPDB(filepath, structure_name="X"):
    parser = PDB.PDBParser()
    return parser.get_structure(structure_name, filepath)


def getPDBChains(fpath, structure_name="X") -> List[str]:

    structure = loadPDB(fpath, structure_name)

    print(structure.id)
    print(structure.child_dict)
    chain_ids = []

    for chain in structure.get_chains():
        com = center_of_mass(chain)
        gcom = center_of_mass(chain, geometric=True)
        print(f"Chain {chain.id} - COM: {com}  GCOM: {gcom}")

        for xchain in structure.get_chains():

            xcom = center_of_mass(xchain)
            xgcom = center_of_mass(xchain, geometric=True)
            print(distance(com, xcom))
            print(distance(gcom, xgcom))

        chain_ids.append(chain.id)

        residues = list(chain.get_residues())
        print(len(residues))
        if not residues:
            exit(1)

    return chain_ids


def MergePDB(pdbinputpaths, pdboutputpath):
    parser = PDB.PDBParser()

    output_struct = parser.get_structure("out", pdbinputpaths[0])
    output_model = output_struct[0]

    for i, pdbinputpath in enumerate(pdbinputpaths[1:]):
        structure = parser.get_structure(str(i), pdbinputpath)

        for model in structure.get_list():
            for chain in structure.get_chains():
                output_model.add(chain.copy())

    out = PDB.PDBIO()
    out.set_structure(output_struct)

    out.save(pdboutputpath, preserve_atom_numbering=False)


def distance(a, b):
    w = [(a[x] - b[x]) ** 2 for x in range(3)]
    s = math.sqrt(sum(w))
    return s


def MoveChains(pdbfile, movchain, tgtchain, target_distance=30):

    parser = PDB.PDBParser()
    structure = parser.get_structure("X", pdbfile)

    structs = [
        [chain for chain in structure.get_chains() if chain.id == x][0]
        for x in [movchain, tgtchain]
    ]

    travelDistance, offset =\
        OffsetCoordinateToBringNearer(*structs, target_distance)

    print("Moving chains:")
    print(travelDistance)
    print(offset)

    for atom in structs[0].get_atoms():
        coordinate = atom.get_coord()
        atom.set_coord([coordinate[i] - offset[i] for i in range(3)])

    w = PDB.PDBIO()
    w.set_structure(structure)
    w.save(pdbfile)


def OffsetCoordinateToBringNearer(a, b, target_distance=18):
    centers_of_mass = [center_of_mass(w) for w in [a, b]]
    d = distance(*centers_of_mass)

    print("Original distance:")
    print(d)

    def get_directions(coms, i):
        w = coms[0][i] - coms[1][i]
        return w / abs(w)

    pointDirections = [get_directions(centers_of_mass, x) for x in range(3)]

    travelDistance = abs(target_distance - d)

    STEP = 0.2
    for i in range(1, int(1e3)):
        W = i * STEP
        wcoord = [W] * 3
        dw = distance([0] * 3, wcoord)

        if dw > travelDistance:
            break
        final_wcoord = wcoord

    final_offset = [(final_wcoord[i] * pointDirections[i]) for i in range(3)]
    return travelDistance, final_offset


def center_of_mass(entity, geometric=False):
    """
    Returns gravitic [default] or geometric center of mass of an Entity.
    Geometric assumes all masses are equal (geometric=True)
    """
    print(entity)
    # Structure, Model, Chain, Residue
    if isinstance(entity, PDB.Entity.Entity):
        atom_list = entity.get_atoms()
    # List of Atoms
    elif hasattr(entity, '__iter__') and [x for x in entity if x.level == 'A']:
        atom_list = entity
    else:
        # Some other weirdo object
        raise ValueError(
            "Center of Mass can only be calculated "
            "from the following objects:\n"
            "Structure, Model, Chain, Residue, list of Atoms.")

    masses = []
    positions = [[], [], []]
    # [ [X1, X2, ..] , [Y1, Y2, ...] , [Z1, Z2, ...] ]

    for atom in atom_list:
        masses.append(atom.mass)
        for i, coord in enumerate(atom.coord.tolist()):
            positions[i].append(coord)

    # If there is a single atom with undefined mass complain loudly.
    if 'ukn' in set(masses) and not geometric:
        raise ValueError(
            "Some Atoms don't have an element assigned.\n"
            "Try adding them manually or calculate the "
            "geometrical center of mass instead.")

    if geometric:
        return [sum(coord_list)/len(masses) for coord_list in positions]
    else:
        w_pos = [[], [], []]
        for atom_index, atom_mass in enumerate(masses):
            w_pos[0].append(positions[0][atom_index]*atom_mass)
            w_pos[1].append(positions[1][atom_index]*atom_mass)
            w_pos[2].append(positions[2][atom_index]*atom_mass)

        return [
            sum(coord_list) / sum(masses)
            for coord_list in w_pos
        ]


def GetResidueIndex(Residue):
    return Residue.get_full_id()[-1][1]


def GetAtomName(Atom):
    ID = Atom.get_full_id()[-1][0]
    Atoms = "HCOFSN"
    for A in Atoms:
        if A in ID:
            return A


def RemoveSolvent(inputpdbpath, outputpdbpath):
    parser = PDB.PDBParser()

    class NotSolvent(PDB.Select):
        def accept_residue(self, res):
            return "SOL" not in res.get_resname()

    output_struct = parser.get_structure("out", inputpdbpath)

    output = PDB.PDBIO()
    output.set_structure(output_struct)
    output.save(outputpdbpath, select=NotSolvent())
