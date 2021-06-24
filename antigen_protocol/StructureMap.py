#!/bin/python

import os
import biostructmap

from .StructureMutator import parse_arguments
biostructmap.seqtools.LOCAL_BLAST = False
biostructmap.seqtools.LOCAL_EXONERATE = False


def GenerateStructure(pdbfile, AlignmentFile, workingDirectory, mutationPositions=None):
    file_prefix = os.path.split(os.path.splitext(pdbfile)[0])[-1]
    # Initialise structure object
    structure = biostructmap.Structure(pdbfile, 'test_pdb_name')

    # The location of known polymorphisms relative to the PDB sequence (we are not
    # providing a reference sequence for this example), for each chain.
    if mutationPositions is not None:
        data = {('F',): mutationPositions}

        results = structure.map(
            data,
            method='snps',
            ref=None,
            radius=15
        )
    else:
        msa_data = biostructmap.SequenceAlignment(AlignmentFile,
                                                  file_format="clustal")
        data = {('F',): msa_data}
        # Map polymorphism data using a radius of 15 Angstrom. Results are returned
        # in a new object.
        reference_seq = {'F': str(msa_data[0].seq)}
        results = structure.map(
            data,
            method='tajimasd',
            ref=reference_seq,
            radius=15,
            map_to_dna=True
        )


    # Use the results object to write data to a local PDB file, with data saved
    # in the B-factor column
    output_pdb = os.path.join(workingDirectory, file_prefix + "aln.pdb")
    results.write_data_to_pdb_b_factor(fileobj=output_pdb)


def main():
    options = parse_arguments()

    if not os.path.isdir(options.WorkingDirectory):
        os.mkdir(options.WorkingDirectory)

    GenerateStructure(
        options.PDBFile,
        options.AlignmentFile,
        options.WorkingDirectory
    )


if __name__ == "__main__":
    main()
