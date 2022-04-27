from typing import List, Dict
import sys
import os
from Bio import SeqIO, SeqUtils

from autogromacs.mdanalysis import process_simulation_name
import pandas as pd

MutationSpecifiers = Dict[str, List[str]]


def is_specifier(fname: str,
                 fractionary: bool,
                 normal: bool) -> bool:
    if fname.startswith("DUMMY"):
        if fname.endswith(".pdb"):
            return False

        if "-" in fname:
            return fractionary

        return normal

    return False


def process_directory(
        ROOT_DIR: str) -> MutationSpecifiers:

    mutation_specifiers = [
        F
        for F in os.listdir(ROOT_DIR)
        if is_specifier(F, True, True)
    ]

    mutation_specifiers = sorted(mutation_specifiers)

    specifiers = {}
    for specifier in mutation_specifiers:
        fpath = os.path.join(ROOT_DIR, specifier)
        with open(fpath) as f:
            W = []
            muts = f.readlines()
            for mut in muts:
                W.append(mut.strip("\n"))
        specifiers[specifier] = W

    return specifiers


def process_file_table(
        mutation_specifiers: MutationSpecifiers,
        output_path: str):
    fn = ['Identificador', 'Substituições']

    data = []
    for mut, mutations in mutation_specifiers.items():
        data.append({
            fn[0]: process_simulation_name(mut),
            fn[1]: "; ".join(mutations)
        })

    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)


def load_from_recipe(recipe_dir):
    pass


def main():
    directories = sys.argv[1:]
    specifiers = {}
    for directory in directories:
        k = process_directory(directory)
        specifiers.update(k)

    process_file_table(specifiers, "muts.csv")
