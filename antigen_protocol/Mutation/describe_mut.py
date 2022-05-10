
from typing import Dict, List
import os
import csv
from .structure_name import process_simulation_name


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
        ROOT_DIR: str,
        load_fractionary: bool,
        load_normal: bool,
        output_path: str) -> None:

    mutation_specifiers = [
        F
        for F in os.listdir(ROOT_DIR)
        if is_specifier(F, load_fractionary, load_normal)
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

    process_file_tex(specifiers, output_path + ".tex")
    process_file_table(specifiers, output_path + ".csv")


K = Dict[str, List[str]]


def process_file_table(mutation_specifiers: K, output_path: str):
    with open(output_path, 'w') as f:
        fn = ['Identificador', 'Substituições']

        out = csv.DictWriter(f, fn)
        out.writeheader()
        for mut, mutations in mutation_specifiers.items():
            out.writerow({
                fn[0]: process_simulation_name(mut),
                fn[1]: "; ".join(mutations)
            })


def process_file_tex(mutation_specifiers: K, output_path: str):
    IO = "\\begin{itemize}"
    IE = "\\end{itemize}"

    output = [IO]

    for specifier_file, mutations in mutation_specifiers.items():
        Name = process_simulation_name(specifier_file)

        output.append("\\item " + Name.replace("#", "\\#"))
        output.append(IO)
        for mut in mutations:
            output.append("\\item " + mut)

        output.append(IE)
    output.append(IE)

    with open(output_path, 'w') as out:
        out.write("\n".join(output))


if __name__ == "__main__":
    process_directory("/home/Science/MD", True, True, "muts")
