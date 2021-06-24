
from pymol import cmd

import Bio.Data.IUPACData as IUPACD

def MutatePymol(input_structure, mutation_list, output_file):
    cmd.load(input_structure)

    cmd.wizard("mutagenesis")
    cmd.do("refresh_wizard")

    for (pos, fromAA, toAA) in mutation_list:

        toAA3lcode = IUPACD.protein_letters_1to3[toAA].upper()
        cmd.get_wizard().set_mode(toAA3lcode)
        cmd.get_wizard().do_select(f"{pos}/")

        cmd.frame(1)
        cmd.get_wizard().apply()

    cmd.save(output_file)
    cmd.delete('all')
