#!/bin/python


import argparse
from ..Wrapper import Rosetta


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", dest="StructureFile", required=True)
    parser.add_argument("--chain", dest="WantedChainName", required=False)
    return parser.parse_args()


def main():
    options = parse_arguments()
    Rosetta.ScorePdb(options.StructureFile)

    
if __name__ == "__main__":
    main()
