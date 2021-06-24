#!python

"""

Parse the csv output from IEDB.org's BepiPred

B-Cell epitope prediction method.

"""
import csv
import argparse

from .. OutputConstructor import plotEpitopeScore


def parse_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", dest="CsvFilePath", required=True)
    parser.add_argument("-o", dest="OutputFilepath")
    return parser.parse_args()


def main():
    options = parse_arguments()
    epitope_data = parse_bepipred(options.CsvFilePath)
    #return epitope_data
    plotEpitopeScore.plotEpitopeScores(epitope_data,
                                       output_filepath=options.OutputFilepath)


# FIXME: Change to pandas.
def parse_bepipred(fpath):
    with open(fpath) as f:
        reader = csv.reader(f)
        header = next(reader, None)

        records = []
        for record in reader:
            out_row = {}
            for r, col_name in enumerate(header):
                out_row[col_name] = record[r]
            records.append(out_row)

    return records


if __name__ == "__main__":
    main()
