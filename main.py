#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        ./fold_u PDB_FILE

    Arguments:
        PDB_FILE                       PDB file

    Options:
        -h, --help                     Show this
"""


# Third-party modules
from Bio.PDB.PDBParser import PDBParser
from datetime import datetime
from docopt import docopt
from schema import Schema, And, Use, SchemaError
import contextlib
import subprocess


def check_args():
    """
        Checks and validates the types of inputs parsed by docopt from command line.
    """
    schema = Schema({
        'PDB_FILE': Use(open, error='PDB_FILE should be readable')
        })
    try:
        schema.validate(ARGUMENTS)
    except SchemaError as err:
        exit(err)




if __name__ == "__main__":

    START_TIME = datetime.now()

    ### Parse command line
    ######################
    ARGUMENTS = docopt(__doc__, version='Protein Flexible Alignment Tool 1.0')
    # Check the types and ranges of the command line arguments parsed by docopt
    check_args()

    PDB_FILE = ARGUMENTS["PDB_FILE"]
    PDB_NAME = PDB_FILE.split(".")[0]

    with contextlib.redirect_stdout(None):
        p = PDBParser()
        structure = p.get_structure(PDB_NAME, PDB_FILE)
        parser = PDBParser(PERMISSIVE=1)
        model = structure[0]
    dssp = subprocess.Popen(["./bin/mkdssp", "-i", PDB_FILE, "-o", "out.dssp"], stdout=subprocess.PIPE).communicate()[0]
    print(dssp)

    print("\nTotal runtime: {} seconds".format(str(datetime.now() - START_TIME)))
