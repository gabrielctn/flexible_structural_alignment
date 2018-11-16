#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        ./fold_u PDB_FILE_1 PDB_FILE_2

    Arguments:
        PDB_FILE_1         Path to the first PDB file to align
        PDB_FILE_2         Path to the second PDB file to align

    Options:
        -h, --help        Show this
"""


# Third-party modules
from Bio.PDB.PDBParser import PDBParser
from datetime import datetime
from docopt import docopt
from schema import Schema, And, Use, SchemaError
import subprocess
import os
import re


def check_args():
    """
        Checks and validates the types of inputs parsed by docopt from command line.
    """
    schema = Schema({
        'PDB_FILE_1': Use(open, error='PDB_FILE_1 should be readable'),
        'PDB_FILE_2': Use(open, error='PDB_FILE_2 should be readable')
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

    PDB_FILE_1 = ARGUMENTS["PDB_FILE_1"]
    PDB_FILE_2 = ARGUMENTS["PDB_FILE_2"]
    PDB_NAME_1 = os.path.splitext(os.path.basename(PDB_FILE_1))[0]
    PDB_NAME_2 = os.path.splitext(os.path.basename(PDB_FILE_2))[0]
    DSSP_FILE_1 = "./data/" + PDB_NAME_1 + ".dss"
    DSSP_FILE_2 = "./data/" + PDB_NAME_2 + ".dss"

    ### Launch TMalign
    ##################
    tm_align_res = subprocess.Popen(["./bin/TMalign", PDB_FILE_1, PDB_FILE_2],
                     stdout=subprocess.PIPE).communicate()[0].decode("UTF-8").split("\n")

    # for line in tm_align_res.split("\n"):
    #     print(line)

    tm_score_regex = re.compile("^TM-score= (\\d+\\.\\d+).*$")
    aligned_len_regex = re.compile("^Aligned length=\\s*(\\d+).*$")

    for line in tm_align_res:
        a = re.search(tm_score_regex, line)
        b = re.search(aligned_len_regex, line)
        if a:
            tm_score = float(a.group(1))
        if b:
            aligned_len = int(b.group(1))
    print(tm_score, aligned_len)

    ### Parse PDB
    #############
    structure = PDBParser(QUIET=True).get_structure(PDB_NAME_1, PDB_FILE_1)

    ### Launch DSSP
    ###############
    subprocess.Popen(["./bin/mkdssp", "-i", PDB_FILE_1, "-o", DSSP_FILE_1],
                     stdout=subprocess.PIPE).communicate()[0]

    ### Launch Peeling
    ##################
    results = subprocess.Popen(["./bin/peeling11_4.1", "-pdb", PDB_FILE_1, "-dssp", DSSP_FILE_1, "-R2", "98",
                      "-ss2", "8", "-lspu", "20", "-mspu", "0", "-d0", "6.0", "-delta", "1.5", "-oss", "0",
                      "-p", "0", "-cp", "0", "-npu", "16"], stdout=subprocess.PIPE).communicate()[0].decode("UTF-8")
    #print(results)











    print("\nTotal runtime: {} seconds".format(str(datetime.now() - START_TIME)))
