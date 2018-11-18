#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        ./flex PDB_FILE_1 PDB_FILE_2

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
from schema import Schema, Use, SchemaError
import subprocess
import os
import re

# Local modules
import src.parse as parse


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


def clean_peeling_outputs():
    """
    Removes unused files generated by Protein Peeling program.
    """
    pattern = "^file_.*$"
    for file in os.listdir("."):
        if re.search(pattern, file):
            os.remove(os.path.join(".", file))


if __name__ == "__main__":

    START_TIME = datetime.now()

    ### Parse command line
    ######################
    ARGUMENTS = docopt(__doc__, version='Protein Flexible Alignment Tool 1.0')
    # Check the types and ranges of the command line arguments parsed by docopt
    check_args()

    ### Definition of filenames and paths
    ##################
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

    ### Parse TMalign
    #################
    tm_score1, tm_score2, rmsd, aligned_len = parse.parse_tm_align(tm_align_res)

    ### Parse PDB
    #############
    # The PDB file is reindexed to fit the Protein Peeling program
    reindexed_pdb = parse.reindex_pdb(1, PDB_FILE_1, True)
    # Write the new reindexed PDB
    NEW_PDB_NAME_1 = PDB_NAME_1 + "_new.pdb"
    NEW_PDB_FILE_1 = "data/" + NEW_PDB_NAME_1
    with open(NEW_PDB_FILE_1, "w") as f_out:
        f_out.write(reindexed_pdb)
    structure = PDBParser(QUIET=True).get_structure(NEW_PDB_NAME_1, NEW_PDB_FILE_1)

    ### Launch DSSP
    ###############
    subprocess.Popen(["./bin/mkdssp", "-i", NEW_PDB_FILE_1, "-o", DSSP_FILE_1],
                     stdout=subprocess.PIPE).communicate()[0]

    ### Launch Peeling
    ##################
    peeling_res = subprocess.Popen(["./bin/peeling11_4.1", "-pdb", NEW_PDB_FILE_1, "-dssp", DSSP_FILE_1,
                                    "-R2", "98", "-ss2", "8", "-lspu", "20", "-mspu", "0", "-d0",
                                    "6.0", "-delta", "1.5", "-oss", "1", "-p", "0", "-cp", "0",
                                    "-npu", "16"],
                                   stdout=subprocess.PIPE).communicate()[0].decode("UTF-8").split("\n")
    clean_peeling_outputs()
    for i in peeling_res:
        print(i)

    ### Parse Peeling
    #################
    peeling_dict = parse.parse_protein_peeling(peeling_res)

    # Write all PU in different PDB files
    for i, iter in enumerate(peeling_dict["PU_bounds"]):
        for j, (start, end) in enumerate(iter):
            out_file = "results/" + PDB_NAME_1 + "_" + str(i+1) + "_pu_" + str(j+1) + ".pdb"
            parse.write_pdb_portion(structure, "A", int(start), int(end), out_file)

    print("\nTotal runtime: {} seconds".format(str(datetime.now() - START_TIME)))