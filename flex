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
import subprocess
import os
from datetime import datetime
from Bio.PDB.PDBParser import PDBParser
from docopt import docopt

# Local modules
import src.parse as parse
import src.pdb as pdb
import src.utils as utils


if __name__ == "__main__":

    START_TIME = datetime.now()

    ### Parse command line
    ######################
    ARGUMENTS = docopt(__doc__, version='Protein Flexible Alignment Tool 1.0')
    # Check the types and ranges of the command line arguments parsed by docopt
    utils.check_args(ARGUMENTS)

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
    TM_ALIGN_RES = subprocess.Popen(["./bin/TMalign", PDB_FILE_1, PDB_FILE_2],
                                    stdout=subprocess.PIPE).communicate()[0]\
                                    .decode("UTF-8").split("\n")

    ### Parse TMalign
    #################
    TM_SCORE1, TM_SCORE2, RMSD, ALIGNED_LEN = parse.parse_tm_align(TM_ALIGN_RES)

    ### Parse PDB
    #############
    # The PDB file is reindexed to fit the Protein Peeling program
    REINDEXED_PDB = pdb.reindex_pdb(1, PDB_FILE_1, True)
    # Write the new reindexed PDB
    NEW_PDB_NAME_1 = PDB_NAME_1 + "_new.pdb"
    NEW_PDB_FILE_1 = "data/" + NEW_PDB_NAME_1
    with open(NEW_PDB_FILE_1, "w") as f_out:
        f_out.write(REINDEXED_PDB)
    STRUCTURE = PDBParser(QUIET=True).get_structure(NEW_PDB_NAME_1, NEW_PDB_FILE_1)

    ### Launch DSSP
    ###############
    OUT, ERR = subprocess.Popen(["./bin/mkdssp", "-i", NEW_PDB_FILE_1, "-o", DSSP_FILE_1],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

    ### Launch Peeling
    ##################
    PEELING_PROCESS_RES = subprocess.Popen(["./bin/peeling11_4.1", "-pdb", NEW_PDB_FILE_1,
                                            "-dssp", DSSP_FILE_1, "-R2", "98", "-ss2", "8",
                                            "-lspu", "20", "-mspu", "0", "-d0", "6.0",
                                            "-delta", "1.5", "-oss", "1", "-p", "0", "-cp", "0",
                                            "-npu", "16"], stdout=subprocess.PIPE).communicate()[0]
    PEELING_RES = PEELING_PROCESS_RES.decode("UTF-8").split("\n")
    utils.clean_peeling_output_files()

    for i in PEELING_RES:
        print(i)

    ### Parse Peeling
    #################
    PEELING_DICT = parse.parse_protein_peeling(PEELING_RES)


    ### Flexible alignment
    ######################
    SCORES = {}


    # After each alignment we remove the aligned portion to avoid bias
    for peeling_level, all_pu in enumerate(PEELING_DICT["PU_BOUNDS"], 1):
        for pu_index, (start, end) in enumerate(all_pu, 1):
            # File name = PDB_NAME_1 + peeling_level + "pu" + PU_index + ".pdb"
            pu_pdb = "results/"+PDB_NAME_1+"_"+str(peeling_level)+"_pu_"+str(pu_index)+".pdb"
            # Write the PDB file for the actual PU
            pdb.write_pdb_portion(STRUCTURE, "A", int(start), int(end), pu_pdb)
            # Calculate TMscore of this PU aligned of the static PDB chain
            if pu_index == 1:
                TM_SCORE_RES = subprocess.Popen(["./bin/TMscore", pu_pdb, PDB_FILE_2],
                                                stdout=subprocess.PIPE).communicate()[0]\
                                                .decode("UTF-8").split("\n")
            else:
                TM_SCORE_RES = subprocess.Popen(["./bin/TMscore", pu_pdb, cropped_pdb],
                                                stdout=subprocess.PIPE).communicate()[0]\
                                                .decode("UTF-8").split("\n")
            SCORES[(peeling_level, pu_index)] = parse.parse_tm_score(TM_SCORE_RES)
            # Remove the aligned PU from the static PDB before continuing
            cropped_pdb = "tmp/"+PDB_NAME_1+"_new"+"_"+str(peeling_level)+"_pu_"+str(pu_index)+".pdb"
            pdb.extract_pu_from_pdb(STRUCTURE, cropped_pdb, int(start), int(end))


    for key, value in SCORES.items():
        print("Peeling Level: ", key[0], " PU index: ", key[1], "TM-score = ", value)

    print("\nTotal runtime: {} seconds".format(str(datetime.now() - START_TIME)))
