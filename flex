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
import operator
import copy
from datetime import datetime
from Bio.PDB.PDBParser import PDBParser
from docopt import docopt

# Local modules
import src.parse as parse
import src.pdb as pdb
import src.utils as utils


def flex_align(PDB_FILE_1, PDB_FILE_2, PDB_NAME_1, PDB_NAME_2, DSSP_FILE_1, DSSP_FILE_2):
    """
    """
    ### Launch TMalign
    ##################
    # TM_ALIGN_RES = subprocess.Popen(["./bin/TMalign", PDB_FILE_1, PDB_FILE_2],
    #                                 stdout=subprocess.PIPE).communicate()[0]\
    #     .decode("UTF-8").split("\n")
    #
    # ### Parse TMalign
    # #################
    # TM_SCORE1, TM_SCORE2, RMSD, ALIGNED_LEN = parse.parse_tm_align(TM_ALIGN_RES)

    # print("\n### TMalign Results ###")
    # print("\nTM-Score 1     : " + str(TM_SCORE1))
    # print("TM-Score 2     : " + str(TM_SCORE2))
    # print("RMSD           : " + str(RMSD))
    # print("Aligned length : " + str(ALIGNED_LEN))
    # print("\n\n\n")

    ### Parse PDB
    #############

    # The PDB files are reindexed to fit the Protein Peeling program and keep only chain A
    REINDEXED_PDB_1 = pdb.reindex_pdb(1, PDB_FILE_1, True)
    REINDEXED_PDB_2 = pdb.reindex_pdb(1, PDB_FILE_2, True)
    # Write the new reindexed PDBs
    NEW_PDB_NAME_1 = PDB_NAME_1 + "_new.pdb"
    NEW_PDB_NAME_2 = PDB_NAME_2 + "_new.pdb"
    NEW_PDB_FILE_1 = "data/" + NEW_PDB_NAME_1
    NEW_PDB_FILE_2 = "data/" + NEW_PDB_NAME_2
    with open(NEW_PDB_FILE_1, "w") as f_out_1, open(NEW_PDB_FILE_2, "w") as f_out_2:
        f_out_1.write(REINDEXED_PDB_1)
        f_out_2.write(REINDEXED_PDB_2)
    # Parse the reindexed PDBs
    STRUCTURE_1 = PDBParser(QUIET=True).get_structure(NEW_PDB_NAME_1, NEW_PDB_FILE_1)
    STRUCTURE_2 = PDBParser(QUIET=True).get_structure(NEW_PDB_NAME_1, NEW_PDB_FILE_1)

    ### Launch DSSP
    # The DSSP file is required for Protein Peeling program
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
    # Clean workspace
    utils.clean_files(dir="data/", pattern="*.dss")
    utils.clean_files(dir=".", pattern="file_*")

    ### Parse Peeling
    #################
    PEELING_DICT = parse.parse_protein_peeling(PEELING_RES)

    ### Flexible alignment
    ######################
    SCORES = {}
    GENERAL_RESULTS = {}
    aligned_pu_pdb = "results/aligned_pu.pdb"
    aligned_max_pu = ""

    # For every peeling level, every PU of the 1st protein are aligned against the 2nd protein iteratively.
    # The PU with the best TM-score obtained is kept. Then the portion of this best alignment is
    # removed in order not to bias next alignments. Then other PUs are aligned, etc.

    # For every peeling level
    for peeling_level, all_pu in enumerate(PEELING_DICT["PU_BOUNDS"]):
        cropped_pdb = NEW_PDB_FILE_2
        all_pu_ref = copy.deepcopy(all_pu)
        # While there are still PUs to align for this peeling level
        while len(all_pu) > 1:
            # Align each PU to the protein
            for pu_index, (start, end) in enumerate(all_pu):
                # File name = PDB_NAME_1 + peeling_level + "pu" + PU_index + ".pdb"
                pu_pdb = "results/{}_{}_pu_{}.pdb".format(PDB_NAME_1, peeling_level+1, pu_index+1)
                # Write the PDB file for the actual PU
                pdb.write_pu_pdb(STRUCTURE_1, "A", int(start), int(end), pu_pdb)
                tm_aligned_file = "tmp/TMaligned_"+str(peeling_level)+"_"+str(pu_index)+".sup"
                # Calculate TMscore of the PU aligned to the static PDB chain
                #print(pu_pdb, cropped_pdb)
                TM_ALIGN_RES = subprocess.Popen(["./bin/TMalign", pu_pdb, cropped_pdb,
                                                 "-o", tm_aligned_file],
                                                stdout=subprocess.PIPE).communicate()[0]\
                                                .decode("UTF-8").split("\n")
                TM_SCORE1, TM_SCORE2, RMSD, ALIGNED_LEN = parse.parse_tm_align(TM_ALIGN_RES)
                SCORES[(peeling_level, pu_index)] = TM_SCORE2
            # Get the (peeling_level, pu_index) associated to the best TM-Score
            key_max_value = max(SCORES.items(), key=operator.itemgetter(1))[0]
            # Get the boundaries of the best PU
            if len(all_pu) == 1:
                best_pu_start = int(all_pu[0][0])
                best_pu_end = int(all_pu[0][1])
            else:
                best_pu_start = int(all_pu[key_max_value[1]][0])
                best_pu_end = int(all_pu[key_max_value[1]][1])
            # Save the coordinates of the best aligned PU
            with open(tm_aligned_file + "_all_atm", "r") as f_in:
                line = f_in.readline()
                # Skip the PDB header, go to the coordinates of the aligned PU
                while line[0:4].strip() != "ATOM":
                    line = f_in.readline()
                # Save the coordinates of the aligned PU
                aligned_max_pu = "tmp/{}_aln_max_pu_{}.pdb".format(peeling_level+1, key_max_value[1])
                with open(aligned_max_pu, "a") as f_out:
                    while line[0:4].strip() != "TER":
                        f_out.write(line)
                        line = f_in.readline()

            # Remove the aligned PU from the static PDB before continuing to align other PUs
            cropped_pdb = "tmp/"+PDB_NAME_2+"_new"+"_"+str(key_max_value[0]+1)+"_pu_"+str(key_max_value[1]+1)+".pdb"
            pdb.remove_aligned_pu_from_pdb(STRUCTURE_2, cropped_pdb, best_pu_start, best_pu_end)
            # Remove the best aligned PU to align the others afterwards
            all_pu.remove((str(best_pu_start), str(best_pu_end)))
        # Write the full PDB of the best aligned PUs by concatenating them in the right order
        aln_max_pu_file = "tmp/{}_aln_max_pu_{}.pdb".format(peeling_level+1, pu_index+1)
        with open(aligned_pu_pdb, "a") as f_out:
            print(all_pu_ref)
            for pu_index in range(1, len(all_pu_ref) + 1):
                aln_max_pu_file = "tmp/{}_aln_max_pu_{}.pdb".format(peeling_level+1, pu_index)
                with open(aln_max_pu_file, "r") as f_in:
                    for line in f_in:
                        f_out.write(line)
        # Launch TMscore between the protein 1 (aligned PUs) and protein 2
        TM_SCORE_RES = subprocess.Popen(["./bin/TMscore", NEW_PDB_FILE_2, aligned_pu_pdb],
                                        stdout=subprocess.PIPE).communicate()[0].decode("UTF-8").split("\n")
        TM_SCORE = parse.parse_tm_score(TM_SCORE_RES)
        GENERAL_RESULTS[PEELING_DICT["NB_PU"][peeling_level]] = TM_SCORE
        print("Number of PUs: {}\t\tTM-Score: {}".format(PEELING_DICT["NB_PU"][peeling_level], TM_SCORE))
        # Clean directories
        utils.clean_files(dir="results", pattern="aligned_pu.pdb")
        utils.clean_files(dir="tmp/", pattern="TM*")

if __name__ == "__main__":

    START_TIME = datetime.now()

    ### Parse command line
    ######################
    ARGUMENTS = docopt(__doc__, version='Protein Flexible Alignment Tool 1.0')
    # Check the types and ranges of the command line arguments parsed by docopt
    utils.check_args(ARGUMENTS)

    ### Filenames and paths
    #######################
    PDB_FILE_1 = ARGUMENTS["PDB_FILE_1"]
    PDB_FILE_2 = ARGUMENTS["PDB_FILE_2"]
    PDB_NAME_1 = os.path.splitext(os.path.basename(PDB_FILE_1))[0]
    PDB_NAME_2 = os.path.splitext(os.path.basename(PDB_FILE_2))[0]
    DSSP_FILE_1 = "./data/" + PDB_NAME_1 + ".dss"
    DSSP_FILE_2 = "./data/" + PDB_NAME_2 + ".dss"

    flex_align(PDB_FILE_1, PDB_FILE_2, PDB_NAME_1, PDB_NAME_2, DSSP_FILE_1, DSSP_FILE_2)
    print("\n")
    #flex_align(PDB_FILE_2, PDB_FILE_1, PDB_NAME_2, PDB_NAME_1, DSSP_FILE_2, DSSP_FILE_1)

    print("\nTotal runtime: {} seconds".format(str(datetime.now() - START_TIME)))
