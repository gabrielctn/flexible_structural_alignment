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
from docopt import docopt
from Bio.PDB.PDBParser import PDBParser

# Local modules
import src.parse as parse
import src.pdb as pdb
import src.utils as utils
import src.benchmarking as benchmarking


def flex_align(PDB_FILE_1, PDB_FILE_2, PDB_NAME_1, PDB_NAME_2, DSSP_FILE_1, DSSP_FILE_2):
    """
        Do the actual flexible alignment !
    """

    ### Clean and Parse PDBs
    ########################
    NEW_PDB_NAME_1, NEW_PDB_FILE_1 = pdb.clean_pdb(PDB_FILE_1, PDB_NAME_1)
    NEW_PDB_NAME_2, NEW_PDB_FILE_2 = pdb.clean_pdb(PDB_FILE_2, PDB_NAME_2)
    # Parse the reindexed PDBs
    STRUCTURE_1 = PDBParser(QUIET=True).get_structure(NEW_PDB_NAME_1, NEW_PDB_FILE_1)
    STRUCTURE_2 = PDBParser(QUIET=True).get_structure(NEW_PDB_NAME_2, NEW_PDB_FILE_2)

    ### Launch DSSP
    ###############
    # The DSSP file is required for Protein Peeling program
    OUT, ERR = subprocess.Popen(["./bin/mkdssp", "-i", NEW_PDB_FILE_1, "-o", DSSP_FILE_1],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

    ### Launch Peeling
    ##################
    PEELING_PROCESS_RES = subprocess.Popen(["./bin/peeling11_4.1", "-pdb", NEW_PDB_FILE_1,
                                            "-dssp", DSSP_FILE_1, "-R2", "98", "-ss2", "8",
                                            "-lspu", "20", "-mspu", "0", "-d0", "6.0",
                                            "-delta", "1.5", "-oss", "0", "-p", "0", "-cp", "0",
                                            "-npu", "16"], stdout=subprocess.PIPE).communicate()[0]
    PEELING_RES = PEELING_PROCESS_RES.decode("UTF-8").split("\n")
    # Clean workspace
    utils.clean_files(dir="data/", pattern="^.*\.dss$")
    utils.clean_files(dir=".", pattern="^file_.*$")

    ### Parse Peeling
    #################
    PEELING_DICT = parse.parse_protein_peeling(PEELING_RES)

    ### Flexible alignment
    ######################
    GENERAL_RESULTS = {}
    GDT = {}

    # Create folders if they don't exist already
    utils.create_dir_safely("tmp")
    utils.create_dir_safely("results")

    # For every peeling level, every PU of the 1st protein are aligned against the 2nd protein iteratively.
    # The PU with the best TM-score obtained is kept. Then the portion of this best alignment is
    # removed in order not to bias next alignments. Then other PUs are aligned, etc.

    # For every peeling level
    for peeling_level, all_pu in enumerate(PEELING_DICT["PU_BOUNDS"]):
        SCORES = {}
        aligned_pu_pdb = ("results/" + PDB_NAME_1 + "_" + PDB_NAME_2 + "_aligned_PUs_level_"
                          + str(peeling_level+1) + ".pdb")
        cropped_pdb = NEW_PDB_FILE_2
        all_pu_ref = copy.deepcopy(all_pu)
        # Write the PDBs of all the PUs of the actual peeling level
        pdb.write_pu_pdb(all_pu, PDB_NAME_1, peeling_level, STRUCTURE_1)
        # While there are still PUs to align for this peeling level
        while len(all_pu) >= 1:
            # Align each PU to the protein 2
            for pu_index, (start, end) in enumerate(all_pu):
                # File name = PDB_NAME_1 + peeling_level + "pu" + PU_index + ".pdb"
                pu_pdb = "results/{}_{}_pu_{}.pdb".format(PDB_NAME_1, peeling_level+1, pu_index+1)
                tm_aligned_file = "tmp/TMaligned_"+str(peeling_level)+"_"+str(pu_index+1)+".sup"
                # Calculate TMscore of the PU aligned to the static PDB chain
                TM_ALIGN_RES = subprocess.Popen(["./bin/TMalign", pu_pdb, cropped_pdb,
                                                 "-o", tm_aligned_file],
                                                stdout=subprocess.PIPE).communicate()[0]\
                    .decode("UTF-8").split("\n")
                TM_SCORE1, TM_SCORE2, RMSD, ALIGNED_LEN = parse.parse_tm_align(TM_ALIGN_RES)
                # Keep the TM-score normalized by the longest protein (PROT 2)
                SCORES[(peeling_level, start, end)] = TM_SCORE2

            # Get the PU (peeling_level, start, end) associated to the best TM-Score
            key_max_value = max(SCORES.items(), key=operator.itemgetter(1))[0]
            best_pu_peel_level = key_max_value[0]
            # Get the boundaries of the best PU
            best_pu_start = int(key_max_value[1])
            best_pu_end = int(key_max_value[2])
            best_pu_index = all_pu_ref.index((str(best_pu_start), str(best_pu_end)))
            # Save the coordinates of the PU best ligned with TM-align
            pdb.save_best_aligned_pu(best_pu_peel_level, best_pu_index)
            # Remove the aligned PU from the static PDB before continuing to align other PUs
            cropped_pdb = ("tmp/"+PDB_NAME_2+"_cropped"+"_"+str(best_pu_peel_level+1)+"_pu_"
                           + str(best_pu_index+1)+".pdb")
            pdb.remove_aligned_pu_from_pdb(STRUCTURE_2, cropped_pdb, best_pu_start, best_pu_end)
            # Remove the best aligned PU to align the others afterwards
            all_pu.remove((str(best_pu_start), str(best_pu_end)))
            del SCORES[key_max_value]
        pdb.concatenate_best_aligned_pus(peeling_level, aligned_pu_pdb, all_pu_ref)
        # Launch TMscore between the protein 1 (aligned PUs) and protein 2
        TM_SCORE_RES = subprocess.Popen(["./bin/TMscore", aligned_pu_pdb, NEW_PDB_FILE_2],
                                        stdout=subprocess.PIPE).communicate()[0].decode("UTF-8")\
            .split("\n")
        TM_SCORE = parse.parse_tm_score(TM_SCORE_RES)
        # Save general results for further benchmarking
        GENERAL_RESULTS[PEELING_DICT["NB_PU"][peeling_level]] = TM_SCORE
        # Optimize alignement between PUs and the reference protein with gdt perl program
        GDT_TM_SCORE = benchmarking.gdt(aligned_pu_pdb, NEW_PDB_FILE_2)
        GDT[PEELING_DICT["NB_PU"][peeling_level]] = GDT_TM_SCORE
        print("{:>13}{:>16}{:>11}{:>9}".format(peeling_level+1, PEELING_DICT["NB_PU"][peeling_level],
                                                TM_SCORE, GDT_TM_SCORE))
        # Clean directory for next alignments
        utils.clean_files(dir="tmp/", pattern="^TM.*$")
    return(GENERAL_RESULTS, GDT)


if __name__ == "__main__":

    START_TIME = datetime.now()

    ### Parse command line
    ######################

    ARGUMENTS = docopt(__doc__, version='Flexible Structural Alignment Tool 1.0')
    # Check the types and ranges of the command line arguments parsed by docopt
    utils.check_args(ARGUMENTS)

    ### Set filenames and paths
    ###########################

    # PDBs
    PDB_FILE_1 = ARGUMENTS["PDB_FILE_1"]
    PDB_FILE_2 = ARGUMENTS["PDB_FILE_2"]
    # PDBs Basenames
    PDB_NAME_1 = os.path.splitext(os.path.basename(PDB_FILE_1))[0]
    PDB_NAME_2 = os.path.splitext(os.path.basename(PDB_FILE_2))[0]
    # DSSP files names
    DSSP_FILE_1 = "./data/" + PDB_NAME_1 + ".dss"
    DSSP_FILE_2 = "./data/" + PDB_NAME_2 + ".dss"

    ### Launch Flexible alignement with Protein Peeling method
    ##########################################################

    print("\n\t Flexible Structural Alignment Tool 1.0\n\n")
    print("\t\t     {} vs {}\n\t\t     {}\n".format(PDB_NAME_1, PDB_NAME_2, 12*"*"))
    print("Peeling level\tNumber of PUs\tTM-Score    TM-score (gdt)")
    ### Launch flexible alignment: PROT1 vs PROT2
    (RESULTS_1, GDT_1) = flex_align(PDB_FILE_1, PDB_FILE_2, PDB_NAME_1, PDB_NAME_2, DSSP_FILE_1, DSSP_FILE_2)

    # Clean workspace
    utils.clean_files(dir="results", pattern="^((?!aligned).)*$")
    utils.clean_directory(dir="tmp")

    print("\n\n\t\t     {} vs {}\n\t\t     {}\n".format(PDB_NAME_2, PDB_NAME_1, 12*"*"))
    print("Peeling level\tNumber of PUs\tTM-Score    TM-score (gdt)")
    ### Launch flexible alignement: PROT2 vs PROT1
    (RESULTS_2, GDT_2) = flex_align(PDB_FILE_2, PDB_FILE_1, PDB_NAME_2, PDB_NAME_1, DSSP_FILE_2, DSSP_FILE_1)

    # Clean workspace
    utils.clean_files(dir="results", pattern="^((?!aligned).)*$")
    utils.clean_directory(dir="tmp")

    ### Launch TMalign for benchmarking
    ###################################
    TM_ALIGN_RES = benchmarking.bench_tmalign(PDB_FILE_1, PDB_FILE_2)

    ### Launch parMATT for benchmarking
    ###################################
    PARMATT_SCORE = benchmarking.bench_parmatt(PDB_FILE_1, PDB_FILE_2, PDB_NAME_1, PDB_NAME_2)

    # Benchmarking results (with and without gdt TM-scores)
    benchmarking.plot_benchmark(RESULTS_1, RESULTS_2, TM_ALIGN_RES, PARMATT_SCORE, PDB_NAME_1, PDB_NAME_2, False)
    benchmarking.plot_benchmark(GDT_1, GDT_2, TM_ALIGN_RES, PARMATT_SCORE, PDB_NAME_1, PDB_NAME_2, True)
    print("\nResult plots are stored in results/")
    # Display runtime
    print("\nTotal runtime: {} seconds".format(str(datetime.now() - START_TIME)))
