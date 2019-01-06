"""
.. module:: Benchmark
   :synopsis: This module implements functions that launch external programs for benchmarking.
"""

# Third-party modules
import os
import subprocess
import re
import logging
import matplotlib.pyplot as plt
from multiprocessing import cpu_count

# Local modules
import src.parse as parse
import src.utils as utils


# Set logging info
logging.basicConfig(filename='flex.log', filemode='w',
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S')


def bench_tmalign(pdb_file_1, pdb_file_2):
    """
        Launch TM-align between the two given proteins to compare TM-scores with peeling method.
        The function prints directly in the terminal the output of TM-align.

        Args:
            pdb_file_1: PDB file of protein 1
            pdb_file_2: PDB file of protein 2
    """
    # Launch TM-align
    TM_ALIGN_RES = subprocess.Popen(["./bin/TMalign", pdb_file_1, pdb_file_2],
                                    stdout=subprocess.PIPE).communicate()[0]\
        .decode("UTF-8").split("\n")

    # Parse TM-align results
    TM_SCORE1, TM_SCORE2, RMSD, ALIGNED_LEN = parse.parse_tm_align(TM_ALIGN_RES)

    # Print TM-align results
    print("\n\n"+"#"*22+" BENCHMARKING "+"#"*22)
    print("\n\t\t    TMalign Results\n\t\t    "+"*"*15)
    print("\n\tTM-Score 1     : " + str(TM_SCORE1))
    print("\tTM-Score 2     : " + str(TM_SCORE2))
    print("\tRMSD           : " + str(RMSD))
    print("\tAligned length : " + str(ALIGNED_LEN))
    print("\n")

    return TM_SCORE2


def bench_parmatt(pdb_file_1, pdb_file_2, pdb_name_1, pdb_name_2):
    """
        parMATT is a parallelized version of MATT, which does local flexible structural alignments.
    """
    # Launch parMatt
    logging.info(subprocess.Popen(["./bin/parMatt", pdb_file_1, pdb_file_2, "-t", str(cpu_count()),
                      "-o", "parMatt.out"], stdout=subprocess.PIPE).communicate()[0].decode("UTF-8").split("\n"))
    # Split the result file of parMatt into 2 different PDBs
    os.system("bin/splitMATT2chains.sh parMatt.out.pdb parMatt.out.txt tmp")
    # Calculate the TM-score between the two aligned structures by parMatt
    PARMATT_RES = subprocess.Popen(["./bin/TMscore", "tmp/"+pdb_name_1+"_A.pdb", "tmp/"+pdb_name_2+"_A.pdb"],
                                   stdout=subprocess.PIPE).communicate()[0].decode("UTF-8").split("\n")
    PARMATT_SCORE = parse.parse_tm_score(PARMATT_RES)
    utils.clean_directory(dir="tmp")
    utils.clean_files(dir=".", pattern="^parMatt.*")
    print("\n\t\t    parMATT Results\n\t\t    "+"*"*15)
    print("\n\tTM-Score     : " + str(PARMATT_SCORE))
    print("\n")
    return PARMATT_SCORE


def gdt(aligned_all_pu_pdb, ref_protein):
    """
        The program gdt.pl determines the pairs of residues of each chain maximizing the TM-score.
        It returns several measures as well.

        Args:
            aligned_all_pu_pdb (str): PDB containing all the aligned PUs
            ref_protein (str): PDB of the protein against which the PUs are aligned (protein 2)

        Returns:
            (float): The maximum TM-score determined by GDT, for the alignment between the PUs and
            the protein 2
    """
    # gdt program needs a file containing the best aligned PUs together with the reference protein.
    with open("gdt.pdb", "w") as f_out:
        with open(aligned_all_pu_pdb) as f_in, open(ref_protein) as f_in2:
            for line in f_in:
                f_out.write(line)
            f_out.write("TER\n")
            for line in f_in2:
                f_out.write(line)

    GDT_RES = subprocess.Popen(["./bin/gdt.pl", "gdt.pdb"],
                               stdout=subprocess.PIPE).communicate()[0].decode("UTF-8").split("\n")
    utils.clean_files(dir=".", pattern="gdt.pdb")
    # Longest protein is the 2nd so we take the score normalized by chain 2
    regex_gdt = re.compile("^TM-score\s+:\s+(\d\.\d+).*Chain 2\)$")
    gdt_tm_score = -1
    for line in GDT_RES:
        score_found = re.search(regex_gdt, line)
        if score_found:
            gdt_tm_score = score_found.group(1)
    return gdt_tm_score


def plot_benchmark(results1, results2, parmatt_score, tm_align_score, prot1, prot2, gdt):
    """
        Plot the number of PU according to the TM-score

        Args:
            gdt (bool): True or False, if the results are with gdt tm-scores
    """
    # sorted by key, returns a list of tuples
    dict_to_tuples_1 = sorted(results1.items())
    dict_to_tuples_2 = sorted(results2.items())
    # unpack a list of pairs into two tuples: (nb_PUs), (tm-scores)
    x1, y1 = zip(*dict_to_tuples_1)
    x2, y2 = zip(*dict_to_tuples_2)
    # TM-align alone results
    x3 = x1 if len(x1) > len(x2) else x2
    y3 = [tm_align_score] * len(x3)
    # parMatt result
    x4 = x1 if len(x1) > len(x2) else x2
    y4 = [parmatt_score] * len(x3)

    fig, ax = plt.subplots()
    if gdt:
        ax.plot(x1, y1, linestyle='solid', marker="*", c='magenta', lw=1, alpha=0.8, label=prot1+" vs "+prot2+" (Peeling + gdt)")
        ax.plot(x2, y2, linestyle='solid', marker="*", c='cyan', lw=1, alpha=0.8, label=prot2+" vs "+prot1+" (Peeling + gdt)")
        ax.invert_yaxis()
    else:
        ax.plot(x1, y1, linestyle='solid', marker="*", c='red', lw=1, alpha=0.8, label=prot1+" vs "+prot2+" (Peeling)")
        ax.plot(x2, y2, linestyle='solid', marker="*", c='blue', lw=1, alpha=0.8, label=prot2+" vs "+prot1+" (Peeling)")
    ax.plot(x3, y3, c='green', label="TM-align (No Peeling)")
    ax.plot(x4, y4, c='black', label="parMATT")
    ax.legend(loc='best', frameon=False)
    ax.set_title('Benchmark of flexible alignment between '+prot1+" and "+prot2)
    ax.set_ylabel('TM-score')
    ax.set_xlabel('Number of Protein Units')
    plt.xticks(x3)

    # Add a second axis below the first one to show the peeling level
    ax2 = ax.twiny()
    # Add some extra space for the second axis at the bottom
    fig.subplots_adjust(bottom=0.3)
    ax2.set_xticks(x3)
    ax2.set_xticklabels(range(1, len(x3)+1))
    # Move twinned axis ticks and label from top to bottom
    ax2.xaxis.set_ticks_position('bottom')
    ax2.xaxis.set_label_position('bottom')
    # Offset the twin axis below the host
    ax2.spines["bottom"].set_position(("axes", -0.2))
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xlabel("Peeling Level")
    if gdt:
        fig.savefig("results/benchmark_results_with_gdt.png", bbox_inches='tight', dpi=150)
    else:
        fig.savefig("results/benchmark_results.png", bbox_inches='tight', dpi=150)
    plt.show()
