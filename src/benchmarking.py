"""
.. module:: Benchmark
   :synopsis: This module implements functions that launch external programs for benchmarking.
"""

# Third-party modules
import os
import subprocess
import re
import operator
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
    logging.info(subprocess.Popen(["./bin/parMatt", pdb_file_1, pdb_file_2, "-o", "parMatt.out"],
                                  stdout=subprocess.PIPE).communicate()[0].decode("UTF-8").split("\n"))
    # Split the result file of parMatt into 2 different PDBs
    os.system("bin/splitMATT2chains.sh parMatt.out.pdb parMatt.out.txt tmp")
    # Calculate the TM-score between the two aligned structures by parMatt
    out_files = os.listdir("tmp")
    PARMATT_RES = subprocess.Popen(["./bin/TMscore", "tmp/"+out_files[0], "tmp/"+out_files[1]],
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


def plot_benchmark(results1, results2, gdt_results1, gdt_results2, tm_align_score, parmatt_score, prot1, prot2):
    """
        Plot the peeling level and number of PUs according to the TM-score
    """
    # TM-align alone results
    x1 = list(results1.keys()) if len(list(results1.keys())) > len(list(results2.keys())) else list(results2.keys())
    y1 = [tm_align_score] * len(list(gdt_results1.keys()))
    # parMatt result
    x2 = list(results1.keys()) if len(list(results1.keys())) > len(list(results2.keys())) else list(results2.keys())
    y2 = [parmatt_score] * len(list(gdt_results1.keys()))

    fig, ax = plt.subplots()
    ax.plot(x1, y1, c='green', label="TM-align (No Peeling)")
    ax.plot(x2, y2, c='black', label="parMATT")
    ax.scatter(list(results1.keys()), list(results1.values()), linestyle='solid', marker="*",
               c='red', lw=1, alpha=0.8, label=prot1+" vs "+prot2+" (Peeling)")
    ax.plot(list(results2.keys()), list(results2.values()), linestyle='solid', marker="*",
            c='blue', lw=1, alpha=0.8, label=prot2+" vs "+prot1+" (Peeling)")
    ax.scatter(list(gdt_results1.keys()), list(gdt_results1.values()), linestyle='solid', marker="*",
               c='magenta', lw=1, alpha=0.8, label=prot1+" vs "+prot2+" (Peeling + gdt)")
    ax.plot(list(gdt_results2.keys()), list(gdt_results2.values()), linestyle='solid', marker="*",
            c='cyan', lw=1, alpha=0.8, label=prot2+" vs "+prot1+" (Peeling + gdt)")
    ax.invert_yaxis()
    ax.legend(loc='best', frameon=False)
    ax.set_title('Benchmark of flexible alignment between '+prot1+" and "+prot2)
    ax.set_ylabel('TM-score')
    ax.set_xlabel('Number of Protein Units')
    plt.xticks(list(results1.keys()))

    # Add a second axis below the first one to show the peeling level
    ax2 = ax.twiny()
    # Add some extra space for the second axis at the bottom
    fig.subplots_adjust(bottom=0.3)
    ax2.set_xticks(list(results1.keys()))
    ax2.set_xticklabels(range(1, len(list(results1.keys()))+1))
    # Move twinned axis ticks and label from top to bottom
    ax2.xaxis.set_ticks_position('bottom')
    ax2.xaxis.set_label_position('bottom')
    # Offset the twin axis below the host
    ax2.spines["bottom"].set_position(("axes", -0.2))
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xlabel("Peeling Level")
    fig.savefig("results/benchmark_results.png", bbox_inches='tight', dpi=150)
    plt.show()
