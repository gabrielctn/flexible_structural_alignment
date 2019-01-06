"""
.. module:: Benchmark
   :synopsis: This module implements functions that launch external programs for benchmarking.
"""

# Third-party modules
import subprocess

# Local modules
import src.parse as parse



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
    print("\n\n"+"#"*13+" BENCHMARKING "+"#"*13)
    print("\n\t    TMalign Results\n\t    "+"*"*15)
    print("\nTM-Score 1     : " + str(TM_SCORE1))
    print("TM-Score 2     : " + str(TM_SCORE2))
    print("RMSD           : " + str(RMSD))
    print("Aligned length : " + str(ALIGNED_LEN))
    print("\n\n\n")
