"""
.. module:: Parse
   :synopsis: This module implements parsing functions.
"""

# Third-party modules
import re


def parse_tm_align(tm_align_res):
    """
    This function parses TMalign's terminal output.
    It retrieves:
    1. TMscore for chain 1 (normalized by the length of chain 1)
    2. TMscore for chain 2 (normalized by the length of chain 2)
    3. RMSD
    4. Alignment length

    Args:
        tm_align_res (list): Terminal's output of TMalign

    Returns:
        Tuple: (tm_score1, tm_score2, rmsd, aligned_length)
    """
    tm_score_regex = re.compile("^TM-score= (\\d+\\.\\d+).*$")
    rmsd_regex = re.compile("^.*RMSD=\\s+(\\d+\\.\\d+).*$")
    aligned_len_regex = re.compile("^Aligned length=\\s*(\\d+).*$")

    tm_score1 = None
    tm_score2 = None
    rmsd = None
    aligned_len = None

    flag = False
    for line in tm_align_res:
        tm_score_found = re.search(tm_score_regex, line)
        rmsd_found = re.search(rmsd_regex, line)
        aligned_len_found = re.search(aligned_len_regex, line)
        if tm_score_found:
            if flag:  # Get 2nd TM-score
                tm_score2 = float(tm_score_found.group(1))
            else:     # Get 1st TM-score
                tm_score1 = float(tm_score_found.group(1))
                flag = True
        if rmsd_found:
            rmsd = float(rmsd_found.group(1))
        if aligned_len_found:
            aligned_len = int(aligned_len_found.group(1))
    # In this case, a reason could be that the first protein from which the PU are generated
    # is smaller than the protein 2 on which they are supposed to be aligned. We juste return -1.
    if tm_score2 is None:
        return -1
    return (tm_score1, tm_score2, rmsd, aligned_len)


def parse_tm_score(tm_score_res):
    """
    This function parses TMscore's terminal output.
    It retrieves the TMscore normalized by the length of chain 2.

    Args:
        tm_score_res (list): Terminal's output of TMscore

    Returns:
        float: The TMscore
    """
    tm_score_regex = re.compile("^TM-score\\s*=\\s*(\\d+\\.\\d+).*$")
    tm_score = None
    for line in tm_score_res:
        tm_score_found = re.search(tm_score_regex, line)
        if tm_score_found:
            tm_score = float(tm_score_found.group(1))
    assert tm_score is not None, "Error: The TMscore could not be parsed from TMscore"
    return tm_score


def parse_protein_peeling(peeling_res):
    """
    This function parses the output of Protein Peeling software.

    Args:
        peeling_res (str): Terminal's output of Protein Peeling 3 software.

    Returns:
        dictionary: keys are => "Number of PUs for each peeling level"
                                "Boundaries of each PU for a given peelling level"
    """
    peeling_dict = {"NB_PU": [], "PU_BOUNDS": []}
    res_line_regex = re.compile("^[^#].*$")
    for line in peeling_res:
        # Looking for lines not starting with "#"
        res_line_found = re.search(res_line_regex, line)
        if res_line_found:
            res_line = res_line_found.group(0).split()
            peeling_dict["NB_PU"].append(int(res_line[4]))
            pu_bounds = [(x, y) for x, y in zip(res_line[5::2], res_line[6::2])]
            # Sort the tuples in place
            pu_bounds.sort(key=lambda tup: int(tup[0]))
            peeling_dict["PU_BOUNDS"].append(pu_bounds)
    return peeling_dict
