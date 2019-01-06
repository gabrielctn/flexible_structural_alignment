"""
.. module:: pdb
   :synopsis: This module implements functions for PDB files manipulation.
"""

# Third-party modules
import os
import logging
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Dice
from Bio.PDB import Select


# Set logging info
logging.basicConfig(filename='flex.log', filemode='w',
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S')


def write_pu_pdb(all_pu, pdb_name_1, peeling_level, structure):
    """
    This function writes down the PDB files of all the PUs of a given peeling level.

    Args:
        all_pu (list of tuples): list of the boundaries of the PUs (start, end)
        pdb_name_1 (str): Name of the PDB of the protein 1
        peeling_level (int): Peeling level
        structure (BIO.PDB structure): Parsed structure of the protein 1
    """
    # Write the PDB files for the PUs of the actual peeling level
    for pu_index, (start, end) in enumerate(all_pu):
        # File name = pdb_name_1 + peeling_level + "pu" + PU_index + ".pdb"
        pu_pdb = "results/{}_{}_pu_{}.pdb".format(pdb_name_1, peeling_level+1, pu_index+1)
        # Write the PDB file for the actual PU
        Dice.extract(structure, "A", int(start), int(end), pu_pdb)


def remove_aligned_pu_from_pdb(structure_in, pdb_out, start, end):
    """
    This function writes a new PDB file based on the "structure_in" one, without the residues
    between start and end.
    In other terms, it extracts the residues of the PU (between start-end indexes) from the native
    structure (the PDB on which the PU was just aligned on).

    Args:
        structure_in (PDBParser Structure): The structure of the native PDB.
        pdb_out (str): Filename of the cropped PDB.
        start (int): Start residue index of the PU (included)
        end (int): End residue index of the PU (included)
    """
    class ResSelect(Select):
        def accept_residue(self, res):
            return False if res.id[1] >= start and res.id[1] <= end else True
    io = PDBIO()
    io.set_structure(structure_in)
    io.save(pdb_out, ResSelect())


def clean_pdb(pdb_file, pdb_name):
    """
        Reindex the PDB given in argument.
        Residue id starting at 1 and keep only chain A.

        Args:
            pdb_file (str): PDB file to clean
            pdb_name (str): Basenale of the PDB file
    """
    # The PDB files are reindexed to fit the Protein Peeling program and keep only chain A
    reindexed_pdb = reindex_pdb(1, pdb_file, True)
    # Write the new reindexed PDBs
    new_pdb_name = pdb_name + "_new.pdb"
    new_pdb_file = "data/" + new_pdb_name
    with open(new_pdb_file, "w") as f_out:
        f_out.write(reindexed_pdb)
    return new_pdb_name, new_pdb_file



def reindex_pdb_by_index(in_file, start_index=1, pdb_txt=""):
    """
        Reindex the residue numbers of PDB given in argument
        Source code: https://zhanglab.ccmb.med.umich.edu/reindex_pdb/reindex_pdb.py

        options:
            in_file (str): PDB file being processed
            start_index(int): Index of first residue
            pdb_txt (str): Text of input PDB to be reindexed
    """
    pdb_txt_reindex = ''
    current_old_index = ''  # residue number in origin PDB
    warn_chain_id = ''  # warning about new chain ID

    for line in pdb_txt.splitlines():
        if len(line) < 27 or (not line.startswith("ATOM  ")
                              and not line.startswith("HETATM") and not line.startswith("TER")):
            pdb_txt_reindex += line+'\n'
            continue
        elif not line[16] in ['A', ' ']:  # alternative location identifier
            continue
        res_seq = line[22:27]  # residue sequence number
        current_chain_id = line[21]  # chain identifier

        if not current_old_index:  # first residue encountered
            current_old_index = res_seq  # residue number in origin PDB
            current_new_index = int(start_index)
            chain_id = current_chain_id
            res_seq_new = str(current_new_index)
            res_seq_new = ' '*(4-len(res_seq_new))+res_seq_new+' '
        elif current_chain_id != chain_id:
            if warn_chain_id != current_chain_id:
                logging.warning("Discarding chain '%s' of %s\n" % (current_chain_id, in_file))
                warn_chain_id = current_chain_id
            continue
        elif res_seq != current_old_index:
            current_new_index += 1
            current_old_index = res_seq
            res_seq_new = str(current_new_index)
            res_seq_new = ' '*(4-len(res_seq_new))+res_seq_new+' '
        pdb_txt_reindex += line[:16]+' '+line[17:22]+res_seq_new+line[27:]+'\n'
    return pdb_txt_reindex


def reindex_pdb(start_index, in_file, clean=True):
    """
        Parse the PDB file "in_file", reindex it according to start index (1) by default,
        and return the text of renumbered PDB

        Args:
            start_index (int): Starting index of the reindexing (1 by default)
            in_file (str): PDB file to reindex
            clean (bool): Have a clean end of file, no flourishes
    """
    f_in = open(in_file, 'rU')
    pdb_txt = ''
    for line in f_in.read().splitlines():
        if line.startswith("END"):
            if clean:
                line = line.replace("ENDMDL", "END   ")
            pdb_txt += line+'\n'
            break
        if (line.startswith("ATOM  ") or line.startswith("TER")
            or (not clean and not line[:6]
                in ["DBREF ", "SEQADV", "MODRES", "HELIX ", "SHEET ", "SSBOND", "SITE  "])):
            pdb_txt += line+'\n'
    f_in.close()

    pdb_txt_reindex = reindex_pdb_by_index(in_file, start_index, pdb_txt)
    return pdb_txt_reindex


def save_best_aligned_pu(best_pu_peel_level, best_pu_index):
    """
        Save the coordinates (PDB) of the PU that was best aligned with TM-align.
    """
    # Save the coordinates of the best aligned PU
    best_tm_aligned_file = "tmp/TMaligned_"+str(best_pu_peel_level)+"_"+str(best_pu_index+1)+".sup"
    with open(best_tm_aligned_file + "_all_atm", "r") as f_in:
        line = f_in.readline()
        # Skip the PDB header, go to the coordinates of the aligned PU
        while line[0:4].strip() != "ATOM":
            line = f_in.readline()
        # Save the coordinates of the aligned PU
        aligned_max_pu = "tmp/{}_aln_max_pu_{}.pdb".format(best_pu_peel_level+1, best_pu_index+1)
        with open(aligned_max_pu, "w") as f_out:
            while line[0:4].strip() != "TER":
                f_out.write(line)
                line = f_in.readline()

def concatenate_best_aligned_pus(peeling_level, aligned_pu_pdb, all_pu_ref):
    """
        Write the full PDB of the best aligned PUs by concatenating them in the right order.

        Args:
            peeling_level (int): Actual peeling level
            aligned_pu_pdb (str): Path to the PDB into which will be written all the best aligned PUs
            all_pu_ref (list of tuples): Full list of the PU bounds of the actual peeling level
    """
    open_mode = "w" if os.path.isfile(aligned_pu_pdb) else "a"
    with open(aligned_pu_pdb, open_mode) as f_out:
        # Loop over the PUs
        for pu_index in range(len(all_pu_ref)):
            # Find PDB of the best PU for this index
            aln_max_pu_file = "tmp/{}_aln_max_pu_{}.pdb".format(peeling_level+1, pu_index+1)
            # Concatenate to the full PDB
            with open(aln_max_pu_file, "r") as f_in:
                for line in f_in:
                    f_out.write(line)
