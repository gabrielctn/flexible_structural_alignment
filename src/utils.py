"""
.. module:: Utils
   :synopsis: This module implements few utility functions.
"""

# Third-party modules
import os
import shutil
import re
from schema import Schema, Use, SchemaError


def check_args(ARGUMENTS):
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


def clean_files(dir, pattern):
    """
    Removes files given a regex-like pattern and a directory.
    """
    for files in os.listdir(dir):
        if re.search(pattern, files):
            os.remove(os.path.join(dir, files))


def clean_directory(dir):
    """
        Delete a directory recursively.
    """
    shutil.rmtree(dir)


def create_dir_safely(dir):
    """
        Create the directory given in argument if not already created.
    """
    if not os.path.exists(dir):
        os.makedirs(dir)
