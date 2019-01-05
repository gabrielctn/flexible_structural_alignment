"""
.. module:: Utils
   :synopsis: This module implements few utility functions.
"""

# Third-party modules
from pathlib import Path
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
    Removes files given a regex-like pattern.
    """
    for p in Path(dir).glob(pattern):
        p.unlink()
