"""This module manages the IO operations

This module uses Biopython package for sequence IO operations using Bio.SeqIO.

More about it:
    1. https://biopython.org/

    2. https://biopython.org/docs/1.75/api/Bio.SeqIO.html
"""

from pathlib import Path

from Bio import SeqIO


def retrieve(directory: str):
    """Generator for retrieving sequences.

    This function retrieves all sequences in FASTA format in given directory.

    Yields:
        seq: Sequences retrieved from directory.

    Raises:
        IsADirectoryError: When given path is not directory.
        OSError: When given path is an empty directory.
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise IsADirectoryError("Given path is not directory.")
    if len(list(directory.iterdir())) == 0:
        raise OSError("Given directory is empty.")
    for file in directory.iterdir():
        for seq in SeqIO.parse(file, "fasta"):
            yield seq
