"""This module manages the blast operations.

This module uses Biopython package for making local blast queries and parsing the output.

More about it:
    1. https://biopython.org/

    2. https://biopython.org/docs/dev/api/Bio.Blast.Applications.html

    3. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2447716/
"""
import os
import time
from pathlib import Path

from Bio.Blast.Applications import NcbiblastnCommandline as BlastCmd

TYPES = ('fasta', 'fa', 'fna')


def make_local_db(host_path: str):
    """Function for making local blast database."""
    host_path = Path(host_path)
    if not host_path.is_dir():
        raise IsADirectoryError('Given path is not a directory.')

    print("Combining files...")
    cmd = "find " + str(host_path.absolute()) + "/ \( -iname \"*.fasta\" -o -iname \"*.fa\" -o -iname \"*.fna\" \) -exec cat {} \; > temp.fasta"
    os.system(cmd)
    BlastCmd(dbtype="nucl", input_file="temp.fasta",
             title="Host_DB", out="host_db")()
    os.remove("temp.fasta")


def query(vir_file: str):
    """Function for making a BLAST query with local database.

    Raises:
        IsADirectoryError: When given path is not a directory.

    Returns:
        list: List containing tuples (query.name, target.name, best score)
    """
    vir_file = Path(vir_file)
    if not vir_file.is_file():
        raise OSError('Given path is not a file.')

    out_file = f"{vir_file.stem}.txt"
    try:
        if vir_file.name.endswith(TYPES):
            BlastCmd(query=vir_file, db="host_db", outfmt="10 qseqid sseqid score",
                     out=out_file, num_alignments=1)()
            with open(out_file, 'r') as fh:
                line = fh.readline().rstrip().split(',')
                result = (line[0], line[1], line[2])
    except (IndexError, ValueError):
        result = None
    finally:
        os.remove(out_file)

    return result
