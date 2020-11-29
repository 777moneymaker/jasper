"""This module manages the blast operations.

This module uses Biopython package for making local blast queries and parsing the output.

More about it:
    1. https://biopython.org/

    2. https://biopython.org/docs/dev/api/Bio.Blast.Applications.html

    3. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2447716/
"""
import os

from Bio import SeqIO
from Bio.Blast import Applications
from jasper import io


def make_local_db(sequences: list):
    with SeqIO.write(sequences, 'temp/temp.fasta', 'fasta'):
        console = Applications.NcbimakeblastdbCommandline(dtype="nucl",
                                                          input_file="temp/temp.fasta",
                                                          title="host_db")
        os.remove("temp/temp.fasta")
        os.rmdir("temp")
    print(console)


if __name__ == "__main__":
    make_local_db(io.retrieve("tests/data/fasta_test_data/seqs.fasta"))