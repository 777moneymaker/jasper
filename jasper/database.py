"""This module manages the blast operations.

This module uses Biopython package for making local blast queries and parsing the output.

More about it:
    1. https://biopython.org/

    2. https://biopython.org/docs/dev/api/Bio.Blast.Applications.html

    3. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2447716/
"""
import os
import threading
from pathlib import Path
import time
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

from jasper import utils

TYPES = ('fasta', 'fa', 'fna')


class Database:
    global_lock = threading.Lock()
    file_contents = []

    def __init__(self, host_path: str, db_name: str = "host_db"):
        """Inits database with given name"""
        self.name = db_name
        self.source_dir = host_path
        self._create()

    def _concat_file(self, host_fl: Path):
        """This function creates command for joining multiple host fasta files.

        Raises:
            IsADirectoryError: When source_dir is not directory.
        Returns:
            (bool): Value indicating if file concatenation succeeded.
        """
        while Database.global_lock.locked():
            continue
        with open(host_fl, 'r') as host_fh:
            Database.global_lock.acquire()
            Database.file_contents.append(host_fh.read())
            Database.global_lock.release()

    def _create(self):
        """Function for making local blast database.

        Raises:
            FileNotFoundError: When given host file is not a file.
        Creates:
            (*.nhr, *.nin, *.nsq): Created database's files in LMBD format.
        """
        start = time.time()
        print("Processing host files in chunks...")
        for host_chunk in utils.chunks(list(Path(self.source_dir).iterdir()), 5):
            threads = []
            for host_fl in host_chunk:
                t = threading.Thread(target=self._concat_file, args=[host_fl])
                threads.append(t)
                t.start()
            [thread.join() for thread in threads]
            with open('temp.fasta', 'a+') as fh:
                fh.write("\n".join(Database.file_contents))
                Database.file_contents.clear()
        end = time.time()
        print(end - start)
        stdout, stderr = NcbimakeblastdbCommandline(input_file="temp.fasta", dbtype="nucl", title="Host_DB", out=f"{self.name}")()
        print(stdout)
        os.remove("temp.fasta")

    def delete_concat_fasta_files(self, db_name: str = "temp.fasta"):
        # TODO: Write a function that will delete file containing concatenated host fasta files.
        pass

    def query(self, vir_file: str):
        """Function for making a BLAST query with local database.

        Args:
            vir_file (str): Virus file that will be used as query.
        Raises:
            FileNotFoundError: When given vir_file is not a file.
        Returns:
            list: List containing tuples (query.name, target.name, best score)
        """
        vir_file = Path(vir_file)
        if not vir_file.is_file():
            raise FileNotFoundError('Given path is not a file.')

        out_file = f"{vir_file.stem}.score.txt"
        try:
            if vir_file.name.lower().endswith(TYPES):
                NcbiblastnCommandline(query=vir_file, db=f"{self.name}", outfmt="10 qseqid sseqid score",
                                      out=out_file, num_alignments=1)()
                with open(out_file, 'r') as fh:
                    line = fh.readline().rstrip().split(',')
                    result = (line[0], line[1], line[2])
        except (IndexError, ValueError):
            result = 0
        finally:
            os.remove(out_file)
        return result
