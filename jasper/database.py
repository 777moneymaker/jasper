"""This module manages the blast operations.

This module uses Biopython package for making local blast queries and parsing the output.

More about it:
    1. https://biopython.org/

    2. https://biopython.org/docs/dev/api/Bio.Blast.Applications.html

    3. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2447716/
"""
from __future__ import annotations

import time
import io
import pandas as pd
from pathlib import Path
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

from jasper import utils

TYPES = ('fasta', 'fa', 'fna')


class Database:
    def __init__(self, source_dir: Path, name: str) -> None:
        """Inits database with given name"""
        if not Path(source_dir).exists():
            raise FileNotFoundError("Given path does not exist.")
        self.source_dir = Path(source_dir)
        self.name = name

    def _read_fasta(self, file: Path) -> str:
        """This function reads single fasta file and returns its content.

        Raises:
            FileNotFoundError: When given path is not file.
        Returns:
            (str): Content of the file.
        """
        if not file.is_file():
            raise FileNotFoundError("Given path is not a file.")

        with open(file, 'r') as fh:
            return "".join(fh.readlines())

    def _repair_fasta(self, file: Path):
        """This function reads single fasta file and returns its content.

        Raises:
            FileNotFoundError: When given path is not file.
        Returns:
            (str): Repaired content of the file.
        """

        if not file.is_file():
            raise FileNotFoundError("Given path is not a file.")

        repaired_content: list = []
        contig: int = 1
        with open(file, 'r') as fh:
            for line in fh:
                if line.startswith(">"):
                    seq_id = f">{file.stem}|{contig}\n"
                    repaired_content.append(seq_id)
                    contig += 1
                else:
                    repaired_content.append(line)
        return "".join(repaired_content)

    def create(self) -> Database:
        """Function for making local blast database.

        Returns:
            (bool): Value indicating if database creation succeeded.

        Creates:
            (*.nhr, *.nin, *.nsq): Created database's files in LMBD format.
        """

        print("Processing source files...")
        try:
            start: float = time.time()

            for i, source_file in enumerate(self.source_dir.iterdir(), 1):
                if not source_file.name.endswith(TYPES):
                    continue
                with open('temp.fasta', 'a+') as fh:
                    fh.write(self._repair_fasta(source_file))
                    print(f'Processed {i} host files', end='\r')

            end: float = time.time()
            print(f"Source files aggregation time: {end - start :.2f}")

            stdout, stderr = NcbimakeblastdbCommandline(input_file="temp.fasta",
                                                        dbtype="nucl",
                                                        title=self.name,
                                                        out=self.name)()
            print(stdout)
            Path("temp.fasta").unlink()
            return self
        except Exception as e:
            print(e)
            exit(0)

    def query(self, query_file: Path, config: dict, blast_format: str = "10 qseqid sseqid score"):
        """Function for making a BLAST query with database.

        Args:
            query_file (str): Query file that will be used as query.
            blast_format (str): Format for blast results.
            config (dict): Configuration for BLAST.
        Raises:
            FileNotFoundError: When given vir_file is not a file.
            ValueError: File has wrong extension.
        Returns:
            tuple or int: Tuple(query.name, target.name, best score) or 0 if query invalid.
        """
        if not query_file.is_file():
            raise FileNotFoundError('Given path is not a file.')
        if not query_file.name.endswith(TYPES):
            raise ValueError("Given file must end with *.fa | *.fna | .*fasta")

        blast_result_file: Path = Path(f"{query_file.stem}.score")
        try:
            NcbiblastnCommandline(query=query_file,
                                  db=f"{self.name}",
                                  outfmt=blast_format,
                                  out=blast_result_file,
                                  num_alignments=1,
                                  num_threads=5,
                                  **config)()
            with open(blast_result_file, 'r') as fh:
                line: list = fh.readline().rstrip().split(',')
                result: tuple = tuple(line)
        except (IndexError, ValueError):
            pass
        finally:
            Path(blast_result_file).unlink()
        return result

    def query_multiple(self, query_dir: Path, config: dict,
                       headers=("Virus", "Host", "Score"), blast_format: str = "10 qseqid sseqid score") -> pd.DataFrame:
        """TODO"""

        if not query_dir.is_dir():
            raise FileNotFoundError('Given path is not a directory.')

        print("Processing target files...")
        start: float = time.time()

        for query_file in query_dir.iterdir():
            with open('temp_vir.fasta', 'a+') as fh:
                if not query_file.name.endswith(TYPES):
                    continue
                fh.write(self._repair_fasta(query_file))  # Repair target seq
        end: float = time.time()
        print(f"Target files aggregation ended. Time: {end - start:.2f}")

        start = time.time()
        print("Quering...")
        stdout, stderr = NcbiblastnCommandline(query="temp_vir.fasta",
                                               db=f"{self.name}",
                                               outfmt=blast_format,
                                               num_threads=5,
                                               num_alignments=1,
                                               **config)()
        end = time.time()
        print(f"Query time: {end - start :.2f}")
        Path("temp_vir.fasta").unlink()
        results_df: pd.DataFrame = pd.read_csv(io.StringIO(stdout), header=None, names=headers)
        return results_df

    def clear_files(self):
        for file in Path(".").iterdir():
            if file.stem == self.name and file.suffix in (".nhr", ".nin", ".nsq"):
                file.unlink()
