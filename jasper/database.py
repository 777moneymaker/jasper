"""This module manages the blast operations.

This module uses Biopython package for making local blast queries and parsing the output.

More about it:
    1. https://biopython.org/

    2. https://biopython.org/docs/dev/api/Bio.Blast.Applications.html

    3. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2447716/
"""
import time
import io
import pandas as pd
from pathlib import Path, PurePath
from __future__ import annotations
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

from jasper import utils

TYPES = ('fasta', 'fa', 'fna')


class Database:
    def __init__(self, config: dict) -> None:
        """Inits database with given name"""
        self.source_dir: Path = Path(config["source_path"])
        self.name: str = config["db_name"]

    @staticmethod
    def _read_fasta(file: Path) -> str:
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

    @staticmethod
    def _repair_fasta(file: Path):
        """This function reads single fasta file and returns its content.

        Raises:
            FileNotFoundError: When given path is not file.
        Returns:
            (str): Repaired content of the file.
        """

        if not file.is_file():
            print(file)
            raise FileNotFoundError("Given path is not a file.")

        repaired_content: list = []
        with open(file, 'r') as fh:
            i: int = 1
            for line in fh:
                if line.startswith(">"):
                    repaired_content.append(f">{file.stem}|{i}\n")
                    i += 1
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

        print("Processing host files...")
        try:
            start: float = time.time()
            for i, host_file in enumerate(self.source_dir.iterdir(), 1):
                if not host_file.name.endswith(TYPES):
                    continue
                with open('temp.fasta', 'a+') as fh:
                    fh.write(self._repair_fasta(host_file))
                print(f'Processed {i} host files', end='\r')
            end: float = time.time()
            print(f"Source files concatenation time: {end - start}")
            stdout, stderr = NcbimakeblastdbCommandline(input_file="temp.fasta",
                                                        dbtype="nucl",
                                                        title="Host_DB",
                                                        out=f"{self.name}")()
            Path("temp.fasta").unlink()

            print(stdout)
            return self
        except Exception as e:
            print(e)
            exit(0)

    def query(self, query_file: Path, config: dict, blast_format: str = "10 qseqid sseqid score"):
        """Function for making a BLAST query with database.

        Args:
            query_file (str): Query file that will be used as query.
            blast_format (str): Format for blast results.
        Raises:
            TypeError: When file is not a Path object.
            FileNotFoundError: When given vir_file is not a file.
            ValueError: File has wrong extension.
        Returns:
            tuple or int: Tuple(query.name, target.name, best score) or 0 if query invalid.
        """
        if not isinstance(query_file, Path):
            raise TypeError("File is not a Path object.")
        if not query_file.is_file():
            raise FileNotFoundError('Given path is not a file.')
        if not query_file.name.endswith(TYPES):
            raise ValueError("Given file must end with *.fa | *.fna | .*fasta")

        out_file: str = f"{query_file.stem}.score.txt"
        try:
            NcbiblastnCommandline(query=query_file,
                                  db=f"{self.name}",
                                  outfmt=blast_format,
                                  out=out_file,
                                  num_alignments=1,
                                  num_threads=5,
                                  **config)()
            with open(out_file, 'r') as fh:
                line: list = fh.readline().rstrip().split(',')
                result: tuple = (line[0], line[1], line[2])
        except (IndexError, ValueError):
            result: int = 0
        finally:
            Path(out_file).unlink()
        return result

    def query_multiple(self, query_dir: str, config: dict, headers=None) -> pd.DataFrame:
        """TODO"""
        if headers is None:
            headers = ["Virus", "Host", "Score"]
        query_dir: Path = Path(query_dir)
        if not query_dir.is_dir():
            raise FileNotFoundError('Given path is not a directory.')

        print("Processing target files...")
        start: float = time.time()
        for query_fl in query_dir.iterdir():
            with open('temp_vir.fasta', 'a+') as fh:
                if not query_fl.name.endswith(TYPES):
                    continue
                # Repair target seq
                fh.write(self._repair_fasta(query_fl))
        end: float = time.time()
        print(f"Target files concatenation ended. Time: {end - start}")

        start = time.time()
        print("Quering...")
        stdout, stderr = NcbiblastnCommandline(query="temp_vir.fasta",
                                               db=f"{self.name}",
                                               outfmt="10 qseqid sseqid score",
                                               num_threads=5,
                                               num_alignments=1,
                                               **config)()
        end = time.time()
        print(f"Query time: {end - start}")
        Path("temp_vir.fasta").unlink()
        if not headers:
            headers = ["Virus", "Host", "Score"]
        results_df: pd.DataFrame = pd.read_csv(io.StringIO(stdout), header=None, names=headers)
        return results_df

    def clear_files(self):
        for file in Path(".").iterdir():
            if file.stem == self.name and file.suffix in (".nhr", ".nin", ".nsq"):
                file.unlink()
