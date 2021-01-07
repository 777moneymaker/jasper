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
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

from jasper import utils

TYPES = ('fasta', 'fa', 'fna')


class Database:
    def __init__(self, config: dict) -> None:
        """Inits database with given name"""
        self.host_dir: Path = Path(config["host_path"])
        self.name: str = config["db_name"]
        self._repair_host_files: bool = config["repair_host_files"]
        self._repair_vir_files: bool = config["repair_host_files"]

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

    def create(self) -> bool:
        """Function for making local blast database.

        Returns:
            (bool): Value indicating if database creation succeeded.

        Creates:
            (*.nhr, *.nin, *.nsq): Created database's files in LMBD format.
        """

        print("Processing host files...")
        try:
            start: float = time.time()
            for i, host_file in enumerate(self.host_dir.iterdir(), 1):
                if not host_file.name.endswith(TYPES):
                    continue
                with open('temp.fasta', 'a+') as fh:
                    if self._repair_host_files:
                        fh.write(self._repair_fasta(host_file))
                    else:
                        fh.write(self._read_fasta(host_file))
                print(f'Processed {i} host files', end='\r')
            end: float = time.time()
            print(f"Host files concatenation time: {end - start}")
            stdout, stderr = NcbimakeblastdbCommandline(input_file="temp.fasta",
                                                        dbtype="nucl",
                                                        title="Host_DB",
                                                        out=f"{self.name}")()
            Path("temp.fasta").unlink()

            print(stdout)
            return True
        except Exception as e:
            print(e)
            return False

    def query(self, vir_file: Path):
        """Function for making a BLAST query with database.

        Args:
            vir_file (str): Virus file that will be used as query.
        Raises:
            TypeError: When file is not a Path object.
            FileNotFoundError: When given vir_file is not a file.
            ValueError: File has wrong extension.
        Returns:
            tuple or int: Tuple(query.name, target.name, best score) or 0 if query invalid.
        """
        if not isinstance(vir_file, Path):
            raise TypeError("File is not a Path object.")
        if not vir_file.is_file():
            raise FileNotFoundError('Given path is not a file.')
        if not vir_file.name.endswith(TYPES):
            raise ValueError("Given file must end with *.fa | *.fna | .*fasta")

        out_file: str = f"{vir_file.stem}.score.txt"
        try:
            NcbiblastnCommandline(query=vir_file,
                                  db=f"{self.name}",
                                  outfmt="10 qseqid sseqid score",
                                  out=out_file,
                                  num_alignments=1,
                                  num_threads=5)()
            with open(out_file, 'r') as fh:
                line: list = fh.readline().rstrip().split(',')
                result: tuple = (line[0], line[1], line[2])
        except (IndexError, ValueError):
            result: int = 0
        finally:
            Path(out_file).unlink()
        return result

    def query_multiple(self, vir_dir: str) -> pd.DataFrame:
        """TODO"""
        vir_dir: Path = Path(vir_dir)
        if not vir_dir.is_dir():
            raise FileNotFoundError('Given path is not a directory.')

        print("Processing virus files...")
        start: float = time.time()
        i = 1
        for vir_fl in vir_dir.iterdir():
            if i == 20: break
            with open('temp_vir.fasta', 'a+') as fh:
                if not vir_fl.name.endswith(TYPES):
                    continue
                if self._repair_host_files:
                    fh.write(self._repair_fasta(vir_fl))
                else:
                    fh.write(self._read_fasta(vir_fl))
                i += 1
        end: float = time.time()
        print(f"Virus files concatenation ended. Time: {end - start}")

        start = time.time()
        print("Quering...")
        stdout, stderr = NcbiblastnCommandline(query="temp_vir.fasta",
                                               db=f"{self.name}",
                                               outfmt="10 qseqid sseqid score",
                                               # out="vir_results.txt",
                                               num_threads=5,
                                               num_alignments=1)()
        end = time.time()
        print(f"Query time: {end - start}")
        Path("temp_vir.fasta").unlink()

        results_df: pd.DataFrame = pd.read_csv(io.StringIO(stdout), header=None, names=["Virus", "Host", "Score"])
        return results_df.groupby(["Virus", "Host"]).max().sort_values(by='Score', ascending=False)
