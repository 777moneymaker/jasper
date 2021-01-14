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
import os
import subprocess
import argparse
import pandas as pd
from pathlib import Path
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

from jasper import utils


class Database:
    def __init__(self, source_dir: Path, name: str) -> None:
        """Inits database with given name"""
        if not Path(source_dir).exists():
            raise FileNotFoundError("Given path does not exist.")
        self.source_dir = Path(source_dir)
        self.name = name

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
        start: float = time.time()
        if Path('blast_input.fasta').exists():
            Path('blast_input.fasta').unlink()
        for i, source_file in enumerate(self.source_dir.iterdir(), 1):
            if i == 100:
                break;
            if not source_file.name.endswith(utils.TYPES):
                continue
            with open('blast_input.fasta', 'a+') as fh:
                fh.write(self._repair_fasta(source_file))
                print(f'Processed {i} host files', end='\r')

        end: float = time.time()
        print(f"Source files aggregation time: {end - start :.2f}")

        try:
            cmd = NcbimakeblastdbCommandline(input_file="blast_input.fasta",
                                             dbtype="nucl",
                                             title=self.name,
                                             out=self.name)
            makeblastdb_output = subprocess.run(str(cmd), capture_output=True, shell=True)
            if makeblastdb_output.stderr:
                if not Path("blast_input.fasta").exists():
                    raise subprocess.SubprocessError(
                        "Makeblastdb returned error. Input file for Makeblastdb does not exist. Check your input.",
                        str(cmd))
                else:
                    raise subprocess.SubprocessError("Makeblastdb returned error. Check your input.",
                                                     makeblastdb_output.stderr.decode())
            else:
                print(makeblastdb_output.stdout.decode())
        except Exception as e:
            print(e)
            exit(0)
        finally:
            if Path("blast_input.fasta").exists():
                Path("blast_input.fasta").unlink()
            return self

    def query(self, query_dir: Path, config: dict,
              headers=("Virus", "Host", "Score"),
              blast_format: str = "10 qseqid sseqid score") -> pd.DataFrame:
        """TODO"""

        if not query_dir.is_dir():
            raise FileNotFoundError('Given path is not a directory.')

        print("Processing target files...")
        start: float = time.time()

        if Path("blast_query.fasta").exists():
            Path("blast_query.fasta").unlink()
        for query_file in query_dir.iterdir():
            with open('blast_query.fasta', 'a+') as fh:
                if not query_file.name.endswith(utils.TYPES):
                    continue
                fh.write(self._repair_fasta(query_file))  # Repair target seq
        end: float = time.time()
        print(f"Target files aggregation ended. Time: {end - start:.2f}")

        start: float = time.time()
        print("Quering...")
        try:
            cmd = NcbiblastnCommandline(query="blast_query.fasta",
                                        db=f"{self.name}",
                                        outfmt=blast_format,
                                        num_threads=os.cpu_count(),
                                        max_target_seqs=1,
                                        **config)
            blastn_output = subprocess.run(str(cmd), capture_output=True, shell=True)
            if blastn_output.stderr and "Examining 5 or more matches" not in blastn_output.stderr.decode():
                if not Path("blast_query.fasta").exists():
                    raise subprocess.SubprocessError("Blastn returned error. Input file for Blastn does not exist. Check your input.", str(cmd))
                else:
                    raise subprocess.SubprocessError("Blastn returned error. Check your input.", blastn_output.stderr.decode())
        except Exception as e:
            print(e)
            exit(0)
        finally:
            if Path("blast_query.fasta").exists():
                Path("blast_query.fasta").unlink()
        end: float = time.time()
        print(f"Query time: {end - start :.2f}")
        print(cmd)
        results_df: pd.DataFrame = pd.read_csv(io.StringIO(blastn_output.stdout.decode()), header=None, names=headers)
        return results_df

    def clear_files(self):
        for file in Path(".").iterdir():
            if file.stem == self.name and file.suffix in (".nhr", ".nin", ".nsq"):
                file.unlink()


def parse_args():
    parser = argparse.ArgumentParser(description="JASPER is a program for bacterial hosts prediction. This module is "
                                                 "performing genome-genome blast query with given data and config.",
                                     epilog="Made by Milosz Chodkowski 2020, PUT Poznan."
                                            "Check my github @ github.com/777moneymaker",
                                     usage="jasper.database [-h] [-vir VIRUS_DIR] [-c BLASTN_CONFIG] (--use_db USE_DB_NAME |  (--create_db CREATE_DB_NAME -hst HOST_DIR))")
    parser.add_argument("-vir", "--virus",
                        required=True,
                        type=str,
                        dest='virus_dir',
                        help='directory containing virus seq files.')
    parser.add_argument("-c", "--config",
                        required=False,
                        type=str,
                        dest='blastn_config',
                        help='File containing megablast config used in genome-genome query.')
    parser.add_argument("--clear",
                        action='store_true',
                        dest='clear_after',
                        help='Specifies if the database files should be deleted after analysis.')
    parser.add_argument('-o', '--output',
                        required=False,
                        type=str,
                        dest="output_file",
                        default="blast_results.csv",
                        help="Output file with final results.")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--use_db",
                       type=str,
                       dest='use_db_name',
                       help='Name for a existing database that will be used.')
    subgroup = group.add_argument_group()
    subgroup.add_argument("--create_db",
                          type=str,
                          dest='create_db_name',
                          help='Name for a new database that will be created.')
    subgroup.add_argument("-hst", "--host",
                          type=str,
                          dest='host_dir',
                          help='directory containing host seq files.')

    args = parser.parse_args()
    if args.use_db_name and any([args.create_db_name, args.host_dir]):
        parser.error("Mutually exclusive argument groups were used. Check usage.")
    if not args.use_db_name and not all([args.create_db_name, args.host_dir]):
        parser.error("You must specify database name and a host files directory. Check usage.")

    return args


if __name__ == '__main__':
    print(utils.LOGO)
    args = parse_args()
    args.blastn_config = utils.parse_config(args.blastn_config, {
        "task": "megablast",
    })
    if args.create_db_name:
        db = Database(Path(args.host_dir), args.create_db_name).create()
    else:
        db = Database(Path('.'), args.use_db_name)
    query_df = db.query(Path(args.virus_dir), config=args.blastn_config)

    query_df['Score'] = query_df['Score'].apply(pd.to_numeric)
    mega_results = query_df.loc[query_df.reset_index().groupby(['Virus'])['Score'].idxmax()]
    genome_results: pd.DataFrame = mega_results.sort_values(by="Score", ascending=False).reset_index(drop=True)
    print("Blastn results (genome-genome query): ", genome_results, sep='\n')

    genome_results.to_csv(args.output_file, index=False)

    if args.clear_after:
        db.clear_files()
