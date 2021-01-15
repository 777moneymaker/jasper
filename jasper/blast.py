#!/usr/bin/env python3
"""This module manages the blast operations.

This module uses Biopython and NCBI-Blast+ for making local blast queries and parsing the output.

More about it:
    1. https://biopython.org/

    2. https://biopython.org/docs/dev/api/Bio.Blast.Applications.html

    3. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2447716/
"""
from __future__ import annotations

import io
import os
import subprocess
import pandas as pd
from pathlib import Path
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

from . import utils


class Database:
    def __init__(self, source_dir: Path, name: str) -> None:
        """Inits database with given name
        Args:
            source_dir (Path): Path to directory containing host files.
            name (str): Name of database.
        Raises:
            TypeError: When source_dir or path is not valid obj type
            FileNotFoundError: When given path does not exist or is not a directory.
        """
        if not isinstance(source_dir, Path):
            raise TypeError("Given object is not a Path object")
        if not isinstance(name, str):
            raise TypeError(f"Name expected to be str")
        if not Path(source_dir).exists():
            raise FileNotFoundError("Given path does not exist.")
        if not Path(source_dir).is_dir():
            raise FileNotFoundError("Given path is not a directory")
        self.source_dir = Path(source_dir)
        self.name = name

    def _repair_fasta(self, file: Path):
        """This function reads single fasta file and returns its content.

        Raises:
            TypeError: When given file is of wrong type.
            FileNotFoundError: When given path is not file.
        Returns:
            (str): Repaired content of the file.
        """
        if not isinstance(file, Path):
            raise TypeError("File is not a Path object.")
        if not file.is_file():
            raise FileNotFoundError("Given path is not a file.")

        repaired_content: list = []
        contig: int = 1
        with open(file, 'r') as fh:
            for line in fh:
                if line.startswith(">"):
                    seq_id = f">{file.stem}|{contig}\n"  # Read header
                    repaired_content.append(seq_id)
                    contig += 1
                else:
                    repaired_content.append(line)
        return "".join(repaired_content)

    def _aggregate(self, directory: Path, outfile: Path):
        """This function aggregates every file in source directory.

        Aggregated file consists of every repaired file in source dir.

        Args:
             directory (Path): Path to directory containing files.
             outfile (Path): Name of the aggregation file to be created.
        Raises:
            TypeError: When given object is of wrong type.
            FileNotFoundError: When given directory does not exist.
            ValueError: When output file is empty.
        """
        if not isinstance(directory, Path):
            raise TypeError("Given object is not Path object.")
        if not isinstance(outfile, Path):
            raise TypeError("Given object is not Path object.")
        if not directory.exists():
            raise FileNotFoundError('Given path does not exist.')
        if not directory.is_dir():
            raise FileNotFoundError('Given path is not a directory.')

        if outfile.exists():
            outfile.unlink()

        for source_file in directory.iterdir():
            if not source_file.name.endswith(utils.TYPES):
                continue
            with open(outfile, 'a+') as fh:
                fh.write(self._repair_fasta(source_file))

        if os.path.getsize(outfile) == 0:  # size == 0
            raise ValueError("Blast input file is empty. Check your input.")

    def create(self) -> tuple:
        """Function for making local blast database.

        This function creates database from files found in source_dir.

        Returns:
            tuple(Database, str): Database object, output from makeblastdb.
        Creates:
            (*.nhr, *.nin, *.nsq): Created database's files in LMBD format.
        Raises:
            SubprocessError: When makeblastdb returns error or when input file does not exist.
        """

        self._aggregate(self.source_dir, Path("blast_input.fasta"))
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
        except Exception:
            raise
        finally:
            if Path("blast_input.fasta").exists():
                Path("blast_input.fasta").unlink()
            return self, makeblastdb_output.stdout.decode()

    def query(self, query_dir: Path, config: dict, blast_format: str, headers: tuple) -> pd.DataFrame:
        """This function queries content of the directory to created database

        Args:
            query_dir (Path): Directory containing query files.
            config (dict): Blast configuration dict.
            blast_format (str): Blast output format.
            headers ( tuple(*str) ): Headers matching blast output for final DataFrame.
        Raises:
            TypeError: When given obj is of wrong type.
            FileNotFoundError: When given path does not exist or when given path is not a directory.
            ValueError: When forbidden blast option was provided.
        Returns:
            (pd.DataFrame): Pandas DataFrame containing query results.
        """

        if not isinstance(query_dir, Path):
            raise TypeError("Given object is not Path object")
        if not query_dir.exists():
            raise FileNotFoundError("Given path does not exist")
        if not query_dir.is_dir():
            raise FileNotFoundError("Given path is not directory")

        if not isinstance(config, dict):
            raise TypeError("Config file is not a dict object")
        if any(kwarg in ('query', 'db', 'outfmt', 'max_target_seqs', 'num_alignments') for kwarg in config.keys()):
            used = filter(lambda k: k in config.keys(), ('query', 'db', 'outfmt', 'max_target_seqs', 'num_alignments'))
            raise ValueError("Given kwargs are not valid in terms of blast usage", list(used))

        self._aggregate(query_dir, Path("blast_query.fasta"))
        try:
            cmd = NcbiblastnCommandline(query="blast_query.fasta",
                                        db=f"{self.name}",
                                        outfmt=blast_format,
                                        max_target_seqs=1,
                                        **config)
            blastn_output = subprocess.run(str(cmd), capture_output=True, shell=True)

            # Error only occurs if it's not this stupid warning.
            if blastn_output.stderr and "Examining 5 or more matches" not in blastn_output.stderr.decode():
                if not Path("blast_query.fasta").exists():
                    raise subprocess.SubprocessError("Blastn returned error. Input file for Blastn does not exist. Check your input.", str(cmd))
                else:
                    raise subprocess.SubprocessError("Blastn returned error. Check your input.", blastn_output.stderr.decode())
        except Exception:
            raise
        finally:
            if Path("blast_query.fasta").exists():
                Path("blast_query.fasta").unlink()
        results_df: pd.DataFrame = pd.read_csv(io.StringIO(blastn_output.stdout.decode()), header=None, names=headers)
        return results_df

    def clear_files(self):
        """Function clears database files if any present."""
        for file in Path(".").iterdir():
            if file.stem == self.name and file.suffix in (".nhr", ".nin", ".nsq"):
                file.unlink()


def main(args):
    """Main function for running module.
    Args:
        args (argparse obj): Arguments from right subparser.
    """
    print(utils.LOGO)

    # Default config if not specified
    args.blastn_config = utils.parse_config(args.blastn_config, {
        "task": "blastn",
        "num_threads": os.cpu_count(),
    })
    print("Starting analysis...")
    if args.create_db_name:
        print("Aggregating files, creating database...")
        db, db_output = Database(Path(args.host_dir), args.create_db_name).create()
        print(db_output)
    else:
        db = Database(Path('.'), args.use_db_name)
    print("Quering...")
    query_df = db.query(Path(args.virus_dir),
                        blast_format="10 qseqid sseqid score",
                        headers=("Virus", "Host", "Score"),
                        config=args.blastn_config)
    if args.clear_after:
        db.clear_files()

    query_df['Score'] = query_df['Score'].apply(pd.to_numeric)
    print("Blastn results (genome-genome query): ", query_df, sep='\n')

    query_df.to_csv(args.output_file, index=False)
    print("Saved files to", args.output_file)


class Expando(object):
    pass


if __name__ == "__main__":
    args = Expando()
    args.virus_dir = Path("example_data/virus")
    args.host_dir = Path("example_data/host")
    args.create_db_name = "host_db"
    args.blastn_config = ""
    args.clear_after = True
    args.output_file = "blast_results.csv"
    main(args)
