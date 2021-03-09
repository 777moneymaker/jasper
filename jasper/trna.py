#!/usr/bin/env python3
"""This module manages the tRNA scanning and query operations.

This module uses Biopython, tRNAscan-SE 2.0 and NCBI+Blast+ software for prediction of host and phage tRNAs and
for creation of tRNA database with tRNA query.

More about it:
    1. https://biopython.org/

    2. http://lowelab.ucsc.edu/tRNAscan-SE/

    3. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2447716/
"""
import io
import os
import subprocess
from pathlib import Path

import pandas as pd
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline

from . import blast
from . import utils


class tRNAScanner:
    def __init__(self, source_dir: Path, target_dir: Path, name: str) -> None:
        """Inits obj with args"""
        if not isinstance(source_dir, Path) or not isinstance(target_dir, Path):
            raise TypeError("Given object is not a Path object")
        if not isinstance(name, str):
            raise TypeError(f"Name expected to be str")
        if not source_dir.exists() or not target_dir.exists():
            raise FileNotFoundError("Given path does not exist.")
        if not source_dir.is_dir() or not target_dir.is_dir():
            raise FileNotFoundError("Given path is not a directory")
        self.source_dir = source_dir
        self.target_dir = target_dir
        self.name = name

    def clear_files(self):
        """Function clears database files if any present."""
        for file in Path(".").iterdir():
            if file.stem == self.name and file.suffix in (".nhr", ".nin", ".nsq"):
                file.unlink()

    def scan(self, source_dir: Path, output_filename: str, threads: int = os.cpu_count()):
        """This function aggregates every host file and every page file into two summary files.

        After that, this function runs the tRNAscan-SE 2.0 and returns the output *.fasta file (Path obj).

        Args:
            source_dir (Path): Path to directory which should be scanned.
            output_filename (Path): Path to file used as tRNAscan output.
            threads (int): Num of threads to run tRNAscan on.
        Returns:
            (Path): Path to output file.
        """
        aggregate = Path("trna_aggregate.fasta")
        for i, source_file in enumerate(source_dir.iterdir(), 1):
            if not source_file.name.endswith(utils.TYPES):
                continue

            with open(aggregate, 'a+') as aggregate_fh:
                aggregate_fh.write(blast.Database.repair_fasta(source_file))
        out_file = Path(output_filename)
        trna_output = self.run_trnascan(aggregate, out_file, threads)
        aggregate.unlink()

        return out_file

    def run_trnascan(self, file: Path, out_file: Path, num_threads: int):
        """This function runs the tRNAscan subprocess

        Args:
            file (Path): Path to file for tRNAscan to use.
            out_file (Path): Path to file used as output.
            num_threads (int): Num of threads to run tRNAscan on.
        """
        res = subprocess.run(['tRNAscan-SE', str(file), '--fasta', str(out_file), '--thread', str(num_threads), '-b'],
                             capture_output=True)
        return res

    def create(self, input_file: Path):
        """Function for making local blast database.

        This function creates database from tRNAs retrieved from tRNAscan-SE.

        Args:
            input_file (Path): Path to file containing DB sequences.
        Returns:
            str: Database output.
        Creates:
            (*.nhr, *.nin, *.nsq): Created database's files in LMBD format.
        Raises:
            SubprocessError: When makeblastdb returns error or when input file does not exist.
        """

        try:
            cmd = NcbimakeblastdbCommandline(input_file=str(input_file),
                                             dbtype="nucl",
                                             title=self.name,
                                             out=self.name)
            cmd()
            makeblastdb_output = subprocess.run(str(cmd), capture_output=True, shell=True)
            if makeblastdb_output.stderr:
                raise subprocess.SubprocessError(f"Makeblastdb returned error: {makeblastdb_output.stderr.decode()}")
        except Exception:
            raise
        finally:
            if input_file.exists():
                input_file.unlink()
        return makeblastdb_output.stdout.decode()

    def query(self, query_file: Path, config: dict, blast_format: str, headers: tuple):
        """This function queries file (aggregated or not) to created database.

        Args:
            query_file (Path): Path to file containing query seqs.
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
        if not isinstance(query_file, Path):
            raise TypeError("Given object is not Path object")
        if not query_file.exists():
            raise FileNotFoundError("Given path does not exist")
        if not query_file.is_file():
            raise FileNotFoundError("Given path is not file")

        if not isinstance(config, dict):
            raise TypeError("Config file is not a dict object")
        if any(kwarg in ('query', 'db', 'outfmt', 'max_target_seqs', 'num_alignments') for kwarg in config.keys()):
            used = filter(lambda k: k in config.keys(), ('query', 'db', 'outfmt', 'max_target_seqs', 'num_alignments'))
            raise ValueError("Given kwargs are not valid in terms of blast usage", list(used))

        try:
            cmd = NcbiblastnCommandline(query=str(query_file),
                                        db=self.name,
                                        outfmt=blast_format,
                                        **config)
            blast_output = subprocess.run(str(cmd), capture_output=True, shell=True)

            # Error only occurs if it's not this stupid warning.
            if blast_output.stderr and "Examining 5 or more matches" not in blast_output.stderr.decode():
                raise subprocess.SubprocessError(f"Blastn returned error: {blast_output.stderr.decode()}")
        except Exception:
            raise
        finally:
            if query_file.exists():
                query_file.unlink()
        results_df: pd.DataFrame = pd.read_csv(io.StringIO(blast_output.stdout.decode()), header=None, names=headers)
        return results_df


def main(args):
    print(utils.LOGO)

    args.blastn_config = utils.parse_config(args.blastn_config, {
        "task": "blastn",
        "num_threads": os.cpu_count(),
    })

    t = tRNAScanner(source_dir=Path(args.host_dir), target_dir=Path(args.virus_dir), name="trna_db")
    print("Aggregation and scanning host files...")
    host = t.scan(t.source_dir, output_filename="host_trnas.fasta", threads=args.num_threads)
    print("Done\n")

    print("Aggregation and scanning phage files...")
    vir = t.scan(t.target_dir, output_filename="phage_trnas.fasta", threads=args.num_threads)
    print("Done\n")

    print("Making a trna host database...")
    db_output = t.create(host)

    print("Quering a phage trnas")
    results_df = t.query(vir, {"task": "blastn"}, blast_format="10 qseqid sseqid score",
                         headers=("Virus", "Host", "Score"))
    print("Done")

    results_df['Virus'] = results_df['Virus'].apply(lambda x: x.split('.')[0])
    results_df['Host'] = results_df['Host'].apply(lambda x: x.split('.')[0])

    results_df['Virus'] = results_df['Virus'].map(lambda x: x.split("|")[0])
    results_df['Host'] = results_df['Host'].map(lambda x: x.split("|")[0])

    results_df = results_df.groupby(['Virus', 'Host'], observed=True).max().reset_index()
    print(results_df)

    results_df.to_csv(args.output_file, index=False)
    print("Saved trnas results")
    if args.clear_after:
        t.clear_files()


if __name__ == '__main__':
    pass
