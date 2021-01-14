"""This module manages the blast operations.

This module uses Biopython package for making local blast queries and parsing the output.

More about it:
    1. https://biopython.org/

    2. https://biopython.org/docs/dev/api/Bio.Blast.Applications.html

    3. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2447716/
"""
from __future__ import annotations

import io
import os
import subprocess
import argparse
import pandas as pd
from pathlib import Path
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

from . import utils


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

    def create(self) -> tuple:
        """Function for making local blast database.

        Returns:
            (bool): Value indicating if database creation succeeded.

        Creates:
            (*.nhr, *.nin, *.nsq): Created database's files in LMBD format.
        """

        if Path('blast_input.fasta').exists():
            Path('blast_input.fasta').unlink()
        for i, source_file in enumerate(self.source_dir.iterdir(), 1):
            if not source_file.name.endswith(utils.TYPES):
                continue
            with open('blast_input.fasta', 'a+') as fh:
                fh.write(self._repair_fasta(source_file))

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

    def query(self, query_dir: Path, config: dict,
              headers=("Virus", "Host", "Score"),
              blast_format: str = "10 qseqid sseqid score") -> pd.DataFrame:
        """TODO"""

        if not query_dir.is_dir():
            raise FileNotFoundError('Given path is not a directory.')

        if Path("blast_query.fasta").exists():
            Path("blast_query.fasta").unlink()
        for query_file in query_dir.iterdir():
            with open('blast_query.fasta', 'a+') as fh:
                if not query_file.name.endswith(utils.TYPES):
                    continue
                fh.write(self._repair_fasta(query_file))  # Repair target seq
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
        except Exception:
            raise
        finally:
            if Path("blast_query.fasta").exists():
                Path("blast_query.fasta").unlink()
        results_df: pd.DataFrame = pd.read_csv(io.StringIO(blastn_output.stdout.decode()), header=None, names=headers)
        return results_df

    def clear_files(self):
        for file in Path(".").iterdir():
            if file.stem == self.name and file.suffix in (".nhr", ".nin", ".nsq"):
                file.unlink()

def main(args):
    print(utils.LOGO)
    args.blastn_config = utils.parse_config(args.blastn_config, {
        "task": "blastn",
    })
    print("Starting analysis...")
    if args.create_db_name:
        print("Aggregating files, creating database...")
        db, db_output = Database(Path(args.host_dir), args.create_db_name).create()
        print(db_output)
    else:
        db = Database(Path('.'), args.use_db_name)
    print("Quering...")
    query_df = db.query(Path(args.virus_dir), config=args.blastn_config)

    query_df['Score'] = query_df['Score'].apply(pd.to_numeric)
    mega_results = query_df.loc[query_df.reset_index().groupby(['Virus'])['Score'].idxmax()]
    genome_results: pd.DataFrame = mega_results.sort_values(by="Score", ascending=False).reset_index(drop=True)
    print("Blastn results (genome-genome query): ", genome_results, sep='\n')
    if args.clear_after:
        db.clear_files()

    genome_results.to_csv(args.output_file, index=False)
    print("Saved files to", args.output_file)
