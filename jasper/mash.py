#!/usr/bin/env python3
"""This module manages the alignment-free operations and analysis.

This module uses Biopython and Mash software for virus's host prediction by using alignment-free methods.

More about it:
    1. https://biopython.org/

    2. https://github.com/soedinglab/WIsH
"""

import io
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from . import blast
from . import utils


class Mash:
    def __init__(self, host_dir: Path, phage_dir: Path):
        """Inits Mash object with 2 directories. One phage directory and one host directory.
        Args:
            host_dir (Path): Path to directory containing host files.
            phage_dir (str): Path to directory containing phage files.
        Raises:
            ValueError: When given path is not a directory.
            FileNotFoundError: When given path does not exist.
            TypeError: When given object is not a Path object.
        """
        if not isinstance(host_dir, Path) or not isinstance(phage_dir, Path):
            raise TypeError("Given object is not a Path object.")
        if not host_dir.exists():
            raise FileNotFoundError("Host directory does not exist.")
        if not host_dir.is_dir():
            raise ValueError("Host directory is not a directory.")
        if not phage_dir.exists():
            raise FileNotFoundError("Phage directory does not exist.")
        if not phage_dir.is_dir():
            raise ValueError("Phage directory is not a directory.")

        self.host_dir = host_dir
        self.phage_dir = phage_dir

    def sketch(self, host_sname: str, phage_sname: str, threads: int = os.cpu_count()):
        """This function uses Mash software to build a sketch for source files.
        For more information, check Mash documentation.

        Args:
            host_sname (str): Prefix for host sketch file.
            phage_sname (str): Prefix for phage sketch file.
            threads (int): Num of threads to use
        Returns:
            tuple(Path, Path): Paths to sketch files.
        """
        host_aggr = Path("host_mash_aggr.txt")
        phage_aggr = Path("phage_mash_aggr.txt")

        host_repaired_folder = Path("host_repaired_mash")
        phage_repaired_folder = Path("phage_repaired_mash")

        if not host_repaired_folder.exists():
            host_repaired_folder.mkdir()
        if not phage_repaired_folder.exists():
            phage_repaired_folder.mkdir()

        for i, host_file in enumerate(self.host_dir.iterdir(), 1):
            repaired = blast.Database.repair_fasta(host_file)
            for record in SeqIO.parse(io.StringIO(repaired), 'fasta'):
                with open(host_repaired_folder / Path(f"{record.id}.fasta"), 'w+') as fh:
                    SeqIO.write(record, fh, 'fasta')

        with open(host_aggr, 'w+') as fh:
            hosts = [str(host.absolute()) for host in host_repaired_folder.iterdir()]
            fh.write("\n".join(hosts))
        host_sketch = self.run_sketch(host_aggr, threads, host_sname)
        host_aggr.unlink()
        shutil.rmtree(host_repaired_folder)

        for i, phage_file in enumerate(self.phage_dir.iterdir(), 1):
            repaired = blast.Database.repair_fasta(phage_file)
            for record in SeqIO.parse(io.StringIO(repaired), 'fasta'):
                with open(phage_repaired_folder / Path(f"{record.id}.fasta"), 'w+') as fh:
                    SeqIO.write(record, fh, 'fasta')
        with open(phage_aggr, 'w+') as fh:
            phages = [str(phage.absolute()) for phage in phage_repaired_folder.iterdir()]
            fh.write("\n".join(phages))
        phage_sketch = self.run_sketch(phage_aggr, threads, phage_sname)
        phage_aggr.unlink()
        shutil.rmtree(phage_repaired_folder)

        return host_sketch, phage_sketch

    def run_sketch(self, file: Path, threads: int, outname: str):
        """Method for running mash sketch module

        Args:
            file (Path): File with paths to fasta files to run sketch on.
            threads (Path): Num of threads to use.
            outname (str): Prefix for sketch file.
        Returns:
            (Path): Path to produced sketch file.
        """
        output = subprocess.run(
            ['mash', 'sketch', '-k', "21", '-s', "100000", '-p', str(threads), '-o', outname, '-l', str(file)],
            capture_output=True)
        outfile = Path(outname + ".msh")
        if output.stderr and "Sketching" not in output.stderr.decode():
            raise subprocess.SubprocessError(output.stderr.decode())
        return outfile

    def run_mash(self, host_sketch: Path, phage_sketch: Path, outfile: Path):
        """Function for running mash prediction

        Args:
            host_sketch (Path): Path to host sketch file.
            phage_sketch (Path): Path to phage sketch file.
            outfile (Path): Path to results_file.
        """
        with open(outfile, 'w+') as fh:
            output = subprocess.run(['mash', 'dist', str(host_sketch), str(phage_sketch)], stdout=fh)
        if output.stderr:
            raise subprocess.SubprocessError(output.stderr.decode())
        return outfile


def main(args):
    """Main method, used with argparse.args"""
    print(utils.LOGO)

    m = Mash(Path(args.host_dir), Path(args.virus_dir))
    print("sketching host and phage files")
    host_sketch, phage_sketch = m.sketch("host_sketch", "phage_sketch", args.threads)
    print("Running mash prediction")
    outfile = m.run_mash(host_sketch, phage_sketch, Path(args.results_file))
    if args.clear_after:
        host_sketch.unlink()
        phage_sketch.unlink()
    df = pd.read_table(outfile, names=("Host", "Virus", "Distance", "P-value", "Hashes"))
    df["Host"] = df["Host"].map(lambda x: x.split("/")[-1].split(".")[0])
    df["Virus"] = df["Virus"].map(lambda x: x.split("/")[-1].split(".")[0])
    df = df.reindex(['Virus', 'Host', 'Distance', 'P-value', 'Hashes'], axis=1)
    df = df.drop(columns=["P-value", "Hashes"])

    df = df[df["Distance"] < 1]
    df["Distance"] = 1 - df["Distance"]

    df['Virus'] = df['Virus'].map(lambda x: x.split("|")[0])
    df['Host'] = df['Host'].map(lambda x: x.split("|")[0])

    # df["MashRank"] = df.groupby(["Virus"])['Distance'].rank(method='dense', ascending=True).astype(int)
    df.rename(columns={'Distance': 'Score'}, inplace=True)
    df.to_csv(outfile, index=False)

    print("Mash results: \n", df)
    print(f"Saved mash results to {str(outfile)}")
