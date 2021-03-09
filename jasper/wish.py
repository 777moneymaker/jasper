#!/usr/bin/env python3
"""This module manages the alignment-free operations and analysis.

This module uses Biopython and WIsH software for virus's host prediction by using alignment-free methods.

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


class Wish:
    def __init__(self, host_dir: Path, phage_dir: Path):
        """Inits Wish object with 2 directories. One phage directory and one host directory.
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

    def build(self, model_dir: Path = Path("afree_model_dir"), threads: int = os.cpu_count()):
        """This function uses WIsH software to build a model for source files.
        For more information, check WIsH documentation.

        Args:
            model_dir (Path): Path to directory in which model will be created.
            threads (int): Number of threads to pass to WIsH.
        Raises:
            SubprocessError: When WIsH returned error.
        """
        # WIsH -c build -g ~/Desktop/JASPER/example_data/host/ -m modelDir -t 3
        if not model_dir.exists():
            model_dir.mkdir()

        temp_input_dir = Path("temp_input_dir")
        if not temp_input_dir.exists():
            temp_input_dir.mkdir()

        for i, host_file in enumerate(self.host_dir.iterdir(), 1):
            repaired = blast.Database.repair_fasta(host_file)
            for record in SeqIO.parse(io.StringIO(repaired), 'fasta'):
                with open(temp_input_dir / Path(f"{record.id}.fasta"), 'w+') as fh:
                    SeqIO.write(record, fh, 'fasta')

        output = subprocess.run(
            ['WIsH', '-c', 'build', '-g', str(temp_input_dir), '-m', str(model_dir), '-t', str(threads)],
            capture_output=True
        )

        shutil.rmtree(temp_input_dir)

        if output.stderr.decode():
            raise subprocess.SubprocessError(output.stderr.decode())

    def predict(self, model_dir: Path = Path("afree_model_dir"), output_dir: Path = Path("afree_output_dir"),
                threads: int = os.cpu_count()):
        """This function uses WIsH software to predict host for a phage by using a-free methods and precomputed model.
        For more information, check WIsH documentation.

        Args:
            model_dir (Path): Path to directory containing model.
            output_dir (Path): Path to directory in which output from WIsH will be placed. Note: it's temporary dir.
            threads (int): Number of threads to pass to WIsH.
        Raises:
            SubprocessError: When WIsH returned error.
            FileNotFoundError: When directory containing model does not exist.
        Returns:
            (pd.DataFrame): Dataframe (in form of matrix) containing log-likelihood results from WIsH.
        """
        # WIsH -c predict -g ~/Desktop/JASPER/example_data/virus/ -m modelDir -r outputResultDir -b
        if not model_dir.exists():
            raise FileNotFoundError("Model directory does not exist.")

        if not output_dir.exists():
            output_dir.mkdir()

        temp_input_dir = Path("temp_input_dir")
        if not temp_input_dir.exists():
            temp_input_dir.mkdir()

        res = pd.DataFrame({'Host': [], 'Virus': [], 'Score': []})
        for phage_file in self.phage_dir.iterdir():
            repaired = blast.Database.repair_fasta(phage_file)
            for record in SeqIO.parse(io.StringIO(repaired), 'fasta'):
                xtracted = temp_input_dir / Path(f"{record.id}.fasta")
                with open(xtracted, 'w+') as fh:
                    SeqIO.write(record, fh, 'fasta')
                output = subprocess.run(
                    ['WIsH', '-c', 'predict', '-g', str(temp_input_dir), '-m', str(model_dir), '-r',
                     str(output_dir),
                     '-b',
                     '-t', str(threads)], capture_output=True)
                if output.stderr.decode():
                    raise subprocess.SubprocessError(output.stderr.decode())
                xtracted.unlink()
                df = pd.read_csv(output_dir / Path("llikelihood.matrix"),
                                 sep='\t', index_col=0)
                tuplz = df.stack().reset_index().agg(tuple, 1).to_list()
                df = pd.DataFrame(tuplz, columns=['Host', 'Virus', 'Score'])
                res = pd.concat([res, df])
        shutil.rmtree(temp_input_dir)
        return res


def main(args):
    """Main function for running module.
    Args:
        args (argparse obj): Arguments from right subparser.
    """
    print(utils.LOGO)
    if args.host_dir is None:
        args.host_dir = "."

    w = Wish(Path(args.host_dir), Path(args.virus_dir))

    if not args.model_dir:
        print("Building model...")
        w.build(threads=args.threads)
        print('Done.', end='\n\n')

    print("Predicting...")
    if args.model_dir:
        df = w.predict(model_dir=Path(args.model_dir), threads=args.threads)
    else:
        df = w.predict(threads=args.threads)
    shutil.rmtree(Path("afree_output_dir"))

    df = df.reindex(['Virus', 'Host', 'Score'], axis=1)

    df['Virus'] = df['Virus'].map(lambda x: x.split("|")[0])
    df['Host'] = df['Host'].map(lambda x: x.split("|")[0])

    df.to_csv(Path(args.results_file), index=False)
    print("Done.")

    if args.clear_after:
        shutil.rmtree(Path("afree_model_dir"))
    print("Saved results to", args.results_file)
