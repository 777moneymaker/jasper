#!/usr/bin/env python3
import io
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from . import utils


class Wish:
    def __init__(self, host_dir: Path, phage_dir: Path):
        if not host_dir.exists():
            raise TypeError("Host directory does not exist.")
        if not host_dir.is_dir():
            raise TypeError("Host directory is not a directory.")
        if not phage_dir.exists():
            raise TypeError("Phage directory does not exist.")
        if not phage_dir.is_dir():
            raise TypeError("Phage directory is not a directory.")

        self.host_dir = host_dir
        self.phage_dir = phage_dir

    def _repair_fasta(self, file: Path):
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

    def build(self, model_dir: Path = Path("temp_model_dir"), threads: int = os.cpu_count()):
        # WIsH -c build -g ~/Desktop/JASPER/example_data/host/ -m modelDir -t 3

        if not model_dir.exists():
            model_dir.mkdir()

        temp_input_dir = Path("temp_input_dir")
        if not temp_input_dir.exists():
            temp_input_dir.mkdir()

        for i, host_file in enumerate(self.host_dir.iterdir(), 1):
            print(f"host concat {i}", end='\r')
            repaired = self._repair_fasta(host_file)
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

    def predict(self, model_dir: Path = Path("temp_model_dir"), output_dir: Path = Path("temp_output_dir"), threads: int = os.cpu_count()):
        # WIsH -c predict -g ~/Desktop/JASPER/example_data/virus/ -m modelDir -r outputResultDir -b
        if not model_dir.exists():
            raise FileNotFoundError("Model directory does not exist.")
        if not output_dir.exists():
            output_dir.mkdir()

        temp_input_dir = Path("temp_input_dir")
        if not temp_input_dir.exists():
            temp_input_dir.mkdir()

        for i, phage_file in enumerate(self.phage_dir.iterdir(), 1):
            print(f"phage concat {i}", end='\r')
            repaired = self._repair_fasta(phage_file)
            for record in SeqIO.parse(io.StringIO(repaired), 'fasta'):
                with open(temp_input_dir / Path(f"{record.id}.fasta"), 'w+') as fh:
                    SeqIO.write(record, fh, 'fasta')

        output = subprocess.run(
            ['WIsH', '-c', 'predict', '-g', str(temp_input_dir), '-m', str(model_dir), '-r', str(output_dir), '-b',
             '-t', str(threads)],
            capture_output=True)
        shutil.rmtree(temp_input_dir)

        if output.stderr.decode():
            raise subprocess.SubprocessError(output.stderr.decode())

        df = pd.read_csv(output_dir / Path("prediction.list"),
                         sep=None, engine='python').iloc[:, 0:3]
        df.columns = pd.Index(['Virus', 'Host', 'Log'])
        return df


def main(args):
    print(utils.LOGO)
    w = Wish(Path(args.host_dir), Path(args.virus_dir))

    if not args.model_dir:
        print("Building model...")
        w.build(threads=args.threads)
        print('Done.')

    print("Predicting...")
    if args.model_dir:
        df = w.predict(model_dir=Path(args.model_dir), threads=args.threads)
    else:
        df = w.predict(threads=args.threads)
    df.to_csv(Path(args.results_file), index=False)
    print("Done.")
    shutil.rmtree(Path("temp_output_dir"))

    if args.clear_after:
        shutil.rmtree(Path("temp_model_dir"))

    print("Saved results to", args.results_file)
