#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline

from . import blast
from . import utils


class tRNAScanner:
    def __init__(self, source_dir: Path, target_dir: Path, name: str) -> None:
        self.source_dir = source_dir
        self.target_dir = target_dir
        self.name = name

    def scan(self, source_dir: Path, mode: str = 'general', threads: int = os.cpu_count()):
        aggregate = Path("trna_aggregate.fasta")
        for i, source_file in enumerate(source_dir.iterdir(), 1):
            if not source_file.name.endswith(utils.TYPES):
                continue

            with open(aggregate, 'a+') as aggregate_fh:
                aggregate_fh.write(blast.Database.repair_fasta(source_file))
            print(f"Performed {i} files aggregation for tRNA-scan", end='\r')

        print()
        print("\nRunning tRNAscan-SE")
        out_file = Path("trnascan_out.fasta")
        trna_output = self.run_trnascan(aggregate, out_file, threads, mode)
        print(trna_output.stdout.decode())
        aggregate.unlink()

        with open(out_file, 'r') as trna_fh:
            trnas = [(str(record.seq), record.name) for record in SeqIO.parse(trna_fh, 'fasta')]

        out_file.unlink()
        return trnas

    def run_trnascan(self, file: Path, out_file: Path, num_threads: int, mode: str):
        scan_mode = {'bacterial': "-B", "phage": "-P", "general": "-G"}[mode]

        res = subprocess.run(['tRNAscan-SE', str(file), '-a', str(out_file), '--thread', str(num_threads), scan_mode],
                             capture_output=True)
        return res

    def create(self):
        '''
        TODO: FIX THAT
        '''

        try:
            cmd = NcbimakeblastdbCommandline(input_file="trna_blast_input.fasta",
                                             dbtype="nucl",
                                             title=self.name,
                                             out=self.name)
            makeblastdb_output = subprocess.run(str(cmd), capture_output=True)  # shell=True
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


if __name__ == '__main__':
    print(utils.LOGO)

    t = tRNAScanner(source_dir=Path("example_data/host"), target_dir=Path("example_data/virus"), name="trna_db")
    print("Aggregation and scanning host files...")
    host_trnas = t.scan(t.source_dir, mode="bacterial")
    print("Done\n")

    print("Aggregation and scanning phage files...")
    vir_trnas = t.scan(t.target_dir, mode="phage")
    print("Done\n")

    print(*vir_trnas, sep="\n")
    print(*host_trnas, sep="\n")
