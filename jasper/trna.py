#!/usr/bin/env python3
import io
import os
import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline

from . import blast
from . import utils


class tRNAScanner:
    def __init__(self, source_dir: Path, target_dir: Path, name: str) -> None:
        self.source_dir = source_dir
        self.target_dir = target_dir
        self.name = name

    def clear_files(self):
        """Function clears database files if any present."""
        for file in Path(".").iterdir():
            if file.stem == self.name and file.suffix in (".nhr", ".nin", ".nsq"):
                file.unlink()

    # max_iter for tests
    def scan(self, source_dir: Path, output_filename: str, mode: str = 'general', threads: int = os.cpu_count(), max_iter: int = 10):
        aggregate = Path("trna_aggregate.fasta")
        for i, source_file in enumerate(source_dir.iterdir(), 1):
            if i == max_iter: break
            if not source_file.name.endswith(utils.TYPES):
                continue

            with open(aggregate, 'a+') as aggregate_fh:
                aggregate_fh.write(blast.Database.repair_fasta(source_file))
            print(f"Performed {i} files aggregation for tRNA-scan", end='\r')

        print("\nRunning tRNAscan-SE")
        out_file = Path(output_filename)
        trna_output = self.run_trnascan(aggregate, out_file, threads, mode)
        aggregate.unlink()

        # with open(out_file, 'r') as trna_fh:
        #     trnas = [(str(record.seq), record.name) for record in SeqIO.parse(trna_fh, 'fasta')]
        #
        # out_file.unlink()
        # return trnas
        return out_file

    def run_trnascan(self, file: Path, out_file: Path, num_threads: int, mode: str):
        scan_mode = {'bacterial': "-B", "phage": "-P", "general": "-G"}[mode]

        res = subprocess.run(['tRNAscan-SE', str(file), '-a', str(out_file), '--thread', str(num_threads), scan_mode],
                             capture_output=True)
        return res


    def create(self, input_file: Path):
        '''
        TODO: TEST THIS
        '''

        try:
            cmd = NcbimakeblastdbCommandline(input_file=str(input_file),
                                             dbtype="nucl",
                                             title=self.name,
                                             out=self.name)
            cmd()
            makeblastdb_output = subprocess.run(str(cmd), capture_output=True, shell=True)
            if makeblastdb_output.stderr:
                if not input_file.exists():
                    raise subprocess.SubprocessError(
                        "Makeblastdb returned error. Input file for Makeblastdb does not exist. Check your input.",
                        str(cmd))
                else:
                    raise subprocess.SubprocessError("Makeblastdb returned error. Check your input.",
                                                     makeblastdb_output.stderr.decode())
        except Exception:
            raise
        finally:
            if input_file.exists():
                input_file.unlink()
        return makeblastdb_output.stdout.decode()

    def query(self, query_file: Path, config: dict, blast_format: str, headers: tuple):
        try:
            cmd = NcbiblastnCommandline(query=str(query_file),
                                        db=self.name,
                                        outfmt=blast_format,
                                        **config)
            blast_output = subprocess.run(str(cmd), capture_output=True, shell=True)

            # Error only occurs if it's not this stupid warning.
            if blast_output.stderr and "Examining 5 or more matches" not in blast_output.stderr.decode():
                if not Path("blast_query.fasta").exists():
                    raise subprocess.SubprocessError(
                        "Blastn returned error. Input file for Blastn does not exist. Check your input.", str(cmd))
                else:
                    raise subprocess.SubprocessError("Blastn returned error. Check your input.",
                                                     blast_output.stderr.decode())
            print(blast_output.stdout)
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
    host = t.scan(t.source_dir, output_filename="host_trnas.fasta", mode="general", threads=6, max_iter=10)
    print("Done\n")

    print("Aggregation and scanning phage files...")
    vir = t.scan(t.target_dir, output_filename="phage_trnas.fasta", mode="general", threads=6, max_iter=60)
    print("Done\n")

    print("Making a trna host database...")
    db_output = t.create(host)
    print(db_output)

    print("Quering a phage trnas")
    results_df = t.query(vir, {"task": "blastn-short"}, blast_format="10 qseqid sseqid score", headers=("Virus", "Host", "Score"))
    print("Done")
    print(results_df)

    results_df.to_csv(args.output_file, index=False)
    print("Saved trnas results")
    t.clear_files()


if __name__ == '__main__':
    main()
