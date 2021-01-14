from __future__ import annotations

import shutil
import subprocess
import time
import pandas as pd
from pathlib import Path
import argparse
from collections import defaultdict

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from . import database
from . import utils


class CrisprFinder(database.Database):
    def __init__(self, source_dir: Path, name: str) -> None:
        super(CrisprFinder, self).__init__(source_dir, name)

    def retrieve_spacers(self) -> CrisprFinder:
        results_dir = Path("crispr_spacers")  # Make directory for PILERCR output.
        if not results_dir.exists():
            results_dir.mkdir()

        for i, host in enumerate(self.source_dir.iterdir(), 1):
            if not host.name.endswith(utils.TYPES):
                continue

            repaired_file = Path(f"{host.stem}.repaired")  # Temp file for repaired seq
            with open(repaired_file, "w+") as repaired_fh:
                repaired_fh.write(self._repair_fasta(host))

            piler_output_file = results_dir / Path(f"{host.stem}.piler")  # Output file for PILERCR.
            piler_output = self.find_crispr_spacers(repaired_file, piler_output_file)  # Find crispr spacers for given file.
            repaired_file.unlink()  # Remove temp file for repaired seq

            # Process the PILERCR output
            with open(piler_output_file, 'r') as res_fh:
                # Read the content of PILERCR file.
                piler_content = res_fh.read()
                # If spacer not found then continue
                if "DETAIL REPORT" not in piler_content:
                    piler_output_file.unlink()
                    continue
                # Get only spacers region
                piler_content = piler_content.split("DETAIL REPORT")[1].split("SUMMARY BY SIMILARITY")[0]
            # Delete PILERCR output file
            piler_output_file.unlink()

            for line in piler_content.splitlines():
                if line.startswith(">"):
                    name = f"{line.lstrip('>').split(' ')[0]}"
                try:
                    line = list(filter(lambda x: x, line.split(" ")))  # Filter non empty strings
                    if len(line) == 7 and any(base in line[6] for base in "ATGC"):
                        # If line had 7 fields and any of 'ATGC' in 7th field
                        seq = SeqRecord(Seq(line[6]), id=name, description="", name=name)  # Create seq
                        spacers_file = results_dir / Path(f"{name}.fasta")  # Result file for spacers
                        with open(spacers_file, 'a+') as final_fh:  # Append seq to file
                            SeqIO.write(seq, final_fh, 'fasta')
                except (ValueError, IndexError):
                    continue
        return self

    def find_crispr_spacers(self, host_file: Path, out_file: Path):
        try:
            return subprocess.run(['pilercr', '-in', str(host_file), '-out', str(out_file)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            raise subprocess.SubprocessError("PilerCR returned error. Check your input", e.output)


def main(args):
    print(utils.LOGO)
    args.short_config = utils.parse_config(args.short_config, {
        "task": "blastn-short",
        "evalue": 1,
        "gapopen": 10,
        "gapextend": 2,
        "penalty": -1,
        "word_size": 7,
        "dust": "no",
    })

    print("Starting analysis...")
    print("Aggregating files, retrieving crispr spacers...")
    crispr_finder = CrisprFinder(Path(args.host_dir), "-").retrieve_spacers()

    if args.create_db_name:
        vir_db, vir_db_output = database.Database(args.virus_dir, args.create_db_name).create()
        print(vir_db_output)
    else:
        vir_db = database.Database(Path('.'), args.use_db_name)

    print("Quering...")
    query_df = vir_db.query(Path("crispr_spacers/"),
                            config=args.short_config,
                            blast_format="10 qseqid sseqid score qlen length mismatch gaps",
                            headers=["Spacer", "Virus", "Score", "Qlen", "Alen", "Mis", "Gap"])
    shutil.rmtree(Path("crispr_spacers/"))

    query_df[["Score", "Qlen", "Alen", "Mis", "Gap"]].apply(pd.to_numeric)
    query_df['Allowed'] = query_df['Qlen'] - (
            query_df['Alen'] - query_df['Mis'] - query_df['Gap'])
    query_df['Allowed'] = query_df['Allowed'].apply(pd.to_numeric)
    short_results = query_df.drop(columns=['Qlen', 'Alen', 'Mis', 'Gap'])

    short_results = short_results[short_results['Allowed'] <= args.allowed_mis].drop(columns='Allowed').reset_index()
    short_results.reindex(["Virus", "Spacer", "Score"], axis=1)
    short_results["Spacer"] = short_results["Spacer"].map(lambda x: x.split("|")[0])

    short_results = short_results.groupby(["Virus", "Spacer"]).sum().reset_index()
    idx = short_results.groupby(['Virus'])['Score'].idxmax()
    crispr_results: pd.DataFrame = short_results.loc[idx].sort_values('Score', ascending=False).reset_index(drop=True)
    if args.clear_after:
        vir_db.clear_files()
    print("blastn-short results (vir_genome-spacers query): ", crispr_results, sep='\n')

    crispr_results.to_csv(args.output_file, index=False)
    print("Saved files to", args.output_file)
