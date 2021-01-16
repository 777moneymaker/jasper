#!/usr/bin/env python3
"""This module manages the crispr parsing and analysis.

This module uses Biopython and Piler-CR for identification of CRISPR spacers and their analysis.

More about it:
    1. https://biopython.org/

    2. https://biopython.org/docs/dev/api/Bio.Blast.Applications.html

    3. http://www.ncbi.nlm.nih.gov/pubmed/17239253
"""

from __future__ import annotations

import os
import shutil
import subprocess
import pandas as pd
from pathlib import Path

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from . import blast
from . import utils


class CrisprFinder(blast.Database):
    def __init__(self, source_dir: Path, name: str) -> None:
        """Inits obj with args"""
        super(CrisprFinder, self).__init__(source_dir, name)

    def _make_output_dir(self, directory: Path):
        """This function creates output dir for PILER-CR

        Args:
            directory (Path): Directory to be created.
        Raises:
            TypeError: When given path is not Path obj.
        """
        if not isinstance(directory, Path):
            raise TypeError("Given object is not Path object.")

        if not directory.exists():
            directory.mkdir()

    def retrieve_spacers(self) -> CrisprFinder:
        """This function retrieves spacers from files in source directory

        This method uses blast module functionality.
        It runs PILER-CR and parses output for each file in source dir.
        Then it creates a virus database and uses retrieved spacers as a blastn-short query.
        Resulting dataframe is populated with blast query results.
        """
        res_dir = Path("crispr_spacers")
        self._make_output_dir(res_dir)

        for host in self.source_dir.iterdir():
            if not host.name.endswith(utils.TYPES):
                continue

            repaired_file = Path(f"{host.stem}.repaired")
            with open(repaired_file, "w+") as repaired_fh:
                repaired_fh.write(self._repair_fasta(host))

            piler_file = res_dir / Path(f"{host.stem}.piler")  # Output file for PILERCR.
            self.find_crispr_spacers(repaired_file, piler_file)  # Find crispr spacers for given file.
            repaired_file.unlink()  # Remove temp file for repaired seq

            # Process the PILERCR output
            with open(piler_file, 'r') as res_fh:
                # Read the content of PILERCR file.
                piler_content = res_fh.read()
                # If spacer not found then continue
                if "DETAIL REPORT" not in piler_content:
                    piler_file.unlink()
                    continue
                # Get only spacers region
                piler_content = piler_content.split("DETAIL REPORT")[1].split("SUMMARY BY SIMILARITY")[0]
            # Delete PILERCR output file
            piler_file.unlink()

            for line in piler_content.splitlines():
                if line.startswith(">"):
                    name = f"{line.lstrip('>').split(' ')[0]}"
                try:
                    line = list(filter(lambda x: x, line.split(" ")))  # Filter non empty strings
                    if len(line) == 7 and any(base in line[6] for base in "ATGC"):
                        # If line had 7 fields and any of 'ATGC' in 7th field
                        seq = SeqRecord(Seq(line[6]), id=name, description="", name=name)  # Create seq
                        spacers_file = res_dir / Path(f"{name}.fasta")  # Result file for spacers
                        with open(spacers_file, 'a+') as final_fh:  # Append seq to file
                            SeqIO.write(seq, final_fh, 'fasta')
                except (ValueError, IndexError):
                    continue
        return self

    def find_crispr_spacers(self, host_file: Path, out_file: Path):
        """This function runs PILER-CR for a given host file and produces output to given file

        Args:
            host_file (Path): Host file to be passed as PILER input.
            out_file (Path): Output file to be created by PILER.
        Raises:
            TypeError: When given obj is of wrong type.
            FileNotFoundError: When given host_file is not file or doesn't exist.
            SubprocessError: When PILER-CR returned error.
        """
        if not isinstance(host_file, Path):
            raise TypeError("Given file is not Path obj.")
        if not isinstance(out_file, Path):
            raise TypeError("Given outfile is not a Path obj.")
        if not host_file.exists():
            raise FileNotFoundError("Host file does not exist.")
        if not host_file.is_file():
            raise FileNotFoundError("Given file is not a file.")

        try:
            subprocess.run(['pilercr', '-in', str(host_file), '-out', str(out_file)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            raise subprocess.SubprocessError("PilerCR returned error. Check your input", e.output)


def main(args):
    """Main function for running module.
    Args:
        args (argparse obj): Arguments from right subparser.
    """
    print(utils.LOGO)
    args.short_config = utils.parse_config(args.short_config, {
        "task": "blastn-short",
        "num_threads": os.cpu_count(),
        "evalue": 1,
        "gapopen": 10,
        "gapextend": 2,
        "penalty": -1,
        "word_size": 7,
        "dust": "no",
    })

    print("Starting analysis...")
    print("Aggregating files, retrieving crispr spacers...")
    finder = CrisprFinder(Path(args.host_dir), "-")
    finder.retrieve_spacers()

    if args.create_db_name:
        vir_db, vir_db_output = blast.Database(Path(args.virus_dir), args.create_db_name).create()
        print(vir_db_output)
    else:
        vir_db = blast.Database(Path('.'), args.use_db_name)

    print("Quering...")
    query_df = vir_db.query(Path("crispr_spacers/"),
                            config=args.short_config,
                            blast_format="10 qseqid sseqid score qlen length mismatch gaps",
                            headers=("Spacer", "Virus", "Score", "Qlen", "Alen", "Mis", "Gap"))
    if args.clear_after:
        vir_db.clear_files()
    shutil.rmtree(Path("crispr_spacers/"))

    query_df[["Score", "Qlen", "Alen", "Mis", "Gap"]].apply(pd.to_numeric)
    query_df['Allowed'] = query_df['Qlen'] - (
            query_df['Alen'] - query_df['Mis'] - query_df['Gap'])
    query_df['Allowed'] = query_df['Allowed'].apply(pd.to_numeric)
    short_results = query_df.drop(columns=['Qlen', 'Alen', 'Mis', 'Gap'])
    short_results = short_results[short_results['Allowed'] <= args.allowed_mis].drop(columns='Allowed').reset_index(drop=True)

    short_results.rename(columns={"Spacer": "Host"}, inplace=True)
    short_results = short_results.reindex(["Virus", "Host", "Score"], axis=1)
    short_results["Host"] = short_results["Host"].map(lambda x: "|".join(x.split("|")[:2]))
    short_results = short_results.reset_index(drop=True)

    # short_results = short_results.groupby(["Virus", "Spacer"]).sum().reset_index()
    # idx = short_results.groupby(['Virus'])['Score'].idxmax()
    # crispr_results: pd.DataFrame = short_results.loc[idx].sort_values('Score', ascending=False).reset_index(drop=True)

    print("blastn-short results (vir_genome-spacers query): ", short_results, sep='\n')

    short_results.to_csv(args.output_file, index=False)
    print("Saved files to", args.output_file)


class Expando(object):
    pass


if __name__ == '__main__':
    args = Expando()
    args.virus_dir = Path("example_data/virus")
    args.host_dir = Path("example_data/host")
    args.create_db_name = "vir_db"
    args.short_config = ""
    args.clear_after = True
    args.output_file = "crispr_results.csv"
    args.allowed_mis = 1
    main(args)
