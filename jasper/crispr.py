from __future__ import annotations

import subprocess
import time
import pandas as pd
from pathlib import Path
from collections import defaultdict

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from jasper import database

from jasper.database import TYPES


class CrisprFinder(database.Database):
    def __init__(self, source_dir: Path, name: str) -> None:
        super(CrisprFinder, self).__init__(source_dir, name)

    def retrieve_spacers(self) -> CrisprFinder:
        results_dir = Path("crispr_spacers")  # Make directory for PILERCR output.
        if not results_dir.exists():
            results_dir.mkdir()

        start: float = time.time()
        for i, host in enumerate(self.source_dir.iterdir(), 1):
            if not host.name.endswith(TYPES):
                continue

            # Remove this line. Only for tests.
            # if i == 100:
            #     break

            # Repair files if user want's to
            repaired_file = Path(f"{host.stem}.repaired")  # Temp file for repaired seq
            with open(repaired_file, "w+") as repaired_fh:
                repaired_fh.write(self._repair_fasta(host))

            piler_output_file = results_dir / Path(f"{host.stem}.piler")  # Output file for PILERCR.
            self.find_crispr_spacers(repaired_file, piler_output_file)  # Find crispr spacers for given file.
            repaired_file.unlink()  # Remove temp file for repaired seq
            print(f'Processed {i} host files for CRISPR identification', end='\r')

            # Process the PILERCR output
            with open(piler_output_file, 'r') as res_fh:
                piler_content = res_fh.read()  # Read the content of PILERCR file.
                if "DETAIL REPORT" not in piler_content:  # If spacer not found then continue
                    piler_output_file.unlink()
                    continue
                piler_content = piler_content.split("DETAIL REPORT")[1].split("SUMMARY BY SIMILARITY")[
                    0]  # Get only spacers region
            piler_output_file.unlink()  # Delete PILERCR output file

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
        end: float = time.time()
        print(f"Crispr finding including repairment: {end - start :.2f}")
        return self

    def find_crispr_spacers(self, host_file: Path, out_file: Path):
        subprocess.run(['pilercr', '-in', str(host_file), '-out', str(out_file)],
                       stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
