import subprocess
import time
import pandas as pd
from pathlib import Path
from collections import defaultdict

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from jasper import database


class Crispr(database.Database):
    def __init__(self, config: dict) -> None:
        super(Crispr, self).__init__(config)

    def create(self) -> dict:
        results = defaultdict(list)
        for i, host_file in enumerate(self.host_dir.iterdir(), 1):
            # TODO: Remove this line. Only for tests.
            if i == 10:
                break

            # Repair files if user want's to
            if self._repair_host_files:
                repaired_content = self._repair_fasta(host_file)
            else:
                repaired_content = self._read_fasta(host_file)

            name = f"{str(host_file)}.repaired{host_file.suffix}"  # Temp file for repaired seq
            with open(name, "w+") as repaired_fh:
                repaired_fh.write(repaired_content)

            res_dir = Path("crispr_spacers")  # Make directory for PILERCR output.
            if not res_dir.exists():
                res_dir.mkdir()

            out_name = res_dir / Path(f"{host_file.stem}.txt")  # Output file for PILERCR.
            host_file = Path(name)  # Host file from host dir.

            self.find_crispr_spacers(host_file, out_name)  # Find crispr spacers for given file.
            Path(name).unlink()  # Remove temp file for repaired seq
            print(f'Processed {i} host files for CRISPR identification', end='\r')

            # Process the PILERCR output
            with open(out_name, 'r') as res_fh:
                content = res_fh.read()  # Read the content of PILERCR file.
                if "DETAIL REPORT" not in content:  # If spacer not found then continue
                    continue

                content = content.split("DETAIL REPORT")[1].split("SUMMARY BY SIMILARITY")[0]  # Get only spacers region
                for line in content.splitlines():
                    if line.startswith(">"):
                        name = f"{line.lstrip('>').split(' ')[0]}"

                    try:
                        line = list(filter(lambda x: x, line.split(" ")))  # Filter non empty strings

                        # If line had 7 fields and any of 'ATGC' in last field
                        if len(line) == 7 and any(base in line[6] for base in "ATGC"):
                            results[name].append(line[6])  # Append spacer
                            res_fasta = res_dir / Path(f"{out_name.stem}.fasta")  # Open handle.
                            seq = SeqRecord(Seq(line[6]), id=name, description="", name=name)  # Create seq
                            with open(res_fasta, 'a+') as final_fh:  # Append seq to file
                                SeqIO.write(seq, final_fh, 'fasta')
                    except (ValueError, IndexError):
                        continue
            out_name.unlink()  # Delete PILERCR output file
        print()
        return dict(results)

    def find_crispr_spacers(self, host_file: Path, out_file: Path):
        result = subprocess.run(['pilercr', '-in', str(host_file), '-out', str(out_file)],
                                stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        return result
