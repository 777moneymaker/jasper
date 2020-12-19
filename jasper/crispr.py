import subprocess
import time
import pandas as pd
from pathlib import Path

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from jasper import database


class Crispr(database.Database):
    def __init__(self, host_path: str, db_name: str = "host_db", repair_host_files: bool = True,
                 repair_vir_files: bool = True):
        super(Crispr, self).__init__(host_path, db_name, repair_host_files, repair_vir_files)

    def create(self):
        results = []
        for i, host_file in enumerate(self.source_dir.iterdir(), 1):
            if i == 20:
                break

            if self._repair_host_files:
                repaired_content = self._repair_fasta(host_file)
            else:
                repaired_content = self._read_fasta(host_file)
            name = f"{str(host_file)}.repaired{host_file.suffix}"
            with open(name, "w+") as repaired_fh:
                repaired_fh.write(repaired_content)
            res_dir = Path("crispr_spacers")
            if not res_dir.exists():
                res_dir.mkdir()
            out_name = res_dir / Path(f"{host_file.stem}.txt")
            host_file = Path(name)
            self.find_crispr_spacers(host_file, out_name)
            Path(name).unlink()
            print(f'Processed {i} host files for CRISPR identification')

            with open(out_name, 'r') as res_fh:
                content = res_fh.read()
                if "DETAIL REPORT" not in content:
                    results.append((out_name.stem, None))
                else:
                    content = content.split("DETAIL REPORT")[1].split("SUMMARY BY SIMILARITY")[0]
                    print(content)
                    i = 1
                    for line in content.splitlines():
                        if line.startswith(">"):
                            name = f"{line.lstrip('>').split('|')[0]}|{i}"
                            i += 1
                        try:
                            line = list(filter(lambda x: x, line.split(" ")))
                            if len(line) == 7 and any(base in line[6] for base in "ATGC"):
                                results.append((name, line[6]))
                                res_fasta = res_dir / Path(f"{out_name.stem}.fasta")
                                seq = SeqRecord(Seq(line[6]), id=name, description="", name=name)
                                with open(res_fasta, 'a+') as final_fh:
                                    SeqIO.write(seq, final_fh, 'fasta')
                        except ValueError:
                            continue
            out_name.unlink()
        print(*results, sep='\n')

    def find_crispr_spacers(self, host_file: Path, out_file: Path):
        result = subprocess.run(['pilercr', '-in', str(host_file), '-out', str(out_file)], stdout=subprocess.DEVNULL)
        return result
