from pathlib import Path
import subprocess
from jasper import blast

from Bio import SeqIO
from jasper import utils


class tRNAScanner(blast.Database):
    def __init__(self, source_dir: Path, name: str) -> None:
        super(tRNAScanner, self).__init__(source_dir, name)

    def scan(self):
        aggregate_file = Path("temp.fasta", "a+")
        for i, source_file in self.source_dir.iterdir():
            if not source_file.name.endswith(utils.TYPES):
                continue

            with open(aggregate_file, 'r') as aggregate_fh:
                aggregate_fh.write(self._repair_fasta(source_file))
            print(f"Performed {i} files aggregation for tRNA-scan", end='\r')

        trna_out = Path("trnas.fasta")
        print("Running tRNAscan-SE")
        self.run_trnascan(aggregate_file, trna_out)
        aggregate_file.unlink()

        with open(trna_out, 'r') as trna_fh:
            trnas = [(record.seq, record.name) for record in SeqIO.parse(trna_fh, 'fasta')]
        trna_out.unlink()
        return trnas

    def run_trnascan(self, file: Path, out_file: Path):
        res = subprocess.run(['tRNAscan-SE', str(file), '-F', str(out_file)],
                             stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        return res
