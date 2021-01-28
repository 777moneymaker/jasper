import unittest
from pathlib import Path

from Bio import SeqIO

from jasper import trna


class tRNATests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(tRNATests, self).__init__(*args, **kwargs)
        self.test_dir = Path(__file__).parent.absolute()

    def test_scanner_init_ok(self):
        try:
            trna.tRNAScanner(self.test_dir / Path("data/test_trnascan_host"),
                             self.test_dir / Path("data/test_trnascan_phage"),
                             "trna_db")
        except:
            self.fail("Valid path caused exception.")

    def test_scanner_init_wrong_path(self):
        with self.assertRaises(FileNotFoundError):
            trna.tRNAScanner(self.test_dir / Path("data/not_existing"),
                             self.test_dir / Path("not_existing"),
                             "trna_db")
            trna.tRNAScanner(self.test_dir / Path("data/blast_query_data/seqs.fasta"),
                             self.test_dir / Path("data/blast_query_data/seqs.fasta"),
                             "trna_db")

    def test_scanner_init_wrong_instance(self):
        with self.assertRaises(TypeError):
            trna.tRNAScanner("not path",
                             self.test_dir / Path("data/test_trnascan_phage"),
                             "trna_db")
            trna.tRNAScanner(self.test_dir / Path("data/test_trnascan_host"),
                             "not path",
                             "trna_db")
            trna.tRNAScanner(self.test_dir / Path("data/test_trnascan_host"),
                             self.test_dir / Path("data/test_trnascan_phage"),
                             Path("not a str"))

    def test_scanner_scans(self):
        t = trna.tRNAScanner(self.test_dir / Path("data/test_trnascan_host"),
                             self.test_dir / Path("data/test_trnascan_phage"),
                             "trna_db")
        host = t.scan(t.source_dir, output_filename="test_host_trnas.fasta")
        phage = t.scan(t.target_dir, output_filename="test_phage_trnas.fasta")
        self.assertTrue(host.exists())
        self.assertTrue(phage.exists())
        self.assertTrue(len(list(seq for seq in SeqIO.parse(host, 'fasta'))) > 0)
        host.unlink()
        phage.unlink()

    def test_scanner_creates(self):
        t = trna.tRNAScanner(self.test_dir / Path("data/test_trnascan_host"),
                             self.test_dir / Path("data/test_trnascan_phage"),
                             "trna_db")
        host = t.scan(t.source_dir, output_filename="test_host_trnas.fasta")
        t.create(host)
        self.assertTrue(all(map(lambda x: Path(x).exists(), ['trna_db.nhr', 'trna_db.nin', 'trna_db.nsq'])))
        t.clear_files()

    def test_scanner_query(self):
        t = trna.tRNAScanner(self.test_dir / Path("data/test_trnascan_host"),
                             self.test_dir / Path("data/test_trnascan_phage"),
                             "trna_db")
        host = t.scan(t.source_dir, output_filename="test_host_trnas.fasta")
        phage = t.scan(t.target_dir, output_filename="test_phage_trnas.fasta")
        t.create(host)
        res = t.query(phage, {"task": "blastn-short"}, blast_format="10 qseqid sseqid score",
                      headers=("Virus", "Host", "Score"))
        self.assertTrue(len(res) == 157)
        t.clear_files()


if __name__ == "__main__":
    unittest.main()
