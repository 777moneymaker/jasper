import subprocess
import unittest
import os
import pandas as pd
from pathlib import Path

from jasper import blast


class BlastTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(BlastTests, self).__init__(*args, **kwargs)
        self.test_data = Path(__file__).parent.absolute()

    def test_database_init_ok(self):
        try:
            blast.Database(source_dir=self.test_data / Path("data"), name="test_db")
        except:
            self.fail("Valid path caused exception.")

    def test_database_init_wrong(self):
        with self.assertRaises(FileNotFoundError):
            blast.Database(source_dir=self.test_data / Path("not_existing_path"), name="test_db")
            blast.Database(source_dir=self.test_data / Path("data/fasta_test_data/seqs.fasta"), name="test_db")

    def test_repair_ok(self):
        db = blast.Database(source_dir=Path(""), name="test_db")
        result = db._repair_fasta(self.test_data / Path("data/fasta_test_data/seqs.fasta"))
        expected = """>seqs|1
ATGCTGATCG
>seqs|2
GACGGTACG"""
        self.assertEqual(result, expected)

    def test_repair_wrong_file(self):
        db = blast.Database(source_dir=Path(""), name="test_db")
        with self.assertRaises(FileNotFoundError):
            db._repair_fasta(self.test_data / Path("data/fasta_test_data"))

    def test_aggregate_ok(self):
        blast.Database(source_dir=self.test_data / Path("data/fasta_test_data"), name="test_db")._aggregate()
        self.assertTrue(Path("blast_input.fasta").exists())
        self.assertTrue(os.path.getsize(Path("blast_input.fasta")) > 0)

    def test_aggregate_empty(self):
        with self.assertRaises(ValueError):
            blast.Database(source_dir=self.test_data / Path("data/fasta_test_empty_data"), name="test_db")._aggregate()

    def test_create(self):
        db, output = blast.Database(source_dir=self.test_data / Path("data/fasta_test_data"), name="test_db").create()
        self.assertTrue(map(lambda p: (self.test_data / Path(p)).exists(), "test_db.nhr", "test_db.nin", "test_db.nsq"))
        db.clear_files()

    def test_query_ok(self):
        db = blast.Database(source_dir=self.test_data / Path("data/fasta_test_data"), name="test_db")
        db.create()
        expected = pd.DataFrame(columns=['Virus', 'Host', 'Score'], data=[['seqs|1', 'seqs|1', 10], ['seqs|2', 'seqs|2', 9]])
        res = db.query(self.test_data / Path("data/fasta_test_data"), config={"task": "blastn-short"}, num_threads=1)
        self.assertTrue(expected.equals(res))

    def test_query_wrong_dir(self):
        with self.assertRaises(FileNotFoundError):
            db = blast.Database(source_dir=self.test_data / Path("data/fasta_test_data"), name="test_db")
            db.create()
            db.query(self.test_data / Path("data/fasta_test_data/seqs.fasta"), config={}, num_threads=1)

    def test_query_wrong_command(self):
        with self.assertRaises(subprocess.SubprocessError):
            db = blast.Database(source_dir=self.test_data / Path("data/fasta_test_data"), name="test_db")
            db.create()
            db.query(self.test_data / Path("data/fasta_test_data/"), blast_format="test_wrong_format", config={"task": "blastn-short"}, num_threads=1)
        db.clear_files()

    def test_clear_files(self):
        db = blast.Database(source_dir=self.test_data / Path("data/fasta_test_data"), name="test_db")
        db.create()
        db.clear_files()
        self.assertTrue(map(lambda p: not (self.test_data / Path(p)).exists(), "test_db.nhr", "test_db.nin", "test_db.nsq"))

    def tearDown(self):
        files = [
            Path("test_db.nhr"),
            Path("test_db.nin"),
            Path("test_db.nsq"),
            Path("blast_input.fasta")
        ]
        for f in files:
            if f.exists():
                f.unlink()


if __name__ == '__main__':
    unittest.main()
