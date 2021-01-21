#!/usr/bin/env python3
import os
import subprocess
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd

from jasper import blast


class BlastTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(BlastTests, self).__init__(*args, **kwargs)
        self.test_dir = Path(__file__).parent.absolute()

    def test_database_init_ok(self):
        try:
            blast.Database(source_dir=self.test_dir / Path("data"), name="test_db")
        except:
            self.fail("Valid path caused exception.")

    def test_database_init_wrong_path(self):
        with self.assertRaises(FileNotFoundError):
            blast.Database(source_dir=self.test_dir / Path("not_existing_path"), name="test_db")
            blast.Database(source_dir=self.test_dir / Path("data/blast_query_data/seqs.fasta"), name="test_db")

    def test_database_init_wrong_instance(self):
        with self.assertRaises(TypeError):
            blast.Database(source_dir="not a path object", name="")
            blast.Database(source_dir=Path(""), name="not a str")

    def test_repair_ok(self):
        db = blast.Database(source_dir=Path(""), name="test_db")
        result = db.repair_fasta(self.test_dir / Path("data/blast_query_data/seqs.fasta"))
        expected = """>seqs|1
ATGCTGATCG
>seqs|2
GACGGTACG"""
        self.assertEqual(result, expected)

    def test_repair_wrong_file(self):
        db = blast.Database(source_dir=Path(""), name="test_db")
        with self.assertRaises(FileNotFoundError):
            db.repair_fasta(self.test_dir / Path("data/fasta_test_data"))

    def test_repair_wrong_type(self):
        db = blast.Database(source_dir=self.test_dir / Path("data/fasta_test_data"), name="test_db")
        with self.assertRaises(TypeError):
            db.repair_fasta("not a path object")

    def test_aggregate_ok(self):
        blast.Database(source_dir=Path(""), name="test_db")._aggregate(self.test_dir / Path("data/fasta_test_data"),
                                                                       Path("blast_input.fasta"))
        self.assertTrue(Path("blast_input.fasta").exists())
        self.assertTrue(os.path.getsize(Path("blast_input.fasta")) > 0)

    def test_aggregate_empty(self):
        with self.assertRaises(ValueError):
            blast.Database(source_dir=Path(""), name="test_db")._aggregate(
                self.test_dir / Path("data/fasta_test_empty_data"), Path("blast_input.fasta"))

    def test_aggregate_wrong_dir(self):
        with self.assertRaises(FileNotFoundError):
            blast.Database(source_dir=Path(""), name="test_db")._aggregate(self.test_dir / Path("not existing"),
                                                                           Path("blast_input.fasta"))
            blast.Database(source_dir=Path(""), name="test_db")._aggregate(self.test_dir / Path(
                "data/blast_query_data/seqs.fasta"), Path("blast_input.fasta"))

    def test_aggregate_wrong_type(self):
        with self.assertRaises(TypeError):
            db = blast.Database(source_dir=Path(""), name="test_db")
            db._aggregate("not a path obj", Path("blast_input.fasta"))
            db._aggregate(Path(""), "not a path obj")

    def test_create(self, ):
        db, output = blast.Database(source_dir=self.test_dir / Path("data/fasta_test_data"), name="test_db").create()
        self.assertTrue(Path("test_db.nhr").exists())
        self.assertTrue(Path("test_db.nhr").exists())
        self.assertTrue(Path("test_db.nsq").exists())
        db.clear_files()

    @patch('jasper.blast.Database._aggregate')
    def test_create_call(self, mock: Mock):
        blast.Database(source_dir=self.test_dir / Path("data/fasta_test_data"), name="test_db").create()
        mock.assert_called_once()

    def test_query_ok(self):
        db = blast.Database(source_dir=self.test_dir / Path("data/blast_query_data"), name="test_db")
        db.create()
        expected = pd.DataFrame(columns=['Virus', 'Host', 'Score'],
                                data=[['seqs|1', 'seqs|1', 10], ['seqs|2', 'seqs|2', 9]])
        res = db.query(self.test_dir / Path("data/blast_query_data"),
                       config={"task": "blastn-short"},
                       blast_format="10 qseqid sseqid score",
                       headers=('Virus', 'Host', 'Score'))
        self.assertIsInstance(res, pd.DataFrame)
        self.assertTrue(expected.equals(res))

    def test_query_wrong_dir(self):
        with self.assertRaises(FileNotFoundError):
            db = blast.Database(source_dir=self.test_dir / Path("data/fasta_test_data"), name="test_db")
            db.create()
            db.query(self.test_dir / Path("data/blast_query_data/seqs.fasta"),
                     config={},
                     blast_format="10 qseqid sseqid score",
                     headers=tuple())

    def test_query_wrong_command(self):
        with self.assertRaises(subprocess.SubprocessError):
            db = blast.Database(source_dir=self.test_dir / Path("data/fasta_test_data"), name="test_db")
            db.create()
            db.query(self.test_dir / Path("data/fasta_test_data/"),
                     blast_format="test_wrong_format",
                     config={"task": "blastn-short"},
                     headers=tuple())
        db.clear_files()

    def test_query_wrong_type(self):
        with self.assertRaises(TypeError):
            db = blast.Database(source_dir=self.test_dir / Path("data/fasta_test_data"), name="test_db")
            db.query("not a path obj",
                     blast_format="",
                     config={"task": "blastn-short"},
                     headers=tuple())
            db.query(self.test_dir / Path("data/fasta_test_data"),
                     config="not a dict", blast_format="",
                     headers=tuple())

    def test_query_wrong_kwargs(self):
        with self.assertRaises(ValueError):
            db = blast.Database(source_dir=self.test_dir / Path("data/fasta_test_data"), name="test_db")
            db.query(self.test_dir / Path("data/fasta_test_data"), blast_format="", headers=tuple(), config={
                "outfmt": "test",
                "query": "test",
            })

    def test_clear_files(self):
        db, tmp = blast.Database(source_dir=self.test_dir / Path("data/fasta_test_data"), name="test_db").create()
        db.clear_files()
        self.assertFalse((self.test_dir / Path("test_db.nhr")).exists())
        self.assertFalse((self.test_dir / Path("test_db.nin")).exists())
        self.assertFalse((self.test_dir / Path("test_db.nsq")).exists())

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
