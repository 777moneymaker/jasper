import shutil
import subprocess
import unittest
from unittest.mock import patch
from jasper import crispr
from tests import utils
from pathlib import Path


class BlastTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(BlastTests, self).__init__(*args, **kwargs)
        self.test_dir = Path(__file__).parent.absolute()

    def test_finder_init_ok(self):
        try:
            finder = crispr.CrisprFinder(source_dir=self.test_dir / Path("data"), name="test_db")
        except:
            self.fail("Valid path caused exception.")

    def test_finder_init_wrong(self):
        with self.assertRaises(FileNotFoundError):
            crispr.CrisprFinder(source_dir=self.test_dir / Path("not_existing_path"), name="test_db")
            crispr.CrisprFinder(source_dir=self.test_dir / Path("data/blast_query_data/seqs.fasta"), name="test_db")

    def test_make_output_dir_ok(self):
        finder = crispr.CrisprFinder(source_dir=self.test_dir / Path("data"), name="test_db")
        finder._make_output_dir(Path("test_dir"))
        self.assertTrue(Path("test_dir").exists())
        self.assertTrue(Path("test_dir").is_dir())
        Path("test_dir").rmdir()

    def test_make_output_dir_wrong_type(self):
        with self.assertRaises(TypeError):
            finder = crispr.CrisprFinder(source_dir=self.test_dir / Path("data"), name="test_db")
            finder._make_output_dir("not a path obj")

    def test_retrieve_spacers_ok(self):
        finder = crispr.CrisprFinder(source_dir=self.test_dir / Path("data/fasta_test_data"), name="test_db")
        finder.retrieve_spacers()
        self.assertTrue(Path("crispr_spacers/").exists())
        self.assertTrue(Path("crispr_spacers/").is_dir())
        self.assertTrue(Path("crispr_spacers/test_genome_spacers|1.fasta").exists())
        self.assertTrue(Path("crispr_spacers/test_genome_spacers|1.fasta").is_file())
        with open(Path("crispr_spacers/test_genome_spacers|1.fasta"), 'r') as fh:
            self.assertEqual(utils.crispr_output, fh.read())
        if Path("crispr_spacers/").exists():
            shutil.rmtree(Path("crispr_spacers/"))

    def test_find_crispr_spacers_ok(self):
        finder = crispr.CrisprFinder(source_dir=self.test_dir / Path("data/fasta_test_data"), name="test_db")
        finder.find_crispr_spacers(self.test_dir / Path("data/fasta_test_data/test_genome_spacers.fasta"), Path("piler_output"))
        self.assertTrue(Path("piler_output").exists())

        with open(Path("piler_output"), 'r') as fh:
            self.assertTrue(">NC_006449.1" in fh.read())
        Path("piler_output").unlink()

    def test_find_crispr_spacers_wrong_type(self):
        with self.assertRaises(TypeError):
            finder = crispr.CrisprFinder(source_dir=self.test_dir / Path("data/fasta_test_data"), name="test_db")
            finder.find_crispr_spacers("not a path obj", self.test_dir / Path("piler_output"))
            finder.find_crispr_spacers(self.test_dir / Path("data/fasta_test_data/test_genome_spacers.fasta"), "not a path obj")

    def test_find_crispr_spacers_wrong_file(self):
        with self.assertRaises(FileNotFoundError):
            finder = crispr.CrisprFinder(source_dir=self.test_dir / Path("data/fasta_test_data"), name="test_db")
            finder.find_crispr_spacers(Path("not_existing_file"), self.test_dir / Path("piler_output"))


if __name__ == '__main__':
    unittest.main()
