import shutil
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
            crispr.CrisprFinder(source_dir=self.test_dir / Path("data/fasta_test_data/seqs.fasta"), name="test_db")

    def test_make_output_dir(self):
        finder = crispr.CrisprFinder(source_dir=self.test_dir / Path("data"), name="test_db")
        finder._make_output_dir(Path("test_dir"))
        self.assertTrue(Path("test_dir").exists())
        self.assertTrue(Path("test_dir").is_dir())
        Path("test_dir").rmdir()

    def test_retrieve_spacers_ok(self):
        finder = crispr.CrisprFinder(source_dir=self.test_dir / Path("data/fasta_test_data"), name="test_db")
        finder.retrieve_spacers()
        self.assertTrue(Path("crispr_spacers/").exists())
        self.assertTrue(Path("crispr_spacers/").is_dir())
        self.assertTrue(Path("crispr_spacers/test_genome_spacers|1.fasta").exists())
        self.assertTrue(Path("crispr_spacers/test_genome_spacers|1.fasta").is_file())
        with open(Path("crispr_spacers/test_genome_spacers|1.fasta"), 'r') as fh:
            self.assertEqual(fh.read(), utils.crispr_output)
        shutil.rmtree(self.test_dir / Path("crispr_spacers/"))
