#!/usr/bin/env python3
import shutil
import unittest
from pathlib import Path

from jasper import wish


class WishTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(WishTests, self).__init__(*args, **kwargs)
        self.test_dir = Path(__file__).parent.absolute()

    def test_wish_init_ok(self):
        try:
            wish.Wish(self.test_dir / Path("data/fasta_test_data"), self.test_dir / Path("data/blast_query_data"))
        except:
            self.fail("Init failed with valid args")

    def test_wish_init_wrong_path(self):
        with self.assertRaises(FileNotFoundError):
            wish.Wish(self.test_dir / Path("not_existing"), self.test_dir / Path("data/fasta_test_data"))
            wish.Wish(self.test_dir / Path("data/fasta_test_data"), self.test_dir / Path("not_existing"))
            wish.Wish(self.test_dir / Path("not_existing"), self.test_dir / Path("not_existing"))

    def test_wish_init_not_dir(self):
        with self.assertRaises(ValueError):
            wish.Wish(self.test_dir / Path("data/test_df_merge1.csv"), self.test_dir / Path("data/test_df_merge1.csv"))
            wish.Wish(self.test_dir / Path("data/fasta_test_data"), self.test_dir / Path("data/test_df_merge1.csv"))
            wish.Wish(self.test_dir / Path("data/test_df_merge1.csv"), self.test_dir / Path("data/fasta_test_data"))

    def test_wish_init_wrong_type(self):
        with self.assertRaises(TypeError):
            wish.Wish("data/test_df_merge1.csv", "data/test_df_merge1.csv")
            wish.Wish(self.test_dir / Path("data/fasta_test_data"), "data/test_df_merge1.csv")
            wish.Wish("data/test_df_merge1.csv", self.test_dir / Path("data/fasta_test_data"))

    def test_builds(self):
        ws = wish.Wish(self.test_dir / Path("data/fasta_test_data"), self.test_dir / Path("data/blast_query_data"))
        ws.build(Path("utest_model_dir"))

        self.assertTrue((Path("utest_model_dir").exists()))
        self.assertTrue((Path("utest_model_dir/test_genome_spacers|1.mm")).exists())

        shutil.rmtree(Path("utest_model_dir"))

    def test_predicts(self):
        ws = wish.Wish(self.test_dir / Path("data/fasta_test_data"), self.test_dir / Path("data/blast_query_data"))
        ws.build(Path("utest_model_dir"))
        res = ws.predict(Path("utest_model_dir"), Path("utest_output_dir"))

        self.assertTrue((Path("utest_output_dir") / Path("prediction.list")).exists())
        self.assertTrue((Path("utest_output_dir") / Path("llikelihood.matrix")).exists())

        shutil.rmtree(Path("utest_output_dir"))
        shutil.rmtree(Path("utest_model_dir"))


if __name__ == '__main__':
    unittest.main()
