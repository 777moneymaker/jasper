#!/usr/bin/env python3
import shutil
import unittest
from pathlib import Path

import pandas as pd
from jasper import afree


class AfreeTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(AfreeTests, self).__init__(*args, **kwargs)
        self.test_dir = Path(__file__).parent.absolute()

    def test_wish_init_ok(self):
        try:
            afree.Wish(self.test_dir / Path("data/fasta_test_data"), self.test_dir / Path("data/blast_query_data"))
        except:
            self.fail("Init failed with valid args")

    def test_wish_init_wrong_path(self):
        with self.assertRaises(FileNotFoundError):
            afree.Wish(self.test_dir / Path("not_existing"), self.test_dir / Path("data/fasta_test_data"))
            afree.Wish(self.test_dir / Path("data/fasta_test_data"), self.test_dir / Path("not_existing"))
            afree.Wish(self.test_dir / Path("not_existing"), self.test_dir / Path("not_existing"))

    def test_wish_init_not_dir(self):
        with self.assertRaises(ValueError):
            afree.Wish(self.test_dir / Path("data/test_df_merge1.csv"), self.test_dir / Path("data/test_df_merge1.csv"))
            afree.Wish(self.test_dir / Path("data/fasta_test_data"), self.test_dir / Path("data/test_df_merge1.csv"))
            afree.Wish(self.test_dir / Path("data/test_df_merge1.csv"), self.test_dir / Path("data/fasta_test_data"))

    def test_wish_init_wrong_type(self):
        with self.assertRaises(TypeError):
            afree.Wish("data/test_df_merge1.csv", "data/test_df_merge1.csv")
            afree.Wish(self.test_dir / Path("data/fasta_test_data"), "data/test_df_merge1.csv")
            afree.Wish("data/test_df_merge1.csv", self.test_dir / Path("data/fasta_test_data"))

    def test_builds(self):
        wish = afree.Wish(self.test_dir / Path("data/fasta_test_data"), self.test_dir / Path("data/blast_query_data"))
        wish.build(Path("utest_model_dir"))

        self.assertTrue((Path("utest_model_dir").exists()))
        self.assertTrue((Path("utest_model_dir/test_genome_spacers|1.mm")).exists())

        shutil.rmtree(Path("utest_model_dir"))

    def test_predicts(self):
        wish = afree.Wish(self.test_dir / Path("data/fasta_test_data"), self.test_dir / Path("data/blast_query_data"))
        wish.build(Path("utest_model_dir"))
        res = wish.predict(Path("utest_model_dir"), Path("utest_output_dir"))
        self.assertTrue((Path("utest_output_dir") / Path("prediction.list")).exists())
        self.assertTrue((Path("utest_output_dir") / Path("llikelihood.matrix")).exists())

        shutil.rmtree(Path("utest_output_dir"))
        shutil.rmtree(Path("utest_model_dir"))

        expected = pd.DataFrame({
            "seqs|1": [-1.55978],
            "seqs|2": [-1.48808],
        })
        expected.index = pd.Index(["test_genome_spacers|1"])

        self.assertTrue(res.equals(expected))


if __name__ == '__main__':
    unittest.main()
