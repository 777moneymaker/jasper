#!/usr/bin/env python3
import shutil
import unittest
from pathlib import Path

import pandas as pd

from jasper import mash

class MashTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(MashTests, self).__init__(*args, **kwargs)
        self.test_dir = Path(__file__).parent.absolute()

    def test_mash_init_ok(self):
        try:
            mash.Mash(self.test_dir / Path("data/fasta_test_data"), self.test_dir / Path("data/blast_query_data"))
        except:
            self.fail("Init failed with valid args")

    def test_mash_init_wrong_path(self):
        with self.assertRaises(FileNotFoundError):
            mash.Mash(self.test_dir / Path("not_existing"), self.test_dir / Path("data/fasta_test_data"))
            mash.Mash(self.test_dir / Path("data/fasta_test_data"), self.test_dir / Path("not_existing"))
            mash.Mash(self.test_dir / Path("not_existing"), self.test_dir / Path("not_existing"))

    def test_mash_init_not_dir(self):
        with self.assertRaises(ValueError):
            mash.Mash(self.test_dir / Path("data/test_df_merge1.csv"), self.test_dir / Path("data/test_df_merge1.csv"))
            mash.Mash(self.test_dir / Path("data/fasta_test_data"), self.test_dir / Path("data/test_df_merge1.csv"))
            mash.Mash(self.test_dir / Path("data/test_df_merge1.csv"), self.test_dir / Path("data/fasta_test_data"))

    def test_mash_init_wrong_type(self):
        with self.assertRaises(TypeError):
            mash.Mash("data/test_df_merge1.csv", "data/test_df_merge1.csv")
            mash.Mash(self.test_dir / Path("data/fasta_test_data"), "data/test_df_merge1.csv")
            mash.Mash("data/test_df_merge1.csv", self.test_dir / Path("data/fasta_test_data"))

    def test_mash_sketches(self):
        ws = mash.Mash(self.test_dir / Path("data/fasta_test_data"), self.test_dir / Path("data/blast_query_data"))
        host_s, phage_s = ws.sketch("test_host_s", "test_phage_s")
        self.assertTrue(host_s.exists())
        self.assertTrue(phage_s.exists())
        host_s.unlink()
        phage_s.unlink()

    def test_mash_predicts(self):
        ws = mash.Mash(self.test_dir / Path("data/test_trnascan_host"), self.test_dir / Path("data/test_trnascan_phage"))
        host_s, phage_s = ws.sketch("test_host_s", "test_phage_s")
        outfile = ws.run_mash(host_s, phage_s, Path("test_outfile.csv"))
        self.assertTrue(outfile.exists())
        df = pd.read_table(outfile, names=("Query", "Target", "Distance", "pval", "hashes"))
        df["Query"] = df["Query"].map(lambda x: x.split("/")[-1].split(".")[0])
        df["Target"] = df["Target"].map(lambda x: x.split("/")[-1].split(".")[0])
        df.drop(columns=["hashes"], inplace=True)
        outfile.unlink()
        self.assertTrue(df.equals(pd.DataFrame(
            {
                "Query": ["NC_000117|1"],
                "Target": ["NC_000902|1"],
                "Distance": [1],
                "pval": [1]
            }
        )))
