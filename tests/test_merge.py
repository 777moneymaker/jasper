#!/usr/bin/env python3
import os
import sys
import unittest
from pathlib import Path

import pandas as pd

from jasper import merge


class Expando(object):
    pass


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


class MergeTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(MergeTests, self).__init__(*args, **kwargs)
        self.test_dir = Path(__file__).parent.absolute()

    def test_merge_ok(self):
        args = Expando()
        args.files = [self.test_dir / Path('data/test_df_merge1.csv'), self.test_dir / Path('data/test_df_merge2.csv')]
        args.output = "test_jasper_results.csv"

        with HiddenPrints():
            merge.main(args)

        self.assertTrue(Path("test_jasper_results.csv").exists())
        res = pd.read_csv(args.output)
        self.assertTrue(len(res) == 9)

        if Path("test_jasper_results.csv").exists():
            Path("test_jasper_results.csv").unlink()


if __name__ == "__main__":
    unittest.main()
