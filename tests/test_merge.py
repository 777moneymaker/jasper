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

    def test_multiply_frames_ok(self):
        frames = [pd.DataFrame({"T": [1, 2, 3], "Score": [1, 2, 3]}),
                  pd.DataFrame({"T": [1, 2, 3], "Score": [1, 2, 3]})]
        res = merge.multiply_frames_by(frames, [2, 2])
        self.assertTrue(res[0].equals(pd.DataFrame({"T": [1, 2, 3], "Score": [2, 4, 6]})))
        self.assertTrue(res[1].equals(pd.DataFrame({"T": [1, 2, 3], "Score": [2, 4, 6]})))

    def test_multiply_frames_wrong_types(self):
        with self.assertRaises(TypeError):
            merge.multiply_frames_by(1, 2)
            merge.multiply_frames_by(pd.DataFrame({"T": [1, 2, 3], "Score": [1, 2, 3]}), 2)
            merge.multiply_frames_by(4, pd.DataFrame({"T": [1, 2, 3], "Score": [1, 2, 3]}))
            merge.multiply_frames_by(
                [pd.DataFrame({"T": [1, 2, 3], "Score": [1, 2, 3]}),
                 pd.DataFrame({"T": [1, 2, 3], "Score": [1, 2, 3]})],
                ['nan', 'nan']
            )

    def test_merge_diff_args(self):
        args = Expando()
        args.files = ['test1', 'test2']
        args.weights = [.6, .2, .2]
        with self.assertRaises(ValueError):
            merge.main(args)

    def test_merge_wrong_sum(self):
        args = Expando()
        args.files = ['test1', 'test2']
        args.weights = [1, 2, 3]
        with self.assertRaises(ValueError):
            merge.main(args)

    def test_merge_ok(self):
        args = Expando()
        args.files = [self.test_dir / Path('data/test_df_merge1.csv'), self.test_dir / Path('data/test_df_merge2.csv')]
        args.weights = [.4, .6]
        args.output = "merged_jasper_results.csv"

        with HiddenPrints():
            merge.main(args)

        self.assertTrue(Path("merged_jasper_results.csv").exists())
        res = pd.read_csv(args.output)
        self.assertTrue(res.equals(pd.DataFrame({
            "Virus": ["v3", "v2", "v1"],
            "Host": ["h3", "h2", "h1"],
            "Score": [7.0, 6.4, 3.4]
        })))

        if Path("merged_jasper_results.csv").exists():
            Path("merged_jasper_results.csv").unlink()


if __name__ == "__main__":
    unittest.main()
