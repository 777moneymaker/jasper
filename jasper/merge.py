#!/usr/bin/env python3
"""This module manages the final merging of results from each single module."""

import numbers
import numpy as np
from itertools import chain
from functools import reduce
import pandas as pd


def main(args):
    """Main function for running module.
    Args:
        args (argparse obj): Arguments from right subparser.
    """

    frames = list(map(pd.read_csv, args.files))
    concat = reduce(lambda left, right: pd.merge(left, right, on=['Virus', 'Host'], how='outer'), frames)
    concat[concat.columns.to_list()[2:]].apply(pd.to_numeric)

    concat = concat[[c for c in concat if c.endswith('Rank') or c in ("Virus", "Host")]]
    concat = concat.reindex(["Virus", "Host", "BlastRank", "CrisprRank", "tRNARank", "WishRank", "MashRank"], axis=1)
    concat.dropna(axis=1, inplace=True, how='all')
    concat.fillna({c: np.inf for c in concat.columns.to_list()[2:]}, inplace=True)

    concat = concat.sort_values(by=concat.columns.to_list()[2:]).reset_index(drop=True)
    concat.to_csv(args.output, index=False)
    print(concat)
    print('Saved final results to', args.output)
