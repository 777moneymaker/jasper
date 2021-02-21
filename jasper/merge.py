#!/usr/bin/env python3
"""This module manages the final merging of results from each single module."""

import numbers
import numpy as np
from itertools import chain
from functools import reduce
from pathlib import Path
import pandas as pd


def rank_frames(frames):
    for i, _ in enumerate(frames):
        col_name = frames[i].name + "Rank"
        frames[i][col_name] = frames[i].groupby(["Virus"], as_index=False)["Score"].rank(method='dense', ascending=False).astype(int)
    return frames


def main(args):
    """Main function for running module.
    Args:
        args (argparse obj): Arguments from right subparser.
    """

    frames = []
    for fl in args.files:
        frame = pd.read_csv(fl)
        frame.name = Path(fl).stem
        frames.append(frame)

    ranked = rank_frames(frames)
    # ranked.to_csv()
    # print(*ranked, sep='\n')

    frames = [frame.groupby(['Virus']).apply(lambda x: x.nlargest(1, ['Score'], keep='all')).reset_index(drop=True) for frame in frames]
    concat = reduce(lambda left, right: pd.merge(left, right, on=['Virus', 'Host'], how='outer'), frames)
    concat[concat.columns.to_list()[2:]].apply(pd.to_numeric)

    concat = concat[[c for c in concat if c.endswith('Rank') or c in ("Virus", "Host")]]
    # concat = concat.reindex(["Virus", "Host", "BlastRank", "CrisprRank", "tRNARank", "WishRank", "MashRank"], axis=1)
    concat.dropna(axis=1, inplace=True, how='all')
    concat.fillna({c: np.inf for c in concat.columns.to_list()[2:]}, inplace=True)

    ### MAGIC HERE ###
    filtered = concat.groupby(['Virus'], as_index=False).apply(lambda x: x.loc[x[x.columns.to_list()[2:]].idxmin()])
    # filtered = filtered.reset_index(drop=True)
    print(filtered)
    ##################

    # Sorting
    concat = concat.sort_values(by=concat.columns.to_list()[2:]).reset_index(drop=True)
    concat.to_csv(args.output, index=False)
    # print(concat)
    print('Saved final results to', args.output)
