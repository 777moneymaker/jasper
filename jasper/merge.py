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


def get_all_min(sdf):
    sdf_min = sdf.min().to_frame().T
    sdf_min.dropna(axis=1, how='all', inplace=True)
    result = pd.concat([pd.merge(sdf, sdf_min[[c]], how='inner') for c in sdf_min.columns.to_list()[2:]])
    result.loc["#"] = ["# Std", ""] + list(result[result.columns.to_list()[2:]].std(ddof=0).round(3))
    return result


def main(args):
    """Main function for running module.
    Args:
        args (argparse obj): Arguments from right subparser.
    """
    pd.options.display.float_format = '{:.3f}'.format

    frames = []
    for fl in args.files:
        frame = pd.read_csv(fl)
        frame.name = Path(fl).stem
        frames.append(frame)
    ranked = rank_frames(frames)

    concat = reduce(lambda left, right: pd.merge(left, right, on=['Virus', 'Host'], how='outer', copy=False), ranked)
    concat[concat.columns.to_list()[2:]].apply(pd.to_numeric)
    concat = concat[[c for c in concat if c.endswith('Rank') or c in ("Virus", "Host")]]
    concat.fillna({c: np.nan for c in concat.columns.to_list()[2:]}, inplace=True)

    filtered = concat.groupby('Virus', as_index=False).apply(get_all_min).reset_index(drop=True)
    filtered.to_csv(args.output, index=False, na_rep="NaN")
    print(filtered)
    print('Saved final results to', args.output)

