#!/usr/bin/env python3
"""This module manages the final merging of results from each single module."""

import pandas as pd
from itertools import chain
import numbers


def multiply_frames_by(frames: list, weights: list):
    """This function multiplies 'Score' column in each DataFrame.

    Each DataFrame's 'Score' column is multiplied by corresponding value in weights list.

    Args:
        frames (list): List of dataframes to be multiplied.
        weights (list): Corresponding list of weights to be used in multiplication.
    Raises:
        TypeError: When given obj is of wrong type.
    Returns:
        (list) List of multiplied DataFrames.
    """
    if not isinstance(frames, list):
        raise TypeError("Given obj is not a list.")
    if not isinstance(weights, list):
        raise TypeError("Given obj is not a list.")
    if not all(isinstance(x, pd.DataFrame) for x in frames):
        raise TypeError("Some of the DataFrames are not DataFrame objects.")
    if not all(isinstance(x, numbers.Number) for x in weights):
        raise TypeError("Some of the weights are not numbers.")

    results = []
    for weight, dframe in zip(weights, frames):
        dframe['Score'] = dframe['Score'] * weight
        results.append(dframe)
    return results


def main(args):
    """Main function for running module.
    Args:
        args (argparse obj): Arguments from right subparser.
    """
    if len(args.files) != len(args.weights):
        raise ValueError("Number of files and weights is different.")
    if sum(args.weights) != 1:
        raise ValueError("Weights don't sum up to 1.")
    frames = list(map(pd.read_csv, args.files))
    frames = multiply_frames_by(frames, args.weights)

    final = pd.DataFrame(columns=['Virus', 'Host', 'Score'])
    final['Virus'] = list(chain(*list(map(lambda df: df['Virus'].tolist(), frames))))
    final['Host'] = list(chain(*list(map(lambda df: df['Host'].tolist(), frames))))
    final['Score'] = list(chain(*list(map(lambda df: df['Score'].tolist(), frames))))
    final['Score'] = final['Score'].apply(pd.to_numeric)

    final = final.groupby(["Virus", "Host"]).sum().reset_index()
    idx = final.groupby(['Virus'])['Score'].idxmax()
    final = final.loc[idx].sort_values('Score', ascending=False).reset_index(drop=True)

    final['Score'] = final['Score'].round(decimals=2)
    final.to_csv(args.output, index=False)
