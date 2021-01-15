import pandas as pd
from itertools import chain


def multiply_frames_by(frames: list, weigths: list):
    results = []
    for weight, dframe in zip(weigths, frames):
        dframe['Score'] = dframe['Score'] * weight
        results.append(dframe)
    return results


def main(args):
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
