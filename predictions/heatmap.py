import numpy as np
import json
from functools import reduce
from pathlib import Path
import pandas as pd
# import seaborn as sns
import matplotlib.pyplot as plt
import itertools

def rank_frames(frames):
    for i, _ in enumerate(frames):
        col_name = frames[i].name + "Rank"
        frames[i][col_name] = frames[i].groupby(["Virus"], as_index=False)["Score"].rank(method='dense', ascending=False).astype(int)
    return frames


def main():
    frames = []
    for fl in ['blast.csv', 'crispr.csv', 'mash.csv', 'wish.csv']:
        frame = pd.read_csv(fl)
        frame.name = Path(fl).stem
        frames.append(frame)

    ranked = rank_frames(frames)

    ranked_fixed = []
    for frame in ranked:
        f = frame.drop(columns=['Score'])
        f.name = frame.name
        ranked_fixed.append(f)

    with open('./true_positives.json') as fh:
        d = json.load(fh)

    filtered = []
    for frame in ranked_fixed:
        df = frame.groupby(['Virus'], as_index=False).apply(lambda grp: grp.loc[grp.iloc[:, -1] == grp.iloc[:, -1].min()]).reset_index(drop=True)
        df.name = frame.name
        filtered.append(df)

    # counts = {}
    # for frame in filtered:
    #     counts[frame.name] = {}
    #     for vir in d.keys():     
    #         if not vir in frame["Virus"]:
    #             counts[frame.name][vir] = 0

    #         df = frame[frame["Virus"] == vir]
    #         if any(hst in d[vir] for hst in df["Host"]):
    #             counts[frame.name][vir] = 1
    #         else:
    #             counts[frame.name][vir] = 0

    upset = {}
    for frame in filtered:
        upset[frame.name] = {}
        for vir in d.keys():     
            if not vir in frame["Virus"]:
                upset[frame.name][vir] = 0

            df = frame[frame["Virus"] == vir]
            if any(hst in d[vir] for hst in df["Host"]):
                upset[frame.name][vir] = 1
            else:
                upset[frame.name][vir] = 0
    df = pd.DataFrame.from_dict(upset).T
    print(df)
    df.to_csv("upset.csv")
    print(sum(df["blast"]))
    exit(0)

    for k in counts:
        counts[k] = dict(sorted(counts[k].items(), key=lambda x: x[0]))
    df = pd.DataFrame.from_dict(counts).transpose()
    df.fillna(0, inplace=True)
    print(df)
    # plt.figure(figsize = (16, 7))
    # ax = sns.heatmap(df, cbar=True, xticklabels=False, cmap='inferno_r')
    # plt.savefig('heat_prediction_sorted_name_value.png', dpi=200)

    for tool in counts.keys():
        score = (list(counts[tool].values()).count(1) / 820) * 100
        print(tool, score)

    combs = list(itertools.combinations(counts.keys(), 2))
    group_counts = {}
    for comb in combs:
        r1 = comb[0]
        r2 = comb[1]
        name = r1 + " & " + r2

        valid = {v1[0]: 1 if v1[1] == v2[1] == 1 else 0 for v1, v2 in zip(sorted(counts[r1].items()), sorted(counts[r2].items()))}
        group_counts[name] = valid
    print()
    for k in group_counts:
        group_counts[k] = dict(sorted(set(group_counts[k].items()), key=lambda x: (-x[1], x[0])))
    df = pd.DataFrame.from_dict(group_counts).transpose()
    print(df)
    plt.figure(figsize = (16, 7))
    ax = sns.heatmap(df, cbar=True, xticklabels=False, cmap='inferno_r')
    plt.savefig('heat_prediction_combs.png', dpi=200)
    
    for tool in group_counts.keys():
        score = (list(group_counts[tool].values()).count(1) / 820) * 100
        print(tool, score)


if __name__ == '__main__':
    main()

