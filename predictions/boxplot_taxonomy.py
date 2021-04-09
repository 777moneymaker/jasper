import numpy as np
import json
from functools import reduce
from pathlib import Path
import pandas as pd
#import seaborn as sns
import matplotlib.pyplot as plt
import itertools
from collections import defaultdict

def rank_frames(frames):
    for i, _ in enumerate(frames):
        col_name = frames[i].name + "Rank"
        frames[i][col_name] = frames[i].groupby(["Virus"], as_index=False)["Score"].rank(method='dense', ascending=False).astype(int)
    return frames


def main():
    frames = []
    # ['blast.csv', 'crispr.csv', 'mash.csv', 'wish.csv']
    for fl in ['blast.csv']:
        frame = pd.read_csv(fl)
        frame.name = Path(fl).stem
        frames.append(frame)

    ranked = rank_frames(frames)

    ranked_fixed = []
    for frame in ranked:
        f = frame.drop(columns=['Score'])
        f.name = frame.name
        ranked_fixed.append(f)

    # filtered = []
    # for frame in ranked_fixed:
    #     df = frame.groupby(['Virus'], as_index=False).apply(lambda grp: grp.loc[grp.iloc[:, -1] == grp.iloc[:, -1].min()]).reset_index(drop=True)
    #     df.name = frame.name
    #     filtered.append(df)

    with open('true_positives.json') as fh:
        d = json.load(fh)

    with open('edwards2016/virus/virus.json') as fh:
        d_vir = json.load(fh)

    with open('edwards2016/host/host.json') as fh:
        d_host = json.load(fh)


    taxonomy = [
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species"]

    # counts = defaultdict(list)
    # for frame in ranked_fixed:
    #     for vir in d.keys():     
    #         df = frame[frame["Virus"] == vir]
    #         for hst in df["Host"]:
    #             if hst in d[vir]:
    #                 counts[frame.name].append(list(df["Host"]).index(hst) + 1)

    counts = defaultdict(list)
    for frame in ranked_fixed:
        for i, vir in enumerate(d.keys(), 1):
            df = frame[frame["Virus"] == vir]

            vir_record = d_vir[vir]
            vir_host_tax = vir_record["host"]["lineage_ranks"][-2]

            #genuses = set([d_host[hst]["lineage_names"][-2] for hst in df["Host"]])
            for hst in df["Host"]:
                hst_record = d_host[hst]
                hst_tax = hst_record["lineage_names"][-2]
                
                if hst_tax == vir_host_tax:
                    counts[frame.name].append(list(df["Host"]).index(hst) + 1)
            print(frame.name, i)

    print(counts)
    labels, data = [*zip(*counts.items())]
    labels, data = counts.keys(), counts.values()
    plt.boxplot(data)
    plt.xticks(range(1, len(labels) + 1), labels)
    plt.savefig('boxplot_genus.png', dpi=200)

if __name__ == '__main__':
    main()