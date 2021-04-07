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

    filtered = []
    for frame in ranked_fixed:
        df = frame.groupby(['Virus'], as_index=False).apply(lambda grp: grp.loc[grp.iloc[:, -1] == grp.iloc[:, -1].min()]).reset_index(drop=True)
        df.name = frame.name
        filtered.append(df)

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

    counts = defaultdict(lambda: defaultdict(list))
    for frame in filtered:
        for vir in d.keys():
            if vir not in set(frame["Virus"]):
                counts[frame.name]["species"].append(False)
                counts[frame.name]["genus"].append(False)
                counts[frame.name]["family"].append(False)
                counts[frame.name]["order"].append(False)
                counts[frame.name]["class"].append(False)
                counts[frame.name]["phylum"].append(False)
            else:
                df = frame[frame["Virus"] == vir]
                vir_record = d_vir[vir]

                for hst in df["Host"]:
                    host_record = d_host[hst]

                    species = host_record["lineage_names"][-1] == vir_record["host"]["lineage_names"][-1]
                    genus = host_record["lineage_names"][-2] == vir_record["host"]["lineage_names"][-2]
                    family = host_record["lineage_names"][-3] == vir_record["host"]["lineage_names"][-3]
                    order = host_record["lineage_names"][-4] == vir_record["host"]["lineage_names"][-4]
                    clas = host_record["lineage_names"][-5] == vir_record["host"]["lineage_names"][-5]
                    phylum = host_record["lineage_names"][-6] == vir_record["host"]["lineage_names"][-6]

                    counts[frame.name]["species"].append(species)
                    counts[frame.name]["genus"].append(genus)
                    counts[frame.name]["family"].append(family)
                    counts[frame.name]["order"].append(order)
                    counts[frame.name]["class"].append(clas)
                    counts[frame.name]["phylum"].append(phylum)

    counts = {frame: {tax: sum(counts[frame][tax]) / len(counts[frame][tax]) * 100 for tax in taxonomy} for frame in counts.keys()}
    
    # for k in counts.keys():
    #     for tax in counts[k].keys():
    #         print(k, tax, f"{sum(counts[k][tax]) / len(counts[k][tax]) * 100 :.2f}")


    # with open("taxonomy_predictions_2.json", "w+") as fh:
    #     json.dump(counts, fh)


    df = pd.DataFrame(counts)
    print(df)

    plt.figure()
    ax = df.T.plot.bar(rot=45)
    # x_offset = -0.03
    # y_offset = 0.02
    # for p in ax.patches:
    #     b = p.get_bbox()
    #     val = f"{b.y1 + b.y0 :.2f}"      
    #     ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset))

    plt.savefig('taxonomy_predictions.png', dpi=200)

if __name__ == '__main__':
    main()