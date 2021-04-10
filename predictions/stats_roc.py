import numpy as np
import json
from functools import reduce
from pathlib import Path
import pandas as pd
#import seaborn as sns
import matplotlib.pyplot as plt
import itertools
from collections import defaultdict
from sklearn import metrics
import matplotlib.gridspec as gridspec

# def name_frames(frames):
#     for i, _ in enumerate(frames):
#         col_name = frames[i].name + "Score"
#         frames[i].rename(columns={"Score": col_name}, inplace=True)
#     return frames

def get_all_max(sdf):
    sdf_max = sdf.max().to_frame().T
    sdf_max.dropna(axis=1, how='all', inplace=True)
    result = pd.concat([pd.merge(sdf, sdf_max[[c]], how='inner') for c in sdf_max.columns.to_list()[2:]])
    result.loc["#"] = ["# Std", ""] + list(result[result.columns.to_list()[2:]].std(ddof=0).round(3))
    return result


def main():
    frames = []
    # ['blast.csv', 'crispr.csv', 'mash.csv', 'wish.csv']
    for fl in ['blast.csv', 'crispr.csv', 'mash.csv', 'wish.csv']:
        frame = pd.read_csv(fl)
        frame.name = Path(fl).stem
        frames.append(frame)


    filtered = []
    for frame in frames:
        df = frame.groupby('Virus', as_index=False).apply(get_all_max).reset_index(drop=True)
        df.name = frame.name
        filtered.append(df)

    with open('true_positives.json') as fh:
        d = json.load(fh)

    with open('edwards2016/virus/virus.json') as fh:
        d_vir = json.load(fh)

    with open('edwards2016/host/host.json') as fh:
        d_host = json.load(fh)

    positives = {(v, h) for v in d.keys() for h in d[v]}
    all_relations = {(v, h) for v in d_vir.keys() for h in d_host.keys()}

    gs = gridspec.GridSpec(2, 2)
    fig = plt.figure(figsize=(12, 18))

    for i, frame in enumerate(frames):
        decisions = []
        values = []
        
        records = {(x[0], x[1]): x[2] for x in list(frame.sort_values(["Score"], ascending=False).to_records(index=False))}

        for j, rel in enumerate(all_relations, 1):
            value = records.get(rel, 0)
            decision = 1 if rel in positives else 0

            decisions.append(decision)
            values.append(value)
            # print(frame.name, j)
        # print(frame.name, f"decs: {sum(decisions)}", len(decisions))

        if i < 2:
            fig.add_subplot(gs[0, i])
        else:
            fig.add_subplot(gs[1, i - 2])


        fallout, recall, thresholds = metrics.roc_curve(decisions, values, drop_intermediate=False)
        # print(fallout, recall, thresholds)
        # print([d for d in decisions if d > 0][:10], [v for v in values if v > 0][:10])
        plt.title(frame.name + " ROC")
        plt.xlabel('Fallout')
        plt.ylabel('Recall')
        plt.plot(fallout, recall, label=f'AUC = {metrics.roc_auc_score(decisions, values) :.2f}')
        plt.legend(loc='lower right')
    plt.savefig('ROC_curve.png', dpi=200)
    print("ROC Done.")

    # # PR CURVE

    gs = gridspec.GridSpec(2, 2)
    fig = plt.figure(figsize=(12, 18))

    for i, frame in enumerate(frames):
        decisions = []
        values = []
        
        records = {(x[0], x[1]): x[2] for x in list(frame.sort_values(["Score"], ascending=False).to_records(index=False))}
        # positive_records = {(x[0], x[1]) for x in records if x in positives}

        for j, rel in enumerate(all_relations, 1):
            value = records.get(rel, 0)
            decision = 1 if rel in positives else 0

            decisions.append(decision)
            values.append(value)

        if i < 2:
            fig.add_subplot(gs[0, i])
        else:
            fig.add_subplot(gs[1, i - 2])

        precision, recall, thresholds = metrics.precision_recall_curve(decisions, values)
        # print(precision, recall, thresholds)
        plt.title(frame.name + " Precision-Recall")
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.plot(recall, precision, label=f'AP = {metrics.average_precision_score(decisions, values) :.2f}')
        plt.legend(loc='lower right')
    
    plt.savefig('PR_curve.png', dpi=200)
    print("PR Done.")

    # F1 SCORE
    # print("F1's:")
    # for i, frame in enumerate(frames):
    #     ground_true = []
    #     predictions = []
        
    #     records = {(x[0], x[1]): x[2] for x in list(frame.sort_values(["Score"], ascending=False).to_records(index=False))}
    #     # positive_records = {(x[0], x[1]) for x in records if x in positives}

    #     for j, rel in enumerate(all_relations, 1):
    #         ground_true_val = 1 if rel in positives else 0
    #         prediction = 1 if rel in records else 0

    #         ground_true.append(ground_true_val)
    #         predictions.append(prediction)
    
    #     f1 = metrics.f1_score(ground_true, predictions)

    #     print(frame.name, f"{f1 :.3f}")

if __name__ == '__main__':
    main()