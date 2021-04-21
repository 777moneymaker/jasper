import pandas as pd
import json
from collections import defaultdict

def frames():
  res = []
  for fl in ['blast.csv', 'crispr.csv', 'wish.csv', 'mash.csv']:
    fh = open(fl)
    df = pd.read_csv(fl)
    fh.close()
    df.name = fl.rstrip(".csv")
    res.append(df)
  return res

def get_all_max(sdf):
  sdf_max = sdf.max().to_frame().T
  sdf_max.dropna(axis=1, how='all', inplace=True)
  result = pd.concat([pd.merge(sdf, sdf_max[[c]], how='inner') for c in sdf_max.columns.to_list()[2:]])
  result.loc["#"] = ["# Std", ""] + list(result[result.columns.to_list()[2:]].std(ddof=0).round(3))
  return result

if __name__ == '__main__':
  # df = pd.read_csv('blast.csv')
  # hosts = set(df.Host)
  # viruses = set(df.Virus)

  with open("true_positives.json") as fh:
    d = json.load(fh)

  all_pos = {v: set(h) for v, h in d.items()}

  results = defaultdict(lambda: defaultdict(int))
  for frame in frames():
    filtered = frame.groupby('Virus', as_index=False).apply(get_all_max).reset_index(drop=True)
    records = {v: set(filtered.loc[filtered.Virus == v].Host) for v in set(frame.Virus)}
    i = 1
    for v, hsts in all_pos.items():
        if v not in records.keys():
          results[frame.name][v] = 0
          continue
        
        if len(hsts.intersection(records[v])) > 0:
          results[frame.name][v] = 1
        else:
          results[frame.name][v] = 0

  final_df.to_csv("upset.csv")
  print(final_df)

