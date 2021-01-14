import argparse
from pathlib import Path
from multiprocessing import Pool
from jasper import database
from jasper import crispr
from jasper import trna
from jasper import utils
import time
import shutil

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord




def parse_args():
    parser = argparse.ArgumentParser(description="JASPER is a program for bacterial hosts prediction",
                                     epilog="Made by Milosz Chodkowski 2020, PUT Poznan."
                                            "Check my github @ github.com/777moneymaker")
    # parser.add_argument("-vir", "--virus",
    #                     required=True,
    #                     type=str,
    #                     dest='virus_dir',
    #                     help='directory containing virus seq files.')
    # parser.add_argument("-hst", "--host",
    #                     required=True,
    #                     type=str,
    #                     dest='host_dir',
    #                     help='directory containing host seq files.')
    parser.add_argument("--allowed",
                        required=False,
                        type=int,
                        default=1,
                        dest='allowed_mis',
                        help='Number of allowed mismatches between vir genome and crispr spacer in blast query.')
    parser.add_argument("-mc", "--mega_config",
                        required=False,
                        type=str,
                        dest='mega_config',
                        help='File containing megablast config used in genome-genome query.')
    parser.add_argument("-sc", "--short_config",
                        required=False,
                        type=str,
                        dest='short_config',
                        help='File containing blastn-short config used in genome-crispr_spacer query.')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    args.short_config = utils.parse_config(args.short_config, {
        "task": "blastn-short",
        "evalue": 1,
        "gapopen": 10,
        "gapextend": 2,
        "penalty": -1,
        "word_size": 7,
        "dust": "no",
    })


    # crispr_finder = crispr.CrisprFinder(Path(args.host_dir), "NoName").retrieve_spacers()
    crispr_finder = crispr.CrisprFinder(Path("example_data/host/"), "NoName").retrieve_spacers()
    # vir_db = database.Database(Path(args.virus_dir), "vir_db").create()
    vir_db = database.Database(Path("example_data/virus"), "vir_db").create()

    short_results = vir_db.query(Path("crispr_spacers/"),
                                 config=args.short_config,
                                 blast_format="10 qseqid sseqid score qlen length mismatch gaps",
                                 headers=["Spacer", "Virus", "Score", "Qlen", "Alen", "Mis", "Gap"])
    vir_db.clear_files()
    shutil.rmtree(Path("crispr_spacers/"))

    short_results[["Score", "Qlen", "Alen", "Mis", "Gap"]].apply(pd.to_numeric)
    short_results['Allowed'] = short_results['Qlen'] - (
                short_results['Alen'] - short_results['Mis'] - short_results['Gap'])
    short_results['Allowed'] = short_results['Allowed'].apply(pd.to_numeric)
    short_results = short_results.drop(columns=['Qlen', 'Alen', 'Mis', 'Gap'])

    short_results = short_results[short_results['Allowed'] <= args.allowed_mis].drop(columns='Allowed').reset_index()
    short_results.reindex(["Virus", "Spacer", "Score"], axis=1)
    short_results["Spacer"] = short_results["Spacer"].map(lambda x: x.split("|")[0])

    short_results = short_results.groupby(["Virus", "Spacer"]).sum().reset_index()
    idx = short_results.groupby(['Virus'])['Score'].idxmax()
    crispr_results = short_results.loc[idx].sort_values('Score', ascending=False).reset_index(drop=True)

    print("blastn-short results (vir_genome-spacers query): ", crispr_results, sep='\n')

# if __name__ == "__main__":
#     trna_scanner = trna.tRNAScanner(Path("example_data/virus"), "no_name")
#     results = trna_scanner.scan()
#     print(results)