import argparse
from pathlib import Path
from multiprocessing import Pool
from jasper import database
from jasper import crispr
import time
import shutil

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

LOGO = r"""
      _   _    ____  ____  _____ ____  
     | | / \  / ___||  _ \| ____|  _ \ 
  _  | |/ _ \ \___ \| |_) |  _| | |_) |
 | |_| / ___ \ ___) |  __/| |___|  _ < 
  \___/_/   \_\____/|_|   |_____|_| \_\    
  
  ========================================
  
            |-----------------|
             \---------------/
              ~-_---------_-~
                ~-_---_-~
                   ~-_
                _-~---~-_
              _-~---------~-_
             /---------------\
            |-----------------|                               
"""


def parse_args():
    parser = argparse.ArgumentParser(description="JASPER is a program for bacterial hosts prediction",
                                     epilog="Made by Milosz Chodkowski 2020, PUT Poznan."
                                            "Check my github @ github.com/777moneymaker")
    parser.add_argument("-vir", "--virus",
                        required=True,
                        dest='virus_dir',
                        help='directory containing virus seq files')
    parser.add_argument("-hst", "--host",
                        required=True,
                        dest='host_dir',
                        help='directory containing host seq files')
    return parser.parse_args()


if __name__ == "__main__":
    print(LOGO)
    # args = parse_args()
    megablast_config = {
        "task": "megablast",
        "evalue": 10,
        "gapopen": 0,
        "gapextend": 2,
        "penalty": -2,
        "word_size": 28,
        "dust": "20 64 1",
    }
    db = database.Database({
        "source_path": "example_data/host",
        "db_name": "my_db",
    })
    db.create()
    query_df = db.query_multiple("example_data/virus", config=megablast_config)
    db.clear_files()

    blast_results = query_df.sort_values(['Virus', 'Score'], ascending=False).groupby(['Virus']).first().reset_index()
    print("Megablast results (genome-genome query): ")
    print(blast_results.sort_values(by="Score", ascending=False))

    crispr_finder = crispr.CrisprFinder({
        "source_path": "example_data/host",
        "db_name": "my_db",
    })
    crispr_finder.retrieve_spacers()

    vir_db = database.Database({
        "source_path": "example_data/virus",
        "db_name": "vir_db",
    })
    vir_db.create()

    spacers_blast_config = {
        "task": "blastn-short",
        "evalue": 1,
        "gapopen": 10,
        "gapextend": 2,
        "penalty": -1,
        "word_size": 7,
        "dust": "no",
    }

    crispr_blast_results = vir_db.query_multiple(Path("crispr_spacers/"), config=spacers_blast_config, headers=["Spacer", "Virus", "Score"])
    vir_db.clear_files()
    shutil.rmtree(Path("crispr_spacers/"))

    crispr_blast_results.reindex(["Virus", "Spacer", "Score"], axis=1)
    crispr_blast_results = crispr_blast_results.sort_values(['Virus', 'Score'], ascending=False).groupby(['Virus']).first().reset_index()

    print("blastn-short results (vir_genome-spacers query): ")
    print(crispr_blast_results)
