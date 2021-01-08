import argparse
from pathlib import Path
from multiprocessing import Pool
from jasper import database
from jasper import crispr
import time

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
    config = {
        "source_path": "example_data/host",
        "db_name": "my_db",
        "repair_source_files": True,
        "repair_target_files": True,
    }
    # db = database.Database(config)
    # db.create()
    # query_df = db.query_multiple("example_data/virus")
    # print(query_df)
    crispr_finder = crispr.CrisprFinder(config)
    crispr_dict = crispr_finder.retrieve_spacers()
    assert not Path("spacers.fasta").exists()
    for key in crispr_dict.keys():
        for i, spacer in enumerate(crispr_dict[key], 1):
            with open("spacers.fasta", "a+") as fh:
                seq_name = f"{key}.spacer|{i}"
                SeqIO.write(SeqRecord(Seq(spacer),
                                      name=seq_name,
                                      id=seq_name,
                                      description=""),
                            fh, 'fasta')
    vir_db = database.Database({
        "source_path": "example_data/virus",
        "db_name": "vir_db",
        "repair_source_files": True,
        "repair_target_files": False,
    })
    vir_db.create()
    crispr_blast_results = vir_db.query(Path("spacers.fasta"))
    print(crispr_blast_results)
