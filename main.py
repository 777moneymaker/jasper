import argparse
from pathlib import Path
from multiprocessing import Pool
from jasper import database
from jasper import crispr
import time

LOGO = r"""
      _   _    ____  ____  _____ ____  
     | | / \  / ___||  _ \| ____|  _ \ 
  _  | |/ _ \ \___ \| |_) |  _| | |_) |
 | |_| / ___ \ ___) |  __/| |___|  _ < 
  \___/_/   \_\____/|_|   |_____|_| \_\                                   
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
    # db = database.Database("example_data/host", "my_db", repair_host_files=True)
    # db.create()
    # query_df = db.query_multiple("example_data/virus")
    # print(query_df)
    crispr_db = crispr.Crispr("example_data/host", db_name="crispr_db", repair_host_files=True)
    crispr_db.create()
