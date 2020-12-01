import argparse
from pathlib import Path
from multiprocessing import Pool
from jasper import database
import time



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
    # args = parse_args()
    # database.concat_fasta_files_cmd("")

    db = database.Database("example_data/host", "my_db")
    print("Quering")
    start = time.time()
    with Pool() as p:
        blast_results = p.map(db.query, list(Path("example_data/virus/").iterdir()))
    end = time.time()
    print(*blast_results, sep='\n')
    print(end - start, len(blast_results))
