import argparse
from pathlib import Path

from jasper import  database


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
    # database.make_local_db("example_data/host/")
    results = []
    for i, vir_fl in enumerate(Path('example_data/virus/').iterdir()):
        res = database.query(str(vir_fl.absolute()))
        results.append(res)
    print(*results, sep='\n')
