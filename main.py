import argparse
from jasper import io


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
    args = parse_args()
    print(*io.retrieve(args.virus_dir), sep='\n')
    print(*io.retrieve(args.host_dir), sep='\n')
