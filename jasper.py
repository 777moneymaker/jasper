import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(prog='jasper',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="")

    subparser = parser.add_subparsers(dest="subparser_name", required=True)

    # Blast parser
    parser_b = subparser.add_parser(name='blast',
                                    description="JASPER is a program for prediction of virus's host. This module is "
                                                "performing genome-genome blast query with given data and config.",
                                    epilog="Made by Milosz Chodkowski 2020, PUT Poznan. "
                                           "Check my github: github.com/777moneymaker",
                                    usage="jasper[.py] blast [-h] --virus VIRUS_DIR [--config BLASTN_CONFIG] --use_db "
                                          "USE_DB_NAME |  (--create_db CREATE_DB_NAME --host HOST_DIR)")
    parser_b._optionals.title = "arguments"
    parser_b.add_argument("--virus",
                          required=True,
                          type=str,
                          dest='virus_dir',
                          help='directory containing virus seq files.')
    parser_b.add_argument("--config",
                          required=False,
                          type=str,
                          dest='blastn_config',
                          help='File containing megablast config used in genome-genome query.')
    parser_b.add_argument("--clear",
                          action='store_true',
                          dest='clear_after',
                          help='Specifies if the database files should be deleted after analysis.')
    parser_b.add_argument('--output',
                          required=False,
                          type=str,
                          dest="output_file",
                          default="blast_results.csv",
                          help="Output file with final results.")


    group_b = parser_b.add_argument_group()
    group_b.add_argument("--use_db",
                         type=str,
                         dest='use_db_name',
                         help='Name for a existing database that will be used.')
    group_b.add_argument("--create_db",
                            type=str,
                            dest='create_db_name',
                            help='Name for a new database that will be created.')
    group_b.add_argument("--host",
                            type=str,
                            dest='host_dir',
                            help='directory containing host seq files.')

    # Crispr parser
    parser_c = subparser.add_parser(name='crispr',
                                    description="JASPER is a program for prediction of virus's host. This module is "
                                                "performing crispr spacers analysis which is "
                                                "vir_genome-spacer blast query with given data and config.",
                                    epilog="Made by Milosz Chodkowski 2020, PUT Poznan. "
                                           "Check my github: github.com/777moneymaker",
                                    usage="jasper crispr[.py] [-h] -host HOST_DIR [--config BLASTN_CONFIG] --use_db "
                                          "USE_DB_NAME |  (--create_db CREATE_DB_NAME --virus VIRUS_DIR)")
    parser_c._optionals.title = "arguments"
    parser_c.add_argument("--host",
                          type=str,
                          required=True,
                          dest='host_dir',
                          help='directory containing host seq files.')
    parser_c.add_argument("--config",
                          required=False,
                          type=str,
                          dest='short_config',
                          help='File containing megablast config used in genome-genome query.')
    parser_c.add_argument("--max_mismatches",
                          required=False,
                          type=int,
                          default=1,
                          dest='allowed_mis',
                          help='Number of allowed mismatches between vir genome and crispr spacer in blast query.')
    parser_c.add_argument("--clear",
                          action='store_true',
                          dest='clear_after',
                          help='Specifies if the database files should be deleted after analysis.')
    parser_c.add_argument('--output',
                          required=False,
                          type=str,
                          dest="output_file",
                          default="crispr_results.csv",
                          help="Output file with final results.")

    group_c = parser_c.add_argument_group()
    group_c.add_argument("--use_db",
                         type=str,
                         dest='use_db_name',
                         help='Name for a existing database that will be used.')
    group_c.add_argument("--create_db",
                            type=str,
                            dest='create_db_name',
                            help='Name for a new database that will be created.')
    group_c.add_argument("--virus",
                            type=str,
                            dest='virus_dir',
                            help='directory containing virus seq files.')

    args = parser.parse_args()
    if args.subparser_name == 'blast':
        if args.use_db_name and any([args.create_db_name, args.host_dir]):
            parser.error("Mutually exclusive argument groups were used. Check usage.")
        if not args.use_db_name and not all([args.create_db_name, args.host_dir]):
            parser.error("You must specify database name and a host files directory. Check usage.")

    if args.subparser_name == 'crispr':
        if args.use_db_name and any([args.create_db_name, args.virus_dir]):
            parser.error("Mutually exclusive argument groups were used. Check usage.")
        if not args.use_db_name and not all([args.create_db_name, args.virus_dir]):
            parser.error("You must specify database name and a virus files directory. Check usage.")
    return args


if __name__ == '__main__':
    args = parse_args()
    subparser_name = args.subparser_name
    if subparser_name == 'blast':
        from jasper import blast
        blast.main(args)
    elif subparser_name == 'crispr':
        from jasper import crispr
        crispr.main(args)
