import argparse
from taxontools.pipelines import *


def args_init():

    # argument parse
    parser = argparse.ArgumentParser(
        prog='TaxonTools',
    )

    subparsers = parser.add_subparsers(
        title='subcommands', dest="subcommand_name")

    # argparse for StoreTaxonDB
    parser_a = subparsers.add_parser('StoreTaxonDB',
                                     help='store a taxon database into a sqlite3 db', description='')
    parser_a.add_argument(
        "taxdump_tar_gz_file", help="Path for taxdump.tar.gz file", type=str)
    parser_a.add_argument(
        "taxdump_db_file", help="Path for sqlite3 file", type=str)

    # argparse for ID2Lineage
    parser_a = subparsers.add_parser('ID2Lineage',
                                     help='get taxonomy lineage by taxon ID', description='')
    parser_a.add_argument(
        "taxdump_db_file", help="Path for taxdump_db_file from StoreTaxonDB", type=str)
    parser_a.add_argument(
        "-l", "--taxon_id_list", help="query taxon ID, you can give me a file path or taxon IDs separate by comma", default=None, type=str)
    parser_a.add_argument(
        "-o", "--output_file", help="Output file name (default as $STDOUT)", default=None, type=str)

    # argparse for Name2Lineage
    parser_a = subparsers.add_parser('Name2Lineage',
                                     help='get taxonomy lineage by taxon name', description='')
    parser_a.add_argument(
        "taxdump_db_file", help="Path for taxdump_db_file from StoreTaxonDB", type=str)
    parser_a.add_argument(
        "-l", "--taxon_name_list", help="query taxon Name, you can give me a file path or taxon names separate by comma", default=None, type=str)
    parser_a.add_argument(
        "-o", "--output_file", help="Output file name (default as $STDOUT)", default=None, type=str)

    # # argparse for RankStat
    # parser_a = subparsers.add_parser('RankStat',
    #                                  help='to stat a rank taxonomy under given taxonmoy, such as return all order under fungi')

    # parser_a.add_argument(
    #     "tax_dir", help="Path for uncompressed taxdump.tar.gz", type=str)
    # parser_a.add_argument("top_taxonomy_id", help="the top taxonomy want to search (taxon_id: 4751 for fungi)",
    #                       type=str)
    # parser_a.add_argument(
    #     "search_rank", help="rank want to search such as order", type=str)
    # parser_a.add_argument(
    #     "-o", "--output_file", help="Output file name (default as $STDOUT)", default=None, type=str)

    # # argparse for CommonTree
    # parser_a = subparsers.add_parser('CommonTree',
    #                                  help='get common tree by leaf taxon ID', description='')
    # parser_a.add_argument(
    #     "tax_db", help="Formatted taxonomy database", type=str)
    # parser_a.add_argument(
    #     "tax_id_list", help="The file containing the tax_id to query", type=str)
    # parser_a.add_argument(
    #     "-o", "--output_file", help="Output file name (default as $STDOUT)", default=None, type=str)
    # parser_a.add_argument("-m", "--mode", help='common tree mode',
    #                       default="NCBI", choices=['NCBI', '1kp'])

    # # argparse for DiamondTaxonAssign
    # parser_a = subparsers.add_parser('DiamondTaxonAssign',
    #                                  help='parse blast results to get protein seq taxon', description='')
    # parser_a.add_argument("blast_results",
    #                       help="blast results need diamond output: qseqid sseqid staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore",
    #                       type=str)
    # parser_a.add_argument(
    #     "tax_dir", help="Path for uncompressed taxdump.tar.gz", type=str)
    # parser_a.add_argument(
    #     "-n", "--num_threads", help="cpu number (default as 10)", default=10, type=int)
    # parser_a.add_argument("-o", "--output_file", help="Output file name (default as taxon.assign.txt)",
    #                       default="taxon.assign.txt", type=str)

    # # argparse for ExcludeTaxon
    # parser_a = subparsers.add_parser('ExcludeTaxon',
    #                                  help='Exclude a taxon id, to get all taxon_id not include given taxon id',
    #                                  description='')
    # parser_a.add_argument(
    #     "taxon_id", help="taxon id which want to be excluded", type=str)
    # parser_a.add_argument(
    #     "tax_dir", help="Path for uncompressed taxdump.tar.gz", type=str)
    # parser_a.add_argument("-r", "--root_taxon_id", help="root for excluded (default as root for all database)",
    #                       default='1', type=str)

    args = parser.parse_args()

    return args


def main():

    args = args_init()
    args_dict = vars(args)

    if args_dict["subcommand_name"] == "StoreTaxonDB":
        StoreTaxonDB_main(args)

    elif args_dict["subcommand_name"] == "ID2Lineage":
        ID2Lineage_main(args)

    elif args_dict["subcommand_name"] == "Name2Lineage":
        Name2Lineage_main(args)

    # elif args_dict["subcommand_name"] == "CommonTree":
    #     CommonTree_main(args)

    # elif args_dict["subcommand_name"] == "RankStat":
    #     RankStat_main(args)

    # elif args_dict["subcommand_name"] == "ExcludeTaxon":
    #     ExcludeTaxon_main(args)

    # elif args_dict["subcommand_name"] == "DiamondTaxonAssign":
    #     DiamondTaxonAssign_main(args)


if __name__ == '__main__':
    main()
