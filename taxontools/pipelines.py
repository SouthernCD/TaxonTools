# import os
import sys
# import copy
# from thefuzz import process
# from toolbiox.api.common.mapping.blast import outfmt6_read_big
# from toolbiox.lib.common.evolution.tree_operate import Phylo, draw_ascii, remove_pass_node, get_newick_string
from toolbiox.lib.common.fileIO import read_list_file
from toolbiox.lib.common.os import is_path
# from toolbiox.lib.common.os import multiprocess_running, is_path
# from toolbiox.lib.common.util import printer_list
from taxontools.src.taxon import build_taxon_database, read_taxon_record_dict_db, store_taxon_record_into_sqlite, read_taxon_name_record_dict_db
# from taxontools.src.taxon import build_taxon_database, blast_taxon_assign, get_common_tree, read_taxon_record_dict_db, store_taxon_record_into_sqlite, get_common_tree_by_1kp, read_taxon_name_record_dict_db


def StoreTaxonDB_main(args):
    taxon_record_dict = build_taxon_database(args.taxdump_tar_gz_file)
    store_taxon_record_into_sqlite(taxon_record_dict, args.taxdump_db_file)


def ID2Lineage_main(args):
    # parse query list
    try:
        if is_path(args.taxon_id_list):
            taxon_id_list = read_list_file(args.taxon_id_list)
        else:
            taxon_id_list = args.taxon_id_list.split(",")
    except:
        raise Exception("Not query list")

    # load database
    taxon_record_dict = read_taxon_record_dict_db(
        args.taxdump_db_file, taxon_id_list)
    ancient_tax_id_list = []
    for tax_id in taxon_record_dict:
        taxon = taxon_record_dict[tax_id]
        ancient_tax_id_list.extend([i[0] for i in taxon.lineage])
    ancient_tax_id_list = list(set(ancient_tax_id_list))
    taxon_record_dict = read_taxon_record_dict_db(
        args.taxdump_db_file, list(set(taxon_id_list + ancient_tax_id_list)))

    # output as a tsv file
    output_file = args.output_file
    f = sys.stdout if output_file is None else open(output_file, 'w')
    f.write("tax_id\tsci_name\trank\tlineage\n")

    for tax_id in taxon_id_list:
        if tax_id in taxon_record_dict:
            taxon = taxon_record_dict[tax_id]
            lineage_str = ""
            for i in taxon.lineage:
                tax_id_tmp, rank_tmp = i
                lineage_str = lineage_str + \
                    taxon_record_dict[tax_id_tmp].sci_name + \
                    " (" + rank_tmp + ");"
            lineage_str = lineage_str.rstrip(";")
            f.write(tax_id + "\t" + taxon.sci_name + "\t" +
                    taxon.rank + "\t" + lineage_str + "\n")
        else:
            lineage_str = "error"
            f.write(tax_id + "\t" + lineage_str + "\n")

    # close
    if output_file is not None:
        f.close()


def Name2Lineage_main(args):
    # parse query list
    try:
        if is_path(args.taxon_name_list):
            taxon_name_list = read_list_file(args.taxon_name_list)
        else:
            taxon_name_list = args.taxon_name_list.split(",")
    except:
        raise Exception("Not query list")

    # load database
    taxon_record_dict = read_taxon_name_record_dict_db(
        args.taxdump_db_file, taxon_name_list)

    all_tax_id_list = []
    for tax_name in taxon_record_dict:
        taxon = taxon_record_dict[tax_name]
        all_tax_id_list.extend([i[0] for i in taxon.lineage])
    all_tax_id_list = list(set(all_tax_id_list))
    all_taxon_record_dict = read_taxon_record_dict_db(
        args.taxdump_db_file, all_tax_id_list)

    # output as a tsv file
    output_file = args.output_file
    f = sys.stdout if output_file is None else open(output_file, 'w')
    f.write("tax_id\tsci_name\trank\tlineage\n")

    for tax_name in taxon_name_list:
        if tax_name in taxon_record_dict:
            taxon = taxon_record_dict[tax_name]
            lineage_str = ""
            for i in taxon.lineage:
                tax_id_tmp, rank_tmp = i
                lineage_str = lineage_str + \
                    all_taxon_record_dict[tax_id_tmp].sci_name + \
                    " (" + rank_tmp + ");"
            lineage_str = lineage_str.rstrip(";")
            f.write(taxon.tax_id + "\t" + taxon.sci_name + "\t" +
                    taxon.rank + "\t" + lineage_str + "\n")
        else:
            lineage_str = "error"
            f.write(taxon.tax_id + "\t" + lineage_str + "\n")

    # close
    if output_file is not None:
        f.close()


if __name__ == '__main__':
    class abc():
        pass

    args = abc()
    args.taxon_id_list = "3702,3701"
    args.output_file = None
    args.taxdump_db_file = "/lustre/home/xuyuxing/Database/NCBI/taxonomy/taxdump.db"

    ID2Lineage_main(args)

    class abc():
        pass

    args = abc()
    args.taxon_name_list = "Solanum lycopersicum,Cuscuta"
    args.output_file = None
    args.taxdump_db_file = "/lustre/home/xuyuxing/Database/NCBI/taxonomy/taxdump.db"

    Name2Lineage_main(args)


    # def ExcludeTaxon_main(args):
    #     tax_record_dict = build_taxon_database(args.tax_dir)
    #     root_taxon = tax_record_dict[args.root_taxon_id]
    #     exc_taxon = tax_record_dict[args.taxon_id]

    #     # get_lineage
    #     root_lineage = root_taxon.get_lineage(tax_record_dict)
    #     exc_lineage = exc_taxon.get_lineage(tax_record_dict)

    #     # diff lineage
    #     diff_lineage = copy.deepcopy(exc_lineage)
    #     for i in root_lineage[:-1]:
    #         diff_lineage.remove(i)

    #     output_list = []
    #     for i in range(len(diff_lineage) - 1):
    #         tax_id, rank_note = diff_lineage[i]
    #         next_tax_id, rank_note = diff_lineage[i + 1]

    #         tax_tmp = tax_record_dict[tax_id]
    #         output_list.extend(
    #             list(set([j for j in tax_tmp.get_sons(tax_record_dict) if j != next_tax_id])))

    #     print(printer_list(output_list, ','))

    # def RankStat_main(args):
    #     # build database
    #     tax_record_dict = build_taxon_database(args.tax_dir, True)

    #     # get list
    #     want_id_list = []
    #     for tax_id in tax_record_dict:
    #         tax_tmp = tax_record_dict[tax_id]
    #         if tax_tmp.rank == args.search_rank and tax_tmp.is_child_of(args.top_taxonomy_id, tax_record_dict):
    #             want_id_list.append(tax_id)

    #     # output
    #     if not args.output_file is None:
    #         sys.stdout = open(args.output_file, 'w')

    #     print("tax_id\tsci_name\trank\tlineage")
    #     for tax_id in want_id_list:
    #         tax_tmp = tax_record_dict[tax_id]
    #         lineage_str = ""
    #         lineage_list = tax_tmp.get_lineage(tax_record_dict)
    #         for j in lineage_list:
    #             tax_id_tmp, rank_tmp = j
    #             lineage_str = lineage_str + \
    #                 tax_record_dict[tax_id_tmp].sci_name + \
    #                 " (" + rank_tmp + ");"
    #         lineage_str = lineage_str.rstrip(";")

    #         print("%s\t%s\t%s\t%s" %
    #               (tax_id, tax_tmp.sci_name, tax_tmp.rank, lineage_str))

    # def CommonTree_main(args):

    #     query_list = read_list_file(args.tax_id_list)

    #     if args.mode == 'NCBI':

    #         tax_record_dict = read_taxon_record_dict_db(
    #             args.tax_db, query_list)
    #         tax_id_list = []
    #         for tax_id in query_list:
    #             for tax_id_tmp, rank_tmp in tax_record_dict[tax_id].get_lineage:
    #                 tax_id_list.append(tax_id_tmp)
    #         tax_id_list = list(set(tax_id_list)) + query_list
    #         tax_record_dict = read_taxon_record_dict_db(
    #             args.tax_db, tax_id_list)

    #         tree = get_common_tree(query_list, tax_record_dict)
    #         tree = remove_pass_node(tree)

    #     elif args.mode == '1kp':
    #         one_kp_tree_file = os.path.dirname(os.path.abspath(
    #             __file__)) + "/one_kp/estimated_species_tree.tree"
    #         acc_to_taxon_map = os.path.dirname(
    #             os.path.abspath(__file__)) + "/one_kp/1kp.acc.taxon.txt"
    #         tree = get_common_tree_by_1kp(
    #             query_list, one_kp_tree_file, acc_to_taxon_map, args.tax_db)

    #     simple_tree = copy.deepcopy(tree)
    #     for clade in simple_tree.find_clades():
    #         clade.confidence = None
    #         clade.branch_length = None

    #     simple_tree_string = get_newick_string(
    #         simple_tree).replace(":0.00000", "")

    #     if args.output_file is None:
    #         draw_ascii(simple_tree)
    #     else:
    #         with open(args.output_file, 'w') as f:
    #             f.write(simple_tree_string)

    # def DiamondTaxonAssign_main(args):
    #     taxon_dict = build_taxon_database(args.tax_dir)

    #     args_list = []
    #     args_id_list = []
    #     for query_record in outfmt6_read_big(args.blast_results,
    #                                          fieldname=["qseqid", "sseqid", "staxids", "pident", "length",
    #                                                     "mismatch",
    #                                                     "gapopen", "qstart", "qend", "sstart", "send", "evalue",
    #                                                     "bitscore"]):

    #         query_name = query_record.qDef

    #         taxon_score_list = []
    #         for hit_tmp in query_record.hit:
    #             for taxon in hit_tmp.Hit_taxon_id:
    #                 if taxon == '':
    #                     continue
    #                 if taxon is None:
    #                     continue
    #                 score = hit_tmp.hsp[0].Hsp_bit_score
    #                 taxon_score_list.append((taxon, score))

    #         taxon_id_list = list(set([i[0] for i in taxon_score_list]))

    #         if taxon_id_list == []:
    #             continue

    #         common_tree = get_common_tree(taxon_id_list, taxon_dict)

    #         # root_taxon, fisher_taxon, q_value = blast_taxon_assign(taxon_score_list, taxon_dict, True)

    #         args_list.append((taxon_score_list, common_tree))
    #         args_id_list.append(query_name)

    #     cmd_result = multiprocess_running(
    #         blast_taxon_assign, args_list, args.num_threads, args_id_list=args_id_list)

    #     with open(args.output_file, 'w') as f:
    #         for query_id in args_id_list:
    #             if cmd_result[query_id]['output'] == "bug":
    #                 print(query_id)
    #             else:
    #                 root_taxon, fisher_taxon, q_value = cmd_result[query_id]['output']
    #                 f.write("%s\t%s\t%s\t%s\t%s\t%.5f\n" % (
    #                     query_id, root_taxon, taxon_dict[root_taxon].sci_name, fisher_taxon,
    #                     taxon_dict[fisher_taxon].sci_name, q_value))
