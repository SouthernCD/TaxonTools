import os
import sys
import copy
from thefuzz import process
from toolbiox.api.common.genome.blast import outfmt6_read_big
from toolbiox.lib.common.evolution.tree_operate import Phylo, draw_ascii, remove_pass_node
from toolbiox.lib.common.fileIO import read_list_file
from toolbiox.lib.common.os import multiprocess_running
from toolbiox.lib.common.util import printer_list
from toolbiox.lib.xuyuxing.evolution.taxonomy import build_taxon_database, blast_taxon_assign, get_common_tree, read_tax_record_dict_db, store_tax_record_into_sqlite, get_common_tree_by_1kp


def ExcludeTaxon_main(args):
    tax_record_dict = build_taxon_database(args.tax_dir)
    root_taxon = tax_record_dict[args.root_taxon_id]
    exc_taxon = tax_record_dict[args.taxon_id]

    # get_lineage
    root_lineage = root_taxon.get_lineage(tax_record_dict)
    exc_lineage = exc_taxon.get_lineage(tax_record_dict)

    # diff lineage
    diff_lineage = copy.deepcopy(exc_lineage)
    for i in root_lineage[:-1]:
        diff_lineage.remove(i)

    output_list = []
    for i in range(len(diff_lineage) - 1):
        tax_id, rank_note = diff_lineage[i]
        next_tax_id, rank_note = diff_lineage[i + 1]

        tax_tmp = tax_record_dict[tax_id]
        output_list.extend(
            list(set([j for j in tax_tmp.get_sons(tax_record_dict) if j != next_tax_id])))

    print(printer_list(output_list, ','))


def RankStat_main(args):
    # build database
    tax_record_dict = build_taxon_database(args.tax_dir, True)

    # get list
    want_id_list = []
    for tax_id in tax_record_dict:
        tax_tmp = tax_record_dict[tax_id]
        if tax_tmp.rank == args.search_rank and tax_tmp.is_child_of(args.top_taxonomy_id, tax_record_dict):
            want_id_list.append(tax_id)

    # output
    if not args.output_file is None:
        sys.stdout = open(args.output_file, 'w')

    print("tax_id\tsci_name\trank\tlineage")
    for tax_id in want_id_list:
        tax_tmp = tax_record_dict[tax_id]
        lineage_str = ""
        lineage_list = tax_tmp.get_lineage(tax_record_dict)
        for j in lineage_list:
            tax_id_tmp, rank_tmp = j
            lineage_str = lineage_str + \
                tax_record_dict[tax_id_tmp].sci_name + " (" + rank_tmp + ");"
        lineage_str = lineage_str.rstrip(";")

        print("%s\t%s\t%s\t%s" %
              (tax_id, tax_tmp.sci_name, tax_tmp.rank, lineage_str))


def Name2Lineage_main(args):
    taxonomy_dir = args.tax_dir
    query_file = args.query_file
    output_file = args.output_file
    fuzzy_search = args.fuzzy_search

    # build database
    tax_record_dict = build_taxon_database(taxonomy_dir, True)

    # make hash table
    tax_record_sciname_hash = {}
    tax_record_othername_hash = {}
    for tax_id in tax_record_dict:
        tax_tmp = tax_record_dict[tax_id]

        if tax_tmp.sci_name not in tax_record_sciname_hash:
            tax_record_sciname_hash[tax_tmp.sci_name] = []
        tax_record_sciname_hash[tax_tmp.sci_name].append(tax_id)

        if hasattr(tax_tmp, 'other_name'):
            for name_tmp in tax_tmp.other_name:
                if name_tmp not in tax_record_othername_hash:
                    tax_record_othername_hash[name_tmp] = []
                tax_record_othername_hash[name_tmp].append(tax_id)

    # read query file
    query_tax_name = read_list_file(query_file)

    # search query in hash table
    subject_list = []
    sci_name_list = list(tax_record_sciname_hash.keys())
    other_name_list = list(tax_record_othername_hash.keys())
    for search_string in query_tax_name:
        other_name_flag = False
        fuzzy_flag = False
        homonym_flag = False
        # search sci_name
        if search_string in sci_name_list:
            search_out = tax_record_sciname_hash[search_string]
        # search other_name
        elif search_string in other_name_list:
            other_name_flag = True
            search_out = tax_record_othername_hash[search_string]
        # fuzzy search all name
        elif fuzzy_search:
            fuzzy_flag = True
            search_out = tax_record_sciname_hash[
                process.extractOne(search_string, sci_name_list + other_name_list)[0]]
        else:
            search_out = []
        # remove older taxid
        search_out = list(set([tax_record_dict[i].tax_id for i in search_out]))

        # check if homonym
        if len(search_out) > 1:
            homonym_flag = True
            print(search_out)
        # add failed info for None
        if len(search_out) == 0:
            subject_list.append((search_string,))
        # add to subject_list
        for tax_id in search_out:
            tax_tmp = tax_record_dict[tax_id]
            subject_list.append((search_string, tax_tmp.tax_id, tax_tmp.sci_name,
                                 tax_tmp.get_lineage(tax_record_dict), other_name_flag, fuzzy_flag, homonym_flag))

    if output_file is None:
        file = sys.stdout
    else:
        file = open(output_file, 'w')

    file.write(
        "query\ttax_id\tsci_name\tlineage\tother_name_flag\tfuzzy_flag\thomonym_flag\n")
    for i in subject_list:
        if len(i) == 1:
            file.write(i[0] + "\terror\n")
        else:
            lineage_str = ""
            lineage_list = i[3]
            for j in lineage_list:
                tax_id_tmp, rank_tmp = j
                lineage_str = lineage_str + \
                    tax_record_dict[tax_id_tmp].sci_name + \
                    " (" + rank_tmp + ");"
            lineage_str = lineage_str.rstrip(";")
            file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                       (i[0], i[1], i[2], lineage_str, i[4], i[5], i[6]))


def CommonTree_main(args):

    query_list = read_list_file(args.tax_id_list)

    if args.mode == 'NCBI':

        tax_record_dict = read_tax_record_dict_db(args.tax_db, query_list)
        tax_id_list = []
        for tax_id in query_list:
            for tax_id_tmp, rank_tmp in tax_record_dict[tax_id].get_lineage:
                tax_id_list.append(tax_id_tmp)
        tax_id_list = list(set(tax_id_list)) + query_list
        tax_record_dict = read_tax_record_dict_db(args.tax_db, tax_id_list)

        tree = get_common_tree(query_list, tax_record_dict)
        tree = remove_pass_node(tree)

    elif args.mode == '1kp':
        one_kp_tree_file = os.path.dirname(os.path.abspath(
            __file__)) + "/one_kp/estimated_species_tree.tree"
        acc_to_taxon_map = os.path.dirname(
            os.path.abspath(__file__)) + "/one_kp/1kp.acc.taxon.txt"
        tree = get_common_tree_by_1kp(
            query_list, one_kp_tree_file, acc_to_taxon_map, args.tax_db)

    if args.output_file is None:
        draw_ascii(tree)
    else:
        with open(args.output_file, 'w') as f:
            Phylo.write(tree, f, 'newick')


def ID2Lineage_main(args):

    taxonomy_dir = args.tax_dir
    query_file = args.query_tax_id_file
    output_file = args.output_file

    # output as a tsv file

    tax_record_dict = build_taxon_database(taxonomy_dir)

    query_tax_id = read_list_file(query_file)
    if query_file is not None:
        query_tax_id = read_list_file(query_file)
    else:
        query_tax_id = tax_record_dict.keys()

    if output_file is not None:
        with open(output_file, 'w') as f:
            for tax_id in query_tax_id:
                if tax_id in tax_record_dict:
                    taxon = tax_record_dict[tax_id]
                    lineage_str = ""
                    lineage_list = taxon.get_lineage(tax_record_dict)
                    for i in lineage_list:
                        tax_id_tmp, rank_tmp = i
                        lineage_str = lineage_str + \
                            tax_record_dict[tax_id_tmp].sci_name + \
                            " (" + rank_tmp + ");"
                    lineage_str = lineage_str.rstrip(";")
                    f.write(tax_id + "\t" + lineage_str + "\n")
                else:
                    lineage_str = "error"
                    f.write(tax_id + "\t" + lineage_str + "\n")
    else:
        for tax_id in query_tax_id:
            taxon = tax_record_dict[tax_id]
            lineage_str = ""
            lineage_list = taxon.get_lineage(tax_record_dict)
            for i in lineage_list:
                tax_id_tmp, rank_tmp = i
                lineage_str = lineage_str + \
                    tax_record_dict[tax_id_tmp].sci_name + \
                    " (" + rank_tmp + ");"
            lineage_str = lineage_str.rstrip(";")
            print(tax_id + "\t" + lineage_str)


def DiamondTaxonAssign_main(args):
    taxon_dict = build_taxon_database(args.tax_dir)

    args_list = []
    args_id_list = []
    for query_record in outfmt6_read_big(args.blast_results,
                                         fieldname=["qseqid", "sseqid", "staxids", "pident", "length",
                                                    "mismatch",
                                                    "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                                                    "bitscore"]):

        query_name = query_record.qDef

        taxon_score_list = []
        for hit_tmp in query_record.hit:
            for taxon in hit_tmp.Hit_taxon_id:
                if taxon == '':
                    continue
                if taxon is None:
                    continue
                score = hit_tmp.hsp[0].Hsp_bit_score
                taxon_score_list.append((taxon, score))

        taxon_id_list = list(set([i[0] for i in taxon_score_list]))

        if taxon_id_list == []:
            continue

        common_tree = get_common_tree(taxon_id_list, taxon_dict)

        # root_taxon, fisher_taxon, q_value = blast_taxon_assign(taxon_score_list, taxon_dict, True)

        args_list.append((taxon_score_list, common_tree))
        args_id_list.append(query_name)

    cmd_result = multiprocess_running(
        blast_taxon_assign, args_list, args.num_threads, args_id_list=args_id_list)

    with open(args.output_file, 'w') as f:
        for query_id in args_id_list:
            if cmd_result[query_id]['output'] == "bug":
                print(query_id)
            else:
                root_taxon, fisher_taxon, q_value = cmd_result[query_id]['output']
                f.write("%s\t%s\t%s\t%s\t%s\t%.5f\n" % (
                    query_id, root_taxon, taxon_dict[root_taxon].sci_name, fisher_taxon,
                    taxon_dict[fisher_taxon].sci_name, q_value))


def StoreTaxonDB_main(args):
    tax_record_dict = build_taxon_database(args.tax_dir)
    store_tax_record_into_sqlite(tax_record_dict, args.tax_db_file)
