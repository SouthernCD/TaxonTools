from itertools import combinations
from mlxtend.evaluate import permutation_test
from scipy import stats
from toolbiox.lib.common.evolution.tree_operate import add_clade_name, lookup_by_names, get_MRCA, remove_pass_node, get_offspring, draw_ascii, get_parent, get_ancestors, get_sons
from toolbiox.lib.common.math.set import merge_subset
from toolbiox.lib.common.util import printer_list
from toolbiox.lib.xuyuxing.math.stats import get_qvalue
import numpy as np


def blast_taxon_assign(taxon_score_list, common_tree, silence=True):
    """
    Give me blast score results of a seq, I will tell you the taxon of this sequence
    :param taxon_blast_score_dict: [("3702",100.00)]
    :param taxon_dict: build_taxon_database(taxon_dir)
    :return:
    """
    if not silence:
        print("# begin analysis:")

    # get the highest score for each taxon
    taxon_hash = {}
    for taxon_id, score in taxon_score_list:
        if taxon_id not in taxon_hash or taxon_hash[taxon_id] <= score:
            taxon_hash[taxon_id] = score

    if not silence:
        print("## get %d taxon\n" % len(taxon_hash))

    # build common_tree
    if not silence:
        print("# build common_tree:")

    taxon_id_list = list(taxon_hash.keys())
    # common_tree = get_common_tree(taxon_id_list, taxon_dict)
    common_tree = remove_pass_node(common_tree)
    common_tree = add_clade_name(common_tree)

    MRCA = common_tree.root.name
    if not silence:
        print("## get common_tree, MRCA: %s\n" % taxon_dict[MRCA].sci_name)

    leaf = common_tree.get_terminals()

    tree_dict = lookup_by_names(common_tree)

    # based on rank and fisher exact test
    if not silence:
        print("# based on rank and fisher exact test:")

    taxon_score_rank = sorted(
        taxon_hash, key=lambda x: taxon_hash[x], reverse=True)
    taxon_score_rank = [i for i in taxon_score_rank if i in tree_dict]

    fisher_p_dict = {}

    for clade in common_tree.find_clades(order='preorder'):
        all_sub_clade_id = [i.name for i in get_offspring(clade) + [clade]]
        all_sub_clade_have_hit = list(
            set(all_sub_clade_id) & set(taxon_id_list))

        top_hit_num = 0
        for i in all_sub_clade_have_hit:
            if taxon_score_rank.index(i) + 1 <= len(all_sub_clade_have_hit):
                top_hit_num += 1

        a = top_hit_num
        b = len(all_sub_clade_have_hit) - a
        c = b
        d = len(taxon_score_rank) - b

        oddsratio, pvalue = stats.fisher_exact([[a, b], [c, d]])
        clade.confidence = pvalue
        # print(pvalue)

        fisher_p_dict[clade.name] = pvalue

    fisher_qvalue_dict = get_qvalue(fisher_p_dict)

    good_clade_name = [
        i for i in sorted(
            fisher_qvalue_dict,
            key=lambda x: fisher_qvalue_dict[x]) if fisher_qvalue_dict[i] <= 0.05]

    if not silence:
        print("# There are %d taxon have q_value <= 0.05" %
              len(good_clade_name))

    # based on score value to say how
    # compared between brother score
    if not silence:
        permutation_p_value = {}
        for clade_name in good_clade_name:
            clade_tmp = tree_dict[clade_name]
            parent_tmp = get_parent(clade_tmp, common_tree)

            if len(parent_tmp) != 0:
                clade_leaf = clade_tmp.get_terminals()
                brother_leaf = [
                    i for i in parent_tmp.get_terminals() if i not in clade_leaf]

                clade_score = [taxon_hash[i.name]
                               for i in clade_leaf if i.name in taxon_hash]
                brother_score = [taxon_hash[i.name]
                                 for i in brother_leaf if i.name in taxon_hash]

                p_value = permutation_test(clade_score, brother_score,
                                           method='approximate',
                                           func=lambda x, y: np.mean(
                                               x) - np.mean(y),
                                           seed=0)

                permutation_p_value[clade_name] = p_value

            else:
                permutation_p_value[clade_name] = 0

        print("Taxon_SciName\tTaxon_id\tFisher\tFDR\tPermutation")
        for i in good_clade_name:
            print(
                "%s\t%s\t%.5e\t%.5e\t%.5f" %
                (taxon_dict[i].sci_name,
                 i,
                 fisher_p_dict[i],
                 fisher_qvalue_dict[i],
                 permutation_p_value[i]))
        # print(taxon_dict[i].sci_name, fisher_p_dict[i], fisher_qvalue_dict[i], permutation_p_value[i])

    # get full_supported_list
    if not silence:
        print("\n# based on node to find fully supported lineage:\n")

    full_supported_list = []
    for clade_name in good_clade_name:
        lineage_list = [i.name for i in get_ancestors(
            clade_name, common_tree)] + [clade_name]
        if set(lineage_list) & set(good_clade_name) == set(lineage_list):
            full_supported_list.append(clade_name)

    # full_supported_list = good_clade_name

    # get best lineage
    good_clade_lineage_list = []
    for clade_name in full_supported_list:
        lineage_list = [i.name for i in get_ancestors(
            clade_name, common_tree)] + [clade_name]
        good_clade_lineage_list.append(lineage_list)

    # if have two good lineage, use vote way
    merge_set_list = merge_subset(good_clade_lineage_list)

    vote_dict = {tuple(i): 0 for i in good_clade_lineage_list if set(
        i) in merge_set_list}

    for i in good_clade_lineage_list:
        for j in vote_dict:
            if set(i) & set(j) == set(i):
                vote_dict[j] += 1

    if not silence:
        num = 0
        for i in sorted(vote_dict, key=lambda x: vote_dict[x], reverse=True):
            print("#( %d ):" % num)
            lineage_txt = printer_list(
                [taxon_dict[j].sci_name for j in i], head="\t", sep=";")
            print(lineage_txt)
            print("Vote num: %d\n" % (vote_dict[i]))
            num = num + 1

    good_lineage_list = list(vote_dict.keys())

    # use score compare
    # find conflict point
    # get full_supported_list
    if not silence:
        print("\n# get final results:")

    while len(good_lineage_list) > 1:
        for comp_pair in combinations(good_lineage_list, 2):
            if not silence:
                print("\tcompare with %s vs %s:" %
                      (comp_pair[0][-1], comp_pair[1][-1]))

            A_clade = tree_dict[comp_pair[0][-1]]
            B_clade = tree_dict[comp_pair[1][-1]]
            mrca_clade = get_MRCA(A_clade, B_clade, common_tree)
            A_conflict_clade = [i for i in get_sons(
                mrca_clade) if i.name in comp_pair[0]][0]
            B_conflict_clade = [i for i in get_sons(
                mrca_clade) if i.name in comp_pair[1]][0]

            A_conflict_score = [taxon_hash[i.name]
                                for i in A_conflict_clade.get_terminals() if i.name in taxon_hash]
            B_conflict_score = [taxon_hash[i.name]
                                for i in B_conflict_clade.get_terminals() if i.name in taxon_hash]

            p_value = permutation_test(A_conflict_score, B_conflict_score,
                                       method='approximate',
                                       func=lambda x, y: np.mean(
                                           x) - np.mean(y),
                                       seed=0)
            if not silence:
                print("\tp_value: %f" % p_value)

            # A > B
            if p_value <= 0.05:
                good_lineage_list.remove(comp_pair[1])
                if not silence:
                    print("\tremove %s" % comp_pair[1][-1])
            # A < B
            elif p_value >= 0.95:
                good_lineage_list.remove(comp_pair[0])
                if not silence:
                    print("\tremove %s" % comp_pair[0][-1])
            else:
                good_lineage_list.remove(comp_pair[0])
                good_lineage_list.remove(comp_pair[1])

                mrca_lineage = tuple([i.name for i in get_ancestors(
                    mrca_clade.name, common_tree)] + [mrca_clade.name])
                good_lineage_list.append(mrca_lineage)
                if not silence:
                    print("\tremove %s and %s, change to %s\n" %
                          (comp_pair[0][-1], comp_pair[1][-1], mrca_clade.name))

            merge_set_list = merge_subset(good_lineage_list)

            good_lineage_list = [
                tuple(i) for i in good_clade_lineage_list if set(i) in merge_set_list]

            break

    if not silence:
        lineage_txt = printer_list([taxon_dict[j].sci_name for j in good_lineage_list[0]],
                                   head="Results:\t", sep="; ")

        print("\n" + lineage_txt + "\n")

        rename_dict = {}
        for i in tree_dict:
            sci_name = taxon_dict[i].sci_name
            if tree_dict[i].is_terminal() and i in taxon_hash:
                sci_name = sci_name + ": %.1f" % taxon_hash[i]

            rename_dict[i] = sci_name

        draw_ascii(common_tree, column_width=160,
                   clade_name=True, rename_dict=rename_dict)

    if good_lineage_list == []:
        recommand_taxon = MRCA
        recommand_taxon_fisher = 1.0
    else:
        recommand_taxon = good_lineage_list[0][-1]
        recommand_taxon_fisher = fisher_qvalue_dict[good_lineage_list[0][-1]]

    return MRCA, recommand_taxon, recommand_taxon_fisher
