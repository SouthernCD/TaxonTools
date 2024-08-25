from yxtree import Clade, Phylo, reroot_by_outgroup_clade, add_clade_name, lookup_by_names, get_MRCA, remove_given_node_from_tree
from yxutil import tsv_file_dict_parse
from .taxon import read_taxon_record_dict_db

# use 1kp
def load_one_kp_tree(one_kp_tree_file, acc_to_taxon_map):
    acc_map_dict = tsv_file_dict_parse(
        acc_to_taxon_map, fieldnames=[
            'acc', 'tax_id', 'lineage'], key_col='acc')

    one_kp_tree = Phylo.read(one_kp_tree_file, 'newick')
    one_kp_tree = add_clade_name(one_kp_tree)
    one_kp_node_dict = lookup_by_names(one_kp_tree)

    outgroup = ['DZPJ', 'RAWF']
    out_clade = get_MRCA(
        one_kp_node_dict[outgroup[0]], one_kp_node_dict[outgroup[1]], one_kp_tree)

    one_kp_rooted_tree, one_kp_rooted_node_dict = reroot_by_outgroup_clade(
        one_kp_tree, one_kp_node_dict, out_clade.name, False)

    return one_kp_rooted_tree, one_kp_rooted_node_dict, acc_map_dict


def extract_common_tree(one_kp_tree, keep_leave_list):

    node_name_list = [i.name for i in one_kp_tree.get_terminals()]
    valid_keep_leave_list = list(set(keep_leave_list) & set(node_name_list))
    removed_node_list = list(set(node_name_list) - set(valid_keep_leave_list))

    removed_tree = remove_given_node_from_tree(one_kp_tree, removed_node_list)

    for c in removed_tree.find_clades():
        if c.branch_length:
            c.branch_length = None

    return removed_tree


def get_best_1kp_represent(input_taxon_id_list, acc_map_dict, taxon_db_file):

    acc_taxon_list = list(set([acc_map_dict[i]['tax_id']
                               for i in acc_map_dict]))
    taxon_dict = read_taxon_record_dict_db(
        taxon_db_file, tax_id_list=input_taxon_id_list + acc_taxon_list)

    rep_dict = {}
    for i in input_taxon_id_list:

        if i in set(acc_taxon_list):
            rep_dict[i] = i
        else:
            it = taxon_dict[i]
            it_lineage = set([j[0] for j in it.get_lineage])

            ca_dict = {}
            for a in acc_taxon_list:
                at = taxon_dict[a]
                at_lineage = set([j[0] for j in at.get_lineage])
                ca_dict[a] = list(it_lineage & at_lineage)

            best_rep_ca_num = len(
                ca_dict[sorted(acc_taxon_list, key=lambda x: len(ca_dict[x]), reverse=True)[0]])
            best_rep = sorted([j for j in acc_taxon_list if len(
                ca_dict[j]) == best_rep_ca_num])[0]
            rep_dict[i] = best_rep

    for i in rep_dict:
        rep_dict[i] = sorted(
            [j for j in acc_map_dict if acc_map_dict[j]['tax_id'] == rep_dict[i]])[0]

    return rep_dict


def get_common_tree_by_1kp(
        input_taxon_id_list,
        one_kp_tree_file,
        acc_to_taxon_map,
        taxon_db_file):

    one_kp_rooted_tree, one_kp_rooted_node_dict, acc_map_dict = load_one_kp_tree(
        one_kp_tree_file, acc_to_taxon_map)
    rep_dict = get_best_1kp_represent(
        input_taxon_id_list, acc_map_dict, taxon_db_file)

    removed_tree = extract_common_tree(
        one_kp_rooted_tree, set([rep_dict[i] for i in rep_dict]))
    removed_tree = add_clade_name(removed_tree)
    removed_tree_dict = lookup_by_names(removed_tree)

    okp2input_dict = {}
    for i in rep_dict:
        okp2input_dict.setdefault(rep_dict[i], []).append(i)

    for okp in okp2input_dict:
        input_it_list = okp2input_dict[okp]
        if len(input_it_list) > 1:
            clade_list = [Clade(name=i, clades=[]) for i in input_it_list]
            removed_tree_dict[okp].clades = clade_list
        else:
            removed_tree_dict[okp].name = input_it_list[0]

    return removed_tree
