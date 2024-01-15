from Bio.Phylo.BaseTree import Tree, Clade
from itertools import combinations
from mlxtend.evaluate import permutation_test
from scipy import stats
from toolbiox.lib.common.evolution.tree_operate import Clade, Phylo, reroot_by_outgroup_clade, add_clade_name, lookup_by_names, get_MRCA, remove_given_node_from_tree, remove_pass_node, get_offspring, draw_ascii, get_parent, get_ancestors, get_sons
from toolbiox.lib.common.fileIO import tsv_file_dict_parse
from toolbiox.lib.common.math.set import merge_subset
from toolbiox.lib.common.util import printer_list
from toolbiox.lib.xuyuxing.base.common_command import log_print
from toolbiox.lib.xuyuxing.math.lcs import lcs
from toolbiox.lib.xuyuxing.math.set_operating import uniqify
from toolbiox.lib.xuyuxing.math.stats import get_qvalue
import copy
import numpy as np
import re
import sqlite3
import tarfile
import time
import toolbiox.lib.common.sqlite_command as sc

# this list is from https://en.wikipedia.org/wiki/Taxonomic_rank
RANK_LIST = [
    "superkingdom",
    "kingdom",
    "subkingdom",
    "superphylum",
    "phylum",
    "subphylum",
    "superclass",
    "class",
    "subclass",
    "infraclass",
    "cohort",
    "subcohort",
    "superorder",
    "order",
    "suborder",
    "infraorder",
    "parvorder",
    "superfamily",
    "family",
    "subfamily",
    "tribe",
    "subtribe",
    "genus",
    "subgenus",
    "section",
    "subsection",
    "series",
    "species group",
    "species subgroup",
    "species",
    "subspecies",
    "varietas",
    "forma",
    "strain",
    "genotype"]

OTHER_RANK_STRING = [
    "no rank",
    "clade",
    "serogroup",
    "biotype",
    "forma specialis",
    "isolate",
    "serotype",
    "morph",
    "pathogroup"]

SQLITE_DB_COLUMN_LIST = ["tax_id", "sci_name", "parent_tax_id",
                         "rank", "merged_new_id", "lineage"]


def check_new_rank_name(taxdump_db_file):
    tax_id_list = [
        i[0] for i in sc.sqlite_select(
            taxdump_db_file,
            'tax_record',
            column_list=['tax_id'])]
    record_dict = read_taxon_record_dict_db(taxdump_db_file, tax_id_list)

    for i in record_dict:
        temp_list = tuple(uniqify(
            [j[1] for j in record_dict[i].get_lineage if j[1] not in OTHER_RANK_STRING]))
        LCS_list = tuple(lcs(temp_list, RANK_LIST)[0])
        if not temp_list == LCS_list:
            print(temp_list)
            raise ValueError


def get_all_taxon_id_from_db_file(taxdump_db_file):
    return [
        i[0] for i in sc.sqlite_select(
            taxdump_db_file,
            'tax_record',
            column_list=['tax_id'])]


class Taxon(object):
    """
    Taxon is a class for each taxon in taxonomy database.
    """

    def __init__(
            self,
            tax_id=None,
            sci_name=None,
            parent_tax_id=None,
            rank=None,
            lineage=None,
            other_name=False):
        self.tax_id = tax_id
        self.sci_name = sci_name
        self.parent_tax_id = parent_tax_id
        self.rank = rank
        self.lineage = lineage
        self.other_name = other_name

    def __str__(self):
        return "ID: %s Sci_Name: %s (%s)" % (
            self.tax_id, self.sci_name, self.rank)

    def get_lineage(self, taxon_dict):
        """
        Get all lineage information from root to given taxon
        :param taxon_dict: The output from func "build_taxon_database"
        :return: a lineage list from root to given taxon
        """
        lineage_list = [(self.tax_id, self.rank)]
        iter_taxon = self
        while not iter_taxon.tax_id == iter_taxon.parent_tax_id:
            iter_taxon = taxon_dict[iter_taxon.parent_tax_id]
            lineage_list.append((iter_taxon.tax_id, iter_taxon.rank))

        self.lineage = lineage_list[::-1]

    def is_a(self, taxon_id):
        if taxon_id in [i[0] for i in self.lineage]:
            return True
        else:
            return False

    def lower_than(self, rank):
        temp_list = tuple([j[1] for j in self.lineage if j[1]
                           != 'no rank' and j[1] != 'clade'])

        if len(temp_list) == 0:  # Other sequence
            return False

        max_index = max(Taxon.rank_list.index(i) for i in temp_list)

        if max_index >= Taxon.rank_list.index(rank):
            return True
        else:
            return False

    def is_child_of(self, taxon_id):
        parent_list = self.lineage[:-1]
        if taxon_id in [i[0] for i in parent_list]:
            return True
        else:
            return False

    def is_ancestor_of(self, target_taxon):
        parent_list = target_taxon.lineage[:-1]
        if self.tax_id in [i[0] for i in parent_list]:
            return True
        else:
            return False

    def get_sons(self, taxdump_db_file):
        sons_id_list = []
        record_dict = read_taxon_record_dict_db(taxdump_db_file)
        for i in record_dict:
            i = record_dict[i]
            if self.tax_id == i.parent_tax_id:
                sons_id_list.append(i.tax_id)
        return sons_id_list

    def get_offspring(self, taxdump_db_file):
        offspring_id_list = []
        record_dict = read_taxon_record_dict_db(taxdump_db_file)
        for i in record_dict:
            i = record_dict[i]
            if self.tax_id in [j[0] for j in i.lineage]:
                if not i.tax_id == self.tax_id:
                    offspring_id_list.append(i.tax_id)
        return offspring_id_list


def build_taxon_database(taxdump_tar_gz_file, full_name_class=False):
    """
    NCBI taxonomy database can be download from: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    :param taxonomy_dir: file path for taxdump.tar.gz
    :return: a dict whose key is tax_id and value is a Taxon with attributes named "tax_id", "sci_name", "parent_tax_id\
    ", "rank" and "lineage" (Big memory)
    """
    tar = tarfile.open(taxdump_tar_gz_file, "r:gz")

    names_dmp_file = tar.extractfile('names.dmp')
    nodes_dmp_file = tar.extractfile('nodes.dmp')

    taxon_record_dict = {}
    for each_line in names_dmp_file.readlines():
        each_line = each_line.strip().decode()
        info = each_line.split("|")
        tax_id = re.sub('\t', '', info[0])
        name = re.sub('\t', '', info[1])
        # unique_name = re.sub('\t', '', info[2])
        name_class = re.sub('\t', '', info[3])

        if tax_id not in taxon_record_dict:
            taxon_record_dict[tax_id] = Taxon()
            taxon_record_dict[tax_id].tax_id = tax_id

        if name_class == "scientific name":
            taxon_record_dict[tax_id].sci_name = name
        elif full_name_class:
            if not hasattr(taxon_record_dict[tax_id], 'other_name'):
                taxon_record_dict[tax_id].other_name = []
            taxon_record_dict[tax_id].other_name.append(name)

    for each_line in nodes_dmp_file.readlines():
        each_line = each_line.strip().decode()
        info = each_line.split("|")
        tax_id = re.sub('\t', '', info[0])
        parent_tax_id = re.sub('\t', '', info[1])
        rank = re.sub('\t', '', info[2])
        taxon_record_dict[tax_id].parent_tax_id = parent_tax_id
        taxon_record_dict[tax_id].rank = rank

    for tax_id in taxon_record_dict:
        taxon = taxon_record_dict[tax_id]
        taxon.get_lineage(taxon_record_dict)
        # taxon.lineage = lineage_list
        for i in taxon.lineage:
            tax_id_tmp, rank_tmp = i
            if not (rank_tmp == 'no rank' or rank_tmp == 'clade'):
                setattr(taxon, rank_tmp, tax_id_tmp)

    merged_dmp_file = tar.extractfile('merged.dmp')

    merged_record = {}

    for each_line in merged_dmp_file.readlines():
        each_line = each_line.strip().decode()
        info = each_line.split("|")
        merged_id = re.sub('\t', '', info[0])
        new_id = re.sub('\t', '', info[1])
        merged_record[merged_id] = new_id

    for i in merged_record:
        taxon_record_dict[i] = copy.deepcopy(
            taxon_record_dict[merged_record[i]])
        taxon_record_dict[i].older = True

    tar.close()

    return taxon_record_dict


def store_taxon_record_into_sqlite(taxon_record_dict, taxdump_db_file):
    log_print("\t\tBegin: convert NCBI taxonomy database to sqlite3 file")
    start_time = time.time()
    sc.init_sql_db(taxdump_db_file, "tax_record", SQLITE_DB_COLUMN_LIST)
    num = 0
    record_tmp_dict = []
    for tax_id in taxon_record_dict:
        record_tmp = []
        record_tmp.append(tax_id)
        taxon_tmp = taxon_record_dict[tax_id]

        if hasattr(taxon_tmp, "sci_name"):
            record_tmp.append(getattr(taxon_tmp, "sci_name").replace("\"", ""))
        else:
            record_tmp.append("")

        if hasattr(taxon_tmp, "parent_tax_id"):
            record_tmp.append(getattr(taxon_tmp, "parent_tax_id"))
        else:
            record_tmp.append("")

        if hasattr(taxon_tmp, "rank"):
            record_tmp.append(getattr(taxon_tmp, "rank"))
        else:
            record_tmp.append("")

        if hasattr(taxon_tmp, "older"):
            record_tmp.append(getattr(taxon_tmp, "tax_id"))
        else:
            record_tmp.append("")

        taxon_tmp.get_lineage(taxon_record_dict)
        lineage_string = ""
        for tax_id_rank, rank in taxon_tmp.lineage:
            lineage_string = lineage_string + "%s:%s;" % (tax_id_rank, rank)
        lineage_string = lineage_string.rstrip(";")
        record_tmp.append(lineage_string)

        num = num + 1
        record_tmp_dict.append(record_tmp)
        # sc.insert_one_record_to_sql_table(tuple(record_tmp), tuple(column_name), db_file, "tax_record")

        if num % 10000 == 0:
            sc.sqlite_write(record_tmp_dict, taxdump_db_file,
                            "tax_record", SQLITE_DB_COLUMN_LIST)
            record_tmp_dict = []

        round_time = time.time()
        if round_time - start_time > 10:
            log_print("\t\t\t%d/%d" % (num, len(taxon_record_dict)))
            start_time = round_time

    if len(record_tmp_dict) > 0:
        sc.sqlite_write(record_tmp_dict, taxdump_db_file,
                        "tax_record", SQLITE_DB_COLUMN_LIST)
        record_tmp_dict = []
        log_print("\t\t\t%d/%d" % (num, len(taxon_record_dict)))

    conn = sqlite3.connect(taxdump_db_file)
    conn.execute("CREATE UNIQUE INDEX tax_id_index on tax_record (tax_id)")
    conn.close()

    log_print("\t\tEnd: convert NCBI taxonomy database to sqlite3 file")


def read_taxon_record_dict_db(taxdump_db_file, tax_id_list=[]):

    if len(tax_id_list) == 0:
        record_tuple_list = sc.sqlite_select(taxdump_db_file, "tax_record")
    else:
        record_tuple_list = sc.sqlite_select_by_a_key(
            taxdump_db_file, "tax_record", "tax_id", tuple(tax_id_list))

    record_dict = {}
    for record_tuple in record_tuple_list:
        record = Taxon()
        record.tax_id = record_tuple[0]
        if not record_tuple[1] == "":
            record.sci_name = record_tuple[1]
        if not record_tuple[2] == "":
            record.parent_tax_id = record_tuple[2]
        if not record_tuple[3] == "":
            record.rank = record_tuple[3]
        if not record_tuple[4] == "":
            record.tax_id = record_tuple[4]
            record.older = True

        record.lineage = [tuple(i.split(":"))
                          for i in record_tuple[5].split(";")]

        for i in record.lineage:
            tax_id_tmp, rank_tmp = i
            if not (rank_tmp == 'no rank' or rank_tmp == 'clade'):
                setattr(record, rank_tmp, tax_id_tmp)
        record_dict[record_tuple[0]] = record

    return record_dict


def read_taxon_name_record_dict_db(taxdump_db_file, tax_name_list=[]):

    record_tuple_list = sc.sqlite_select_by_a_key(
        taxdump_db_file, "tax_record", "sci_name", tuple(tax_name_list))

    record_dict = {}
    for record_tuple in record_tuple_list:
        record = Taxon()
        record.tax_id = record_tuple[0]
        if not record_tuple[1] == "":
            record.sci_name = record_tuple[1]
        if not record_tuple[2] == "":
            record.parent_tax_id = record_tuple[2]
        if not record_tuple[3] == "":
            record.rank = record_tuple[3]
        if not record_tuple[4] == "":
            record.tax_id = record_tuple[4]
            record.older = True

        record.lineage = [tuple(i.split(":"))
                          for i in record_tuple[5].split(";")]

        lineage_list = record.lineage
        # taxon.lineage = lineage_list
        for i in lineage_list:
            tax_id_tmp, rank_tmp = i
            if not (rank_tmp == 'no rank' or rank_tmp == 'clade'):
                setattr(record, rank_tmp, tax_id_tmp)
        record_dict[record_tuple[1]] = record

    return record_dict


def get_MRCA_taxon_id(taxon_list):
    if len(taxon_list) == 1:
        return taxon_list[0]

    lineage_sum_list = []
    for tax_tmp in taxon_list:
        lineage_sum_list.append(tax_tmp.lineage)

    range_num = min([len(i) for i in lineage_sum_list])

    MRCA = None
    for i in range(0, range_num):
        tmp_taxon_list = [j[i] for j in lineage_sum_list]
        if len(set(tmp_taxon_list)) == 1:
            MRCA = tmp_taxon_list[0]

    return MRCA[0]


def get_common_tree(taxon_id_list, taxdump_db_file):
    taxon_dict = read_taxon_record_dict_db(
        taxdump_db_file, tax_id_list=tax_id_list)
    MRCA_id = get_MRCA_taxon_id([taxon_dict[i] for i in taxon_dict])

    clade_MRCA = Clade(branch_length=None, name=MRCA_id,
                       clades=None, confidence=None, color=None, width=None)
    clade_dict = {
        MRCA_id: clade_MRCA
    }

    for taxon_id in taxon_id_list:
        lineage = [i[0] for i in taxon_dict[taxon_id].lineage]

        for taxon_id_tmp in lineage:
            taxon_id_tmp_index = lineage.index(taxon_id_tmp)
            # old than MRCA
            if taxon_id_tmp_index < lineage.index(MRCA_id):
                continue
            # leaf taxon
            if taxon_id_tmp_index + 1 >= len(lineage):
                continue
            # good now add son to dict
            my_clade = clade_dict[taxon_id_tmp]
            son_taxon_id = lineage[taxon_id_tmp_index + 1]
            if my_clade.clades is None:
                my_clade.clades = []
            already_have_son = [i.name for i in my_clade.clades]
            if son_taxon_id not in already_have_son:
                son_clade = Clade(name=son_taxon_id)
                my_clade.clades.append(son_clade)
                clade_dict[son_taxon_id] = son_clade

    tree = Tree.from_clade(clade_MRCA)

    return tree


if __name__ == '__main__':
    taxdump_tar_gz_file = "/lustre/home/xuyuxing/Database/NCBI/taxonomy/taxdump.tar.gz"
    taxdump_db_file = "/lustre/home/xuyuxing/Database/NCBI/taxonomy/taxdump.db"

    # convert NCBI taxonomy database to sqlite3 file
    taxon_record_dict = build_taxon_database(taxdump_tar_gz_file)
    store_taxon_record_into_sqlite(taxon_record_dict, taxdump_db_file)

    # extract taxon record by taxon ID
    tax_id = '3702'
    taxon_dict = read_taxon_record_dict_db(
        taxdump_db_file, tax_id_list=[tax_id])
    vars(taxon_dict[tax_id])

    # extract taxon record by taxon name
    tax_name = 'Arabidopsis thaliana'
    taxon_dict = read_taxon_name_record_dict_db(
        taxdump_db_file, tax_name_list=[tax_name])
    vars(taxon_dict[tax_name])

    # get MRCA of a taxon list and common tree
    tax_id_list = ['3702', '4128']
    taxon_dict = read_taxon_record_dict_db(
        taxdump_db_file, tax_id_list=tax_id_list)
    MRCA_id = get_MRCA_taxon_id([taxon_dict[i] for i in taxon_dict])
    print(read_taxon_record_dict_db(
        taxdump_db_file, tax_id_list=[MRCA_id])[MRCA_id])

    # get common tree
    tax_id_list = ['3702', '4128', '4081']
    common_tree = get_common_tree(tax_id_list, taxdump_db_file)
