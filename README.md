# TaxonTools
TaxonTools can help you use [NCBI Taxonomy database](https://www.ncbi.nlm.nih.gov/taxonomy)

## Installation
```
pip install taxontools
```

## Usage

**StoreTaxonDB**: convert NCBI Taxonomy database to sqlite format

```
usage: TaxonTools StoreTaxonDB [-h] taxdump_tar_gz_file taxdump_db_file

positional arguments:
  taxdump_tar_gz_file  Path for taxdump.tar.gz file
  taxdump_db_file      Path for sqlite3 file

optional arguments:
  -h, --help           show this help message and exit
```

**ID2Lineage**: extract taxonomy information and lineage information from NCBI Taxonomy database by taxon ID

```
usage: TaxonTools ID2Lineage [-h] [-l TAXON_ID_LIST] [-o OUTPUT_FILE] taxdump_db_file

positional arguments:
  taxdump_db_file       Path for taxdump_db_file from StoreTaxonDB

optional arguments:
  -h, --help            show this help message and exit
  -l TAXON_ID_LIST, --taxon_id_list TAXON_ID_LIST
                        query taxon ID, you can give me a file path or taxon IDs separate by comma
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output file name (default as $STDOUT)
```

**Name2Lineage**: extract taxonomy information and lineage information from NCBI Taxonomy database by taxon Name

```
usage: TaxonTools Name2Lineage [-h] [-l TAXON_NAME_LIST] [-o OUTPUT_FILE] taxdump_db_file

positional arguments:
  taxdump_db_file       Path for taxdump_db_file from StoreTaxonDB

optional arguments:
  -h, --help            show this help message and exit
  -l TAXON_NAME_LIST, --taxon_name_list TAXON_NAME_LIST
                        query taxon Name, you can give me a file path or taxon names separate by comma
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output file name (default as $STDOUT)
```

## Example
1. Download NCBI Taxonomy database
```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
```

2. Convert NCBI Taxonomy database to sqlite format
```
TaxonTools StoreTaxonDB taxdump.tar.gz taxdump.db
```

3. Extract taxonomy information and lineage information from NCBI Taxonomy database by taxon ID
```
# query taxon ID 9606,10090,10116
TaxonTools ID2Lineage -l 9606,10090,10116 -o taxon_lineage.txt taxdump.db
# query taxon ID from file
TaxonTools ID2Lineage -l taxon_id_list.txt -o taxon_lineage.txt taxdump.db
```

4. Extract taxonomy information and lineage information from NCBI Taxonomy database by taxon Name
```
# query taxon ID Arabidopsis thaliana, Mus musculus, Rattus norvegicus
TaxonTools Name2Lineage -l "Arabidopsis\ thaliana,Mus\ musculus,Rattus\ norvegicus" -o taxon_lineage.txt taxdump.db
# query taxon ID from file
TaxonTools Name2Lineage -l taxon_name_list.txt -o taxon_lineage.txt taxdump.db
```

5. Use in python
- convert NCBI taxonomy database to sqlite3 file
```
from taxontools import build_taxon_database, store_taxon_record_into_sqlite

taxon_record_dict = build_taxon_database(taxdump_tar_gz_file)
store_taxon_record_into_sqlite(taxon_record_dict, taxdump_db_file)
```
- extract taxon record by taxon ID
```
from taxontools import read_taxon_record_dict_db

tax_id = '3702'
taxon_dict = read_taxon_record_dict_db(
    taxdump_db_file, tax_id_list=[tax_id])
vars(taxon_dict[tax_id])
```
- extract taxon record by taxon name
```
from taxontools import read_taxon_name_record_dict_db

tax_name = 'Arabidopsis thaliana'
taxon_dict = read_taxon_name_record_dict_db(
    taxdump_db_file, tax_name_list=[tax_name])
vars(taxon_dict[tax_name])
```
- get MRCA of a taxon list and common tree
```
from taxontools import get_MRCA_taxon_id

tax_id_list = ['3702', '4128']
taxon_dict = read_taxon_record_dict_db(
    taxdump_db_file, tax_id_list=tax_id_list)
MRCA_id = get_MRCA_taxon_id([taxon_dict[i] for i in taxon_dict])
print(read_taxon_record_dict_db(
    taxdump_db_file, tax_id_list=[MRCA_id])[MRCA_id])
```
- get common tree
```
from taxontools import get_common_tree
from Bio import Phylo

tax_id_list = ['3702', '4128', '4081']
common_tree = get_common_tree(tax_id_list, taxdump_db_file)
Phylo.draw_ascii(common_tree)
```