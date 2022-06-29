# TaxonTools
TaxonTools can help you use [NCBI Taxonomy database](https://www.ncbi.nlm.nih.gov/taxonomy)

## Installation

```
pip install taxontools
```

## usage

Download NCBI Taxonomy database
```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
```

Convert NCBI database to sqlite format
```
tar -zxvf taxdump.tar.gz
TaxonTools StoreTaxonDB ./taxdump taxdump.db
```

# TODO

```

    # --------------------------------------------
    # command detail

    if args_dict["subcommand_name"] == "ID2Lineage":
        from toolbiox.src.xuyuxing.tools.taxontools import ID2Lineage_main
        ID2Lineage_main(args)

    elif args_dict["subcommand_name"] == "StoreTaxonDB":
        from toolbiox.lib.xuyuxing.evolution.taxonomy import build_taxon_database, store_tax_record_into_sqlite
        tax_record_dict = build_taxon_database(args.tax_dir)
        store_tax_record_into_sqlite(tax_record_dict, args.tax_db_file)

    elif args_dict["subcommand_name"] == "CommonTree":
        from toolbiox.src.xuyuxing.tools.taxontools import CommonTree_main
        CommonTree_main(args)

    elif args_dict["subcommand_name"] == "Name2Lineage":
        from toolbiox.src.xuyuxing.tools.taxontools import Name2Lineage_main
        Name2Lineage_main(args)

    elif args_dict["subcommand_name"] == "RankStat":
        from toolbiox.src.xuyuxing.tools.taxontools import RankStat_main
        RankStat_main(args)

    elif args_dict["subcommand_name"] == "ExcludeTaxon":
        """
        class abc():
            pass

        args = abc()
        args.tax_dir = '/lustre/home/xuyuxing/Database/genome2020/genome/info/NCBI/taxonomy'
        args.root_taxon_id = '73496'
        args.taxon_id = '40553'
        """
        from toolbiox.src.xuyuxing.tools.taxontools import ExcludeTaxon_main
        ExcludeTaxon_main(args)

    elif args_dict["subcommand_name"] == "DiamondTaxonAssign":
        """
        class abc():
            pass

        args = abc()
        args.blast_results = '/lustre/home/xuyuxing/Database/Cuscuta/Cau/genomev1.1/gene_taxon/Cuscuta.pt.v1.1.fasta.bls'
        args.tax_dir = '/lustre/home/xuyuxing/Database/NCBI/nr/2020/taxdmp'
        args.output_file = '/lustre/home/xuyuxing/Database/Cuscuta/Cau/genomev1.1/gene_taxon/Cuscuta.pt.v1.1.fasta.taxon'
        args.num_threads = 56
        """
        from toolbiox.src.xuyuxing.tools.taxontools import DiamondTaxonAssign_main
        DiamondTaxonAssign_main(args)

```
