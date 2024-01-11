def get_versions():
    return versions[0]["number"]


versions = [
    {
        "number": "0.2.1",
        "features": [
            "1. update README.md",
        ],
    },
    {
        "number": "0.2.0",
        "features": [
            "1. total rewrite",
            "2. only ID2Lineage and Name2Lineage are available",
        ],
    },
    {
        "number": "0.1.1",
        "features": [
            "1. Remove branch length and confidence in common tree",
        ],
    },
    {
        "number": "0.1.0",
        "features": [
            "1. Separate TaxonTools from ToolBiox",
        ],
    },
]
