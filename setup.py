# coding utf8
import setuptools
from taxontools.versions import get_versions

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

setuptools.setup(
    name="taxontools",
    version=get_versions(),
    author="Yuxing Xu",
    author_email="xuyuxing@mail.kib.ac.cn",
    description="A toolkit for analyzing Taxonomy",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url="https://github.com/SouthernCD/taxontools",

    entry_points={
        "console_scripts": ["TaxonTools = taxontools.cli:main"]
    },    

    package_data={
        "taxontools": ['one_kp/*'],
    },

    packages=setuptools.find_packages(),

    install_requires=[
        "yxtree",
        "yxutil",
        "yxmath",
        "mlxtend>=0.17.2",
        "scipy>=1.4.1",
        "numpy>=1.18.1",
        "biopython<=1.80",        
        # "thefuzz",
    ],

    python_requires='>=3.5',
)
