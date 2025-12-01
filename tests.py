import configuration
import fasttarget
import os
import pandas as pd
import argparse
import multiprocessing
from ftscripts import files, structures, pathways, offtargets, genome, essentiality, metadata
from datetime import datetime
import logging
import sys
from ftscripts.logger import logger 
import shutil

def run_test(test_config, base_path):
    try:
        print("Running test script...")
        result = fasttarget.main(test_config, base_path)
        print(f"Test script completed.")
        
        return result
        
    except Exception as e:
        print(f"An error occurred: {e}")
        raise


if __name__ == "__main__":

    base_path = os.path.dirname(os.path.abspath(__file__))
    test_path = os.path.join(base_path, 'organism', 'test')
    
    config_dict = {
        "organism": {
            "name": "test",
            "tax_id": 2104,
            "strain_taxid": 272634,
            "gbk_file": f'{test_path}/test.gbk'
        },
        "cpus": None,
        "metabolism-PathwayTools": {
            "enabled": True,
            "sbml_file": f'{test_path}/test.sbml',
            "chokepoint_file": f'{test_path}/test_chokepoints.txt',
            "smarttable_file": f'{test_path}/test_genes_smarttable.tsv'
        },
        "metabolism-SBML": {
            "enabled": False,
            "sbml_file": f'{test_path}/test.sbml',
            "filter_file": f'{test_path}/test_ubiquitous.txt'
        },
        "structures": {
            "enabled": True
        },
        "core": {
            "enabled": True,
            "roary": True,
            "corecruncher": True,
            "min_identity": 95,
            "min_core_freq": 99
        },
        "offtarget": {
            "enabled": True,
            "human": True,
            "microbiome": True,
            "microbiome_identity_filter": 50,
            "microbiome_coverage_filter": 70,
            "foldseek_human": True
        },
        "deg": {
            "enabled": True,
            "deg_identity_filter": 50,
            "deg_coverage_filter": 70
        },
        "psortb": {
            "enabled": True,
            "gram_type": "n"
        },
        "metadata": {
            "enabled": True,
            "meta_tables": [f'{test_path}/test_meta_table1.txt']
        }
    }

    test_config = configuration.Config(config_dict)
    run_test(test_config, base_path)