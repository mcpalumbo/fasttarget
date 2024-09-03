import yaml
import os
import argparse

class Config:
    def __init__(self, config):
        self.organism = config['organism']
        self.cpus = config['cpus']
        self.structures = config['structures'] if config['structures']['enabled'] else False
        self.metabolism = config['metabolism'] if config['metabolism']['enabled'] else False
        self.core = config['core'] if config['core']['enabled'] else False
        self.metadata = config['metadata'] if config['metadata']['enabled'] else False
        self.offtarget = config['offtarget'] if config['offtarget']['enabled'] else False
        self.deg = config['deg'] if config['deg']['enabled'] else False
        self.psortb = config['psortb'] if config['psortb']['enabled'] else False

def load_config(config_path):
    """
    Load configuration from a YAML file.

    :param config_path: Path to the configuration file.

    :return: Configuration dictionary.
    """
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def validate_config(config):
    """
    Validate the configuration dictionary.
    
    :param config: Configuration dictionary.
    """
    required_keys = [
        'organism', 'structures', 'metabolism', 
        'core', 'metadata', 'offtarget', 'deg', 'psortb'
    ]
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required config key: {key}")

def print_config(config):
    """
    Print the configuration to the console.
    
    :param config: Configuration object.
    """

    print("Configuration:")
       
    print("\n---- Organism information ----")
    print(f"Organism Name: {config.organism['name']}")
    print(f"Tax ID: {config.organism['tax_id']}")
    print(f"GBK File: {config.organism['gbk_file']}")
    
    if config.cpus:
        print(f"CPUS: {config.cpus}")
        
    print("\n---- Structures ----")
    if config.structures:
        print(f"Proteome Uniprot: {config.structures['proteome_uniprot']}")
    else:
        print(f"Structures Enabled: {config.structures}")
    
    print("\n---- Metabolism ----")
    if config.metabolism:
        print(f"SBML File: {config.metabolism['sbml_file']}")
        print(f"Chokepoint File: {config.metabolism['chokepoint_file']}")
        print(f"Smarttable File: {config.metabolism['smarttable_file']}")
    else:
        print(f"Metabolism Enabled: {config.metabolism}")

    print("\n---- Core ----")
    
    if config.core:
        print(f"Roary Enabled: {config.core['roary']}")       
        print(f"CoreCruncher Enabled: {config.core['corecruncher']}")
        if config.core['corecruncher']:
            print(f"CoreCruncher Script: {config.core['corecruncher_script']}")
    else:
        print(f"Core Enabled: {config.core}")
      
    print("\n---- Offtarget ----")
    
    if config.offtarget:
        print(f"Human Offtarget Enabled: {config.offtarget['human']}")
        print(f"Microbiome Offtarget Enabled: {config.offtarget['microbiome']}")
        if config.offtarget['microbiome']:
            print(f"Microbiome Identity Filter: {config.offtarget['microbiome_identity_filter']}")
            print(f"Microbiome Coverage Filter: {config.offtarget['microbiome_coverage_filter']}")
    else:
        print(f"Offtarget Enabled: {config.offtarget}")
        
    print("\n---- DEG ----")
    
    if config.deg:
        print(f"DEG Identity Filter: {config.deg['deg_identity_filter']}")
        print(f"DEG Coverage Filter: {config.deg['deg_coverage_filter']}")
    else:
        print(f"DEG Enabled: {config.deg}")

    
    print("\n---- Localization ----")
    if config.psortb:
        print(f"Gram Type: {config.psortb['gram_type']}")
    else:
        print(f"PSortB Enabled: {config.psortb}")

    print("\n---- Metadata ----")
    if config.metadata:
        for table in config.metadata['meta_tables']:
            print(f"Meta Table: {table}")
    else:
        print(f"Metadata Enabled: {config.metadata}")


def get_config(config_path):
    """
    Load and validate the configuration.

    :param config_path: Path to the configuration file.
    
    :return: Configuration object.
    """

    base_path = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(base_path, config_path)
    config = load_config(config_path)
    validate_config(config)
    return Config(config)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some configuration.')
    parser.add_argument('--config_file', type=str, default='config.yml', help='Path to the configuration file')
    args = parser.parse_args()
    
    config = get_config(args.config_file)
    print_config(config)
    print('Please modify config.yml to change data.')
