import yaml
import os
import argparse

class Config:
    def __init__(self, config):
        self.organism = config['organism']
        self.cpus = config['cpus']
        self.structures = config['structures'] if config['structures']['enabled'] else False
        self.metabolism_pathwaytools = config['metabolism-PathwayTools'] if config['metabolism-PathwayTools']['enabled'] else False
        self.metabolism_sbml = config['metabolism-SBML'] if config['metabolism-SBML']['enabled'] else False
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
        'organism', 'structures', 'metabolism-PathwayTools', 'metabolism-SBML',
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
    print(f"Species Tax ID: {config.organism['tax_id']}")
    print(f"Strain Tax ID: {config.organism['strain_taxid']}")
    print(f"GBK File: {config.organism['gbk_file']}")
    
    if config.cpus:
        print(f"CPUS: {config.cpus}")
        
    print("\n---- Structures ----")
    if config.structures:
        print("Structures will be used in the analysis.")
    else:
        print(f"Structures Enabled: {config.structures}")
    
    print("\n---- Metabolism ----")
    if config.metabolism_pathwaytools:
        print("PathwayTools Metabolism analysis will be used.")
        print(f"SBML File: {config.metabolism_pathwaytools['sbml_file']}")
        print(f"Chokepoint File: {config.metabolism_pathwaytools['chokepoint_file']}")
        print(f"Smarttable File: {config.metabolism_pathwaytools['smarttable_file']}")
    else:
        print(f"Metabolism PathwayTools Enabled: {config.metabolism_pathwaytools}")

    if config.metabolism_sbml:
        print("SBML Metabolism analysis will be used.")
        print(f"SBML File: {config.metabolism_sbml['sbml_file']}")
        print(f"Filter File: {config.metabolism_sbml['filter_file']}")
    else:
        print(f"Metabolism SBML Enabled: {config.metabolism_sbml}")

    print("\n---- Core ----")
    
    if config.core:
        print(f"Roary Enabled: {config.core['roary']}")       
        print(f"CoreCruncher Enabled: {config.core['corecruncher']}")
        print(f"Minimum Identity: {config.core['min_identity']}%")
        print(f"Minimum Core Frequency: {config.core['min_core_freq']}%")
    else:
        print(f"Core Enabled: {config.core}")
      
    print("\n---- Offtarget ----")
    
    if config.offtarget:
        print(f"Human Offtarget Enabled: {config.offtarget['human']}")
        print(f"Microbiome Offtarget Enabled: {config.offtarget['microbiome']}")
        if config.offtarget['microbiome']:
            print(f"Microbiome Identity Filter: {config.offtarget['microbiome_identity_filter']}")
            print(f"Microbiome Coverage Filter: {config.offtarget['microbiome_coverage_filter']}")
        if config.structures:
            print(f"Foldseek Human Offtarget Enabled: {config.offtarget['foldseek_human']}")
        else:
            print('Structures must be enabled to use Foldseek')
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

    # Add summary section
    print("\n" + "="*50)
    print("CONFIGURATION SUMMARY")
    print("="*50)
    enabled_modules = []
    if config.metabolism_pathwaytools or config.metabolism_sbml: enabled_modules.append("Metabolism")
    if config.structures: enabled_modules.append("Structures") 
    if config.core: enabled_modules.append("Core Analysis")
    if config.offtarget: enabled_modules.append("Off-target")
    if config.deg: enabled_modules.append("DEG")
    if config.psortb: enabled_modules.append("Localization")
    if config.metadata: enabled_modules.append("Metadata")
    
    print(f"Organism: {config.organism['name']}")
    print(f"Enabled modules ({len(enabled_modules)}): {', '.join(enabled_modules)}")
    print("="*50)


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
    parser = argparse.ArgumentParser(description='Validate and display FastTarget configuration.')
    parser.add_argument('--config_file', type=str, default='config.yml', help='Path to the configuration file')
    args = parser.parse_args()
    
    try:
        config = get_config(args.config_file)
        print_config(config)
        print('\n✅ Configuration is valid!')
        print('To modify settings, edit config.yml and run this script again.')
        print('To run the pipeline: python fasttarget.py --config_file config.yml')
        exit(0)
    except Exception as e:
        print(f'\n❌ Configuration error: {e}')
        print('Please fix config.yml and try again.')
        exit(1)
