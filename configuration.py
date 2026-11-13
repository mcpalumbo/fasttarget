import yaml
import os
import argparse
from ftscripts.microbiome_catalogues import MGNIFY_CATALOGUES

def validate_file_path(filepath, name, required_extensions=None):
    """
    Validate that a file path exists and optionally has the correct extension.
    
    :param filepath: Path to validate
    :param name: Human-readable name for error messages
    :param required_extensions: List of valid extensions (e.g., ['.gbk', '.gb'])
    :return: True if valid, error message string if invalid
    """
    if not isinstance(filepath, str) or not filepath.strip():
        return f"{name} must be a non-empty string"
    
    if not os.path.exists(filepath):
        return f"{name} not found: {filepath}"
    
    if not os.path.isfile(filepath):
        return f"{name} is not a file: {filepath}"
    
    if required_extensions:
        if not any(filepath.lower().endswith(ext) for ext in required_extensions):
            ext_str = ", ".join(required_extensions)
            return f"{name} should have one of these extensions: {ext_str} (got: {filepath})"
    
    return None  # Valid

class Config:
    """
    Validated configuration, and the single source of default values.

    Defaults live here and nowhere else. They used to be applied both here and
    inside validate_config, which mutated the dictionary it was validating, so
    the two copies could drift apart and the validated object depended on
    validation having run first.

    Optional sections and keys are read with .get() so that a config omitting
    them produces a working object rather than a bare KeyError. 'cpus' and
    'colabfold' in particular are not required keys, so reading them directly
    crashed on a config that had passed validation.
    """

    # Defaults applied to a section when it is enabled but the key is absent.
    SECTION_DEFAULTS = {
        'structures': {'pocket_full_mode': False, 'uniparc_fallback': False},
        'colabfold': {'amber': False, 'gpu': False, 'colabfold_run_all': False},
        'core': {'panaroo': False, 'panaroo_mode': 'strict'},
    }

    @staticmethod
    def _section(config, name, defaults=None):
        """
        Return a copy of an enabled config section, or False when it is disabled.

        :param config: Configuration dictionary.
        :param name: Section name.
        :param defaults: Optional mapping of default values for missing keys.
        :return: Dict for an enabled section, False otherwise.
        """
        section = config.get(name) or {}
        if not isinstance(section, dict) or not section.get('enabled', False):
            return False

        section = section.copy()
        for key, value in (defaults or {}).items():
            section.setdefault(key, value)
        return section

    def __init__(self, config):
        self.organism = config['organism']
        self.cpus = config.get('cpus')
        self.container_engine = config.get('container_engine', 'docker')  # Default to docker if not specified

        self.structures = self._section(config, 'structures', self.SECTION_DEFAULTS['structures'])
        self.colabfold = self._section(config, 'colabfold', self.SECTION_DEFAULTS['colabfold'])
        self.metabolism_pathwaytools = self._section(config, 'metabolism-PathwayTools')
        self.metabolism_sbml = self._section(config, 'metabolism-SBML')
        self.core = self._section(config, 'core', self.SECTION_DEFAULTS['core'])
        self.metadata = self._section(config, 'metadata')
        self.offtarget = self._section(config, 'offtarget')
        self.deg = self._section(config, 'deg')
        self.psortb = self._section(config, 'psortb')
        ligands_config = config.get('ligands', {})
        self.ligands = ligands_config.copy() if ligands_config.get('enabled', False) else False
        alphafill_config = self.ligands.get('alphafill', {}) if self.ligands else {}
        self.alphafill = alphafill_config.copy() if alphafill_config.get('enabled', False) else False
        for source in ('chembl', 'bindingdb'):
            section = self.ligands.get(source, {}) if self.ligands else {}
            setattr(self, source, section.copy() if section.get('enabled', False) else False)

def load_config(config_path):
    """
    Load configuration from a YAML file.

    :param config_path: Path to the configuration file.

    :return: Configuration dictionary.
    """
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def as_integer(value):
    """
    Interpret a configuration value as an integer, or report that it is not one.

    int() was used directly, which truncates: a tax_id of 1.5 became 1, the NCBI
    root node, and the run queried that instead of failing. A value is accepted
    only when nothing is lost, so 2104, '2104' and 2104.0 are the same number and
    2104.5 is not a number of a taxon.

    :param value: Value from the configuration.
    :return: The integer, or None when the value is not an exact integer.
    """
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value) if value.is_integer() else None
    if isinstance(value, str):
        text = value.strip()
        try:
            number = float(text)
        except ValueError:
            return None
        return int(number) if number.is_integer() else None
    return None


def validate_config(config):
    """
    Validate the configuration dictionary with comprehensive checks.
    
    :param config: Configuration dictionary.
    """
    errors = []
    
    # Check required top-level keys
    required_keys = [
        'organism', 'structures', 'metabolism-PathwayTools', 'metabolism-SBML',
        'core', 'metadata', 'offtarget', 'deg', 'psortb'
    ]
    for key in required_keys:
        if key not in config:
            errors.append(f"Missing required config key: '{key}'")
    
    # If critical keys missing, raise early
    if errors:
        raise ValueError("\n".join(errors))
    
    # Validate organism section
    if 'organism' in config:
        org = config['organism']
        
        # Check required organism fields
        required_org_fields = ['name', 'tax_id', 'strain_taxid', 'gbk_file']
        for field in required_org_fields:
            if field not in org:
                errors.append(f"Missing required organism field: '{field}'")
        
        # Validate organism name (no spaces, reasonable length)
        if 'name' in org:
            name = org['name']
            if not isinstance(name, str) or not name.strip():
                errors.append("Organism name must be a non-empty string")
            elif ' ' in name:
                errors.append(f"Organism name '{name}' should not contain spaces")
            elif len(name) > 50:
                errors.append(f"Organism name '{name}' is too long (max 50 characters)")
        
        # Validate tax IDs (must be positive integers)
        for tax_field in ['tax_id', 'strain_taxid']:
            if tax_field in org:
                tax_val = as_integer(org[tax_field])
                if tax_val is None:
                    errors.append(f"Organism {tax_field} must be an integer (got {org[tax_field]})")
                elif tax_val <= 0:
                    errors.append(f"Organism {tax_field} must be a positive integer (got {tax_val})")
        
        # Validate GBK file exists
        if 'gbk_file' in org:
            gbk_path = org['gbk_file']
            if not isinstance(gbk_path, str) or not gbk_path.strip():
                errors.append("GBK file path must be a non-empty string")
            elif not os.path.exists(gbk_path):
                errors.append(f"GBK file not found: {gbk_path}")
            elif not gbk_path.lower().endswith(('.gbk', '.gbff', '.gb', '.genbank')):
                errors.append(f"GBK file should have .gbk, .gbff, .gb, or .genbank extension: {gbk_path}")
    
    # Validate CPUs
    if 'cpus' in config and config['cpus'] is not None:
        cpus = as_integer(config['cpus'])
        if cpus is None:
            errors.append(f"CPUs must be an integer or None (got {config['cpus']})")
        else:
            available = os.cpu_count()
            if cpus < 1:
                errors.append(f"CPUs must be at least 1 (got {cpus})")
            elif available and cpus > available:
                # Oversubscription is suboptimal, not invalid: it is a legitimate
                # choice on a shared or hyperthreaded machine, so warn rather than
                # refuse to run.
                print(
                    f"⚠️  Warning: cpus ({cpus}) exceeds the {available} CPUs detected. "
                    f"This may slow the run down rather than speed it up."
                )

    #Validate container engine
    if 'container_engine' in config:
        engine = config['container_engine']
        valid_engines = ['docker', 'singularity', 'apptainer']
        if engine not in valid_engines:
            errors.append(f"container_engine must be one of {valid_engines} (got '{engine}')")

    # Validate structures section
    if 'structures' in config:
        if not isinstance(config['structures'], dict):
            errors.append("'structures' must be a dictionary with 'enabled' key")
        elif 'enabled' not in config['structures']:
            errors.append("'structures' section missing 'enabled' key")
        elif not isinstance(config['structures']['enabled'], bool):
            errors.append(f"structures.enabled must be boolean (got {config['structures']['enabled']})")
        else:
            if config['structures']['enabled']:
                # proteome_uniprot is accepted for backwards compatibility but no
                # longer required: nothing in the pipeline reads it. The structures
                # module queries UniProt by organism.tax_id and organism.strain_taxid.
                # Requiring it only made users invent a value with no effect.
                # Present or absent, it is accepted and ignored.


                # Optional booleans: validate when present, never fill in defaults
                # here. Config owns the defaults, so validation stays read-only.
                for bool_field in ('pocket_full_mode', 'uniparc_fallback'):
                    if bool_field in config['structures'] and not isinstance(
                        config['structures'][bool_field], bool
                    ):
                        errors.append(
                            f"structures.{bool_field} must be boolean "
                            f"(got {config['structures'][bool_field]})"
                        )
    
    #Validate colabfold section
    if 'colabfold' in config:
        if not isinstance(config['colabfold'], dict):
            errors.append("'colabfold' must be a dictionary with 'enabled' key")
        elif 'enabled' not in config['colabfold']:
            errors.append("'colabfold' section missing 'enabled' key")
        elif not isinstance(config['colabfold']['enabled'], bool):
            errors.append(f"colabfold.enabled must be boolean (got {config['colabfold']['enabled']})")
        else:
            if config['colabfold']['enabled']:
                # Optional booleans: validate when present, never fill in defaults
                # here. Config owns the defaults, so validation stays read-only.
                for bool_field in ('amber', 'gpu', 'colabfold_run_all'):
                    if bool_field in config['colabfold'] and not isinstance(
                        config['colabfold'][bool_field], bool
                    ):
                        errors.append(f"colabfold.{bool_field} must be boolean (got {config['colabfold'][bool_field]})")

    # Ligands is optional for older configurations and independent of pockets.
    if 'alphafill' in config:
        errors.append("Move top-level alphafill settings under ligands.alphafill")
    if 'ligands' in config:
        ligands_section = config['ligands']
        if not isinstance(ligands_section, dict):
            errors.append("'ligands' must be a dictionary with 'enabled' key")
        elif 'enabled' not in ligands_section:
            errors.append("'ligands' section missing 'enabled' key")
        elif not isinstance(ligands_section['enabled'], bool):
            errors.append("ligands.enabled must be boolean")
        elif 'alphafill' in ligands_section:
            section = ligands_section['alphafill']
            if not isinstance(section, dict):
                errors.append("'ligands.alphafill' must be a dictionary with 'enabled' key")
            elif 'enabled' not in section:
                errors.append("'ligands.alphafill' section missing 'enabled' key")
            elif not isinstance(section['enabled'], bool):
                errors.append("ligands.alphafill.enabled must be boolean")
            elif ligands_section['enabled'] and section['enabled']:
                from ftscripts.alphafill import validate_options
                try:
                    section.update(validate_options(section))
                except (ValueError, TypeError, OverflowError) as error:
                    errors.append(f"ligands.alphafill: {error}")

    if isinstance(config.get('ligands'), dict):
        parent = config['ligands']
        taxonomy = parent.get('taxonomy_dir')
        if taxonomy is not None and (not isinstance(taxonomy, str) or not taxonomy.strip()):
            errors.append('ligands.taxonomy_dir must be a non-empty path or null')
        from ftscripts.ligand_inputs import validate_options as validate_ligand_options
        for source in ('chembl', 'bindingdb'):
            if source not in parent:
                continue
            section = parent[source]
            if not isinstance(section, dict) or not isinstance(section.get('enabled'), bool):
                errors.append(f'ligands.{source} must contain a boolean enabled key')
            elif parent.get('enabled') and section['enabled']:
                try:
                    section.update(validate_ligand_options(source, section))
                except (ValueError, TypeError) as error:
                    errors.append(f'ligands.{source}: {error}')

    # Validate metabolism-PathwayTools section
    if 'metabolism-PathwayTools' in config and config['metabolism-PathwayTools'].get('enabled'):
        met = config['metabolism-PathwayTools']
        required_met_files = ['sbml_file', 'chokepoint_file', 'smarttable_file']
        for field in required_met_files:
            if field not in met:
                errors.append(f"metabolism-PathwayTools missing required field: '{field}'")
            else:
                filepath = met[field]
                if not os.path.exists(filepath):
                    errors.append(f"metabolism-PathwayTools {field} not found: {filepath}")
        
        # Validate optional curated_ubiquitous_file
        if 'curated_ubiquitous_file' in met and met['curated_ubiquitous_file']:
            ubiq_path = met['curated_ubiquitous_file']
            if not os.path.exists(ubiq_path):
                errors.append(f"metabolism-PathwayTools curated_ubiquitous_file not found: {ubiq_path}")
    
    # Validate metabolism-SBML section
    if 'metabolism-SBML' in config and config['metabolism-SBML'].get('enabled'):
        met_sbml = config['metabolism-SBML']
        if 'sbml_file' not in met_sbml:
            errors.append("metabolism-SBML missing required field: 'sbml_file'")
        else:
            sbml_path = met_sbml['sbml_file']
            if not os.path.exists(sbml_path):
                errors.append(f"metabolism-SBML sbml_file not found: {sbml_path}")
        
        # Filter file is optional, but if provided, validate it
        if 'filter_file' in met_sbml and met_sbml['filter_file']:
            filter_path = met_sbml['filter_file']
            if not os.path.exists(filter_path):
                errors.append(f"metabolism-SBML filter_file not found: {filter_path}")
    
    for section in ('metabolism-PathwayTools', 'metabolism-SBML'):
        settings = config.get(section, {}).get('network_metrics', {})
        if not isinstance(settings, dict):
            errors.append(f'{section}.network_metrics must be a mapping')
            continue
        repetitions, seed = settings.get('random_networks', 100), settings.get('seed', 42)
        alpha = settings.get('pagerank_alpha', 0.85)
        if isinstance(repetitions, bool) or not isinstance(repetitions, int) or repetitions < 1:
            errors.append(f'{section}.random_networks must be a positive integer')
        if isinstance(seed, bool) or not isinstance(seed, int):
            errors.append(f'{section}.seed must be an integer')
        if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) or not 0 < alpha < 1:
            errors.append(f'{section}.pagerank_alpha must be in (0,1)')

    # Validate core genome section
    if 'core' in config and config['core'].get('enabled'):
        core = config['core']
        from ftscripts.dereplication import options as dereplication_options
        try:
            dereplication_options(core)
        except ValueError as error:
            errors.append(str(error))
        
        # Check at least one core method is enabled
        if not core.get('roary') and not core.get('corecruncher') and not core.get('panaroo'):
            errors.append("Core analysis enabled but none of 'roary', 'corecruncher' or 'panaroo' is enabled")
        if core.get('panaroo_mode', 'strict') not in ('strict', 'moderate', 'sensitive'):
            errors.append(f"core.panaroo_mode must be strict, moderate or sensitive (got {core.get('panaroo_mode')})")
        
        # Validate identity and frequency (0-100%)
        if 'min_identity' in core:
            try:
                identity = float(core['min_identity'])
                if not (0 <= identity <= 100):
                    errors.append(f"core.min_identity must be between 0 and 100 (got {identity})")
            except (ValueError, TypeError):
                errors.append(f"core.min_identity must be a number (got {core['min_identity']})")
        
        if 'min_core_freq' in core:
            try:
                freq = float(core['min_core_freq'])
                if not (0 <= freq <= 100):
                    errors.append(f"core.min_core_freq must be between 0 and 100 (got {freq})")
            except (ValueError, TypeError):
                errors.append(f"core.min_core_freq must be a number (got {core['min_core_freq']})")
    
    # Validate offtarget section
    if 'offtarget' in config and config['offtarget'].get('enabled'):
        offt = config['offtarget']
        
        if offt.get('microbiome'):
            threads_per_genome = offt.get('microbiome_threads_per_genome', 4)
            if (
                isinstance(threads_per_genome, bool)
                or not isinstance(threads_per_genome, int)
                or threads_per_genome < 1
            ):
                errors.append(
                    "offtarget.microbiome_threads_per_genome must be a "
                    "positive integer"
                )

            catalogue_configs = offt.get("microbiome_catalogues")
            if not isinstance(catalogue_configs, list) or not catalogue_configs:
                errors.append("offtarget.microbiome_catalogues must be a non-empty list")
                catalogue_configs = []

            seen_catalogues = set()
            for index, catalogue in enumerate(catalogue_configs):
                prefix = f"offtarget.microbiome_catalogues[{index}]"
                if not isinstance(catalogue, dict):
                    errors.append(f"{prefix} must be a mapping")
                    continue

                name = catalogue.get('name')
                if name not in MGNIFY_CATALOGUES:
                    errors.append(
                        f"{prefix}.name must be one of "
                        f"{sorted(MGNIFY_CATALOGUES)} (got {name!r})"
                    )
                elif name in seen_catalogues:
                    errors.append(f"Duplicate microbiome catalogue: {name}")
                else:
                    seen_catalogues.add(name)

                for field in ('identity_filter', 'coverage_filter'):
                    if field not in catalogue:
                        errors.append(f"{prefix}.{field} is required")
                        continue
                    try:
                        value = float(catalogue[field])
                        if not (0 <= value <= 100):
                            errors.append(
                                f"{prefix}.{field} must be between 0 and 100 "
                                f"(got {value})"
                            )
                    except (ValueError, TypeError):
                        errors.append(f"{prefix}.{field} must be a number")
        
        # Validate foldseek_human requires structures
        if offt.get('foldseek_human') and not config.get('structures', {}).get('enabled'):
            errors.append("offtarget.foldseek_human requires structures to be enabled")
    
    # Validate DEG section
    if 'deg' in config and config['deg'].get('enabled'):
        deg = config['deg']
        
        if 'deg_identity_filter' in deg:
            try:
                identity = float(deg['deg_identity_filter'])
                if not (0 <= identity <= 100):
                    errors.append(f"deg.deg_identity_filter must be between 0 and 100 (got {identity})")
            except (ValueError, TypeError):
                errors.append(f"deg.deg_identity_filter must be a number (got {deg['deg_identity_filter']})")
        
        if 'deg_coverage_filter' in deg:
            try:
                coverage = float(deg['deg_coverage_filter'])
                if not (0 <= coverage <= 100):
                    errors.append(f"deg.deg_coverage_filter must be between 0 and 100 (got {coverage})")
            except (ValueError, TypeError):
                errors.append(f"deg.deg_coverage_filter must be a number (got {deg['deg_coverage_filter']})")
    
    # Validate PSortB section
    if 'psortb' in config and config['psortb'].get('enabled'):
        psortb = config['psortb']
        
        if 'gram_type' not in psortb:
            errors.append("psortb missing required field: 'gram_type'")
        else:
            gram = psortb['gram_type']
            valid_gram_types = ['n', 'p', 'a']
            if gram not in valid_gram_types:
                errors.append(f"psortb.gram_type must be one of {valid_gram_types} (got '{gram}')")
    
    # Validate metadata section
    if 'metadata' in config and config['metadata'].get('enabled'):
        meta = config['metadata']
        
        if 'meta_tables' not in meta:
            errors.append("metadata missing required field: 'meta_tables'")
        else:
            meta_tables = meta['meta_tables']
            if not isinstance(meta_tables, list):
                errors.append(f"metadata.meta_tables must be a list (got {type(meta_tables).__name__})")
            elif len(meta_tables) == 0:
                errors.append("metadata.meta_tables is empty but metadata is enabled")
            else:
                for i, table_path in enumerate(meta_tables):
                    if not isinstance(table_path, str):
                        errors.append(f"metadata.meta_tables[{i}] must be a string path")
                    elif not os.path.exists(table_path):
                        errors.append(f"metadata table not found: {table_path}")
                    elif not table_path.lower().endswith(('.tsv', '.csv', '.txt')):
                        errors.append(f"metadata table should be .tsv, .csv, or .txt file: {table_path}")
    
    # Raise all validation errors together
    if errors:
        error_msg = "\n❌ Configuration validation failed:\n\n" + "\n".join(f"  • {err}" for err in errors)
        raise ValueError(error_msg)

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
        print(f"Pocket Full Mode: {config.structures['pocket_full_mode']}")
        print(f"UniParc Fallback: {config.structures['uniparc_fallback']}")
    else:
        print(f"Structures Enabled: {config.structures}")

    print("\n---- ColabFold ----")
    if config.colabfold:
        print("ColabFold options:")
        print(f"  Amber Refinement: {config.colabfold['amber']}")
        print(f"  GPU Acceleration: {config.colabfold['gpu']}")
        print(f"  Run All Models: {config.colabfold['colabfold_run_all']}")
    else:
        print(f"ColabFold Enabled: {config.colabfold}")
    
    for source in ("chembl", "bindingdb"):
        print(f"{source} Enabled: {bool(getattr(config, source, False))}")
    print("\n---- Ligands ----")
    print(f"Ligands Enabled: {bool(getattr(config, 'ligands', False))}")
    if getattr(config, 'alphafill', False):
        print("AlphaFill Enabled: True")
        for key, value in config.alphafill.items():
            if key != 'enabled':
                print(f"  {key}: {value}")
    else:
        print("AlphaFill Enabled: False")

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
        print(f"Panaroo Enabled: {config.core['panaroo']} (clean mode {config.core['panaroo_mode']})")
        print(f"Minimum Identity: {config.core['min_identity']}%")
        print(f"Minimum Core Frequency: {config.core['min_core_freq']}%")
    else:
        print(f"Core Enabled: {config.core}")
      
    print("\n---- Offtarget ----")
    
    if config.offtarget:
        print(f"Human Offtarget Enabled: {config.offtarget['human']}")
        print(f"Microbiome Offtarget Enabled: {config.offtarget['microbiome']}")
        if config.offtarget['microbiome']:
            print(
                "Microbiome Threads per Genome: "
                f"{config.offtarget.get('microbiome_threads_per_genome', 4)}"
            )
            for catalogue in config.offtarget["microbiome_catalogues"]:
                print(
                    f"Microbiome Catalogue: {catalogue['name']} "
                    f"(identity={catalogue['identity_filter']}%, "
                    f"coverage={catalogue['coverage_filter']}%)"
                )
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
    if config.colabfold: enabled_modules.append("ColabFold")
    if getattr(config, "ligands", False):
        enabled_modules.append("Ligands (AlphaFill)" if config.alphafill else "Ligands")
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

    # Handle config_path: if absolute, use it; if relative, check current dir first
    if os.path.isabs(config_path):
        # Absolute path provided by user
        final_config_path = config_path
    elif os.path.exists(config_path):
        # Relative path exists in current working directory
        final_config_path = os.path.abspath(config_path)
    else:
        # Fallback: look in the script's directory (for default config.yml)
        base_path = os.path.dirname(os.path.abspath(__file__))
        final_config_path = os.path.join(base_path, config_path)
    
    config = load_config(final_config_path)
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
