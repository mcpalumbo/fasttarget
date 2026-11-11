import os
import configuration
import pandas as pd
import argparse
import multiprocessing
from ftscripts import files, structures, pathways, offtargets, genome, essentiality, metadata, alphafill
from datetime import datetime
import logging
import sys
from ftscripts.logger import logger 
import shutil

# Steps that failed during the current run.
#
# Sub-steps of a module are independent (off-target runs human, microbiome and
# Foldseek searches separately), so a failure must not abort the siblings that
# could still succeed. Instead each failure is recorded here, the run finishes
# everything it can, and the process exits non-zero at the end. Without this the
# pipeline wrote a partial table and exited 0, which is indistinguishable from
# a clean run.
_failed_steps = []


def reset_failed_steps():
    """
    Clear the recorded failures. Called at the start of each run.
    """
    _failed_steps.clear()


def failed_steps():
    """
    Names of the steps that failed during the current run.

    :return: List of step names, in the order they failed.
    """
    return list(_failed_steps)


def record_failure(step, error):
    """
    Record a failed step, log the traceback and warn on stdout.

    :param step: Human-readable name of the step that failed.
    :param error: The exception that was raised.
    """
    _failed_steps.append(step)
    if sys.exc_info()[0] is not None:
        logging.exception(f'Error in {step}: {error}')
    else:
        # Called outside an except block: there is no traceback to attach.
        logging.error(f'Error in {step}: {error}')
    print(f'ERROR: {step} failed - {error}')


def print_stylized(title, width=80):
    """
    Print a stylized title to the console.
    """
    dash_line = '-' * width
    asterisk_line = '*' * width
    print(dash_line)
    print(f'{title.center(width)}')
    print(asterisk_line)

def prepare_genome_files(config, output_path):
    """
    Prepare genome files and organism subfolders.
    :param config: Configuration object.
    :param output_path: Path to the output directory.
    """

    # Organism data
    print_stylized('GENOME')

    organism_name = config.organism['name']
    tax_id = config.organism['tax_id']
    strain_taxid = config.organism['strain_taxid']
    gbk_file = config.organism['gbk_file']
    
    print(f'Organism name: {organism_name}')
    print(f'Species Tax ID: {tax_id}')
    print(f'Strain Tax ID: {strain_taxid}')
    print(f'Genome file: {gbk_file}')

    # Create organism subfolders
    files.create_organism_subfolders(output_path, organism_name)
    logging.info(f'Organism subfolders created in {output_path}/{organism_name}')

    # Create organism genome files (gbk, gff3 and fasta)
    files_dir = os.path.join(output_path, organism_name, 'genome')
    if not os.path.exists(files_dir):
        os.makedirs(files_dir, exist_ok=True)
        logging.info(f'Genome directory created in {files_dir}')
        
    genome.ref_genome_files(gbk_file, files_dir, organism_name, container_engine=config.container_engine)

    logging.info(f'Genome files created in {output_path}/{organism_name}')

def metabolic_module(config, output_path):
    """
    Run metabolic analysis based on configuration.
    :param config: Configuration object.
    :param output_path: Path to the output directory.
    :return: List of resulting DataFrames.
    """

    organism_name = config.organism['name']
    module_tables = []

    # Run METABOLIC ANALYSIS
    if config.metabolism_pathwaytools:
        try:
            print_stylized('Metabolic analysis - Pathway Tools')
            print('Using Pathway Tools output files')

            logging.info('Starting metabolic analysis')

            sbml_file =  config.metabolism_pathwaytools['sbml_file']
            chokepoint_file = config.metabolism_pathwaytools['chokepoint_file']
            smarttable_file = config.metabolism_pathwaytools['smarttable_file']
            curated_ubiquitous_file = config.metabolism_pathwaytools.get('curated_ubiquitous_file', None)
            
            # Handle empty string or None for curated_ubiquitous_file
            if curated_ubiquitous_file == "" or curated_ubiquitous_file is None:
                curated_ubiquitous_file = None

            logging.info(f'SBML file: {sbml_file}')
            logging.info(f'Chokepoint file: {chokepoint_file}')
            logging.info(f'Smarttable file: {smarttable_file}')
            if curated_ubiquitous_file:
                logging.info(f'Curated ubiquitous file: {curated_ubiquitous_file}')
            else:
                logging.info('Curated ubiquitous file: Not provided (will use auto-generated)')

            # Parse metabolic files, make network and calculate centrality
            df_centrality, df_edges, producing_df, consuming_df, both_df = pathways.run_metabolism_ptools (output_path, organism_name, sbml_file, chokepoint_file, smarttable_file, curated_ubiquitous_file)
            module_tables.append(df_centrality)
            module_tables.append(df_edges)
            module_tables.append(producing_df)
            module_tables.append(consuming_df)
            module_tables.append(both_df)

            logging.info('Metabolic analysis finished')

        except Exception as e:
            record_failure('metabolic analysis (Pathway Tools)', e)
    else:
        logging.info('Metabolic analysis with Pathway Tools files not enabled')
    

    if config.metabolism_sbml:
        try:
            print_stylized('Metabolic analysis - Custom SBML file')
            print('Using SBML file and MetaGraphTools')
            
            logging.info('Starting metabolic analysis')

            sbml_file =  config.metabolism_sbml['sbml_file']
            filter_file = config.metabolism_sbml.get('filter_file', None)
            
            # Handle empty string or None for filter_file
            if filter_file == "" or filter_file is None:
                filter_file = None

            logging.info(f'SBML file: {sbml_file}')
            if filter_file:
                logging.info(f'Filter file: {filter_file}')
            else:
                logging.info('Filter file: Not provided (will use default frequency filter)')

            # Parse metabolic files, make network and calculate centrality
            mgt_bc_df, mgt_degree_df, mgt_consumption_df, mgt_production_df = pathways.run_metabolism_sbml (output_path, organism_name, sbml_file, filter_file, container_engine=config.container_engine)

            module_tables.append(mgt_bc_df)
            module_tables.append(mgt_degree_df)
            module_tables.append(mgt_consumption_df)
            module_tables.append(mgt_production_df)

            logging.info('Metabolic analysis from SBML with MetaGraphTools finished')

        except Exception as e:
            record_failure('metabolic analysis (SBML)', e)
    else:
        logging.info('Metabolic analysis with SBML file not enabled')
    
    return module_tables
    
def structure_module(config, databases_path, output_path, cpus):
    """
    Run structure analysis based on configuration.
    :param config: Configuration object.
    :param databases_path: Path to the databases directory.
    :param output_path: Path to the output directory.
    :param cpus: Number of CPUs to use.
    :return: List of resulting DataFrames.
    """

    module_tables = []
    
    # Run STRUCTURES
    organism_name = config.organism['name']
    tax_id = config.organism['tax_id']
    strain_taxid = config.organism['strain_taxid']

    if config.structures:
        try:
            print_stylized('STRUCTURES')
            logging.info('Starting structures analysis')

            logging.info(f'Species Tax ID: {tax_id}')
            logging.info(f'Strain Tax ID: {strain_taxid}')

            pocket_full_mode = config.structures.get('pocket_full_mode', False)
            logging.info(f'Pocket full mode: {pocket_full_mode}')

            #Colabfold options
            amber_option = False
            gpu_option = False
            colabfold_all_models_option = False
            if config.colabfold:
                amber_option = config.colabfold.get('amber', False)
                gpu_option = config.colabfold.get('gpu', False)
                colabfold_all_models_option = config.colabfold.get('colabfold_run_all', False)
            logging.info(f'ColabFold Amber option: {amber_option}')
            logging.info(f'ColabFold GPU option: {gpu_option}')
            logging.info(f'ColabFold Run All Models option: {colabfold_all_models_option}')

            uniparc_fallback_option = bool(config.structures.get('uniparc_fallback', False)) if config.structures else False
            logging.info(f'UniParc exact-identity fallback option: {uniparc_fallback_option}')


            # Run complete structure pipeline: UniProt mapping, structure download, and pocket detection
            df_structures = structures.pipeline_structures(output_path, 
                                                            organism_name, 
                                                            tax_id, 
                                                            strain_taxid, 
                                                            cpus=cpus,
                                                            resolution_cutoff=3.5, 
                                                            coverage_cutoff=40.0,
                                                            container_engine=config.container_engine, 
                                                            full_mode=pocket_full_mode, 
                                                            amber_option=amber_option, 
                                                            gpu_option=gpu_option,
                                                            colabfold=config.colabfold,
                                                            colabfold_all_models=colabfold_all_models_option,
                                                            uniparc_fallback=uniparc_fallback_option)
            logging.info('Structures analysis finished')
            module_tables.append(df_structures)

        except Exception as e:
            record_failure('structures analysis', e)
    else:
        logging.info('Structures analysis not enabled')
    
    return module_tables

def prepare_alphafill_mapping(config, output_path, cpus):
    """Reuse the structure mapping or build only the UniProt mapping stages."""
    organism_name = config.organism['name']
    tax_id = config.organism['tax_id']
    strain_taxid = config.organism['strain_taxid']
    loci = list(dict.fromkeys(metadata.ref_gbk_locus(output_path, organism_name)))
    if not loci:
        raise ValueError('AlphaFill requires genome CDS features with locus_tag and translation')
    uniprot_dir = os.path.join(output_path, organism_name, 'structures', 'uniprot_files')
    mapping_path = os.path.join(uniprot_dir, f'uniprot_{organism_name}_id_mapping.json')
    if files.file_check(mapping_path):
        logging.info('AlphaFill reusing UniProt mapping: %s', mapping_path)
        mapping = files.json_to_dict(mapping_path)
    else:
        logging.info('Preparing UniProt mapping for AlphaFill without pocket analysis')
        structures.download_species_uniprot_data(output_path, organism_name, tax_id)
        species_file = os.path.join(uniprot_dir, f'uniprot_specie_taxid_{tax_id}_data.tsv')
        structures.parse_uniprot_species_data(species_file, tax_id, strain_taxid)
        structures.cluster_uniprot_specie(output_path, organism_name, tax_id)
        structures.create_uniprot_blast_db(output_path, organism_name, tax_id, strain_taxid)
        structures.uniprot_proteome_blast(output_path, organism_name, tax_id, strain_taxid, cpus=cpus)
        uniparc_fallback_option = bool(config.structures.get('uniparc_fallback', False)) if config.structures else False
        mapping = structures.uniprot_proteome_mapping(output_path, organism_name, tax_id, strain_taxid, uniparc_fallback=uniparc_fallback_option)
    if not isinstance(mapping, dict):
        raise ValueError(f'Invalid UniProt mapping for AlphaFill: {mapping_path}')
    genome_mapping = {gene: mapping.get(gene) for gene in loci}
    for gene, values in genome_mapping.items():
        if values is None:
            continue
        values = values if isinstance(values, list) else [values]
        if any(value is not None and not isinstance(value, str) for value in values):
            raise ValueError(f'Invalid UniProt accession value in mapping for {gene}')
    return genome_mapping


def alphafill_module(config, output_path, cpus):
    """Return only the one-row-per-locus summary for the final FastTarget table."""
    options = getattr(config, 'alphafill', False)
    if not options:
        logging.info('AlphaFill analysis not enabled')
        return []
    logging.info('Starting AlphaFill ligand annotation')
    mapping = prepare_alphafill_mapping(config, output_path, cpus)
    summary = alphafill.run_alphafill(output_path, config.organism['name'], mapping, options)
    if summary['gene'].duplicated().any():
        raise ValueError('AlphaFill summary must have exactly one row per locus_tag')
    logging.info('AlphaFill finished for %s loci; details in ligands/alphafill', len(summary))
    return [summary]


def ligands_module(config, output_path, cpus, databases_path=None):
    """Run independent ligand sources and retain one summary row per locus."""
    if not getattr(config, 'ligands', False):
        return []
    from pathlib import Path
    import json
    from ftscripts import chembl, bindingdb
    from ftscripts.ligand_inputs import resolve_database, local_path, validate_options
    print_stylized('LIGANDS')
    try:
        tables = alphafill_module(config, output_path, cpus)
    except Exception:
        if getattr(config, 'alphafill', False) and config.alphafill.get('strict_errors', False):
            raise
        record_failure('AlphaFill ligand analysis', error)
        tables = []
    databases_path = Path(databases_path or Path(__file__).resolve().parent / 'databases')
    organism = config.organism['name']
    base = Path(output_path) / organism
    for source, runner in (('chembl', chembl.run_chembl), ('bindingdb', bindingdb.run_bindingdb)):
        options = getattr(config, source, False)
        if not options:
            continue
        destination = base / 'ligands' / source
        try:
            settings = validate_options(source, options)
            database = resolve_database(source, settings.pop('database'), databases_path)
            settings.pop('enabled', None)
            settings.pop('strict_errors')
            taxonomy = local_path(config.ligands.get('taxonomy_dir') or 'ncbi_taxonomy', databases_path)
            logging.info('Starting %s with local database %s', source, database)
            result = runner(database, config.organism['tax_id'], destination,
                            proteins_fasta=base / 'genome' / f'{organism}.faa',
                            taxonomy_dir=taxonomy, threads=cpus, **settings)
            summary = result['summary']
            genes = list(dict.fromkeys(metadata.ref_gbk_locus(output_path, organism)))
            expected = ['gene', f'{source}_status', f'{source}_target_count',
                        f'{source}_compounds', f'{source}_compound_count']
            if list(summary.columns) != expected or summary['gene'].duplicated().any() or set(summary['gene']) != set(genes):
                raise ValueError(f'{source} summary must contain exactly one row per reference locus and the five summary columns')
            destination.mkdir(parents=True, exist_ok=True)
            (destination / 'pipeline_status.json').write_text(json.dumps({'status': 'completed'}) + '\n')
            tables.append(summary)
        except Exception as error:
            record_failure(f'{source} ligand analysis', error)
            destination.mkdir(parents=True, exist_ok=True)
            (destination / 'pipeline_status.json').write_text(json.dumps({'status': 'failed', 'error': str(error)}) + '\n')
            genes = list(dict.fromkeys(metadata.ref_gbk_locus(output_path, organism)))
            summary = pd.DataFrame({'gene': genes, f'{source}_status': 'error',
                                    f'{source}_target_count': pd.NA, f'{source}_compounds': pd.NA,
                                    f'{source}_compound_count': pd.NA})
            for extension, separator in (('tsv', '\t'), ('csv', ',')):
                summary.to_csv(destination / f'{source}_summary.{extension}', sep=separator, index=False)
            if options.get('strict_errors', False):
                raise
            tables.append(summary)
    return tables


def conservation_module(config, databases_path, output_path, cpus):
    """
    Run conservation analysis based on configuration.
    :param config: Configuration object.
    :param databases_path: Path to the databases directory.
    :param output_path: Path to the output directory.
    :param cpus: Number of CPUs to use.
    :return: List of resulting DataFrames.
    """

    module_tables = []
    
    # Run CORE ANALYSIS
    organism_name = config.organism['name']
    tax_id = config.organism['tax_id']

    if config.core:
        try:
            print_stylized('CORE ANALYSIS')
            
            logging.info('Starting core analysis')

            # Read accession list if provided
            accession_list = None
            if config.core.get('accession_file'):
                accession_file = config.core['accession_file']
                if accession_file and accession_file != 'null' and os.path.exists(accession_file):
                    print(f'Reading accession list from: {accession_file}')
                    with open(accession_file, 'r') as f:
                        accession_list = [line.strip() for line in f if line.strip() and not line.startswith('#')]
                    print(f'Found {len(accession_list)} accessions to download')
                    logging.info(f'Using accession list with {len(accession_list)} genomes')
                else:
                    logging.warning(f'Accession file not found or invalid: {accession_file}')

            # Download complete NCBI genomes from organism tax id or accession list
            if accession_list:
                print('----- 1. Downloading specific genomes from accession list -----')
            else:
                print('----- 1. Downloading tax_id genomes from NCBI -----')
            genome.core_download_genomes_ncbi(output_path, organism_name, tax_id, accession_list=accession_list)
            genome.core_download_missing_accessions(output_path, organism_name, tax_id, accession_list=accession_list)
            logging.info('Genomes downloaded')

            # Keep genomes with human as host. Check presence of .gff and .faa files for each strain
            print('----- 1. Selecting genomes -----')
            # If accession_list was provided, don't filter by host; otherwise filter for human host
            filter_by_host = accession_list is None or len(accession_list) == 0
            genome.core_check_files(output_path, organism_name, container_engine=config.container_engine, filter_by_host=filter_by_host)
            logging.info('Genomes filtered')
            print('----- 1. Finished -----')

            min_identity = config.core['min_identity']
            min_core_freq = config.core['min_core_freq']

            if config.core['roary']:
                try:
                    #Run roary
                    print('----- 2. Running Roary -----')
                    logging.info('Starting Roary analysis')
                    genome.core_genome_programs(output_path, organism_name, min_core_freq, min_identity, cpus, program_list=['roary'], container_engine=config.container_engine)
                    # Parse output
                    print('----- 2. Parsing Roary results -----')
                    df_roary = genome.roary_output(output_path, organism_name, core_threshold=min_core_freq/100)
                    module_tables.append(df_roary)
                    logging.info('Roary analysis finished')
                    print('----- 2. Finished -----')
                except Exception as e:
                    record_failure('Roary analysis', e)
            else:
                logging.info('Roary not enabled')

            if config.core['corecruncher']:
                try:
                    # Run CoreCruncher
                    print('----- 2. Running CoreCruncher -----')
                    logging.info('Starting CoreCruncher analysis')
                    genome.core_genome_programs(output_path, organism_name, min_core_freq, min_identity, cpus, program_list=['corecruncher'], container_engine=config.container_engine)
                    # Parse output
                    print('----- 2. Parsing CoreCruncher results -----')
                    df_cc = genome.corecruncher_output(output_path, organism_name, core_threshold=min_core_freq)
                    module_tables.append(df_cc)
                    logging.info('CoreCruncher analysis finished')
                    print('----- 2. Finished -----')
                except Exception as e:
                    record_failure('CoreCruncher analysis', e)
            else:
                logging.info('CoreCruncher not enabled')
        except Exception as e:
            record_failure('core analysis', e)
    else:
        logging.info('Core analysis not enabled')
    
    return module_tables

def offtarget_module(config, databases_path, output_path, cpus):
    """
    Run offtarget analysis based on configuration.
    :param config: Configuration object.
    :param databases_path: Path to the databases directory.
    :param output_path: Path to the output directory.
    :param cpus: Number of CPUs to use.
    :return: List of resulting DataFrames.
    """

    module_tables = []
    organism_name = config.organism['name']
    
    # Run OFFTARGETS
    if config.offtarget:
        try:
            offtarget_path = os.path.join(output_path, organism_name, 'offtarget')

            if config.offtarget['human']:
                try:
                    print_stylized('HUMAN OFFTARGET')

                    human_blast_output = os.path.join(offtarget_path, 'human_offtarget_blast.tsv')
                    
                    if not files.step_completed(human_blast_output):
                        # Run blastp search
                        print('-----  Blastp search -----')
                        offtargets.human_offtarget_blast(databases_path, output_path, organism_name, cpus)
                        logging.info('Human offtarget blast search finished')
                    else:
                        logging.info('Blast with human already done')
                        print('Blast output file found')
                        print(human_blast_output)

                    # Parse results
                    df_human = offtargets.human_offtarget_parse(output_path, organism_name)
                    module_tables.append(df_human)
                    logging.info('Human offtarget analysis finished')
                    print('----- Finished -----')
                except Exception as e:
                    record_failure('human offtarget analysis', e)
            else:
                logging.info('Human offtarget analysis not enabled')

            if config.offtarget['microbiome']:
                try:

                    print_stylized('MICROBIOME OFFTARGET')

                    threads_per_genome = config.offtarget.get(
                        'microbiome_threads_per_genome',
                        4,
                    )
                    for catalogue in config.offtarget["microbiome_catalogues"]:
                        catalogue_name = catalogue['name']
                        identity_filter = float(catalogue['identity_filter'])
                        coverage_filter = float(catalogue['coverage_filter'])
                        print(f"----- {catalogue_name}: DIAMOND search -----")
                        offtargets.microbiome_offtarget_blast_species(
                            databases_path,
                            output_path,
                            organism_name,
                            catalogue_name,
                            identity_filter,
                            coverage_filter,
                            cpus,
                            threads_per_genome,
                        )
                        logging.info(
                            "%s microbiome search finished with identity=%s and coverage=%s",
                            catalogue_name,
                            identity_filter,
                            coverage_filter,
                        )
                        result_tables = offtargets.microbiome_species_parse(
                            databases_path,
                            output_path,
                            organism_name,
                            catalogue_name,
                            identity_filter,
                            coverage_filter,
                            export_taxonomy_details=True,
                        )
                        module_tables.extend(result_tables)
                        if catalogue_name == "human-gut":
                            try:
                                module_tables.append(offtargets.microbiome_silva_core_parse(
                                    databases_path, output_path, organism_name,
                                ))
                            except Exception as error:
                                record_failure('Silva core annotation', error)
                            try:
                                module_tables.append(offtargets.microbiome_silva_prevalence_parse(
                                    databases_path, output_path, organism_name,
                                ))
                            except Exception as error:
                                record_failure('Silva prevalence annotation', error)
                    logging.info('Microbiome offtarget analysis finished')
                    print('----- Finished -----')
                            
                except Exception as e:
                    record_failure('microbiome offtarget analysis', e)
            else:
                logging.info('Microbiome offtarget analysis not enabled')
            
            available = {}

            if config.structures and config.offtarget['foldseek_human']:
                # config.structures only says structures were requested, not that
                # they were produced. Searching an empty or partial structure set
                # yields output that looks successful while silently covering
                # fewer genes, so require real reference structures first.
                reference_structures = structures.get_all_reference_structures(
                    output_path, organism_name, path_mode=True
                )
                available = {
                    locus_tag: path for locus_tag, path in (reference_structures or {}).items() if path
                }
                total_loci = len(reference_structures or {})

                if not available:
                    skipped = (
                        'Foldseek human offtarget skipped: the structures module left no '
                        'reference structures to search.'
                    )
                    record_failure(
                        'Foldseek human offtarget analysis (no reference structures)',
                        RuntimeError(skipped),
                    )
                elif len(available) < total_loci:
                    # Partial coverage is legitimate (a gene may have no structure at
                    # all), but it must be visible rather than inferred from counts.
                    logging.warning(
                        f'Foldseek will cover {len(available)} of {total_loci} genes; '
                        f'{total_loci - len(available)} have no reference structure.'
                    )
                    print(
                        f'WARNING: Foldseek will cover {len(available)} of {total_loci} genes; '
                        f'{total_loci - len(available)} have no reference structure.'
                    )

            if (config.structures and config.offtarget['foldseek_human']
                    and available):
                try:
                    print_stylized('FOLDSEEK HUMAN OFFTARGET')
                   
                    # Run foldseek against human structures
                    foldseek_mapping = offtargets.run_foldseek_human_structures (databases_path, output_path, organism_name, container_engine=config.container_engine)
                    logging.info('Foldseek human offtarget search finished')
                    
                    # Parse results
                    results_foldseek_dict = offtargets.foldseek_human_parser (output_path, organism_name, foldseek_mapping)
                    mapped_dict_foldseek = offtargets.merge_foldseek_data (output_path, organism_name)
                    final_foldseek_df = offtargets.final_foldseek_structure_table (output_path, organism_name, mapped_dict_foldseek)
                    module_tables.append(final_foldseek_df)

                    if config.colabfold and config.colabfold.get('colabfold_run_all', False):
                        foldseek_colab_mapping = offtargets.run_foldseek_human_colabfold_structures(
                            databases_path,
                            output_path,
                            organism_name,
                            container_engine=config.container_engine
                        )
                        logging.info('Foldseek human offtarget search for ColabFold structures finished')

                        results_foldseek_colab_dict = offtargets.foldseek_human_colabfold_parser(
                            output_path,
                            organism_name,
                            foldseek_colab_mapping
                        )
                        mapped_dict_foldseek_colab = offtargets.merge_foldseek_colabfold_data(
                            output_path,
                            organism_name
                        )
                        final_foldseek_colab_df = offtargets.final_foldseek_colabfold_structure_table(
                            output_path,
                            organism_name,
                            mapped_dict_foldseek_colab
                        )
                        module_tables.append(final_foldseek_colab_df)

                    logging.info('Foldseek human offtarget analysis finished')
                    print('----- Finished -----')
                except Exception as e:
                    record_failure('Foldseek human offtarget analysis', e)
            else:
                logging.info('Foldseek human offtarget analysis not enabled')           

        except Exception as e:
            record_failure('offtarget analysis', e)
    else:
        logging.info('Offtarget analysis not enabled')
    
    return module_tables

def essentiality_module(config, databases_path, output_path, cpus):
    """
    Run essentiality analysis based on configuration.
    :param config: Configuration object.
    :param databases_path: Path to the databases directory.
    :param output_path: Path to the output directory.
    :param cpus: Number of CPUs to use.
    :return: List of resulting DataFrames.
    """

    module_tables = []
    organism_name = config.organism['name']

    # Run ESSENTIALITY
    if config.deg:
        essentiality_path = os.path.join(output_path, organism_name, 'essentiality')
        try:
            print_stylized('ESSENTIALITY')

            deg_blast_output = os.path.join(essentiality_path, 'deg_blast.tsv')
            
            if not files.step_completed(deg_blast_output):
                # Run blastp search
                print('----- Blastp search -----')
                essentiality.essential_deg_blast(databases_path, output_path, organism_name, cpus)
                logging.info('DEG blast search finished')
            else:
                logging.info('Blast with DEG already done')
                print('Blast output file found')
                print(deg_blast_output)      

            # Parse results
            deg_identity_filter = config.deg['deg_identity_filter']
            deg_coverage_filter = config.deg['deg_coverage_filter']
            logging.info(f'DEG identity filter: {deg_identity_filter}')
            logging.info(f'DEG coverage filter: {deg_coverage_filter}')
            df_deg = essentiality.deg_parse(output_path, organism_name, deg_identity_filter, deg_coverage_filter)
            module_tables.append(df_deg)
            logging.info('DEG analysis finished')
            print('----- Finished -----')

        except Exception as e:
            record_failure('essentiality analysis', e)
    else:
        logging.info('Essentiality analysis not enabled')
    
    return module_tables

def localization_module(config, output_path):
    """
    Run localization analysis based on configuration.
    :param config: Configuration object.
    :param output_path: Path to the output directory.
    :return: List of resulting DataFrames.
    """

    module_tables = []
    organism_name = config.organism['name']
    # Run LOCALIZATION
    if config.psortb:
        try:
            print_stylized('LOCALIZATION')

            gram_type = config.psortb['gram_type']
            logging.info(f'Gram type: {gram_type}')
            
            #Run psortb
            print('----- Running psort -----')
            df_psort = genome.localization_prediction(output_path, organism_name, gram_type, container_engine=config.container_engine)
            module_tables.append(df_psort)
            logging.info('Psortb analysis finished')
            print('----- Finished -----')
        except Exception as e:
            record_failure('localization analysis', e)
    else:
        logging.info('Localization analysis not enabled')
    
    return module_tables

def metadata_loading(config, output_path):
    """
    Load metadata tables based on configuration.
    :param config: Configuration object.
    :param output_path: Path to the output directory.
    :return: List of resulting DataFrames.
    """

    module_tables = []
    organism_name = config.organism['name']

    # Load METADATA
    if config.metadata:
        try:
            print_stylized('METADATA')
            for table in config.metadata['meta_tables']:
                print(f'----- Loading metadata table: {table} -----')
                shutil.copy(table, os.path.join(output_path, organism_name, 'metadata'))
                with open(table, 'r') as file:
                    first_line = file.readline()
                    if '\t' in first_line:
                        sep = '\t'
                    elif ',' in first_line:
                        sep = ','
                    elif ';' in first_line:
                        sep = ';'
                    else:
                        raise ValueError('Invalid file format. Only CSV and TSV metadata files are supported.')

                df_meta = pd.read_csv(table, header=0, sep=sep)
                module_tables.append(df_meta)
                logging.info(f'Metadata table {table} loaded')
                print('----- Finished -----')
        except Exception as e:
            record_failure('metadata analysis', e)
    else:
        logging.info('Metadata analysis not enabled')
    
    return module_tables

def merge_final_tables(config, output_path, tables):
    """
    Merge final result tables and save to output.
    :param config: Configuration object.
    :param output_path: Path to the output directory.
    :param tables: List of DataFrames to merge.
    :return: Combined DataFrame with results.
    """

    organism_name = config.organism['name']

    tables = [table for table in tables if table is not None]

    # Anchor every merge on the reference gene list from the annotation.
    #
    # This used to happen only when a ligand source was configured. Without the
    # anchor the gene set of the output is whichever module's table happens to be
    # first in the list, and because every join below is a left join, any gene
    # that module dropped disappears from every other module's columns too. That
    # made the row set of the results table depend on unrelated config and on
    # module ordering. The annotated CDS set is the only correct anchor: it is the
    # user's own gene list, so the output always has exactly one row per gene.
    reference_genes = list(dict.fromkeys(metadata.ref_gbk_locus(output_path, organism_name)))

    if reference_genes:
        for table in tables:
            if 'gene' not in getattr(table, 'columns', []):
                continue
            unknown = set(table['gene']) - set(reference_genes)
            if unknown:
                # Not fatal, but it means a module emitted genes the annotation
                # does not contain; those rows are dropped by the join below.
                logging.warning(
                    f'{len(unknown)} gene(s) in a module table are absent from the '
                    f'reference annotation and will be dropped, e.g. '
                    f'{sorted(unknown)[:5]}'
                )
        tables = [pd.DataFrame({'gene': reference_genes}), *tables]
    else:
        logging.error(
            'No reference genes found in the annotation; falling back to merging on '
            'the first module table, whose gene set may be incomplete.'
        )

    # Merge dfs
    print_stylized('RESULTS')
    current_date = datetime.now().strftime('%Y-%m-%d-%H-%M')
    results_path = os.path.join(output_path, organism_name, f'{organism_name}_results_{current_date}')
    results_table_path = os.path.join(results_path, f'{organism_name}_results_table.tsv')

    if not os.path.exists(results_path):
        os.makedirs(results_path, exist_ok=True)
        logging.info(f'Results directory created in {results_path}')
    

    # The anchor is not a result: a run where every module failed must not look
    # like a successful run that produced a gene list.
    module_table_count = len(tables) - (1 if reference_genes else 0)

    if module_table_count == 0:
        logging.error('----- Error: No module produced a results table. -----')
        print('ERROR: no module produced a results table; nothing to merge.')
        results = None
    elif len(tables) > 1:
        combined_df = tables[0]
        for df in tables[1:]:
            if df is not None:
                combined_df = pd.merge(combined_df, df, on='gene', how='left')
        
        # Add gene_name and product information
        combined_df = metadata.add_gene_product_info(combined_df, output_path, organism_name)
        
        combined_df.to_csv(results_table_path, sep='\t', index=False)
        
        print(f'Final FastTarget results saved in {results_table_path}.')
        logging.info(f'Final FastTarget results saved.')
        
        results = combined_df

        # Create metadata tables for Target Pathogen
        print('\n')
        print('Creating metadata tables for Target Pathogen')
        metadata.tables_for_TP(organism_name, results_path)
        logging.info('Tables for Target Pathogen created')

    elif len(tables) == 1:
        tables[0] = metadata.add_gene_product_info(tables[0], output_path, organism_name)
        tables[0].to_csv(results_table_path, sep='\t', index=False)
        print(f'Final FastTarget results saved in {results_table_path}.')
        logging.info(f'Final FastTarget results saved.')
        results = tables[0]

        # Create metadata tables for Target Pathogen
        print('\n')
        print('Creating metadata tables for Target Pathogen')
        metadata.tables_for_TP(organism_name, results_path)
        logging.info('Tables for Target Pathogen created')
    else:
        logging.error('----- Error: No final DataFrame data. -----')
        results = None
    
    return results

def main(config, databases_path, output_path):
    """
    Main function to run FastTarget.
    
    :param config: Configuration object.
    :param databases_path: Path to the databases directory.
    :param output_path: Path to the output directory.

    :return: DataFrame with the results.
    """
    results = None

    # Prepare genome files
    prepare_genome_files(config, output_path)
    
    # Number of CPUS
    if isinstance(config.cpus, int):
        cpus = config.cpus
    else:
        cpus = multiprocessing.cpu_count()

    logging.info(f'CPUS: {cpus}')

    tables = []

    reset_failed_steps()

    # Run metabolic module
    try:
        logging.info('Starting metabolic module')
        metabolic_tables = metabolic_module(config, output_path)
        if metabolic_tables:
            tables.extend(metabolic_tables)
            logging.info(f'Metabolic module completed - {len(metabolic_tables)} table(s) generated')
        else:
            logging.info('Metabolic module completed - no tables generated')
    except Exception as e:
        record_failure('metabolic module', e)

    # Run structure module
    try:
        logging.info('Starting structure module')
        structure_tables = structure_module(config, databases_path, output_path, cpus)
        if structure_tables:
            tables.extend(structure_tables)
            logging.info(f'Structure module completed - {len(structure_tables)} table(s) generated')
        else:
            logging.info('Structure module completed - no tables generated')
    except Exception as e:
        record_failure('structure module', e)

    # Ligands can run on its own, or reuse the mapping made by structures.
    try:
        ligand_tables = ligands_module(config, output_path, cpus, databases_path=databases_path)
        tables.extend(ligand_tables)
    except Exception as e:
        if any(getattr(config, source, False) and getattr(config, source).get('strict_errors', False)
               for source in ('alphafill', 'chembl', 'bindingdb')):
            raise
        record_failure('ligands module', e)

    # Run conservation module
    try:
        logging.info('Starting conservation module')
        conservation_tables = conservation_module(config, databases_path, output_path, cpus)
        if conservation_tables:
            tables.extend(conservation_tables)
            logging.info(f'Conservation module completed - {len(conservation_tables)} table(s) generated')
        else:
            logging.info('Conservation module completed - no tables generated')
    except Exception as e:
        record_failure('conservation module', e)

    # Run offtarget module
    try:
        logging.info('Starting offtarget module')
        offtarget_tables = offtarget_module(config, databases_path, output_path, cpus)
        if offtarget_tables:
            tables.extend(offtarget_tables)
            logging.info(f'Offtarget module completed - {len(offtarget_tables)} table(s) generated')
        else:
            logging.info('Offtarget module completed - no tables generated')
    except Exception as e:
        record_failure('offtarget module', e)

    # Run essentiality module
    try:
        logging.info('Starting essentiality module')
        essentiality_tables = essentiality_module(config, databases_path, output_path, cpus)
        if essentiality_tables:
            tables.extend(essentiality_tables)
            logging.info(f'Essentiality module completed - {len(essentiality_tables)} table(s) generated')
        else:
            logging.info('Essentiality module completed - no tables generated')
    except Exception as e:
        record_failure('essentiality module', e)

    # Run localization module
    try:
        logging.info('Starting localization module')
        localization_tables = localization_module(config, output_path)
        if localization_tables:
            tables.extend(localization_tables)
            logging.info(f'Localization module completed - {len(localization_tables)} table(s) generated')
        else:
            logging.info('Localization module completed - no tables generated')
    except Exception as e:
        record_failure('localization module', e)

    # Load metadata
    try:
        logging.info('Starting metadata loading')
        metadata_tables = metadata_loading(config, output_path)
        if metadata_tables:
            tables.extend(metadata_tables)
            logging.info(f'Metadata loading completed - {len(metadata_tables)} table(s) loaded')
        else:
            logging.info('Metadata loading completed - no tables loaded')
    except Exception as e:
        record_failure('metadata loading', e)

    # Merge final tables
    logging.info(f'Total tables collected: {len(tables)}')

    results = merge_final_tables(config, output_path, tables)


    errors = failed_steps()
    if errors:
        # The tables above are still written, but the run is not a success: say so
        # loudly here and let the caller exit non-zero.
        logging.error(f'Steps with errors ({len(errors)}): {", ".join(errors)}')
        print('------------------------------------- FINISHED WITH ERRORS ----------------------')
        for step in errors:
            print(f'  FAILED: {step}')
        print(f'Results are incomplete: {len(errors)} step(s) failed. See the log for tracebacks.')
    else:
        print('------------------------------------- FINISHED ----------------------------------------')

    return results

if __name__ == "__main__":
    
    base_path = os.path.dirname(os.path.abspath(__file__))

    databases_default_path = os.path.join(base_path, 'databases')
    output_default_path = os.path.join(base_path, 'organism')

    parser = argparse.ArgumentParser(description='FastTarget script')
    parser.add_argument('--config_file', type=str, default='config.yml', help='Path to the configuration file')
    parser.add_argument('--databases_path', type=str, default=databases_default_path, help='Path to the databases directory')
    parser.add_argument('--output_path', type=str, default=output_default_path, help='Path to the output directory')
    
    args = parser.parse_args()

    config = configuration.get_config(args.config_file)


    main(config, args.databases_path, args.output_path)

    # Exit non-zero when any step failed, so wrappers and CI can tell a partial
    # run from a clean one. Outputs produced by the steps that did succeed are
    # already written at this point.
    if failed_steps():
        sys.exit(1)
