import os
import configuration
import pandas as pd
import argparse
import multiprocessing
from ftscripts import files, structures, pathways, offtargets, genome, essentiality, metadata
from datetime import datetime
import logging
import sys
from ftscripts.logger import logger 
import shutil

def print_stylized(title, width=80):
    """
    Print a stylized title to the console.
    """
    dash_line = '-' * width
    asterisk_line = '*' * width
    print(dash_line)
    print(f'{title.center(width)}')
    print(asterisk_line)

def main(config, base_path):
    """
    Main function to run FastTarget.
    
    :param config: Configuration object.
    :param base_path: Base path of the FastTarget repository folder.

    :return: DataFrame with the results.
    """
    results = None

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
    files.create_organism_subfolders(base_path, organism_name)
    logging.info(f'Organism subfolders created in {base_path}/organism/{organism_name}')

    # Create organism genome files (gbk, gff3 and fasta)
    genome.ref_genome_files(gbk_file, base_path, organism_name)
    logging.info(f'Genome files created in {base_path}/organism/{organism_name}')

    # Number of CPUS
    if isinstance(config.cpus, int):
        cpus = config.cpus
    else:
        cpus = multiprocessing.cpu_count()

    logging.info(f'CPUS: {cpus}')

    tables = []

    # Run METABOLIC ANALYSIS
    if config.metabolism:
        try:
            print_stylized('METABOLIC ANALYSIS')

            logging.info('Starting metabolic analysis')

            sbml_file =  config.metabolism['sbml_file']
            chokepoint_file = config.metabolism['chokepoint_file']
            smarttable_file = config.metabolism['smarttable_file']

            logging.info(f'SBML file: {sbml_file}')
            logging.info(f'Chokepoint file: {chokepoint_file}')
            logging.info(f'Smarttable file: {smarttable_file}')

            # Parse metabolic files, make network and calculate centrality
            df_centrality, df_edges, producing_df, consuming_df, both_df = pathways.run_metabolism (base_path, organism_name, sbml_file, chokepoint_file, smarttable_file)
            tables.append(df_centrality)
            tables.append(df_edges)
            tables.append(producing_df)
            tables.append(consuming_df)
            tables.append(both_df)

            logging.info('Metabolic analysis finished')

        except Exception as e:
            logging.error(f'Error in metabolic analysis: {e}')        
    else:
        logging.info('Metabolic analysis not enabled')

    # Run STRUCTURES
    if config.structures:
        try:
            print_stylized('STRUCTURES')
            logging.info('Starting structures analysis')

            logging.info(f'Species Tax ID: {tax_id}')
            logging.info(f'Strain Tax ID: {strain_taxid}')
            
            # Run complete structure pipeline: UniProt mapping, structure download, and pocket detection
            df_structures = structures.pipeline_structures(base_path, organism_name, tax_id, strain_taxid, cpus=cpus)
            logging.info('Structures analysis finished')
            tables.append(df_structures)

        except Exception as e:
            logging.error(f'Error in structures analysis: {e}')
    else:
        logging.info('Structures analysis not enabled')

    # Run CORE ANALYSIS
    if config.core:
        try:
            print_stylized('CORE ANALYSIS')
            
            logging.info('Starting core analysis')

            # Download complete NCBI genomes from organism tax id
            print('----- 1. Downloading tax_id genomes from NCBI -----')   
            genome.core_download_genomes_ncbi(base_path, organism_name, tax_id)
            genome.core_download_missing_accessions(base_path, organism_name, tax_id)
            logging.info('Genomes downloaded')

            # Keep genomes with human as host. Check presence of .gff and .faa files for each strain
            print('----- 1. Selecting genomes -----')
            genome.core_check_files(base_path, organism_name)
            logging.info('Genomes filtered')
            print('----- 1. Finished -----')

            min_identity = config.core['min_identity']
            min_core_freq = config.core['min_core_freq']

            if config.core['roary']:
                try:
                    #Run roary
                    print('----- 2. Running Roary -----')
                    logging.info('Starting Roary analysis')
                    genome.core_genome_programs(base_path, organism_name, min_core_freq, min_identity, cpus, program_list=['roary'])
                    # Parse output
                    print('----- 2. Parsing Roary results -----')
                    df_roary = genome.roary_output(base_path, organism_name)
                    tables.append(df_roary)
                    logging.info('Roary analysis finished')
                    print('----- 2. Finished -----')
                except Exception as e:
                    logging.error(f'Error in Roary analysis: {e}')
            else:
                logging.info('Roary not enabled')

            if config.core['corecruncher']:
                try:
                    # Run CoreCruncher
                    print('----- 2. Running CoreCruncher -----')
                    logging.info('Starting CoreCruncher analysis')
                    genome.core_genome_programs(base_path, organism_name, min_core_freq, min_identity, cpus, program_list=['corecruncher'])
                    # Parse output
                    print('----- 2. Parsing CoreCruncher results -----')
                    df_cc = genome.corecruncher_output(base_path, organism_name)
                    tables.append(df_cc)
                    logging.info('CoreCruncher analysis finished')
                    print('----- 2. Finished -----')
                except Exception as e:
                    logging.error(f'Error in CoreCruncher analysis: {e}')
            else:
                logging.info('CoreCruncher not enabled')
        except Exception as e:
            logging.error(f'Error in core analysis: {e}')
    else:
        logging.info('Core analysis not enabled')

    # Run OFFTARGETS
    if config.offtarget:
        try:
            offtarget_path = os.path.join(base_path, 'organism', f'{organism_name}', 'offtarget')

            if config.offtarget['human']:
                try:
                    print_stylized('HUMAN OFFTARGET')

                    human_blast_output = os.path.join(offtarget_path, 'human_offtarget_blast.tsv')
                    
                    if not files.file_check(human_blast_output):
                        # Run blastp search
                        print('-----  Blastp search -----')
                        offtargets.human_offtarget_blast(base_path, organism_name, cpus)
                        logging.info('Human offtarget blast search finished')
                    else:
                        logging.info('Blast with human already done')
                        print('Blast output file found')
                        print(human_blast_output)

                    # Parse results
                    df_human = offtargets.human_offtarget_parse(base_path, organism_name)
                    tables.append(df_human)
                    logging.info('Human offtarget analysis finished')
                    print('----- Finished -----')
                except Exception as e:
                    logging.error(f'Error in human offtarget analysis: {e}')
            else:
                logging.info('Human offtarget analysis not enabled')

            if config.offtarget['microbiome']:
                try:

                    print_stylized('MICROBIOME OFFTARGET')

                    # Run blastp search
                    print('----- Blastp search -----')
                    offtargets.microbiome_offtarget_blast_species(base_path, organism_name, cpus)
                    logging.info('Microbiome offtarget blast search finished')

                    # Parse results
                    microbiome_identity_filter = config.offtarget['microbiome_identity_filter']
                    microbiome_coverage_filter = config.offtarget['microbiome_coverage_filter']
                    logging.info(f'Microbiome identity filter: {microbiome_identity_filter}')
                    logging.info(f'Microbiome coverage filter: {microbiome_coverage_filter}')
                    df_microbiome = offtargets.microbiome_species_parse(base_path, organism_name, microbiome_identity_filter, microbiome_coverage_filter)
                    tables.append(df_microbiome)
                    logging.info('Microbiome offtarget analysis finished')
                    print('----- Finished -----')
                            
                except Exception as e:
                    logging.error(f'Error in microbiome offtarget analysis: {e}')
            else:
                logging.info('Microbiome offtarget analysis not enabled')
            
            if config.structures and config.offtarget['foldseek_human']:
                try:
                    print_stylized('FOLDSEEK HUMAN OFFTARGET')

                    # Get annotation information from UniProt proteome and links IDs with the genome locus_tags of the organism
                    uniprot_proteome_annotations, id_equivalences = structures.uniprot_proteome(base_path, organism_name, proteome_uniprot, cpus=8)
                    # Run foldseek against human structures
                    offtargets.run_foldseek_human_structures (base_path, organism_name)
                    logging.info('Foldseek human offtarget search finished')
                    # Parse results
                    results_foldseek_dict = offtargets.foldseek_human_parser (base_path, organism_name)
                    mapped_dict_foldseek = offtargets.merge_foldseek_data (base_path, organism_name, id_equivalences, uniprot_proteome_annotations)
                    final_foldseek_df = offtargets.final_foldseek_structure_table (base_path, organism_name, mapped_dict_foldseek)
                    tables.append(final_foldseek_df)
                    logging.info('Foldseek human offtarget analysis finished')
                    print('----- Finished -----')
                except Exception as e:
                    logging.error(f'Error in foldseek human offtarget analysis: {e}')
            else:
                logging.info('Foldseek human offtarget analysis not enabled')           

        except Exception as e:
            logging.error(f'Error in offtarget analysis: {e}')
    else:
        logging.info('Offtarget analysis not enabled')

    # Run ESSENTIALITY
    if config.deg:
        essentiality_path = os.path.join(base_path, 'organism', f'{organism_name}', 'essentiality')
        try:
            print_stylized('ESSENTIALITY')

            deg_blast_output = os.path.join(essentiality_path, 'deg_blast.tsv')
            
            if not files.file_check(deg_blast_output):
                # Run blastp search
                print('----- Blastp search -----')
                essentiality.essential_deg_blast(base_path, organism_name, cpus)
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
            df_deg = essentiality.deg_parse(base_path, organism_name, deg_identity_filter, deg_coverage_filter)
            tables.append(df_deg)
            logging.info('DEG analysis finished')
            print('----- Finished -----')

        except Exception as e:
            logging.error(f'Error in essentiality analysis: {e}')
    else:
        logging.info('Essentiality analysis not enabled')

    # Run LOCALIZATION
    if config.psortb:
        try:
            print_stylized('LOCALIZATION')

            gram_type = config.psortb['gram_type']
            logging.info(f'Gram type: {gram_type}')
            
            #Run psortb
            print('----- Running psort -----')
            df_psort = genome.localization_prediction(base_path, organism_name, gram_type)
            tables.append(df_psort)
            logging.info('Psortb analysis finished')
            print('----- Finished -----')
        except Exception as e:
            logging.error(f'Error in localization analysis: {e}')
    else:
        logging.info('Localization analysis not enabled')

    # Load METADATA
    if config.metadata:
        try:
            print_stylized('METADATA')
            for table in config.metadata['meta_tables']:
                print(f'----- Loading metadata table: {table} -----')
                shutil.copy(table, os.path.join(base_path, 'organism', organism_name, 'metadata'))
                with open(table, 'r') as file:
                    first_line = file.readline()
                    if '\t' in first_line:
                        sep = '\t'
                    elif ',' in first_line:
                        sep = ','
                    else:
                        raise ValueError('Invalid file format. Only CSV and TSV metadata files are supported.')

                df_meta = pd.read_csv(table, header=0, sep=sep)
                tables.append(df_meta)
                logging.info(f'Metadata table {table} loaded')
                print('----- Finished -----')
        except Exception as e:
            logging.error(f'Error in metadata analysis: {e}')
    else:
        logging.info('Metadata analysis not enabled')

    # Merge dfs
    print_stylized('RESULTS')
    current_date = datetime.now().strftime('%Y-%m-%d-%H-%M')
    results_path = os.path.join(base_path, 'organism', organism_name, f'{organism_name}_results_{current_date}')
    results_table_path = os.path.join(base_path, 'organism', organism_name, results_path, f'{organism_name}_results_table.tsv')

    if not os.path.exists(results_path):
        os.makedirs(results_path, exist_ok=True)
        logging.info(f'Results directory created in {results_path}')
    

    if len(tables) > 1:
        combined_df = tables[0]
        for df in tables[1:]:
            if df is not None:
                combined_df = pd.merge(combined_df, df, on='gene')
        
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

    print('------------------------------------- FINISHED ----------------------------------------')
    
    return results 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='FastTarget script')
    parser.add_argument('--config_file', type=str, default='config.yml', help='Path to the configuration file')
    args = parser.parse_args()

    config = configuration.get_config(args.config_file)
    base_path = os.path.dirname(os.path.abspath(__file__))

    main(config, base_path)
