from ftscripts import programs, metadata, files
import os
import shutil
import sys
import re
import numpy
import pandas as pd
import logging
import hashlib
import networkx as nx
from networkx.algorithms.centrality.betweenness import betweenness_centrality
from networkx.algorithms.components.connected import connected_components
from contextlib import contextmanager


@contextmanager
def metabolism_step(number, title):
    """
    Announce a numbered step and report whether it finished or failed.

    Every step used to print 'Finished' unconditionally, so a log read as a list
    of completed steps even where one of them had failed. The failure is now
    named at the step that caused it, and then re-raised for the caller to
    record.

    :param number: Step number, as shown in the log.
    :param title: What the step does.
    """
    print(f'---------- {number}. {title} ----------')
    try:
        yield
    except Exception as error:
        print(f'---------- {number}. FAILED: {title} ----------')
        logging.error(f'Metabolism step {number} ({title}) failed: {error}')
        raise
    print(f'---------- {number}. Finished ----------')


# Functions to process Pathway Tools output files

def process_sbml(output_path, organism_name, sbml_file):
    """
    Process the SBML file to generate the Ubiquitous Compounds file inside the metabolism directory.
    This function uses the `run_ubiquitous` function from the `programs` module.

    :param output_path: Path where output will be stored.
    :param organism_name: Name of the organism.
    :param sbml_file: Path to the SBML file
    """
    
    metabolism_dir = os.path.join(output_path, organism_name, 'metabolism')
    sbml_path = os.path.join(metabolism_dir, f'{organism_name}.sbml')

    if not os.path.exists(sbml_file):
        raise FileNotFoundError(f"SBML file '{sbml_file}' not found.")

    shutil.copyfile(sbml_file, sbml_path)

    print(f'Processing {sbml_path}')

    programs.run_ubiquitous(sbml_path, metabolism_dir)

    # The message used to be printed whether or not the file had been written, so
    # a failed run still reported a saved file.
    ubiquitous_path = os.path.join(metabolism_dir, 'ubiquitous_compounds.txt')
    if not files.file_check(ubiquitous_path):
        raise FileNotFoundError(
            f'{ubiquitous_path} was not written, although generating it reported no error.'
        )
    compounds = files.file_to_list(ubiquitous_path)
    print(f'Ubiquitous compounds: {len(compounds)} saved to {ubiquitous_path}')

def ubiquitous_checker(output_path, organism_name, curated_ubiquitous_file=None):
    """
    Process the ubiquitous compounds file for metabolic network generation.
    If a curated file is provided, it will be copied to the expected location.
    If not provided, the auto-generated file will be used with a warning.

    :param output_path: Path where output will be stored.
    :param organism_name: Name of the organism.
    :param curated_ubiquitous_file: Optional path to a pre-curated ubiquitous compounds file.
    """

    metabolism_dir = os.path.join(output_path, organism_name, 'metabolism')
    ubiquitous_file = os.path.join(metabolism_dir, f'{organism_name}_ubiquitous.txt')
    auto_generated_file = os.path.join(metabolism_dir, 'ubiquitous_compounds.txt')
    
    if curated_ubiquitous_file:
        # User provided a curated file
        if os.path.exists(curated_ubiquitous_file):
            shutil.copy(curated_ubiquitous_file, ubiquitous_file)
            print(f"Using curated ubiquitous compounds file: {curated_ubiquitous_file}")
            logging.info(f"Copied curated ubiquitous file to {ubiquitous_file}")
        else:
            logging.error(f"Curated ubiquitous file not found: {curated_ubiquitous_file}")
            raise FileNotFoundError(f"Curated ubiquitous file not found: {curated_ubiquitous_file}")
    
    if not curated_ubiquitous_file or not os.path.exists(curated_ubiquitous_file):
        # No curated file provided, use auto-generated one
        if os.path.exists(auto_generated_file):
            shutil.copy(auto_generated_file, ubiquitous_file)
            print(f"No curated ubiquitous file provided. Using auto-generated file.")
            logging.warning("Using auto-generated ubiquitous compounds file without manual curation. "
                          "For better results, consider providing a curated file via the config.")
        else:
            logging.error(f"Auto-generated ubiquitous file not found: {auto_generated_file}")
            raise FileNotFoundError(f"Auto-generated ubiquitous file not found: {auto_generated_file}")

def generate_sif(output_path, organism_name):
    """
    Generate the network.sif file from the SBML file and the Ubiquitous Compounds file.
    This function uses the `run_sbml_to_sif` function from the `programs` module.
    
    :param output_path: Path where output will be stored.
    :param organism_name: Name of the organism.
    """

    metabolism_dir = os.path.join(output_path, organism_name, 'metabolism')
    sbml_path = os.path.join(metabolism_dir, f'{organism_name}.sbml')
    ubiquitous_file = os.path.join(metabolism_dir, f'{organism_name}_ubiquitous.txt')

    sif_path = os.path.join(metabolism_dir, 'network.sif')
    previous_sif = os.path.getmtime(sif_path) if os.path.exists(sif_path) else None

    programs.run_sbml_to_sif(sbml_path, ubiquitous_file, metabolism_dir)

    # The message used to be printed even when generation had failed, and step 5
    # then read whatever network.sif happened to be in the directory: a stale file
    # from an earlier run produced centrality values for the current one.
    if not files.file_check(sif_path):
        raise FileNotFoundError(f'{sif_path} was not written by the SIF generation step.')
    if previous_sif is not None and os.path.getmtime(sif_path) == previous_sif:
        raise RuntimeError(
            f'{sif_path} was not rewritten by the SIF generation step, so it is left '
            f'over from an earlier run and does not describe the current SBML file.'
        )
    print(f'Sif file generated in {metabolism_dir}.')

def process_chokepoint_ptools(output_path, organism_name, chokepoint_file):
    """
    Process the chokepoint file generated by Pathway Tools.
    This function reads the chokepoint file, extracts the producing and consuming chokepoints, and returns the sets of reactions.
    
    :param output_path: Path where output will be stored.
    :param organism_name: Name of the organism.
    :param chokepoint_file: Path to the chokepoint file.

    :return: Sets of producing chokepoints.
    :return: Sets of consuming chokepoints.
    :return: Sets of both producing and consuming chokepoints.
    """

    metabolism_dir = os.path.join(output_path, organism_name, 'metabolism')
    chokepoint_path = os.path.join(metabolism_dir, f'{organism_name}_chokepoints.txt')

    # The error used to be logged and the copy from a previous run read instead,
    # so a wrong path in the config produced the previous organism's chokepoints
    # without a word.
    if not os.path.exists(chokepoint_file):
        raise FileNotFoundError(f"Chokepoint file '{chokepoint_file}' not found.")

    shutil.copy(chokepoint_file, chokepoint_path)

    with open(chokepoint_path, 'r') as file:
        text = file.read()

    # Define regular expressions for finding sections
    producing_pattern = r"Reactions that are chokepoints on the producing side \(\d+\):\n((?:.*\n)*)\nReactions that are chokepoints on the consuming side"
    consuming_pattern = r"Reactions that are chokepoints on the consuming side \(\d+\):\n((?:.*\n)*)\nReactions that are not chokepoints, but are candidates"

    # Extract the producing side reactions
    producing_match = re.search(producing_pattern, text)
    producing_reactions = producing_match.group(1).strip().split('\n') if producing_match else []

    # Extract the consuming side reactions
    consuming_match = re.search(consuming_pattern, text)
    consuming_reactions = consuming_match.group(1).strip().split('\n') if consuming_match else []

    # ''.split('\n') is [''], so an empty section used to yield one chokepoint
    # named '', and the counts printed below reported it as a real one.
    producing_list = [
        reaction.split('  ', 1)[0] for reaction in producing_reactions
        if reaction.strip()
    ]
    consuming_list = [
        reaction.split('  ', 1)[0] for reaction in consuming_reactions
        if reaction.strip()
    ]

    producing_set = set(producing_list)
    consuming_set = set(consuming_list)
    
    print('Total producing chokepoints')
    print(len(producing_set))
    print('Total consuming chokepoints')
    print(len(consuming_set))
    
    both_set = producing_set & consuming_set
    producing_only_set = producing_set - both_set
    consuming_only_set = consuming_set - both_set
    
    print('Only producing chokepoints')
    print(len(producing_only_set))
    print('Only consuming chokepoints')
    print(len(consuming_only_set))
    print('Both producing and consuming chokepoints')
    print(len(both_set))
    
    return producing_only_set, consuming_only_set, both_set

def process_genes_smarttable(output_path, organism_name, smarttable_file):
    """
    This function reads the Smarttable file, extracts the genes and their associated reactions, and returns a dictionary with the genes as keys and the reactions as values.
    
    :param output_path: Path where output will be stored.
    :param organism_name: Name of the organism.
    :param smarttable_file: Path to the Smarttable file.
    
    :return: Dictionary with the genes as keys and the reactions as values.
    """

    metabolism_dir = os.path.join(output_path, organism_name, 'metabolism')
    smarttable_path = os.path.join(metabolism_dir, f'{organism_name}_smarttable.txt')

    # As for the chokepoint file: reading the previous run's copy is worse than
    # stopping, because the gene-to-reaction mapping would be the wrong organism's.
    if not os.path.exists(smarttable_file):
        raise FileNotFoundError(f"Smarttable file '{smarttable_file}' not found.")

    shutil.copy(smarttable_file, smarttable_path)

    df_genes = pd.read_csv(smarttable_path, delimiter='\t', usecols=['ID', 'Reactions'])
    df_genes = df_genes.dropna(subset=['Reactions'], how='all')
    reaction_dict = {row['ID']: [reaction.strip() for reaction in row['Reactions'].split(',')] for index, row in df_genes.iterrows()}

    return reaction_dict

def write_sif(graph, file_path):
    """
    Write a networkx graph to a SIF file.

    :param graph: Networkx graph.
    :param file_path: Path to the SIF file.
    """
    with open(file_path, 'w') as file:
        for edge in graph.edges():
            file.write(f"{edge[0]}\tinteracts\t{edge[1]}\n")

def read_sif(output_path, organism_name):
    """
    This function reads the network.sif file, generates the principal component of the network, and calculates the betweenness centrality and node degrees.

    :param output_path: Path where output will be stored.
    :param organism_name: Name of the organism.

    :return: The complete network.
    :return: The principal component of the network.
    :return: Dictionary with the betweenness centrality values.
    :return: Dictionary with the node degrees.
    """

    metabolism_dir = os.path.join(output_path, organism_name, 'metabolism')
    sif_path = os.path.join(metabolism_dir, 'network.sif')

    G = nx.DiGraph()
    with open(sif_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 3:
                source, interaction, target = parts
                G.add_edge(source, target, interaction=interaction)
                
    # A network with no edges has no components, and max() over nothing raised
    # ValueError from here, which reached the caller as an unexplained failure of
    # the whole metabolic module.
    if G.number_of_edges() == 0:
        raise ValueError(
            f'{sif_path} contains no usable interactions, so there is no metabolic '
            f'network to analyse. Check that the SIF generation step produced a '
            f'network from the SBML file.'
        )

    largest_component = max(nx.weakly_connected_components(G), key=len)
    subgraph = G.subgraph(largest_component).copy()
    write_sif(subgraph, os.path.join(metabolism_dir, f'{organism_name}_principal_component.sif'))

    betweenness_centrality = nx.betweenness_centrality(subgraph)
    # Every centrality is zero whenever no node lies on a shortest path between
    # two others, which is the case for any component small enough: a two-node
    # component, a self-loop, a star. That is a legitimate result, but dividing
    # by it raised ZeroDivisionError and lost the whole module.
    max_value = max(betweenness_centrality.values(), default=0)
    if max_value > 0:
        normalized_betweenness_centrality = {
            key: value / max_value for key, value in betweenness_centrality.items()
        }
    else:
        print(
            f'No node of the principal component ({subgraph.number_of_nodes()} nodes) lies '
            f'between two others, so every betweenness centrality is zero.'
        )
        normalized_betweenness_centrality = dict.fromkeys(betweenness_centrality, 0.0)

    node_degrees = dict(subgraph.degree())
    
    return G, subgraph, normalized_betweenness_centrality, node_degrees   

def mapping_metabolism(locus_reaction_dict, reaction_property_dict):
    """
    This function maps the betweenness centrality values to the genes.

    :param locus_reaction_dict: Dictionary with the genes as keys and the reactions as values.
    :param reaction_property_dict: Dictionary with the reactions as keys and the betweenness centrality values as values.

    :return: Dictionary with the genes as keys and the betweenness centrality values as values.
    """
    
    locustag_to_values = {}

    for gene, rxns in locus_reaction_dict.items():
        values = [reaction_property_dict[rxn] for rxn in rxns if rxn in reaction_property_dict]
        locustag_to_values[gene] = values
        
    locustag_max = {}
    for gene, values in locustag_to_values.items():
        if len(values)>0:
            locustag_max[gene] = max(values)
        else:
            locustag_max[gene] = 0   
    
    return locustag_max

def mapping_chokepoints(locus_reaction_dict, producing_ck, consuming_ck, both_ck):
    """
    This function maps the chokepoints to the genes.

    :param locus_reaction_dict: Dictionary with the genes as keys and the reactions as values.
    :param producing_ck: Set of producing chokepoints.
    :param consuming_ck: Set of consuming chokepoints.
    :param both_ck: Set of both producing and consuming chokepoints.

    :return: Dictionaries with the genes as keys and the chokepoints as values. Returns three dictionaries: producing chokepoints, consuming chokepoints, and both chokepoints.
    """
    
    producing_chokepoint = {}
    consuming_chokepoint = {}
    both_chokepoint = {}
        
    for gene, rxns in locus_reaction_dict.items():
        for rxn in rxns:
            if rxn in producing_ck:
                producing_chokepoint[gene] = rxn
            elif rxn in consuming_ck:
                consuming_chokepoint[gene] = rxn
            elif rxn in both_ck:
                both_chokepoint[gene] = rxn

    return producing_chokepoint, consuming_chokepoint, both_chokepoint

def run_metabolism_ptools (output_path, organism_name, sbml_file, chokepoint_file, smarttable_file, curated_ubiquitous_file=None, network_settings=None):
    """
    Runs the metabolism pipeline. 
    It processes the SBML file, generates the Ubiquitous Compounds file, the network SIF file, and the gene-reaction mapping. 
    Returns DataFrames with locus_tag as key and the values of betweenness centrality, node degrees, and chokepoints.
    The DataFrames are created using the functions `metadata_table_with_values` and `metadata_table_with_values` from the `metadata` module.   

    :param output_path: Path where the output will be stored.
    :param organism_name: Name of the organism.
    :param sbml_file: Path to the SBML file.
    :param chokepoint_file: Path to the chokepoint file.
    :param smarttable_file: Path to the Smarttable file.
    :param curated_ubiquitous_file: Optional path to a pre-curated ubiquitous compounds file.

    :return: DataFrame with the betweenness centrality values.
    :return: DataFrame with the node degrees.
    :return: DataFrame with the chokepoints.
    """

    metabolism_dir = os.path.join(output_path, organism_name, 'metabolism')

    # Generating the ubiquitous compounds file is the one step whose failure can be
    # survived, and only because step 2 can fall back on a curated file. Without
    # one there is nothing for step 2 to use, so the failure stands.
    try:
        with metabolism_step(1, 'Processing SBML file'):
            process_sbml(output_path, organism_name, sbml_file)
    except Exception as error:
        if not curated_ubiquitous_file:
            raise
        print(
            f'Continuing with the curated ubiquitous compounds file '
            f'{curated_ubiquitous_file}, since generating one failed: {error}'
        )

    with metabolism_step(2, 'Processing Ubiquitous compounds file'):
        ubiquitous_checker(output_path, organism_name, curated_ubiquitous_file)

    with metabolism_step(3, 'Generating network SIF file'):
        generate_sif(output_path, organism_name)

    with metabolism_step(4, 'Map proteins and reactions'):
        locus_reaction_dict = process_genes_smarttable(output_path, organism_name, smarttable_file)

    with metabolism_step(5, 'Calculate Centrality'):
        graph, subgraph, reaction_centrality_dict, reaction_edges_dict = read_sif(output_path, organism_name)
        bc_dict = mapping_metabolism(locus_reaction_dict, reaction_centrality_dict)
        edge_dict = mapping_metabolism(locus_reaction_dict, reaction_edges_dict)

    with metabolism_step(6, 'Obtaining Choke-points'):
        producing_set, consuming_set, both_set = process_chokepoint_ptools(output_path, organism_name, chokepoint_file)
        producing_chokepoint, consuming_chokepoint, both_chokepoint = mapping_chokepoints(locus_reaction_dict, producing_set, consuming_set, both_set)

    print(f'---------- 7. Generating results files ----------')
    bcentrality_df = metadata.metadata_table_with_values(output_path, organism_name, bc_dict, "PTOOLS_betweenness_centrality", metabolism_dir, 0)
    edges_df = metadata.metadata_table_with_values(output_path, organism_name, edge_dict, "PTOOLS_edges", metabolism_dir, 0)
    producing_df = metadata.metadata_table_with_values(output_path, organism_name, producing_chokepoint, "PTOOLS_producing_chokepoints", metabolism_dir, 'None')
    consuming_df = metadata.metadata_table_with_values(output_path, organism_name, consuming_chokepoint, "PTOOLS_consuming_chokepoints", metabolism_dir, 'None')
    both_df = metadata.metadata_table_with_values(output_path, organism_name, both_chokepoint, "PTOOLS_both_chokepoints", metabolism_dir, 'None')
    print(f'---------- 7. Finished ----------')

    bcentrality_df = add_network_gene_metrics(bcentrality_df, subgraph, locus_reaction_dict, metabolism_dir, 'PTOOLS', network_settings)
    bcentrality_df.to_csv(os.path.join(metabolism_dir, 'PTOOLS_betweenness_centrality.tsv'), sep='\t', index=False)
    return bcentrality_df, edges_df, producing_df, consuming_df, both_df

# Functions for process an external SBML with MetaGraphTools (docker)

def process_external_sbml(sbml_file, output_path, filter_file=None, container_engine='docker'):

    """
    Process an external SBML file using MetaGraphTools container image to generate the Ubiquitous Compounds file.

    :param sbml_file: Path to the SBML file.
    :param output_path: Path to the output directory.
    :param filter_file: Optional path to the filter file for ubiquitous compounds.
    :param container_engine: Container engine to use ('docker' or 'singularity').
    """

    if not os.path.exists(output_path):
        os.makedirs(output_path, exist_ok=True)

    # Move SBML file to output directory
    sbml_basename = os.path.basename(sbml_file)
    sbml_file_destination = os.path.join(output_path, sbml_basename)
    if not files.file_check(sbml_file_destination):
        if files.file_check(sbml_file):
            shutil.copy2(sbml_file, sbml_file_destination)
        else:
            logging.error(f"SBML file '{sbml_file}' not found.")
            return
    else:
        print(f"SBML file '{sbml_file_destination}' already exists in the output directory.")
    
    #Move filter file to output directory if provided
    if filter_file:
        filter_basename = os.path.basename(filter_file)
        filter_file_destination = os.path.join(output_path, filter_basename)
        if not files.file_check(filter_file_destination):
            if files.file_check(filter_file):
                shutil.copy2(filter_file, filter_file_destination)
            else:
                logging.warning(f"Filter file '{filter_file}' not found.")
                print(f'Proceeding without filter file.')
                filter_file_destination = None
        else:
            print(f"Filter file '{filter_file_destination}' already exists in the output directory.")
    else:
        filter_file_destination = None

    try:    
        programs.run_metagraphtools(output_path, sbml_file_destination, filter_file_destination, chokepoints=True, graph=True, container_engine=container_engine)
        print(f'MetaGraphTools processing completed. Results are in {output_path}.')
        try:
            programs.change_permission_user_dir(output_path)
        except Exception as e:
            logging.exception(f"Error changing file permissions: {e}")
    except Exception as e:
        logging.exception(f"Error running MetaGraphTools: {e}")

def parse_mgt_results(mgt_results_dir):
    """
    Parse the results generated by MetaGraphTools.

    :param mgt_results_dir: Path to the MetaGraphTools results directory.

    :return: Tuple of dictionaries (consumption_chokepoint_dict, production_chokepoint_dict, betweenness_centrality_dict, degree_dict).
    
    """
    programs.change_permission_user_dir(mgt_results_dir)

    chokepoint_file = os.path.join(mgt_results_dir, 'chokepoint_genes.tsv')
    betweenness_centrality_file = os.path.join(mgt_results_dir, 'betweenness_centrality.tsv')
    
    consumption_chokepoint_dict = {}
    production_chokepoint_dict = {}
    betweenness_centrality_dict = {}
    degree_dict = {}
    
    # Parse chokepoint genes file
    if os.path.exists(chokepoint_file):
        try:
            chokepoint_df = pd.read_csv(chokepoint_file, sep='\t', usecols=['Gene', 'Consumption_Chokepoint', 'Production_Chokepoint'])
            
            # Convert boolean columns
            for col in ['Consumption_Chokepoint', 'Production_Chokepoint']:
                if col in chokepoint_df.columns:
                    chokepoint_df[col] = chokepoint_df[col].astype(bool)
            
            # Create dictionaries with Gene as key and boolean values
            for _, row in chokepoint_df.iterrows():
                gene = row['Gene']
                if row['Consumption_Chokepoint']:
                    consumption_chokepoint_dict[gene] = True
                if row['Production_Chokepoint']:
                    production_chokepoint_dict[gene] = True
            
            print(f'Parsed {len(chokepoint_df)} chokepoint gene entries.')
            print(f'Consumption chokepoints: {len(consumption_chokepoint_dict)}')
            print(f'Production chokepoints: {len(production_chokepoint_dict)}')
            
        except Exception as e:
            logging.exception(f'Error parsing chokepoint genes file: {e}')
    else:
        logging.error(f'Chokepoint genes file not found: {chokepoint_file}')
    
    # Parse betweenness centrality file
    if os.path.exists(betweenness_centrality_file):
        try:
            betweenness_df = pd.read_csv(betweenness_centrality_file, sep='\t', usecols=['Gene', 'Betweenness_Centrality', 'Degree'])
            
            # Ensure numeric columns are proper floats/ints
            if 'Betweenness_Centrality' in betweenness_df.columns:
                betweenness_df['Betweenness_Centrality'] = betweenness_df['Betweenness_Centrality'].astype(float)
            if 'Degree' in betweenness_df.columns:
                betweenness_df['Degree'] = betweenness_df['Degree'].astype(int)
            
            # Create dictionaries with Gene as key

            betweenness_centrality_dict = betweenness_df.groupby('Gene')['Betweenness_Centrality'].max().to_dict()
            degree_dict = betweenness_df.groupby('Gene')['Degree'].max().to_dict()

            print(f'Parsed {len(betweenness_df)} betweenness centrality entries.')
            
        except Exception as e:
            logging.exception(f'Error parsing betweenness centrality file: {e}')
    else:
        logging.error(f'Betweenness centrality file not found: {betweenness_centrality_file}')
    
    return consumption_chokepoint_dict, production_chokepoint_dict, betweenness_centrality_dict, degree_dict

def pick_mgt_results_dir(metabolism_dir):
    """
    Pick the MetaGraphTools results directory from the metabolism directory.
    Removes incomplete or older duplicate directories, keeping only the most recent valid one.
    
    :param metabolism_dir: Path to the metabolism directory.
    :return: Path to the MetaGraphTools results directory if found, None otherwise.
    """

    def check_mgt_output_files(mgt_results_dir):
        """
        Check if all expected MetaGraphTools output files are present.
        :param mgt_results_dir: Path to the MetaGraphTools results directory.
        :return: True if all expected files are present, False otherwise.
        """
        expected_files = [
            "all_metabolite_frequencies.tsv",
            "betweenness_centrality.tsv",
            "chokepoint_genes.tsv",
            "chokepoint_reactions.tsv",
            "metabolic_network_filter.html",
            "metabolic_network_filter.sif",
            "metabolic_network.html",
            "metabolic_network.sif",
            "metabolite_reactant_product.tsv"
        ]

        missing_files = []
        for filename in expected_files:
            file_path = os.path.join(mgt_results_dir, filename)
            if not os.path.exists(file_path):
                missing_files.append(filename)
        
        if missing_files:
            logging.warning(f'Missing MetaGraphTools output files in {mgt_results_dir}: {", ".join(missing_files)}')
            return False
        return True

    # Find all MGT_results directories
    mgt_dirs = [
        os.path.join(metabolism_dir, d)
        for d in os.listdir(metabolism_dir)
        if d.startswith("MGT_results") and os.path.isdir(os.path.join(metabolism_dir, d))
    ]

    if not mgt_dirs:
        print("No MGT_results folders found.")
        return None

    # Sort by modification time (oldest → newest)
    mgt_dirs.sort(key=lambda d: os.path.getmtime(d))

    # Check each folder and keep track of valid ones
    valid_dirs = []
    for d in mgt_dirs:
        if check_mgt_output_files(d):
            valid_dirs.append(d)
        else:
            try:
                programs.change_permission_user_dir(d)
                shutil.rmtree(d)
                print(f"Removed incomplete MGT_results folder: {d}")
            except Exception as e:
                logging.exception(f"Error removing folder {d}: {e}")

    if not valid_dirs:
        print("No valid MGT_results folders found.")
        return None

    # If multiple valid folders exist, keep only the newest
    if len(valid_dirs) > 1:
        for d in valid_dirs[:-1]:  # Remove all but the last (most recent)
            try:
                shutil.rmtree(d)
                print(f"Removed older valid MGT_results folder: {d}")
            except Exception as e:
                logging.exception(f"Error removing folder {d}: {e}")

    return valid_dirs[-1]  # Return the most recent valid directory

def run_metabolism_sbml (output_path, organism_name, sbml_file, filter_file, container_engine='docker', network_settings=None):

    """
    Runs the metabolism pipeline using an external SBML file with MetaGraphTools.
    It processes the SBML file, generates the Ubiquitous Compounds file, and
    parses the MetaGraphTools results to return DataFrames with locus_tag as key and the values of betweenness centrality,
    node degrees, and chokepoints.
    The DataFrames are created using the functions `metadata_table_with_values` from the `metadata
    module.
    :param output_path: Path where the output will be stored.
    :param organism_name: Name of the organism.
    :param sbml_file: Path to the SBML file.
    :param filter_file: Path to the filter file for ubiquitous compounds.
    :param container_engine: Container engine to use ('docker' or 'singularity').
    :return: DataFrame with the betweenness centrality values.
    :return: DataFrame with the node degrees.
    :return: DataFrame with the consuming chokepoints.
    :return: DataFrame with the producing chokepoints.
    """

    metabolism_dir = os.path.join(output_path, organism_name, 'metabolism')

    print(f'---------- Processing external SBML file with MetaGraphTools ----------')

    mgt_results_dir = pick_mgt_results_dir(metabolism_dir)

    if not mgt_results_dir:
        try:
            process_external_sbml(sbml_file, metabolism_dir, filter_file, container_engine=container_engine)
            # After running MetaGraphTools, try to pick the results directory again
            mgt_results_dir = pick_mgt_results_dir(metabolism_dir)
        except Exception as e:
            logging.exception(f'Error processing external SBML file: {e}')
            return None, None, None, None

    # Parse results if we have a valid MGT results directory
    if mgt_results_dir:
        consumption_chokepoint_dict, production_chokepoint_dict, betweenness_centrality_dict, degree_dict = parse_mgt_results(mgt_results_dir)

        bcentrality_df = metadata.metadata_table_with_values(output_path, organism_name, betweenness_centrality_dict, "MGT_betweenness_centrality", metabolism_dir, 0)
        degree_df = metadata.metadata_table_with_values(output_path, organism_name, degree_dict, "MGT_edges", metabolism_dir, 0)
        consumption_df = metadata.metadata_table_with_values(output_path, organism_name, consumption_chokepoint_dict, "MGT_consuming_chokepoints", metabolism_dir, 'False')
        production_df = metadata.metadata_table_with_values(output_path, organism_name, production_chokepoint_dict, "MGT_producing_chokepoints", metabolism_dir, 'False')

        print(f'---------- Finished processing external SBML file ----------')

        graph = nx.DiGraph()
        with open(os.path.join(mgt_results_dir, 'metabolic_network_filter.sif')) as handle:
            for line in handle:
                fields = line.split()
                if len(fields) == 3:
                    graph.add_edge(fields[0], fields[2])
                elif fields:
                    raise ValueError('Malformed filtered MGT SIF')
        associations = pd.read_csv(os.path.join(mgt_results_dir, 'betweenness_centrality.tsv'),
                                   sep='\t', dtype=str, keep_default_na=False)
        mapping = associations.groupby('Gene')['Reaction'].agg(list).to_dict()
        bcentrality_df = add_network_gene_metrics(bcentrality_df, graph, mapping, metabolism_dir, 'MGT', network_settings)
        bcentrality_df.to_csv(os.path.join(metabolism_dir, 'MGT_betweenness_centrality.tsv'), sep='\t', index=False)
        return bcentrality_df, degree_df, consumption_df, production_df
    else:
        logging.error('MetaGraphTools did not produce valid results.')
        return None, None, None, None




def extended_network_metrics(graph, destination, prefix, settings=None):
    """Node metrics and descriptive global comparisons to directed G(n,m)."""
    import random
    settings = settings or {}
    repetitions = settings.get('random_networks', 100)
    seed = settings.get('seed', 42)
    alpha = settings.get('pagerank_alpha', 0.85)
    if isinstance(repetitions, bool) or not isinstance(repetitions, int) or repetitions < 1:
        raise ValueError('random_networks must be a positive integer')
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 < alpha < 1:
        raise ValueError('seed must be an integer and pagerank_alpha must be in (0,1)')
    graph = nx.DiGraph(graph)
    n, m = graph.number_of_nodes(), graph.number_of_edges()
    loops = nx.number_of_selfloops(graph)

    def measure(network):
        ranks = nx.pagerank(network, alpha=alpha, max_iter=1000, tol=1e-10) if len(network) else {}
        indegree, outdegree = dict(network.in_degree()), dict(network.out_degree())
        summary = {'pagerank_max': max(ranks.values(), default=0.0),
                   'pagerank_std': float(numpy.std(list(ranks.values()))) if ranks else 0.0,
                   'in_degree_std': float(numpy.std(list(indegree.values()))) if indegree else 0.0,
                   'out_degree_std': float(numpy.std(list(outdegree.values()))) if outdegree else 0.0,
                   'weak_components': nx.number_weakly_connected_components(network),
                   'strong_components': nx.number_strongly_connected_components(network)}
        return {'pagerank': ranks, 'in_degree': indegree, 'out_degree': outdegree}, summary

    metrics, observed = measure(graph)
    rng = random.Random(seed)
    null_rows = []
    for iteration in range(repetitions):
        null = nx.gnm_random_graph(n, m - loops, seed=rng, directed=True)
        null.add_edges_from((node, node) for node in rng.sample(range(n), loops))
        _, summary = measure(null)
        null_rows.append({'iteration': iteration, **summary})
    nulls = pd.DataFrame(null_rows)
    comparisons = []
    for metric, value in observed.items():
        values = nulls[metric]
        mean, spread = float(values.mean()), float(values.std(ddof=0))
        comparisons.append({'metric': metric, 'observed': value, 'null_mean': mean,
                            'null_std': spread, 'difference': value - mean,
                            'z_score': (value - mean) / spread if spread > 0 else None,
                            'null_min': float(values.min()), 'null_max': float(values.max())})
    os.makedirs(destination, exist_ok=True)
    pd.DataFrame(metrics).rename_axis('reaction').to_csv(os.path.join(destination, prefix + '_reaction_metrics.tsv'), sep='\t')
    nulls.to_csv(os.path.join(destination, prefix + '_random_networks.tsv'), sep='\t', index=False)
    pd.DataFrame(comparisons).to_csv(os.path.join(destination, prefix + '_network_comparison.tsv'), sep='\t', index=False)
    files.atomic_write_json({'nodes': n, 'edges': m, 'self_loops': loops,
                             'random_networks': repetitions, 'seed': seed, 'pagerank_alpha': alpha,
                             'networkx_version': nx.__version__, 'model': 'directed G(n,m), fixed self-loop count',
                             'pagerank_policy': 'uniform personalization and dangling redistribution',
                             'interpretation': 'descriptive global comparison, no node-level significance',
                             'edges_sha256': hashlib.sha256(
                                 repr(sorted((str(u), str(v)) for u,v in graph.edges())).encode()).hexdigest()},
                            os.path.join(destination, prefix + '_network_parameters.json'))
    return metrics


def add_network_gene_metrics(table, graph, mapping, destination, prefix, settings=None):
    metrics = extended_network_metrics(graph, destination, prefix, settings)
    for name, values in metrics.items():
        # Match the established maximum-over-associated-reactions aggregation;
        # leave genes outside the analyzed graph unavailable rather than zero.
        mapped = {gene: max(values[r] for r in reactions if r in values)
                  for gene, reactions in mapping.items() if any(r in values for r in reactions)}
        table[prefix + '_' + name] = table['gene'].map(mapped)
    return table
