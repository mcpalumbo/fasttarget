
from ftscripts import programs, metadata, files, structures
import os
import pandas as pd
import multiprocessing
import glob
from tqdm import tqdm

def human_offtarget_blast (base_path, organism_name, cpus=multiprocessing.cpu_count()):

    """
    Runs NCBI BLASTP against the human proteome.
    This function uses the `run_blastp` function from the `programs` module.
    It uses the HUMAN_DB database created by the `index_db_blast_human` function from the `databases` module.
    The blast output is saved in the 'offtarget' folder of the organism directory.

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.
    :param cpus: Number of threads (CPUs) to use in the blast search.
    
    """

    #Database files
    databases_path = os.path.join(base_path, 'databases')
    humanprot_index_path = os.path.join(databases_path, 'HUMAN_DB')

    #Organism files
    organism_path = os.path.join(base_path, f'organism/{organism_name}')
    organism_prot_seq_path = os.path.join(organism_path, f'genome/{organism_name}.faa')

    offtarget_path = os.path.join(organism_path, 'offtarget')
    blast_output_path = os.path.join(offtarget_path, 'human_offtarget_blast.tsv')

    programs.run_blastp(
        blastdb= humanprot_index_path,
        query= organism_prot_seq_path,
        output=blast_output_path,
        evalue= '1e-5',
        outfmt= '6 std qcovhsp qcovs',
        cpus=cpus
    )

def microbiome_offtarget_blast (base_path, organism_name, cpus=multiprocessing.cpu_count()):

    """
    Runs ncbi blastp against microbiome proteome.
    This function uses the `run_blastp` function from the `programs` module.
    It uses the MICROBIOME_DB database created by the `index_db_blast_microbiome` function from the `databases` module.
    The blast output is saved in the 'offtarget' folder of the organism directory.

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.
    :param cpus: Number of threads (CPUs) to use in the blast search.
    
    """
    
    #Database files
    databases_path = os.path.join(base_path, 'databases')
    microbiome_index_path = os.path.join(databases_path, 'MICROBIOME_DB')

    #Organism files
    organism_path = os.path.join(base_path, f'organism/{organism_name}')
    organism_prot_seq_path = os.path.join(organism_path, f'genome/{organism_name}.faa')

    offtarget_path = os.path.join(organism_path, 'offtarget')
    blast_output_path = os.path.join(offtarget_path, 'microbiome_offtarget_blast.tsv')

    programs.run_blastp(
        blastdb= microbiome_index_path,
        query= organism_prot_seq_path,
        output=blast_output_path,
        evalue= '1e-5',
        outfmt= '6 std qcovhsp qcovs',
        max_target_seqs = '1000',
        cpus=cpus
    )

def read_blast_output(file_path):
    """
    Read BLASTP output file and return a pandas DataFrame.

    :param file_path: File path of blast output file.

    :return: Pandas DataFrame with blast output.
    """

    if os.path.exists(file_path):
        blast_output_df = pd.read_csv(file_path, sep='\t', header=None)

        # Columns used for blastp
        blast_output_df.columns = [
        "qseqid",   # query or source (gene) sequence id
        "sseqid",   # subject or target (reference genome) sequence id
        "pident",   # percentage of identical positions
        "length",   # alignment length (sequence overlap)
        "mismatch", # number of mismatches
        "gapopen",  # number of gap openings
        "qstart",   # start of alignment in query
        "qend",     # end of alignment in query
        "sstart",   # start of alignment in subject
        "send",     # end of alignment in subject
        "evalue",   # expect value
        "bitscore", # bit score
        "qcovhsp",  # Query Coverage hsp
        "qcovs"     # Query Coverage full
        ]

        return blast_output_df
    
    else:
        print(f'File {file_path} not found.')

def human_offtarget_parse (base_path, organism_name):

    """
    Parse NCBI BLASTP results against human proteome, stored in the file 'human_offtarget_blast.tsv'.
    Obtains the hit with the highest percentage of identity for each locus_tag.
    Returns a dictionary with locus_tag as key and highest percentage of identity as value, 
    and a DataFrame with all locus_tags from the genome and their respective values.
    The DataFrame is created using the `metadata_table_with_values` function from the `metadata` module.

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.
    
    :return: Dictionary with locus_tag as key and highest percentage of identity value.
    :return: DataFrame with all locus_tags from the genome and their respective values.
    """

    offtarget_path = os.path.join(base_path, 'organism', f'{organism_name}', 'offtarget')
    human_blast_output = os.path.join(offtarget_path, 'human_offtarget_blast.tsv')
    human_results = os.path.join(offtarget_path, 'human_offtarget.tsv')

    if not files.file_check(human_results):
        blast_output_df = read_blast_output(human_blast_output)

        highest_pident_values = {}

        for index,row in blast_output_df.iterrows():
            qseqid = row['qseqid']
            pident = row['pident']

            if qseqid not in highest_pident_values or pident > highest_pident_values[qseqid]:
                highest_pident_values[qseqid] = pident

        df_human = metadata.metadata_table_with_values(base_path, organism_name, highest_pident_values, 
                                            'human_offtarget', offtarget_path, 'no_hit')
    else:
        print('Human offtarget analysis already done, output file found')
        print(human_results)
        df_human = pd.read_csv(human_results, sep='\t', header=0)

    return df_human

def microbiome_offtarget_parse (base_path, organism_name, identity_filter, coverage_filter):

    """
    Parse NCBI BLASTP results against microbiome proteome, stored in the file 'microbiome_offtarget_blast.tsv'.
    Filters results based on identity and coverage thresholds, then counts occurrences of each locus_tag.
    Returns a dictionary with normalized counts of each locus_tag that has at least one hit with UHGP90, 
    and a DataFrame with all locus_tags from the genome and their respective values.
    The DataFrame is created using the `metadata_table_with_values` function from the `metadata` module.

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.
    :param identity_filter: Percentage identity filter value. Keeps results above this value in the pident column.
    :param coverage_filter: Query coverage filter value. Keeps results above this value in the qcovs column.
    
    :return: Dictionary with normalized counts of locus_tags that have at least one hit with UHGP90.
    :return: DataFrame with all locus_tags from the genome and their respective normalized values.
    """

    offtarget_path = os.path.join(base_path, 'organism', f'{organism_name}', 'offtarget')
    microbiome_blast_output = os.path.join(offtarget_path, 'microbiome_offtarget_blast.tsv')
    microbiome_results = os.path.join(offtarget_path, 'gut_microbiome_offtarget.tsv')

    if not files.file_check(microbiome_results):

        blast_output_df = read_blast_output(microbiome_blast_output)

        #Filter % identity and coverage
        filtered_df = blast_output_df[(blast_output_df['pident'] > identity_filter) 
                                    & (blast_output_df['qcovs'] > coverage_filter)]

        value_counts = filtered_df['qseqid'].value_counts()
        max_count = value_counts.max()
        norm_counts = value_counts / max_count

        normalized_counts_dict = norm_counts.to_dict()

        df_microbiome = metadata.metadata_table_with_values(base_path, organism_name, normalized_counts_dict, 
                                            'gut_microbiome_offtarget', offtarget_path, 'no_hit')

    else:
        print('Microbiome offtarget analysis already done, output file found')
        print(microbiome_results)
        df_microbiome = pd.read_csv(microbiome_results, sep='\t', header=0)
        
    return df_microbiome

def run_foldseek_human_structures (base_path, organism_name):

    """
    Runs Foldseek easy-search for human proteome PDB and AlphaFold structures.
    Returns a tsv file in databases/human_structures/PDB_files and databases/human_structures/AlphaFold_files for each structure.

    :param base_path =  Base path of fasttarget folder.

    """

    # Human databases of PDB and AlphaFold structures
    db_human_PDB_path = os.path.join(base_path, 'databases', 'human_structures', 'PDB_files', 'DB_foldseek')
    db_human_AF_path = os.path.join(base_path, 'databases', 'human_structures', 'AlphaFold_files', 'DB_foldseek')

    # Organism files of PDB and AlphaFold structures
    structure_dir = os.path.join(base_path, 'organism', f'{organism_name}', 'structures')

    structures_PDB_path = os.path.join(structure_dir, 'PDB_files')
    structures_AF_path = os.path.join(structure_dir, 'AlphaFold_files')

    PDB_files = glob.glob(os.path.join(structures_PDB_path, '*.pdb'))
    AF_files = glob.glob(os.path.join(structures_AF_path, '*.pdb'))

    #Run Foldseek PDB structures
    for pdb in PDB_files:
        pdb_name = os.path.basename(pdb)
        #Search with human PDB database
        try:
            programs.run_foldseek_search(structures_PDB_path, db_human_PDB_path, 'DB_human_PDB', pdb_name)
        except Exception as e:
            print(f'Error running Foldseek for PDB structure {pdb_name} with human PDB database: {e}')
        #Search with human AlphaFold database
        try:
            programs.run_foldseek_search(structures_PDB_path, db_human_AF_path, 'DB_human_AF', pdb_name)
        except Exception as e:
            print(f'Error running Foldseek for PDB structure {pdb_name} with human AlphaFold database: {e}')

    #Run Foldseek AlphaFold structures
    for AF in AF_files:
        AF_name = os.path.basename(AF)
        #Search with human PDB database
        try:
            programs.run_foldseek_search(structures_AF_path, db_human_PDB_path, 'DB_human_PDB', AF_name)
        except Exception as e:
            print(f'Error running Foldseek for AlphaFold structure {AF_name} with human PDB database: {e}')
        #Search with human AlphaFold database
        try:
            programs.run_foldseek_search(structures_AF_path, db_human_AF_path, 'DB_human_AF', AF_name)
        except Exception as e:
            print(f'Error running Foldseek for AlphaFold structure {AF_name} with human AlphaFold database: {e}')

    #Check results
    foldseek_pdb_files = glob.glob(os.path.join(structures_PDB_path, 'foldseek_results', '*foldseek_results.tsv'))
    foldseek_AF_files = glob.glob(os.path.join(structures_AF_path, 'foldseek_results', '*foldseek_results.tsv'))

    print(f'Of a total of {len(PDB_files)} PDB structures: number of Foldseek PDB results:', len(foldseek_pdb_files)/2)
    print(f'Of a total of {len(AF_files)} AlphaFold structures: number of Foldseek AlphaFold results:', len(foldseek_AF_files)/2)

    if len(PDB_files) == len(foldseek_pdb_files)/2:
        print('Foldseek PDB search completed successfully.')
    else:
        print('Missing PDB Foldseek results.')

    if len(AF_files) == len(foldseek_AF_files)/2:
        print('Foldseek AlphaFold search completed successfully.')
    else:
        print('Missing AlphaFold Foldseek results.')

def foldseek_human_parser (base_path, organism_name):

    """
    Parses Foldseek results for human proteome PDB and AlphaFold structures.
    Returns a dictionary with the results.

    :param base_path =  Base path of fasttarget folder.
    
    :return: results_foldseek_dict = Dictionary with Foldseek results of search against human proteome structures.
    """

    # Organism files of PDB and AlphaFold structures
    structure_dir = os.path.join(base_path, 'organism', f'{organism_name}', 'structures')

    if not files.file_check(os.path.join(structure_dir, 'human_foldseek_dict.json')):
  
        # PDB structures
        foldseek_pdb_results = glob.glob(os.path.join(structure_dir, 'PDB_files', 'foldseek_results', '*foldseek_results.tsv'))

        results_foldseek_dict = {}

        queries_pdb = {}

        #Parse PDB results
        for pdb_res in foldseek_pdb_results:
          
            file_name = os.path.splitext(os.path.basename(pdb_res))[0]
            query = file_name.split('_')[1].replace('.pdb','')

            if query not in queries_pdb:
                queries_pdb[query] = []
            
            queries_pdb[query].append(pdb_res)


        for query, files_query in queries_pdb.items():
            dfs = []
            for file in files_query:
                df = pd.read_csv(file, sep='\t', usecols=['query', 'target', 'rmsd', 'prob', 'pident'])
                if not df.empty:
                    dfs.append(df)
            
            if len(dfs) > 0:
                dfs_combined = pd.concat(dfs, ignore_index=True)
                best_row = dfs_combined.sort_values(by='prob', ascending=False).iloc[0]
                pdb_id = best_row['query'].split('_')[1]

                if pdb_id == query:
                        
                    results_foldseek_dict[query] = {
                            'target_foldseek': best_row['target'],
                            'rmsd_foldseek': best_row['rmsd'],
                            'prob_foldseek': best_row['prob'],
                            'pident_foldseek': best_row['pident']
                        }
            else:
                results_foldseek_dict[query] = {
                            'target_foldseek': None,
                            'rmsd_foldseek': None,
                            'prob_foldseek': None,
                            'pident_foldseek': None
                        }

        # AlphaFold structures
        structures_AF_path = os.path.join(structure_dir, 'AlphaFold_files')
        foldseek_results_AF_path = os.path.join(structures_AF_path, 'foldseek_results')

        foldseek_AF_results = glob.glob(os.path.join(foldseek_results_AF_path, '*foldseek_results.tsv'))

        queries_AF = {}

        #Parse AlphaFold results
        for AF_res in foldseek_AF_results:
             
            file_name = os.path.splitext(os.path.basename(AF_res))[0]
            query = file_name.split('_')[1].replace('.pdb','')

            if query not in queries_AF:
                queries_AF[query] = []

            queries_AF[query].append(AF_res)
        
        for query, files_query in queries_AF.items():
            dfs = []
            for file in files_query:
                df = pd.read_csv(file, sep='\t', usecols=['query', 'target', 'rmsd', 'prob', 'pident'])
                if not df.empty:
                    dfs.append(df)
            
            if len(dfs) > 0:
                dfs_combined = pd.concat(dfs, ignore_index=True)
                best_row = dfs_combined.sort_values(by='prob', ascending=False).iloc[0]
                af_id = best_row['query'].split('_')[1]

                if af_id == query and not query in results_foldseek_dict.keys():
                        
                    results_foldseek_dict[query] = {
                            'target_foldseek': best_row['target'],
                            'rmsd_foldseek': best_row['rmsd'],
                            'prob_foldseek': best_row['prob'],
                            'pident_foldseek': best_row['pident']
                        }
            elif len(dfs) == 0 and not query in results_foldseek_dict.keys():
                results_foldseek_dict[query] = {
                        'target_foldseek': None,
                        'rmsd_foldseek': None,
                        'prob_foldseek': None,
                        'pident_foldseek': None
                    }

        files.dict_to_json(structure_dir, 'human_foldseek_dict.json', results_foldseek_dict)
    else:
        results_foldseek_dict = files.json_to_dict(os.path.join(structure_dir, 'human_foldseek_dict.json'))
    
    return results_foldseek_dict

def merge_foldseek_data (base_path, organism_name, id_equivalences, uniprot_proteome_annotations):
    """
    Merge all the foldseek results in a single dictionary. For each locus_tag, the dictionary contains the human structure (PDB or AlphaFold) with the highest probability.
    Saves the dictionary in a .json file named using the organism name followed by '_final_foldseek_results.json' in the 'structures' directory.
    Returns a dictionary with the merged data.

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.
    :param id_equivalences: Dictionary with locus_tag and uniprot_id.
    :param uniprot_proteome_annotations: Dictionary with annotations for the UniProt proteome ID.

    :return: Dictionary with the merged data.
    """

    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    foldseek_res_file = os.path.join(structure_dir, 'human_foldseek_dict.json')
    foldseek_mapped_file = os.path.join(structure_dir, f'{organism_name}_final_foldseek_results.json')

    if not files.file_check(foldseek_mapped_file):
        if files.file_check(foldseek_res_file):  

            results_foldseek_dict = files.json_to_dict(foldseek_res_file)

            mapped_dict = {}

            # Parse the dictionaries and map the data
            for locus_tag, uniprot_ids in id_equivalences.items():
                mapped_dict[locus_tag] = {  'gene': locus_tag,
                                            'structure': None,
                                            'target': None,
                                            'rmsd': None,
                                            'prob': 0,
                                            'pident': 0
                                            }

                for uniprot_id in uniprot_ids:

                    proteome_data = uniprot_proteome_annotations.get(uniprot_id, {})

                    pdb_ids = []
                    pdb_ids_dict = proteome_data.get('PDB_id', [])
                    if pdb_ids_dict:
                        for entry in pdb_ids_dict:
                            pdb_ids.append(entry['ID'])
                        pdb_ids = list(set(pdb_ids))

                    if pdb_ids:
                        for pdb_id in pdb_ids:
                            if pdb_id in results_foldseek_dict:
                                if results_foldseek_dict[pdb_id]['prob_foldseek']:
                                    if mapped_dict[locus_tag]['prob'] < results_foldseek_dict[pdb_id]['prob_foldseek']:
                                        mapped_dict[locus_tag] = {
                                                    'gene': locus_tag,
                                                    'structure': pdb_id,
                                                    'target': results_foldseek_dict[pdb_id]['target_foldseek'],
                                                    'rmsd': results_foldseek_dict[pdb_id]['rmsd_foldseek'],
                                                    'prob': results_foldseek_dict[pdb_id]['prob_foldseek'],
                                                    'pident': results_foldseek_dict[pdb_id]['pident_foldseek']
                                                }
                                    elif mapped_dict[locus_tag]['prob'] == results_foldseek_dict[pdb_id]['prob_foldseek'] and mapped_dict[locus_tag]['pident'] < results_foldseek_dict[pdb_id]['pident_foldseek']:
                                        mapped_dict[locus_tag] = {
                                                    'gene': locus_tag,
                                                    'structure': pdb_id,
                                                    'target': results_foldseek_dict[pdb_id]['target_foldseek'],
                                                    'rmsd': results_foldseek_dict[pdb_id]['rmsd_foldseek'],
                                                    'prob': results_foldseek_dict[pdb_id]['prob_foldseek'],
                                                    'pident': results_foldseek_dict[pdb_id]['pident_foldseek']
                                                }
                        
                    if uniprot_id in results_foldseek_dict:
                        if results_foldseek_dict[uniprot_id]['prob_foldseek']:
                            if mapped_dict[locus_tag]['prob'] < results_foldseek_dict[uniprot_id]['prob_foldseek']:
                                mapped_dict[locus_tag] = {
                                            'gene': locus_tag,
                                            'structure': uniprot_id,
                                            'target': results_foldseek_dict[uniprot_id]['target_foldseek'],
                                            'rmsd': results_foldseek_dict[uniprot_id]['rmsd_foldseek'],
                                            'prob': results_foldseek_dict[uniprot_id]['prob_foldseek'],
                                            'pident': results_foldseek_dict[uniprot_id]['pident_foldseek']
                                        }

            files.dict_to_json(structure_dir, f'{organism_name}_final_foldseek_results.json', mapped_dict)                    
            print(f'Foldseek results saved to {organism_name}_final_foldseek_results.json.')
        else:
            print(f'File {foldseek_res_file} not found.')
    else:
        mapped_dict = files.json_to_dict(foldseek_mapped_file)
        print(f'Foldseek results in {foldseek_mapped_file}.')

    return mapped_dict

def final_foldseek_structure_table (base_path, organism_name, mapped_dict):
    """
    Create a final table with the merge dictionary of the Foldseek results. 
    Saves the table in a .tsv file named using the organism name followed by '_final_foldseek_results.tsv'in the 'structures' directory.
    Returns a DataFrame with the final results.

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.
    :param mapped_dict: Dictionary with the merged data.

    :return: DataFrame with the final results.
    """
    
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')

    all_locus_tags = metadata.ref_gbk_locus(base_path, organism_name)
    
    # Create a DataFrame with final results
    final_foldseek_file = os.path.join(structure_dir, f'{organism_name}_final_foldseek_results.tsv')

    if not files.file_check(final_foldseek_file):

        rows = []

        for locus_tag in all_locus_tags:
            if locus_tag in mapped_dict:
                rows.append(mapped_dict[locus_tag])
            else:
                rows.append({'gene': locus_tag, 'structure': 'No hit', 'target': None, 'rmsd': None, 'prob': None, 'pident': None })

        final_foldseek_df = pd.DataFrame(rows).rename(columns={
            'structure': 'FS_organism_structure_query',
            'target': 'FS_human_structure_hit',
            'rmsd': 'FS_rmsd',
            'prob': 'FS_prob',
            'pident': 'FS_pident'
        })

        final_foldseek_df.to_csv(final_foldseek_file, sep='\t', index=False)
   
        print(f'Foldseek final results saved to {final_foldseek_file}.')
    
    else:
        print(f'Foldseek final results in {final_foldseek_file}.')
        final_foldseek_df = pd.read_csv(final_foldseek_file, sep='\t')

    return final_foldseek_df