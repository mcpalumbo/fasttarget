import os
import sys
import glob
import requests
import xml.etree.ElementTree as ET
import urllib
import pandas as pd
import json
import tqdm
import shutil
import multiprocessing
from ftscripts import programs, files, offtargets, metadata
import re
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed


def uniprot_protein_annotations(uniprot_id):
    """
    Return a dictionary with annotations for a UniProt ID: Dataset, Protein name, Refseq ID, Sequence, 
    Version, EC, Locus_tag and PDB/Alphafold IDs. 

    :param uniprot_id: UniProt ID.

    :return: Dictionary with annotations for the UniProt ID.
    """

    uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}?format=xml"

    response = requests.get(uniprot_url)
    if response.status_code == 200:
        uniprot_xml = response.content
        root = ET.fromstring(uniprot_xml)
        namespaces = {'up': 'http://uniprot.org/uniprot'}

        # Extract accession number
        accession = root.find('up:entry/up:accession', namespaces).text
        version = root.find('up:entry', namespaces).attrib['version']
        dataset = root.find('up:entry', namespaces).attrib['dataset']

        if accession == uniprot_id:
            result = {accession: {'Version': version, 'Dataset': dataset}}

            # Extract protein names

            recommended_name = root.findall('up:entry/up:protein/up:recommendedName/up:fullName', namespaces)
            submitted_name = root.findall('up:entry/up:protein/up:submittedName/up:fullName', namespaces)

            if recommended_name:
                result[accession]['Protein_name'] = [protein_name.text for protein_name in recommended_name][0]
            else:
                if submitted_name:
                    result[accession]['Protein_name'] = [protein_name.text for protein_name in submitted_name][0]
                else:
                    result[accession]['Protein_name'] = None

            # Extract gene names
            ordered_locus = root.findall('up:entry/up:gene/up:name[@type="ordered locus"]', namespaces)
            if ordered_locus:
                result[accession]['Locus_tag'] = [gene_name.text for gene_name in ordered_locus]
            else:
                result[accession]['Locus_tag'] = None

            # Find the catalytic activity comment
            catalytic_activity_comment = root.find('.//up:comment[@type="catalytic activity"]', namespaces)
            if catalytic_activity_comment is not None:
                # Find the reaction element within the catalytic activity comment
                reaction = catalytic_activity_comment.find('up:reaction', namespaces)
                if reaction is not None:
                    # Find the dbReference element with type 'EC'
                    ec_reference = reaction.find('up:dbReference[@type="EC"]', namespaces)
                    if ec_reference is not None:
                        ec_number = ec_reference.get('id')
                        result[accession]['EC_number'] = ec_number
                    else:
                        result[accession]['EC_number'] = None
                else:
                    result[accession]['EC_number'] = None
            else:
                result[accession]['EC_number'] = None

            # Find Refseq references
            refseq = root.findall('.//up:dbReference[@type="RefSeq"]', namespaces)
            if refseq:
                for refseq_id in refseq:
                    id = refseq_id.get('id')
                    if id.startswith('WP_'):
                        result[accession]['Refseq_ProtID_nd'] = id
                    else:
                        result[accession]['Refseq_ProtID'] = id
            else:
                result[accession]['Refseq_ProtID_nd'] = None
                result[accession]['Refseq_ProtID'] = None

            # Find Alphafold ids
            alphafold = root.findall('.//up:dbReference[@type="AlphaFoldDB"]', namespaces)
            if alphafold:
                for alphafold_id in alphafold:
                    id = alphafold_id.get('id')
                    result[accession]['AlphaFoldDB'] = id
            else:
                result[accession]['AlphaFoldDB'] = None

            # Find PDB references
            pdb_references = root.findall('.//up:dbReference[@type="PDB"]', namespaces)
            if pdb_references:
                pdb_list = []
                for pdb_ref in pdb_references:
                    pdb_id = pdb_ref.get('id')
                    properties = pdb_ref.findall('up:property', namespaces)
                    pdb_info = {"ID": pdb_id}
                    for prop in properties:
                        prop_type = prop.get('type')
                        prop_value = prop.get('value')
                        pdb_info[prop_type] = prop_value
                    pdb_list.append(pdb_info)
                result[accession]['PDB_id'] = pdb_list
            else:
                result[accession]['PDB_id'] = None

            # Extract sequence
            sequence_element = root.find('up:entry/up:sequence', namespaces)
            if sequence_element is not None:
                sequence = sequence_element.text
                result[accession]['Sequence'] = sequence
            else:
                result[accession]['Sequence'] = None
        else:
            result = None
    else:
        result = None
    
    return result

def uniprot_proteome_ids(proteome_id):
    """
    Get the list of UniProt IDs from a proteome ID.

    :param proteome_id: UniProt proteome ID.

    :return: List of UniProt IDs.
    """
    
    url = f"https://rest.uniprot.org/uniprotkb/stream?format=list&query=%28%28proteome%3A{proteome_id}%29%29"
    data = urllib.request.urlopen(url).read()
    string = data.decode('UTF-8')
    uniprot_list = list(string.split('\n'))
    uniprot_list = list(filter(None, uniprot_list))
    return uniprot_list

def uniprot_proteome_save_annotations (base_path, organism_name, proteome_id, cpus=multiprocessing.cpu_count()):
    """
    Return a dictionary with UniProt annotations for a proteome ID.
    Create a .json file with this data.
    
    :param base_path: Base path.
    :param organism_name: Name of organism.
    :param proteome_id: Uniprot proteome ID for the organism.
    :param cpus: Number of threads (CPUs) to use.

    :return: Dictionary with annotations for the UniProt proteome ID.
    """
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    proteome_json = os.path.join(structure_dir, f"{proteome_id}_ann.json")

    if not files.file_check(proteome_json):
        print(f'Getting data from UniProt proteome {proteome_id}.')

        uniprot_id_list = uniprot_proteome_ids(proteome_id)
        dict_proteome = {}

        def fetch_annotation(n):
            dict_protein = uniprot_protein_annotations(n)
            single_id = str(list(dict_protein.keys())[0])
            return single_id, dict_protein[single_id]

        with ThreadPoolExecutor(max_workers=cpus) as executor:  
            future_to_protein = {executor.submit(fetch_annotation, n): n for n in uniprot_id_list}

            for future in tqdm.tqdm(as_completed(future_to_protein), total=len(uniprot_id_list), desc="Processing", unit="protein"):
                single_id, protein_data = future.result()
                dict_proteome[single_id] = protein_data

        print(f'UniProt proteome {proteome_id} data retrieved.')
        files.dict_to_json(structure_dir, f"{proteome_id}_ann.json", dict_proteome)
        print(f'UniProt proteome data saved to "{proteome_json}".')

    else:
        print(f'UniProt proteome data in "{proteome_json}".')
        dict_proteome = files.json_to_dict(proteome_json)

    return dict_proteome

def uniprot_proteome_structure_ids (base_path, organism_name, proteome_id):

    """
    Parse a dictionary with annotations for a UniProt proteme ID, located in structures dir.
    Returns two lists with IDs of AlphaFold and PDB structures.
    
    :param base_path: Base path.
    :param organism_name: Name of organism.
    :param proteome_id: Uniprot proteome ID for the organism.

    :return: Lists with AlphaFold IDs.
    :return: Lists with PDB IDs.
    """

    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    file_json = os.path.join(structure_dir, f"{proteome_id}_ann.json")
    
    proteome_ann_dict = files.json_to_dict(file_json)

    alphafold_ids = []
    pdb_ids = []

    for uniprot_id in proteome_ann_dict.keys():

        if proteome_ann_dict[uniprot_id]['AlphaFoldDB']:
            alphafold_ids.append(proteome_ann_dict[uniprot_id]['AlphaFoldDB'])
        
        if proteome_ann_dict[uniprot_id]['PDB_id']:
            for pdb_entry in proteome_ann_dict[uniprot_id]['PDB_id']:
                pdb_ids.append(pdb_entry['ID'])

    alphafold_ids_uniq = list(set(alphafold_ids))
    pdb_ids_uniq = list(set(pdb_ids))

    return alphafold_ids_uniq, pdb_ids_uniq

def uniprot_proteome_sequences (base_path, organism_name, proteome_id):

    """
    Parse a dictionary with annotations for a UniProt proteme ID, located in structures dir.
    Returns a fasta file with sequences of the proteome_id.
    
    :param base_path: Base path.
    :param organism_name: Name of organism.
    :param proteome_id: Uniprot proteome ID for the organism.
    """

    print(f'Getting sequences from uniprot proteome {proteome_id}.')

    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    proteome_json = os.path.join(structure_dir, f"{proteome_id}_ann.json")
    proteome_fasta = os.path.join(structure_dir, f"{proteome_id}.faa")

    if not files.file_check(proteome_fasta):
        if files.file_check(proteome_json):
            proteome_ann_dict = files.json_to_dict(proteome_json)
            if proteome_ann_dict:
                with open(proteome_fasta, 'w') as fasta_file:
                    for protein_id, data in proteome_ann_dict.items():
                        sequence = data['Sequence']
                        fasta_file.write(f">{protein_id}\n")
                        for i in range(0, len(sequence), 60):
                            fasta_file.write(sequence[i:i+60] + '\n')
                print(f'Uniprot proteome {proteome_id} sequences saved to {proteome_fasta}.')
            else:
                print(f'{proteome_json} cannot be parsed.')
        else:
            print(f'Uniprot proteome data not found in "{proteome_json}.')
    else:
        print(f'Uniprot proteome {proteome_id} sequences in {proteome_fasta}.')
       
def uniprot_proteome_index (base_path, organism_name, proteome_id ):
    """
    Runs NCBI makeblastdb against uniprot proteome.
    This function uses the `run_makeblastdb` function from the `programs` module.
    Database is saved in the structures directory.

    :param base_path: Base path.
    :param organism_name: Name of organism.
    :param proteome_id: Uniprot proteome ID.    
    """
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    blast_output_path = os.path.join(structure_dir, f'{proteome_id}_blast.tsv')

    #Index
    proteome_fasta = os.path.join(structure_dir, f'{proteome_id}.faa')
    proteome_index_path = os.path.join(structure_dir, f'{proteome_id}')

    if not files.file_check(blast_output_path):
        print(f'Indexing uniprot proteome  {proteome_id}')
        try:
            programs.run_makeblastdb(
            input= proteome_fasta,
            output= proteome_index_path,
            title= f'{proteome_id}',
            dbtype= 'prot'
            )
            print(f'Index finished for uniprot proteome  {proteome_id}')
        except Exception as e:
            print(f"Failed to run makeblastdb to file {proteome_fasta}: {e}")
   
def uniprot_proteome_blast (base_path, organism_name, proteome_id, cpus=multiprocessing.cpu_count()):
    """
    Runs NCBI Blastp against uniprot proteome.
    This function uses the `run_blastp` function from the `programs` module. 
    Results are saved in a .tsv file in the structures directory, named using the Uniprot proteome ID followed by '_blast.tsv'.

    :param base_path: Base path.
    :param organism_name: Name of organism.
    :param proteome_id: Uniprot proteome ID.
    :param cpus: Number of threads (CPUs) to use in blast search.
    """

    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    proteome_index_path = os.path.join(structure_dir, f'{proteome_id}')

    organism_prot_seq_path = os.path.join(base_path, 'organism', f'{organism_name}','genome', f'{organism_name}.faa')
    blast_output_path = os.path.join(structure_dir, f'{proteome_id}_blast.tsv')

    if not files.file_check(blast_output_path):
        print(f'Runing blastp for {organism_name} and uniprot proteome {proteome_id}')

        try:
            programs.run_blastp(
                blastdb= proteome_index_path,
                query= organism_prot_seq_path,
                output=blast_output_path,
                evalue= '1e-5',
                outfmt= '6 std qcovhsp qcovs',
                max_target_seqs = '5',
                cpus=cpus
            )
            print(f'Blastp finished for uniprot proteome  {proteome_id}.')
            print(f'Blastp results saved in  {blast_output_path}.')

        except Exception as e:
            print(f"Failed to run blastp to file {proteome_index_path}: {e}")
    else:
        print(f'Blastp results in  {blast_output_path}.')

def uniprot_proteome_blast_parse (base_path, organism_name, proteome_id):

    """
    Parse NCBI blastp results against uniprot proteome.
    Returns a dictionary with each locus_tag and the identical hit in proteome. 
    This dictionary is saved in a .json file in the structures directory, named using the Uniprot proteome ID followed by '_ids.json'.

    :param base_path: Base path.
    :param organism_name: Name of organism.
    :param proteome_id: Uniprot proteome ID.

    :return: Dictionary with locus_tag and uniprot_id.
    """

    print(f'Reading blastp results for {proteome_id}.')

    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    blast_output_path = os.path.join(structure_dir, f'{proteome_id}_blast.tsv')
    proteome_ids_file = os.path.join(structure_dir, f'{proteome_id}_ids.json')

    if files.file_check(blast_output_path):
        if not files.file_check(proteome_ids_file):

            blast_output_df = offtargets.read_blast_output(blast_output_path)

            #Filter % identity and coverage - identical proteins
            filtered_df = blast_output_df[(blast_output_df['pident'] > 95) 
                                        & (blast_output_df['qcovs'] > 50)]

            id_equivalence = {}

            for qseqid, group in filtered_df.groupby('qseqid'):

                # Find the maximum identity value in the group
                max_identity = group['pident'].max()
                # Filter rows that have the maximum identity
                best_identity = group[group['pident'] == max_identity]

                # Find the maximum coverage value in the group
                max_coverage = best_identity['qcovs'].max()
                # Filter rows that have the maximum coverage
                best_results = best_identity[best_identity['qcovs'] == max_coverage]
        
                # Store all sseqid values with the maximum identity and coverage
                id_equivalence[qseqid] = best_results['sseqid'].tolist()
        
                if len(id_equivalence[qseqid]) > 1:
                    print(f'Locus_tag {qseqid} with several uniprot ids: {id_equivalence[qseqid]}.')

            files.dict_to_json(structure_dir, f'{proteome_id}_ids.json', id_equivalence)
            print(f'IDs file saved in {proteome_ids_file}.')

        else:
            id_equivalence = files.json_to_dict(proteome_ids_file)
            print(f'IDs file in {proteome_ids_file}.')
    
        return id_equivalence
    
    else:
        print(f'{blast_output_path} not found.')   

def get_structure_PDB (output_path, PDB_id):
    """
    Download structure from PDB. Returns True if successful.
    
    :param output_path: Output directory path.
    :param PDB_id: ID of PDB structure.

    :return: True if successful.
    """
    res = False
    file_path = os.path.join(output_path, f"PDB_{PDB_id}.pdb")

    pdb_url = f"https://files.rcsb.org/download/{PDB_id}.pdb"

    max_retries = 0

    while max_retries < 3:
        try:
            response = requests.get(pdb_url)
            if response.status_code == 200:
                with open(file_path, 'wb') as file:
                    file.write(response.content)
                print(f"Download {PDB_id} finished.")
                max_retries += 3
                res = True
            else:
                max_retries += 1
                if max_retries == 3:
                    print(f"Failed to download {PDB_id}.")
        except requests.exceptions.RequestException as e:
            print(f'Failed to download {PDB_id}, attempt number {max_retries}: {e}')
            max_retries += 1
        
    return res

def get_structure_alphafold(output_path, uniprot_id):
    """
    Download structure from AlphaFold. Returns True if successful.
    
    :param output_path: Output directory path.
    :param uniprot_id: ID of UniProt with AlphaFold structure.

    :return: True if successful.
    """

    res = False
    file_path = os.path.join(output_path, f"AF_{uniprot_id}.pdb")

    alphafold_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"

    max_retries = 0

    while max_retries < 3:
        try: 
            response = requests.get(alphafold_url)
            if response.status_code == 200:
                with open(file_path, 'wb') as file:
                    file.write(response.content)
                print(f"Download AlphaFold model {uniprot_id} finished.")
                max_retries += 3
                res = True
            else:
                max_retries += 1
                if max_retries == 3:
                    print(f"Failed to download AlphaFold prediction for {uniprot_id}.")
        except requests.exceptions.RequestException as e:
            print(f"Failed to download AlphaFold prediction for {uniprot_id}, attempt number {max_retries}: {e}")
            max_retries += 1
   
    return res

def get_structures (base_path, organism_name, proteome_id, cpus=multiprocessing.cpu_count()):
    """
    Download PDB and AlphaFold structures from a list of IDs. 
    This function creates two directories in the 'structures' directory: 'PDB_files' and 'AlphaFold_files'.
    Returns lists of failed IDs.
    
    :param base_path: Base path.
    :param organism_name: Name of organism.
    :param proteome_id: Uniprot proteome ID for the organism.
    :param cpus: Number of threads (CPUs) to use.

    :return: Lists of failed PDB IDs.
    :return: Lists of failed AlphaFold IDs.
    """

    alphafold_ids, pdb_ids = uniprot_proteome_structure_ids(base_path, organism_name, proteome_id)
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')

    if os.path.exists(structure_dir):
        output_dir_PDB = os.path.join(structure_dir, 'PDB_files')
        os.makedirs(output_dir_PDB, exist_ok=True)

        output_dir_AF = os.path.join(structure_dir, 'AlphaFold_files')
        os.makedirs(output_dir_AF, exist_ok=True)

        failed_pdb = []
        failed_AF = []

        def download_pdb(pdb_id):
            """
            Download PDB structure. Returns PDB ID if failed.
            
            :param pdb_id: ID of PDB structure.
            
            :return: PDB ID if failed.
            """

            file_path_pdb = os.path.join(output_dir_PDB, f'PDB_{pdb_id}.pdb')
            if not files.file_check(file_path_pdb):
                pdb_res = get_structure_PDB(output_dir_PDB, pdb_id)
                if not pdb_res:
                    return pdb_id
            return None

        def download_af(af_id):
            """
            Download AlphaFold structure. Returns AlphaFold ID if failed.

            :param af_id: ID of AlphaFold structure.

            :return: AlphaFold ID if failed.
            """

            file_path_af = os.path.join(output_dir_AF, f'AF_{af_id}.pdb')
            if not files.file_check(file_path_af):
                af_res = get_structure_alphafold(output_dir_AF, af_id)
                if not af_res:
                    return af_id
            return None

        # Parallel download of PDB structures
        print('----Downloading PDB structures----')
        with ThreadPoolExecutor(max_workers=cpus) as executor:
            futures_pdb = {executor.submit(download_pdb, pdb_id): pdb_id for pdb_id in pdb_ids}
            for future in tqdm.tqdm(as_completed(futures_pdb), total=len(pdb_ids), desc="Downloading from PDB", unit="model"):
                result = future.result()
                if result:
                    failed_pdb.append(result)
        print('----Download finished!----')

        fail_file_pdb = os.path.join(output_dir_PDB, 'Failed_PDB_ids.txt')
        files.list_to_file(fail_file_pdb, failed_pdb)
        print(f"The IDs that could not be downloaded are saved in {fail_file_pdb}")

        # Parallel download of AlphaFold structures
        print('----Downloading AlphaFold structures----')
        with ThreadPoolExecutor(max_workers=cpus) as executor:
            futures_af = {executor.submit(download_af, af_id): af_id for af_id in alphafold_ids}
            for future in tqdm.tqdm(as_completed(futures_af), total=len(alphafold_ids), desc="Downloading from AlphaFold", unit="model"):
                result = future.result()
                if result:
                    failed_AF.append(result)
        print('----Download finished!----')

        fail_file_af = os.path.join(output_dir_AF, 'Failed_AlphaFold_ids.txt')
        files.list_to_file(fail_file_af, failed_AF)
        print(f"The IDs that could not be downloaded are saved in {fail_file_af}")

        return failed_pdb, failed_AF

    else:
        print(f"The directory '{structure_dir}' not found.", file=sys.stderr)

def check_complete_structures (base_path, organism_name, proteome_id):

    """
    Check if all PDB and AlphaFold structures from list of IDs has been downloaded or are in the list of failed IDs. 
    Returns True if the download was successful.
    
    :param base_path: Base path.
    :param organism_name: Name of organism.
    :param proteome_id: Uniprot proteome ID for the organism.

    :return: True if the download was successful.
    """

    alphafold_ids, pdb_ids = uniprot_proteome_structure_ids (base_path, organism_name, proteome_id)
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')

    res = False
    
    if os.path.exists(structure_dir):

        #PDB
        output_dir_PDB = os.path.join(structure_dir,'PDB_files')
        fail_file_pdb = os.path.join(output_dir_PDB,'Failed_PDB_ids.txt')
        pdb_files = glob.glob(os.path.join(output_dir_PDB, '*.pdb'))
        
        #AlphaFold
        output_dir_AF = os.path.join(structure_dir,'AlphaFold_files')
        fail_file_af = os.path.join(output_dir_AF,'Failed_AlphaFold_ids.txt')
        AF_files = glob.glob(os.path.join(output_dir_AF, '*.pdb'))

        if os.path.isfile(fail_file_pdb) and os.path.isfile(fail_file_af):              
            failed_pdb =  files.file_to_list(fail_file_pdb)
            failed_AF = files.file_to_list(fail_file_af)

            PDB_complete = len(pdb_ids) == len(pdb_files) + len(failed_pdb)
            AF_complete = len(alphafold_ids) == len(AF_files) + len(failed_AF)

            if PDB_complete and AF_complete:
                res = True
                print(f'All structures downloaded correctly in {output_dir_PDB} and {output_dir_AF}.')
    else:
        print(f"The directory '{structure_dir}' not found.", file=sys.stderr)
        
    return res

def fpocket_for_structure(pdb, output_path, fpocket_results):
    """
    Run Fpocket for a PDB structure.
    This function uses the `run_fpocket` function from the `programs` module.
    Results are saved in a directory named 'fpocket_results' in the output path.

    :param pdb: Path to PDB structure.
    :param output_path: Output directory path.
    :param fpocket_results: Directory to save fpocket results.
    """

    fpocket_outdir = os.path.basename(os.path.splitext(pdb)[0]) + "_out"
    results_path = os.path.join(fpocket_results, fpocket_outdir)

    if not os.path.exists(results_path):
        print(f'Running Fpocket for {pdb}')
        programs.run_fpocket(output_path, pdb)  # Run the fpocket Docker command
        shutil.move(os.path.join(output_path, fpocket_outdir), fpocket_results)
        print(f'FPocket prediction for {pdb} completed.')
    print(f'{pdb} already processed.')

def pockets_finder (output_path, cpus=multiprocessing.cpu_count()):
    '''
    Run Fpocket for all PDB structures in a directory.
    Results are saved in a directory named 'fpocket_results' in the output path.

    :param output_path:
    :param cpus: Number of threads (CPUs) to use
    '''

    if os.path.exists(output_path):
        fpocket_results = os.path.join(output_path, 'fpocket_results')
        os.makedirs(fpocket_results, exist_ok=True)

        pdb_files = glob.glob(os.path.join(output_path, '*.pdb'))

        # Use ProcessPoolExecutor for parallel processing of CPU-bound tasks
        with ProcessPoolExecutor(max_workers=cpus) as executor:
            futures = {executor.submit(fpocket_for_structure, pdb, output_path, fpocket_results): pdb for pdb in pdb_files}

            for future in tqdm.tqdm(as_completed(futures), total=len(pdb_files), desc="Processing", unit="structure"):
                result = future.result()

    else:
        print(f"The directory '{output_path}' not found.", file=sys.stderr)

def pockets_parse_output(content):
    """
    Parse the output of Fpocket.
    Returns a dictionary with the pockets data.

    :param content: Content of the fpocket output file.

    :return: Dictionary with the data
    """

    data = {}
    current_pocket = None

    pocket_re = re.compile(r"^Pocket (\d+) :")
    key_value_re = re.compile(r"\t(.+?) :\s+(.+)")

    for line in content.split('\n'):
        pocket_match = pocket_re.match(line)
        key_value_match = key_value_re.match(line)

        if pocket_match:
            current_pocket = f"Pocket {pocket_match.group(1)}"
            data[current_pocket] = {}
        elif key_value_match and current_pocket:
            key = key_value_match.group(1).strip()
            value = key_value_match.group(2).strip()
            try:
                if '.' in value:
                    value = float(value)
                else:
                    value = int(value)
            except ValueError:
                pass
            data[current_pocket][key] = value

    return data

def pockets_data_to_dict(directory):
    """
    Parse the output of Fpocket for all .pdb structures in a directory.
    Saves the data in a file named 'all_pockets.json'.
    Returns a dictionary with the pockets data.

    :param directory: Directory with the Fpocket results for each structure.

    :return: Dictionary with the pockets data.
    """

    if os.path.exists(directory):

        pockets_file = os.path.join(directory, 'all_pockets.json')

        if not files.file_check(pockets_file):
            print('Parsing pockets results.')
            pockets_dict = {}
            
            for root, dirs, list_files in os.walk(directory):
                for file in list_files:
                    if file.endswith('_info.txt'):
                        file_path = os.path.join(root, file)
                        protein_ID = file.split("_")[1]
                        with open(file_path, 'r') as f:
                            content = f.read()
                            if content:
                                pockets_dict[protein_ID] = pockets_parse_output(content)
                            else:
                                pockets_dict[protein_ID] = 'No_pockets'

            files.dict_to_json(directory, 'all_pockets.json', pockets_dict)
            print(f'All pockets data saved to {pockets_file}.')
        
        else:
            pockets_dict = files.json_to_dict(pockets_file)
            print(f'All pockets data in {pockets_file}.')
        
        return pockets_dict
    
    else:
        print(f'{directory} not found.')

def pockets_filter(pockets_dict, out_dir):
    """
    Filter the pockets data to keep only the one with the highest Druggability Score.
    Saves the filtered data in a file named 'best_pockets_DS.json'.
    Returns a dictionary with the filtered data.

    :param pockets_dict: Dictionary with the pockets data.
    :param out_dir: Output directory path.

    :return: Dictionary with the best pockets data.
    """

    best_pockets_file = os.path.join(out_dir, 'best_pockets_DS.json')
    
    if not files.file_check(best_pockets_file):
        
        if not pockets_dict or not isinstance(pockets_dict, dict):
            raise ValueError("pockets_dict must be a dictionary with the pockets data.")

        print('Obtaining best druggability score for each structure.')

        DS_dict = {}
        for prots, pockets in pockets_dict.items():
            if pockets_dict[prots] == 'No_pockets':
                DS_dict[prots] = 'No_pockets'
            else:
                n = 0
                max_pocket = None
                for pocket, props in pockets.items():
                    for prop, value in props.items():
                        if prop == 'Druggability Score':
                            if value > n:
                                n = value
                                max_pocket = pocket

                DS_dict[prots] = {}
                DS_dict[prots]['maxDS'] = n
                DS_dict[prots]['pocket'] = max_pocket

        files.dict_to_json(out_dir, 'best_pockets_DS.json', DS_dict)
        print(f'Best pockets data saved to {best_pockets_file}.')

    else:
        DS_dict = files.json_to_dict(best_pockets_file)
        print(f'Best pockets data in {best_pockets_file}.')

    return DS_dict

def merge_structure_data (base_path, organism_name, id_equivalences, uniprot_proteome_annotations, best_pockets_dict_PDB, best_pockets_dict_AF):
    """
    Merge all the structure data in a single dictionary.
    Saves the dictionary in a .json file named using the organism name followed by '_structure_data.json' in the 'structures' directory.
    Returns a dictionary with the merged data.

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.
    :param id_equivalences: Dictionary with locus_tag and uniprot_id.
    :param uniprot_proteome_annotations: Dictionary with annotations for the UniProt proteome ID.
    :param best_pockets_dict_PDB: Dictionary with the best pockets (highest DS) data for PDB structures.
    :param best_pockets_dict_AF: Dictionary with the best pockets data (highest DS) for AlphaFold structures.

    :return: Dictionary with the merged data.
    """

    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    all_struct_data_file = os.path.join(structure_dir, f'{organism_name}_structure_data.json')
    
    if not files.file_check(all_struct_data_file):
        
        pdbs_not_found = []
        af_not_found = []

        mapped_dict = {}

        # Parse the dictionaries and map the data
        for locus_tag, uniprot_ids in id_equivalences.items():
            mapped_dict[locus_tag] = {}

            for uniprot_id in uniprot_ids:
                mapped_dict[locus_tag][uniprot_id] = {'PDB': {}, 'AF': {}}

                # Process PDB data
                proteome_data = uniprot_proteome_annotations.get(uniprot_id, {})
                pdb_ids = proteome_data.get('PDB_id', [])
                

                if pdb_ids:
                    for pdb_entry in pdb_ids:
                        pdb_id = pdb_entry['ID']
                        pdb_info = best_pockets_dict_PDB.get(pdb_id, 'No_PDB')
                        if pdb_info == 'No_PDB':
                            print(f'Pockets results not found for PDB {pdb_id}.')
                            pdbs_not_found.append(pdb_id)
                        elif pdb_info == 'No_pockets':
                            mapped_dict[locus_tag][uniprot_id]['PDB'][pdb_id] = 'No_pockets'
                        else:
                            mapped_dict[locus_tag][uniprot_id]['PDB'][pdb_id] = {
                                'maxDS': pdb_info['maxDS'],
                                'pocket': pdb_info['pocket']
                            }

                # Process AlphaFold (AF) data
                af_info = best_pockets_dict_AF.get(uniprot_id, 'No_AF')
                if af_info == 'No_AF':
                    print(f'Pockets results not found for AlphaFold {uniprot_id}.')
                    af_not_found.append(uniprot_id)
                elif af_info == 'No_pockets':
                    mapped_dict[locus_tag][uniprot_id]['AF'] = 'No_pockets'
                else:
                    mapped_dict[locus_tag][uniprot_id]['AF'] = {
                        'maxDS': af_info['maxDS'],
                        'pocket': af_info['pocket']
                }

        print(f'Pocket results not found for {len(set(pdbs_not_found))} PDB structures: {set(pdbs_not_found)}')
        print(f'Check if .pdb files were not downloaded in  {structure_dir}/PDB_files/Failed_PDB_ids.txt')

        print(f'Pocket results not found for {len(set(af_not_found))} AlphaFold structures: {set(af_not_found)}')
        print(f'If a UniProt ID has an AF model, check whether it has been downloaded in {structure_dir}/AlphaFold_files/Failed_AlphaFold_ids.txt')
              
        # Save the parsed data
        files.dict_to_json(structure_dir, f'{organism_name}_structure_data.json', mapped_dict)
        print(f'All pockets results saved to {all_struct_data_file}.')

    else:
        print(f'All pockets results in {all_struct_data_file}.')
        mapped_dict = files.json_to_dict(all_struct_data_file)

    return mapped_dict

def final_structure_table (base_path, organism_name, mapped_dict):
    """
    Create a final table with the structure data. 
    Saves the table in a .tsv file named using the organism name followed by '_structure_results.tsv'in the 'structures' directory.
    Returns a DataFrame with the final results.

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.
    :param mapped_dict: Dictionary with the merged data.

    :return: DataFrame with the final results.
    """
    
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')

    all_locus_tags = metadata.ref_gbk_locus(base_path, organism_name)

    # Create a DataFrame with final results
    final_structure_file = os.path.join(structure_dir, f'{organism_name}_structure_results.tsv')

    if not files.file_check(final_structure_file):

        rows = []

        for locus_tag in all_locus_tags:
            if locus_tag in mapped_dict:
                max_DS = 0
                structure_max_DS = 'No_structure'
                best_pocket = 'No_structure'
                final_uniprot_id = ''

                for uniprot_id, structures in mapped_dict[locus_tag].items():
                    if structures['PDB']:
                        for pdb_id, pdb_data in structures['PDB'].items():
                            if pdb_data == 'No_pockets':
                                if best_pocket == 'No_structure':
                                    final_uniprot_id = uniprot_id
                                    structure_max_DS = f'PDB_{pdb_id}'
                                    best_pocket = 'No_pockets'
                            else:
                                if pdb_data['maxDS'] > max_DS or 'PDB' not in structure_max_DS:
                                    max_DS = pdb_data['maxDS']
                                    best_pocket = pdb_data['pocket']
                                    structure_max_DS = f'PDB_{pdb_id}'
                                    final_uniprot_id = uniprot_id
                    elif structures['AF']:
                        if structures['AF'] == 'No_pockets':
                            if best_pocket == 'No_structure':
                                final_uniprot_id = uniprot_id
                                structure_max_DS = f'AF_{uniprot_id}'
                                best_pocket = 'No_pockets'
                        else:
                            if structures['AF']['maxDS'] > max_DS and 'PDB' not in structure_max_DS:
                                max_DS = structures['AF']['maxDS']
                                best_pocket = structures['AF']['pocket']
                                structure_max_DS = f'AF_{uniprot_id}'
                                final_uniprot_id = uniprot_id
                    else:
                        if final_uniprot_id == '':
                            final_uniprot_id = uniprot_id

                rows.append([locus_tag, final_uniprot_id, max_DS if max_DS > 0 else None, best_pocket, structure_max_DS])
            else:
                rows.append([locus_tag, None, None, 'No_structure', 'No_structure'])

        final_pockets_df = pd.DataFrame(rows, columns=['gene', 'uniprot', 'druggability_score', 'pocket', 'structure'])

        final_pockets_df.to_csv(final_structure_file, sep='\t', index=False)
   
        print(f'Structure final results saved to {final_structure_file}.')

        return final_pockets_df
    
    else:
        print(f'Structure final results in {final_structure_file}.')     

def uniprot_proteome(base_path, organism_name, proteome_id, cpus=multiprocessing.cpu_count()):

    """
    First step to process the structure data.
    Gets annotation data from a UniProt proteome ID and saves it in a .json file in the structures directory.
    Maps each locus_tag to the equivalent uniprot_id and saves it in a .json file in the structures directory.
    Returns the annotation dictionary and the id equivalences dictionary.

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.
    :param proteome_id: Uniprot proteome ID for the organism.
    :param cpus: Number of threads (CPUs) to use.
    
    :return: Dictionary with annotations for the UniProt proteome ID.
    :return: Dictionary with equivalent locus_tag and uniprot_id.
    """

    print(f'---------- 1. Processing uniprot proteome {proteome_id} data ----------')
    
    uniprot_proteome_annotations = uniprot_proteome_save_annotations (base_path, organism_name, proteome_id, cpus)

    uniprot_proteome_sequences (base_path, organism_name, proteome_id)
    

    uniprot_proteome_index(base_path, organism_name, proteome_id)
    uniprot_proteome_blast(base_path, organism_name, proteome_id, cpus)
    
    id_equivalences = uniprot_proteome_blast_parse(base_path, organism_name, proteome_id)

    print(f'---------- 1. Finished ----------')

    return uniprot_proteome_annotations, id_equivalences

def structures (base_path, organism_name, proteome_id, cpus=multiprocessing.cpu_count()):
    """
    Second step to process the structure data.
    Downloads PDB and AlphaFold structures from a UniProt proteome ID.
    Saves the structures in the 'PDB_files' and 'AlphaFold_files' directories in the 'structures' directory.

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.
    :param proteome_id: Uniprot proteome ID for the organism.
    :param cpus: Number of threads (CPUs) to use.
    """

    print(f'---------- 2. Obtaining  uniprot proteome {proteome_id} structures ----------')

    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    output_dir_PDB = os.path.join(structure_dir,'PDB_files')
    output_dir_AF = os.path.join(structure_dir,'AlphaFold_files')

    status = check_complete_structures (base_path, organism_name, proteome_id)
    while not status:
        get_structures (base_path, organism_name, proteome_id, cpus)
        status = check_complete_structures (base_path, organism_name, proteome_id)
   
    print(f'All PDB structures are found in {output_dir_PDB}')
    print(f'All AlphaFold structures are found in {output_dir_AF}')

    print(f'---------- 2. Finished ----------')

def pockets (base_path, organism_name, id_equivalences, uniprot_proteome_annotations, cpus=multiprocessing.cpu_count()):
    """
    Third step to process the structure data.
    Runs Fpocket for PDB and AlphaFold structures.
    Parses the Fpocket results and filters the data to keep only the pockets with the highest Druggability Score.
    Merges all the structure data in a single dictionary.
    Creates a final table with the structure data.

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.
    :param id_equivalences: Dictionary with equivalent locus_tag and uniprot_id.
    :param uniprot_proteome_annotations: Dictionary with annotations for the UniProt proteome ID.
    :param cpus: Number of threads (CPUs) to use.

    :return: DataFrame with the final structure results.
    """

    print(f'---------- 3. Obtaining pockets data ----------')

    # Directories
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    output_dir_PDB = os.path.join(structure_dir,'PDB_files')
    output_dir_AF = os.path.join(structure_dir,'AlphaFold_files')
    fpocket_results_PDB = os.path.join(output_dir_PDB,'fpocket_results')
    fpocket_results_AF = os.path.join(output_dir_AF,'fpocket_results')
    
    print(f'----- 3.1. PDB structures -----')
    pockets_finder (output_dir_PDB, cpus)
    pockets_dict_PDB = pockets_data_to_dict(fpocket_results_PDB)
    best_pockets_dict_PDB = pockets_filter(pockets_dict_PDB, output_dir_PDB)


    print(f'----- 3.2. AlphaFold models -----')
    pockets_finder (output_dir_AF, cpus)  
    pockets_dict_AF = pockets_data_to_dict(fpocket_results_AF)    
    best_pockets_dict_AF = pockets_filter(pockets_dict_AF, output_dir_AF)

    print(f'----- 3.3. Processing results -----')
    mapped_dict = merge_structure_data (base_path, organism_name, id_equivalences, uniprot_proteome_annotations, best_pockets_dict_PDB, best_pockets_dict_AF)
    final_df = final_structure_table (base_path, organism_name, mapped_dict)

    print(f'---------- 3. Finished ----------')

    return final_df

def FPocket_models(directory):
    """
    Run Fpocket for all .pdb in a directory.
    Results are saved inside each fpocket output directory as all_pockets.json.

    :param directory: Directory with the models.
    """

    if os.path.exists(directory):
        for protein in tqdm.tqdm(os.listdir(directory), desc="Processing", unit="protein"):
            protein_path = os.path.join(directory, protein)
            for template in os.listdir(protein_path):
                template_path = os.path.join(protein_path, template)
                for model_file in os.listdir(template_path):
                    model_pdb_path = os.path.join(template_path, model_file)
                    if os.path.isfile(model_pdb_path) and model_file.endswith('.pdb'):
                        fpocket_outdir = os.path.splitext(model_pdb_path)[0] + "_out"
                        if not os.path.exists(fpocket_outdir):
                            print(f'Run Fpocket for {model_pdb_path}')
                            programs.run_fpocket(template_path, model_pdb_path)
                            print(f'FPocket prediction for {model_file} in {fpocket_outdir}.')
                            pockets_data_to_dict(fpocket_outdir)
                            pockets_filter(f'{fpocket_outdir}/all_pockets.json')
    else:
        print(f"The directory '{directory}' not found.", file=sys.stderr)