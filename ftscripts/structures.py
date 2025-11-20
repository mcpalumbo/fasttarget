import os
import sys
import urllib.request
from tqdm import tqdm
import pandas as pd
import re
import tempfile
from Bio import SeqIO
import multiprocessing
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select
from ftscripts import programs, files, metadata, logger
import databases
import glob
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
import requests
import xml.etree.ElementTree as ET
import time
from functools import wraps  
from urllib3.exceptions import SSLError
from requests.exceptions import (
    SSLError as RequestsSSLError, 
    ConnectionError, 
    Timeout
)


## ------------------- UNIPROT PROTEOME FUNCTIONS ------------------- ##

def parse_pdb_refs(pdb_string):
    """
    Parse PDB references from metadata string
    E.g., "1ABC;2DEF;3GHI" -> ["1ABC", "2DEF", "3GHI"]
    """
    if pd.isna(pdb_string) or pdb_string == "":
        return []
    else:
        return [pdb.strip() for pdb in str(pdb_string).split(';') if pdb.strip()]



@databases.retry_with_backoff(max_retries=10, initial_delay=2, backoff_factor=2, max_delay=80)
def fetch_uniprot_xml(uniprot_id):
    """
    Fetch UniProt XML data with automatic retry on network/SSL errors.
    
    Uses the retry_with_backoff decorator from databases module for robust
    error handling with exponential backoff.
    
    :param uniprot_id: UniProt accession ID.
    :return: XML content as bytes.
    :raises: requests.exceptions.RequestException if all retries fail.
    """
    uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}?format=xml"
    
    # Add User-Agent header to identify the tool
    headers = {
        'User-Agent': 'FastTarget/1.0',
        'Accept': 'application/xml'
    }
    
    response = requests.get(uniprot_url, headers=headers, timeout=30)
    response.raise_for_status()  # Raise exception for bad status codes
    
    return response.content

@databases.retry_with_backoff(max_retries=10, initial_delay=2, backoff_factor=2, max_delay=80)
def fetch_uniprot_batch_xml(uniprot_ids):
    """
    Fetch UniProt XML data for MULTIPLE UniProt IDs in a single batch request.
    
    Uses UniProt's batch endpoint: /uniprotkb/accessions?accessions=ID1,ID2,ID3
    Returns a single XML document containing multiple <entry> elements.
    
    :param uniprot_ids: List of UniProt accession IDs (max 500 recommended).
    :return: XML content as bytes containing all entries.
    :raises: requests.exceptions.RequestException if all retries fail.
    """
    if not uniprot_ids:
        return b'<?xml version="1.0" encoding="UTF-8"?><uniprot xmlns="http://uniprot.org/uniprot"></uniprot>'
    
    # Join IDs with commas (UniProt batch endpoint format)
    ids_string = ",".join(uniprot_ids)
    uniprot_url = f"https://rest.uniprot.org/uniprotkb/accessions?accessions={ids_string}&format=xml"
    
    headers = {
        'User-Agent': 'FastTarget/1.0',
        'Accept': 'application/xml'
    }
    
    logger.logger.info(f"Fetching batch of {len(uniprot_ids)} UniProt entries...")
    response = requests.get(uniprot_url, headers=headers, timeout=60)
    response.raise_for_status()
    
    return response.content

def parse_uniprot_entry_xml(entry_element):
    """
    Parse a single UniProt <entry> XML element into annotation dictionary.
    
    Extracted from uniprot_protein_annotations() to allow reuse for both
    single and batch queries.
    
    :param entry_element: ElementTree Element representing a single <entry>.
    :return: Dictionary with annotations {uniprot_id: {annotation_data}}, or None if parsing fails.
    """
    try:
        namespaces = {'up': 'http://uniprot.org/uniprot'}
        
        # Extract accession number
        accession_elem = entry_element.find('up:accession', namespaces)
        if accession_elem is None:
            logger.logger.warning("No accession found in entry")
            return None
        
        accession = accession_elem.text
        version = entry_element.attrib.get('version')
        dataset = entry_element.attrib.get('dataset')
        
        result = {accession: {'Version': version, 'Dataset': dataset}}
        
        # Extract protein names
        recommended_name = entry_element.findall('up:protein/up:recommendedName/up:fullName', namespaces)
        submitted_name = entry_element.findall('up:protein/up:submittedName/up:fullName', namespaces)
        
        if recommended_name:
            result[accession]['Protein_name'] = recommended_name[0].text
        elif submitted_name:
            result[accession]['Protein_name'] = submitted_name[0].text
        else:
            result[accession]['Protein_name'] = None
        
        # Extract gene names (locus tags)
        ordered_locus = entry_element.findall('up:gene/up:name[@type="ordered locus"]', namespaces)
        if ordered_locus:
            result[accession]['Locus_tag'] = [gene_name.text for gene_name in ordered_locus]
        else:
            result[accession]['Locus_tag'] = None
        
        # Extract EC number
        catalytic_activity_comment = entry_element.find('.//up:comment[@type="catalytic activity"]', namespaces)
        if catalytic_activity_comment is not None:
            reaction = catalytic_activity_comment.find('up:reaction', namespaces)
            if reaction is not None:
                ec_reference = reaction.find('up:dbReference[@type="EC"]', namespaces)
                result[accession]['EC_number'] = ec_reference.get('id') if ec_reference is not None else None
            else:
                result[accession]['EC_number'] = None
        else:
            result[accession]['EC_number'] = None
        
        # Extract RefSeq IDs
        refseq = entry_element.findall('.//up:dbReference[@type="RefSeq"]', namespaces)
        result[accession]['Refseq_ProtID_nd'] = None
        result[accession]['Refseq_ProtID'] = None
        if refseq:
            for refseq_id in refseq:
                id_val = refseq_id.get('id')
                if id_val.startswith('WP_'):
                    result[accession]['Refseq_ProtID_nd'] = id_val
                else:
                    result[accession]['Refseq_ProtID'] = id_val
        
        # Extract AlphaFold IDs
        alphafold = entry_element.findall('.//up:dbReference[@type="AlphaFoldDB"]', namespaces)
        result[accession]['AlphaFoldDB'] = alphafold[0].get('id') if alphafold else None
        
        # Extract PDB references
        pdb_references = entry_element.findall('.//up:dbReference[@type="PDB"]', namespaces)
        if pdb_references:
            pdb_list = []
            for pdb_ref in pdb_references:
                pdb_id = pdb_ref.get('id')
                properties = pdb_ref.findall('up:property', namespaces)
                pdb_info = {"ID": pdb_id}
                for prop in properties:
                    pdb_info[prop.get('type')] = prop.get('value')
                pdb_list.append(pdb_info)
            result[accession]['PDB_id'] = pdb_list
        else:
            result[accession]['PDB_id'] = None
        
        # Extract sequence
        sequence_element = entry_element.find('up:sequence', namespaces)
        result[accession]['Sequence'] = sequence_element.text if sequence_element is not None else None
        
        return result
        
    except Exception as e:
        logger.logger.error(f"Failed to parse UniProt entry: {e}")
        return None

def uniprot_protein_annotations(uniprot_id):
    """
    Return a dictionary with annotations for a UniProt ID: Dataset, Protein name, Refseq ID, Sequence, 
    Version, EC, Locus_tag and PDB/Alphafold IDs. 

    :param uniprot_id: UniProt ID.
    :return: Dictionary with annotations for the UniProt ID, or None if failed.
    """
    
    try:
        # Use retry-enabled fetch function
        uniprot_xml = fetch_uniprot_xml(uniprot_id)
        
        root = ET.fromstring(uniprot_xml)
        namespaces = {'up': 'http://uniprot.org/uniprot'}

        # Extract accession number
        accession = root.find('up:entry/up:accession', namespaces).text
        version = root.find('up:entry', namespaces).attrib['version']
        dataset = root.find('up:entry', namespaces).attrib['dataset']

        if accession != uniprot_id:
            logger.logger.warning(f"Accession mismatch: requested {uniprot_id}, got {accession}")
            return None

        result = {accession: {'Version': version, 'Dataset': dataset}}

        # Extract protein names
        recommended_name = root.findall('up:entry/up:protein/up:recommendedName/up:fullName', namespaces)
        submitted_name = root.findall('up:entry/up:protein/up:submittedName/up:fullName', namespaces)

        if recommended_name:
            result[accession]['Protein_name'] = [protein_name.text for protein_name in recommended_name][0]
        elif submitted_name:
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
            reaction = catalytic_activity_comment.find('up:reaction', namespaces)
            if reaction is not None:
                ec_reference = reaction.find('up:dbReference[@type="EC"]', namespaces)
                if ec_reference is not None:
                    result[accession]['EC_number'] = ec_reference.get('id')
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
            
        return result
        
    except Exception as e:
        logger.logger.error(f"Failed to fetch annotations for {uniprot_id} after all retries: {e}")
        return None

def uniprot_protein_annotations_batch(uniprot_ids, batch_size=500):
    """
    Return a dictionary with annotations for MULTIPLE UniProt IDs using batch queries.
    
    Fetches data from UniProt REST API in batches (default 500 IDs per request).
    Returns combined annotations for all requested IDs.
    
    :param uniprot_ids: List of UniProt accession IDs.
    :param batch_size: Number of IDs to fetch per request (default 500, max ~1000).
    :return: Dictionary with annotations {uniprot_id: {annotation_data}} or empty dict if failed.
    """
    if not uniprot_ids:
        return {}
    
    all_annotations = {}
    total_batches = (len(uniprot_ids) + batch_size - 1) // batch_size
    
    logger.logger.info(f"Fetching annotations for {len(uniprot_ids)} UniProt IDs in {total_batches} batch(es)...")
    
    for i in range(0, len(uniprot_ids), batch_size):
        batch = uniprot_ids[i:i + batch_size]
        batch_num = (i // batch_size) + 1
        
        logger.logger.info(f"Processing batch {batch_num}/{total_batches} ({len(batch)} entries)...")
        
        try:
            # Fetch batch XML
            batch_xml = fetch_uniprot_batch_xml(batch)
            
            # Parse batch XML
            root = ET.fromstring(batch_xml)
            namespaces = {'up': 'http://uniprot.org/uniprot'}
            
            # Find all entry elements in batch response
            entries = root.findall('up:entry', namespaces)
            logger.logger.info(f"  Found {len(entries)} entries in batch response")
            
            # Parse each entry
            for entry in entries:
                entry_result = parse_uniprot_entry_xml(entry)
                if entry_result:
                    all_annotations.update(entry_result)
            
            # Log batch progress
            logger.logger.info(f"  ✓ Batch {batch_num} complete: {len(all_annotations)}/{len(uniprot_ids)} total annotations retrieved")
            
            # Small delay between batches to be respectful to API
            if i + batch_size < len(uniprot_ids):
                time.sleep(2)
                
        except Exception as e:
            logger.logger.error(f"Failed to process batch {batch_num}: {e}")
            continue
    
    # Summary
    success_rate = (len(all_annotations) / len(uniprot_ids)) * 100 if uniprot_ids else 0
    logger.logger.info(f"Batch fetch complete: {len(all_annotations)}/{len(uniprot_ids)} annotations retrieved ({success_rate:.1f}%)")
    
    # Log missing IDs
    missing_ids = set(uniprot_ids) - set(all_annotations.keys())
    if missing_ids:
        logger.logger.warning(f"Missing annotations for {len(missing_ids)} IDs: {list(missing_ids)[:10]}{'...' if len(missing_ids) > 10 else ''}")
    
    return all_annotations

def download_species_uniprot_data(base_path, organism_name, specie_taxid):
    """
    Download all UniProt data for a taxonomic ID (species-level)
    Including sequences, strain information and structure availability
    
    :param tax_id: NCBI Taxonomy ID
    :param output_path: Directory to save the data
    """
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    output_path = os.path.join(structure_dir, 'uniprot_files')

    os.makedirs(output_path, exist_ok=True)
       
    # Download data from Uniprot, including strain and structure information
    uniptot_url = f"https://rest.uniprot.org/uniprotkb/stream?format=tsv&query=taxonomy_id:{specie_taxid}&fields=accession,id,reviewed,annotation_score,protein_name,organism_name,organism_id,length,xref_pdb,xref_alphafolddb,sequence"
    uniptot_file = os.path.join(output_path, f"uniprot_specie_taxid_{specie_taxid}_data.tsv")
    
    if not files.file_check(uniptot_file):
        print(f"Downloading metadata for Tax ID {specie_taxid}...")
        databases.download_with_progress(uniptot_url, uniptot_file)
    else:
        print(f"UniProt data for Tax ID {specie_taxid} already exists at {uniptot_file}. Skipping download.")

def parse_uniprot_species_data(uniptot_file, specie_taxid, strain_taxid):
    """
    Parse downloaded UniProt data into FASTA files for:
    1) Strain-specific proteins
    2) Proteins with PDB structures for the species
    3) Rest of species proteins (no PDB, not strain)

    :param uniptot_file: Path to the downloaded UniProt TSV file with species data
    :param specie_taxid: Species Taxonomy ID (e.g. 287 for P. aeruginosa)
    :param strain_taxid: Strain Taxonomy ID (e.g. 208964 for PAO1)
    """

    if files.file_check(uniptot_file) is False:
        print(f"Error: File {uniptot_file} does not exist or is empty.")
        raise FileNotFoundError(f'{uniptot_file} not found.')
    else:
        # Parse UniProt species data
        # Keep Entry and PDB columns as string to prevent scientific notation
        uniprot_df = pd.read_csv(uniptot_file, sep='\t', dtype={'Entry': str, 'PDB': str})

        # 1) Obtain strain specific data
        strain_df = uniprot_df[uniprot_df['Organism (ID)'] == strain_taxid]

        fasta_file = os.path.join(
            os.path.dirname(uniptot_file),
            f"uniprot_strain_taxid_{strain_taxid}.faa")
        with open(fasta_file, 'w') as fasta_out:
            for _, row in strain_df.iterrows():
                fasta_out.write(f">{row['Entry']}\n{row['Sequence']}\n")
        print(f"UniProt strain-specific FASTA file created: {fasta_file}")
        # 2) Obtain proteins with PDB structure data
        structure_df = uniprot_df[(uniprot_df['PDB'].notna())] 
        fasta_file_struct = os.path.join(
            os.path.dirname(uniptot_file),
            f"uniprot_PDB_structures_specie_taxid_{specie_taxid}.faa")
        with open(fasta_file_struct, 'w') as fasta_out:
            for _, row in structure_df.iterrows():
                fasta_out.write(f">{row['Entry']}\n{row['Sequence']}\n")
        print(f"UniProt PDB structure FASTA file created: {fasta_file_struct}")
        # 3) Obtain the rest of UniProt proteins for the species
        rest_mask = (~uniprot_df['Entry'].isin(strain_df['Entry']) &
                        ~uniprot_df['Entry'].isin(structure_df['Entry']))
        rest_df = uniprot_df[rest_mask]
        # Order rest_df by reviewed status
        rest_df = rest_df.sort_values(
            by=['Reviewed', 'Annotation'],
            ascending=[True, False])

        fasta_file_species_rest = os.path.join(
            os.path.dirname(uniptot_file),
            f"uniprot_species_taxid_{specie_taxid}_rest.faa")
        with open(fasta_file_species_rest, 'w') as fasta_out:
            for _, row in rest_df.iterrows():
                fasta_out.write(f">{row['Entry']}\n{row['Sequence']}\n")
        print(f"Uniprot species REST FASTA file created: {fasta_file_species_rest}")

def cluster_uniprot_specie(base_path, organism_name, specie_taxid):
    """
    Cluster uniprot proteins at 100% identity and a minimum alignment coverage of 90% using CD-HIT.
    This function uses the `run_cd_hit` function from the `programs` module.
    Clustered sequences are saved in a new .faa file.

    :param fasta_proteins: Path to the fasta file with specie-level proteins.
    """
    uniprot_dir = os.path.join(base_path, 'organism', organism_name, 'structures', 'uniprot_files')
    uniprot_species_rest_faa = os.path.join(uniprot_dir, f"uniprot_species_taxid_{specie_taxid}_rest.faa")

    try:
        programs.run_cd_hit(
            input_fasta= uniprot_species_rest_faa,
            output_fasta= uniprot_species_rest_faa.replace('.faa', '_cdhit100.faa'),
            identity= 1.0,
            aln_coverage_short= 0.9,
            aln_coverage_long= 0.9,
            use_global_seq_identity= True,
            accurate_mode= True,
            cpus= multiprocessing.cpu_count()
        )
    except Exception as e:
        print(f"Failed to run CD-HIT on file {uniprot_species_rest_faa}: {e}")

def create_uniprot_blast_db (base_path, organism_name, specie_taxid, strain_taxid):
    
    """
    Runs NCBI makeblastdb against uniprot proteome. Do this for:
    1) Strain-specific proteome
    2) PDB structure proteome (specie-level)
    3) Rest of proteome (specie-level, clustered with CD-HIT)

    This function uses the `run_makeblastdb` function from the `programs` module.
    Database is saved in the structures directory.

    :param base_path: Base path.
    :param organism_name: Name of organism.
    :param proteome_id: Uniprot proteome ID.
    """

    uniprot_dir = os.path.join(base_path, 'organism', organism_name, 'structures', 'uniprot_files')

    # For strain-specific database
    uniprot_strain_faa = os.path.join(uniprot_dir, f"uniprot_strain_taxid_{strain_taxid}.faa")
    uniprot_strain_index_path = os.path.join(uniprot_dir, f'uniprot_strain_taxid_{strain_taxid}')
    
    blast_output_strain_path= os.path.join(uniprot_dir, f'uniprot_strain_taxid_{strain_taxid}_blast.tsv')

    #Index
    if not files.file_check(blast_output_strain_path):
        print(f'Indexing uniprot strain-specific proteome for taxid {strain_taxid}')
    
        try:
            programs.run_makeblastdb(
            input= uniprot_strain_faa,
            output= uniprot_strain_index_path,
            title= f'uniprot_strain_taxid_{strain_taxid}',
            dbtype= 'prot'
            )
            print(f'Index built for strain uniprot proteome: uniprot_strain_taxid_{strain_taxid}')
        except Exception as e:
            print(f"Failed to run makeblastdb to file {uniprot_strain_faa}: {e}")
    else:
        print(f'Index already exists for strain uniprot proteome: uniprot_strain_taxid_{strain_taxid}')

    # For PDB structure database
    uniprot_pdb_faa = os.path.join(uniprot_dir, f"uniprot_PDB_structures_specie_taxid_{specie_taxid}.faa")
    uniprot_pdb_index_path = os.path.join(uniprot_dir, f'uniprot_PDB_structures_specie_taxid_{specie_taxid}')
    
    blast_output_pdb_path= os.path.join(uniprot_dir, f'uniprot_PDB_structures_specie_taxid_{specie_taxid}_blast.tsv')

    #Index
    if not files.file_check(blast_output_pdb_path):
        print(f'Indexing uniprot PDB structure proteome for taxid {specie_taxid}')
    
        try:
            programs.run_makeblastdb(
            input= uniprot_pdb_faa,
            output= uniprot_pdb_index_path,
            title= f'uniprot_PDB_structures_specie_taxid_{specie_taxid}',
            dbtype= 'prot'
            )
            print(f'Index built for PDB uniprot proteome: uniprot_PDB_structures_specie_taxid_{specie_taxid}')
        except Exception as e:
            print(f"Failed to run makeblastdb to file {uniprot_pdb_faa}: {e}")
    else:
        print(f'Index already exists for PDB uniprot proteome: uniprot_PDB_structures_specie_taxid_{specie_taxid}')
    
    # For rest of uniprot protein (specie-level)
    uniprot_species_rest_faa = os.path.join(uniprot_dir, f"uniprot_species_taxid_{specie_taxid}_rest_cdhit100.faa")
    uniprot_species_rest_index_path = os.path.join(uniprot_dir, f'uniprot_species_taxid_{specie_taxid}_rest_cdhit100')

    blast_output_species_rest_path= os.path.join(uniprot_dir, f'uniprot_species_taxid_{specie_taxid}_rest_cdhit100_blast.tsv')

    #Index
    if not files.file_check(blast_output_species_rest_path):
        print(f'Indexing uniprot species-level rest proteome for taxid {specie_taxid}')

        try:
            programs.run_makeblastdb(
            input= uniprot_species_rest_faa,
            output= uniprot_species_rest_index_path,
            title= f'uniprot_species_taxid_{specie_taxid}_rest',
            dbtype= 'prot'
            )
            print(f'Index finished for species-level rest uniprot proteome: uniprot_species_taxid_{specie_taxid}_rest')
        except Exception as e:
            print(f"Failed to run makeblastdb to file {uniprot_species_rest_faa}: {e}")
    else:
        print(f'Index already exists for species-level rest uniprot proteome: uniprot_species_taxid_{specie_taxid}_rest')   

def uniprot_proteome_blast (base_path, organism_name, specie_taxid, strain_taxid, cpus=multiprocessing.cpu_count()):
    """
    Runs NCBI Blastp against uniprot proteome databases. Do this for:
    1) Strain-specific proteome
    2) PDB structure proteome (specie-level)
    3) Rest of proteome (specie-level, clustered with CD-HIT)
    This function uses the `run_blastp` function from the `programs` module. 
    Results are saved in a .tsv file in the structures/uniprot_files directory, with separate files for each database.

    :param base_path: Base path.
    :param organism_name: Name of organism.
    :param specie_taxid: Species Taxonomy ID.
    :param strain_taxid: Strain Taxonomy ID.
    :param cpus: Number of threads (CPUs) to use in blast search.
    """

    uniprot_dir = os.path.join(base_path, 'organism', organism_name, 'structures', 'uniprot_files')

    # Query: proteins from the organism genome
    organism_prot_seq_path = os.path.join(base_path, 'organism', f'{organism_name}','genome', f'{organism_name}.faa')

    # For strain-specific database
    uniprot_strain_index_path = os.path.join(uniprot_dir, f'uniprot_strain_taxid_{strain_taxid}')
    blast_output_strain_path= os.path.join(uniprot_dir, f'uniprot_strain_taxid_{strain_taxid}_blast.tsv')
    
    if not files.file_check(blast_output_strain_path):
        print(f'Runing blastp for {organism_name} and uniprot strain-specific proteome taxid {strain_taxid}')

        try:
            programs.run_blastp(
                blastdb= uniprot_strain_index_path,
                query= organism_prot_seq_path,
                output=blast_output_strain_path,
                evalue= '1e-5',
                outfmt= '6 std qcovhsp qcovs',
                max_target_seqs = '5',
                cpus=cpus
            )
            print(f'Blastp finished for uniprot strain-specific proteome taxid {strain_taxid}')
            print(f'Blastp results saved in {blast_output_strain_path}.')

        except Exception as e:
            print(f"Failed to run blastp to file {uniprot_strain_index_path}: {e}")
    else:
        print(f'Blastp results in {blast_output_strain_path}.')
    
    # For PDB structure database
    uniprot_pdb_index_path = os.path.join(uniprot_dir, f'uniprot_PDB_structures_specie_taxid_{specie_taxid}')
    blast_output_pdb_path= os.path.join(uniprot_dir, f'uniprot_PDB_structures_specie_taxid_{specie_taxid}_blast.tsv')

    if not files.file_check(blast_output_pdb_path):
        print(f'Runing blastp for {organism_name} and uniprot PDB structure proteome for taxid {specie_taxid}')

        try:
            programs.run_blastp(
                blastdb= uniprot_pdb_index_path,
                query= organism_prot_seq_path,
                output=blast_output_pdb_path,
                evalue= '1e-5',
                outfmt= '6 std qcovhsp qcovs',
                max_target_seqs = '5',
                cpus=cpus
            )
            print(f'Blastp finished for uniprot PDB structure proteome for taxid {specie_taxid}')
            print(f'Blastp results saved in {blast_output_pdb_path}.')

        except Exception as e:
            print(f"Failed to run blastp to file {uniprot_pdb_index_path}: {e}")
    else:
        print(f'Blastp results in {blast_output_pdb_path}.')

    # For rest of uniprot protein (specie-level)
    uniprot_species_rest_index_path = os.path.join(uniprot_dir, f'uniprot_species_taxid_{specie_taxid}_rest_cdhit100')
    blast_output_species_rest_path= os.path.join(uniprot_dir, f'uniprot_species_taxid_{specie_taxid}_rest_cdhit100_blast.tsv')


    if not files.file_check(blast_output_species_rest_path):
        print(f'Runing blastp for {organism_name} and uniprot Specie structure proteome for taxid {specie_taxid}')

        try:
            programs.run_blastp(
                blastdb= uniprot_species_rest_index_path,
                query= organism_prot_seq_path,
                output=blast_output_species_rest_path,
                evalue= '1e-5',
                outfmt= '6 std qcovhsp qcovs',
                max_target_seqs = '5',
                cpus=cpus
            )
            print(f'Blastp finished for uniprot Specie structure proteome for taxid {specie_taxid}')
            print(f'Blastp results saved in {blast_output_species_rest_path}.')

        except Exception as e:
            print(f"Failed to run blastp to file {uniprot_species_rest_index_path}: {e}")
    else:
        print(f'Blastp results in {blast_output_species_rest_path}.')

def parse_best_result_blast (file, identity_cutoff=95, coverage_cutoff=90, all_hits=False):
    """
    Parse blast output and return best hit per query based on identity and coverage.

    :param file: Path to blast output file.
    :param identity_cutoff: Minimum percentage identity to consider a hit. Default is 95.
    :param coverage_cutoff: Minimum query coverage to consider a hit. Default is 90.
    :return: Dictionary with query id as key and best subject id as value.
    """

    print(f'Reading blastp results for {file}.')

    blast_output_df = files.read_blast_output(file)

    # Apply filters
    blast_output_filtered_df = blast_output_df[(blast_output_df['pident'] >= identity_cutoff) & 
                                                (blast_output_df['qcovs'] >= coverage_cutoff)]

    best_hits = {}

    for qseqid, group in blast_output_filtered_df.groupby('qseqid'):

        # Find the maximum identity value in the group
        max_identity = group['pident'].max()
        # Filter rows that have the maximum identity
        best_identity = group[group['pident'] == max_identity]

        # Find the maximum coverage value in the group
        max_coverage = best_identity['qcovs'].max()
        # Filter rows that have the maximum coverage
        best_results = best_identity[best_identity['qcovs'] == max_coverage]

        if all_hits:
            best_hits[qseqid] = best_results['sseqid'].tolist()
        else:
            best_hits[qseqid] = best_results['sseqid'].iloc[0]

    return best_hits

def uniprot_proteome_mapping (base_path, organism_name, specie_taxid, strain_taxid):

    """
    Map locus_tags to uniprot ids using blastp results for specie specific proteome.
    Follows this priorrity to asign uniprot ids:
    1) Look for proteins that have an structures in PDB (specie-level)
    2) Map proteins to strain-specific uniprot proteome
    3) Map proteins to rest of specie-level uniprot proteome (cross-strain search)

    :param base_path: Base path.
    :param organism_name: Name of organism.
    :param strain_taxid: Strain Taxonomy ID.
    :return: Dictionary with locus_tag as key and list of uniprot ids as values.

    """
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    uniprot_dir = os.path.join(structure_dir, 'uniprot_files')
    
    # Blastp search result files
    # PDB
    blast_output_pdb_path= os.path.join(uniprot_dir, f'uniprot_PDB_structures_specie_taxid_{specie_taxid}_blast.tsv')
    # Strain specific 
    blast_output_strain_path= os.path.join(uniprot_dir, f'uniprot_strain_taxid_{strain_taxid}_blast.tsv')
    # Rest of Uniprot 
    blast_output_species_rest_path= os.path.join(uniprot_dir, f'uniprot_species_taxid_{specie_taxid}_rest_cdhit100_blast.tsv')

    # File with mapping IDs 
    proteome_ids_file = os.path.join(uniprot_dir, f'uniprot_{organism_name}_id_mapping.json')

    if not files.file_check(proteome_ids_file):
           
        # 1) Map with IDs with PDB structures
        if files.file_check(blast_output_pdb_path):
            parse_pdb = parse_best_result_blast (blast_output_pdb_path, 
                                    identity_cutoff=95, 
                                    coverage_cutoff=90, 
                                    all_hits=True)
        else:
            raise FileNotFoundError(f'{blast_output_pdb_path} not found.')
        
        # 2) Map with IDs within the strain
        if files.file_check(blast_output_strain_path):
            parse_strain = parse_best_result_blast (blast_output_strain_path, 
                                    identity_cutoff=95, 
                                    coverage_cutoff=90, 
                                    all_hits=False)
        else:
            raise FileNotFoundError(f'{blast_output_strain_path} not found.')
        
        # 3) Map with IDs within the specie
        if files.file_check(blast_output_species_rest_path):
            parse_rest = parse_best_result_blast (blast_output_species_rest_path, 
                                    identity_cutoff=95, 
                                    coverage_cutoff=90, 
                                    all_hits=False)
        else:
            raise FileNotFoundError(f'{blast_output_species_rest_path} not found.')

        # Merge mapping results

        map_results = {}

        all_locus_tags = set(metadata.ref_gbk_locus(base_path, organism_name))
        for locus_tag in all_locus_tags:
            if locus_tag in parse_pdb:
                map_results[locus_tag] = parse_pdb[locus_tag]
            elif locus_tag in parse_strain:
                map_results[locus_tag] = parse_strain[locus_tag]
            elif locus_tag in parse_rest:
                map_results[locus_tag] = parse_rest[locus_tag]
            else:
                map_results[locus_tag] = None
                print(f'No UniProt ID found for {locus_tag}')
        
        files.dict_to_json(uniprot_dir, f'uniprot_{organism_name}_id_mapping.json', map_results)
    else:
        map_results = files.json_to_dict(proteome_ids_file)
        print(f'IDs file in {proteome_ids_file}.')
    
    return map_results

def parse_resolution(res_str):
    """
    Converts resolution string to float (in Angstroms).
    """
    if not res_str:
        return None
    try:
        return float(str(res_str).split()[0])
    except Exception:
        return None

def parse_chain_string(chains_str):
    """
    Converts chain string from UniProt into list of (chain_id, start, end) tuples.
    
    Handles multiple formats from UniProt PDB chain annotations:
    - Single chain: "A=1-139" → [('A', 1, 139)]
    - Homooligomer (identical chains): "A/B/C=1-139" → [('A', 1, 139)]  # Only first chain
    - Heterooligomer (different subunits): "A=24-193, B=217-762" → [('A', 24, 193), ('B', 217, 762)]  # ALL chains
    - Mixed format: "A/B=27-512, C/P=528-536" → [('A', 27, 512), ('C', 528, 536)]  # First from each group
    - Chain without range: "A" → [('A', None, None)]
    
    Strategy:
    - For chains separated by "/" (e.g., "A/B/C"): These are identical copies, select FIRST only
    - For chains separated by "," (e.g., "A=..., B=..."): These are different subunits, keep ALL
    
    :param chains_str: Chain string from UniProt (e.g., "A/B=1-100, C=101-200")
    :return: List of (chain_id, start, end) tuples representing unique functional chains.
    """
    if not chains_str:
        return []

    result = []
    
    # Split by comma first to handle different subunits (HETEROOLIGOMER)
    subunits = [s.strip() for s in str(chains_str).split(',')]
    
    for subunit in subunits:
        if not subunit:
            continue
            
        # Check if there's an equals sign (chain=range format)
        if '=' in subunit:
            chain_part, range_part = subunit.split('=', 1)
            chain_part = chain_part.strip()
            range_part = range_part.strip()
            
            # Handle multiple chains with same range (A/B/C=1-139) - HOMOOLIGOMER
            # These are identical - select FIRST chain only
            if '/' in chain_part:
                chains = [c.strip() for c in chain_part.split('/')]
                chain_id = chains[0]  # Select first chain from homooligomer
                print(f"    Homooligomer detected: '{chain_part}' → selecting chain {chain_id}")
            else:
                chain_id = chain_part
            
            # Parse range
            if '-' in range_part:
                try:
                    start_str, end_str = range_part.split('-', 1)
                    start = int(start_str.strip())
                    end = int(end_str.strip())
                except ValueError:
                    print(f"    Warning: Could not parse range '{range_part}', using None")
                    start = end = None
            else:
                # Single position
                try:
                    start = int(range_part.strip())
                    end = start  # Single residue
                except ValueError:
                    print(f"    Warning: Could not parse position '{range_part}', using None")
                    start = end = None
        else:
            # No range specified, just chain ID (may have slashes)
            if '/' in subunit:
                chains = [c.strip() for c in subunit.split('/')]
                chain_id = chains[0]
                print(f"    Multiple chains without range: '{subunit}' → selecting chain {chain_id}")
            else:
                chain_id = subunit.strip()
            start = end = None
        
        result.append((chain_id, start, end))
    
    return result

def compute_coverage(start, end, seq_len):
    """
    Calculates coverage percentage given start, end, and sequence length.
    coverage = (end - start + 1) / seq_len * 100
    """
    if start is None or end is None or seq_len <= 0:
        return None
    covered = max(0, end - start + 1)
    return (covered / seq_len) * 100.0

def collect_structures_for_uniprot(uniprot_id, uni_info, resolution_cutoff = 3.5):
    """
    Given a Uniprot ID and its info dict from Uniprot,
    collects all structure entries (PDB and AlphaFold) and returns
    a list of dicts with structure information.
    
    For heterooligomers (multiple chains separated by commas), stores ALL chain IDs
    in the 'chain' field as a semicolon-separated string.
    
    :param uniprot_id: Uniprot accession ID.
    :param uni_info: Dict with Uniprot info.
    :param resolution_cutoff: Resolution (Å) above which AlphaFold is also added alongside PDB.
    :return: List of dicts with structure information.
    """
    seq = uni_info.get("Sequence", "")
    seq_len = len(seq)

    rows = []

    # --- PDBs ---
    pdb_entries = uni_info.get("PDB_id")
    has_good_resolution_pdb = False  # Track if any PDB passes resolution cutoff
    
    if pdb_entries:

        best_by_id = {}

        for entry in pdb_entries:
            if not isinstance(entry, dict):
                continue
            pdb_id = entry.get("ID")
            if not pdb_id:
                continue

            existing = best_by_id.get(pdb_id, {})

            merged = dict(existing)
            merged.update({k: v for k, v in entry.items() if v is not None})
            best_by_id[pdb_id] = merged

        for pdb_id, entry in best_by_id.items():
            method = entry.get("method")
            resolution_str = entry.get("resolution")
            chains_str = entry.get("chains")

            resolution = parse_resolution(resolution_str)
            chain_info_list = parse_chain_string(chains_str)

            # Check if this PDB has good resolution
            if resolution is not None and resolution <= resolution_cutoff:
                has_good_resolution_pdb = True

            # For heterooligomers: store ALL chains as semicolon-separated string
            if chain_info_list:
                if len(chain_info_list) > 1:
                    # Heterooligomer - multiple different subunits
                    chain_ids = [c[0] for c in chain_info_list]
                    chain_id_str = ';'.join(chain_ids)
                    print(f"    Heterooligomer detected for {pdb_id}: chains {chain_id_str}")
                    # Use range from all chain for coverage calculation
                    start, end = chain_info_list[0][1], chain_info_list[len(chain_info_list)-1][2]
                else:
                    # Single chain or homooligomer
                    chain_id_str = chain_info_list[0][0]
                    start, end = chain_info_list[0][1], chain_info_list[0][2]
            else:
                chain_id_str = None
                start = end = None

            coverage = compute_coverage(start, end, seq_len)

            rows.append({
                "uniprot_id": uniprot_id,
                "structure_type": "PDB",
                "structure_id": pdb_id,
                "method": method,
                "resolution": resolution,
                "chain": chain_id_str,  # Now can be "A;B;C" for heterooligomers
                "residue_range": f"{start}-{end}" if start is not None and end is not None else None,
                "coverage": coverage,
                "sequence_length": seq_len,
                "is_reference": False,
            })

        # Add AlphaFold only if NO PDB passed the resolution cutoff
        if not has_good_resolution_pdb:
            af_id = uni_info.get("AlphaFoldDB")
            if af_id:
                rows.append({
                    "uniprot_id": uniprot_id,
                    "structure_type": "AlphaFold",
                    "structure_id": af_id,
                    "method": "AlphaFold",
                    "resolution": None,
                    "chain": "A",
                    "residue_range": f"1-{seq_len}" if seq_len > 0 else None,
                    "coverage": 100.0 if seq_len > 0 else None,
                    "sequence_length": seq_len,
                    "is_reference": False,
                })
    else:
        # --- AlphaFold ---
        af_id = uni_info.get("AlphaFoldDB")
        if af_id:
            rows.append({
                "uniprot_id": uniprot_id,
                "structure_type": "AlphaFold",
                "structure_id": af_id,
                "method": "AlphaFold",
                "resolution": None,
                "chain": "A",
                "residue_range": f"1-{seq_len}" if seq_len > 0 else None,
                "coverage": 100.0 if seq_len > 0 else None,
                "sequence_length": seq_len,
                "is_reference": False,
            })
        else:
            # No AlphaFoldDB ID in UniProt, but we can still try to download from AlphaFold
            # using the UniProt ID directly (fallback mechanism)
            print(f"  [INFO] No AlphaFoldDB ID in UniProt for {uniprot_id}, will attempt direct AlphaFold download")
            rows.append({
                "uniprot_id": uniprot_id,
                "structure_type": "AlphaFold",
                "structure_id": uniprot_id,  # Use UniProt ID as structure ID for fallback
                "method": "AlphaFold",
                "resolution": None,
                "chain": "A",
                "residue_range": f"1-{seq_len}" if seq_len > 0 else None,
                "coverage": 100.0 if seq_len > 0 else None,
                "sequence_length": seq_len,
                "is_reference": False,
            })

    return rows

def select_reference_structure(structs, resolution_cutoff= 3.5):
    """
    Given a list of structure dicts (from collect_structures_for_uniprot),
    selects the best structure to use as reference according to criteria:
      1) Prefer PDB X-ray structures, best resolution, best coverage.
      2) If no X-ray, prefer PDB EM structures, best resolution, best coverage.
         If best EM res > resolution_cutoff and AlphaFold exists, prefer AlphaFold.
      3) If no PDB, use AlphaFold if exists.
      4) Otherwise, None.
    Returns the index of the selected structure in the list, or None.
    :param structs: List of structure dicts.
    :return: Index of selected structure, or None.
    """
    if not structs:
        return None

    def method_class(s):
        m = (s.get("method") or "").lower()
        if "x-ray" in m or "xray" in m or "x-ray diffraction" in m:
            return "xray"
        if "em" in m or "electron microscopy" in m:
            return "em"
        if "alphafold" in m:
            return "alphafold"
        return "other"

    pdb_xray_idx = []
    pdb_em_idx = []
    alphafold_idx = []

    for i, s in enumerate(structs):
        mclass = method_class(s)
        stype = s.get("structure_type")
        if stype == "PDB":
            if mclass == "xray":
                pdb_xray_idx.append(i)
            elif mclass == "em":
                pdb_em_idx.append(i)
        elif stype == "AlphaFold":
            alphafold_idx.append(i)

    def sort_key(idx):
        s = structs[idx]
        res = s.get("resolution")
        cov = s.get("coverage")
        if res is None:
            res = float("inf")
        if cov is None:
            cov = -1.0
        return (res, -cov)

    # 1) If there are X-ray structures, pick the best
    if pdb_xray_idx:
        return sorted(pdb_xray_idx, key=sort_key)[0]

    # 2) If no X-ray, look at EM
    if pdb_em_idx:
        # best EM
        best_em_idx = sorted(pdb_em_idx, key=sort_key)[0]
        best_em = structs[best_em_idx]
        best_em_res = best_em.get("resolution")

        # if best EM res > resolution_cutoff and AlphaFold exists, prefer AlphaFold
        if best_em_res is not None and best_em_res > resolution_cutoff and alphafold_idx:
            return sorted(alphafold_idx, key=sort_key)[0]
        else:
            return best_em_idx

    # 3) If no PDB, use AlphaFold if exists
    if alphafold_idx:
        return sorted(alphafold_idx, key=sort_key)[0]

    # 4) Otherwise, None
    return None

def create_subfolder_structures(base_path, organism_name):
    """
    Create structures subfolder within structure folder.
    Create one folder per locus_tag and one per UniProt ID.

    :param base_path: Base path.
    :param organism_name: Name of organism.
    """
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    proteome_ids_file = os.path.join(structure_dir, 'uniprot_files', f'uniprot_{organism_name}_id_mapping.json')

    if files.file_check(proteome_ids_file):
        print(f'Creating subfolders for structures in {structure_dir}...')
        map_results = files.json_to_dict(proteome_ids_file)
        
        for locus_tag in map_results.keys():
            locus_tag_dir = os.path.join(structure_dir, locus_tag)
            os.makedirs(locus_tag_dir, exist_ok=True)

            uniprot_ids = map_results[locus_tag]
            if uniprot_ids:
                if isinstance(uniprot_ids, list):
                    for uniprot_id in uniprot_ids:
                        uniprot_dir = os.path.join(locus_tag_dir, uniprot_id)
                        os.makedirs(uniprot_dir, exist_ok=True)
                else:
                    uniprot_dir = os.path.join(locus_tag_dir, uniprot_ids)
                    os.makedirs(uniprot_dir, exist_ok=True)
    else:
        raise FileNotFoundError(f'{proteome_ids_file} not found.')

def create_summary_structure_table(batch_annotations, mapping_dict, locus_tag, resolution_cutoff = 3.5):
    """
    Create a summary table with structure information for a locus_tag.
    Only ONE structure will be marked as reference across all UniProt IDs.
    
    Priority for reference selection:
    1. Best structure from first UniProt ID with PDB structures
    2. If no PDB in first UniProt, check AlphaFold
    3. If first UniProt has no structures, check next UniProt IDs
    
    :param mapping_dict: Dictionary mapping locus_tags to UniProt IDs.
    :param batch_annotations: Dictionary with UniProt annotations for all UniProt IDs.
    :param locus_tag: Locus tag to create summary for.
    :return: DataFrame with structure summary.
    """
    
    summary_data = []
    uniprot_ids = mapping_dict.get(locus_tag)

    if uniprot_ids:
        # Normalize to list
        if not isinstance(uniprot_ids, list):
            uniprot_ids = [uniprot_ids]
        
        # Collect all structures from all UniProt IDs
        for uniprot_id in uniprot_ids:    
            if uniprot_id in batch_annotations:
                uni_info = batch_annotations[uniprot_id]
                structures_info = collect_structures_for_uniprot(uniprot_id, uni_info, resolution_cutoff)
                for struct in structures_info:
                    summary_data.append(struct)
            else:
                print(f'  Warning: No annotations found for UniProt ID {uniprot_id}')
        

    if summary_data:
        # Select ONE reference structure across ALL UniProt IDs for this locus_tag
        ref_idx = select_reference_structure(summary_data, resolution_cutoff)
        if ref_idx is not None:
            summary_data[ref_idx]["is_reference"] = True
            print(f'  Selected reference: {summary_data[ref_idx]["structure_id"]} (UniProt: {summary_data[ref_idx]["uniprot_id"]})')
        else:
            print(f'  Warning: No reference structure selected for {locus_tag}.')

        summary_df = pd.DataFrame(summary_data)
    else:
        print(f'  No structure data found for {locus_tag}, skipping summary table creation.')
        summary_df = pd.DataFrame()

    return summary_df

def create_summary_structure_file(base_path, organism_name, resolution_cutoff= 3.5):

    """
    Create a TSV summary file with structure information for each locus_tag.
    The table includes:
    - uniprot_id
    - AlphaFold ID
    - PDB IDs

    :param base_path: Base path.
    :param organism_name: Name of organism.

    """
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    uniprot_files_dir = os.path.join(structure_dir, 'uniprot_files')
    proteome_ids_file = os.path.join(uniprot_files_dir, f'uniprot_{organism_name}_id_mapping.json')

    if not files.file_check(proteome_ids_file):
        raise FileNotFoundError(f'{proteome_ids_file} not found.')

    map_results = files.json_to_dict(proteome_ids_file)

    print(f'Creating structure summary table for each locus_tag in {structure_dir}...')

    all_uniprot_ids = []
    for uniprot_ids in map_results.values():
        if isinstance(uniprot_ids, list):
            all_uniprot_ids.extend(uniprot_ids)
        elif uniprot_ids:
            all_uniprot_ids.append(uniprot_ids)

    # Fetch all annotations at once
    if not files.file_check(os.path.join(structure_dir, 'uniprot_files', f'uniprot_{organism_name}_annotations.json')):
        print(f'Fetching UniProt annotations for {len(set(all_uniprot_ids))} unique UniProt IDs...')
        batch_annotations = uniprot_protein_annotations_batch(list(set(all_uniprot_ids)))
        files.dict_to_json(uniprot_files_dir, f'uniprot_{organism_name}_annotations.json', batch_annotations)
    else:
        batch_annotations = files.json_to_dict(os.path.join(uniprot_files_dir, f'uniprot_{organism_name}_annotations.json'))
        print(f'UniProt annotations file found for {len(batch_annotations)} UniProt IDs.')

    for locus_tag in tqdm(map_results.keys(), desc='Locus tags'):
        
        summary_table_path = os.path.join(structure_dir, locus_tag, f"{locus_tag}_structure_summary.tsv")

        if not files.file_check(summary_table_path):
            
            summary_df = create_summary_structure_table(batch_annotations, map_results, locus_tag, resolution_cutoff)
            if not summary_df.empty:
                # Ensure structure_id is stored as string to prevent scientific notation issues
                summary_df['structure_id'] = summary_df['structure_id'].astype(str)
                summary_df.to_csv(summary_table_path, sep='\t', index=False)
            else:
                print(f'No structure data found for {locus_tag}, summary table not created.')

        else:
            print(f'Structure summary table already exists for {locus_tag}.')

##  ------------------- Download structures functions ------------------- ## 
def get_structure_PDB (output_path, PDB_id):
    """
    Download structure from PDB. Returns True if successful.
    
    :param output_path: Output directory path.
    :param PDB_id: ID of PDB structure.

    :return: True if successful, False otherwise.
    """
    res = False
    file_path = os.path.join(output_path, f"PDB_{PDB_id}.pdb")

    pdb_url = f"https://files.rcsb.org/download/{PDB_id}.pdb"

    for attempt in range(3):
        try:
            response = requests.get(pdb_url, timeout=30)
            if response.status_code == 200:
                # Verify we got content
                if len(response.content) > 0:
                    with open(file_path, 'wb') as file:
                        file.write(response.content)
                    # Verify file was written successfully
                    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
                        print(f"Downloaded {PDB_id}.pdb.")
                        res = True
                        break
                    else:
                        print(f"Warning: File {file_path} not written correctly, retrying...")
                else:
                    print(f"Warning: Empty content received for {PDB_id}.pdb")
            elif response.status_code == 404:
                # 404 means PDB format doesn't exist, no point retrying
                print(f"PDB format not available for {PDB_id} (404).")
                break
            else:
                print(f"Failed to download .pdb for {PDB_id}, status code: {response.status_code}")
        except requests.exceptions.RequestException as e:
            print(f'Failed to download .pdb for {PDB_id}, attempt {attempt + 1}/3: {e}')
        
        if attempt == 2 and not res:
            print(f"Failed to download .pdb for {PDB_id} after 3 attempts.")
    
    return res

def get_structure_CIF (output_path, PDB_id):
    """
    Download structure from PDB in CIF format. Returns True if successful.
    
    :param output_path: Output directory path.
    :param PDB_id: ID of PDB structure.

    :return: True if successful, False otherwise.
    """
    res = False
    file_path = os.path.join(output_path, f"PDB_{PDB_id}.cif")

    pdb_url = f"https://files.rcsb.org/download/{PDB_id}.cif"

    for attempt in range(3):
        try:
            response = requests.get(pdb_url, timeout=30)
            if response.status_code == 200:
                # Verify we got content
                if len(response.content) > 0:
                    with open(file_path, 'wb') as file:
                        file.write(response.content)
                    # Verify file was written successfully
                    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
                        print(f"Downloaded {PDB_id}.cif.")
                        res = True
                        break
                    else:
                        print(f"Warning: File {file_path} not written correctly, retrying...")
                else:
                    print(f"Warning: Empty content received for {PDB_id}.cif")
            elif response.status_code == 404:
                # 404 means structure doesn't exist at all
                print(f"Structure {PDB_id} not found (404).")
                break
            else:
                print(f"Failed to download .cif for {PDB_id}, status code: {response.status_code}")
        except requests.exceptions.RequestException as e:
            print(f'Failed to download .cif for {PDB_id}, attempt {attempt + 1}/3: {e}')
        
        if attempt == 2 and not res:
            print(f"Failed to download .cif for {PDB_id} after 3 attempts.")
    
    return res

def get_structure_alphafold(output_path, uniprot_id):
    """
    Download structure from AlphaFold. Returns True if successful.
    
    :param output_path: Output directory path.
    :param uniprot_id: ID of UniProt with AlphaFold structure.

    :return: True if successful, False otherwise.
    """

    res = False
    file_path = os.path.join(output_path, f"AF_{uniprot_id}.pdb")

    alphafold_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v6.pdb"

    for attempt in range(3):
        try: 
            response = requests.get(alphafold_url, timeout=30)
            if response.status_code == 200:
                # Verify we got content
                if len(response.content) > 0:
                    with open(file_path, 'wb') as file:
                        file.write(response.content)
                    # Verify file was written successfully
                    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
                        print(f"Downloaded AlphaFold model {uniprot_id}.")
                        res = True
                        break
                    else:
                        print(f"Warning: File {file_path} not written correctly, retrying...")
                else:
                    print(f"Warning: Empty content received for AlphaFold {uniprot_id}")
            elif response.status_code == 404:
                # 404 means AlphaFold prediction doesn't exist for this UniProt ID
                print(f"AlphaFold prediction not available for {uniprot_id} (404).")
                break
            else:
                print(f"Failed to download AlphaFold for {uniprot_id}, status code: {response.status_code}")
        except requests.exceptions.RequestException as e:
            print(f"Failed to download AlphaFold prediction for {uniprot_id}, attempt {attempt + 1}/3: {e}")
        
        if attempt == 2 and not res:
            print(f"Failed to download AlphaFold prediction for {uniprot_id} after 3 attempts.")
   
    return res

def download_structures(base_path, organism_name):

    """
    Download structures for each locus_tag based on the summary table.
    
    Priority: Download .pdb if available, otherwise download .cif (for cryo-EM structures).
    
    :param base_path: Base path.
    :param organism_name: Name of organism.

    """
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    
    all_locus_tags = metadata.ref_gbk_locus(base_path, organism_name)

    for locus_tag in tqdm(all_locus_tags, desc='Locus tags'):

        summary_table_path = os.path.join(structure_dir, locus_tag, f"{locus_tag}_structure_summary.tsv")

        if files.file_check(summary_table_path):
            # Read structure_id as string to prevent scientific notation (e.g., 3E59 -> 3e+59)
            summary_df = pd.read_csv(summary_table_path, sep='\t', dtype={'structure_id': str})

            uniprot_ids = summary_df['uniprot_id'].unique()

            for uniprot_id in uniprot_ids:
                uniprot_dir = os.path.join(structure_dir, locus_tag, uniprot_id)

                struct_rows = summary_df[summary_df['uniprot_id'] == uniprot_id]

                for _, row in struct_rows.iterrows():
                    struct_type = row['structure_type']
                    struct_id = row['structure_id']

                    if struct_type == 'PDB':
                        # Define file paths
                        pdb_file_path = os.path.join(uniprot_dir, f"PDB_{struct_id}.pdb")
                        cif_file_path = os.path.join(uniprot_dir, f"PDB_{struct_id}.cif")
                        
                        # Check if we already have a valid structure file
                        # Use try-except to handle potential race conditions or filesystem issues
                        try:
                            pdb_exists = files.file_check(pdb_file_path) and os.path.getsize(pdb_file_path) > 0
                        except OSError:
                            pdb_exists = False
                        
                        try:
                            cif_exists = files.file_check(cif_file_path) and os.path.getsize(cif_file_path) > 0
                        except OSError:
                            cif_exists = False
                        
                        if pdb_exists or cif_exists:
                            continue
                        
                        # Try PDB format first
                        success_pdb = get_structure_PDB(uniprot_dir, struct_id)
                        
                        # If PDB fails, try CIF (common for cryo-EM structures)
                        if not success_pdb:
                            print(f"  PDB format not available for {struct_id}, trying CIF format...")
                            success_cif = get_structure_CIF(uniprot_dir, struct_id)
                            
                            if not success_cif:
                                print(f"  Warning: Could not download {struct_id} in either PDB or CIF format.")

                    elif struct_type == 'AlphaFold':
                        # Check and download AlphaFold structure
                        # Note: This handles both cases:
                        # 1. AlphaFoldDB ID from UniProt API
                        # 2. Fallback attempt using UniProt ID directly when no AlphaFoldDB ID was found
                        af_file_path = os.path.join(uniprot_dir, f"AF_{uniprot_id}.pdb")
                        
                        # Check if we already have a valid AlphaFold file
                        try:
                            af_exists = files.file_check(af_file_path) and os.path.getsize(af_file_path) > 0
                        except OSError:
                            af_exists = False
                        
                        if not af_exists:
                            # Download using UniProt ID (works for both explicit and fallback cases)
                            success_af = get_structure_alphafold(uniprot_dir, uniprot_id)
                            if not success_af:
                                print(f"  Warning: Could not download AlphaFold model for {uniprot_id}.")
        else:
            print(f'Structure summary table not found for {locus_tag}, skipping download check.')

def extract_chain_from_pdb(pdb_file, chain_ids, output_file):
    """
    Extracts one or more specific chains from a PDB file using Biopython.
    
    For heterooligomers, pass multiple chain IDs to extract all functional subunits
    into a single output file.

    :param pdb_file: Path to the input PDB file.
    :param chain_ids: Single chain ID (str) or list of chain IDs to extract.
    :param output_file: Path to the output PDB file where the extracted chain(s) will be saved.
    """

    # Normalize to list
    if isinstance(chain_ids, str):
        if ';' in chain_ids:
            chain_ids = [c.strip() for c in chain_ids.split(';')]
        else:
            chain_ids = [chain_ids]
    
    class MultiChainSelect(Select):
        """
        Biopython Select subclass to filter multiple chains.
        
        Accepts chains whose IDs match any in the provided list.
        Used for extracting heterooligomers (multiple different subunits).
        """
        def __init__(self, chain_ids: list):
            super().__init__()
            self.chain_ids = [str(c).strip() for c in chain_ids]

        def accept_chain(self, chain) -> bool:
            """Return True if chain ID is in the list of requested chains."""
            return chain.id in self.chain_ids
    
    parser = PDBParser(QUIET=True)

    structure_id = Path(pdb_file).stem or "structure"
    structure = parser.get_structure(structure_id, pdb_file)
    
    # Verify that requested chains exist in the structure
    available_chains = [chain.id for model in structure for chain in model]
    missing_chains = [c for c in chain_ids if c not in available_chains]
    
    if missing_chains:
        print(f"    Warning: Chain(s) {missing_chains} not found in {Path(pdb_file).name}")
        print(f"    Available chains: {available_chains}")
        # Continue with available chains only
        chain_ids = [c for c in chain_ids if c in available_chains]
        if not chain_ids:
            raise ValueError(f"None of the requested chains found in {pdb_file}")

    io = PDBIO()
    io.set_structure(structure)

    selector = MultiChainSelect(chain_ids)
    io.save(output_file, select=selector)
    
    # Verify output file was created and has content
    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
        raise IOError(f"Failed to create output file {output_file}")
    
    chains_str = ','.join(chain_ids)
    print(f"    Extracted chain(s) {chains_str} from {Path(pdb_file).name}")

def extract_chain_from_cif(cif_file, chain_ids, output_file):
    """
    Extracts one or more specific chains from a CIF file using Biopython.
    
    This function is similar to extract_chain_from_pdb but works with mmCIF format files,
    which are commonly used for cryo-EM structures.
    
    For heterooligomers, pass multiple chain IDs to extract all functional subunits
    into a single output file.

    :param cif_file: Path to the input CIF file.
    :param chain_ids: Single chain ID (str) or list of chain IDs to extract.
    :param output_file: Path to the output CIF file where the extracted chain(s) will be saved.
    """
    from Bio.PDB import MMCIFParser
    
    # Normalize to list
    if isinstance(chain_ids, str):
        if ';' in chain_ids:
            chain_ids = [c.strip() for c in chain_ids.split(';')]
        else:
            chain_ids = [chain_ids]
    
    class MultiChainSelect(Select):
        """
        Biopython Select subclass to filter multiple chains.
        
        Accepts chains whose IDs match any in the provided list.
        Used for extracting heterooligomers (multiple different subunits).
        """
        def __init__(self, chain_ids: list):
            super().__init__()
            self.chain_ids = [str(c).strip() for c in chain_ids]

        def accept_chain(self, chain) -> bool:
            """Return True if chain ID is in the list of requested chains."""
            return chain.id in self.chain_ids
    
    parser = MMCIFParser(QUIET=True)

    structure_id = Path(cif_file).stem or "structure"
    structure = parser.get_structure(structure_id, cif_file)
    
    # Verify that requested chains exist in the structure
    available_chains = [chain.id for model in structure for chain in model]
    missing_chains = [c for c in chain_ids if c not in available_chains]
    
    if missing_chains:
        print(f"    Warning: Chain(s) {missing_chains} not found in {Path(cif_file).name}")
        print(f"    Available chains: {available_chains}")
        # Continue with available chains only
        chain_ids = [c for c in chain_ids if c in available_chains]
        if not chain_ids:
            raise ValueError(f"None of the requested chains found in {cif_file}")

    io = PDBIO()  # PDBIO can save in PDB format even from CIF input
    io.set_structure(structure)

    selector = MultiChainSelect(chain_ids)
    io.save(output_file, select=selector)
    
    # Verify output file was created and has content
    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
        raise IOError(f"Failed to create output file {output_file}")
    
    chains_str = ','.join(chain_ids)
    print(f"    Extracted chain(s) {chains_str} from {Path(cif_file).name} (CIF)")

def get_chain_reference_structure(base_path, organism_name):
    """
    For each locus_tag, get the reference structure chain(s) and save to a separate PDB file.
    
    Handles both:
    - Single chains (homooligomers or monomers): extracts one chain
    - Multiple chains (heterooligomers): extracts ALL functional subunits
    
    IMPORTANT: Only ONE structure should be marked as reference per locus_tag.
    This function extracts the chain(s) from that single reference structure.

    :param base_path: Base path.
    :param organism_name: Name of organism.
    """
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    
    all_locus_tags = metadata.ref_gbk_locus(base_path, organism_name)
    
    for locus_tag in all_locus_tags:
        locus_dir = os.path.join(structure_dir, locus_tag)
        
        if not os.path.isdir(locus_dir):
            continue
            
        summary_table_path = os.path.join(locus_dir, f"{locus_tag}_structure_summary.tsv")

        if files.file_check(summary_table_path):
            # Read structure_id as string to prevent scientific notation (e.g., 3E59 -> 3e+59)
            summary_df = pd.read_csv(summary_table_path, sep='\t', dtype={'structure_id': str})

            ref_rows = summary_df[summary_df['is_reference'] == True]

            if ref_rows.empty:
                print(f'  Warning: No reference structure found for locus_tag {locus_tag}.')
                continue
            
            if len(ref_rows) > 1:
                print(f'  ERROR: Multiple reference structures found for {locus_tag}. There should be only ONE.')
                print(f'  Found {len(ref_rows)} reference structures. Using first one.')
                ref_row = ref_rows.iloc[0]
            else:
                ref_row = ref_rows.iloc[0]
            
            struct_type = ref_row['structure_type']
            struct_id = ref_row['structure_id']
            chain_field = ref_row['chain']
            uniprot_id = ref_row['uniprot_id']
            
            # Parse chain field - can be single "A" or multiple "A;B;C" for heterooligomers
            if not chain_field or pd.isna(chain_field):
                chain_ids = ['A']
                print(f'  Warning: No chain ID specified for {struct_id} in {locus_tag}, defaulting to chain A.')
            elif ';' in str(chain_field):
                # Heterooligomer - multiple chains
                chain_ids = [c.strip() for c in str(chain_field).split(';')]
                print(f'  Heterooligomer detected for {locus_tag}: extracting chains {chain_ids}')
            else:
                # Single chain
                chain_ids = [str(chain_field)]
            
            uniprot_dir = os.path.join(locus_dir, uniprot_id)
            
            # Check if uniprot_dir exists
            if not os.path.isdir(uniprot_dir):
                print(f'  ERROR: UniProt directory {uniprot_dir} not found for {locus_tag}.')
                continue

            # Only extract chain(s) for PDB structures (AlphaFold doesn't need extraction)
            if struct_type == 'PDB':
                pdb_file_path = os.path.join(uniprot_dir, f"PDB_{struct_id}.pdb")
                cif_file_path = os.path.join(uniprot_dir, f"PDB_{struct_id}.cif")
                
                # Create filename with chain notation
                # Limit filename length for many chains
                if len(chain_ids) > 1:
                    if len(chain_ids) <= 5:  # Reasonable limit
                        chains_suffix = '_'.join(chain_ids)
                    else:
                        # For structures with many chains, use compact notation
                        chains_suffix = f"{chain_ids[0]}_to_{chain_ids[-1]}_multi"
                else:
                    chains_suffix = chain_ids[0]
                    
                output_chain_file = os.path.join(uniprot_dir, f"PDB_{struct_id}_{chains_suffix}_ref.pdb")
                
                # Try PDB format first, then CIF format
                # Verify files exist AND have content (not empty/corrupt)
                try:
                    pdb_valid = files.file_check(pdb_file_path) and os.path.getsize(pdb_file_path) > 0
                except OSError:
                    pdb_valid = False
                
                try:
                    cif_valid = files.file_check(cif_file_path) and os.path.getsize(cif_file_path) > 0
                except OSError:
                    cif_valid = False
                
                if pdb_valid:
                    if not files.file_check(output_chain_file):
                        print(f'  Extracting chain(s) {chain_ids} from {struct_id} for {locus_tag}')
                        try:
                            extract_chain_from_pdb(pdb_file_path, chain_ids, output_chain_file)
                        except Exception as e:
                            print(f'  ERROR: Failed to extract chain(s) {chain_ids} from {pdb_file_path}: {e}')
                    else:
                        print(f'  Chain file already exists: {output_chain_file}')
                elif cif_valid:
                    # Use CIF file if PDB not available (common for cryo-EM)
                    if not files.file_check(output_chain_file):
                        print(f'  Extracting chain(s) {chain_ids} from {struct_id} CIF for {locus_tag}')
                        try:
                            extract_chain_from_cif(cif_file_path, chain_ids, output_chain_file)
                        except Exception as e:
                            print(f'  ERROR: Failed to extract chain(s) {chain_ids} from {cif_file_path}: {e}')
                    else:
                        print(f'  Chain file already exists: {output_chain_file}')
                else:
                    print(f'  ERROR: Neither PDB nor CIF file found (or files are empty) for {struct_id} in locus_tag {locus_tag}.')
                    if files.file_check(pdb_file_path):
                        try:
                            size = os.path.getsize(pdb_file_path)
                            print(f'    Note: {pdb_file_path} exists but is empty (size: {size} bytes)')
                        except OSError:
                            print(f'    Note: {pdb_file_path} exists but size cannot be determined')
                    if files.file_check(cif_file_path):
                        try:
                            size = os.path.getsize(cif_file_path)
                            print(f'    Note: {cif_file_path} exists but is empty (size: {size} bytes)')
                        except OSError:
                            print(f'    Note: {cif_file_path} exists but size cannot be determined')
            
            elif struct_type == 'AlphaFold':
                af_file_path = os.path.join(uniprot_dir, f"AF_{uniprot_id}.pdb")
                if not files.file_check(af_file_path):
                    print(f'  ERROR: AlphaFold file {af_file_path} not found for locus_tag {locus_tag}.')
        else:
            print(f'  Structure summary table not found for {locus_tag}, skipping chain extraction.')

##  ------------------- FPocket functions ------------------- ##
def FPocket_models(directory):
    """
    Run Fpocket for all .pdb in a directory.
    Results are saved inside each fpocket output directory as all_pockets.json.

    :param directory: Directory with the models.
    """

    if os.path.exists(directory):
        for protein in tqdm(os.listdir(directory), desc="Processing", unit="protein"):
            protein_path = os.path.join(directory, protein)
            for template in os.listdir(protein_path):
                template_path = os.path.join(protein_path, template)
                for model_file in os.listdir(template_path):
                    model_pdb_path = os.path.join(template_path, model_file)
                    if os.path.isfile(model_pdb_path) and model_file.endswith('.pdb'):
                        fpocket_outdir = os.path.splitext(model_pdb_path)[0] + "_fpocket"
                        if not os.path.exists(fpocket_outdir):
                            print(f'Run Fpocket for {model_pdb_path}')
                            programs.run_fpocket(template_path, model_pdb_path)
                            print(f'FPocket prediction for {model_file} in {fpocket_outdir}.')
                            pockets_dict = pockets_data_to_dict(fpocket_outdir)
                            pockets_filter(pockets_dict)
    else:
        print(f"The directory '{directory}' not found.", file=sys.stderr)

def fpocket_for_structure(pdb, output_path, pockets_dir):
    """
    Run Fpocket for a single PDB structure.

    This function uses the `run_fpocket` function from the `programs` module.
    Results are generated by Fpocket in a directory named '<pdb_basename>_fpocketc'
    inside `output_path`, and then moved into `pockets_dir`.

    :param pdb: Path to the PDB file.
    :param output_path: Directory where Fpocket will write its output.
    :param pockets_dir: Directory where the final pockets results will be stored.
    """

    fpocket_outdir = os.path.basename(os.path.splitext(pdb)[0]) + "_fpocket"
    results_path = os.path.join(pockets_dir, fpocket_outdir)

    if not os.path.exists(results_path):
        print(f'Running Fpocket for {pdb}')
        programs.run_fpocket(output_path, pdb)

        shutil.move(os.path.join(output_path, fpocket_outdir), pockets_dir)

        print(f'FPocket prediction for {pdb} completed.')
    else:
        print(f'Fpocket results for {pdb} already present.')


def find_structures_for_locus(locus_dir):
    """
    Find ONLY the reference structure for a given locus_tag directory.
    
    Selection rules:
    1) Look for PDB reference structures (*_ref.pdb)
       (These are PDB chains extracted by extract_chain_from_pdb)
    2) If no '_ref.pdb' found, use AlphaFold models: files starting with 'AF_'
       (AlphaFold files don't need chain extraction)
    3) If nothing found, return None
    
    :param locus_dir: Path to the locus_tag directory.
    :return: Path to reference PDB file, or None.
    """
    # 1) First priority: Reference structures (*_ref.pdb)
    # Exclude folder pockets
    ref_files = glob.glob(os.path.join(locus_dir, '**', '*_ref.pdb'), recursive=True)
    ref_files = [f for f in ref_files if 'pockets' not in f.split(os.sep)]
    ref_files = sorted(set(ref_files))
    
    if len(ref_files) > 1:
        print(f"  ERROR: Found {len(ref_files)} '*_ref.pdb' files in {locus_dir}. There should be only ONE reference.")
        print(f"  Files found: {ref_files}")
        print(f"  Using first file: {ref_files[0]}")
        return ref_files[0]
    elif len(ref_files) == 1:
        print(f"  Using PDB reference: {ref_files[0]}")
        return ref_files[0]
    
    # 2) Second priority: AlphaFold models (if no PDB ref exists)
    # Exclude folder pockets
    af_pdbs = glob.glob(os.path.join(locus_dir, '**', 'AF_*.pdb'), recursive=True)
    af_pdbs = [f for f in af_pdbs if 'pockets' not in f.split(os.sep)]
    af_pdbs = sorted(set(af_pdbs))
    
    if len(af_pdbs) > 1:
        print(f"  Warning: Found {len(af_pdbs)} AlphaFold models in {locus_dir}. Using first one.")
        print(f"  Using AlphaFold: {af_pdbs[0]}")
        return af_pdbs[0]
    elif len(af_pdbs) == 1:
        #print(f"  Using AlphaFold model: {af_pdbs[0]}")
        return af_pdbs[0]
    
    # 3) No reference structure found
    print(f"  Warning: No reference structure (*_ref.pdb) or AlphaFold (AF_*.pdb) found in {locus_dir}. Skipping FPocket.")
    return None

def pockets_finder_for_locus(locus_dir):
    """
    Run Fpocket for the selected structures in a locus_tag directory.

    This function:
      - Creates a 'pockets' directory inside the locus_tag directory.
      - Selects the structures to run Fpocket on, according to:
            * Reference PDB chains '*_ref.pdb' if any.
            * Else, AlphaFold PDB files (heuristic by name).
            * Else, all '.pdb' as a last fallback.
      - Runs Fpocket in parallel using ProcessPoolExecutor.

    :param locus_dir: Path to the locus_tag directory.

    """

    if os.path.exists(locus_dir):

        pockets_dir = os.path.join(locus_dir, 'pockets')
        os.makedirs(pockets_dir, exist_ok=True)

        # Select structure to process
        pdb_file = find_structures_for_locus(locus_dir)
        if not pdb_file:
            print(f"No structures to process in {locus_dir}")
            return

        print(f"Running Fpocket on {pdb_file} structures in {locus_dir}")
        
        pdb_parent_dir = os.path.dirname(pdb_file)
        fpocket_for_structure(pdb_file, pdb_parent_dir, pockets_dir)

    else:
        print(f"The directory '{locus_dir}' was not found.", file=sys.stderr)


def pockets_finder_for_all_loci(base_path, organism_name):
    """
    Run Fpocket for all locus_tag directories under the 'structures' folder.

    :param base_path: Base path.
    :param organism_name: Name of organism.

    """

    structures_dir = os.path.join(base_path, 'organism', organism_name, "structures")
    
    if not os.path.exists(structures_dir):
        print(f"The directory '{structures_dir}' was not found.", file=sys.stderr)
        return

    all_locus = metadata.ref_gbk_locus(base_path, organism_name)

    for locus_tag in tqdm(all_locus, desc='Locus tags'):
        locus_dir = os.path.join(structures_dir, locus_tag)
        if os.path.exists(locus_dir):
            pockets_finder_for_locus(locus_dir)
        else:
            print(f"Locus directory '{locus_dir}' does not exist, skipping.")

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
    if os.path.exists(directory):
        pockets_file = os.path.join(directory, 'all_pockets.json')

        if not files.file_check(pockets_file):
            print('Parsing pockets results.')
            pockets_dict = {}
            
            for root, dirs, list_files in os.walk(directory):
                for file in list_files:
                    if file.endswith('_info.txt'):
                        file_path = os.path.join(root, file)
                        
                        # Extract structure ID from folder name
                        # E.g., "PDB_1ABC_A_ref_fpocket" -> "1ABC"
                        # or "AF_P12345_fpocket" -> "P12345"
                        folder_name = os.path.basename(root)
                        if folder_name.endswith('_fpocket'):
                            folder_name = folder_name[:-8]  # Remove '_fpocket'
                        
                        # Parse structure ID
                        parts = folder_name.split('_')
                        if parts[0] == 'PDB' and len(parts) >= 2:
                            protein_ID = parts[1]  # Get PDB code (e.g., "1ABC")
                        elif parts[0] == 'AF' and len(parts) >= 2:
                            protein_ID = '_'.join(parts[1:])  # Get full UniProt ID (e.g., "P12345" or "P12345-F1")
                        else:
                            protein_ID = folder_name
                        
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
        return {}

def pockets_filter(pockets_dict):
    """
    Filter the pockets data to keep only the one with the highest Druggability Score.
    Returns a dictionary with the filtered data.

    :param pockets_dict: Dictionary with the pockets data.

    :return: Dictionary with the best pockets data. 
    """

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

    return DS_dict

##  ------------------- P2Rank functions ------------------- ##
def p2rank_for_structure(pdb, output_path, p2rank_dir, cpus, alphafold=False):
    """
    Run P2Rank for a single PDB structure.

    This function uses the `run_p2rank` function from the `programs` module.
    Results are generated by P2Rank in a directory inside `output_path`, 
    and then moved into `p2rank_dir`.

    :param pdb: Path to the PDB file.
    :param output_path: Directory where P2Rank will write its output.
    :param p2rank_dir: Directory where the final P2Rank results will be stored.
    :param cpus: Number of CPUs/threads to use.
    :param alphafold: Boolean, if True adds '-c alphafold' to the command.
    """
    pdb_basename = os.path.basename(os.path.splitext(pdb)[0])
    
    otput_folder = pdb_basename + "_p2rank"
    
    p2rank_output_dir_tmp = os.path.join(output_path, otput_folder)
    p2rank_output_dir_final = os.path.join(p2rank_dir, otput_folder)

    if os.path.exists(p2rank_output_dir_tmp):
        shutil.rmtree(p2rank_output_dir_tmp)

    if not os.path.exists(p2rank_output_dir_final) or not os.listdir(p2rank_output_dir_final):       

        try:
    
            print(f'Running P2Rank for {pdb}')
            # Run the p2rank Docker command
            programs.run_p2rank(output_path, pdb, cpus, alphafold=alphafold)

            # Move the output directory to the final destination
            shutil.move(p2rank_output_dir_tmp, p2rank_dir)
            print(f'P2Rank prediction for {pdb} completed.')
        except Exception as e:
            print(f'Error running P2Rank for {pdb}: {e}')
    else:
        print(f'P2Rank results directory for {pdb} already present.')


def p2rank_finder_for_locus(locus_dir, cpus):
    """
    Run P2Rank for the selected structures in a locus_tag directory.

    This function:
      - Creates a 'pockets' directory inside the locus_tag directory.
      - Selects the structures to run P2Rank on, according to:
            * Reference PDB chains '*_ref.pdb' if any.
            * Else, AlphaFold PDB files (heuristic by name).
            * Else, all '.pdb' as a last fallback.
      - Runs P2Rank with appropriate parameters.

    :param locus_dir: Path to the locus_tag directory.
    :param cpus: Number of CPUs/threads to use.
    """

    if os.path.exists(locus_dir):

        p2rank_dir = os.path.join(locus_dir, 'pockets')
        os.makedirs(p2rank_dir, exist_ok=True)

        # Select structure to process
        pdb_file = find_structures_for_locus(locus_dir)
        if not pdb_file:
            print(f"No structures to process in {locus_dir}")
            return

        # Determine if it's an AlphaFold structure
        is_alphafold = 'AF_' in os.path.basename(pdb_file)
        
        print(f"Running P2Rank on {pdb_file} in {locus_dir} (AlphaFold={is_alphafold})")
        
        pdb_parent_dir = os.path.dirname(pdb_file)
        p2rank_for_structure(pdb_file, pdb_parent_dir, p2rank_dir, cpus, alphafold=is_alphafold)

    else:
        print(f"The directory '{locus_dir}' was not found.", file=sys.stderr)

def p2rank_finder_for_all_loci(base_path, organism_name, cpus):
    """
    Run P2Rank for all locus_tag directories under the 'structures' folder.

    :param base_path: Base path.
    :param organism_name: Name of organism.
    :param cpus: Number of CPUs/threads to use.
    """

    structures_dir = os.path.join(base_path, 'organism', organism_name, "structures")
    
    if not os.path.exists(structures_dir):
        print(f"The directory '{structures_dir}' was not found.", file=sys.stderr)
        return

    all_locus = metadata.ref_gbk_locus(base_path, organism_name)

    for locus_tag in tqdm(all_locus, desc='Locus tags'):
        locus_dir = os.path.join(structures_dir, locus_tag)
        if os.path.exists(locus_dir):
            p2rank_finder_for_locus(locus_dir, cpus)
        else:
            print(f"Locus directory '{locus_dir}' does not exist, skipping.")


def p2rank_parse_predictions(csv_file):
    """
    Parse P2Rank predictions CSV file to extract pocket information.
    
    The CSV file contains columns including 'name' (pocket number) and 'probability'
    (probability of being a ligand-binding site).
    
    :param csv_file: Path to the *_predictions.csv file from P2Rank.
    :return: Dictionary with pocket information, or 'No_pockets' if no data found.
    """
    if not os.path.exists(csv_file):
        return 'No_pockets'
    
    try:
        # Use skipinitialspace=True to strip leading whitespace after commas
        df = pd.read_csv(csv_file, sep=',', skipinitialspace=True)
        
        # Strip whitespace from column names
        df.columns = df.columns.str.strip()
        
        # Check if required columns exist
        if 'name' not in df.columns or 'probability' not in df.columns:
            print(f"Warning: Required columns 'name' and/or 'probability' not found in {csv_file}")
            return 'No_pockets'
        
        # Check if there are any pockets
        if len(df) == 0:
            return 'No_pockets'
        
        # Create dictionary with pocket data
        pockets_dict = {}
        for _, row in df.iterrows():
            pocket_name = f"Pocket {row['name']}"
            pockets_dict[pocket_name] = {
                'probability': float(row['probability'])
            }
        
        return pockets_dict
        
    except Exception as e:
        print(f"Error parsing P2Rank predictions file {csv_file}: {e}")
        return 'No_pockets'


def p2rank_data_to_dict(directory):
    """
    Parse P2Rank results from a directory to extract pocket predictions.
    
    This function processes P2Rank output directories and creates a dictionary
     with pocket information.
    
    :param directory: Path to the directory containing P2Rank results.
    :return: Dictionary with P2Rank pocket data for each structure.
    """
    if not os.path.exists(directory):
        print(f'Error: {directory} not found.')
        return {}
    
    p2rank_file = os.path.join(directory, 'all_p2rank_pockets.json')
    
    if not files.file_check(p2rank_file):
        print('Parsing P2Rank results.')
        p2rank_dict = {}
        
        # Find all P2Rank output directories (ending with _p2rank)
        for item in os.listdir(directory):
            item_path = os.path.join(directory, item)
            
            if os.path.isdir(item_path) and item.endswith('_p2rank'):
                # Extract structure ID from folder name
                # E.g., "PDB_1ABC_A_ref_p2rank" -> "1ABC"
                # or "AF_P12345_p2rank" -> "P12345"
                folder_name = item.replace('_p2rank', '')  # Remove '_p2rank' (7 characters)
                
                # Parse structure ID
                parts = folder_name.split('_')
                if parts[0] == 'PDB' and len(parts) >= 2:
                    protein_ID = parts[1]  # Get PDB code (e.g., "1ABC")
                elif parts[0] == 'AF' and len(parts) >= 2:
                    protein_ID = '_'.join(parts[1:])  # Get full UniProt ID (e.g., "P12345" or "P12345-F1")
                else:
                    protein_ID = folder_name
                
                # Find the predictions CSV file
                predictions_file = None
                for file in os.listdir(item_path):
                    if file.endswith('_predictions.csv'):
                        predictions_file = os.path.join(item_path, file)
                        break
                
                if predictions_file:
                    p2rank_dict[protein_ID] = p2rank_parse_predictions(predictions_file)
                else:
                    print(f"Warning: No predictions CSV file found in {item_path}")
                    p2rank_dict[protein_ID] = 'No_pockets'
        
        files.dict_to_json(directory, 'all_p2rank_pockets.json', p2rank_dict)
        print(f'All P2Rank pocket data saved to {p2rank_file}.')
    else:
        p2rank_dict = files.json_to_dict(p2rank_file)
        print(f'All P2Rank pocket data loaded from {p2rank_file}.')
    
    return p2rank_dict


def p2rank_filter(p2rank_dict):
    """
    Filter P2Rank pockets to keep only the one with the highest probability.
    
    This function returns a dictionary
    with the best pocket (highest probability) for each structure.
    
    :param p2rank_dict: Dictionary with P2Rank pocket data.
    :return: Dictionary with the best pocket data (highest probability).
    """
    best_pockets_dict = {}
    
    for protein_id, pockets in p2rank_dict.items():
        if pockets == 'No_pockets':
            best_pockets_dict[protein_id] = 'No_pockets'
        else:
            max_prob = 0
            max_pocket = None
            
            for pocket_name, props in pockets.items():
                probability = props.get('probability', 0)
                if probability > max_prob:
                    max_prob = probability
                    max_pocket = pocket_name
            
            best_pockets_dict[protein_id] = {
                'max_probability': max_prob,
                'pocket': max_pocket
            }
    
    return best_pockets_dict


def merge_structure_data (base_path, organism_name):
    """
    Merge all the structure data in a single dictionary.
    Saves the dictionary in a .json file named using the organism name followed by '_structure_data.json' in the 'structures' directory.
    Returns a dictionary with the merged data.

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.

    :return: Dictionary with the merged data.
    """

    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    merged_file = os.path.join(structure_dir, f'{organism_name}_structure_data.json')

    all_locus_tags = metadata.ref_gbk_locus(base_path, organism_name)

    if not files.file_check(merged_file):

        print('Merging all structure data.')

        merged_dict = {}

        for locus_tag in tqdm(all_locus_tags, desc='Locus tags'):

            locus_dir = os.path.join(structure_dir, locus_tag)
            
            # Skip if locus directory doesn't exist
            if not os.path.isdir(locus_dir):
                continue

            # Initialize entry for this locus_tag
            if locus_tag not in merged_dict:
                merged_dict[locus_tag] = {}

            # ========== FPocket Data ==========
            pockets_dir = os.path.join(locus_dir, 'pockets')
            
            if os.path.isdir(pockets_dir):
                best_pockets_file = os.path.join(pockets_dir, f'best_DS_fpocket_{locus_tag}.json')

                # Obtain best pockets data
                if not files.file_check(best_pockets_file):

                    # Fpocket
                    all_fpocket_folders = [d for d in os.listdir(pockets_dir) if d.endswith('_fpocket') and os.path.isdir(os.path.join(pockets_dir, d))]

                    if len(all_fpocket_folders) == 1:
                        fpocket_folder = os.path.join(pockets_dir, all_fpocket_folders[0])
                        pockets_dict = pockets_data_to_dict(fpocket_folder)
                        best_pockets_dict = pockets_filter(pockets_dict)
                        files.dict_to_json(pockets_dir, f'best_DS_fpocket_{locus_tag}.json', best_pockets_dict)

                        merged_dict[locus_tag]['fpocket_best_pockets'] = best_pockets_dict

                    elif len(all_fpocket_folders) > 1:
                        print(f'Multiple Fpocket output folders found in {pockets_dir}, skipping best pockets extraction for {locus_tag}.')
                    else:
                        print(f'No Fpocket output folder found in {pockets_dir}, skipping best pockets extraction for {locus_tag}.')
                
                else:
                    best_pockets_dict = files.json_to_dict(best_pockets_file)
                    merged_dict[locus_tag]['fpocket_best_pockets'] = best_pockets_dict

            # ========== P2Rank Data ==========
            p2rank_dir = os.path.join(locus_dir, 'pockets')
            
            if os.path.isdir(p2rank_dir):
                best_p2rank_file = os.path.join(p2rank_dir, f'best_prob_p2rank_{locus_tag}.json')

                # Obtain best P2Rank pockets data
                if not files.file_check(best_p2rank_file):

                    # P2Rank
                    all_p2rank_folders = [d for d in os.listdir(p2rank_dir) if d.endswith('_p2rank') and os.path.isdir(os.path.join(p2rank_dir, d))]

                    if len(all_p2rank_folders) >= 1:
                        # Process P2Rank results
                        p2rank_dict = p2rank_data_to_dict(p2rank_dir)
                        best_p2rank_dict = p2rank_filter(p2rank_dict)
                        files.dict_to_json(p2rank_dir, f'best_prob_p2rank_{locus_tag}.json', best_p2rank_dict)

                        merged_dict[locus_tag]['p2rank_best_pockets'] = best_p2rank_dict

                    else:
                        print(f'No P2Rank output folder found in {p2rank_dir}, skipping P2Rank extraction for {locus_tag}.')
                
                else:
                    best_p2rank_dict = files.json_to_dict(best_p2rank_file)
                    merged_dict[locus_tag]['p2rank_best_pockets'] = best_p2rank_dict

        files.dict_to_json(structure_dir, f'{organism_name}_structure_data.json', merged_dict)
        return merged_dict
    
    else:
        merged_dict = files.json_to_dict(merged_file)
        print(f'Merged structure data in {merged_file}.')
        return merged_dict

def final_structure_table(base_path, organism_name):
    """
    Create a final summary table with structure and pocket information for each locus_tag.
    The table columns are ['gene', 'uniprot', 'druggability_score', 'pocket', 'structure']
    
    :param base_path: Base path.
    :param organism_name: Name of organism.
    :return: DataFrame with the final summary table (ONE row per locus_tag).
    """
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    merged_file = os.path.join(structure_dir, f'{organism_name}_structure_data.json')
    final_table_file = os.path.join(structure_dir, f'{organism_name}_final_structure_summary.tsv')

    if not files.file_check(final_table_file):
        
        if not files.file_check(merged_file):
            merged_dict = merge_structure_data(base_path, organism_name)
        else:
            merged_dict = files.json_to_dict(merged_file)

        rows = []
        all_locus_tags = metadata.ref_gbk_locus(base_path, organism_name)

        for locus_tag in all_locus_tags:
            # Get UniProt ID and structure ID from REFERENCE structure only
            structure_summary_path = os.path.join(
                structure_dir, locus_tag, f"{locus_tag}_structure_summary.tsv"
            )
            
            uniprot_id = None
            structure_id = None
            structure_type = None
            
            if files.file_check(structure_summary_path):
                # Read structure_id as string to prevent scientific notation (e.g., 3E59 -> 3e+59)
                struct_df = pd.read_csv(structure_summary_path, sep='\t', dtype={'structure_id': str})
                ref_rows = struct_df[struct_df['is_reference'] == True]
                
                if not ref_rows.empty:
                    if len(ref_rows) > 1:
                        print(f'  Warning: Multiple reference structures for {locus_tag}, using first.')
                    ref_row = ref_rows.iloc[0]
                    uniprot_id = ref_row['uniprot_id']
                    structure_id = ref_row['structure_id']
                    structure_type = ref_row['structure_type']
            
            # Get pocket data for the reference structure
            data = merged_dict.get(locus_tag, {})
            
            # ========== FPocket Data ==========
            fpocket_data = data.get('fpocket_best_pockets', {})
            
            druggability_score = None
            fpocket_pocket = None
            
            if fpocket_data and structure_id:
                # For AlphaFold structures, use UniProt ID as key
                # For PDB structures, use PDB code as key
                if structure_type == 'AlphaFold':
                    pocket_key = uniprot_id
                else:
                    pocket_key = structure_id
                
                pocket_info = fpocket_data.get(pocket_key)
                if pocket_info:
                    if pocket_info != 'No_pockets':
                        druggability_score = pocket_info.get('maxDS')
                        fpocket_pocket = pocket_info.get('pocket')
                    else:
                        fpocket_pocket = 'No_pockets'
            
            # ========== P2Rank Data ==========
            p2rank_data = data.get('p2rank_best_pockets', {})
            
            p2rank_probability = None
            p2rank_pocket = None
            
            if p2rank_data and structure_id:
                # For AlphaFold structures, use UniProt ID as key
                # For PDB structures, use PDB code as key
                if structure_type == 'AlphaFold':
                    pocket_key = uniprot_id
                else:
                    pocket_key = structure_id
                
                pocket_info = p2rank_data.get(pocket_key)
                if pocket_info:
                    if pocket_info != 'No_pockets':
                        p2rank_probability = pocket_info.get('max_probability')
                        p2rank_pocket = pocket_info.get('pocket')
                    else:
                        p2rank_pocket = 'No_pockets'
            
            # Add ONE row per locus_tag
            rows.append({
                'gene': locus_tag,
                'uniprot': uniprot_id,
                'structure': structure_id,
                'druggability_score': druggability_score,
                'fpocket_pocket': fpocket_pocket,
                'p2rank_probability': p2rank_probability,
                'p2rank_pocket': p2rank_pocket
            })

        final_df = pd.DataFrame(rows)
        # Ensure structure column is stored as string to prevent scientific notation issues
        if 'structure' in final_df.columns:
            final_df['structure'] = final_df['structure'].astype(str)
        final_df.to_csv(final_table_file, sep='\t', index=False)
        print(f'\nFinal structure summary table saved to {final_table_file}')
        print(f'  Total genes: {len(final_df)}')
        print(f'  Genes with structures: {final_df["structure"].notna().sum()}')
        print(f'  Genes with FPocket pockets: {final_df["fpocket_pocket"].notna().sum()}')
        print(f'  Genes with P2Rank pockets: {final_df["p2rank_pocket"].notna().sum()}')
        return final_df
    
    else:
        # Read structure column as string to prevent scientific notation (e.g., 3E59 -> 3e+59)
        final_df = pd.read_csv(final_table_file, sep='\t', dtype={'structure': str})
        print(f'Final structure summary table loaded from {final_table_file}')
        return final_df
    
def pipeline_structures(base_path, organism_name, specie_taxid, strain_taxid, cpus=multiprocessing.cpu_count(), resolution_cutoff = 3.5):
    """
    Complete pipeline to obtain and process structures for drug target identification.
    
    This pipeline:
    1. Downloads and processes UniProt proteome data
    2. Maps organism genes to UniProt IDs via BLAST
    3. Downloads PDB and AlphaFold structures
    4. Runs FPocket and P2Rank for druggable pocket detection
    5. Generates final summary tables
    
    :param base_path: Base path to project directory.
    :param organism_name: Name of organism.
    :param specie_taxid: Species-level NCBI Taxonomy ID (e.g., 287 for P. aeruginosa).
    :param strain_taxid: Strain-level NCBI Taxonomy ID (e.g., 208964 for PAO1).
    :param cpus: Number of CPU cores to use. Default is all available cores.
    
    :return: DataFrame with final structure summary table.
    """
    
    print(f'\n{"="*80}')
    print(f'FASTTARGET STRUCTURE ANALYSIS PIPELINE')
    print(f'Organism: {organism_name}')
    print(f'Species TaxID: {specie_taxid} | Strain TaxID: {strain_taxid}')
    print(f'Using {cpus} CPU cores')
    print(f'{"="*80}\n')
    
    # ========== STAGE 1: UniProt Proteome Acquisition ==========
    print(f'\n{"─"*80}')
    print(f'STAGE 1: UNIPROT PROTEOME ACQUISITION AND MAPPING')
    print(f'{"─"*80}')
    
    try:
        print(f'\n[1.1] Downloading UniProt species data (TaxID: {specie_taxid})...')
        download_species_uniprot_data(base_path, organism_name, specie_taxid)
        
        print(f'\n[1.2] Parsing UniProt data into FASTA files...')
        uniprot_dir = os.path.join(base_path, 'organism', organism_name, 'structures', 'uniprot_files')
        uniprot_file = os.path.join(uniprot_dir, f"uniprot_specie_taxid_{specie_taxid}_data.tsv")
        parse_uniprot_species_data(uniprot_file, specie_taxid, strain_taxid)
        
        print(f'\n[1.3] Clustering species proteome with CD-HIT...')
        cluster_uniprot_specie(base_path, organism_name, specie_taxid)
        
        print(f'\n[1.4] Creating BLAST databases...')
        create_uniprot_blast_db(base_path, organism_name, specie_taxid, strain_taxid)
        
        print(f'\n[1.5] Running BLAST searches against UniProt databases...')
        uniprot_proteome_blast(base_path, organism_name, specie_taxid, strain_taxid, cpus=cpus)
        
        print(f'\n[1.6] Mapping organism genes to UniProt IDs...')
        mapping_dict = uniprot_proteome_mapping(base_path, organism_name, specie_taxid, strain_taxid)
        print(f'    ✓ Mapped {len(mapping_dict)} genes to UniProt IDs')
        
    except Exception as e:
        print(f'\n    ✗ ERROR in Stage 1: {e}')
        raise
    
    # ========== STAGE 2: Structure Acquisition ==========
    print(f'\n{"─"*80}')
    print(f'STAGE 2: STRUCTURE DOWNLOAD AND ORGANIZATION')
    print(f'{"─"*80}')
    
    try:
        print(f'\n[2.1] Creating directory structure for each gene...')
        create_subfolder_structures(base_path, organism_name)
        
        print(f'\n[2.2] Generating structure summary tables...')
        create_summary_structure_file(base_path, organism_name, resolution_cutoff=resolution_cutoff)
        
        print(f'\n[2.3] Downloading PDB and AlphaFold structures...')
        download_structures(base_path, organism_name)
        
        print(f'\n[2.4] Extracting reference structure chains...')
        get_chain_reference_structure(base_path, organism_name)
        
    except Exception as e:
        print(f'\n    ✗ ERROR in Stage 2: {e}')
        raise
    
    # ========== STAGE 3: Pocket Detection ==========
    print(f'\n{"─"*80}')
    print(f'STAGE 3: DRUGGABLE POCKET DETECTION')
    print(f'{"─"*80}')
    
    try:
        print(f'\n[3.1] Running FPocket for all structures...')
        structures_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
        pockets_finder_for_all_loci(base_path, organism_name)
        
        print(f'\n[3.2] Running P2Rank for all structures...')
        p2rank_finder_for_all_loci(base_path, organism_name, cpus)
        
        print(f'\n[3.3] Merging structure and pocket data...')
        merged_data = merge_structure_data(base_path, organism_name)
        print(f'    ✓ Processed {len(merged_data)} genes')
        
        print(f'\n[3.4] Creating final summary table...')
        final_df = final_structure_table(base_path, organism_name)
        
    except Exception as e:
        print(f'\n    ✗ ERROR in Stage 3: {e}')
        raise
    
    # ========== Pipeline Complete ==========
    print(f'\n{"="*80}')
    print(f'Pipeline finished')
    print(f'{"="*80}')
    
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    final_table_path = os.path.join(structure_dir, f'{organism_name}_final_structure_summary.tsv')
    
    print(f'\nResults saved to: {final_table_path}')
    print(f'Structure data: {os.path.join(structure_dir, f"{organism_name}_structure_data.json")}')
    
    return final_df


def get_reference_structure_path(base_path, organism_name, locus_tag):
    """
    Get the file path of the reference structure for a given locus_tag.
    
    This is a convenience wrapper around find_structures_for_locus() that takes
    base_path, organism_name, and locus_tag as parameters.
    
    Selection priority:
    1) PDB reference structures (*_ref.pdb) - extracted chains from PDB
    2) AlphaFold models (AF_*.pdb) - full AlphaFold predictions
    3) None if no reference structure found
    
    :param base_path: Base path to project directory.
    :param organism_name: Name of organism.
    :param locus_tag: Locus tag identifier.
    
    :return: Absolute path to reference structure file, or None if not found.
    """
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    locus_dir = os.path.join(structure_dir, locus_tag)
    
    # Check if locus directory exists
    if not os.path.isdir(locus_dir):
        return None
    
    # Use existing function to find reference structure
    reference_path = find_structures_for_locus(locus_dir)
    
    return reference_path


def get_all_reference_structures(base_path, organism_name, path_mode=True):
    """
    Get a dictionary mapping all locus_tags to their reference structure paths.
    
    This function iterates through all locus_tag directories and finds their
    reference structures, returning a complete mapping.
    
    :param base_path: Base path to project directory.
    :param organism_name: Name of organism.
    
    :return: Dictionary with locus_tag as key and reference structure path as value.
             Locus tags without structures will have None as value.
    """
    structure_dir = os.path.join(base_path, 'organism', organism_name, 'structures')
    
    # Get all locus_tags from genome
    all_locus_tags = metadata.ref_gbk_locus(base_path, organism_name)
    
    reference_structures = {}
    
    for locus_tag in all_locus_tags:
        if path_mode:
            ref_path = get_reference_structure_path(base_path, organism_name, locus_tag)
            reference_structures[locus_tag] = ref_path
        else:
            #get the name of the reference structure only
            ref_path = get_reference_structure_path(base_path, organism_name, locus_tag)
            if ref_path:
                ref_name = os.path.basename(ref_path).split('.')[0]
                reference_structures[locus_tag] = ref_name
            else:
                reference_structures[locus_tag] = None

    return reference_structures

