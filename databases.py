#!/usr/bin/env python3
"""
FastTarget Databases Script

This script contains functions to download and prepare various biological databases,
including the human proteome, microbiome species catalogue, and DEG database.
It also includes functions to create BLAST and FOLDSEEK databases.

"""


import urllib.request
import gzip
import tarfile
import shutil
import os
import time
from ftscripts import programs,files, structures,logger
import tqdm
import requests
import multiprocessing
import urllib.request
from bs4 import BeautifulSoup
import re
import glob
import json
import zlib
from urllib.parse import urlparse, parse_qs, urlencode
from requests.adapters import HTTPAdapter, Retry
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import socket
import argparse
from functools import wraps  
from urllib3.exceptions import SSLError
from requests.exceptions import (
    SSLError as RequestsSSLError, 
    ConnectionError, 
    Timeout
)
import pandas as pd

def retry_on_timeout(func, *args, max_retries=3, delay=5, **kwargs):
    """
    Retry a function call if a timeout occurs.

    :param func: The function to call.
    :param max_retries: Maximum number of retries.
    :param delay: Delay between retries (in seconds).
    :param args: Positional arguments to pass to the function.
    :param kwargs: Keyword arguments to pass to the function.
    :return: The result of the function call.
    """
    attempt = 0
    while attempt < max_retries:
        try:
            return func(*args, **kwargs)
        except (requests.exceptions.Timeout, socket.timeout) as e:
            attempt += 1
            print(f'Timeout occurred: {e}. Retrying {attempt}/{max_retries}...')
            time.sleep(delay)  # Add a delay before retrying
    raise Exception(f'Failed after {max_retries} retries due to timeout.')

def retry_with_backoff(max_retries=5, initial_delay=2, backoff_factor=2, max_delay=60):
    """
    Enhanced retry decorator with exponential backoff for network/SSL errors.
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            delay = initial_delay
            last_exception = None
            
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                except (SSLError, RequestsSSLError, ConnectionError, Timeout, 
                        requests.exceptions.ChunkedEncodingError) as e:
                    last_exception = e
                    if attempt < max_retries - 1:
                        sleep_time = min(delay, max_delay)
                        logger.logger.warning(
                            f"Attempt {attempt + 1}/{max_retries} failed for {func.__name__}: {str(e)}. "
                            f"Retrying in {sleep_time}s..."
                        )
                        time.sleep(sleep_time)
                        delay *= backoff_factor
                    else:
                        logger.logger.error(
                            f"All {max_retries} attempts failed for {func.__name__}: {str(e)}"
                        )
            
            raise last_exception
        return wrapper
    return decorator

def download_with_progress(url, filepath, max_retries=3, backoff_factor=2):
    """
    Download a file from a URL with a progress indicator and retry logic.

    :param url: The URL to download from.
    :param filepath: The local path where the file should be saved.
    :param max_retries: Maximum number of attempts (default: 3).
    :param backoff_factor: Multiplier for exponential backoff (default: 2).
    """

    def reporthook(count, block_size, total_size):
        downloaded_mb = count * block_size / (1024 * 1024)
        if total_size > 0:
            total_mb = total_size / (1024 * 1024)
            pct = (downloaded_mb / total_mb) * 100
            print(f"\rDownloaded {downloaded_mb:.2f} / {total_mb:.2f} MB ({pct:.1f}%)",
                  end="", flush=True)
        else:
            print(f"\rDownloaded {downloaded_mb:.2f} MB",
                  end="", flush=True)

    for attempt in range(1, max_retries + 1):
        print(f"\nAttempt {attempt}/{max_retries}...")
        try:
            urllib.request.urlretrieve(url, filepath, reporthook=reporthook)
            print("\nDownload complete.")
            return True
        except Exception as e:
            print(f"\nError on attempt {attempt}: {e}")
            if attempt < max_retries:
                wait = backoff_factor ** (attempt - 1)
                print(f"Retrying in {wait} seconds...")
                time.sleep(wait)
            else:
                print("All retry attempts failed.")
                return False

def download_with_wget(url, filepath, max_retries=5):
    """
    Download a file using wget with progress, resume capability, and retry logic.
    
    :param url: The URL to download from.
    :param filepath: The local path where the file should be saved.
    :param max_retries: Number of retry attempts if download fails.
    """
    attempt = 0
    success = False
    
    while attempt < max_retries and not success:
        try:
            # Use wget with resume and progress reporting
            command = f'wget -c --progress=dot:mega -O {filepath} "{url}"'
            programs.run_bash_command(command)
            print('Download complete.')
            success = True  # Mark success if wget completes
        except Exception as e:
            attempt += 1
            print(f'Error downloading {url}: {e}')
            if attempt < max_retries:
                print(f'Retrying... ({attempt}/{max_retries})')
                time.sleep(2)  # Optional: wait before retrying
            else:
                print('Max retries reached. Download failed.')


def simple_concat_faa(input_dir: str, output_path: str):
    """
    Concatenate all .faa files from a directory into a single FASTA file.

    :param input_dir : Directory containing .faa files.
    :param output_path : Path to the output concatenated FASTA file.
    
    """
    faa_files = sorted(glob.glob(os.path.join(input_dir, "**", "*.faa"), recursive=True))

    if not faa_files:
        raise FileNotFoundError(f"No .faa files found in {input_dir}")

    with open(output_path, "w") as outfile:
        for fname in tqdm.tqdm(faa_files, desc="Concatenating .faa files", unit="file"):
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)


def batch_uniprot_mapping(source, dest, ids, max_retries=5, sleep_time=5):
    """
    Maps a list of UniProt IDs from a source database to a destination database using the UniProt API.
    More information: https://www.uniprot.org/help/id_mapping

    :param source (str): The source database (e.g., 'UniProtKB_AC-ID').
    :param dest (str): The destination database (e.g., 'PDB').
    :param ids (list): A list of UniProt IDs to map.
    :param max_retries (int): Number of retries in case of request failure.
    :param sleep_time (int): Number of seconds to wait between retries.
    
    :return A dictionary where keys are UniProt IDs and values are lists of mapped IDs from the destination database.
    """

    POLLING_INTERVAL = sleep_time
    API_URL = "https://rest.uniprot.org"
    
    retries = Retry(total=max_retries, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))


    def check_response(response):
        try:
            response.raise_for_status()
        except requests.HTTPError:
            print(response.json())
            raise


    def submit_id_mapping(from_db, to_db, ids):
        request = requests.post(
            f"{API_URL}/idmapping/run",
            data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
        )
        check_response(request)
        return request.json()["jobId"]


    def get_next_link(headers):
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)


    def check_id_mapping_results_ready(job_id):
        while True:
            request = session.get(f"{API_URL}/idmapping/status/{job_id}")
            check_response(request)
            j = request.json()
            if "jobStatus" in j:
                if j["jobStatus"] in ("NEW", "RUNNING"):
                    print(f"Retrying in {POLLING_INTERVAL}s")
                    time.sleep(POLLING_INTERVAL)
                else:
                    raise Exception(j["jobStatus"])
            else:
                return bool(j["results"] or j["failedIds"])


    def get_batch(batch_response, file_format, compressed):
        batch_url = get_next_link(batch_response.headers)
        while batch_url:
            batch_response = session.get(batch_url)
            batch_response.raise_for_status()
            yield decode_results(batch_response, file_format, compressed)
            batch_url = get_next_link(batch_response.headers)


    def combine_batches(all_results, batch_results, file_format):
        if file_format == "json":
            for key in ("results", "failedIds"):
                if key in batch_results and batch_results[key]:
                    all_results[key] += batch_results[key]
        elif file_format == "tsv":
            return all_results + batch_results[1:]
        else:
            return all_results + batch_results
        return all_results


    def get_id_mapping_results_link(job_id):
        url = f"{API_URL}/idmapping/details/{job_id}"
        request = session.get(url)
        check_response(request)
        return request.json()["redirectURL"]


    def decode_results(response, file_format, compressed):
        if compressed:
            decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
            if file_format == "json":
                j = json.loads(decompressed.decode("utf-8"))
                return j
            elif file_format == "tsv":
                return [line for line in decompressed.decode("utf-8").split("\n") if line]
            elif file_format == "xlsx":
                return [decompressed]
            elif file_format == "xml":
                return [decompressed.decode("utf-8")]
            else:
                return decompressed.decode("utf-8")
        elif file_format == "json":
            return response.json()
        elif file_format == "tsv":
            return [line for line in response.text.split("\n") if line]
        elif file_format == "xlsx":
            return [response.content]
        elif file_format == "xml":
            return [response.text]
        return response.text


    def get_xml_namespace(element):
        m = re.match(r"\{(.*)\}", element.tag)
        return m.groups()[0] if m else ""


    def merge_xml_results(xml_results):
        merged_root = ElementTree.fromstring(xml_results[0])
        for result in xml_results[1:]:
            root = ElementTree.fromstring(result)
            for child in root.findall("{http://uniprot.org/uniprot}entry"):
                merged_root.insert(-1, child)
        ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
        return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


    def print_progress_batches(batch_index, size, total):
        n_fetched = min((batch_index + 1) * size, total)
        print(f"Fetched: {n_fetched} / {total}")


    def get_id_mapping_results_search(url):
        parsed = urlparse(url)
        query = parse_qs(parsed.query)
        file_format = query["format"][0] if "format" in query else "json"
        if "size" in query:
            size = int(query["size"][0])
        else:
            size = 500
            query["size"] = size
        compressed = (
            query["compressed"][0].lower() == "true" if "compressed" in query else False
        )
        parsed = parsed._replace(query=urlencode(query, doseq=True))
        url = parsed.geturl()
        request = session.get(url)
        check_response(request)
        results = decode_results(request, file_format, compressed)
        total = int(request.headers["x-total-results"])
        print_progress_batches(0, size, total)
        for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
            results = combine_batches(results, batch, file_format)
            print_progress_batches(i, size, total)
        if file_format == "xml":
            return merge_xml_results(results)
        return results


    def get_id_mapping_results_stream(url):
        if "/stream/" not in url:
            url = url.replace("/results/", "/results/stream/")
        request = session.get(url)
        check_response(request)
        parsed = urlparse(url)
        query = parse_qs(parsed.query)
        file_format = query["format"][0] if "format" in query else "json"
        compressed = (
            query["compressed"][0].lower() == "true" if "compressed" in query else False
        )
        return decode_results(request, file_format, compressed)
                                                        
    # Submit mapping from UniProt to PDB
    job_id = submit_id_mapping(from_db=source, to_db=dest, ids=ids)
    print(f"Job ID: {job_id}")

    # Check if results are ready
    if check_id_mapping_results_ready(job_id):
        # Get the link to fetch results
        link = get_id_mapping_results_link(job_id)
        print(f"Results link: {link}")
        
        # Retrieve results
        results = get_id_mapping_results_search(link)
        print(f"Number of results: {len(results)}")

    # Parse the results into a dictionary
    dict_proteome = {}

    ids_not_mapped = list(set(results["failedIds"]))
    ids_mapped = results["results"]

    for result in ids_mapped:
        uniprot_id = result['from']
        mapped_id = result['to']
        if uniprot_id in dict_proteome:
            dict_proteome[uniprot_id].append(mapped_id)
        else:
            dict_proteome[uniprot_id] = [mapped_id]

    print(f"Number of IDs mapped: {len(ids_mapped)}")
    print(f"Number of IDs not mapped: {len(ids_not_mapped)}")

    return dict_proteome, ids_not_mapped

def download_DEG(base_path):

    """
    Downloads DEG Bacteria database. 
    The fasta file is located in "databases" folder.
    https://tubic.org/deg/public/index.php
    DOI: 10.1093/nar/gkaa917

    :param base_path =  Path to the directory where 'databases' folder is located.
    :return: download_status = True if download was successful, False otherwise.
    """
    databases_path = os.path.join(base_path, 'databases')
    deg_path = os.path.join(databases_path, 'DEG10.aa.gz')
   
    url = 'http://tubic.org/deg/public/download/DEG10.aa.gz'
    
    download_status = False

    print('Downloading DEG database.')
    download_with_wget(url, deg_path)

    if files.file_check(deg_path):
        print('Finished.')
    
        try: 
            with gzip.open(deg_path, 'rb') as f_in:
                deg_fasta = os.path.splitext(deg_path)[0] + ".faa"
                with open(deg_fasta, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    print(f'File {deg_fasta} downloaded and decompressed successfully.')
            os.remove(deg_path)
            download_status = True
        except Exception as e:
            print(f"An error occurred while decompressing the DEG database: {e}")
    else:
        print(f'Failed to download DEG database from {url}.')
        print('Please check if DEG database is available at the source URL.'
              ' You can also try to download it manually and place it in the databases folder.')
    return download_status

def extract_microbiome_species_ids(base_path):
    """
    Extracts unique species IDs from the Human Gut Microbiome Species Catalogue metadata file.
    :param base_path =  Path to the directory where 'databases' folder is located.
    :return: species_ids = List of unique species IDs.
    """

    species_path = os.path.join(base_path, 'databases', 'species_catalogue')
    meta_path = os.path.join(species_path, 'genomes-all_metadata.tsv')

    if  os.path.exists(meta_path):
        df = pd.read_csv(meta_path, sep='\t')
        species_ids = df['Species_rep'].unique().tolist()
        return species_ids
    else:
        raise FileNotFoundError(f"Metadata file not found: {meta_path}")


def download_microbiome_species_catalogue(base_path):

    """
    Downloads the Human Gut Microbiome Species Catalogue from the EBI Metagenomics (MGnify) database, 
    Each genome has its own folder with a .faa file inside 'genome'.
    Files are saved into 'databases/species_catalogue'.
    https://www.ebi.ac.uk/metagenomics/genome-catalogues/human-gut-v2-0-2
    DOI: 10.1093/nar/gkz1035

    :param base_path =  Path to the directory where 'databases' folder is located.

    """

    databases_path = os.path.join(base_path, 'databases')  
    species_path = os.path.join(databases_path, 'species_catalogue')
    os.makedirs(species_path, exist_ok=True)

    # Download README file
    info_path = os.path.join(species_path, 'README_uhgp_v2.0.2.txt')
    info_url = 'ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/README_v2.0.2.txt'
    download_with_wget(info_url, info_path)

    # Download metadata file
    meta_path = os.path.join(species_path, 'genomes-all_metadata.tsv')
    meta_url = 'ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/genomes-all_metadata.tsv'
    download_with_wget(meta_url, meta_path)

    # Download species catalogue

    check_file = os.path.join(species_path, "download_check.txt")

    # Check if download was ssuccessfully completed before
    
    success_file = False
    if os.path.exists(check_file):
        with open(check_file, "a+") as f:
            f.seek(0)
            content = f.read()
            if "successfully" in content:
                success_file = True
            else:
                os.remove(check_file)

    success_check = check_microbiome_species_catalogue_download(base_path)

    if success_file and success_check:
        print("Species catalogue download was already completed successfully.")

    else:

        # Parse the main directory to get species folders        
        base_url = "http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/species_catalogue/"

        r = requests.get(base_url)
        r.raise_for_status()
        soup = BeautifulSoup(r.text, "html.parser")

        # Obtain folder links 
        species_links = [a["href"] for a in soup.find_all("a") if a["href"].startswith("MGYG")]
        expected_files = []
        genome_list = []

        for sp in species_links:
            sp_url = base_url + sp
            r2 = requests.get(sp_url)
            r2.raise_for_status()
            soup2 = BeautifulSoup(r2.text, "html.parser")

            genome_links = [a["href"] for a in soup2.find_all("a") if a["href"].startswith("MGYG")]

            for genome in genome_links:
                genome_list.append(genome[:-1])

                faa_url = f"{sp_url}{genome}genome/{genome[:-1]}.faa"
                faa_ftp_url = faa_url.replace("http://", "ftp://")

                genome_path = os.path.join(species_path, genome[:-1])
                os.makedirs(genome_path, exist_ok=True)

                faa_name = os.path.basename(faa_url)
                faa_out = os.path.join(genome_path, faa_name)
                expected_files.append(faa_name)

                if not files.file_check(faa_out):
                    print(f"Downloading {faa_name}...")
                    download_with_wget(faa_ftp_url, faa_out)
                else:
                    print(f"{faa_name} already exists.")
        
        # Check for missing files
        downloaded_files = [os.path.basename(f) for f in glob.glob(os.path.join(species_path, "**", "*.faa"), recursive=True)]
        missing_files = sorted(set(expected_files) - set(downloaded_files))
        
        with open(check_file, "w") as f:
            f.write(f"Expected files: {len(expected_files)}\n")
            f.write(f"Downloaded files: {len(downloaded_files)}\n")
            f.write(f"Missing files: {len(missing_files)}\n\n")
            if missing_files:
                f.write("List of missing files:\n")
                for mf in missing_files:
                    f.write(mf + "\n")
            else:
                f.write("All files downloaded successfully!\n")

def check_microbiome_species_catalogue_download(base_path):
    """
    Checks if all .faa files from the species catalogue have been downloaded.
    The concatenated fasta file is located in "databases/species_catalogue" folder.

    :param base_path =  Path to the directory where 'databases' folder is located.

    """

    species_path = os.path.join(base_path, 'databases', 'species_catalogue')

    species_ids = extract_microbiome_species_ids(base_path)

    for id in species_ids:
        faa_file = os.path.join(species_path, id, f"{id}.faa")
        if not files.file_check(faa_file):
            print(f"Missing file: {faa_file}")
            return False

def concat_microbiome_species_catalogue(base_path):
    """
    Concatenates all .faa files from the species catalogue into a single FASTA file.
    The concatenated fasta file is located in "databases/species_catalogue" folder.

    :param base_path =  Path to the directory where 'databases' folder is located.

    """

    species_path = os.path.join(base_path, 'databases', 'species_catalogue')
    output_faa = os.path.join(species_path, 'species_catalogue.faa')

    if not os.path.exists(output_faa):
        print('Concatenating species catalogue .faa files...')
        simple_concat_faa(species_path, output_faa)
        print(f'File {output_faa} created successfully.')
    else:
        print(f'{output_faa} already exists.')

    
def download_microbiome_protein_catalogue(base_path):

    """
    Downloads Human-gut protein catalogue from EBI Metagenomics (MGnify) clustered at 90% aa identity. 
    The fasta file is located in "databases" folder.
    https://www.ebi.ac.uk/metagenomics/genome-catalogues/human-gut-v2-0-2
    DOI: 10.1093/nar/gkz1035

    :param base_path =  Path to the directory where 'databases' folder is located.

    """

    databases_path = os.path.join(base_path, 'databases')
    extracted_dir = os.path.join(databases_path, 'uhgp-90')
    faa_file_path = os.path.join(extracted_dir, 'uhgp-90.faa')

    info_path = os.path.join(databases_path, 'README_uhgp_v2.0.2.txt')
    info_url = 'ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/README_v2.0.2.txt'
    download_with_wget(info_url, info_path)

    if not files.file_check(os.path.join(databases_path, 'uhgp-90.faa')):
        uhgp_path = os.path.join(databases_path, 'uhgp-90.tar.gz')
        uhgp_url = 'ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/protein_catalogue/uhgp-90.tar.gz'
        
        print('Downloading UHGP database.')
        download_with_wget(uhgp_url, uhgp_path)
        print('Finished.')

        with tarfile.open(uhgp_path, 'r:gz') as tar:
            tar.extractall(path=databases_path)

        if os.path.exists(faa_file_path):
            shutil.move(faa_file_path, databases_path)
        
        shutil.rmtree(extracted_dir)

        print(f'File {uhgp_path} downloaded and decompressed successfully')
    else:
        print(f'{faa_file_path} already exists.')

def download_human_PDB (base_path, cpus=multiprocessing.cpu_count()):
    """
    Downloads Human PDB structures.
    The PDB files are located in "databases/human_structures/PDB_files" folder.
    The list of PDB IDs is saved in "databases/human_structures/PDB_files/Human_PDB_ids.txt".
    The IDs that could not be downloaded are saved in "databases/human_structures/PDB_files/Failed_Human_PDB_ids.txt".

    :param base_path =  Path to the directory where 'fasttarget' folder is located.
    :param cpus = Number of CPUs to use for downloading.

    :return: all_pdb_ids = List with all PDB IDs.
    :return: failed_pdb = List with PDB IDs that failed to download.

    """

    human_structures_pdb_path = os.path.join(base_path, 'databases', 'human_structures', 'PDB_files')
    os.makedirs(human_structures_pdb_path, exist_ok=True)
    os.makedirs(os.path.join(human_structures_pdb_path, 'DB_foldseek'), exist_ok=True)

    human_pdb_ids_path = os.path.join(human_structures_pdb_path, 'Human_PDB_ids.txt')
    human_ids_without_pdb_path = os.path.join(human_structures_pdb_path, 'Human_prots_without_PDB.txt')

    def download_pdb(pdb_id):
        """
        Download PDB structure. Returns PDB ID if failed.
        
        :param pdb_id: ID of PDB structure.
        
        :return: PDB ID if failed.
        """

        file_path_pdb = os.path.join(human_structures_pdb_path, f'PDB_{pdb_id}.pdb')
        file_path_cif = os.path.join(human_structures_pdb_path, f'PDB_{pdb_id}.cif')

        if not files.file_check(file_path_pdb) and not files.file_check(file_path_cif):
            pdb_res = structures.get_structure_PDB(human_structures_pdb_path, pdb_id)
            if not pdb_res:
                cif_res = structures.get_structure_CIF(human_structures_pdb_path, pdb_id)
                if not cif_res:
                    print(f'Failed to download PDB {pdb_id}.')
                    return pdb_id
                    
        else:
            print(f'PDB {pdb_id} already exists.')
        return None

    if not files.file_check(human_pdb_ids_path):

        human_uniprot_ids = structures.uniprot_proteome_ids('UP000005640')
        print('Number of Human Uniprot IDs:', len(human_uniprot_ids))

        print('----Mapping Human Uniprot IDs to PDB----')
        result_dict, ids_not_mapped= batch_uniprot_mapping('UniProtKB_AC-ID', 'PDB', human_uniprot_ids)

        all_pdb_ids = []
        for values in result_dict.values():
            all_pdb_ids.extend(values)

        all_pdb_ids = list(set(all_pdb_ids))

        # Save the list of PDB IDs to a text file
        files.list_to_file(human_pdb_ids_path, all_pdb_ids)      
        files.list_to_file(human_ids_without_pdb_path, ids_not_mapped)    
    else:
        all_pdb_ids = files.file_to_list(human_pdb_ids_path)
        all_pdb_ids = list(set(all_pdb_ids))

    print('Number of Human PDB structures:', len(all_pdb_ids))

    # Parallel download of PDB structures
    failed_pdb = []
    print('----Downloading Human PDB structures----')
    with ThreadPoolExecutor(max_workers=cpus) as executor:
        futures_pdb = {executor.submit(download_pdb, pdb_id): pdb_id for pdb_id in all_pdb_ids}
        for future in tqdm.tqdm(as_completed(futures_pdb), total=len(all_pdb_ids), desc="Downloading from PDB", unit="pdb"):
            result = future.result()
            if result:
                failed_pdb.append(result)
    print('Number of PDB structures that could not be downloaded:', len(failed_pdb))
    
    # Last try for failed PDB IDs
    print('----Trying to download failed PDB IDs----')
    final_failed_pdb = []
    if failed_pdb:
        for fail_id in tqdm.tqdm(failed_pdb, desc="Downloading failed PDB IDs", unit="pdb"):
            res = download_pdb(fail_id)
            if res:
                final_failed_pdb.append(res)
    
    print('----Download finished!----')

    fail_file_pdb = os.path.join(human_structures_pdb_path, 'Failed_Human_PDB_ids.txt')
    files.list_to_file(fail_file_pdb, final_failed_pdb)
    print('Number of PDB structures that could not be downloaded:', len(final_failed_pdb))
    print(f"The IDs that could not be downloaded are saved in {fail_file_pdb}")  

    return all_pdb_ids, final_failed_pdb


def download_human_AF (base_path):
    """
    Downloads Human AlphaFold structures.
    The PDB files are located in "databases/human_structures/AlphaFold_files" folder.
    The list of AlphaFold IDs is saved in "databases/human_structures/AlphaFold_files/Human_AF_ids.txt".

    :param base_path =  Path to the directory where 'fasttarget' folder is located.

    :return: AF_models = List with AlphaFold IDs.

    """
    human_structures_AF_path = os.path.join(base_path, 'databases', 'human_structures', 'AlphaFold_files')
    os.makedirs(human_structures_AF_path, exist_ok=True)
    os.makedirs(os.path.join(human_structures_AF_path, 'DB_foldseek'), exist_ok=True)

    af_ids_file = os.path.join(human_structures_AF_path, 'Human_AF_ids.txt')

    if not files.file_check(af_ids_file):

        human_AF_path = os.path.join(human_structures_AF_path, 'UP000005640_9606_HUMAN_v6.tar')

        if not os.path.exists(human_AF_path):       
            url_human_AF = 'ftp://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v6.tar'
            try:
                print('Downloading AlphaFold human database with progress bar.')
                download_with_wget(url_human_AF, human_AF_path)
                print('AlphaFold human database downloaded.')
            except Exception as e:
                print(f'Error downloading AlphaFold files: {e}')
                # remove incomplete file if exists
                if os.path.exists(human_AF_path):
                    os.remove(human_AF_path)

            try:
                with tarfile.open(human_AF_path, 'r') as tar:
                    tar.extractall(path=human_structures_AF_path)
            except Exception as e:
                print(f'Error decompressing AlphaFold files: {e}')

        AF_models = []

        for filename in os.listdir(human_structures_AF_path):
            file_path = os.path.join(human_structures_AF_path, filename)

            # Remove cif files
            if filename.endswith('.cif.gz'):
                os.remove(file_path)
            
            # Uncompress gz files
            elif filename.endswith('.pdb.gz'):
                with gzip.open(file_path, 'rb') as f_in:
                    pdb_filename = file_path[:-3]
                    with open(pdb_filename, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                    f_name = os.path.basename(pdb_filename)
                    AF_models.append(f_name)
                os.remove(file_path)

        files.list_to_file(af_ids_file, AF_models)

    else:
        AF_models = files.file_to_list(af_ids_file)
        print('AlphaFold human database already exists.')

    print('Number of AlphaFold structures:', len(AF_models))
    
    return AF_models
    

def download_human(base_path, cpus=multiprocessing.cpu_count()):
    """
    Downloads Human proteome form Uniprot (UP000005640). 
    Obtains fasta sequences and structures from AlphaFold and PDB.
    The fasta file and annotations are located in "databases" folder.
    The structures are located in "databases/human_structures" folder.

    :param base_path =  Path to the directory where 'fasttarget' repository is located.

    :return: human_ann_dict = Dictionary with Uniprot annotations for human proteome.
    :return: human_failed_pdb = List with PDB IDs that failed to download.
    :return: human_failed_AF = List with AlphaFold IDs that failed to download.

    """
    #Download fasta sequences
    databases_path = os.path.join(base_path, 'databases')
    humanprot_path = os.path.join(databases_path, 'human_uniprot_UP000005640.faa')

    if not os.path.exists(humanprot_path):
    
        url = 'ftp://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28proteome%3AUP000005640%29%29'
        
        print('Downloading human proteome: fasta sequences.')
        download_with_wget(url, humanprot_path)
        print('Finished.')

        print(f'File {humanprot_path} downloaded successfully')
    else:
        print(f'{humanprot_path} already exists.')

    #Download structures
    human_structures_path = os.path.join(databases_path, 'human_structures')

    try:
        print('Downloading human PDB structures.')
        human_pdb_ids, human_failed_pdb = retry_on_timeout(download_human_PDB, base_path, cpus)
        print('PDB structures downloaded.')
    except Exception as e:
        print(f'Error downloading human PDB structures: {e}')

    try:
        print('Downloading AlphaFold structures.')
        human_alphafold_ids = retry_on_timeout(download_human_AF, base_path)
        print('AlphaFold structures downloaded.')
    except Exception as e:
        print(f'Error downloading AlphaFold structures: {e}')

    #Check if all structures were downloaded

    human_downloaded_pdb_files = glob.glob(os.path.join(human_structures_path, 'PDB_files', '*.pdb'))
    human_downloaded_cif_files = glob.glob(os.path.join(human_structures_path, 'PDB_files', '*.cif'))

    human_downloaded_AF_files = glob.glob(os.path.join(human_structures_path, 'AlphaFold_files', '*.pdb'))

    check = False
    PDB_complete = len(human_pdb_ids) == len(human_downloaded_pdb_files) + len(human_downloaded_cif_files) + len(human_failed_pdb)
    AF_complete = len(human_alphafold_ids) == len(human_downloaded_AF_files)

    if PDB_complete and AF_complete:
        print(f'All structures downloaded correctly in {human_structures_path}.')
        check = True
        # Save the list of downloaded PDB and AlphaFold files to a text file
        downloaded_files_path = os.path.join(human_structures_path, 'human_structures.txt')
        with open(downloaded_files_path, 'w') as f:
            for pdb_file in human_downloaded_pdb_files:
                f.write(f'{pdb_file}\n')
            for af_file in human_downloaded_AF_files:
                f.write(f'{af_file}\n')
        print(f'List of downloaded files saved to {downloaded_files_path}')
    else:
        print(f'PDB complete: {PDB_complete}, AlphaFold complete: {AF_complete}')
        print(f'Number of downloaded PDB files: {len(human_downloaded_pdb_files) + len(human_downloaded_cif_files)} / total: {len(human_pdb_ids)}')
        print(f'Number of downloaded AlphaFold files: {len(human_downloaded_AF_files)} / total: {len(human_alphafold_ids)}')

        raise Exception(f'Error downloading structures in {human_structures_path}. Check the failed files.')

    return check

def index_db_blast_human (base_path):

    """
    Makes a BLAST db for human proteome. 

    :param base_path =  Base path of fasttarget folder.

    """

    databases_path = os.path.join(base_path, 'databases')

    #Human
    humanprot_path = os.path.join(databases_path, 'human_uniprot_UP000005640.faa')
    humanprot_index_path = os.path.join(databases_path, 'HUMAN_DB')

    try:
        programs.run_makeblastdb(
        input= humanprot_path,
        output= humanprot_index_path,
        title= 'HUMAN_DB',
        dbtype= 'prot'
        )
    except Exception as e:
        print('Error indexing Human database:', e)

def index_db_foldseek_human_structures (base_path):

    """
    Makes a FOLDSEEK db for human proteome PDB and AlphaFold structures. 

    :param base_path =  Base path of fasttarget folder.

    """

    human_structures_path = os.path.join(base_path, 'databases', 'human_structures')
    human_PDB_path = os.path.join(human_structures_path, 'PDB_files')
    human_AF_path = os.path.join(human_structures_path, 'AlphaFold_files')

    #Index PDB structures
    try:
        programs.run_foldseek_create_index_db(human_PDB_path, 'DB_human_PDB')
        print('Foldseek PDB database created successfully.')
    except Exception as e:
        print('Error indexing human PDB structures for Foldseek:', e)

    #Index AlphaFold structures
    try:
        programs.run_foldseek_create_index_db(human_AF_path, 'DB_human_AF')
        print('Foldseek AF database created successfully.')
    except Exception as e:
        print('Error indexing human AlphaFold structures for Foldseek.', e)


def index_db_blast_microbiome_protein_catalogue (base_path):

    """
    Makes a BLAST db for gut microbiome protein catalogue. 

    :param base_path =  Base path of fasttarget folder.

    """

    databases_path = os.path.join(base_path, 'databases')

    #Gut Microbiome
    microbiome_path = os.path.join(databases_path, 'uhgp-90.faa')
    microbiome_index_path = os.path.join(databases_path, 'MICROBIOME_DB')

    try:
        programs.run_makeblastdb(
        input= microbiome_path,
        output= microbiome_index_path,
        title= 'MICROBIOME_DB',
        dbtype= 'prot'
        )
    except Exception as e:
        print('Error indexing Microbiome database:', e)


def index_db_blast_microbiome_species_catalogue (base_path, specific_file=None):

    """
    Makes a Diamond BLAST db for gut microbiome species catalogue. 

    :param base_path =  Base path of fasttarget folder.
    :param specific_file =  Specific .faa file to index. If None, all .faa files in the
    species_catalogue folder will be indexed.

    """

    species_path = os.path.join(base_path, 'databases', 'species_catalogue')
    faa_files = glob.glob(os.path.join(species_path, "**", "*.faa"), recursive=True)

    if specific_file is None:

        for faa_file in faa_files:
            if not files.file_check(faa_file):
                raise FileNotFoundError(f"File {faa_file} does not exist.")
            else:
                try:
                    programs.run_makediamonddb(
                    input= faa_file,
                    output= faa_file.replace(".faa", "_DB")
                    )
                except Exception as e:
                    print(f'Error indexing {faa_file} for Diamond DB:', e)
    else:
        if not files.file_check(specific_file):
            raise FileNotFoundError(f"File {specific_file} does not exist.")
        else:
            try:
                programs.run_makediamonddb(
                input= specific_file,
                output= specific_file.replace(".faa", "_DB")
                )
            except Exception as e:
                print(f'Error indexing {specific_file} for Diamond DB:', e)

def index_db_blast_deg (base_path):

    """
    Makes a BLAST db for DEG. 

    :param base_path =  Base path of fasttarget folder.

    """

    databases_path = os.path.join(base_path, 'databases')

    #DEG
    deg_path = os.path.join(databases_path, 'DEG10.aa.faa')
    deg_index_path = os.path.join(databases_path, 'DEG_DB')

    try:
        programs.run_makeblastdb(
        input= deg_path,
        output= deg_index_path,
        title= 'DEG_DB',
        dbtype= 'prot'
        )
    except Exception as e:
        print('Error indexing DEG database:', e)



def download_and_index_human(databases_path, force=False):
    """
    Downloads and indexes Human proteome.

    :param databases_path =  Path to the 'databases' folder.
    :param force = Force re-download even if files exist.

    """
    
    # Download and index Human proteome
    print('----- 1. Downloading and indexing human proteome -----')
    
    humanprot_path = os.path.join(databases_path, 'human_uniprot_UP000005640.faa')
    human_structures_path = os.path.join(databases_path, 'human_structures')
    downloaded_files_path = os.path.join(human_structures_path, 'human_structures.txt')

    if force or not files.file_check(humanprot_path) or not files.file_check(downloaded_files_path):
        try:
            check = download_human(base_path)
            print('Human proteome downloaded')

            if check:
                print('Download complete.')
            else:
                print('Error downloading human proteome. Check the failed files.')
            # Index human proteome
            index_db_blast_human(base_path)
            print('Human blast db created')
            index_db_foldseek_human_structures(base_path)
            print('Human foldseek db created')
            print('Human proteome downloaded and indexed')
        except Exception as e:
            print(f'Error downloading human proteome: {e}')  
    else:
        print('Human proteome already exists.')
        if force or not files.file_check(os.path.join(databases_path, 'HUMAN_DB.phr')):
            print('Indexing human proteome')
            index_db_blast_human(base_path)
            print('Human proteome indexed')
        else:
            print('Human proteome already indexed')

        human_PDB_db = os.path.join(human_structures_path, 'PDB_files', 'DB_foldseek', 'DB_human_PDB')
        human_AF_db = os.path.join(human_structures_path, 'AlphaFold_files', 'DB_foldseek', 'DB_human_AF')

        if force or not os.path.exists(human_PDB_db) or not os.path.exists(human_AF_db):
            print('Indexing human proteome for Foldseek')
            index_db_foldseek_human_structures(base_path)
            print('Human foldseek db created')
        else:
            print('Human foldseek db already created')

    print('----- 1. Finished -----')

def download_and_index_microbiome(databases_path):
    """
    Downloads and indexes Microbiome database.

    :param databases_path =  Path to the 'databases' folder.

    """

    # Download and index Microbiome database
    print('----- 2. Downloading and indexing microbiome database -----')
    species_path = os.path.join(databases_path, 'species_catalogue')
    check_file = os.path.join(species_path, 'download_check.txt')

    if not files.file_check(check_file):
        download_microbiome_species_catalogue (base_path)
        index_db_blast_microbiome_species_catalogue (base_path)
        print('Microbiome database downloaded and indexed')
    else:
        with open(check_file) as f:
            content = f.read()
            if "List of missing files:" in content:
                print('Some files are missing in the species catalogue. Please check download_check.txt file.')                     
                download_microbiome_species_catalogue (base_path)
                index_db_blast_microbiome_species_catalogue (base_path)
            else:
                print('Microbiome sequences already exists.')

        # Check if all .faa files have been indexed
        check = check_microbiome_species_catalogue_downloaded(databases_path)
        if len(check) > 0:
            for faa_path in check:
                index_db_blast_microbiome_species_catalogue(base_path, specific_file=faa_path)

    print('----- 2. Finished -----')

def check_microbiome_species_catalogue_downloaded(databases_path):
    """
    Check if all species in the microbiome species catalogue have been indexed.
    
    :param databases_path =  Path to the 'databases' folder.
    """

    species_path = os.path.join(databases_path, 'species_catalogue')
    required_extensions = {".dmnd"}

    check = []
    
    for root, dirs, allfiles in os.walk(species_path):
        if root == species_path:
            continue

        files_set = set(os.path.splitext(f)[1] for f in allfiles)

        if required_extensions.issubset(files_set):
            continue

        faa_files = [f for f in allfiles if f.endswith(".faa")]

        if not faa_files:
            print(f"⚠️ No .faa or index files found in {root}")
        else:
            faa_path = os.path.join(root, faa_files[0])
            check.append(faa_path)
            print(f"❌ No diamond blast index files in {root} for {faa_path}")
    
    return check


def download_and_index_deg(databases_path):
    """
    Downloads and indexes DEG database.

    :param databases_path =  Path to the 'databases' folder.

    """

    # Download and index DEG database
    print('----- 3. Downloading and indexing DEG database -----')
    if not files.file_check(os.path.join(databases_path, 'DEG10.aa.faa')):
        download_status = download_DEG (base_path)
        if download_status:
            index_db_blast_deg (base_path)
            print('DEG database downloaded and indexed')
        else:
            print('Failed to download DEG database.')
    else:
        print('DEG database already exists.')
        if not files.file_check(os.path.join(databases_path, 'DEG_DB.phr')):
            print('Indexing DEG database')
            index_db_blast_deg (base_path)
            print('DEG database indexed')
        else:
            print('DEG database already indexed')
    print('----- 3. Finished -----')



def main_download(base_path, selected_databases, force=False, parallel=4):
    """Main download function with enhanced options"""
    print(f"📥 Downloading databases: {', '.join(selected_databases)}")
    if force:
        print("🔄 Force mode: Will re-download existing databases")
    
    databases_path = os.path.join(base_path, 'databases')
    if not os.path.exists(databases_path):
        os.makedirs(databases_path)
    
    success_count = 0
    for db_name in selected_databases:
        try:
            print(f"\n{'='*20} {db_name.upper()} {'='*20}")
            
            if db_name == 'human':
                download_and_index_human(databases_path)
                success_count += 1
                    
            elif db_name == 'microbiome':
                download_and_index_microbiome(databases_path)
                success_count += 1
                    
            elif db_name == 'deg':
                download_and_index_deg(databases_path)
                success_count += 1
                    
        except Exception as e:
            print(f"❌ Failed to download {db_name}: {e}")
    
    print(f"\n🎯 Summary: {success_count}/{len(selected_databases)} databases completed successfully")


if __name__ == '__main__':
    base_path = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(
        description="FastTarget Database Management",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Examples:
        python databases.py --download all          # Download all databases
        python databases.py --download human        # Download only human proteome
        python databases.py --update microbiome     # Update microbiome database
        python databases.py --verify                # Check all databases
        python databases.py --status                # Show database status
        python databases.py --clean                 # Clean incomplete downloads
                """
        )
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--download', choices=['human', 'microbiome', 'deg', 'all'], 
                      help="Download specified database(s)")
    group.add_argument('--update', choices=['human', 'microbiome', 'deg', 'all'],
                      help="Update specified database(s)")
    group.add_argument('--verify', action='store_true',
                      help="Verify all databases are complete and indexed")
    parser.add_argument('--force', action='store_true',
                       help="Force re-download even if database exists")
    parser.add_argument('--parallel', type=int, default=4,
                       help="Number of parallel downloads (default: 4)")

    args = parser.parse_args()

    # Determine databases to process
    if hasattr(args, 'download') and args.download:
        databases_to_process = ['human', 'microbiome', 'deg'] if args.download == 'all' else [args.download]
        main_download(base_path, databases_to_process, force=args.force, parallel=args.parallel)
    elif hasattr(args, 'update') and args.update:
        databases_to_process = ['human', 'microbiome', 'deg'] if args.update == 'all' else [args.update]
        main_update(base_path, databases_to_process, parallel=args.parallel)
    elif args.verify:
        verify_databases(base_path)