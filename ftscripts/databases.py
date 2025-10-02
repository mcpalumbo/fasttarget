import urllib.request
import gzip
import tarfile
import shutil
import os
import time
from ftscripts import programs,files, structures
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

def download_with_progress(url, filepath):
    """
    Download a file from a URL with a progress indicator.

    :param url: The URL to download from.
    :param filepath: The local path where the file should be saved.
    """
    def reporthook(count, block_size, total_size):
        # If total_size is 0 or unknown, show the amount downloaded in MB
        downloaded_mb = count * block_size / (1024 * 1024)  # Convert to MB
        print(f'\rDownloaded {downloaded_mb:.2f} MB', end='', flush=True)

    try:
        # Make the download request and use the reporthook to show progress
        urllib.request.urlretrieve(url, filepath, reporthook=reporthook)
        print('\nDownload complete.')
    except Exception as e:
        print(f'Error downloading {url}: {e}')

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

    """
    databases_path = os.path.join(base_path, 'databases')
    deg_path = os.path.join(databases_path, 'DEG10.aa.gz')
   
    url = 'http://tubic.org/deg/public/download/DEG10.aa.gz'
    
    print('Downloading DEG database.')
    download_with_wget(url, deg_path)
    print('Finished.')

    with gzip.open(deg_path, 'rb') as f_in:
        deg_fasta = os.path.splitext(deg_path)[0] + ".faa"
        with open(deg_fasta, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            print(f'File {deg_fasta} downloaded and decompressed successfully.')
    os.remove(deg_path)        

def download_microbiome(base_path):

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
    info_url = 'https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/README_v2.0.2.txt'
    download_with_wget(info_url, info_path)

    if not files.file_check(os.path.join(databases_path, 'uhgp-90.faa')):
        uhgp_path = os.path.join(databases_path, 'uhgp-90.tar.gz')
        uhgp_url = 'https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/protein_catalogue/uhgp-90.tar.gz'
        
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
                    return pdb_id
                    print(f'PDB {pdb_id} could not be downloaded in PDB or CIF format.')
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

        human_AF_path = os.path.join(human_structures_AF_path, 'UP000005640_9606_HUMAN_v4.tar')

        if not os.path.exists(human_AF_path):       
            url_human_AF = 'https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar'
            download_with_wget(url_human_AF, human_AF_path)

            print('AlphaFold human database downloaded.')
            print('Decompressing AlphaFold files.')
        else:
            with tarfile.open(human_AF_path, 'r') as tar:
                #tar.extractall(path=human_structures_AF_path)
                print('AlphaFold human database already exists.')

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
    
        url = 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28proteome%3AUP000005640%29%29'
        
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
        print(f'Error downloading structures in {human_structures_path}. Check the failed files.')

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

    programs.run_makeblastdb(
    input= humanprot_path,
    output= humanprot_index_path,
    title= 'HUMAN_DB',
    dbtype= 'prot'
    )

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


def index_db_blast_microbiome (base_path):

    """
    Makes a BLAST db for gut microbiome. 

    :param base_path =  Base path of fasttarget folder.

    """

    databases_path = os.path.join(base_path, 'databases')

    #Gut Microbiome
    microbiome_path = os.path.join(databases_path, 'uhgp-90.faa')
    microbiome_index_path = os.path.join(databases_path, 'MICROBIOME_DB')

    programs.run_makeblastdb(
    input= microbiome_path,
    output= microbiome_index_path,
    title= 'MICROBIOME_DB',
    dbtype= 'prot'
    )

def index_db_blast_deg (base_path):

    """
    Makes a BLAST db for DEG. 

    :param base_path =  Base path of fasttarget folder.

    """

    databases_path = os.path.join(base_path, 'databases')

    #DEG
    deg_path = os.path.join(databases_path, 'DEG10.aa.faa')
    deg_index_path = os.path.join(databases_path, 'DEG_DB')

    programs.run_makeblastdb(
    input= deg_path,
    output= deg_index_path,
    title= 'DEG_DB',
    dbtype= 'prot'
    )

def main(base_path):
    """
    Downloads and indexes Human proteome, Microbiome database and DEG database.

    :param base_path =  Path to the repository folder.

    """

    databases_path = os.path.join(base_path, 'databases')
    if not os.path.exists(databases_path):
        os.makedirs(databases_path)
    
    # Download and index Human proteome
    print('----- 1. Downloading and indexing human proteome -----')
    humanprot_path = os.path.join(databases_path, 'human_uniprot_UP000005640.faa')
    human_structures_path = os.path.join(databases_path, 'human_structures')
    downloaded_files_path = os.path.join(human_structures_path, 'human_structures.txt')

    if not files.file_check(os.path.join(databases_path, 'human_uniprot_UP000005640.faa')) or not files.file_check(downloaded_files_path):
        try:
            check = download_human (base_path)
            print('Human proteome downloaded')

            if check:
                print('Download complete.')
            else:
                print('Error downloading human proteome. Check the failed files.')
            # Index human proteome
            index_db_blast_human (base_path)
            print('Human blast db created')
            index_db_foldseek_human_structures (base_path)
            print('Human foldseek db created')
            print('Human proteome downloaded and indexed')
        except Exception as e:
            print(f'Error downloading human proteome: {e}')  
    else:
        print('Human proteome already exists.')
        if not files.file_check(os.path.join(databases_path, 'HUMAN_DB.phr')):
            print('Indexing human proteome')
            index_db_blast_human (base_path)
            print('Human proteome indexed')
        else:
            print('Human proteome already indexed')

        human_PDB_db = os.path.join(human_structures_path, 'PDB_files', 'DB_foldseek', 'DB_human_PDB')
        human_AF_db = os.path.join(human_structures_path, 'AlphaFold_files', 'DB_foldseek', 'DB_human_AF')

        if not os.path.exists(human_PDB_db) or not os.path.exists(human_AF_db):
            print('Indexing human proteome for Foldseek')
            index_db_foldseek_human_structures (base_path)
            print('Human foldseek db created')
        else:
            print('Human foldseek db already created')

    print('----- 1. Finished -----')


    # Download and index Microbiome database
    print('----- 2. Downloading and indexing microbiome database -----')
    if not files.file_check(os.path.join(databases_path, 'uhgp-90.faa')):
        download_microbiome (base_path)
        index_db_blast_microbiome (base_path)
        print('Microbiome database downloaded and indexed')
    else:
        print('Microbiome database already exists.')
        if not files.file_check(os.path.join(databases_path, 'MICROBIOME_DB.00.phr')):
            print('Indexing microbiome database')
            index_db_blast_microbiome (base_path)
            print('Microbiome database indexed')
        else:
            print('Microbiome database already indexed')
    print('----- 2. Finished -----')

    # Download and index DEG database
    print('----- 3. Downloading and indexing DEG database -----')
    if not files.file_check(os.path.join(databases_path, 'DEG10.aa.faa')):
        download_DEG (base_path)
        index_db_blast_deg (base_path)
        print('DEG database downloaded and indexed')
    else:
        print('DEG database already exists.')
        if not files.file_check(os.path.join(databases_path, 'DEG_DB.phr')):
            print('Indexing DEG database')
            index_db_blast_deg (base_path)
            print('DEG database indexed')
        else:
            print('DEG database already indexed')
    print('----- 3. Finished -----')

    print('All databases downloaded and indexed successfully.')

if __name__ == '__main__':
    base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    main(base_path)

