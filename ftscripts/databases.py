import urllib.request
import gzip
import tarfile
import shutil
import os
from ftscripts import programs,files

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
    urllib.request.urlretrieve(url, deg_path)
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
    urllib.request.urlretrieve(info_url, info_path)

    if not files.file_check(os.path.join(databases_path, 'uhgp-90.faa')):
        uhgp_path = os.path.join(databases_path, 'uhgp-90.tar.gz')
        uhgp_url = 'https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/protein_catalogue/uhgp-90.tar.gz'
        
        print('Downloading UHGP database.')
        urllib.request.urlretrieve(uhgp_url, uhgp_path)
        print('Finished.')

        with tarfile.open(uhgp_path, 'r:gz') as tar:
            tar.extractall(path=databases_path)

        if os.path.exists(faa_file_path):
            shutil.move(faa_file_path, databases_path)
        
        shutil.rmtree(extracted_dir)

        print(f'File {uhgp_path} downloaded and decompressed successfully')
    else:
        print(f'{faa_file_path} already exists.')

def download_human(base_path):
    """
    Downloads Human proteome form Uniprot (UP000005640). 
    The fasta file is located in "databases" folder.

    :param base_path =  Path to the directory where 'databases' folder is located.

    """
    databases_path = os.path.join(base_path, 'databases')
    humanprot_path = os.path.join(databases_path, 'human_uniprot_UP000005640.faa')
    
    url = 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28proteome%3AUP000005640%29%29'
    
    print('Downloading human proteome.')
    urllib.request.urlretrieve(url, humanprot_path)
    print('Finished.')

    print(f'File {humanprot_path} downloaded successfully')

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
    if not files.file_check(os.path.join(databases_path, 'human_uniprot_UP000005640.faa')):
        download_human (base_path)
        index_db_blast_human (base_path)
        print('Human proteome downloaded and indexed')
    else:
        print('Human proteome already exists.')
        if not files.file_check(os.path.join(databases_path, 'HUMAN_DB.phr')):
            print('Indexing human proteome')
            index_db_blast_human (base_path)
            print('Human proteome indexed')
        else:
            print('Human proteome already indexed')
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

