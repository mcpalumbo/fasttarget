import subprocess
import docker
import os
import sys
import multiprocessing
from ftscripts import structures, files
import SNDG
import json
import sys
import glob
import pwd
import grp
from pathlib import Path
import shutil

def change_permission_user_file(file_path):
    """
    Change the permissions of a file to the current user.
    
    :param file_path: The file path.
    """

    username = os.getenv('SUDO_USER') if os.getenv('SUDO_USER') else os.getlogin()
    groupname = username
    uid = pwd.getpwnam(username).pw_uid
    gid = grp.getgrnam(groupname).gr_gid

    try:
        os.chmod(file_path, 0o777)
        os.chown(file_path, uid, gid)
        print(f"Permissions changed for file {file_path}.")
    except Exception as e:
        print(f"Failed to change permissions for file {file_path}: {e}")

def change_permission_user_dir(directory_path):
    """
    Change the permissions of a directory and its contents to the current user.
    
    :param directory_path: The directory path.
    """
    
    username = os.getenv('SUDO_USER') if os.getenv('SUDO_USER') else os.getlogin()
    groupname = username
    uid = pwd.getpwnam(username).pw_uid
    gid = grp.getgrnam(groupname).gr_gid

    for root, dirs, files in os.walk(directory_path):
        for dir_name in dirs:
            dir_path = os.path.join(root, dir_name)
            try:
                os.chmod(dir_path, 0o777)
                os.chown(dir_path, uid, gid)
                #print(f"Permissions changed for directory {dir_path}.")
            except Exception as e:
                print(f"Failed to change permissions for directory {dir_path}: {e}")
        
        for file_name in files:
            file_path = os.path.join(root, file_name)
            try:
                os.chmod(file_path, 0o777)
                os.chown(file_path, uid, gid)
                #print(f"Permissions changed for file {file_path}.")
            except Exception as e:
                print(f"Failed to change permissions for file {file_path}: {e}")

def load_config(base_path):
    """
    Load the configuration 'config.json' file.

    :param base_path: Base path where the repository data is stored.
    """
    
    config_path = f'{base_path}/config.json'
    with config_path.open() as file:
        config = json.load(file)
    return config

def run_bash_command(command):
    """
    Run a bash command.

    :param command: The command to run.
    """
    try:
        print(f'Running: {command}')
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        stdout, stderr = process.communicate()
        print(f"STDOUT:\n{stdout}")  # Print stdout regardless of success/failure
        print(f"STDERR:\n{stderr}")

        if process.returncode == 0:
            print("Command executed successfully.")
            if stdout:
                print(f"STDOUT:\n{stdout}")
            if stderr:
                print(f"STDERR:\n{stderr}")
        else:
            print(f"Command failed with return code {process.returncode}.")
            if stderr:
                print(f"STDERR:\n{stderr}")
            raise Exception(f"Command failed with return code {process.returncode}.")

    except Exception as e:
        print(f"An error occurred: {str(e)}")

def run_docker_container(work_dir, bind_dir, image_name, command, env_vars=None):
    """
    Run a docker container.

    :param work_dir: Working directory path.
    :param bind_dir: Binding directory path in container.
    :param image_name: The image to run.
    :param command: The command to run.
    :param env_vars: Dictionary of environment variables to set in the container.
    """
    client = docker.from_env()

    volumes = {work_dir: {'bind': bind_dir, 'mode': 'rw'}}
    user_str = f"{os.getuid()}:{os.getgid()}"

    try:
        print(f'Running docker image {image_name}, command: {command}')
        container = client.containers.run(
            image_name,
            command,
            volumes=volumes,
            working_dir= bind_dir,
            user=user_str,
            remove=True,
            stdout=True,
            stderr=True,
            environment=env_vars 
        )
        print(container.decode('utf-8'))
    except docker.errors.ContainerError as e:
        print(f"Error running container: {e}")

def run_fpocket(work_dir, pdb_file):

    """
    Run FPocket, using the docker image 'ezequieljsosa/fpocket'.
    
    :param work_dir:  Working directory path.
    :param pdb_file:  Structure file (.pdb) path.

    """

    if os.path.exists(pdb_file):
        FPOCKET_image = "ezequieljsosa/fpocket"
        FPOCKET_command = f"fpocket -f {pdb_file}"

        run_docker_container(work_dir, work_dir, FPOCKET_image, FPOCKET_command)
    else:
        print(f"The file '{pdb_file}' not found.", file=sys.stderr)

def run_fpocket2(work_dir, pdb_file):

    """
    Run FPocket, using the docker image 'ezequieljsosa/fpocket'.
    
    :param work_dir:  Working directory path.
    :param pdb_file:  Structure file (.pdb) path.
    """

    if os.path.exists(pdb_file):
        user = os.getuid()

        FPOCKET_BIN = f"docker run -v {work_dir}:{work_dir} -w {work_dir} --user {user}:{user} --rm ezequieljsosa/fpocket fpocket"
        FPOCKET_COMMAND = f"{FPOCKET_BIN} -f {pdb_file}"

        run_bash_command(FPOCKET_COMMAND)
    else:
        print(f"The file '{pdb_file}' not found.", file=sys.stderr)

def run_blastp(blastdb, query, output, evalue='1e-5', max_hsps='1', outfmt='6', max_target_seqs='500',cpus=multiprocessing.cpu_count()):

    """
    Runs Protein-Protein BLAST command line.
    
    :param blastdb: BLAST database name (full path).
    :param query: Query fasta (protein) file path.
    :param output: Output file path.
    :param evalue: Expect value (E) for saving hits. Deafult 1e-5.
    :param max_hsps: Maximum number of HSPs (alignments) to keep for any single query-subject pair. Dafault 1.
    :param outfmt: Output format. Default 6 (tabular).
    :param max_target_seqs: Number of aligned sequences to keep.
    :param cpus: Number of threads (CPUs) to use in blast search.
    
    """
    if files.file_check(query):

        blastp_command = f'blastp -evalue {evalue} -max_hsps {max_hsps} -outfmt "{outfmt}" -db {blastdb} -query {query} -num_threads {cpus} -max_target_seqs {max_target_seqs} > {output}'
        run_bash_command(blastp_command)
    else:
        print(f"Query file '{query}' not found.", file=sys.stderr)

def run_makeblastdb(input, output, title, dbtype):

    """
    Creates a BLAST database from command line.
    
    :param input: Input fasta file path.
    :param output: Output file path.
    :param title: Name of database.
    :param dbtype:  Molecule type of target db. Values: 'nucl' or 'prot'.
        
    """
    
    if files.file_check(input):
            makeblast_command = f'makeblastdb -in {input} -title {title} -out {output} -parse_seqids -dbtype {dbtype}'
            run_bash_command(makeblast_command)
    else:
        print(f"Blast database file '{input}' not found.", file=sys.stderr)

def run_genbank2gff3(input, output):

    """
    Runs script bp_genbank2gff3.pl from BioPerl.
    Convert a GenBank file to a gff3 for Roary.

    :param input: Input gbk file path.
    :param output: Output directory path.
        
    """
    work_dir = os.path.dirname(input)
    bind_dir = '/data'
    image_name = 'docker-bioperl'
    command = f'bp_genbank2gff3 {os.path.basename(input)}'

    if files.file_check(input):
        try:
            run_docker_container(
                work_dir=work_dir,
                bind_dir=bind_dir,
                image_name=image_name,
                command=command
            )
            file_output = f'{input}.gff'
            file_final = os.path.splitext(input)[0]+'.gff'
            if os.path.exists(file_output):
                print('bp_genbank2gff3.pl executed successfully.')
                if os.path.exists(output):
                    shutil.move(file_output, file_final)
                    if not os.path.exists(os.path.join(output, file_final)):
                        shutil.copy(file_final, output)
                        print(f'Gff3 file saved in {output}')
                else:
                    print(f"Directory '{output}' not found.", file=sys.stderr)
        except Exception as e:
            print(f'Error running bp_genbank2gff3.pl: {e}')
            print(f"An error occurred: {e}")
    else:
        print(f"GenBank file '{input}' not found.", file=sys.stderr)

def run_roary(work_dir:str, input:str, output:str, cpus=multiprocessing.cpu_count()):

    """
    Runs the docker image sangerpathogens/roary, a pan genome pipeline. Default options.
    More info: https://sanger-pathogens.github.io/Roary/

    :param work_dir: Directory where gff and roary output folders are.
    :param input: Directory where gff3 files are found.
    :param output: Output directory path.
    :param cpus: Number of threads.
    """

    if os.path.exists(input):

        result = subprocess.run([f'ls {input}/*.gff'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if result.returncode != 0:
            raise Exception(f"Error listing .gff files: {result.stderr.decode('utf-8')}")
        gff_files = result.stdout.decode('utf-8').strip().split()
        gff_files_str = " ".join(gff_files)

        ROARY_image = "sangerpathogens/roary"
        ROARY_command = f"roary -p {cpus} -f {output} {gff_files_str}"

        run_docker_container(work_dir, work_dir, ROARY_image, ROARY_command)
    else:
        print(f"Directory '{input}' not found.", file=sys.stderr)

def run_core_cruncher(ccruncher_script:str, input:str, output:str, reference:str):
    """
    Runs CoreCruncher, a core genome tool. Default options.
    More info: https://github.com/lbobay/CoreCruncher

    :param ccruncher_script: Path of corecruncher_master.py script.
    :param input: Path to a folder containing the genomes to analyze.
    :param output: CoreCruncher will create an output folder containg all the results of the analysis.
    :param reference: Pivot genome, specify the name of the file.
    
    """
   
    if os.path.exists(ccruncher_script):
        if os.path.exists(input):
            if os.path.exists(output):
            
                ccruncher_command = f'{ccruncher_script} -in {input} -out {output} -ref {reference}'
                run_bash_command(ccruncher_command)

            else:
                print(f"Directory '{output}' not found.", file=sys.stderr)
        else:
            print(f"Directory '{input}' not found.", file=sys.stderr)
    else:
        print(f"Script '{ccruncher_script}' not found.", file=sys.stderr)

def run_panx(panx_script, input, species_name, cpus=multiprocessing.cpu_count()):
    """
    Runs pan-genome-analysis, a pan-genome tool. Default options.
    More info: https://github.com/neherlab/pan-genome-analysis

    :param panx_script: Path of panX.py script.
    :param input: Path to a folder containing the genomes to analyze. It is the run directory.
    :param species_name: specie name. Used as prefix for some temporary folders (e.g.: P_aeruginosa).
    :param cpus: Number of threads.
    
    """
   
    if os.path.exists(panx_script):
        if os.path.exists(input):            
            panx_command = f'{panx_script} -fn {input}/ -sl {species_name} -t {cpus}'
            run_bash_command(panx_command)

        else:
            print(f"Directory '{input}' not found.", file=sys.stderr)
    else:
        print(f"Script '{panx_script}' not found.", file=sys.stderr)

def run_unzip(input_file, output_dir):
    
    """
    Unzip a file in a directory.
    
    :param input_file: File to unzip.
    :param output_dir: Directory where to unzip the file.
    """

    if os.path.exists(input_file):
        if os.path.exists(output_dir):            
            unzip_command = f'unzip {input_file} -d {output_dir}'
            run_bash_command(unzip_command)
            print(f'Unzip {input_file} in {output_dir}')
        else:
            print(f"Directory '{output_dir}' not found.", file=sys.stderr)
    else:
        print(f"File '{input_file}' not found.", file=sys.stderr)

def run_ncbi_datasets(tax_id, organism_name, output_dir):

    """
    Downloads genomes from ncbi for any taxonomy ID, in gbff format.
    Assembly level complete and annotated by RefSeq.
    More info: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

    :param tax_id: Taxonomy ID (eg. Pseudomonas aeruginosa taxonomy ID: 287).
    :param output_dir: Output directory path.
    :param organism_name: Name given to directories to download (eg. PAO).
    
    """

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")

    dehydrated_file =  os.path.join(output_dir, f'{organism_name}_dehydrated.zip')
    dehydrated_dir =  os.path.join(output_dir, f'{organism_name}_dataset')
    checkpoint_file =  os.path.join(output_dir, f'{organism_name}_checkpoint_ncbi_datasets.txt')

    if not files.file_check(checkpoint_file):
        if not os.path.exists(dehydrated_file):
            datasets_command = f"datasets download genome taxon {tax_id} --assembly-level complete --annotated --assembly-source 'RefSeq' --include gbff,gff3 --exclude-atypical --dehydrated  --filename {dehydrated_file}"
            run_bash_command(datasets_command)
            print(f'Download a dehydrated data package in {dehydrated_file}')
        else:
            print(f'Dehydrated data package in {dehydrated_file}')
        
        if not os.path.exists(dehydrated_dir):
            os.makedirs(dehydrated_dir)
            print(f"Created directory: {dehydrated_dir}")
            run_unzip(dehydrated_file, dehydrated_dir)
            print(f'Unzip a dehydrated data package in {dehydrated_dir}')
        else:
            print(f'Dehydrated data package in {dehydrated_dir}')
        
        rehydrate_command = f'datasets rehydrate --directory {dehydrated_dir}'
        run_bash_command(rehydrate_command)
        print(f'Rehydrated complete')
        
        print('NCBI download complete')

        with open(checkpoint_file, 'w') as f:
            f.write("Download complete: " + str(dehydrated_dir))
            f.close()
    else:
        print(f"Checkpoint file '{checkpoint_file}' already exists.", file=sys.stderr)
        print(f"NCBI datasets download already completed.", file=sys.stderr)

def run_ncbi_accession(accession, output_dir):

    """
    Downloads genomes from ncbi for any accession ID in gbff format.
    More info: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

    :param accession: Accession ID.
    :param output_dir: Output directory path.
    
    """

    if os.path.exists(output_dir):

        file_download = os.path.join(output_dir, f'{accession}.zip')

        datasets_command = f'datasets download genome accession {accession} --include gbff,gff3 --filename {file_download}'
        run_bash_command(datasets_command)
        run_unzip(file_download, output_dir)
        run_bash_command(f'rm {file_download} {output_dir}/README.md')

        old_accesion_dir = os.path.join(output_dir, accession)
        if os.path.exists(old_accesion_dir):
            run_bash_command(f'rm -r {old_accesion_dir}')

        run_bash_command(f'mv {output_dir}/ncbi_dataset/data/{accession} {output_dir}')
        run_bash_command(f'rm -r {output_dir}/ncbi_dataset')

    else:
        print(f"Directory '{output_dir}' not found.", file=sys.stderr)
    
    print(f'NCBI {accession} download complete in {output_dir }')

def run_ubiquitous(sbml_file, out_dir):
    
    """
    Generate a file with the ubiquitous compounds from a SBML file.
    
    :param sbml_file: SBML file path.
    :param out_dir: Output directory path.
    """
    
    if os.path.exists(sbml_file):
        if os.path.exists(out_dir):
            try:
                ubiq_command = f'python3 -m "SNDG.Network.SBMLProcessor" -i {sbml_file} -o {out_dir}/'
                run_bash_command(ubiq_command)
                print(f'Ubiquitous compounds file generated.')
            except Exception as e:
                print(f"An error occurred: {e}")
        else:
            print(f"Directory '{out_dir}' not found.", file=sys.stderr)
    else:
        print(f"SBML file '{sbml_file}' not found.", file=sys.stderr)
  
def run_sbml_to_sif(sbml_file, ubiquitous_file, out_dir):
    
    """
    Generate a SIF file from a SBML file.
    
    :param sbml_file: SBML file path.
    :param ubiquitous_file: Ubiquitous compounds file path.
    :param out_dir: Output directory path.
    """

    if os.path.exists(sbml_file):
        if os.path.exists(ubiquitous_file):
            if os.path.exists(out_dir):
                try:
                    sif_command = f'python3 -m "SNDG.Network.SBMLProcessor" -i {sbml_file} -o {out_dir}/ -f {ubiquitous_file}'
                    run_bash_command(sif_command)
                    print(f'Sif file generated.')
                except Exception as e:
                    print(f"An error occurred: {e}")
            else:
                print(f"Directory '{out_dir}' not found.", file=sys.stderr)
        else:
            print(f"Ubiquitous compounds file '{ubiquitous_file}' not found.", file=sys.stderr)
    else:
        print(f"SBML file '{sbml_file}' not found.", file=sys.stderr) 

def run_psort(input, organism_type, output_dir, output_format='terse'):
    """
    Runs PSORTb, a tool for predicting subcellular localization for a given set of protein sequences.
    More info: https://hub.docker.com/r/brinkmanlab/psortb_commandline

    :param input: Path to fasta file with protein sequences.
    :param organism_type: Type of organism, it can be Gram negative/positive bacteria or archaea. Only can take these values: n, p or a.
    :param output_dir:  Path of where to save results files.
    :param output_format: Format of output files. Value can be normal, terse or long. Default:terse.
    
    """

    valid_type = ['a','n','p']
    if organism_type not in valid_type:
        raise ValueError("output_dir must be one of %r." % valid_type)

    valid_format = ['normal','terse','long']
    if output_format not in valid_format:
        raise ValueError("output_format must be one of %r." % valid_format)   
   
    if os.path.exists(input):
        if os.path.exists(output_dir):          

            shutil.copy(input, output_dir)
            file_name = os.path.basename(input)

            image_name = 'brinkmanlab/psortb_commandline:1.0.2'
            psortb_command = f'/usr/local/psortb/bin/psort -{organism_type} -o {output_format} -i /tmp/results/{file_name}'
            env_vars = {'MOUNT': output_dir}

            run_docker_container(output_dir, '/tmp/results', image_name, psortb_command, env_vars)

            os.remove(os.path.join(output_dir, file_name))

            print(f"Psort results in '{output_dir}'.", file=sys.stderr)
        else:
            print(f"Directory '{output_dir}' not found.", file=sys.stderr)
    else:
        print(f"Directory '{input}' not found.", file=sys.stderr)

