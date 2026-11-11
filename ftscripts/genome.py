from Bio import SeqIO
from Bio.SeqFeature import ExactPosition
import os
import shutil
import tqdm
from ftscripts import programs, metadata, files
import multiprocessing
import glob
import csv
import pandas as pd
from BCBio import GFF
import sys
import subprocess
import logging
import json
import tempfile


def gbk_to_gff3(gbk_file, gff_dir):

    """
    Convert a GenBank file to GFF3 format.

    This function converts a GenBank file to GFF3 format and saves the output to a file in the specified directory.

    :param gbk_file: GenBank file path.
    :param gff_dir: Directory where the GFF3 file will be saved.
    """

    def fix_genbank_file(input_file, output_file):
        """
        Fix a GenBank file to be compatible with GFF3 format.

        This function fixes a GenBank file to be compatible with GFF3 format.
        It rejects locations/reading frames that cannot be safely converted.
        It also changes the type of features in pseudogenes.

        :param input_file: Input GenBank file path.
        :param output_file: Output GenBank file path.
        """

        
        fixed_records = []
        with open(input_file, "r") as in_handle, open(output_file, "w") as out_handle:
            for record in SeqIO.parse(in_handle, "genbank"):
                for feature in record.features:
                    if 'pseudo' in feature.qualifiers:
                        if feature.type == "gene":
                            feature.type = 'pseudogenic'
                        if feature.type == "CDS":
                            feature.type = 'pseudogene'
                    if feature.location is None:
                        raise ValueError(
                            f"Missing feature location in {gb_record.id}: {feature.type}"
                        )
                    location = feature.location
                    if (len(location.parts) != 1
                            or not isinstance(location.start, ExactPosition)
                            or not isinstance(location.end, ExactPosition)
                            or location.ref is not None
                            or location.ref_db is not None):
                        raise ValueError(
                            f"Python GFF fallback cannot safely represent location "
                            f"{location} in {gb_record.id}; use the BioPerl converter."
                        )
                    if "codon_start" in feature.qualifiers:
                        values = feature.qualifiers["codon_start"]
                        if len(values) != 1 or values[0] not in ('1', '2', '3'):
                            raise ValueError(
                                f"Invalid codon_start in {gb_record.id}: {values}"
                            )
                fixed_records.append(record)

            SeqIO.write(fixed_records, out_handle, "genbank")
            print(f"Fixed GenBank file saved to {output_file}")

    def write_gff3(gff_file, gb_records):
        """
        Write a list of GenBank records to a GFF3 file.

        :param gff_file: Output GFF3 file path.
        :param gb_records: List of GenBank records.

        """
        with open(gff_file, "w") as output_handle:
            GFF.write(gb_records, output_handle)

            # Add FASTA sequence at the end of the GFF file
            output_handle.write("##FASTA\n")
            for gb_record in gb_records:
                
                output_handle.write(f">{gb_record.id}\n")
                output_handle.write(str(gb_record.seq) + "\n")
                
                # Add CDS sequences
                for feature in gb_record.features:
                    if feature.type == "CDS":
                        if 'translation' in feature.qualifiers:
                            #cds_id = feature.qualifiers['locus_tag'][0]
                            cds_seq = feature.qualifiers['translation'][0]
                            cds_id = feature.qualifiers.get("locus_tag", ["unknown"])[0]

                            output_handle.write(f">{cds_id}.p01\n")
                            output_handle.write(str(cds_seq) + "\n")
        print(f"GFF3 file saved to {gff_file}")

    if os.path.exists(gff_dir):
        # Fix: Remove .gbk extension before adding .gff
        base_name = os.path.splitext(os.path.basename(gbk_file))[0]
        gff_file = os.path.join(gff_dir, f'{base_name}.gff')
        if os.path.exists(gbk_file):
            try: 
                fixed_gbk = os.path.join(os.path.dirname(gbk_file), f'fix_{os.path.basename(gbk_file)}')
                fix_genbank_file(gbk_file, fixed_gbk)
                records = []
                with open(fixed_gbk) as input_handle:
                    for gb_record in SeqIO.parse(input_handle, "genbank"):
                        records.append(gb_record)

                write_gff3(gff_file, records)
            except Exception as e:
                logging.exception(f"An error occurred during GBK to GFF3 conversion: {e}")
                raise  # Re-raise to let caller know it failed
        else:
            logging.error(f"GBK file '{gbk_file}' not found.")
    else:
        logging.error(f"Gff dir '{gff_dir}' not found.")

def add_sequences_to_gff3(gff_file, gbk_file, output_file=None):
    """Replace the FASTA section atomically, optionally into a separate GFF.

    Preserve the existing genome plus translated-protein (.p01) convention.
    Supplying output_file leaves the downloaded source unchanged.
    """
    destination = os.path.abspath(output_file or gff_file)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir=os.path.dirname(destination),
                                         prefix='.gff-', delete=False) as handle:
            temporary = handle.name
            with open(gff_file) as source:
                last_line = ''
                for line in source:
                    if line.strip() == '##FASTA':
                        break
                    handle.write(line)
                    last_line = line
                if last_line and not last_line.endswith('\n'):
                    handle.write('\n')
            handle.write('##FASTA\n')
            identifiers = set()

            def write_sequence(identifier, sequence):
                if (not identifier or any(c.isspace() for c in identifier)
                        or '>' in identifier or identifier in identifiers):
                    raise ValueError(f'Invalid or duplicate GFF FASTA ID: {identifier!r}')
                if not sequence or any(c.isspace() for c in sequence) or '>' in sequence:
                    raise ValueError(f'Invalid GFF FASTA sequence for {identifier}')
                identifiers.add(identifier)
                handle.write(f'>{identifier}\n{sequence}\n')

            with open(gbk_file) as source:
                for record in SeqIO.parse(source, 'genbank'):
                    write_sequence(record.id, str(record.seq))
                    for feature in record.features:
                        if feature.type != 'CDS' or 'translation' not in feature.qualifiers:
                            continue
                        loci = feature.qualifiers.get('locus_tag', [])
                        translations = feature.qualifiers['translation']
                        if len(loci) != 1 or not loci[0] or len(translations) != 1:
                            raise ValueError(f'Invalid translated CDS in {record.id}')
                        write_sequence(f'{loci[0]}.p01', translations[0])
            if not identifiers:
                raise ValueError(f'No GenBank sequences found in {gbk_file}')
        os.replace(temporary, destination)
        temporary = None
    finally:
        if temporary is not None:
            os.unlink(temporary)
    logging.info('Wrote GFF with one FASTA section: %s', destination)

def gbk_locus_strain_host(gbk_file):
    """
    Extracts the locus tag, strain, and host from a GenBank file.

    :param gbk_file: GenBank file path.
    
    :return: Locus tag list, strain, and host.
    """
    
    locus_tag_list = []
    strain = None
    host = None
    for record in SeqIO.parse(gbk_file, "genbank"):
        locus_tag = record.name
        locus_tag_list.append(locus_tag)
        print(locus_tag)
        for feature in record.features:
            if feature.type == "source":
                strain = feature.qualifiers.get("strain", ["unknown"])[0]
                strain = strain.replace(" ", "_")
                strain = strain.replace("/", "_")
                strain = strain.replace("(", "_")
                strain = strain.replace(")", "_")
                strain = strain.replace(":", "_")
                strain = strain.replace(";", "_")

                if 'host' in feature.qualifiers:
                    host = feature.qualifiers['host'][0]

    return locus_tag_list, strain, host

def gbk_to_fasta(gbk_file, output_file_fna=None, output_file_faa=None, output_file_ffn=None):
    """
    Generates the .faa, .fna, and/or .ffn files for a GenBank file.

    :param gbk_file: GenBank file path.
    :param output_file_fna: Output fna file path or None to skip.
    :param output_file_faa: Output faa file path or None to skip.
    :param output_file_ffn: Output ffn file path or None to skip.
    """

    # FAA and FFN describe the same universe: CDS with annotated translation.
    # Validate before opening any requested destination.
    if not any((output_file_fna, output_file_faa, output_file_ffn)):
        return
    records = list(SeqIO.parse(gbk_file, "genbank"))
    if not records:
        raise ValueError(f"No GenBank records found in {gbk_file}")
    record_ids = set()
    locus_ids = set()
    protein_records = []
    nucleotide_records = []
    skipped = 0

    def checked_id(values, label):
        if (len(values) != 1 or not values[0]
                or any(char.isspace() for char in values[0]) or '>' in values[0]):
            raise ValueError(f"Invalid or missing {label}: {values}")
        return values[0]

    def fasta_entry(identifier, sequence):
        return f">{identifier}\n" + "\n".join(
            sequence[i:i + 60] for i in range(0, len(sequence), 60)) + "\n"

    for record in records:
        record_id = checked_id([record.id], 'record ID')
        if record_id in record_ids:
            raise ValueError(f"Duplicate GenBank record ID: {record_id}")
        record_ids.add(record_id)
        if not (output_file_faa or output_file_ffn):
            continue
        for feature in record.features:
            if feature.type != 'CDS':
                continue
            if 'translation' not in feature.qualifiers:
                skipped += 1
                continue
            locus = checked_id(feature.qualifiers.get('locus_tag', []),
                               f'CDS locus_tag in {record.id}')
            if locus in locus_ids:
                raise ValueError(f"Duplicate translated CDS locus_tag: {locus}")
            locus_ids.add(locus)
            translations = feature.qualifiers['translation']
            if (len(translations) != 1 or not translations[0]
                    or any(char.isspace() for char in translations[0])
                    or '>' in translations[0]):
                raise ValueError(f"Invalid translation for CDS {locus}")
            if output_file_faa:
                protein_records.append(fasta_entry(locus, translations[0]))
            if output_file_ffn:
                if feature.location is None:
                    raise ValueError(f"Missing location for CDS {locus}")
                sequence = str(feature.extract(record.seq))
                if not sequence:
                    raise ValueError(f"Empty nucleotide sequence for CDS {locus}")
                nucleotide_records.append(fasta_entry(locus, sequence))
    if skipped:
        logging.warning('%s: excluded %d CDS without annotated translation from FAA/FFN',
                        gbk_file, skipped)
    if (output_file_faa or output_file_ffn) and not locus_ids:
        raise ValueError(f"No translated CDS available in {gbk_file}")
    if output_file_fna:
        with open(output_file_fna, 'w') as handle:
            SeqIO.write(records, handle, 'fasta')
    for destination, entries in ((output_file_faa, protein_records),
                                 (output_file_ffn, nucleotide_records)):
        if destination:
            with open(destination, 'w') as handle:
                handle.writelines(entries)

def ref_genome_files(gbk_file, output_path, organism_name, container_engine='docker'):
    """Prepare a reference generation, or reuse only byte-verified artifacts.

    Incompatible/legacy generations require a new organism output directory:
    rebuilding these files alone cannot invalidate all downstream caches.
    """
    if not os.path.isfile(gbk_file):
        raise FileNotFoundError(f"Reference GenBank not found: {gbk_file}")
    os.makedirs(output_path, exist_ok=True)
    names = {ext: f'{organism_name}.{ext}' for ext in ('gbk', 'fna', 'faa', 'ffn', 'gff')}
    manifest_path = os.path.join(output_path, f'{organism_name}.reference.json')
    source_hash = files.sha256_file(gbk_file)
    expected = {'schema': 3, 'source_sha256': source_hash,
                'container_engine': container_engine}
    existing = any(os.path.lexists(os.path.join(output_path, name)) for name in names.values())
    if existing or os.path.lexists(manifest_path):
        try:
            with open(manifest_path) as handle:
                manifest = json.load(handle)
            valid = isinstance(manifest, dict) and all(manifest.get(k) == v for k, v in expected.items())
            valid = valid and manifest.get('artifacts') == {
                name: files.sha256_file(os.path.join(output_path, name))
                for name in names.values()
            }
        except (OSError, ValueError):
            valid = False
        if not valid:
            raise RuntimeError(
                'Reference artifacts are unverified or incompatible with the current GenBank. '
                'Use a new organism output directory to avoid reusing stale downstream results. '
                f'Existing files were preserved: {output_path}'
            )
        logging.info('Reusing verified reference artifacts in %s', output_path)
        return

    # Publish the completion manifest last; interrupted publication cannot be reused.
    with tempfile.TemporaryDirectory(prefix='.reference-', dir=output_path) as staging:
        staged = {ext: os.path.join(staging, name) for ext, name in names.items()}
        shutil.copyfile(gbk_file, staged['gbk'])
        if files.sha256_file(staged['gbk']) != source_hash:
            raise RuntimeError('Reference GenBank changed during preparation')
        gbk_to_fasta(staged['gbk'], output_file_fna=staged['fna'],
                     output_file_faa=staged['faa'], output_file_ffn=staged['ffn'])
        try:
            programs.run_genbank2gff3(staged['gbk'], staging, container_engine=container_engine)
            if not files.file_check(staged['gff']):
                raise RuntimeError('GenBank conversion produced no nonempty GFF')
        except Exception as error:
            logging.warning('Container conversion failed; trying Python fallback: %s', error)
            if os.path.exists(staged['gff']):
                os.remove(staged['gff'])
            gbk_to_gff3(staged['gbk'], staging)
        if not all(files.file_check(path) for path in staged.values()):
            raise RuntimeError('Reference preparation produced missing or empty artifacts')
        manifest = dict(expected)
        manifest['artifacts'] = {names[ext]: files.sha256_file(path) for ext, path in staged.items()}
        for ext, path in staged.items():
            os.replace(path, os.path.join(output_path, names[ext]))
        files.atomic_write_json(manifest, manifest_path)

def id_to_locustag_gff(file_path:str, ids, fix=False):
    """

    Parses a .gff file, from a list of IDs retrieves a dictionary with the locus_tag for each ID.

    :param file_path: Path to the gff file.
    :param ids: List of IDs.
    :param fix: If True, work with IDs of roary index.

    :return: Dictionary with the locus_tag for each ID.

    """
    locus_tags = {}
    with open(file_path, 'r') as file:
        for rec in GFF.parse(file):
            for feature in rec.features:
                if 'ID' in feature.qualifiers:
                    feature_id = feature.qualifiers['ID'][0]
                    for id in ids:
                        if fix:
                            if feature_id == id and feature_id not in locus_tags:
                                locus_tags[feature_id] = feature.qualifiers.get('locus_tag', [feature_id])[0]
                        else:
                            if feature_id in id and feature_id not in locus_tags:
                                locus_tags[feature_id] = feature.qualifiers.get('locus_tag', [feature_id])[0]
    
    list_locus = list(locus_tags.values())
    return list_locus

def core_download_genomes_ncbi(output_path, organism_name, tax_id, accession_list=None):
    """
    Download genomes from NCBI.
    Can download either by taxonomy ID (all genomes) or from a specific list of accessions.
    This function uses the `run_ncbi_datasets` or `run_ncbi_datasets_accessions` function from the `programs` module.

    :param output_path: Directory of the oraganism output.
    :param organism_name: Name of the organism.
    :param tax_id: Taxonomy ID from ncbi.
    :param accession_list: Optional list of specific assembly accessions to download. If provided, only these genomes will be downloaded.
    """

    conservation_dir = os.path.join(output_path, organism_name, 'conservation')
    
    if accession_list and len(accession_list) > 0:
        print(f'Downloading {len(accession_list)} specific genomes from accession list')
        programs.run_ncbi_datasets_accessions(accession_list, organism_name, conservation_dir)
    else:
        print(f'Downloading all genomes from taxonomy ID: {tax_id}')
        programs.run_ncbi_datasets(tax_id, organism_name, conservation_dir)

def core_check_genomes_ncbi(output_path, organism_name):
    """

    Check genomes downloaded from NCBI. 
    Parse assembly_data_report.jsonl and look for .gbff and .gff files not downloaded.

    :param output_path: Directory of the oraganism output.
    :param organism_name: Name of the organism.

    :return: List of missing genomes.
    """

    conservation_dir = os.path.join(output_path, organism_name, 'conservation')
    ncbi_download_data = os.path.join(conservation_dir, f'{organism_name}_dataset', 'ncbi_dataset', 'data')

    assembly_json = os.path.join(ncbi_download_data,'assembly_data_report.jsonl')
    assembly_dict = files.jsonl_to_dict(assembly_json)

    required_files = ['*.gbff', '*.gff']
    all_accessions = []
    missing_files_accessions = []

    for item in assembly_dict:
        accession = item['accession']
        all_accessions.append(accession)
        
        accession_dir = os.path.join(ncbi_download_data, accession)
        
        for acc in required_files:
            if os.path.exists(accession_dir):
                files_gff = glob.glob(os.path.join(accession_dir, acc))
                if files_gff:
                    print(f'{acc} file found for {accession}')
                else:
                    missing_files_accessions.append(accession)
                    print(f'{acc} file from {accession} is missing.')
            else:
                missing_files_accessions.append(accession)
                print(f'{accession} is missing.')

    if len(missing_files_accessions) == 0:
        print(f'All genomes downloaded correctly.')
    else:
        print(f'Missing genomes: {missing_files_accessions}.')

    return missing_files_accessions
        
def core_download_missing_accessions(output_path, organism_name, tax_id, accession_list=None):
    """
    Check for missing genome files and download them if necessary.

    This function checks if there are any missing genome files for the given organism.
    If there are missing files, it attempts to download them from NCBI.

    :param output_path: Directory of the oraganism output.
    :param organism_name: Name of the organism.
    :param tax_id: Taxonomy ID of the organism.
    :param accession_list: Optional list of specific assembly accessions. If provided, only these genomes are expected.
    """

    conservation_dir = os.path.join(output_path, organism_name, 'conservation')
    ncbi_download_data = os.path.join(conservation_dir, f'{organism_name}_dataset', 'ncbi_dataset', 'data')
    checkpoint_file =  os.path.join(ncbi_download_data, f'{organism_name}_checkpoint_check_datasets.txt')

    if not files.file_check(checkpoint_file):
        # Check for missing genome files
        missing_files_accessions = core_check_genomes_ncbi(output_path, organism_name)
        if len(missing_files_accessions) > 0:
            # Download missing genome files
            core_download_genomes_ncbi(output_path, organism_name, tax_id, accession_list=accession_list)

        # Re-check for missing genome files after download attempt
        missing_files_accessions = core_check_genomes_ncbi(output_path, organism_name)
        if len(missing_files_accessions) > 0:
            print(f'Missing genomes: {missing_files_accessions}.')
            print(len(missing_files_accessions)) 
            # Uncomment the following lines to download each missing accession
            # for accession in missing_files_accessions:
            #     programs.run_ncbi_accession(accession, ncbi_download_data)
            #     print(f'{accession} downloaded.')
        
        with open(checkpoint_file, 'w') as f:
            f.write("Check core genomes complete. Missing genomes:"+ str(set(missing_files_accessions)))
            f.close()
    else:
        print('Check already performed.')

def core_files(output_path, organism_name, container_engine='docker', filter_by_host=True):
    """
    Select genomes for core genome analysis and generate .faa and .gff3 files.

    This function processes the genomes of the given organism.
    If filter_by_host=True, it selects only those with Human as host.
    If filter_by_host=False, it includes all genomes that don't duplicate loci or strain.
    It generates the necessary .faa and .gff3 files for each genome.

    :param output_path: Directory of the oraganism output.
    :param organism_name: Name of the organism.
    :param container_engine: Container engine to use ('docker' or 'singularity').
    :param filter_by_host: If True, filter by 'Homo sapiens' host. If False, use all genomes.

    :return: List of core genomes IDs.
    """

    conservation_dir = os.path.join(output_path, organism_name, 'conservation')
    ncbi_download_data = os.path.join(conservation_dir, f'{organism_name}_dataset', 'ncbi_dataset', 'data')
           
    #Output directory for Roary, for input takes gff files
    gff_dir = os.path.join(conservation_dir, 'roary_output','gff')
    os.makedirs(gff_dir, exist_ok=True)
    #Output directory for CoreCruncher, for input takes faa files
    fasta_dir = os.path.join(conservation_dir, 'corecruncher_output','faa')
    os.makedirs(fasta_dir, exist_ok=True)
    
    print('Conservation gff and faa folders created.')
    core_genomes = [organism_name]

    #Reference genome
    ref_gbk = os.path.join(output_path, organism_name, 'genome', f'{organism_name}.gbk')
    ref_gff = os.path.join(output_path, organism_name, 'genome', f'{organism_name}.gff')
    ref_faa = os.path.join(output_path, organism_name, 'genome', f'{organism_name}.faa')

    if not os.path.exists(os.path.join(gff_dir, f'{organism_name}.gff')):
        shutil.copy(ref_gff, gff_dir)
    if not os.path.exists(os.path.join(fasta_dir, f'{organism_name}.faa')):
        shutil.copy(ref_faa, fasta_dir)

    ref_locus_list, ref_strain, ref_host = gbk_locus_strain_host(ref_gbk)

    #Core genomes from NCBI
    for root, dirs, files_names in os.walk(ncbi_download_data):
        for i, dir_name in enumerate(dirs):
            assembly_id = str(dir_name).split(".")[0]

            print(f'Genome {i} from {len(dirs)}: {dir_name}')
            gbff_pattern = os.path.join(root, dir_name, "*.gbff")
            gbff_files = glob.glob(gbff_pattern)
            gbk_pattern = os.path.join(root, dir_name, "*.gbk")
            gbk_files = glob.glob(gbk_pattern)
            genome_files = gbff_files if len(gbff_files) == 1 else gbk_files

            gff_pattern = os.path.join(root, dir_name, "*.gff")
            gff_files = glob.glob(gff_pattern)
    
            if len(genome_files) == 1:

                locus_list, strain, host = gbk_locus_strain_host(genome_files[0])

                # Check if genome should be included based on filtering criteria
                include_genome = not set(locus_list) & set(ref_locus_list) and strain != ref_strain
                
                if filter_by_host and include_genome:
                    # When filtering by host, only include 'Homo sapiens'
                    include_genome = host == 'Homo sapiens'
                    if not include_genome:
                        print(f'{dir_name}: Host not Human: {host}, for {assembly_id}, {locus_list}, strain {strain}.')
                
                if include_genome:
                    core_genomes.append(f'{assembly_id}')  
                    
                    new_gbk = os.path.join(root, dir_name, f'{assembly_id}.gbk')
                    new_faa = os.path.join(fasta_dir, f'{assembly_id}.faa')
                    new_gff = os.path.join(gff_dir, f'{assembly_id}.gff')

                    if not os.path.exists(new_gbk):
                        if len(gbff_files) == 1:
                            shutil.move(gbff_files[0], new_gbk)
                        else:
                            shutil.copy(genome_files[0], new_gbk)
                    else:
                        print(f'{new_gbk} already exists.')

                    if not os.path.exists(new_faa):
                        gbk_to_fasta(new_gbk, output_file_faa=new_faa)
                    else:
                        print(f'{new_faa} already exists.')

                    if not os.path.exists(new_gff):
                        try:      
                            programs.run_genbank2gff3(new_gbk, gff_dir, container_engine=container_engine)
                            print(f"Processed {assembly_id}.gff")
                        except Exception as e:
                            logging.exception(f"Error in run_genbank2gff3: {e}")
                            if len(gff_files) == 1:
                                add_sequences_to_gff3(gff_files[0], new_gbk, output_file=new_gff)
                                print(f"Processed {assembly_id}.gff")
                            else:
                                logging.error(
                                    f"Expected exactly one GFF fallback file for {assembly_id}, found {len(gff_files)}."
                                )
                    else:
                        print(f"{assembly_id}.gff already exists.")
                        shutil.rmtree(os.path.join(root, dir_name))
                else:
                    shutil.rmtree(os.path.join(root, dir_name))
            elif len(gbff_files) == 0 and len(gbk_files) == 0:
                logging.error(f"No .gbff or .gbk file found for {assembly_id} in {os.path.join(root, dir_name)}")
            else:
                print(f'More than one genome files.')

    print("Files processed successfully.")
    files.list_to_file(os.path.join(conservation_dir, 'core_genomes_IDs.txt'), core_genomes)
        
    return core_genomes


def core_check_files(output_path, organism_name, container_engine='docker', filter_by_host=True):
    """

    Check for missing files in the core genomes analysis.

    This function checks if the necessary files for the core genomes analysis are present.
    It checks for the presence of the .gff and .faa files.
    If pending files are found, it runs the `core_files` function to generate the missing files.

    :param output_path: Directory of the oraganism output.
    :param organism_name: Name of the organism.
    :param container_engine: Container engine to use ('docker' or 'singularity').
    :param filter_by_host: If True, filter genomes by 'Homo sapiens' host. If False, use all genomes.

    """
    conservation_dir = os.path.join(output_path, organism_name, 'conservation')  
    genomes_file = os.path.join(conservation_dir,'core_genomes_IDs.txt')
    fasta_dir = os.path.join(conservation_dir, 'corecruncher_output','faa')
    gff_dir = os.path.join(conservation_dir, 'roary_output','gff')

    def find_pending_files(genomes_list):

        pending_files = []

        for genome in tqdm.tqdm(genomes_list, initial=1):
            gff_file = os.path.join(gff_dir, f"{genome}.gff")
            faa_files = os.path.join(fasta_dir, f"{genome}.faa")

            if not os.path.exists(gff_file):
                pending_files.append(gff_file)
                print(f'{gff_file} not found.')

            if not os.path.exists(faa_files):
                pending_files.append(faa_files)
                print(f'{faa_files} not found.')
        
        if len(pending_files) == 0:
            print('All core files correct.')
        else:
            print('Error: Some core files are missing.')
            print(f'Number of pending files: {len(pending_files)}')
        
        return pending_files

    # Previously this returned the inner function itself on the failure path,
    # because the function name was rebound by the assignment below. Worse, the
    # core programs then ran against missing inputs. Without a genome list there
    # is nothing to check, so fail with a message that says that.
    if not files.file_check(genomes_file):
        core_files(output_path, organism_name, container_engine=container_engine, filter_by_host=filter_by_host)

        if not files.file_check(genomes_file):
            raise FileNotFoundError(
                f'No core genomes file found or empty: {genomes_file}. '
                f'Cannot verify the core genome input files.'
            )

    genomes_list = files.file_to_list(genomes_file)
    # The caller only logged 'Genomes filtered', with no indication of how many
    # survived, which is the number the core genome is computed over.
    print(f'Genomes selected for the core genome analysis: {len(genomes_list)}')
    pending = find_pending_files(genomes_list)

    if pending:
        core_files(output_path, organism_name, container_engine=container_engine, filter_by_host=filter_by_host)
        pending = find_pending_files(genomes_list)
        if pending:
            logging.warning(
                f'{len(pending)} core input file(s) are still missing after retrying, '
                f'e.g. {pending[:3]}'
            )

    return pending

def report_core_table(table, column, program, core_percent, table_path):

    """
    Report the size of a core genome read back from a cached table.

    :param table: Table with one boolean column per gene.
    :param column: Name of the boolean column.
    :param program: Name of the program that produced it.
    :param core_percent: Frequency threshold the table was built with.
    :param table_path: Path of the table, used to read the stored genome count.
    """

    manifest = files.read_manifest(table_path) or {}
    genomes = manifest.get('extra', {}).get('genomes_analyzed')
    genome_note = f' over {genomes} genomes' if genomes else ''

    if column in table.columns:
        core_genes = int(table[column].astype(str).str.lower().eq('true').sum())
        print(
            f'{program} core genes (> {core_percent:g}% strains){genome_note}: '
            f'{core_genes} of {len(table)} genes (from cache)'
        )
    else:
        print(f'{program} table read from cache: {len(table)} genes, column {column} absent')

def _check_core_run_parameters(output_file, parameters, program, stale_dir):

    """
    Refuse to reuse a core genome run that used different parameters.

    Roary and CoreCruncher take hours to days, so their outputs are never
    recomputed automatically: a changed parameter stops the run and names the
    directory to remove, instead of either silently reusing the old clustering or
    deleting expensive results on its own.

    A parameter the manifest does not mention is unknown rather than wrong, which
    is the case for an output built by an earlier version of FastTarget: those are
    reused with a warning. Only a recorded value that differs stops the run.

    :param output_file: Path of the output file produced by the program.
    :param parameters: Mapping of parameter name to the value for this run.
    :param program: Name of the program, used in the messages.
    :param stale_dir: Directory to remove in order to force a re-run.
    :raises RuntimeError: If the recorded parameters differ from the current ones.
    """

    unrecorded = files.manifest_unrecorded(output_file, parameters)
    if unrecorded:
        known = 'nothing' if files.read_manifest(output_file) is None else 'only some parameters'
        print(
            f'Warning: {output_file} records {known} about how it was built '
            f'({", ".join(unrecorded)} unknown). Reusing it; remove {stale_dir} to '
            f'rebuild it with core frequency {parameters["min_core_freq"]}% and '
            f'identity {parameters["min_identity"]}% over the current genomes.'
        )

    differences = files.manifest_differences(output_file, parameters)
    if differences:
        # The genome set is not a config value to put back, so say what actually
        # changed rather than offering a remedy that does not apply.
        if any(difference.startswith('genome_') for difference in differences):
            remedy = (
                f'The set of genomes changed, so the old clustering no longer '
                f'describes it. Remove {stale_dir} to re-run {program} over the '
                f'current genomes, or restore the previous core.accession_file.'
            )
        else:
            remedy = (
                f'Remove {stale_dir} to re-run {program}, or restore the previous '
                f'values in the config.'
            )
        raise RuntimeError(
            f'{program} output in {stale_dir} was built with different parameters '
            f'({"; ".join(differences)}). {remedy}'
        )

def core_genome_programs(output_path, organism_name, core_threshold=99, identity=95, cpus=multiprocessing.cpu_count(), program_list=['roary','corecruncher'], container_engine='docker'):

    """
    Run Roary and CoreCruncher programs for the core genomes analysis of the given organism.

    This function runs Roary and CoreCruncher programs for the core genomes of the given organism to identify the core genes.
    It generates the necessary input files and runs the programs.

    This function uses the `run_roary` and `run_core_cruncher` functions from the `programs` module.

    :param output_path: Directory of the oraganism output.
    :param organism_name: Name of the organism.
    :param cpus: Number of CPUs to use.
    :param program_list: List of programs to run, can be 'roary', 'corecruncher', or both.
    :param container_engine: Container engine to use ('docker' or 'singularity').
    
    """

    conservation_dir = os.path.join(output_path, organism_name, 'conservation')
    gff_dir = os.path.join(conservation_dir, 'roary_output','gff')
    fasta_dir = os.path.join(conservation_dir, 'corecruncher_output','faa')

    def core_parameters_for(input_dir, pattern):
        """
        Parameters identifying one core genome run: the thresholds and the input set.

        The genome set comes from core.accession_file or organism.tax_id, so it is
        as much an input as the thresholds are: adding genomes to the list and
        re-running used to reuse a clustering computed over the old set. The names
        are reduced to a digest, since there can be thousands of them.
        """
        genomes = glob.glob(os.path.join(input_dir, pattern))
        return {
            'min_core_freq': core_threshold,
            'min_identity': identity,
            'genome_count': len(genomes),
            'genome_set': files.fileset_fingerprint(genomes),
        }

    if 'roary' in program_list:
        roary_out_dir = os.path.join(conservation_dir, 'roary_output')
        results_dir = os.path.join(roary_out_dir, 'results')
        genes_csv_file = os.path.join(results_dir, 'gene_presence_absence.csv')
        if os.path.exists(gff_dir):
            core_parameters = core_parameters_for(gff_dir, '*.gff')
            if not files.file_check(genes_csv_file):
                print('Running Roary')
                # These used to be hardcoded to 99/95, so core.min_core_freq and
                # core.min_identity from the config had no effect on the Roary run.
                print(f'Roary parameters: core frequency {core_threshold}%, identity {identity}%')
                programs.run_roary(roary_out_dir, gff_dir, results_dir, core_threshold=core_threshold, identity=identity, cpus=cpus, container_engine=container_engine)
                files.write_manifest(genes_csv_file, core_parameters, extra={'program': 'roary'})
                print('Finished')
            else:
                _check_core_run_parameters(genes_csv_file, core_parameters, 'Roary', results_dir)
                print(f'Roary output file already exists: {genes_csv_file}.')
        else:
            print(f'{gff_dir} not found.')

    if 'corecruncher' in program_list:
        ccruncher_out_dir = os.path.join(conservation_dir, 'corecruncher_output')
        cc_output_file = os.path.join(ccruncher_out_dir, 'families_core.txt')

        if os.path.exists(fasta_dir):
            core_parameters = core_parameters_for(fasta_dir, '*.faa')
            if not files.file_check(cc_output_file):
                reference_file = f'{organism_name}.faa'
                print('Running CoreCruncher')
                # Same here: the config values were received and then ignored.
                print(f'CoreCruncher parameters: core frequency {core_threshold}%, identity {identity}%')
                programs.run_core_cruncher(ccruncher_out_dir, reference_file, core_threshold=core_threshold, identity=identity, container_engine=container_engine)
                files.write_manifest(cc_output_file, core_parameters, extra={'program': 'corecruncher'})
                print('Finished')
            else:
                _check_core_run_parameters(cc_output_file, core_parameters, 'CoreCruncher', ccruncher_out_dir)
                print(f'CoreCruncher output file already exists: {cc_output_file}.')
        else:
            print(f'{fasta_dir} not found.')

def roary_output(output_path, organism_name, core_threshold=0.99):
    """    
    This function processes the output files generated by Roary for the specified organism.
    It reads the 'gene_presence_absence.csv' file to extract the core locus tags and generates a metadata table with these IDs.

    :param output_path: Directory of the oraganism output.
    :param organism_name: Name of the organism.
    :core_threshold: Threshold to consider a gene as core (default is 0.99 for 99%).
    
    :return: A list containing the core locus tags and a true/false table.
    """
    conservation_dir = os.path.join(output_path, organism_name, 'conservation')
    roary_out_dir = os.path.join(conservation_dir, 'roary_output')
    gff_dir = os.path.join(roary_out_dir,'gff')
    results_dir = os.path.join(roary_out_dir, 'results')
    roary_results_table = os.path.join(conservation_dir, 'core_roary.tsv')

    # Each missing input used to fall through to `return roary_table` with the
    # variable never assigned, surfacing as UnboundLocalError instead of naming
    # the file that was absent.
    if not os.path.exists(results_dir):
        raise FileNotFoundError(f'No Roary output found in {roary_out_dir}.')

    fixed_input_dir = os.path.join(results_dir, 'fixed_input_files')
    genes_csv_file = os.path.join(results_dir, 'gene_presence_absence.csv')

    if not os.path.exists(genes_csv_file):
        raise FileNotFoundError(
            f'Roary did not produce "gene_presence_absence.csv" in {results_dir}.'
        )

    # The threshold is applied here, at parse time, so a cached table built with a
    # different one does not describe the requested core genome. Rebuilding it only
    # costs a pandas filter, unlike the Roary run itself.
    table_parameters = {'core_threshold': core_threshold}

    if files.manifest_matches(roary_results_table, table_parameters, label='core_roary.tsv'):
        print(f'Roary output file already exists: {roary_results_table}.')
        cached_table = pd.read_csv(roary_results_table, sep='\t', index_col=0, header=0)
        # Reporting the counts only when the table is computed meant a resumed run
        # said nothing about the size of the core genome, which is the number the
        # run exists to produce.
        report_core_table(cached_table, 'core_roary', 'Roary', core_threshold * 100, roary_results_table)
        return cached_table

    dtypes = {'No. isolates': 'int', f'{organism_name}': 'str'}
    df = pd.read_csv(genes_csv_file, usecols= ['No. isolates', f'{organism_name}'], dtype=dtypes)
    max_isolates = df['No. isolates'].max()
    print(f'Total genomes analyzed by Roary: {max_isolates}')
    max_isolates_rows = df[df['No. isolates'] >= int(max_isolates)*core_threshold]
    print(f'Roary total core genes (> {core_threshold*100}% strains): {len(max_isolates_rows)}')
    gbk_ids = max_isolates_rows[f'{organism_name}'].dropna().tolist()
    print(f'Roary {organism_name} core genes: {len(gbk_ids)}')

    gff_file = os.path.join(gff_dir, f'{organism_name}.gff')
    gff_file_fixed = os.path.join(fixed_input_dir, f'{organism_name}.gff')

    if os.path.exists(gff_file_fixed):
        print(f'Reading {gff_file_fixed}')
        core_locus_tag = id_to_locustag_gff(gff_file_fixed, gbk_ids, fix=True)
    elif os.path.exists(gff_file):
        print(f'Reading {gff_file}')
        core_locus_tag = id_to_locustag_gff(gff_file, gbk_ids)
    else:
        # Without the GFF the Roary IDs cannot be translated to locus tags, so
        # there is no way to build the table. This used to leave core_locus_tag
        # unassigned and fail on the next line with no indication of the cause.
        raise FileNotFoundError(
            f'{organism_name}.gff not found in {fixed_input_dir} or {gff_dir}; '
            f'cannot map Roary gene IDs to locus tags.'
        )

    roary_table = metadata.metadata_table_bool(output_path, organism_name, core_locus_tag, 'core_roary', conservation_dir)
    files.write_manifest(
        roary_results_table,
        table_parameters,
        extra={'genomes_analyzed': int(max_isolates), 'core_genes': len(gbk_ids)},
    )

    return roary_table

def corecruncher_output(output_path, organism_name, core_threshold=None):
    
    """
    Process the output from CoreCruncher for the given organism.

    This function processes the output file 'families_core.txt' generated by CoreCruncher for the specified organism.
    It reads the file to extract the core locus tags and generates a metadata table with these IDs.

    :param output_path: Directory of the oraganism output.
    :param organism_name: Name of the organism.
    
    :return: A list containing the core locus tags and a true/false table.
    """
        
    conservation_dir = os.path.join(output_path, organism_name, 'conservation')
    cc_output_file = os.path.join(conservation_dir, 'corecruncher_output', 'families_core.txt')
    ref_genome = f'{organism_name}.faa'
    cc_results_table = os.path.join(conservation_dir, 'core_corecruncher.tsv')

    core_locus_tag = []
    core_total = 0

    if not files.file_check(cc_results_table):
        # A missing output file used to be logged and then ignored, leaving
        # core_locus_tag empty and producing a table that declares no gene core.
        # That is wrong data rather than a failure, so refuse instead.
        if not os.path.exists(cc_output_file):
            raise FileNotFoundError(
                f"CoreCruncher output file '{cc_output_file}' not found; "
                f'cannot determine the core genome.'
            )

        with open(cc_output_file, 'r') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            for row in reader:
                core_total += 1
                for value in row:
                    if isinstance(value, str) and value.startswith(ref_genome):
                        core_locus_tag.append(value.split('&')[1])

        # CoreCruncher already applied the frequency threshold when it ran, so this
        # is only a label. It used to say 90% unconditionally, which did not match
        # the value the run actually used.
        threshold_label = f'{core_threshold}%' if core_threshold is not None else 'the configured frequency'
        print(f'CoreCruncher total core genes (> {threshold_label} strains): {core_total}')
        print(f'CoreCruncher {organism_name} core genes: {len(core_locus_tag)}')

        corecruncher_table = metadata.metadata_table_bool(output_path, organism_name, core_locus_tag, 'core_corecruncher',conservation_dir)
        # Unlike core_roary.tsv, this parse applies no threshold of its own: the
        # frequency was applied by CoreCruncher when it ran, and a changed value is
        # refused at that gate. The manifest is provenance, not a cache key.
        files.write_manifest(
            cc_results_table,
            {'core_threshold': core_threshold},
            extra={'core_families': core_total, 'core_genes': len(core_locus_tag)},
        )
    else:
        print(f'CoreCruncher output file already exists: {cc_results_table}.')
        corecruncher_table = pd.read_csv(cc_results_table, sep='\t', index_col=0, header=0)
        threshold_percent = core_threshold if core_threshold is not None else float('nan')
        report_core_table(
            corecruncher_table, 'core_corecruncher', 'CoreCruncher',
            threshold_percent, cc_results_table,
        )

    return corecruncher_table

GRAM_TYPE_NAMES = {'n': 'Gram-negative', 'p': 'Gram-positive', 'a': 'archaea'}


def _check_psortb_gram_type(psortb_file, organism_type):

    """
    Refuse to reuse a PSORTb result produced for a different organism type.

    Unlike a threshold, the gram type is not a value to sweep: a Gram-negative
    prediction reused for a Gram-positive organism is simply the wrong answer.
    The file is not removed automatically, both because PSORTb takes a while and
    because leaving two results behind would trip the "leave only one" check.

    :param psortb_file: Path of the existing PSORTb output file.
    :param organism_type: Gram type for this run ('n', 'p' or 'a').
    :raises RuntimeError: If the file was produced for another gram type.
    """

    requested = GRAM_TYPE_NAMES.get(organism_type, organism_type)

    if files.manifest_unrecorded(psortb_file, {'gram_type': organism_type}):
        print(
            f'Warning: {psortb_file} does not record the gram type it was produced '
            f'with. Reusing it as a {requested} prediction; remove it to re-run PSORTb.'
        )
        return

    differences = files.manifest_differences(psortb_file, {'gram_type': organism_type})
    if differences:
        recorded = (files.read_manifest(psortb_file) or {}).get('parameters', {}).get('gram_type')
        raise RuntimeError(
            f'{psortb_file} was produced for '
            f'{GRAM_TYPE_NAMES.get(recorded, recorded)}, but psortb.gram_type is now '
            f'{organism_type} ({requested}). Remove that file to re-run PSORTb, or '
            f'restore the previous gram type in the config.'
        )

def localization_prediction(output_path, organism_name, organism_type, container_engine='docker'):

    """
    Predicts the subcellular localization of proteins using PSORTb.
    
    This function predicts the subcellular localization of proteins for the given organism using PSORTb.
    It first checks if a localization prediction has already been performed and saved to disk.
    If not, it runs PSORTb to predict the localization and saves the results to disk.
    This function uses the `run_psort` function from the `programs` module.

    :param output_path: Directory of the oraganism output.
    :param organism_name: Name of the organism.
    :param organism_type: Type of organism, can be 'a' (archaea), 'n' (gram-negative bacteria), or 'p' (gram-positive bacteria).
    :param container_engine: Container engine to use ('docker' or 'singularity').

    :return: A dictionary with the localization predictions for each protein.
    
    """
    
    valid_type = ['a','n','p']
    if organism_type not in valid_type:
        raise ValueError("output_dir must be one of %r." % valid_type)
    
    faa_path = os.path.join(output_path, organism_name, 'genome', f'{organism_name}.faa')
    localization_dir = os.path.join(output_path, organism_name, 'localization')

    # PSORTb names its own output after the gram type it was run with, for example
    # 20261111021226_psortb_gramneg.txt. The gram type is in the file name but not
    # in this pattern, so a result for the wrong organism type used to satisfy the
    # gate: switching psortb.gram_type from 'n' to 'p' over an existing output
    # directory reused the Gram-negative predictions without a word. The run is
    # recorded in a manifest instead of being inferred from the name, which varies
    # by PSORTb version.
    file_pattern = os.path.join(localization_dir, "*_psortb_*.txt")
    all_files = sorted(glob.glob(file_pattern))

    if not all_files:
        programs.run_psort(faa_path, organism_type, localization_dir, output_format='terse', container_engine=container_engine)
        # After running, find the newly created file
        all_files = sorted(glob.glob(file_pattern))
        if not all_files:
            raise FileNotFoundError(f"PSORTb ran but no output file found in {localization_dir}")
        files.write_manifest(all_files[0], {'gram_type': organism_type})
    elif len(all_files) == 1:
        _check_psortb_gram_type(all_files[0], organism_type)

    if len(all_files) > 1:
        raise FileNotFoundError(f"Multiple psort result files found in {localization_dir}. Leave only one.")

    # Process PSORTb results
    print(f"Processing PSORTb results from {all_files[0]}")

    psort_results = os.path.join(localization_dir, 'psortb_localization.tsv')
    psort_parameters = {'gram_type': organism_type}

    if not files.manifest_matches(psort_results, psort_parameters, label='psortb_localization.tsv'):
        # Parse the PSORTb output file
        file_path = all_files[0]
        df = pd.read_csv(file_path, sep=r'\s+', usecols=[0, 1])

        result_dict = df.set_index('SeqID')['Localization'].to_dict()

        psort_df = metadata.metadata_table_with_values(output_path, organism_name, result_dict,'psortb_localization',localization_dir, 'Unknown')
        files.write_manifest(
            psort_results,
            psort_parameters,
            extra={'source': os.path.basename(file_path), 'proteins': len(result_dict)},
        )

        return psort_df
    else:
        # Load previously processed results
        print(f"Loading previously processed PSORTb results from {psort_results}")
        psort_df = pd.read_csv(psort_results, sep='\t', index_col=0, header=0)
        if 'psortb_localization' in psort_df.columns:
            counts = psort_df['psortb_localization'].value_counts()
            summary = ', '.join(f'{name} {count}' for name, count in counts.items())
            print(f'Localizations for {len(psort_df)} genes (from cache): {summary}')

        return psort_df
        
    
