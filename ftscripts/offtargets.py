from ftscripts import programs, metadata
import os
import pandas as pd
import multiprocessing

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

def essential_deg_blast(base_path, organism_name, cpus=multiprocessing.cpu_count()):

    """
    Runs NCBI BLASTP against DEG bacterial proteome.
    This function uses the `run_blastp` function from the `programs` module.
    It uses the DEG_DB database created by the `index_db_blast_deg` function from the `databases` module.
    The blast output is saved in the 'offtarget' folder of the organism directory.

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.
    :param cpus: Number of threads (CPUs) to use in blast search.
    
    """

    #Database files
    databases_path = os.path.join(base_path, 'databases')
    deg_index_path = os.path.join(databases_path, 'DEG_DB')

    #Organism files
    organism_path = os.path.join(base_path, f'organism/{organism_name}')
    organism_prot_seq_path = os.path.join(organism_path, f'genome/{organism_name}.faa')

    offtarget_path = os.path.join(organism_path, 'offtarget')
    blast_output_path = os.path.join(offtarget_path, 'deg_blast.tsv')

    programs.run_blastp(
        blastdb= deg_index_path,
        query= organism_prot_seq_path,
        output=blast_output_path,
        evalue= '1e-5',
        outfmt= '6 std qcovhsp qcovs',
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

    blast_output_df = read_blast_output(human_blast_output)

    highest_pident_values = {}

    for index,row in blast_output_df.iterrows():
        qseqid = row['qseqid']
        pident = row['pident']

        if qseqid not in highest_pident_values or pident > highest_pident_values[qseqid]:
             highest_pident_values[qseqid] = pident

    df_human = metadata.metadata_table_with_values(base_path, organism_name, highest_pident_values, 
                                        'human_offtarget', offtarget_path, 'no_hit')

    return highest_pident_values, df_human

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

    return normalized_counts_dict, df_microbiome

def deg_parse (base_path, organism_name, identity_filter, coverage_filter):

    """
    Parse NCBI BLASTP results against DEG bacterial proteome, stored in the file 'deg_blast.tsv'.
    Filters results based on identity and coverage thresholds, and obtains the proteins that has at least one hit with DEG.
    Returns a list with locus_tags that have at least one hit with DEG, and a DataFrame with all locus_tags from the genome and a TRUE/FALSE value.
    The DataFrame is created using the `metadata_table_bool` function from the `metadata` module.
    
    :param base_path:  Base path where the repository data is stored.
    :param organism_name: Name of organism.
    :param identity_filter: Value of % of Identity filter. Keeps results above this value in column pident.
    :param coverage_filter: Value of query coverage filter. Keeps results above this value in column qcovs.
    
    :return: List with locus_tags that have at least one hit in DEG.
    :return: DataFrame with all locus_tags from genome and a TRUE/FALSE value.
    """

    offtarget_path = os.path.join(base_path, 'organism', f'{organism_name}', 'offtarget')
    deg_blast_output = os.path.join(offtarget_path, 'deg_blast.tsv')

    blast_output_df = read_blast_output(deg_blast_output)

    #Filter % identity and coverage
    filtered_df = blast_output_df[(blast_output_df['pident'] > identity_filter) 
                                  & (blast_output_df['qcovs'] > coverage_filter)]

    unique_values = filtered_df['qseqid'].unique()
    deg_hits = unique_values.tolist()

    df_deg = metadata.metadata_table_bool(base_path, organism_name, deg_hits, 
                                        'hit_in_deg', offtarget_path)

    return deg_hits, df_deg

    
