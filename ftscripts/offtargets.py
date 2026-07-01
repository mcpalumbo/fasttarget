
from ftscripts import programs, metadata, files, structures
from ftscripts.microbiome_catalogues import (
    catalogue_column_prefix,
    catalogue_species_path,
    get_catalogue,
)
import os
import json
import pandas as pd
import multiprocessing
import glob
from tqdm import tqdm
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass

MICROBIOME_BLAST_COLUMNS = 13


@dataclass(frozen=True)
class GenomeSearchResult:
    genome_id: str
    status: str
    output_path: str
    error: str = None


def _format_filter_value(value):
    return f"{float(value):g}"


def _microbiome_result_suffix(identity_filter, coverage_filter):
    identity = _format_filter_value(identity_filter)
    coverage = _format_filter_value(coverage_filter)
    return f"_offtarget_id{identity}_cov{coverage}.tsv"


def _validate_microbiome_blast_output(output_path):
    with open(output_path, "r", encoding="utf-8") as output_file:
        for line_number, line in enumerate(output_file, start=1):
            if len(line.rstrip("\n").split("\t")) != MICROBIOME_BLAST_COLUMNS:
                raise ValueError(
                    f"Invalid DIAMOND output in {output_path} at line {line_number}: "
                    f"expected {MICROBIOME_BLAST_COLUMNS} columns."
                )


def _microbiome_catalogue_status(species_path):
    genome_dirs = {
        entry for entry in os.listdir(species_path)
        if os.path.isdir(os.path.join(species_path, entry))
    }
    indexed_genomes = {
        genome for genome in genome_dirs
        if os.path.isfile(os.path.join(species_path, genome, f"{genome}_DB.dmnd"))
    }
    return genome_dirs, indexed_genomes


def _available_cpus():
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0))
    return multiprocessing.cpu_count()


def search_one_genome(
    genome_id,
    genome_db,
    query_faa,
    output_path,
    identity_filter,
    coverage_filter,
    threads=4,
):
    """
    Runs and validates one DIAMOND search against a representative genome.

    :return: GenomeSearchResult with status success, skipped, or error.
    """

    temporary_output_path = f"{output_path}.tmp"

    try:
        if os.path.exists(output_path):
            try:
                _validate_microbiome_blast_output(output_path)
                return GenomeSearchResult(genome_id, "skipped", output_path)
            except ValueError:
                os.remove(output_path)

        if os.path.exists(temporary_output_path):
            os.remove(temporary_output_path)

        programs.run_diamond_blastp(
            blastdb=genome_db,
            query=query_faa,
            output=temporary_output_path,
            evalue="1e-5",
            outfmt=(
                "6 qseqid sseqid pident length mismatch gapopen qstart qend "
                "sstart send evalue bitscore qcovhsp"
            ),
            cpus=threads,
            identity=identity_filter,
            query_cover=coverage_filter,
            max_target_seqs=1,
        )
        _validate_microbiome_blast_output(temporary_output_path)
        os.replace(temporary_output_path, output_path)
        return GenomeSearchResult(genome_id, "success", output_path)
    except Exception as error:
        if os.path.exists(temporary_output_path):
            os.remove(temporary_output_path)
        return GenomeSearchResult(
            genome_id,
            "error",
            output_path,
            str(error),
        )


def _warn_incomplete_microbiome_catalogue(
    catalogue_name,
    expected_genomes,
    genome_dirs,
    indexed_genomes,
):
    if len(genome_dirs) != expected_genomes:
        logging.warning(
            "The %s catalogue contains %d genome directories; %d were expected.",
            catalogue_name,
            len(genome_dirs),
            expected_genomes,
        )
        print(
            f"Warning: The {catalogue_name} catalogue contains {len(genome_dirs)} "
            f"genome directories; {expected_genomes} were expected."
        )
    if len(indexed_genomes) != expected_genomes:
        logging.warning(
            "Only %d of %d expected %s genomes have a DIAMOND index.",
            len(indexed_genomes),
            expected_genomes,
            catalogue_name,
        )
        print(
            f"Warning: Only {len(indexed_genomes)} of {expected_genomes} expected "
            f"{catalogue_name} genomes have a DIAMOND index."
        )


def human_offtarget_blast (databases_path, output_path, organism_name, cpus=multiprocessing.cpu_count()):

    """
    Runs NCBI BLASTP against the human proteome.
    This function uses the `run_blastp` function from the `programs` module.
    It uses the HUMAN_DB database created by the `index_db_blast_human` function from the `databases` module.
    The blast output is saved in the 'offtarget' folder of the organism directory.

    :param databases_path: Directory where HUMAN databases are stored.
    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param cpus: Number of threads (CPUs) to use in the blast search.
    
    """

    #Database files
    humanprot_index_path = os.path.join(databases_path, 'HUMAN_DB')

    #Organism files
    organism_path = os.path.join(output_path, organism_name)
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

def microbiome_offtarget_blast_species(
    databases_path,
    output_path,
    organism_name,
    catalogue_name,
    identity_filter,
    coverage_filter,
    cpus=multiprocessing.cpu_count(),
    threads_per_genome=4,
):
    """
    Runs Diamond BLASTP of the organism proteome against each genome in the microbiome species catalogue.
    Each representative genome is stored under
    `databases/microbiomes/<catalogue>/species_catalogue`, with its DIAMOND index.

    For each genome, the BLAST output is stored in the 'offtarget' folder of the organism directory,
    with one result file per genome. The process can be resumed if interrupted.

    :param databases_path: Path where the species catalogue databases are stored.
    :param output_path: Path of the organism output.
    :param organism_name: Name of the organism (folder name under 'organism').
    :param catalogue_name: Name of a supported MGnify catalogue.
    :param identity_filter: Identity threshold associated with the result files.
    :param coverage_filter: Query coverage threshold associated with the result files.
    :param cpus: Total CPU budget for concurrent searches.
    :param threads_per_genome: Maximum DIAMOND threads assigned to each genome.
    """

    catalogue = get_catalogue(catalogue_name)
    expected_genomes = catalogue["number_of_species"]
    species_databases_path = catalogue_species_path(databases_path, catalogue_name)

    # Path to organism proteome (.faa file)
    organism_path = os.path.join(output_path, organism_name)
    organism_prot_seq_path = os.path.join(organism_path, "genome", f"{organism_name}.faa")

    # Output folder
    offtarget_path = os.path.join(
        organism_path,
        "offtarget",
        "microbiomes",
        catalogue_name,
        "species_blast_results",
    )
    os.makedirs(offtarget_path, exist_ok=True)

    genome_dirs, indexed_genomes = _microbiome_catalogue_status(species_databases_path)
    _warn_incomplete_microbiome_catalogue(
        catalogue_name,
        expected_genomes,
        genome_dirs,
        indexed_genomes,
    )

    if not indexed_genomes:
        raise RuntimeError(
            f"No indexed microbiome genomes were found in {species_databases_path}."
        )

    result_suffix = _microbiome_result_suffix(identity_filter, coverage_filter)
    available_cpus = _available_cpus()
    cpu_budget = max(1, min(int(cpus), available_cpus))
    genome_threads = max(1, min(int(threads_per_genome), cpu_budget))
    parallel_genomes = max(1, cpu_budget // genome_threads)

    print(
        f"{catalogue_name}: {parallel_genomes} concurrent searches, "
        f"{genome_threads} DIAMOND threads per genome "
        f"({cpu_budget} CPUs available to this module)."
    )

    search_arguments = []
    for genome_dir in sorted(indexed_genomes):
        genome_path = os.path.join(species_databases_path, genome_dir)
        blast_output_path = os.path.join(offtarget_path, f"{genome_dir}{result_suffix}")
        genome_db = os.path.join(genome_path, f'{genome_dir}_DB')
        search_arguments.append((genome_dir, genome_db, blast_output_path))

    results = []
    with ThreadPoolExecutor(max_workers=parallel_genomes) as executor:
        futures = {
            executor.submit(
                search_one_genome,
                genome_id,
                genome_db,
                organism_prot_seq_path,
                blast_output_path,
                identity_filter,
                coverage_filter,
                genome_threads,
            ): genome_id
            for genome_id, genome_db, blast_output_path in search_arguments
        }
        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc=f"Searching {catalogue_name}",
        ):
            results.append(future.result())

    status_counts = {
        status: sum(result.status == status for result in results)
        for status in ("success", "skipped", "error")
    }
    print(
        f"{catalogue_name} search summary: {status_counts['success']} completed, "
        f"{status_counts['skipped']} skipped, {status_counts['error']} failed."
    )

    failed = [result for result in results if result.status == "error"]
    if failed:
        examples = "; ".join(
            f"{result.genome_id}: {result.error}"
            for result in failed[:5]
        )
        raise RuntimeError(
            f"{len(failed)} {catalogue_name} DIAMOND searches failed. {examples}"
        )

    if len(results) != len(indexed_genomes):
        raise RuntimeError(
            f"{catalogue_name} search accounting mismatch: received "
            f"{len(results)} results for {len(indexed_genomes)} indexed genomes."
        )

    return results

def microbiome_offtarget_blast_allproteins (databases_path, output_path, organism_name, cpus=multiprocessing.cpu_count()):

    """
    Runs ncbi blastp against microbiome proteome.
    This function uses the `run_blastp` function from the `programs` module.
    It uses the MICROBIOME_DB database created by the `index_db_blast_microbiome` function from the `databases` module.
    The blast output is saved in the 'offtarget' folder of the organism directory.

    :param databases_path: Directory where MICROBIOME database is stored.
    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param cpus: Number of threads (CPUs) to use in the blast search.
    
    """
    
    #Database files
    species_databases_path = os.path.join(databases_path,  'species_catalogue')
    microbiome_index_path = os.path.join(species_databases_path, 'MICROBIOME_DB')

    #Organism files
    organism_path = os.path.join(output_path, organism_name)
    organism_prot_seq_path = os.path.join(organism_path, f'genome/{organism_name}.faa')

    offtarget_path = os.path.join(organism_path, 'offtarget')
    blast_output_path = os.path.join(offtarget_path, 'microbiome_offtarget_blast.tsv')

    programs.run_diamond_blastp(
        blastdb= microbiome_index_path,
        query= organism_prot_seq_path,
        output=blast_output_path,
        evalue= '1e-5',
        outfmt= '6 std qcovhsp qcovs',
        cpus=cpus
    )

def human_offtarget_parse (output_path, organism_name):

    """
    Parse NCBI BLASTP results against human proteome, stored in the file 'human_offtarget_blast.tsv'.
    Obtains the hit with the highest percentage of identity for each locus_tag.
    Returns a dictionary with locus_tag as key and highest percentage of identity as value, 
    and a DataFrame with all locus_tags from the genome and their respective values.
    The DataFrame is created using the `metadata_table_with_values` function from the `metadata` module.

    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    
    :return: Dictionary with locus_tag as key and highest percentage of identity value.
    :return: DataFrame with all locus_tags from the genome and their respective values.
    """

    offtarget_path = os.path.join(output_path, organism_name, 'offtarget')
    human_blast_output = os.path.join(offtarget_path, 'human_offtarget_blast.tsv')
    human_results = os.path.join(offtarget_path, 'human_offtarget.tsv')

    if not files.file_check(human_results):
        blast_output_df = files.read_blast_output(human_blast_output)

        highest_pident_values = {}

        for index,row in blast_output_df.iterrows():
            qseqid = row['qseqid']
            pident = row['pident']

            if qseqid not in highest_pident_values or pident > highest_pident_values[qseqid]:
                highest_pident_values[qseqid] = pident

        df_human = metadata.metadata_table_with_values(output_path, organism_name, highest_pident_values, 
                                            'human_offtarget', offtarget_path, 'no_hit')
    else:
        print('Human offtarget analysis already done, output file found')
        print(human_results)
        df_human = pd.read_csv(human_results, sep='\t', header=0)

    return df_human

def microbiome_species_parse(
    databases_path,
    output_path,
    organism_name,
    catalogue_name,
    identity_filter,
    coverage_filter,
):
    """
    Parse Diamond BLASTP results against all genomes in the microbiome species catalogue.
    Each genome has its own BLAST output file under 'offtarget' folder of the organism.

    For each protein (qseqid) of the organism, this function determines in which genomes
    it has at least one hit passing the identity and coverage filters.

    :param output_path: Directory of the organism output.
    :param databases_path: Base path where MICROBIOME database is stored.
    :param organism_name: Name of the organism.
    :param catalogue_name: Name of a supported MGnify catalogue.
    :param identity_filter: Minimum percentage identity accepted in the pident column.
    :param coverage_filter: Minimum query coverage accepted in the qcovhsp column.

    Returns:
        - df_microbiome_norm: DataFrame with one row per protein and a column with normalized counts
        - df_microbiome_counts: DataFrame with one row per protein and a column with number of genomes with hits
        - df_microbiome_total_genomes: DataFrame with one row per protein and a column with total number of genomes analyzed
    """

    offtarget_path = os.path.join(
        output_path,
        organism_name,
        "offtarget",
        "microbiomes",
        catalogue_name,
        "species_blast_results",
    )

    catalogue = get_catalogue(catalogue_name)
    expected_genomes = catalogue["number_of_species"]
    species_path = catalogue_species_path(databases_path, catalogue_name)
    genome_dirs, indexed_genomes = _microbiome_catalogue_status(species_path)
    _warn_incomplete_microbiome_catalogue(
        catalogue_name,
        expected_genomes,
        genome_dirs,
        indexed_genomes,
    )

    result_suffix = _microbiome_result_suffix(identity_filter, coverage_filter)
    genome_files = {
        f.removesuffix(result_suffix): f
        for f in os.listdir(offtarget_path)
        if f.endswith(result_suffix)
    }
    searched_genomes = sorted(indexed_genomes.intersection(genome_files))

    if not searched_genomes:
        raise RuntimeError(
            f"No completed microbiome searches were found in {offtarget_path}. "
            "Run microbiome_offtarget_blast_species first."
        )

    if len(searched_genomes) < len(indexed_genomes):
        print(
            f"Warning: Only {len(searched_genomes)} of {len(indexed_genomes)} indexed "
            "microbiome genomes have search results."
        )
    if len(searched_genomes) < expected_genomes:
        print(
            f"Warning: Microbiome scores will be normalized using the "
            f"{len(searched_genomes)} genomes actually analyzed instead of the "
            f"{expected_genomes} genomes expected for {catalogue_name}."
        )

    print("Parsing microbiome species BLAST results...")
    print(f"Catalogue: {catalogue_name}")
    print(f"Expected genomes: {expected_genomes}")
    print(f"Genome directories found: {len(genome_dirs)}")
    print(f"Indexed genomes found: {len(indexed_genomes)}")
    print(f"Genomes analyzed: {len(searched_genomes)}")

    protein_hits = {}
    for genome_name in tqdm(searched_genomes, desc="Parsing microbiome species BLAST results"):
        blast_output_path = os.path.join(offtarget_path, genome_files[genome_name])

        if os.stat(blast_output_path).st_size == 0:
            continue

        df = pd.read_csv(blast_output_path, sep="\t", header=None)
        df.columns = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "qcovhsp"
        ]

        filtered_df = df[
            (df["pident"] >= identity_filter)
            & (df["qcovhsp"] >= coverage_filter)
        ]
        for qseqid in filtered_df["qseqid"].unique():
            protein_hits.setdefault(qseqid, set()).add(genome_name)

    protein_hits = {
        protein: sorted(genomes)
        for protein, genomes in protein_hits.items()
    }
    analyzed_genomes = len(searched_genomes)
    protein_hit_counts = {
        protein: len(genomes) / analyzed_genomes
        for protein, genomes in protein_hits.items()
    }
    protein_hit_totals = {
        protein: len(genomes)
        for protein, genomes in protein_hits.items()
    }
    protein_total_genomes = {
        protein: analyzed_genomes
        for protein in protein_hits
    }

    column_prefix = catalogue_column_prefix(catalogue_name)
    metadata.metadata_table_with_values(
        output_path,
        organism_name,
        protein_hits,
        f'{column_prefix}_offtarget',
        offtarget_path,
        'no_hit',
    )
    df_microbiome_norm = metadata.metadata_table_with_values(
        output_path,
        organism_name,
        protein_hit_counts,
        f'{column_prefix}_offtarget_norm',
        offtarget_path,
        0,
    )
    df_hit_totals = metadata.metadata_table_with_values(
        output_path,
        organism_name,
        protein_hit_totals,
        f'{column_prefix}_offtarget_counts',
        offtarget_path,
        0,
    )
    df_total_genomes = metadata.metadata_table_with_values(
        output_path,
        organism_name,
        protein_total_genomes,
        f'{column_prefix}_genomes_analyzed',
        offtarget_path,
        analyzed_genomes,
    )

    return df_microbiome_norm, df_hit_totals, df_total_genomes


def microbiome_protein_clusters_parse (output_path, organism_name, identity_filter, coverage_filter):

    """
    Parse NCBI BLASTP results against microbiome proteome, stored in the file 'microbiome_offtarget_blast.tsv'.
    Filters results based on identity and coverage thresholds, then counts occurrences of each locus_tag.
    Returns a dictionary with normalized counts of each locus_tag that has at least one hit with UHGP90, 
    and a DataFrame with all locus_tags from the genome and their respective values.
    The DataFrame is created using the `metadata_table_with_values` function from the `metadata` module.

    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param identity_filter: Percentage identity filter value. Keeps results above this value in the pident column.
    :param coverage_filter: Query coverage filter value. Keeps results above this value in the qcovs column.
    
    :return: Dictionary with normalized counts of locus_tags that have at least one hit with UHGP90.
    :return: DataFrame with all locus_tags from the genome and their respective normalized values.
    """

    offtarget_path = os.path.join(output_path, organism_name, 'offtarget')
    microbiome_blast_output = os.path.join(offtarget_path, 'microbiome_offtarget_blast.tsv')
    microbiome_results = os.path.join(offtarget_path, 'gut_microbiome_offtarget.tsv')

    if not files.file_check(microbiome_results):

        blast_output_df = files.read_blast_output(microbiome_blast_output)

        #Filter % identity and coverage
        filtered_df = blast_output_df[(blast_output_df['pident'] > identity_filter) 
                                    & (blast_output_df['qcovs'] > coverage_filter)]

        value_counts = filtered_df['qseqid'].value_counts()
        max_count = value_counts.max()
        norm_counts = value_counts / max_count

        normalized_counts_dict = norm_counts.to_dict()

        df_microbiome = metadata.metadata_table_with_values(output_path, organism_name, normalized_counts_dict, 
                                            'gut_microbiome_offtarget', offtarget_path, 'no_hit')

    else:
        print('Microbiome offtarget analysis already done, output file found')
        print(microbiome_results)
        df_microbiome = pd.read_csv(microbiome_results, sep='\t', header=0)
        
    return df_microbiome


def run_foldseek_human_structures (databases_path, output_path, organism_name, container_engine='docker'):

    """
    Runs Foldseek easy-search against unified human reference structures database.
    Uses ONLY the reference structures for each locus_tag.
    
    Reference structures are obtained via structures.get_all_reference_structures():
    - PDB reference structures: PDB_{uniprot}_{pdb_id}_{chain}.pdb (extracted chains)
    - AlphaFold models: AF_{uniprot}.pdb (full predictions)
    
    Each reference structure is searched against the unified human reference database
    containing both PDB and AlphaFold structures.

    :param output_path: Directory of the organism output.
    :param databases_path: Directory where FOLDSEEK structure databases are stored.
    :param organism_name: Name of the organism.
    :param container_engine: Container engine to use ('docker' or 'singularity').
    :return: Dictionary with locus_tag as key and path to foldseek results file as value.
    """
    foldseek_results_mapping = {}

    print(f'\n{"="*80}')
    print('FOLDSEEK HUMAN OFFTARGET ANALYSIS')
    print(f'{"="*80}\n')

    # Human unified reference database (PDB + AlphaFold structures)
    # DB is created in DB_foldseek subdirectory by programs.run_foldseek_create_index_db
    db_human_path = os.path.join(databases_path, 'human_structures', 'DB_foldseek')

    # Get all reference structures using helper function
    print('Getting reference structures for all locus_tags...')
    reference_dict = structures.get_all_reference_structures(output_path, organism_name, path_mode=True)
    
    # Filter out None values (locus_tags without structures)
    reference_dict = {k: v for k, v in reference_dict.items() if v is not None}
    
    print(f'Found {len(reference_dict)} locus_tags with reference structures')
    
    # Offtarget path - create directory before any early returns
    offtarget_path = os.path.join(output_path, organism_name, 'offtarget')
    foldseek_results_path = os.path.join(offtarget_path, 'foldseek_results')

    if not os.path.exists(foldseek_results_path):
        os.makedirs(foldseek_results_path, exist_ok=True)
    
    if not reference_dict:
        print('ERROR: No reference structures found. Make sure to run structures pipeline first.')
        return {}

    # Run Foldseek for each reference structure
    success_count = 0
    error_count = 0
    skipped_count = 0
    
    for locus_tag, struct_path in reference_dict.items():
        struct_name = os.path.basename(struct_path)
        struct_dir = os.path.dirname(struct_path)
        
        print(f'\n[{locus_tag}] Processing {struct_name}')
        
        # Search against unified human reference database (PDB + AlphaFold)
        try:
            programs.run_foldseek_search(struct_dir, db_human_path, 'DB_human_reference', struct_name, foldseek_results_path, container_engine=container_engine)
            
            # Validate output file exists before reporting success
            struct_basename = struct_name.split('.')[0]
            result_file = os.path.join(foldseek_results_path, f'{struct_basename}_output_foldseek', f'{struct_basename}_vs_DB_human_reference_foldseek_results.tsv')
            
            if os.path.exists(result_file) and os.path.getsize(result_file) > 0:
                print(f'  ✓ Foldseek search completed')
                success_count += 1
                foldseek_results_mapping[locus_tag] = [result_file]
            else:
                print(f'  ✗ Foldseek search failed: output file not found or empty')
                logging.error(f'Foldseek output file missing or empty: {result_file}')
                error_count += 1
                foldseek_results_mapping[locus_tag] = []

        except Exception as e:
            print(f'  ✗ Foldseek search failed')
            logging.exception(f'Error running Foldseek search: {e}')
            error_count += 1
            foldseek_results_mapping[locus_tag] = []
    
    print(f'\n{"="*80}')
    print('FOLDSEEK SUMMARY')
    print(f'{"="*80}')
    print(f'Total locus_tags with structures: {len(reference_dict)}')
    print(f'Successfully processed: {success_count}')
    print(f'Errors: {error_count}')
    print(f'{"="*80}\n')

    return foldseek_results_mapping

def run_foldseek_human_colabfold_structures(databases_path, output_path, organism_name, container_engine='docker'):

    """
    Runs Foldseek easy-search against unified human reference structures database
    using ColabFold structures (CB_*.pdb) for each locus_tag.

    Raw Foldseek results are stored in the shared 'foldseek_results' folder so
    existing CB_* queries can be reused if they were already generated by the
    default analysis.

    :param output_path: Directory of the organism output.
    :param databases_path: Directory where FOLDSEEK structure databases are stored.
    :param organism_name: Name of the organism.
    :param container_engine: Container engine to use ('docker' or 'singularity').
    :return: Dictionary with locus_tag as key and path to foldseek results file as value.
    """
    foldseek_results_mapping = {}

    print(f'\n{"="*80}')
    print('FOLDSEEK HUMAN OFFTARGET ANALYSIS (COLABFOLD)')
    print(f'{"="*80}\n')

    db_human_path = os.path.join(databases_path, 'human_structures', 'DB_foldseek')
    organism_structures_path = os.path.join(output_path, organism_name, 'structures')
    all_locus_tags = metadata.ref_gbk_locus(output_path, organism_name)

    offtarget_path = os.path.join(output_path, organism_name, 'offtarget')
    foldseek_results_path = os.path.join(offtarget_path, 'foldseek_results')

    if not os.path.exists(foldseek_results_path):
        os.makedirs(foldseek_results_path, exist_ok=True)

    success_count = 0
    error_count = 0
    skipped_count = 0

    for locus_tag in all_locus_tags:
        locus_dir = os.path.join(organism_structures_path, locus_tag)
        cb_structures = structures.find_colabfold_for_locus(locus_dir)

        if not cb_structures:
            print(f'\n[{locus_tag}] No ColabFold structure found, skipping')
            skipped_count += 1
            foldseek_results_mapping[locus_tag] = []
            continue

        struct_path = cb_structures[0]
        struct_name = os.path.basename(struct_path)
        struct_dir = os.path.dirname(struct_path)

        print(f'\n[{locus_tag}] Processing {struct_name}')

        try:
            programs.run_foldseek_search(
                struct_dir,
                db_human_path,
                'DB_human_reference',
                struct_name,
                foldseek_results_path,
                container_engine=container_engine
            )

            struct_basename = struct_name.split('.')[0]
            result_file = os.path.join(
                foldseek_results_path,
                f'{struct_basename}_output_foldseek',
                f'{struct_basename}_vs_DB_human_reference_foldseek_results.tsv'
            )

            if os.path.exists(result_file) and os.path.getsize(result_file) > 0:
                print('  ✓ Foldseek search completed')
                success_count += 1
                foldseek_results_mapping[locus_tag] = [result_file]
            else:
                print('  ✗ Foldseek search failed: output file not found or empty')
                logging.error(f'Foldseek output file missing or empty: {result_file}')
                error_count += 1
                foldseek_results_mapping[locus_tag] = []

        except Exception as e:
            print('  ✗ Foldseek search failed')
            logging.exception(f'Error running Foldseek search with ColabFold structure: {e}')
            error_count += 1
            foldseek_results_mapping[locus_tag] = []

    print(f'\n{"="*80}')
    print('FOLDSEEK COLABFOLD SUMMARY')
    print(f'{"="*80}')
    print(f'Total locus_tags in genome: {len(all_locus_tags)}')
    print(f'Successfully processed: {success_count}')
    print(f'Skipped (no ColabFold structure): {skipped_count}')
    print(f'Errors: {error_count}')
    print(f'{"="*80}\n')

    return foldseek_results_mapping

def foldseek_human_parser (output_path, organism_name, map_foldseek):

    """
    Parses Foldseek results for reference structures and maps them to locus_tags.
    
    For each locus_tag, selects the best match (highest TM-score) from the unified 
    human reference database search.

    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param map_foldseek: Dictionary mapping locus_tag to foldseek result file (single file per locus_tag).
    
    :return: Dictionary with locus_tag as key and best foldseek match.
    """

    offtargets_dir = os.path.join(output_path, organism_name, 'offtarget')
    foldseek_results_path = os.path.join(offtargets_dir, 'foldseek_results')

    foldseek_dict_file = os.path.join(foldseek_results_path, 'human_foldseek_dict.json')

    if not files.file_check(foldseek_dict_file):
        
        print('\nParsing Foldseek results...')
        results_foldseek_dict = {}
        
        # Guard against empty mapping (no structures found)
        if not map_foldseek:
            print('No Foldseek results to parse (empty mapping)')
            with open(foldseek_dict_file, 'w') as f:
                json.dump({}, f)
            return {}
        
        # Parse results for each locus_tag
        for locus_tag, result_files in map_foldseek.items():
            
            # Skip if no result files (e.g., foldseek search failed)
            if not result_files:
                results_foldseek_dict[locus_tag] = {
                    'query_structure': None,
                    'target_foldseek': None,
                    'alnlen_foldseek': None,
                    'qcov_foldseek': None,
                    'tcov_foldseek': None,
                    'lddt_foldseek': None,
                    'qtmscore_foldseek': None,
                    'ttmscore_foldseek': None,
                    'alntmscore_foldseek': None,
                    'rmsd_foldseek': None,
                    'prob_foldseek': None,
                    'pident_foldseek': None,
                    'evalue_foldseek': None
                }
                print(f'  {locus_tag}: Foldseek search failed, no results')
                continue
            
            # Process the single result file for this locus_tag
            result_file = result_files[0]
            
            if not files.file_check(result_file):
                print(f'  {locus_tag}: Result file not found: {result_file}')
                results_foldseek_dict[locus_tag] = {
                    'query_structure': None,
                    'target_foldseek': None,
                    'alnlen_foldseek': None,
                    'qcov_foldseek': None,
                    'tcov_foldseek': None,
                    'lddt_foldseek': None,
                    'qtmscore_foldseek': None,
                    'ttmscore_foldseek': None,
                    'alntmscore_foldseek': None,
                    'rmsd_foldseek': None,
                    'prob_foldseek': None,
                    'pident_foldseek': None,
                    'evalue_foldseek': None
                }
                continue
            
            try:
                df = pd.read_csv(result_file, sep='\t', usecols=['query', 'target', 'alnlen', 'qcov', 'tcov', 'lddt', 'qtmscore', 'ttmscore', 'alntmscore', 'rmsd', 'prob', 'pident', 'evalue'])
                
                if not df.empty:
                    # Sort by TM score (max of query and target TM-scores)
                    df['max_tmscore'] = df[['qtmscore', 'ttmscore']].max(axis=1)
                    best_row = df.sort_values(by='max_tmscore', ascending=False).iloc[0]
                    
                    results_foldseek_dict[locus_tag] = {
                        'query_structure': best_row['query'],
                        'target_foldseek': best_row['target'],
                        'alnlen_foldseek': best_row['alnlen'],
                        'qcov_foldseek': best_row['qcov'],
                        'tcov_foldseek': best_row['tcov'],
                        'lddt_foldseek': best_row['lddt'],
                        'qtmscore_foldseek': best_row['qtmscore'],
                        'ttmscore_foldseek': best_row['ttmscore'],
                        'alntmscore_foldseek': best_row['alntmscore'],
                        'rmsd_foldseek': best_row['rmsd'],
                        'prob_foldseek': best_row['prob'],
                        'pident_foldseek': best_row['pident'],
                        'evalue_foldseek': best_row['evalue']
                    }
                    print(f'  {locus_tag}: Best match = {best_row["target"]} (TM-score={best_row["max_tmscore"]:.3f})')
                else:
                    # Empty result file
                    results_foldseek_dict[locus_tag] = {
                        'query_structure': None,
                        'target_foldseek': None,
                        'alnlen_foldseek': None,
                        'qcov_foldseek': None,
                        'tcov_foldseek': None,
                        'lddt_foldseek': None,
                        'qtmscore_foldseek': None,
                        'ttmscore_foldseek': None,
                        'alntmscore_foldseek': None,
                        'rmsd_foldseek': None,
                        'prob_foldseek': None,
                        'pident_foldseek': None,
                        'evalue_foldseek': None
                    }
                    print(f'  {locus_tag}: No hits found in foldseek results')
                    
            except Exception as e:
                logging.exception(f'Could not read {result_file}: {e}')
                results_foldseek_dict[locus_tag] = {
                    'query_structure': None,
                    'target_foldseek': None,
                    'alnlen_foldseek': None,
                    'qcov_foldseek': None,
                    'tcov_foldseek': None,
                    'lddt_foldseek': None,
                    'qtmscore_foldseek': None,
                    'ttmscore_foldseek': None,
                    'alntmscore_foldseek': None,
                    'rmsd_foldseek': None,
                    'prob_foldseek': None,
                    'pident_foldseek': None,
                    'evalue_foldseek': None
                }
        
        # Save results
        files.dict_to_json(foldseek_results_path, 'human_foldseek_dict.json', results_foldseek_dict)
        print(f'\nFoldseek results saved to {foldseek_dict_file}')
        print(f'Total locus_tags with matches: {len([v for v in results_foldseek_dict.values() if v["target_foldseek"] is not None])}')
    else:
        print(f'Loading existing foldseek results from {foldseek_dict_file}')
        results_foldseek_dict = files.json_to_dict(foldseek_dict_file)
    
    return results_foldseek_dict

def foldseek_human_colabfold_parser(output_path, organism_name, map_foldseek):

    """
    Parses Foldseek results for ColabFold structures and maps them to locus_tags.

    For each locus_tag, selects the best match (highest TM-score) from the unified
    human reference database search using the ColabFold structure as query.

    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param map_foldseek: Dictionary mapping locus_tag to foldseek result file.

    :return: Dictionary with locus_tag as key and best foldseek match.
    """

    offtargets_dir = os.path.join(output_path, organism_name, 'offtarget')
    foldseek_results_path = os.path.join(offtargets_dir, 'foldseek_results')

    foldseek_dict_file = os.path.join(foldseek_results_path, 'human_foldseek_colabfold_dict.json')

    if not files.file_check(foldseek_dict_file):

        print('\nParsing Foldseek ColabFold results...')
        results_foldseek_dict = {}

        if not map_foldseek:
            print('No Foldseek ColabFold results to parse (empty mapping)')
            with open(foldseek_dict_file, 'w') as f:
                json.dump({}, f)
            return {}

        for locus_tag, result_files in map_foldseek.items():

            if not result_files:
                results_foldseek_dict[locus_tag] = {
                    'query_structure': None,
                    'target_foldseek': None,
                    'alnlen_foldseek': None,
                    'qcov_foldseek': None,
                    'tcov_foldseek': None,
                    'lddt_foldseek': None,
                    'qtmscore_foldseek': None,
                    'ttmscore_foldseek': None,
                    'alntmscore_foldseek': None,
                    'rmsd_foldseek': None,
                    'prob_foldseek': None,
                    'pident_foldseek': None,
                    'evalue_foldseek': None
                }
                print(f'  {locus_tag}: Foldseek ColabFold search failed, no results')
                continue

            result_file = result_files[0]

            if not files.file_check(result_file):
                print(f'  {locus_tag}: Result file not found: {result_file}')
                results_foldseek_dict[locus_tag] = {
                    'query_structure': None,
                    'target_foldseek': None,
                    'alnlen_foldseek': None,
                    'qcov_foldseek': None,
                    'tcov_foldseek': None,
                    'lddt_foldseek': None,
                    'qtmscore_foldseek': None,
                    'ttmscore_foldseek': None,
                    'alntmscore_foldseek': None,
                    'rmsd_foldseek': None,
                    'prob_foldseek': None,
                    'pident_foldseek': None,
                    'evalue_foldseek': None
                }
                continue

            try:
                df = pd.read_csv(result_file, sep='\t', usecols=['query', 'target', 'alnlen', 'qcov', 'tcov', 'lddt', 'qtmscore', 'ttmscore', 'alntmscore', 'rmsd', 'prob', 'pident', 'evalue'])

                if not df.empty:
                    df['max_tmscore'] = df[['qtmscore', 'ttmscore']].max(axis=1)
                    best_row = df.sort_values(by='max_tmscore', ascending=False).iloc[0]

                    results_foldseek_dict[locus_tag] = {
                        'query_structure': best_row['query'],
                        'target_foldseek': best_row['target'],
                        'alnlen_foldseek': best_row['alnlen'],
                        'qcov_foldseek': best_row['qcov'],
                        'tcov_foldseek': best_row['tcov'],
                        'lddt_foldseek': best_row['lddt'],
                        'qtmscore_foldseek': best_row['qtmscore'],
                        'ttmscore_foldseek': best_row['ttmscore'],
                        'alntmscore_foldseek': best_row['alntmscore'],
                        'rmsd_foldseek': best_row['rmsd'],
                        'prob_foldseek': best_row['prob'],
                        'pident_foldseek': best_row['pident'],
                        'evalue_foldseek': best_row['evalue']
                    }
                    print(f'  {locus_tag}: Best ColabFold match = {best_row["target"]} (TM-score={best_row["max_tmscore"]:.3f})')
                else:
                    results_foldseek_dict[locus_tag] = {
                        'query_structure': None,
                        'target_foldseek': None,
                        'alnlen_foldseek': None,
                        'qcov_foldseek': None,
                        'tcov_foldseek': None,
                        'lddt_foldseek': None,
                        'qtmscore_foldseek': None,
                        'ttmscore_foldseek': None,
                        'alntmscore_foldseek': None,
                        'rmsd_foldseek': None,
                        'prob_foldseek': None,
                        'pident_foldseek': None,
                        'evalue_foldseek': None
                    }
                    print(f'  {locus_tag}: No hits found in Foldseek ColabFold results')

            except Exception as e:
                logging.exception(f'Could not read {result_file}: {e}')
                results_foldseek_dict[locus_tag] = {
                    'query_structure': None,
                    'target_foldseek': None,
                    'alnlen_foldseek': None,
                    'qcov_foldseek': None,
                    'tcov_foldseek': None,
                    'lddt_foldseek': None,
                    'qtmscore_foldseek': None,
                    'ttmscore_foldseek': None,
                    'alntmscore_foldseek': None,
                    'rmsd_foldseek': None,
                    'prob_foldseek': None,
                    'pident_foldseek': None,
                    'evalue_foldseek': None
                }

        files.dict_to_json(foldseek_results_path, 'human_foldseek_colabfold_dict.json', results_foldseek_dict)
        print(f'\nFoldseek ColabFold results saved to {foldseek_dict_file}')
        print(f'Total locus_tags with ColabFold matches: {len([v for v in results_foldseek_dict.values() if v["target_foldseek"] is not None])}')
    else:
        print(f'Loading existing Foldseek ColabFold results from {foldseek_dict_file}')
        results_foldseek_dict = files.json_to_dict(foldseek_dict_file)

    return results_foldseek_dict

def merge_foldseek_data (output_path, organism_name):
    """
    Format foldseek results for final output table.
    
    With the new structure organization, foldseek results are already keyed by locus_tag,
    so this function primarily reformats the data for the final table.
    
    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param id_equivalences: Dictionary with locus_tag and uniprot_id (for compatibility).
    :param uniprot_proteome_annotations: Dictionary with annotations (for compatibility).

    :return: Dictionary with the merged data formatted for final table.
    """

    offtargets_dir = os.path.join(output_path, organism_name, 'offtarget')
    foldseek_results_path = os.path.join(offtargets_dir, 'foldseek_results')
    foldseek_res_file = os.path.join(foldseek_results_path, 'human_foldseek_dict.json')
    foldseek_mapped_file = os.path.join(offtargets_dir, f'{organism_name}_final_foldseek_results.json')

    if not files.file_check(foldseek_mapped_file):
        if files.file_check(foldseek_res_file):  

            results_foldseek_dict = files.json_to_dict(foldseek_res_file)

            mapped_dict = {}

            # Results are already organized by locus_tag, just reformat
            for locus_tag, foldseek_data in results_foldseek_dict.items():
                mapped_dict[locus_tag] = {
                    'gene': locus_tag,
                    'query_structure': foldseek_data.get('query_structure'),
                    'structure': foldseek_data.get('query_structure'),  # The reference structure used
                    'target': foldseek_data.get('target_foldseek'),
                    'alnlen': foldseek_data.get('alnlen_foldseek'),
                    'qcov': foldseek_data.get('qcov_foldseek'),
                    'tcov': foldseek_data.get('tcov_foldseek'),
                    'lddt': foldseek_data.get('lddt_foldseek'),
                    'qtmscore': foldseek_data.get('qtmscore_foldseek'),
                    'ttmscore': foldseek_data.get('ttmscore_foldseek'),
                    'alntmscore': foldseek_data.get('alntmscore_foldseek'),
                    'rmsd': foldseek_data.get('rmsd_foldseek'),
                    'prob': foldseek_data.get('prob_foldseek'),
                    'pident': foldseek_data.get('pident_foldseek'),
                    'evalue': foldseek_data.get('evalue_foldseek')
                }

            files.dict_to_json(offtargets_dir, f'{organism_name}_final_foldseek_results.json', mapped_dict)
            print(f'\nFoldseek data merged and saved to {foldseek_mapped_file}')
            print(f'Total genes with foldseek results: {len([v for v in mapped_dict.values() if v["target"] is not None])}')
        else:
            print(f'File {foldseek_res_file} not found.')
            mapped_dict = {}
    else:
        mapped_dict = files.json_to_dict(foldseek_mapped_file)
        print(f'Foldseek results in {foldseek_mapped_file}.')

    return mapped_dict

def merge_foldseek_colabfold_data(output_path, organism_name):
    """
    Format Foldseek ColabFold results for final output table.

    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.

    :return: Dictionary with the merged data formatted for the final table.
    """

    offtargets_dir = os.path.join(output_path, organism_name, 'offtarget')
    foldseek_results_path = os.path.join(offtargets_dir, 'foldseek_results')
    foldseek_res_file = os.path.join(foldseek_results_path, 'human_foldseek_colabfold_dict.json')
    foldseek_mapped_file = os.path.join(offtargets_dir, f'{organism_name}_final_foldseek_colabfold_results.json')

    if not files.file_check(foldseek_mapped_file):
        if files.file_check(foldseek_res_file):

            results_foldseek_dict = files.json_to_dict(foldseek_res_file)

            mapped_dict = {}

            for locus_tag, foldseek_data in results_foldseek_dict.items():
                mapped_dict[locus_tag] = {
                    'gene': locus_tag,
                    'query_structure': foldseek_data.get('query_structure'),
                    'structure': foldseek_data.get('query_structure'),
                    'target': foldseek_data.get('target_foldseek'),
                    'alnlen': foldseek_data.get('alnlen_foldseek'),
                    'qcov': foldseek_data.get('qcov_foldseek'),
                    'tcov': foldseek_data.get('tcov_foldseek'),
                    'lddt': foldseek_data.get('lddt_foldseek'),
                    'qtmscore': foldseek_data.get('qtmscore_foldseek'),
                    'ttmscore': foldseek_data.get('ttmscore_foldseek'),
                    'alntmscore': foldseek_data.get('alntmscore_foldseek'),
                    'rmsd': foldseek_data.get('rmsd_foldseek'),
                    'prob': foldseek_data.get('prob_foldseek'),
                    'pident': foldseek_data.get('pident_foldseek'),
                    'evalue': foldseek_data.get('evalue_foldseek')
                }

            files.dict_to_json(offtargets_dir, f'{organism_name}_final_foldseek_colabfold_results.json', mapped_dict)
            print(f'\nFoldseek ColabFold data merged and saved to {foldseek_mapped_file}')
            print(f'Total genes with Foldseek ColabFold results: {len([v for v in mapped_dict.values() if v["target"] is not None])}')
        else:
            print(f'File {foldseek_res_file} not found.')
            mapped_dict = {}
    else:
        mapped_dict = files.json_to_dict(foldseek_mapped_file)
        print(f'Foldseek ColabFold results in {foldseek_mapped_file}.')

    return mapped_dict

def final_foldseek_structure_table (output_path, organism_name, mapped_dict):
    """
    Create a final table with the merge dictionary of the Foldseek results. 
    Saves the table in a .tsv file named using the organism name followed by '_final_foldseek_results.tsv'in the 'structures' directory.
    Returns a DataFrame with the final results.

    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param mapped_dict: Dictionary with the merged data.

    :return: DataFrame with the final results.
    """
    
    offtargets_dir = os.path.join(output_path, organism_name, 'offtarget')

    all_locus_tags = metadata.ref_gbk_locus(output_path, organism_name)
    
    # Create a DataFrame with final results
    final_foldseek_file = os.path.join(offtargets_dir, f'{organism_name}_final_foldseek_results.tsv')

    if not files.file_check(final_foldseek_file):

        rows = []

        for locus_tag in all_locus_tags:
            if locus_tag in mapped_dict:
                rows.append(mapped_dict[locus_tag])
            else:
                rows.append({'gene': locus_tag, 'query_structure': None, 'structure': 'No hit', 'target': None, 'alnlen': None, 'qcov': None, 'tcov': None, 'lddt': None, 'qtmscore': None, 'ttmscore': None, 'alntmscore': None, 'rmsd': None, 'prob': None, 'pident': None, 'evalue': None })

        final_foldseek_df = pd.DataFrame(rows).rename(columns={
            'query_structure': 'FS_query_structure',
            'structure': 'FS_organism_structure_query',
            'target': 'FS_human_structure_hit',
            'alnlen': 'FS_alnlen',
            'qcov': 'FS_qcov',
            'tcov': 'FS_tcov',
            'lddt': 'FS_lddt',
            'qtmscore': 'FS_qtmscore',
            'ttmscore': 'FS_ttmscore',
            'alntmscore': 'FS_alntmscore',
            'rmsd': 'FS_rmsd',
            'prob': 'FS_prob',
            'pident': 'FS_pident',
            'evalue': 'FS_evalue'
        })

        final_foldseek_df.to_csv(final_foldseek_file, sep='\t', index=False)
   
        print(f'Foldseek final results saved to {final_foldseek_file}.')
    
    else:
        print(f'Foldseek final results in {final_foldseek_file}.')
        final_foldseek_df = pd.read_csv(final_foldseek_file, sep='\t')

    return final_foldseek_df

def final_foldseek_colabfold_structure_table(output_path, organism_name, mapped_dict):
    """
    Create a final table with the merged dictionary of the Foldseek ColabFold results.
    Saves the table in a .tsv file named using the organism name followed by
    '_final_foldseek_colabfold_results.tsv'.

    :param output_path: Directory of the organism output.
    :param organism_name: Name of the organism.
    :param mapped_dict: Dictionary with the merged data.

    :return: DataFrame with the final results.
    """

    offtargets_dir = os.path.join(output_path, organism_name, 'offtarget')
    all_locus_tags = metadata.ref_gbk_locus(output_path, organism_name)
    final_foldseek_file = os.path.join(offtargets_dir, f'{organism_name}_final_foldseek_colabfold_results.tsv')

    if not files.file_check(final_foldseek_file):

        rows = []

        for locus_tag in all_locus_tags:
            if locus_tag in mapped_dict:
                rows.append(mapped_dict[locus_tag])
            else:
                rows.append({'gene': locus_tag, 'query_structure': None, 'structure': 'No hit', 'target': None, 'alnlen': None, 'qcov': None, 'tcov': None, 'lddt': None, 'qtmscore': None, 'ttmscore': None, 'alntmscore': None, 'rmsd': None, 'prob': None, 'pident': None, 'evalue': None})

        final_foldseek_df = pd.DataFrame(rows).rename(columns={
            'query_structure': 'FS_CB_query_structure',
            'structure': 'FS_CB_organism_structure_query',
            'target': 'FS_CB_human_structure_hit',
            'alnlen': 'FS_CB_alnlen',
            'qcov': 'FS_CB_qcov',
            'tcov': 'FS_CB_tcov',
            'lddt': 'FS_CB_lddt',
            'qtmscore': 'FS_CB_qtmscore',
            'ttmscore': 'FS_CB_ttmscore',
            'alntmscore': 'FS_CB_alntmscore',
            'rmsd': 'FS_CB_rmsd',
            'prob': 'FS_CB_prob',
            'pident': 'FS_CB_pident',
            'evalue': 'FS_CB_evalue'
        })

        final_foldseek_df.to_csv(final_foldseek_file, sep='\t', index=False)

        print(f'Foldseek ColabFold final results saved to {final_foldseek_file}.')

    else:
        print(f'Foldseek ColabFold final results in {final_foldseek_file}.')
        final_foldseek_df = pd.read_csv(final_foldseek_file, sep='\t')

    return final_foldseek_df
