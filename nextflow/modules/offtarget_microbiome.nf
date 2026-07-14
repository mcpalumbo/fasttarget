#!/usr/bin/env nextflow

/*
 * Microbiome off-target scatter-gather workflow.
 *
 * Representative genomes are balanced into shards by FASTA size. Each shard
 * runs sequentially with a fixed CPU allocation, while Nextflow executes
 * multiple shards concurrently.
 */

process PREPARE_MICROBIOME_SHARDS {
    tag "${catalogue_name}"
    label 'low_resources'

    input:
    tuple val(catalogue_name), val(identity_filter), val(coverage_filter)
    val databases_path
    val shard_size

    output:
    tuple val(catalogue_name), val(identity_filter), val(coverage_filter),
        path("shards/*.txt"), emit: shards

    script:
    def base_path = workflow.projectDir.parent
    """
    #!/usr/bin/env python3
    import os
    import sys

    sys.path.insert(0, '${base_path}')

    from ftscripts import offtargets
    from ftscripts.microbiome_catalogues import catalogue_species_path

    species_path = catalogue_species_path('${databases_path}', '${catalogue_name}')
    shard_paths = offtargets.create_microbiome_shards(
        species_path,
        os.path.join(os.getcwd(), 'shards'),
        ${shard_size},
    )
    print(
        f"Prepared {len(shard_paths)} shards for ${catalogue_name} "
        f"with a target size of ${shard_size} genomes."
    )
    """

    stub:
    """
    mkdir -p shards
    echo "MGYG000000001" > shards/shard_0001.txt
    """
}


process MICROBIOME_SHARD_SEARCH {
    tag "${catalogue_name}:${shard_file.simpleName}"
    label 'microbiome_shard'
    cache false
    cpus { threads_per_genome as int }
    maxForks params.microbiome_max_forks
    memory params.microbiome_search_memory
    time params.microbiome_search_time
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(catalogue_name), val(identity_filter), val(coverage_filter),
        path(shard_file)
    path query_faa
    val organism_name
    val output_path
    val databases_path
    val threads_per_genome

    output:
    tuple val(catalogue_name), val(identity_filter), val(coverage_filter),
        path("completed/${catalogue_name}_${shard_file.simpleName}.done"),
        emit: completed

    script:
    def base_path = workflow.projectDir.parent
    """
    #!/usr/bin/env python3
    import os
    import sys

    sys.path.insert(0, '${base_path}')

    from ftscripts import offtargets
    from ftscripts.microbiome_catalogues import catalogue_species_path

    catalogue_name = '${catalogue_name}'
    species_path = catalogue_species_path('${databases_path}', catalogue_name)
    results_path = os.path.join(
        '${output_path}',
        '${organism_name}',
        'offtarget',
        'microbiomes',
        catalogue_name,
        'species_blast_results',
    )
    os.makedirs(results_path, exist_ok=True)
    offtargets.validate_microbiome_search_environment(
        '${query_faa}',
        results_path,
    )

    suffix = offtargets.microbiome_result_suffix(
        ${identity_filter},
        ${coverage_filter},
    )
    with open('${shard_file}', 'r', encoding='utf-8') as shard_handle:
        genome_ids = [line.strip() for line in shard_handle if line.strip()]

    results = []
    consolidated = offtargets.is_microbiome_consolidation_compatible(
        '${output_path}',
        '${organism_name}',
        catalogue_name,
        ${identity_filter},
        ${coverage_filter},
        verify_checksum=False,
    )
    if not consolidated:
        for genome_id in genome_ids:
            genome_db = os.path.join(species_path, genome_id, f'{genome_id}_DB')
            result = offtargets.search_one_genome(
                genome_id,
                genome_db,
                '${query_faa}',
                os.path.join(results_path, f'{genome_id}{suffix}'),
                ${identity_filter},
                ${coverage_filter},
                ${threads_per_genome},
            )
            results.append(result)
            if result.status == 'system_error':
                raise RuntimeError(
                    f'Systemic DIAMOND search failure for {genome_id}: {result.error}'
                )

    failed = [result for result in results if result.status == 'error']
    if failed:
        examples = '; '.join(
            f'{result.genome_id}: {result.error}' for result in failed[:5]
        )
        raise RuntimeError(
            f'{len(failed)} searches failed in ${shard_file.simpleName}: {examples}'
        )

    os.makedirs('completed', exist_ok=True)
    marker = os.path.join(
        'completed',
        '${catalogue_name}_${shard_file.simpleName}.done',
    )
    with open(marker, 'w', encoding='utf-8') as marker_file:
        marker_file.write(
            f'completed={len(results)}\\n'
            f'skipped={sum(r.status == "skipped" for r in results)}\\n'
        )
    """

    stub:
    """
    mkdir -p completed
    echo "completed=1" > completed/${catalogue_name}_${shard_file.simpleName}.done
    """
}


process PARSE_MICROBIOME_RESULTS {
    tag "${catalogue_name}"
    label 'microbiome_parse'
    memory params.microbiome_parse_memory
    time params.microbiome_parse_time
    publishDir "${output_path}", mode: 'copy',
        pattern: "${organism_name}/offtarget/microbiomes/${catalogue_name}/**"

    input:
    tuple val(catalogue_name), val(identity_filter), val(coverage_filter),
        path(completion_markers)
    val organism_name
    val output_path
    val databases_path
    path genome_gbk

    output:
    path "${organism_name}/offtarget/microbiomes/${catalogue_name}/species_blast_results/*_offtarget_norm.tsv",
        emit: normalized_table
    path "${organism_name}/offtarget/microbiomes/${catalogue_name}/species_blast_results/*_offtarget_counts.tsv",
        emit: counts_table
    path "${organism_name}/offtarget/microbiomes/${catalogue_name}/species_blast_results/*_genomes_analyzed.tsv",
        emit: genomes_analyzed_table
    path "${organism_name}/offtarget/microbiomes/${catalogue_name}/species_blast_results/${catalogue_name}_offtarget_hits.parquet",
        emit: consolidated_hits
    path "${organism_name}/offtarget/microbiomes/${catalogue_name}/species_blast_results/${catalogue_name}_offtarget_manifest.json",
        emit: consolidated_manifest
    val organism_name, emit: organism_name

    script:
    def base_path = workflow.projectDir.parent
    """
    #!/usr/bin/env python3
    import os
    import glob
    import shutil
    import sys

    sys.path.insert(0, '${base_path}')

    from ftscripts import offtargets
    from ftscripts.microbiome_catalogues import catalogue_column_prefix

    work_dir = os.getcwd()
    genome_dir = os.path.join(work_dir, '${organism_name}', 'genome')
    os.makedirs(genome_dir, exist_ok=True)
    shutil.copy2(
        '${genome_gbk}',
        os.path.join(genome_dir, '${organism_name}.gbk'),
    )

    offtargets.consolidate_microbiome_hits(
        '${databases_path}',
        '${output_path}',
        '${organism_name}',
        '${catalogue_name}',
        ${identity_filter},
        ${coverage_filter},
    )

    tables = offtargets.microbiome_species_parse(
        '${databases_path}',
        '${output_path}',
        '${organism_name}',
        '${catalogue_name}',
        ${identity_filter},
        ${coverage_filter},
        genome_output_path=work_dir,
    )

    prefix = catalogue_column_prefix('${catalogue_name}')
    persistent_results = os.path.join(
        '${output_path}',
        '${organism_name}',
        'offtarget',
        'microbiomes',
        '${catalogue_name}',
        'species_blast_results',
    )
    local_results = os.path.join(
        '${organism_name}',
        'offtarget',
        'microbiomes',
        '${catalogue_name}',
        'species_blast_results',
    )
    os.makedirs(local_results, exist_ok=True)

    result_tables = []
    for pattern in (
        f'{prefix}*_offtarget_norm.tsv',
        f'{prefix}*_offtarget_counts.tsv',
        f'{prefix}_genomes_analyzed.tsv',
    ):
        result_tables.extend(glob.glob(os.path.join(persistent_results, pattern)))

    for source in sorted(set(result_tables)):
        shutil.copy2(source, os.path.join(local_results, os.path.basename(source)))

    for consolidated_name in (
        '${catalogue_name}_offtarget_hits.parquet',
        '${catalogue_name}_offtarget_manifest.json',
    ):
        source = os.path.join(persistent_results, consolidated_name)
        shutil.copy2(source, os.path.join(local_results, consolidated_name))

    print('Parsed ${catalogue_name} microbiome results.')
    """

    stub:
    """
    results=${organism_name}/offtarget/microbiomes/${catalogue_name}/species_blast_results
    mkdir -p "\$results"
    echo -e "gene\tstub_offtarget_norm\ngene1\t0.15" > "\$results/stub_offtarget_norm.tsv"
    echo -e "gene\tstub_offtarget_counts\ngene1\t5" > "\$results/stub_offtarget_counts.tsv"
    echo -e "gene\tstub_genomes_analyzed\ngene1\t1" > "\$results/stub_genomes_analyzed.tsv"
    touch "\$results/${catalogue_name}_offtarget_hits.parquet"
    echo '{}' > "\$results/${catalogue_name}_offtarget_manifest.json"
    """
}
