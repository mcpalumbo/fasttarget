"""Isolated dRep selection; input datasets are read-only sources."""
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
from importlib.metadata import version
import json
import multiprocessing
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

from Bio import SeqIO
import pandas as pd


def digest(path):
    with open(path, 'rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


# ---------------------------------------------------------------------------
# Persistent per-genome cache (conservation/genome_cache/)
#
# Reading 1000+ GenBank files is the slow part of the selection and its result
# only depends on the GenBank bytes, so it is computed once per gbff and reused
# by every later run (e.g. when only the dRep thresholds change):
#   genomes.tsv    one row per genome: gbk/fna sha256, sequence digest, hosts and
#                  translated-CDS count (the protein QC input). Human-readable.
#   fna/<acc>.fna  FASTA written from the GenBank (dRep input).
#   prepared/<acc>-<sha12>/  gff/faa/... from genome.ref_genome_files.
# A row is reused only if the gbff sha256 matches and the cached FASTA still
# matches its recorded sha256; anything else is recomputed. Files whose size and
# mtime are unchanged are not re-hashed (the disk read dominated cached runs).
# ---------------------------------------------------------------------------
CACHE_SCHEMA = 1
CACHE_COLUMNS = ['accession', 'gbk_sha256', 'fna_sha256', 'sequence_digest',
                 'hosts', 'protein_count', 'n_records', 'schema',
                 'gbk_stat', 'fna_stat']


def _stat(path):
    """Quick fingerprint (size:mtime_ns). Equal fingerprint -> the recorded sha256 is
    trusted without re-reading the file (make/rsync-style check); otherwise it is hashed."""
    st = os.stat(path)
    return f'{st.st_size}:{st.st_mtime_ns}'


def _verified_sha(path, cached_sha, cached_stat):
    if cached_sha and cached_stat == _stat(path):
        return cached_sha
    return digest(path)
CACHE_FLUSH_EVERY = 50


def _translated_cds_count(records):
    return sum(feature.type == 'CDS' and bool(feature.qualifiers.get('translation', [''])[0])
               for record in records for feature in record.features)


def _sequence_digest(records):
    sequences = sorted(hashlib.sha256(str(r.seq).upper().encode()).hexdigest() for r in records)
    return hashlib.sha256('\n'.join(sequences).encode()).hexdigest()


def _read_genbank(path):
    records = list(SeqIO.parse(path, 'genbank'))
    if not records or any(not len(record.seq) for record in records):
        raise ValueError(f'Empty genome: {path}')
    return records


def _scan_genome(job):
    """Worker: one GenBank -> cache row (+ FASTA). Reuses a verified cached row."""
    accession, gbk, fna, cached = job
    cached = cached or {}
    gbk_sha = _verified_sha(gbk, cached.get('gbk_sha256'), cached.get('gbk_stat'))
    if (cached.get('gbk_sha256') == gbk_sha and str(cached.get('schema')) == str(CACHE_SCHEMA)
            and os.path.isfile(fna)
            and _verified_sha(fna, cached.get('fna_sha256'), cached.get('fna_stat')) == cached.get('fna_sha256')):
        row = dict(cached)
        row.update(gbk_stat=_stat(gbk), fna_stat=_stat(fna))
        return row, True
    records = _read_genbank(gbk)
    hosts = sorted({h for r in records for f in r.features if f.type == 'source'
                    for h in f.qualifiers.get('host', [])})
    tmp = f'{fna}.tmp-{os.getpid()}'
    SeqIO.write(records, tmp, 'fasta')
    os.replace(tmp, fna)
    return {'accession': accession, 'gbk_sha256': gbk_sha, 'fna_sha256': digest(fna),
            'sequence_digest': _sequence_digest(records), 'hosts': json.dumps(hosts),
            'protein_count': _translated_cds_count(records), 'n_records': len(records),
            'schema': CACHE_SCHEMA, 'gbk_stat': _stat(gbk), 'fna_stat': _stat(fna)}, False


def _load_cache(table):
    if not table.is_file():
        return {}
    with table.open(newline='') as handle:
        return {row['accession']: row for row in csv.DictReader(handle, delimiter='\t')}


def _write_cache(table, rows):
    tmp = table.with_name(table.name + '.tmp')
    with tmp.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=CACHE_COLUMNS, delimiter='\t')
        writer.writeheader()
        for accession in sorted(rows):
            writer.writerow({k: rows[accession].get(k, '') for k in CACHE_COLUMNS})
    os.replace(tmp, table)


def scan_genomes(folders, cache, cpus):
    """Return {accession: (cache_row, gbk_path)} for every dataset folder, in parallel.

    The cache table is rewritten every CACHE_FLUSH_EVERY genomes, so an interrupted
    run keeps the work already done.
    """
    cache = Path(cache)
    (cache / 'fna').mkdir(parents=True, exist_ok=True)
    table = cache / 'genomes.tsv'
    rows = _load_cache(table)
    jobs, gbks = [], {}
    for folder in folders:
        paths = sorted(folder.glob('*.gbff')) or sorted(folder.glob('*.gbk'))
        if len(paths) != 1:
            raise ValueError(f'Expected one GenBank in {folder}, found {len(paths)}')
        gbks[folder.name] = str(paths[0].resolve())
        jobs.append((folder.name, gbks[folder.name], str(cache / 'fna' / f'{folder.name}.fna'),
                     rows.get(folder.name)))
    reused = done = 0
    with multiprocessing.get_context('spawn').Pool(max(1, cpus)) as pool:
        for row, hit in pool.imap_unordered(_scan_genome, jobs, chunksize=4):
            rows[row['accession']] = row
            reused += hit
            done += 1
            if not hit and done % CACHE_FLUSH_EVERY == 0:
                _write_cache(table, rows)
    _write_cache(table, rows)
    print(f'Genome cache: {reused} reused, {len(jobs) - reused} scanned ({table})')
    return {name: (rows[name], gbks[name]) for name in gbks}


def _link_or_copy(source, target):
    try:
        os.link(source, target)
    except OSError:
        shutil.copyfile(source, target)


# ---------------------------------------------------------------------------
# dRep result cache (genome_cache/drep/<key>.json)
#
# dRep is deterministic given the input FASTA bytes, its clustering parameters and
# the tool versions, so its outcome (secondary cluster + centrality per genome) is
# stored under a key derived from exactly those. A later run with the same genomes
# and parameters (e.g. after Roary failed) skips dRep; changing any of them -- a
# threshold, a genome, a dRep/fastANI/mash version -- gives a new key and reruns it.
# Paths and the thread count are excluded from the key: they do not change the result.
# ---------------------------------------------------------------------------
def _drep_parameters(command):
    params, skip = [], False
    for i, arg in enumerate(command[3:], start=3):   # drop 'dRep dereplicate <workdir>'
        if skip:
            skip = False
            continue
        if arg in ('-g', '-p'):                      # genome list path / threads
            skip = True
            continue
        params.append(arg)
    return params


def drep_cache_key(candidates, command, versions):
    payload = {'genomes': {name: info['fna_sha256'] for name, info in sorted(candidates.items())},
               'parameters': _drep_parameters(command), 'versions': versions}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def load_drep_result(cache, candidates, command, versions):
    key = drep_cache_key(candidates, command, versions)
    path = Path(cache) / 'drep' / f'{key}.json'
    if not path.is_file():
        return None
    try:
        stored = json.loads(path.read_text())
        clusters, scores = stored['clusters'], {k: float(v) for k, v in stored['centrality'].items()}
    except (OSError, ValueError, KeyError):
        return None
    if (stored.get('key') != key or set(clusters) != set(candidates) or set(scores) != set(candidates)
            or any(pd.isna(v) or not 0 <= v <= 1 for v in scores.values())):
        return None
    return clusters, scores, key


def save_drep_result(cache, candidates, command, versions, clusters, scores):
    key = drep_cache_key(candidates, command, versions)
    folder = Path(cache) / 'drep'
    folder.mkdir(parents=True, exist_ok=True)
    stored = {'key': key, 'parameters': _drep_parameters(command), 'versions': versions,
              'genomes': {name: info['fna_sha256'] for name, info in sorted(candidates.items())},
              'clusters': {k: str(v) for k, v in clusters.items()},
              'centrality': {k: float(v) for k, v in scores.items()}}
    tmp = folder / f'{key}.json.tmp'
    tmp.write_text(json.dumps(stored, indent=1, sort_keys=True) + '\n')
    os.replace(tmp, folder / f'{key}.json')
    return key


def options(core):
    value = core.get('dereplication', {})
    if not isinstance(value, dict):
        raise ValueError('core.dereplication must be a mapping')
    enabled = value.get('enabled', False)
    threshold = value.get('min_genomes', 300)
    if not isinstance(enabled, bool):
        raise ValueError('core.dereplication.enabled must be boolean')
    if isinstance(threshold, bool) or not isinstance(threshold, int) or threshold < 0:
        raise ValueError('core.dereplication.min_genomes must be a nonnegative integer')
    quality = core.get('protein_qc', {})
    if not isinstance(quality, dict):
        raise ValueError('core.protein_qc must be a mapping')
    warn, exclude = quality.get('warn_ratio', 0.8), quality.get('exclude_ratio', 0.5)
    if (any(isinstance(x, bool) or not isinstance(x, (int, float)) for x in (warn, exclude))
            or not 0 <= exclude <= warn <= 1):
        raise ValueError('protein_qc requires 0 <= exclude_ratio <= warn_ratio <= 1')
    return enabled, threshold


def select(dataset, destination, min_genomes=300, cpus=4, reference=None, human_only=False, quality=None,
           enabled=True, cache=None):
    """Run above threshold, counting downloaded assemblies before host selection.

    Returns a manifest. Destination must be new; never resume partial dRep output.
    Native dRep/FastANI coverage uses aligned query fragments (larger method).
    Per-genome parsing (FASTA, hosts, translated-CDS counts) comes from the persistent
    genome cache (default: <dataset>/../../../genome_cache) and runs in parallel.
    """
    if isinstance(min_genomes, bool) or not isinstance(min_genomes, int) or min_genomes < 0:
        raise ValueError("min_genomes must be a nonnegative integer")
    if isinstance(cpus, bool) or not isinstance(cpus, int) or cpus < 1:
        raise ValueError("cpus must be a positive integer")
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    inputs = destination / 'inputs'
    inputs.mkdir()
    candidates, excluded = {}, {}
    cache = Path(cache) if cache else Path(dataset).resolve().parents[2] / 'genome_cache'

    ref_digest = _sequence_digest(_read_genbank(reference)) if reference else None
    folders = sorted(p for p in Path(dataset).iterdir() if p.is_dir() and p.name.startswith(('GCF_', 'GCA_')))
    scanned = scan_genomes(folders, cache, cpus)
    counts = {}
    for folder in folders:
        row, gbk = scanned[folder.name]
        if row['sequence_digest'] == ref_digest:
            excluded[folder.name] = 'identical_reference_sequences'
            continue
        if human_only and 'Homo sapiens' not in json.loads(row['hosts']):
            excluded[folder.name] = 'host_filter'
            continue
        fasta = inputs / f'{folder.name}.fna'
        _link_or_copy(cache / 'fna' / f'{folder.name}.fna', fasta)
        candidates[fasta.name] = {'accession': folder.name, 'gbk': gbk,
                                  'gbk_sha256': row['gbk_sha256'], 'fna_sha256': row['fna_sha256']}
        counts[folder.name] = int(row['protein_count'])
    from ftscripts import genome
    quality = quality or {}
    qc_excluded = genome.conservation_protein_qc(
        {info['accession']: info['gbk'] for info in candidates.values()},
        str(destination / 'protein_qc.tsv'), reference,
        quality.get('warn_ratio', 0.8), quality.get('exclude_ratio', 0.5), counts=counts)
    for name in list(candidates):
        if candidates[name]['accession'] in qc_excluded:
            excluded[candidates[name]['accession']] = 'protein_count_below_median_cutoff'
            del candidates[name]
    genome.conservation_annotation_report(str(dataset), str(destination / 'annotation_provenance.tsv'))
    active = enabled and len(folders) > min_genomes and len(candidates) > 1
    clusters = {name: name for name in candidates}
    scores = {name: 0.0 for name in candidates}
    command = None
    versions = {}
    drep_key, drep_reused = None, False
    if active:
        versions['dRep'] = version('drep')
        if versions['dRep'] != '3.6.2':
            raise RuntimeError('This integration requires dRep 3.6.2')
        if not shutil.which('dRep'):
            raise RuntimeError('dRep not installed in the active environment')
        for exe, args in [('fastANI', ['--version']), ('mash', ['--version'])]:
            location = shutil.which(exe)
            if not location:
                raise RuntimeError(f'{exe} not installed in the active environment')
            result = subprocess.run([location, *args], capture_output=True, text=True, check=True)
            versions[exe] = (result.stdout + result.stderr).strip()
        listing = destination / 'genomes.txt'
        listing.write_text(''.join(str(inputs / name) + '\n' for name in sorted(candidates)))
        command = ['dRep', 'dereplicate', str(destination / 'drep'), '-g', str(listing),
                   # -pa 0.998: MASH primary clusters at 2x the secondary ANI distance, so
                   # fastANI only compares near-identical genomes instead of all-vs-all
                   # (a single-species cohort otherwise collapses into one primary cluster).
                   '-p', str(cpus), '--S_algorithm', 'fastANI', '-pa', '0.998', '-sa', '0.999', '-nc', '0.90',
                   '-cm', 'larger', '--clusterAlg', 'average', '--ignoreGenomeQuality', '-l', '0',
                   '-comW', '0', '-conW', '0', '-strW', '0', '-N50W', '0', '-sizeW', '0', '-centW', '1']
        cached = load_drep_result(cache, candidates, command, versions)
        if cached is not None:
            clusters, scores, drep_key = cached
            drep_reused = True
            print(f'dRep result reused from cache ({drep_key}): same genomes, parameters and versions')
        else:
            with (destination / 'drep.log').open('w') as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
            tables = destination / 'drep' / 'data_tables'
            cdb = pd.read_csv(tables / 'Cdb.csv')
            gdb = pd.read_csv(tables / 'genomeInformation.csv')
            if cdb.genome.duplicated().any() or set(cdb.genome) != set(candidates):
                raise ValueError('dRep cluster membership differs from requested genomes')
            clusters = dict(zip(cdb.genome, cdb.secondary_cluster.astype(str)))
            if gdb.genome.duplicated().any() or set(gdb.genome) != set(candidates):
                raise ValueError('Incomplete dRep centrality table')
            scores = dict(zip(gdb.genome, gdb.centrality))
            if any(pd.isna(v) or not 0 <= v <= 1 for v in scores.values()):
                raise ValueError('Invalid dRep centrality values')
            drep_key = save_drep_result(cache, candidates, command, versions, clusters, scores)
    representatives = []
    for cluster in sorted(set(clusters.values())):
        members = [name for name in candidates if clusters[name] == cluster]
        representatives.append(min(members, key=lambda name: (-scores[name], name)))
    manifest = {'downloaded_count': len(folders), 'eligible_count': len(candidates),
                'min_genomes': min_genomes, 'applied': active, 'versions': versions,
                'command': command, 'drep_cache_key': drep_key, 'drep_reused': drep_reused,
                'excluded': excluded, 'candidates': candidates,
                'policy': {'ani': 0.999, 'coverage': 0.90, 'linkage': 'average',
                           'backend': 'fastANI', 'primary_ani': 0.998,
                           'coverage_semantics': 'dRep native directional fragment fraction > 0.90 before reciprocal ANI averaging',
                           'representative': 'maximum dRep within-cluster mean ANI; accession lexical tie-break'},
                'clusters': clusters, 'centrality': scores, 'representatives': representatives,
                'qc_status': 'host, readable sequence and relative translated-CDS count; contamination filtering out of scope'}
    (destination / 'selection.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def _sync_inputs(folder, sources, link):
    """Make folder hold exactly {name: source}, as hardlinks (link=True) or copies.

    Files already identical to their source (same inode) are left alone, files of
    genomes no longer in the cohort are removed, so the folder is the cohort.
    """
    folder.mkdir(parents=True, exist_ok=True)
    for stale in folder.iterdir():
        if stale.name not in sources:
            stale.unlink()
    for name, source in sources.items():
        target = folder / name
        if link and target.exists() and os.path.samefile(target, source):
            continue
        tmp = folder / f'.{name}.tmp'
        if tmp.exists():
            tmp.unlink()
        if link:
            _link_or_copy(source, tmp)
        else:
            shutil.copyfile(source, tmp)
        os.replace(tmp, target)


def prepare(output_path, organism, core, cpus, human_only, engine):
    """Select the cohort and stage the Roary/CoreCruncher inputs in <organism>/conservation.

    Everything lives in the organism's own conservation folder (selection/,
    roary_output/, corecruncher_output/, cohort.json), so results are where a user
    looks for them and a finished Roary or CoreCruncher run is reused by the next
    run when the cohort and parameters are the same. A changed cohort changes the
    content fingerprint of the staged inputs, which core_genome_programs refuses to
    mix with an old clustering. Returns output_path (unchanged).
    """
    from ftscripts import genome
    enabled, threshold = options(core)
    base = Path(output_path).resolve() / organism
    conservation = base / 'conservation'
    reference = base / 'genome' / f'{organism}.gbk'
    cache = conservation / 'genome_cache'
    # select() needs a fresh directory (dRep must never resume partial output); it is
    # built aside and swapped in, so a failed selection keeps the previous one.
    staging = Path(tempfile.mkdtemp(prefix='.selection-', dir=conservation))
    try:
        selection = select(conservation / f'{organism}_dataset' / 'ncbi_dataset' / 'data',
                           staging / 'selection', threshold, cpus, reference, human_only,
                           core.get('protein_qc'), enabled, cache=cache)
        final_selection = conservation / 'selection'
        if final_selection.exists():
            shutil.rmtree(final_selection)
        os.replace(staging / 'selection', final_selection)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    shutil.copyfile(final_selection / 'protein_qc.tsv', conservation / 'protein_qc.tsv')

    gffs = conservation / 'roary_output' / 'gff'
    faas = conservation / 'corecruncher_output' / 'faa'
    cohort = [organism]
    infos = [selection['candidates'][name] for name in selection['representatives']]
    for info in infos:
        if digest(info['gbk']) != info['gbk_sha256']:
            raise RuntimeError(f"Comparator changed after selection: {info['accession']}")
        if info['accession'] == organism:
            raise ValueError('Comparator accession collides with reference name')

    def prepared_dir(info):
        # Persistent and keyed by the GenBank hash: a changed gbff gets a new directory,
        # so ref_genome_files never sees artifacts from other bytes.
        accession = info['accession']
        prepared = cache / 'prepared' / f"{accession}-{info['gbk_sha256'][:12]}"
        try:
            genome.ref_genome_files(info['gbk'], str(prepared), accession, container_engine=engine)
        except RuntimeError:
            # Interrupted publication left unverifiable files in our own cache: rebuild.
            shutil.rmtree(prepared, ignore_errors=True)
            genome.ref_genome_files(info['gbk'], str(prepared), accession, container_engine=engine)
        return accession, prepared

    with ThreadPoolExecutor(max_workers=max(1, cpus)) as executor:
        prepared_dirs = dict(executor.map(prepared_dir, infos))
    gff_sources = {f'{organism}.gff': base / 'genome' / f'{organism}.gff'}
    faa_sources = {f'{organism}.faa': base / 'genome' / f'{organism}.faa'}
    for info in infos:
        accession = info['accession']
        gff_sources[f'{accession}.gff'] = prepared_dirs[accession] / f'{accession}.gff'
        faa_sources[f'{accession}.faa'] = prepared_dirs[accession] / f'{accession}.faa'
        cohort.append(accession)
    # Roary only reads its GFFs (it writes fixed copies under results/), so they are
    # hardlinked from the verified caches; protein FASTAs are small and copied.
    _sync_inputs(gffs, gff_sources, link=True)
    _sync_inputs(faas, faa_sources, link=False)
    for folder, suffix in [(gffs, '.gff'), (faas, '.faa')]:
        if {p.stem for p in folder.glob('*' + suffix)} != set(cohort):
            raise ValueError('Prepared tools have different cohorts')
    selection['reference'] = {'name': organism, 'gbk_sha256': digest(reference)}
    selection['cohort'] = cohort
    selection['N'] = len(cohort)
    selection['artifacts'] = {str(p.relative_to(conservation)): digest(p)
                              for folder in [gffs, faas] for p in folder.iterdir()}
    (conservation / 'cohort.json').write_text(json.dumps(selection, indent=2) + '\n')
    (conservation / 'core_genomes_IDs.txt').write_text('\n'.join(cohort) + '\n')
    print(f'Shared core cohort: {len(cohort)} genomes; inputs in {conservation}')
    return str(Path(output_path).resolve())
