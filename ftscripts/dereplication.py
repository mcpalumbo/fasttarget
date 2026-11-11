"""Isolated dRep selection; input datasets are read-only sources."""
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from Bio import SeqIO
import pandas as pd


def digest(path):
    with open(path, 'rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


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
    return enabled, threshold


def select(dataset, destination, min_genomes=300, cpus=4, reference=None, human_only=False):
    """Run above threshold, counting downloaded assemblies before host selection.

    Returns a manifest. Destination must be new; never resume partial dRep output.
    Native dRep/FastANI coverage uses aligned query fragments (larger method).
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

    def read(path):
        records = list(SeqIO.parse(path, 'genbank'))
        if not records or any(not len(record.seq) for record in records):
            raise ValueError(f'Empty genome: {path}')
        sequences = sorted(hashlib.sha256(str(r.seq).upper().encode()).hexdigest() for r in records)
        return records, sequences

    ref_sequences = read(reference)[1] if reference else None
    folders = sorted(p for p in Path(dataset).iterdir() if p.is_dir() and p.name.startswith(('GCF_', 'GCA_')))
    for folder in folders:
        paths = sorted(folder.glob('*.gbff')) or sorted(folder.glob('*.gbk'))
        if len(paths) != 1:
            raise ValueError(f'Expected one GenBank in {folder}, found {len(paths)}')
        path = paths[0]
        records, sequences = read(path)
        if sequences == ref_sequences:
            excluded[folder.name] = 'identical_reference_sequences'
            continue
        hosts = {h for r in records for f in r.features if f.type == 'source'
                 for h in f.qualifiers.get('host', [])}
        if human_only and 'Homo sapiens' not in hosts:
            excluded[folder.name] = 'host_filter'
            continue
        fasta = inputs / f'{folder.name}.fna'
        SeqIO.write(records, fasta, 'fasta')
        candidates[fasta.name] = {'accession': folder.name, 'gbk': str(path.resolve()),
                                  'gbk_sha256': digest(path), 'fna_sha256': digest(fasta)}
    active = len(folders) > min_genomes and len(candidates) > 1
    clusters = {name: name for name in candidates}
    scores = {name: 0.0 for name in candidates}
    command = None
    versions = {}
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
                   '-p', str(cpus), '--S_algorithm', 'fastANI', '-sa', '0.999', '-nc', '0.90',
                   '-cm', 'larger', '--clusterAlg', 'average', '--ignoreGenomeQuality', '-l', '0',
                   '-comW', '0', '-conW', '0', '-strW', '0', '-N50W', '0', '-sizeW', '0', '-centW', '1']
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
    representatives = []
    for cluster in sorted(set(clusters.values())):
        members = [name for name in candidates if clusters[name] == cluster]
        representatives.append(min(members, key=lambda name: (-scores[name], name)))
    manifest = {'downloaded_count': len(folders), 'eligible_count': len(candidates),
                'min_genomes': min_genomes, 'applied': active, 'versions': versions,
                'command': command, 'excluded': excluded, 'candidates': candidates,
                'policy': {'ani': 0.999, 'coverage': 0.90, 'linkage': 'average',
                           'backend': 'fastANI', 'primary_ani': 0.90,
                           'coverage_semantics': 'dRep native directional fragment fraction > 0.90 before reciprocal ANI averaging',
                           'representative': 'maximum dRep within-cluster mean ANI; accession lexical tie-break'},
                'clusters': clusters, 'centrality': scores, 'representatives': representatives,
                'qc_status': 'host and readable sequence only; contamination/low-CDS policy pending'}
    (destination / 'selection.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def prepare(output_path, organism, core, cpus, human_only, engine):
    """Prepare a fresh shared cohort workspace for the Python core entry point."""
    from ftscripts import genome
    enabled, threshold = options(core)
    if not enabled:
        return output_path
    base = Path(output_path).resolve() / organism
    runs = base / 'conservation' / 'cohort_runs'
    runs.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='cohort-', dir=runs))
    reference = base / 'genome' / f'{organism}.gbk'
    selection = select(base / 'conservation' / f'{organism}_dataset' / 'ncbi_dataset' / 'data',
                       work / 'selection', threshold, cpus, reference, human_only)
    target = work / organism
    shutil.copytree(base / 'genome', target / 'genome')
    conservation = target / 'conservation'
    gffs = conservation / 'roary_output' / 'gff'
    faas = conservation / 'corecruncher_output' / 'faa'
    gffs.mkdir(parents=True)
    faas.mkdir(parents=True)
    for suffix, folder in [('gff', gffs), ('faa', faas)]:
        shutil.copyfile(target / 'genome' / f'{organism}.{suffix}', folder / f'{organism}.{suffix}')
    cohort = [organism]
    for name in selection['representatives']:
        info = selection['candidates'][name]
        accession = info['accession']
        if digest(info['gbk']) != info['gbk_sha256']:
            raise RuntimeError(f'Comparator changed after selection: {accession}')
        if accession == organism:
            raise ValueError('Comparator accession collides with reference name')
        prepared = conservation / 'prepared' / accession
        prepared.mkdir(parents=True)
        genome.ref_genome_files(info['gbk'], str(prepared), accession, container_engine=engine)
        shutil.copyfile(prepared / f'{accession}.gff', gffs / f'{accession}.gff')
        shutil.copyfile(prepared / f'{accession}.faa', faas / f'{accession}.faa')
        cohort.append(accession)
    for folder, suffix in [(gffs, '.gff'), (faas, '.faa')]:
        if {p.stem for p in folder.glob('*' + suffix)} != set(cohort):
            raise ValueError('Prepared tools have different cohorts')
    selection['reference'] = {'name': organism, 'gbk_sha256': digest(reference)}
    selection['cohort'] = cohort
    selection['N'] = len(cohort)
    selection['artifacts'] = {str(p.relative_to(work)): digest(p)
                              for folder in [gffs, faas] for p in folder.iterdir()}
    (conservation / 'cohort.json').write_text(json.dumps(selection, indent=2) + '\n')
    (conservation / 'core_genomes_IDs.txt').write_text('\n'.join(cohort) + '\n')
    print(f'Shared core cohort: {len(cohort)} genomes; workspace: {work}')
    return str(work)
