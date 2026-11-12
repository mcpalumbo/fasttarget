"""Persistent per-genome cache used by the core selection (ftscripts.dereplication)."""
import csv

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from ftscripts import dereplication, genome


def _write_genbank(folder, accession, n_translated, n_untranslated=0, host='Homo sapiens', seq='ATG' * 40):
    folder.mkdir(parents=True)
    record = SeqRecord(Seq(seq), id=accession, name=accession[:16], annotations={'molecule_type': 'DNA'})
    record.features.append(SeqFeature(FeatureLocation(0, len(seq)), type='source',
                                      qualifiers={'host': [host]}))
    for i in range(n_translated):
        record.features.append(SeqFeature(FeatureLocation(0, 9), type='CDS',
                                          qualifiers={'locus_tag': [f'{accession}_{i}'], 'translation': ['MMM']}))
    for i in range(n_untranslated):
        record.features.append(SeqFeature(FeatureLocation(0, 9), type='CDS',
                                          qualifiers={'locus_tag': [f'{accession}_u{i}']}))
    path = folder / f'{accession}.gbff'
    SeqIO.write([record], path, 'genbank')
    return path


def _dataset(tmp_path):
    data = tmp_path / 'data'
    _write_genbank(data / 'GCF_1', 'GCF_1', 5, n_untranslated=2)
    _write_genbank(data / 'GCF_2', 'GCF_2', 3, host='Sus scrofa', seq='ATGC' * 30)
    return data, sorted(p for p in data.iterdir())


def _rows(cache):
    with (cache / 'genomes.tsv').open() as handle:
        return {r['accession']: r for r in csv.DictReader(handle, delimiter='\t')}


def test_scan_counts_translated_cds_and_writes_a_readable_table(tmp_path):
    _, folders = _dataset(tmp_path)
    cache = tmp_path / 'genome_cache'
    result = dereplication.scan_genomes(folders, cache, cpus=2)
    assert int(result['GCF_1'][0]['protein_count']) == 5      # untranslated CDS do not count
    assert int(result['GCF_2'][0]['protein_count']) == 3
    rows = _rows(cache)
    assert set(rows) == {'GCF_1', 'GCF_2'}
    assert (cache / 'fna' / 'GCF_1.fna').is_file()


def test_second_scan_reuses_rows_and_a_changed_genbank_is_rescanned(tmp_path, capsys):
    data, folders = _dataset(tmp_path)
    cache = tmp_path / 'genome_cache'
    dereplication.scan_genomes(folders, cache, cpus=2)
    capsys.readouterr()
    dereplication.scan_genomes(folders, cache, cpus=2)
    assert '2 reused, 0 scanned' in capsys.readouterr().out

    # Rewrite one GenBank with a different CDS count: its row must not be reused.
    (data / 'GCF_2' / 'GCF_2.gbff').unlink()
    (data / 'GCF_2').rmdir()
    _write_genbank(data / 'GCF_2', 'GCF_2', 7, host='Sus scrofa', seq='ATGC' * 30)
    result = dereplication.scan_genomes(sorted(data.iterdir()), cache, cpus=2)
    assert '1 reused, 1 scanned' in capsys.readouterr().out
    assert int(result['GCF_2'][0]['protein_count']) == 7


def test_a_tampered_cached_fasta_is_regenerated(tmp_path, capsys):
    _, folders = _dataset(tmp_path)
    cache = tmp_path / 'genome_cache'
    dereplication.scan_genomes(folders, cache, cpus=1)
    (cache / 'fna' / 'GCF_1.fna').write_text('>broken\nNNNN\n')
    capsys.readouterr()
    dereplication.scan_genomes(folders, cache, cpus=1)
    assert '1 reused, 1 scanned' in capsys.readouterr().out
    assert 'broken' not in (cache / 'fna' / 'GCF_1.fna').read_text()


def test_protein_qc_uses_precomputed_counts_without_parsing(tmp_path):
    destination = tmp_path / 'protein_qc.tsv'
    missing = str(tmp_path / 'does_not_exist.gbff')   # would raise if it were parsed
    excluded = genome.conservation_protein_qc({'A': missing, 'B': missing, 'C': missing},
                                              str(destination), counts={'A': 100, 'B': 100, 'C': 10})
    assert excluded == {'C'}


def _command(workdir='/w1', listing='/w1/g.txt', cpus='8', sa='0.999'):
    return ['dRep', 'dereplicate', workdir, '-g', listing, '-p', cpus, '--S_algorithm', 'fastANI',
            '-pa', '0.998', '-sa', sa, '-nc', '0.90']


CANDIDATES = {'GCF_1.fna': {'fna_sha256': 'a'}, 'GCF_2.fna': {'fna_sha256': 'b'}}
VERSIONS = {'dRep': '3.6.2', 'fastANI': '1.34', 'mash': '2.3'}


def test_drep_key_ignores_paths_and_threads_but_not_thresholds_genomes_or_versions():
    key = dereplication.drep_cache_key(CANDIDATES, _command(), VERSIONS)
    assert key == dereplication.drep_cache_key(CANDIDATES, _command('/w2', '/w2/g.txt', '16'), VERSIONS)
    assert key != dereplication.drep_cache_key(CANDIDATES, _command(sa='0.995'), VERSIONS)
    assert key != dereplication.drep_cache_key({**CANDIDATES, 'GCF_2.fna': {'fna_sha256': 'c'}}, _command(), VERSIONS)
    assert key != dereplication.drep_cache_key(CANDIDATES, _command(), {**VERSIONS, 'dRep': '3.6.3'})


def test_drep_result_round_trip_and_rejects_other_genome_sets(tmp_path):
    clusters, scores = {'GCF_1.fna': '1_1', 'GCF_2.fna': '1_1'}, {'GCF_1.fna': 0.9995, 'GCF_2.fna': 0.9991}
    assert dereplication.load_drep_result(tmp_path, CANDIDATES, _command(), VERSIONS) is None
    dereplication.save_drep_result(tmp_path, CANDIDATES, _command(), VERSIONS, clusters, scores)
    loaded = dereplication.load_drep_result(tmp_path, CANDIDATES, _command('/other', '/other/g', '2'), VERSIONS)
    assert loaded[0] == clusters and loaded[1] == scores
    fewer = {'GCF_1.fna': CANDIDATES['GCF_1.fna']}
    assert dereplication.load_drep_result(tmp_path, fewer, _command(), VERSIONS) is None


def test_sync_inputs_holds_exactly_the_cohort_and_relinks_changed_sources(tmp_path):
    src = tmp_path / 'src'
    src.mkdir()
    for name in ('A.gff', 'B.gff', 'C.gff'):
        (src / name).write_text(name)
    folder = tmp_path / 'gff'
    dereplication._sync_inputs(folder, {n: src / n for n in ('A.gff', 'B.gff')}, link=True)
    assert sorted(p.name for p in folder.iterdir()) == ['A.gff', 'B.gff']
    assert (folder / 'A.gff').samefile(src / 'A.gff')

    # Next cohort drops B, adds C, and A now comes from a different (rebuilt) source.
    (tmp_path / 'A2.gff').write_text('rebuilt')
    dereplication._sync_inputs(folder, {'A.gff': tmp_path / 'A2.gff', 'C.gff': src / 'C.gff'}, link=True)
    assert sorted(p.name for p in folder.iterdir()) == ['A.gff', 'C.gff']
    assert (folder / 'A.gff').read_text() == 'rebuilt'

    copies = tmp_path / 'faa'
    dereplication._sync_inputs(copies, {'A.gff': src / 'A.gff'}, link=False)
    assert not (copies / 'A.gff').samefile(src / 'A.gff')
