"""Panaroo integration and per-protein presence for CoreCruncher."""
import csv
import os

import pytest

from conftest import GENES, ROARY_METADATA_COLUMNS, write_corecruncher_inputs
from ftscripts import genome

BIOPERL_CDS = 'chr\tGenBank\tCDS\t1\t9\t.\t+\t1\tID=A.p01\n'


def test_panaroo_gff_shifts_bioperl_codon_start_to_a_gff3_phase(tmp_path):
    source = tmp_path / 'in.gff'
    source.write_text('##gff-version 3\n# conversion-by bp_genbank2gff3.pl\n'
                      'chr\tGenBank\tgene\t1\t9\t.\t+\t1\tID=A\n' + BIOPERL_CDS +
                      '##FASTA\n>chr\nATGATGATG\n')
    target = tmp_path / 'out.gff'
    genome.panaroo_gff(str(source), str(target))
    text = target.read_text()
    assert 'chr\tGenBank\tCDS\t1\t9\t.\t+\t0\tID=A.p01\n' in text
    assert 'chr\tGenBank\tgene\t1\t9\t.\t+\t1\tID=A\n' in text      # only CDS carry a phase
    assert target.stat().st_size == source.stat().st_size


def test_panaroo_gff_leaves_files_from_other_converters_alone(tmp_path):
    source = tmp_path / 'in.gff'
    source.write_text('##gff-version 3\n' + BIOPERL_CDS + '##FASTA\n>chr\nATG\n')
    target = tmp_path / 'out.gff'
    genome.panaroo_gff(str(source), str(target))
    assert target.read_text() == source.read_text()


def test_prepare_panaroo_inputs_mirrors_the_cohort_and_rewrites_rebuilt_genomes(tmp_path):
    gff, out = tmp_path / 'gff', tmp_path / 'pan'
    gff.mkdir()
    for name in ('A', 'B'):
        (gff / f'{name}.gff').write_text(f'##gff-version 3\n# {name}\n')
    assert genome.prepare_panaroo_inputs(str(gff), str(out)) == ['A.gff', 'B.gff']
    first = (out / 'A.gff').stat().st_mtime_ns
    genome.prepare_panaroo_inputs(str(gff), str(out))
    assert (out / 'A.gff').stat().st_mtime_ns == first           # unchanged source: reused

    (gff / 'B.gff').unlink()
    (gff / 'A.gff').unlink()
    (gff / 'A.gff').write_text('##gff-version 3\n# rebuilt\n')   # a different file
    assert genome.prepare_panaroo_inputs(str(gff), str(out)) == ['A.gff']
    assert sorted(p.name for p in out.glob('*.gff')) == ['A.gff']
    assert 'rebuilt' in (out / 'A.gff').read_text()


def _write_panaroo(output_path, organism_name, rows, genomes=4):
    conservation = os.path.join(output_path, organism_name, 'conservation')
    results = os.path.join(conservation, 'panaroo_output', 'results')
    gff_dir = os.path.join(conservation, 'panaroo_output', 'gff')
    os.makedirs(results)
    os.makedirs(gff_dir)
    samples = [organism_name] + [f'genome_{i:02d}' for i in range(1, genomes)]
    for sample in samples:
        with open(os.path.join(gff_dir, f'{sample}.gff'), 'w') as handle:
            handle.write('##gff-version 3\n')
            if sample == organism_name:
                for i, gene in enumerate(GENES, start=1):
                    handle.write(f'chr\tsynthetic\tCDS\t{i * 100}\t{i * 100 + 59}\t.\t+\t0\t'
                                 f'ID={gene}.p01;locus_tag={gene}\n')
    with open(os.path.join(results, 'gene_presence_absence_roary.csv'), 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(ROARY_METADATA_COLUMNS + samples)
        for n, (reference_cell, isolates) in enumerate(rows):
            cells = [reference_cell] + [f'g{i}_{n}' if i < isolates else '' for i in range(1, genomes)]
            writer.writerow([f'group_{n}', '', 'p', isolates, isolates, 1.0,
                             1, 1, '', '', '', 180, 180, 180] + cells)


def test_panaroo_table_handles_merged_fragments_and_refound_genes(organism):
    output_path, organism_name = organism
    g1, g2, g3, g4, g5 = GENES
    _write_panaroo(output_path, organism_name, [
        (f'{g1}.p01', 4),                       # in every genome
        (f'{g2}.p01;{g3}.p01', 4),              # two fragments Panaroo merged into one gene
        (f'{g4}.p01', 2),                       # half of the genomes
        (f'{organism_name}_refound_7', 4),      # re-found in unannotated reference sequence
        (f'{g5}.p01', 3),
    ])
    table = genome.panaroo_output(output_path, organism_name, core_threshold=0.99).set_index('gene')
    assert table.loc[g1, 'core_panaroo'] and table.loc[g2, 'core_panaroo'] and table.loc[g3, 'core_panaroo']
    assert not table.loc[g4, 'core_panaroo']
    assert table.loc[g4, 'panaroo_conservation_pct'] == 50
    assert table.loc[g5, 'panaroo_genomes_present'] == 3
    assert set(table['panaroo_occupancy_status']) == {'evaluated'}
    assert os.path.isfile(os.path.join(output_path, organism_name, 'conservation', 'core_panaroo.tsv'))


def test_corecruncher_reports_presence_for_every_protein_not_only_core(organism):
    output_path, organism_name = organism
    g1, g2 = GENES[:2]
    write_corecruncher_inputs(output_path, organism_name, [g1], genomes=10, hits={g2: 4})
    table = genome.corecruncher_output(output_path, organism_name, core_threshold=99).set_index('gene')
    assert table.loc[g1, 'corecruncher_hit_genomes'] == 10 and table.loc[g1, 'corecruncher_hit_pct'] == 100
    assert table.loc[g2, 'corecruncher_hit_genomes'] == 5          # 4 comparators + the reference
    assert table.loc[g2, 'corecruncher_hit_pct'] == 50
    assert not table.loc[g2, 'core_corecruncher']


def test_corecruncher_hit_occupancy_applies_identity_and_length_filters(organism):
    output_path, organism_name = organism
    write_corecruncher_inputs(output_path, organism_name, [], genomes=3)
    cc = os.path.join(output_path, organism_name, 'conservation', 'corecruncher_output')
    cohort = {f'{organism_name}.faa', 'genome_01.faa', 'genome_02.faa'}
    assert set(genome.corecruncher_hit_occupancy(cc, f'{organism_name}.faa', cohort, 95).values()) == {3}
    assert set(genome.corecruncher_hit_occupancy(cc, f'{organism_name}.faa', cohort, 99.5).values()) == {1}
    os.remove(os.path.join(cc, 'CC', f'genome_02.faa-{organism_name}.faa'))
    with pytest.raises(FileNotFoundError):
        genome.corecruncher_hit_occupancy(cc, f'{organism_name}.faa', cohort, 95)
