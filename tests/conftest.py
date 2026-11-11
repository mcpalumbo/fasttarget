"""Synthetic fixtures for the fast test suite.

These build a tiny organism output directory from scratch: five genes, a few
hundred residues in total, no containers, no network and no databases. They
exercise the parsing, caching and reporting logic, which is where every bug
found during the pipeline review lived.

This suite does not replace the twenty-protein test case. That one runs the real
tools on real data and is the only thing that can catch an environment problem,
a container argument or an external API change. These tests catch regressions in
the logic, in seconds, on any machine.
"""

import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

ORGANISM = 'synth'

# Five genes is enough for every case the suite needs: one that every module
# resolves, one that no module resolves, and three to make counts meaningful.
GENES = ['SYN_0001', 'SYN_0002', 'SYN_0003', 'SYN_0004', 'SYN_0005']

# Short, distinct, valid protein sequences.
SEQUENCES = {
    'SYN_0001': 'MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCK',
    'SYN_0002': 'MSNQVLQNLSSTVNNLSHTVSSLSQSVSALSHSVSQQSAMMNHLSQSVSHLSQSVSSLSQ',
    'SYN_0003': 'MAHHHHHHVGTGSNDDDDKSPDPMTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIE',
    'SYN_0004': 'MGSSHHHHHHSSGLVPRGSHMRGSHHHHHHGMASMTGGQQMGRDLYDDDDKDRWGSELEQ',
    'SYN_0005': 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG',
}


def _genbank_record():
    """
    A minimal GenBank record with one CDS per gene.

    Built with Biopython rather than written by hand, so the file is well formed:
    a malformed fixture parses with warnings and can mask a real problem.
    """
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    features = []
    position = 0
    for index, gene in enumerate(GENES, start=1):
        sequence = SEQUENCES[gene]
        length = len(sequence) * 3
        features.append(SeqFeature(
            FeatureLocation(position, position + length, strand=1),
            type='CDS',
            qualifiers={
                'locus_tag': [gene],
                'gene': [f'g{index}'],
                'product': [f'synthetic protein {index}'],
                'translation': [sequence],
            },
        ))
        position += length

    record = SeqRecord(
        Seq('ATGCATGCAT' * ((position // 10) + 1)),
        id='SYNTH001',
        name='SYNTH001',
        description='Synthetic organism for FastTarget tests',
        annotations={'molecule_type': 'DNA', 'topology': 'circular', 'organism': 'synthetic construct'},
        features=features,
    )
    return record


@pytest.fixture
def organism(tmp_path):
    """
    A synthetic organism output directory.

    :return: Tuple of (output_path, organism_name).
    """
    from ftscripts import files

    output_path = str(tmp_path)
    files.create_organism_subfolders(output_path, ORGANISM)

    from Bio import SeqIO

    genome_dir = os.path.join(output_path, ORGANISM, 'genome')
    with open(os.path.join(genome_dir, f'{ORGANISM}.gbk'), 'w') as handle:
        SeqIO.write(_genbank_record(), handle, 'genbank')

    faa_path = os.path.join(genome_dir, f'{ORGANISM}.faa')
    with open(faa_path, 'w') as handle:
        for gene in GENES:
            handle.write(f'>{gene}\n{SEQUENCES[gene]}\n')

    return output_path, ORGANISM


def write_blast_output(path, rows):
    """
    Write a BLAST tabular file with the 14 columns the pipeline reads.

    :param path: Destination path.
    :param rows: Iterable of (qseqid, sseqid, pident, qcovs) tuples.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as handle:
        for qseqid, sseqid, pident, qcovs in rows:
            handle.write(
                '\t'.join(str(value) for value in [
                    qseqid, sseqid, pident, 100, 0, 0, 1, 100, 1, 100,
                    '1e-50', 200.0, qcovs, qcovs,
                ]) + '\n'
            )


def write_roary_inputs(output_path, organism_name, presence, genomes=10):
    """
    Write the Roary outputs the parser reads: a presence table and a GFF.

    :param output_path: Output directory.
    :param organism_name: Organism name.
    :param presence: Mapping of gene to the number of isolates it appears in.
    :param genomes: Total number of isolates.
    """
    conservation = os.path.join(output_path, organism_name, 'conservation')
    results_dir = os.path.join(conservation, 'roary_output', 'results')
    gff_dir = os.path.join(conservation, 'roary_output', 'gff')
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(gff_dir, exist_ok=True)

    with open(os.path.join(results_dir, 'gene_presence_absence.csv'), 'w') as handle:
        handle.write(f'"No. isolates","{organism_name}"\n')
        # One row per gene, plus a row carrying the isolate maximum so the parser
        # computes the same denominator a real run would.
        handle.write(f'{genomes},"{GENES[0]}"\n')
        for gene, isolates in presence.items():
            handle.write(f'{isolates},"{gene}"\n')

    with open(os.path.join(gff_dir, f'{organism_name}.gff'), 'w') as handle:
        handle.write('##gff-version 3\n')
        for index, gene in enumerate(GENES, start=1):
            start = index * 100
            handle.write(
                f'{organism_name}_chr\tsynthetic\tCDS\t{start}\t{start + 59}\t.\t+\t0\t'
                f'ID={gene};locus_tag={gene}\n'
            )


def write_corecruncher_inputs(output_path, organism_name, core_genes):
    """
    Write the families_core.txt the CoreCruncher parser reads.

    :param output_path: Output directory.
    :param organism_name: Organism name.
    :param core_genes: Genes that belong to a core family.
    """
    ccruncher = os.path.join(output_path, organism_name, 'conservation', 'corecruncher_output')
    os.makedirs(os.path.join(ccruncher, 'faa'), exist_ok=True)

    with open(os.path.join(ccruncher, 'families_core.txt'), 'w') as handle:
        for gene in core_genes:
            # The parser keys on '<reference>.faa&<locus_tag>'.
            handle.write(f'{organism_name}.faa&{gene}\tother.faa&X_{gene}\n')
