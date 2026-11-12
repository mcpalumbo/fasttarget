"""Roary needs every GFF feature seqid to have a sequence in the ##FASTA section."""
import pytest

from ftscripts import genome

REGION = 'NC_1\tGenBank\tDNA\t1\t9\t.\t.\t1\tID=NC_1;Alias=NC_1.1\n'


def _gff(tmp_path, feature_seqid, fasta_header='NC_1'):
    path = tmp_path / 'ref.gff'
    path.write_text(f'##gff-version 3\n##sequence-region NC_1 1 9\n{REGION}'
                    f'{feature_seqid}\tGenBank\tCDS\t1\t9\t.\t+\t0\tID=A_1.p01;locus_tag=A_1\n'
                    f'##FASTA\n>{fasta_header}\nATGATGATG\n>A_1.p01\nMMM\n')
    return path


def test_versioned_feature_seqid_is_rewritten_to_the_fasta_record(tmp_path):
    # Shape of the Pseudomonas Genome DB PA14 GenBank after bp_genbank2gff3.
    path = _gff(tmp_path, 'NC_1.1')
    assert genome.reconcile_gff_seqids(str(path)) == {'NC_1.1': 'NC_1'}
    text = path.read_text()
    assert 'NC_1\tGenBank\tCDS\t1\t9' in text and 'NC_1.1\t' not in text
    assert REGION in text and '>NC_1\nATGATGATG\n>A_1.p01\nMMM\n' in text
    assert genome.reconcile_gff_seqids(str(path)) == {}


def test_consistent_gff_is_left_untouched(tmp_path):
    path = _gff(tmp_path, 'NC_1')
    before = path.read_bytes()
    assert genome.reconcile_gff_seqids(str(path)) == {}
    assert path.read_bytes() == before


def test_unrelated_fasta_record_is_an_error(tmp_path):
    path = _gff(tmp_path, 'NC_1.1', fasta_header='OTHER_9')
    with pytest.raises(ValueError, match='no sequence'):
        genome.reconcile_gff_seqids(str(path))
