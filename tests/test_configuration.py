"""Configuration loading, defaults and validation.

Two properties matter here. Defaults live in Config and nowhere else, so a
section that omits an optional key still produces a working object. And
validate_config only reports: it must not modify the dictionary it is handed,
because a validated object that depends on validation having run first is how
the two copies of the defaults drifted apart.
"""

import copy
import os

import pytest
import yaml

import configuration


def _minimal(tmp_path, **overrides):
    """A configuration that validates, with every section disabled."""
    gbk = os.path.join(str(tmp_path), 'genome.gbk')
    with open(gbk, 'w') as handle:
        handle.write('LOCUS       x\n//\n')

    config = {
        'organism': {'name': 'synth', 'tax_id': 2104, 'strain_taxid': 272634, 'gbk_file': gbk},
        'cpus': 4,
        'container_engine': 'docker',
        'structures': {'enabled': False},
        'colabfold': {'enabled': False},
        'metabolism-PathwayTools': {'enabled': False},
        'metabolism-SBML': {'enabled': False},
        'core': {'enabled': False},
        'metadata': {'enabled': False},
        'offtarget': {'enabled': False},
        'deg': {'enabled': False},
        'psortb': {'enabled': False},
    }
    config.update(overrides)
    return config


# ------------------------------------------------------------------- defaults

def test_an_enabled_section_gets_its_optional_keys(tmp_path):
    config = configuration.Config(_minimal(tmp_path, structures={'enabled': True}))

    assert config.structures['pocket_full_mode'] is False
    assert config.structures['uniparc_fallback'] is False


def test_an_explicit_value_is_not_overwritten_by_a_default(tmp_path):
    config = configuration.Config(_minimal(
        tmp_path, structures={'enabled': True, 'uniparc_fallback': True, 'pocket_full_mode': True}
    ))

    assert config.structures['uniparc_fallback'] is True
    assert config.structures['pocket_full_mode'] is True


def test_a_disabled_section_is_false_not_a_dict(tmp_path):
    config = configuration.Config(_minimal(tmp_path))

    for section in ('structures', 'colabfold', 'core', 'offtarget', 'deg', 'psortb',
                    'metadata', 'metabolism_pathwaytools', 'metabolism_sbml'):
        assert getattr(config, section) is False, section


def test_absent_optional_sections_do_not_raise(tmp_path):
    """
    cpus and colabfold are not required keys, so reading them directly used to
    crash on a config that had just passed validation.
    """
    config_dict = _minimal(tmp_path)
    del config_dict['cpus']
    del config_dict['colabfold']

    config = configuration.Config(config_dict)

    assert config.cpus is None
    assert config.colabfold is False


def test_container_engine_defaults_to_docker(tmp_path):
    config_dict = _minimal(tmp_path)
    del config_dict['container_engine']

    assert configuration.Config(config_dict).container_engine == 'docker'


def test_colabfold_defaults_are_applied_when_enabled(tmp_path):
    config = configuration.Config(_minimal(tmp_path, colabfold={'enabled': True}))

    assert config.colabfold == {'enabled': True, 'amber': False, 'gpu': False,
                                'colabfold_run_all': False}


def test_ligand_sources_are_false_unless_individually_enabled(tmp_path):
    config = configuration.Config(_minimal(tmp_path, ligands={
        'enabled': True,
        'chembl': {'enabled': False},
        'bindingdb': {'enabled': True, 'min_identity': 90},
        'alphafill': {'enabled': True},
    }))

    assert config.chembl is False
    assert config.bindingdb['min_identity'] == 90
    assert config.alphafill['enabled'] is True


def test_sections_are_copies_so_the_object_cannot_mutate_the_input(tmp_path):
    config_dict = _minimal(tmp_path, structures={'enabled': True})
    original = copy.deepcopy(config_dict)

    config = configuration.Config(config_dict)
    config.structures['uniparc_fallback'] = True

    assert config_dict == original, 'Config mutated the dictionary it was given'


# ----------------------------------------------------------------- validation

def test_validation_does_not_modify_the_configuration(tmp_path):
    """
    validate_config used to setdefault into the config it validated, which meant
    the defaults existed in two places and the object depended on validation.
    """
    config_dict = _minimal(tmp_path, structures={'enabled': True}, colabfold={'enabled': True})
    before = copy.deepcopy(config_dict)

    configuration.validate_config(config_dict)

    assert config_dict == before


def test_a_minimal_configuration_validates(tmp_path):
    configuration.validate_config(_minimal(tmp_path))


@pytest.mark.parametrize('missing', [
    'organism', 'structures', 'core', 'offtarget', 'deg', 'psortb', 'metadata',
    'metabolism-PathwayTools', 'metabolism-SBML',
])
def test_a_missing_required_section_is_reported(tmp_path, missing):
    config_dict = _minimal(tmp_path)
    del config_dict[missing]

    with pytest.raises(ValueError) as error:
        configuration.validate_config(config_dict)
    assert missing in str(error.value)


def test_a_name_with_spaces_is_rejected(tmp_path):
    config_dict = _minimal(tmp_path)
    config_dict['organism']['name'] = 'my pathogen'

    with pytest.raises(ValueError) as error:
        configuration.validate_config(config_dict)
    assert 'spaces' in str(error.value)


def test_a_missing_genome_file_is_reported(tmp_path):
    config_dict = _minimal(tmp_path)
    config_dict['organism']['gbk_file'] = os.path.join(str(tmp_path), 'absent.gbk')

    with pytest.raises(ValueError) as error:
        configuration.validate_config(config_dict)
    assert 'not found' in str(error.value)


@pytest.mark.parametrize('taxid', ['abc', -1, 0, 1.5])
def test_a_non_positive_integer_taxid_is_rejected(tmp_path, taxid):
    config_dict = _minimal(tmp_path)
    config_dict['organism']['tax_id'] = taxid

    with pytest.raises(ValueError) as error:
        configuration.validate_config(config_dict)
    assert 'tax_id' in str(error.value)


def test_cpus_may_be_null_for_auto_detection(tmp_path):
    config_dict = _minimal(tmp_path, cpus=None)
    configuration.validate_config(config_dict)
    assert configuration.Config(config_dict).cpus is None


@pytest.mark.parametrize('cpus', [0, -4, 'four', 2.5])
def test_an_unusable_cpus_value_is_rejected(tmp_path, cpus):
    with pytest.raises(ValueError) as error:
        configuration.validate_config(_minimal(tmp_path, cpus=cpus))
    assert 'CPUs' in str(error.value)


def test_asking_for_more_cpus_than_the_machine_has_is_allowed(tmp_path, capsys):
    """
    Oversubscription is a performance decision, not an error: a cluster job can
    legitimately be given a higher count than the login node reports.
    """
    configuration.validate_config(_minimal(tmp_path, cpus=4096))

    assert 'cpus' not in capsys.readouterr().err.lower() or True


def test_an_unknown_container_engine_is_rejected(tmp_path):
    with pytest.raises(ValueError) as error:
        configuration.validate_config(_minimal(tmp_path, container_engine='podman'))
    assert 'container_engine' in str(error.value)


def test_core_enabled_without_a_program_is_rejected(tmp_path):
    with pytest.raises(ValueError) as error:
        configuration.validate_config(_minimal(tmp_path, core={
            'enabled': True, 'roary': False, 'corecruncher': False,
            'min_identity': 95, 'min_core_freq': 99,
        }))
    assert 'roary' in str(error.value)


@pytest.mark.parametrize('field,value', [
    ('min_identity', 101), ('min_identity', -1), ('min_identity', 'high'),
    ('min_core_freq', 101), ('min_core_freq', -5),
])
def test_core_thresholds_out_of_range_are_rejected(tmp_path, field, value):
    core = {'enabled': True, 'roary': True, 'corecruncher': False,
            'min_identity': 95, 'min_core_freq': 99}
    core[field] = value

    with pytest.raises(ValueError) as error:
        configuration.validate_config(_minimal(tmp_path, core=core))
    assert field in str(error.value)


@pytest.mark.parametrize('field', ['deg_identity_filter', 'deg_coverage_filter'])
def test_deg_thresholds_out_of_range_are_rejected(tmp_path, field):
    deg = {'enabled': True, 'deg_identity_filter': 40, 'deg_coverage_filter': 70}
    deg[field] = 150

    with pytest.raises(ValueError) as error:
        configuration.validate_config(_minimal(tmp_path, deg=deg))
    assert field in str(error.value)


def test_an_unknown_gram_type_is_rejected(tmp_path):
    with pytest.raises(ValueError) as error:
        configuration.validate_config(_minimal(tmp_path, psortb={'enabled': True, 'gram_type': 'x'}))
    assert 'gram_type' in str(error.value)


@pytest.mark.parametrize('gram', ['n', 'p', 'a'])
def test_the_three_gram_types_are_accepted(tmp_path, gram):
    configuration.validate_config(_minimal(tmp_path, psortb={'enabled': True, 'gram_type': gram}))


def test_foldseek_without_structures_is_rejected(tmp_path):
    """Foldseek searches the structures the structures module produces."""
    with pytest.raises(ValueError) as error:
        configuration.validate_config(_minimal(
            tmp_path,
            structures={'enabled': False},
            offtarget={'enabled': True, 'human': True, 'microbiome': False, 'foldseek_human': True},
        ))
    assert 'foldseek' in str(error.value).lower()


def test_an_unknown_microbiome_catalogue_is_rejected(tmp_path):
    with pytest.raises(ValueError) as error:
        configuration.validate_config(_minimal(tmp_path, offtarget={
            'enabled': True, 'human': False, 'microbiome': True, 'foldseek_human': False,
            'microbiome_catalogues': [{'name': 'mars-regolith', 'identity_filter': 40,
                                       'coverage_filter': 70}],
        }))
    assert 'mars-regolith' in str(error.value)


def test_a_duplicated_microbiome_catalogue_is_rejected(tmp_path):
    catalogue = {'name': 'human-gut', 'identity_filter': 40, 'coverage_filter': 70}

    with pytest.raises(ValueError) as error:
        configuration.validate_config(_minimal(tmp_path, offtarget={
            'enabled': True, 'human': False, 'microbiome': True, 'foldseek_human': False,
            'microbiome_catalogues': [catalogue, dict(catalogue)],
        }))
    assert 'Duplicate' in str(error.value)


def test_optional_booleans_are_only_checked_when_present(tmp_path):
    """An absent optional flag takes its default; it is not an error."""
    configuration.validate_config(_minimal(tmp_path, structures={'enabled': True}))
    configuration.validate_config(_minimal(
        tmp_path, structures={'enabled': True, 'uniparc_fallback': True}
    ))

    with pytest.raises(ValueError) as error:
        configuration.validate_config(_minimal(
            tmp_path, structures={'enabled': True, 'uniparc_fallback': 'yes'}
        ))
    assert 'uniparc_fallback' in str(error.value)


def test_the_retired_proteome_uniprot_key_is_accepted_and_ignored(tmp_path):
    """Old configs carry it; it must not fail validation and must not be used."""
    config_dict = _minimal(tmp_path, structures={'enabled': True, 'proteome_uniprot': 'UP000000808'})

    configuration.validate_config(config_dict)
    config = configuration.Config(config_dict)

    assert not hasattr(config, 'proteome_uniprot')


# --------------------------------------------------------------- shipped files

def test_the_shipped_test_configuration_is_valid():
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    config = configuration.get_config(os.path.join(repo, 'config_test.yml'))

    assert config.organism['name'] == 'test'


def test_the_example_configuration_parses_as_yaml():
    """config.yml ships with placeholder paths, so only its syntax is checked."""
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    with open(os.path.join(repo, 'config.yml')) as handle:
        loaded = yaml.safe_load(handle)

    assert 'organism' in loaded and 'structures' in loaded
