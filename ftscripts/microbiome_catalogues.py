import os


MGNIFY_CATALOGUES = {
    "human-gut": {
        "number_of_species": 4744,
        "ftp_site": "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2",
    },
    "human-oral": {
        "number_of_species": 452,
        "ftp_site": "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-oral/v1.0.1",
    },
    "human-skin": {
        "number_of_species": 579,
        "ftp_site": "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-skin/v1.0",
    },
    "human-vaginal": {
        "number_of_species": 280,
        "ftp_site": "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-vaginal/v1.0",
    },
    "marine": {
        "number_of_species": 13223,
        "ftp_site": "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/marine/v2.0",
    },
    "soil": {
        "number_of_species": 19472,
        "ftp_site": "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/soil/v1.0",
    },
    "barley-rhizosphere": {
        "number_of_species": 86,
        "ftp_site": "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/barley-rhizosphere/v2.0",
    },
    "maize-rhizosphere": {
        "number_of_species": 336,
        "ftp_site": "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/maize-rhizosphere/v1.0",
    },
    "tomato-rhizosphere": {
        "number_of_species": 579,
        "ftp_site": "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/tomato-rhizosphere/v1.0",
    },
}


def get_catalogue(catalogue_name):
    try:
        return MGNIFY_CATALOGUES[catalogue_name]
    except KeyError as error:
        supported = ", ".join(sorted(MGNIFY_CATALOGUES))
        raise ValueError(
            f"Unsupported microbiome catalogue '{catalogue_name}'. "
            f"Supported catalogues: {supported}."
        ) from error


def catalogue_species_path(databases_path, catalogue_name, allow_legacy=True):
    catalogue_path = os.path.join(
        databases_path,
        "microbiomes",
        catalogue_name,
        "species_catalogue",
    )
    legacy_path = os.path.join(databases_path, "species_catalogue")
    if (
        allow_legacy
        and catalogue_name == "human-gut"
        and not os.path.isdir(catalogue_path)
        and os.path.isdir(legacy_path)
    ):
        return legacy_path
    return catalogue_path


def catalogue_column_prefix(catalogue_name):
    return catalogue_name.replace("-", "_")


def configured_catalogues(offtarget_config):
    catalogue_configs = offtarget_config.get("microbiome_catalogues")
    if catalogue_configs is None:
        catalogue_configs = [{
            "name": "human-gut",
            "identity_filter": offtarget_config.get("microbiome_identity_filter", 40),
            "coverage_filter": offtarget_config.get("microbiome_coverage_filter", 70),
        }]
    return catalogue_configs
