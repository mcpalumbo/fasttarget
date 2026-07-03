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
    """
    Returns the configuration of a supported MGnify catalogue.

    :param catalogue_name: Name of the microbiome catalogue.
    :return: Dictionary containing the catalogue configuration.
    :raises ValueError: If the catalogue is not supported.
    """
    try:
        return MGNIFY_CATALOGUES[catalogue_name]
    except KeyError as error:
        supported = ", ".join(sorted(MGNIFY_CATALOGUES))
        raise ValueError(
            f"Unsupported microbiome catalogue '{catalogue_name}'. "
            f"Supported catalogues: {supported}."
        ) from error


def catalogue_species_path(databases_path, catalogue_name):
    """
    Returns the species catalogue directory for a microbiome catalogue.

    :param databases_path: Path to the databases folder.
    :param catalogue_name: Name of the microbiome catalogue.
    :return: Path to the species catalogue directory.
    """
    return os.path.join(
        databases_path,
        "microbiomes",
        catalogue_name,
        "species_catalogue",
    )


def catalogue_column_prefix(catalogue_name):
    """
    Creates a column-safe prefix from a microbiome catalogue name.

    :param catalogue_name: Name of the microbiome catalogue.
    :return: Catalogue name with hyphens replaced by underscores.
    """
    return catalogue_name.replace("-", "_")
