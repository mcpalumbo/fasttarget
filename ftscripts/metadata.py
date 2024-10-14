from Bio import SeqIO
import pandas as pd
import os
from ftscripts import files

def ref_gbk_locus(base_path, organism_name):
    """
    Returns a list of locus_tags from the reference genome.

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.

    :return: List of locus_tags.
    """
    ref_gbk = os.path.join(base_path, 'organism', organism_name, 'genome', f'{organism_name}.gbk')

    locus_tags = []

    for record in SeqIO.parse(ref_gbk, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if 'translation' in feature.qualifiers and 'locus_tag' in feature.qualifiers:
                    locus_tags.append(feature.qualifiers["locus_tag"][0])
    
    return locus_tags

def metadata_table_bool(base_path, organism_name, locus_tag_true:str, property:str, out_dir:str):
    """
    Makes a metadata table. Each locus_tag has a boolean value for a property.    

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.
    :param locus_tag_true: List of locus_tag TRUE for a property.
    :param property: Name of the property.
    :param out_dir : Output directory path.

    :return: Metadata table. Each locus_tag has TRUE/FALSE value for a property.
    """

    locus_tags = ref_gbk_locus(base_path, organism_name)

    data = {
        "gene": locus_tags,
        property: [locus_tag in locus_tag_true for locus_tag in locus_tags]
    }
    
    metadata_table = pd.DataFrame(data)
    metadata_table.to_csv(f"{out_dir}/{property}.csv", index=False)
    metadata_table.to_csv(f"{out_dir}/{property}.tsv", index=False, sep='\t')

    print(f"{out_dir}/{property}.csv and .tsv have been created.")

    return metadata_table

def metadata_table_with_values(base_path, organism_name, values_dict:str, property:str, out_dir:str, neg_value=None):
    """
    Makes a metadata table. Each locus_tag has a numerical or categorical value for a property.

    :param base_path: Base path where the repository data is stored.
    :param organism_name: Name of the organism.
    :param values_dict: Dictionary with locus_tag and value.
    :param property: Name of the property.
    :param neg_value: The value assigned to the `locus_tag` if it lacks the specified property. Defaults to None.
    :param out_dir : Output directory path.

    :return: Metadata table. Each locus_tag has a numerical or categorical value for a property.
    """

    locus_tags = ref_gbk_locus(base_path, organism_name)

    data = {"gene": [], property: []}
    for locus_tag in locus_tags:
        if locus_tag in values_dict.keys():
            data["gene"].append(locus_tag)
            data[property].append(values_dict[locus_tag])
        else:
            data["gene"].append(locus_tag)
            data[property].append(neg_value)            
    
    metadata_table = pd.DataFrame(data)
    metadata_table.to_csv(f"{out_dir}/{property}.csv", index=False)
    metadata_table.to_csv(f"{out_dir}/{property}.tsv", index=False, sep='\t')

    print(f"{out_dir}/{property}.csv and .tsv have been created.")

    return metadata_table

def tables_for_TP(organism_name, results_path):
    """
    Generates separate metadata tables for each property in the results table.
    This tables can be imported as metadata in Target Pathogen.
    They are saved in the 'tables_for_TP' directory.

    :param organism_name: Name of the organism
    :param results_path: Path to the results directory.
    """

    TP_metadata_path = os.path.join(results_path, 'tables_for_TP')

    results_table_path = os.path.join(results_path, f'{organism_name}_results_table.tsv')

    os.makedirs(TP_metadata_path, exist_ok=True)
    print(f"{TP_metadata_path} has been created.")

    if files.file_check(results_table_path):
        
        results_table = pd.read_csv(results_table_path, sep='\t', header=0)

        # Generate separate metadata tables for each property

        for column in results_table.columns[1:]:
            sub_table = results_table[['gene', column]]
            sub_table.to_csv(f"{TP_metadata_path}/{column}.tsv", index=False, sep='\t')
            print(f"{TP_metadata_path}/{column}.tsv has been created.")
    
    else:
        print(f"{results_table_path} not found.")


