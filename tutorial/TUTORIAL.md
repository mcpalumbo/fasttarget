# Tutorial: How to Obtain Data for FastTarget

This tutorial will guide you through the steps to obtain the necessary data for using this pipeline.

## 1. Getting the GenBank (gbk) File

To get the GenBank file (`.gbk`) for your organism, you have a couple of options:

- **NCBI Database**: Download the GenBank file directly from the NCBI website. Visit the [NCBI Genome Database](https://www.ncbi.nlm.nih.gov/genome/) and search for your organism. Once you find the desired genome, you can download the GenBank file from the "Download" section.

- **Annotation Programs**: If you're using an annotation tool like Prokka or PGAP, you can generate the GenBank file as part of the annotation process. Follow the documentation of the respective tool to obtain the GenBank file.

## 2. Obtaining Metabolic Files

To obtain metabolic files, you can use the Pathway Tools program. Here's a brief guide:

1. **Install Pathway Tools**: Download and install Pathway Tools from the [Pathway Tools website](https://bioinformatics.ai.sri.com/ptools/).

2. **Generate Metabolic Data**: Use Pathway Tools to generate metabolic files for your organism. This usually involves loading your GenBank file into Pathway Tools and running the metabolic pathway analysis.

3. **Save and Export**: Once the analysis is complete, export the metabolic files in the required format for use in this pipeline.

## 3. Finding the Proteome ID

To find the proteome ID for your organism:

1. **Visit UniProt**: Go to the [UniProt website](https://www.uniprot.org/).

2. **Change Database**: Before searching, make sure to select "Proteomes" from the dropdown menu in the search bar.

3. **Search for Your Organism**: Enter the name of your organism in the search bar.

4. **Obtain Proteome ID**: Locate the proteome ID in the search results that matches your organism. 

## 4. Obtaining the Taxonomy ID

To find the taxonomy ID for your species:

1. **Visit NCBI**: Go to the [NCBI Taxonomy Database](https://www.ncbi.nlm.nih.gov/taxonomy).

2. **Search for Your Species**: Enter the name of your species in the search bar.

3. **Find the Tax ID**: Locate the taxonomy ID associated with your species.
