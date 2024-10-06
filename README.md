# FastTarget Setup Instructions

## Configuration Setup

### Getting Started

To begin, make sure you have the latest version of the repository. You can do this by either downloading the repository or running a `git pull` if you have already cloned it. Then, navigate to the repository folder:

```bash
git clone git@github.com:mcpalumbo/fasttarget.git
cd fasttarget
```

### Conda Environment Setup and Docker Images

To set up the conda environment for this project, follow these steps:

1. **Create the conda environment:**
   ```bash
   conda env create -f requirements.yml
   ```

After successfully creating the environment, activate it:

2. **Activate the conda environment:**
   ```bash
   conda activate fasttarget
   ```

This project requires the use of Docker images to run specific bioinformatics tools.

To download these images, execute the following command:

3. **Set up Docker images:**
   ```bash
   bash setup_docker.sh
   ```

**Note:** This pipeline is designed to run Docker commands without `sudo`.
To set up Docker so that `sudo` is not required, follow these steps:

1. **Create the Docker group** (if it doesn't already exist):

   ```bash
   sudo groupadd docker
   ```

2. **Add your user to the Docker group**:
   
   Replace `your-username` with your actual username:

   ```bash
   sudo usermod -aG docker your-username
   ```

3. **Log out and log back in** for the changes to take effect.

After completing these steps, you will be able to run Docker commands without needing to prepend `sudo`.

##### Warning ⚠️

If you have Docker Desktop installed, this solution may not work and could result in errors, certain configurations with Docker Desktop can lead to issues.

## Databases

This repository downloads and uses several key databases for analysis:

1. **Human Proteome from UniProt (UP000005640)**: The complete *Homo sapiens* proteome will be downloaded from UniProt, including all annotated proteins for the organism.

2. **Database of Essential Genes (DEG)**: This database contains essential bacterial genes and can be obtained from [DEG](http://origin.tubic.org/deg/public/index.php/download).

3. **Unified Human Gastrointestinal Protein (UHGP) catalog (UHGP-90)**: The [UHGP](https://doi.org/10.1038/s41587-020-0603-3) catalog clustered at 90% amino acid sequence identity will be downloaded.  **Note:** This file is over 10 GB in size, so ensure you have enough disk space available.

During pipeline execution, if structural data calculations are enabled, proteins from **PDB** and **AlphaFold** will also be downloaded, which may take up several GB depending on the size of the organism being analyzed.


## Test Case

To test the pipeline, we provide the small genome of *Mycoplasma pneumoniae*.
This dataset includes the GenBank (GBK) file, the proteome ID for structural data, and the SBML file for metabolism analysis. 
You can find the test dataset in the `organism/test` folder.

1. **Run the pipeline with the test dataset:**
   ```bash
   python tests.py
   ```

      #### Example Files

      In the `organism/test` folder, you can find example files that are used in this pipeline. These files include the GenBank (GBK) file, the SBML file for metabolism analysis, and others. These examples can serve as a reference for setting up your own dataset.


## Editing the `config.yml` File

The `config.yml` file is the central configuration file for this repository. It allows you to specify various settings related to your organism, CPU usage, structural data, metabolism, core genome analysis, offtarget analysis, DEG (Database of Essential Genes) analysis, localization, and metadata. Follow these steps to correctly fill out this file:

1. **Organism Information:**
   - `organism.name`: Enter a short alias/name for your genome without spaces. Example: `PAO`.
   - `organism.tax_id`: Fill in the NCBI Taxonomy ID of your organism's species. Example: 287 for Pseudomonas aeruginosa. This is used to obtain the genomes for calculating the core genome.
   - `organism.gbk_file`: Provide the path to the GenBank (GBK) file of your genome.

2. **CPU Usage Preferences:**
   - `cpus`: Specify the number of CPUs to be used for running this pipeline. Set to `None` if not specified.

3. **Structures:**
   - `structures`: Set to `True` if structural data is to be used; otherwise, set to `False`.
   - `proteome_uniprot`: If `structures` is enabled, provide the proteome ID from UniProt.

4. **Core Genome Analysis:**
   - `core`: Set to `True` if core genome analysis is required; otherwise, set to `False`.
   - `roary`: Enable Roary for core genome analysis by setting this to `True`.
   - `corecruncher`: Enable CoreCruncher by setting this to `True`.

5. **Metabolism:**
   - `metabolism`: Set to `True` if metabolic data is to be analyzed; otherwise, set to `False`.
   - Provide paths to your SBML file, chokepoint file, and smarttable file if metabolism is enabled.

   **Note:** To get help on how to generate these files, please read the tutorials in the `tutorial` folder.

6. **Offtarget Analysis:**
   - `offtarget.human`: Set to `True` to enable human offtarget analysis.
   - `offtarget.microbiome`: Set to `True` to enable microbiome offtarget analysis.
   - If microbiome offtarget analysis is enabled, specify identity and coverage filters.

7. **DEG Analysis:**
   - `deg`: Set to `True` to enable DEG analysis.
   - Specify identity and coverage filters if DEG analysis is enabled.

8. **Localization:**
   - `psortb`: Set to `True` to enable localization analysis. Specify the gram type (`n`, `p`, or `a`).

9. **Metadata:**
   - `metadata`: Set to `True` if you have additional metadata tables to include.
   - Provide a list of paths to your metadata tables under `meta_tables`.

### Validating and Viewing the Configuration

After editing the `config.yml` file, you can use the `configuration.py` script to validate and print the configuration. This will help you ensure that everything is loaded correctly.

#### Steps:

1. **Run the Configuration Script:**
   - In your terminal, navigate to the repository directory and run the following command:
     ```bash
     python configuration.py
     ```

2. **Validation:**
   - The script will first validate that all the required keys are present in the `config.yml` file. If any required key is missing, an error will be raised.

3. **View the Loaded Configuration:**
   - After validation, the script will print the loaded configuration details, allowing you to verify that everything is set up correctly.

## Running the FastTarget Pipeline

Once you have set up the configuration file, you can run the FastTarget pipeline using the `fasttarget.py` script.

### Steps:

1. **Run the FastTarget Script:**
   - In your terminal, navigate to the repository directory and run the following command:
     ```bash
     python fasttarget.py
     ```

   **Note:** If you want to run the pipeline with your own config file, you can use the following command:
   ```bash
   python fasttarget.py --config_file mi_config.yml
   ```


The results are stored in a file called `results_table.tsv` in the folder of your organism. This file contains the following columns:

- `gene`: Locus_tag ID.
- `betweenness_centrality`: Measure of the protein's centrality within the network.
- `edges`: Number of edges connected to the protein in the network.
- `chokepoints`: Indication of whether the protein is a chokepoint.
- `uniprot`: Identifier of the protein in the UniProt database.
- `druggability_score`: Druggability score assigned by FPocket to the best binding pocket.
- `pocket`: Best binding pocket of the protein.
- `structure`: Structure containing the best pocket (could be an AlphaFold model or a PDB crystal structure).
- `core_roary`: Indication of the protein's presence in the core genome according to Roary.
- `core_corecruncher`: Indication of the protein's presence in the core genome according to CoreCruncher.
- `human_offtarget`: BLASTp hit in the human proteome.
- `gut_microbiome_offtarget`: This score is calculated based on the number of hits in the microbiota and is normalized. A value of 1 indicates that there are at least 1000 BLASTp hits.
- `hit_in_deg`: BLASTp hit in the Database of Essential Genes (DEG).
- `psortb_localization`: Predicted localization of the protein based on PSORTb analysis.
- `metadata`: Each metadata table should have its own column.


## Visualization on Target Pathogen


The data generated in this pipeline can be loaded as metadata in Target Pathogen. 

You can find each metadata table to upload to Target Pathogen in `tables_for_TP` inside the folder of your organism.

This is a web interface that allows for the integration of multi-omics data to identify attractive targets in pathogens. 
One of its main features is the visualization of genome structures along with their druggable pockets, and the ability to customize filters and scoring functions to prioritize targets. 

It is available on Docker to generate this interface locally at:
https://github.com/sndg-arg/targetpathogenweb

Visit the Target Pathogen website at [http://target.sbg.qb.fcen.uba.ar/patho/](http://target.sbg.qb.fcen.uba.ar/patho/).
