# FastTarget Setup Instructions

## Configuration Setup


### Getting Started

To begin, ensure that you have the latest version of the repository. You can do this by either downloading the repository or running a `git pull` if you've already cloned it. Then, navigate into the repository folder:

```bash
git clone git@github.com:mcpalumbo/fasttarget.git
cd fasttarget
```

### Conda Environment Setup and Docker images

To set up the conda environment for this project, follow the steps below.

Use the provided environment.yml file to create an environment with all the required dependencies:

1. **Create the conda environment:**
    ```bash
    conda env create -f requirements.yml
    ```

After the environment is successfully created, activate it:

2. **Activate the conda environment:**
    ```bash
    conda activate fasttarget
    ```

This project requires the use of Docker images to run specific bioinformatics tools. To download these images, execute the following command:

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

After completing these steps, you'll be able to run Docker commands without needing to prepend `sudo`.

### Test case

To test the pipeline, we provide the small genome of *Mycoplasma pneumoniae*.
This dataset includes the GenBank (GBK) file, the proteome ID for structural data, and the SBML file for metabolism analysis. 
You can find the test dataset in the `organism/test` folder.

To run the pipeline with the test dataset:

    ```bash
    python tests.py
    ```

### Editing the `config.yml` File

The `config.yml` file is the central configuration file for this repository. It allows you to specify various settings related to your organism, CPU usage, structural data, metabolism, core genome analysis, offtarget analysis, DEG (Database of Essential Genes) analysis, localization, and metadata. Below are the steps to correctly fill out this file:

1. **Organism Information:**
   - `organism.name`: Enter a short alias/name for your genome without spaces. Example: `PAO`.
   - `organism.tax_id`: Fill with the number of NCBI Taxonomy ID of the specie of your organism. Example: 287 for Pseudomonas aeruginosa. It is used to obtain the genomes to calculate the core genome.
   - `organism.gbk_file`: Provide the path to the GenBank (GBK) file of your genome.

2. **CPU Usage Preferences:**
   - `cpus`: Specify the number of CPUs to be used for running this pipeline. Set to `None` if not specified.

3. **Structures:**
   - `structures`: Set to `True` if structural data is to be used; otherwise, set to `False`.
   - `proteome_uniprot`: If `structures` is enabled, provide the proteome ID from UniProt.

4. **Core Genome Analysis:**
   - `core`: Set to `True` if core genome analysis is required; otherwise, set to `False`.
   - `roary`: Enable Roary for core genome analysis by setting this to `True`.
   - `corecruncher`: Enable CoreCruncher by setting this to `True`. Provide the path to the CoreCruncher script under `corecruncher_script`.

5. **Metabolism:**
   - `metabolism`: Set to `True` if metabolic data is to be analyzed; otherwise, set to `False`.
   - Provide paths to your SBML file, chokepoint file, and smarttable file if metabolism is enabled.

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
   - In your terminal, in the repo directory you can find the file `configuration.py`. Run with the following command:
     ```bash
     python configuration.py
     ```

2. **Validation:**
   - The script will first validate that all the required keys are present in the `config.yml` file. If any required key is missing, an error will be raised.

3. **View the Loaded Configuration:**
   - After validation, the script will print the loaded configuration details, allowing you to verify that everything is set up correctly.


### Running the FastTarget Pipeline

Once you have set up the configuration file, you can run the FastTarget pipeline using the `fasttarget.py` script.

#### Steps:

1. **Run the FastTarget Script:**
   - In your terminal, run the main script with following command:
     ```bash
     python fasttarget.py
     ```

