# Gene Strand Information Script

## Environment Setup

Before running the script, it's essential to set up the correct environment using conda. Follow the steps below to create and activate a conda environment:

1. Open a terminal or command prompt.
2. Run the following command to create a new conda environment:
   ```
   conda create --name myenv python=3.7
   ```
3. Activate the conda environment by running:
   ```
   conda activate myenv
   ```
4. Install necessary packages using the following command:
   ```
   conda install biopython
   ```

## Script Explanation

This script is written in Python and utilizes the Biopython library to fetch gene information from the NCBI Gene database. Specifically, it retrieves the strand information of a given gene.

### Function: `get_gene_strand_info(gene)`

- **Input**: A string representing the gene name.
- **Output**: A dictionary containing the gene name and its corresponding strand information.

### How the Function Works

1. The function starts by searching the NCBI Gene database for the input gene name.
2. It then fetches detailed gene information in XML format.
3. The XML data is parsed to extract the strand information, located within the `<Seq-interval_strand>` tag.
4. The function returns the gene name along with the strand information ("plus" or "minus"). If the strand information is not available, it returns "Information not available."

### Execution

Run the script by executing the following command in the terminal:
```
python your_script_name.py
```

Note: Make sure the conda environment is activated, and you have set your email address in the `Entrez.email` field in the script.

