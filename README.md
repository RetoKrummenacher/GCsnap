# GCsnap2.0 Desktop

This is the implementation of GCsnap2.0 as a desktop application.

GCsnap is a flexible Python-based tool that allows for the interactive comparison of the genomic contexts of protein-coding genes from any genome at any taxonomic level, integrating them with functional and structural information for any of the genes shown. 

![figure1](https://github.com/JoanaMPereira/GCsnap/blob/master/examples/Fig1.png)

By connecting the output to different protein databases, the user can navigate through the different genomic contexts from a simple interactive platform, facilitating the further analysis of the contexts found. 

GCsnap is not limited to a single input format, can preform batch jobs and accepts protein classification maps. 

All information is stored in detailed, human and machine-readable files, and customable publication-ready figures.


Thank you for using and showing interest on GCsnap!

## Dependencies

GCsnap2.0 Desktop was written in Python 3.11. It was tested on Python 3.11. It requires mostly core Python modules and some external packages: 
  - Biopython
  - Bokeh
  - Matplotlib
  - Networkx 
  - PaCMAP (requires Numba which requires NumPy)
  - Scikit-learn
  - Pandas
  - Rich
  - Urllib3 and Requests
  - Jinja2

For detailed requirements including working versions, check ```pyproject.toml```.

Additionally, GCsnap relies on a local installation of MMseqs or for Windows users, the static binary. See below for installation details.
## Installation

### Installing from Source

Clone the repository with git, create a Conda environment and install:

**Linux/MacOS**
```
# To download
git clone https://github.com/RetoKrummenacher/GCsnap
cd GCsnap

# Change to new developpment branch 
git checkout gcsnap2desktop

# To install
conda create -n GCsnap -c conda-forge -c bioconda gcc=14.1 python=3.11 mmseqs2
conda activate GCsnap
pip install .
```

**Windows**  
There is no MMseqs2 installation candidate, just create Conda environment without mmseqs2 and get the static binary executable from here: https://mmseqs.com/latest/. Download, extract and store locally.
When running GCsnap, pass the path to the executable (i.e., mmseqs.bat) via the ```--mmseqs-executable-path``` argument.
```
# To download
git clone https://github.com/RetoKrummenacher/GCsnap
cd GCsnap

# Change to new developpment branch 
git checkout gcsnap2desktop

# To install
conda create -n GCsnap -c conda-forge -c bioconda gcc=14.1 python=3.11
conda activate GCsnap
pip install .
```


## Allowed inputs

GCsnap takes as main input a list of sequence identifiers, which can be in **Entrez, UniprotKB, UniRef, GeneID, and ENSEMBLE ID formats, or a mix**. These identifiers can be given as:
  - a text file, where each is in a different line
  - a fasta file, where the sequence header starts with the sequence identifier
  - a sequences cluster file in CLANS format
  - direct input in the terminal as a space-separated list
  
## Usage

In its most simple mode of usage, GCsnap only requires a list of sequence identifiers. 

**Required** arguments are:
```
  --targets: which can be a list of sequence identifiers, a text file, a fasta file, a clans file or a list of files
```
**Optional** arguments allow for the tweaking of GCsnap behaviour. Default values for arguments are taken from the  ```config.yaml```. They can be changed there directly or pass via the CLI, e.g., --ncpu 2.
A list of all possible arguments and their current default value can be show in the terminal via:
```  
  GCsnap --help 
```

The most relevant arguments are:
```  
  --n-cpu: the number of cores of zour CPU to be used when processing.
  --n-flanking5: the number of flanking genes to be taken on the 5' side.
  --n-flanking3: the number of flanking genes to be taken on the 3' side.
  --ncbi-user-email: it may be required to access the NCBI databases. It is not used for anything else.
  --ncbi-api-key: the key for NCBI API, which allows for up to 10 queries per second to NCBI databases. Can be obtained after   obtaing an NCBI account.
  --get-taxonomy: set to false if no taxonomy is to be collected.
  --annotate-TM: set to true to annotate the presence of transmembrane segments and signal peptides.
  --annotation-TM-mode: the mode to use to collect transmembrane and signal peptide annotations (phobius, tmhmm or uniprot).
  --clans-pattern: a set of patterns in CLANS groups names that define different groups to be considered as independent jobs.
  --operon-cluster-advanced: set to true to have a more comprehensive analysis/summary of the genomic contexts found. Ideal for very large input sets.  
```

For Windows users:
```  
  --mmseqs-executable-path: path to MMseqs executable (i.e., mmseqs.bat) if not installed in Conda environment.
```
