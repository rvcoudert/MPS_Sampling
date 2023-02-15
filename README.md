# MPS-Sampling

> An easy-to-use workflow for the selection of representative genomes from large genome databases.

## ⚙️ Installation

### 1) Install SnakeMake

SnakeMake installation method is available on the [official website](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). It needs to be installed through [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) in a bash environment.

### 2) Set the Conda Environment

Install MMseq2 in the [Conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) dedicated to MPS-sampling Pipeline.

## 🔧 Usage

### 1) Prepare the Input

See the [inputs of the example](https://github.com/rvcoudert/MPS_Sampling/tree/main/data/example/input).
- *genome_index.csv*: A `csv` file with one row per genome, one column for the primary key "genomeAccession" and one column for each priority tag.
- *fasta/protein_family/protein_family.fasta*: A `fasta` file for each protein family.

### 2) Launch

`snakemake --cores 1 -s ~/GIT/MPS-sampling_pipeline/snakemake/Snakefile -d /zfs/RVC/RUN/MPS-sampling/article/ --latency-wait 60000  --config verbose=TRUE`


## 📖 Publication


[Coudert RV, Charrier JP, Jauffrit F, Flandrois JP, Brochier-Armanet C. Multi Proteins Similarity based sampling to select evolutionary significant representative genomes from large databases. BMC Bioinformatics, doi: XXX (2023)](https://github.com/rvcoudert/MPS_Sampling).

## 📲 Contact


Feel free to contact me for any comment, concern or discussion : 
[Mail](rv.coudert@gmail.com)
[LinkedIn](https://www.linkedin.com/in/rvcoudert/)
