# 🧬 MPS-Sampling

> An easy-to-use workflow for the selection of representative genomes from large genome databases.

## ⚙️ Installation

### 💻 Install a Conda-based Python3 distribution

It could be :
- [Miniforge](https://github.com/conda-forge/miniforge), a new name project for Mambaforge (recommanded)
- or [Conda](https://docs.conda.io/) (not recommanded)

Miniforge can be installed for your own system with the following code :

```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

`$(uname -m)` and `$(uname -m)` will provide the required informations for choosing the adequate version of Miniforge, e.g. respectively `Linux` and `x86_64`.

### 🐍 Install Snakemake

[Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) can be installed with the following code :
```
conda activate base
mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal
```

### 📇 Import MPS-Sampling

Simply import the code from the GitHub repository: `git clone https://github.com/rvcoudert/MPS_Sampling`

### 🎉 Finish

Activate the Conda environment for Snakemake, preferably in a dedicated [screen](https://linux.die.net/man/1/screen).

```
screen -S MPS_Sampling
conda activate snakemake
```

Congratulations, MPS-Sampling is ready to go !

## 🔧 Usage

### 1) Prepare the input

The data directory needs to contain the input data in a subdirectory called `input` and containing :
- *genome_index.csv*: A `csv` file with one row per genome and two columns, one for the primary key "genomeAccession" and another one for the priority score.
- *fasta/protein_family/protein_family.fasta*: A `fasta` file for each protein family.

See the organization of the [example_1](https://github.com/rvcoudert/MPS_Sampling/tree/main/data/example_1).

### 2) Launch

```
snakemake --use-conda --cores [NUMBER_OF_CORES] -s [PATH_TO_SNAKEFILE] -d [PATH_TO_DATA]
```

```
Usage: snakemake  [--use-conda (necessary to use conda env)]
                  [--cores NUMBER OF CORES]
                  [-s|--snakefile PATH_TO_SNAKEFILE]
                  [-d|--directory PATH_TO_DATA]
                  [--config input_directory=input (optional)
                            LinclustTemp_directory=LinclustTemp (optional)
                            output_directory=output (optional)
                            e_values=0.00001 (optional)
                            cov_modes=0 (optional)
                            min_covs=0.8 (optional)
                            min_seq_ids=0.6 (optional)
                            deltas=0.5 (optional)]
                       
 Options:
  -s, --snakefile                   Path to the SnakeFile of MPS-Sampling, located at [MPS-sampling_pipeline/snakemake/Snakefile](https://github.com/rvcoudert/MPS_Sampling/blob/main/MPS-sampling_pipeline/snakemake/Snakefile).
  -d, --directory                   Path to the dataset to process.
  --config  input_directory=        Set the input directory. (default: "input")
            LinclustTemp_directory= Set the temporary directory for Linclust files. (default: "LinclustTemp")
            output_directory=       Set the output directory. (default: "output")
            e_values=               Set the e-value for Linclust (step 1). It can be a set, e.g. [0.0001,0.0001]. (default: 0.00001)
            cov_modes=              Set the coverage mode for Linclust (step 1). It can be a set, e.g. [0,1,2]. (default: 0)
            min_covs=               Set the minimum coverage for Linclust (step 1). It can be a set, e.g. [0.7,0.8,0.9]. (default: 0.8)
            min_seq_ids=            Set the minimum sequence identity for Linclust (step 1). It can be a set, e.g. [0.4,0.6,0.8]. (default: 0.6)
            deltas=                 Set the minimum similarity Δ for the hierarchical clustering with complete-linkage (step 4). It can be a set, e.g. [0.4,0.5,0.6]. (default: 0.5)
    
```

The examples can analyzed by MPS-Sampling with :

```
snakemake --use-conda --cores 1 -s MPS-sampling_pipeline/snakemake/Snakefile -d data/example_1
snakemake --use-conda --cores 1 -s MPS-sampling_pipeline/snakemake/Snakefile -d data/example_2
snakemake --use-conda --cores 1 -s MPS-sampling_pipeline/snakemake/Snakefile -d data/example_3
snakemake --use-conda --cores 1 -s MPS-sampling_pipeline/snakemake/Snakefile -d data/example_4
snakemake --use-conda --cores 1 -s MPS-sampling_pipeline/snakemake/Snakefile -d data/example_5
snakemake --use-conda --cores 1 -s MPS-sampling_pipeline/snakemake/Snakefile -d data/example_bio
```

## 💾 Data

### 👉 Examples

The directory [`data`](https://github.com/rvcoudert/MPS_Sampling/tree/main/data) proposes the same examples as presented in the article.
<br />
The sequences from 4 protein families enables the reduction of 10 genomes to 5 genomes.
<br />
There are also several versions of the same example, withthe same data but in another order or with another labels, as described in the [ReadMe.md](https://github.com/rvcoudert/MPS_Sampling/blob/main/data/ReadMe.md).


### 🎥 Reproducibility

The data about the 178,203 genomes analyzed in the publication are available [online on figshare](https://figshare.com/articles/dataset/Article_-_MPS-Sampling/24552160).

```
snakemake --use-conda --cores 1 -s [PATH_TO_SNAKEFILE] -d [PATH_TO_DATA] --config deltas=[1,0.9,0.8,0.7,0.6,0.5,0.4]
```

## 🐍 Workflow


### MPS-Sampling


0. Check integrity and compute input stats.
1. Build Lin-clusters whithin protein families.
2. Build the Lin-combinations.
    1. Build the Lin-clustering matrix.
    2. Dereplicate the rows.
    3. Build the Lin-combination matrix.
3. Build MPS-clusters.
    1. Compute the similarity matrix.
    2. Build MPS-clusters.
4. Select MPS-representatives.

### 🎨 Illustrations

#### 💡 Complete workflow

Below the complete workflow of MPS-Sampling.

<p align="center"><img src="https://github.com/rvcoudert/MPS_Sampling/blob/main/illustrations/1.workflow.png" alt="Complete WorkFlow of MPS-Sampling" width=80% /></p>


#### 📈  DAG

Below the directed acyclic graph (DAG) describing the different steps of the Snakemake pipeline of MPS-Sampling.

<p align="center"><img src="https://github.com/rvcoudert/MPS_Sampling/blob/main/illustrations/2.DAG.png" alt="Snakemake DAG (directed acyclic graph)" width=80%/></p>


#### 📋 ERD

Below the entity relationship diagram (ERD) describing the different tables processed by MPS-Sampling.

<p align="center"><img src="https://github.com/rvcoudert/MPS_Sampling/blob/main/illustrations/3.ERD_.png" alt="Tables and ERD (entity relationship diagram)" width=60%/></p>


## 📖 Publication


[(Coudert, R. V., Charrier, J. P., Jauffrit, F., Flandrois, J. P., & Brochier-Armanet, C. (2025). Multi-proteins similarity-based sampling to select representative genomes from large databases. BMC bioinformatics, 26(1), 121.)](https://doi.org/10.1186/s12859-025-06095-3)


## 📲 Contact


Feel free to contact me for any comment, concern or discussion :
<br/>
[Mail](rv.coudert@gmail.com)
<br/>
[LinkedIn](https://www.linkedin.com/in/rvcoudert/)
