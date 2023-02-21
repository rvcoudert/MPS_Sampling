# 🧬 MPS-Sampling

> An easy-to-use workflow for the selection of representative genomes from large genome databases.

## ⚙️ Installation

### 🐍 Install SnakeMake

SnakeMake installation method is available on the [official website](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### 📇 Import MPS-Sampling

Simply import the code from the GitHub repository: `git clone https://github.com/rvcoudert/MPS_Sampling`

### 🎉 Finish

Congratulations, MPS-Sampling is ready to go !

## 💾 Data

### 👉 Pedagogical example

The directory [`example/article`](https://github.com/rvcoudert/MPS_Sampling/tree/main/example/article) proposes the same example as presented in the article.
<br />
The sequences from 4 protein families enables the reduction of 10 genomes to 5 genomes.

### 🎥 Reproducibility

The data about the 158,027 genomes analyzed in the publication are available [online](https://figshare.com/articles/dataset/MPS-Sampling/22133120) on [figshare](https://figshare.com/).

## 🔧 Usage

### 1) Prepare the input

See the [input of the example](https://github.com/rvcoudert/MPS_Sampling/tree/main/example/article/input).
- *genome_index.csv*: A `csv` file with one row per genome, one column for the primary key "genomeAccession" and one column for each priority tag.
- *fasta/protein_family/protein_family.fasta*: A `fasta` file for each protein family.

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
                            min_nb_Linclusters=2 (optional)
                            deltas=0.5 (optional)]
                       
 Options:
  -s, --snakefile                   Path to the SnakeFile of MPS-Sampling.
  -d, --directory                   Path to the dataset to process.
  --config  input_directory=        Set the input directory. (default: "input")
            LinclustTemp_directory= Set the temporary directory for Linclust files. (default: "LinclustTemp")
            output_directory=       Set the output directory. (default: "output")
            e_values=               Set the e-value for Linclust (step 1). It can be a set, e.g. {0.0001, 0.0001}. (default: 0.00001)
            cov_modes=              Set the coverage mode for Linclust (step 1). It can be a set, e.g. {0, 1, 2}. (default: 0)
            min_covs=               Set the minimum coverage for Linclust (step 1). It can be a set, e.g. {0.7, 0.8, 0.9}. (default: 0.8)
            min_seq_ids=            Set the minimum sequence identity for Linclust (step 1). It can be a set, e.g. {0.4, 0.6, 0.8}. (default: 0.6)
            min_nb_Linclusters=     Set the minimum number of common Lin-clusters identity for pre-connection (step 3). It can be a set, e.g. {10, 20, 30}. (default: 2)
            deltas=                 Set the minimum similarity Δ for the hierarchical clustering with complete-linkage (step 4). It can be a set, e.g. {0.4, 0.5, 0.6}. (default: 0.5)
    
```

## 🐍 Workflow


### MPS-Sampling


0. Check integrity and compute input stats.
1. Build Lin-clusters whithin protein families.
2. Build the Lin-combination matrix.
    1. Build Lin-cluster matrix.
    2. Rereplicate the rows.
    3. Build the Lin-combination matrix.
3. Build pre-connected components.
4. Build MPS-clusters.
    1. Compute similarity submatrices.
    2. Build MPS-clusters.
5. Select MPS-representatives.

### 🎨 Illustrations

#### 💡 Complete workflow

Below the complete workflow of MPS-Sampling.
- First line: Genome level.
- Second line: Encoding levels.

<p align="center"><img src="https://github.com/rvcoudert/MPS_Sampling/blob/main/Illustrations/Fig.3_workflow_complete.png" alt="Complete WorkFlow of MPS-Sampling" width=80% /></p>


#### 📈  DAG

Below the directed acyclic graph (DAG) describing the different steps of the SnakeMake pipeline of MPS-Sampling.

<p align="center"><img src="https://github.com/rvcoudert/MPS_Sampling/blob/main/Illustrations/Fig.S1_Flowchart.png" alt="SnakeMake DAG (directed acyclic graph)" width=80%/></p>


#### 📋 ERD

Below the entity relationship diagram (ERD) describing the different tables processed by MPS-Sampling.

<p align="center"><img src="https://github.com/rvcoudert/MPS_Sampling/blob/main/Illustrations/Fig.S2_Tables.png" alt="Tables and ERD (entity relationship diagram)" width=60%/></p>


## 📖 Publication


[Coudert RV, Charrier JP, Jauffrit F, Flandrois JP, Brochier-Armanet C. Multi Proteins Similarity based sampling to select evolutionary significant representative genomes from large databases. BMC Bioinformatics, doi: XXX (2023)](https://github.com/rvcoudert/MPS_Sampling).

## 📲 Contact


Feel free to contact me for any comment, concern or discussion :
<br/>
[Mail](rv.coudert@gmail.com)
<br/>
[LinkedIn](https://www.linkedin.com/in/rvcoudert/)
