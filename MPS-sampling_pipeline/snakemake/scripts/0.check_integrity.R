# Copyright or © or Copr.
# Rémi-Vinh COUDERT, Céline BROCHIER-ARMANET (2023-02-14)
# 
# rv.coudert@gmail.com
# celine.brochier-armanet@univ-lyon1.fr
# 
# This software is a computer program whose purpose is to [describe
# functionalities and technical features of your software].
# 
# This software is governed by the [CeCILL|CeCILL-B|CeCILL-C] license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the [CeCILL|CeCILL-B|CeCILL-C]
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
# 
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
# 
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
# 
# The fact that you are presently reading this means that you have had
# knowledge of the [CeCILL|CeCILL-B|CeCILL-C] license and that you accept its terms.

# 0.check_integrity.R
# Check data integrity and compute basic stats.


# save.image("~/0.check_integrity.RDS")
# stop()
# load("~/0.check_integrity.RDS")


##### Input #####


library(magrittr)
print("sessionInfo()")
print(sessionInfo())
print(".libPaths()")
print(.libPaths())
print("loadedNamespeces()")
print(loadedNamespaces())
cat("\n")


cores <- snakemake@params[["cores"]]
cat("Parallelization on", cores, "cores.\n")

if (cores > 1) {
  doParallel::registerDoParallel(cores = cores)
  parallel <- TRUE
} else {
  parallel <- FALSE
}


# Get the list of the fasta files.
fastaFiles_list <- snakemake@input[["fasta_files"]]
# Load the genome index.
genome_index <- snakemake@input[["genome_index"]] %>%
  data.table::fread()


##### Load fasta contents #####


# Remove extension from file name.
names(fastaFiles_list) <- fastaFiles_list %>%
  as.character() %>%
  basename() %>%
  stringr::str_remove(".fasta$")
# Open fasta files.
cat("Open", length(fastaFiles_list), "fasta files.")
cat("\n")
fastaContents_list <- fastaFiles_list %>%
  plyr::llply(Biostrings::readAAStringSet,
              .progress = "text",
              .parallel = parallel)
names(fastaContents_list) <- names(fastaFiles_list)


##### Check duplicates #####


# Check dupliated sequence names into fasta files.
nDuplicates_perFile <- fastaContents_list %>%
  plyr::llply(function(fastaContent) {
    fastaContent %>%
      names() %>%
      duplicated() %>%
      sum()
  }) %>% unlist()

if (sum(nDuplicates_perFile) > 0)
  stop("Duplications in the sequence names in the following fasta files:\n",
       nDuplicates_perFile %>%
         .[. > 0] %>%
         names() %>%
         paste(sep = " ", collapse = " "))

# Check duplicated genome names.
genomeList_index <- genome_index$genomeAccession %>%
  sort()

if (sum(duplicated(genomeList_index)) > 0)
  stop("Some genome names are duplicated.")


##### Input stats #####


# From each fasta content, extract the names of all protein sequences.
fastaContents_sequenceNames <- fastaContents_list %>%
  plyr::llply(names)
# Tidy all sequence names into a list.
genomeList_fasta <- fastaContents_sequenceNames %>%
  # Remove the protein names to avoid merging them with the assembly accession.
  unname() %>%
  # Tidy the list of lists into a simple list.
  unlist() %>%
  # Remove duplicates.
  unique() %>%
  # Order  it.
  sort()


##### Check integrity #####


## Check genome accessions are equal in genome index and fasta files.
# Get genome list from the genome index.
# genomeList_index <- genome_index$genomeAccession %>%
#   sort()
# Check genome list integrity between genome index and fasta files.
identical_genomeLists <- identical(genomeList_index,
                                   genomeList_fasta)
# If the two lists are not identical, then stop.
if (!identical_genomeLists)
  stop("Genome list are different from genome index and fasta files.")


##### Input stats #####


nFasta <- length(fastaFiles_list)
nGenomes <- nrow(genome_index)

input_stats <- data.table::data.table(
  feature = c("nFasta", "nGenomes"),
  value = c(nFasta, nGenomes))


##### Output #####


data.table::fwrite(
  x = input_stats,
  file = snakemake@output[["input_stats"]],
  sep = ";")

