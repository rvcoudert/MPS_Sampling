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

