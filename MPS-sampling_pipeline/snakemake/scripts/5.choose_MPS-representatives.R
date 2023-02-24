# 5.choose_MPS-representatives.R
# Select one MPS-representative genome per MPS-cluster.


# save.image("~/5.choose_MPS-representatives.RDS")
# stop()
# load("~/5.choose_MPS-representatives.RDS")


##### Input #####


library(magrittr)
print("sessionInfo()")
print(sessionInfo())
print(".libPaths()")
print(.libPaths())
print("loadedNamespeces()")
print(loadedNamespaces())
cat("\n")

genome_index <- snakemake@input[["genome_index"]] %>%
  data.table::fread()
genome_cascading_links <- snakemake@input[["genome_cascading_links"]] %>%
  data.table::fread()
Lincombination_matrix <- snakemake@input[["Lincombination_matrix"]] %>%
  data.table::fread()


##### Pre-processing #####


##### __Functions #####


# Replace with relative frequency
#
#
# From a Lin-combination matrix,
#   replace each Lin-custer with its (vertical) relative frequency.
replace_with_relativeFrequency <- function(Lincombination_matrix) {
  # List Lin-combination columns.
  Lincombination_columns <- colnames(Lincombination_matrix)

  # Replace each Lin-custer by its (vertical) relative requency.
  Lincluster_frequency_DT <- Lincombination_matrix[, lapply(
    .SD,
    function(Lincombination_list) {
      # DEBUG #
      # Lincombination_list <- Lincombination_matrix$uL1
      # Compute the frequencies of each different Lin-cluster.
      Lincluster_table <- as.data.frame(table(
        Lincombination_list))
      # Compute the relative frequency for each Lin-cluster of the list.
      Lincluster_table$Freq[match(Lincombination_list,
                                  Lincluster_table[, 1])] /
        #   divided by the length of the list.
        length(Lincombination_list)
    }),
    , .SDcols = Lincombination_columns]

  # Export the result.
  return(Lincluster_frequency_DT)
}


# Compute relative frequency per Lin-combination
#
#
# From a Lin-combination matrix,
#   compute the relative frequency for each Lin-combination
#   of their Lin-clusters.
compute_relativeFrequency_perLincombination <- function(Lincombination_matrix) {
  # Replace each Lin-cluster with its relative frequency.
  Lincluster_frequency_DT <- replace_with_relativeFrequency(
    Lincombination_matrix = Lincombination_matrix)

  # Sum all the relative frequencies per Lin-combination.
  sum_relFrequency <- rowSums(Lincluster_frequency_DT, na.rm = TRUE)

  # Export the result.
  return(sum_relFrequency)
}


##### __Genomc full index #####


list_tags <- colnames(genome_index) %>%
  .[!. %in% "genomeAccession"]
genome_fullIndex <- data.table::merge.data.table(
  x = genome_index,
  y = genome_cascading_links,
  by = "genomeAccession",
  all = TRUE,
  sort = TRUE)


##### __Pseudo-randomize #####


# Compute genome hash from genome accession.
genome_fullIndex[
  , genomeHash := genomeAccession %>%
    plyr::llply(rlang::hash) %>% stringr::str_sub(-6) %>% strtoi(16L)]


##### __Completeness #####


# Count number of Lin-clusters per Lin-combination.
nLinclusters_perLincombination <- Lincombination_matrix[
  , .(Lincombination_ID = .I,
      nLinclusters = ncol(Lincombination_matrix) -
        Reduce("+", lapply(.SD, is.na)))]
# Add this information to the genome full index.
genome_fullIndex <- data.table::merge.data.table(
  x = genome_fullIndex,
  y = nLinclusters_perLincombination,
  by = "Lincombination_ID",
  all = TRUE)


##### MPS-rep #####


# Prepare the Lin-combination links.
Lincomb_links <- genome_cascading_links %>%
  .[, .(Lincombination_ID, preconnected_component_ID, MPS_cluster_ID) ] %>%
  unique()
# Count the number of genomes per MPS_cluster.
nGenomes_perMPS_cluster <- genome_cascading_links %>%
  .[, .(nGenomes = .N), by = MPS_cluster_ID]
# List the MPS_clusters with more than one associated genomes.
rich_MPS_clusters <- nGenomes_perMPS_cluster[
  nGenomes > 1, MPS_cluster_ID] %>%
  magrittr::set_names(., .)
# Count number of MPS_clusters to process.
nb_rich_MPS_clusters <- length(rich_MPS_clusters)
cat("Choose MPS-reprentatives for",
    format(x = nb_rich_MPS_clusters, big.mark = " "),
    "rich MPS_clusters.\n")
# Select MPS-representatives for MPS-cluster with at least several genomes.
MPS_representatives <- rich_MPS_clusters %>%
  plyr::ldply(function(current_MPS_cluster) {
    # DEBUG #
    # current_MPS_cluster <- "1_1"

    # Get the list of the associated Lin-combinations.
    list_Lincombinations <- Lincomb_links[MPS_cluster_ID == current_MPS_cluster,
                                          Lincombination_ID]

    # Get the Lin-combination submatrix.
    Lincombination_submatrix <- Lincombination_matrix %>%
      .[row.names(.) %in% list_Lincombinations]

    # Compute the cumulated relative frequencies.
    relativeFrequencies <- Lincombination_submatrix %>%
      compute_relativeFrequency_perLincombination()

    # Link the Lincombination_ID with the relative frequency.
    relativeFrequencies_DT <- data.table::data.table(
      Lincombination_ID = list_Lincombinations,
      relativeFrequency = relativeFrequencies)

    # Get associated genome sub-index.
    genome_subIndex <- genome_fullIndex[
      Lincombination_ID %in% list_Lincombinations]

    # Add the relative frequency.
    genome_subIndex <- data.table::merge.data.table(
      x = genome_subIndex,
      y = relativeFrequencies_DT,
      by = "Lincombination_ID")


    # 1) Choose according to priority tags.
    # 2) Choose according to data coverage.
    # 3) Choose according to centrality.
    # 4) Choose according to genome hash.
    genome_subIndex <- genome_subIndex[
      order(-nLinclusters,
            -relativeFrequency,
            genomeHash)]
    genome_selection <- genome_subIndex[0]

    # Try each priority tag until find at list one genome.
    i <- 1
    nTags <- length(list_tags)
    while (i <= nTags & nrow(genome_selection) == 0) {
      genome_selection <- genome_subIndex[get(list_tags[[i]]) == TRUE]
      i <- i + 1
    }
    rm(i)

    # If no genome has been selected, then select all the genomes.
    if (nrow(genome_selection) == 0)
      genome_selection <- genome_subIndex

    # Select the best genome for this tag.
    genome_selection <- genome_selection[1, genomeAccession]

    data.table::data.table(MPS_representative = genome_selection)
  }, .progress = "text", .id = "MPS_cluster_ID")
# Convert into a <data.table>.
MPS_representatives <- data.table::as.data.table(MPS_representatives)
# Prepare the simple Lin-combinations.
simple_MPS_clusters_DT <- genome_cascading_links[
  !MPS_cluster_ID %in% rich_MPS_clusters,
  .(MPS_cluster_ID,
    MPS_representative = genomeAccession)]
# Add the simple MPS-clusters.
#   If there is no rich MPS_cluster, then there are only simple MPS-clusters.
if (length(rich_MPS_clusters) == 0) {
  MPS_representatives <- simple_MPS_clusters_DT
} else {
  MPS_representatives <- rbind(MPS_representatives,
                               simple_MPS_clusters_DT)
}
MPS_representatives[, MPS_cluster_ID := as.character(MPS_cluster_ID)]
MPS_representatives <- MPS_representatives[order(MPS_cluster_ID)]
cat("\nDone.\n")


##### MPS-links #####


# Build MPS Links.
MPS_links <- data.table::merge.data.table(
  x = genome_cascading_links[, .(genomeAccession,
                                 MPS_cluster_ID)],
  y = MPS_representatives[, .(MPS_cluster_ID,
                              MPS_representative)],
  by = "MPS_cluster_ID")


##### Output #####


data.table::fwrite(
  x = MPS_representatives[, .(MPS_representative)],
  file = snakemake@output[["MPS_representatives"]],
  sep = ";")
data.table::fwrite(
  x = MPS_links[, .(genomeAccession, MPS_representative)],
  file = snakemake@output[["MPS_links"]],
  sep = ";")
