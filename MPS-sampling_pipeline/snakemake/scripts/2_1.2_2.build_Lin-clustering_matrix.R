# 2_1.2_2.build_Lin-clustering_matrix.R
# Label and merge Lin-clusters within protein families to build a single Lin-clustering matrix.


# save.image("~/2_1.2_2.build_Lin-clustering_matrix.RDS")
# stop()
# load("~/2_1.2_2.build_Lin-clustering_matrix.RDS")


##### Input #####


library(magrittr)
print("sessionInfo()")
print(sessionInfo())
print(".libPaths()")
print(.libPaths())
print("loadedNamespeces()")
print(loadedNamespaces())
cat("\n")


ascending_order <- snakemake@params[["ascending_order"]]
cores <- snakemake@params[["cores"]]
cat("Parallelization on", cores, "cores.\n")


if (cores > 1) {
  doParallel::registerDoParallel(cores = cores)
  parallel <- TRUE
} else {
  parallel <- FALSE
}


# Get the list of the Lin-cluster files.
LinclusterFiles_List <- snakemake@input[["Lincluster_files"]]


##### Compute #####


## Compute the Lin-clustering matrix.
nFamilies <- length(LinclusterFiles_List)
cat(paste0("Merge Lin-clusters from ", nFamilies, " cluster files.\n"))
# Name the Lin-cluster files according to the captured protein family name.
names(LinclusterFiles_List) <- LinclusterFiles_List %>%
  as.character() %>%
  basename() %>%
  stringr::str_remove(".fasta_cluster.tsv$")
# Read all Lin-cluster files and number the Lin-clusters.
Lincluster_DTs <- LinclusterFiles_List %>%
  plyr::llply(function(LinclusterFile) {
    # DEBUG #
    # LinclusterFile <- LinclusterFiles_List[[1]]
    # Read each file, without header and using the tabulation as separator.
    Lincluster_DT <- LinclusterFile %>%
      data.table::fread(header = FALSE, sep = "\t")
    # Rename the <data.table> columns with "center" and "sequence" respectively.
    colnames(Lincluster_DT) <- c("centroid", "sequence")
    # Enumerate all Lin-clusters, by numbering each centroid sequence.
    Lincluster_DT[, centroid := centroid %>% factor() %>% as.integer()]
    # Compute the frequency of each Lin-cluster.
    clusterFrequencies <- Lincluster_DT$centroid %>%
      table() %>%
      sort(decreasing = TRUE)
    # Generate the new order.
    cluster_newOrder <- 1:length(clusterFrequencies) %>%
      magrittr::set_names(names(clusterFrequencies))
    # Assign the <Lincluster_ID> according to the most frequent Lin-clusters.
    Lincluster_DT$Lincluster_ID <- cluster_newOrder[
      match(Lincluster_DT$centroid, names(cluster_newOrder))] %>% as.integer()
    # Re-order according to the <Lincluster_ID>
    #   and rename the "sequence" column into "genomeAccession".
    Lincluster_DT <- Lincluster_DT[
      order(Lincluster_ID),
      .(genomeAccession = sequence, Lincluster_ID)]
    # Export the <data.table> of the Lin-cluster.
    return(Lincluster_DT)
  }, .parallel = parallel)
# Rename the second column according to the protein family name.
Lincluster_DTs <- Lincluster_DTs %>%
  purrr::imap(
    function(Lincluster_DT, proteinFamilyName) {
      colnames(Lincluster_DT)[[2]] <- proteinFamilyName
      return(Lincluster_DT)
    })
# Tidy all Lin-clusters into a wide <data.table>.
Linclustering_matrix <- Reduce(
  function(DT_1, DT_2)
    data.table::merge.data.table(DT_1,
                                 DT_2,
                                 by = "genomeAccession",
                                 all = TRUE),
  Lincluster_DTs)


##### Order #####


## (i) Order the Lin-clustering matrix according to hashing values.


# List the protein families and the genomes.
proteinFamily_list <- colnames(Linclustering_matrix) %>%
  .[2:(nFamilies + 1)]
genomeAccession_list <- Linclustering_matrix$genomeAccession

# Compute the hashing values for protein families and genomes.
proteinFamily_hashValues <- proteinFamily_list %>%
  plyr::llply(rlang::hash) %>%
  unlist()
genomeAccession_hashValues <- genomeAccession_list %>%
  plyr::llply(rlang::hash) %>%
  unlist()

# if (any(duplicated(genomeAccession_hashValues)))
#   stop("Overlapping hashing values from genome names.")
# if (any(duplicated(proteinFamily_hashValues)))
#   stop("Overlapping hashing values from protein family names.")

# Re-order protein families and genomes according to hashing values.
proteinFamily_match <- match(x = sort(proteinFamily_hashValues),
                             table = proteinFamily_hashValues)
genomeAccession_match <- match(x = sort(genomeAccession_hashValues),
                               table = genomeAccession_hashValues)

# Generate the new list of protein families and genomes.
proteinFamily_newList <- proteinFamily_list[proteinFamily_match]
genomeAccession_newList <- genomeAccession_list[genomeAccession_match]

# Order the Lin-clustering matrix according to hashing values.
Linclustering_matrix %>%
  data.table::setcolorder(c("genomeAccession", proteinFamily_newList))
Linclustering_matrix <- Linclustering_matrix[genomeAccession_match]


## (ii) Order the Lin-clustering matrix according to Lin-clustering.


# Order the columns of the <data.table> according to the number of Lin-clusters.
nLinclusters_perProteinFamily <- lapply(
  Linclustering_matrix[, -1],
  max,
  na.rm = TRUE) %>%
  unlist() %>%
  sort(decreasing = !ascending_order)
proteinFamily_newList_2 <- names(nLinclusters_perProteinFamily)
Linclustering_matrix %>%
  data.table::setcolorder(c("genomeAccession", proteinFamily_newList_2))


# Order the rows of the <data.table> according to the Lin-cluster numbering.
Linclustering_matrix <- Linclustering_matrix[
  do.call(order, list(cols = Linclustering_matrix[, -1],
                      decreasing = !ascending_order))]
# Order the <data.table> according to "genomeAccession".
# Linclustering_matrix <- Linclustering_matrix[order(genomeAccession), ]
cat("Done.\n")


##### Output #####


data.table::fwrite(
  x = Linclustering_matrix,
  file = snakemake@output[["Linclustering_matrix"]],
  sep = ";")
