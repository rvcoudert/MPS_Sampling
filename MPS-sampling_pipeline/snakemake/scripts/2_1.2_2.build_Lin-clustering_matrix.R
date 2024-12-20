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
