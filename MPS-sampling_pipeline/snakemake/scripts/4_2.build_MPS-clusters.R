# 4_2.build_MPS-clusters.R
# Build MPS-clusters through hierarchical clustering with complete-linkage.


# save.image("~/4_2.build_MPS-clusters.RDS")
# stop()
# load("~/4_2.build_MPS-clusters.RDS")


##### Input #####


library(magrittr)
print("sessionInfo()")
print(sessionInfo())
print(".libPaths()")
print(.libPaths())
print("loadedNamespeces()")
print(loadedNamespaces())
cat("\n")

Lincomb_precomp_links <- snakemake@input[["Lincomb_precomp_links"]] %>%
  data.table::fread()
similarity_submatrices_directory <- snakemake@params[["similarity_submatrices_directory"]]
delta <- snakemake@params[["delta"]] %>%
  as.numeric()


##### Lincomb_MPS_cluster_links ######


# Count the number of Lin-ombinations per pre-connected component.
nb_Lincombs_per_preconnected_component <- Lincomb_precomp_links %>%
  .[, .(nb_Lincombinations = .N), preconnected_component_ID] %>%
  .[order(preconnected_component_ID)]

# Build MPS-clusters for each pre-connected component.
MPS_clusterings <- nb_Lincombs_per_preconnected_component %>%
  # Remove pre-connected component with only one Lin-combination.
  .[nb_Lincombinations > 1, preconnected_component_ID] %>%
  # Remove the NA pre-connected component,
  #   ie the Lin-combinations with no pre-connected component.
  .[!is.na(.)] %>%
  # Name the list for <plyr::ldply>.
  magrittr::set_names(., .) %>%
  plyr::ldply(function(selected_precomponent_ID) {
    # DEBUG #
    # selected_precomponent_ID <- "1"
    # Display a progression message.
    cat(paste0(
      "Compute the hierarchical clustering with complete-linkage",
      " for the pre-connected component : ",
      selected_precomponent_ID,
      " (size : ",
      nb_Lincombs_per_preconnected_component %>%
        .[preconnected_component_ID == selected_precomponent_ID, nb_Lincombinations] %>%
        format(big.mark = " "),
      ").",
      "\n"))
    # Get the file name of the associated similarity matrix.
    filename <- paste0(similarity_submatrices_directory,
                       selected_precomponent_ID,
                       ".RDS")
    # Load the similarity submatrix of the associated pre-connected component.
    similarity_submatrix <- readRDS(file = filename)
    # Compute hierarchical clustering.
    my_hclust <- hclust(d = as.dist(1 - similarity_submatrix),
                        method = "complete")
    my_cutree <- cutree(tree = my_hclust,
                        h = 1 - delta)
    # Tidy the hierarchical clustering.
    MPS_clustering <- data.table::data.table(
      Lincombination_ID = as.integer(names(my_cutree)),
      MPS_cluster_ID = paste0(selected_precomponent_ID,
                              "_",
                              my_cutree))
    # Return the expected MPS-cluster.
    return(MPS_clustering)
  }, .id = "preconnected_component_ID")
# If there is no hierarchical clustering.
if (nrow(MPS_clusterings) == 0) {
  Lincomb_MPS_cluster_links <- Lincomb_precomp_links[, .(Lincombination_ID,
                                                         MPS_cluster_ID = NA)]
} else {
  # Convert into a <data.table>.
  MPS_clusterings <- MPS_clusterings %>%
    data.table::as.data.table()
  # Re-order the columns.
  Lincomb_MPS_cluster_links <- MPS_clusterings[
    , .(Lincombination_ID, MPS_cluster_ID)]
  # Add missing combination_ID.
  Lincomb_MPS_cluster_links <- data.table::merge.data.table(
    x = Lincomb_precomp_links[, .(Lincombination_ID)],
    y = Lincomb_MPS_cluster_links,
    by = "Lincombination_ID",
    all.x = TRUE,
    order = TRUE)
}
# Order the rows.
Lincomb_MPS_cluster_links <- Lincomb_MPS_cluster_links[
  order(Lincombination_ID)]


##### Output #####


data.table::fwrite(
  x = Lincomb_MPS_cluster_links,
  file = snakemake@output[["Lincomb_MPS_cluster_links"]],
  sep = ";")

