# 3_1a.compute_dissimilarity_matrix.R
# Compute dissimilarity matrix.


# save.image("~/3_1a.compute_dissimilarity_matrix.RDS")
# stop()
# load("~/3_1a.compute_dissimilarity_matrix.RDS")


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


Lincombination_matrix <- snakemake@input[["Lincombination_matrix"]] %>%
  data.table::fread()


##### Compute #####


## Count the number of Linclusters per Lincombination.
print("__Computing a mask of each Lincluster__")
# Count the number of missing values per row.
nb_Linclusters_per_Lincombination <- Lincombination_matrix %>%
  is.na() %>% rowSums() %>% `*`(-1) %>% `+`(ncol(Lincombination_matrix))


# For each protein family, generate a mask for each Lincluster.
list_list_masks <- colnames(Lincombination_matrix) %>%
  # For a given protein family, generate a mask for each involved Lincluster.
  plyr::llply(function(column_name) {
    # DEBUG #
    # column_name <- "bl34"
    # Get the content of the focused protein family.
    column_content <- Lincombination_matrix[, get(column_name)]
    # For each Lincluster, generate the associated mask.
    # (Column-list.)
    list_masks <- (1:max(column_content, na.rm = TRUE)) %>%
      plyr::llply(function(Linclust_ID) {column_content == Linclust_ID})
    # Name this list.
    names(list_masks) <- 1:max(column_content, na.rm = TRUE)
    # Return this list.
    return(list_masks)
  })
# Name the generated list of lists.
names(list_list_masks) <- colnames(Lincombination_matrix)


print(Sys.time())
list_list_masks %>% object.size() %>% print(units = "Mb")
gc(reset = TRUE)


print("__Computing pairwise dissimilarities__")
# For each Lincombination, compute the pairwise number of common Linclusters
#   against all other Lincombination.
compute_dissimilarity <- function(Lincombination_ID) {
  # DEBUG #
  # Lincombination_ID <- 5

  ## Prepare the Dice denominator.
  dice_denominator <- nb_Linclusters_per_Lincombination +
    nb_Linclusters_per_Lincombination[[Lincombination_ID]]

  ## Prepare the Dice nominator.
  # Get the content of the focused Lincombination.
  Lincombination_content <- Lincombination_matrix[Lincombination_ID, ] %>%
    as.list()
  # Get the list of masks of the focused Lincombination.
  # (Row-list.)
  list_masks <- Lincombination_content %>%
    purrr::imap(function(Linclust_ID, column_name)
      list_list_masks[[column_name]][[Linclust_ID]])
  # Sum the all the masks of the focused Lincombination,
  dice_numerator <- do.call(rbind, list_masks) %>%
    colSums(na.rm = TRUE) %>%
    # multiply by two to obtain the Dice numerator,
    magrittr::multiply_by(2)
  # Get the Dice dissimilarity.
  1 - dice_numerator / dice_denominator
}

dissimilarity_matrix <- Lincombination_matrix[
  , as.list(compute_dissimilarity(.I)),
  by = 1:nrow(Lincombination_matrix)]

rm(list_list_masks)
rm(Lincombination_matrix)

dissimilarity_matrix <- as.data.frame(dissimilarity_matrix)

dissimilarity_matrix$nrow <- NULL

# dissimilarity_matrix <- as.dist(dissimilarity_matrix)


##### Output #####


print("__Saving the dissimilarity matrix__")
saveRDS(object = dissimilarity_matrix,
        file = snakemake@output[["dissimilarity_matrix_raw"]])



