# 4_1.compute_similarity_submatrices.R
# Compute similarity submatrix for each pre-connected component.


# save.image("~/4_1.compute_similarity_submatrices.R")
# stop()
# load("~/4_1.compute_similarity_submatrices.R")


##### Input #####


library(magrittr)

Lincombination_matrix <- snakemake@input[["Lincombination_matrix"]] %>%
  data.table::fread()
Lincomb_precomp_links <- snakemake@input[["Lincomb_precomp_links"]] %>%
  data.table::fread()
similarity_submatrices_directory <- snakemake@params[["similarity_submatrices_directory"]]


#### Similarity submatrices #####


# Add a first column with the Lin-combination ID.
Lincombination_matrix <- cbind(
  "Lincombination_ID" = 1:nrow(Lincombination_matrix),
  Lincombination_matrix)

# Count the number of Lin-ombinations per pre-connected component.
nb_Lincombs_per_preconnected_component <- Lincomb_precomp_links %>%
  .[, .(nb_Lincombinations = .N), preconnected_component_ID] %>%
  .[order(preconnected_component_ID)]

# If it not exists, then create the similarity matrices directory.
if (!dir.exists(similarity_submatrices_directory))
  dir.create(similarity_submatrices_directory, recursive = TRUE)

# Compute the similarity submatrices.
nb_Lincombs_per_preconnected_component %>%
  # Remove pre-connected component with only one Lin-combination.
  .[nb_Lincombinations > 1, preconnected_component_ID] %>%
  # Remove the NA pre-connected component,
  #   ie the Lin-combinations with no pre-connected component.
  .[!is.na(.)] %>%
  # Compute the similarity submatrix for each pre-connected component.
  plyr::l_ply(function(selected_precomponent_ID) {
    # DEBUG #
    # selected_precomponent_ID <- 1
    # Prepare the file name of the similarity submatrix.
    filename <- paste0(similarity_submatrices_directory,
                       selected_precomponent_ID,
                       ".RDS")
    # If the file of the similarity matrix does not exist,
    #   then compute and save it.
    if (!file.exists(filename)) {
      # Display a progression message.
      cat(paste0(
        "Compute the similarity matrix of the pre-connected component : ",
        selected_precomponent_ID,
        " (",
        nb_Lincombs_per_preconnected_component %>%
          .[preconnected_component_ID == selected_precomponent_ID,
            nb_Lincombinations],
        " Lin-combinations).",
        "\n"))
      # List the associated Lin-combinations.
      list_Lincombinations <- Lincomb_precomp_links %>%
        .[preconnected_component_ID == selected_precomponent_ID,
          Lincombination_ID]
      # Filter the associated Lin-combination matrix.
      Lincombination_subMatrix <- Lincombination_matrix %>%
        .[Lincombination_ID %in% list_Lincombinations]
      # Remove the first column of "Lincombination_ID".
      transaction_subMatrix <- Lincombination_subMatrix %>%
        .[, 2:ncol(.)]
      # Convert each Linclustering into a factor.
      transaction_subMatrix <- transaction_subMatrix %>%
        .[, lapply(.SD, as.factor), .SDcols = colnames(.)]
      # Convert the subMatrix into transactions.
      transaction_subMatrix <- transaction_subMatrix %>%
        arules::transactions()
      # Compute the dissimilarity.
      selected_similMatrix <- transaction_subMatrix %>%
        arules::dissimilarity(method = "dice")
      # Convert into dissimilarity into similarity.
      selected_similMatrix <- 1 - selected_similMatrix
      # Add the Lincombination_ID.
      attr(selected_similMatrix, "Labels") <- list_Lincombinations
      # Save the computed similarity matrix.
      selected_similMatrix %>%
        saveRDS(file = paste0(
          similarity_submatrices_directory,
          selected_precomponent_ID,
          ".RDS"))
    }
  })


##### List similarity submatrices #####


list_similarity_submatrices <- nb_Lincombs_per_preconnected_component[
  , .(preconnected_component_ID)]


##### Output #####


data.table::fwrite(
  x = list_similarity_submatrices,
  file = snakemake@output[["list_similarity_submatrices"]],
  sep = ";")

