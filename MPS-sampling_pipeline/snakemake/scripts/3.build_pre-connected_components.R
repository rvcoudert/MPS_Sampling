# 3.build_pre-connected_components.R
# Build pre-connected components.


# save.image("~/3.build_pre-connected_components.RDS")
# stop()
# load("~/3.build_pre-connected_components.RDS")


##### Input #####


library(magrittr)
print("sessionInfo()")
print(sessionInfo())
print(".libPaths()")
print(.libPaths())
print("loadedNamespeces()")
print(loadedNamespaces())
cat("\n")

Lincombination_matrix <- snakemake@input[["Lincombination_matrix"]] %>%
  data.table::fread()
min_nb_Linclusters <- as.integer(snakemake@params[["min_nb_Lincluster"]])


##### Functions #####


# Count common Lin-clusters
#
#
# From a Lin-combination matrix,
# count the number of common Lin-clusters with each other Lin-combination.
#
#
# The first column of the Lin-combination matrix
# needs to be the Lin-combination_ID.
count_common_Linclusters <- function(
  Lincombination_matrix,
  target_Lincombinations,
  ignore_target_Lincombinations = TRUE,
  ignored_Lincombinations = NULL) {

  ## () List the concerned protein families for the Lin-combinations.
  # List the protein family names of the Lin-combination matrix;
  #   exclude the first column which is supposed to be the Lin-combination_ID.
  theColnames <- colnames(Lincombination_matrix)
  column_ID <- theColnames[[1]]
  proteinFamilyColnames <- theColnames[2:length(theColnames)]
  names(proteinFamilyColnames) <- proteinFamilyColnames

  ## () Get the target Lin-clusters for each protein family.
  # Filter the specific Lin-combination_ID
  # according to the given Lin-combination_IDs.
  target_Lincombination_matrix <- Lincombination_matrix[
    get(column_ID) %in% target_Lincombinations]
  # If no Lin-combination is associated to the Lin-combination_ID,
  #   then return an empty list with a warning.
  if (nrow(target_Lincombination_matrix) == 0) {
    warning("Lin-combination_ID ",
            "\"",
            target_Lincombinations,
            "\"",
            " is missing in the first column called ",
            "\"",
            column_ID,
            "\"",
            ".")}
  # List the target Lin-clusters for each protein family,
  # according to the target Lin-combinations.
  target_Linclusters <- plyr::llply(
    proteinFamilyColnames,
    function(LincombinationColname) {
      # DEBUG #
      # LincombinationColname <- proteinFamilyColnames[[1]]

      # Get the target Lin-clusters.
      target_Lincluster <- target_Lincombination_matrix[
        , get(LincombinationColname)]
      target_Lincluster <- unique(target_Lincluster)
      target_Lincluster <- target_Lincluster[!is.na(target_Lincluster)]
      return(target_Lincluster)
    })

  ## () Filter the remaining Lin-combinations.
  # If asked, then remove the target Lin-combinations
  #   and the ignored Lin-combinations.
  # Else, only remove the ignored Lin-combinations.
  if (ignore_target_Lincombinations) {
    remaining_Lincombination_matrix <- Lincombination_matrix[
      !get(column_ID) %in% c(target_Lincombinations,
                             ignored_Lincombinations)]
  } else {
    remaining_Lincombination_matrix <- Lincombination_matrix[!get(column_ID) %in%
                                                               ignored_Lincombinations]
  }


  ## () Vertically count the number of common Lin-clusters.
  # For each protein family, list the common Lin-clusters.
  list_common_Linclusters <- plyr::llply(
    proteinFamilyColnames,
    function(LincombinationColname) {
      # DEBUG #
      # LincombinationColname <- proteinFamilyColnames[[1]]

      # Check which Lin-combinations shared the target Lin-clusters.
      concerned_Lincombinations <- remaining_Lincombination_matrix[
        get(LincombinationColname) %in%
          target_Linclusters[[LincombinationColname]],
        1]
      concerned_Lincombinations <- unlist(concerned_Lincombinations, use.names = FALSE)
      return(concerned_Lincombinations)
    })
  # Flatten the list of column lists of common Lin-clusters.
  list_common_Linclusters <- unlist(list_common_Linclusters, use.names = FALSE)
  # Then count the number of common Lin-clusters per Lin-combination.
  list_common_Linclusters <- table(list_common_Linclusters)
  # Format the results.
  list_common_Linclusters <- as.list(list_common_Linclusters)
  list_common_Linclusters <- unlist(list_common_Linclusters)
  # Export the expected results.
  return(list_common_Linclusters)
}


# Get isolated Lin-combinations
#
#
# From a Lin-combination matrix, list the isolated Lin-combinations.
#
# An isolated Lin-combination is a Lin-combination
# which does not contain enough information
# to be connected to another one,
# depending on a minimum number of shared Lin-clusters.
#
# It happens if the number of connectable Lin-clusters
# (not Lin-singletons and not NA) is less than the minimum of shared Lin-clusters.
get_isolated_Lincombinations <- function(
  Lincombination_matrix,
  min_nb_Linclusters = ceiling(ncol(Lincombination_matrix) / 2),
  verbose = FALSE) {
  # Remove the first column, which supposed to be the "Lin-combination_ID".
  if (colnames(Lincombination_matrix)[[1]] == "Lincombination_ID") {
    Lincombination_matrix <- Lincombination_matrix %>%
      .[, 2:ncol(.)]
  } else {
    warning("Lin-combination_ID not found as first column.",
            "\nNo column has been removed.")
  }

  # Count missing values per Lin-combination.
  count_NA <- Lincombination_matrix %>%
    .[, (nb_obs = rowSums(is.na(Lincombination_matrix)))]

  # For each protein family, put a 0 if it's a Lin-singleton.
  count_Linsingletons_DT <- Lincombination_matrix %>%
    data.table::copy() %>%
    .[
      ,
      lapply(.SD,
             function(Linclustering) {
               # DEBUG #
               # Linclustering <- Lincombination_matrix[, 2]
               # List all the Lin-singletons of
               #  the current protein family Lin-clustering.
               list_Linsingletons <- Linclustering %>%
                 table() %>%
                 .[. == 1] %>%
                 names() %>%
                 as.integer() %>%
                 .[!is.na(.)]
               # Put a 0 for these Lin-singletons.
               Linclustering[Linclustering %in% list_Linsingletons] <- 0
               # Return the corrected protein clustering.
               return(Linclustering)
             }),
      .SDcols = colnames(.)]

  # Count number of Lin-singletons per row.
  count_Linsingletons_perRow <- count_Linsingletons_DT %>%
    .[, (num_obs = rowSums(count_Linsingletons_DT == 0,
                           na.rm = TRUE))]

  # Remove Lin-combinations according to the threshold.
  nb_cols <- ncol(Lincombination_matrix)

  # List the isolated Lin-combinations.
  isolated_Lincombinations <- which((count_NA + count_Linsingletons_perRow) >
                                      (nb_cols - min_nb_Linclusters))

  # If asked, then print a message
  # about the isolated Lin-combinations that have been found.
  if (verbose) {
    cat(isolated_Lincombinations %>% length() %>% format(big.mark = " "),
        "over the",
        Lincombination_matrix %>% nrow() %>% format(big.mark = " "),
        "Lin-combinations have more than",
        (nb_cols - min_nb_Linclusters) %>% format(big.mark = " "),
        "unusable values\nover the",
        Lincombination_matrix %>% ncol() %>% format(big.mark = " "),
        "protein families (missing values or Lin-singletons)",
        "\nand don't have enough information",
        "to be linked to another Lin-combination.",
        "\n")
  }

  # Export the expected results.
  return(isolated_Lincombinations)
}


# Get pre-connected Lin-combinations
#
#
# From a Lin-combination matrix and the ID of some target Lin-combinations,
# filter the pre-connected Lin-combinations
# depending on a minimum number of shared Lin-clusters.
#
#
# The first column of the Lin-combination matrix
# needs to be the Lin-combination_ID.
get_preconnected_Lincombinations <- function(
  Lincombination_matrix,
  target_Lincombinations,
  min_nb_Linclusters = ceiling(ncol(Lincombination_matrix) / 2),
  ignored_Lincombinations = NULL) {

  # Count common Lin-clusters between the target and the rest.
  list_common_Linclusters <- count_common_Linclusters(
    Lincombination_matrix = Lincombination_matrix,
    target_Lincombinations = target_Lincombinations,
    ignored_Lincombinations = ignored_Lincombinations)
  # Keep only the Lin-combinations above the given threshold of common elements.
  preconnected_Lincombinations <- names(list_common_Linclusters[
    list_common_Linclusters >= min_nb_Linclusters])
  # Convert the index of the pre-connected Lin-combinations as an integer.
  preconnected_Lincombinations <- as.integer(preconnected_Lincombinations)
  # Add the target Lin-combinations.
  preconnected_Lincombinations <- c(target_Lincombinations,
                                    preconnected_Lincombinations)
  # Export the expected results.
  return(preconnected_Lincombinations)
}


# Get pre-connected component
#
#
# From a Lin-combination matrix and the ID of some target Lin-combinations,
# get the pre-connected component of theses Lin-combinations
# depending on a minimum number of shared Lin-clusters.
#
#
# The first column of the Lin-combination matrix
# needs to be the Lin-combination_ID.
get_preconnected_component <- function(
  Lincombination_matrix,
  target_Lincombinations,
  min_nb_Linclusters = ceiling(ncol(Lincombination_matrix) / 2),
  ignored_Lincombinations = NULL) {

  # Initialize the old pre-connected Lin-combinations.
  old_preconnected_Lincombinations <- target_Lincombinations

  # Initialize the new pre-connected Lin-combinations.
  new_preconnected_Lincombinations <- get_preconnected_Lincombinations(
    Lincombination_matrix = Lincombination_matrix,
    target_Lincombinations = old_preconnected_Lincombinations,
    min_nb_Linclusters = min_nb_Linclusters,
    ignored_Lincombinations = ignored_Lincombinations)

  # Count the number of pre-connected Lin-combinations.
  n_old <- length(old_preconnected_Lincombinations)
  n_new <- length(new_preconnected_Lincombinations)

  # If new pre-connected Lin-combinations have been found,
  #   then start a new cycle.
  # If no new pre-connected Lin-combination has been found,
  #   then end the process.
  while (n_old != n_new) {
    # Update the old pre-connected Lin-combinations,
    #   ie the new ones become the old ones.
    old_preconnected_Lincombinations <- new_preconnected_Lincombinations
    # Get the new pre-connected Lin-combinations.
    new_preconnected_Lincombinations <- get_preconnected_Lincombinations(
      target_Lincombinations = old_preconnected_Lincombinations,
      Lincombination_matrix = Lincombination_matrix,
      min_nb_Linclusters = min_nb_Linclusters,
      ignored_Lincombinations = ignored_Lincombinations)
    # Count the number of pre-connected Lin-combinations.
    n_old <- length(old_preconnected_Lincombinations)
    n_new <- length(new_preconnected_Lincombinations)
  }

  # Export the results.
  return(old_preconnected_Lincombinations)
}


# Split a Lin-combination matrix into pre-connected components
#
#
# From a Lin-combination matrix,
# split it into pre-connected components.
#
# Begin with the first Lin-combination, get its pre-connected component,
# then do it again with the next first remaining Lin-combination,
# iterate the process until all the Lin-combinations
# have their own pre-connected component.
#
#
# The first column of the Lin-combination matrix
# needs to be the Lin-combination_ID.
split_into_preconnected_components <- function(
  Lincombination_matrix,
  min_nb_Linclusters = ceiling(ncol(Lincombination_matrix) / 2),
  ignored_Lincombinations = NULL,
  long_format = TRUE,
  verbose = FALSE) {

  ## () Pre-process.
  # Get the name of the row ID (first column).
  theColnames <- colnames(Lincombination_matrix)
  column_ID <- theColnames[[1]]

  ## () Initialize the loop.
  # Initialize the list of pre-connected components.
  list_component <- list()
  # List the isolated Lin-combinations.
  isolated_Lincombinations <- get_isolated_Lincombinations(
    Lincombination_matrix = Lincombination_matrix,
    min_nb_Linclusters = min_nb_Linclusters,
    verbose = verbose)

  # Initialize the remaining Lin-combination matrix to split.
  remaining_Lincombination_matrix <- Lincombination_matrix[
    !get(column_ID) %in% ignored_Lincombinations]
  remaining_Lincombination_matrix <- Lincombination_matrix[
    !get(column_ID) %in% isolated_Lincombinations]

  # If asked, then print a progression message.
  if (verbose) {
    cat("\n")
    cat("Remaining Lin-combinations  : ")
    cat(format(x = nrow(remaining_Lincombination_matrix), big.mark = " "))
    cat("\n")
  }

  ## () Start the loop.
  # While at least one Lin-combination is remaining.
  while (nrow(remaining_Lincombination_matrix) > 0) {
    # Get the first remaining Lin-combination.
    first_remaining_Lincombination <- remaining_Lincombination_matrix[
      1, get(column_ID)]

    # Get the next component.
    current_component <- get_preconnected_component(
      Lincombination_matrix = remaining_Lincombination_matrix,
      target_Lincombinations = first_remaining_Lincombination,
      min_nb_Linclusters = min_nb_Linclusters)

    # Add this component to the list.
    list_component <- append(list_component, list(current_component))

    # Update the remaining Data Table.
    remaining_Lincombination_matrix <- remaining_Lincombination_matrix[
      !get(column_ID) %in% current_component]

    # If asked, then print a progression message.
    if (verbose) {
      cat("Remaining Lin-combinations  : ")
      cat(format(x = nrow(remaining_Lincombination_matrix), big.mark = " "))
      cat("\n")
    }
  }

  # Export the results
  return(list_component)
}



#### Lincomb_precomp_links #####


# Add a first column with the Lin-combination ID.
Lincombination_matrix <- cbind(
  "Lincombination_ID" = 1:nrow(Lincombination_matrix),
  Lincombination_matrix)

# List the pre-connected components.
preconnected_components <- split_into_preconnected_components(
    Lincombination_matrix = Lincombination_matrix,
    min_nb_Linclusters = min_nb_Linclusters,
    verbose = TRUE)

# Convert into long format.
Lincomb_precomp_links <- preconnected_components %>%
  magrittr::set_names(1:length(.)) %>%
  plyr::ldply(function(list_Lincombinations) {
    data.table::data.table(Lincombination_ID = list_Lincombinations)
  }, .progress = "text", .id = "preconnected_component_ID")
Lincomb_precomp_links <- data.table::as.data.table(Lincomb_precomp_links)

# Order the columns.
Lincomb_precomp_links <- Lincomb_precomp_links[, .(Lincombination_ID,
                                             preconnected_component_ID)]
# Add missing Lin-combinations.
Lincomb_precomp_links <- data.table::merge.data.table(
  x = Lincombination_matrix[, .(Lincombination_ID)],
  y = Lincomb_precomp_links,
  by = "Lincombination_ID",
  all.x = TRUE,
  sort = TRUE)
# Order the rows.
Lincomb_precomp_links <- Lincomb_precomp_links[order(Lincombination_ID)]


##### Output #####


data.table::fwrite(
  x = Lincomb_precomp_links,
  file = snakemake@output[["Lincomb_precomp_links"]],
  sep = ";")

