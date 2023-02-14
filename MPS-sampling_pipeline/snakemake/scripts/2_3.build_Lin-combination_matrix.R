# 2_3.build_Lin-combination_matrix.R
# Remove redundant Lin-combinations to build a Lin-combination matrix.


# save.image("~/2_3.build_Lin-combination_matrix.RDS")
# stop()
# load("~/2_3.build_Lin-combination_matrix.RDS")


##### Input #####


library(magrittr)
print(sessionInfo())

Linclustering_matrix <- snakemake@input[["Linclustering_matrix"]] %>%
  data.table::fread()


#### Lin-combination matrix #####


# List the names of the protein family columns.
proteinFamily_columns <- colnames(Linclustering_matrix[, -1])
# Compute number of genomes per Lin-combination.
Lincombination_stats <- Linclustering_matrix %>%
  .[, .(nb_genomes = .N),
    by = proteinFamily_columns]
# Select only the protein family columns.
Lincombination_matrix <- Lincombination_stats %>%
  .[, ..proteinFamily_columns]


##### Gen_Lincomb_links #####


# Compute the column <Lincombination_ID>
#   by counting according to the number of genomes of each Lin-combination.
Lincombination_ID_column <- Lincombination_stats[, nb_genomes] %>%
  magrittr::set_names(1:length(.)) %>%
  purrr::imap(function(the_nb_genomes, the_Lincombination_ID) {
    the_Lincombination_ID %>%
      rep(the_nb_genomes) %>%
      as.integer()
  }) %>% unname() %>% unlist()
gen_Lincomb_links <- cbind(
  Linclustering_matrix[, .(genomeAccession)],
  Lincombination_ID = Lincombination_ID_column)


##### Output #####


data.table::fwrite(
  x = Lincombination_matrix,
  file = snakemake@output[["Lincombination_matrix"]],
  sep = ";")
data.table::fwrite(
  x = gen_Lincomb_links,
  file = snakemake@output[["gen_Lincomb_links"]],
  sep = ";")

