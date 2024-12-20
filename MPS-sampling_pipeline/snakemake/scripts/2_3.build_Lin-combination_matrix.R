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

# 2_3.build_Lin-combination_matrix.R
# Remove redundant Lin-combinations to build a Lin-combination matrix.


# save.image("~/2_3.build_Lin-combination_matrix.RDS")
# stop()
# load("~/2_3.build_Lin-combination_matrix.RDS")


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

