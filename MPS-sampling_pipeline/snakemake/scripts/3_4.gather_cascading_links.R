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

# 3_4.gather_cascading_links.R
# Build MPS-clusters through hierarchical clustering with complete-linkage.


# save.image("~/3_4.gather_cascading_links.RDS")
# stop()
# load("~/3_4.gather_cascading_links.RDS")


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


gen_Lincomb_links <- snakemake@input[["gen_Lincomb_links"]] %>%
  data.table::fread()
Lincomb_MPS_cluster_links <- snakemake@input[["Lincomb_MPS_cluster_links"]] %>%
  data.table::fread()


##### Compute ######


# Simply merge all <data.tables>.
genome_cascading_links <- data.table::merge.data.table(
  x = gen_Lincomb_links,
  y = Lincomb_MPS_cluster_links,
  by = "Lincombination_ID",
  all = TRUE,
  order = TRUE)
# Re-order the columns.
genome_cascading_links <- genome_cascading_links[
  , .(genomeAccession,
      Lincombination_ID,
      MPS_cluster_ID)]


##### Output #####


data.table::fwrite(
  x = genome_cascading_links,
  file = snakemake@output[["genome_cascading_links"]],
  sep = ";")
