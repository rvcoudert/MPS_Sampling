# 3_3.build_MPS-clusters.R
# Build MPS-clusters through hierarchical clustering with complete-linkage.


# save.image("~/3_3.build_MPS-clusters.RDS")
# stop()
# load("~/3_3.build_MPS-clusters.RDS")


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


hierachical_clustering <- snakemake@input[["hierachical_clustering"]] %>%
  readRDS()
delta <- snakemake@params[["delta"]] %>%
  as.numeric()


##### Compute ######


# Compute hierarchical clustering for the specified delta.
my_cutree <- cutree(tree = hierachical_clustering,
                    h = 1 - delta)
# Tidy the hierarchical clustering.
Lincomb_MPS_cluster_links <- data.table::data.table(
  Lincombination_ID = 1:length(my_cutree),
  MPS_cluster_ID = my_cutree)


##### Output #####


data.table::fwrite(
  x = Lincomb_MPS_cluster_links,
  file = snakemake@output[["Lincomb_MPS_cluster_links"]],
  sep = ";")

