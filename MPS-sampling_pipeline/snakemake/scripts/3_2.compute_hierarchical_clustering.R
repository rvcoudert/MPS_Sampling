# 3_2.compute_hierarchical_clustering.R
# Compute hierarchical tree.


# save.image("~/3_2.compute_hierarchical_clustering.RDS")
# stop()
# load("~/3_2.compute_hierarchical_clustering.RDS")


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

dissimilarity_matrix <- snakemake@input[["dissimilarity_matrix"]] %>%
  readRDS()


##### Compute #####


hierarchical_tree <- fastcluster::hclust(dissimilarity_matrix)


##### Output #####


saveRDS(object = hierarchical_tree,
        file = snakemake@output[["hierachical_clustering"]])
