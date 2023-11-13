# 3_1b.convert_dissimilarity_matrix.R
# Convert dissimilarity matrix.


# save.image("~/3_1b.convert_dissimilarity_matrix.RDS")
# stop()
# load("~/3_1b.convert_dissimilarity_matrix.RDS")


##### Input #####


library(magrittr)
print("sessionInfo()")
print(sessionInfo())
print(".libPaths()")
print(.libPaths())
print("loadedNamespeces()")
print(loadedNamespaces())
cat("\n")


# cores <- snakemake@params[["cores"]]
# cat("Parallelization on", cores, "cores.\n")
#
#
# if (cores > 1) {
#   doParallel::registerDoParallel(cores = cores)
#   parallel <- TRUE
# } else {
#   parallel <- FALSE
# }

dissimilarity_matrix <- snakemake@input[["dissimilarity_matrix_raw"]] %>%
  readRDS()


##### Compute #####


dissimilarity_matrix <- as.dist(dissimilarity_matrix)


##### Output #####


saveRDS(object = dissimilarity_matrix,
        file = snakemake@output[["dissimilarity_matrix"]])


print(Sys.time())
